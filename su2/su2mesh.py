"""gmsh -> SU2 native ASCII mesh writer (avoids the CGNS ADF/HDF5 mismatch)."""
# gmsh element type -> SU2 element type. The 3-D codes matter for
# boundary-layer meshes: the 2-D BoundaryLayer field emits quads, and
# extruding those gives PRISMS (gmsh type 6), not tetrahedra. Omitting type 6
# silently drops every near-wall cell.
GMSH2SU2 = {1: 3, 2: 5, 3: 9, 4: 10, 5: 12, 6: 13, 7: 14}
# SU2: 3=line, 5=tri, 9=quad, 10=tet, 12=hex, 13=prism, 14=pyramid
# gmsh and SU2 disagree on the prism node order: gmsh lists the two triangles
# as (0,1,2,3,4,5), SU2 wants the top triangle first, (3,4,5,0,1,2).
SU2_REORDER = {}


def write_cgns_mesh(gmsh_mod, path):
    """Write CGNS (HDF5) via meshio, so gmsh owns the connectivity.

    SU2 8.x reads HDF5 CGNS natively. Routing the mesher output through meshio
    avoids hand-writing element connectivity entirely - the prism ordering and
    the quad-face association are produced by gmsh and passed through
    unchanged, which is what SU2's face matching needs.
    """
    import meshio
    g = gmsh_mod
    node_tags, node_coords, _ = g.model.mesh.getNodes()
    tag2idx = {int(t): i for i, t in enumerate(node_tags)}
    points = node_coords.reshape(-1, 3)

    cells, cell_data = [], {}
    for dim in (3, 2, 1):
        for ptag in [p for d, p in g.model.getPhysicalGroups() if d == dim]:
            name = g.model.getPhysicalName(dim, ptag)
            for dd, ent in g.model.getEntities(dim):
                if (dd, ent) not in {(dim, e) for e in
                                     g.model.getEntitiesForPhysicalGroup(dim, ptag)}:
                    continue
                etypes, _, enodes = g.model.mesh.getElements(dd, ent)
                for et, conn in zip(etypes, enodes):
                    _, _, _, nn, _, _ = g.model.mesh.getElementProperties(et)
                    blocks = [conn[i:i + nn] for i in range(0, len(conn), nn)]
                    arr = [[tag2idx[int(x)] for x in b] for b in blocks]
                    if not arr:
                        continue
                    # gmsh -> meshio type names
                    typ = {1: 'line', 2: 'triangle', 3: 'quad',
                           4: 'tetra', 5: 'hexahedron', 6: 'wedge',
                           7: 'pyramid'}[int(et)]
                    cells.append((typ, arr))
                    cell_data.setdefault('gmsh:physical', []).append(
                        [ptag] * len(arr))
                    cell_data.setdefault('gmsh:geometrical', []).append(
                        [ent] * len(arr))
    m = meshio.Mesh(points, cells,
                    cell_data={'gmsh:physical': cell_data['gmsh:physical'],
                               'gmsh:geometrical':
                                   cell_data['gmsh:geometrical']})
    meshio.write(path, m, file_format='cgns')
    return len(points), sum(len(c[1]) for c in cells)


def write_msh_mesh(gmsh_mod, path):
    """Write a gmsh MSH 2.2 file instead of hand-rolling SU2 connectivity.

    SU2 reads MSH 2.2 natively (format detected from the extension), so the
    prism connectivity lives in gmsh's own numbering. This sidesteps the
    quad-face ordering problem that made hand-written prism meshes fail with
    "surface element does not have an associated volume element".

    Physical tags are required for the boundary elements; the first tag of an
    element is its physical group id, the second its geometric entity.
    """
    g = gmsh_mod
    node_tags, node_coords, _ = g.model.mesh.getNodes()
    tag2idx = {int(t): i for i, t in enumerate(node_tags)}

    phys_of = {}
    for dim, ptag in g.model.getPhysicalGroups():
        for ent in g.model.getEntitiesForPhysicalGroup(dim, ptag):
            phys_of[(dim, ent)] = ptag

    elems = []          # (gmsh type, phys tag, geom entity, [node tags])
    for dim in (1, 2, 3):
        for dd, ent in g.model.getEntities(dim):
            ptag = phys_of.get((dd, ent))
            if ptag is None:
                continue
            etypes, _, enodes = g.model.mesh.getElements(dd, ent)
            for et, conn in zip(etypes, enodes):
                _, _, _, nn, _, _ = g.model.mesh.getElementProperties(et)
                for i in range(0, len(conn), nn):
                    elems.append((int(et), ptag, ent,
                                  [int(x) for x in conn[i:i + nn]]))

    with open(path, 'w') as fh:
        fh.write('$MeshFormat\n2.2 0 8\n$EndMeshFormat\n')
        fh.write(f'$Nodes\n{len(node_tags)}\n')
        for t in node_tags:
            i = tag2idx[int(t)]
            x, y, z = node_coords[3 * i:3 * i + 3]
            fh.write(f'{int(t)} {x:.12e} {y:.12e} {z:.12e}\n')
        fh.write('$EndNodes\n')
        fh.write(f'$Elements\n{len(elems)}\n')
        for n, (et, ptag, ent, ids) in enumerate(elems, start=1):
            fh.write(f'{n} {et} 2 {ptag} {ent} ' +
                     ' '.join(str(x) for x in ids) + '\n')
        fh.write('$EndElements\n')
    return len(node_tags), len(elems)


def write_su2_mesh(gmsh_mod, path, to_su2_frame=False):
    """Export the current gmsh model to a .su2 ASCII mesh.

    Physical groups must already be defined; their names become SU2 markers.
    The highest-dimensional group(s) form the bulk elements.

    to_su2_frame: rotate the nodes from the wing frame (x chord, y span,
    z up) into SU2's frame (x chord, y up, z span) with the proper rotation
    (x, y, z) -> (x, -z, y). SU2 rotates the freestream by AOA in its x-y
    plane, so the span must lie along z and lift along y (ONERA M6
    convention). A plain y/z exchange would be an improper transform
    (determinant -1) and invert every element; the inverse rotation
    (x, z, -y) is proper but turns the wing upside down, reversing the lift.
    """
    g = gmsh_mod
    node_tags, node_coords, _ = g.model.mesh.getNodes()
    tag2idx = {int(t): i for i, t in enumerate(node_tags)}
    coord_of = {int(t): node_coords[3 * i:3 * i + 3]
                for i, t in enumerate(node_tags)}
    ndim = int(g.model.getDimension())

    pg = g.model.getPhysicalGroups()   # [(dim, tag), ...]
    all_dims = [d for d, _ in pg]
    bulk_dims = [d for d in set(all_dims) if d == ndim]

    bulk = []
    seen_nodes = None
    for d in bulk_dims:
        for (dd, pt) in pg:
            if dd != d:
                continue
            for ent in g.model.getEntitiesForPhysicalGroup(dd, pt):
                etypes, _, enodes = g.model.mesh.getElements(dd, ent)
                for et, conn in zip(etypes, enodes):
                    _, _, _, nn, _, _ = g.model.mesh.getElementProperties(et)
                    et2 = GMSH2SU2.get(int(et))
                    if et2 is None:
                        continue
                    for i in range(0, len(conn), nn):
                        ids = [tag2idx[int(x)] for x in conn[i:i + nn]]
                        re_order = SU2_REORDER.get(et2)
                        if re_order:
                            ids = [ids[j] for j in re_order]
                        bulk.append((et2, ids))

    markers = []
    for (dd, pt) in pg:
        if dd >= ndim or dd != ndim - 1:
            continue          # SU2 markers must be (ndim-1)-D faces
        name = g.model.getPhysicalName(dd, pt)
        conns, etype = [], None
        for ent in g.model.getEntitiesForPhysicalGroup(dd, pt):
            etypes, _, enodes = g.model.mesh.getElements(dd, ent)
            for et, conn in zip(etypes, enodes):
                _, _, _, nn, _, _ = g.model.mesh.getElementProperties(et)
                etype = GMSH2SU2.get(int(et))
                if etype is None:
                    continue
                for i in range(0, len(conn), nn):
                    conns.append([tag2idx[int(x)] for x in conn[i:i + nn]])
        if conns:
            markers.append((name, etype, conns))

    # keep only nodes referenced by any element (SU2 rejects orphans)
    used = set()
    for _, nodes in bulk:
        used.update(nodes)
    for _, _, conns in markers:
        for nodes in conns:
            used.update(nodes)
    order = sorted(used)
    remap = {old: new for new, old in enumerate(order)}
    npoin = len(order)

    with open(path, 'w') as fh:
        fh.write(f'NDIME= {ndim}\n')
        fh.write(f'NELEM= {len(bulk)}\n')
        for i, (et, nodes) in enumerate(bulk):
            fh.write(f'{et} ' + ' '.join(str(remap[x]) for x in nodes) + f' {i}\n')
        fh.write(f'NPOIN= {npoin}\n')
        tags_by_idx = {v: k for k, v in tag2idx.items()}
        for new, old in enumerate(order):
            c = coord_of[tags_by_idx[old]]
            if to_su2_frame:
                c = (c[0], -c[2], c[1])
            fh.write(' '.join(f'{v:.12e}' for v in c[:ndim]) + f' {new}\n')
        fh.write(f'NMARK= {len(markers)}\n')
        for name, et, conns in markers:
            fh.write(f'MARKER_TAG= {name}\n')
            fh.write(f'MARKER_ELEMS= {len(conns)}\n')
            for nodes in conns:
                fh.write(f'{et} ' + ' '.join(str(remap[x]) for x in nodes) + '\n')
    return npoin, len(bulk), len(markers)
