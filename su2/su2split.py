"""Split an internal surface into two coincident faces for SU2's actuator disk.

SU2's VARIABLE_LOAD actuator disk needs the inlet and outlet as two distinct
markers whose nodes are geometrically coincident. Two constraints make this
awkward to obtain from gmsh:

  * gmsh's ``fragment()`` always produces a conformal mesh, so the shared
    interface is a single merged face - ``occ.copy()`` of the face gets merged
    as well;
  * SU2 maps every node to a radial station through
    ``r_ = |P - CENTER| / RADIUS`` (full 3-D distance, not the perpendicular
    distance from the axis), so a physical gap between the faces corrupts the
    radial distribution and can push rim nodes past the last station.

This module solves it in the mesh instead: nodes on the interface are
duplicated, every volume element on the +x side is re-indexed onto the
duplicates, and the surface is written twice - once with the original nodes
(``tag_in``) and once with the duplicates (``tag_out``).
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from su2mesh import GMSH2SU2, SU2_REORDER   # noqa: E402


def write_su2_mesh_doubled(gmsh_mod, path, split_tag, x_plane=0.0,
                           tag_in='disk_in', tag_out='disk_out',
                           exclude=('disk',)):
    """Write the gmsh model with ``split_tag`` doubled into two markers.

    Returns (npoin, nelem, nmark, n_split_tri, n_split_nodes).
    """
    g = gmsh_mod
    node_tags, node_coords, _ = g.model.mesh.getNodes()
    coord = {int(t): list(node_coords[3 * i:3 * i + 3])
             for i, t in enumerate(node_tags)}

    # ---- the surface to split (triangles) -------------------------------
    split_tris = []
    for et, conn in zip(*g.model.mesh.getElements(2, split_tag)[::2]):
        _, _, _, nn, _, _ = g.model.mesh.getElementProperties(et)
        if nn != 3:
            continue
        for i in range(0, len(conn), nn):
            split_tris.append([int(x) for x in conn[i:i + nn]])
    split_nodes = {n for tri in split_tris for n in tri}

    # ---- volume elements, tagged by their side of the plane -------------
    bulk = []
    for dim, pt in g.model.getPhysicalGroups():
        if dim != 3:
            continue
        for ent in g.model.getEntitiesForPhysicalGroup(dim, pt):
            ets, _, enodes = g.model.mesh.getElements(dim, ent)
            for et, conn in zip(ets, enodes):
                _, _, _, nn, _, _ = g.model.mesh.getElementProperties(et)
                et2 = GMSH2SU2.get(int(et))
                if et2 is None:
                    continue
                for i in range(0, len(conn), nn):
                    ids = [int(x) for x in conn[i:i + nn]]
                    cx = sum(coord[j][0] for j in ids) / nn
                    bulk.append((et2, ids, cx > x_plane))

    # ---- duplicate the interface nodes ----------------------------------
    dup = {}
    next_tag = max(coord) + 1
    for n in sorted(split_nodes):
        dup[n] = next_tag
        coord[next_tag] = list(coord[n])
        next_tag += 1

    remap_bulk = []
    for et2, ids, plus_side in bulk:
        if plus_side:
            ids = [dup.get(j, j) for j in ids]
        re_order = SU2_REORDER.get(et2)
        if re_order:
            ids = [ids[j] for j in re_order]
        remap_bulk.append((et2, ids))

    # ---- boundary markers (the split surface is handled separately) -----
    markers = []
    for dim, pt in g.model.getPhysicalGroups():
        if dim != 2:
            continue
        name = g.model.getPhysicalName(dim, pt)
        if name in exclude:
            continue
        conns, etype = [], None
        for ent in g.model.getEntitiesForPhysicalGroup(dim, pt):
            ets, _, enodes = g.model.mesh.getElements(dim, ent)
            for et, conn in zip(ets, enodes):
                _, _, _, nn, _, _ = g.model.mesh.getElementProperties(et)
                etype = GMSH2SU2.get(int(et))
                if etype is None:
                    continue
                for i in range(0, len(conn), nn):
                    conns.append([int(x) for x in conn[i:i + nn]])
        if conns:
            markers.append((name, etype, conns))
    # the two copies must face OPPOSITE ways: the interface winding follows
    # the -x side volume, so the +x copy is reversed. Emitting both with the
    # same winding makes SU2 treat them as co-oriented and it cannot tell
    # which one is upstream of the disk.
    markers.append((tag_in, 5, [list(t) for t in split_tris]))
    markers.append((tag_out, 5,
                    [[dup[n] for n in (t[0], t[2], t[1])]
                     for t in split_tris]))

    # ---- keep only referenced nodes -------------------------------------
    used = set()
    for _, ids in remap_bulk:
        used.update(ids)
    for _, _, conns in markers:
        for ids in conns:
            used.update(ids)
    order = sorted(used)
    idx = {old: new for new, old in enumerate(order)}

    with open(path, 'w') as fh:
        fh.write('NDIME= 3\n')
        fh.write(f'NELEM= {len(remap_bulk)}\n')
        for i, (et2, ids) in enumerate(remap_bulk):
            fh.write(f'{et2} ' + ' '.join(str(idx[j]) for j in ids)
                     + f' {i}\n')
        fh.write(f'NPOIN= {len(order)}\n')
        for new, old in enumerate(order):
            c = coord[old]
            fh.write(' '.join(f'{v:.12e}' for v in c) + f' {new}\n')
        fh.write(f'NMARK= {len(markers)}\n')
        for name, et2, conns in markers:
            fh.write(f'MARKER_TAG= {name}\n')
            fh.write(f'MARKER_ELEMS= {len(conns)}\n')
            for ids in conns:
                fh.write(f'{et2} ' + ' '.join(str(idx[j]) for j in ids) + '\n')
    return (len(order), len(remap_bulk), len(markers), len(split_tris),
            len(split_nodes))
