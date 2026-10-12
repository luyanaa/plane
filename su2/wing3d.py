"""Layer 12: 3-D half-wing mesh in SU2's frame, built by extrusion.

SU2 rotates the freestream by AOA inside its own x-y plane and reports
CL = CFz / CD = CFx, so the model is built as

    x  chordwise (flow direction)
    y  spanwise   (symmetry plane at y = 0)
    z  vertical   (lift direction)

The section is meshed as a 2-D O-grid and extruded along y, which is far more
robust than boolean solids. The deflected flap makes the contour C0-kinked at
the hinge, so the outline is split into separate curve segments there; a
single spline through a kink stalls the 2-D mesher.

Used for the layer-12 prop-wing-flap closure (k_pw). The wing is modelled at
constant chord equal to the MAC: k_pw is a ratio of lift increments normalised
by one reference area, so the taper idealisation cancels out.
"""
from math import cos, pi, radians, sin
from pathlib import Path

import numpy as np

MAC = 0.6788            # D mean aerodynamic chord (layer 1)
FLAP_FRAC = 0.25        # flap chord fraction (layer 3)
Y_HALF = 0.5 * 6.317    # D half span
FLAP_SPAN = 0.48 * Y_HALF


def read_section(dat, n=110, te_cut=0.995):
    raw = []
    for ln in Path(dat).read_text().splitlines():
        ln = ln.strip()
        if not ln or ln[0] in '#%':
            continue
        try:
            raw.append([float(v) for v in ln.split()])
        except ValueError:
            continue
    arr = np.array(raw)
    x, t = arr[:, 0], arr[:, 1]
    i_le = int(np.argmin(x))
    xu, tu = x[:i_le + 1][::-1], t[:i_le + 1][::-1]
    xl, tl = x[i_le:], t[i_le:]
    s2 = 0.5 * (1.0 - np.cos(np.linspace(0.0, pi, n)))

    def rs(xs, ts):
        s = np.linspace(0.0, 1.0, len(xs))
        m = xs <= te_cut
        return np.interp(s2, s[m], xs[m]), np.interp(s2, s[m], ts[m])

    a, b = rs(xu, tu)
    c, d = rs(xl, tl)
    return (a, b), (c, d)


def deflected_section(dat, chord, delta_deg, flap_frac=FLAP_FRAC):
    """Upper/lower contours in metres with the aft flap rotated about the
    hinge (TE-down positive)."""
    (xu, zu), (xl, zl) = read_section(dat)
    xh = 1.0 - flap_frac
    # TE-down is a positive deflection; flipping the section's thickness
    # convention here (thickness along -z) keeps the trailing edge moving
    # downward for delta > 0
    d = -radians(delta_deg)

    def rot(x, t):
        x = np.asarray(x, dtype=float)
        t = np.asarray(t, dtype=float)
        dx = x - xh
        aft = x >= xh
        return (np.where(aft, xh + dx * cos(d) - t * sin(d), x) * chord,
                np.where(aft, dx * sin(d) + t * cos(d), t) * chord)

    xu2, zu2 = rot(xu, zu)
    xl2, zl2 = rot(xl, zl)
    return (xu2, zu2), (xl2, zl2)


def _dedupe(pts):
    out = []
    for x, t in pts:
        if out and abs(x - out[-1][0]) < 1e-10 and abs(t - out[-1][1]) < 1e-10:
            continue
        out.append((float(x), float(t)))
    return out


def _contour_curves(geo, dat, chord, delta, flap_frac=FLAP_FRAC):
    """Split the closed section outline at the hinge and build each piece as
    its own OCC curve. Returns the ordered curve tags forming a closed loop.

    The flap hinge puts a C0 kink in the contour; meshing a single spline
    through it stalls the 2-D mesher, so the closed sequence is cut into the
    two arcs between the upper and lower hinge points.
    """
    (xu, zu), (xl, zl) = deflected_section(dat, chord, delta, flap_frac)
    xh = (1.0 - flap_frac) * chord
    _ = xh  # hinge station, kept for reference
    seq = _dedupe(list(zip(xu[::-1], zu[::-1])) +
                  list(zip(xl[1:], zl[1:])))
    if abs(seq[-1][0] - seq[0][0]) < 1e-12 and abs(seq[-1][1] - seq[0][1]) < 1e-12:
        seq = seq[:-1]

    # the sequence runs TE -> upper -> LE -> lower -> TE; split each surface at
    # the point nearest the hinge so the C0 kink is a curve endpoint
    npts = len(seq)
    i_le = int(np.argmin([x for x, _ in seq]))
    def nearest(lo, hi):
        sub = range(lo, hi)
        return min(sub, key=lambda i: abs(seq[i][0] - xh))
    i_up = nearest(0, i_le + 1) if i_le > 1 else 0
    i_lo = nearest(i_le, npts) if npts - i_le > 1 else npts - 1
    tags = [geo.addPoint(float(x), 0.0, float(t)) for x, t in seq]
    if i_up <= 0 or i_lo <= i_le or i_lo >= npts - 1:
        segs = [tags + [tags[0]]]
    else:
        segs = [tags[i_up:i_lo + 1], tags[i_lo:] + tags[:i_up + 1]]
    curves = []
    for seg in segs:
        if len(seg) >= 2:
            curves.append(geo.addBSpline(seg))
    return curves


def _wire_from_curves(geo, dat, chord, delta, y, want_tags=None,
                      flap_frac=FLAP_FRAC):
    """Build one closed section wire at span station y (2-D curves in the
    x-z plane). Mirrors _contour_curves but lifts the points to y."""
    (xu, zu), (xl, zl) = deflected_section(dat, chord, delta, flap_frac)
    xh = (1.0 - flap_frac) * chord
    seq = _dedupe(list(zip(xu[::-1], zu[::-1])) +
                  list(zip(xl[1:], zl[1:])))
    if abs(seq[-1][0] - seq[0][0]) < 1e-12 and abs(seq[-1][1] - seq[0][1]) < 1e-12:
        seq = seq[:-1]

    npts = len(seq)
    i_le = int(np.argmin([x for x, _ in seq]))

    def nearest(lo, hi):
        return min(range(lo, hi), key=lambda i: abs(seq[i][0] - xh))

    i_up = nearest(0, i_le + 1) if i_le > 1 else 0
    i_lo = nearest(i_le, npts) if npts - i_le > 1 else npts - 1
    tags = [geo.addPoint(float(x), float(y), float(t)) for x, t in seq]
    if i_up <= 0 or i_lo <= i_le or i_lo >= npts - 1:
        segs = [tags + [tags[0]]]
    else:
        segs = [tags[i_up:i_lo + 1], tags[i_lo:] + tags[:i_up + 1]]
    cur = [geo.addBSpline(seg) for seg in segs if len(seg) >= 2]
    loop = geo.addCurveLoop(cur)
    if want_tags is not None:
        want_tags.extend(cur)
    return loop


def build(out_su2, dat, delta_deg=0.0, span=Y_HALF, far=8.0, n_span=24,
          wall_size=0.012, far_size=1.2, chord=MAC, disk=None, verbose=False,
          bl_layers=0, bl_h0=1.5e-5, bl_growth=1.2, bl_thick=0.05):
    """Build and write the extruded half-wing mesh. Returns marker info.

    disk: optional dict(x=station, y=engine span station, z=axis height,
    radius, gap) - two coincident circular faces are embedded in the volume
    as the SU2 VARIABLE_LOAD inlet/outlet pair.
    """
    import gmsh
    gmsh.initialize()
    gmsh.option.setNumber('General.Terminal', 1 if verbose else 0)
    gmsh.model.add('halfwing')
    g = gmsh.model
    geo = g.occ

    curves = _contour_curves(geo, dat, chord, delta_deg)
    loop_af = geo.addCurveLoop(curves)

    nf = 72
    cx, cz = 0.5 * chord, 0.0
    p_far = [geo.addPoint(cx + far * cos(2 * pi * i / nf), 0.0,
                          cz + far * sin(2 * pi * i / nf)) for i in range(nf)]
    arcs = [geo.addLine(p_far[i], p_far[(i + 1) % nf]) for i in range(nf)]
    loop_far = geo.addCurveLoop(arcs[::-1])
    surf = geo.addPlaneSurface([loop_far, loop_af])
    geo.synchronize()


    gmsh.option.setNumber('Mesh.MeshSizeExtendFromBoundary', 0)
    gmsh.option.setNumber('Mesh.MeshSizeFromPoints', 0)
    gmsh.option.setNumber('Mesh.MeshSizeFromCurvature', 0)
    gmsh.option.setNumber('Mesh.MeshSizeMin', wall_size)
    gmsh.option.setNumber('Mesh.MeshSizeMax', far_size)
    gmsh.option.setNumber('Mesh.Algorithm', 5)
    # graded sizing in the section plane
    fd = g.mesh.field.add('Distance')
    g.mesh.field.setNumbers(fd, 'CurvesList', curves)
    g.mesh.field.setNumber(fd, 'Sampling', 80)
    ft = g.mesh.field.add('Threshold')
    g.mesh.field.setNumber(ft, 'InField', fd)
    g.mesh.field.setNumber(ft, 'SizeMin', wall_size)
    g.mesh.field.setNumber(ft, 'SizeMax', far_size)
    g.mesh.field.setNumber(ft, 'DistMin', 0.01)
    g.mesh.field.setNumber(ft, 'DistMax', 2.0)
    if bl_layers:
        # Boundary-layer extrusion in the SECTION plane: gmsh's BoundaryLayer
        # field is documented for 2-D, and here that is exactly what is
        # needed - the near-wall quads become hexahedral prisms once the
        # O-grid is extruded along the span. First height follows
        # y1 = y+ nu / u_tau, i.e. tens of microns, NOT millimetres; the
        # millimetre figure is a tangential surface size.
        blf = g.mesh.field.add('BoundaryLayer')
        g.mesh.field.setNumbers(blf, 'CurvesList', curves)
        g.mesh.field.setNumber(blf, 'Size', bl_h0)
        g.mesh.field.setNumber(blf, 'Ratio', bl_growth)
        g.mesh.field.setNumber(blf, 'Quads', 1)
        g.mesh.field.setNumber(blf, 'Thickness', bl_thick)
        g.mesh.field.setAsBoundaryLayer(blf)
    g.mesh.generate(2)

    # The spanwise direction must be discretised explicitly: OCC extrude with
    # no numElements produces a single 3.15 m cell across the span. Passing
    # numElements=[n_span] keeps the 2-D O-grid (so the farfield surface and
    # the section contours both survive) and layers it along the span; the
    # size field then only governs the in-section element size.
    out = geo.extrude([(2, surf)], 0.0, span, 0.0, [n_span],
                     recombine=bool(bl_layers))
    geo.synchronize()

    disk_face = None
    if disk:
        _ = g.occ
        d = g.occ.addDisk(disk['x'], disk['y'], disk['z'], disk['radius'],
                          disk['radius'])
        g.occ.rotate([(2, d)], disk['x'], disk['y'], disk['z'],
                     0.0, 1.0, 0.0, pi / 2.0)
        g.occ.synchronize()
        disk_face = d
        # the disc is far finer than the background size; without a matching
        # local size the volume mesher emits overlapping facets on it
        fd_d = g.mesh.field.add('Distance')
        g.mesh.field.setNumbers(fd_d, 'SurfacesList', [d])
        g.mesh.field.setNumber(fd_d, 'Sampling', 60)
        ft_d = g.mesh.field.add('Threshold')
        g.mesh.field.setNumber(ft_d, 'InField', fd_d)
        g.mesh.field.setNumber(ft_d, 'SizeMin', disk['radius'] / 5.0)
        g.mesh.field.setNumber(ft_d, 'SizeMax', far_size)
        g.mesh.field.setNumber(ft_d, 'DistMin', 0.0)
        g.mesh.field.setNumber(ft_d, 'DistMax', 1.5)
        bg = g.mesh.field.add('Min')
        g.mesh.field.setNumbers(bg, 'FieldsList', [ft, ft_d])
        g.mesh.field.setAsBackgroundMesh(bg)
        for _, vt in g.getEntities(3):
            g.mesh.embed(2, [d], 3, vt)

    g.mesh.generate(3)

    groups = {}
    cx_f, cz_f = 0.5 * chord, 0.0
    for d, t in g.getEntities(2):
        bx0, by0, bz0, bx1, by1, bz1 = g.getBoundingBox(2, t)
        if abs(by0) < 1e-6 and abs(by1) < 1e-6:
            groups.setdefault('symmetry', []).append(t)
        elif abs(by0 - span) < 1e-6 and abs(by1 - span) < 1e-6:
            groups.setdefault('tip', []).append(t)
        else:
            r = ((0.5 * (bx0 + bx1) - cx_f) ** 2 +
                 (0.5 * (bz0 + bz1) - cz_f) ** 2) ** 0.5
            groups.setdefault('wing' if r < 0.5 * far else 'farfield',
                              []).append(t)
    for name, ts in groups.items():
        g.addPhysicalGroup(2, ts, name=name)
    if disk_face is not None:
        g.addPhysicalGroup(2, [disk_face], name='disk')
    g.addPhysicalGroup(3, [t for _, t in g.getEntities(3)], name='fluid')
    geo.synchronize()

    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from su2mesh import write_su2_mesh, write_cgns_mesh
    if bl_layers:
        # hand-written prism connectivity fails SU2's face association, so
        # boundary-layer meshes go out as CGNS with gmsh's own connectivity
        npoin, nelem = write_cgns_mesh(gmsh, out_su2)
        nmark = len(gmsh.model.getPhysicalGroups(2))
    else:
        npoin, nelem, nmark = write_su2_mesh(gmsh, out_su2)
    gmsh.finalize()
    return {'npoin': npoin, 'nelem': nelem, 'nmark': nmark,
            'groups': {k: len(v) for k, v in groups.items()}}


if __name__ == '__main__':
    root = Path('/home/yanlu/Documents/plane')
    print(build(str(root / 'su2/cases/HW.su2'),
                str(root / 'data/airfoils/sd7062.dat'), delta_deg=20.0,
                verbose=True))
