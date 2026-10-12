"""Layer 12: D half-model mesh (inner wing + deflected flap + tractor nacelle +
actuator disk) for the SU2 prop-wing-flap closure.

The mesh is built directly in SU2's aerodynamic frame (the ONERA M6
convention), because SU2 rotates the freestream by AOA within its own x-y
plane and therefore expects lift along +y and the span along z:

    x  chordwise, wing LE at x = 0, flow +x
    y  vertical, wing chord plane at y = 0, lift direction
    z  spanwise, symmetry plane at z = 0, modelled panel to z = Y_TIP

Building in this frame avoids any post-hoc coordinate permutation, which is
either improper (y/z exchange inverts every element) or flips the wing.

The actuator disk is an internal baffle pair (two discs DISK_T apart); SU2's
VARIABLE_LOAD double-surface BC wants two coincident faces. The nacelle is a
body of revolution on the tractor axis, ahead of the wing LE.
"""
from math import cos, pi, radians, sin
from pathlib import Path

import numpy as np

# ---- D geometry (layers 1/4) ----------------------------------------------
C_ROOT = 0.831          # root chord, m
TAPER = 0.60
B_FULL = 6.317          # full span, m
Y_ENG = 0.75            # engine/nacelle spanwise station (model z)
D_PROP = 0.53           # propeller diameter, m
R_PROP = D_PROP / 2.0
Y_TIP = B_FULL / 2.0    # model the full half-span (symmetry plane at z = 0)
FLAP_SPAN = 0.48 * B_FULL / 2.0   # flap covers 0..FLAP_SPAN (layer 3)
FLAP_FRAC = 0.25        # flap chord fraction (layer 3)
NAC_R = 0.085           # nacelle max radius
NAC_Z = 0.055           # nacelle axis height above the chord plane (model y)
NAC_X0, NAC_X1 = -0.45, -0.05   # nacelle extent (ahead of the LE)
DISK_X = -0.52          # actuator disk station (ahead of the nacelle)
DISK_T = 0.05           # disk inlet/outlet separation, m (must exceed the
                        # wall mesh size or the mesher stalls on slivers)


def chord_at(y):
    return C_ROOT * (1.0 - 2.0 * (1.0 - TAPER) * y / B_FULL)


def read_section(dat, n=140, te_cut=0.99):
    """Return (xu, yu), (xl, yl) cosine-resampled LE -> TE, unit chord.

    The second coordinate is the section thickness direction.
    """
    raw = []
    for ln in Path(dat).read_text().splitlines():
        ln = ln.strip()
        if not ln or ln[0] in '#%':
            continue
        try:
            raw.append([float(v) for v in ln.split()])
        except ValueError:
            continue                      # airfoil name header
    arr = np.array(raw)
    x, z = arr[:, 0], arr[:, 1]
    i_le = int(np.argmin(x))
    xu, zu = x[:i_le + 1][::-1], z[:i_le + 1][::-1]
    xl, zl = x[i_le:], z[i_le:]
    s2 = 0.5 * (1.0 - np.cos(np.linspace(0.0, pi, n)))

    def rs(xs, zs):
        s = np.linspace(0.0, 1.0, len(xs))
        m = xs <= te_cut
        return np.interp(s2, s[m], xs[m]), np.interp(s2, s[m], zs[m])

    a, b = rs(xu, zu)
    c, d = rs(xl, zl)
    return (a, b), (c, d)


def deflected_section(dat, chord, delta_deg, flap_frac=FLAP_FRAC):
    """Section in metres with the aft flap_frac rotated about the hinge,
    TE-down positive. Returns (xc, yt) arrays, upper then lower."""
    (xu, zu), (xl, zl) = read_section(dat)
    xh = 1.0 - flap_frac
    d = radians(delta_deg)

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


def _section_wire(g, dat, chord, span_z, delta, flap_frac=FLAP_FRAC):
    """Closed OCC wire for one span station, built in the SU2 frame:
    (x, y, z) = (chordwise, vertical, spanwise)."""
    (xu, zu), (xl, zl) = deflected_section(dat, chord, delta, flap_frac)
    xs = list(xu[::-1]) + list(xl[1:])       # TE -> upper -> LE -> lower -> TE
    ys = list(zu[::-1]) + list(zl[1:])

    pts = []
    for x, y in zip(xs, ys):
        if pts and abs(x - pts[-1][0]) < 1e-9 and abs(y - pts[-1][1]) < 1e-9:
            continue
        pts.append((x, y))
    if abs(pts[-1][0] - pts[0][0]) < 1e-9 and abs(pts[-1][1] - pts[0][1]) < 1e-9:
        pts.pop()
    tags = [g.occ.addPoint(float(x), float(span_z), float(y)) for x, y in pts]
    su = g.occ.addSpline(tags)
    te = g.occ.addLine(tags[-1], tags[0])
    return g.occ.addCurveLoop([su, te])


def build(out_su2, dat, delta_deg=20.0, z_tip=Y_TIP, far=4.0,
          mesh_size=0.05, far_size=0.5, wall_size=None, disk=True,
          nacelle=True, verbose=False):
    """Build the half-model mesh and write it. Returns marker info."""
    import gmsh
    gmsh.initialize()
    gmsh.option.setNumber('General.Terminal', 1 if verbose else 0)
    gmsh.model.add('D_half')
    g = gmsh.model
    occ = g.occ

    # ---- wing + flap solid (span along z) ----------------------------
    # three spanwise stations: root, flap end (deflected), tip (clean). The
    # flap covers only the inboard FLAP_SPAN, as in layer 3/11.
    z_flap = min(FLAP_SPAN, z_tip)
    st = [(0.0, delta_deg), (z_flap, delta_deg), (z_tip, 0.0)]
    loops = [_section_wire(g, dat, chord_at(z), z, d) for z, d in st]
    wing = occ.addThruSections(loops, makeSolid=True, makeRuled=False)
    sol = [(3, wing[0][1])]

    # ---- nacelle (ahead of the LE, at z = Y_ENG, height y = NAC_Z) ---
    if nacelle:
        nac = occ.addCylinder(NAC_X0, Y_ENG, NAC_Z,
                              NAC_X1 - NAC_X0, 0.0, 0.0, NAC_R)
        sol.append((3, nac))

    # ---- actuator disk: a thin solid carved out of the fluid ---------
    # Cutting a solid slab leaves two flat faces DISK_T apart - exactly the
    # inlet/outlet pair SU2's ACTDISK_DOUBLE_SURFACE expects. A 2-D baffle
    # instead stalls the mesher on slivers in the gap.
    if disk:
        dsk = occ.addCylinder(DISK_X, Y_ENG, NAC_Z, DISK_T, 0.0, 0.0, R_PROP)
        sol.append((3, dsk))

    # ---- fluid box with the body carved out -------------------------
    x0, x1 = -2.5, 3.5
    y0, y1 = 0.0, far
    z0, z1 = -far * 0.6, far * 0.6
    box = occ.addBox(x0, y0, z0, x1 - x0, y1 - y0, z1 - z0)

    # boolean subtract: fragment() does not reliably carve a cavity when a
    # solid lies fully inside another, and without the cavity the wing faces
    # become internal baffles the flow passes straight through (CL ~ 0).
    carved, _ = occ.cut([(3, box)], sol, removeTool=True)
    occ.synchronize()
    fluid = [t for d, t in carved if d == 3]



    # ---- classify surfaces -------------------------------------------
    eps = 1e-6
    groups = {}
    for d, t in g.getEntities(2):
        bx0, by0, bz0, bx1, by1, bz1 = occ.getBoundingBox(2, t)
        cx, cy, cz = (bx0 + bx1) / 2, (by0 + by1) / 2, (bz0 + bz1) / 2
        if disk and abs(bx1 - bx0) < 1e-4 and abs(bz1 - bz0 - D_PROP) < 1e-3 \
                and abs(cy - Y_ENG) < 1e-3 and abs(cz - NAC_Z) < 1e-3:
            groups.setdefault('disk_in' if abs(cx - DISK_X) <
                              abs(cx - (DISK_X + DISK_T)) else 'disk_out',
                              []).append(t)
        elif abs(by0 - y0) < eps and abs(by1 - y0) < eps:
            groups.setdefault('symmetry', []).append(t)
        elif (abs(bx0 - x0) < eps or abs(bx1 - x1) < eps or
              abs(by1 - y1) < eps or abs(bz0 - z0) < eps or
              abs(bz1 - z1) < eps):
            groups.setdefault('farfield', []).append(t)
        elif nacelle and abs(cx - (NAC_X0 + NAC_X1) / 2) < 0.25 and \
                abs(cy - Y_ENG) < 0.12 and abs(cz - NAC_Z) < 0.12 and \
                max(bx1 - bx0, by1 - by0, bz1 - bz0) < 0.45:
            groups.setdefault('nacelle', []).append(t)
        else:
            groups.setdefault('wing', []).append(t)

    for name, ts in groups.items():
        if ts:
            g.addPhysicalGroup(2, ts, name=name)
    g.addPhysicalGroup(3, fluid, name='fluid')
    occ.synchronize()

    # ---- graded size field + mesh ------------------------------------
    body = [t for name in ('wing', 'nacelle')
            for t in groups.get(name, [])]
    f_dist = gmsh.model.mesh.field.add('Distance')
    gmsh.model.mesh.field.setNumbers(f_dist, 'SurfacesList', body)
    gmsh.model.mesh.field.setNumber(f_dist, 'Sampling', 40)
    wall_size = wall_size or mesh_size
    f_th = gmsh.model.mesh.field.add('Threshold')
    gmsh.model.mesh.field.setNumber(f_th, 'InField', f_dist)
    gmsh.model.mesh.field.setNumber(f_th, 'SizeMin', wall_size)
    gmsh.model.mesh.field.setNumber(f_th, 'SizeMax', far_size)
    gmsh.model.mesh.field.setNumber(f_th, 'DistMin', 0.005)
    gmsh.model.mesh.field.setNumber(f_th, 'DistMax', 1.2)
    gmsh.model.mesh.field.setAsBackgroundMesh(f_th)
    gmsh.option.setNumber('Mesh.MeshSizeExtendFromBoundary', 0)
    gmsh.option.setNumber('Mesh.MeshSizeFromPoints', 0)
    gmsh.option.setNumber('Mesh.MeshSizeFromCurvature', 14)
    gmsh.option.setNumber('Mesh.Algorithm3D', 1)
    gmsh.model.mesh.generate(3)

    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from su2mesh import write_su2_mesh
    npoin, nelem, nmark = write_su2_mesh(gmsh, out_su2)
    gmsh.finalize()
    return {'npoin': npoin, 'nelem': nelem, 'nmark': nmark,
            'groups': {k: len(v) for k, v in groups.items()}}


if __name__ == '__main__':
    root = Path('/home/yanlu/Documents/plane')
    info = build(str(root / 'su2/cases/D_half.su2'),
                 str(root / 'data/airfoils/sd7062.dat'),
                 delta_deg=20.0, verbose=False)
    print(info)
