"""SU2 mesh generation helpers (gmsh -> CGNS).

2-D airfoil mesh for the layer-12 smoke/V&V cases, and the 3-D
prop-wing-flap half-model builder used for the D/C closure study.
"""
from pathlib import Path
from math import cos, pi, radians, sin

import numpy as np


def read_airfoil(path):
    """Read a .dat/Selig airfoil, return (x, y) ordered TE -> upper -> LE -> lower -> TE."""
    pts = []
    for line in Path(path).read_text().splitlines():
        s = line.strip()
        if not s or s.startswith('#') or s.startswith('%'):
            continue
        parts = s.split()
        if len(parts) >= 2:
            try:
                pts.append((float(parts[0]), float(parts[1])))
            except ValueError:
                continue
    p = np.array(pts)
    return p[:, 0], p[:, 1]


def rebuild_airfoil(x, y, n=200, te_cut=0.99):
    """Split into upper/lower surfaces with a blunt TE at te_cut*chord.

    Returns (xu, yu, xl, yl): upper from TE_blunt -> LE -> lower to TE_blunt,
    each cosine-resampled. The TE face is the straight segment between the
    two TE points (avoids the degenerate sharp-TE line that hangs gmsh).
    """
    i_le = int(np.argmin(x))
    xu, yu = x[:i_le + 1][::-1], y[:i_le + 1][::-1]     # LE -> TE (upper)
    xl, yl = x[i_le:], y[i_le:]                          # LE -> TE (lower)
    # resample on a normalized arc coordinate, then cut at te_cut
    def resample(xs, ys):
        s = np.linspace(0.0, 1.0, len(xs))
        m = xs <= te_cut
        s2 = 0.5 * (1.0 - np.cos(np.linspace(0.0, pi, n)))
        xr = np.interp(s2, s[m], xs[m])
        yr = np.interp(s2, s[m], ys[m])
        return xr, yr
    xu_r, yu_r = resample(xu, yu)
    xl_r, yl_r = resample(xl, yl)
    return xu_r, yu_r, xl_r, yl_r


def mesh_airfoil_2d(af_dat, out_cgns, chord=1.0, far=25.0, n_surf=220,
                    bl_layers=32, bl_h0=1e-5, bl_growth=1.15, bl_thick=0.06,
                    surf_size=0.002, far_size=0.8, ramp_chords=4.0,
                    verbose=False, generate=True, mesh_algo=5):
    """O-grid airfoil mesh with a boundary layer; writes .cgns.

    Upper and lower surfaces are separate splines joined by a blunt TE face.
    """
    import gmsh
    xu, yu, xl, yl = rebuild_airfoil(*read_airfoil(af_dat), n=n_surf)
    gmsh.initialize()
    gmsh.option.setNumber('General.Verbosity', 3 if verbose else 0)
    gmsh.option.setNumber('Geometry.OCCAutoFix', 0)
    gmsh.model.add('airfoil')
    geo = gmsh.model.geo

    # shared leading-edge point, then upper/lower point chains
    le = geo.addPoint(xu[0] * chord, yu[0] * chord, 0.0)
    pu = [geo.addPoint(x * chord, y * chord, 0.0) for x, y in zip(xu[1:], yu[1:])]
    pl = [geo.addPoint(x * chord, y * chord, 0.0) for x, y in zip(xl[1:], yl[1:])]
    # pu / pl run LE->TE; build the loop as
    #   su: LE -> TE_upper,  te: TE_upper -> TE_lower,  sl: TE_lower -> LE
    su = geo.addSpline([le] + pu)
    sl = geo.addSpline(pl[::-1] + [le])
    te = geo.addLine(pu[-1], pl[-1])
    loop_af = geo.addCurveLoop([su, te, sl])

    # farfield circle
    nf = 60
    fp = [geo.addPoint(0.5 * chord + far * chord * cos(2 * pi * i / nf),
                       far * chord * sin(2 * pi * i / nf), 0.0) for i in range(nf)]
    arcs = [geo.addLine(fp[i], fp[(i + 1) % nf]) for i in range(nf)]

    loop_far = geo.addCurveLoop(arcs[::-1])
    surf = geo.addPlaneSurface([loop_far, loop_af])
    geo.synchronize()

    # size field: fine at the wall, growing to the farfield
    f_dist = gmsh.model.mesh.field.add('Distance')
    gmsh.model.mesh.field.setNumbers(f_dist, 'CurvesList', [su, sl, te])
    gmsh.model.mesh.field.setNumber(f_dist, 'Sampling', 200)
    f_th = gmsh.model.mesh.field.add('Threshold')
    gmsh.model.mesh.field.setNumber(f_th, 'InField', f_dist)
    gmsh.model.mesh.field.setNumber(f_th, 'SizeMin', surf_size * chord)
    gmsh.model.mesh.field.setNumber(f_th, 'SizeMax', far_size * chord)
    gmsh.model.mesh.field.setNumber(f_th, 'DistMin', 0.005 * chord)
    gmsh.model.mesh.field.setNumber(f_th, 'DistMax', ramp_chords * chord)

    fields = [f_th]
    if bl_layers > 0:
        blf = gmsh.model.mesh.field.add('BoundaryLayer')
        gmsh.model.mesh.field.setNumbers(blf, 'CurvesList', [su, sl, te])
        gmsh.model.mesh.field.setNumber(blf, 'Size', bl_h0 * chord)
        gmsh.model.mesh.field.setNumber(blf, 'Ratio', bl_growth)
        gmsh.model.mesh.field.setNumber(blf, 'Quads', 1)
        gmsh.model.mesh.field.setNumber(blf, 'Thickness', bl_thick * chord)
        gmsh.model.mesh.field.setAsBoundaryLayer(blf)   # applies *in addition*
        fields.append(blf)
    bg = gmsh.model.mesh.field.add('Min')
    gmsh.model.mesh.field.setNumbers(bg, 'FieldsList', fields)
    gmsh.model.mesh.field.setAsBackgroundMesh(bg)
    gmsh.option.setNumber('Mesh.MeshSizeExtendFromBoundary', 0)
    gmsh.option.setNumber('Mesh.MeshSizeFromPoints', 0)
    gmsh.option.setNumber('Mesh.MeshSizeFromCurvature', 0)

    if generate:
        gmsh.option.setNumber('Mesh.Algorithm', mesh_algo)
    gmsh.model.mesh.generate(2)
    gmsh.model.addPhysicalGroup(1, [su, te, sl], name='airfoil')
    gmsh.model.addPhysicalGroup(1, arcs, name='farfield')
    gmsh.model.addPhysicalGroup(2, [surf], name='fluid')
    geo.synchronize()
    nn = len(gmsh.model.mesh.getNodes()[0])
    ne = sum(len(e) for e in gmsh.model.mesh.getElements()[1])
    if generate:
        import sys as _sys
        _sys.path.insert(0, str(Path(__file__).resolve().parent))
        from su2mesh import write_su2_mesh
        write_su2_mesh(gmsh, out_cgns)
    gmsh.finalize()
    return nn, ne
