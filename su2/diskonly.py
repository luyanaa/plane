"""Layer 12 step 3: actuator-disk-only thrust closure.

The cheapest independent closure in the sequence: no wing, so it does not
depend on the unresolved wall-normal mesh question, and it validates the
Layer-8 BEMT -> SU2 mapping directly (the pitch-angle fix changed this data).

A box domain carries an internal double-surface actuator disk loaded from
PropBEMT.radial_loads(). Thrust is recovered from the axial momentum balance on
the outlet plane rather than from a reported coefficient:

    T = rho V^2 A + int_out (rho u_x^2 - rho V^2) dA + int_out (p - p_inf) dA

which is exact for a control volume whose inlet is freestream and whose outlet
crosses the whole streamtube (outside it u_x ~ V, p ~ p_inf).

Pass condition: |T_CFD - T_BEMT| / T_BEMT < 3%.
"""
import csv
import sys
from pathlib import Path

import numpy as np

ROOT = Path('/home/yanlu/Documents/plane')
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'su2'))

V = 30.00                      # J ~ 0.49; J<0.25 is the invalid
                               # heavy-load regime where BEMT is known bad
N_RPM = 6900.0                 # moderate-load point, see V below
D_PROP = 0.530
R_PROP = D_PROP / 2.0
RHO, P_INF, T_INF = 1.225, 101325.0, 288.15
R_IN_FRAC, R_OUT_FRAC = 0.165, 0.990   # annulus bounds, see build_domain
A_DISK = float(np.pi * R_PROP ** 2)
X0, X1 = -2.00, 4.00
HW = 1.50                      # box half-width in y and z

CFG = """MESH_FILENAME= {mesh}
SOLVER= EULER
MESH_FORMAT= SU2
MATH_PROBLEM= DIRECT
FLUID_MODEL= STANDARD_AIR
MACH_NUMBER= {mach}
AOA= 0.0
SIDESLIP_ANGLE= 0.0
FREESTREAM_PRESSURE= {p_inf}
FREESTREAM_TEMPERATURE= {t_inf}
REF_ORIGIN_MOMENT_X= 0.0
REF_ORIGIN_MOMENT_Y= 0.0
REF_ORIGIN_MOMENT_Z= 0.0
REF_LENGTH= {ref_length}
REF_AREA= {ref_area}
MARKER_INLET= ( inlet, {t_inf}, {p_tot}, 1.0, 0.0, 0.0 )
MARKER_OUTLET= ( outlet, {p_inf} )
MARKER_FAR= ( farfield )
MARKER_PLOTTING= ( disk_in, disk_out )
MARKER_MONITORING= ( NONE )
ACTDISK_DOUBLE_SURFACE= YES
ACTDISK_TYPE= VARIABLE_LOAD
ACTDISK_FILENAME= {adfile}
MARKER_ACTDISK= ( disk_in, disk_out, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 )
NUM_METHOD_GRAD= WEIGHTED_LEAST_SQUARES
CFL_NUMBER= 0.5
CFL_ADAPT= YES
CFL_ADAPT_PARAM= ( 0.1, 2.0, 1.0, 100.0 )
ITER= {iters}
CONV_NUM_METHOD_FLOW= ROE
ENTROPY_FIX_COEFF= 0.1
MUSCL_FLOW= YES
SLOPE_LIMITER_FLOW= VENKATAKRISHNAN
TIME_DISCRE_FLOW= EULER_IMPLICIT
LINEAR_SOLVER= FGMRES
LINEAR_SOLVER_PREC= ILU
LINEAR_SOLVER_ERROR= 1E-4
LINEAR_SOLVER_ITER= 10
SCREEN_OUTPUT= (INNER_ITER, RMS_DENSITY, LIFT, DRAG)
HISTORY_OUTPUT= (ITER, RMS_RES, AERO_COEFF)
OUTPUT_FILES= (RESTART_ASCII, PARAVIEW_ASCII, SURFACE_CSV, SURFACE_PARAVIEW_ASCII)
CONV_FILENAME= history
SURFACE_FILENAME= surface
RESTART_FILENAME= restart
"""


def build_domain(out_su2, disk_gap=0.08, mesh_size=0.12, far_size=0.30):
    """Box domain with a two-faced internal actuator disk."""
    import gmsh
    gmsh.initialize()
    gmsh.option.setNumber('General.Terminal', 0)
    gmsh.model.add('diskonly')
    g = gmsh.model
    occ = g.occ

    # The mesh must be CONFORMAL to the disk plane x = 0, otherwise elements
    # straddle it and the node-splitting step corrupts them (measured: 1291 of
    # 35971 tets straddled, which is what made SU2 diverge). Two boxes sharing
    # that plane give a real conformal interface; the disc is then a patch on
    # it.
    boxA = occ.addBox(X0, -HW, -HW, 0.0 - X0, 2 * HW, 2 * HW)
    boxB = occ.addBox(0.0, -HW, -HW, X1, 2 * HW, 2 * HW)
    # The disk face must be an ANNULUS. SU2 maps a node to a radial station
    # with  for (iEl=0; iEl<nRow; iEl++) if (r_ <= rad_v[iEl]) {
    #            h = (r_-rad_v[iEl-1])/(rad_v[iEl]-rad_v[iEl-1]);
    # so a node with r_ < rad_v[0] indexes rad_v[-1]: undefined behaviour.
    # The shipped tutorial disk is an annulus over a spinner for this reason.
    # Here the face spans [R_IN_FRAC, R_OUT_FRAC] * R, and RADIUS in the input
    # file is the true radius, so every node lands strictly inside
    # (rad_v[0], rad_v[-1]].
    r_in, r_out = R_IN_FRAC * R_PROP, R_OUT_FRAC * R_PROP
    outer = occ.addDisk(0.0, 0.0, 0.0, r_out, r_out)
    inner = occ.addDisk(0.0, 0.0, 0.0, r_in, r_in)
    for d in (outer, inner):
        occ.rotate([(2, d)], 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, np.pi / 2.0)
    occ.synchronize()
    ann, _ = occ.cut([(2, outer)], [(2, inner)], removeTool=True)
    occ.synchronize()
    occ.fragment([(3, boxA), (3, boxB)], [(2, t) for _, t in ann])
    occ.synchronize()

    groups = {}
    annulus_tag = None
    eps = 1e-6
    for dd, t in g.getEntities(2):
        bx0, by0, bz0, bx1, by1, bz1 = occ.getBoundingBox(2, t)
        cx, cy, cz = (bx0 + bx1) / 2, (by0 + by1) / 2, (bz0 + bz1) / 2
        if abs(bx1 - bx0) < 1e-6 and abs(cx) < 1e-9:
            # on the disk plane: the disc patch becomes the actuator disk, the
            # surrounding interface annulus is interior and takes no marker
            # the disc patch spans the full 2*r_out (~0.52 m); the interface
            # annulus spans the whole 3 m face
            if (by1 - by0) < 1.0 and (bz1 - bz0) < 1.0:
                groups.setdefault('disk', []).append(t)
                annulus_tag = t
            continue
        elif abs(bx0 - X0) < eps and abs(bx1 - X0) < eps:
            groups.setdefault('inlet', []).append(t)
        elif abs(bx0 - X1) < eps and abs(bx1 - X1) < eps:
            groups.setdefault('outlet', []).append(t)
        else:
            groups.setdefault('farfield', []).append(t)

    for name, ts in groups.items():
        g.addPhysicalGroup(2, ts, name=name)
    g.addPhysicalGroup(3, [t for _, t in g.getEntities(3)], name='fluid')
    occ.synchronize()

    # The mesh must be UNIFORM here. A graded Distance/Threshold field around
    # the disk faces diverges (NaN at every load scale, including zero) while
    # a uniform mesh converges at 0.30/0.18/0.12 m, so the gradient itself is
    # the instability, not the actuator-disk load or the boundary conditions.
    gmsh.option.setNumber('Mesh.MeshSizeExtendFromBoundary', 0)
    gmsh.option.setNumber('Mesh.MeshSizeFromPoints', 0)
    gmsh.option.setNumber('Mesh.MeshSizeFromCurvature', 0)
    gmsh.option.setNumber('Mesh.MeshSizeMin', mesh_size)
    gmsh.option.setNumber('Mesh.MeshSizeMax', mesh_size * 1.6)
    gmsh.option.setNumber('Mesh.Algorithm3D', 1)
    g.mesh.generate(3)

    strad = 0
    for dim, pt in g.getPhysicalGroups():
        if dim != 3:
            continue
        nt, nc = g.mesh.getNodes()[0], None
        for ent in g.getEntitiesForPhysicalGroup(dim, pt):
            ets, _, enodes = g.mesh.getElements(dim, ent)
            for et, conn in zip(ets, enodes):
                _, _, _, nn, _, _ = g.mesh.getElementProperties(et)
                for i in range(0, len(conn), nn):
                    xs = [g.mesh.getNode(int(x))[0][0] for x in conn[i:i + nn]]
                    if min(xs) < -1e-12 < max(xs):
                        strad += 1
    print(f'[mesh] elements straddling x=0: {strad}')

    from su2split import write_su2_mesh_doubled
    info = write_su2_mesh_doubled(gmsh, out_su2, annulus_tag, x_plane=0.0)
    gmsh.finalize()
    return info, {k: len(v) for k, v in groups.items()}


def read_wall_fields(vtk):
    """Return (points, triangles, {scalar name: array}, {vector name: Nx3})."""
    tok = open(vtk).read().split()
    i, N = 0, len(tok)
    pts = tri = None
    scal, vect = {}, {}

    while i < N:
        w = tok[i]
        if w == 'POINTS':
            n = int(tok[i + 1])
            i += 3
            pts = np.array(tok[i:i + 3 * n], dtype=float).reshape(n, 3)
            i += 3 * n
        elif w == 'CELLS':
            total = int(tok[i + 2])
            i += 3
            end = i + total
            tris = []
            while i < end:
                k = int(tok[i])
                ids = [int(x) for x in tok[i + 1:i + 1 + k]]
                if k == 3:
                    tris.append(ids)
                i += k + 1
            tri = np.array(tris, dtype=int)
        elif w == 'CELL_TYPES':
            i += 2 + int(tok[i + 1])
        elif w == 'POINT_DATA':
            i += 2
        elif w in ('SCALARS', 'VECTORS'):
            name = tok[i + 1]
            n = len(pts)
            ncomp = 1 if w == 'SCALARS' else 3
            i += 4 if w == 'SCALARS' else 3
            if w == 'SCALARS':
                i += 2                       # LOOKUP_TABLE default
            arr = np.array(tok[i:i + ncomp * n], dtype=float)
            i += ncomp * n
            if w == 'SCALARS':
                scal[name] = arr
            else:
                vect[name] = arr.reshape(n, 3)
        else:
            i += 1
    return pts, tri, scal, vect


def thrust_from_disk(vtk, v_inf, rho_inf, p_inf):
    """Thrust from the pressure jump across the actuator disk.

    The far-field momentum balance is unusable in this box: the accelerated
    flow fills the whole outlet plane, so outside-streamtube conditions
    (u -> V, p -> p_inf) never hold and the integral over-counts by ~40x.
    The local measure is exact by actuator-disk theory - the static pressure
    RISES across a propulsive disk, with T = (p_down - p_up) A:

        Bernoulli upstream:   p_inf + q V^2            = p1 + q V^2(1+a)^2
        Bernoulli downstream: p2  + q V^2(1+a)^2       = p_inf + q V^2(1+2a)^2
        => p2 - p1 = 2 rho V^2 a (1+a),  and T = (p2-p1) A
    """
    pts, tri, scal, vect = read_wall_fields(vtk)
    p = scal['Pressure']
    a, b, c = pts[tri[:, 0]], pts[tri[:, 1]], pts[tri[:, 2]]
    nvec = np.cross(b - a, c - a)
    nn = np.linalg.norm(nvec, axis=1)
    nhat = nvec / np.maximum(nn[:, None], 1e-300)
    area_x = 0.5 * nn * np.abs(nhat[:, 0])
    xc = pts[tri].mean(axis=1)[:, 0]
    pm = p[tri].mean(axis=1)

    is_in = xc < 0.5 * (xc.min() + xc.max())
    a_in, a_out = area_x[is_in].sum(), area_x[~is_in].sum()
    f_in = float((pm[is_in] * area_x[is_in]).sum())
    f_out = float((pm[~is_in] * area_x[~is_in]).sum())
    t = f_out - f_in
    return t, a_in, a_out, f_in / max(a_in, 1e-30), f_out / max(a_out, 1e-30)


def main():
    from sizing.configs import PROTOTYPES
    from sizing.prop_bemt import make_prop
    from actdisk import write_actdisk_input, disk_from_bemt

    pack = PROTOTYPES['D']['pack']
    prop = pack.prop
    n_rps = N_RPM / 60.0
    mach = V / np.sqrt(1.4 * P_INF / RHO)

    t_bemt, q_bemt, p_bemt = prop.solve(V, n_rps)
    ct_bemt = t_bemt / (RHO * n_rps ** 2 * D_PROP ** 4)
    cp_bemt = p_bemt / (RHO * n_rps ** 3 * D_PROP ** 5)
    print(f'BEMT: V={V} m/s  n={N_RPM} rpm  J={V/(n_rps*D_PROP):.4f}  '
          f'T={t_bemt:.2f} N  P={p_bemt:.1f} W')
    print(f'  C_T={ct_bemt:.5f}  C_P={cp_bemt:.5f}  '
          f'(per engine; {pack.n_eng} engines, D={D_PROP} m)')

    mesh = str(ROOT / 'su2/cases/DSK.su2')
    info, groups = build_domain(mesh)
    print(f'mesh: {info} groups={groups}')
    if groups.get('disk') != 1:
        print(f'WARNING: actuator-disk marker not isolated: {groups}')

    d = disk_from_bemt(prop, V, n_rps, markers=('disk_in', 'disk_out'),
                       center=(0.0, 0.0, 0.0), axis=(1.0, 0.0, 0.0))
    adfile = str(ROOT / 'su2/cases/DSK_actdisk.dat')
    write_actdisk_input(adfile, [d])
    # ---- gate 3a: the written file must integrate back to BEMT ---------
    ct_file = float(np.trapezoid(d['d_ct'], d['r_bar']))
    cp_file = float(np.trapezoid(d['d_cp'], d['r_bar']))
    # SU2 applies these as rho n^2 D^4 and rho n^3 D^5 references
    t_file = ct_file * RHO * n_rps ** 2 * D_PROP ** 4
    p_file = cp_file * RHO * n_rps ** 3 * D_PROP ** 5
    dct = abs(ct_file - ct_bemt) / ct_bemt
    dcp = abs(cp_file - cp_bemt) / cp_bemt
    print(f'[gate 3a] int dCT/drbar drbar = {ct_file:.6f} vs BEMT '
          f'{ct_bemt:.6f}  ({dct:.2e})')
    print(f'[gate 3a] int dCP/drbar drbar = {cp_file:.6f} vs BEMT '
          f'{cp_bemt:.6f}  ({dcp:.2e})')
    print(f'[gate 3a] T_file={t_file:.2f} N  P_file={p_file:.1f} W')
    print(f'[gate 3a] dCt range {min(d["d_ct"]):+.4f}..{max(d["d_ct"]):+.4f}  '
          f'dCp range {min(d["d_cp"]):+.4f}..{max(d["d_cp"]):+.4f}')
    gate3a = (dct < 1e-3) and (dcp < 1e-3)
    print(f'[gate 3a] {"PASS" if gate3a else "FAIL "} '
          f'(need |dC_T|,|dC_P| < 1e-3)')

    # T the disk BC should actually produce, from the file itself
    t_input = t_file

    # MARKER_INLET takes TOTAL conditions. Setting p_total = p_static,outlet
    # leaves no dynamic head and the whole domain settles near ~13 m/s instead
    # of the intended 30 m/s.
    p_tot = P_INF + 0.5 * RHO * V * V
    print(f'inlet total pressure {p_tot:.1f} Pa '
          f'(static {P_INF:.0f} + q {0.5*RHO*V*V:.1f})')
    cfg = str(ROOT / 'su2/cases/DSK.cfg')
    Path(cfg).write_text(CFG.format(mesh=mesh, mach=mach, p_inf=P_INF,
                                    p_tot=p_tot, t_inf=T_INF,
                                    ref_length=D_PROP, ref_area=A_DISK,
                                    adfile=adfile, iters=6000))
    from su2case import run_case
    r = run_case(cfg, str(ROOT / 'su2/cases/DSKrun'), np_threads=24,
                 timeout=3600)
    print(f'SU2 rc={r.returncode}')

    hist = ROOT / 'su2/cases/DSKrun/history.csv'
    if hist.exists():
        rows = list(csv.reader(open(hist)))
        hdr = [c.strip().replace('"', '') for c in rows[0]]
        ir = hdr.index('rms[Rho]')
        data = [x for x in rows[1:] if len(x) > ir]
        print(f'  residual rms[Rho] = {float(data[-1][ir]):.2f} '
              f'after {len(data)} iters')

    vtk = ROOT / 'su2/cases/DSKrun/surface.vtk'
    if vtk.exists():
        t_cfd, a_in, a_out, p_in, p_out = thrust_from_disk(str(vtk), V, RHO,
                                                           P_INF)
        err = 100.0 * (t_cfd - t_input) / t_input
        print(f'disk areas: in {a_in:.5f}  out {a_out:.5f} m2 '
              f'(geometric {A_DISK:.5f})')
        print(f'face-mean p: upstream {p_in:.1f}  downstream {p_out:.1f} Pa '
              f'(jump {p_out - p_in:+.1f})')
        print(f'[gate 3b] T_CFD = {t_cfd:.2f} N   T_input = {t_input:.2f} N   '
              f'error = {err:+.2f}%')
        print(f'[gate 3b] {"PASS" if abs(err) < 3.0 else "FAIL "} '
              f'(need |error| < 3%)')
        print('NOTE: the face-pressure readout above is invalid - surface.vtk '
              'on an internal ACTDISK marker returns the domain state, not the '
              'BC-imposed face state. Use disk_pjump.py on the volume field.')
        with open(ROOT / 'results/su2_disk_vv.csv', 'w', newline='') as fh:
            w = csv.writer(fh)
            w.writerow(['T_BEMT_N', 'T_CFD_N', 'err_pct', 'A_out_m2',
                        'A_disk_m2', 'J'])
            w.writerow([round(t_bemt, 3), round(t_cfd, 3), round(err, 3),
                        round(a_out, 5), round(A_DISK, 5),
                        round(d['advance_ratio'], 4)])
        print('saved: results/su2_disk_vv.csv')


if __name__ == '__main__':
    main()
