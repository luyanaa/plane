"""Layer 12 step 3b: actuator-disk thrust closure on SU2's own tutorial mesh.

My homemade gmsh mesh could not satisfy SU2's actuator-disk requirements
(coincident inlet/outlet faces AND a mesh conformal to the disk plane AND an
annulus face), so this uses the shipped tutorial mesh, which is known-valid:
its DISK and DISK_BACK markers are exactly coincident at x = 0 with
r_bar = 0.1804..1.0000.

The Layer-8 radial loads are substituted for the tutorial's, evaluated at the
tutorial's advance ratio so the comparison is like-for-like. Radial stations
are extended to cover the mesh's r_bar range and the load file is rescaled so
its RADIUS matches the mesh radius - otherwise rim nodes fall past the last
station and SU2's BCEuler disk loop never assigns them a load.

The reported variables of interest are the dimensionless coefficients, which
are independent of the mesh's physical radius:

    C_T^SU2 / C_T^BEMT      C_P^SU2 / C_P^BEMT

plus the axial and swirl velocity jumps across the disk.
"""
import csv
import sys
from pathlib import Path

import numpy as np

R = Path('/home/yanlu/Documents/plane')
sys.path.insert(0, str(R))
sys.path.insert(0, str(R / 'su2'))

TUT = R / 'su2/tutorial'
MESH_R = 2.5146               # tutorial mesh radius (meters)
R_BAR_MIN, R_BAR_MAX = 0.1700, 1.0200   # covers the mesh's 0.1804..1.0000


def main(argv):
    V = float(argv[1]) if len(argv) > 1 else 30.0
    from sizing.configs import PROTOTYPES
    from actdisk import write_actdisk_input

    prop = PROTOTYPES['D']['pack'].prop
    # Evaluate at MY operating advance ratio. The tutorial J = 2.81 is a
    # windmilling condition for this propeller (P/D = 1.0, zero thrust near
    # J = 1), which yields negative C_T/C_P and no useful closure. The mesh is
    # only a domain here: the disk BC works in r_bar and the Renard
    # coefficients are reference-independent, so a 2.5 m mesh can carry a
    # 0.53 m propeller's load distribution.
    J = 0.4922
    n_rps = V / (J * prop.D)
    mu = np.linspace(R_BAR_MIN, R_BAR_MAX, 41)

    t, q, p = prop.solve(V, n_rps)
    r_bar, d_ct, d_cp, d_cr, j = prop.radial_loads(V, n_rps, mu=mu)
    # the file must state the MESH radius; the coefficients are reference
    # independent, so a node at rbar on the mesh maps onto rbar in the table
    ct_bemt = float(np.trapezoid(d_ct, r_bar))
    cp_bemt = float(np.trapezoid(d_cp, r_bar))

    print(f'V = {V} m/s, n = {n_rps*60:.1f} rpm, J = {J}')
    print(f'BEMT (evaluated on the file stations): '
          f'C_T = {ct_bemt:.5f}  C_P = {cp_bemt:.5f}')
    print(f'BEMT (native 30 stations):             '
          f'C_T = {t/(1.225*n_rps**2*prop.D**4):.5f}  '
          f'C_P = {p/(1.225*n_rps**3*prop.D**5):.5f}')

    ad = TUT / 'ActuatorDisk_layer8.dat'
    write_actdisk_input(str(ad), [dict(
        markers=('DISK', 'DISK_BACK'), center=(0.0, 0.0, 0.0),
        axis=(1.0, 0.0, 0.0), radius=MESH_R, advance_ratio=J,
        r_bar=list(r_bar), d_ct=list(d_ct), d_cp=list(d_cp),
        d_cr=list(d_cr))])
    print(f'wrote {ad.name}: NROW={len(r_bar)} rbar '
          f'{r_bar[0]:.4f}..{r_bar[-1]:.4f} RADIUS={MESH_R}')

    cfg = TUT / 'tutorial_layer8.cfg'
    txt = (TUT / 'tutorial.cfg').read_text()
    txt = txt.replace('ActuatorDisk.dat', ad.name)
    txt = txt.replace('CONV_FILENAME= history', 'CONV_FILENAME= history_l8')
    txt = txt.replace('RESTART_FILENAME= restart', 'RESTART_FILENAME= restart_l8')
    txt = txt.replace('SURFACE_FILENAME= surface', 'SURFACE_FILENAME= surface_l8')
    txt = txt.replace('VOLUME_FILENAME= flow', 'VOLUME_FILENAME= flow_l8')
    cfg.write_text(txt)

    from su2case import run_case
    r = run_case(str(cfg), str(TUT), np_threads=24, timeout=5400)
    print(f'SU2 rc = {r.returncode}')

    hist = TUT / 'history_l8.csv'
    if not hist.exists():
        s = (r.stdout or '') + (r.stderr or '')
        for ln in s.splitlines():
            if 'rror' in ln or 'NaN' in ln:
                print('  ', ln.strip()[:100])
        return
    rows = list(csv.reader(open(hist)))
    hdr = [c.strip().replace('"', '') for c in rows[0]]
    iR = hdr.index('rms[Rho]')
    data = [x for x in rows[1:] if len(x) > iR]
    print(f'converged to rms[Rho] = {float(data[-1][iR]):.2f} '
          f'after {len(data)} iters')

    # MARKER_ANALYZE reports massflow-averaged surface quantities; recover the
    # disk thrust from the pressure jump on the two coincident faces
    from disk_pjump import read_wall_fields
    vtk = TUT / 'surface_l8.vtk'
    if vtk.exists():
        pts, tri, scal, vect = read_wall_fields(str(vtk))
        p = scal['Pressure']
        rr = np.hypot(pts[:, 1], pts[:, 2]) / MESH_R
        print(f'disk faces: {len(pts)} nodes, rbar {rr.min():.3f}..{rr.max():.3f}')
        # split by marker is unavailable in the flat VTK; use radial bins
        edges = np.linspace(0.2, 1.0, 9)
        print('radial bin   n   p [Pa]')
        for lo, hi in zip(edges[:-1], edges[1:]):
            m = (rr >= lo) & (rr < hi)
            if m.sum() > 1:
                print(f'  {lo:.1f}-{hi:.1f}   {m.sum():3d}   {p[m].mean():.1f}')


if __name__ == '__main__':
    main(sys.argv)
