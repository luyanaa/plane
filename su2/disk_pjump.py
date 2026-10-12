"""Gate 3b: thrust from the pressure jump across the actuator disk.

Measured from the VOLUME solution at two planes just outside the disk faces,
because surface.vtk on an internal ACTDISK marker returns the domain state and
not the BC-imposed fictitious face state (the readout there gave a 15 kPa drop,
contradicting CEulerSolver::BC_ActDisk_VariableLoad, which sets
Pressure_out = Pressure_in + Fa).

    T = (p_down - p_up) * A_disk

is exact by actuator-disk theory: the static pressure rises across a
propulsive disk, with p2 - p1 = 2 rho V^2 a (1+a).
"""
import sys
from pathlib import Path

import numpy as np

R = Path('/home/yanlu/Documents/plane')
sys.path.insert(0, str(R / 'su2'))

from diskonly import read_wall_fields, V, RHO, R_PROP, A_DISK   # noqa: E402


def plane_mean(pts, p, x0, tol, r_max):
    m = (np.abs(pts[:, 0] - x0) < tol) & (np.hypot(pts[:, 1], pts[:, 2]) < r_max)
    if m.sum() == 0:
        return None, 0
    return float(p[m].mean()), int(m.sum())


def main():
    vtk = R / 'su2/cases/DSKrun/vol_solution.vtk'
    if not vtk.exists():
        print('no volume output; enable PARAVIEW_ASCII')
        return
    pts, tri, scal, vect = read_wall_fields(str(vtk))
    p = scal['Pressure']
    print(f'volume: {len(pts)} points, p {p.min():.0f}..{p.max():.0f} Pa')

    up, nu = plane_mean(pts, p, -0.06, 0.07, 0.9 * R_PROP)
    dn, nd = plane_mean(pts, p, 0.14, 0.07, 0.9 * R_PROP)
    if up is None or dn is None:
        import collections
        near = pts[np.abs(pts[:, 0]) < 0.4, 0]
        print('x stations near disk:', sorted(set(np.round(near, 3)))[:12])
        return
    print(f'upstream  plane x=-0.06: n={nu:4d}  p={up:.1f} Pa')
    print(f'downstream plane x=+0.14: n={nd:4d}  p={dn:.1f} Pa')
    t = (dn - up) * A_DISK
    print(f'pressure jump = {dn - up:+.1f} Pa   (theory: T/A = '
          f'{100.31 / A_DISK:+.1f} Pa mean)')
    print(f'T_CFD = {t:.2f} N   T_input = 100.31 N   '
          f'error = {100 * (t - 100.31) / 100.31:+.2f}%')

    # induced velocity check: u_x just downstream vs freestream
    ux = vect['Momentum'][:, 0] / scal['Density']
    mu = (np.abs(pts[:, 0] + 0.06) < 0.07) & (np.hypot(pts[:, 1], pts[:, 2]) < 0.9 * R_PROP)
    md = (np.abs(pts[:, 0] - 0.14) < 0.07) & (np.hypot(pts[:, 1], pts[:, 2]) < 0.9 * R_PROP)
    if mu.sum() and md.sum():
        u1, u2 = ux[mu].mean(), ux[md].mean()
        print(f'u_x: upstream {u1:.2f}  downstream {u2:.2f}  V_inf {V:.2f} m/s')
        print(f'   a = (u2-V)/V = {(u2 - V) / V:+.4f}')


if __name__ == '__main__':
    main()
