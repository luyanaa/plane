"""Layer 12 step 1: independent force regression from the SU2 surface VTK.

Integrates the wall pressure over the triangulated surface,

    F_body = -oint p n dA

with n the outward (fluid-side) normal, then non-dimensionalises with q_inf
and the reference area. Agreement with SU2's reported CL validates the frame,
the marker set and REF_AREA together - it is the check that catches the classic
"half-model forces against full-wing reference area" factor-of-two error, which
no residual or convergence criterion would ever reveal.

The pressure-only integral reproduces CL to high accuracy at this Reynolds
number (the viscous normal-stress contribution to lift is O(0.1%)). CD needs
the skin-friction term as well, so only CL is used as the pass criterion;
CD is reported for information (and is expected to be *lower* than SU2's).

Usage:  python su2/force_check.py <surface.vtk> <v_inf> <p_inf> <rho_inf> <S_ref>
"""
import sys

import numpy as np


def read_su2_vtk(path):
    """Return (points Nx3, triangles Mx3, pressure N) from an SU2 ASCII VTK."""
    tok = open(path).read().split()
    i, N = 0, len(tok)
    pts = tri = pressure = None

    while i < N:
        w = tok[i]
        if w == 'POINTS':
            n = int(tok[i + 1])
            i += 3                                  # skip <n> <type>
            pts = np.array(tok[i:i + 3 * n], dtype=float).reshape(n, 3)
            i += 3 * n
        elif w == 'CELLS':
            nc, total = int(tok[i + 1]), int(tok[i + 2])
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
            nc = int(tok[i + 1])
            i += 2 + nc
        elif w == 'POINT_DATA':
            i += 2
        elif w == 'SCALARS':
            name = tok[i + 1]
            i += 4                                  # SCALARS <n> <t> 1
            i += 2                                  # LOOKUP_TABLE default
            n = len(pts)
            vals = np.array(tok[i:i + n], dtype=float)
            i += n
            if name == 'Pressure':
                pressure = vals
        elif w == 'VECTORS':
            i += 3                                  # VECTORS <n> <t>
            i += 3 * len(pts)
        else:
            i += 1
    return pts, tri, pressure


def main():
    vtk, v_inf, p_inf, rho_inf, s_ref = sys.argv[1:6]
    v_inf, p_inf, rho_inf, s_ref = map(float, (v_inf, p_inf, rho_inf, s_ref))
    pts, tri, p = read_su2_vtk(vtk)

    a, b, c = pts[tri[:, 0]], pts[tri[:, 1]], pts[tri[:, 2]]
    nvec = np.cross(b - a, c - a)
    a2 = np.linalg.norm(nvec, axis=1)
    nhat = nvec / np.maximum(a2[:, None], 1e-300)
    pm = p[tri].mean(axis=1)
    # outward-normal integral of p n dA, then the force on the body
    f = -(nhat * (pm * 0.5 * a2)[:, None]).sum(axis=0)

    q = 0.5 * rho_inf * v_inf * v_inf
    print(f'points {len(pts)}  triangles {len(tri)}  '
          f'p {p.min():.0f}..{p.max():.0f} Pa')
    print(f'F = ({f[0]:+.2f}, {f[1]:+.2f}, {f[2]:+.2f}) N')
    print(f'q_inf = {q:.2f} Pa   S_ref = {s_ref:.4f} m2')
    print(f'CL_manual = {f[2] / (q * s_ref):+.5f}')
    print(f'CD_manual = {f[0] / (q * s_ref):+.5f}  (pressure only)')


if __name__ == '__main__':
    main()
