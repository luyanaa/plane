"""OF3R closure: rotorDiskSource vs the Layer-8 BEMT.

Rotordisk emits both axial and tangential blade force, so the readout has to be
a momentum/angular-momentum balance, not a pressure integral (the wing's
`forces.py` route does not apply - there is no wall here).

    T_wake = rho int u_x (u_x - V) dA          axial momentum excess
    Q_wake = rho int r u_theta u_x dA          angular momentum flux
    P      = Q * Omega                         shaft power

Acceptance (user's):
    cruise : |dT| < 5 %  AND  |dQ| (or |dP|) < 5-8 %
    TO     : wider, but record a genuine divergence if T/P differ by > 10 %
             rather than tuning the blade to force agreement

Layer-8 carries a J < 0.25 validity warning and every point here is at or below
J = 0.42, so at TO the Layer-8 radial loads are NOT ground truth: the rotorDisk
model is an INDEPENDENT blade-element formulation and only the integrated
quantities are compared.
"""
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np

R = Path('/home/yanlu/Documents/plane')
sys.path.insert(0, str(R / 'foam'))

import disk3d                                              # noqa: E402
import wake as W                                           # noqa: E402

RHO = 1.225
SHELL = 0.030
X_STATION = 1.40


def load_bemt():
    return json.loads((R / 'results/prop_points_D.json').read_text())


def measure(case, Rc):
    dirs = [p for p in (case / 'VTK').glob(f'{case.name}_*') if p.is_dir()]
    vtu = sorted(dirs)[-1] / 'internal.vtu'
    pts, conn, cell = W.read_vtu(vtu)
    cen, vol = W.cell_geometry(pts, conn)
    U = cell['U']
    m = np.abs(cen[:, 0] - X_STATION) < 0.5 * SHELL
    rr = np.hypot(cen[m, 1], cen[m, 2])
    k = rr <= Rc
    dA = vol[m][k] / SHELL
    ux = U[m][k, 0]
    # tangential component: azimuthal unit vector (-y, z)/r in the y-z plane
    y, z = cen[m][k, 1], cen[m][k, 2]
    r = np.hypot(y, z)
    ut = np.where(r > 1e-9, (-U[m][k, 1] * y + U[m][k, 2] * z) / r, 0.0)
    return dict(npix=int(k.sum()), ux=ux, ut=ut, r=r, dA=dA)


def main():
    point = sys.argv[1] if len(sys.argv) > 1 else 'cruise'
    case = R / f'foam/cases/rotor_{point}'
    if not case.is_absolute():
        case = R / case
    bem = load_bemt()[point]
    V, n_rps = bem['V'], bem['n_rps']
    om = 2.0 * math.pi * n_rps
    Rc = 1.30 * disk3d.R_PROP

    mm = measure(case, Rc)
    T = float(RHO * np.sum(mm['ux'] * (mm['ux'] - V) * mm['dA']))
    Q = float(RHO * np.sum(mm['r'] * mm['ut'] * mm['ux'] * mm['dA']))
    P = Q * om

    print(f'--- OF3R {point}: V={V:.3f} m/s  n={n_rps * 60:.1f} rpm  '
          f'J={bem["J"]:.4f}  Omega={om:.1f} rad/s')
    print(f'    station x={X_STATION} m ({X_STATION / disk3d.D_PROP:.2f} D), '
          f'Rc={Rc / disk3d.R_PROP:.2f} R, {mm["npix"]} cells')
    print()
    print(f'{"quantity":>10} {"CFD":>10} {"BEMT":>10} {"diff":>9}')
    for name, a, b in (('T [N]', T, bem['T']), ('Q [Nm]', Q, bem['Q']),
                       ('P [W]', P, bem['P'])):
        d = 100.0 * (a - b) / b if b else float('nan')
        print(f'{name:>10} {a:10.3f} {b:10.3f} {d:+8.2f}%')

    dT = abs(T - bem['T']) / bem['T'] * 100
    dQ = abs(Q - bem['Q']) / bem['Q'] * 100
    dP = abs(P - bem['P']) / bem['P'] * 100
    print()
    if point == 'cruise':
        ok = dT < 5.0 and min(dQ, dP) < 8.0
        print(f'gate (cruise): |dT|<5% -> {dT:.2f}%, |dQ|<8% -> {dQ:.2f}%  '
              f'=> {"PASS" if ok else "FAIL"}')
    else:
        print(f'{point}: |dT|={dT:.2f}%  |dP|={dP:.2f}%'
              + ('  -> formulations DIVERGE at heavy loading (>10%): record, '
                 'do not tune' if max(dT, dP) > 10.0 else ''))

    # tangential inflow angle at the flap-relevant radius, the quantity the
    # scalar q_slip/q_inf ROM cannot represent
    r_b = mm['r']
    sel = (r_b > 0.5 * disk3d.R_PROP) & (r_b < 1.05 * disk3d.R_PROP)
    print()
    print(f'mean swirl angle in r/R 0.5-1.05: '
          f'{math.degrees(math.atan(np.mean(mm["ut"][sel]) / np.mean(mm["ux"][sel]))):.2f} deg'
          f'   (u_theta mean {np.mean(mm["ut"][sel]):.3f} m/s)')

    out = R / f'results/of_rotor_{point}.json'
    out.write_text(json.dumps(dict(point=point, V=V, n_rpm=n_rps * 60,
                                   T=T, Q=Q, P=P, T_bemt=bem['T'],
                                   Q_bemt=bem['Q'], P_bemt=bem['P'],
                                   dT_pct=dT, dQ_pct=dQ, dP_pct=dP),
                              indent=2))
    print('saved', out)


if __name__ == '__main__':
    main()
