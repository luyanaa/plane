"""Prototype D operating points for OF3R, with corrected units.

`Aircraft._prop_full(V, rho)` returns `(thrust, n)` where n is in REVOLUTIONS
PER SECOND (it bisects over 5..190 rps against the engine power curve), not
rpm. Reading it as rpm gave J = 11.99 at takeoff; as rps it gives J = 0.1998,
which matches results/prop_operating.csv exactly.

The three points also sit at three different speeds - takeoff is at
v_r = 1.1 V_sTO = 16.62 m/s, not the 14.07 m/s used in an earlier draft:

    takeoff   V = 16.622 m/s   n = 157.0 rps = 9420 rpm   J = 0.1998
    climb     V = 17.717 m/s   n = 157.1 rps = 9426 rpm   J = 0.2128
    cruise    V = 21.989 m/s   n =  99.8 rps = 5990 rpm   J = 0.4156

Usage:  python scripts/prop_points.py
"""
import json
import math
import sys
from pathlib import Path

import numpy as np

R = Path('/home/yanlu/Documents/plane')
sys.path.insert(0, str(R))

from sizing.atmosphere import RHO0                          # noqa: E402
from sizing.configs import MISSION, PROTOTYPES              # noqa: E402
from sizing.mission import Aircraft                         # noqa: E402
from sizing.prop_bemt import PropBEMT                       # noqa: E402

KEY = 'D'


def main():
    frame, pack = PROTOTYPES[KEY]['frame'], PROTOTYPES[KEY]['pack']
    ac = Aircraft(frame, pack, MISSION)
    prop = pack.prop
    met = ac.metrics(RHO0)
    A = math.pi * (prop.D / 2.0) ** 2

    v_r = met['V_s_TO'] * 1.1
    _, n_to = ac._prop_full(v_r, RHO0)          # (thrust, n[rps])
    _, n_cl = ac._prop_full(met['V_best_ROC'], RHO0)

    pts = {
        'TO':     dict(V=v_r, n_rps=n_to),
        'climb':  dict(V=met['V_best_ROC'], n_rps=n_cl),
        'cruise': dict(V=met['V_cruise'], n_rps=met['n_cruise'] / 60.0),
    }

    print(f'prototype {KEY}: D={prop.D:.4f} m  B={prop.B}  P/D={prop.P_over_D}  '
          f'n_eng={pack.n_eng}  MTOW={frame.MTOW_kg:.1f} kg')
    print(f'V_s_TO={met["V_s_TO"]:.3f}  V_ROC={met["V_best_ROC"]:.3f}  '
          f'V_cruise={met["V_cruise"]:.3f} m/s')
    print()
    hdr = (f'{"point":>7} {"V[m/s]":>8} {"n[rpm]":>9} {"J":>7} {"T[N]":>8} '
           f'{"Q[Nm]":>8} {"P[W]":>8} {"Ct_std":>7} {"a_mom":>6} '
           f'{"u_th":>6} {"q_sw%":>6}')
    print(hdr)
    out = {}
    for k, d in pts.items():
        V, n = d['V'], d['n_rps']
        p = PropBEMT(D=prop.D, P_over_D=prop.P_over_D, B=prop.B)
        T, Q, P = p.solve(V, n, iters=120)
        J = V / (n * prop.D)
        ct_std = T / (0.5 * RHO0 * A * V ** 2)
        a = (math.sqrt(1.0 + ct_std) - 1.0) / 2.0
        # swirl
        _, _, _, _, dQdr, mu, r = p._solve_elements(V, n, iters=120)
        mu = np.asarray(mu, float)
        r = mu * p.R
        om = 2.0 * math.pi * n
        ap = dQdr / (4.0 * math.pi * RHO0 * r ** 3 * V * om)
        ut = ap * om * r
        g = mu > 0.3
        q_sw = 0.5 * RHO0 * np.mean(ut[g] ** 2) / (0.5 * RHO0 * V ** 2)
        print(f'{k:>7} {V:8.3f} {n * 60:9.1f} {J:7.4f} {T:8.2f} {Q:8.4f} '
              f'{P:8.1f} {ct_std:7.3f} {a:6.3f} {np.mean(ut[g]):6.2f} '
              f'{q_sw * 100:5.1f}')
        out[k] = dict(V=V, n_rps=n, n_rpm=n * 60, J=J, T=T, Q=Q, P=P,
                      Ct_std=ct_std, a_mom=a, u_theta_mean=float(np.mean(ut[g])),
                      q_swirl_frac=q_sw)

    dest = R / 'results/prop_points_D.json'
    dest.write_text(json.dumps(out, indent=2))
    print('\nsaved', dest)
    print('\nLayer-8 BEMT validity note: the Layer-8 model carries a J < 0.25')
    print('warning, and all three points sit at or below J = 0.42, with TO and')
    print('climb deep inside it. Layer-8 radial loads at TO are therefore NOT')
    print('ground truth for rotorDiskSource; only the integrated T/P and the')
    print('momentum-theory consistency checks are used as constraints.')


if __name__ == '__main__':
    main()
