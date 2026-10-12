"""Layer 9: 6-DoF trim, dynamic modes, and engine-out V_MC (prototypes A-D).

Uses sizing/sixdof.py (nonlinear rigid-body EOM with the VSPAERO aero tables,
engine map + BEMT thrust grid). Emits results/sixdof_trim.csv,
results/sixdof_modes.csv, results/sixdof_vmc.csv.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from sizing.atmosphere import RHO0
from sizing.configs import PROTOTYPES
from sizing.sixdof import SixDoF

OUT = Path(__file__).resolve().parents[1] / 'results'


def state_at_trim(m, z):
    alpha, de, th, thr = z
    return np.array([m.trim_V * np.cos(alpha), 0.0, m.trim_V * np.sin(alpha),
                     0, 0, 0, 0, th, 0, 0, 0, 500.0])


def linearize(m, x0, ctrl, states, eps=1e-4):
    n = len(states)
    a = np.zeros((n, n))
    for i, si in enumerate(states):
        for sgn in (-1, 1):
            xp = x0.copy()
            xp[si] += sgn * eps
            a[:, i] += sgn * m.derivs(xp, ctrl)[[s for s in states]] / (2 * eps)
    return a


def mode_report(name, a):
    out = []
    ev = np.linalg.eigvals(a)
    used = set()
    for i, lam in enumerate(ev):
        if i in used or lam.real > 0:
            continue
        if abs(lam.imag) > 1e-6:
            j = int(np.argmin(abs(ev - lam.conj())))
            used |= {i, j}
            w = abs(lam.imag) / (2 * np.pi)
            z = -lam.real / abs(lam)
            out.append((name, 'osc', w, z, 1.0 / abs(lam.real)))
        else:
            out.append((name, 'real', 0.0, 1.0,
                        -1.0 / lam.real if abs(lam.real) > 1e-8 else float('inf')))
    return out


def main():
    rows_t, rows_m = [], []
    for key, cfg in PROTOTYPES.items():
        m = SixDoF(key, cfg['frame'], cfg['pack'])
        met = m.met
        trims = [('cruise', met['V_cruise']), ('approach', 1.3 * met['V_s_clean']),
                 ('V_roc', met['V_best_ROC'])]
        for tag, v in trims:
            z, res = m.trim_level(v)
            m.trim_V = v
            rows_t.append({'prototype': key, 'case': tag, 'V': v,
                           'alpha_deg': np.degrees(z[0]), 'theta_deg': np.degrees(z[2]),
                           'elevator_deg': np.degrees(z[1]), 'throttle': z[3],
                           'residual': res.max()})
        # modes around cruise trim
        z, _ = m.trim_level(met['V_cruise'])
        m.trim_V = met['V_cruise']
        x0 = state_at_trim(m, z)
        ctrl = [z[1], 0.0, 0.0, z[3], z[3]]
        a_long = linearize(m, x0, ctrl, [0, 2, 4, 7])
        a_lat = linearize(m, x0, ctrl, [1, 3, 5, 6])
        for r in mode_report('long', a_long) + mode_report('lat', a_lat):
            rows_m.append({'prototype': key, 'axis': r[0], 'type': r[1],
                           'freq_hz': r[2], 'damping': r[3], 'tau_s': r[4]})
        print(f'{key}: trim cruise alpha={np.degrees(z[0]):.1f} deg, '
              f'de={np.degrees(z[1]):.1f} deg, throttle={z[3]:.2f}')

    # engine-out V_MC for the twins (C, D)
    rows_v = []
    for key in ('C', 'D'):
        m = SixDoF(key, PROTOTYPES[key]['frame'], PROTOTYPES[key]['pack'])
        v_s = m.met['V_s_land']
        v_r = m.met['V_s_TO'] * 1.1
        vmc = None
        for v in np.linspace(v_s * 1.05, m.met['V_cruise'] * 1.5, 12):
            z, _ = m.trim_level(v)
            m.trim_V = v
            x0 = state_at_trim(m, z)
            ctrl = [z[1], 0.0, 1.0, 0.0, 1.0]   # max rudder, left off, right full
            x = x0.copy()
            dt, t, max_bank, max_yaw, h0, uncontrolled = 0.02, 0.0, 0.0, 0.0, x0[11], False
            while t < 20.0 and not uncontrolled:
                k1 = m.derivs(x, ctrl)
                k2 = m.derivs(x + 0.5 * dt * k1, ctrl)
                k3 = m.derivs(x + 0.5 * dt * k2, ctrl)
                k4 = m.derivs(x + dt * k3, ctrl)
                x += dt / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
                bank = abs(np.degrees(x[6]))
                max_bank = max(max_bank, bank)
                max_yaw = max(max_yaw, abs(np.degrees(x[8])))
                if not np.isfinite(x).all() or bank > 85.0:
                    uncontrolled = True
                t += dt
            if not uncontrolled and max_bank <= 15.0 and vmc is None:
                vmc = v
            rows_v.append({'prototype': key, 'V': v, 'V_over_VR': v / v_r,
                           'max_bank_deg': max_bank if not uncontrolled else float('nan'),
                           'max_yaw_deg': max_yaw if not uncontrolled else float('nan'),
                           'alt_loss_m': (h0 - x[11]) if not uncontrolled else float('nan'),
                           'uncontrolled': uncontrolled})
        print(f'{key}: V_MC ~ {vmc:.1f} m/s = {vmc / v_r:.2f} V_R'
              if vmc else f'{key}: uncontrolled up to 1.5x V_cr')
    pd.DataFrame(rows_t).to_csv(OUT / 'sixdof_trim.csv', index=False)
    pd.DataFrame(rows_m).to_csv(OUT / 'sixdof_modes.csv', index=False)
    pd.DataFrame(rows_v).to_csv(OUT / 'sixdof_vmc.csv', index=False)
    print('\n=== modes ===')
    print(pd.DataFrame(rows_m).to_string(index=False, float_format=lambda x: f'{x:8.3f}'))
    print('\n=== engine-out (C/D) ===')
    print(pd.DataFrame(rows_v).to_string(index=False, float_format=lambda x: f'{x:8.2f}'))
    print(f'\nsaved: {OUT}/sixdof_trim.csv, sixdof_modes.csv, sixdof_vmc.csv')


if __name__ == '__main__':
    main()
