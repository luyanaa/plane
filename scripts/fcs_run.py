"""Layer 10: FCS test battery on the 6-DoF model.

Scenarios (per prototype): cruise hold, 8 m/s crosswind, gust + turbulence,
engine-out (twins), GPS loss, airspeed bias. Metrics per scenario with
pass/fail thresholds. Emits results/fcs_metrics.csv + fcs_*.png histories.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from sizing.atmosphere import RHO0
from sizing.autopilot import Autopilot
from sizing.configs import PROTOTYPES
from sizing.sixdof import SixDoF

OUT = Path(__file__).resolve().parents[1] / 'results'
DT = 0.02


def gust_u(t, start, dur, amp):
    """1-cos discrete gust profile."""
    return amp * 0.5 * (1 - np.cos(2 * np.pi * (t - start) / dur)) \
        if start <= t <= start + dur else 0.0


def sim(m, ap, t_end, scenario):
    """Run the FCS loop; returns recorded arrays + summary dict."""
    x = ap._x0.copy()
    rec = {'t': [], 'h': [], 'V': [], 'psi': [], 'phi': [], 'nz': [], 'de': [],
           'da': [], 'dr': [], 'thr': [], 'beta': []}
    turb = 0.0
    rng = np.random.default_rng(1)
    for step in range(int(t_end / DT)):
        t = step * DT
        # wind
        w = np.zeros(3)
        if scenario == 'crosswind':
            w[1] = 8.0                                # 8 m/s from the right
        elif scenario == 'gust':
            w[0] = gust_u(t, 5.0, 2.0, 10.0) + turb
            turb += DT / 0.5 * (2.0 * rng.standard_normal() - turb)
        m.set_wind(w)
        # sensors
        gps_ok = not (scenario == 'gps_loss' and t > 20.0)
        v_bias = 3.0 if scenario == 'airspeed_bias' else 0.0
        ctrl = ap.update(x, DT, gps_ok=gps_ok, V_bias=v_bias)
        if scenario == 'engine_out' and t > 5.0:
            ctrl[3] = 0.0                            # left engine off
        # RK4
        k1 = m.derivs(x, ctrl)
        k2 = m.derivs(x + 0.5 * DT * k1, ctrl)
        k3 = m.derivs(x + 0.5 * DT * k2, ctrl)
        k4 = m.derivs(x + DT * k3, ctrl)
        x += DT / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
        ua, va, wa = m.air_relative(x)
        V = np.sqrt(ua * ua + va * va + wa * wa)
        nz = (-m.derivs(x, ctrl)[2] + 9.80665) / 9.80665 if V > 1 else 0
        rec['t'].append(t); rec['h'].append(x[11]); rec['V'].append(V)
        rec['psi'].append(np.degrees(x[8])); rec['phi'].append(np.degrees(x[6]))
        rec['nz'].append(nz); rec['beta'].append(np.degrees(np.arcsin(va / V)))
        rec['de'].append(np.degrees(ctrl[0])); rec['da'].append(np.degrees(ctrl[1]))
        rec['dr'].append(np.degrees(ctrl[2])); rec['thr'].append(ctrl[3])
    return rec


def main():
    rows = []
    histories = {}
    for key, cfg in PROTOTYPES.items():
        m = SixDoF(key, cfg['frame'], cfg['pack'])
        v_cr = m.met['V_cruise']
        z, _ = m.trim_level(v_cr)
        m.trim_V = v_cr
        x0 = np.array([v_cr * np.cos(z[0]), 0.0, v_cr * np.sin(z[0]),
                       0, 0, 0, 0, z[0], 0, 0, 0, 500.0])
        for scen in ('cruise_hold', 'crosswind', 'gust', 'gps_loss', 'airspeed_bias'):
            if scen == 'engine_out' and m.pack.n_eng != 2:
                continue
            ap = Autopilot(m, v_cr)
            ap._x0 = x0.copy()
            ap._h_prev = None
            rec = sim(m, ap, 40.0, scen)
            histories[f'{key}_{scen}'] = rec
            dh = max(abs(np.array(rec['h']) - 500.0))
            dpsi = max(abs(np.array(rec['psi'])))
            dV = max(abs(np.array(rec['V']) - v_cr))
            nzmax = max(abs(np.array(rec['nz'])))
            dmx = np.max([np.max(np.abs(np.array(rec[k]))) for k in ('de', 'da', 'dr')])
            sat = np.mean([abs(v) > 24.5 for v in rec['de']])
            beta_ss = rec['beta'][-1]
            ok = dh < 25 and dpsi < 10 and dV < 2.5 and nzmax < 3.0
            rows.append({'prototype': key, 'scenario': scen,
                         'max_dh_m': dh, 'max_dpsi_deg': dpsi, 'max_dV': dV,
                         'max_nz': nzmax, 'max_de_da_dr': dmx, 'sat_frac': sat,
                         'beta_ss_deg': beta_ss, 'PASS': ok})
        # engine-out with FCS
        if m.pack.n_eng == 2:
            ap = Autopilot(m, v_cr)
            ap._x0 = x0.copy()
            ap._h_prev = None
            rec = sim(m, ap, 40.0, 'engine_out')
            histories[f'{key}_engine_out'] = rec
            dh = max(abs(np.array(rec['h']) - 500.0))
            dpsi = max(abs(np.array(rec['psi'])))
            bank = max(abs(np.array(rec['phi'])))
            ok = dpsi < 90 and bank < 45 and dh < 200
            rows.append({'prototype': key, 'scenario': 'engine_out',
                         'max_dh_m': dh, 'max_dpsi_deg': dpsi,
                         'max_dV': max(abs(np.array(rec['V']) - v_cr)),
                         'max_nz': max(abs(np.array(rec['nz']))),
                         'max_de_da_dr': max(abs(np.array(rec['dr']))),
                         'sat_frac': np.mean([abs(v) > 24.5 for v in rec['dr']]),
                         'beta_ss_deg': rec['beta'][-1], 'PASS': ok})
    df = pd.DataFrame(rows)
    df.to_csv(OUT / 'fcs_metrics.csv', index=False)
    pd.set_option('display.width', 220)
    print(df.to_string(index=False, float_format=lambda x: f'{x:7.2f}'))
    print(f'\nsaved: {OUT}/fcs_metrics.csv')


if __name__ == '__main__':
    main()
