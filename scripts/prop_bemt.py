"""Layer 8: BEMT propeller curves, operating points, mission closure comparison.

Per prototype:
  - C_T(J), C_P(J), eta_p(J) curves (fixed reference RPM)
  - operating points: takeoff (full throttle at V_R), climb (full at best-ROC V),
    cruise (from the coupled engine-map+BEMT closure in metrics())
  - endurance/range with BEMT eta_p(J) vs the eta_p=0.75 baseline

Emits results/prop_curves.png, results/prop_operating.csv.
"""
import sys
from dataclasses import replace
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from sizing.atmosphere import RHO0
from sizing.configs import MISSION, PROTOTYPES
from sizing.mission import Aircraft

OUT = Path(__file__).resolve().parents[1] / 'results'


def main():
    rows = []
    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    for ax, (key, cfg) in zip(axes.ravel(), PROTOTYPES.items()):
        frame, pack = cfg['frame'], cfg['pack']
        ac = Aircraft(frame, pack, MISSION)
        prop = pack.prop
        js, ct, cp, eta = prop.curves()
        ax.plot(js, ct, label='C_T')
        ax.plot(js, cp, label='C_P')
        ax.plot(js, eta, label='eta_p')
        ax.axvspan(0, 0.25, alpha=0.12, color='gray')
        ax.text(0.08, 0.85, 'J<0.25:\nmomentum\nmodel zone', fontsize=7, ha='center')
        ax.set_title(f'{key} prop D={prop.D:.2f} m, P/D={prop.P_over_D:.2f}, B={prop.B}')
        ax.set_xlabel('J')
        ax.grid(alpha=0.3)
        ax.legend()
        # operating points via the mission closure
        met = ac.metrics(RHO0)
        old = Aircraft(frame, replace(pack, prop=None), MISSION).metrics(RHO0)
        v_r = met['V_s_TO'] * 1.1
        t_full, n_full = ac._prop_full(v_r, RHO0)
        j_to = v_r / (n_full * prop.D)
        v_cl = met['V_best_ROC']
        t_cl, n_cl = ac._prop_full(v_cl, RHO0)
        j_cl = v_cl / (n_cl * prop.D)
        for jj, tag, c in ((j_to, 'TO', 'g'), (j_cl, 'climb', 'm'), (met['J_cruise'], 'cruise', 'r')):
            ax.axvline(jj, ls='--', c=c, lw=0.8)
            ax.annotate(tag, (jj, 0.02), rotation=90, fontsize=7, c=c)
        rows.append({'prototype': key, 'D': prop.D, 'P_over_D': prop.P_over_D,
                     'J_TO': j_to, 'J_climb': j_cl, 'J_cruise': met['J_cruise'],
                     'eta_cruise': met['eta_p_cruise'], 'n_cruise': met['n_cruise'],
                     'throttle': met['throttle_frac'],
                     'T_static': ac._prop_full(0.5, RHO0)[0],
                     'E_new_h': met['endurance_h'], 'E_old_h': old['endurance_h'],
                     'R_new_km': met['range_km'], 'R_old_km': old['range_km'],
                     'E_delta_pct': 100 * (met['endurance_h'] / old['endurance_h'] - 1)})
    fig.tight_layout()
    fig.savefig(OUT / 'prop_curves.png', dpi=120)
    plt.close(fig)
    df = pd.DataFrame(rows).set_index('prototype')
    df.to_csv(OUT / 'prop_operating.csv')
    pd.set_option('display.width', 220)
    print(df.to_string(float_format=lambda x: f'{x:8.3f}'))
    print(f'\nsaved: {OUT}/prop_curves.png, prop_operating.csv')


if __name__ == '__main__':
    main()
