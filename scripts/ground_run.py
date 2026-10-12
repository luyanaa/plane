"""Layer 11: takeoff/landing distance per flap state; is the big flap worth it?

For each prototype, compute s_TO and s_landing for cruise/takeoff/landing flap
settings, with and without propwash flap enhancement, and break out the
contributions. Emits results/ground_run.csv (+ png).
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from sizing.configs import PROTOTYPES
from sizing.ground import GroundRun, FLAP_DB, AIRFOIL_OF, FLAP_SPAN_FRAC

OUT = Path(__file__).resolve().parents[1] / 'results'


def main():
    rows, detail = [], {}
    for key, cfg in PROTOTYPES.items():
        gr = GroundRun(key, cfg['frame'], cfg['pack'])
        d = {}
        for state in ('clean', 'takeoff', 'landing'):
            to = gr.takeoff(state=state)
            ld = gr.landing(state=state)
            d[state] = (to, ld)
            q_to = gr.q_slip_ratio(to.V_r, gr.thrust(to.V_r, 1.0) * gr.pack.n_eng)
            thr_ld = 0.25
            q_ld = gr.q_slip_ratio(ld.V_approach,
                                   gr.thrust(ld.V_approach, thr_ld) * gr.pack.n_eng)
            rows.append({'prototype': key, 'state': state,
                         'clmax_eff_TO': gr.clmax_effective(state, True, to.V_r, 1.0),
                         'q_ratio_TO': q_to,
                         'clmax_eff_LDG': gr.clmax_effective(state, True, ld.V_approach, thr_ld),
                         'q_ratio_LDG': q_ld,
                         'V_s_TO': to.V_s, 'V_R': to.V_r,
                         's_ground_m': to.s_ground, 's_air_m': to.s_air,
                         's_TO_m': to.s_total, 'accel_g': to.accel_g,
                         'V_s_land': ld.V_s, 'V_app': ld.V_approach,
                         's_land_air_m': ld.s_air, 's_land_ground_m': ld.s_ground,
                         's_land_total_m': ld.s_total})
        # propwash contribution: takeoff with vs without the flap boost
        gr_noboost = GroundRun(key, cfg['frame'], cfg['pack'])
        gr_noboost.boost_frac = 0.0      # same geometry, no slipstream
        to_nb = gr_noboost.takeoff(state='takeoff')
        to_b = gr.takeoff(state='takeoff')
        detail[key] = (to_b.s_total, to_nb.s_total,
                       d['clean'][0].s_total, d['landing'][0].s_total)
        print(f'{key}: s_TO(clean)={d["clean"][0].s_total:6.1f} '
              f'(TO flap)={to_b.s_total:6.1f} (land flap)={d["landing"][0].s_total:6.1f} m | '
              f'propwash boost saves {to_nb.s_total - to_b.s_total:5.1f} m | '
              f's_land={d["landing"][1].s_total:6.1f} m')

    df = pd.DataFrame(rows)
    df.to_csv(OUT / 'ground_run.csv', index=False)

    # chart: s_TO and s_landing by flap state
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    keys = list(PROTOTYPES)
    states = ['clean', 'takeoff', 'landing']
    w = 0.26
    for ax, metric, title in ((axes[0], 's_TO_m', 'takeoff distance'),
                              (axes[1], 's_land_total_m', 'landing distance')):
        for i, st in enumerate(states):
            vals = [df[(df.prototype == k) & (df.state == st)][metric].iloc[0] for k in keys]
            ax.bar(np.arange(len(keys)) + (i - 1) * w, vals, w, label=st)
        ax.set_xticks(np.arange(len(keys)))
        ax.set_xticklabels(keys)
        ax.set_ylabel('m')
        ax.set_title(title)
        ax.legend()
        ax.grid(axis='y', alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / 'ground_run.png', dpi=120)
    print(f'\nsaved: {OUT}/ground_run.csv, ground_run.png')

    # the question: is the big flap worth it?
    print('\n=== is the landing flap worth it? (s_land_total vs clean) ===')
    for k in keys:
        cl = df[(df.prototype == k) & (df.state == 'clean')].iloc[0]
        la = df[(df.prototype == k) & (df.state == 'landing')].iloc[0]
        d_cl = cl.s_land_total_m
        d_la = la.s_land_total_m
        print(f'  {k}: clean {d_cl:6.1f} m -> landing flap {d_la:6.1f} m '
              f'({100*(d_la/d_cl-1):+.0f}%), V_app {cl.V_app:4.1f} -> {la.V_app:4.1f} m/s')


if __name__ == '__main__':
    main()
