"""Four digital prototypes A-D, side-by-side layer-1 comparison.

Saves results/prototypes_table.csv + results/prototypes_compare.png.
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from sizing.configs import MISSION, PROTOTYPES
from sizing.mission import Aircraft

OUT = Path(__file__).resolve().parents[1] / 'results'

ROWS = ['MTOW_kg', 'M_empty_kg', 'M_fuel_kg', 'fuel_frac', 'W_S_kgm2', 'b_m',
        'V_s_land', 'V_cruise', 'C_L_cruise', 'L_D_max', 'P_cruise_kW',
        'throttle_frac', 'ROC', 's_TO_m', 'T0_W', 'endurance_h', 'range_km', 'feasible']

UNIT = {'MTOW_kg': 'kg', 'M_empty_kg': 'kg', 'M_fuel_kg': 'kg', 'fuel_frac': '',
        'W_S_kgm2': 'kg/m^2', 'b_m': 'm', 'V_s_land': 'm/s', 'V_cruise': 'm/s',
        'C_L_cruise': '', 'L_D_max': '', 'P_cruise_kW': 'kW', 'throttle_frac': '',
        'ROC': 'm/s', 's_TO_m': 'm', 'T0_W': '', 'endurance_h': 'h',
        'range_km': 'km', 'feasible': ''}

GROUPS = [
    ('Range / endurance', ['range_km', 'endurance_h']),
    ('Takeoff & stall', ['s_TO_m', 'V_s_land']),
    ('Climb / power / fuel', ['ROC', 'throttle_frac', 'fuel_frac']),
]


def main():
    OUT.mkdir(exist_ok=True)
    data, meta = {}, {}
    for key, cfg in PROTOTYPES.items():
        ac = Aircraft(cfg['frame'], cfg['pack'], MISSION)
        met = ac.metrics()
        data[key] = {r: met[r] for r in ROWS}
        meta[key] = (cfg['name'], cfg['layout'])
    df = pd.DataFrame(data).T
    df.index.name = 'prototype'
    df.to_csv(OUT / 'prototypes_table.csv')

    print('=' * 100)
    for key in PROTOTYPES:
        print(f"{key}: {meta[key][0]}  [{meta[key][1]}]")
    print('=' * 100)
    print(df.to_string(float_format=lambda x: f'{x:,.2f}'))

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    keys = list(PROTOTYPES)
    for ax, (gtitle, rows) in zip(axes, GROUPS):
        vals = np.array([df.loc[k, rows].astype(float).values for k in keys])
        norm = vals / np.max(vals, axis=0, keepdims=True)
        x = np.arange(len(rows))
        w = 0.2
        for i, k in enumerate(keys):
            ax.bar(x + (i - 1.5) * w, norm[i], w, label=k)
        ax.set_title(gtitle)
        ax.set_xticks(x, [f'{r}\n[{UNIT[r]}]' for r in rows])
        ax.set_ylim(0, 1.15)
        ax.legend()
    fig.tight_layout()
    fig.savefig(OUT / 'prototypes_compare.png', dpi=130)
    print(f'\nsaved: {OUT}/prototypes_table.csv, {OUT}/prototypes_compare.png')


if __name__ == '__main__':
    main()
