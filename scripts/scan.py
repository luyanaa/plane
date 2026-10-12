"""Design-space sweep: S x AR x MTOW for each engine pack (layer 1).

Saves results/scan_<pack>.csv + results/scan_all.csv, prints per-pack
feasibility and top-5 by range / endurance / short takeoff, and writes
results/scan_pareto.png and results/scan_feasibility_mtow100.png.
"""
import itertools
import sys
from dataclasses import replace
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from sizing.configs import MISSION, PACKS
from sizing.mission import Aircraft, Airframe

OUT = Path(__file__).resolve().parents[1] / 'results'

S_GRID = np.arange(2.5, 4.51, 0.1)
AR_GRID = np.arange(8.0, 14.01, 0.5)
MTOW_GRID = np.arange(70.0, 140.01, 5.0)

COLS = ['b_m', 'W_S_kgm2', 'M_fuel_kg', 'fuel_frac', 'V_s_land', 'V_cruise',
        'C_L_cruise', 'L_D_max', 'P_cruise_kW', 'throttle_frac', 'ROC',
        's_TO_m', 'T0_W', 'range_km', 'endurance_h', 'feasible']


def main():
    OUT.mkdir(exist_ok=True)
    frames = []
    for pack in PACKS:
        pack = replace(pack, map=None, prop=None)  # scan uses the generic model (fast)
        rows = []
        for s, ar, m in itertools.product(S_GRID, AR_GRID, MTOW_GRID):
            ac = Aircraft(Airframe(S=s, AR=ar, MTOW_kg=m), pack, MISSION)
            met = ac.metrics()
            rows.append({'pack': pack.name, 'S': s, 'AR': ar, 'MTOW': m,
                         **{k: met[k] for k in COLS}})
        df = pd.DataFrame(rows)
        df.to_csv(OUT / f'scan_{pack.name}.csv', index=False)
        frames.append(df)

        feas = df[df.feasible]
        print(f'\n=== pack {pack.name}: {len(feas)}/{len(df)} feasible '
              f'({100.0 * len(feas) / len(df):.0f}%) ===')
        if len(feas):
            t5 = feas.nlargest(5, 'range_km')
            print('top-5 range (km):      ',
                  t5[['S', 'AR', 'MTOW', 'range_km', 'endurance_h']].round(1).values.tolist())
            t5 = feas.nlargest(5, 'endurance_h')
            print('top-5 endurance (h):   ',
                  t5[['S', 'AR', 'MTOW', 'endurance_h', 'range_km']].round(1).values.tolist())
            t5 = feas.nsmallest(5, 's_TO_m')
            print('top-5 short TO (m):    ',
                  t5[['S', 'AR', 'MTOW', 's_TO_m', 'ROC', 'V_s_land']].round(1).values.tolist())

    all_df = pd.concat(frames, ignore_index=True)
    all_df.to_csv(OUT / 'scan_all.csv', index=False)
    print(f'\nsaved: {OUT}/scan_all.csv + scan_<pack>.csv')

    _plot_pareto()
    _plot_feasibility()


def _plot_pareto():
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    for ax, pack in zip(axes.ravel(), PACKS):
        df = pd.read_csv(OUT / f'scan_{pack.name}.csv')
        feas = df[df.feasible]
        sc = ax.scatter(feas.s_TO_m, feas.range_km, c=feas.MTOW, s=8, cmap='viridis')
        ax.set_title(pack.name)
        ax.set_xlabel('s_TO [m]')
        ax.set_ylabel('Range [km]')
        ax.axvline(MISSION.s_TO_max, ls='--', c='r', lw=0.8)
        fig.colorbar(sc, ax=ax, label='MTOW [kg]')
    fig.tight_layout()
    fig.savefig(OUT / 'scan_pareto.png', dpi=130)
    print(f'saved: {OUT}/scan_pareto.png')


def _plot_feasibility():
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    for ax, pack in zip(axes.ravel(), PACKS):
        df = pd.read_csv(OUT / f'scan_{pack.name}.csv')
        sub = df[df.MTOW == 100.0]
        piv = sub.pivot(index='S', columns='AR', values='feasible')
        ax.pcolormesh(piv.columns.values, piv.index.values, piv.values,
                      cmap='RdYlGn', vmin=0, vmax=1)
        ax.set_title(f'{pack.name} @ MTOW=100 kg (green=feasible)')
        ax.set_xlabel('AR')
        ax.set_ylabel('S [m^2]')
    fig.tight_layout()
    fig.savefig(OUT / 'scan_feasibility_mtow100.png', dpi=130)
    print(f'saved: {OUT}/scan_feasibility_mtow100.png')


if __name__ == '__main__':
    main()
