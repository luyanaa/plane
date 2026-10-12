"""Layer 7: engine map visualization + operating points + mission comparison.

Answers the design question: single 150cc vs twin 150cc - which has the
longer range/endurance once the 2-D engine map replaces the flat 400 g/kWh
+ throttle-penalty model?

Emits results/engine_map_150cc.png (torque + BSFC surface + op points),
results/engine_operating_points.csv and the mission comparison table.
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from sizing.atmosphere import RHO0
from sizing.configs import MISSION, PROTOTYPES, PACK_150, PACK_200, PACK_2X150
from sizing.engine_map import map_150
from sizing.mission import Aircraft

OUT = Path(__file__).resolve().parents[1] / 'results'


def plot_map(emap, ops, title):
    n = np.linspace(1500, 10500, 200)
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    axes[0].plot(n, [emap.T_max(nn) for nn in n])
    axes[0].axvline(emap.n_r_rpm, ls='--', c='r', lw=0.8, label='rated RPM')
    axes[0].set_xlabel('RPM')
    axes[0].set_ylabel('T_max [N m]')
    axes[0].set_title('torque curve')
    axes[0].legend()
    # BSFC surface vs (n, load)
    nr = np.linspace(2000, 9000, 60)
    lr = np.linspace(0.1, 1.0, 40)
    NR, LR = np.meshgrid(nr, lr)
    BS = np.array([[emap.bsfc(nn, ll) for nn in nr] for ll in lr])
    cf = axes[1].contourf(NR, LR, BS, levels=24, cmap='viridis')
    fig.colorbar(cf, ax=axes[1], label='BSFC [g/kWh]')
    axes[1].set_xlabel('RPM')
    axes[1].set_ylabel('load fraction')
    axes[1].set_title('BSFC surface')
    # operating points overlaid on BSFC
    ax = axes[2]
    for tag, n_, load_, bsfc_, p_kw in ops:
        ax.scatter(n_, load_, s=60, marker='o')
        ax.annotate(f'{tag}\n{bsfc_:.0f} g/kWh', (n_, load_),
                    textcoords='offset points', xytext=(6, 6), fontsize=7)
    ax.set_xlim(2000, 9000)
    ax.set_ylim(0.1, 1.0)
    ax.set_xlabel('RPM')
    ax.set_ylabel('load')
    ax.set_title('operating points (cruise/loiter)')
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(OUT / 'engine_map.png', dpi=120)
    plt.close(fig)


def main():
    ops = []
    rows = []
    for key, cfg in PROTOTYPES.items():
        ac = Aircraft(cfg['frame'], cfg['pack'], MISSION)
        met = ac.metrics(RHO0)
        p = cfg['pack']
        p_cr = met['P_cruise_kW'] / p.n_eng
        # loiter point at the end-of-loiter weight (W1)
        m1 = ac.frame.MTOW_kg - met['fuel_climb_kg'] - (ac.m_fuel - met['fuel_climb_kg']
                                                        - MISSION.reserve_frac * ac.m_fuel)
        v1 = ac.polar.V_LDmax(m1, RHO0)
        p_lo = m1 * 9.80665 * v1 / (ac.polar.L_D_max * p.eta_p_cruise) / 1000.0 / p.n_eng
        op_cr = p.map.operating_point(p_cr)
        op_lo = p.map.operating_point(p_lo)
        ops += [(f'{key}-cr', op_cr['n_rpm'], op_cr['load'], op_cr['bsfc_g_kwh'], p_cr),
                (f'{key}-lo', op_lo['n_rpm'], op_lo['load'], op_lo['bsfc_g_kwh'], p_lo)]
        rows.append({'prototype': key, 'engines': p.n_eng, 'P_cruise_kW': met['P_cruise_kW'],
                     'cr_RPM': op_cr['n_rpm'], 'cr_load': op_cr['load'],
                     'cr_BSFC': op_cr['bsfc_g_kwh'],
                     'lo_RPM': op_lo['n_rpm'], 'lo_load': op_lo['load'],
                     'lo_BSFC': op_lo['bsfc_g_kwh'],
                     'fuel_flow_cr_gh': met['fuel_flow_cruise_gh'],
                     'endurance_h': met['endurance_h'], 'range_km': met['range_km'],
                     'throttle_frac': met['throttle_frac']})
    df = pd.DataFrame(rows).set_index('prototype')
    df.to_csv(OUT / 'engine_operating_points.csv')
    pd.set_option('display.width', 220)
    print(df.to_string(float_format=lambda x: f'{x:8.2f}'))
    plot_map(map_150(), ops, '150cc 2-stroke map (placeholder, bench data pending)')

    # the question: single 150 vs twin 150 (same total power, same payload)
    a = df.loc['A']
    c = df.loc['C']
    print(f"\nSINGLE 150 (A):  endurance {a.endurance_h:.1f} h, range {a.range_km:.0f} km "
          f"@ BSFC {a.cr_BSFC:.0f} g/kWh (load {a.cr_load:.0%})")
    print(f"TWIN 150 (C):    endurance {c.endurance_h:.1f} h, range {c.range_km:.0f} km "
          f"@ BSFC {c.cr_BSFC:.0f} g/kWh (load {c.cr_load:.0%})")
    print(f"ratio A/C endurance: {a.endurance_h / c.endurance_h:.2f}x")
    print(f'\nsaved: {OUT}/engine_map.png, engine_operating_points.csv')


if __name__ == '__main__':
    main()
