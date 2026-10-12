"""Layer 3: XFOIL plain-flap 2-D sweep -> first aerodynamic flap lookup table.

Sweeps flap chord (20-35%) x deflection (0-40 deg) x transition mode at
Re = 8e5 for each candidate airfoil. The flapped geometry is built in Python
(rotate the part x >= 1 - flap_chord about the hinge by +delta, TE down),
then run through XFOIL with the same batch/VACC/Xvfb machinery as layer 2.

Emits:
  results/flap_polars.csv      raw polars
  results/flap_summary.csv     C_Lmax, dC_Lmax, alpha_stall, CD@CL0.8, Cm per config
  results/flap_lookup.csv      three states: cruise 0 / takeoff 10 / landing 30
  results/flap_<airfoil>.png   dC_Lmax vs delta + C_L(alpha) families

Usage: .venv/bin/python scripts/flap_scan.py [xfoil_binary]
"""
import os
import subprocess
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor
from math import cos, radians, sin
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts.airfoil_scan import (WORK as L2_WORK, DISPLAY_NUM, ensure_airfoils,
                                  ensure_display, parse_polar)

OUT = Path(__file__).resolve().parents[1] / 'results'
WORK = Path(tempfile.gettempdir()) / 'flap_scan'

AIRFOILS = ['naca2412', 'sd7062']      # A (endurance) and D (STOL) candidates
FLAP_CHORDS = [0.20, 0.25, 0.30, 0.35]
DEFLECTIONS = [0, 5, 10, 15, 20, 25, 30, 35, 40]
RE = 8e5
MODES = {'free': None, 'forced5': 0.05}
ALPHAS = list(range(-6, 17))           # -6..16 deg

# three flight states for the lookup table
STATES = [('cruise', 0.0), ('takeoff', 10.0), ('landing', 30.0)]


def run_flap_polar(xfoil, src_dat, hinge_x, delta, re, xtr, polar_path, log_path):
    """XFOIL built-in GDES-FLAP (deflect TE flap) + polar, VACC retry.

    GDES exit is a blank line; 'eXec' writes the buffer back to current.
    """
    lines = [f'LOAD {src_dat}']
    if delta > 0:
        lines += ['GDES', 'Flap', f'{hinge_x}', '0', f'{delta}', 'eXec', '']
    lines += ['PANE', 'OPER', 'ITER 300', f'VISC {re:.0f}', 'M 0.0']
    if xtr is not None:
        lines += ['VPAR', 'XTR', f'{xtr}', f'{xtr}', '']
    lines += ['PACC', str(polar_path), '']
    lines += [f'ALFA {a}' for a in ALPHAS]
    lines += ['PCLR', 'QUIT']
    env = {**os.environ, 'DISPLAY': os.environ.get('DISPLAY', f':{DISPLAY_NUM}')}
    for attempt, extra in enumerate(([], ['VPAR', 'VACC', '0.3', ''])):
        batch = lines[:]
        if extra:
            idx = batch.index('M 0.0') + 1
            batch[idx:idx] = extra
        try:
            with open(log_path, 'w') as lf:
                subprocess.run([xfoil], input='\n'.join(batch) + '\n', text=True,
                               stdout=lf, stderr=subprocess.STDOUT, timeout=45, env=env)
            break
        except subprocess.TimeoutExpired:
            if attempt == 0:
                continue
    return polar_path


def one_job(args):
    xfoil, name, chord, delta, mode, xtr, src_dat = args
    tag = f'{name}_c{chord:.2f}_d{delta:02d}_{mode}'
    pol = WORK / f'{tag}.pol'
    log = WORK / f'{tag}.log'
    pol.unlink(missing_ok=True)
    run_flap_polar(xfoil, src_dat, 1.0 - chord, delta, RE, xtr, pol, log)
    rows = parse_polar(pol)
    return [{'airfoil': name, 'flap_chord': chord, 'delta': delta, 'mode': mode,
             'alpha': r[0], 'CL': r[1], 'CD': r[2], 'CDp': r[3], 'Cm': r[4],
             'TopXtr': r[5], 'BotXtr': r[6]} for r in rows]


def summarize(df):
    rows = []
    for (name, chord, delta, mode), g in df.groupby(['airfoil', 'flap_chord', 'delta', 'mode']):
        g = g.sort_values('alpha')
        if len(g) < 4:
            continue
        g = g[g.CD >= 0.003]                      # below laminar friction floor = bad
        if len(g) < 4:
            continue
        imax = g.CL.idxmax()
        clmax = g.CL.max()
        cd_cruise = np.interp(0.8, g.CL.values, g.CD.values) if g.CL.min() < 0.8 < g.CL.max() else np.nan
        cm_cruise = np.interp(0.8, g.CL.values, g.Cm.values) if g.CL.min() < 0.8 < g.CL.max() else np.nan
        rows.append({'airfoil': name, 'flap_chord': chord, 'delta': delta, 'mode': mode,
                     'C_Lmax': clmax, 'alpha_stall': g.loc[imax, 'alpha'],
                     'CD_at_CL08': cd_cruise, 'Cm_at_CL08': cm_cruise,
                     'n_points': len(g)})
    s = pd.DataFrame(rows)
    # dC_Lmax vs the clean (delta=0) case, per airfoil/chord/mode
    clean = s[s.delta == 0][['airfoil', 'flap_chord', 'mode', 'C_Lmax']]
    clean = clean.rename(columns={'C_Lmax': 'C_Lmax_clean'})
    s = s.merge(clean, on=['airfoil', 'flap_chord', 'mode'], how='left')
    s['dC_Lmax'] = s.C_Lmax - s.C_Lmax_clean
    return s


def lookup_table(summ):
    rows = []
    for name in AIRFOILS:
        for mode in MODES:
            for state, delta in STATES:
                row = summ[(summ.airfoil == name) & (summ['mode'] == mode)
                           & (summ.delta == delta)]
                if len(row):
                    # pick the chord with max C_Lmax for TO/landing, min CD@CL0.8 for cruise
                    if state == 'cruise':
                        best = row.loc[row.CD_at_CL08.idxmin()]
                    else:
                        best = row.loc[row.C_Lmax.idxmax()]
                    rows.append({'airfoil': name, 'mode': mode, 'state': state,
                                 'delta': delta, 'flap_chord': best.flap_chord,
                                 'C_Lmax': best.C_Lmax, 'alpha_stall': best.alpha_stall,
                                 'CD_at_CL08': best.CD_at_CL08, 'Cm_at_CL08': best.Cm_at_CL08})
    return pd.DataFrame(rows)


def plot_airfoil(df, summ, name):
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
    sub = summ[(summ.airfoil == name) & (summ['mode'] == 'free')]
    for chord in FLAP_CHORDS:
        g = sub[sub.flap_chord == chord].sort_values('delta')
        axes[0].plot(g.delta, g.dC_Lmax, marker='o', label=f'flap {chord:.0%}')
    axes[0].set_title(f'{name}: dC_Lmax vs deflection (free transition, Re={RE:.0e})')
    axes[0].set_xlabel('delta [deg]')
    axes[0].set_ylabel('dC_Lmax')
    axes[0].grid(alpha=0.3)
    axes[0].legend()
    for delta in (0, 10, 20, 30, 40):
        g = df[(df.airfoil == name) & (df['mode'] == 'free')
               & (df.flap_chord == 0.25) & (df.delta == delta)].sort_values('alpha')
        axes[1].plot(g.alpha, g.CL, label=f'delta={delta}')
    axes[1].set_title(f'{name}: C_L(alpha), 25% flap, free transition')
    axes[1].set_xlabel('alpha [deg]')
    axes[1].set_ylabel('C_L')
    axes[1].grid(alpha=0.3)
    axes[1].legend()
    fig.tight_layout()
    fig.savefig(OUT / f'flap_{name}.png', dpi=120)
    plt.close(fig)


def main():
    xfoil = sys.argv[1] if len(sys.argv) > 1 else 'xfoil'
    os.environ['DISPLAY'] = ensure_display()
    OUT.mkdir(exist_ok=True)
    WORK.mkdir(exist_ok=True)
    names = ensure_airfoils()          # populates layer-2 WORK with source dats
    jobs = [(xfoil, n, c, d, mode, xtr, str(L2_WORK / f'{n}.dat'))
            for n in AIRFOILS if n in names
            for c in FLAP_CHORDS for d in DEFLECTIONS for mode, xtr in MODES.items()]
    print(f'{len(jobs)} flap runs: {len(AIRFOILS)} airfoils x {len(FLAP_CHORDS)} chords '
          f'x {len(DEFLECTIONS)} deflections x {len(MODES)} modes @ Re {RE:.0e}')
    all_rows = []
    with ProcessPoolExecutor(max_workers=min(24, 48)) as ex:
        for i, rows in enumerate(ex.map(one_job, jobs)):
            all_rows += rows
            if (i + 1) % 50 == 0:
                print(f'  {i + 1}/{len(jobs)} done', flush=True)
    df = pd.DataFrame(all_rows)
    df.to_csv(OUT / 'flap_polars.csv', index=False)
    summ = summarize(df)
    summ.to_csv(OUT / 'flap_summary.csv', index=False)
    look = lookup_table(summ)
    look.to_csv(OUT / 'flap_lookup.csv', index=False)
    for name in AIRFOILS:
        plot_airfoil(df, summ, name)
    pd.set_option('display.width', 220)
    print('\n=== three-state lookup (chord chosen per state) ===')
    print(look.to_string(index=False, float_format=lambda x: f'{x:.3f}'))
    print('\n=== dC_Lmax summary (free) ===')
    f = summ[summ['mode'] == 'free'].pivot_table(index=['airfoil', 'flap_chord'],
                                                 columns='delta', values='dC_Lmax')
    print(f.round(3).to_string())
    print(f'\nsaved: {OUT}/flap_polars.csv, flap_summary.csv, flap_lookup.csv, flap_*.png')


if __name__ == '__main__':
    main()
