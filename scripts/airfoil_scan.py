"""Layer 2: XFOIL airfoil polar sweep at UAV Reynolds numbers.

Per candidate airfoil, sweep Re = {0.4, 0.6, 0.8, 1.0, 1.5}e6 x transition
mode {free, forced x/c=0.05, forced x/c=0.10}. Emits:
  results/airfoil_polars.csv            all raw polars
  results/airfoil_summary.csv           C_Lmax, min CD, max L/D, cruise plateau
  results/airfoil_<name>.png            3x3 panel figure per airfoil

Usage:  .venv/bin/python scripts/airfoil_scan.py [xfoil_binary]
"""
import csv
import os
import subprocess
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor
from math import cos, pi, sin, sqrt
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

OUT = Path(__file__).resolve().parents[1] / 'results'
WORK = Path(tempfile.gettempdir()) / 'xfoil_scan'

DISPLAY_NUM = 99
XVFB_CANDIDATES = [
    os.environ.get('XVFB', ''),
    '/tmp/xvfb_local/usr/bin/Xvfb',
    '/usr/bin/Xvfb',
]

RES = [0.4e6, 0.6e6, 0.8e6, 1.0e6, 1.5e6]
MODES = {'free': None, 'forced5': 0.05, 'forced10': 0.10}
ALPHAS = list(range(-6, 17))            # -6..16 deg

NACA4 = ['naca2412', 'naca4412', 'naca2415']          # generated locally
UIUC = ['sd7037', 'sd7062', 'e197', 'ag36', 'fx63137']  # from data/airfoils/ (synced)
AIRFOIL_DIR = Path(__file__).resolve().parents[1] / 'data' / 'airfoils'


# ---------------------------------------------------------------- airfoils
def naca4_coords(name, n=120):
    """Analytic NACA 4-digit coordinates (cosine spacing)."""
    m = int(name[4]) / 100.0
    p = int(name[5]) / 10.0
    t = int(name[6:8]) / 100.0
    x = 0.5 * (1.0 - np.cos(np.linspace(0.0, pi, n)))
    yt = 5 * t * (0.2969 * np.sqrt(x) - 0.1260 * x - 0.3516 * x ** 2
                  + 0.2843 * x ** 3 - 0.1015 * x ** 4)
    yc = np.where(x < p, m / p ** 2 * (2 * p * x - x ** 2),
                  m / (1 - p) ** 2 * ((1 - 2 * p) + 2 * p * x - x ** 2))
    dyc = np.where(x < p, 2 * m / p ** 2 * (p - x), 2 * m / (1 - p) ** 2 * (p - x))
    th = np.arctan(dyc)
    xu, yu = x - yt * np.sin(th), yc + yt * np.cos(th)
    xl, yl = x + yt * np.sin(th), yc - yt * np.cos(th)
    coords = np.vstack([np.hstack([xu[::-1], xl[1:]]), np.hstack([yu[::-1], yl[1:]])]).T
    return coords


def fetch_uiuc(name, dest):
    """Best-effort download from the UIUC airfoil database."""
    url = f'https://m-selig.ae.illinois.edu/ads/coord/{name}.dat'
    try:
        import urllib.request
        data = urllib.request.urlopen(url, timeout=20).read().decode()
        lines = [l for l in data.splitlines() if l.strip() and not l.startswith(('#', '%'))]
        pts = [tuple(map(float, l.split())) for l in lines]
        if len(pts) > 10:
            with open(dest, 'w') as fh:
                fh.write(f'{name}\n' + '\n'.join(f'{x: .6f} {y: .6f}' for x, y in pts) + '\n')
            return True
    except Exception:
        pass
    return False


def ensure_airfoils():
    WORK.mkdir(parents=True, exist_ok=True)
    ok = []
    for name in NACA4:
        path = WORK / f'{name}.dat'
        np.savetxt(path, naca4_coords(name), fmt='% .6f', header=name, comments='')
        ok.append(name)
    for name in UIUC:
        src = AIRFOIL_DIR / f'{name}.dat'
        dest = WORK / f'{name}.dat'
        if src.exists():
            import shutil
            shutil.copy(src, dest)
            ok.append(name)
        elif dest.exists() or fetch_uiuc(name, dest):
            ok.append(name)
        else:
            print(f'[warn] airfoil {name} unavailable, skipping')
    return ok


# ---------------------------------------------------------------- xfoil run
def run_polar(xfoil, dat, re, mode, xtr, alpha_list, polar_path, log_path):
    cmd = [xfoil]
    lines = [f'LOAD {dat}', 'PANE', 'OPER', 'ITER 300', f'VISC {re:.0f}', 'M 0.0']
    if xtr is not None:
        lines += ['VPAR', 'XTR', f'{xtr}', f'{xtr}', '']
    lines += ['PACC', str(polar_path), '']
    lines += [f'ALFA {a}' for a in alpha_list]
    lines += ['PCLR', 'QUIT']
    env = {**os.environ, 'DISPLAY': os.environ.get('DISPLAY', f':{DISPLAY_NUM}')}
    for attempt, extra in enumerate(([], ['VPAR', 'VACC', '0.3', ''])):
        batch = lines[:]
        if extra:
            # viscous acceleration helps BL oscillation (SD70xx low-Re)
            idx = batch.index('M 0.0') + 1
            batch[idx:idx] = extra
        try:
            with open(log_path, 'w') as lf:
                subprocess.run(cmd, input='\n'.join(batch) + '\n', text=True,
                               stdout=lf, stderr=subprocess.STDOUT, timeout=45, env=env)
            break
        except subprocess.TimeoutExpired:
            if attempt == 0:
                continue  # retry with VACC
            # keep whatever polar points were already written and move on
    return polar_path


def parse_polar(path):
    rows = []
    if not Path(path).exists():
        return rows
    with open(path) as fh:
        lines = fh.read().splitlines()
    for i, line in enumerate(lines):
        if line.split() and line.split()[0].lower() == 'alpha':
            for data in lines[i + 2:]:
                if not data.strip():
                    break
                try:
                    vals = [float(v) for v in data.split()]
                    if len(vals) >= 7:
                        rows.append(vals[:7])
                except ValueError:
                    break
            break
    return rows  # [alpha, CL, CD, CDp, Cm, TopXtr, BotXtr]


def one_job(args):
    xfoil, name, re, mode, xtr, dat = args
    tag = f'{name}_Re{re:.1e}_{mode}'
    polar = WORK / f'{tag}.pol'
    log = WORK / f'{tag}.log'
    polar.unlink(missing_ok=True)
    run_polar(xfoil, dat, re, mode, xtr, ALPHAS, polar, log)
    rows = parse_polar(polar)
    return [{'airfoil': name, 'Re': re, 'mode': mode, 'alpha': r[0], 'CL': r[1],
             'CD': r[2], 'CDp': r[3], 'Cm': r[4], 'TopXtr': r[5], 'BotXtr': r[6]}
            for r in rows]


# ---------------------------------------------------------------- summary
def summarize(df):
    rows = []
    for (name, re, mode), g in df.groupby(['airfoil', 'Re', 'mode']):
        g = g.sort_values('alpha')
        if len(g) < 5:
            continue
        g = g[g.CD >= 0.003]  # below laminar skin-friction floor = bad point
        if len(g) < 5:
            continue
        g['LD'] = g.CL / g.CD
        clmax = g.CL.max()
        cdmin = g.CD.min()
        ldmax = g.LD.max()
        row_ldmax = g.loc[g.LD.idxmax()]
        ld90 = 0.9 * ldmax
        plateau = g[(g.LD >= ld90) & (g.CL.between(0.4, 0.9))]
        rows.append({'airfoil': name, 'Re': re, 'mode': mode,
                     'C_Lmax': clmax, 'CD_min': cdmin, 'LD_max': ldmax,
                     'CL_at_LDmax': row_ldmax.CL, 'alpha_at_LDmax': row_ldmax.alpha,
                     'plateau_CL_min': plateau.CL.min() if len(plateau) else np.nan,
                     'plateau_CL_max': plateau.CL.max() if len(plateau) else np.nan,
                     'plateau_width': (plateau.CL.max() - plateau.CL.min())
                                      if len(plateau) else 0.0})
    return pd.DataFrame(rows)


def plot_airfoil(df, name):
    fig, axes = plt.subplots(2, 3, figsize=(15, 9))
    for ax, re in zip(axes.ravel(), RES):
        g = df[(df.airfoil == name) & (df.Re == re)]
        for mode, sub in g.groupby('mode'):
            sub = sub.sort_values('alpha')
            ax.plot(sub.alpha, sub.CL, label=mode)
            ax.set_title(f'{name}  Re={re:.1e}')
            ax.set_xlabel('alpha [deg]')
            ax.set_ylabel('C_L')
            ax.legend(fontsize=7)
    for ax, re in zip(axes.ravel()[3:], RES):
        g = df[(df.airfoil == name) & (df.Re == re)]
        for mode, sub in g.groupby('mode'):
            sub = sub.sort_values('CL')
            ax.plot(sub.CL, sub.CL / sub.CD, label=mode)
            ax.set_xlabel('C_L')
            ax.set_ylabel('L/D')
            ax.axvspan(0.4, 0.9, alpha=0.15, color='g', label='cruise plateau')
            ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(OUT / f'airfoil_{name}.png', dpi=120)
    plt.close(fig)


def ensure_display():
    """Ubuntu xfoil needs an X server even in batch; start Xvfb if needed."""
    if os.environ.get('DISPLAY'):
        return os.environ['DISPLAY']
    import socket
    import time
    x_socket = Path(f'/tmp/.X11-unix/X{DISPLAY_NUM}')
    if x_socket.exists():
        return f':{DISPLAY_NUM}'
    for cand in XVFB_CANDIDATES:
        if cand and Path(cand).exists():
            subprocess.Popen([cand, f':{DISPLAY_NUM}', '-screen', '0', '1024x768x24'],
                             stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            print(f'[xfoil] started {cand} on :{DISPLAY_NUM}')
            for _ in range(30):
                time.sleep(0.2)
                if x_socket.exists():
                    return f':{DISPLAY_NUM}'
            print('[xfoil] WARNING: Xvfb did not become ready')
            return f':{DISPLAY_NUM}'
    print('[xfoil] WARNING: no Xvfb found; xfoil batch runs will abort '
          '(Ubuntu build requires a display). Set DISPLAY or XVFB.')
    return f':{DISPLAY_NUM}'


# ---------------------------------------------------------------- main
def main():
    xfoil = sys.argv[1] if len(sys.argv) > 1 else 'xfoil'
    os.environ['DISPLAY'] = ensure_display()
    OUT.mkdir(exist_ok=True)
    names = ensure_airfoils()
    jobs = [(xfoil, n, re, mode, xtr, str(WORK / f'{n}.dat'))
            for n in names for re in RES for mode, xtr in MODES.items()]
    print(f'{len(jobs)} XFOIL runs on {len(names)} airfoils '
          f'(Re {RES[0]:.0e}-{RES[-1]:.0e}, modes {list(MODES)})')
    all_rows = []
    with ProcessPoolExecutor(max_workers=min(24, 48)) as ex:
        for i, rows in enumerate(ex.map(one_job, jobs)):
            all_rows += rows
            if (i + 1) % 100 == 0:
                print(f'  {i + 1}/{len(jobs)} done', flush=True)
    df = pd.DataFrame(all_rows)
    df.to_csv(OUT / 'airfoil_polars.csv', index=False)
    summ = summarize(df)
    summ.to_csv(OUT / 'airfoil_summary.csv', index=False)
    for name in names:
        plot_airfoil(df, name)
    pd.set_option('display.width', 220)
    print(summ.sort_values(['airfoil', 'Re', 'mode']).to_string(index=False,
          float_format=lambda x: f'{x:.3f}'))
    print(f'\nsaved: {OUT}/airfoil_polars.csv, airfoil_summary.csv, airfoil_*.png')


if __name__ == '__main__':
    main()
