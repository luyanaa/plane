"""Layer 5: VSPAERO whole-aircraft aero database + stability derivatives (A-D).

Per prototype (vsp_models/<key>.vsp3 from layer 4):
  - set AnalysisRef: Sref=S, cref=MAC, bref=b, Xcg from component CG
  - alpha sweep -4..16 deg at beta=0   -> C_L(a), C_D(a), C_m(a), C_Lalpha, C_malpha
  - beta sweep 0..8 deg at alpha=4     -> C_Ybeta, C_lbeta, C_nbeta
  - checks: C_malpha<0 (longitudinal), C_nbeta>0 (directional), SM 5-15% MAC

Emits results/vspaero_db.csv, results/vspaero_derivs.csv, results/vspaero_*.png
"""
import multiprocessing as mp
import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import openvsp as vsp

from sizing.configs import MISSION, PROTOTYPES
from sizing.geometry import build_geometry
from sizing.weight import WeightModel

OUT = Path(__file__).resolve().parents[1] / 'results'
MODELS = Path(__file__).resolve().parents[1] / 'vsp_models'
WORK = Path('/tmp/vspaero_l5')

ALPHA_RANGE = (-4.0, 16.0, 11)     # start, end, npts
BETA_RANGE = (0.0, 8.0, 3)         # at alpha = 4
BETA_ALPHA = 4.0
MACH = 0.05
WAKE_ITERS = 10


# ---------------------------------------------------------------- CG model
def cg_x(frame, pack, geo):
    """Component CG from the layer-1 weight breakdown (placeholder, layer 13 refines)."""
    _, items = WeightModel().empty_mass(geo, pack.n_eng, pack.P_kw_per_eng)
    x_wing = 0.33 * frame.l_fus
    b_h = float(np.sqrt(geo.S_h * 4.0))
    c_h = geo.S_h / b_h
    x_htp = frame.l_fus - 0.35
    x_eng = 0.30 if pack.n_eng == 1 else x_wing - 0.05
    pos = {
        'wing': x_wing + 0.30 * geo.MAC,
        'fuselage': 0.45 * frame.l_fus,
        'HTP': x_htp + 0.30 * c_h,
        'VTP': x_htp + 0.30 * c_h,
        'landing gear': 0.55 * frame.l_fus,
        'engine+mount': x_eng,
        'fuel system': 0.40 * frame.l_fus,
        'propellers': x_eng,
        'avionics/FCS': 0.50 * frame.l_fus,
        'servos/wiring': 0.55 * frame.l_fus,
        'payload bay': 0.42 * frame.l_fus,
        'electrical': 0.55 * frame.l_fus,
    }
    m_empty = sum(items.values())
    m_fuel = frame.MTOW_kg - m_empty - MISSION.payload_kg
    x_cg = (sum(items[k] * pos[k] for k in items)
            + m_fuel * 0.42 * frame.l_fus
            + MISSION.payload_kg * 0.42 * frame.l_fus) / frame.MTOW_kg
    return x_cg


# ---------------------------------------------------------------- sweep run
def parse_polar(path):
    rows = []
    for line in Path(path).read_text().splitlines():
        parts = line.split()
        if len(parts) > 24:
            try:
                vals = [float(p) for p in parts[:24]]
            except ValueError:
                continue
            rows.append(vals)
    return rows  # [Beta, Mach, AoA, Re, CLo, CLi, CLtot, CDo, CDi, CDtot,
    #             CSo, CSi, CStot, L/D, E, CMox..z, CMxtot, CMytot, CMztot]


def run_sweep(key):
    model = MODELS / f'{key}.vsp3'
    wdir = WORK / key
    wdir.mkdir(parents=True, exist_ok=True)
    import shutil
    shutil.copy(model, wdir / model.name)
    os.chdir(wdir)

    frame, pack = PROTOTYPES[key]['frame'], PROTOTYPES[key]['pack']
    geo = build_geometry(frame.S, frame.AR, frame.taper, frame.t_over_c,
                         frame.l_fus, frame.d_fus, n_nacelles=frame.n_nacelles)
    x_cg = cg_x(frame, pack, geo)

    vsp.ReadVSPFile(str(wdir / model.name))
    vsp.SetAnalysisInputDefaults('VSPAEROComputeGeometry')
    vsp.ExecAnalysis('VSPAEROComputeGeometry')

    def setd(nm, val):
        if nm in vsp.GetAnalysisInputNames('VSPAEROSweep'):
            vsp.SetDoubleAnalysisInput('VSPAEROSweep', nm, [val])

    def seti(nm, val):
        if nm in vsp.GetAnalysisInputNames('VSPAEROSweep'):
            vsp.SetIntAnalysisInput('VSPAEROSweep', nm, [val])

    vsp.SetAnalysisInputDefaults('VSPAEROSweep')
    # reference values
    setd('Sref', frame.S)
    setd('cref', geo.MAC)
    setd('bref', geo.b)
    setd('Xcg', x_cg)
    setd('Ycg', 0.0)
    setd('Zcg', 0.0)
    seti('UseCGModeFlag', 0)
    setd('MachStart', MACH)
    setd('MachEnd', MACH)
    seti('MachNpts', 1)
    seti('WakeNumIter', WAKE_ITERS)
    setd('Rho', 1.225)

    # alpha sweep
    setd('AlphaStart', ALPHA_RANGE[0])
    setd('AlphaEnd', ALPHA_RANGE[1])
    seti('AlphaNpts', ALPHA_RANGE[2])
    setd('BetaStart', 0.0)
    setd('BetaEnd', 0.0)
    seti('BetaNpts', 1)
    vsp.ExecAnalysis('VSPAEROSweep')
    a_rows = parse_polar(wdir / f'{key}.polar')

    # beta sweep at fixed alpha
    setd('AlphaStart', BETA_ALPHA)
    setd('AlphaEnd', BETA_ALPHA)
    seti('AlphaNpts', 1)
    setd('BetaStart', BETA_RANGE[0])
    setd('BetaEnd', BETA_RANGE[1])
    seti('BetaNpts', BETA_RANGE[2])
    vsp.ExecAnalysis('VSPAEROSweep')
    b_rows = parse_polar(wdir / f'{key}.polar')
    return key, x_cg, a_rows, b_rows


# ---------------------------------------------------------------- derivatives
def derivs(key, x_cg, a_rows, b_rows):
    a = pd.DataFrame(a_rows, columns=['Beta', 'Mach', 'AoA', 'Re', 'CLo', 'CLi',
                                      'CLtot', 'CDo', 'CDi', 'CDtot', 'CSo', 'CSi',
                                      'CStot', 'LD', 'E', 'CMox', 'CMoy', 'CMoz',
                                      'CMix', 'CMiy', 'CMiz', 'CMxtot', 'CMytot', 'CMztot'])
    a = a.sort_values('AoA')
    a = a[a.AoA <= 10.0]                       # linear region
    p = np.polyfit(a.AoA.values, a.CLtot.values, 1)
    q = np.polyfit(a.AoA.values, a.CMytot.values, 1)
    cl_a, cm_a = p[0], q[0]                    # per degree
    cl_a_r, cm_a_r = cl_a * 57.2958, cm_a * 57.2958   # per radian
    sm = -cm_a_r / cl_a_r
    ldmax = (a.CLtot / a.CDtot).max()
    i_ld = (a.CLtot / a.CDtot).idxmax()

    b = pd.DataFrame(b_rows, columns=['Beta', 'Mach', 'AoA', 'Re', 'CLo', 'CLi',
                                      'CLtot', 'CDo', 'CDi', 'CDtot', 'CSo', 'CSi',
                                      'CStot', 'LD', 'E', 'CMox', 'CMoy', 'CMoz',
                                      'CMix', 'CMiy', 'CMiz', 'CMxtot', 'CMytot', 'CMztot'])
    b = b.sort_values('Beta')
    r = np.polyfit(b.Beta.values, b.CStot.values, 1)
    l = np.polyfit(b.Beta.values, b.CMxtot.values, 1)
    n = np.polyfit(b.Beta.values, b.CMztot.values, 1)
    cy_b, cl_b, cn_b = r[0] * 57.2958, l[0] * 57.2958, n[0] * 57.2958  # per radian

    return {'prototype': key, 'x_cg': x_cg, 'CL_alpha': cl_a_r, 'Cm_alpha': cm_a_r,
            'SM': sm, 'L_D_max': ldmax, 'alpha_LDmax': a.loc[i_ld, 'AoA'],
            'CY_beta': cy_b, 'Cl_beta': cl_b, 'Cn_beta': cn_b,
            'Cm_alpha_ok': cm_a_r < 0, 'Cn_beta_ok': cn_b > 0,
            'SM_ok': 0.05 <= sm <= 0.15,
            'a_pts': len(a), 'b_pts': len(b)}


# ---------------------------------------------------------------- plots
def plot_polars(db, keys):
    fig, axes = plt.subplots(4, 3, figsize=(14, 15))
    for i, key in enumerate(keys):
        a = db[key]
        a = a.sort_values('AoA')
        axes[i, 0].plot(a.AoA, a.CLtot, 'o-')
        axes[i, 0].set_ylabel('C_L')
        axes[i, 1].plot(a.AoA, a.CDtot, 'o-')
        axes[i, 1].set_ylabel('C_D')
        axes[i, 2].plot(a.AoA, a.CMytot, 'o-')
        axes[i, 2].set_ylabel('C_m')
        for j in range(3):
            axes[i, j].set_title(f'{key}' if j == 1 else '')
            axes[i, j].set_xlabel('alpha [deg]')
            axes[i, j].grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / 'vspaero_polars.png', dpi=120)
    plt.close(fig)


def plot_derivs(d):
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))
    keys = list(d.keys())
    labels = {'SM': 'static margin [MAC]', 'Cn_beta': 'C_n_beta [1/rad]',
              'Cm_alpha': 'C_m_alpha [1/rad]'}
    for ax, (col, lab) in zip(axes, labels.items()):
        vals = [d[k][col] for k in keys]
        bars = ax.bar(keys, vals, color=['g' if d[k][f'{col}_ok' if col != 'Cm_alpha' else 'Cm_alpha_ok'] else 'r'
                                         for k in keys])
        ax.axhline(0, color='k', lw=0.8)
        ax.set_title(lab)
        for bar, v in zip(bars, vals):
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height(),
                    f'{v:.3f}', ha='center', va='bottom', fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / 'vspaero_derivs.png', dpi=120)
    plt.close(fig)


# ---------------------------------------------------------------- main
def main():
    import os
    os.environ['DISPLAY'] = ':99'
    OUT.mkdir(exist_ok=True)
    WORK.mkdir(exist_ok=True)
    keys = list(PROTOTYPES)
    # one process per prototype (each runs its own openvsp/vspaero)
    with mp.Pool(min(4, len(keys))) as pool:
        results = pool.map(run_sweep, keys)
    db = {k: pd.DataFrame(a, columns=['Beta', 'Mach', 'AoA', 'Re', 'CLo', 'CLi',
                                      'CLtot', 'CDo', 'CDi', 'CDtot', 'CSo', 'CSi',
                                      'CStot', 'LD', 'E', 'CMox', 'CMoy', 'CMoz',
                                      'CMix', 'CMiy', 'CMiz', 'CMxtot', 'CMytot', 'CMztot'])
          for k, _, a, _ in results}
    rows = {k: derivs(k, x, a, b) for k, x, a, b in results}
    pd.concat([pd.DataFrame([v]) for v in rows.values()], ignore_index=True).to_csv(
        OUT / 'vspaero_derivs.csv', index=False)
    pd.concat([db[k].assign(prototype=k) for k in keys], ignore_index=True).to_csv(
        OUT / 'vspaero_db.csv', index=False)
    plot_polars(db, keys)
    plot_derivs(rows)
    d = pd.DataFrame(rows).T
    pd.set_option('display.width', 200)
    print(d[['x_cg', 'CL_alpha', 'Cm_alpha', 'SM', 'L_D_max', 'alpha_LDmax',
             'CY_beta', 'Cl_beta', 'Cn_beta', 'Cm_alpha_ok', 'Cn_beta_ok', 'SM_ok',
             'a_pts', 'b_pts']].to_string(float_format=lambda x: f'{x:8.3f}'))
    print(f'\nsaved: {OUT}/vspaero_db.csv, vspaero_derivs.csv, vspaero_*.png')


if __name__ == '__main__':
    main()
