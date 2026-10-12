"""Layer 6: propeller actuator-disk model + engine-out yaw (prototypes A-D).

VSPAERO 3.53's Python API does not expose rotor/actuator-disk operating
points (only GUI/CLI), so the disk is modeled here with momentum theory
(which is what an actuator disk IS), consuming the layer-1 engine packs and
the layer-5 VSPAERO baseline (Cn_beta). Quantifies:

  - slipstream velocity/dynamic-pressure ratio at the tail (single) and
    inner wing (twin) at takeoff and cruise
  - elevator-authority boost from propwash (low-speed pitch)
  - thrust-line pitch moment
  - engine-out yaw moment  Cn_EO = T_max*y_eng/(q*S*b)   (T_L=0, T_R=Tmax)
  - rudder authority and V_MC estimate (control-limited)

Emits results/propwash.csv and a comparison table.
"""
import sys
from math import pi, sqrt
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from sizing.atmosphere import RHO0
from sizing.configs import MISSION, PROTOTYPES
from sizing.geometry import build_geometry
from sizing.mission import Aircraft

OUT = Path(__file__).resolve().parents[1] / 'results'

RUDDER_CHORD_FRAC = 0.25
RUDDER_DEFLECTION_DEG = 25.0
RUDDER_TAU = 0.6           # flap effectiveness of rudder on the VTP
VTP_CL_ALPHA = 3.0         # 1/rad, AR~1.8 fin (VLM estimate)
JET_GROWTH_DEG = 4.0       # slipstream cone half-angle
FAR_WAKE_FRAC = 0.9        # tail sits ~90% of the way to the far wake


def thrust_at(V, ac):
    """Full-throttle thrust (N) at speed V: layer-1 linear model."""
    p, m = ac.pack, ac.metrics(RHO0)
    v_cr = m['V_cruise']
    t0 = p.static_thrust(RHO0)
    t_cr = p.eta_p_cruise * p.P_kw * 1000.0 / v_cr
    return t0 + (t_cr - t0) * (V / v_cr) if V < v_cr else p.eta_p_cruise * p.P_kw * 1000.0 / V


def disk_induced(V, T, d_prop):
    """Actuator-disk induced velocity v_i at the disk (momentum theory)."""
    a = pi * (d_prop / 2.0) ** 2
    return (-V + sqrt(V * V + 2.0 * T / (RHO0 * a))) / 2.0


def slipstream(V, T, d_prop, dist):
    """Slipstream velocity and diameter at distance `dist` behind the disk."""
    v_i = disk_induced(V, T, d_prop)
    v_s = V + 2.0 * FAR_WAKE_FRAC * v_i
    d_jet = d_prop * sqrt((V + v_i) / (V + 2.0 * v_i)) + 2.0 * dist * np.tan(
        np.radians(JET_GROWTH_DEG))
    return v_i, v_s, d_jet


def main():
    rows = []
    for key, cfg in PROTOTYPES.items():
        frame, pack = cfg['frame'], cfg['pack']
        ac = Aircraft(frame, pack, MISSION)
        geo = build_geometry(frame.S, frame.AR, frame.taper, frame.t_over_c,
                             frame.l_fus, frame.d_fus, n_nacelles=frame.n_nacelles)
        met = ac.metrics(RHO0)
        v_r = met['V_s_TO'] * 1.1
        v_cr = met['V_cruise']

        r = {'prototype': key, 'S': frame.S, 'b': geo.b, 'MAC': geo.MAC,
             'n_eng': pack.n_eng, 'D_prop': pack.D_prop, 'V_s_TO': met['V_s_TO'],
             'V_R': v_r, 'V_cruise': v_cr}

        # ---- per-engine thrust at V_R and cruise
        t_r = thrust_at(v_r, ac) / pack.n_eng
        t_c = thrust_at(v_cr, ac) / pack.n_eng
        r['T_VR_per_eng'] = t_r
        r['T_cruise_per_eng'] = t_c

        # ---- propwash: tail for singles, inner wing for twins
        if pack.n_eng == 1:
            # central tractor: disk at x~0.55, HTP at x_htp
            dist = (frame.l_fus - 0.35) - 0.55
            b_h = sqrt(geo.S_h * 4.0)
            for tag, V, T in (('TO', v_r, t_r), ('CR', v_cr, t_c)):
                v_i, v_s, d_jet = slipstream(V, T, pack.D_prop, dist)
                frac = min(d_jet / b_h, 1.0)
                q_ratio = (v_s / V) ** 2
                r[f'{tag}_v_i'] = v_i
                r[f'{tag}_q_tail'] = q_ratio
                r[f'{tag}_tail_frac'] = frac
                r[f'{tag}_tail_qboost'] = (q_ratio - 1.0) * frac  # effective q gain
        else:
            # wing-root tractors: jet over inner wing, flap span y~0.45-1.4
            d_flap_in, d_flap_out = 0.45, 1.40
            for tag, V, T in (('TO', v_r, t_r), ('CR', v_cr, t_c)):
                v_i, v_s, d_jet = slipstream(V, T, pack.D_prop, 0.2)
                y_c = 0.75
                cover = (min(y_c + d_jet / 2, d_flap_out)
                         - max(y_c - d_jet / 2, d_flap_in))
                frac = max(cover / (d_flap_out - d_flap_in), 0.0)
                q_ratio = (v_s / V) ** 2
                r[f'{tag}_v_i'] = v_i
                r[f'{tag}_q_wing'] = q_ratio
                r[f'{tag}_flap_frac'] = frac
                r[f'{tag}_flap_qboost'] = (q_ratio - 1.0) * frac

        # ---- thrust-line pitch moment at takeoff (z_T below CG ~0.05 m)
        z_t = 0.05
        q_to = 0.5 * RHO0 * v_r ** 2
        r['dCm_thrust_TO'] = -pack.n_eng * t_r * z_t / (q_to * frame.S * geo.MAC)

        # ---- engine-out yaw (twin only): T_L=0, T_R=Tmax at V_R
        y_eng = 0.75
        if pack.n_eng == 2:
            cn_eo = t_r * y_eng / (q_to * frame.S * geo.b)
            # rudder authority: dCn_r = C_Lalpha_v * tau * (S_v/S)*(l_v/b) * delta
            cn_dr = VTP_CL_ALPHA * RUDDER_TAU * (geo.S_v / frame.S) * (1.15 / geo.b)
            dcn_r = cn_dr * np.radians(RUDDER_DEFLECTION_DEG)
            # control-limited V_MC: T*y = q*S*b*dCn_r
            v_mc = sqrt(2.0 * t_r * y_eng / (RHO0 * frame.S * geo.b * dcn_r))
            r['Cn_EO_VR'] = cn_eo
            r['Cn_dr'] = cn_dr
            r['dCn_r_max'] = dcn_r
            r['V_MC'] = v_mc
            r['V_MC_over_VR'] = v_mc / v_r
            r['V_MC_over_Vs'] = v_mc / met['V_s_land']
        rows.append(r)

    df = pd.DataFrame(rows).set_index('prototype')
    df.to_csv(OUT / 'propwash.csv')
    pd.set_option('display.width', 220)
    show = ['S', 'b', 'V_s_TO', 'V_R', 'V_cruise', 'T_VR_per_eng', 'T_cruise_per_eng']
    for k in df.columns:
        if any(t in k for t in ('q_', 'frac', 'boost', 'v_i')):
            show.append(k)
    show += [c for c in df.columns if c not in show]
    print(df[show].to_string(float_format=lambda x: f'{x:8.2f}'))
    print(f'\nsaved: {OUT}/propwash.csv')


if __name__ == '__main__':
    main()
