"""Layer 9: JSBSim 6-DoF models and analysis for prototypes A-D.

Builds JSBSim aircraft XMLs (aero tables from the VSPAERO DB, analytic
control derivatives, mass/inertia from layer-1 weights, piston engine from
the layer-7 map, propeller CT/CP tables from the layer-8 BEMT), then:
  - trims: level cruise / climb (pitch) / approach
  - dynamic modes: linearized state matrix eigen-decomposition
    (phugoid, short period, Dutch roll, roll subsidence, spiral)
  - twin engine-out: left engine off at a speed sweep with max rudder
    -> yaw rate, bank, altitude loss, control-limited V_MC

Emits results/jsbsim_trim.csv, jsbsim_modes.csv, jsbsim_vmc.csv.
"""
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import jsbsim

from sizing.atmosphere import G, RHO0
from sizing.configs import MISSION, PROTOTYPES
from sizing.geometry import build_geometry
from sizing.mission import Aircraft
from sizing.weight import WeightModel

OUT = Path(__file__).resolve().parents[1] / 'results'
JROOT = Path(__file__).resolve().parents[1] / 'jsbsim_root'

ELEV_MAX_DEG, AIL_MAX_DEG, RUD_MAX_DEG = 25.0, 20.0, 25.0
TAU_ELEV, TAU_RUD = 0.5, 0.6
CL_ALPHA_H, CL_ALPHA_V = 4.5, 3.0
BSFC_LB_HPHR = 0.66          # 400 g/kWh


# ---------------------------------------------------------------- helpers
def cg_x(frame, pack, geo):
    """Same component CG as layer 5."""
    _, items = WeightModel().empty_mass(geo, pack.n_eng, pack.P_kw_per_eng)
    x_wing = 0.33 * frame.l_fus
    b_h = float(np.sqrt(geo.S_h * 4.0))
    c_h = geo.S_h / b_h
    x_htp = frame.l_fus - 0.35
    x_eng = 0.30 if pack.n_eng == 1 else x_wing - 0.05
    pos = {'wing': x_wing + 0.30 * geo.MAC, 'fuselage': 0.45 * frame.l_fus,
           'HTP': x_htp + 0.30 * c_h, 'VTP': x_htp + 0.30 * c_h,
           'landing gear': 0.55 * frame.l_fus, 'engine+mount': x_eng,
           'fuel system': 0.40 * frame.l_fus, 'propellers': x_eng,
           'avionics/FCS': 0.50 * frame.l_fus, 'servos/wiring': 0.55 * frame.l_fus,
           'payload bay': 0.42 * frame.l_fus, 'electrical': 0.55 * frame.l_fus}
    m_e = sum(items.values())
    m_f = frame.MTOW_kg - m_e - MISSION.payload_kg
    return (sum(items[k] * pos[k] for k in items) + (m_f + MISSION.payload_kg)
            * 0.42 * frame.l_fus) / frame.MTOW_kg


def inertias(frame, pack, geo):
    m = frame.MTOW_kg
    _, items = WeightModel().empty_mass(geo, pack.n_eng, pack.P_kw_per_eng)
    ixx = items['wing'] * (geo.b / 2) ** 2 / 3.0 + 2.0
    iyy = 0.10 * m * frame.l_fus ** 2
    izz = ixx + iyy
    return ixx, iyy, izz


def control_derivs(frame, geo):
    sh_s = geo.S_h / frame.S
    sv_s = geo.S_v / frame.S
    cmde = -sh_s * (1.15 / geo.MAC) * CL_ALPHA_H * TAU_ELEV
    clde = sh_s * CL_ALPHA_H * TAU_ELEV
    cndr = sv_s * (1.15 / geo.b) * CL_ALPHA_V * TAU_RUD
    cydr = sv_s * CL_ALPHA_V * TAU_RUD
    return {'cmde': cmde, 'clde': clde, 'cndr': cndr, 'cydr': cydr,
            'clda': -0.15, 'cnda': 0.01}


def vsp_aero(key):
    """CL/CD/Cm vs alpha (rad), CD shifted to the layer-1 profile drag."""
    db = pd.read_csv(OUT / 'vspaero_db.csv')
    d = db[db.prototype == key].sort_values('AoA')
    a_deg = d.AoA.values
    cl = d.CLtot.values
    cd = d.CDtot.values
    cm = d.CMytot.values
    # extend to +-20 deg with the linear slope
    p_cl = np.polyfit(a_deg, cl, 1)
    p_cm = np.polyfit(a_deg, cm, 1)
    a_ext = np.concatenate([[-20.0], a_deg, [20.0]])
    cl = np.concatenate([[p_cl[0] * -20 + p_cl[1]], cl, [p_cl[0] * 20 + p_cl[1]]])
    cm = np.concatenate([[p_cm[0] * -20 + p_cm[1]], cm, [p_cm[0] * 20 + p_cm[1]]])
    cd0_shift = layer1_cd0(key) - cd[0]
    cd = np.concatenate([[cd[0] + cd0_shift], cd + cd0_shift])
    return np.radians(a_ext), cl, cd, cm


def layer1_cd0(key):
    from sizing.configs import PROTOTYPES as P
    f, p = P[key]['frame'], P[key]['pack']
    g = build_geometry(f.S, f.AR, f.taper, f.t_over_c, f.l_fus, f.d_fus,
                       n_nacelles=f.n_nacelles)
    return 0.009 * g.S_wet / f.S


def lateral_derivs(key):
    dr = pd.read_csv(OUT / 'vspaero_derivs.csv')
    r = dr[dr.prototype == key].iloc[0]
    return r.CY_beta, r.Cl_beta, r.Cn_beta


def table_xml(x, y, indent='          '):
    pairs = ' '.join(f'{a:.5f} {b:.6f}' for a, b in zip(x, y))
    return f'{indent}<tableData>{pairs}</tableData>\n'


# ---------------------------------------------------------------- builders
def build_engine(key, pack):
    emap = pack.map
    n = np.linspace(1800, 10500, 9)
    p = [emap.P_max_kw(nn) * 1.341 for nn in n]      # kW -> hp
    rows = ' '.join(f'{nn:.0f} {pp:.2f}' for nn, pp in zip(n, p))
    xml = f'''<?xml version="1.0"?>
<piston_engine name="eng{key}">
  <minmp unit="INHG">10.0</minmp>
  <maxmp unit="INHG">28.5</maxmp>
  <maxhp>{pack.P_kw * 1.341:.1f}</maxhp>
  <bsfc>{BSFC_LB_HPHR}</bsfc>
  <displacement unit="IN3">{150 if pack.P_kw_per_eng <= 11 else 222:.0f}</displacement>
  <cycles>2.0</cycles>
  <idlerpm>1800</idlerpm>
  <maxrpm>10500</maxrpm>
  <maxthrottle>1.0</maxthrottle>
  <minthrottle>0.05</minthrottle>
</piston_engine>
'''
    (JROOT / 'engine').mkdir(parents=True, exist_ok=True)
    (JROOT / 'engine' / f'eng{key}.xml').write_text(xml)


def build_prop(key, pack):
    js, ct, cp, _ = pack.prop.curves()
    ct_tab = table_xml(js, ct)
    cp_tab = table_xml(js, cp)
    d_in = pack.D_prop * 39.3701
    pitch_in = pack.prop.P_over_D * d_in
    xml = f'''<?xml version="1.0"?>
<propeller name="prop{key}">
  <ixx>0.05</ixx>
  <diameter unit="IN">{d_in:.2f}</diameter>
  <numblades>{pack.prop.B}</numblades>
  <minpitch>{pitch_in:.2f}</minpitch>
  <maxpitch>{pitch_in:.2f}</maxpitch>
  <table name="C_THRUST" type="internal">
{ct_tab}  </table>
  <table name="C_POWER" type="internal">
{cp_tab}  </table>
</propeller>
'''
    (JROOT / 'engine' / f'prop{key}.xml').write_text(xml)


def build_aircraft(key, cfg):
    frame, pack = cfg['frame'], cfg['pack']
    geo = build_geometry(frame.S, frame.AR, frame.taper, frame.t_over_c,
                         frame.l_fus, frame.d_fus, n_nacelles=frame.n_nacelles)
    x_cg = cg_x(frame, pack, geo)
    ixx, iyy, izz = inertias(frame, pack, geo)
    cd = control_derivs(frame, geo)
    a_rad, cl, cdd, cm = vsp_aero(key)
    cyb, clb, cnb = lateral_derivs(key)
    m = frame.MTOW_kg
    el = f'{np.radians(ELEV_MAX_DEG):.6f}'
    al = f'{np.radians(AIL_MAX_DEG):.6f}'
    rd = f'{np.radians(RUD_MAX_DEG):.6f}'

    def fn(name, expr, axis=None):
        inner = f'  <function name="aero/coefficient/{name}">\n{expr}  </function>\n'
        return inner

    cl_expr = f'''    <product>
      <sum>
        <table><independentVar>alpha</independentVar>{table_xml(a_rad, cl, '          ')}        </table>
        <product><property>fcs/elevator-pos-rad</property><value>{cd['clde']:.5f}</value></product>
      </sum>
      <property>aero/qbar-area</property>
    </product>
'''
    cd_expr = f'''    <product>
      <table><independentVar>alpha</independentVar>{table_xml(a_rad, cdd, '        ')}      </table>
      <property>aero/qbar-area</property>
    </product>
'''
    cm_expr = f'''    <product>
      <sum>
        <table><independentVar>alpha</independentVar>{table_xml(a_rad, cm, '          ')}        </table>
        <product><property>fcs/elevator-pos-rad</property><value>{cd['cmde']:.5f}</value></product>
      </sum>
      <property>aero/qbar-area</property>
      <property>metrics/cbarw-ft</property>
    </product>
'''
    cy_expr = f'''    <product>
      <sum>
        <product><property>aero/beta</property><value>{cyb:.5f}</value></product>
        <product><property>fcs/rudder-pos-rad</property><value>{cd['cydr']:.5f}</value></product>
      </sum>
      <property>aero/qbar-area</property>
    </product>
'''
    cll_expr = f'''    <product>
      <sum>
        <product><property>aero/beta</property><value>{clb:.5f}</value></product>
        <product><property>fcs/aileron-pos-rad</property><value>{cd['clda']:.5f}</value></product>
      </sum>
      <property>aero/qbar-area</property>
      <property>metrics/bw-ft</property>
    </product>
'''
    cn_expr = f'''    <product>
      <sum>
        <product><property>aero/beta</property><value>{cnb:.5f}</value></product>
        <product><property>fcs/rudder-pos-rad</property><value>{cd['cndr']:.5f}</value></product>
        <product><property>fcs/aileron-pos-rad</property><value>{cd['cnda']:.5f}</value></product>
      </sum>
      <property>aero/qbar-area</property>
      <property>metrics/bw-ft</property>
    </product>
'''
    if pack.n_eng == 1:
        prop = f'''    <engine file="eng{key}">
      <name>ENG</name>
      <throttle>fcs/throttle-cmd-norm</throttle>
      <location unit="M"><x>0.55</x><y>0.0</y><z>0.0</z></location>
      <thruster file="prop{key}">
        <location unit="M"><x>0.80</x><y>0.0</y><z>0.0</z></location>
      </thruster>
    </engine>
'''
    else:
        prop = ''
        for sgn in ('L', 'R'):
            y = -0.75 if sgn == 'L' else 0.75
            prop += f'''    <engine file="eng{key}">
      <name>ENG_{sgn}</name>
      <throttle>fcs/throttle-cmd-norm</throttle>
      <location unit="M"><x>0.74</x><y>{y}</y><z>0.10</z></location>
      <thruster file="prop{key}">
        <location unit="M"><x>0.99</x><y>{y}</y><z>0.10</z></location>
      </thruster>
    </engine>
'''
    prop += f'''    <tank type="FUEL">
      <capacity unit="LBS">{max(frame.MTOW_kg * 0.3 * 2.2046, 5):.1f}</capacity>
      <contents unit="LBS">{max(frame.MTOW_kg * 0.3 * 2.2046, 5):.1f}</contents>
      <center unit="M"><x>{0.42 * frame.l_fus}</x><y>0.0</y><z>0.0</z></center>
    </tank>
'''
    xml = f'''<?xml version="1.0"?>
<fdm_config name="uav{key}" version="2.0" release="LAYER9">
  <metrics>
    <wingarea unit="M2">{frame.S:.4f}</wingarea>
    <wingspan unit="M">{geo.b:.4f}</wingspan>
    <chord unit="M">{geo.MAC:.4f}</chord>
    <htailarea unit="M2">{geo.S_h:.4f}</htailarea>
    <htailarm unit="M">1.15</htailarm>
    <vtailarea unit="M2">{geo.S_v:.4f}</vtailarea>
    <vtailarm unit="M">1.15</vtailarm>
    <location name="AERORP" unit="M"><x>{x_cg:.4f}</x><y>0.0</y><z>0.05</z></location>
    <location name="EYEPOINT" unit="M"><x>0.5</x><y>0.0</y><z>0.10</z></location>
    <location name="VRP" unit="M"><x>{x_cg:.4f}</x><y>0.0</y><z>0.05</z></location>
  </metrics>
  <mass_balance>
    <ixx unit="KG*M2">{ixx:.3f}</ixx>
    <iyy unit="KG*M2">{iyy:.3f}</iyy>
    <izz unit="KG*M2">{izz:.3f}</izz>
    <emptywt unit="N">{m * G:.1f}</emptywt>
    <location name="CG" unit="M"><x>{x_cg:.4f}</x><y>0.0</y><z>0.05</z></location>
  </mass_balance>
  <propulsion>
{prop}  </propulsion>
  <flight_control>
    <channel name="pitch">
      <gain name="elevator gain" input="fcs/elevator-cmd-norm" gain="{el}" output="fcs/elevator-pos-rad"/>
    </channel>
    <channel name="roll">
      <gain name="aileron gain" input="fcs/aileron-cmd-norm" gain="{al}" output="fcs/aileron-pos-rad"/>
    </channel>
    <channel name="yaw">
      <gain name="rudder gain" input="fcs/rudder-cmd-norm" gain="{rd}" output="fcs/rudder-pos-rad"/>
    </channel>
  </flight_control>
  <aerodynamics>
    <axis name="LIFT">{fn('CL', cl_expr)}</axis>
    <axis name="DRAG">{fn('CD', cd_expr)}</axis>
    <axis name="PITCH">{fn('CM', cm_expr)}</axis>
    <axis name="SIDE">{fn('CY', cy_expr)}</axis>
    <axis name="ROLL">{fn('CLl', cll_expr)}</axis>
    <axis name="YAW">{fn('CN', cn_expr)}</axis>
  </aerodynamics>
</fdm_config>
'''
    (JROOT / 'aircraft' / f'uav{key}').mkdir(parents=True, exist_ok=True)
    (JROOT / 'aircraft' / f'uav{key}' / f'uav{key}.xml').write_text(xml)


def build_all():
    JROOT.mkdir(exist_ok=True)
    for key, cfg in PROTOTYPES.items():
        build_engine(key, cfg['pack'])
        build_prop(key, cfg['pack'])
        build_aircraft(key, cfg)
    print(f'models written to {JROOT}')


# ---------------------------------------------------------------- analysis
def make_fdm(key):
    fdm = jsbsim.FGFDMExec(str(JROOT))
    fdm.load_model(f'uav{key}')
    return fdm


def trim_case(key, v, gamma_deg=0.0, h=500.0):
    fdm = make_fdm(key)
    fdm['ic/h-sl-meters'] = h
    fdm['ic/vtrue-fps'] = v / 0.3048
    fdm['ic/gamma-deg'] = gamma_deg
    fdm['ic/alpha-deg'] = 4.0
    fdm.run_ic()
    status = fdm.do_trim(0)          # tLongitudinal: steady level at IC speed
    return fdm, status


def state_vec(fdm):
    return np.array([fdm['velocities/v-wind-u'], fdm['velocities/v-wind-w'],
                     fdm['velocities/q-rad_sec'], fdm['attitude/theta-rad']])


def lat_state(fdm):
    return np.array([fdm['velocities/v-wind-v'], fdm['velocities/p-rad_sec'],
                     fdm['velocities/r-rad_sec'], fdm['attitude/phi-rad']])


def set_state(fdm, x, which):
    names = (['velocities/v-wind-u', 'velocities/v-wind-w', 'velocities/q-rad_sec',
              'attitude/theta-rad'] if which == 'long'
             else ['velocities/v-wind-v', 'velocities/p-rad_sec',
                   'velocities/r-rad_sec', 'attitude/phi-rad'])
    for n, val in zip(names, x):
        fdm[n] = val


def jacobian(fdm, which, eps=1e-5, dt=2e-4):
    x0 = state_vec(fdm) if which == 'long' else lat_state(fdm)
    n = len(x0)
    a = np.zeros((n, n))
    fdm['simulation/dt'] = dt
    for i in range(n):
        xp = x0.copy()
        xp[i] += eps
        set_state(fdm, xp, which)
        fdm.run()
        x_plus = state_vec(fdm) if which == 'long' else lat_state(fdm)
        xp = x0.copy()
        xp[i] -= eps
        set_state(fdm, xp, which)
        fdm.run()
        x_minus = state_vec(fdm) if which == 'long' else lat_state(fdm)
        a[:, i] = (x_plus - x_minus) / (2.0 * eps * dt)
    set_state(fdm, x0, which)
    return a


def modes_of(a):
    """Return list of (label, freq_hz, damp, tau) for each eigen-pair/real."""
    out = []
    ev = np.linalg.eigvals(a)
    used = set()
    for i, lam in enumerate(ev):
        if i in used:
            continue
        if abs(lam.imag) > 1e-8:
            j = int(np.argmin(abs(ev - lam.conj())))
            used |= {i, j}
            w = abs(lam.imag) / (2 * np.pi)
            z = -lam.real / abs(lam)
            out.append(('complex', w, z, 1.0 / max(abs(lam.real), 1e-9)))
        else:
            out.append(('real', 0.0, 1.0, -1.0 / lam.real if lam.real < 0 else np.inf))
    return out


def main():
    build_all()
    rows_t, rows_m = [], []
    for key in PROTOTYPES:
        frame, pack = PROTOTYPES[key]['frame'], PROTOTYPES[key]['pack']
        ac = Aircraft(frame, pack, MISSION)
        met = ac.metrics(RHO0)
        # trims
        for tag, v in (('cruise', met['V_cruise']), ('approach', 1.3 * met['V_s_clean'])):
            fdm, st = trim_case(key, v)
            rows_t.append({'prototype': key, 'case': tag, 'V': v,
                           'trim_ok': st,
                           'alpha': np.degrees(fdm['aero/alpha']),
                           'elevator': np.degrees(fdm['fcs/elevator-pos-rad']),
                           'throttle': fdm['propulsion/total-thrust'] / (pack.P_kw * 1000.0),
                           'CL': fdm['aero/cl-total'], 'Cm': fdm['aero/cm-total']})
        fdm, st = trim_case(key, met['V_best_ROC'], gamma_deg=3.0)
        rows_t.append({'prototype': key, 'case': 'climb3deg', 'V': met['V_best_ROC'],
                       'trim_ok': st, 'alpha': np.degrees(fdm['aero/alpha']),
                       'elevator': np.degrees(fdm['fcs/elevator-pos-rad']),
                       'throttle': fdm['propulsion/total-thrust'] / (pack.P_kw * 1000.0),
                       'CL': fdm['aero/cl-total'], 'Cm': fdm['aero/cm-total']})
        # modes around cruise trim
        fdm, _ = trim_case(key, met['V_cruise'])
        a_long = jacobian(fdm, 'long')
        a_lat = jacobian(fdm, 'lat')
        ml = modes_of(a_long)
        mlat = modes_of(a_lat)
        for lab, w, z, tau in ml:
            rows_m.append({'prototype': key, 'axis': 'long', 'mode': lab,
                           'freq_hz': w, 'damping': z, 'tau_s': tau})
        for lab, w, z, tau in mlat:
            rows_m.append({'prototype': key, 'axis': 'lat', 'mode': lab,
                           'freq_hz': w, 'damping': z, 'tau_s': tau})
    pd.DataFrame(rows_t).to_csv(OUT / 'jsbsim_trim.csv', index=False)
    pd.DataFrame(rows_m).to_csv(OUT / 'jsbsim_modes.csv', index=False)
    print('\n=== trim ===')
    print(pd.DataFrame(rows_t).to_string(index=False, float_format=lambda x: f'{x:8.3f}'))
    print('\n=== modes ===')
    print(pd.DataFrame(rows_m).to_string(index=False, float_format=lambda x: f'{x:8.4f}'))


if __name__ == '__main__':
    main()
