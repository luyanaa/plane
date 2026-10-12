"""Layer 4: parametric OpenVSP geometry for prototypes A-D.

Mirrors the layer-1 sizing configs (sizing/configs.py) and the tail-volume
geometry (sizing/geometry.py): fuselage, wing, HTP, VTP, engine nacelles,
propeller disks. Saves vsp_models/<key>.vsp3 and verifies wing span/area and
tail areas against the sizing values.

Run headless on tinybox:  .venv/bin/python scripts/vsp_geometry.py
"""
import os
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import openvsp as vsp

from sizing.configs import PROTOTYPES
from sizing.geometry import build_geometry

OUT = Path(__file__).resolve().parents[1] / 'vsp_models'

# tail planforms (AR choices typical for this class)
HTP_AR = 4.0
VTP_AR = 1.8
TAIL_TAPER = 0.6

# prototype -> XFOIL-verified airfoil file in data/airfoils/
AIRFOIL = {'A': 'naca2412', 'B': 'sd7062', 'C': 'naca2412', 'D': 'sd7062'}
AIRFOIL_DIR = Path(__file__).resolve().parents[1] / 'data' / 'airfoils'


def set_parm(gid, name, group, value):
    """Set ALL parms of (name, group) on the geom (flat pid list)."""
    hits = [pid for pid in vsp.GetGeomParmIDs(gid)
            if vsp.GetParmName(pid) == name and vsp.GetParmGroupName(pid) == group]
    if not hits:
        raise RuntimeError(f'parm not found: {name}/{group} on {vsp.GetGeomName(gid)}')
    for pid in hits:
        vsp.SetParmVal(pid, value)


def set_xsec(xsec, name, value):
    vsp.SetParmVal(vsp.FindParm(xsec, name, 'XSec'), value)


def add_wing(parent, area, ar, taper, x, z, sweep=0.0, dihedral=0.0, airfoil=None):
    """Wing from area/AR/taper.

    VSP 3.53: the LAST XSec drives the planform (Span = half-span per side,
    doubled by symmetry); the first XSec is collapsed to a stub.
    """
    b = float(np.sqrt(area * ar))
    c_root = 2.0 * area / (b * (1.0 + taper))
    c_tip = c_root * taper
    gid = vsp.AddGeom('WING') if not parent else vsp.AddGeom('WING', parent)
    surf = vsp.GetXSecSurf(gid, 0)
    x0, x1 = vsp.GetXSec(surf, 0), vsp.GetXSec(surf, 1)
    # main planform on the last section
    set_xsec(x1, 'Span', 0.5 * b)
    set_xsec(x1, 'Root_Chord', c_root)
    set_xsec(x1, 'Tip_Chord', c_tip)
    set_xsec(x1, 'Sweep', sweep)
    set_xsec(x1, 'Dihedral', dihedral)
    # collapsed stub at the root (0.15 m, non-degenerate for VSPAERO meshing)
    set_xsec(x0, 'Span', 0.15)
    set_xsec(x0, 'Root_Chord', c_root)
    set_xsec(x0, 'Tip_Chord', c_root)
    set_parm(gid, 'X_Rel_Location', 'XForm', x)
    set_parm(gid, 'Z_Rel_Location', 'XForm', z)
    if airfoil:
        af = AIRFOIL_DIR / f'{airfoil}.dat'
        if af.exists():
            surf = vsp.GetXSecSurf(gid, 0)
            vsp.ChangeXSecShape(surf, 0, vsp.XS_FILE_AIRFOIL)
            vsp.ChangeXSecShape(surf, 1, vsp.XS_FILE_AIRFOIL)
            vsp.ReadFileAirfoil(vsp.GetXSec(surf, 0), str(af))
            vsp.ReadFileAirfoil(vsp.GetXSec(surf, 1), str(af))
        else:
            print(f'  [warn] airfoil file {af} missing, using default')
    return gid, b, c_root


def build_prototype(key, cfg):
    frame, pack = cfg['frame'], cfg['pack']
    props_as_pods = os.environ.get('VSP_NO_PROPS') == '1'
    no_wing_pods = os.environ.get('VSP_NO_WING_PODS') == '1'
    g = build_geometry(frame.S, frame.AR, frame.taper, frame.t_over_c,
                       frame.l_fus, frame.d_fus, n_nacelles=frame.n_nacelles)
    print(f'--- {key} {cfg["name"]}  [{cfg["layout"]}]')

    fus = vsp.AddGeom('FUSELAGE')
    set_parm(fus, 'Length', 'Design', frame.l_fus)
    try:
        set_parm(fus, 'Ellipse_Width', 'XSecCurve', frame.d_fus)
        set_parm(fus, 'Ellipse_Height', 'XSecCurve', frame.d_fus)
    except RuntimeError:
        print('  [warn] fuselage diameter not settable, using default body')

    # wing: high wing, x ~ 33% of fuselage
    x_wing = 0.33 * frame.l_fus
    wing, b, c_root = add_wing(fus, frame.S, frame.AR, frame.taper,
                               x_wing, 0.10 * frame.d_fus,
                               airfoil=AIRFOIL[key])

    # HTP / VTP from tail-volume areas, at tail arm l_h=1.25 / l_v=1.15
    x_htp = frame.l_fus - 0.35
    z_h = 0.30 * frame.d_fus          # clear of the body sheet (z=0)
    z_v = 0.36 * frame.d_fus
    htp, b_h, _ = add_wing(fus, g.S_h, HTP_AR, TAIL_TAPER, x_htp, z_h)
    vtp, b_v, _ = add_wing(fus, g.S_v, VTP_AR, TAIL_TAPER, x_htp, z_v)
    set_parm(vtp, 'Y_Rel_Location', 'XForm', 0.0)
    set_parm(vtp, 'Z_Rotation', 'XForm', 90.0)   # vertical fin

    # engines + propeller disks (spinner pods when VSP_NO_PROPS=1, for solver smoke)
    n_eng = pack.n_eng
    if n_eng == 1:
        # central tractor: cowl pod at nose
        pod = vsp.AddGeom('POD', fus)
        set_parm(pod, 'Length', 'Design', 0.45)
        set_parm(pod, 'X_Rel_Location', 'XForm', 0.30)
        if props_as_pods:
            spinner = vsp.AddGeom('POD', pod)
            set_parm(spinner, 'Length', 'Design', 0.15)
            set_parm(spinner, 'X_Rel_Location', 'XForm', 0.4)
        else:
            prop = vsp.AddGeom('PROP', pod)
            set_parm(prop, 'Diameter', 'Design', pack.D_prop)
            set_parm(prop, 'X_Rel_Location', 'XForm', 0.25)
    else:
        # wing-root tractors at y = +-0.75 m (inboard, outside flap span)
        if no_wing_pods:
            print('  [sweep mode] wing pods omitted (thin sheets break VLM)')
        else:
            for sign in (-1, 1):
                pod = vsp.AddGeom('POD', wing)
                set_parm(pod, 'Length', 'Design', 0.60)
                set_parm(pod, 'X_Rel_Location', 'XForm', -0.05)
                set_parm(pod, 'Y_Rel_Location', 'XForm', sign * 0.75)
                set_parm(pod, 'Z_Rel_Location', 'XForm', 0.10)  # above wing plane (tractor)
                if props_as_pods:
                    spinner = vsp.AddGeom('POD', pod)
                    set_parm(spinner, 'Length', 'Design', 0.15)
                    set_parm(spinner, 'X_Rel_Location', 'XForm', 0.5)
                else:
                    prop = vsp.AddGeom('PROP', pod)
                    set_parm(prop, 'Diameter', 'Design', pack.D_prop)
                    set_parm(prop, 'X_Rel_Location', 'XForm', 0.55)

    vsp.Update()
    file = OUT / f'{key}.vsp3'
    vsp.WriteVSPFile(str(file))

    # verify against sizing
    tot_span = vsp.GetParmVal(vsp.FindParm(wing, 'TotalSpan', 'WingGeom'))
    tot_area = vsp.GetParmVal(vsp.FindParm(wing, 'TotalArea', 'WingGeom'))

    def absx(gid):
        pids = [p for p in vsp.GetGeomParmIDs(gid) if vsp.GetParmName(p) == 'X_Location']
        return vsp.GetParmVal(pids[0]) if pids else float('nan')

    return {'file': file, 'b_target': g.b, 'b_vsp': tot_span,
            'S_target': frame.S, 'S_vsp': tot_area,
            'S_h_target': g.S_h, 'S_v_target': g.S_v,
            'x_wing': absx(wing), 'x_htp': absx(htp), 'x_vtp': absx(vtp)}


def main():
    OUT.mkdir(exist_ok=True)
    rows = {}
    for key, cfg in PROTOTYPES.items():
        rows[key] = build_prototype(key, cfg)
    print('\n=== verification (VSP readback vs layer-1 sizing) ===')
    print(f"{'proto':6s} {'b target':>9s} {'b VSP':>8s} {'S target':>9s} {'S VSP':>8s} "
          f"{'S_h':>6s} {'S_v':>6s} {'x_wing':>7s} {'x_htp':>7s} {'x_vtp':>7s}")
    ok = True
    for key, r in rows.items():
        db = abs(r['b_vsp'] - r['b_target']) / r['b_target']
        ds = abs(r['S_vsp'] - r['S_target']) / r['S_target']
        flag = '' if db < 0.01 and ds < 0.01 else '  <-- MISMATCH'
        if flag:
            ok = False
        print(f"{key:6s} {r['b_target']:9.3f} {r['b_vsp']:8.3f} "
              f"{r['S_target']:9.3f} {r['S_vsp']:8.3f} "
              f"{r['S_h_target']:6.3f} {r['S_v_target']:6.3f} "
              f"{r['x_wing']:7.3f} {r['x_htp']:7.3f} {r['x_vtp']:7.3f}{flag}")
    print(f'\nsaved: {OUT}/A.vsp3, B.vsp3, C.vsp3, D.vsp3  |  all match: {ok}')


if __name__ == '__main__':
    main()
