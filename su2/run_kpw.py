"""Layer 12: D prop-wing-flap closure (k_pw).

Runs the D half-wing model (sd7062, 25% flap, tractor actuator disk) with SU2
and computes the powered-lift closure factor

    k_pw = dCL_CFD(powered) / [ dCL_flap * (q_ratio - 1) ]

where dCL_flap is the unpowered flap increment at the same incidence and
q_ratio is the layer-11 momentum-disk slipstream ratio at the case speed.
k_pw ~ 1 means the simplified q-boost model is right; 0.7-0.9 means mildly
optimistic; 0.4-0.6 means D is over-sized.

Case TO2: 1.2 V_s, 20 deg flap, full throttle (layer-11 D takeoff state).
Reference area is the modelled half wing (constant MAC x half span).
"""
import csv
import sys
from pathlib import Path

import numpy as np

ROOT = Path('/home/yanlu/Documents/plane')
BASE = ROOT / 'su2'
OUT = ROOT / 'results'
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(BASE))

CASE = dict(name='TO2', V=14.07, flap=20.0, n_rpm=6900.0, aoa=6.0, mach=0.041)
DISK = dict(x=-0.62, y=0.75, z=0.055, radius=0.265)
# half-wing reference area = MAC * half span (constant-chord idealisation)
S_REF = 0.6788 * 0.5 * 6.317


def latest(hist):
    rows = list(csv.reader(open(hist)))
    hdr = [c.strip().replace('"', '') for c in rows[0]]
    iCL, iCD, iN = hdr.index('CL'), hdr.index('CD'), hdr.index('rms[Rho]')
    data = [x for x in rows[1:] if len(x) > max(iCL, iCD, iN)
            and x[0].strip() == '0']
    x = data[-1]
    return float(x[iCL]), float(x[iCD]), float(x[iN]), len(data)


def main():
    from sizing.configs import PROTOTYPES, MISSION      # noqa: E402
    from sizing.mission import Aircraft                 # noqa: E402
    from sizing.ground import GroundRun                 # noqa: E402
    from wing3d import build                            # noqa: E402
    from su2case import write_cfg_3d, run_case          # noqa: E402
    from actdisk import write_actdisk_input, disk_from_bemt  # noqa: E402

    cfg_proto = PROTOTYPES['D']
    ac = Aircraft(cfg_proto['frame'], cfg_proto['pack'], MISSION)
    gr = GroundRun('D', cfg_proto['frame'], cfg_proto['pack'])
    prop = cfg_proto['pack'].prop
    V = CASE['V']
    n_rps = CASE['n_rpm'] / 60.0
    print(f'reference area (half wing) = {S_REF:.4f} m2')

    # ---- meshes ------------------------------------------------------
    # SU2 requires every mesh marker to appear in the config, so the prop-off
    # cases need a mesh without the disk baffle
    meshes = {}
    combos = (('flap_on', CASE['flap'], True),
              ('flap_off', CASE['flap'], False),
              ('clean_off', 0.0, False))
    for key, delta, disk_on in combos:
        path = BASE / f'cases/D12_{key}.su2'
        info = build(str(path), str(ROOT / 'data/airfoils/sd7062.dat'),
                     delta_deg=delta, far=8.0, n_span=20, wall_size=0.012,
                     far_size=1.2, disk=DISK if disk_on else None)
        meshes[key] = path
        print(f'{key}: {info["npoin"]} nodes, {info["nelem"]} elems, '
              f'{info["nmark"]} markers {info["groups"]}')

    # ---- actuator-disk input -----------------------------------------
    disk = disk_from_bemt(prop, V, n_rps,
                          markers=('disk',),
                          center=(DISK['x'], DISK['y'], DISK['z']),
                          axis=(1.0, 0.0, 0.0))
    adfile = BASE / 'cases/D12_actdisk.dat'
    write_actdisk_input(str(adfile), [disk])
    # prop-off: the mesh still carries the disk markers, so declare the same
    # boundary with zero load (a frictionless internal surface) rather than
    # leaving an undeclared marker in the config
    zeros = dict(disk)
    zeros['d_ct'] = [0.0] * len(disk['r_bar'])
    zeros['d_cp'] = [0.0] * len(disk['r_bar'])
    zero_file = BASE / 'cases/D12_actdisk_zero.dat'
    write_actdisk_input(str(zero_file), [zeros])
    print(f'disk: J={disk["advance_ratio"]:.3f} '
          f'sum dCT={sum(disk["d_ct"]):.4f} R={disk["radius"]:.3f}')
    q_ratio = gr.q_slip_ratio(V, gr.thrust(V, 1.0) * cfg_proto['pack'].n_eng)

    # ---- runs ---------------------------------------------------------
    rows = []
    for key, path in meshes.items():
        disk_on = key.endswith('_on')
        name = f'V12_{CASE["name"]}_{key}'
        ad = ('disk',) if disk_on else None
        if True:
            cfg = write_cfg_3d(
                str(BASE / f'cases/{name}.cfg'), str(path), ('wing',),
                sym='symmetry', far='farfield', plot=('wing',), actdisk=ad,
                mach=CASE['mach'], aoa=CASE['aoa'],
                re=V * ac.geo.MAC / 1.46e-5, ref_length=ac.geo.MAC,
                ref_area=S_REF, x_ref=0.25, iters=3000, cfl=5.0,
                actdisk_file=str(adfile),
                extra_walls=('tip',))
            r = run_case(cfg, str(BASE / f'cases/{name}run'), np_threads=24,
                         timeout=7200)
            hist = BASE / f'cases/{name}run' / 'history.csv'
            if not hist.exists():
                print(f'{name}: no history (rc={r.returncode})')
                continue
            cl, cd, rms, nit = latest(hist)
            rows.append(dict(name=name, tag=key, disk=int(disk_on),
                             CL=round(cl, 4), CD=round(cd, 5),
                             rms=round(rms, 3), iters=nit))
            print(f'{name}: CL={cl:.4f} CD={cd:.5f} rms={rms:.2f} it={nit}')

    # ---- k_pw ---------------------------------------------------------
    get = {(r_['tag'], r_['disk']): r_ for r_ in rows}
    if (('flap_on', 1) in get and ('flap_off', 0) in get
            and ('clean_off', 0) in get):
        d_on = get[('flap_on', 1)]['CL'] - get[('flap_off', 0)]['CL']
        d_flap = get[('flap_off', 0)]['CL'] - get[('clean_off', 0)]['CL']
        d_model = d_flap * (q_ratio - 1.0)
        k_pw = d_on / d_model if d_model else float('nan')
        print(f'\nq_ratio={q_ratio:.2f}')
        print(f'dCL_powered={d_on:+.4f}  dCL_flap={d_flap:+.4f}  '
              f'dCL_model={d_model:+.4f}')
        print(f'k_pw = {k_pw:.3f}')
        rows.append(dict(name='K_PW', tag='-', disk=-1, CL=round(k_pw, 4),
                         CD=round(d_model, 4), rms=round(q_ratio, 4),
                         iters=0))

    OUT.mkdir(exist_ok=True)
    with open(OUT / 'su2_kpw.csv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f'saved: {OUT}/su2_kpw.csv')


if __name__ == '__main__':
    main()
