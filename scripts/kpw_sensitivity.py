"""Layer 12 / Layer 11: k_pw sensitivity band (no CFD required).

Layer 12 CFD has cleared the reference-area/frame chain (0.000%), the
actuator-disk input integration (1e-16) and - via Euler - the
geometry/topology chain, but no converged wall-resolved 3-D viscous case
exists yet. So k_pw is reported as a SENSITIVITY BAND over 0.4-1.0 rather
than as a nominal CFD calibration.

The engineering question this answers immediately: if the true k_pw turns out
to be only 0.4, does D still beat A/B/C on takeoff and landing? If yes,
Layer 12 only calibrates the size of D's advantage; if D only works for
k_pw > 0.7, then CFD is a genuine design gate.

Emits results/kpw_sensitivity.csv and prints the stored-distance crossover.
"""
import csv
import sys
from pathlib import Path

R = Path('/home/yanlu/Documents/plane')
sys.path.insert(0, str(R))

from sizing.configs import PROTOTYPES, MISSION     # noqa: E402
from sizing.ground import GroundRun, K_PW_BAND     # noqa: E402


def run(key, k_pw):
    cfg = PROTOTYPES[key]
    gr = GroundRun(key, cfg['frame'], cfg['pack'])
    gr.k_pw = k_pw
    to = gr.takeoff()
    ld = gr.landing()
    return to, ld


def main():
    rows = []
    for key in PROTOTYPES:
        for k in K_PW_BAND:
            to, ld = run(key, k)
            rows.append(dict(prototype=key, k_pw=k,
                             s_TO_m=round(to.s_total, 1),
                             V_s_TO=round(to.V_s, 2),
                             s_land_m=round(ld.s_total, 1),
                             V_s_land=round(ld.V_s, 2),
                             V_app=round(ld.V_approach, 2)))
            print(f'{key} k_pw={k:.1f}  s_TO={to.s_total:6.1f} m  '
                  f's_land={ld.s_total:6.1f} m  V_s_TO={to.V_s:5.2f} m/s',
                  flush=True)

    out = R / 'results'
    out.mkdir(exist_ok=True)
    with open(out / 'kpw_sensitivity.csv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    # ---- the design question: is D's STOL advantage robust to k_pw? ----
    print('\n=== D vs best of A/B/C (STOL advantage) ===')
    print('k_pw   s_TO_D   best_other   D gain   s_land_D  best_other   D gain')
    for k in K_PW_BAND:
        d = [r for r in rows if r['prototype'] == 'D' and r['k_pw'] == k][0]
        others_to = min(r['s_TO_m'] for r in rows
                        if r['prototype'] != 'D' and r['k_pw'] == k)
        others_ld = min(r['s_land_m'] for r in rows
                        if r['prototype'] != 'D' and r['k_pw'] == k)
        print(f'{k:.1f}    {d["s_TO_m"]:6.1f}   {others_to:8.1f}   '
              f'{others_to - d["s_TO_m"]:6.1f}   {d["s_land_m"]:7.1f}   '
              f'{others_ld:8.1f}   {others_ld - d["s_land_m"]:6.1f}')

    print(f'\nsaved: {out}/kpw_sensitivity.csv')


if __name__ == '__main__':
    main()
