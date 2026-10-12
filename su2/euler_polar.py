"""Layer 12 step 2b: Euler angle-of-attack sweep on the SAME tetra mesh.

Cheapest discriminator for the C_L(0) < 0 anomaly. Reuses the existing
8497-node clean-wing mesh and simply turns viscosity off, so any change is
attributable to the viscous near-wall treatment rather than geometry, frame,
reference area or volume topology.

If Euler gives C_L(0) > 0 and a slope near the finite-wing oracle
(0.0835-0.0881 /deg), then geometry/frame/reference/topology are all cleared
and the remaining error is wall resolution. If Euler also gives C_L(0) < 0,
the fault is in the CFD surface facets themselves and needs a camber-line
refit from the exported boundary nodes.
"""
import csv
import sys
from pathlib import Path

import numpy as np

R = Path('/home/yanlu/Documents/plane')
sys.path.insert(0, str(R / 'su2'))

from su2case import run_case                     # noqa: E402

MESH = R / 'su2/cases/V13_clean.su2'
ALPHAS = (-4.0, -2.0, 0.0, 2.0, 4.0, 6.0)
MAC = 0.6788
S_HALF = 2.1440
RE = 14.07 * MAC / 1.46e-5
MACH = 0.041

A0_PER_DEG = 0.109
AR_3D = (2 * 0.5 * 6.317) ** 2 / (2 * 2.1440)


def euler_cfg(path, aoa, iters=4000):
    txt = f"""MESH_FILENAME= {MESH}
SOLVER= EULER
MESH_FORMAT= SU2
MATH_PROBLEM= DIRECT
FLUID_MODEL= STANDARD_AIR
MACH_NUMBER= {MACH}
AOA= {aoa}
SIDESLIP_ANGLE= 0.0
FREESTREAM_PRESSURE= 101325.0
FREESTREAM_TEMPERATURE= 288.15
REF_ORIGIN_MOMENT_X= 0.25
REF_ORIGIN_MOMENT_Y= 0.0
REF_ORIGIN_MOMENT_Z= 0.0
REF_LENGTH= {MAC}
REF_AREA= {S_HALF}
MARKER_SYM= ( symmetry )
MARKER_EULER= ( wing, tip )
MARKER_FAR= ( farfield )
MARKER_PLOTTING= ( wing )
MARKER_MONITORING= ( wing )
NUM_METHOD_GRAD= WEIGHTED_LEAST_SQUARES
CFL_NUMBER= 5.0
CFL_ADAPT= YES
CFL_ADAPT_PARAM= ( 0.1, 2.0, 1.0, 100.0 )
ITER= {iters}
CONV_NUM_METHOD_FLOW= ROE
ENTROPY_FIX_COEFF= 0.1
MUSCL_FLOW= YES
SLOPE_LIMITER_FLOW= VENKATAKRISHNAN
TIME_DISCRE_FLOW= EULER_IMPLICIT
LINEAR_SOLVER= FGMRES
LINEAR_SOLVER_PREC= ILU
LINEAR_SOLVER_ERROR= 1E-4
LINEAR_SOLVER_ITER= 10
HISTORY_OUTPUT= (ITER, RMS_RES, AERO_COEFF)
CONV_FILENAME= history
"""
    Path(path).write_text(txt)
    return path


def latest(hist):
    rows = list(csv.reader(open(hist)))
    hdr = [c.strip().replace('"', '') for c in rows[0]]
    iC, iD, iN = hdr.index('CL'), hdr.index('CD'), hdr.index('rms[Rho]')
    d = [x for x in rows[1:] if len(x) > max(iC, iD, iN)
         and x[0].strip() == '0']
    return float(d[-1][iC]), float(d[-1][iD]), float(d[-1][iN])


def main():
    a0 = A0_PER_DEG * 180.0 / np.pi
    oracle = {e: a0 / (1 + a0 / (np.pi * e * AR_3D)) * np.pi / 180.0
              for e in (0.7, 0.8, 0.9)}
    lo, hi = min(oracle.values()), max(oracle.values())
    print(f'oracle C_Lalpha = {lo:.4f}..{hi:.4f} /deg  (AR_3D={AR_3D:.2f})')

    rows = []
    for aoa in ALPHAS:
        name = f'E2_a{int(aoa):+d}'
        cfg = euler_cfg(R / f'su2/cases/{name}.cfg', aoa)
        run_case(str(cfg), str(R / f'su2/cases/{name}run'), np_threads=24,
                 timeout=3600)
        cl, cd, rms = latest(R / f'su2/cases/{name}run' / 'history.csv')
        rows.append(dict(aoa=aoa, CL=round(cl, 4), CD=round(cd, 5),
                         rms=round(rms, 3)))
        print(f'  Euler alpha={aoa:+5.1f}  CL={cl:+.4f}  CD={cd:.5f}  '
              f'rms={rms:.2f}', flush=True)

    aa = np.array([r['aoa'] for r in rows])
    cc = np.array([r['CL'] for r in rows])
    slope = float(np.polyfit(aa, cc, 1)[0])
    cl0 = float(np.interp(0.0, aa, cc))
    print(f'\n[step2b] Euler C_Lalpha = {slope:.4f} /deg  vs oracle '
          f'{lo:.4f}..{hi:.4f}  ({100*(slope/hi-1):+.0f}..'
          f'{100*(slope/lo-1):+.0f}%)')
    print(f'[step2b] Euler C_L(0) = {cl0:+.4f}   (viscous run gave -0.0869; '
          f'oracle +0.43)')
    verdict = ('geometry/frame/reference/topology CLEARED - remaining error '
               'is wall resolution' if cl0 > 0.2 and slope > 0.075 else
               'geometry or surface facets still suspect')
    print(f'[step2b] VERDICT: {verdict}')

    with open(R / 'results/su2_euler_polar.csv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print('saved: results/su2_euler_polar.csv')


if __name__ == '__main__':
    main()
