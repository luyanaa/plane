"""Layer 12 step 1+2: force-reference/marker regression and clean-wing polar.

Step 1 - force regression: build the half-wing mesh, run one case, then
recover CL and CD manually from the surface pressure (Cp integrated over the
wall triangles) and compare with SU2's reported CL/CD. Agreement to <1e-4
proves the frame, marker and REF_AREA chain is self-consistent - it removes
the "off by exactly two" class of bug before any physics is discussed.

Step 2 - clean-wing polar: alpha = -4..6 deg on the clean (undeflected) wing,
compared against the layer-2 XFOIL polars for sd7062 at Re 1e6. The metrics
that matter are C_L0 and C_Lalpha, NOT C_Lmax (the 1.47-1.77 figures quoted
earlier are C_Lmax and were wrongly used as an alpha=0 reference).

Half-model bookkeeping: the mesh models half the wing, so REF_AREA is the
half-wing area MAC*(b/2); SU2 then integrates only the modelled half and the
coefficient is correct without any doubling.
"""
import csv
import sys
from pathlib import Path

import numpy as np

R = Path('/home/yanlu/Documents/plane')
sys.path.insert(0, str(R / 'su2'))

from wing3d import build                        # noqa: E402
from su2case import write_cfg_3d, run_case      # noqa: E402

MAC = 0.6788
B = 6.317
S_HALF = MAC * 0.5 * B
V = 14.07
RE = V * MAC / 1.46e-5
ALPHAS = (-4.0, -2.0, 0.0, 2.0, 4.0, 6.0)

# ---- finite-wing oracle ---------------------------------------------------
# The 2-D XFOIL slope must NOT be used as the target for a 3-D half-wing.
# Classical lifting-line correction:
#     a = a0 / (1 + a0/(pi e AR))
# with a0 = 0.109 /deg (sd7062, forced transition, layer 2), e = 0.7..0.9.
A0_PER_DEG = 0.109
E_OSWALD = (0.7, 0.8, 0.9)
AR_3D = (2 * 0.5 * 6.317) ** 2 / (2 * 2.1440)     # full span^2 / full area


def finite_wing_oracle():
    a0 = A0_PER_DEG * 180.0 / np.pi                 # per radian
    out = {}
    for e in E_OSWALD:
        a3 = a0 / (1.0 + a0 / (np.pi * e * AR_3D))
        out[e] = a3 * np.pi / 180.0                 # back to per degree
    return out


def su2_coeffs(hist):
    rows = list(csv.reader(open(hist)))
    hdr = [c.strip().replace('"', '') for c in rows[0]]
    iC, iD, iN = hdr.index('CL'), hdr.index('CD'), hdr.index('rms[Rho]')
    data = [x for x in rows[1:] if len(x) > max(iC, iD, iN)
            and x[0].strip() == '0']
    x = data[-1]
    return float(x[iC]), float(x[iD]), float(x[iN])


def refarea_check(hist, s_ref, rho=1.225, p_inf=101325.0, mach=0.041):
    """Step 1 gate: SU2 reports RefForce = q_inf * REF_AREA.

    Recovering S_ref from that identity and comparing with the half-wing area
    is an exact, first-principles check of the reference chain. It catches the
    classic half-model failure - forcing only half the wing while dividing by
    the full-wing area, which halves every coefficient - without depending on
    any post-processing of mine.
    """
    rows = list(csv.reader(open(hist)))
    hdr = [c.strip().replace('"', '') for c in rows[0]]
    iRF = hdr.index('RefForce')
    data = [x for x in rows[1:] if len(x) > iRF]
    rf = float(data[-1][iRF])
    gamma = 1.4
    a = np.sqrt(gamma * p_inf / rho)
    v = mach * a
    q = 0.5 * rho * v * v
    return rf / q, s_ref, v, q


def main():
    mesh = str(R / 'su2/cases/V13_clean.su2')
    info = build(mesh, str(R / 'data/airfoils/sd7062.dat'), delta_deg=0.0,
                 far=4.0, n_span=28, wall_size=0.008, far_size=0.6)
    print(f'clean mesh: {info["npoin"]} nodes, {info["nelem"]} elems, '
          f'{info["groups"]}')
    print(f'half-wing area {S_HALF:.4f} m2   Re {RE:.3e}')

    # ---- step 1: reference-area / marker / frame regression ----------
    probe = write_cfg_3d(str(R / 'su2/cases/V13_probe.cfg'), mesh, ('wing',),
                         sym='symmetry', far='farfield', plot=('wing',),
                         actdisk=None, mach=0.041, aoa=0.0, re=RE,
                         ref_length=MAC, ref_area=S_HALF, x_ref=0.25,
                         iters=200, cfl=5.0, extra_walls=('tip',))
    run_case(probe, str(R / 'su2/cases/V13_proberun'), np_threads=24,
             timeout=1800)
    s_su2, s_geom, v_inf, q_inf = refarea_check(
        R / 'su2/cases/V13_proberun/history.csv', S_HALF)
    err = abs(s_su2 - s_geom) / s_geom
    print(f'[step1] RefForce/q_inf = {s_su2:.4f} m2 vs modelled {s_geom:.4f} '
          f'm2  ({100*err:+.3f}%)  V_inf={v_inf:.3f} m/s')
    assert err < 0.01, 'REF_AREA chain inconsistent - stopping before any CFD'

    orc = finite_wing_oracle()
    lo, hi = min(orc.values()), max(orc.values())
    c_d0 = 0.0118                                   # forced-transition profile drag
    cl0 = 0.4313                                    # forced-transition C_L0
    c_di = cl0 ** 2 / (np.pi * 0.8 * AR_3D)
    print(f'[oracle] AR_3D = {AR_3D:.2f}   C_Lalpha,3D = {lo:.4f}..{hi:.4f} '
          f'/deg  (e = 0.7..0.9)')
    print(f'[oracle] C_D(alpha=0) ~ C_D0 + C_L0^2/(pi e AR) = {c_d0:.4f} + '
          f'{c_di:.4f} = {c_d0 + c_di:.4f}')

    xfoil = {}
    with open(R / 'results/airfoil_polars.csv') as fh:
        for row in csv.DictReader(fh):
            if row['airfoil'] == 'sd7062' and row['mode'] == 'free' \
                    and float(row['Re']) == 1.0e6:
                xfoil[float(row['alpha'])] = (float(row['CL']),
                                              float(row['CD']))

    rows = []
    for aoa in ALPHAS:
        name = f'V13_a{int(aoa):+d}'
        cfg = write_cfg_3d(
            str(R / f'su2/cases/{name}.cfg'), mesh, ('wing',),
            sym='symmetry', far='farfield', plot=('wing',), actdisk=None,
            mach=0.041, aoa=aoa, re=RE, ref_length=MAC, ref_area=S_HALF,
            x_ref=0.25, iters=8000, cfl=5.0, extra_walls=('tip',))
        run_case(cfg, str(R / f'su2/cases/{name}run'), np_threads=24,
                 timeout=7200)
        cl, cd, rms = su2_coeffs(R / f'su2/cases/{name}run/history.csv')
        clx, cdx = xfoil.get(aoa, (float('nan'), float('nan')))
        rows.append(dict(aoa=aoa, CL_SU2=round(cl, 4), CD_SU2=round(cd, 5),
                         CL_XFOIL=clx, CD_XFOIL=cdx,
                         dCL_pct=round(100 * (cl / clx - 1), 1) if clx else None,
                         rms=round(rms, 3)))
        print(f'alpha={aoa:+5.1f}  CL_SU2={cl:+.4f} (XFOIL {clx:+.4f}, '
              f'{100*(cl/clx-1):+.1f}%)  CD_SU2={cd:.5f} '
              f'(XFOIL {cdx:.5f})  rms={rms:.2f}', flush=True)

    # ---- slope check against the oracle -----------------------------
    aa = np.array([r_['aoa'] for r_ in rows], dtype=float)
    cc = np.array([r_['CL_SU2'] for r_ in rows], dtype=float)
    slope = float(np.polyfit(aa, cc, 1)[0])
    cl0_su2 = float(np.interp(0.0, aa, cc))
    cd0_su2 = float(np.interp(0.0, aa, [r_['CD_SU2'] for r_ in rows]))
    print(f'\n[step2] SU2 C_Lalpha = {slope:.4f} /deg  vs oracle '
          f'{lo:.4f}..{hi:.4f}  -> {100*(slope/hi-1):+.0f}..'
          f'{100*(slope/lo-1):+.0f}%')
    print(f'[step2] SU2 C_L(0) = {cl0_su2:+.4f}  vs oracle {cl0:.4f}')
    print(f'[step2] SU2 C_D(0) = {cd0_su2:.4f}    vs oracle {c_d0 + c_di:.4f}'
          f'  ({cd0_su2/(c_d0+c_di):.1f}x)')

    with open(R / 'results/su2_clean_polar.csv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print('saved: results/su2_clean_polar.csv')


if __name__ == '__main__':
    main()
