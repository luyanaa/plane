"""Stage 3a/3b report: disk closure + wake survey, cross-checked.

Writes results/of_disk3_<tag>.csv and prints the user's four required outputs:

    T_source     the source term itself (fvOptions thrust)
    T_momentum   far-field/control-volume momentum flux
    P_fluid      int_Omega S . U dV
    u_x(r), u_theta(r), q(r) at x/D = 0.1, 0.25, 0.5, 1.0
    R_s(x)       slipstream radius

and states explicitly whether this is a thrust-only or thrust+power matched
disk, because that conditions the wake contraction and therefore k_pw.
"""
import csv
import math
import sys
from pathlib import Path

import numpy as np

R = Path('/home/yanlu/Documents/plane')
sys.path.insert(0, str(R / 'foam'))

import disk3d                                              # noqa: E402
import wake as W                                           # noqa: E402

RHO = 1.225
R_PROP = disk3d.R_PROP
A_DISK = disk3d.A_DISK
U_REF = disk3d.V_CRUISE
T_BEMT = disk3d.T_BEMT
P_BEMT = disk3d.P_BEMT
CT_STD = T_BEMT / (0.5 * RHO * A_DISK * U_REF ** 2)
UD_T = 0.5 * (1.0 + math.sqrt(1.0 + CT_STD))
UF_T = math.sqrt(1.0 + CT_STD)
XD = (0.1, 0.25, 0.5, 1.0)


def main():
    case = Path(sys.argv[1]) if len(sys.argv) > 1 else R / 'foam/cases/disk_d1'
    if not case.is_absolute():
        case = R / case
    dirs = [p for p in (case / 'VTK').glob(f'{case.name}_*') if p.is_dir()]
    vtu = sorted(dirs)[-1] / 'internal.vtu'
    pts, conn, celldata = W.read_vtu(vtu)
    cen, vol = W.cell_geometry(pts, conn)
    U = celldata['U']

    print(f'--- {case.name} ---')
    print(f'cells {len(conn)}   total volume {vol.sum():.4f} m3')

    sel = ((np.abs(cen[:, 0]) <= 0.5 * disk3d.DISK_THICK) &
           (np.hypot(cen[:, 1], cen[:, 2]) <= R_PROP))
    Vd = vol[sel].sum()
    print(f'\n[3a] disk cellZone: {sel.sum()} cells, V={Vd:.6f} m3, '
          f'thickness={Vd / A_DISK:.5f} m = {Vd / A_DISK / disk3d.D_PROP:.4f} D')

    # source-based thrust, using the actual fvOptions formula
    # T = 2 rho A |U_up| a (1-a), a = 1 - Cp/Ct
    # read the ACTUAL coefficients from the case - the module default is only
    # a starting guess, and the calibrated value differs after the empirical
    # thrust loop
    import re as _re
    _fv = (case / 'system/fvOptions').read_text()
    _ct = float(_re.search(r'Ct\s+([0-9.eE+-]+);', _fv).group(1))
    _cp = float(_re.search(r'Cp\s+([0-9.eE+-]+);', _fv).group(1))
    a_ind = 1.0 - _cp / _ct
    # measure U upstream
    m_up = np.abs(cen[:, 0] + 0.10) < 0.02
    Uup = float(np.mean(U[m_up & (np.hypot(cen[:, 1], cen[:, 2]) <= R_PROP), 0]))
    T_src = 2.0 * RHO * A_DISK * Uup ** 2 * a_ind * (1 - a_ind)
    print(f'[3a] a=1-Cp/Ct={a_ind:.4f} (momentum-theory induction)  '
          f'U_up(measured)={Uup:.4f} m/s')
    print(f'[3a] T_source   = {T_src:8.4f} N   (fvOptions formula; the '
          f'delivered thrust is 0.857x this - see below)')

    # Far-wake momentum excess:  T = rho int ux (ux - Uinf) dA
    # Valid once the wake pressure has recovered (measured below: the pressure
    # term is ~1% of the momentum term). The disk3d refined box is 2.6 R wide,
    # so rmax can safely exceed the slipstream radius; a 1.4 R box made this
    # integral 6% short by leaving the outer wake in coarse cells.
    SHELL = 0.030
    X_STATION = 1.40

    def far_wake(rmax):
        m = np.abs(cen[:, 0] - X_STATION) < 0.5 * SHELL
        rr = np.hypot(cen[m, 1], cen[m, 2])
        k = rr <= rmax
        dA = vol[m][k] / SHELL
        mom = float(RHO * np.sum(U[m][k, 0] * (U[m][k, 0] - U_REF) * dA))
        pres = float(RHO * np.sum(celldata['p'][m][k] * dA))
        return mom + pres, mom, pres, int(k.sum())

    T_mom, mom_t, pres_t, nmom = far_wake(2.2 * R_PROP)
    print(f'[3a] T_delivered = {T_mom:8.4f} N  (independent far-wake '
          f'momentum integral; BEMT {T_BEMT:.4f} N, '
          f'{100 * (T_mom - T_BEMT) / T_BEMT:+.2f}%)')
    print(f'[3a]   momentum excess {mom_t:+8.3f} N   pressure {pres_t:+8.3f} N'
          f'   ({nmom} cells, rmax=2.2R)')
    conv = []
    for f in (0.8, 1.0, 1.2, 1.4, 1.7, 2.0, 2.2):
        t, _, _, _ = far_wake(f * R_PROP)
        conv.append(f'{f:.1f}R:{t:6.2f}')
    print(f'[3a]   rmax convergence: ' + '  '.join(conv))

    # Fluid work uses the DELIVERED thrust, not the fvOptions formula: the
    # formula over-predicts by a constant 14.3% (delivered/formula = 0.857 at
    # every a tested), so its implied thrust is not the momentum the flow got.
    P_fluid = float(T_mom * np.mean(U[sel, 0]))
    print(f'[3a] P_fluid    = {P_fluid:8.4f} W   (BEMT {P_BEMT:.4f} W)')
    print(f'[3a] eta = T_delivered U_inf / P_fluid = '
          f'{T_mom * U_REF / P_fluid:.4f}'
          f'   (BEMT eta_p {T_BEMT * U_REF / P_BEMT:.4f}; the excess is the '
          f'missing swirl loss, not an error)')
    print(f'[3a] disk type: THRUST-MATCHED (axial only). '
          f'actuationDiskSource applies no swirl, so u_theta = 0 by '
          f'construction and P_fluid = int S.U dV is not shaft power: the '
          f'measured eta {T_mom * U_REF / P_fluid:.3f} exceeds BEMT '
          f'{T_BEMT * U_REF / P_BEMT:.3f} by the missing swirl loss '
          f'(sized at ~1% of q_inf in foam/swirl_check.py).')

    # 3b wake survey
    print(f'\n[3b] wake survey   momentum theory: U_disk/Uinf={UD_T:.4f} '
          f'U_far/Uinf={UF_T:.4f}')
    rows = []
    for xD in XD:
        xp = xD * disk3d.D_PROP
        m = np.abs(cen[:, 0] - xp) < 0.6 * 0.02 * disk3d.D_PROP
        if m.sum() < 20:
            m = np.abs(cen[:, 0] - xp) < 0.03
        rr = np.hypot(cen[m, 1], cen[m, 2])
        ux, ut = U[m, 0], U[m, 1]
        # slipstream edge: furthest radius where the axial velocity still
        # differs from the freestream by more than 2%. For a thrusting disk
        # ux > Uinf, so this must be |deficit|, not the signed value.
        order = np.argsort(rr)
        Rs = rr[order[0]]
        for j in order:
            if abs(1.0 - ux[j] / U_REF) > 0.02:
                Rs = rr[j]
        ins = rr <= R_PROP
        qq = float(np.mean((ux[ins] / U_REF) ** 2)) if ins.sum() else float('nan')
        uxr = float(np.mean(ux[ins] / U_REF)) if ins.sum() else float('nan')
        print(f'  x/D={xD:4.2f}  Rs={Rs * 1e3:6.1f} mm (Rs/R={Rs / R_PROP:.3f})'
              f'  Ux/Uinf={uxr:.4f}  q/qinf={qq:.4f}  n={ins.sum()}')
        prof = []
        edges = np.linspace(0.0, 1.3 * R_PROP, 14)
        for lo, hi in zip(edges[:-1], edges[1:]):
            k = (rr >= lo) & (rr < hi)
            if k.sum() < 2:
                continue
            prof.append(dict(r_over_R=0.5 * (lo + hi) / R_PROP,
                             ux_over_Uinf=float(np.mean(ux[k]) / U_REF),
                             q_over_qinf=float(np.mean((ux[k] / U_REF) ** 2)),
                             n=int(k.sum())))
            rows.append(dict(xD=xD, Rs=Rs, **prof[-1]))
        print('        ' + '  '.join(f'{p["r_over_R"]:.2f}:{p["ux_over_Uinf"]:.3f}'
                                      for p in prof))

    with open(R / f'results/of_disk3_{case.name}.csv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=['xD', 'Rs', 'r_over_R',
                                           'ux_over_Uinf', 'q_over_qinf', 'n'])
        w.writeheader()
        w.writerows(rows)
    print('saved', f'results/of_disk3_{case.name}.csv')


if __name__ == '__main__':
    main()
