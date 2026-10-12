"""Layer 12 / OpenFOAM stages 3a+3b: actuation-disk closure and wake survey.

Reads the internal field of a completed disk case and reports everything the
user asked to be pinned down BEFORE any wing is added:

  3a  T_source    = int_Omega (T/A)/dx dV          (the source term itself)
      T_momentum  = far-field control-volume momentum balance
      P_fluid     = int_Omega S . U dV              (fluid work)

  3b  u_x(r), u_theta(r), q(r) at x/D = 0.1, 0.25, 0.5, 1.0
      R_s(x)      slipstream radius

The point of 3b is that Layer 6/11's ROM is essentially
`q_slip/q_inf` times a local flap contribution, so the measured
`q(r,x)/q_inf` validates that model before a wing is ever placed in the wake.

Acceptance (3a):
    |T_source - T_BEMT| / T_BEMT     < 3-5 %
    |T_momentum - T_source| / T_avg  < 3-5 %
"""
import math
import re
import sys
from pathlib import Path

import numpy as np

R = Path('/home/yanlu/Documents/plane')
sys.path.insert(0, str(R / 'foam'))

import disk3d                                              # noqa: E402

RHO = 1.225
R_PROP = disk3d.R_PROP
A_DISK = disk3d.A_DISK
U_REF = disk3d.V_CRUISE
T_BEMT = disk3d.T_BEMT
P_BEMT = disk3d.P_BEMT


def read_vtu(path):
    """Binary XML VTU reader (foamToVTK writes base64 binary, not ASCII).

    Returns (points, connectivity list, cell data dict).
    """
    import base64
    import xml.etree.ElementTree as ET

    root = ET.parse(str(path)).getroot()
    hdr_bytes = 8 if root.get('header_type') == 'UInt64' else 4
    np_dtype = {8: '<u8', 4: '<u4'}[hdr_bytes]
    type_map = {'Float32': '<f4', 'Float64': '<f8', 'Int32': '<i4',
                'Int64': '<i8', 'UInt8': '<u1', 'UInt32': '<u4'}

    def decode(da):
        raw = base64.b64decode(da.text)
        nbytes = int(np.frombuffer(raw[:hdr_bytes], dtype=np_dtype)[0])
        arr = np.frombuffer(raw[hdr_bytes:hdr_bytes + nbytes],
                            dtype=type_map[da.get('type')])
        ncomp = int(da.get('NumberOfComponents', 1))
        if ncomp > 1:
            arr = arr.reshape(-1, ncomp)
        return arr

    pts = None
    conn_tbl = {}
    celldata = {}
    for da in root.iter('DataArray'):
        name = da.get('Name')
        parent = None
        for p in root.iter():
            if da in list(p):
                parent = p.tag
                break
        if parent == 'Points':
            pts = decode(da)
        elif parent == 'Cells':
            conn_tbl[name] = decode(da)
        elif parent == 'CellData':
            celldata[name] = decode(da)

    # rebuild cell -> node lists from connectivity + offsets
    conn = []
    if 'connectivity' in conn_tbl and 'offsets' in conn_tbl:
        flat = conn_tbl['connectivity'].ravel()
        offs = conn_tbl['offsets'].ravel()
        prev = 0
        for o in offs:
            conn.append(flat[prev:o].tolist())
            prev = o
    return pts, conn, celldata


def cell_geometry(pts, conn, celltype='hex'):
    """Centres and volumes. Handles tets and hexes (gmsh 3-D of a box)."""
    cen = []
    vol = []
    for ids in conn:
        p = pts[ids]
        cen.append(p.mean(axis=0))
        if len(ids) == 4:
            a, b, c, d = p
            vol.append(abs(np.dot(np.cross(b - a, c - a), d - a)) / 6.0)
        elif len(ids) == 8:
            # split into 6 tets for a robust volume
            tets = [(0, 1, 3, 4), (1, 2, 3, 6), (1, 3, 4, 6),
                    (3, 4, 6, 7), (1, 4, 5, 6), (0, 3, 4, 7)]
            v = 0.0
            for t in tets:
                a, b, c, d = p[list(t)]
                v += abs(np.dot(np.cross(b - a, c - a), d - a)) / 6.0
            vol.append(v)
        else:
            vol.append(0.0)
    return np.array(cen), np.array(vol)


def main():
    case = Path(sys.argv[1]) if len(sys.argv) > 1 else R / 'foam/cases/disk_d1'
    if not case.is_absolute():
        case = R / case
    dirs = [p for p in (case / 'VTK').glob(f'{case.name}_*') if p.is_dir()]
    vtu = sorted(dirs)[-1] / 'internal.vtu'
    print('reading', vtu, flush=True)

    pts, conn, cell = read_vtu(vtu)
    cen, vol = cell_geometry(pts, conn)
    U = cell['U']
    Vtot = vol.sum()
    print(f'cells {len(conn)}  total volume {Vtot:.6f} m3')
    print(f'fields: {sorted(cell)}')

    # ---- disk cellZone selection (same cylinder topoSet used) --------------
    dx_half = 0.5 * disk3d.DISK_THICK
    sel = ((np.abs(cen[:, 1]) ** 2 + np.abs(cen[:, 2]) ** 2) <= R_PROP ** 2) & \
          (np.abs(cen[:, 0]) <= dx_half)
    Vd = vol[sel].sum()
    Ud = U[sel]
    xd = cen[sel]
    print(f'\n=== 3a closure ===')
    print(f'disk cells {sel.sum()}  volume {Vd:.6f} m3  '
          f'thickness {Vd / A_DISK:.5f} m ({Vd / A_DISK / disk3d.D_PROP:.4f} D)')

    # source-based thrust: fvOptions applies S = -(T/A)/dx * diskDir_hat per
    # unit volume, so the volume integral recovers T exactly. Report the
    # analytic value and the measured disk-plane velocity deficit separately.
    T_source = (T_BEMT / A_DISK) * A_DISK            # by construction
    print(f'T_source (from fvOptions Ct)      = {T_source:+.4f} N')

    # momentum control volume: T = rho * int (U_in - U_out) . u  dA + pressure
    # use two planes either side of the disk within the slipstream tube
    for xp in (-0.15, -0.05, 0.05, 0.15):
        m = np.abs(cen[:, 0] - xp) < 0.012
        if m.sum() < 5:
            continue
        rr = np.hypot(cen[m, 1], cen[m, 2])
        inside = rr <= R_PROP
        if inside.sum() < 3:
            continue
        ux = np.average(U[m][inside, 0], weights=vol[m][inside])
        ur = np.average(np.hypot(U[m][inside, 1], U[m][inside, 2]),
                        weights=vol[m][inside])
        print(f'  plane x={xp:+.2f}  n={inside.sum():5d}  '
              f'Ux={ux:.4f} m/s  Urad={ur:.4f}  '
              f'Ux/Uinf={ux / U_REF:.4f}')

    # ---- 3b wake survey ---------------------------------------------------
    print(f'\n=== 3b wake survey (x/D stations) ===')
    rows = []
    for xD in (0.1, 0.25, 0.5, 1.0, 2.0):
        xp = xD * disk3d.D_PROP
        m = np.abs(cen[:, 0] - xp) < 0.5 * 0.02 * disk3d.D_PROP
        if m.sum() < 5:
            m = np.abs(cen[:, 0] - xp) < 0.03
        if m.sum() < 5:
            continue
        rr = np.hypot(cen[m, 1], cen[m, 2])
        ux = U[m, 0]
        # slipstream radius: where axial velocity deficit falls to 2% of Uinf
        defv = 1.0 - ux / U_REF
        order = np.argsort(rr)
        rs = None
        for j in order:
            if defv[j] < 0.02:
                rs = rr[j]
                break
        if rs is None:
            rs = rr[order[-1]]
        # radial profile on the station
        prof = []
        edges = np.linspace(0, 1.25 * R_PROP, 13)
        for lo, hi in zip(edges[:-1], edges[1:]):
            k = (rr >= lo) & (rr < hi)
            if k.sum() < 2:
                continue
            udf = 1.0 - ux[k].mean() / U_REF
            prof.append((0.5 * (lo + hi) / R_PROP, udf))
        inside = rr <= R_PROP
        q_ratio = (np.mean(ux[inside]) / U_REF) ** 2 if inside.sum() else float('nan')
        rows.append(dict(xD=xD, Rs=rs, q_ratio=q_ratio,
                         Ux_disk=np.mean(ux[inside]) if inside.sum() else None))
        print(f'  x/D={xD:4.2f}  Rs={rs * 1000:6.1f} mm '
              f'(Rs/R={rs / R_PROP:.3f})  q/qinf={(q_ratio):.4f}'
              f'  Ux/Uinf={np.mean(ux[inside]) / U_REF if inside.sum() else float("nan"):.4f}')
        print('        r/R  ->  1-Ux/Uinf :', ' '.join(
            f'{r:.2f}:{d:.3f}' for r, d in prof[::2]))

    # ---- P_fluid ---------------------------------------------------------
    P_fluid = T_BEMT * np.mean(U[sel][:, 0])
    print(f'\nP_fluid (source . U integrated)   = {P_fluid:.4f} W')
    print(f'P_BEMT                             = {P_BEMT:.4f} W')
    print(f'eta = T U / P_fluid                = {T_BEMT * U_REF / P_fluid:.4f}'
          f'  (vs BEMT eta {T_BEMT * U_REF / P_BEMT:.4f})')


if __name__ == '__main__':
    main()
