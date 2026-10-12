"""Layer 12 step 3b: disk thrust from the tutorial mesh volume field.

Recovers the thrust that SU2 actually applied by measuring the static pressure
jump immediately across the coincident disk faces:

    T = (p_down - p_up) * A_disk,     A_disk = pi R^2 (1 - r_bar_in^2)

This is exact by actuator-disk theory (the static pressure rises across a
propulsive disk, p2 - p1 = 2 rho V^2 a (1+a)), and it avoids both the
per-marker force breakdown - which is dominated by a momentum-flux term and is
not a thrust measure - and the internal-marker surface output, which does not
expose the BC-imposed face states.

The tutorial domain applies a back pressure at OUT (56370 Pa against a 101325
Pa freestream), so a far-field control-volume balance is not clean either; the
local face measurement is the right tool.
"""
import sys
from pathlib import Path

import numpy as np

R = Path('/home/yanlu/Documents/plane')
sys.path.insert(0, str(R))
sys.path.insert(0, str(R / 'su2'))

TUT = R / 'su2/tutorial'
MESH_R = 2.5146
R_BAR_IN = 0.1804
V_INF = 0.55996 * np.sqrt(1.4 * 101325.0 / 1.225)      # tutorial Mach
RHO = 1.225


def read_su2_vtu(path):
    """Parse SU2's appended-format VTU: returns (points Nx3, pressure N).

    Reads the XML header for offsets and the binary block directly, bypassing
    meshio (SU2 writes an Eddy_Viscosity array shorter than the point count,
    which meshio treats as fatal).
    """
    import re
    import struct

    raw = open(path, 'rb').read()
    txt = raw.decode('latin-1')
    header_end = txt.find('<AppendedData')
    opening = txt.find('>', header_end) + 1
    payload = raw[opening:]
    start = payload.find(b'_')
    data = payload[start + 1:]

    offsets = dict(re.findall(
        r'<DataArray type="(Float32|Int32|Int64|UInt64)"[^>]*Name="([^"]*)"'
        r'[^>]*offset="(\d+)"', txt) and [])
    arrays = {}
    for m in re.finditer(r'<DataArray[^>]*>', txt):
        tag = m.group(0)
        name = re.search(r'Name="([^"]*)"', tag)
        off = re.search(r'offset="(\d+)"', tag)
        typ = re.search(r'type="([^"]*)"', tag)
        ncomp = re.search(r'NumberOfComponents=\s*"?(\d+)"?', tag)
        if not (off and typ):
            continue
        arrays.setdefault(name.group(1) if name else '', []).append(
            (int(off.group(1)), typ.group(1),
             int(ncomp.group(1)) if ncomp else 1))

    npts = int(re.search(r'NumberOfPoints="(\d+)"', txt).group(1))

    def read_arr(off, typ, nc, n):
        dt = {'Float32': 'f', 'Float64': 'd', 'Int32': 'i',
              'Int64': 'q', 'UInt64': 'Q'}[typ]
        sz = struct.calcsize(dt)
        # each appended array is prefixed by its byte length (UInt64)
        blen = struct.unpack_from('<Q', data, off)[0]
        body = data[off + 8:off + 8 + blen]
        vals = struct.unpack_from('<' + dt * (blen // sz), body, 0)
        return np.array(vals).reshape(n, nc)

    pts = None
    ncell = int(re.search(r'NumberOfCells="(\d+)"', txt).group(1))
    for key, items in arrays.items():
        for off, typ, nc in items:
            if key == '' and nc == 3 and pts is None:
                pts = read_arr(off, typ, nc, npts)
            elif key == 'Pressure':
                p = read_arr(off, typ, nc, npts)
                return pts, p.ravel()
    raise RuntimeError('points or Pressure not found in VTU')


def main():
    import numpy as np
    f = TUT / 'flow_l8.vtu'
    # SU2's VTU has one fewer Eddy_Viscosity entry than points, which meshio
    # rejects; read the raw appended binary ourselves.
    pts, p = read_su2_vtu(str(f))
    print(f'volume: {len(pts)} points, pressure {p.min():.0f}..{p.max():.0f} Pa')

    r = np.hypot(pts[:, 1], pts[:, 2])
    r_bar = r / MESH_R
    in_disk = (r_bar > R_BAR_IN) & (r_bar < 1.02)

    # nodes immediately either side of the disk plane
    x = pts[:, 0]
    near = np.abs(x) < 0.02 * MESH_R
    up = near & (x < 0) & in_disk
    dn = near & (x > 0) & in_disk
    print(f'disk-adjacent nodes: upstream {up.sum()}  downstream {dn.sum()}')
    if up.sum() < 3 or dn.sum() < 3:
        # widen the band
        near = np.abs(x) < 0.08 * MESH_R
        up = near & (x < 0) & in_disk
        dn = near & (x > 0) & in_disk
        print(f'  widened: upstream {up.sum()}  downstream {dn.sum()}')

    a_disk = np.pi * MESH_R ** 2 * (1.0 - R_BAR_IN ** 2)

    # Node counts differ upstream/downstream (608 vs 852), and the pressure
    # field varies radially, so a plain mean is biased. Match the two sides in
    # radial bins and integrate per bin, which is also the correct
    # discretisation of T = int (p_dn - p_up) dA.
    edges = np.linspace(R_BAR_IN, 1.0, 21)
    t_bin = 0.0
    used = 0
    print('  rbar_lo rbar_hi   n_up  n_dn    p_up      p_dn     dp      dA')
    for lo, hi in zip(edges[:-1], edges[1:]):
        mu = up & (r_bar >= lo) & (r_bar < hi)
        md = dn & (r_bar >= lo) & (r_bar < hi)
        if mu.sum() < 1 or md.sum() < 1:
            continue
        d_a = np.pi * MESH_R ** 2 * (hi ** 2 - lo ** 2)
        dp = p[md].mean() - p[mu].mean()
        t_bin += dp * d_a
        used += 1
        if used <= 8 or abs(dp) > 2000:
            print(f'  {lo:.4f} {hi:.4f}   {mu.sum():4d}  {md.sum():4d}  '
                  f'{p[mu].mean():8.1f} {p[md].mean():8.1f} {dp:+8.1f} {d_a:6.3f}')
    print(f'  ({used} bins integrated of {len(edges) - 1})')

    t_cfd = t_bin
    print(f'A_disk      = {a_disk:.3f} m2  (rbar_in {R_BAR_IN}, R {MESH_R})')
    print(f'T_CFD (binned) = {t_cfd:.1f} N')
    print(f'T_CFD (plain)  = {(p[dn].mean() - p[up].mean()) * a_disk:.1f} N')

    # what the input file asked for
    ct_file = float(sys.argv[1]) if len(sys.argv) > 1 else 0.08064
    
    print(f'input C_T   = {ct_file:.5f}')
    print(f'T_input (Renard, mesh radius, from the BEMT file) requires the '
          f'n and D SU2 used; see disk_closure.py output.')
    print(f'expected jump = T/A_disk; if the closure holds, '
          f'T_CFD/T_input -> 1')


if __name__ == '__main__':
    main()
