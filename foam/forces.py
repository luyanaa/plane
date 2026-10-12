"""Own force integration for OpenFOAM results (no function objects needed).

This Ubuntu OpenFOAM v1912 package ships an incomplete runtime layout
(FOAM_ETC and FOAM_LIBBIN are empty and /usr/share/openfoam/bin is absent), so
every function object - `forces`, `yPlus` - falls back to the dynamicCode path
and dies writing a `sha1` stream. The solver, gmshToFoam and foamToVTK all
work, so forces are integrated here from the exported patch surfaces.

    F = sum_faces  p_face * A_face * n_face

The patch normal points outward from the fluid, i.e. into the body, which is
exactly the direction in which the fluid pushes the wall. simpleFoam stores
kinematic pressure (p/rho), hence the density factor.
"""
import base64
import re
import struct

import numpy as np

_DT = {'Float32': ('f', 4), 'Float64': ('d', 8),
       'Int32': ('i', 4), 'Int64': ('q', 8),
       'UInt32': ('I', 4), 'UInt64': ('Q', 8),
       'UInt8': ('B', 1), 'Int8': ('b', 1)}


def _decode(block_text, typ):
    body = block_text.strip()
    buf = base64.b64decode(body)
    dt, sz = _DT[typ]
    if len(buf) < 8:
        return np.zeros(0)
    n = (len(buf) - 8) // sz
    return np.array(struct.unpack_from('<' + dt * n, buf, 8))


def read_vtk_xml(path):
    """Parse a VTK XML file (VTP or VTU) with base64 binary arrays.

    Returns points, polys, and the point/cell data dictionaries. Sections are
    located by absolute offsets so the point and cell blocks cannot be
    confused.
    """
    txt = open(path, 'rb').read().decode('latin-1')

    def section(name):
        m = re.search(rf'<{name}>(.*?)</{name}>', txt, re.S)
        return (m.start(1), m.end(1)) if m else None

    pd_span = section('PointData')
    cd_span = section('CellData')

    point_data, cell_data = {}, {}
    for m in re.finditer(r'<DataArray[^>]*>.*?</DataArray>', txt, re.S):
        if pd_span and pd_span[0] <= m.start() < pd_span[1]:
            dst = point_data
        elif cd_span and cd_span[0] <= m.start() < cd_span[1]:
            dst = cell_data
        else:
            continue
        head = m.group(0)[:m.group(0).find('>')]
        name = re.search(r"Name='([^']*)'", head)
        typ = re.search(r"type='([^']*)'", head)
        ncomp = re.search(r"NumberOfComponents='(\d+)'", head)
        if not (name and typ):
            continue
        body = m.group(0)[m.group(0).find('>') + 1:
                          m.group(0).rfind('</DataArray>')]
        vals = _decode(body, typ.group(1))
        nc = int(ncomp.group(1)) if ncomp else 1
        dst[name.group(1)] = vals.reshape(-1, nc) if nc > 1 else vals

    pts = None
    mp = re.search(r'<Points>.*?</Points>', txt, re.S)
    if mp:
        for m in re.finditer(r'<DataArray[^>]*>.*?</DataArray>',
                             mp.group(0), re.S):
            head = m.group(0)[:m.group(0).find('>')]
            typ = re.search(r"type='([^']*)'", head)
            body = m.group(0)[m.group(0).find('>') + 1:
                              m.group(0).rfind('</DataArray>')]
            vals = _decode(body, typ.group(1))
            pts = vals.reshape(-1, 3)
            break

    def conn_offsets():
        c = o = None
        for m in re.finditer(r'<DataArray[^>]*>.*?</DataArray>', txt, re.S):
            head = m.group(0)[:m.group(0).find('>')]
            name = re.search(r"Name='([^']*)'", head)
            typ = re.search(r"type='([^']*)'", head)
            if not (name and typ):
                continue
            body = m.group(0)[m.group(0).find('>') + 1:
                              m.group(0).rfind('</DataArray>')]
            if name.group(1) == 'connectivity':
                c = _decode(body, typ.group(1)).astype(int)
            elif name.group(1) == 'offsets':
                o = _decode(body, typ.group(1)).astype(int)
        return c, o

    c, o = conn_offsets()
    polys = None
    if c is not None and o is not None:
        polys = [c[o[i - 1] if i else 0:o[i]] for i in range(len(o))]

    return dict(points=pts, polys=polys, point_data=point_data,
                cell_data=cell_data)


def poly_area_normals(pts, polys):
    """Newell area vectors: magnitude = polygon area, direction = winding."""
    out = np.empty((len(polys), 3))
    for i, poly in enumerate(polys):
        p = pts[poly]
        n = np.zeros(3)
        for j in range(len(p)):
            n += np.cross(p[j], p[(j + 1) % len(p)])
        out[i] = 0.5 * n
    return out


def face_pressure(d, polys):
    """Pressure per face, from point data if present else cell data."""
    for src in (d['point_data'], d['cell_data']):
        if 'p' in src:
            p = np.asarray(src['p'])
            p = p[:, 0] if p.ndim > 1 else p
            if len(p) == len(d['points']):
                return np.array([p[poly].mean() for poly in polys])
            return p
    raise RuntimeError('no pressure field in the patch file')


def patch_forces(vtp, rho):
    """Pressure force on the body from one wall patch (kinematic p -> *rho)."""
    d = read_vtk_xml(vtp)
    av = poly_area_normals(d['points'], d['polys'])
    p = face_pressure(d, d['polys'])
    f = (p * rho)[:, None] * av
    return f.sum(axis=0), d


def coefs(f, rho, u_inf, chord, alpha_deg):
    """Wind-axis CL, CD for a 2-D case (reference area = chord x 1 m)."""
    import math
    a = math.radians(alpha_deg)
    drag = f[0] * math.cos(a) + f[2] * math.sin(a)
    lift = -f[0] * math.sin(a) + f[2] * math.cos(a)
    q = 0.5 * rho * u_inf * u_inf * chord
    return lift / q, drag / q
