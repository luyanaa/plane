"""Post-processing shared by the body-fitted 3-D wing drivers.

Extracted from the abandoned snappyHexMesh route (see README) so the working
body-fitted drivers do not import a dead module.
"""
import math
import sys
from pathlib import Path

import numpy as np

R = Path('/home/yanlu/Documents/plane')
sys.path.insert(0, str(R / 'foam'))

import ofcase                                              # noqa: E402

S_HALF = 2.1440              # half-wing reference area
B_HALF = 3.1585
C_ROOT = 0.8310
TAPER = 0.60


def chord_at(y):
    return C_ROOT * (1.0 - (1.0 - TAPER) * y / B_HALF)


def latest_vtp(case):
    case = Path(case)
    dirs = [p for p in (case / 'VTK').glob(f'{case.name}_*') if p.is_dir()]
    if not dirs:
        return None
    v = sorted(dirs)[-1] / 'boundary' / 'wing.vtp'
    return v if v.exists() else None


def alpha_result(case, alpha_deg):
    """CL, CD from our own pressure integration of the wing patch."""
    import forces as F
    vtp = latest_vtp(case)
    if not vtp:
        return None
    f, _ = F.patch_forces(str(vtp), ofcase.RHO)
    a = math.radians(alpha_deg)
    drag = f[0] * math.cos(a) + f[2] * math.sin(a)
    lift = -f[0] * math.sin(a) + f[2] * math.cos(a)
    q = 0.5 * ofcase.RHO * ofcase.U_INF ** 2 * S_HALF
    return dict(alpha=alpha_deg, CL=lift / q, CD=drag / q, vtp=str(vtp))


def span_loading(vtp, alpha_deg, u_inf=None):
    """L'(y)/L and c_l(y) from the wing patch pressure.

    Faces sit on DISCRETE span stations (the extrusion divisions), so uniform
    y bins alternately capture and miss whole stations and produce a spurious
    zig-zag. Group by the actual station instead.
    """
    import forces as F
    d = F.read_vtk_xml(str(vtp))
    pts = d['points']
    av = F.poly_area_normals(pts, d['polys'])
    p = F.face_pressure(d, d['polys'])
    f = (p * ofcase.RHO)[:, None] * av

    yc = np.array([pts[poly][:, 1].mean() for poly in d['polys']])
    a = math.radians(alpha_deg)
    lift_face = -f[:, 0] * math.sin(a) + f[:, 2] * math.cos(a)

    q = 0.5 * ofcase.RHO * ofcase.U_INF ** 2
    key = np.round(yc / 1e-4).astype(int)
    uk = np.unique(key)
    ys = np.array([yc[key == k].mean() for k in uk])
    Ls = np.array([lift_face[key == k].sum() for k in uk])
    out = []
    for i in range(len(uk)):
        lo = 0.5 * (ys[i - 1] + ys[i]) if i > 0 else ys[0]
        hi = 0.5 * (ys[i] + ys[i + 1]) if i < len(uk) - 1 else ys[-1]
        dy = max(hi - lo, 1e-6)
        out.append(dict(y=float(ys[i]), Lprime=float(Ls[i] / dy),
                        cl=float((Ls[i] / dy) / (q * chord_at(ys[i]))),
                        n=int((key == uk[i]).sum())))
    if out:
        tot = float(np.trapezoid([o['Lprime'] for o in out],
                                 [o['y'] for o in out]))
        for o in out:
            o['Ln'] = o['Lprime'] / tot if tot else 0.0
    return out
