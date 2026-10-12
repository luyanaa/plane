"""Layer 12 / OpenFOAM stage 3R: rotorDiskSource closure.

Why this exists: the k_pw design point is prototype D's TAKEOFF state,
J_TO = 0.1998 at 5990 rpm, i.e. Ct_std = 7.21 - a heavily loaded disk. At that
loading the swirl is first-order (q_swirl/q_inf ~ 27%, swirl flow angle
~14 deg, comparable to the 20 deg flap deflection), and swirl acts
azimuthally, u_theta(r, psi), so it redistributes the flap's spanwise loading
rather than shifting a scalar incidence. A scalar q_slip/q_inf ROM therefore
lacks a degree of freedom there, and the axial-only actuationDiskSource cannot
supply it.

rotorDiskSource is the right tool: it is a cell-based momentum source over a
cylindrical region that computes the blade force from LOCAL blade kinematics
(local axial and tangential velocity + blade twist), producing both axial and
tangential (swirl) loading.

Formulation, per the user's instruction:
    fixed blade geometry + fixed RPM + local inflow
NOT `targetForce` trim - that would tune the blade characteristics to hit a
target thrust and would erase the very difference between the Layer-8 and
rotorDisk blade-element formulations that this stage exists to test.

Layer-8 blade geometry is mapped verbatim:
    blade  data ( (profile1 (radius[m] twist[deg] chord[m])) ... )
    profile data ( (AOA[deg] Cd Cl) ... )      <- Cd BEFORE Cl
"""
import math
import shutil
import sys
from pathlib import Path

import numpy as np

R = Path('/home/yanlu/Documents/plane')
sys.path.insert(0, str(R))
sys.path.insert(0, str(R / 'foam'))

import disk3d as disk_d                                    # noqa: E402
import ofcase                                              # noqa: E402

from sizing.prop_bemt import PropBEMT                      # noqa: E402

D_PROP = 0.530
R_PROP = D_PROP / 2.0
RHO = 1.225
N_BLADES = 2
P_OVER_D = 1.0
TIP_EFFECT = 0.98
U_REF = 14.07

# Prototype D operating points, from scripts/prop_points.py. Each point has
# its OWN speed and rpm - takeoff is at v_r = 1.1 V_sTO = 16.62 m/s and
# 9420 rpm, NOT the 14.07 m/s / 132.9 rps used in an earlier draft.
POINTS_FILE = R / 'results/prop_points_D.json'


def load_points():
    import json
    return json.loads(POINTS_FILE.read_text())


def layer8_blade(n_stations=26):
    """Radial blade geometry from the Layer-8 BEMT, mapped verbatim."""
    p = PropBEMT(D=D_PROP, P_over_D=P_OVER_D)
    mu = np.linspace(p.r_hub_frac + 0.01, 0.99, n_stations)
    rows = [(float(m * p.R), float(p.twist(m)), float(p.chord(m))) for m in mu]
    return p, rows


def layer8_profile(n_lin=31, alpha_max=30.0, alpha_clamp=90.0):
    """(AOA deg, Cd, Cl) table from the Layer-8 linear polar.

    Cl = cl_alpha (alpha - alpha0),  Cd = cd0 + cd_k Cl^2
    The linear region is reproduced exactly; beyond +-alpha_max the last value
    is held out to the clamp angle, which is standard practice and keeps the
    solver stable without inventing polar behaviour Layer-8 does not have.
    """
    p = PropBEMT(D=D_PROP, P_over_D=P_OVER_D)
    alphas = np.linspace(-alpha_max, alpha_max, n_lin)
    rows = []
    for a in alphas:
        cl = p.cl_alpha * math.radians(a - p.alpha0_deg)
        cd = p.cd0 + p.cd_k * cl * cl
        rows.append((float(a), float(cd), float(cl)))
    # clamp ends
    rows = ([(float(-alpha_clamp), rows[0][1], rows[0][2])] + rows +
            [(float(alpha_clamp), rows[-1][1], rows[-1][2])])
    return rows


FVOPT_HEAD = """\
FoamFile {{ version 2.0; format ascii; class dictionary; object fvOptions; }}
rotorDisk
{{
    type            rotorDisk;
    active          true;
    selectionMode   cellZone;
    cellZone        propDiskZone;
    fields          (U);

    nBlades         {nblades};
    tipEffect       {tip};

    // local inflow: the blade force comes from the LOCAL flow the solver sees,
    // not from a prescribed freestream, which is what makes this a genuine
    // independent blade-element model rather than a re-labelled disk
    inletFlowType   local;

    geometryMode    specified;
    origin          (0 0 0);
    axis            (1 0 0);
    refDirection    (0 0 1);
    rpm             {rpm};

    // fixed blade geometry, NOT targetForce: targetForce would tune the blade
    // to hit a thrust target and hide any formulation difference from Layer-8
    trimModel       fixedTrim;
    rhoRef          {rho};
    rhoInf          {rho};

    fixedTrimCoeffs
    {{
        theta0      0;
        theta1c     0;
        theta1s     0;
    }}

    flapCoeffs
    {{
        beta0       0;
        beta1c      0;
        beta1s      0;
        beta2s      0;
    }}

    blade
    {{
        data
        (
{blade}
        );
    }}

    profiles
    {{
        p1
        {{
            type    lookup;
            data
            (
{profile}
            );
        }}
    }}
}}
"""


def write_rotor_fvoptions(case, rpm):
    _, rows = layer8_blade()
    prof = layer8_profile()
    blade = '\n'.join('            (p1 (%.5f %.4f %.6f))' % (r, t, c)
                      for r, t, c in rows)
    profile = '\n'.join('                (%.1f %.6f %.6f)' % (a, cd, cl)
                        for a, cd, cl in prof)
    (case / 'system/fvOptions').write_text(
        FVOPT_HEAD.format(nblades=N_BLADES, tip=TIP_EFFECT, rpm=rpm, rho=RHO,
                          blade=blade, profile=profile))


# rotorDiskSource distributes blade force over the cellZone volume, so the
# zone must be at least as thick as the blade chord or the force is
# over-concentrated and the explicit source destabilises the solve. The blade
# chord here is 0.024-0.042 m, so 0.12 m (0.23 D) gives ~3x chord.
ZONE_HALF_THICK = 0.06

TOPOSET = """\
FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      topoSetDict;
}}

actions
(
    {{
        name    propDiskZone;
        type    cellSet;
        action  new;
        source  cylinderToCell;
        origin  (0 0 0);
        axis    (1 0 0);
        radius  {radius};
        p1      ({x0} 0 0);
        p2      ({x1} 0 0);
    }}
);
"""


def recut_zone(case, half=ZONE_HALF_THICK):
    """Re-cut propDiskZone thicker than the blade chord (see note above)."""
    (case / 'system/topoSetDict').write_text(
        TOPOSET.format(radius=disk_d.R_PROP, x0=-half, x1=half))
    ofcase.run_cmd(['topoSet'], case, log=case / 'log.topoSet')
    ofcase.run_cmd(['setsToZones', '-noFlipMap'], case,
                   log=case / 'log.setsToZones')


def set_relax(case, u=0.5, kom=0.5):
    """A strong explicit momentum source needs heavier under-relaxation than
    the default 0.9, otherwise the explicit coupling diverges."""
    f = case / 'system/fvSolution'
    t = f.read_text()
    import re as _re
    t = _re.sub(r'(?m)^        U               [0-9.]+;',
                f'        U               {u};', t)
    t = _re.sub(r'(?m)^        "\(k\|omega\)"     [0-9.]+;',
                f'        "(k|omega)"     {kom};', t)
    f.write_text(t)


def make_case(point, tag=None, base=R / 'foam/cases/disk_d1'):
    """Reuse the validated disk mesh; only the source and freestream change."""
    pts = load_points()
    if point not in pts:
        raise SystemExit(f'unknown point {point}; have {list(pts)}')
    P = pts[point]
    tag = tag or point
    case = R / f'foam/cases/rotor_{tag}'
    if case.exists():
        shutil.rmtree(case)
    shutil.copytree(base, case)
    for d in case.glob('[0-9]*'):
        if d.is_dir() and d.name != '0':
            shutil.rmtree(d)
    for f in case.glob('log.*'):
        f.unlink()

    cd = case / 'system/controlDict'
    t = cd.read_text()
    t = t.replace('endTime         1500;', 'endTime         2500;')
    t = t.replace('writeInterval   1500;', 'writeInterval   2500;')
    cd.write_text(t)

    ofcase.U_INF = P['V']
    for name in ('U', 'k', 'omega', 'nut', 'p'):
        f = case / '0' / name
        txt = f.read_text()
        txt = txt.replace('uniform (14.07 0 0)', f'uniform ({P["V"]:.6f} 0 0)')
        txt = txt.replace('freestreamValue uniform (14.07 0 0);',
                          f'freestreamValue uniform ({P["V"]:.6f} 0 0);')
        f.write_text(txt)
    import re as _re
    k_inf = 1.5 * (P['V'] * 1e-3) ** 2
    om_inf = k_inf ** 0.5 / (0.09 ** 0.25 * 0.1 * 0.6788)
    for name, val in (('k', k_inf), ('omega', om_inf)):
        f = case / '0' / name
        txt = _re.sub(r'uniform [0-9.eE+-]+;', f'uniform {val};',
                      f.read_text(), count=1)
        f.write_text(txt)

    recut_zone(case)
    set_relax(case)
    write_rotor_fvoptions(case, P['n_rpm'])
    return case, P


def main():
    point = sys.argv[1] if len(sys.argv) > 1 else 'cruise'
    case, P = make_case(point)
    print(f'[{point}] V={P["V"]:.3f} m/s  n={P["n_rpm"]:.1f} rpm  '
          f'J={P["J"]:.4f}  BEMT T={P["T"]:.2f} N  P={P["P"]:.1f} W  '
          f'Ct_std={P["Ct_std"]:.3f}  -> {case.name}', flush=True)
    rc = ofcase.run_cmd(['simpleFoam'], case, log=case / 'log.simpleFoam',
                        timeout=14400)
    print(f'[{point}] simpleFoam rc={rc}', flush=True)
    if rc != 0:
        for ln in (case / 'log.simpleFoam').read_text().splitlines()[-25:]:
            print('   ', ln[:130])
        return
    ofcase.run_cmd(['foamToVTK', '-latestTime'], case, log=case / 'log.vtk')


if __name__ == '__main__':
    main()
