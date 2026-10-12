"""Layer 12 / OpenFOAM stage 3a: actuator disk as a VOLUMETRIC momentum source.

Design decision (user-directed): the disk is a thin cellZone carrying a
volumetric momentum source, NOT an internal physical surface. A real actuator
disk is a zero-thickness pressure/momentum jump; modelling it as a wall or
baffle would introduce no-slip, blockage, separation and artificial pressure
drag that have no place in the ROM calibration.

    disk = thin cellZone + volumetric source      (this)
    disk = internal physical surface              (rejected)

This also sidesteps SU2's coincident front/back face requirement entirely.

The cellZone is cut from a plain box mesh with topoSet, so no thin-cylinder
geometry ever enters the mesher - the cylinder exists only as a cell selector.

Closure targets (stage 3a), all cross-checked:
    T_source    = int_Omega S_x dV            (the source term itself)
    T_momentum  = far-field control-volume momentum balance
    T_BEMT      = 16.18 N                     (layer-8 BEMT at J = 0.509)
    |T_source - T_BEMT| / T < 3-5 %
    |T_momentum - T_source| / T < 3-5 %

    P_fluid     = int_Omega S . U dV          (fluid work, not shaft power)

Stage 3 records explicitly whether this is a thrust-only or thrust+power
matched disk, because that determines the wake contraction and therefore k_pw.
"""
import csv
import math
import shutil
import sys
from pathlib import Path

import numpy as np

R = Path('/home/yanlu/Documents/plane')
sys.path.insert(0, str(R / 'foam'))

import ofcase                                              # noqa: E402

# --- layer-8 BEMT reference, cruise point (per engine) --------------------
D_PROP = 0.530
R_PROP = D_PROP / 2.0
A_DISK = math.pi * R_PROP ** 2
J_CRUISE = 0.509
V_CRUISE = 14.07
N_RPS = V_CRUISE / (J_CRUISE * D_PROP)      # 52.16
T_BEMT = 16.1806
P_BEMT = 313.0747
CT_BEMT = 0.061540
CP_BEMT = 0.043076

# disk cellZone thickness: 0.02-0.05 D per the standard volumetric-source
# practice (the source is per-unit-volume, so thickness only sets the smearing)
DISK_THICK = 0.028 * D_PROP                  # 0.01484 m

ALPHAS = [0.0, 4.0, 8.0]


def build_mesh(out_msh, cell_disk=0.0050, cell_near=0.030, cell_far=0.25,
               verbose=False):
    """Box mesh with a refined disk+wake region; no disk geometry."""
    import gmsh
    gmsh.initialize()
    gmsh.option.setNumber('General.Terminal', 1 if verbose else 0)
    gmsh.model.add('disk')
    g = gmsh.model
    geo = g.occ

    x0, x1 = -2.0 * D_PROP, 5.0 * D_PROP
    hw = 2.5 * D_PROP
    geo.addBox(x0, -hw, -hw, x1 - x0, 2 * hw, 2 * hw)
    geo.synchronize()

    gmsh.option.setNumber('Mesh.MeshSizeExtendFromBoundary', 0)
    gmsh.option.setNumber('Mesh.MeshSizeFromPoints', 0)
    gmsh.option.setNumber('Mesh.MeshSizeFromCurvature', 0)

    def box_field(xmin, xmax, h, size):
        f = g.mesh.field.add('Box')
        g.mesh.field.setNumber(f, 'VIn', size)
        g.mesh.field.setNumber(f, 'VOut', 1e30)
        g.mesh.field.setNumber(f, 'XMin', xmin)
        g.mesh.field.setNumber(f, 'XMax', xmax)
        g.mesh.field.setNumber(f, 'YMin', -h)
        g.mesh.field.setNumber(f, 'YMax', h)
        g.mesh.field.setNumber(f, 'ZMin', -h)
        g.mesh.field.setNumber(f, 'ZMax', h)
        return f

    # thin slab for the disk itself, then a coarser slab over the wake
    f_disk = box_field(-0.5 * DISK_THICK, 0.5 * DISK_THICK, 1.15 * R_PROP,
                       cell_disk)
    # The refined near box must be WIDER THAN THE SLIPSTREAM: at x/D = 1 the
    # wake radius is ~1.4 R, so a 1.4 R box leaves the outer wake in 0.25 m
    # coarse cells and any momentum integral over it comes out short. 2.6 R
    # leaves clear margin for expansion.
    f_near = box_field(-0.2 * D_PROP, 3.0 * D_PROP, 2.6 * R_PROP, cell_near)
    f_back = g.mesh.field.add('Box')
    g.mesh.field.setNumber(f_back, 'VIn', cell_far)
    g.mesh.field.setNumber(f_back, 'VOut', 1e30)
    g.mesh.field.setNumber(f_back, 'XMin', x0)
    g.mesh.field.setNumber(f_back, 'XMax', x1)
    g.mesh.field.setNumber(f_back, 'YMin', -hw)
    g.mesh.field.setNumber(f_back, 'YMax', hw)
    g.mesh.field.setNumber(f_back, 'ZMin', -hw)
    g.mesh.field.setNumber(f_back, 'ZMax', hw)
    f_min = g.mesh.field.add('Min')
    g.mesh.field.setNumbers(f_min, 'FieldsList', [f_disk, f_near, f_back])
    g.mesh.field.setAsBackgroundMesh(f_min)

    g.mesh.generate(3)

    # farfield boundary; the disk cellZone is cut later with topoSet
    surf = [t for _, t in g.getEntities(2)]
    g.addPhysicalGroup(2, surf, name='farfield')
    g.addPhysicalGroup(3, [t for _, t in g.getEntities(3)], name='internal')
    geo.synchronize()

    gmsh.option.setNumber('Mesh.MshFileVersion', 2.2)
    gmsh.write(str(out_msh))
    info = dict(nodes=len(g.mesh.getNodes()[0]),
                elems=sum(len(c) for c in g.mesh.getElements(3)[1]))
    gmsh.finalize()
    return info


FVOPTIONS = """\
FoamFile {{ version 2.0; format ascii; class dictionary; object fvOptions; }}
disk1
{{
    type            {dtype};
    active          true;
    selectionMode   cellZone;
    cellZone        propDiskZone;
    fields          (U);

    // surface-normal vector of the disk pointing UPSTREAM. The prop is a
    // tractor, so the freestream arrives along +x and the normal is -x.
    diskDir         (-1 0 0);
    diskArea        {A_DISK};

    // inflow monitoring: the source needs the mean incoming velocity. Place
    // the probe just upstream of the disk (freestream is +x, so upstream is
    // at negative x).
    monitorMethod   points;
    monitorCoeffs
    {{
        points      ((-0.10 0 0));
    }}
    upstreamPoint   (-0.10 0 0);

    // ---- thrust/power coefficients ----
    // OpenFOAM's Froude variant defines a = 1 - Cp/Ct with TURBINE
    // coefficients Cr = T/(0.5 rho A U^2). The BEMT uses the PROPULSOR
    // normalisation Ct = T/(rho n^2 D^4), which is a different quantity:
    //     Ct_propulsor 0.06154  ->  Ct_std = T/(0.5 rho A Uinf^2) = 0.6047
    //     Cp_propulsor 0.043076 ->  Cp_std = P/(0.5 rho A Uinf^3) = 0.8312
    // With the standard values Froude gives a = 1 - Cp/Ct = -0.375 (negative
    // thrust), so Froude is wrong for a propeller; variableScaling uses Ct
    // directly and reproduces T = 16.18 N exactly.
    variant         {variant};
    Cp              {CP};
    Ct              {CT};

    // ---- radial variant: per-radius thrust/power coefficient profiles ----
    // (present only when dtype is radialActuationDiskSource; the profiles are
    // written by a helper because they must match the BEMT stations)
{radial}
    // keep the source as a source, not a sink
    sink            true;
    writeFileStart  0;
    writeFileEnd    100000;
}}
"""


# standard-normalisation coefficients (see the fvOptions note in disk3d.py)
CT_STD = T_BEMT / (0.5 * 1.225 * A_DISK * V_CRUISE ** 2)
CP_STD = P_BEMT / (0.5 * 1.225 * A_DISK * V_CRUISE ** 3)


# U_up the fvOptions monitor reads for this mesh. It is NOT the freestream:
# the probe sits in the disk's induction region and reads ~7% high, so the
# coefficients are calibrated against the measured value rather than U_inf.
# (T_source error was exactly the U_up^2 error until this was accounted for.)
U_UP_MEASURED = 14.9081


def froude_coeffs(a_override=None, u_up=None):
    """Coefficients for variants/Froude.

    Froude computes T = 2 rho A U_up^2 a (1-a) with a = 1 - Cp/Ct. Two caveats
    learned the hard way:

    1. The a(1-a) form is the TURBINE one (decelerating disk). Used on a
       propeller it still accelerates the flow, but the quantitative link
       between the coefficient and the delivered thrust is NOT trustworthy:
       the source was measured to deliver only ~86% of the thrust its own
       coefficient formula implies (13.52 N of 15.77 N).

    2. The upstreamPoint monitor sits INSIDE the induction region and reads
       ~6% above the freestream (14.91 vs 14.07 m/s).

    So this returns a starting guess only; the real calibration is empirical,
    against the independent far-wake momentum integral (see calibrate()).
    """
    u = u_up if u_up is not None else U_UP_MEASURED
    if a_override is not None:
        a = a_override
    else:
        disc = 1.0 - 2.0 * T_BEMT / (1.225 * A_DISK * u ** 2)
        a = 0.5 * (1.0 - math.sqrt(max(disc, 0.0)))
    return 1.0, 1.0 - a, a


def moment_ratio(a):
    """The T/(2 rho A U^2) the Froude formula would need, given induction a."""
    return a * (1.0 - a)


def write_fvoptions(case, dtype='actuationDiskSource', variant='Froude'):
    radial = ''
    ct, cp, a = froude_coeffs()
    (case / 'system/fvOptions').write_text(
        FVOPTIONS.format(dtype=dtype, A_DISK=A_DISK, CP=cp, CT=ct,
                         variant=variant, radial=radial))
    u = U_UP_MEASURED
    print(f'  Froude: Ct={ct:.5f} Cp={cp:.5f} -> a=1-Cp/Ct={a:.4f}  '
          f'T=2 rho A U_up^2 a(1-a) at U_up={u:.4f} = '
          f'{2 * 1.225 * A_DISK * u ** 2 * a * (1 - a):.4f} N '
          f'(BEMT {T_BEMT:.4f} N)', flush=True)


TOPO_DICT = """\
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


def make_case(cell_disk=0.0050, cell_near=0.034, cell_far=0.25, alpha=0.0,
              dtype='actuationDiskSource', rebuild=False, tag='d1'):
    case = R / f'foam/cases/disk_{tag}'
    if rebuild and case.exists():
        shutil.rmtree(case)
    if (case / 'constant/polyMesh/boundary').exists():
        return case
    if case.exists():
        shutil.rmtree(case)
    case.mkdir(parents=True)

    ofcase.U_INF = V_CRUISE
    mesh = R / f'foam/mesh/disk_{tag}.msh'
    mesh.parent.mkdir(parents=True, exist_ok=True)
    print('mesh:', build_mesh(mesh, cell_disk, cell_near, cell_far), flush=True)
    ofcase.write_case(case, mesh,
                      dict(walls=[], symmetry=None, farfield='farfield'),
                      converged_iters=1500, write_interval=1500)
    ofcase.mesh_to_foam(case, mesh, {})
    ofcase.run_cmd(['checkMesh'], case, log=case / 'log.checkMesh')

    (case / 'system/topoSetDict').write_text(
        TOPO_DICT.format(radius=R_PROP, x0=-0.5 * DISK_THICK,
                         x1=0.5 * DISK_THICK))
    ofcase.run_cmd(['topoSet'], case, log=case / 'log.topoSet')
    ofcase.run_cmd(['setsToZones', '-noFlipMap'], case,
                   log=case / 'log.setsToZones')
    write_fvoptions(case, dtype)

    # rotate the freestream for the sweep
    ux = V_CRUISE * math.cos(math.radians(alpha))
    uz = V_CRUISE * math.sin(math.radians(alpha))
    for name in ('U', 'k', 'omega', 'nut', 'p'):
        f = case / '0' / name
        t = f.read_text()
        t = t.replace(f'uniform ({V_CRUISE} 0 0)',
                      f'uniform ({ux:.6f} 0 {uz:.6f})')
        f.write_text(t)
    return case


def main():
    tag = sys.argv[1] if len(sys.argv) > 1 else 'd1'
    case = make_case(tag=tag, rebuild=True)
    print(case, flush=True)
    rc = ofcase.run_cmd(['simpleFoam'], case, log=case / 'log.simpleFoam',
                        timeout=14400)
    print('simpleFoam rc=', rc, flush=True)
    if rc != 0:
        t = (case / 'log.simpleFoam').read_text().splitlines()
        for ln in t[-30:]:
            print('   ', ln[:120])
        return
    ofcase.run_cmd(['foamToVTK', '-latestTime'], case, log=case / 'log.vtk')
    # fvOptions writes its own thrust report
    pp = case / 'postProcessing'
    if pp.exists():
        for f in sorted(pp.rglob('*')):
            if f.is_file():
                print('fvOptions output:', f.relative_to(case))
                print(f.read_text()[-600:])


def calibrate(a_target, iters=900, tag='cal'):
    """Rewrite fvOptions for induction `a_target` and restart the solve.

    Restarting is legitimate only from a WELL-DEVELOPED wake: an earlier
    attempt restarted from a decelerated (wrong-sign) field and 500 iterations
    were not enough to overwrite it. This starts from the converged accelerated
    field, so it only has to adjust the magnitude.
    """
    case = R / f'foam/cases/disk_{tag}'
    if case.exists():
        shutil.rmtree(case)
    shutil.copytree(R / 'foam/cases/disk_d1', case)
    latest = max((d.name for d in case.glob('[0-9]*')
                  if d.is_dir() and d.name != '0'), key=int)
    cd = case / 'system/controlDict'
    t = cd.read_text()
    t = t.replace('startFrom       startTime;', 'startFrom       latestTime;')
    t = t.replace('endTime         1500;',
                  f'endTime         {int(latest) + iters};')
    t = t.replace('writeInterval   1500;', f'writeInterval   {iters};')
    cd.write_text(t)
    ct, cp, a = froude_coeffs(a_override=a_target)
    s2 = (case / 'system/fvOptions').read_text()
    import re as _re
    s2 = _re.sub(r'Ct              [^;]+;', f'Ct              {ct};', s2)
    s2 = _re.sub(r'Cp              [^;]+;', f'Cp              {cp};', s2)
    (case / 'system/fvOptions').write_text(s2)
    print(f'  calibrate: a={a} -> Cp/Ct={cp / ct:.4f}', flush=True)
    rc = ofcase.run_cmd(['simpleFoam'], case, log=case / 'log.simpleFoam',
                        timeout=14400)
    ofcase.run_cmd(['foamToVTK', '-latestTime'], case, log=case / 'log.vtk')
    return case, rc


if __name__ == '__main__':
    if len(sys.argv) > 2 and sys.argv[1] == 'cal':
        case, rc = calibrate(float(sys.argv[2]))
        print('rc=', rc, case)
    else:
        main()
