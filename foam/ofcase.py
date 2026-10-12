"""OpenFOAM case driver (layer 12, second and independent closure route).

OpenFOAM is the 3-D viscous / boundary-layer / blown-flap route; SU2 stays the
geometry, inviscid and interface-regression route. The decisive advantage here
is that OpenFOAM accepts general polyhedra, so the prism and coincident-face
constraints that make SU2 actuator disks painful do not exist.

This module builds an extruded O-grid wing mesh in gmsh (2-D BoundaryLayer
field, which is documented to work in 2-D, then recombined into hexes/prisms),
writes gmsh's own MSH 2.2, converts with gmshToFoam, and runs a steady
incompressible RANS solve.
"""
import os
import shutil
from math import cos, pi, sin
import subprocess
from pathlib import Path

R = Path('/home/yanlu/Documents/plane')
FOAM = Path('/home/yanlu/Documents/plane/foam')
OPENFOAM_BASHRC = '/usr/share/openfoam/etc/bashrc'

# --- flow state (layer-11 TO2, sea level) ---------------------------------
U_INF = 14.07           # m/s
NU = 1.46e-5            # m^2/s
RHO = 1.225
CHORD = 0.6788
SPAN = 0.5 * 6.317      # half span
RE = U_INF * CHORD / NU


def foam_env():
    """Environment for OpenFOAM binaries."""
    env = dict(os.environ)
    out = subprocess.run(
        ['bash', '-lc', f'source {OPENFOAM_BASHRC} && env'],
        capture_output=True, text=True)
    for line in out.stdout.splitlines():
        if '=' in line:
            k, v = line.split('=', 1)
            env[k] = v
    return env


def run_cmd(cmd, cwd, log=None, timeout=7200):
    env = foam_env()
    # OpenFOAM resolves the case from PWD. Inheriting PWD from the login shell
    # while running in a different cwd makes it read/write the wrong case and
    # trip the sha1 IOstream check.
    env['PWD'] = str(Path(cwd).resolve())
    with open(os.devnull, 'wb') as devnull:
        if log:
            fh = open(log, 'w')
        else:
            fh = devnull
        try:
            p = subprocess.run(cmd, cwd=str(cwd), env=env, stdout=fh,
                               stderr=subprocess.STDOUT, timeout=timeout)
        finally:
            if log:
                fh.close()
    return p.returncode


# ---------------------------------------------------------------- mesh
def build_wing_msh(out_msh, delta_deg=0.0, chord=CHORD, span=SPAN,
                   far=8.0, n_span=24, wall_size=0.012, far_size=0.35,
                   bl_layers=20, bl_h0=8e-6, bl_growth=1.2, bl_thick=0.02,
                   wing_name='wing', sym_name='symmetry', tip_name='tip',
                   far_name='farfield', split_tip=False, verbose=False):
    """Extruded O-grid half wing, written as gmsh MSH 2.2.

    split_tip separates the two ends of the extrusion by whether they touch the
    airfoil contour: the wing tip becomes a `wall` while the rest of that plane
    stays `farfield`. The 2-D stage wants the whole cap as one `empty` patch,
    so this is off by default.
    """
    import sys
    # load su2/wing3d.py by path: a plain `import wing3d` is ambiguous because
    # foam/ was on sys.path first and used to hold a same-named module
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        'su2_wing3d', str(R / 'su2' / 'wing3d.py'))
    W = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(W)
    sys.path.insert(0, str(R / 'su2'))
    import gmsh

    gmsh.initialize()
    gmsh.option.setNumber('General.Terminal', 1 if verbose else 0)
    gmsh.model.add('wing')
    g = gmsh.model
    geo = g.occ

    curves = W._contour_curves(geo, str(R / 'data/airfoils/sd7062.dat'),
                               chord, delta_deg)
    loop_af = geo.addCurveLoop(curves)
    nf = 72
    cx = 0.5 * chord
    pf = [geo.addPoint(cx + far * cos(2 * pi * i / nf), 0.0,
                       far * sin(2 * pi * i / nf))
          for i in range(nf)]
    arcs = [geo.addLine(pf[i], pf[(i + 1) % nf]) for i in range(nf)]
    loop_far = geo.addCurveLoop(arcs[::-1])
    surf = geo.addPlaneSurface([loop_far, loop_af])
    geo.synchronize()

    gmsh.option.setNumber('Mesh.MeshSizeExtendFromBoundary', 0)
    gmsh.option.setNumber('Mesh.MeshSizeFromPoints', 0)
    gmsh.option.setNumber('Mesh.MeshSizeFromCurvature', 0)
    gmsh.option.setNumber('Mesh.Algorithm', 5)

    fd = g.mesh.field.add('Distance')
    g.mesh.field.setNumbers(fd, 'CurvesList', curves)
    g.mesh.field.setNumber(fd, 'Sampling', 80)
    ft = g.mesh.field.add('Threshold')
    g.mesh.field.setNumber(ft, 'InField', fd)
    g.mesh.field.setNumber(ft, 'SizeMin', wall_size)
    g.mesh.field.setNumber(ft, 'SizeMax', far_size)
    g.mesh.field.setNumber(ft, 'DistMin', 0.01)
    g.mesh.field.setNumber(ft, 'DistMax', 2.0)
    g.mesh.field.setAsBackgroundMesh(ft)

    if bl_layers:
        blf = g.mesh.field.add('BoundaryLayer')
        g.mesh.field.setNumbers(blf, 'CurvesList', curves)
        g.mesh.field.setNumber(blf, 'Size', bl_h0)
        g.mesh.field.setNumber(blf, 'Ratio', bl_growth)
        g.mesh.field.setNumber(blf, 'Quads', 1)
        g.mesh.field.setNumber(blf, 'Thickness', bl_thick)
        g.mesh.field.setAsBoundaryLayer(blf)

    g.mesh.generate(2)
    out = geo.extrude([(2, surf)], 0.0, span, 0.0, [n_span],
                      recombine=bool(bl_layers))
    geo.synchronize()
    g.mesh.generate(3)

    # classify: caps by y, walls by boundary curves, rest is farfield
    cont = set(curves)
    groups = {}
    for dd, t in g.getEntities(2):
        bb = g.getBoundingBox(2, t)
        if abs(bb[1]) < 1e-6 and abs(bb[4]) < 1e-6:
            if split_tip:
                # root plane: only the part away from the airfoil is external
                bnd = {c for _, c in g.getBoundary([(2, t)])}
                if not (bnd & cont):
                    groups.setdefault(far_name, []).append(t)
                    continue
            groups.setdefault(sym_name, []).append(t)
        elif abs(bb[1] - span) < 1e-6 and abs(bb[4] - span) < 1e-6:
            if split_tip:
                bnd = {c for _, c in g.getBoundary([(2, t)])}
                # wing tip is a solid surface; the rest of the plane is the
                # outer far-field boundary
                groups.setdefault(tip_name if (bnd & cont) else far_name,
                                  []).append(t)
            else:
                groups.setdefault(tip_name, []).append(t)
        else:
            bnd = {c for _, c in g.getBoundary([(2, t)])}
            groups.setdefault(wing_name if (bnd & cont) else far_name,
                              []).append(t)
    for name, ts in groups.items():
        if ts:
            g.addPhysicalGroup(2, ts, name=name)
    g.addPhysicalGroup(3, [t for _, t in g.getEntities(3)], name='internal')
    geo.synchronize()

    gmsh.option.setNumber('Mesh.MshFileVersion', 2.2)
    gmsh.write(str(out_msh))
    counts = {k: len(v) for k, v in groups.items()}
    nn = len(g.mesh.getNodes()[0])
    ne = sum(len(c) for c in g.mesh.getElements(3)[1])
    gmsh.finalize()
    return dict(nodes=nn, elems=ne, groups=counts)


# ---------------------------------------------------------------- case
def write_case(case, mesh_msh, patch_map, converged_iters=4000,
               turbulence='kOmegaSST', write_interval=500):
    """Write a steady incompressible RANS case around an existing mesh."""
    import sys
    sys.path.insert(0, str(R / 'su2'))
    from math import cos, pi, sin

    case = Path(case)
    for sub in ('0', 'constant', 'system'):
        (case / sub).mkdir(parents=True, exist_ok=True)

    (case / 'system' / 'controlDict').write_text(f"""\
FoamFile {{ version 2.0; format ascii; class dictionary;
            object controlDict; }}
application     simpleFoam;
startFrom       startTime;
startTime       0;
stopAt          endTime;
endTime         {converged_iters};
deltaT          1;
writeControl    timeStep;
writeInterval   {write_interval};
purgeWrite      0;
writeFormat     ascii;
writePrecision  8;
writeCompression off;
timeFormat      general;
timePrecision   6;
runTimeModifiable true;
// NOTE: no functions block. This Ubuntu v1912 package has an incomplete
// runtime layout (FOAM_ETC/FOAM_LIBBIN empty, /usr/share/openfoam/bin absent),
// so any function object falls back to dynamicCode and dies writing a sha1
// stream. Forces are integrated from foamToVTK output in foam/forces.py.
""")

    (case / 'system' / 'fvSchemes').write_text("""\
FoamFile { version 2.0; format ascii; class dictionary; object fvSchemes; }
ddtSchemes      { default steadyState; }
gradSchemes     { default cellLimited Gauss linear 1; }
divSchemes
{
    default         none;
    div(phi,U)      bounded Gauss linearUpwind grad(U);
    div(phi,k)      bounded Gauss upwind;
    div(phi,omega)  bounded Gauss upwind;
    div((nuEff*dev2(T(grad(U))))) Gauss linear;
}
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes   { default corrected; }
wallDist        { method meshWave; }
""")

    (case / 'system' / 'fvSolution').write_text("""\
FoamFile { version 2.0; format ascii; class dictionary; object fvSolution; }
solvers
{
    p
    {
        solver          GAMG;
        smoother        GaussSeidel;
        tolerance       1e-7;
        relTol          0.05;
    }
    "(U|k|omega)"
    {
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-8;
        relTol          0.05;
    }
}
SIMPLE
{
    nNonOrthogonalCorrectors 2;
    consistent      yes;
    residualControl
    {
        // lift is fixed by the pressure field, so letting p converge tighter
        // than default is what actually terminates the run
        p               1e-3;
        U               1e-5;
        "(k|omega)"     1e-5;
    }
}
relaxationFactors
{
    equations
    {
        U               0.9;
        "(k|omega)"     0.9;
    }
}
""")

    # ---- 0/ fields -------------------------------------------------
    k_inf = 1.5 * (U_INF * 0.001) ** 2       # TI ~ 0.1%
    omega_inf = k_inf ** 0.5 / (0.09 ** 0.25 * 0.1 * CHORD)
    nu_t_inf = k_inf / omega_inf
    far_fields = dict(k=k_inf, omega=omega_inf, nut=nu_t_inf)

    for name, txt in _field_files(far_fields, patch_map).items():
        (case / '0' / name).write_text(txt)

    (case / 'constant' / 'transportProperties').write_text(f"""\
FoamFile {{ version 2.0; format ascii; class dictionary;
            object transportProperties; }}
transportModel  Newtonian;
nu              [0 2 -1 0 0 0 0] {NU};
""")
    (case / 'constant' / 'turbulenceProperties').write_text(f"""\
FoamFile {{ version 2.0; format ascii; class dictionary;
            object turbulenceProperties; }}
simulationType  RAS;
RAS
{{
    RASModel        {turbulence};
    turbulence      on;
    printCoeffs     on;
}}
""")
    # copy the mesh in and convert
    shutil.copy(mesh_msh, case / 'mesh.msh')
    return case


def mesh_to_foam(case, mesh_msh, patch_types):
    """Copy the MSH into the case and convert with gmshToFoam.

    patch_types maps patch name -> OpenFOAM type for the *inlet/outlet* hint;
    gmshToFoam derives patch names from the MSH physical names.
    """
    case = Path(case)
    shutil.copy(mesh_msh, case / 'mesh.msh')
    rc = run_cmd(['gmshToFoam', 'mesh.msh'], case, log=case / 'gmshToFoam.log')
    return rc


def _field_files(far, pm):
    """Build 0/U, 0/p, 0/k, 0/omega, 0/nut with one block per patch."""
    walls = list(pm['walls'])
    sym = pm.get('symmetry')
    empty = list(pm.get('empty', []))
    farfield = pm['farfield']

    def wall_block(kind, value=None):
        v = f' value uniform {value};' if value is not None else ''
        return ''.join(f'    {w} {{ type {kind};{v} }}\n' for w in walls)

    def sym_block(kind='symmetryPlane'):
        out = f'    {sym} {{ type {kind}; }}\n' if sym else ''
        out += ''.join(f'    {e} {{ type empty; }}\n' for e in empty)
        return out

    U = f"""FoamFile {{ version 2.0; format ascii; class volVectorField;
            object U; }}
dimensions      [0 1 -1 0 0 0 0];
internalField   uniform ({U_INF} 0 0);
boundaryField
{{
{wall_block('noSlip')}{sym_block()}\
    {farfield} {{ type freestream; freestreamValue uniform ({U_INF} 0 0); }}
}}
"""
    p = f"""FoamFile {{ version 2.0; format ascii; class volScalarField;
            object p; }}
dimensions      [0 2 -2 0 0 0 0];
internalField   uniform 0;
boundaryField
{{
{wall_block('zeroGradient')}{sym_block()}\
    {farfield} {{ type freestreamPressure; freestreamValue uniform 0; }}
}}
"""
    k = f"""FoamFile {{ version 2.0; format ascii; class volScalarField;
            object k; }}
dimensions      [0 2 -2 0 0 0 0];
internalField   uniform {far['k']};
boundaryField
{{
{wall_block('kqRWallFunction', far['k'])}{sym_block()}\
    {farfield} {{ type inletOutlet; inletValue uniform {far['k']};
                  value uniform {far['k']}; }}
}}
"""
    omega = f"""FoamFile {{ version 2.0; format ascii; class volScalarField;
            object omega; }}
dimensions      [0 0 -1 0 0 0 0];
internalField   uniform {far['omega']};
boundaryField
{{
{wall_block('omegaWallFunction', far['omega'])}{sym_block()}\
    {farfield} {{ type inletOutlet; inletValue uniform {far['omega']};
                  value uniform {far['omega']}; }}
}}
"""
    nut = f"""FoamFile {{ version 2.0; format ascii; class volScalarField;
            object nut; }}
dimensions      [0 2 -1 0 0 0 0];
internalField   uniform {far['nut']};
boundaryField
{{
{wall_block('nutkWallFunction', 0)}{sym_block()}\
    {farfield} {{ type calculated; value uniform {far['nut']}; }}
}}
"""
    return {'U': U, 'p': p, 'k': k, 'omega': omega, 'nut': nut}
