"""OpenFOAM clean-wing gate, stage 1: 2-D (infinite wing) validation.

The extruded O-grid is solved with `empty` front/back patches, which makes it
a genuine 2-D case in OpenFOAM. That removes all finite-wing geometry
ambiguity and tests the parts that matter together - gmshToFoam, the
kOmegaSST setup, wall resolution and the forces function object - against the
Layer-2 XFOIL polars, which are an independent reference already in the repo.

Re is set to 1e6 to match the layer-2 XFOIL Re = 1e6 series exactly, so the
comparison is like-for-like. SST has no transition model, so the correct
reference is the forced-transition (fully turbulent) polar.

Stage 2 (finite wing, snappyHexMesh, CLalpha vs the 0.0835-0.0881 oracle)
follows once this passes.
"""
import csv
import math
import sys
from pathlib import Path

R = Path('/home/yanlu/Documents/plane')
sys.path.insert(0, str(R / 'foam'))

import ofcase                                            # noqa: E402

RE_TARGET = 1.0e6
CHORD = 0.6788
NU = 1.46e-5
U = RE_TARGET * NU / CHORD          # 21.5 m/s
ALPHAS = (-4.0, -2.0, 0.0, 2.0, 4.0, 6.0)

# OF free-stream vector for each alpha is set through the farfield patch, so a
# single mesh can serve every angle by rotating the inlet vector.
BASE = ofcase


def set_patch_types(bnd_path, types):
    """Rewrite the `type` entry of named patches in constant/polyMesh/boundary.

    gmshToFoam maps each MSH physical group to a `patch` regardless of what it
    represents, so wall functions and 2-D `empty` patches must be retyped.
    """
    import re
    txt = bnd_path.read_text()
    for name, ptype in types.items():
        txt, n = re.subn(
            r'(\n    ' + re.escape(name) + r'\n    \{\n        type\s+)\w+;',
            r'\g<1>' + ptype + ';', txt)
        if n == 0:
            raise RuntimeError(f'patch {name} not found in boundary')
    bnd_path.write_text(txt)


def angle_case(case, alpha, mesh_msh):
    """Write a 2-D case at angle alpha by rotating the free-stream vector."""
    import math
    ux = U * math.cos(math.radians(alpha))
    uz = U * math.sin(math.radians(alpha))
    ofcase.U_INF = U
    ofcase.write_case(case, mesh_msh,
                      dict(walls=['wing'], symmetry=None,
                           farfield='farfield', empty=['sym', 'tip']),
                      converged_iters=3000, write_interval=3000)
    # gmshToFoam needs the case dictionaries present (it reads controlDict) and
    # it is what writes constant/polyMesh/boundary, so the conversion must sit
    # between write_case and the patch typing below
    ofcase.mesh_to_foam(case, mesh_msh, {})
    # gmshToFoam types every physical group as `patch`. Wall functions require
    # the wing to be a `wall`, and a 2-D case wants `empty` on both span ends.
    set_patch_types(case / 'constant' / 'polyMesh' / 'boundary',
                    {'wing': 'wall', 'symmetry': 'empty', 'tip': 'empty'})

    for name in ('U', 'p', 'k', 'omega', 'nut'):
        f = case / '0' / name
        t = f.read_text()
        if name == 'U':
            t = t.replace(f'uniform ({U} 0 0)', f'uniform ({ux:.6f} 0 {uz:.6f})')
            t = t.replace(f'freestreamValue uniform ({U} 0 0);',
                          f'freestreamValue uniform ({ux:.6f} 0 {uz:.6f});')
        t = t.replace('    symmetry { type symmetryPlane; }',
                      '    symmetry { type empty; }')
        f.write_text(t)
    return case


# The extruded mesh is one CELL thick but SPAN metres long (the standard
# OpenFOAM way to run 2-D), so the wing patch area is perimeter x SPAN and the
# reference area must be chord x SPAN. Using chord x 1 m inflates every
# coefficient by the span factor.
S_REF_2D = CHORD * ofcase.SPAN


def latest_forces(case, alpha):
    """CL, CD from our own integration of the foamToVTK wing patch."""
    import forces as F
    vtkroot = case / 'VTK'
    # foamToVTK also writes <case>_<time>.vtm beside the directory of the same
    # stem, and sorted() puts the file last, so filter to directories
    dirs = [p for p in vtkroot.glob(f'{case.name}_*') if p.is_dir()]
    if not dirs:
        return None
    vtp = sorted(dirs)[-1] / 'boundary' / 'wing.vtp'
    if not vtp.exists():
        return None
    f, _ = F.patch_forces(str(vtp), ofcase.RHO)
    a = math.radians(alpha)
    drag = f[0] * math.cos(a) + f[2] * math.sin(a)
    lift = -f[0] * math.sin(a) + f[2] * math.cos(a)
    q = 0.5 * ofcase.RHO * U * U * S_REF_2D
    return lift / q, drag / q


def main():
    mesh = R / 'foam/mesh/wing2d.msh'
    mesh.parent.mkdir(parents=True, exist_ok=True)
    info = ofcase.build_wing_msh(mesh, delta_deg=0.0, far=8.0, n_span=1,
                                 wall_size=0.012, far_size=0.3,
                                 bl_layers=22, bl_h0=8e-6, bl_growth=1.15,
                                 bl_thick=0.015)
    print('mesh:', info, flush=True)

    xfoil = {}
    with open(R / 'results/airfoil_polars.csv') as fh:
        for row in csv.DictReader(fh):
            if row['airfoil'] == 'sd7062' and row['mode'] == 'forced5' \
                    and float(row['Re']) == 1.0e6:
                xfoil[float(row['alpha'])] = (float(row['CL']),
                                              float(row['CD']))
    print(f'Re = {RE_TARGET:.1e}  U = {U:.2f} m/s  (XFOIL forced5 reference)')

    out = []
    for a in ALPHAS:
        case = R / f'foam/cases/w2d_a{int(a):+d}'
        if case.exists():
            import shutil
            shutil.rmtree(case)
        case.mkdir(parents=True)
        angle_case(case, a, mesh)
        ofcase.U_INF = U
        rc = ofcase.run_cmd(['simpleFoam'], case, log=case / 'solve.log',
                            timeout=7200)
        ofcase.run_cmd(['foamToVTK', '-latestTime'], case,
                       log=case / 'vtk.log', timeout=1800)
        res = latest_forces(case, a)
        if res is None:
            print(f'alpha={a:+.0f}: rc={rc} no forces', flush=True)
            continue
        cl, cd = res
        clx, cdx = xfoil.get(a, (float('nan'), float('nan')))
        out.append(dict(alpha=a, CL_OF=round(cl, 4), CD_OF=round(cd, 5),
                        CL_XFOIL=clx, CD_XFOIL=cdx,
                        dCL_pct=round(100 * (cl / clx - 1), 1) if clx else None))
        print(f'alpha={a:+5.1f}  CL_OF={cl:+.4f} (XFOIL {clx:+.4f}, '
              f'{100*(cl/clx-1):+.1f}%)  CD_OF={cd:.5f} '
              f'(XFOIL {cdx:.5f})', flush=True)

    if out:
        import numpy as np
        aa = np.array([r['alpha'] for r in out])
        cc = np.array([r['CL_OF'] for r in out])
        slope = float(np.polyfit(aa, cc, 1)[0])
        print(f'\n[stage1 2-D] C_Lalpha = {slope:.4f} /deg  '
              f'(2-D XFOIL oracle ~0.109-0.110)')
        with open(R / 'results/of_wing2d.csv', 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=list(out[0]))
            w.writeheader()
            w.writerows(out)
        print('saved: results/of_wing2d.csv')


if __name__ == '__main__':
    main()
