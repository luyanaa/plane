"""Layer 12 / OpenFOAM stage 2 (body-fitted route): 3-D finite wing.

The snappyHexMesh route cannot mesh this wing: the sd7062 trailing edge is
0.47 mm thick and the layer algorithm truncates the prism stack to zero there
and then discards the whole stack, on every level and for every layer setting
tried (see README). gmsh's BoundaryLayer field handles it because it is
applied to the 2-D O-grid *before* extrusion, so the prisms are swept along a
well-shaped surface mesh - the same machinery that already passed stage 1.

Mesh: extruded O-grid, symmetry plane at the root, airfoil cap of the top face
treated as a wall (the wing tip). Coarse surface (~10 mm) with 20 boundary
layers, 8 um first cell (y+ < 1).

Oracle (3-D, finite-wing slope):
    C_L0      = 0.330 .. 0.349
    C_Lalpha  = 0.0835 .. 0.0881 /deg
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
import wing_report as wing3d                                       # noqa: E402
from wing2d import set_patch_types                         # noqa: E402

ALPHAS = [-4.0, -2.0, 0.0, 2.0, 4.0, 6.0]
MESH = R / 'foam/mesh/wing3d.msh'


# Grid levels for the convergence study. The WALL TREATMENT is deliberately
# held fixed (first cell, growth, layer count) so a change in the answer comes
# from streamwise/spanwise discretisation and not from a different near-wall
# model; only the surface and span divisions vary.
LEVEL_MESH = dict(
    c=dict(n_span=20, wall_size=0.018, far=5.0, bl_layers=20, bl_h0=8e-6,
           bl_thick=0.02),
    g=dict(n_span=36, wall_size=0.010, far=6.0, bl_layers=20, bl_h0=8e-6,
           bl_thick=0.02),
    f=dict(n_span=52, wall_size=0.007, far=7.0, bl_layers=20, bl_h0=8e-6,
           bl_thick=0.02),
)
_last = {}


def build_mesh(tag='g', verbose=False):
    p = dict(LEVEL_MESH[tag[0]])
    p['far_size'] = 0.30
    p['verbose'] = verbose
    _last.clear()
    _last.update(p)
    return ofcase.build_wing_msh(
        MESH, far=p['far'], n_span=p['n_span'], wall_size=p['wall_size'],
        far_size=p['far_size'], bl_layers=p['bl_layers'], bl_h0=p['bl_h0'],
        bl_growth=1.2, bl_thick=p['bl_thick'], split_tip=True,
        verbose=verbose)


def make_case(level_tag, alpha, rebuild=False):
    case = R / f'foam/cases/w3b_{level_tag}_a{alpha:+.0f}'
    if rebuild and case.exists():
        shutil.rmtree(case)
    if (case / 'constant/polyMesh/boundary').exists():
        return case
    if case.exists():
        shutil.rmtree(case)
    case.mkdir(parents=True)
    ofcase.write_case(case, MESH,
                      dict(walls=['wing'], symmetry='symmetry',
                           farfield='farfield'),
                      converged_iters=1200, write_interval=1200)
    ofcase.mesh_to_foam(case, MESH, {})
    set_patch_types(case / 'constant/polyMesh/boundary',
                    {'wing': 'wall', 'symmetry': 'symmetryPlane'})
    # rotate the freestream (geometry stays put) so one mesh serves every angle
    ux = ofcase.U_INF * math.cos(math.radians(alpha))
    uz = ofcase.U_INF * math.sin(math.radians(alpha))
    for name in ('U', 'k', 'omega', 'nut', 'p'):
        f = case / '0' / name
        t = f.read_text()
        t = t.replace(f'uniform ({ofcase.U_INF} 0 0)',
                      f'uniform ({ux:.6f} 0 {uz:.6f})')
        f.write_text(t)
    return case


def main():
    tag = sys.argv[1] if len(sys.argv) > 1 else 'g20'
    print('mesh:', build_mesh(tag), flush=True)
    rows = []
    for a in ALPHAS:
        case = make_case(tag, a)
        ofcase.U_INF = 14.07
        rc = ofcase.run_cmd(['simpleFoam'], case, log=case / 'log.simpleFoam',
                            timeout=14400)
        ofcase.run_cmd(['foamToVTK', '-latestTime'], case, log=case / 'log.vtk')
        r = wing3d.alpha_result(case, a)
        if r is None:
            print(f'[{tag}] alpha={a:+.0f} no VTK', flush=True)
            continue
        rows.append(r)
        print(f'[{tag}] alpha={a:+.0f}  CL={r["CL"]:+.4f}  CD={r["CD"]:.5f}',
              flush=True)
        # keep the span loading at the two angles that matter for k_pw
        if abs(a) < 1e-6:
            load = wing3d.span_loading(r['vtp'], a)
            with open(R / f'results/of_w3b_{tag}_span.csv', 'w',
                      newline='') as fh:
                w = csv.DictWriter(fh, fieldnames=['y', 'Lprime', 'cl', 'Ln',
                                                   'n'])
                w.writeheader()
                w.writerows(load)
            for o in load[::3]:
                print(f'    y={o["y"]:.2f} cl={o["cl"]:+.3f} Ln={o["Ln"]:.4f}',
                      flush=True)

    if len(rows) >= 3:
        a = np.array([r['alpha'] for r in rows], float)
        cl = np.array([r['CL'] for r in rows], float)
        lin = np.abs(a) <= 2.0                       # exclude stall at a>=4
        A = np.vstack([a[lin], np.ones(lin.sum())]).T
        slope, cl0 = np.linalg.lstsq(A, cl[lin], rcond=None)[0]
        print(f'[{tag}] linear-range fit |alpha|<=2:  CLalpha={slope:.5f} '
              f'/deg  CL0={cl0:+.4f}', flush=True)
    with open(R / f'results/of_w3b_{tag}.csv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=['alpha', 'CL', 'CD'])
        w.writeheader()
        w.writerows([dict(alpha=r['alpha'], CL=r['CL'], CD=r['CD'])
                     for r in rows])


if __name__ == '__main__':
    main()
