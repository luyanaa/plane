"""Smoke test: OpenFOAM wing mesh build + gmshToFoam conversion.

The decisive question is whether OpenFOAM accepts the extruded O-grid with
boundary-layer prisms that SU2 could not take. OpenFOAM handles general
polyhedra, so it should.
"""
import sys
from pathlib import Path

R = Path('/home/yanlu/Documents/plane')
sys.path.insert(0, str(R / 'foam'))
sys.path.insert(0, str(R / 'su2'))

from ofcase import build_wing_msh, run_cmd, write_case    # noqa: E402


def main():
    msh = R / 'foam/mesh/smoke.msh'
    msh.parent.mkdir(parents=True, exist_ok=True)
    info = build_wing_msh(msh, delta_deg=0.0, far=6.0, n_span=12,
                          wall_size=0.02, far_size=0.4, bl_layers=15,
                          bl_h0=1e-5, bl_thick=0.012)
    print('mesh:', info)
    print('msh size: %.1f MB' % (msh.stat().st_size / 1e6))

    import shutil
    case = R / 'foam/cases/smoke'
    if case.exists():
        shutil.rmtree(case)
    case.mkdir(parents=True)
    # OpenFOAM needs a valid case (system/controlDict) before gmshToFoam
    write_case(case, msh, dict(walls=['wing', 'tip'], symmetry='symmetry',
                               farfield='farfield'), converged_iters=10)
    from ofcase import mesh_to_foam
    rc = mesh_to_foam(case, msh, {})
    print('gmshToFoam rc =', rc)
    log = (case / 'gmshToFoam.log').read_text()
    for ln in log.splitlines():
        if any(k in ln for k in ('Creating', 'patch', 'error', 'Error',
                                 'cells', 'points')):
            print('  ', ln.strip()[:110])
    print('\nboundary file:')
    for f in sorted((case / 'constant' / 'polyMesh').glob('*')):
        print('  ', f.name, f.stat().st_size)


if __name__ == '__main__':
    main()
