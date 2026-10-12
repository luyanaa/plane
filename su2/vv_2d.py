"""Layer 12: 2-D airfoil mesh resolution study + V&V vs XFOIL (layer 2).

Sweeps near-field resolution (no gmsh BoundaryLayer field - it produces
degenerate cells here) and reports convergence of CL/CD for naca2412 at
Re = 1e6 against the layer-2 XFOIL polars.
"""
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from meshgen import mesh_airfoil_2d          # noqa: E402
from su2case import write_cfg_2d, run_case   # noqa: E402

ROOT = Path('/home/yanlu/Documents/plane')
BASE = ROOT / 'su2'
OUT = ROOT / 'results'

XFOIL = {0.0: (0.236, 0.00600), 4.0: (0.711, 0.00790)}

MESHES = {
    'm1': dict(bl_layers=0, far=50.0, n_surf=260, surf_size=0.0010,
               far_size=2.0, ramp_chords=8.0),
    'm2': dict(bl_layers=0, far=50.0, n_surf=320, surf_size=0.0004,
               far_size=1.2, ramp_chords=12.0),
}


def latest(hist):
    rows = list(csv.reader(open(hist)))
    hdr = [c.strip().replace('"', '') for c in rows[0]]
    iCL, iCD, iN = hdr.index('CL'), hdr.index('CD'), hdr.index('rms[Rho]')
    data = [x for x in rows[1:] if len(x) > max(iCL, iCD, iN) and x[0].strip() == '0']
    x = data[-1]
    return float(x[iCL]), float(x[iCD]), float(x[iN])


def main():
    OUT.mkdir(exist_ok=True)
    rows = []
    for tag, kw in MESHES.items():
        mesh = BASE / f'cases/{tag}.su2'
        nn, ne = mesh_airfoil_2d(str(ROOT / 'data/airfoils/naca2412.dat'),
                                 str(mesh), **kw)
        print(f'{tag}: {nn} nodes, {ne} elems')
        for aoa in (0.0, 4.0):
            cfg = write_cfg_2d(str(BASE / f'cases/{tag}_{int(aoa)}.cfg'), str(mesh),
                               mach=0.05, aoa=aoa, re=1.0e6, iters=8000,
                               cfl=2.0, trans='none')
            run_case(cfg, str(BASE / f'cases/{tag}_{int(aoa)}run'),
                     np_threads=16, timeout=3600)
            cl, cd, rms = latest(BASE / f'cases/{tag}_{int(aoa)}run' / 'history.csv')
            clx, cdx = XFOIL[aoa]
            rows.append(dict(mesh=tag, nodes=nn, elems=ne, aoa=aoa,
                             CL=round(cl, 4), CD=round(cd, 5),
                             CL_xfoil=clx, CD_xfoil=cdx,
                             dCL_pct=round(100 * (cl / clx - 1), 1),
                             dCD_pct=round(100 * (cd / cdx - 1), 1),
                             rms=round(rms, 3)))
            print(f'  {tag} aoa={aoa:4.1f}  CL={cl:.4f} ({100*(cl/clx-1):+.1f}%)'
                  f'  CD={cd:.5f} ({100*(cd/cdx-1):+.1f}%)  rms={rms:.2f}')
    with open(OUT / 'su2_2d_vv.csv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f'saved: {OUT}/su2_2d_vv.csv')


if __name__ == '__main__':
    main()
