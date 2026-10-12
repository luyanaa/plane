"""Layer 12: grid-convergence check for the 3-D half-wing model.

Runs the flapped TO2 case at three mesh resolutions and reports CL/CD, so the
k_pw result can be quoted with a measured discretisation error.
"""
import csv
import sys
from pathlib import Path

R = Path('/home/yanlu/Documents/plane')
sys.path.insert(0, str(R / 'su2'))

from wing3d import build                        # noqa: E402
from su2case import write_cfg_3d, run_case      # noqa: E402

LEVELS = {
    'coarse': dict(n_span=16, wall_size=0.020, far_size=1.0),
    'medium': dict(n_span=24, wall_size=0.010, far_size=0.8),
    'fine': dict(n_span=34, wall_size=0.005, far_size=0.6),
}


def latest(hist):
    rows = list(csv.reader(open(hist)))
    hdr = [c.strip().replace('"', '') for c in rows[0]]
    iCL, iCD, iN = hdr.index('CL'), hdr.index('CD'), hdr.index('rms[Rho]')
    data = [x for x in rows[1:] if len(x) > max(iCL, iCD, iN)
            and x[0].strip() == '0']
    return float(data[-1][iCL]), float(data[-1][iCD]), float(data[-1][iN])


def main():
    out = R / 'results'
    out.mkdir(exist_ok=True)
    rows = []
    for name, kw in LEVELS.items():
        mesh = str(R / f'su2/cases/GC_{name}.su2')
        info = build(mesh, str(R / 'data/airfoils/sd7062.dat'),
                     delta_deg=20.0, far=8.0, **kw)
        cfg = write_cfg_3d(str(R / f'su2/cases/GC_{name}.cfg'), mesh,
                           ('wing',), sym='symmetry', far='farfield',
                           plot=('wing',), actdisk=None, mach=0.041, aoa=6.0,
                           re=14.07 * 0.6788 / 1.46e-5, ref_length=0.6788,
                           ref_area=0.6788 * 0.5 * 6.317, x_ref=0.25,
                           iters=6000, cfl=5.0, extra_walls=('tip',))
        run_case(cfg, str(R / f'su2/cases/GC_{name}run'), np_threads=24,
                 timeout=7200)
        cl, cd, rms = latest(R / f'su2/cases/GC_{name}run' / 'history.csv')
        rows.append(dict(level=name, nodes=info['npoin'], elems=info['nelem'],
                         CL=round(cl, 4), CD=round(cd, 5), rms=round(rms, 3)))
        print(f'{name:7s} {info["npoin"]:6d} nodes  CL={cl:+.4f} '
              f'CD={cd:.5f} rms={rms:.2f}', flush=True)
    with open(out / 'su2_gridconv.csv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f'saved: {out}/su2_gridconv.csv')


if __name__ == '__main__':
    main()
