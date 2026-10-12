"""Print the normalised spanwise loading from a completed 3-D case.

    python foam/span_report.py <case-dir>
"""
import sys
from pathlib import Path

import numpy as np

R = Path('/home/yanlu/Documents/plane')
sys.path.insert(0, str(R / 'foam'))
import wing_report as wing3d                                       # noqa: E402


def main():
    case = Path(sys.argv[1])
    if not case.is_absolute():
        case = R / case
    dirs = [p for p in (case / 'VTK').glob(f'{case.name}_*') if p.is_dir()]
    vtp = sorted(dirs)[-1] / 'boundary' / 'wing.vtp'
    print(case.name, '->', vtp)
    load = wing3d.span_loading(str(vtp), 0.0)
    ys = [o['y'] for o in load]
    ls = [o['Lprime'] for o in load]
    for o in load:
        print('y=%.3f  cl=%+.3f  Ln=%.4f  n=%d'
              % (o['y'], o['cl'], o['Ln'], o['n']))
    print('total L = %.2f N/m (integral)' % float(np.trapezoid(ls, ys)))


if __name__ == '__main__':
    main()
