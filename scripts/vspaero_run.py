"""VSPAERO sweep runner on a .vsp3 model (layer-4/5 bridge).

Usage:
  .venv/bin/python scripts/vspaero_run.py MODEL.vsp3 OUT_PREFIX [alpha_start alpha_end alpha_npts]

Writes OUT_PREFIX.polar / OUT_PREFIX.history. Note: reference values
(Sref/Cref/Bref) come from the model's AnalysisRef; unset models use
VSPAERO defaults (Sref=100) - coefficients are scaled accordingly.
Propeller geoms are skipped here (static VLM diverges); layer 6 handles
rotor/actuator-disk cases.
"""
import sys
from pathlib import Path

import openvsp as vsp


def main():
    model = Path(sys.argv[1])
    prefix = sys.argv[2]
    a0, a1, an = 0.0, 8.0, 5
    if len(sys.argv) > 5:
        a0, a1, an = map(float, sys.argv[3:6])
        an = int(an)

    vsp.ReadVSPFile(str(model))
    vsp.SetAnalysisInputDefaults('VSPAEROComputeGeometry')
    vsp.ExecAnalysis('VSPAEROComputeGeometry')

    vsp.SetAnalysisInputDefaults('VSPAEROSweep')
    names = vsp.GetAnalysisInputNames('VSPAEROSweep')

    def setd(nm, val):
        if nm in names:
            vsp.SetDoubleAnalysisInput('VSPAEROSweep', nm, [val])

    def seti(nm, val):
        if nm in names:
            vsp.SetIntAnalysisInput('VSPAEROSweep', nm, [val])

    setd('AlphaStart', a0)
    setd('AlphaEnd', a1)
    seti('AlphaNpts', an)
    setd('MachStart', 0.05)
    setd('MachEnd', 0.05)
    seti('MachNpts', 1)
    setd('BetaStart', 0.0)
    setd('BetaEnd', 0.0)
    seti('BetaNpts', 1)
    seti('WakeNumIter', 10)
    vsp.ExecAnalysis('VSPAEROSweep')
    print(f'DONE: {prefix}.polar')


if __name__ == '__main__':
    main()
