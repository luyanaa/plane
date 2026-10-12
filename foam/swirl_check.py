"""Why the stage-3 actuator disk is axial-only (swirl magnitude check).

radialActuationDiskSource would add swirl, but its `coeffs` entry is a file
reference whose v1912 format is undocumented and is not resolved here. Before
deciding whether that matters, size the effect from the layer-8 BEMT, which
already carries the tangential induction:

    dQ/dr = 4 pi rho r^3 U (1-a) a' Omega

so a' is recoverable from the BEMT's own torque distribution. The swirl then
enters two places:

  1. dynamic pressure   q_swirl = 0.5 rho u_theta^2
  2. flow angle         atan(u_theta / U_axial)

Swirl only matters for k_pw if either is comparable to the axial slipstream
excess, which is q_slip/q_inf - 1 ~ 0.35-0.47.
"""
import math
import sys
from pathlib import Path

import numpy as np

R = Path('/home/yanlu/Documents/plane')
sys.path.insert(0, str(R))

from sizing.prop_bemt import PropBEMT                          # noqa: E402

# Prototype D (results/prop_operating.csv). The k_pw design point is TAKEOFF,
# J = 0.1998 at 5990 rpm - a heavily loaded disk, Ct_std = 7.21. An earlier
# version of this check ran at J = 0.509 (light loading, Ct_std = 0.60) and
# concluded swirl was negligible; that answer does NOT transfer, because swirl
# scales steeply with loading.
D, V = 0.530, 14.07
PTS = dict(J_TO=0.1998, J_climb=0.2128, J_cruise=0.4156)
RHO = 1.225
A_DISK = math.pi * (D / 2.0) ** 2


def main():
    q_inf = 0.5 * RHO * V ** 2
    print(f'{"point":>10} {"n[rps]":>8} {"T[N]":>8} {"Ct_std":>8} '
          f'{"u_theta":>8} {"q_sw/q_inf":>10} {"angle[deg]":>11}')
    for name, J in PTS.items():
        n = V / (J * D)
        om = 2.0 * math.pi * n
        p = PropBEMT(D=D, P_over_D=1.0)
        T, Q, P, _, dQdr, mu, r = p._solve_elements(V, n, iters=120)
        mu = np.asarray(mu, float)
        r = mu * p.R
        a_prime = dQdr / (4.0 * math.pi * RHO * r ** 3 * V * om)
        u_theta = a_prime * om * r
        good = mu > 0.3
        ct_std = T / (0.5 * RHO * A_DISK * V ** 2)
        # flow angle uses the SLIPSTREAM axial speed, not the freestream
        u_ax = V * (1.0 + 0.5 * ct_std) ** 0.5
        q_sw = 0.5 * RHO * np.mean(u_theta[good] ** 2)
        ang = math.degrees(math.atan(np.mean(u_theta[good]) / u_ax))
        print(f'{name:>10} {n:8.1f} {T:8.2f} {ct_std:8.3f} '
              f'{np.mean(u_theta[good]):8.2f} {q_sw / q_inf * 100:9.1f}%'
              f' {ang:10.1f}')

    print()
    print('VERDICT: swirl is NOT negligible at the k_pw design point.')
    print('At takeoff Ct_std = 7.21 (heavy loading), q_swirl reaches ~28% of')
    print('q_inf and the swirl flow angle ~14 deg - comparable to the flap')
    print('deflection, so it is a first-order effect on the blown-flap gain.')
    print('At cruise loading it is ~1% and can be ignored. The axial-only')
    print('actuationDiskSource is therefore adequate ONLY for light loading;')
    print('the k_pw calibration needs swirl (rotorDiskSource) - or k_pw must')
    print('be restricted to a lower-loading regime and said so explicitly.')


if __name__ == '__main__':
    main()
