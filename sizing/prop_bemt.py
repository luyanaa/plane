"""Layer 8: propeller BEMT - C_T(J), C_P(J), eta_p(J).

Blade-element momentum theory with Prandtl tip/hub loss, constant-pitch
geometry (twist from P/D), linear 2-D polar. Curves are n-independent to
first order (inviscid scaling), so they are generated once per prop on a
J-grid and interpolated by the mission model. Static thrust (J->0) uses a
small freestream with under-relaxed induction iteration.
"""
from dataclasses import dataclass
from math import acos, atan, atan2, cos, exp, pi, sin, sqrt, tan

import numpy as np

from .atmosphere import RHO0


@dataclass(frozen=True)
class PropBEMT:
    D: float               # diameter, m
    B: int = 2             # blade count
    P_over_D: float = 0.85 # geometric pitch ratio (design variable)
    n_rpm: float = 7500.0  # reference RPM for curve generation
    r_hub_frac: float = 0.15
    c_ref_frac: float = 0.08   # max chord / D
    alpha0_deg: float = 0.0    # zero-lift angle (symmetric polar placeholder)
    cl_alpha: float = 5.7      # 2-D lift slope, 1/rad
    cd0: float = 0.008
    cd_k: float = 0.015        # CD = cd0 + cd_k*CL^2

    @property
    def R(self):
        return self.D / 2.0

    def chord(self, mu):
        """Chord at radial station mu = r/R (tapered, max ~0.75R)."""
        return self.D * self.c_ref_frac * (0.55 + 0.45 * mu) ** 1.2 / 1.0

    def twist(self, mu, r_deg=True):
        """Geometric twist (deg). Constant-pitch blade: beta(r) = atan(P/(2 pi r)),
        which in terms of the pitch ratio and r = mu R = mu D/2 becomes
        beta = atan( (P/D) / (pi mu) )."""
        beta = atan(self.P_over_D / (pi * mu))
        return beta * (180.0 / pi) if r_deg else beta

    # ------------------------------------------------------------ solver
    def solve(self, V, n_rps, a0=0.1, a0p=0.0, iters=60):
        """BEMT solution at flight speed V and rotation rate n (rev/s).

        Returns (T, Q, P_shaft) with 30 blade elements.
        """
        T, Q, P, *_ = self._solve_elements(V, n_rps, iters, a0, a0p)
        return T, Q, P

    def _solve_elements(self, V, n_rps, iters=60, a0=0.1, a0p=0.0, mu=None):
        """Full BEMT state: (T, Q, P, dT/dr, dQ/dr, mu, r)."""
        if mu is None:
            mu = np.linspace(self.r_hub_frac + 0.01, 0.99, 30)
        else:
            mu = np.asarray(mu, dtype=float)
        r = mu * self.R
        c = np.array([self.chord(m) for m in mu])
        beta = np.array([self.twist(m) for m in mu])
        omega = 2.0 * pi * n_rps
        a = np.full_like(mu, a0)
        ap = np.full_like(mu, a0p)
        sig = self.B * c / (2.0 * pi * r)
        for _ in range(iters):
            phi = np.arctan2(V * (1.0 + a), omega * r * (1.0 - ap))
            alpha = np.radians(beta) - phi
            cl = self.cl_alpha * (alpha - np.radians(self.alpha0_deg))
            cd = self.cd0 + self.cd_k * cl ** 2
            # Prandtl tip + hub loss
            f_tip = self.B / 2.0 * (1.0 - mu) / np.maximum(np.sin(phi), 1e-6)
            f_hub = self.B / 2.0 * (mu - self.r_hub_frac) / np.maximum(np.sin(phi), 1e-6)
            f = (2.0 / pi) * np.arccos(np.exp(-np.maximum(f_tip, 0))) \
                * (2.0 / pi) * np.arccos(np.exp(-np.maximum(f_hub, 0)))
            cn = cl * np.cos(phi) - cd * np.sin(phi)
            ct = cl * np.sin(phi) + cd * np.cos(phi)
            # axial:  a = X/(1-X), X = sig*cn/(4 f sin^2 phi); heavy-load clamp
            x = sig * cn / (4.0 * f * np.sin(phi) ** 2 + 1e-9)
            a_new = np.where(x < 0.9, x / (1.0 - x), 4.0)
            # tangential: a' = sig*ct/(4 f sin cos + sig*ct)
            ap_new = sig * ct / (4.0 * f * np.sin(phi) * np.cos(phi) + sig * ct)
            a = 0.6 * a + 0.4 * a_new
            ap = 0.6 * ap + 0.4 * ap_new
        # integrate forces
        v_rel2 = (V * (1.0 + a)) ** 2 + (omega * r * (1.0 - ap)) ** 2
        v_rel = np.sqrt(v_rel2)
        phi = np.arctan2(V * (1.0 + a), omega * r * (1.0 - ap))
        alpha = np.radians(beta) - phi
        cl = self.cl_alpha * (alpha - np.radians(self.alpha0_deg))
        cd = self.cd0 + self.cd_k * cl ** 2
        dT = 0.5 * RHO0 * v_rel2 * c * (cl * np.cos(phi) - cd * np.sin(phi))
        dQ = 0.5 * RHO0 * v_rel2 * c * r * (cl * np.sin(phi) + cd * np.cos(phi))
        T = float(self.B * np.trapezoid(dT, r))
        Q = float(self.B * np.trapezoid(dQ, r))
        return T, Q, Q * omega, self.B * dT, self.B * dQ, mu, r

    # ------------------------------------------------------------ curves
    def curves(self, j_max=1.2, n_pts=25, n_rpm=None):
        """C_T(J), C_P(J), eta_p(J) at fixed reference RPM."""
        n_rpm = n_rpm or self.n_rpm
        n = n_rpm / 60.0
        js = np.linspace(0.0, j_max, n_pts)
        ct, cp, eta = [], [], []
        for j in js:
            v = j * n * self.D
            if j < 0.02:
                v = 0.5                      # static: tiny freestream
            t, q, p = self.solve(v, n)
            ct.append(t / (RHO0 * n ** 2 * self.D ** 4))
            cp.append(p / (RHO0 * n ** 3 * self.D ** 5))
        ct, cp = np.array(ct), np.array(cp)
        eta = np.where(cp > 0, js * ct / np.maximum(cp, 1e-9), 0.0)
        return js, ct, cp, eta

    def point(self, V, n_rps):
        """Non-dimensional coefficients at speed V and rate n."""
        t, q, p = self.solve(V, n_rps)
        ct = t / (RHO0 * n_rps ** 2 * self.D ** 4)
        cp = p / (RHO0 * n_rps ** 3 * self.D ** 5)
        j = V / (n_rps * self.D)
        return ct, cp, j * ct / cp if cp > 0 else 0.0

    def thrust_at_speed(self, V, n_rps):
        return self.solve(V, n_rps)[0]

    def radial_loads(self, V, n_rps, iters=60, mu=None):
        """Radial load distribution for an SU2 VARIABLE_LOAD actuator disk.

        Returns (r_bar, dCT_drbar, dCP_drbar, dCr_drbar, J) where the
        coefficients use the Renard convention SU2 expects:
        dC_T/drbar = R * (dT/dr) / (rho n^2 D^4)   [force ref rho n^2 D^4]
        dC_P/drbar = omega * R * (dQ/dr) / (rho n^3 D^5)  [power ref rho n^3 D^5]
        dC_R/drbar = 0 (radial force; swirl is carried by dC_P in this model).
        """
        T, Q, P, dT, dQ, mu, r = self._solve_elements(V, n_rps, iters, mu=mu)
        D, R, n = self.D, self.R, n_rps
        omega = 2.0 * pi * n
        dct = R * dT / (RHO0 * n ** 2 * D ** 4)
        dcp = omega * R * dQ / (RHO0 * n ** 3 * D ** 5)
        j = V / (n * D) if n > 0 else 0.0
        return mu, dct, dcp, np.zeros_like(mu), j


def make_prop(pack):
    """Prop geometry per engine pack (design choice: P/D by cruise J target)."""
    p_d = 0.9 if pack.D_prop >= 0.6 else 1.0   # singles: lower pitch, twins: high-J
    return PropBEMT(D=pack.D_prop, B=2, P_over_D=p_d, n_rpm=7500.0)
