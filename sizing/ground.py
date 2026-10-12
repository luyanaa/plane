"""Layer 11: takeoff / landing ground-run model with flap and propwash effects.

Ground roll:  m dV/dt = T(V) - D - mu (W - L)
Air phase:    accelerate V_R -> V_2, climb to obstacle
Landing:      air phase (approach at 1.3 V_s, flare) + ground phase with braking

Flap effects (layer-3 lookup): C_Lmax and CD increments per flap state.
Propwash (layer-6): slipstream dynamic-pressure boost on the flap-covered span
  increases the effective C_L on that span; applied to the inner-flap fraction
  for the twins and to the tail for the singles.
"""
from dataclasses import dataclass, field
from math import cos, pi, radians, sin, sqrt

import numpy as np

from .atmosphere import G, RHO0, isa
from .configs import MISSION, PROTOTYPES
from .engine_map import EngineMap
from .geometry import build_geometry
from .mission import Aircraft
from .prop_bemt import PropBEMT
from .propulsion import EnginePack

# layer-3 flap lookup (XFOIL, Re 8e5): (delta, chord) -> (C_Lmax, dCD)
FLAP_DB = {
    'naca2412': {'clean': (1.47, 0.000),
                 'takeoff': (1.66, 0.002),
                 'landing': (1.87, 0.012)},
    'sd7062': {'clean': (1.69, 0.000),
               'takeoff': (1.92, 0.002),
               'landing': (2.15, 0.014)},
}
AIRFOIL_OF = {'A': 'naca2412', 'B': 'sd7062', 'C': 'naca2412', 'D': 'sd7062'}
# layer-6 flap-span q-boost: slipstream dynamic-pressure ratio on the covered
# span (layer 6 computed q_slip/q = 7.7-9.0 at takeoff; used here on the flap
# lift increment)
FLAP_QBOOST = {'A': 1.0, 'B': 1.0, 'C': 3.12, 'D': 3.84}
FLAP_SPAN_FRAC = {'A': 0.0, 'B': 0.0, 'C': 0.47, 'D': 0.48}
# slipstream velocity at the flap as a fraction of the far-wake value
# (layer 6 used 0.9 to the far wake; at the flap, closer to the disk, the
# jet has not fully contracted - 0.75 is the working value)
K_WAKE = 0.75

# Slipstream-to-flap effectiveness: fraction of the momentum-disk dynamic
# pressure ratio that actually shows up as flap lift increment. Accounts for
# slipstream contraction past the disk, swirl-induced local alpha change,
# spanwise mixing and separation at high deflection.
#
# LAYER 12 STATUS: this is NOT yet a calibrated value. The CFD ladder has
# cleared the reference-area/frame chain (0.000%) and the actuator-disk input
# integration (1e-16), and Euler has cleared the geometry/topology chain, but
# no converged wall-resolved 3-D viscous case exists yet, so k_pw remains a
# sensitivity parameter spanning 0.4-1.0 rather than a nominal.
K_PW = 0.5
K_PW_BAND = (0.4, 0.5, 0.7, 1.0)


@dataclass
class GroundResult:
    s_ground: float = 0.0
    s_air: float = 0.0
    s_total: float = 0.0
    V_r: float = 0.0
    V_2: float = 0.0
    V_s: float = 0.0
    V_approach: float = 0.0
    V_touchdown: float = 0.0
    t_ground: float = 0.0
    accel_g: float = 0.0
    notes: list = field(default_factory=list)


class GroundRun:
    def __init__(self, key, frame, pack, spec=MISSION):
        self.key, self.frame, self.pack, self.spec = key, frame, pack, spec
        self.geo = build_geometry(frame.S, frame.AR, frame.taper, frame.t_over_c,
                                  frame.l_fus, frame.d_fus, n_nacelles=frame.n_nacelles)
        self.ac = Aircraft(frame, pack, spec)
        self.af = AIRFOIL_OF[key]
        met = self.ac.metrics(RHO0)
        self._v_cr = met['V_cruise']
        self._cd0 = self.ac.polar.C_D0
        self._k = self.ac.polar.k
        self._t0_eng = pack.static_thrust(RHO0) / pack.n_eng
        self._t_cr_eng = pack.eta_p_cruise * pack.P_kw_per_eng * 1000.0 / self._v_cr
        self.boost_frac = FLAP_SPAN_FRAC[key]   # instance-overridable
        # slipstream-to-flap effectiveness. Layer 12 CFD has not yet closed a
        # credible 3-D viscous mesh, so this stays a sensitivity parameter and
        # is NOT a calibrated nominal. See kpw_sensitivity.py.
        self.k_pw = K_PW

    def thrust(self, V, throttle=1.0):
        """Full-throttle thrust (momentum-disk linear model, cached params)."""
        if V < self._v_cr:
            t_eng = self._t0_eng + (self._t_cr_eng - self._t0_eng) * (V / self._v_cr)
        else:
            t_eng = self._t_cr_eng * self._v_cr / V
        return t_eng * self.pack.n_eng * throttle

    def drag(self, V, cl, state):
        q = 0.5 * RHO0 * V * V
        cd = self._cd0 + self._k * cl * cl + self.cd_extra(state)
        return q * self.frame.S * cd
    def clmax(self, state):
        return FLAP_DB[self.af][state][0]

    def cd_extra(self, state):
        return FLAP_DB[self.af][state][1]

    def q_slip_ratio(self, V, T_total, throttle=1.0):
        """Slipstream-to-freestream dynamic pressure ratio at the flap.

        Momentum disk: T = 2 rho A (V + v_i) v_i  ->  v_i(V,T).
        Slipstream velocity near the flap: v_s = V + 2 k_wake v_i.
        Returns q_slip/q_infinity (=1 when unpowered).
        """
        if self.pack.n_eng == 0 or V < 1e-3:
            return 1.0
        T = max(T_total, 0.0)
        if T <= 1e-6:
            return 1.0
        a_disk = pi * (self.pack.D_prop / 2.0) ** 2
        v_i = (-V + sqrt(V * V + 2.0 * T / (RHO0 * a_disk))) / 2.0
        v_s = V + 2.0 * K_WAKE * v_i
        return (v_s / V) ** 2

    def clmax_effective(self, state, powered=True, V=None, throttle=1.0):
        """C_Lmax with propwash boost on the flap-covered span.

        When V and the operating throttle are supplied, the boost comes from
        the momentum-disk slipstream at that (V, T) - so a low-power approach
        gets a correspondingly smaller boost.
        """
        base = self.clmax(state)
        frac = self.boost_frac
        if not powered or frac == 0.0:
            return base
        k = self.k_pw
        if V is not None:
            T = self.thrust(V, throttle) * self.pack.n_eng
            q_ratio = 1.0 + k * (self.q_slip_ratio(V, T) - 1.0)
        else:
            q_ratio = 1.0 + k * (FLAP_QBOOST[self.key] - 1.0)
        clean = FLAP_DB[self.af]['clean'][0]
        inc = base - clean
        return clean + inc * (frac * q_ratio + (1.0 - frac))

    # ------------------------------------------------------------ takeoff
    def takeoff(self, state='takeoff', mu=0.05, obstacle=15.0, dt=0.01,
                throttle=1.0):
        f, s = self.frame, self.spec
        m = f.MTOW_kg
        W = m * G
        # C_Lmax with the slipstream actually produced at V_R (not a constant);
        # evaluated once at the clean-speed reference (no self-consistent
        # iteration - that would feed back positively at low speed)
        v_s_clean = sqrt(2.0 * W / (RHO0 * f.S * self.clmax(state)))
        v_r0 = 1.1 * v_s_clean
        clmax = self.clmax_effective(state, powered=True, V=v_r0, throttle=throttle)
        v_s = sqrt(2.0 * W / (RHO0 * f.S * clmax))
        v_r = 1.1 * v_s
        v_2 = 1.2 * v_s
        r = GroundResult(V_s=v_s, V_r=v_r, V_2=v_2)
        r.notes.append(f'clmax_eff {clmax:.3f} at V_R {v_r:.1f} m/s, thr {throttle:.2f}')
        # ground roll (integrate in time)
        V, t, x, it = 0.05, 0.0, 0.0, 0
        while V < v_r and t < 120.0 and it < 20000:
            it += 1
            cl = clmax * min((V / v_r) ** 2, 1.0)
            q = 0.5 * RHO0 * V * V
            L = q * f.S * cl
            D = self.drag(V, cl, state)
            T = self.thrust(V, throttle)
            a = (T - D - mu * (W - L)) / m
            if a <= 0:
                r.notes.append('stalled ground roll (T <= D+mu(W-L))')
                break
            V += a * dt
            x += V * dt
            t += dt
        r.s_ground, r.t_ground = x, t
        r.accel_g = (self.thrust(v_r, throttle) - self.drag(v_r, clmax, state)) / (m * G)
        # air phase: accelerate to V_2 climbing to obstacle
        h, Va, xa, it = 0.0, v_r, x, 0
        cl = clmax * 0.85
        while h < obstacle and Va < v_2 * 1.4 and it < 20000:
            it += 1
            D = self.drag(Va, cl, state)
            T = self.thrust(Va, throttle)
            gam = np.arcsin(np.clip((T - D) / W, 0.02, 0.35))
            a = (T - D) / m
            Va += a * dt
            h += Va * sin(gam) * dt
            xa += Va * cos(gam) * dt
            if Va < 1.0:
                r.notes.append('air phase decelerating')
                break
        r.s_air = xa - x
        r.s_total = xa
        return r

    # ------------------------------------------------------------ landing
    def landing(self, state='landing', mu_brake=0.35, obstacle=15.0, dt=0.01,
                glideslope_deg=5.0, flare_h=3.0, approach_throttle=0.25):
        """Approach (glideslope at v_app) -> flare -> ground braking roll.

        The approach is powered at `approach_throttle`, so the flap gets only
        the slipstream that throttle actually produces - not the takeoff
        value. After touchdown spoilers deploy (ground C_L = 0.3 x C_Lmax).
        """
        f = self.frame
        m = f.MTOW_kg
        W = m * G
        v_s_clean = sqrt(2.0 * W / (RHO0 * f.S * self.clmax(state)))
        v_app0 = 1.3 * v_s_clean
        clmax = self.clmax_effective(state, powered=True, V=v_app0,
                                     throttle=approach_throttle)
        v_s = sqrt(2.0 * W / (RHO0 * f.S * clmax))
        v_app = 1.3 * v_s
        v_td = 1.15 * v_s
        r = GroundResult(V_s=v_s, V_approach=v_app, V_touchdown=v_td)
        r.notes.append(f'clmax_eff {clmax:.3f} at V_app {v_app:.1f} m/s, '
                       f'thr {approach_throttle:.2f}')
        gam = radians(glideslope_deg)
        # circular-arc flare at n_flare: constant-speed pull from the
        # glideslope to level, wheel height = arc height drop
        n_flare = 1.15
        R_f = v_app * v_app / (G * (n_flare - 1.0))
        h_flare = R_f * (1.0 - cos(gam))
        s_flare = R_f * sin(gam)
        if h_flare > obstacle:                      # very flat glide path
            h_flare, s_flare = obstacle, obstacle / max(np.tan(gam), 1e-3)
        s_glide = (obstacle - h_flare) / max(np.tan(gam), 1e-3)
        r.s_air = s_glide + s_flare
        r.notes.append(f'glide {s_glide:.0f} m + flare {s_flare:.0f} m '
                       f'(h_f {h_flare:.2f} m, n {n_flare:.2f})')
        # 3) ground braking roll with spoilers deployed
        cl_g = 0.3 * clmax
        Vg, xg, it = v_td, 0.0, 0
        while Vg > 1.0 and it < 20000:
            it += 1
            q = 0.5 * RHO0 * Vg * Vg
            L = q * f.S * cl_g
            D = self.drag(Vg, cl_g, state)
            a = (-D - mu_brake * (W - L)) / m
            if a >= 0:
                break
            Vg += a * dt
            xg += Vg * dt
        r.s_ground = xg
        r.s_total = r.s_air + xg
        return r


def make_ground(cfg, key):
    return GroundRun(key, cfg['frame'], cfg['pack'])
