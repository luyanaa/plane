"""Layer-1 sizing: stall, cruise, climb, takeoff, range, endurance.

All metrics at MTOW / sea level unless noted. Range and endurance use the
closed-form max-L/D tracking profile (C_L held at C_L*; V ~ sqrt(m)), which is
the exact solution of this model and matches Breguet in the constant-V limit.
"""
from dataclasses import dataclass
from math import asin, cos, log, pi, sqrt

import numpy as np

from .aero import Polar
from .atmosphere import G, RHO0, isa
from .geometry import build_geometry
from .propulsion import EnginePack
from .weight import WeightModel


@dataclass(frozen=True)
class Airframe:
    S: float
    AR: float
    MTOW_kg: float
    taper: float = 0.60
    t_over_c: float = 0.13
    l_fus: float = 2.4
    d_fus: float = 0.32
    n_nacelles: int = 0
    # layer-2/3 XFOIL-informed C_Lmax (None -> MissionSpec defaults)
    C_Lmax_clean: float | None = None
    C_Lmax_TO: float | None = None
    C_Lmax_land: float | None = None


@dataclass(frozen=True)
class MissionSpec:
    payload_kg: float = 12.0
    h_cruise_m: float = 500.0
    reserve_frac: float = 0.10
    obstacle_m: float = 15.0
    V_R_factor: float = 1.1
    V_2_factor: float = 1.2
    mu_roll: float = 0.05
    C_Lmax_clean: float = 1.35
    C_Lmax_TO: float = 1.80
    C_Lmax_land: float = 2.30
    e_oswald: float = 0.80
    C_fe: float = 0.009
    V_s_land_max: float = 17.0
    ROC_min: float = 4.0
    fuel_frac_min: float = 0.03
    fuel_frac_max: float = 0.33
    s_TO_max: float = 250.0


class Aircraft:
    def __init__(self, frame, pack, spec):
        self.frame, self.pack, self.spec = frame, pack, spec
        self.geo = build_geometry(frame.S, frame.AR, frame.taper, frame.t_over_c,
                                  frame.l_fus, frame.d_fus, n_nacelles=frame.n_nacelles)
        self.m_empty, self.weight_breakdown = WeightModel().empty_mass(
            self.geo, pack.n_eng, pack.P_kw_per_eng)
        self.m_fuel = frame.MTOW_kg - self.m_empty - spec.payload_kg
        self.polar = Polar(frame.S, frame.AR, self.geo.S_wet,
                           C_fe=spec.C_fe, e=spec.e_oswald,
                           C_Lmax_clean=frame.C_Lmax_clean or spec.C_Lmax_clean,
                           C_Lmax_TO=frame.C_Lmax_TO or spec.C_Lmax_TO,
                           C_Lmax_land=frame.C_Lmax_land or spec.C_Lmax_land)
        self.W = frame.MTOW_kg * G

    # ------------------------------------------------------------- metrics
    def metrics(self, rho=None):
        rho = RHO0 if rho is None else rho
        p, s, f = self.pack, self.spec, self.frame
        W, S = self.W, self.geo.S
        c = p.bsfc_kg_per_Ws()
        m = {}

        m['b_m'] = self.geo.b
        m['MAC_m'] = self.geo.MAC
        m['MTOW_kg'] = f.MTOW_kg
        m['W_S_kgm2'] = f.MTOW_kg / S
        m['M_empty_kg'] = self.m_empty
        m['M_fuel_kg'] = max(self.m_fuel, 0.0)
        m['fuel_frac'] = self.m_fuel / f.MTOW_kg

        # stall speeds at MTOW, sea level (uses frame/XFOIL C_Lmax overrides)
        m['V_s_clean'] = self.polar.stall(f.MTOW_kg, rho, self.polar.C_Lmax_clean)
        m['V_s_TO'] = self.polar.stall(f.MTOW_kg, rho, self.polar.C_Lmax_TO)
        m['V_s_land'] = self.polar.stall(f.MTOW_kg, rho, self.polar.C_Lmax_land)

        # cruise point (max L/D) at MTOW; BEMT prop closes eta_p(J)
        cl_star, ld = self.polar.C_L_star, self.polar.L_D_max
        m['C_L_cruise'] = cl_star
        m['L_D_max'] = ld
        v_cr = self.polar.V_LDmax(f.MTOW_kg, rho)
        m['V_cruise'] = v_cr
        if p.prop and p.map:
            eta_p = p.eta_p_cruise
            for _ in range(6):
                p_cr = W * v_cr / (ld * eta_p)
                pt = self._prop_point(v_cr, p_cr, rho)
                eta_p = pt['eta']
            m['P_cruise_kW'] = p_cr / 1000.0
            m['eta_p_cruise'] = eta_p
            m['J_cruise'] = pt['J']
            m['n_cruise'] = pt['n_rpm']
            n_rpm = pt['n_rpm']
            load = min(p_cr / p.n_eng / (p.map.P_max_kw(n_rpm) * 1000.0), 1.0)
            m['throttle_frac'] = load
            m['fuel_flow_cruise_gh'] = p.map.bsfc(n_rpm, load) * (p_cr / 1000.0)
        elif p.map:
            p_cr = W * v_cr / (ld * p.eta_p_cruise)
            m['P_cruise_kW'] = p_cr / 1000.0
            op = p.map.operating_point(p_cr / 1000.0 / p.n_eng)
            m['throttle_frac'] = op['load']
            m['fuel_flow_cruise_gh'] = op['bsfc_g_kwh'] * (p_cr / 1000.0)
        else:
            p_cr = W * v_cr / (ld * p.eta_p_cruise)
            m['P_cruise_kW'] = p_cr / 1000.0
            m['throttle_frac'] = p_cr / (p.P_kw * 1000.0)
            m['fuel_flow_cruise_gh'] = p.bsfc_eff_kg_per_Ws(m['throttle_frac']) * p_cr * 3600.0 * 1000.0

        # best climb rate, sea level, full power (excess-power sweep)
        v_grid = np.linspace(1.1 * m['V_s_clean'], 1.6 * v_cr, 60)
        c_l = W / (0.5 * rho * S * v_grid ** 2)
        d_vec = 0.5 * rho * S * v_grid ** 2 * (self.polar.C_D0 + self.polar.k * c_l ** 2)
        roc = (p.power_avail(rho) * p.eta_p_cruise - d_vec * v_grid) / W
        m['ROC'] = float(np.max(roc))
        m['V_best_ROC'] = float(v_grid[int(np.argmax(roc))])

        # climb segment: sea-level ROC, average available power over the climb
        rho_h = isa(s.h_cruise_m)[0]
        p_avg = p.power_avail(0.5 * (rho + rho_h))
        m['t_climb_s'] = s.h_cruise_m / m['ROC']
        m['fuel_climb_kg'] = c * p_avg * m['t_climb_s']

        # takeoff
        s_g, s_air, t0_w = self._takeoff(rho)
        m['s_G_m'] = s_g
        m['s_air_m'] = s_air
        m['s_TO_m'] = s_g + s_air
        m['T0_W'] = t0_w
        m['P_W_kW_kg'] = p.P_kw / f.MTOW_kg

        # range & endurance: closed-form max-L/D tracking, or numerical
        # integration with the layer-7 engine map when one is attached
        fuel_cr = self.m_fuel - m['fuel_climb_kg'] - s.reserve_frac * self.m_fuel
        m0 = f.MTOW_kg - m['fuel_climb_kg']
        if fuel_cr > 1e-6:
            m1 = m0 - fuel_cr
            if p.map:
                t_e, dist = self._integrate_mission(m0, m1, rho)
                m['endurance_h'] = t_e / 3600.0
                m['range_km'] = dist / 1000.0
            else:
                u = sqrt(2.0 * G / (rho * S * cl_star))           # V = u*sqrt(m)
                m_avg = 0.5 * (m0 + m1)
                frac_avg = (m_avg * G * u * sqrt(m_avg) / (ld * p.eta_p_cruise)) / (p.P_kw * 1000.0)
                c_eff = p.bsfc_eff_kg_per_Ws(frac_avg)
                k = c_eff * G * u / (p.eta_p_cruise * ld)         # dm/dt = -k*m^1.5
                m['range_km'] = (p.eta_p_cruise * ld / (c_eff * G)) * log(m0 / m1) / 1000.0
                m['endurance_h'] = (2.0 / (k * sqrt(m0))) * (sqrt(m0 / m1) - 1.0) / 3600.0
        else:
            m['range_km'] = 0.0
            m['endurance_h'] = 0.0

        # design-constraint gates
        m['feasible'] = (
            m['V_s_land'] <= s.V_s_land_max
            and m['ROC'] >= s.ROC_min
            and s.fuel_frac_min <= m['fuel_frac'] <= s.fuel_frac_max
            and np.isfinite(m['s_TO_m']) and m['s_TO_m'] <= s.s_TO_max
            and m['range_km'] > 0.0
        )
        return m

    # ------------------------------------------------------------- mission
    def _bisect(self, f, lo, hi, tol=1e-4):
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            if f(mid) == 0.0:
                return mid
            if f(lo) * f(mid) < 0:
                hi = mid
            else:
                lo = mid
            if hi - lo < tol:
                break
        return 0.5 * (lo + hi)

    def _prop_point(self, V, p_req_w, rho):
        """Solve n so the BEMT prop absorbs p_req at speed V; returns dict."""
        p = self.pack
        prop = p.prop

        def f(n):
            _, _, pw = prop.solve(V, n)
            return pw - p_req_w

        n = self._bisect(f, 5.0, 190.0)
        t, q, pw = prop.solve(V, n)
        ct = t / (rho * n ** 2 * prop.D ** 4)
        cp = pw / (rho * n ** 3 * prop.D ** 5)
        j = V / (n * prop.D)
        return {'n_rps': n, 'n_rpm': n * 60.0, 'T': t, 'eta': j * ct / cp if cp > 0 else 0.0,
                'J': j, 'CT': ct, 'CP': cp}

    def _prop_full(self, V, rho):
        """Full-throttle point: BEMT absorption = engine torque-curve power."""
        p = self.pack

        def f(n):
            _, _, pw = p.prop.solve(V, n)
            return pw - p.map.P_max_kw(n * 60.0) * 1000.0 * (rho / RHO0)

        n = self._bisect(f, 5.0, 190.0)
        t, q, pw = p.prop.solve(V, n)
        return t, n

    def _integrate_mission(self, m0, m1, rho):
        """Numerical cruise/loiter on the max-L/D profile with engine map.

        Returns (time_s, distance_m). Step 120 s is accurate to <1% here.
        """
        p, polar = self.pack, self.polar
        cl_star, ld = polar.C_L_star, polar.L_D_max
        u = sqrt(2.0 * G / (rho * self.geo.S * cl_star))
        eta = p.eta_p_cruise
        m, t, dist, dt = m0, 0.0, 0.0, 120.0
        while m > m1:
            v = u * sqrt(m)
            pw = m * G * v / (ld * eta)                      # W, all engines
            if p.prop:
                pt = self._prop_point(v, pw, rho)
                eta = pt['eta']
                pw = m * G * v / (ld * eta)
                n_rpm, t_shaft = pt['n_rpm'], pw / (2.0 * pi * pt['n_rps'])
            else:
                n_rpm = p.map.operating_point(pw / 1000.0 / p.n_eng)['n_rpm']
                t_shaft = pw / (2.0 * pi * n_rpm / 60.0)
            load = min(pw / p.n_eng / (p.map.P_max_kw(n_rpm) * 1000.0), 1.0)
            bsfc = p.map.bsfc(n_rpm, load)
            gph = bsfc * (pw / 1000.0)                       # total g/h
            m -= gph / 3.6e6 * dt                            # kg
            t += dt
            dist += v * dt
        return t, dist

    # ------------------------------------------------------------- takeoff
    def _takeoff(self, rho):
        """Ground roll + air distance to obstacle (layer-11 model refines this)."""
        p, s, f = self.pack, self.spec, self.frame
        W, m_ac, S = self.W, f.MTOW_kg, self.geo.S
        v_s_to = self.polar.stall(m_ac, rho, self.polar.C_Lmax_TO)
        v_r = s.V_R_factor * v_s_to
        v_2 = s.V_2_factor * v_s_to
        v_cr = self.polar.V_LDmax(m_ac, rho)
        p_tot = p.P_kw * 1000.0
        t0 = p.static_thrust(rho)
        # takeoff stays on the momentum-disk linear model: BEMT is not valid
        # for J < 0.25 (heavy-load regime, no-stall polar over-predicts)
        t_full_cr = p.eta_p_cruise * p_tot / v_cr

        def thrust(v):
            lin = t0 + (t_full_cr - t0) * (v / v_cr)
            return np.where(v < v_cr, lin, p.eta_p_cruise * p_tot / np.maximum(v, 1e-6))

        # ground roll: s = integral of V/a dV (a = T - D - mu(W-L))/m)
        vg = np.linspace(0.05, v_r, 240)
        q = 0.5 * rho * vg ** 2
        c_lr = s.C_Lmax_TO * np.minimum(1.0, (vg / v_r) ** 2)
        lift = q * S * c_lr
        drag = q * S * (self.polar.C_D0 + self.polar.k * c_lr ** 2)
        acc = (thrust(vg) - drag - s.mu_roll * (W - lift)) / m_ac
        if acc[0] <= 0.0:
            return float('inf'), 0.0, t0 / W
        s_g = float(np.trapezoid(vg / acc, vg))

        # air phase: accelerate V_R -> V_2, then climb at V_2 to obstacle height
        va = np.linspace(v_r, v_2, 160)
        qa = 0.5 * rho * va ** 2
        drag_a = qa * S * (self.polar.C_D0 + self.polar.k * s.C_Lmax_TO ** 2)
        acc_a = (thrust(va) - drag_a) / m_ac
        gam = np.arcsin(np.clip((thrust(va) - drag_a) / W, -1.0, 1.0))
        mid = 0.5 * ((va / acc_a)[1:] + (va / acc_a)[:-1])
        gmid = 0.5 * (gam[1:] + gam[:-1])
        h_acc = np.cumsum(mid * np.sin(gmid) * np.diff(va))
        s_acc = np.cumsum(mid * np.cos(gmid) * np.diff(va))
        idx = np.searchsorted(h_acc, s.obstacle_m)
        if idx < len(h_acc):
            if idx == 0:
                s_air = s_acc[0] * (s.obstacle_m / h_acc[0]) if h_acc[0] > 0 else float('inf')
            else:
                s_air = float(np.interp(s.obstacle_m, h_acc[idx - 1:idx + 1], s_acc[idx - 1:idx + 1]))
        elif gam[-1] > 0.02:
            t_rem = (s.obstacle_m - h_acc[-1]) / (v_2 * np.sin(gam[-1]))
            s_air = float(s_acc[-1] + v_2 * np.cos(gam[-1]) * t_rem)
        else:
            s_air = float('inf')
        return s_g, s_air, t0 / W
