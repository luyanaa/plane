"""Layer 7: 2-D engine map - T(RPM, throttle) and BSFC(RPM, T).

Placeholder surface calibrated to public rated points of the 150/200-250cc
2-stroke UAV engine class (e.g. 3W-157/DLE-170: 11 kW @ 7500 RPM; DLE-222:
16.5 kW @ 7500 RPM). Replace with bench data when available (the user's plan
explicitly: rated point first, sensible surface, bench-map later).

Physics:
  - torque curve: flat 0.5-1.05 x n_r, droop above, rise below (2-stroke shape)
  - throttle scales torque at fixed RPM: T = throttle * T_max(RPM)
  - propeller matching: n = n_r * (P/P_r)^(1/3)  (fixed-pitch prop, P ~ n^3)
  - BSFC(RPM, T): bsfc_r * f_load(T/T_max) * f_rpm(n/n_r), with
    f_load = a + (1-a)/load  (part-throttle penalty), f_rpm = 1 + c*(x-1)^2
"""
from dataclasses import dataclass
from math import pi


@dataclass(frozen=True)
class EngineMap:
    P_r_kw: float          # rated shaft power, kW
    n_r_rpm: float         # rated RPM
    bsfc_r_g_kwh: float    # rated BSFC
    n_idle_rpm: float = 1800.0
    n_max_rpm: float = 10500.0
    bsfc_load_a: float = 0.45     # f_load = a + (1-a)/load
    bsfc_rpm_c: float = 0.8       # f_rpm = 1 + c*(x-1)^2

    @property
    def T_r_Nm(self):
        return self.P_r_kw * 1000.0 / (2.0 * pi * self.n_r_rpm / 60.0)

    # ------------------------------------------------------------ torque
    def T_max(self, n_rpm):
        """Max torque (N m) at RPM - 2-stroke platform shape."""
        x = n_rpm / self.n_r_rpm
        if x < 0.5:
            tau = 0.55 + 0.9 * x
        elif x <= 1.05:
            tau = 1.0
        else:
            tau = max(1.0 - 2.2 * (x - 1.05) ** 2, 0.2)
        return self.T_r_Nm * tau

    def P_max_kw(self, n_rpm):
        return self.T_max(n_rpm) * (2.0 * pi * n_rpm / 60.0) / 1000.0

    # ------------------------------------------------------------ bsfc
    def bsfc(self, n_rpm, t_load):
        """BSFC (g/kWh) at RPM and torque load fraction (0..1)."""
        x = n_rpm / self.n_r_rpm
        load = min(max(t_load, 0.05), 1.0)
        f_load = self.bsfc_load_a + (1.0 - self.bsfc_load_a) / load
        f_rpm = 1.0 + self.bsfc_rpm_c * (x - 1.0) ** 2
        return self.bsfc_r_g_kwh * f_load * f_rpm

    # ------------------------------------------------------------ ops
    def operating_point(self, p_kw_per_eng):
        """Best operating point for required per-engine power P (kW).

        Returns dict(n_rpm, T_Nm, load, bsfc_g_kwh, throttle, saturated).
        """
        p = max(p_kw_per_eng, 1e-6)
        n = self.n_r_rpm * (p / self.P_r_kw) ** (1.0 / 3.0)
        n = min(max(n, self.n_idle_rpm), self.n_max_rpm)
        t_req = p * 1000.0 / (2.0 * pi * n / 60.0)
        t_max = self.T_max(n)
        saturated = t_req > t_max
        load = min(t_req / t_max, 1.0)
        return {'n_rpm': n, 'T_Nm': t_req, 'load': load,
                'bsfc_g_kwh': self.bsfc(n, load), 'throttle': load,
                'saturated': saturated}


def map_150():
    """150cc class: 11 kW @ 7500 RPM (3W-157 / DLE-170 family)."""
    return EngineMap(P_r_kw=11.0, n_r_rpm=7500.0, bsfc_r_g_kwh=400.0)


def map_200():
    """200-250cc class: 16.5 kW @ 7500 RPM (DLE-222 family)."""
    return EngineMap(P_r_kw=16.5, n_r_rpm=7500.0, bsfc_r_g_kwh=400.0)
