"""Engine pack and propeller models (layer-1 placeholders).

- Engine: rated shaft power, density-corrected availability, fixed BSFC.
  Layers 7 (2-D engine map) and 8 (BEMT prop) replace these.
- Static thrust: ideal actuator disk with a figure-of-merit eta_0.
"""
from dataclasses import dataclass
from math import pi
from typing import Optional

from .atmosphere import RHO0
from .engine_map import EngineMap
from .prop_bemt import PropBEMT


@dataclass(frozen=True)
class EnginePack:
    name: str
    n_eng: int
    P_kw_per_eng: float        # rated shaft power per engine, kW
    D_prop: float              # m
    bsfc_g_per_kwh: float = 400.0   # placeholder; layer-7 map replaces
    eta_p_cruise: float = 0.75      # placeholder; layer-8 BEMT replaces
    eta_0_static: float = 0.55      # prop figure of merit at J=0
    bsfc_pt_a: float = 0.55         # part-throttle SFC penalty model:
    bsfc_pt_b: float = 0.45         #   eff = bsfc*(a + b/throttle_frac)
    map: Optional[EngineMap] = None   # layer-7 2-D engine map (if set, used)
    prop: Optional[PropBEMT] = None   # layer-8 BEMT prop (if set, closes eta/J)

    @property
    def P_kw(self):
        return self.n_eng * self.P_kw_per_eng

    def bsfc_kg_per_Ws(self):
        """g/kWh -> kg/(W s)."""
        return self.bsfc_g_per_kwh * 1e-3 / 3.6e6

    def bsfc_eff_kg_per_Ws(self, throttle_frac):
        """SFC at part throttle; carbureted 2-strokes run rich off WOT.

        Placeholder until the layer-7 2-D map (BSFC(RPM, T)) lands.
        """
        frac = min(max(throttle_frac, 0.05), 1.0)
        return self.bsfc_kg_per_Ws() * (self.bsfc_pt_a + self.bsfc_pt_b / frac)

    def static_thrust(self, rho=RHO0):
        """Total static thrust, N (ideal actuator disk, one engine)."""
        a_d = pi * (self.D_prop / 2.0) ** 2
        p_shaft = self.P_kw_per_eng * 1000.0
        t0 = (2.0 * rho * a_d * (self.eta_0_static * p_shaft) ** 2) ** (1.0 / 3.0)
        return self.n_eng * t0

    def power_avail(self, rho):
        """Installed shaft power at density rho, W (~ density ratio)."""
        return self.P_kw * 1000.0 * (rho / RHO0)
