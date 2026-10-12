"""Layer-1 drag model: C_D = C_D0 + C_L^2/(pi e AR).

C_D0 = C_fe * S_wet / S with a flat-plate equivalent skin friction C_fe.
Placeholder C_Lmax values until layers 2-3 (airfoil + flap 2-D) refine them.
"""
from dataclasses import dataclass
from math import pi, sqrt

from .atmosphere import G


@dataclass
class Polar:
    S: float
    AR: float
    S_wet: float
    C_fe: float = 0.009
    e: float = 0.80
    C_Lmax_clean: float = 1.35
    C_Lmax_TO: float = 1.80
    C_Lmax_land: float = 2.30

    def __post_init__(self):
        self.C_D0 = self.C_fe * self.S_wet / self.S
        self.k = 1.0 / (pi * self.e * self.AR)

    def C_D(self, c_l):
        return self.C_D0 + self.k * c_l * c_l

    @property
    def C_L_star(self):
        """C_L at best L/D (target cruise region is C_L = 0.4-0.9)."""
        return sqrt(self.C_D0 / self.k)

    @property
    def L_D_max(self):
        return 0.5 * sqrt(pi * self.e * self.AR / self.C_D0)

    def V_LDmax(self, m_kg, rho):
        return sqrt(2.0 * m_kg * G / (rho * self.S * self.C_L_star))

    def stall(self, m_kg, rho, c_lmax):
        return sqrt(2.0 * m_kg * G / (rho * self.S * c_lmax))
