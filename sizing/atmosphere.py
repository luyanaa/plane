"""ISA atmosphere (0-11 km) and physical constants."""
from math import exp, sqrt

G = 9.80665          # m/s^2
RHO0 = 1.225         # kg/m^3, sea level
_R = 287.0529        # J/(kg K)
_L = 0.0065          # K/m
_T0 = 288.15         # K
_P0 = 101325.0       # Pa


def isa(h_m):
    """Return (rho, T, p, a) at geometric altitude h_m."""
    if h_m <= 11000.0:
        t = _T0 - _L * h_m
        p = _P0 * (t / _T0) ** (G / (_R * _L))
    else:
        t = 216.65
        p = 22632.06 * exp(-G * (h_m - 11000.0) / (_R * t))
    rho = p / (_R * t)
    return rho, t, p, sqrt(1.4 * _R * t)
