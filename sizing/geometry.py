"""Parametric geometry for layer-1 sizing.

Fuselage/empennage are sized by volume coefficients; wetted area feeds the
C_D0 estimate. Placeholder until OpenVSP geometry (layer 4) replaces it.
"""
from dataclasses import dataclass
from math import pi, sqrt


@dataclass(frozen=True)
class Geometry:
    S: float        # wing area, m^2
    AR: float
    b: float        # span, m
    c_root: float   # m
    MAC: float      # m
    S_h: float      # HTP area, m^2
    S_v: float      # VTP area, m^2
    S_wet: float    # total wetted area, m^2
    taper: float
    t_over_c: float


def build_geometry(S, AR, taper=0.60, t_over_c=0.13, l_fus=2.4, d_fus=0.32,
                   V_h=0.65, V_v=0.035, l_h=1.25, l_v=1.15, n_nacelles=0):
    """Tail areas from volume coefficients:
        S_h = V_h * S * MAC / l_h,   S_v = V_v * S * b / l_v
    Wetted area = 2.05*S (wing) + 2*(S_h+S_v) + fuselage skin + nacelles.
    """
    b = sqrt(S * AR)
    c_root = 2.0 * S / (b * (1.0 + taper))
    mac = (2.0 / 3.0) * c_root * (1.0 + taper + taper ** 2) / (1.0 + taper)
    s_h = V_h * S * mac / l_h
    s_v = V_v * S * b / l_v
    s_wet = (2.05 * S + 2.0 * (s_h + s_v)
             + 0.82 * pi * d_fus * l_fus + n_nacelles * 0.45)
    return Geometry(S, AR, b, c_root, mac, s_h, s_v, s_wet, taper, t_over_c)
