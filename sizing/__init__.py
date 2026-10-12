"""Layer-1 parametric sizing for the UAV preliminary-design pipeline.

Pipeline roadmap (later layers get their own packages under this repo):

  sizing/       layer 1  parametric sizing, design-space scan (this package)
  airfoil/      layer 2  XFOIL/XFLR5 2-D polars at Re = 0.4-1.5e6
                layer 3  flap 2-D sweep (chord, deflection, slot/gap/overlap)
  vsp/          layers 4-6  OpenVSP parametric geometry, VSPAERO aero DB,
                            actuator-disk propwash / engine-out checks
  engine/       layer 7  2-D engine maps: T(RPM, throttle), BSFC(RPM, T)
  prop/         layer 8  BEMT propeller: C_T(J), C_P(J), eta_p(J)
  jsbsim/       layers 9-11 6-DoF trim/dynamics, FCS, ground-run model
  cfd/          layer 12 targeted SU2/OpenFOAM (nacelle, cooling, flap, junction)
  structures/   layer 13 spar / engine mount / gear FEA from L(y) loads
"""
__version__ = "0.1.0"

from .mission import Aircraft, Airframe, MissionSpec
from .propulsion import EnginePack
