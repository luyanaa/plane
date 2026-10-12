"""Design configurations: engine packs, mission spec, four digital prototypes.

Engine sizing assumptions (rated shaft power for the 150/200-250cc class):
  150cc    ~ 11.0 kW   (e.g. 3W-150 / DLE-170 class)
  200-250cc ~ 16.5 kW  (e.g. DLE-222 class)
Layer 7 will replace these single points with 2-D maps.
"""
from dataclasses import replace

from .engine_map import map_150, map_200
from .mission import Airframe, MissionSpec
from .prop_bemt import make_prop
from .propulsion import EnginePack

PACK_150 = EnginePack('single150cc', 1, 11.0, 0.66, map=map_150())
PACK_200 = EnginePack('single200cc', 1, 16.5, 0.71, map=map_200())
PACK_2X150 = EnginePack('twin150cc', 2, 11.0, 0.51, map=map_150())
PACK_2X150_STOL = EnginePack('twin150ccSTOL', 2, 11.0, 0.53, map=map_150())

PACK_150 = replace(PACK_150, prop=make_prop(PACK_150))
PACK_200 = replace(PACK_200, prop=make_prop(PACK_200))
PACK_2X150 = replace(PACK_2X150, prop=make_prop(PACK_2X150))
PACK_2X150_STOL = replace(PACK_2X150_STOL, prop=make_prop(PACK_2X150_STOL))

PACKS = [PACK_150, PACK_200, PACK_2X150, PACK_2X150_STOL]

# Layer-2/3 XFOIL-informed C_Lmax (Re 8e5; plain flap):
#   naca2412: clean 1.47 / TO 10deg-30% 1.66 / land 30deg-30% 1.87
#   sd7062:   clean 1.69 / TO 10deg-20% 1.92 / land 30deg-25% 2.15
MISSION = MissionSpec()

PROTOTYPES = {
    'A': dict(
        name='A single150 long-endurance',
        frame=Airframe(S=3.2, AR=12.0, MTOW_kg=95.0,
                       C_Lmax_clean=1.47, C_Lmax_TO=1.66, C_Lmax_land=1.87),
        pack=PACK_150,
        layout='high wing, central tractor, single-slotted flap (XFOIL: plain 30% 10/30deg)'),
    'B': dict(
        name='B single200 heavy-lift',
        frame=Airframe(S=3.7, AR=10.5, MTOW_kg=112.0,
                       C_Lmax_clean=1.69, C_Lmax_TO=1.92, C_Lmax_land=2.15),
        pack=PACK_200,
        layout='high wing, central tractor, large flap (sd7062)'),
    'C': dict(
        name='C twin150 cruise',
        frame=Airframe(S=3.4, AR=11.0, MTOW_kg=105.0, n_nacelles=2,
                       C_Lmax_clean=1.47, C_Lmax_TO=1.66, C_Lmax_land=1.87),
        pack=PACK_2X150,
        layout='high wing, wing-root tractors (naca2412)'),
    'D': dict(
        name='D twin150 STOL',
        frame=Airframe(S=4.2, AR=9.5, MTOW_kg=115.0, n_nacelles=2,
                       C_Lmax_clean=1.69, C_Lmax_TO=1.92, C_Lmax_land=2.15),
        pack=PACK_2X150_STOL,
        layout='high wing, root tractors, large flap, propwash coupling (sd7062)'),
}
