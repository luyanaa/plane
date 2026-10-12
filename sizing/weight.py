"""Simple empty-weight buildup.

Placeholder for layer-13 FEA-refined structure mass. Component densities were
tuned so the empty fraction lands ~0.50-0.60 of MTOW for the 70-140 kg class,
consistent with composite endurance UAVs.
"""
from dataclasses import dataclass


@dataclass
class WeightModel:
    wing_kg_m2: float = 7.0
    fuselage_kg: float = 8.0
    tail_kg_m2: float = 2.0
    landing_gear_kg: float = 3.5
    engine_kg_per_kw: float = 0.45
    engine_fixed_kg: float = 1.0
    fuel_system_kg: float = 1.2
    prop_kg: float = 0.40
    avionics_kg: float = 2.2
    servos_kg: float = 1.6
    payload_bay_kg: float = 1.4
    electrical_kg: float = 1.5
    misc_fraction: float = 0.10

    def empty_mass(self, geo, n_eng, p_kw_per_eng):
        items = {
            'wing': self.wing_kg_m2 * geo.S,
            'fuselage': self.fuselage_kg,
            'HTP': self.tail_kg_m2 * geo.S_h,
            'VTP': self.tail_kg_m2 * geo.S_v,
            'landing gear': self.landing_gear_kg,
            'engine+mount': n_eng * (self.engine_kg_per_kw * p_kw_per_eng + self.engine_fixed_kg),
            'fuel system': self.fuel_system_kg,
            'propellers': self.prop_kg * n_eng,
            'avionics/FCS': self.avionics_kg,
            'servos/wiring': self.servos_kg,
            'payload bay': self.payload_bay_kg,
            'electrical': self.electrical_kg,
        }
        total = (1.0 + self.misc_fraction) * sum(items.values())
        return total, items
