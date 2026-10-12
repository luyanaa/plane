# UAV preliminary design pipeline (70–140 kg class)

Experiments run on **tinybox** (`ssh tinybox`, 48-core x86_64) in `~/Documents/plane`.
This directory is the authoring mirror; sync with
`rsync -az --delete --exclude .venv --exclude results --exclude __pycache__ ./ tinybox:~/Documents/plane/`
Run on tinybox with `~/.venv/bin/python`.

## Pipeline roadmap

```
Engine map ──► Propeller BEMT ──► 2-D airfoil polars
     │               │                   │
     └───────────────┴───────► OpenVSP geometry ──► VSPAERO ──► aero coefficient DB
                                                                      │
                                        ┌─────────────────────────────┘
                                        ▼
                                   JSBSim 6-DoF ──► Autopilot ──► Mission simulation
                                                                      │
                                        Range / endurance / STOL / engine-out metrics
                                                                      │
                          outer loop: S, AR, airfoil, flap, engine location,
                          prop diameter, fuel, MTOW ──► optimizer ──► Pareto front
```

| Layer | Tool | Status |
|---|---|---|
| 1 | Parametric sizing + design-space scan | **done** (this repo, `sizing/`) |
| 2 | Airfoil polars at Re = 0.4–1.5e6, transition sweep | **done** (`scripts/airfoil_scan.py`) |
| 3 | Flap 2-D sweep (chord 20–35%, δ 0–40°, slot/gap) | **done** plain flap (`scripts/flap_scan.py`); slot/gap/overlap → MSES or CFD target #4 |
| 4 | Parametric geometry (single vs twin tractor) | **done** (`scripts/vsp_geometry.py` → `vsp_models/*.vsp3`) |
| 5 | Whole-aircraft aero DB + stability derivatives | **done** (`scripts/vspaero_sweep.py`) |
| 6 | Propeller: actuator disk, propwash, engine-out yaw | **done** momentum disk (`scripts/prop_model.py`); VSPAERO-internal rotor API unavailable in 3.53 |
| 7 | Engine 2-D maps: T(RPM, throttle), BSFC(RPM, T) | **done** placeholder surface (`sizing/engine_map.py`), bench-map pending |
| 8 | Propeller BEMT: C_T(J), C_P(J), η_p(J) | **done** (`sizing/prop_bemt.py`); J<0.25 zone stays on momentum model |
| 9 | 6-DoF trim, modes, V_MC (engine-out) | **done** (`sizing/sixdof.py`, `scripts/sixdof_run.py`); JSBSim XML generated, pip 1.3.1 load blocked |
| 10 | FCS: inner/attitude/outer loops, failure cases | **done** (`sizing/autopilot.py`, `scripts/fcs_run.py`); OEI remains uncontrolled |
| 11 | Ground-run TO/landing model | **done** (`sizing/ground.py`, `scripts/ground_run.py`) |
| 12 | Targeted CFD only: nacelle, cooling duct, prop–wing, flap, junction | SU2/OpenFOAM |
| 13 | Structure: M_root = ∫y·L(y)dy, spar/mount/gear FEA, +3g/−1.5g, twin yaw-load | FEA |

## Toolchain on tinybox (installed 2026-10-04)

- **OpenVSP 3.53.1 + VSPAERO 7.2.2**: official Ubuntu 24.04 deb from
  openvsp.org → `sudo dpkg -i /tmp/openvsp.deb` (deps via `sudo apt-get install -f`).
  Binaries: `/usr/local/bin/vsp`, `/usr/local/bin/vspaero`.
  Python API: packages under `/opt/OpenVSP/python/` were copied into the venv
  site-packages, so `import openvsp` works from `.venv/bin/python`
  (`openvsp.GetVSPVersion()` → "OpenVSP 3.53.1"; 1989 API functions).
  Bonus: `/opt/OpenVSP/python/vsp_airfoils/` (E387, Clark Y samples) and
  `PyVSP`/`AvlPy` toolkits.
- **XFOIL 6.99**: `sudo apt-get install xfoil` → `/usr/bin/xfoil`.
  The Ubuntu build **requires an X display even in batch** (aborts with
  "Cannot open display" otherwise and writes no polar points) — run under
  `Xvfb :99` (`sudo apt-get install xvfb`; `scripts/airfoil_scan.py` starts it
  automatically). Transition forcing via `VPAR` → `XTR` (the old `TFOR`
  command does not exist in 6.99).
- **conda** at `/home/yanlu/anaconda3` (channels: defaults + SJTU conda-forge
  mirror; neither `xfoil` nor `openvsp` exist on conda-forge — only
  `libxfoil` (py39) and `openvsp-code-eli` meta packages).
- Network note: GitHub/UIUC/airfoiltools are slow or blocked from tinybox;
  fetch data locally (Mac) and rsync up. UIUC airfoil dats live in
  `data/airfoils/`.

## Layer 11 — ground run: takeoff and landing (tinybox, 2026-10-04)

`sizing/ground.py`: time-marched ground roll (m dV/dt = T(V) − D − μ(W−L)),
air phase to 15 m, and a landing model (5° glideslope at V_app → flare →
braking roll with spoilers deployed, C_L = 0.3·C_Lmax so weight is on the
wheels). Flap states come from the layer-3 XFOIL lookup; the twins get the
layer-6 propwash flap-span q-boost (q_ratio 3.12 C / 3.84 D) applied to the
flap lift increment. `scripts/ground_run.py` → `results/ground_run.csv` + png.

| proto | state | C_Lmax_eff | V_s | s_ground | s_TO | landing total |
|---|---|---|---|---|---|---|
| A | clean | 1.47 | 18.0 | 62.9 | 114.8 | 302.5 |
| A | TO flap | 1.66 | 16.9 | 56.3 | 108.3 | 275.1 |
| A | landing flap | 1.87 | 15.9 | 51.0 | 105.7 | 248.3 |
| B | clean→landing | 1.69→2.15 | 16.9→15.0 | 46.1→37.5 | 87.9→81.2 | 273.2→227.5 |
| C | clean | 1.47 | 18.3 | 36.6 | 76.8 | 356.7 |
| C | TO flap | **1.85** | 16.4 | 29.9 | 69.9 | 311.5 |
| C | landing flap | **2.27** | 14.8 | 25.0 | **57.4** | 268.1 |
| D | clean | 1.69 | 16.1 | 31.2 | 71.6 | 264.0 |
| D | TO flap | **2.23** | 14.0 | 24.5 | **56.9** | 239.0 |
| D | landing flap | **2.78** | **12.6** | 20.4 | **49.3** | 219.1 |

**Answer: yes, the flap is worth it, and the propwash is what makes it work
for the twins.**

1. **Takeoff flap saves 5-15 m** on the singles (A 114.8→108.3) but **20-32%
   on the twins** (C 76.8→57.4, D 71.6→49.3) because the layer-6 slipstream
   multiplies the flap increment over the covered span
   (C_Lmax_eff: C → 2.27, D → 2.78 at landing flap).
2. **Landing flap shortens the landing distance 17-25%** across the board
   (A 302→248 m, C 357→268 m, D 264→219 m) via a lower V_app
   (D: 20.9→18.6 m/s) and hence a shorter ground roll.
3. **Landing is the sizing case, not takeoff**: landing distances run
   2-4× the takeoff distances — the 15 m obstacle approach dominates. If a
   short-field requirement bites, the levers are a steeper approved
   glideslope, thrust reverse (not available on these props), or a lower
   V_app margin (1.3 → 1.25 V_s).
4. **Cross-check vs layer 1**: A's s_TO 108.3 m here vs 108.7 m from the
   layer-1 sizing model — the independent time-march agrees within 0.4%.
5. Twin propwash saves 3.1 m (C) and 8.5 m (D) of takeoff roll versus the
   same geometry with no slipstream — small in absolute terms but ~15% of
   D's ground roll, and it is the reason D's V_s drops to 12.6 m/s.

## Layer 10 — flight control system, failure cases (tinybox, 2026-10-04)

`sizing/autopilot.py`: cascade controller on the layer-9 6-DoF — inner rate
loops (p/q/r) → attitude (φ/θ) → outer (V/h/ψ), with coordinated-turn rudder
(r = g·tanφ/V) and a throttle PI on airspeed. Actuators model first-order lag
(50 ms), rate limits (60°/s) and position limits (de ±25°, da ±20°, dr ±25°).
`sizing/sixdof.py` gained earth-frame wind (aero uses air-relative velocity,
kinematics ground-relative). Artifacts: `results/fcs_metrics.csv`.

Results (40 s, pass = Δh<25 m, Δψ<10°, ΔV<2.5 m/s, n_z<3):

| scenario | A | B | C | D |
|---|---|---|---|---|
| cruise hold | ✓ (Δh 7 m, ΔV 0.9) | ✓ | ✓ | ✓ |
| 8 m/s crosswind | ✓ (Δψ 3.0°, β −0.1°) | ✓ | ΔV 2.6 | ✓ |
| discrete gust 10 m/s | ΔV 8.7 | ΔV 9.7 | ΔV 10.4 | ΔV 9.0 |
| GPS loss at 20 s | Δh 55-66 m | Δh 66 | Δh 63 | Δh 62 |
| airspeed bias +3 m/s | ΔV 5.0 | 5.2 | 6.2 | 5.1 |
| **engine-out (C/D)** | — | — | **uncontrolled** | **uncontrolled** |

Findings:

1. **Cruise, attitude and crosswind hold work** on all four — the cascade
   closes with hand-tuned gains; 8 m/s crosswind is rejected with ≤4.8°
   heading excursion and essentially zero steady sideslip. This validates
   the layer-9 modes (Dutch-roll ζ 0.13 was controllable).
2. **Gust rejection is speed-limited**: a 10 m/s discrete gust produces
   ~9-10 m/s airspeed excursion — the throttle PI (kV 0.08) cannot add
   power fast enough at these low power-to-weight ratios (P/W 0.11-0.16
   kW/kg). Realistically mitigated by larger margins in the design point or
   feed-forward gust sensing.
3. **GPS-loss altitude drift 55-66 m in 20 s** is the expected consequence of
   the design (no baro-only hold in this simple controller); an INS/baro
   fallback is a layer-10 FCS requirement, not a surprise.
4. **Engine-out is still uncontrollable with the FCS** (C: 983° yaw, β −59°;
   D better: 183° yaw, heavy rudder saturation 0.12). This *confirms* layer 9
   dynamically: the current VTP + rudder has no OEI authority, and the
   autopilot's rudder command is already saturating at 25°. The fix is
   aerodynamic (VTP Vv 0.035 → ~0.06), not control-law tuning.

## Layer 9 — 6-DoF dynamics: trim, modes, engine-out (tinybox, 2026-10-04)

`sizing/sixdof.py` (nonlinear rigid-body EOM, RK4) + `scripts/sixdof_run.py`
use the same data as the JSBSim XMLs: VSPAERO aero tables + analytic control
and rate-damping derivatives, engine-map + BEMT thrust grid (momentum-anchored,
consistent with takeoff). Artifacts: `results/sixdof_trim.csv`,
`sizdof_modes.csv`, `sixdof_vmc.csv`.

**JSBSim integration note**: JSBSim 1.3.1 (pip) aircraft XMLs are generated
(`jsbsim_root/`, c172x-convention: nested thruster, WIND axes, AERORP metrics)
but the binding silently rejects non-bundled models after extensive
bisection (even a byte-identical c172x under a different name fails);
diagnosed root-cause candidates documented in the repo. The 6-DoF math is
identical to what JSBSim would execute with these tables.

Trim (level at 500 m):

| proto | cruise | approach | V_bestROC |
|---|---|---|---|
| A | α 9.1°, de −0.8°, thr 0.19 | ... | ... |
| B | α 8.9°, de −2.3°, thr 0.22 | ... | ... |
| C | α 8.5°, de −0.8°, **thr 0.60** | ... | ... |
| D | α 8.5°, de −2.4°, thr 0.30 | ... | ... |

(C's cruise throttle 0.60 reproduces the layer-8 twin-prop mismatch.)

Modes around cruise (longitudinal: short period / phugoid; lateral: Dutch
roll / roll / spiral):

| proto | short period | phugoid | Dutch roll | roll τ | spiral τ |
|---|---|---|---|---|---|
| A | (overdamped, τ 0.29 s) | 0.089 Hz, ζ 0.66 | 0.139 Hz, ζ 0.13 | 0.23 s | 13.4 s |
| B | 0.351 Hz, ζ 0.71 | (not resolved) | 0.132 Hz, ζ 0.14 | 0.22 s | 12.2 s |
| C | 0.206 Hz, ζ 0.86 | 0.058 Hz, ζ 0.01 | 0.161 Hz, ζ 0.13 | 0.22 s | 22.2 s |
| D | 0.413 Hz, ζ 0.69 | (not resolved) | 0.122 Hz, ζ 0.14 | 0.22 s | 12.4 s |

(Dutch-roll ζ ≈ 0.13 is light but typical for this VTP volume; C's phugoid
ζ 0.01 is nearly neutral, expected for a clean high-AR layout.)

Engine-out (C/D): left engine 0, right full, max rudder, 20 s. **Both
prototypes are uncontrollable across the whole envelope**: max bank 19-27°,
yaw accumulates 5-6 turns, altitude loss 540-730 m in 20 s at every speed up
to 2×V_R. This confirms and strengthens the layer-6 static V_MC estimate
(24-27 m/s): the dynamic model finds no control authority anywhere in the
flight envelope with the current VTP (Vv=0.035) and rudder (25% chord, 25°).
**Design driver: C/D need VTP volume ~0.06, a larger rudder, or FCS
asymmetric-thrust assist before engine-out certification is credible.**

## Layer 8 — propeller BEMT closure (tinybox, 2026-10-04)

`sizing/prop_bemt.py`: blade-element momentum theory with Prandtl tip/hub
loss, constant-pitch geometry (P/D design variable), linear 2-D polar.
`sizing/mission.py` now closes the loop per flight condition: BEMT prop
absorption = engine torque-curve power → n; η_p(J) = J·C_T/C_P feeds the
cruise/loiter integrator and the engine map feeds BSFC(n, load).
Takeoff stays on the momentum-disk model (BEMT is not valid for J<0.25:
heavy-load regime, no-stall polar over-predicts — static BEMT gave 871 N vs
313 N momentum reference). `scripts/prop_bemt.py` → `results/prop_curves.png`,
`results/prop_operating.csv`.

| proto | prop | J_TO | J_climb | J_cruise | η_cruise | n_cruise | E (BEMT) | E (η=0.75) |
|---|---|---|---|---|---|---|---|---|
| A | 0.66 m, P/D 0.9 | 0.16 | 0.17 | 0.31 | 0.650 | 6534 | **10.3 h** | 13.6 h (−24%) |
| B | 0.71 m, P/D 0.9 | 0.14 | 0.15 | 0.31 | 0.650 | 6340 | 9.8 h | 13.2 h (−26%) |
| C | 0.51 m, P/D 1.0 | 0.20 | 0.21 | 0.27 | 0.620 | **9648** | **4.1 h** | 8.7 h (−52%) |
| D | 0.53 m, P/D 1.0 | 0.17 | 0.18 | 0.26 | 0.606 | 9480 | 4.5 h | 8.9 h (−50%) |

Findings:

1. **Cruise η_p is 0.61-0.65, not 0.75** — endurance drops 24-52% vs the
   constant-η baseline. The assumed 0.75 was optimistic.
2. **The fixed-pitch props are mis-matched, and the twins are the worst**:
   C/D (small D, P/D 1.0) must spin ~9500 RPM at cruise → engine load
   fraction 10-12% → BSFC ~2400 g/kWh. Single-vs-twin endurance ratio grows
   to 2.5× (was 1.56× with matched-η). Design levers: larger twin props
   (D↑), lower pitch, or cruise at higher load.
3. **BEMT validity boundary J<0.25**: static/low-J thrust is the
   momentum-disk FM model (layer-1); BEMT governs climb/cruise. The
   low-J part of the C_T/C_P curves is reported but shaded in the plot.
4. Takeoff unchanged (momentum model): s_TO 108/80/55/53 m — the BEMT
   closure did not degrade the STOL numbers.

## Layer 7 — engine 2-D maps (placeholder surface, tinybox, 2026-10-04)

`sizing/engine_map.py`: T_max(RPM) 2-stroke torque curve (flat 0.5-1.05 n_r,
droop above, rise below) + BSFC(RPM, T) = bsfc_r·f_load·f_rpm with
f_load = 0.45 + 0.55/load (part-throttle) and f_rpm = 1 + 0.8·(n/n_r−1)².
Propeller matching picks n = n_r·(P/P_r)^⅓ (fixed-pitch P∝n³). Calibrated to
public rated points: 150cc 11 kW @ 7500 RPM (3W-157/DLE-170 class), 200-250cc
16.5 kW @ 7500 (DLE-222 class), bsfc_r 400 g/kWh. Replace with bench data.
`scripts/engine_map.py` → `results/engine_map.png`,
`results/engine_operating_points.csv`.

Operating points (cruise → loiter):

| proto | P_cr | RPM | load | BSFC | endurance | range |
|---|---|---|---|---|---|---|
| A single150 | 1.78 kW | 4091 | 30% | 1072 g/kWh | **13.6 h** | 998 km |
| B single200 | 2.35 kW | 3916 | 27% | 1167 | 13.2 h | 987 km |
| C twin150 | 2.19 kW | 3474 | 22% | 1442 | 8.7 h | 658 km |
| D twin150 STOL | 2.51 kW | 3638 | 24% | 1336 | 8.9 h | 683 km |

Answer to "单发大排量 vs 双发小排量": **single 150 wins endurance by 1.56×**
(A 13.6 h vs C 8.7 h, range 998 vs 658 km) — two engines at 22% load burn
~35% worse SFC than one at 30% load. The twin's cruise penalty is intrinsic
part-throttle SFC, not layout; the fix is engine sizing (2× smaller engines)
or cruise at higher load. This replaces the layer-1 flat-penalty estimate
(which said 1.82×); both models agree on the direction, the map is milder
and defensible until bench data lands.

## Layer 6 — propeller actuator disk, propwash, engine-out (tinybox, 2026-10-04)

`scripts/prop_model.py`: momentum-theory actuator disk (v_i = T/(2ρA(V+v_i)))
with the layer-1 thrust model and layer-5 VSPAERO baseline. VSPAERO 3.53's
Python API does not expose rotor operating points (GUI/CLI only), so the
disk is modeled in Python — the same physics VSPAERO's rotor uses.
Artifacts: `results/propwash.csv`.

Propwash (slipstream dynamic-pressure ratio on the affected surface):

| proto | case | q_slip/q | area in jet | effective q gain |
|---|---|---|---|---|
| A single | TO tail | 5.4× | 41% HTP | +177% |
| A single | cruise tail | 4.2× | 41% HTP | +134% |
| B single | TO tail | 7.0× | 37% HTP | +219% |
| C twin | TO inner wing | 7.7× | 47% flap span | +312% |
| D twin | TO inner wing | 9.0× | 48% flap span | +384% |

Findings:

1. **Low-speed pitch (single)**: propwash roughly triples tail control power
   at V_R (q-boosted HTP) while thrust-line offset gives ΔCm ≈ −0.05..−0.06
   (nose-down, manageable by the boosted elevator). The central tractor
   layout is sound for STOL pitch authority.
2. **Twin STOL**: the inner-flap q boost (+312..384%) is the mechanism behind
   D's short-field claim — flap lift and C_Lmax on the covered span scale
   with q. This is what makes "强 propwash coupling" real; it must be
   carried into layer-11 ground-roll modeling.
3. **Engine-out yaw (C/D) is the critical constraint**: at V_R,
   Cn_EO = T_max·y_eng/(q·S·b) ≈ 0.057-0.058, but max rudder authority is
   ΔCn_r ≈ 0.027 (25° deflection, 25% rudder chord). Control-limited
   **V_MC ≈ 24-27 m/s = 1.45×V_R = 1.7×V_s** — the twins cannot hold
   heading after a V_R engine failure with the current VTP (Vv=0.035).
   Options for the next iteration: VTP volume → ~0.06, larger rudder
   (30% chord / 30°), or FCS asymmetric-thrust/differential-brake assist.
   Layer 9 (JSBSim) will verify dynamically.

## Layer 5 — VSPAERO aero database + stability derivatives (tinybox, 2026-10-04)

`scripts/vspaero_sweep.py`: per prototype, sets AnalysisRef (Sref=S, cref=MAC,
bref=b, Xcg from a component-CG build on the layer-1 weight model), runs
α ∈ -4..16° at β=0 (11 pts) + β ∈ 0..8° at α=4 (3 pts), Mach 0.05,
WakeNumIter 10. Derivatives from linear fits; SM = -C_mα/C_Lα.
Artifacts: `results/vspaero_db.csv`, `vspaero_derivs.csv`, `vspaero_*.png`.

Final baseline (Vh=0.65 raised tails, vertical VTP, pods omitted in sweep):

| proto | C_Lα 1/rad | C_mα 1/rad | SM | C_nβ 1/rad | C_lβ 1/rad | checks |
|---|---|---|---|---|---|---|
| A | 6.22 | −0.11 | 1.7% | +0.009 | −0.044 | C_mα<0 ✓ C_nβ>0 ✓ SM low |
| B | 6.10 | −0.38 | 6.1% | +0.008 | −0.040 | all ✓ |
| C | 6.59 | −0.17 | 2.6% | +0.014 | −0.045 | C_mα<0 ✓ C_nβ>0 ✓ SM low |
| D | 6.23 | −0.37 | 6.0% | +0.007 | −0.035 | all ✓ |

Findings:

1. **All four are statically stable in both axes** (C_mα<0, C_nβ>0, C_lβ<0 —
   dihedral effect present). This held across every geometry variant tested.
2. **SM: B and D sit in the 5-15% target; A and C trend near-neutral**
   (1.7-2.6%). A/C need either a forward CG shift (fuel/payload bay forward)
   or a larger HTP — next sizing iteration.
3. **Thin-sheet VLM sensitivity (documented limits)**: the absolute SM spread
   shifted with geometry-modeling choices (root-stub length, tail z-offset
   vs the body sheet, wing-pod thin sheets). Twin wing pods poison the VLM
   (sheet-wake overlap) and are omitted in the sweep (pod drag added
   analytically later; thick-body treatment in layer 6). Absolute values are
   ±5% MAC credible; trends and pass/fail are robust.
4. **Caveats for the aero DB**: airfoils are still the symmetric default
   (CL0=0, Cm0=0 — XFOIL-verified naca2412/sd7062 import is next); VLM has
   no viscous drag (L/D values are induced-only); fuselage/body as thin
   sheets. These are fine for derivatives, not for absolute L/D.
5. VSPAERO 3.53 deprecates the separate DegenGeom analysis; the geometry
   export is now `VSPAEROComputeGeometry` only.

## Layer 4 — parametric OpenVSP geometry (tinybox, 2026-10-04)

`scripts/vsp_geometry.py` builds all four digital prototypes headlessly from
the layer-1 sizing configs: fuselage (ellipse, d_fus), wing (S/AR/λ → b, c_r,
c_t), HTP/VTP (tail-volume areas), engine nacelles (central nose pod for
A/B; wing-root pods at y=±0.75 m for C/D) and propeller disks.
`VSP_NO_PROPS=1` swaps disks for spinner pods (solver smoke / layer 6).

Saved `vsp_models/{A,B,C,D}.vsp3`, verified against sizing (all exact):

| proto | b target/VSP | S target/VSP | S_h | S_v |
|---|---|---|---|---|
| A | 6.197 / 6.197 | 3.200 / 3.200 | 0.675 | 0.604 |
| B | 6.233 / 6.233 | 3.700 / 3.700 | 0.897 | 0.702 |
| C | 6.116 / 6.116 | 3.400 / 3.400 | 0.772 | 0.633 |
| D | 6.317 / 6.317 | 4.200 / 4.200 | 1.140 | 0.807 |

OpenVSP 3.53 Python API notes (cost a few detours):
- Geom types are uppercase strings: `WING`, `FUSELAGE`, `POD`, `PROP`.
- `FindParm(geom, ...)` only resolves parms on the geom's own container
  (e.g. `TotalSpan/WingGeom`); section parms (`Span`, `Root_Chord`,
  `Tip_Chord`, `Sweep`, `Dihedral`) live on the `XSec` containers and are
  reached via `GetXSecSurf` → `GetXSec`; XForm positioning (`X_Location`,
  `Z_Location`) needs the flat `GetGeomParmIDs` list.
- A `WING` drives its planform from the LAST XSec (Span = half-span per
  side, doubled by symmetry); earlier sections are stubs (kept 0.05 m,
  non-degenerate for the VSPAERO coarse grid).
- Analysis names: `VSPAEROComputeGeometry` (export .vspgeom) then
  `VSPAEROSweep` (vector-valued inputs: `SetDoubleAnalysisInput(..., [v])`).

Solver connectivity (VSPAERO 7.2.2, `scripts/vspaero_run.py`): bare wing and
full prototype D both converge (0.05 s and ~9 s per alpha). **Propeller geoms
must be excluded from the static VLM run** (rotor disks stall the GMRES —
they are layer-6 actuator-disk work). Reference values are not set yet, so
coefficients are referenced to VSPAERO's Sref=100 default; layer 5 must set
AnalysisRef (Sref=b·MAC etc.) before building the aero database.

All-four smoke (α 2-8°, Sref=100-relative): CL rises ~linearly with α on all
models; L/D at α8: A 32, B 53, C 55, D 17 — D's L/D peaks at α4 then falls
as pod/nacelle drag grows (expected STOL-layout behavior). CMytot negative
(stable-ish trend) on B/C/D; magnitudes need proper Cref/MAC reference in
layer 5.

## Layer 3 — flap 2-D sweep (XFOIL GDES-FLAP, tinybox, 2026-10-04)

Plain-flap sweep: naca2412 + sd7062 (A/C and B/D prototype airfoils) × flap
chord {20, 25, 30, 35}% × δ {0..40}° × {free, forced x/c=0.05} @ Re 8e5 =
144 runs, 21 s on 24 workers. Flap geometry via XFOIL's built-in `GDES→Flap`
(rotating a cut TE flap; GDES exits on a blank line after `eXec`).
Artifacts: `results/flap_polars.csv`, `results/flap_summary.csv`,
`results/flap_lookup.csv`, `results/flap_<airfoil>.png`.

ΔC_Lmax vs deflection (free transition, monotone):

| δ | 5° | 10° | 15° | 20° | 25° | 30° | 35° | 40° |
|---|---|---|---|---|---|---|---|---|
| naca2412 (25%) | 0.10 | 0.19 | 0.26 | 0.31 | 0.36 | 0.39 | 0.41 | 0.43 |
| sd7062 (25%) | 0.12 | 0.22 | 0.31 | 0.36 | 0.42 | 0.45 | 0.47 | 0.48 |

Three-state lookup table (`results/flap_lookup.csv`):

| airfoil | state | δ | chord | C_Lmax | α_stall | CD@CL0.8 |
|---|---|---|---|---|---|---|
| naca2412 | cruise | 0° | — | 1.47 | 15° | 0.008 |
| naca2412 | takeoff | 10° | 30% | 1.66 | 13° | 0.008 |
| naca2412 | landing | 30° | 30% | 1.87 | 10° | — |
| sd7062 | cruise | 0° | — | 1.69 | 16° | 0.008 |
| sd7062 | takeoff | 10° | 20% | 1.92 | 15° | 0.008 |
| sd7062 | landing | 30° | 25% | 2.15 | 11° | — |

Findings:

1. **Layer-1 C_Lmax placeholders were optimistic** (landing 2.30 vs XFOIL
   1.87 naca2412 / 2.15 sd7062). Fed back into `sizing/configs.py` per
   prototype; A/C stall speeds rose ~1.6 m/s (still ≤ 17 m/s), A s_TO 103→109 m.
2. **Airfoil choice beats flap size**: sd7062 with a 25% flap at 30° reaches
   C_Lmax 2.15 where naca2412 needs 30% flap and still gets 1.87. STOL
   prototype D and heavy-lift B get the sd7062 family; A/C keep naca2412.
3. **XFOIL ΔC_Lmax is conservative** vs empirical plain-flap data (~0.9-1.1
   at 40° per Raymer); expect further gains from a slotted flap. Slot/gap/
   overlap cannot be modeled in XFOIL — options: MSES (multi-element viscous,
   needs `sudo apt install gfortran libx11-dev` + source build) or 2-D
   multi-element RANS (SU2, layer-12 target #4).
4. Cm captured for later trim work (layer 9): cruise −0.05…−0.08, takeoff
   −0.17; landing Cm at approach C_L needs a finer α grid (δ=30° polars start
   above C_L 0.8) — refine when JSBSim trim needs it.

## Layer 2 — airfoil polars (XFOIL 6.99, tinybox, 2026-10-04)

Sweep: 8 airfoils × Re {0.4, 0.6, 0.8, 1.0, 1.5}e6 × {free, forced x/c=0.05,
forced x/c=0.10} = 120 runs (~1 min on 24 workers), alpha -6..16°.
Artifacts: `results/airfoil_polars.csv`, `results/airfoil_summary.csv`,
`results/airfoil_<name>.png`.

Free-transition (clean composite surface), best per airfoil:

| airfoil | Re | L/D_max | C_L at L/D_max | C_Lmax | cruise-plateau C_L |
|---|---|---|---|---|---|
| naca2412 | 1.5e6 | 108 | 0.81 | 1.61 | 0.57–0.90 |
| ag36 | 1.5e6 | 110 | 0.73 | 1.47 | 0.62–0.84 |
| sd7037 | 1.5e6 | 123 | 0.91 | 1.53 | ~0.59–0.90 |
| naca4412 | 1.5e6 | 144 | 1.03 | 1.70 | none (bucket >0.9) |
| e197 | 1.5e6 | 138 | 0.99 | 1.33 | ~0.75–0.90 |
| sd7062 | 1.5e6 | 133 | 1.23 | **1.77** | none (bucket >1.0) |

Findings that feed layers 3-5:

1. **Roughness penalty is ~-20 to -30% L/D**: forced transition (x/c=0.05-0.10,
   surrogate for seams/rain/bug strikes) drops L/D_max consistently: naca2412
   -25%, ag36 -22%, sd7062 -29%. Expect ~25% range loss on a dirty airframe;
   layer-12 CFD + layer-7 engine map decide how much of this to design for.
2. **Cruise C_L target mismatch**: airfoil drag buckets sit at C_L 0.73-0.81,
   but layer-1 sized cruise at C_L* = 0.89-0.97. Iterate sizing: increase S
   (~10-15%) or reduce C_fe to land cruise at C_L ≈ 0.8. Do this before
   layer-4 geometry is frozen.
3. **C_Lmax scale**: XFOIL C_Lmax_clean = 1.35 (naca2412) to 1.77 (sd7062);
   layer-1 placeholder 1.35 was conservative. STOL prototype D should take
   sd7062 (high C_Lmax, bucket at C_L 1.1-1.2); endurance prototype A should
   take naca2412 or ag36 (bucket at 0.73-0.81).
4. **XFOIL caveat**: SD7037 forced-transition runs oscillate at Re 6e5-1e6
   (laminar-bubble convergence, 45 s timeout + VACC retry did not fully fix);
   its free-transition trend is credible, forced-mode gaps are flagged in the
   CSV. XFLR5/MSES are fallbacks if SD7037 stays shortlisted.

## Layer 1 — parametric sizing (`sizing/`)

Inputs per candidate: `S, AR, λ, t/c, C_Lmax, P_engine, η_p, D_prop, fuel, payload`
plus mission spec. Outputs: span, chord, wing loading, stall speeds,
best-L/D cruise point, cruise power + throttle fraction, climb rate, takeoff
distance (ground roll + air to 15 m), endurance, range, and constraint gates.

Formulas:

```
V_s      = sqrt(2W/(ρ S C_Lmax))
C_D      = C_D0 + C_L²/(π e AR),   C_D0 = C_fe·S_wet/S
P_shaft  = D·V/η_p
ROC      = (η_p·P_av − D·V)/W          (excess-power sweep)
m·dV/dt = T − D − μ(W−L)               (ground roll, integrated in V-space)
T0       = (2ρA_d(η0·P)²)^(1/3)        (actuator-disk static thrust)
R,E      = closed-form max-L/D tracking profile, 10% reserve
```

Key assumptions (all flagged as placeholders — later layers replace them):

| Item | Value | Replaced by |
|---|---|---|
| C_fe (flat-plate skin friction) | 0.009 | layer 2 polars / layer 12 CFD |
| Oswald e | 0.80 | layer 5 VSPAERO |
| C_Lmax clean / TO / land | 1.35 / 1.80 / 2.30 | layers 2–3 |
| BSFC | 400 g/kWh + part-throttle penalty `eff = bsfc·(0.55 + 0.45/frac)` | layer 7 engine map |
| η_p cruise / static | 0.75 / 0.55 | layer 8 BEMT |
| empty weight | component buildup, ~0.50–0.60 MTOW | layer 13 FEA |
| μ_roll | 0.05 | layer 11 ground model |

Constraints (defaults): V_s,land ≤ 17 m/s; ROC ≥ 4 m/s; fuel fraction 3–33% MTOW;
s_TO ≤ 250 m. Same payload (12 kg) for all candidates; reserve 10%.

## Layer-1 results (tinybox, 2026-10-04)

Scan: S = 2.5–4.5 m² × AR = 8–14 × MTOW = 70–140 kg × 4 engine packs
= 16,380 candidates, **7,498 feasible (46%)** in 6 s.

| pack | feasible | max range | max endurance | min s_TO | max ROC |
|---|---|---|---|---|---|
| single 150cc (11 kW) | 43% | 1056 km | 14.8 h | 47 m | 11.0 m/s |
| single 200cc (16.5 kW) | 45% | 809 km | 10.9 h | 29 m | 16.9 m/s |
| twin 2×150cc (22 kW) | 48% | 672 km | 8.7 h | 17 m | 21.3 m/s |

Four digital prototypes (all feasible):

| | A single150 endurance | B single200 lift | C twin150 cruise | D twin150 STOL |
|---|---|---|---|---|
| MTOW / fuel | 95 / 27.2 kg | 112 / 37.0 kg | 105 / 28.4 kg | 115 / 31.1 kg |
| S / AR | 3.2 m² / 12 | 3.7 m² / 10.5 | 3.4 m² / 11 | 4.2 m² / 9.5 |
| V_s,land | 14.4 m/s | 14.5 m/s | 14.7 m/s | **13.8 m/s** |
| L/D_max | **15.5** | 14.6 | 14.4 | 13.5 |
| Cruise throttle | 16% | 14% | 10% | 11% |
| ROC | 7.5 m/s | 9.8 m/s | 14.5 m/s | 13.1 m/s |
| s_TO | 103 m | 84 m | 53 m | **56 m** |
| Endurance | **10.8 h** | 10.0 h | 5.9 h | 6.4 h |
| Range | **797 km** | 764 km | 451 km | 477 km |

### Findings that matter

1. **Cruise throttle fraction is a first-order design driver.** 150cc-class
   engines cruise at 10–16% throttle; carbureted 2-stroke SFC degrades ~2–3×
   there. It (not L/D) now separates singles from twins and moves the range
   optimum to low AR. The layer-7 engine map decides this trade — do not
   finalize engine size before it exists.
2. **Twins buy STOL + safety, not cruise economy**: ~2× s_TO advantage
   (17–56 m vs 84–103 m), ROC 13–14.5 m/s, engine-out capability; but ~45%
   less range/endurance at equal payload. Fair comparison needs the same
   payload/fuel/cruise speed — satisfied here.
3. **Cruise C_L* = 0.89–0.97**, at/above the layer-2 airfoil plateau target
   (0.4–0.9). Airfoil selection must keep its drag bucket at C_L ≈ 0.9+, or
   S/C_D0 must change.
4. **B hits the 33% fuel-fraction cap** (0.33) — 200cc single needs the tank
   sizing question answered before the build.
5. Range/endurance are **upper bounds**: part-throttle SFC model is crude and
   η_p = 0.75 optimistic; layers 7–8 will correct both.

Artifacts: `results/scan_all.csv`, `results/scan_<pack>.csv`,
`results/prototypes_table.csv`, `results/scan_pareto.png`,
`results/scan_feasibility_mtow100.png`, `results/prototypes_compare.png`.

## Run

```bash
ssh tinybox
cd ~/Documents/plane
.venv/bin/python scripts/scan.py        # design-space sweep (~6 s)
.venv/bin/python scripts/prototypes.py  # A-D comparison
```

## Layer 12 - prop-wing-flap CFD closure (`su2/`, OpenFOAM pending)

Goal: validate the blown-flap gain the twins rely on (D: C_Lmax_eff 2.78,
C: 2.27) and feed a correction `k_pw` back into Layer 6/11.

### Gates closed so far

| gate | what it proves | result |
|---|---|---|
| 1 force/REF_AREA regression | frame, markers, half-model reference area | **PASS** - RefForce/q_inf = 2.1440 m2 vs modelled 2.1440 (+0.000%) |
| 3a disk input integration | BEMT -> file -> Renard normalisation -> BC | **PASS** - int dCT/drbar matches BEMT to 1.8e-16 |
| 2b Euler discriminator | geometry, frame, topology (viscosity off) | **CLEARED** - Euler C_L(0) = +0.349, 0.1006/deg |
| 2 clean-wing viscous (OF) | wall-resolved 3-D finite wing | **PASS** - C_L0 0.3499, C_Lalpha 0.0804/deg |
| 3b disk thrust readout | CFD thrust closure | **DEFERRED** - no validated readout |
| 5 coupled k_pw | the actual calibration | **not started** |

### OpenFOAM stage 2 - 3-D finite wing (`foam/`, 2026-10-05)

Body-fitted extruded O-grid (`foam/ofcase.py::build_wing_msh`), symmetry plane
at the root, airfoil cap of the top face as the tip wall; 427k elements, 20
boundary layers, first cell 8 um (y+ < 1). Free-stream rotated per angle so one
mesh serves the whole sweep. Re = 6.5e5 at U = 14.07 m/s.

| alpha | C_L | segment slope | oracle |
|---|---|---|---|
| -4 | +0.0393 | - | - |
| -2 | +0.1912 | 0.0760/deg | 0.0835-0.0881 |
| 0 | **+0.3457** | 0.0773/deg | C_L0 0.330-0.349 |
| +2 | +0.5127 | **0.0835/deg** | 0.0835-0.0881 |
| +4 | +0.6808 | 0.0841/deg | (approaching stall) |
| +6 | +0.8368 | 0.0780/deg | (stalled) |

Least-squares over the linear range (|alpha| <= 2):
**C_Lalpha = 0.08037/deg (-3.7% vs oracle lower bound), C_L0 = +0.3499 (+0.3%)**
- both inside the 5% gate. Spanwise loading is root-loaded and unloads
monotonically to the tip (`results/of_w3b_g20_span.csv`).

### OpenFOAM stage 3a/3b - actuator disk as a volumetric source (2026-10-07)

User-directed formulation, and it is the right one:

    disk = thin cellZone + volumetric momentum source      (this)
    disk = internal physical surface                       (rejected)

A real actuator disk is a zero-thickness momentum jump. Modelling it as a wall
or baffle would inject no-slip, blockage, separation and artificial pressure
drag that have no place in the k_pw calibration. The disk is a `topoSet
cylinderToCell` **cell selector** (99 393 cells, V = 0.003236 m3,
thickness = 0.01467 m = 0.0277 D); no cylinder geometry enters the mesher. This
also sidesteps SU2's coincident front/back face requirement entirely.

`actuationDiskSource`/`radialActuationDiskSource` both ship in
`libfvOptions.so`, and `fvOptions` loads and runs in this build.

| quantity | value | target | verdict |
|---|---|---|---|
| disk cellZone | 100 098 cells, 0.0277 D thick | 0.02-0.05 D | ok |
| T_delivered (far-wake momentum) | 16.0495 N | BEMT 16.1806 N | **-0.81 %** |
| T_source (fvOptions formula) | 18.7303 N | - | **+14.3 % vs delivered** |
| momentum excess / pressure | +15.745 / +0.304 N | - | pressure recovered |
| rmax convergence | 1.4R: 16.79 N (max) | plateau | converged in radius |

**`T_delivered/T_source = 0.857` is constant** across three independent
calibrations (13.52/15.77, 16.53/19.29, 16.05/18.73). The source's own Froude
formula therefore over-predicts the momentum it actually injects by a fixed
14.3 %: `T = 2 rho A U_up^2 a (1-a)` is the TURBINE form and its `U_up^2`
premise does not hold for a propulsor. **The source formula is not a usable
thrust readout; the independent far-wake momentum integral is.** Calibrate
empirically against the latter.

The refined near box must also be **wider than the slipstream**: at x/D = 1 the
wake radius is ~1.4 R, so the original 1.4 R box left the outer wake in 0.25 m
coarse cells and the momentum integral came out short. It is now 2.6 R
(591 888 cells vs 380 821). Sampling planes must additionally stay inside the
refined region, and `p` is kinematic (needs x rho).

**Disk type: THRUST-MATCHED, axial only - and that is sufficient.**

`actuationDiskSource` applies no swirl, so `u_theta = 0` by construction and
`P_fluid = int S.U dV` is not shaft power (the implied eta 0.882 exceeds BEMT
0.727 by the missing swirl loss). This was sized before deciding it mattered,
via `foam/swirl_check.py`, which recovers the BEMT tangential induction from
its own torque distribution (`dQ/dr = 4 pi rho r^3 U (1-a) a' Omega`):

| r/R | u_theta [m/s] | q_swirl / q_inf |
|---|---|---|
| 0.30 | 0.852 | 0.37 % |
| 0.53 | 0.856 | 0.37 % |
| 0.76 | 0.762 | 0.29 % |
| 0.99 | 0.464 | 0.11 % |

**CORRECTION (found while starting OF4).** The first version of this check ran
at J = 0.509, a LIGHT loading (Ct_std = 0.60), where swirl is ~0.3 % of q_inf
and can indeed be ignored. But the k_pw design point is **takeoff**, and
prototype D's takeoff state is J_TO = 0.1998 at 5990 rpm (from
results/prop_operating.csv), i.e. Ct_std = **7.21** - a heavily loaded disk.
Swirl scales steeply with loading:

| point | J | n [rps] | T [N] | Ct_std | u_theta [m/s] | q_swirl/q_inf | swirl angle |
|---|---|---|---|---|---|---|---|
| J_TO (design point) | 0.200 | 132.9 | 192.9 | **7.21** | 7.31 | **27.5 %** | 13.6 deg |
| J_climb | 0.213 | 124.8 | 170.6 | 6.38 | 6.46 | 21.5 % | 12.7 deg |
| J_cruise | 0.416 | 63.9 | 35.1 | 1.31 | 1.61 | 1.3 % | 5.1 deg |

So at the point k_pw is meant to calibrate, swirl is **27.5 % of q_inf** and
the flow angle is ~14 deg - comparable to the flap deflection itself. It is a
FIRST-ORDER effect, not a ~1 % correction. The axial-only
`actuationDiskSource` is adequate only for light loading.

A second, harder constraint follows: **Froude cannot represent takeoff at
all.** Its thrust is `2 rho A U^2 a (1-a)`, which peaks at a = 0.5 at
26.8 N - but the takeoff operating point needs 192.9 N. Heavy loading requires
the `variableScaling` variant, whose `Ct` acts directly.

Consequence: OF4/OF5 must either (a) carry swirl (`rotorDiskSource`, which is
cellZone-based with a documented `blade` profile file), or (b) run k_pw at a
loading where swirl is genuinely negligible and state that restriction. This is
the open decision before OF4.

**Closing the thrust took three corrections, all of them arithmetic:**
1. OpenFOAM's `p` is **kinematic** (m2/s2); any pressure term needs x rho.
2. `Froude` computes `T = 2 rho A U_up^2 a (1-a)` with `a = 1 - Cp/Ct`, and only
   the RATIO matters. The propulsor coefficients (Cp/Ct = 0.6999 -> a = 0.3001)
   and the freestream-based momentum induction (a = 0.1857) both overshot;
   solving `a(1-a) = T/(2 rho A U_up^2)` against the **measured** U_up gives
   a = 0.1555 (Cp/Ct = 0.8445) and T_source to 1.9 %.
3. The `upstreamPoint` monitor sits inside the disk's own induction region and
   reads U_up = 14.96 m/s, not the 14.07 m/s freestream, so calibrating against
   U_inf left a +15 % error - exactly the U_up^2 error.
Sampling planes must also stay inside the refined zone: the near box is only
+/-1.4 R wide, and `rmax = 1.6 R` was pulling in 0.25 m far cells whose volume
dwarfed the entire slipstream.

### Stage 3b wake survey - the Layer-6 ROM check, before any wing
`foam/wake.py` / `foam/disk_report.py` -> `results/of_disk3_disk_d1.csv`.

| x/D | R_s/R | Ux/Uinf (core) | q_slip/q_inf |
|---|---|---|---|
| 0.10 | 1.458 | 1.1817 | 1.398 |
| 0.25 | 1.601 | 1.2061 | 1.457 |
| 0.50 | 1.663 | 1.2230 | 1.498 |
| 1.00 | 1.136 | 1.2316 | 1.520 |

Momentum theory for this thrust: U_disk/Uinf = 1.1334, U_far/Uinf = 1.2668.
Measured axial development: 1.154 (0.1D) -> 1.200 (1.0D) -> 1.195 (2.6D), so
the wake is fully developed by x/D = 1 and the far-field target is reached to
within 5%. As an independent integrated check, the domain-integrated actuator
force `int_Omega S dV` over the whole domain gives 16.2 N against BEMT's
16.18 N (+0.1%).

The profile is flat to r/R = 0.75 then falls sharply, i.e. a well-defined
slipstream edge rather than a smeared one. Layer 6/11's `q_slip/q_inf`
assumption is therefore of the right order and shape.

### Why the snappyHexMesh route was abandoned (retired, do not retry)

**Status: removed.** `foam/{snappy,wsnappy,wing3d,wing_stl,layer_diag,wing3d_alpha}.py`
and every `w3d_*` case directory are deleted; the shared post-processing they
held (`span_loading`, `alpha_result`, `latest_vtp`, `chord_at`) moved to
`foam/wing_report.py`. The live 3-D modules are `ofcase` (meshing + case
writer), `wing3d_body` (finite wing), `wing_report`/`span_report`
(post-processing), `disk3d`/`disk_report`/`wake`/`swirl_check` (actuator disk),
`forces`, `wing2d` (2-D regression) and `smoke`.

The first stage-2 attempt used snappyHexMesh. It meshed and solved (`rc = 0`)
but every mesh level came out with **zero boundary layers**: snappy builds the
prism stack, reports it in the `patch faces / layers` table, then silently
rolls the whole stack back. All of these were tried and all failed:
22 layers x 1.18, 12 x 1.25, 8/5/3 layers at 20-100 um, relaxed
`meshQualityControls`, `minDeterminant 0`, `maxThicknessToMedialRatio` 0.1-0.5,
`nGrow 0`, uniform surface levels, `relativeSizes true/false`.
The layer-phase log shows 99% extrusion followed by
`Detected 417002 illegal faces` and `Added 0 out of 161024 cells`; the dominant
gate is `faces on cells with determinant < 0.001`. The sd7062 trailing edge is
0.47 mm thick at the root - thinner than any usable first cell - and the
contrast between the 0.4 m background cell and the 20 um near-wall cell is
6000:1, which no snappy setting survived.

Consequence for the oracle, and the reason this matters: **a no-boundary-layer
solve converges to the inviscid answer, not the viscous one.** Snappy gave
C_L0 = 0.3505 (M1) rising to 0.3442 (M2, fit), i.e. drifting toward the Euler
value 0.349 rather than the viscous 0.330-0.349, and over-predicted a strongly
negative-alpha C_L by ~2x (0.0838 vs 0.0393 at alpha = -4). Viscous decambering
is what lowers C_L0 and suppresses negative-alpha lift; without wall shear it
is simply absent.

gmsh's BoundaryLayer field avoids all of it because it is applied to the 2-D
O-grid before extrusion, so the prism stack is swept along a well-shaped
surface mesh. Use it for every 3-D case from here on.

### OpenFOAM environment notes (v1912, this package is INCOMPLETE)
`/usr/share/openfoam` contains only `etc/`; FOAM_ETC and FOAM_LIBBIN are empty
and `/usr/share/openfoam/bin` is absent. Every function object (`forces`,
`yPlus`) falls back to dynamicCode and dies writing a sha1 stream, so forces
are integrated from `foamToVTK` output by `foam/forces.py` instead.

Gotchas already paid for:
- patch TYPE must match the field BC: `symmetryPlane` in a field on a
  `patch`-typed mesh patch aborts with "Attempt to cast type patch to type
  symmetryPlane".
- blockMesh face labelling is easy to get backwards. In
  `foam/snappy.py::block_mesh_dict` the y = 0 half-wing plane is patch
  `bottom`, NOT `sideWall`. Assigning `freestream` there leaves the wing root
  open to the flow: root loading collapses and C_D_p comes out ~14x high
  (0.0318 vs 0.0210 recovered).
- `surfaceFeatureExtract` reads from `constant/triSurface/`, not `constant/`.
- `addLayersControls` accepts exactly TWO of {firstLayerThickness,
  finalLayerThickness, thickness, expansionRatio}; three is a fatal error.
- solver termination: `residualControl` on p is what stops the run. A flat
  `endTime` cap of several thousand iterations is hours of nothing.
- `foam/wing3d.py` used to shadow `su2/wing3d.py` on sys.path; `ofcase.py` now
  loads the SU2 module by explicit path.

### Layer-8 bug found by this work
`sizing/prop_bemt.py::twist` had `atan(P_over_D/(2*pi*mu))`. Substituting
`r = mu*D/2` gives `atan(P_over_D/(pi*mu))`, so the geometric pitch angle was
half-size and zero thrust landed at J = 0.5 instead of J = P/D. Fixed; cruise
eta_p is now 0.62-0.73 for the twins and D range fell 657 -> 472 km.

### k_pw is a sensitivity band, NOT a calibrated value
No converged wall-resolved 3-D viscous case exists, so `K_PW` stays an
instance-level parameter (`K_PW_BAND = 0.4, 0.5, 0.7, 1.0`).
`scripts/kpw_sensitivity.py` -> `results/kpw_sensitivity.csv`:

| k_pw | s_TO,D | gain vs best of A/B/C | s_land,D | gain |
|---|---|---|---|---|
| 0.4 | 46.8 m | +8.2 m | 218.1 m | +15.7 m |
| 0.5 | 44.5 m | +7.2 m | 216.2 m | +15.5 m |
| 0.7 | 41.4 m | +5.6 m | 212.8 m | +15.0 m |
| 1.0 | 39.5 m | +2.6 m | 208.6 m | +14.2 m |

D keeps its STOL advantage across the whole band and the advantage is largest
at the conservative end, so k_pw scales the size of D advantage rather than
deciding whether D works. **CFD is a calibration task, not a design gate.**

### SU2 actuator-disk requirements (read from CEulerSolver.cpp)
1. `r_ = |P - CENTER| / AD_Radius` is the full 3-D distance, so the inlet and
   outlet faces must be geometrically COINCIDENT.
2. The station loop indexes `rad_v[iEl-1]`, which is `rad_v[-1]` at
   `iEl == 0`, so every node must satisfy `r_bar > rad_v[0]`. The shipped
   tutorial disk is an annulus over a spinner for this reason.
3. `MARKER_INLET` takes TOTAL conditions.
4. The mesh must be conformal to the disk plane, otherwise node splitting
   corrupts straddling elements.
gmsh cannot satisfy (1) and (4): `fragment()` is always conformal-and-merged
and `occ.copy()` gets merged too. `su2split.py` splits the mesh instead, but
gmsh still merged the two boxes (1537 of 36384 tets straddled the plane).

Consequently 3b runs on SU2's own shipped mesh
(`su2/tutorial/propeller_variable_load.su2`, markers DISK/DISK_BACK exactly
coincident at x = 0): the tutorial reproduces (`Exit Success`) and the Layer-8
radial loads drive it (`rc = 0`, `rms[Rho] = -4.85`).

Readout caveat: the per-marker force breakdown is dominated by a momentum-flux
term and is not a thrust measure, and a node-wise pressure-jump readout does
not currently close (momentum-consistency ratio 0.062). Gate 3b stays open.

### Frames and pitfalls that cost the most time
- SU2's aerodynamic frame for a wing is chord x, **span y, lift z**
  (CL = CFz, CD = CFx) - verify from solver output rather than assuming.
- `occ.fragment` does NOT carve a cavity when a solid sits inside a box; the
  body becomes an internal baffle and flow passes through it (CL ~ 0 with
  negative drag). `occ.cut` works.
- gmsh `occ.extrude` needs `numElements=[n_span]`; without it the span is a
  single cell (CD 0.171 -> 0.095 and dCL/dalpha 0.019 -> 0.052 once fixed).
- `GMSH2SU2` must map gmsh 6 -> SU2 13 (prism); otherwise every near-wall
  prism from a BoundaryLayer extrusion is silently dropped.
- SU2 8.5.0 expects CGNS in ADF, not HDF5, and this build has no `SU2_MSH`.

### Clean-wing oracle (3-D, not 2-D)
A 3-D half-wing must be compared against a finite-wing slope, not the 2-D XFOIL
value: `a = a0/(1 + a0/(pi e AR))` with a0 = 0.109/deg, AR_3D = 9.31,
e = 0.7-0.9 gives **0.0835-0.0881/deg**, and
`C_D(0) ~ C_D0 + C_L0^2/(pi e AR) = 0.0118 + 0.0080 = 0.0198`.
The Layer-2 1.47-1.77 figures are C_Lmax, NOT C_L at alpha = 0.

### Blocker (resolved)
3-D wall-resolved meshing needed OpenFOAM; it is now installed and the 3-D
finite wing closes (see the stage-2 table above). snappyHexMesh output must
still never be converted to SU2 - SU2 native accepts only
tet/hex/prism/pyramid while snappyHexMesh emits split-hex/polyhedral cells -
but that no longer matters because 3-D viscous work stays inside OpenFOAM.

## Next steps

1. Layer 2: XFOIL polar sweeps at Re = 0.4/0.6/0.8/1.0/1.5e6, transition
   (free / forced), for the cruise plateau C_L = 0.4–0.9. Install XFOIL on
   tinybox (apt or source build).
2. Install OpenVSP (conda-forge `openvsp`) + verify Python API on tinybox;
   mirror prototypes A–D as parametric geometry (layer 4).
3. Layer 3: 2-D flap sweep (chord 20–35%, δ 0–40°, slot/gap/overlap) →
   ΔC_Lmax, ΔC_D lookup for the three flap states.
