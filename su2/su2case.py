"""SU2 case setup and run wrapper (layer 12).

Generates .cfg files for the 2-D airfoil V&V and the 3-D prop-wing-flap
half-model, runs SU2_CFD, and parses the history/surface output.
"""
import subprocess
from pathlib import Path

SU2_BIN = Path.home() / 'opt' / 'su2' / 'bin'

BASE_2D = """
SOLVER= RANS
MESH_FORMAT= SU2
KIND_TURB_MODEL= SST
KIND_TRANS_MODEL= {trans}
MATH_PROBLEM= DIRECT
RESTART_SOL= NO

FLUID_MODEL= STANDARD_AIR
VISCOSITY_MODEL= SUTHERLAND
CONDUCTIVITY_MODEL= CONSTANT_PRANDTL

MACH_NUMBER= {mach}
AOA= {aoa}
FREESTREAM_PRESSURE= {p_inf}
FREESTREAM_TEMPERATURE= {t_inf}
REYNOLDS_NUMBER= {re}
REYNOLDS_LENGTH= 1.0
INIT_OPTION= TD_CONDITIONS

REF_ORIGIN_MOMENT_X= 0.25
REF_ORIGIN_MOMENT_Y= 0.0
REF_ORIGIN_MOMENT_Z= 0.0
REF_LENGTH= 1.0
REF_AREA= 1.0

MARKER_HEATFLUX= ( airfoil, 0.0 )
MARKER_FAR= ( farfield )
MARKER_PLOTTING= ( airfoil )
MARKER_MONITORING= ( airfoil )

NUM_METHOD_GRAD= WEIGHTED_LEAST_SQUARES
CFL_NUMBER= {cfl}
CFL_ADAPT= NO
ITER= {iters}
CONV_NUM_METHOD_FLOW= ROE
ENTROPY_FIX_COEFF= 0.1
MUSCL_FLOW= YES
SLOPE_LIMITER_FLOW= VENKATAKRISHNAN
TIME_DISCRE_FLOW= EULER_IMPLICIT

CONV_NUM_METHOD_TURB= SCALAR_UPWIND
MUSCL_TURB= NO
TIME_DISCRE_TURB= EULER_IMPLICIT
CFL_REDUCTION_TURB= 1.0

LINEAR_SOLVER= FGMRES
LINEAR_SOLVER_PREC= ILU
LINEAR_SOLVER_ERROR= 1E-4
LINEAR_SOLVER_ITER= 10
LINEAR_SOLVER_ILU_FILL_IN= 0

SCREEN_OUTPUT= (INNER_ITER, RMS_DENSITY, RMS_ENERGY, LIFT, DRAG)
HISTORY_OUTPUT= (ITER, RMS_RES, AERO_COEFF)
OUTPUT_FILES= (RESTART_ASCII, PARAVIEW_ASCII, SURFACE_CSV)

CONV_FILENAME= history
VOLUME_FILENAME= flow
SURFACE_FILENAME= surface
RESTART_FILENAME= restart
"""

CONV_TRANS = {
    'none': 'NONE',
    'lm': 'LM',
}

BASE_3D = """
SOLVER= RANS
MESH_FORMAT= SU2
KIND_TURB_MODEL= SST
KIND_TRANS_MODEL= {trans}
MATH_PROBLEM= DIRECT
FLUID_MODEL= STANDARD_AIR
VISCOSITY_MODEL= SUTHERLAND
CONDUCTIVITY_MODEL= CONSTANT_PRANDTL
MACH_NUMBER= {mach}
AOA= {aoa}
SIDESLIP_ANGLE= 0.0
FREESTREAM_PRESSURE= {p_inf}
FREESTREAM_TEMPERATURE= {t_inf}
REYNOLDS_NUMBER= {re}
REYNOLDS_LENGTH= {ref_length}
INIT_OPTION= TD_CONDITIONS
REF_ORIGIN_MOMENT_X= {x_ref}
REF_ORIGIN_MOMENT_Y= 0.0
REF_ORIGIN_MOMENT_Z= 0.0
REF_LENGTH= {ref_length}
REF_AREA= {ref_area}
MARKER_SYM= ( {sym} )
{wall_markers}
MARKER_FAR= ( {far} )
MARKER_PLOTTING= ( {plot} )
MARKER_MONITORING= ( {plot} )
{actdisk_block}
NUM_METHOD_GRAD= WEIGHTED_LEAST_SQUARES
CFL_NUMBER= {cfl}
CFL_ADAPT= YES
CFL_ADAPT_PARAM= ( 0.1, 2.0, 1.0, 100.0 )
ITER= {iters}
CONV_NUM_METHOD_FLOW= ROE
ENTROPY_FIX_COEFF= 0.1
MUSCL_FLOW= YES
SLOPE_LIMITER_FLOW= VENKATAKRISHNAN
TIME_DISCRE_FLOW= EULER_IMPLICIT
CONV_NUM_METHOD_TURB= SCALAR_UPWIND
MUSCL_TURB= NO
TIME_DISCRE_TURB= EULER_IMPLICIT
CFL_REDUCTION_TURB= 1.0
LINEAR_SOLVER= FGMRES
LINEAR_SOLVER_PREC= ILU
LINEAR_SOLVER_ERROR= 1E-4
LINEAR_SOLVER_ITER= 10
LINEAR_SOLVER_ILU_FILL_IN= 0
SCREEN_OUTPUT= (INNER_ITER, RMS_DENSITY, RMS_ENERGY, LIFT, DRAG)
HISTORY_OUTPUT= (ITER, RMS_RES, AERO_COEFF)
OUTPUT_FILES= (RESTART_ASCII, PARAVIEW_ASCII, SURFACE_CSV, SURFACE_PARAVIEW_ASCII)
CONV_FILENAME= history
VOLUME_FILENAME= flow
SURFACE_FILENAME= surface
RESTART_FILENAME= restart
"""


def write_cfg_3d(path, mesh, walls, sym='symmetry', far='farfield',
                 plot=('wing',), actdisk=None, mach=0.1, aoa=4.0, re=1.0e6,
                 p_inf=101325.0, t_inf=288.15, iters=4000, cfl=3.0,
                 trans='none', ref_length=1.0, ref_area=1.0, x_ref=0.25,
                 actdisk_file='actdiskinput.dat', extra_walls=()):
    """3-D half-model config. `walls` = iterable of wall marker names.
    `actdisk` = (inlet_marker, outlet_marker) or None for prop-off.
    `extra_walls` are additional viscous walls (e.g. the tip cap of a
    half-wing model, which the image plane does not cover)."""
    walls = tuple(walls) + tuple(extra_walls)
    wall_block = 'MARKER_HEATFLUX= ( ' + ', '.join(
        f'{w}, 0.0' for w in walls) + ' )'
    if actdisk:
        double = len(actdisk) > 1
        entries = list(actdisk) + ['0.0'] * (8 - len(actdisk))
        ad = (f'ACTDISK_DOUBLE_SURFACE= {"YES" if double else "NO"}\n'
              'ACTDISK_TYPE= VARIABLE_LOAD\n'
              f'ACTDISK_FILENAME= {actdisk_file}\n'
              'MARKER_ACTDISK= ( ' + ', '.join(entries) + ' )')
    else:
        ad = ''
    Path(path).write_text(f'MESH_FILENAME= {mesh}\n' + BASE_3D.format(
        trans=CONV_TRANS[trans], mach=mach, aoa=aoa, p_inf=p_inf, t_inf=t_inf,
        re=re, cfl=cfl, iters=iters, sym=sym, far=far,
        plot=', '.join(plot), actdisk_block=ad, wall_markers=wall_block,
        ref_length=ref_length, ref_area=ref_area, x_ref=x_ref))
    return path


def write_cfg_2d(path, mesh, mach=0.1, aoa=4.0, re=1.0e6, p_inf=101325.0,
                 t_inf=288.15, iters=3000, cfl=5.0, trans='none'):
    Path(path).write_text(BASE_2D.format(
        trans=CONV_TRANS[trans], mach=mach, aoa=aoa, p_inf=p_inf,
        t_inf=t_inf, re=re, cfl=cfl, iters=iters,
        mesh=mesh).replace('{mesh}', str(mesh)))
    # inject the mesh filename
    s = Path(path).read_text()
    s = f'MESH_FILENAME= {mesh}\n' + s
    Path(path).write_text(s)
    return path


def run_case(cfg, wdir, np_threads=8, timeout=3600):
    wdir = Path(wdir)
    wdir.mkdir(parents=True, exist_ok=True)
    env = {'OMP_NUM_THREADS': str(np_threads)}
    import os
    e = {**os.environ, **env}
    r = subprocess.run([str(SU2_BIN / 'SU2_CFD'), str(cfg)], cwd=str(wdir),
                       capture_output=True, text=True, timeout=timeout, env=e)
    return r


def parse_history(wdir, name='history.csv'):
    """Return (iters, cl, cd, cm) from the SU2 history file."""
    import csv
    p = Path(wdir) / name
    if not p.exists():
        p = Path(wdir) / 'history.csv'
        if not p.exists():
            return None
    rows = []
    with open(p) as fh:
        rd = csv.DictReader(fh)
        for row in rd:
            rows.append(row)
    return rows
