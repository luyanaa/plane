"""Layer 12: SU2 VARIABLE_LOAD actuator-disk input file writer.

SU2 reads one block per disk (CEulerSolver::ReadActDisk_InputFile). Each block

    MARKER_ACTDISK= <inlet_marker> <outlet_marker>
    CENTER= x y z
    AXIS= ax ay az
    RADIUS= R
    ADV_RATIO= J
    NROW= n
    <one dummy header line, skipped by SU2>
    r_bar dCT/dr_bar dCP/dr_bar dCR/dr_bar     (n lines)

The coefficients use the Renard convention: force ref = rho n^2 D^4,
power ref = rho n^3 D^5, and the radial station is r_bar = r/R.

SU2 converts these to force per unit disk area. Substituting the Renard
definitions and J = V/(n D) into the source expressions collapses both to

    Fa = (dT/dr) / (2 pi r)        (axial)
    Ft = (dQ/dr) / (2 pi r^2)      (tangential, swirl)

i.e. the load depends only on the physical blade-element force gradients, so
the BEMT distributions can be written directly.
"""
from pathlib import Path


def write_actdisk_input(path, disks, header='# rs=r/R  dCT/drs  dCP/drs  dCR/drs'):
    """Write the VARIABLE_LOAD input file.

    disks: iterable of dicts with keys
        markers (name_in, name_out), center (x, y, z), axis (ax, ay, az),
        radius, advance_ratio, r_bar, d_ct, d_cp, d_cr.
    """
    lines = []
    for d in disks:
        names = list(d['markers'])
        ni = names[0]
        no = names[1] if len(names) > 1 else names[0]
        cx, cy, cz = d['center']
        ax, ay, az = d['axis']
        lines.append(f'MARKER_ACTDISK= {ni} {no}')
        lines.append(f'CENTER= {cx:.10g} {cy:.10g} {cz:.10g}')
        lines.append(f'AXIS= {ax:.10g} {ay:.10g} {az:.10g}')
        lines.append(f'RADIUS= {d["radius"]:.10g}')
        lines.append(f'ADV_RATIO= {d["advance_ratio"]:.10g}')
        r_bar = d['r_bar']
        lines.append(f'NROW= {len(r_bar)}')
        lines.append(header)
        for i in range(len(r_bar)):
            lines.append('%.6f %.8e %.8e %.8e' % (
                r_bar[i], d['d_ct'][i], d['d_cp'][i], d['d_cr'][i]))
    Path(path).write_text('\n'.join(lines) + '\n')
    return path


def disk_from_bemt(prop, V, n_rps, markers, center, axis):
    """Build one disk dict from a layer-8 PropBEMT radial-load solution."""
    r_bar, dct, dcp, dcr, j = prop.radial_loads(V, n_rps)
    return dict(markers=markers, center=center, axis=axis,
                radius=prop.R, advance_ratio=j,
                r_bar=r_bar, d_ct=dct, d_cp=dcp, d_cr=dcr)
