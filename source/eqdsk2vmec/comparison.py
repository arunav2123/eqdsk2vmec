# Written by Arunav Kumar, MIT Plasma Science and Fusion center, 10th May, 2026
"""Compare GEQDSK data with a parsed, written VMEC input."""
from pathlib import Path
import json
import f90nml
import numpy as np
from scipy.spatial import cKDTree


def reconstruct_boundary(namelist, theta):
    """Evaluate retained m modes from f90nml's Fortran-indexed boundary arrays."""
    if namelist['ntor'] != 0:
        raise ValueError('This comparison supports axisymmetric inputs only')
    r, z = np.zeros_like(theta), np.zeros_like(theta)
    for name, target, trig in [('rbc', r, np.cos), ('rbs', r, np.sin),
                                ('zbc', z, np.cos), ('zbs', z, np.sin)]:
        if name not in namelist or (name in ('rbs', 'zbc') and not namelist['lasym']):
            continue
        n_start, m_start = namelist.start_index[name]
        for offset, row in enumerate(namelist[name]):
            m = offset + m_start
            if row is None or m >= namelist['mpol']:
                continue
            n_index = -n_start
            if 0 <= n_index < len(row) and row[n_index] is not None:
                target += row[n_index] * trig(m*theta)
    return r, z


def plot_comparison(gdata, input_path, converted, *, plot_path=None, show=False):
    """Plot serialized geometry and profile knots; save PNG and numerical JSON metrics.

    Profile lines connect serialized spline knots, not a VMEC-computed solution.
    The geometric distance is to a densely sampled curve (8192 points).
    """
    import matplotlib.pyplot as plt
    nml = f90nml.read(input_path)['indata']
    theta = np.linspace(0, 2*np.pi, 8192, endpoint=False)
    r, z = reconstruct_boundary(nml, theta)
    source = np.column_stack([gdata['xbndry'], gdata['zbndry']])
    distances = cKDTree(np.column_stack([r, z])).query(source)[0]
    s_p, pressure = np.array(nml['am_aux_s']), np.array(nml['am_aux_f'])
    s_q, iota = np.array(nml['ai_aux_s']), np.array(nml['ai_aux_f'])
    indices = converted.profile_indices
    if len(s_p) != len(indices) or len(s_q) != len(indices):
        raise ValueError('Serialized profile length differs from selected source knots')
    q_recovered = converted.q_sign / iota
    metrics = dict(boundary_rms_distance_m=float(np.sqrt(np.mean(distances**2))),
                   boundary_max_distance_m=float(np.max(distances)),
                   pressure_max_knot_error_Pa=float(np.max(np.abs(pressure-gdata['sp'][indices]))),
                   q_max_knot_error=float(np.max(np.abs(q_recovered-gdata['qpsi'][indices]))),
                   flux_coordinate_max_knot_error=float(np.max(np.abs(s_p-converted.tflux[indices]))),
                   profile_knot_count=len(indices),
                   pressure_linear_diagnostic_max_error_Pa=float(np.max(np.abs(np.interp(converted.tflux, s_p, pressure)-gdata['sp']))),
                   q_linear_diagnostic_max_error=float(np.max(np.abs(np.interp(converted.tflux, s_q, q_recovered)-gdata['qpsi']))),
                   phiedge_Wb=float(nml['phiedge']), lasym=bool(nml['lasym']),
                   mpol=int(nml['mpol']), cocos_assumed=int(converted.cocos),
                   symmetry_removed_max_m=converted.symmetry_removed_max_m)
    fig = plt.figure(figsize=(12, 8), layout='constrained')
    grid = fig.add_gridspec(3, 2, width_ratios=[1, 1.25])
    ax = fig.add_subplot(grid[:, 0])
    ax.contour(gdata['xgrid'], gdata['zgrid'], gdata['psixz'].T, 24,
               colors='0.82', linewidths=.6)
    ax.plot(gdata['xlim'], gdata['zlim'], color='0.5', linewidth=1, label='Limiter')
    ax.plot(*np.vstack([source, source[0]]).T, color='#1764ab', linewidth=2, label='GEQDSK LCFS')
    ax.plot(np.r_[r,r[0]], np.r_[z,z[0]], '--', color='#e07018', linewidth=1.6, label='input boundary')
    ax.plot(gdata['xaxis'], gdata['zaxis'], '+', color='#1764ab', markersize=10, label='GEQDSK axis')
    ax.set(xlabel='R [m]', ylabel='Z [m]')
    ax.legend(fontsize=8)
    for row, original, sx, values, title in [
        (0, gdata['sp'], s_p, pressure, 'Pressure [Pa]'),
        (1, gdata['qpsi'], s_q, q_recovered, ' q')]:
        ap = fig.add_subplot(grid[row, 1])
        ap.plot(converted.tflux, original, color='#1764ab', linewidth=2, label='GEQDSK samples')
        ap.plot(sx, values, '--', color='#e07018', marker='.', markersize=3,
                linewidth=1, label=' spline knots')
        ap.set(xlabel='Normalized toroidal flux s', ylabel=title)
        ap.grid(alpha=.2)
        ap.legend(fontsize=8)
    ae = fig.add_subplot(grid[2, 1])
    ae.plot(np.arange(len(distances)), distances*1000, color='#8d3675')
    ae.set(xlabel='GEQDSK boundary point index', ylabel='Distance to written boundary [mm]',
           title=f'RMS {metrics["boundary_rms_distance_m"]*1000:.2f} mm; max {max(distances)*1000:.2f} mm')
    ae.grid(alpha=.2)
    fig.suptitle(f'{Path(input_path).name} | LASYM={nml["lasym"]}, MPOL={nml["mpol"]}, COCOS={converted.cocos}\n')
    if plot_path is not None:
        plot_path = Path(plot_path)
        fig.savefig(plot_path, dpi=170)
        plot_path.with_suffix('.json').write_text(json.dumps(metrics, indent=2)+'\n')
    if show:
        plt.show()
    plt.close(fig)
    return metrics
