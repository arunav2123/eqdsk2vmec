# Written by Arunav Kumar, MIT Plasma Science and Fusion center, 10th May, 2026
"""Numerical conversion of axisymmetric GEQDSK equilibria (no VMEC solve)."""
from types import SimpleNamespace
import warnings
import numpy as np
from scipy.interpolate import PchipInterpolator


def endpoint_polynomial(s, values, degree=9):
    """Least-squares polynomial with both endpoints imposed exactly, ascending order."""
    s, values = np.asarray(s), np.asarray(values)
    degree = min(degree, len(s) - 1)
    base = values[0] * (1 - s) + values[-1] * s
    if degree < 2:
        return np.array([values[0], values[-1] - values[0]])
    # p(s) = linear endpoint interpolant + s(1-s) h(s).
    design = s[:, None] * (1 - s[:, None]) * np.polynomial.polynomial.polyvander(s, degree - 2)
    h = np.linalg.lstsq(design, values - base, rcond=None)[0]
    coeff = np.zeros(degree + 1)
    coeff[:2] = [values[0], values[-1] - values[0]]
    coeff[1:-1] += h
    coeff[2:] -= h
    return coeff


def eqdisk2vmec_inputfile(filename, data, *, mpol=12, lasym=False, cocos=None):
    """Return boundary modes and profiles. MPOL retains m=0,...,MPOL-1.

    COCOS 1..8 and 11..18 are supported for signed q. With cocos=None,
    assume COCOS 5 (legacy reciprocal-q convention), emitting a warning.
    This function performs no plotting and never runs VMEC.
    """
    if not isinstance(mpol, (int, np.integer)) or mpol < 2:
        raise ValueError('mpol must be an integer >= 2')
    if cocos is None:
        warnings.warn('Input COCOS is unspecified; assuming COCOS 5. Verify with the equilibrium producer.', UserWarning)
        cocos = 5
    if cocos not in (*range(1, 9), *range(11, 19)):
        raise ValueError('cocos must be 1..8 or 11..18')
    base = cocos % 10
    toroidal_sign = 1 if base % 2 else -1
    q_sign = -1 if base in (1, 2, 7, 8) else 1
    flux_factor = 1.0 if cocos >= 10 else 2 * np.pi
    out = SimpleNamespace(mpol=mpol, lasym=bool(lasym), cocos=cocos, q_sign=q_sign)
    for key in ('sf', 'sp', 'qpsi'):
        values = np.asarray(data[key], dtype=float)
        if values.ndim != 1 or len(values) < 2 or not np.all(np.isfinite(values)):
            raise ValueError(f'{key} must be a finite 1D profile with at least two points')
    q = np.asarray(data['qpsi'], dtype=float)
    if len(q) != len(data['sp']) or len(q) != len(data['sf']):
        raise ValueError('Profile lengths must match')
    if np.any(q == 0) or np.any(np.sign(q) != np.sign(q[0])):
        raise ValueError('q must be nonzero and have one sign for nested toroidal-flux surfaces')
    delta = data['psilim'] - data['psiaxis']
    if not np.isfinite(delta) or delta == 0 or not np.isfinite(data['btor']) or data['btor'] == 0:
        raise ValueError('Nonzero finite poloidal flux span and toroidal field are required')
    # Integrate in normalized poloidal flux: always increasing, even if psi decreases.
    u = np.linspace(0, 1, len(q))
    primitive = PchipInterpolator(u, q).antiderivative()
    integral = primitive(u) - primitive(0)
    out.tflux = integral / integral[-1]
    out.phiedge = abs(delta * integral[-1]) * flux_factor * np.sign(data['btor']) * toroidal_sign
    out.iotaf = q_sign / q
    out.curtor = float(data['totcur']) * toroidal_sign
    out.raxis, out.zaxis = float(data['xaxis']), float(data['zaxis'])
    # VMEC2000 vparams limits each auxiliary spline to 101 knots. Keep source
    # samples, selecting additional knots where linear interpolation errs most.
    selected = [0, len(q)-1]
    profiles = np.array([data['sp'], out.iotaf], dtype=float)
    scales = np.maximum(np.max(np.abs(profiles), axis=1), np.finfo(float).tiny)
    while len(selected) < min(101, len(q)):
        selected.sort()
        approximation = np.array([np.interp(out.tflux, out.tflux[selected], profile[selected])
                                  for profile in profiles])
        error = np.max(np.abs(profiles-approximation)/scales[:, None], axis=0)
        error[selected] = -1
        selected.append(int(np.argmax(error)))
    out.profile_indices = np.array(sorted(selected))
    out.am_aux_s = out.tflux[out.profile_indices]
    out.am_aux_f = profiles[0, out.profile_indices]
    out.ai_aux_s = out.am_aux_s.copy()
    out.ai_aux_f = profiles[1, out.profile_indices]
    out.am = endpoint_polynomial(out.am_aux_s, out.am_aux_f)
    out.ai = endpoint_polynomial(out.ai_aux_s, out.ai_aux_f)
    # No current profile is inferred from p' + FF': these have incompatible units.
    out.ac = out.ac_aux_s = out.ac_aux_f = np.array([])

    r, z = np.asarray(data['xbndry']), np.asarray(data['zbndry'])
    if r.ndim != 1 or z.shape != r.shape or len(r) < 2 * mpol + 1 or not np.all(np.isfinite([r, z])):
        raise ValueError('Boundary needs finite paired coordinates and at least 2*mpol+1 points')
    if np.allclose([r[0], z[0]], [r[-1], z[-1]], rtol=0, atol=1e-10):
        r, z = r[:-1], z[:-1]
    center_z = out.zaxis if lasym else 0.0
    angle = np.mod(np.arctan2(z - center_z, r - out.raxis), 2 * np.pi)
    # A polar graph is appropriate for ordinary tokamak boundaries, but not arbitrary contours.
    steps = (np.diff(np.r_[angle, angle[0]]) + np.pi) % (2*np.pi) - np.pi
    nonzero = steps[np.abs(steps) > 1e-10]
    if len(nonzero) == 0 or not (np.all(nonzero > 0) or np.all(nonzero < 0)):
        raise ValueError('Boundary must be ordered and star-shaped about the chosen center')
    order = np.argsort(angle)
    angle, r, z = angle[order], r[order], z[order]
    keep = np.r_[True, np.diff(angle) > 1e-10]
    angle, r, z = angle[keep], r[keep], z[keep]
    if len(angle) < 2 * mpol + 1:
        raise ValueError('Too few distinct boundary angles for requested mpol')
    theta = np.linspace(0, 2*np.pi, max(2048, 16*mpol), endpoint=False)
    # Resample at actual geometric angles, not boundary-point indices.
    rr = np.interp(theta, angle, r, period=2*np.pi)
    zz = np.interp(theta, angle, z, period=2*np.pi)
    rc, zc = np.fft.rfft(rr)/len(rr), np.fft.rfft(zz)/len(zz)
    out.refou = (2*rc[:mpol].real)[:, None]
    out.refou2 = (-2*rc[:mpol].imag)[:, None]
    out.zefou = (-2*zc[:mpol].imag)[:, None]
    out.zefou2 = (2*zc[:mpol].real)[:, None]
    out.refou[0] *= .5
    out.zefou2[0] *= .5
    if not lasym:
        modes = theta[:, None] * np.arange(mpol)
        removed_r = np.sin(modes) @ out.refou2[:, 0]
        removed_z = np.cos(modes) @ out.zefou2[:, 0]
        out.symmetry_removed_max_m = float(np.max(np.hypot(removed_r, removed_z)))
        if out.symmetry_removed_max_m > 1e-3:
            warnings.warn(f'LASYM=0 removes asymmetric boundary components up to {out.symmetry_removed_max_m:.4g} m; inspect the comparison plot or use lasym=True.', UserWarning)
        out.refou2[:] = 0
        out.zefou2[:] = 0
    else:
        out.symmetry_removed_max_m = 0.0
    return out
