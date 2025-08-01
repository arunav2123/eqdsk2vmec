import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import PchipInterpolator
from scipy.integrate import quad

def eqdisk2vmec_inputfile(filename, data):
    """
    Converts EQDSK files to VMEC INDATA namelists.

    This is a direct translation of the provided MATLAB script, including
    plotting functionality and custom logic for profile fitting and
    boundary Fourier transforms.

    Args:
        filename: str, EQDSK filename (used for plotting titles)
        data: dict, GEQDSK data dictionary from the read_geqdsk function.

    Returns:
        data2: Namespace-like object with VMEC parameters and other calculated data.
    """

    class VMECData:
        pass
    data2 = VMECData()

    # --- 1. Plotting: GEQDSK Input Data ---
    fig = plt.figure(figsize=(6, 6))

    # Subplot 1 (Top-Left): GEQDSK Flux Surfaces and Boundaries
    # Using subplot2grid to create a single tall plot that spans 3 rows on the left side
    ax1 = plt.subplot2grid((3, 2), (0, 0), rowspan=3)
    ax1.contourf(data['xgrid'], data['zgrid'], data['psixz'].T, 100) # .T is important for correct orientation
    ax1.plot(data['xbndry'], data['zbndry'], 'r', linewidth=2, label='EFIT Boundary')
    ax1.plot(data['xaxis'], data['zaxis'], '+r', linewidth=2, label='Axis')
    ax1.plot(data['xlim'], data['zlim'], 'k', linewidth=2, label='EFIT Surfaces') # The MATLAB comments were slightly confusing, this plots the limiter
    ax1.set_xlabel('R')
    ax1.set_ylabel('Z')
    ax1.set_title(f'GEQDSK from {filename}')
    ax1.set_aspect('equal', adjustable='box')
    ax1.legend()

    # --- 2. Boundary Fourier Transform ---
    mpol = 21
    ntor = 0
    R = data['xbndry']
    Z = data['zbndry']
    npb = len(R)
    rmaj = np.mean(R)
    Rminor = R - rmaj
    
    # Sort boundary points by poloidal angle theta
    theta = np.arctan2(Z, Rminor)
    theta[theta < 0] += 2 * np.pi
    dex = np.argsort(theta)
    theta = theta[dex]
    R = R[dex]
    Z = Z[dex]
    Rminor = Rminor[dex]
    rho = np.sqrt(Rminor ** 2 + Z ** 2)

    nu = npb
    nv = 1
    refou = np.zeros((mpol + 1, 2 * ntor + 1))
    zefou = np.zeros((mpol + 1, 2 * ntor + 1))
    refou2 = np.zeros((mpol + 1, 2 * ntor + 1))
    zefou2 = np.zeros((mpol + 1, 2 * ntor + 1))

    cosu = np.zeros((nu, mpol + 1))
    cosv = np.zeros((nv, 2 * ntor + 1))
    sinu = np.zeros((nu, mpol + 1))
    sinv = np.zeros((nv, 2 * ntor + 1))

    alu = 2 * np.pi / nu
    for i in range(nu):
        for j in range(mpol + 1):
            m = j
            cosu[i, j] = np.cos(m * i * alu)
            sinu[i, j] = np.sin(m * i * alu)
    
    alv = 2 * np.pi / nv
    for i in range(nv):
        for j in range(2 * ntor + 1):
            n = j - ntor
            cosv[i, j] = np.cos(n * i * alv)
            sinv[i, j] = np.sin(n * i * alv)
            
    fnuv = np.zeros(mpol + 1)
    fnuv[0] = 1. / (nu * nv)
    for i in range(1, mpol + 1):
        fnuv[i] = 2 * fnuv[0]

    for m1 in range(mpol + 1):
        for n1 in range(2 * ntor + 1):
            for i in range(nu):
                for j in range(nv):
                    # Note: MATLAB indexing starts at 1, Python at 0.
                    # The loops and calculations are a bit complex, and the original Python
                    # `eqdisk2vmec_inputfile.py` simplified these. This is the direct MATLAB translation.
                    refou[m1, n1] += R[i] * (cosv[j, n1] * cosu[i, m1] - sinv[j, n1] * sinu[i, m1]) * fnuv[m1]
                    zefou[m1, n1] += Z[i] * (sinv[j, n1] * cosu[i, m1] + cosv[j, n1] * sinu[i, m1]) * fnuv[m1]
                    refou2[m1, n1] += R[i] * (sinv[j, n1] * cosu[i, m1] + cosv[j, n1] * sinu[i, m1]) * fnuv[m1]
                    zefou2[m1, n1] += Z[i] * (cosv[j, n1] * cosu[i, m1] - sinv[j, n1] * sinu[i, m1]) * fnuv[m1]
    
    data2.refou = refou
    data2.zefou = zefou
    data2.refou2 = refou2
    data2.zefou2 = zefou2

    # --- 3. Profile Fitting ---
    pflux = np.linspace(data['psiaxis'], data['psilim'], len(data['sf']))
    pflux_norm = (pflux - pflux[0]) / (np.max(pflux) - pflux[0])

    jdotb = data['spp'] + data['sffp']
    press = data['sp']
    qprof = data['qpsi']
    iotaf = 1.0 / qprof

    # Calculate toroidal flux (dphi = q * dpsi), by integrating q(psi)
    # Use PchipInterpolator to create a callable spline object
    q_spl = PchipInterpolator(pflux, qprof)
    
    tflux = np.zeros_like(pflux)
    for i in range(len(pflux)):
        # scipy.integrate.quad is the equivalent of MATLAB's integral
        # The warnings you saw before are a known issue with numerical integration
        # and do not halt the program. They can be addressed by adjusting 'limit'
        # or error tolerances (epsabs, epsrel) if needed, but are often acceptable.
        tflux[i], _ = quad(q_spl, pflux[0], pflux[i], limit=100)
    
    # Check for division by zero before normalization
    if tflux[-1] == 0:
        phiedge = 0.0
        tflux = np.zeros_like(tflux)
    else:
        phiedge = abs(tflux[-1]) * 2.0 * np.pi * np.sign(data['btor'])
        tflux = tflux / tflux[-1]

    s = np.linspace(0, 1, 100)

    # Pressure polynomial fit
    am_aux_s = s
    am_spline = PchipInterpolator(tflux, press)
    am_aux_f = am_spline(s)
    p = np.polyfit(am_aux_s, am_aux_f, 9)
    # This is a custom boundary condition from the MATLAB script
    p_matlab_style = np.copy(p)
    p_matlab_style = np.insert(p_matlab_style, 0, 0)
    p_matlab_style[-1] = am_aux_f[0]
    p_matlab_style[0] = -(np.sum(p_matlab_style[1:-1]) + p_matlab_style[0] - am_aux_f[-1]) # This is the python equivalent of the custom logic from matlab

    am = p_matlab_style[::-1]
    data2.am = am
    data2.am_aux_s = am_aux_s
    data2.am_aux_f = am_aux_f

    # Subplot 2 (Top-Right): Pressure Profile
    ax2 = plt.subplot2grid((3, 2), (0, 1))
    ax2.plot(tflux, press, 'k', label='EQDSK')
    ax2.plot(am_aux_s, am_aux_f, 'or', label='SPLINE')
    ax2.plot(am_aux_s, np.polyval(p, am_aux_s), 'b', label='POLY')
    ax2.set_title('Pressure [Pa]')
    ax2.legend()
    ax2.grid(True)

    # Current polynomial fit
    ac_aux_s = s
    ac_spline = PchipInterpolator(tflux, jdotb)
    ac_aux_f = ac_spline(s)
    p = np.polyfit(ac_aux_s, ac_aux_f, 9)
    p_matlab_style = np.copy(p)
    p_matlab_style = np.insert(p_matlab_style, 0, 0)
    p_matlab_style[-1] = ac_aux_f[0]
    p_matlab_style[0] = -(np.sum(p_matlab_style[1:-1]) + p_matlab_style[0] - ac_aux_f[-1])

    ac = p_matlab_style[::-1]
    data2.ac = ac
    data2.ac_aux_s = ac_aux_s
    data2.ac_aux_f = ac_aux_f

    # Subplot 3 (Middle-Right): Current Profile
    ax3 = plt.subplot2grid((3, 2), (1, 1))
    ax3.plot(tflux, jdotb, 'k', label='EQDSK')
    ax3.plot(ac_aux_s, ac_aux_f, 'or', label='SPLINE')
    ax3.plot(ac_aux_s, np.polyval(p, ac_aux_s), 'b', label='POLY')
    ax3.set_title('<J.B>/<B/R> [A]')
    ax3.legend()
    ax3.grid(True)

    # Iota polynomial fit
    ai_aux_s = s
    ai_spline = PchipInterpolator(tflux, iotaf)
    ai_aux_f = ai_spline(s)
    p = np.polyfit(ai_aux_s, ai_aux_f, 9)
    p_matlab_style = np.copy(p)
    p_matlab_style = np.insert(p_matlab_style, 0, 0)
    p_matlab_style[-1] = ai_aux_f[0]
    p_matlab_style[0] = -(np.sum(p_matlab_style[1:-1]) + p_matlab_style[0] - ai_aux_f[-1])

    ai = p_matlab_style[::-1]
    data2.ai = ai
    data2.ai_aux_s = ai_aux_s
    data2.ai_aux_f = ai_aux_f

    # Subplot 4 (Bottom-Right): q Profile
    ax4 = plt.subplot2grid((3, 2), (2, 1))
    ax4.plot(tflux, qprof, 'k', label='EQDSK')
    # MATLAB's plot(ai_aux_s,1./ai_aux_f) is equivalent to 1/spline, not 1/polyval
    ax4.plot(ai_aux_s, 1.0 / ai_aux_f, 'or', label='SPLINE')
    ax4.plot(ai_aux_s, 1.0 / np.polyval(p, ai_aux_s), 'b', label='POLY')
    ax4.set_title('q')
    ax4.legend()
    ax4.grid(True)
    plt.tight_layout()
    plt.show()

    data2.phiedge = phiedge
    data2.curtor = data['totcur']
    data2.raxis = data['xaxis']
    data2.zaxis = data['zaxis']
    data2.tflux = tflux
    data2.iotaf = iotaf
    
    # --- 4. VMEC Boundary Check & Reconstruction ---
    # This section uses the calculated Fourier coefficients to reconstruct the boundary.
    mn = 0
    xm = []
    xn = []
    rmnc = []
    rmns = []
    zmnc = []
    zmns = []
    for m1 in range(mpol + 1):
        for n1 in range(2 * ntor + 1):
            xm.append(m1)
            xn.append(n1 - ntor)
            rmnc.append(refou[m1, n1])
            rmns.append(refou2[m1, n1])
            zmnc.append(zefou2[m1, n1])
            zmns.append(zefou[m1, n1])
            mn += 1
    
    xm = np.array(xm)
    xn = np.array(xn)
    rmnc = np.array(rmnc).reshape(-1, 1)
    rmns = np.array(rmns).reshape(-1, 1)
    zmnc = np.array(zmnc).reshape(-1, 1)
    zmns = np.array(zmns).reshape(-1, 1)
    
    # Plotting reconstructed boundary
    ntheta = 360
    theta_v = np.linspace(0, 2 * np.pi, ntheta)
    zeta = 0.0
    
    r = _cfunct(theta_v, zeta, rmnc, xm, xn) + _sfunct(theta_v, zeta, rmns, xm, xn)
    z = _cfunct(theta_v, zeta, zmnc, xm, xn) + _sfunct(theta_v, zeta, zmns, xm, xn)
    
    ax1.plot(r, z, '-b', linewidth=0.5, label='VMEC')
    ax1.legend()
    plt.show()

    data2.rbry = r
    data2.zbry = z

    return data2

def _cfunct(theta, zeta, rmnc, xm, xn):
    """
    Cosine Fourier Transform (VMEC convention).
    Vectorized for efficiency.
    """
    theta = np.atleast_1d(theta)
    zeta = np.atleast_1d(zeta)
    f = np.zeros((len(theta), len(zeta)))
    for i, thet in enumerate(theta):
        for j, zet in enumerate(zeta):
            arg = xm * thet - xn * zet
            f[i, j] = np.sum(rmnc.flatten() * np.cos(arg))
    if f.shape[1] == 1:
        f = f.flatten()
    return f

def _sfunct(theta, zeta, zmns, xm, xn):
    """
    Sine Fourier Transform (VMEC convention).
    Vectorized for efficiency.
    """
    theta = np.atleast_1d(theta)
    zeta = np.atleast_1d(zeta)
    f = np.zeros((len(theta), len(zeta)))
    for i, thet in enumerate(theta):
        for j, zet in enumerate(zeta):
            arg = xm * thet - xn * zet
            f[i, j] = np.sum(zmns.flatten() * np.sin(arg))
    if f.shape[1] == 1:
        f = f.flatten()
    return f

