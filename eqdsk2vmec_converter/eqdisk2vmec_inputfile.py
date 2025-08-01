import numpy as np
from scipy.interpolate import pchip_interpolate
from scipy.integrate import quad

def eqdisk2vmec_inputfile(filename, data):
    """
    Converts EQDSK files to VMEC INDATA namelists.
    Args:
        filename: str, EQDSK filename
        data: dict, GEQDSK data dictionary
    Returns:
        data2: Namespace-like object with VMEC parameters
    """

    class VMECData:
        pass
    data2 = VMECData()

    # 1. Geometry and boundary
    xgrid = np.array(data['xgrid'])
    zgrid = np.array(data['zgrid'])
    psixz = np.array(data['psixz']).T  # Transposed for plotting
    xbndry = np.array(data['xbndry'])
    zbndry = np.array(data['zbndry'])
    xaxis = np.array(data['xaxis'])
    zaxis = np.array(data['zaxis'])
    xlim = np.array(data['xlim'])
    zlim = np.array(data['zlim'])

    mpol = 21
    ntor = 0
    R = xbndry
    Z = zbndry
    npb = len(R)
    rmaj = np.mean(R)
    Rminor = R - rmaj
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
    fnuv[1:] = 2 * fnuv[0]
    for m1 in range(mpol + 1):
        for n1 in range(2 * ntor + 1):
            for i in range(nu):
                for j in range(nv):
                    refou[m1, n1] += R[i] * (cosv[j, n1] * cosu[i, m1] - sinv[j, n1] * sinu[i, m1]) * fnuv[m1]
                    zefou[m1, n1] += Z[i] * (sinv[j, n1] * cosu[i, m1] + cosv[j, n1] * sinu[i, m1]) * fnuv[m1]
                    refou2[m1, n1] += R[i] * (sinv[j, n1] * cosu[i, m1] + cosv[j, n1] * sinu[i, m1]) * fnuv[m1]
                    zefou2[m1, n1] += Z[i] * (cosv[j, n1] * cosu[i, m1] - sinv[j, n1] * sinu[i, m1]) * fnuv[m1]
    data2.refou = refou
    data2.zefou = zefou
    data2.refou2 = refou2
    data2.zefou2 = zefou2

    # 2. Profiles
    sf = np.array(data['sf'])
    psiaxis = data['psiaxis']
    psilim = data['psilim']
    pflux = np.linspace(psiaxis, psilim, len(sf))
    pflux_norm = pflux - pflux[0]
    pflux_norm = pflux_norm / np.max(pflux_norm)

    spp = np.array(data['spp'])
    sffp = np.array(data['sffp'])
    press = np.array(data['sp'])
    qprof = np.array(data['qpsi'])
    iotaf = 1.0 / qprof
    jdotb = spp + sffp

    # dphi = q * dpsi, integrate using pchip spline
    q_spl = lambda x: np.interp(x, pflux, qprof)
    tflux = np.zeros(len(pflux))
    for i in range(len(pflux)):
        tflux[i], _ = quad(q_spl, pflux[0], pflux[i])
    btor = data['btor'] if 'btor' in data else 1.0
    phiedge = abs(tflux[-1]) * 2 * np.pi * np.sign(btor)
    tflux = tflux / tflux[-1]

    s = np.linspace(0, 1, 100)

    # Pressure polynomial fit
    am_aux_s = s
    am_aux_f = pchip_interpolate(tflux, press, s)
    p = np.polyfit(am_aux_s, am_aux_f, 9)
    p[-1] = am_aux_f[0]
    p = np.concatenate(([-np.sum(p[:10]) + am_aux_f[-1]], p))
    am = p[::-1]

    # Current polynomial fit
    ac_aux_s = s
    ac_aux_f = pchip_interpolate(tflux, jdotb, s)
    p = np.polyfit(ac_aux_s, ac_aux_f, 9)
    p[-1] = ac_aux_f[0]
    p = np.concatenate(([-np.sum(p[:10]) + ac_aux_f[-1]], p))
    ac = p[::-1]

    # Iota polynomial fit
    ai_aux_s = s
    ai_aux_f = pchip_interpolate(tflux, iotaf, s)
    p = np.polyfit(ai_aux_s, ai_aux_f, 9)
    p[-1] = ai_aux_f[0]
    p = np.concatenate(([-np.sum(p[:10]) + ai_aux_f[-1]], p))
    ai = p[::-1]

    data2.phiedge = phiedge
    data2.curtor = data['totcur']
    data2.raxis = xaxis
    data2.zaxis = zaxis
    data2.am = am
    data2.am_aux_s = am_aux_s
    data2.am_aux_f = am_aux_f
    data2.ac = ac
    data2.ac_aux_s = ac_aux_s
    data2.ac_aux_f = ac_aux_f
    data2.ai = ai
    data2.ai_aux_s = ai_aux_s
    data2.ai_aux_f = ai_aux_f
    data2.tflux = tflux
    data2.iotaf = iotaf

    # 3. VMEC boundary check
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
    ntheta = 360
    theta_v = np.linspace(0, 2 * np.pi, ntheta)
    zeta = 0
    r = cfunct(theta_v, zeta, rmnc, xm, xn)
    z = cfunct(theta_v, zeta, zmnc, xm, xn) + sfunct(theta_v, zeta, zmns, xm, xn)
    data2.rbry = r
    data2.zbry = z

    return data2

def cfunct(theta, zeta, rmnc, xm, xn):
    """
    Cosine Fourier Transform (VMEC convention)
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

def sfunct(theta, zeta, zmns, xm, xn):
    """
    Sine Fourier Transform (VMEC convention)
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