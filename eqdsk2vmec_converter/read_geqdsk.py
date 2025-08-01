#from freeqdsk import read_geqdsk as freeqdsk_read_geqdsk
import numpy as np

def read_geqdsk(filename):
    """
    Reads a GEQDSK file and returns a dictionary with all relevant fields.
    Faithful to the MATLAB read_geqdsk.m logic.
    """
    efit_data = {}
    with open(filename, 'r') as f:
        header = f.readline()  # Read and skip header
        # Parse header fields
        temp = [int(x) for x in header[49:].split()[:3]]
        efit_data['ipest'] = temp[0]
        nx = temp[1]
        nz = temp[2]
        efit_data['nx'] = nx
        efit_data['nz'] = nz

        def fread_float(count=1):
            vals = []
            while len(vals) < count:
                vals.extend([float(x) for x in f.readline().split()])
            return vals if count > 1 else vals[0]

        efit_data['xdim'] = fread_float()
        efit_data['zdim'] = fread_float()
        efit_data['zc'] = fread_float()
        efit_data['redge'] = fread_float()
        efit_data['zmid'] = fread_float()
        efit_data['xaxis'] = fread_float()
        efit_data['zaxis'] = fread_float()
        efit_data['psiaxis'] = fread_float()
        efit_data['psilim'] = fread_float()
        efit_data['btor'] = fread_float()
        efit_data['totcur'] = fread_float()
        efit_data['psimx'] = [fread_float(), fread_float()]
        efit_data['xax'] = [fread_float(), fread_float()]
        efit_data['zax'] = [fread_float(), fread_float()]
        efit_data['psisep'] = fread_float()
        efit_data['xsep'] = fread_float()
        efit_data['zsep'] = fread_float()
        efit_data['sf'] = np.array(fread_float(nx))
        efit_data['sp'] = np.array(fread_float(nx))
        efit_data['sffp'] = np.array(fread_float(nx))
        efit_data['spp'] = np.array(fread_float(nx))
        # Read psixz as a (nx, nz) Fortran-order array
        psixz_flat = fread_float(nx * nz)
        efit_data['psixz'] = np.array(psixz_flat).reshape((nx, nz), order='F')
        efit_data['qpsi'] = np.array(fread_float(nx))
        efit_data['nbndry'] = int(fread_float())
        efit_data['nlim'] = int(fread_float())
        bdry_temp = fread_float(efit_data['nbndry'] * 2)
        efit_data['xbndry'] = np.array(bdry_temp[::2])
        efit_data['zbndry'] = np.array(bdry_temp[1::2])
        lim_temp = fread_float(efit_data['nlim'] * 2)
        efit_data['xlim'] = np.array(lim_temp[::2])
        efit_data['zlim'] = np.array(lim_temp[1::2])
        # Optional: rotation and mass density (if present)
        try:
            efit_data['kvtor'] = int(f.readline().split()[0])
            efit_data['rvtor'] = fread_float()
            efit_data['nmass'] = int(f.readline().split()[0])
            if efit_data['kvtor'] > 0:
                efit_data['pressw'] = np.array(fread_float(nx))
                efit_data['pwprim'] = np.array(fread_float(nx))
            if efit_data['nmass'] > 0:
                efit_data['rho0'] = np.array(fread_float(nx))
        except Exception:
            pass  # These fields are optional
    efit_data['datatype'] = 'EFIT_G'
    efit_data['xgrid'] = efit_data['redge'] + efit_data['xdim'] * np.arange(nx) / (nx - 1)
    efit_data['zgrid'] = efit_data['zdim'] * np.arange(nz) / (nz - 1) - efit_data['zdim'] / 2.0
    return efit_data