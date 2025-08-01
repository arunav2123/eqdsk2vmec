import numpy as np
import datetime

class VMECNamelist:
    """VMEC namelist data structure"""

    def __init__(self):
        self.datatype = 'VMEC_input'
        self.delt = 1.0
        self.niter = 20000
        self.nstep = 200
        self.tcon0 = 1.0
        self.ns_array = [16, 32, 64, 128]
        self.ftol_array = [1e-30, 1e-30, 1e-30, 1e-12]
        self.niter_array = [1000, 2000, 4000, 20000]
        self.lasym = 0
        self.nfp = 1
        self.mpol = 5
        self.ntor = 0
        self.ntheta = 16
        self.nzeta = 0
        self.phiedge = 1.0
        self.lfreeb = 0
        self.mgrid_file = ''
        self.extcur = [0, 0, 0]
        self.nvacskip = 6
        self.gamma = 0.0
        self.bloat = 1.0
        self.spres_ped = 1.0
        self.pres_scale = 1.0
        self.pmass_type = 'akima_spline'
        self.am = []
        self.am_aux_s = []
        self.am_aux_f = []
        self.pcurr_type = 'akima_spline_ip'
        self.curtor = 0.0
        self.ncurr = 0
        self.ac = []
        self.ac_aux_s = []
        self.ac_aux_f = []
        self.piota_type = 'akima_spline'
        self.ai = []
        self.ai_aux_s = []
        self.ai_aux_f = []
        self.raxis = []
        self.zaxis = []
        self.raxis_cc = []
        self.raxis_cs = []
        self.zaxis_cc = []
        self.zaxis_cs = []
        self.rbc = np.array([])
        self.zbs = np.array([])
        self.rbs = np.array([])
        self.zbc = np.array([])

def vmec_namelist_init(namelist_name):
    """Initialize VMEC namelist"""
    return VMECNamelist()

def write_namelist_int(f, name, val):
    f.write(f"  {name.upper()} = {int(val)}\n")

def write_namelist_flt(f, name, val):
    f.write(f"  {name.upper()} = {float(val):21.14E}\n")

def write_namelist_boo(f, name, val):
    s = 'T' if bool(val) else 'F'
    f.write(f"  {name.upper()} = {s}\n")

def write_namelist_str(f, name, val):
    f.write(f"  {name.upper()} = '{val}'\n")

def write_namelist_vec(f, name, arr, typ=None, per_line=8):
    arr = np.atleast_1d(arr)
    if typ == 'int':
        f.write(f"  {name.upper()} =")
        for i, v in enumerate(arr):
            f.write(f" {int(v)}")
            if (i+1) % per_line == 0 and (i+1) < len(arr):
                f.write("\n    ")
        f.write("\n")
    elif typ == 'flt' or (typ is None and np.issubdtype(arr.dtype, np.floating)):
        f.write(f"  {name.upper()} =")
        for i, v in enumerate(arr):
            f.write(f" {float(v):.12E}")
            if (i+1) % per_line == 0 and (i+1) < len(arr):
                f.write("\n    ")
        f.write("\n")
    else:
        f.write(f"  {name.upper()} =")
        for i, v in enumerate(arr):
            f.write(f" {v}")
            if (i+1) % per_line == 0 and (i+1) < len(arr):
                f.write("\n    ")
        f.write("\n")

def write_vmec_input(filename, data):
    """
    Write a VMEC input structure to a file.
    """
    if not hasattr(data, 'datatype') or data.datatype != 'VMEC_input':
        print('Error: Not VMEC input data!')
        return -1

    with open(filename, 'w') as f:
        f.write("&INDATA\n")
        f.write("!----- Runtime Parameters -----\n")
        write_namelist_flt(f, 'delt', getattr(data, 'delt', 1.0))
        write_namelist_int(f, 'niter', getattr(data, 'niter', 20000))
        write_namelist_int(f, 'nstep', getattr(data, 'nstep', 200))
        write_namelist_flt(f, 'tcon0', getattr(data, 'tcon0', 1.0))

        ns = len(getattr(data, 'ns_array', []))
        if ns > 0:
            write_namelist_vec(f, 'NS_ARRAY', data.ns_array, typ='int')
        niter_array = getattr(data, 'niter_array', [])
        if len(niter_array) == ns:
            write_namelist_vec(f, 'NITER_ARRAY', data.niter_array, typ='int')
        ftol_array = getattr(data, 'ftol_array', [])
        if len(ftol_array) == ns:
            write_namelist_vec(f, 'FTOL_ARRAY', data.ftol_array, typ='flt')
        f.write("!----- Grid Parameters -----\n")
        write_namelist_boo(f, 'lasym', getattr(data, 'lasym', 0))
        write_namelist_int(f, 'nfp', getattr(data, 'nfp', 1))
        write_namelist_int(f, 'mpol', getattr(data, 'mpol', 5))
        write_namelist_int(f, 'ntor', getattr(data, 'ntor', 0))
        if hasattr(data, 'ntheta'):
            write_namelist_int(f, 'ntheta', getattr(data, 'ntheta', 16))
        if hasattr(data, 'nzeta'):
            write_namelist_int(f, 'nzeta', getattr(data, 'nzeta', 0))
        write_namelist_flt(f, 'phiedge', getattr(data, 'phiedge', 1.0))

        f.write("!----- Free Boundary Parameters -----\n")
        write_namelist_boo(f, 'lfreeb', getattr(data, 'lfreeb', 0))
        if getattr(data, 'lfreeb', 0):
            write_namelist_str(f, 'mgrid_file', getattr(data, 'mgrid_file', ''))
            write_namelist_vec(f, 'extcur', getattr(data, 'extcur', []), typ='flt')
            write_namelist_int(f, 'nvacskip', getattr(data, 'nvacskip', 6))
        else:
            if hasattr(data, 'extcur') and np.any(np.array(data.extcur) != 0):
                write_namelist_vec(f, 'extcur', data.extcur, typ='flt')
            if hasattr(data, 'mgrid_file') and data.mgrid_file:
                write_namelist_str(f, 'mgrid_file', data.mgrid_file)

        f.write("!----- Pressure Parameters -----\n")
        write_namelist_flt(f, 'gamma', getattr(data, 'gamma', 0.0))
        write_namelist_flt(f, 'bloat', getattr(data, 'bloat', 1.0))
        write_namelist_flt(f, 'spres_ped', getattr(data, 'spres_ped', 1.0))
        write_namelist_str(f, 'pmass_type', getattr(data, 'pmass_type', 'akima_spline'))
        pres_scale = getattr(data, 'pres_scale', None)
        if pres_scale is not None:
            if isinstance(pres_scale, (list, np.ndarray)) and len(pres_scale) > 1:
                write_namelist_vec(f, 'pres_scale', pres_scale, typ='flt')
            else:
                write_namelist_flt(f, 'pres_scale', pres_scale)
        # am, am_aux_s, am_aux_f
        if hasattr(data, 'am'):
            write_namelist_vec(f, 'am', data.am, typ='flt')
        if hasattr(data, 'am_aux_s') and hasattr(data, 'am_aux_f'):
            dex = len(data.am_aux_s)
            write_namelist_vec(f, 'AM_AUX_S', data.am_aux_s[:dex], typ='flt')
            write_namelist_vec(f, 'AM_AUX_F', data.am_aux_f[:dex], typ='flt')

        # FLOW/ANISOTROPY Parameters (optional)
        if hasattr(data, 'bcrit'):
            f.write("!----- FLOW/ANISOTROPY Parameters -----\n")
            write_namelist_flt(f, 'bcrit', getattr(data, 'bcrit', 0.0))
            write_namelist_str(f, 'ph_type', getattr(data, 'ph_type', 'power_series'))
            if hasattr(data, 'ah'):
                write_namelist_vec(f, 'ah', data.ah, typ='flt')
            if hasattr(data, 'ah_aux_s') and hasattr(data, 'ah_aux_f'):
                dex = len(data.ah_aux_s)
                write_namelist_vec(f, 'AH_AUX_S', data.ah_aux_s[:dex], typ='flt')
                write_namelist_vec(f, 'AH_AUX_F', data.ah_aux_f[:dex], typ='flt')
            write_namelist_str(f, 'pt_type', getattr(data, 'pt_type', 'power_series'))
            if hasattr(data, 'at'):
                write_namelist_vec(f, 'at', data.at, typ='flt')
            if hasattr(data, 'at_aux_s') and hasattr(data, 'at_aux_f'):
                dex = len(data.at_aux_s)
                write_namelist_vec(f, 'AT_AUX_S', data.at_aux_s[:dex], typ='flt')
                write_namelist_vec(f, 'AT_AUX_F', data.at_aux_f[:dex], typ='flt')

        f.write("!----- Current/Iota Parameters -----\n")
        write_namelist_flt(f, 'curtor', getattr(data, 'curtor', 0.0))
        write_namelist_int(f, 'ncurr', getattr(data, 'ncurr', 0))
        if hasattr(data, 'ac_form'):
            write_namelist_int(f, 'AC_FORM', getattr(data, 'ac_form', 0))
        write_namelist_str(f, 'piota_type', getattr(data, 'piota_type', 'akima_spline'))
        if hasattr(data, 'ai'):
            write_namelist_vec(f, 'ai', data.ai, typ='flt')
        if hasattr(data, 'ai_aux_s') and hasattr(data, 'ai_aux_f'):
            dex = len(data.ai_aux_s)
            write_namelist_vec(f, 'AI_AUX_S', data.ai_aux_s[:dex], typ='flt')
            write_namelist_vec(f, 'AI_AUX_F', data.ai_aux_f[:dex], typ='flt')
        write_namelist_str(f, 'pcurr_type', getattr(data, 'pcurr_type', 'akima_spline_ip'))
        if hasattr(data, 'ac'):
            write_namelist_vec(f, 'ac', data.ac, typ='flt')
        if hasattr(data, 'ac_aux_s') and hasattr(data, 'ac_aux_f'):
            dex = len(data.ac_aux_s)
            write_namelist_vec(f, 'AC_AUX_S', data.ac_aux_s[:dex], typ='flt')
            write_namelist_vec(f, 'AC_AUX_F', data.ac_aux_f[:dex], typ='flt')

        f.write("!----- Axis Parameters -----\n")
        if hasattr(data, 'raxis') and np.size(getattr(data, 'raxis', [])) > 0:
            write_namelist_vec(f, 'RAXIS', data.raxis, typ='flt')
        elif hasattr(data, 'raxis_cc') and np.size(getattr(data, 'raaxis_cc', [])) > 0:
            write_namelist_vec(f, 'RAXIS_CC', data.raxis_cc, typ='flt')
            if getattr(data, 'lasym', 0):
                write_namelist_vec(f, 'RAXIS_CS', data.raxis_cs, typ='flt')
        if hasattr(data, 'zaxis') and np.size(getattr(data, 'zaxis', [])) > 0:
            write_namelist_vec(f, 'ZAXIS', data.zaxis, typ='flt')
        elif hasattr(data, 'zaxis_cc') and np.size(getattr(data, 'zaxis_cc', [])) > 0:
            write_namelist_vec(f, 'ZAXIS_CC', data.zaxis_cc, typ='flt')
            if getattr(data, 'lasym', 0):
                write_namelist_vec(f, 'ZAXIS_CS', data.zaxis_cs, typ='flt')

        f.write("!----- Boundary Parameters -----\n")
        # Write boundary Fourier coefficients as in MATLAB
        # 2D arrays, shape: (2*ntor+1, mpol+1) after transpose
        for arr_name in ['rbc', 'zbs', 'rbs', 'zbc']:
            arr = getattr(data, arr_name, None)
            if arr is not None and arr.size > 0:
                arr = np.atleast_2d(arr)
                setattr(data, arr_name, arr)

        mpol = getattr(data, 'mpol', 5)
        ntor = getattr(data, 'ntor', 0)
        rbc = getattr(data, 'rbc', np.zeros((2*ntor+1, mpol+1)))
        zbs = getattr(data, 'zbs', np.zeros((2*ntor+1, mpol+1)))
        rbs = getattr(data, 'rbs', np.zeros((2*ntor+1, mpol+1)))
        zbc = getattr(data, 'zbc', np.zeros((2*ntor+1, mpol+1)))
        lasym = getattr(data, 'lasym', 0)

        for i in range(2*ntor+1):
            for j in range(mpol+1):
                # Only print nonzero boundary coefficients for efficiency/readability
                if (rbc[i, j] != 0.0) or (zbs[i, j] != 0.0):
                    f.write(f"  RBC({i-ntor},{j}) = {rbc[i,j]:17.10e}  ZBS({i-ntor},{j}) = {zbs[i,j]:17.10e}\n")
                    if lasym:
                        f.write(f"    RBS({i-ntor},{j}) = {rbs[i,j]:17.10e}  ZBC({i-ntor},{j}) = {zbc[i,j]:17.10e}\n")

        import datetime
        f.write(f"!----- Created by write_vmec_input {datetime.datetime.now().isoformat()} -----\n")
        f.write("/\n")
    return 1
