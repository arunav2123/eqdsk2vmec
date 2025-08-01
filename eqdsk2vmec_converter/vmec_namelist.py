import numpy as np

class VMECNamelist:
    """VMEC namelist data structure"""
    
    def __init__(self):
        # Default VMEC parameters
        self.delt = 1.0
        self.niter = 20000
        self.tcon0 = 1.0
        self.ns_array = [16, 32, 64, 128]
        self.ftol_array = [1e-30, 1e-30, 1e-30, 1e-12]
        self.niter_array = [1000, 2000, 4000, 20000]
        self.lasym = 0
        self.nfp = 1
        self.mpol = 5
        self.ntor = 0
        self.nstep = 200
        self.ntheta = 16
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

def write_vmec_input(filename, vmec_input):
    """Write VMEC input file"""
    print(f"Writing VMEC input file: {filename}")
    
    with open(filename, 'w') as f:
        f.write("&INDATA\n")
        
        # Write scalar parameters
        f.write(f"  DELT = {vmec_input.delt}\n")
        f.write(f"  NITER = {vmec_input.niter}\n")
        f.write(f"  TCON0 = {vmec_input.tcon0}\n")
        f.write(f"  LASYM = {vmec_input.lasym}\n")
        f.write(f"  NFP = {vmec_input.nfp}\n")
        f.write(f"  MPOL = {vmec_input.mpol}\n")
        f.write(f"  NTOR = {vmec_input.ntor}\n")
        f.write(f"  NSTEP = {vmec_input.nstep}\n")
        f.write(f"  NTHETA = {vmec_input.ntheta}\n")
        f.write(f"  PHIEDGE = {vmec_input.phiedge}\n")
        f.write(f"  LFREEB = {vmec_input.lfreeb}\n")
        f.write(f"  MGRID_FILE = '{vmec_input.mgrid_file}'\n")
        f.write(f"  NVACSKIP = {vmec_input.nvacskip}\n")
        f.write(f"  GAMMA = {vmec_input.gamma}\n")
        f.write(f"  BLOAT = {vmec_input.bloat}\n")
        f.write(f"  SPRES_PED = {vmec_input.spres_ped}\n")
        f.write(f"  PRES_SCALE = {vmec_input.pres_scale}\n")
        f.write(f"  PMASS_TYPE = '{vmec_input.pmass_type}'\n")
        f.write(f"  PCURR_TYPE = '{vmec_input.pcurr_type}'\n")
        f.write(f"  CURTOR = {vmec_input.curtor}\n")
        f.write(f"  NCURR = {vmec_input.ncurr}\n")
        f.write(f"  PIOTA_TYPE = '{vmec_input.piota_type}'\n")
        
        # Write array parameters
        if vmec_input.ns_array:
            f.write(f"  NS_ARRAY = {' '.join(map(str, vmec_input.ns_array))}\n")
        if vmec_input.ftol_array:
            f.write(f"  FTOL_ARRAY = {' '.join([f'{x:.2e}' for x in vmec_input.ftol_array])}\n")
        if vmec_input.niter_array:
            f.write(f"  NITER_ARRAY = {' '.join(map(str, vmec_input.niter_array))}\n")
        if vmec_input.extcur:
            f.write(f"  EXTCUR = {' '.join(map(str, vmec_input.extcur))}\n")
        
        # Write profile arrays
        if vmec_input.am:
            f.write(f"  AM = {' '.join(map(str, vmec_input.am))}\n")
        if vmec_input.am_aux_s:
            f.write(f"  AM_AUX_S = {' '.join(map(str, vmec_input.am_aux_s))}\n")
        if vmec_input.am_aux_f:
            f.write(f"  AM_AUX_F = {' '.join(map(str, vmec_input.am_aux_f))}\n")
        if vmec_input.ac:
            f.write(f"  AC = {' '.join(map(str, vmec_input.ac))}\n")
        if vmec_input.ac_aux_s:
            f.write(f"  AC_AUX_S = {' '.join(map(str, vmec_input.ac_aux_s))}\n")
        if vmec_input.ac_aux_f:
            f.write(f"  AC_AUX_F = {' '.join(map(str, vmec_input.ac_aux_f))}\n")
        if vmec_input.ai:
            f.write(f"  AI = {' '.join(map(str, vmec_input.ai))}\n")
        if vmec_input.ai_aux_s:
            f.write(f"  AI_AUX_S = {' '.join(map(str, vmec_input.ai_aux_s))}\n")
        if vmec_input.ai_aux_f:
            f.write(f"  AI_AUX_F = {' '.join(map(str, vmec_input.ai_aux_f))}\n")
        
        # Write axis coefficients
        if vmec_input.raxis_cc:
            f.write(f"  RAXIS_CC = {' '.join(map(str, vmec_input.raxis_cc))}\n")
        if vmec_input.raxis_cs:
            f.write(f"  RAXIS_CS = {' '.join(map(str, vmec_input.raxis_cs))}\n")
        if vmec_input.zaxis_cc:
            f.write(f"  ZAXIS_CC = {' '.join(map(str, vmec_input.zaxis_cc))}\n")
        if vmec_input.zaxis_cs:
            f.write(f"  ZAXIS_CS = {' '.join(map(str, vmec_input.zaxis_cs))}\n")
        
        # Write boundary coefficients
        if vmec_input.rbc.size > 0:
            for i in range(vmec_input.rbc.shape[0]):
                for j in range(vmec_input.rbc.shape[1]):
                    if vmec_input.rbc[i, j] != 0:
                        f.write(f"  RBC({i},{j}) = {vmec_input.rbc[i, j]}\n")
        
        if vmec_input.zbs.size > 0:
            for i in range(vmec_input.zbs.shape[0]):
                for j in range(vmec_input.zbs.shape[1]):
                    if vmec_input.zbs[i, j] != 0:
                        f.write(f"  ZBS({i},{j}) = {vmec_input.zbs[i, j]}\n")
        
        if vmec_input.rbs.size > 0:
            for i in range(vmec_input.rbs.shape[0]):
                for j in range(vmec_input.rbs.shape[1]):
                    if vmec_input.rbs[i, j] != 0:
                        f.write(f"  RBS({i},{j}) = {vmec_input.rbs[i, j]}\n")
        
        if vmec_input.zbc.size > 0:
            for i in range(vmec_input.zbc.shape[0]):
                for j in range(vmec_input.zbc.shape[1]):
                    if vmec_input.zbc[i, j] != 0:
                        f.write(f"  ZBC({i},{j}) = {vmec_input.zbc[i, j]}\n")
        
        f.write("/\n")
    
    print(f"VMEC input file written successfully: {filename}")

