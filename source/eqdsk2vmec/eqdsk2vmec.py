import os
import numpy as np
from .read_geqdsk import read_geqdsk
from .eqdisk2vmec_inputfile import eqdisk2vmec_inputfile
from .vmec_namelist import vmec_namelist_init, write_vmec_input

def convert_eqdsk_to_vmec(filename):
    """
    Converts an EQDSK file to a VMEC input file.
    
    Args:
        filename (str): The path to the EQDSK file.
    """
    print(f"Starting conversion for {filename}")
    
    # Read GEQDSK data
    gdata = read_geqdsk(filename)
    
    # Convert EQDSK data to VMEC input data structure
    data = eqdisk2vmec_inputfile(filename, gdata)
    
    # Initialize VMEC namelist
    vmec_input = vmec_namelist_init("indata")
    
    # Populate VMEC namelist with data from conversion
    vmec_input.delt = 1.0
    vmec_input.niter = 1000
    vmec_input.tcon0 = 1.0
    vmec_input.ns_array = [16, 32, 64, 128]
    vmec_input.ftol_array = [1e-06, 1e-08, 1e-10, 1e-12]
    vmec_input.niter_array = [1000, 2000, 4000, 20000]
    vmec_input.lasym = 0
    vmec_input.nfp = 1
    # Set mpol and ntor based on the MATLAB script
    vmec_input.mpol = 12
    vmec_input.ntor = 0
    vmec_input.nstep = 200
    vmec_input.ntheta = 2 * vmec_input.mpol + 6
    vmec_input.phiedge = data.phiedge
    vmec_input.lfreeb = 0
    vmec_input.mgrid_file = ''
    vmec_input.extcur = [0, 0, 0]
    vmec_input.nvacskip = 6
    vmec_input.gamma = 0.0
    vmec_input.bloat = 1.0
    vmec_input.spres_ped = 1.0
    vmec_input.pres_scale = 1.0
    vmec_input.pmass_type = 'akima_spline'
    vmec_input.am = data.am
    vmec_input.am_aux_s = data.am_aux_s
    vmec_input.am_aux_f = data.am_aux_f
    vmec_input.pcurr_type = 'akima_spline_ip'
    vmec_input.curtor = data.curtor
    vmec_input.ncurr = 0
    vmec_input.ac = data.ac
    vmec_input.ac_aux_s = data.ac_aux_s
    vmec_input.ac_aux_f = data.ac_aux_f
    vmec_input.piota_type = 'akima_spline'
    vmec_input.ai = data.ai
    vmec_input.ai_aux_s = data.ai_aux_s
    vmec_input.ai_aux_f = data.ai_aux_f

    vmec_input.raxis_cc = np.atleast_1d(gdata['xaxis']) # Using gdata's xaxis
    vmec_input.zaxis_cc = np.atleast_1d(gdata['zaxis']) # Using gdata's zaxis
    vmec_input.raxis_cs = np.atleast_1d(data.zaxis) 
    vmec_input.zaxis_cs = np.atleast_1d(data.zaxis) 

    vmec_input.rbc = data.refou.T
    vmec_input.zbs = data.zefou.T
    vmec_input.rbs = data.refou2.T
    vmec_input.zbc = data.zefou2.T
    
    # Construct new filename
    base_filename = os.path.splitext(os.path.basename(filename))[0]
    new_file = f"input.{base_filename}"
    
    # Write VMEC input file
    write_vmec_input(new_file, vmec_input)
    
    print(f"Conversion complete. VMEC input file saved as {new_file}")

if __name__ == '__main__':
    pass