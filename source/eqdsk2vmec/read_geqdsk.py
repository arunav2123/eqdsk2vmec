# Written by Arunav Kumar, MIT Plasma Science and Fusion center, 10th May, 2026
import numpy as np

def read_geqdsk(filename):
    """
    Reads a GEQDSK file and returns a dictionary of equilibrium data.

    This version uses dedicated functions for parsing different line formats,
    providing more robust and accurate reading.

    Args:
        filename (str): The path to the GEQDSK file.

    Returns:
        dict: A dictionary containing the parsed equilibrium data.
    """
    efit_data = {}
    rotate = True

    with open(filename, 'r') as f:
        # --- Helper functions for reading specific formats ---
        def _read_floats_fixed_width(file_handle, count, width=16):
            """Reads a specified number of floats from fixed-width fields."""
            values = []
            while len(values) < count:
                line = file_handle.readline()
                if not line:
                    raise ValueError(f"Truncated GEQDSK: expected {count} values, got {len(values)}")
                line = line.rstrip('\n\r')
                pos = 0
                while pos < len(line) and len(values) < count:
                    chunk = line[pos : pos + width]
                    if chunk.strip():
                        try:
                            value = float(chunk.replace('D', 'E').replace('d', 'e'))
                        except ValueError as exc:
                            raise ValueError(f'Invalid GEQDSK numeric field: {chunk!r}') from exc
                        if not np.isfinite(value):
                            raise ValueError('Non-finite GEQDSK numeric field')
                        values.append(value)
                    pos += width
            return np.array(values)

        def _read_integers_on_line(file_handle, count):
            """Reads a single line and parses integers from it."""
            line = file_handle.readline()
            if not line:
                return np.array([])
            
            try:
                values = [int(line[i*5:(i+1)*5]) for i in range(count)]
            except ValueError:
                values = [int(p) for p in line.split()]
            if len(values) != count or any(v < 0 for v in values):
                raise ValueError('Invalid boundary/limiter counts')
            return np.array(values)
        
        # --- Main parsing logic ---
        
        # Read header line
        header = f.readline()
        if not header.strip():
            raise ValueError("Invalid GEQDSK header format: too short.")
        
        header_parts = header.split()[-3:]
        if len(header_parts) < 3:
            raise ValueError("Insufficient data in GEQDSK header for ipest, nx, nz.")
        
        try:
            efit_data['ipest'] = int(header_parts[0])
            nx = int(header_parts[1])
            nz = int(header_parts[2])
        except ValueError:
            raise ValueError("Failed to parse integers (ipest, nx, nz) from GEQDSK header.")

        if nx < 2 or nz < 2:
            raise ValueError('GEQDSK grid dimensions must be >= 2')
        efit_data['nx'] = nx
        efit_data['nz'] = nz
        
        # Read the next 20 floats in one block
        initial_params = _read_floats_fixed_width(f, 20)
        
        # Assign values based on the order
        efit_data['xdim'] = initial_params[0]
        efit_data['zdim'] = initial_params[1]
        efit_data['rcentr'] = initial_params[2]
        efit_data['redge'] = initial_params[3]
        efit_data['zmid'] = initial_params[4]
        efit_data['xaxis'] = initial_params[5]
        efit_data['zaxis'] = initial_params[6]
        efit_data['psiaxis'] = initial_params[7]
        efit_data['psilim'] = initial_params[8]
        efit_data['btor'] = initial_params[9]
        efit_data['totcur'] = initial_params[10]

        # Read the remaining grids
        efit_data['sf'] = _read_floats_fixed_width(f, nx)
        efit_data['sp'] = _read_floats_fixed_width(f, nx)
        efit_data['sffp'] = _read_floats_fixed_width(f, nx)
        efit_data['spp'] = _read_floats_fixed_width(f, nx)
        
        psixz_flat = _read_floats_fixed_width(f, nx * nz)
        efit_data['psixz'] = psixz_flat.reshape((nx, nz), order='F')

        efit_data['qpsi'] = _read_floats_fixed_width(f, nx)

        # Read boundary and limiter dimensions
        bdry_lim_counts = _read_integers_on_line(f, 2)
        if len(bdry_lim_counts) != 2:
            raise ValueError("Failed to read boundary and limiter counts.")
        
        efit_data['nbndry'] = int(bdry_lim_counts[0])
        efit_data['nlim'] = int(bdry_lim_counts[1])
        
        # Read boundary coordinates
        if efit_data['nbndry'] > 0:
            bdry_coords_flat = _read_floats_fixed_width(f, efit_data['nbndry'] * 2)
            if len(bdry_coords_flat) != efit_data['nbndry'] * 2:
                raise ValueError(f"Expected {efit_data['nbndry'] * 2} boundary coordinates, but read {len(bdry_coords_flat)}.")
            efit_data['xbndry'] = bdry_coords_flat[::2]
            efit_data['zbndry'] = bdry_coords_flat[1::2]
        else:
            efit_data['xbndry'] = np.array([])
            efit_data['zbndry'] = np.array([])

        # Read limiter coordinates
        if efit_data['nlim'] > 0:
            lim_coords_flat = _read_floats_fixed_width(f, efit_data['nlim'] * 2)
            if len(lim_coords_flat) != efit_data['nlim'] * 2:
                raise ValueError(f"Expected {efit_data['nlim'] * 2} limiter coordinates, but read {len(lim_coords_flat)}.")
            efit_data['xlim'] = lim_coords_flat[::2]
            efit_data['zlim'] = lim_coords_flat[1::2]
        else:
            efit_data['xlim'] = np.array([])
            efit_data['zlim'] = np.array([])

        # Optional: rotation and mass density
        if rotate:
            # kvtor: often a single integer on its own line
            try:
                kvtor_val = int(f.readline().strip().split()[0])
            except (ValueError, IndexError):
                kvtor_val = 0
            efit_data['kvtor'] = kvtor_val

            if efit_data['kvtor'] > 0:
                efit_data['rvtor'] = _read_floats_fixed_width(f, 1)[0]
                efit_data['pressw'] = _read_floats_fixed_width(f, nx)
                efit_data['pwprim'] = _read_floats_fixed_width(f, nx)
            else:
                efit_data['rvtor'] = None
                efit_data['pressw'] = np.array([])
                efit_data['pwprim'] = np.array([])

            try:
                nmass_val = int(f.readline().strip().split()[0])
            except (ValueError, IndexError):
                nmass_val = 0
            efit_data['nmass'] = nmass_val

            if efit_data['nmass'] > 0:
                efit_data['rho0'] = _read_floats_fixed_width(f, nx)
            else:
                efit_data['rho0'] = np.array([])
        else:
            efit_data['kvtor'] = 0
            efit_data['rvtor'] = None
            efit_data['nmass'] = 0
            efit_data['pressw'] = np.array([])
            efit_data['pwprim'] = np.array([])
            efit_data['rho0'] = np.array([])

    efit_data['datatype'] = 'EFIT_G'
    efit_data['xgrid'] = efit_data['redge'] + efit_data['xdim'] * np.arange(nx) / (nx - 1)
    efit_data['zgrid'] = efit_data['zmid'] - efit_data['zdim'] / 2.0 + efit_data['zdim'] * np.arange(nz) / (nz - 1)
    
    return efit_data
