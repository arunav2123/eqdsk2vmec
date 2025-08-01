# GEQDSK to VMEC Converter

A Python package to convert GEQDSK files to VMEC input format, translated from MATLAB code.

## Overview

This package provides functionality to read GEQDSK (Equilibrium Disk) files and convert them to VMEC (Variational Moments Equilibrium Code) input format. It is based on the translation of MATLAB code that performs similar conversions for tokamak plasma equilibrium data.

## Features

- Read GEQDSK files 
- Convert GEQDSK data to VMEC input parameters
- Generate VMEC namelist input files
- Support for various VMEC configuration parameters

## Installation

### From source

```bash
git clone <repository-url>
cd eqdsk2vmec_converter
pip install -e .
```

### Dependencies

- numpy
- freeqdsk

## Usage

### Basic Usage

```python
from eqdsk2vmec_converter.eqdsk2vmec import convert_eqdsk_to_vmec

# Convert an EQDSK file to VMEC input
convert_eqdsk_to_vmec('g02895743.geq')
```

This will create a VMEC input file named `input.g02895743` in the current directory.

### Advanced Usage

You can also use the individual components:

```python
from eqdsk2vmec_converter.read_geqdsk import read_geqdsk
from eqdsk2vmec_converter.eqdisk2vmec_inputfile import eqdisk2vmec_inputfile
from eqdsk2vmec_converter.vmec_namelist import vmec_namelist_init, write_vmec_input

# Read EQDSK file
gdata = read_geqdsk('your_file.geq')

# Convert to VMEC data structure
data = eqdisk2vmec_inputfile('your_file.geq', gdata)

# Initialize VMEC namelist and customize parameters
vmec_input = vmec_namelist_init('indata')
vmec_input.delt = 1.0
vmec_input.niter = 20000
# ... set other parameters

# Write VMEC input file
write_vmec_input('input.your_file', vmec_input)
```

## File Structure

```
eqdsk2vmec_converter/
├── eqdsk2vmec_converter/
│   ├── __init__.py
│   ├── read_geqdsk.py          # EQDSK file reader
│   ├── eqdisk2vmec_inputfile.py # EQDSK to VMEC data converter
│   ├── vmec_namelist.py        # VMEC namelist handling
│   └── eqdsk2vmec.py          # Main conversion function
├── tests/
│   └── __init__.py
├── setup.py
└── README.md
```

## VMEC Parameters

The converter sets the following VMEC parameters based on the original MATLAB code:

- `delt = 1.0`
- `niter = 20000`
- `tcon0 = 1.0`
- `ns_array = [16, 32, 64, 128]`
- `ftol_array = [1e-30, 1e-30, 1e-30, 1e-12]`
- `niter_array = [1000, 2000, 4000, 20000]`
- `lasym = 0`
- `nfp = 1`
- `ntor = 0`
- `nstep = 200`
- `lfreeb = 0`
- `nvacskip = 6`
- `gamma = 0.0`
- `bloat = 1.0`
- `spres_ped = 1.0`
- `pres_scale = 1.0`
- Various profile types set to 'akima_spline'

## Notes

- This is a direct translation from MATLAB code and may require further refinement for specific use cases
- The conversion of boundary coefficients (Fourier modes) is currently simplified and may need enhancement
- Profile data conversion assumes specific mappings from EQDSK to VMEC parameters

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the MIT License.

