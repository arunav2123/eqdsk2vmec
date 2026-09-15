# Written by Arunav Kumar, MIT Plasma Science and Fusion center, 10th May, 2026
"""
EQDSK to VMEC Converter

A Python package to convert EQDSK files to VMEC input format.
"""

from .eqdsk2vmec import convert_eqdsk_to_vmec
from .read_geqdsk import read_geqdsk
from .eqdisk2vmec_inputfile import eqdisk2vmec_inputfile
from .vmec_namelist import vmec_namelist_init, write_vmec_input

__version__ = "0.1.0"
__all__ = [
    "convert_eqdsk_to_vmec",
    "read_geqdsk",
    "eqdisk2vmec_inputfile", 
    "vmec_namelist_init",
    "write_vmec_input"
]

