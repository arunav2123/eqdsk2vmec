# Written by Arunav Kumar, MIT Plasma Science and Fusion center, 10th May, 2026
"""Public GEQDSK-to-VMEC conversion entry point."""
from pathlib import Path
import numpy as np
from .read_geqdsk import read_geqdsk
from .eqdisk2vmec_inputfile import eqdisk2vmec_inputfile
from .vmec_namelist import vmec_namelist_init, write_vmec_input


def convert_eqdsk_to_vmec(filename, *, output_path=None, plot_path=None,
                          show=False, mpol=12, lasym=False, cocos=None):
    """Write a fixed-boundary, prescribed-iota input; optionally compare it visually.

    Returns the output Path. LASYM=False explicitly symmetrizes about Z=0.
    Specify the producer's COCOS convention; None warns and assumes COCOS 5.
    No VMEC executable is invoked.
    """
    gdata = read_geqdsk(filename)
    data = eqdisk2vmec_inputfile(filename, gdata, mpol=mpol, lasym=lasym, cocos=cocos)
    vmec_input = vmec_namelist_init('indata')
    vmec_input.ftol_array = [1e-6, 1e-8, 1e-10, 1e-12]
    vmec_input.lasym = lasym
    vmec_input.mpol = mpol
    vmec_input.ntheta = 2*mpol + 6
    vmec_input.phiedge = data.phiedge
    vmec_input.curtor = data.curtor
    vmec_input.ncurr = 0
    for name in ('am', 'am_aux_s', 'am_aux_f', 'ai', 'ai_aux_s', 'ai_aux_f'):
        setattr(vmec_input, name, getattr(data, name))
    vmec_input.raxis_cc = np.array([data.raxis])
    vmec_input.raxis_cs = np.zeros(1)
    vmec_input.zaxis_cs = np.zeros(1)
    vmec_input.zaxis_cc = np.array([data.zaxis if lasym else 0.0])
    for name, source in [('rbc','refou'), ('zbs','zefou'), ('rbs','refou2'), ('zbc','zefou2')]:
        setattr(vmec_input, name, getattr(data, source).T)
    vmec_input.comment = (f'Input COCOS={data.cocos}; NCURR=0: current is not independently constrained. '
                          f'LASYM={bool(lasym)}; removed asymmetry max={data.symmetry_removed_max_m:.6g} m.')
    output = Path(output_path) if output_path is not None else Path(f'input.{Path(filename).stem}')
    write_vmec_input(output, vmec_input)
    if plot_path is not None or show:
        from .comparison import plot_comparison
        plot_comparison(gdata, output, data, plot_path=plot_path, show=show)
    return output
