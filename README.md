# GEQDSK to VMEC converter

Convert axisymmetric tokamak GEQDSK equilibria into fixed-boundary VMEC INDATA
namelists. The converter prescribes pressure and rotational transform (`NCURR=0`)
and can compare the **written input** with the GEQDSK data. It never runs VMEC.

## Install and run

Requires Python 3.9+, NumPy, SciPy, Matplotlib, and f90nml:

```bash
pip install -e .
python example_usage.py Device/KSTAR/g02895743.geq \
  --output input.g02895743 --plot comparison.png
```

The omitted coordinate convention produces an explicit warning; specify
`--cocos N` when the producer's COCOS index is known. Use `--show` to display
interactively or `MPLBACKEND=Agg` for headless plotting.

```python
from eqdsk2vmec import convert_eqdsk_to_vmec

output = convert_eqdsk_to_vmec(
    'Device/KSTAR/g02895743.geq',
    output_path='input.g02895743',
    plot_path='comparison.png',
    mpol=12,
    lasym=False,
    # cocos=5,  # Set only after checking the source convention.
)
```

The function returns the output `Path`. Plotting is optional and no GUI opens by
default. A plot also writes an adjacent JSON file with comparison metrics.

## Geometry and symmetry

`NTOR=0` expresses tokamak axisymmetry. `LASYM=False` remains the default, explicitly
projecting the boundary onto an **up–down symmetric shape about Z=0**. Axisymmetry
alone does not imply up–down symmetry. The converter warns when discarded
asymmetric components exceed 1 mm. Use `lasym=True` / `--lasym` to preserve them.
The symmetric initial axis has Z=0; the asymmetric case retains the source Z axis.

Boundary points are interpreted as an ordered polar graph about (R_axis, 0) in
symmetric mode or (R_axis, Z_axis) in asymmetric mode. Non-star-shaped boundaries
are rejected. Nonuniform samples are periodically resampled in geometric angle
before the Fourier transform. `MPOL=12` retains precisely m=0 through 11 in both
fitting and writing; increase `--mpol` to investigate truncation error, subject to
sufficient distinct input points. Sharp separatrix tips cannot be represented
exactly by a finite smooth Fourier series.

## Profiles and coordinate conventions

Toroidal flux is integrated analytically from a PCHIP representation of q on
normalized poloidal flux. Both increasing and decreasing physical poloidal flux
are supported. Zero flux span, zero/sign-changing q, nonfinite data, and malformed
files are rejected rather than silently propagated.

COCOS indices 1–8 and 11–18 are supported for **signed q**. Conversion accounts for
q handedness, the toroidal direction of B/current, and whether psi is in Wb or
Wb/rad. With `cocos=None`, the legacy reciprocal-q convention is made explicit as
an assumption of COCOS 5, with a warning and a comment in the namelist. GEQDSK
files do not reliably encode COCOS; the supplied fixtures' conventions have not
been independently established. For sources that replace signed q with |q|,
restore its convention-correct sign before conversion; the code cannot infer it
from positive q alone.

Pressure and iota are written as Akima spline knots at normalized toroidal flux.
At most 101 original profile samples are retained to respect VMEC2000's
`ndatafmax=101`. For larger profiles, knots are selected by repeatedly adding the
point with the largest scaled pressure/iota linear-interpolation error. Axis and
edge are always retained. Polynomial coefficients remain available with both
endpoints constrained, but are inactive under the default spline profile types.

The invalid sum `p' + FF'` has been removed. No toroidal current profile is inferred.
`CURTOR` records the convention-adjusted source total current, but with `NCURR=0`
it is **not an independent constraint**. A caller selecting `NCURR=1` in the
low-level writer must supply a physically derived current profile.

## Comparison plots

See [generated comparisons](artifacts/comparisons/README.md).
The plot reads the namelist back with f90nml and reconstructs only the modes that
will be retained, honoring LASYM and MPOL. It overlays the GEQDSK boundary and
pressure/q samples with the written boundary and spline knots, and shows the
boundary distance errors. The gray interior contours come **only from GEQDSK**;
there is no converted interior equilibrium until VMEC is solved separately.

Profile lines connect knots; they are not VMEC-evaluated splines. JSON metrics
separate serialization errors at retained knots from linear-interpolation
errors evaluated at all source samples. Boundary distances are one-way nearest
Euclidean distances from original points to an 8192-point reconstructed curve,
not a Hausdorff distance or a force-balance error.

## Tests

```bash
MPLBACKEND=Agg python -m unittest discover -v
```

Tests cover analytic flux integration, decreasing flux, COCOS transformations,
polynomial endpoints, nonuniform boundary sampling, asymmetric geometry, parsed
namelist coefficients, malformed files, and both supplied equilibria. These tests
validate conversion and serialization, not VMEC convergence or force balance.

## References

- [VMEC input parameters](https://princetonuniversity.github.io/STELLOPT/VMEC%20Input%20Namelist%20(v8.47).html)
- [FreeQDSK format and conventions](https://freeqdsk.readthedocs.io/en/stable/geqdsk.html)
- [Sauter and Medvedev: COCOS, Table I](https://www.epfl.ch/research/domains/swiss-plasma-center/wp-content/uploads/2018/10/Sauter_COCOS_Tokamak_Coordinate_Conventions.pdf)
- [VMEC sign conventions](https://terpconnect.umd.edu/~mattland/assets/notes/vmec_signs.pdf)
- [VMEC2000 array limits](https://github.com/PrincetonUniversity/STELLOPT/blob/develop/LIBSTELL/Sources/Modules/vparams.f)
