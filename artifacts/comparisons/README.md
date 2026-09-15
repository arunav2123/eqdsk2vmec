# GEQDSK versus written VMEC inputs

These plots were generated without running VMEC. Main files retain the requested
LASYM=False default; `_asymmetric` files show LASYM=True for comparison. Both use
MPOL=12 and assume COCOS=5 because the fixture conventions are undocumented.
The original GEQDSK fixtures under Device are retained. Obsolete pre-fix VMEC
inputs were removed; the current converted inputs are stored alongside these plots.

| Fixture | LASYM | Boundary RMS [mm] | Boundary max [mm] | Plot |
|---|---|---:|---:|---|
| Jacob-DIII-D | False | 32.81 | 166.50 | [g147634.png](g147634.png) |
| Jacob-DIII-D | True | 15.21 | 102.97 | [g147634_asymmetric.png](g147634_asymmetric.png) |
| KSTAR | False | 23.18 | 102.25 | [g02895743.png](g02895743.png) |
| KSTAR | True | 7.46 | 51.08 | [g02895743_asymmetric.png](g02895743_asymmetric.png) |

Pressure and q are compared in the source convention. Their retained-knot errors
are serialization errors; KSTAR is reduced from 129 to 101 knots, while DIII-D
retains all 65. Adjacent JSON files include linear-interpolation diagnostic errors
at every source sample. These are not VMEC spline or equilibrium errors.

The boundary metric measures each original point's distance to the densely sampled
written curve. Sharp lower tips retain Fourier truncation error even for LASYM=True.
LASYM=False additionally forces up–down symmetry about Z=0, altering the source shape.
Increasing MPOL cannot restore symmetry components omitted by LASYM=False.

Reproduce from the repository root after installing the package:

```bash
MPLBACKEND=Agg python example_usage.py Device/KSTAR/g02895743.geq \
  --output artifacts/comparisons/input.g02895743 \
  --plot artifacts/comparisons/g02895743.png
```

For the second fixture, substitute `Device/Jacob-DIII-D/g147634.geq` and the
corresponding output names. Add `--lasym` and an `_asymmetric` suffix to reproduce
the comparison variants. No `--cocos` is supplied intentionally: the warning records
that COCOS 5 is an assumption, not verified fixture metadata.
