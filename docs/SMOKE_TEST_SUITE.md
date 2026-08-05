# Short regression test suite

## Purpose

The `tests/` directory contains reduced versions of every physical validation
case documented for LBFAST. They provide a quick build-and-run regression
suite for development. These short runs do not replace the quantitative
benchmarks: they check that each compile-time configuration builds, starts,
produces finite diagnostics, and follows the expected initial physical trend.

## Cases and checks

| Test | Reduced domain | Steps | Pass condition |
|---|---:|---:|---|
| CSF Laplace | `32^3`, `R=8` | 20 | Positive pressure jump and effective surface tension; finite measured radius. |
| Planar capillary wave | `32 x 64 x 8` | 20 | Conserved phase mass and CSF Fourier projection within a broad regression bound. |
| Lamb oscillation | `40^3` | 20 | Finite ellipsoidal interface dimensions and bounded deformation. |
| Taylor--Green | `32^3` | 30 | At least three finite samples and decreasing kinetic energy. |
| Forced Poiseuille | `16 x 8 x 16` | 100 | Positive, symmetric force-generated flow and zero velocity on solid planes. |

All dimensions are multiples of the default `8 x 8 x 8` GPU tile. Field,
restart, and visualization output is disabled in the reduced inputs.

## Running the suite

From the repository root:

```bash
python3 tests/run_smoke_tests.py
```

The default compiler target is `nvfortran` for compute capability 8.0. These
can be changed without editing the script:

```bash
python3 tests/run_smoke_tests.py --gpu-cc 90
python3 tests/run_smoke_tests.py --make-target nvfortran --case taylor_green
python3 tests/run_smoke_tests.py --case laplace_csf --case poiseuille_force
```

The cases require different compile-time selectors, so the script performs a
clean build for every selected test. It saves the complete initial
`defines.h`, restores it even when a test fails, and finally runs `make clean`
and removes the generated `main.x`, so that no executable remains inconsistent
with the restored macros.

The console summary has the form

```text
PASS  laplace_csf          positive Laplace jump; sigma_eff=...
PASS  capillary_wave       mass conserved; CSF projection=...
PASS  lamb_oscillation     finite interface; deformation=...
PASS  taylor_green         energy decays; log(E/E0)=...
PASS  poiseuille_force     force generates symmetric flow; umax=...

5 passed, 0 failed
```

Individual inputs, build and run logs, compact diagnostics, and the aggregate
`smoke_test_results.csv` are written below `tests/results/`. This directory is
ignored by Git. Use `--keep-results` to preserve existing case directories
when rerunning a subset.

## Interpretation

The tolerances are intentionally broad enough for short transient runs and
small diffuse-interface geometries. A passing smoke test establishes basic
build and runtime integrity, not quantitative agreement with the reference
period, viscosity, or stationary pressure jump. Use the full drivers linked
from the main README for scientific validation.
