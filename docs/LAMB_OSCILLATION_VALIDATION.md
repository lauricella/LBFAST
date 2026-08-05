# Ellipsoidal-droplet oscillation validation

## Purpose

This benchmark validates the dynamic response of the two-component LBFAST
solver by measuring the quadrupolar oscillation period of a freely evolving
prolate droplet. The numerical period is compared with the viscous
Miller--Scriven prediction used by the built-in `LAMBTEST` diagnostic.

The benchmark exercises the conservative Allen--Cahn phase equation, the
variable-density hydrodynamic model, and the continuum-surface-force (CSF)
coupling. The static surface-tension response is assessed separately by the
Laplace validation.

## Configuration

The tracked `lamb.inp` case uses:

| Parameter | Value |
|---|---:|
| Domain | `80 x 80 x 80` |
| Center | `(40, 40, 40)` |
| Initial semi-axes | `22, 22, 30` |
| Equivalent radius | `24.3962011465` |
| Interface width | `4` |
| Cahn number `W/Re` | `0.1639599533` |
| Surface tension | `2.27e-3` |
| Inner/outer density | `4 / 1` |
| Inner/outer kinematic viscosity | `3.333e-3 / 3.333e-3` |
| Allen--Cahn mobility parameter `tau_diff` | `0.1` |
| Time steps | `20000` |
| Diagnostic interval | `10` steps |

The initialization preserves the reference prolate aspect ratio `11:11:15`.
The `radius` input parameter specifies the equivalent radius

```text
Re = (ax ay az)^(1/3),
```

so the tracked value gives the resolved semi-axes `22,22,30`. Keeping `W=4`
while doubling the droplet resolution gives `W/Re=0.164`.

Full raw-field output is disabled because the period measurement only requires
the compact `lamb.dat` time series.

## Theoretical comparison

For the quadrupolar mode, the inviscid angular frequency is

```text
omega_0 = sqrt[24 sigma / (Re^3 (2 rho_out + 3 rho_in))].
```

The built-in diagnostic applies the Miller--Scriven viscous correction using
the dynamic viscosities of both phases. It reports the corrected angular
frequency and period at startup, locates maxima of the measured axial interface
position, and obtains the numerical period from consecutive maxima.

## Current result

The `80^3` double-precision D3Q27 high-order CSF run produced:

| Quantity | Value |
|---|---:|
| Theoretical period | `12772.6392608318` |
| Numerical period | `13600` |
| Period error | `+6.477602%` |
| Theoretical angular frequency | `4.919254e-4` |
| Numerical angular frequency | `4.619989e-4` |
| Angular-frequency error | `-6.083544%` |
| Detected maxima | steps `90` and `13690` |

This measurement contains one complete peak-to-peak interval. A longer run can
be requested to average several periods. The committed machine-readable result
is in `docs/lamb_oscillation_results.csv`.

## Automated build, run, and analysis

From the repository root, run:

```bash
python3 run_lamb_validation.py
```

The script configures `defines.h` with the baseline double-precision settings
expected by `compile.sh`: `LATTICE=27`, `HIGHORDER`, `TWOCOMPONENT`,
`DENSRATIO`, `PRC=8`, `STRPRC=8`, and `noMIXEDPRC`. It also enables `CSF`,
`DOBENCHMARK`, and `LAMBTEST`, disables the other benchmark selectors, performs
a clean CC 8.0 build, and runs the case in `lamb_validation_run/`.

The output directory contains the copied input, `run.log`, `lamb.dat`, the
theoretical time series, the plot, and `lamb_validation_summary.csv`.

Useful alternatives are:

```bash
# Average more oscillations.
python3 run_lamb_validation.py --nsteps 40000 --force

# Reuse an executable compiled with the required macros.
python3 run_lamb_validation.py --skip-build

# Analyse an existing output directory without compiling or running.
python3 run_lamb_validation.py --skip-build --skip-run --output lamb_validation_run

# Select another NVIDIA compute capability.
python3 run_lamb_validation.py --gpu-cc 90 --force
```

The script intentionally leaves `defines.h` configured for the Lamb benchmark,
so the visible macro configuration remains consistent with the executable.
