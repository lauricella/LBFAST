# Taylor--Green vortex validation on a 128-cubed lattice

## Purpose

This benchmark validates the effective viscosity of the single-component
LBFAST solver through the viscous decay of a three-dimensional Taylor--Green
vortex. The compact `128 x 128 x 128` configuration is intended for a single
NVIDIA A30 GPU and runs entirely in double precision.

The test is configured independently from the large reference input:
`taylorgreen.inp` remains unchanged, while this benchmark uses
`taylorgreen128.inp`.

## Initial condition

In a fully periodic cubic domain, LBFAST initializes

```text
x = 2 pi (i-1) / Nx,
y = 2 pi (j-1) / Ny,
z = 2 pi (k-1) / Nz,

u_x =  U0 sin(x) cos(y) cos(z),
u_y = -U0 cos(x) sin(y) cos(z),
u_z =  0,

p = (rho U0^2 / 16) [cos(2x) + cos(2y)] [2 + cos(2z)].
```

The implementation uses global periodic coordinates, including under domain
decomposition. The pressure is converted to the pressure-like lattice variable
by division by `rho cs^2`. The velocity field is analytically divergence-free.

## Configuration

The tracked `taylorgreen128.inp` case uses:

| Parameter | Value |
|---|---:|
| Domain | `128 x 128 x 128` |
| Time steps | `2000` |
| Energy sampling interval | `10` steps |
| Density | `1` |
| Initial velocity `U0` | `0.04` |
| Kinematic viscosity | `0.04074` |
| Characteristic length `L=N/(2 pi)` | `20.3718` |
| Reynolds number `U0 L / nu` | `20.00` |
| Lattice | D3Q27 high-order |
| Precision | double precision |

The lower Reynolds number follows from retaining `U0` and the viscosity while
reducing the lattice size. This provides a well-resolved viscous-decay test in
a small GPU-memory footprint.

## Viscous-decay measurement

For this mode, the normalized kinetic energy follows

```text
E(t) / E(0) = exp[-6 nu (2 pi / N)^2 t].
```

The numerical viscosity is therefore obtained from a linear least-squares fit
to

```text
log[E(t) / E(0)] = slope t + intercept,
nu_num = -slope / [6 (2 pi / N)^2].
```

The default analysis uses samples satisfying
`log(E/E0) >= -0.25`. On this grid that gives 41 samples from the 10-step
diagnostic series, ending at step 410. Restricting the fit to the early decay
avoids giving excessive weight to late-time deviations from the analytical
exponential.

## Current result

The validation run on an NVIDIA A30 produced:

| Quantity | Result |
|---|---:|
| Theoretical viscosity | `0.0407400000` |
| Numerical viscosity | `0.04141210127401` |
| Relative viscosity error | `+1.649733%` |
| Fit slope | `-5.987124858286e-4` |
| Fit intercept | `3.943066331990e-4` |
| RMS error in fitted log-energy | `4.832086838738e-4` |
| Fit samples | `41` |
| Last fitted step | `410` |
| Reported GPU allocation | `0.9135 GB` |

The machine-readable comparison is stored in
`docs/taylorgreen128_results.csv`.

## Automated build, run, and analysis

From the repository root, run:

```bash
python3 run_taylorgreen_validation.py
```

The script configures `defines.h` for a single-component, double-precision,
D3Q27 high-order build. It enables `TAYLORGREEN`, explicitly disables the
two-component, density-ratio, mixed-precision, and CSF options, performs a
clean CC 8.0 build, and runs the case in `taylorgreen_validation_run/`.

The run directory contains the copied input, solver log, numerical and
theoretical energy histories, plot, and `taylorgreen_validation_summary.csv`.

Useful alternatives are:

```bash
# Replace an existing run.
python3 run_taylorgreen_validation.py --force

# Reuse an executable compiled with the required macros.
python3 run_taylorgreen_validation.py --skip-build

# Analyse an existing result without compiling or running.
python3 run_taylorgreen_validation.py \
  --skip-build --skip-run \
  --output taylorgreen_validation_run

# Select a different fit threshold or GPU compute capability.
python3 run_taylorgreen_validation.py --max-log-decay -0.5 --gpu-cc 90 --force
```

After configuring and building, the script leaves `defines.h` consistent with
the generated Taylor--Green executable.
