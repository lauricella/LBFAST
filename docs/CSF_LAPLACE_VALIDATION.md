# CSF Laplace validation

## Purpose

This document describes the implementation and validation of the geometric
continuum-surface-force (CSF) model in the two-component LBFAST solver. The
model is selected at compile time with `CSF` in `defines.h` and is validated
against the three-dimensional Laplace pressure jump.

## CSF formulation

LBFAST uses `phi=1` inside the droplet. The phase-field gradient therefore
points inward. The regularized normal, curvature, and capillary force are

```text
n         = grad(phi) / sqrt(|grad(phi)|^2 + epsilon_n^2),
curvature = -div(n),
F_sigma   = sigma curvature grad(phi),
```

where `epsilon_n=1e-12` in lattice units. This sign convention produces a
positive pressure jump inside a spherical droplet.

Both the phase gradient and the divergence of the normal use the isotropic
D3Q27 26-neighbour stencil. The normal halo exchange already used by the
conservative Allen-Cahn solver is completed before curvature is evaluated, so
the CSF geometry is consistent across tiled and MPI boundaries.

### Memory layout

Curvature is not stored in a persistent auxiliary field. A dedicated CUDA
kernel evaluates `-div(n)` in registers, immediately forms the three components
of `sigma curvature grad(phi)`, and writes them into the existing `forces`
array. The eight hydrodynamic-moment kernels then read this local CSF force,
add the remaining body-force, pressure, viscosity, and density-ratio terms,
and overwrite `forces` with the total force used by the hydrodynamic update.

Consequently, the CSF implementation does not increase `nlocauxfields` or the
persistent GPU-memory footprint. For the `64^3` validation case the reported
GPU allocation remained `0.4975 GB`. The CSF kernel used 24 kB of shared
memory, 54 registers, and no register spills when compiled for NVIDIA CC 8.0.

## Static Laplace test

The baseline validation input is `laplace.inp`:

| Parameter | Value |
|---|---:|
| Domain | `64 x 64 x 64` |
| Radius | `16` |
| Center | `(32, 32, 32)` |
| Interface width | `4` |
| Surface tension | `0.03` |
| Kinematic viscosities | `0.1`, `0.1` |
| Densities | `1`, `1` |
| Time steps | `10000` |
| Sampling interval | `100` steps |

For a three-dimensional spherical droplet, the theoretical pressure jump is

```text
Delta p = 2 sigma / R = 0.00375.
```

LBFAST stores a pressure-like lattice variable. The Laplace diagnostic converts
it to physical pressure according to

```text
p = rho(phi) cs^2 pstar,
```

measures the internal value at the droplet center and the external value at
`(1,32,32)`, and computes

```text
sigma_eff = Delta p R_measured / 2.
```

The initial condition preloads the theoretical spherical pressure jump. The
step-zero sample is therefore not an independent capillary validation. Results
below use the stationary window from step 5000 through step 10000.

### Single-radius result

| `Delta p` | `sigma_eff` | Relative surface-tension error |
|---:|---:|---:|
| `0.00380989` | `0.0305051` | `+1.68%` |

The result is stationary by step 5000; averages over steps 5000--10000 and
7500--10000 agree to much better than the measured model error.

## Radius sweep

The script `run_laplace_csf_sweep.py` generates and runs three isolated cases
with constant interface width `W=4` and constant relative box size `N=4R`:

| Radius | Domain | Numerical `Delta p` | Exact `Delta p` | Relative error |
|---:|---:|---:|---:|---:|
| 16 | `64^3` | `3.809886420e-3` | `3.750000000e-3` | `+1.5970%` |
| 24 | `96^3` | `2.517192400e-3` | `2.500000000e-3` | `+0.6877%` |
| 32 | `128^3` | `1.882137790e-3` | `1.875000000e-3` | `+0.3807%` |

Each value is averaged over steps 5000--10000. A least-squares fit through the
origin to

```text
Delta p = slope / R
```

gives

```text
numerical slope = 6.070741217821e-2,
exact slope     = 6.000000000000e-2,
slope error     = +1.179020%,
sigma from fit  = 3.035370608910e-2.
```

The decreasing pointwise error demonstrates convergence toward the Laplace law
as curvature is reduced. The through-origin fit provides the effective surface
tension independently of any single-radius measurement.

## Reproduction

By default, the sweep script configures the required CSF Laplace macros in
`defines.h`, performs a clean build equivalent to

```bash
make clean
make nvfortran GPUCC=80
```

Run the baseline case with

```bash
./main.x laplace.inp 2>&1 | tee laplace_csf_run.log
```

Run or repeat the radius sweep with

```bash
python3 run_laplace_csf_sweep.py
python3 run_laplace_csf_sweep.py --force
```

The GPU compute capability can be selected, for example, with `--gpu-cc 90`.
An already compiled executable can be reused explicitly with `--skip-build`.
In that mode the user is responsible for ensuring that it was built with the
correct `defines.h` configuration.

The script is compatible with Python 3.6. It creates one directory per case
under `laplace_csf_sweep/`, preserves each generated input, `laplace.dat`, and
run log, and writes the aggregate files

```text
laplace_csf_sweep/laplace_csf_sweep_summary.csv
laplace_csf_sweep/laplace_csf_fit.txt
```

Generated binaries, object files, logs, and numerical output directories are
validation artifacts and are not part of the source commit.
