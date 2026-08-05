# Planar Capillary-Wave Validation

## Purpose

This benchmark validates the dynamic coupling between the hydrodynamic solver, the conservative Allen–Cahn phase equation, and the continuum-surface-force (CSF) model. It complements the static Laplace test: Laplace verifies the equilibrium pressure jump, whereas the planar wave verifies the time-dependent capillary response.

The setup is a canonical equal-density capillary-wave problem. A diffuse liquid slab has two planar interfaces displaced by the same sinusoidal mode. The resulting sinuous mode oscillates under surface tension.

## Canonical configuration

The tracked input file is `capillary_wave.inp`:

| Parameter | Value |
|---|---:|
| Domain | `96 x 288 x 8` |
| Wave mode | `1` |
| Initial amplitude | `0.96` |
| Slab half-thickness | `72` |
| Interface width | `4` |
| Surface tension | `0.001` |
| Phase mobility | `0.05` |
| Kinematic viscosity | `0.0025` in both phases |
| Density | `1` in both phases |
| Diagnostic interval | `20` steps |

For this benchmark only, the generic input variables `uwall` and `radius` represent the initial wave amplitude and slab half-thickness, respectively.

The initial phase field is

```text
eta(x) = A cos(k x),                 k = 2 pi / Lx
phi(x,y) = 0.5 [1 + tanh(2 (h - |y-yc-eta(x)|) / W)].
```

For two equal-density fluids, the finite-depth inviscid angular frequency used by the diagnostic is

```text
omega_theory = sqrt[sigma k^3 tanh(k h) / (rho_1 + rho_2)].
```

## Automated build and run

From the repository root, run:

```bash
python3 analyse_capillary_wave.py
```

The script performs the complete workflow:

1. configures `defines.h` for `CAPILLARYWAVE`, `CSF`, `TWOCOMPONENT`, `DENSRATIO`, `D3Q27`, `HIGHORDER`, `PRINTPHI`, and `DOBENCHMARK`;
2. disables other benchmark selectors;
3. runs `make clean`;
4. builds with `make nvfortran GPUCC=80`;
5. runs the simulation in `capillary_wave_run/`;
6. measures the numerical frequency from interpolated zero crossings;
7. checks mass conservation, phase-field kinematics, and the initial CSF Fourier projection;
8. writes a one-row CSV summary.

The script intentionally leaves `defines.h` configured for this benchmark. This makes the compiled executable and the visible source configuration consistent.

Useful alternatives are:

```bash
# Reuse an executable that is already compiled with the correct macros.
python3 analyse_capillary_wave.py --skip-build

# Analyse an existing output directory without compiling or running.
python3 analyse_capillary_wave.py --skip-build --skip-run --output capillary_wave_run

# Select another NVIDIA compute capability and replace existing results.
python3 analyse_capillary_wave.py --gpu-cc 90 --force
```

## Outputs

The output directory contains:

- `capillary_wave.dat`: instantaneous amplitude, mass, theoretical frequency, amplitude derivative, and interface-normal velocity mode;
- `capillary_force_projection.dat`: initial numerical and theoretical CSF Fourier coefficients;
- `run.log`: solver output;
- `capillary_wave_summary.csv`: final machine-readable comparison.

Important CSV fields include `omega_theory`, `omega_numerical`, `frequency_error_percent`, `max_mass_error`, `kinematic_gain_deta_dt_over_vn`, `kinematic_relative_rms_error`, and `force_ratio`.

## Validation result

For the canonical case:

| Quantity | Result |
|---|---:|
| Theoretical angular frequency | `3.7438023e-4` |
| LBFAST angular frequency | `3.4236607e-4` |
| Frequency error | `-8.55%` |
| `d eta/dt` to normal-velocity gain | `0.99743` |
| Kinematic relative RMS error | `0.516%` |
| Numerical/theoretical CSF projection | `0.99704` |
| Maximum relative mass error | `8.7e-10` |

The force and kinematic checks show that the remaining frequency error is a numerical error of this diffuse-interface configuration rather than a loss of effective surface tension or an incorrect phase-advection speed.

## Phase-field kinematic check

The conservative Allen–Cahn equation contains the advection term

```text
partial_t phi + u dot grad(phi) = diffusion + interface compression.
```

The diagnostic compares the measured amplitude derivative with the Fourier mode of the normal fluid velocity at the interface. Their least-squares gain should remain close to one, providing a direct regression check of phase-field advection.
