# LBFAST v1.0.2: Changes Since v1.0.1

This document summarizes the changes included in LBFAST `v1.0.2`, released on
6 August 2026, relative to the `v1.0.1` release.

## Release scope

LBFAST `v1.0.2` extends the two-component capillary model, adds reproducible
physical validation workflows, improves MPI/OpenACC device initialization, and
documents the supported build and benchmark configurations. The comparison
baseline is tag `v1.0.1` at commit `95f1631`.

## Capillary and phase-field model

### Geometric continuum-surface-force model

A geometric continuum-surface-force (CSF) formulation is now available under
the compile-time macro `CSF`. It computes the interface normal from the phase
field, evaluates curvature as the negative divergence of that normal, and
applies the localized force

```text
f_sigma = sigma kappa n |grad(phi)|.
```

The curvature is evaluated on the GPU using the isotropic D3Q27 stencil. The
normal components are reused from the existing auxiliary fields, while
curvature is computed on the fly in shared memory; no persistent curvature
array is allocated. The force is stored in the existing force array and is
consumed by every active moment kernel, including the interior and all MPI halo
variants.

When `CSF` is disabled, the existing Jacqmin chemical-potential formulation
remains available. The default `defines.h` configuration now enables `CSF`.

### Interface transport and density-ratio correction

The conservative phase-field advection is applied consistently in all eight
update kernels: the full-domain kernel, the interior kernel, and the six halo
face kernels. The density-interface correction profile used by the moment
kernels is now

```text
sqrt(max(0, 4 phi (1 - phi))).
```

The same expression is present in the full-domain, interior, and six halo
moment kernels.

## New physical validation cases

Five documented quantitative benchmarks are included.

| Benchmark | Main purpose | Current full-suite result |
|---|---|---|
| Static Laplace sweep | Effective surface tension and `Delta p proportional to 1/R` | PASS; slope error `1.179%`, maximum point error `1.597%` |
| Planar capillary wave | Dynamic capillary response and force projection | PASS; frequency error `8.464%`, force ratio `0.997039` |
| Lamb droplet oscillation | Miller--Scriven oscillation period | PASS; period error `6.478%` |
| Taylor--Green vortex | Single-component viscous decay in FP64 | PASS; viscosity error `1.650%` |
| Force-driven Poiseuille flow | Single-component body force and wall placement | PASS; relative L2 error `0.340%`, centerline error `0.254%` |

Each benchmark has a dedicated input file, Python driver, Markdown guide, and
machine-readable CSV reference result. The Python drivers configure
`defines.h`, perform a clean build, run in an isolated output directory, apply
the quantitative comparison, and restore the original build configuration.

The new benchmark macros are:

- `CAPILLARYWAVE` for the planar capillary wave;
- `LAMBTEST` for the oscillating ellipsoidal droplet;
- `POISEUILLESTARTREST` to initialize the forced Poiseuille test from rest.

The existing `LAPLACE`, `TAYLORGREEN`, and `POISEUILLE` selectors are used by
their corresponding validation drivers.

## Automated test suites

Two repository-level test suites are now provided:

```bash
python3 tests/run_smoke_tests.py
python3 tests/run_full_validation.py
```

The smoke suite performs short build-and-run regression checks for all five
documented cases. The full suite runs their quantitative production
configurations and writes a combined CSV report to
`tests/full_results/full_validation_results.csv`.

Both suites accept `--mpi-procs 2` and `--mpi-procs 4`. For MPI runs, the
drivers select a tile-compatible decomposition supported by each benchmark.
Incompatible decompositions are rejected before execution.

The complete single-GPU suite was run on physical GPU 3 before preparing this
document and completed with `5 passed, 0 failed`.

## MPI/OpenACC device initialization

GPU selection now occurs immediately after MPI determines the node-local rank
and before the OpenACC runtime queries or allocates on a device. Each local MPI
rank selects

```text
device = node_local_rank modulo visible_device_count.
```

This ordering prevents nonzero ranks from first creating an unused context on
device 0 and then a second context on their assigned GPU. A four-rank test with
four visible GPUs now creates exactly four LBFAST GPU processes, one per GPU.

The non-MPI path explicitly initializes the local rank and local size to zero
and one, respectively, and follows the same device-selection logic.

## Build and default configuration

`compile.sh` now disables `CSF` while producing single-component executables
and restores it when returning to the default two-component configuration.
This keeps all generated executable variants consistent with their intended
physical model.

The repository default is a two-component, high-density-ratio, high-order
D3Q27 build with CSF enabled. Benchmark selectors remain disabled in the base
configuration, so it is compatible with `compile.sh` and with the standard
`test512*.inp` production inputs.

The build matrix and the macro transitions performed by `compile.sh` are
documented in [`COMPILE_SCRIPT.md`](COMPILE_SCRIPT.md).

## Documentation and machine-readable results

The release adds dedicated guides for:

- the CSF Laplace-pressure validation;
- the planar capillary-wave validation;
- the Lamb oscillating-droplet validation;
- the `128^3` Taylor--Green validation;
- the force-driven Poiseuille validation;
- the smoke suite and the complete validation suite;
- operation of `compile.sh`.

Reference results are tracked as CSV files in `docs/`. Runtime logs, plots,
executables, and generated run directories remain untracked artifacts.

## Upgrade notes

Users updating from `v1.0.1` should review the following points:

1. `CSF` is enabled in the default `defines.h`. Define `noCSF` to retain the
   Jacqmin capillary-force path.
2. Recompile after changing any physical-model or benchmark macro.
3. Use a tile-compatible domain and MPI decomposition. The validation drivers
   check this automatically.
4. When selecting a physical GPU with `CUDA_VISIBLE_DEVICES`, the visible
   device is renumbered from zero inside the process. For example,
   `CUDA_VISIBLE_DEVICES=3` exposes physical GPU 3 as OpenACC device 0.
5. Generated benchmark results are written below the selected output directory
   and are not intended to be committed.

## Release checklist

Before publishing the `v1.0.2` tag:

- run the smoke suite from a clean build;
- run the complete suite on at least one supported NVIDIA GPU;
- run at least one two- or four-rank MPI validation;
- confirm that `defines.h` is restored to the documented default;
- confirm that this document matches the final release commit;
- create the `v1.0.2` tag from the validated `main` commit;
- archive the release and record its DOI, if applicable.
