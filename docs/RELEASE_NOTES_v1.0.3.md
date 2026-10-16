# LBFAST v1.0.3: Changes Since v1.0.2

This document summarizes the changes included in LBFAST `v1.0.3`, released on
16 October 2026, relative to the `v1.0.2` release.

## Release scope

LBFAST `v1.0.3` corrects the off-equilibrium Hermite coefficients of orders
four to six in the high-order D3Q27 reconstruction, switches the fused
collide-and-stream kernels from push to pull streaming, adds the citation of
the LBFAST publication and a guide for reproducing the scaling benchmarks, and
regenerates the reference validation results. The comparison baseline is tag
`v1.0.2` at commit `5047275`.

## High-order D3Q27 off-equilibrium coefficients

This change affects builds with `#define LATTICE 27` and `#define HIGHORDER`,
which is the repository default. D3Q15, D3Q19, and lower-order D3Q27 builds
are not affected.

The high-order reconstruction computes the off-equilibrium Hermite
coefficients of orders three to six from the second-order coefficient
`a1^(2)` and the velocity through the recursive relation of Malaspinas
(arXiv:1505.06900, Eq. 31 and Appendix B):

```text
a1^(n)_{k1...kn} = a1^(n-1)_{k1...k(n-1)} u_kn
                 + ( u_k1 ... u_k(n-2) a1^(2)_{k(n-1) kn} + perm(kn) ),
```

where the bracket is summed over the `n - 1` index pairs that contain the last
index `kn`. In `v1.0.2`, the coefficients of orders four to six summed the
bracket over all index pairs, so the pairs that do not contain `kn` were
counted twice. For example, the `xxyy` coefficient was

```text
v1.0.2:  u^2 a1_yy + 6 u v a1_xy + 2 v^2 a1_xx
v1.0.3:  u^2 a1_yy + 4 u v a1_xy +   v^2 a1_xx
```

The ten off-equilibrium coefficients of orders four to six (`xxyy`, `xxzz`,
`yyzz`, `xxyz`, `xyyz`, `xyzz`, `xxyyz`, `xxyzz`, `xyyzz`, and `xxyyzz`) are
corrected in the following routines:

- `fused_LB_kernel` and `fused_LB_kernel_bb` in `lb_cuda_fused_module.f90`;
- `bcs_mesoscopic_hfields` in `boundary_cds_module.f90`.

The Hermite basis, the equilibrium coefficients, and the third-order
off-equilibrium coefficients were already correct and are unchanged. The
corrected expressions agree with the explicit forms of Malaspinas, Appendix B,
and with an independent Chapman--Enskog projection.

The error removed by this correction scales as `u^(n-2) a1^(2)`, so it is small
at low Mach number. In the complete validation suite, the Laplace, Lamb, and
Poiseuille results are identical to those of `v1.0.2` at all printed digits,
while the capillary-wave frequency and the Taylor--Green viscosity change by
`2e-10` and `4e-10` in relative terms. High-order D3Q27 simulations with larger
velocities or stronger non-equilibrium can differ more from `v1.0.2`.

## Pull streaming in the fused kernels

The fused collide-and-stream kernels now use pull streaming in shared memory.
In `v1.0.2`, each thread reconstructed a post-collision population and wrote
it to the shifted position of the shared-memory buffer (push). In `v1.0.3`,
each thread writes the population at its own position; after the barrier,
each interior node of the tile gathers the population arriving from
`x - c_i` and adds it to the moment accumulators. This is the pull
formulation described in the LBFAST publication.

In `fused_LB_kernel_bb`, used when `BOUNCE_BACK` is defined, the bounce-back
reflection is applied during the gather: when the upstream node is solid, the
opposite population of the same node is used. This removes the separate
write-back step and one barrier for each pair of opposite directions.

The change does not alter the arithmetic. In double precision, the final
moment fields are bitwise identical to those of the push kernel with the same
coefficients for the high-order D3Q27 Taylor--Green vortex on one and two MPI
ranks, for the D3Q19 Taylor--Green vortex, for the two-component Lamb droplet,
and for a Poiseuille flow with `BOUNCE_BACK`. All quantities reported by the
complete validation suite are also identical.

The table reports the mean fused-kernel time per step on one NVIDIA A30 for
`256^3` lattices, averaged over two runs. The one-component cases are a
Taylor--Green vortex at `Re = 80`, except the `BOUNCE_BACK` case, which is a
force-driven channel flow with walls normal to `x`. The two-component case is
the static droplet of the default configuration.

| Case | Lattice and precision | Push (ms) | Pull (ms) | Change |
|---|---|---:|---:|---:|
| One-component | D3Q19, FP64 | 23.3 | 22.2 | -4.6% |
| One-component | D3Q27, FP64 | 31.2 | 29.8 | -4.6% |
| One-component | High-order D3Q27, FP64 | 82.0 | 83.4 | +1.7% |
| One-component | High-order D3Q27, FP32 | 34.9 | 32.1 | -7.9% |
| Two-component | High-order D3Q27, FP64 | 108.3 | 109.1 | +0.7% |
| One-component, `BOUNCE_BACK` | High-order D3Q27, FP64 | 119.0 | 115.0 | -3.3% |

The strong- and weak-scaling tables in the README were measured before this
release and have not been repeated.

## Citation

The README now cites the LBFAST publication, with a BibTeX entry:

> Lauricella, M., Montessori, A., Amati, G., Spiga, F., Tiribocchi, A.,
> Bernaschi, M., & Succi, S. (2026). *LBFAST: A Lightweight Moment-Represented
> Lattice Boltzmann Solver for Multi-GPU Architectures*. Procedia Computer
> Science, 286, 34--47. https://doi.org/10.1016/j.procs.2026.08.017

`CITATION.cff` lists the same article as `preferred-citation`.

## Documentation and reference results

- The new guide [`SCALING_BENCHMARKS.md`](SCALING_BENCHMARKS.md) describes the
  inputs, executables, job scripts, MPI decompositions, initial conditions, and
  log metrics used for the strong-scaling, weak-scaling, and energy
  measurements.
- The README description of the weak-scaling mode is corrected: each rank
  receives a cubic subdomain of side `MYSIDE`, a compile-time macro with
  default value `512`, and the domain size in the input file is overridden.
- The README now states that `CSF` is enabled by default, lists D3Q15 among the
  supported lattices, and documents the `CAPILLARYWAVE` and
  `POISEUILLESTARTREST` selectors.
- The reference CSV files in `docs/` are regenerated from the complete
  validation suite of this release, and the tables of the validation guides
  are aligned with them. In `v1.0.2`, the Laplace CSV had been produced with
  an earlier state of the code, and the capillary-wave CSV with a shorter run
  (`1332` samples instead of `4001`). The GPU memory columns are now taken from
  the device memory reported in the run log.

## Validation status

The following runs were performed on NVIDIA A30 GPUs with the release code:

| Suite | Ranks | Result |
|---|---|---|
| Smoke suite | 1 | `5 passed, 0 failed` |
| Complete suite | 1 | `5 passed, 0 failed` |
| Complete suite | 2 | `5 passed, 0 failed` (decomposition `1x1x2`; `1x2x1` for the capillary wave) |

Complete-suite results on one GPU:

| Benchmark | Result |
|---|---|
| Static Laplace sweep | PASS; slope error `1.179%`, maximum point error `1.597%` |
| Planar capillary wave | PASS; frequency error `8.464%`, force ratio `0.997039` |
| Lamb droplet oscillation | PASS; period error `6.478%` |
| Taylor--Green vortex | PASS; viscosity error `1.650%` |
| Force-driven Poiseuille flow | PASS; relative L2 error `0.340%`, centerline error `0.254%` |

## Upgrade notes

Users updating from `v1.0.2` should review the following points:

1. High-order D3Q27 results change. At low Mach number the change is
   negligible, but runs that are compared with `v1.0.2` at higher velocities
   should be repeated.
2. Input files, `defines.h` macros, Make targets, and output and restart
   formats are unchanged. Recompile after updating.
3. Pull streaming replaces push streaming and cannot be disabled; it does not
   change the results.
4. Executables built from earlier versions, including those produced by
   `compile.sh`, must be rebuilt to include the corrections.

## Release checklist

Before publishing the `v1.0.3` tag:

- run the smoke suite from a clean build;
- run the complete suite on at least one supported NVIDIA GPU;
- run at least one two- or four-rank MPI validation;
- confirm that `defines.h` is restored to the documented default;
- confirm that this document matches the final release commit;
- create the `v1.0.3` tag from the validated `main` commit;
- archive the release and record its DOI, if applicable.
