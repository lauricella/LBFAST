# LBFAST v1.1.0: Changes Since v1.0.3

This document summarizes the changes included in LBFAST `v1.1.0`, released on
17 October 2026, relative to the `v1.0.3` release.

## Release scope

LBFAST `v1.1.0` replaces the fused collide-and-stream kernels with new kernels
that give the same results to round-off and are 2.3 to 3.3 times faster for
D3Q15, D3Q19, and D3Q27, and 2.8 to 4.5 times faster for high-order D3Q27, on
an NVIDIA A30. The new kernels require more shared memory per thread block. In
single-component flows without internal obstacles, the fused kernels also apply
the force update of the moments, so the separate moment kernel is no longer
called. In two-component flows, the interface kernels use the thread
organization of the fused kernels. These two changes give bitwise identical
results in double precision. On an A30, the whole time step is 2.0 to 3.8 times
faster than in `v1.0.3`, depending on the case. The release also adds a
single-GPU Taylor--Green performance case, translates the remaining Italian
comments into English, and updates the documentation. The comparison baseline
is tag `v1.0.3` at commit `4f3d255`.

## Fused collide-and-stream kernels

At every time step, one fused kernel performs the collision, the streaming,
and the computation of the new moments. The kernels of `v1.0.3` used blocks of
`10 x 10 x 10` threads: every thread, including the halo threads,
reconstructed all the populations, and each direction needed two barriers. The
new kernels use one thread per interior node of an `8 x 8 x 8` tile; the first
488 threads also process one halo node, so each node of the tile and its halo
is processed once. The values needed by the neighbors are stored in dynamic
shared memory and gathered (pull) after a barrier.

| Lattice | Kernels | Values per node in shared memory | Barriers |
|---|---|---:|---:|
| D3Q15, D3Q19, D3Q27 | `fused_LB_kernel`, `fused_LB_kernel_bb` | 10 | 1 |
| High-order D3Q27 | `fused_LB_kernel_h`, `fused_LB_kernel_bb_h` | 19 | 3 |

- **Second-order lattices.** The post-collision populations are a linear
  combination of ten post-collision Hermite coefficients, which include the
  equilibrium, the relaxed non-equilibrium part, and half of the Guo forcing
  term. The kernels store the ten coefficients and gather them with stages
  written out for each lattice, with literal components and weights.
- **High-order D3Q27.** The kernels first gather the ten second-order
  coefficients, then compute the seventeen coefficients of orders three to six
  once per node and gather them linearly, with the elements of `Minv_mrt` as
  exact fractions.
- The gather stages are in `fused2_stages.h` and `fusedh_stages.h`, generated
  by `gen_fused2_stages.py` from the lattice definitions in `vars_module.f90`.
- With `BOUNCE_BACK`, a value pulled from a solid node is replaced by the value
  of the node itself, with the sign of the reflected direction.
- `fused_LB_cuda` selects the kernel from the lattice and from `BOUNCE_BACK`.
  No macro is needed, and the kernels of `v1.0.3` have been removed.

The formulation, the implementation, and the verification are described in
[`FUSED_KERNELS.md`](FUSED_KERNELS.md).

## Moment update in single-component flows

Before every collision, `moments_LB_kernel` adds half the force per unit
density to the velocity and the force term to the second-order moment of every
fluid node. In single-component flows without internal obstacles the force is
the uniform body force, and the update depends only on the moments of the node.

- `defines.h` defines `FUSED_MOMENTS` when `TWOCOMPONENT`,
  `INTERNAL_OBSTACLES`, `MULTIHIT`, and `ELASTIC_FORCE` are all undefined. No
  user setting is needed. `POISEUILLE` and `TWOPOISEUILLE` define
  `INTERNAL_OBSTACLES` and keep the separate moment kernel.
- The time loop then skips `moments_LB_cuda` and the exchange of the force
  arrays. The fused kernels apply the update when they read the moments of a
  fluid node, with the same operations and rounding as `moments_LB_kernel`.
- The stored moments are the moments before the update; `copy_print` applies
  the update to the printed velocity.

In double precision the interior moments, the kinetic energy, and the probe
values are bitwise identical to those obtained with the separate moment kernel,
on one and two MPI ranks, with `MIXEDPRC`, in single precision, and with the
`BOUNCE_BACK` kernels. On the `256^3` Taylor--Green vortex, the time step is
1.51 times faster for D3Q19, 1.34 to 1.36 times for D3Q15 and D3Q27, and 1.14 to
1.17 times for high-order D3Q27.

## Interface kernels

In two-component flows, `compute_norm_interface_kernel` (with its interior and
face variants), `compute_div_theta_n_kernel`, and `compute_csf_force_kernel`
compute finite-difference stencils of the phase field and of the interface
normal before the moment update. They used blocks of `10 x 10 x 10` threads, of
which 488 stopped after loading the halo. They now use one thread per interior
node of an `8 x 8 x 8` tile, as the fused kernels; each thread loads up to two
nodes of the tile and its halo. With `8 x 8 x 8` tiles, the divergence and CSF
kernels are compiled with launch bounds for four blocks per multiprocessor.

The arithmetic and the order of the operations in the time step are unchanged,
and all moments are bitwise identical in double precision, with and without
`ASYNCMPI`, on one rank and on two ranks split along each axis. On the `256^3`
two-component droplet, these kernels take 7.6 ms per step instead of 10.9 ms,
and the time step is 1.07 times faster with high-order D3Q27 and 1.12 times
faster with D3Q19.

## Hardware requirements

The new kernels need the following shared memory per thread block:

| Kernels | Double precision and `MIXEDPRC` | Single precision |
|---|---:|---:|
| D3Q15, D3Q19, D3Q27 | 80 KB | 40 KB |
| High-order D3Q27 | 152 KB | 76 KB |

At the first launch, `fused_LB_cuda` requests this memory with
`cudaFuncSetAttribute`; if the GPU cannot provide it, the run stops with an
error message. The NVIDIA A30, A100, and H100 provide enough shared memory for
all kernels. A GPU with less than 152 KB of shared memory per block, such as
the V100, cannot run the high-order kernels in double precision. The
high-order kernels require `8 x 8 x 8` tiles, the default.

## Performance

Mean fused-kernel time per step and throughput of the whole time step on one
NVIDIA A30 for `256^3` lattices, averaged over two runs. The single-component
cases are the Taylor--Green vortex of `taylorgreen256.inp` with `nsteps = 500`,
the two-component cases the droplet of `bench.inp` on a `256^3` lattice with
`nsteps = 300`, and the `BOUNCE_BACK` case a force-driven channel:

| Case | Fused kernel, `v1.0.3` (ms) | Fused kernel, `v1.1.0` (ms) | Speed-up | GLUPS, `v1.0.3` | GLUPS, `v1.1.0` | Speed-up of the time step |
|---|---:|---:|---:|---:|---:|---:|
| Single-component D3Q19, FP64 | 22.2 | 7.0 | 3.2 | 0.599 | 1.933 | 3.2 |
| Single-component D3Q27, FP64 | 29.7 | 9.4 | 3.2 | 0.472 | 1.509 | 3.2 |
| Single-component D3Q27, FP32 | 15.1 | 5.8 | 2.6 | 0.886 | 2.298 | 2.6 |
| Two-component D3Q27, FP64 | 35.4 | 10.7 | 3.3 | 0.297 | 0.587 | 2.0 |
| Single-component high-order D3Q27, FP64 | 82.7 | 23.5 | 3.5 | 0.190 | 0.667 | 3.5 |
| Single-component high-order D3Q27, FP32 | 32.1 | 11.8 | 2.7 | 0.467 | 1.257 | 2.7 |
| Two-component high-order D3Q27, FP64 (default configuration) | 108.5 | 24.5 | 4.4 | 0.129 | 0.396 | 3.1 |
| Single-component high-order D3Q27, FP64, `BOUNCE_BACK` | 114.4 | 26.9 | 4.3 | 0.141 | 0.538 | 3.8 |

The single-component cases without `BOUNCE_BACK` include the moment update in
the fused kernels, and the two-component cases include the new interface
kernels. The `BOUNCE_BACK` channel uses `POISEUILLE`, which keeps the separate
moment kernel. More cases, including D3Q15, are in
[`FUSED_KERNELS.md`](FUSED_KERNELS.md) and
[`TAYLOR_GREEN_256_PERFORMANCE.md`](TAYLOR_GREEN_256_PERFORMANCE.md). The
strong- and weak-scaling tables in the README were measured with earlier
versions and have not been repeated.

## Verification

In double precision, the final moments agree with those of the `v1.0.3`
kernels to round-off, with relative differences up to `8e-13`, for D3Q15,
D3Q19, D3Q27, and high-order D3Q27, for single- and two-component flows, with
and without `BOUNCE_BACK`, with `SMAGORINSKI`, and on one and two MPI ranks.
With `MIXEDPRC` the differences are those of the single-precision storage.

The validation cases run at low Mach number, where the high-order terms are
small. A Taylor--Green case with `U0 = 0.15` (Mach 0.26) tests them directly:
the high-order and second-order reconstructions differ there by `1e-2`, while
the new high-order kernels agree with the `v1.0.3` kernels to `1.6e-14`.

In single precision, the gather sums more terms per output moment than the
`v1.0.3` kernels. The deviations of the moment fields from the
double-precision solution are up to 2.6 times those of `v1.0.3`, below `1.5e-4`
in relative terms, and the kinetic energy of the Taylor--Green vortex differs
from the `v1.0.3` single-precision result by `1e-5` to `1e-4`, much less than
the single-precision error itself, about `1.4e-3`.

| Suite | Ranks | Result |
|---|---|---|
| Smoke suite | 1 | `5 passed, 0 failed` |
| Complete suite | 1 | `5 passed, 0 failed` |
| Complete suite | 2 | `5 passed, 0 failed` (decomposition `1x1x2`; `1x2x1` for the capillary wave) |
| `compile.sh` executable matrix | -- | 24 of 24 executables built, no errors |

In the complete suite, which uses high-order D3Q27, the Laplace, Lamb, and
Poiseuille results are identical to those of `v1.0.3` at all printed digits,
and the capillary-wave and Taylor--Green summaries agree to `2e-11` and `1e-9`.
The Taylor--Green and capillary-wave reference CSV files in `docs/` were
regenerated, and the Taylor--Green guide was aligned with them. The moment
update in the fused kernels and the new interface kernels leave all the
summaries of the smoke and complete suites unchanged.

## Other changes

- `taylorgreen256.inp` runs the Taylor--Green vortex at the Reynolds number of
  `taylorgreen.inp` on a `256^3` lattice that fits on one GPU.
  [`TAYLOR_GREEN_256_PERFORMANCE.md`](TAYLOR_GREEN_256_PERFORMANCE.md) describes
  the build, the timing outputs, and reference A30 results.
- The comments and messages in the input files, sources, and job scripts are
  now in English, including the message printed when a restart file is read.
- The README states how an AI agent was used, and describes the fused kernels.
- [`SCALING_BENCHMARKS.md`](SCALING_BENCHMARKS.md) and
  [`COMPILE_SCRIPT.md`](COMPILE_SCRIPT.md) state which kernels the executables
  use and that the published scaling figures predate them.
- `compile.sh` is tracked as an executable file, so it can be run as
  `./compile.sh`, as described in [`COMPILE_SCRIPT.md`](COMPILE_SCRIPT.md).
- Comments in English were added throughout the source. They explain the tiled
  data layout, the time-step sequence, the force computation, the phase-field
  update, the boundary conditions, the halo exchange, the output and restart
  files, and the validation drivers. No code line was changed.

## Upgrade notes

Users updating from `v1.0.3` should review the following points:

1. Rebuild all executables, including those produced by `compile.sh`.
2. Check that the GPU provides the shared memory listed above; otherwise the
   run stops at the first time step.
3. Input files, `defines.h` macros, Make targets, and output and restart
   formats are unchanged. `FUSED_MOMENTS` is set automatically in `defines.h`
   and should not be defined by hand.
4. Results agree with `v1.0.3` to round-off in double and mixed precision; in
   single precision the round-off is slightly larger.

## Release checklist

Before publishing the `v1.1.0` tag:

- run the smoke suite from a clean build;
- run the complete suite on at least one supported NVIDIA GPU;
- run at least one two- or four-rank MPI validation;
- build the complete `compile.sh` executable matrix;
- confirm that `defines.h` is restored to the documented default;
- confirm that this document matches the final release commit;
- create the `v1.1.0` tag from the validated `main` commit;
- archive the release and record its DOI.
