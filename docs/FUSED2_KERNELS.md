# Second-order fused kernels (`FUSED2`)

## Purpose

For D3Q19 and D3Q27 without `HIGHORDER`, LBFAST performs the fused
collide-and-stream step with the kernels `fused_LB_kernel_2` and
`fused_LB_kernel_bb_2`, selected by the `FUSED2` macro. They give the same
results as the general kernels `fused_LB_kernel` and `fused_LB_kernel_bb`, to
round-off, and run about twice as fast. This document describes their
formulation, implementation, verification, and performance.

## Formulation

In the second-order model, the post-collision population of direction `l` at
node `x` is a linear combination of ten post-collision Hermite coefficients:

```text
f*_l(x) = p_l [ A0 + c_l.A1 / cs^2 + H2(c_l):A2 / (2 cs^4) ],

A0 = p*,
A1 = u + F/2,
A2 = u u + (1 - omega) Pi_neq + (u F + F u)/2,
```

where `p_l` are the lattice weights, `Pi_neq = a2 - u u` is the non-equilibrium
part of the second-order coefficient, `F` is the force per unit density, and
`omega` is the local relaxation frequency. This is the sum of the equilibrium,
the relaxed non-equilibrium population, and half of the Guo forcing term that
the general kernels evaluate direction by direction.

Since `H2(c):A2 = sum_ab c_a c_b A2_ab - cs^2 tr(A2)`, the trace term is moved
into the zeroth-order coefficient:

```text
f*_l = p_l [ A0' + c_l.A1 / cs^2 + sum_a c_la^2 A2_aa / (2 cs^4)
                 + sum_(a<b) c_la c_lb A2_ab / cs^4 ],
A0' = A0 - tr(A2) / (2 cs^2).
```

With this form, a diagonal coefficient `A2_aa` contributes only to the
directions with `c_la` different from zero.

The kernels stream in pull form. The moments of node `x` after streaming are

```text
m_j(x) = sum_l T_j(c_l) f*_l(x - c_l)
       = sum_k sum_l T_j(c_l) p_l B_k(c_l) A_k(x - c_l),
```

where `B_k` is the basis function multiplying `A_k` above and `T_j` takes the
values `1`, `c_a`, and `c_a c_b`. The second-order moments `c_a^2 - cs^2` are
obtained at the end from the sums of `c_a^2 f` and of `f`. Therefore each
coefficient can be streamed on its own: for every `k`, the coefficient field
`A_k` is pulled from all directions and added to the ten output moments.

The directions are processed in opposite pairs `(l, l+1)`. For the two values
`A_k(x - c_l)` and `A_k(x + c_l)`, only their sum and difference enter the
output moments.

With `BOUNCE_BACK`, a population pulled from a solid node is the reflected
population of the node itself, `f_l(x) = f*_(l+1)(x)`. In stage `k`, the value
pulled from the solid node is therefore replaced by `A_k(x)`, with a minus sign
for the first-order coefficients.

## Implementation

| Item | Location |
|---|---|
| Kernels | `fused_LB_kernel_2` (compiled without `BOUNCE_BACK`) and `fused_LB_kernel_bb_2` (compiled with `BOUNCE_BACK`) in `lb_cuda_fused_module.f90`. |
| Gather stages | `fused2_stages.h`, included by both kernels. |
| Stage generator | `gen_fused2_stages.py`. |
| Kernel selection | `fused_LB_cuda` in `lb_cuda_driver_module.f90`. |
| Macro | `FUSED2` in `defines.h`. |

The kernels have the same arguments, thread blocks, and prologue as the
general kernels. Each block of `10 x 10 x 10` threads covers an `8 x 8 x 8`
tile and its halo. Every thread reads the moments of its node and computes the
ten post-collision coefficients once. Each of the ten stages then

1. stores one coefficient in shared memory,
2. waits at a barrier, and
3. lets the interior threads add the contribution of every direction pair to
   their output moments.

The stages alternate between two shared buffers of `10^3` values, 16 KB in
double precision, so one barrier per stage is sufficient. The general kernels
need two barriers per direction, and every thread, including the halo threads,
reconstructs all the populations.

The ten stages are written out in `fused2_stages.h` for D3Q19 and D3Q27, with
the lattice components and weights as literals. A loop over the directions
would not be equivalent in cost. nvfortran unrolls such a loop, but it does
not replace elements of parameter arrays by their values, and multiplications
by zero components cannot be removed under IEEE arithmetic. The generator
reads the lattice vectors, their order, and the weight names from
`vars_module.f90`. Run it from the repository root after any change to the
D3Q19 or D3Q27 lattice definitions:

```bash
python3 gen_fused2_stages.py
```

Compiled for compute capability 8.0 with the `maxregcount:64` limit, the
kernels use the following resources:

| Kernel | FP64 instructions per thread | Barriers | Registers | Spills |
|---|---:|---:|---:|---:|
| `fused_LB_kernel`, D3Q19 FP64 | about 1200 | 36 | 64 | 20--28 bytes |
| `fused_LB_kernel_2`, D3Q19 FP64 | 374 | 10 | 62 | none |
| `fused_LB_kernel_bb_2`, D3Q19 FP64 | 377 | 10 | 64 | 4 bytes |
| `fused_LB_kernel`, D3Q27 FP64 | about 1730 | 52 | 64 | 20--28 bytes |
| `fused_LB_kernel_2`, D3Q27 FP64 | 854 | 10 | 64 | none |
| `fused_LB_kernel_bb_2`, D3Q27 FP64 | 857 | 10 | 64 | 8 bytes |

The instruction counts of the general kernels are per thread and include the
reconstruction of all populations, which every thread of the block performs.
The new kernels reconstruct no population: the halo threads only compute the
ten coefficients.

## Selection

`FUSED2` is defined in the default `defines.h`. It takes effect only for
`LATTICE 19` and `LATTICE 27` without `HIGHORDER`. For the high-order D3Q27
reconstruction and for D3Q15, the general kernels are used, and the build
prints the message `FUSED2: not applicable to this lattice`. Define `noFUSED2`
to use the general kernels for every lattice. In the executables produced by
`compile.sh`, `FUSED2` therefore applies to the `_19` and `_27` builds.

## Verification

The final moments of the `FUSED2` kernels were compared with those of the
general kernels, built with `noFUSED2`, in double precision unless stated
otherwise. Only round-off differences were found:

| Case | Lattice | Steps | Largest relative difference |
|---|---|---:|---:|
| Taylor--Green `128^3` | D3Q19 | 2000 | `1.3e-13` |
| Taylor--Green `128^3`, two MPI ranks | D3Q19 | 2000 | `1.3e-13` |
| Taylor--Green `128^3` | D3Q27 | 2000 | `1.1e-13` |
| Taylor--Green `128^3`, `SMAGORINSKI` | D3Q19 | 2000 | `1.1e-13` |
| Lamb droplet, two-component `80^3` | D3Q19 | 2000 | `6.2e-13` |
| Lamb droplet, two-component `80^3` | D3Q27 | 2000 | `7.8e-13` |
| Force-driven Poiseuille `32 x 8 x 32` | D3Q19 | 10000 | `8.1e-14` |
| Force-driven Poiseuille, `BOUNCE_BACK` | D3Q19 and D3Q27 | 10000 | `1.0e-13` |
| Taylor--Green `128^3`, single-precision storage (`MIXEDPRC`) | D3Q19 | 2000 | `9.4e-6` |
| Taylor--Green `128^3`, single precision | D3Q19 | 2000 | `6.8e-5` |

The relative difference is the largest absolute difference of each moment
field divided by the largest magnitude of that field. In the Poiseuille cases
it refers to the nonzero fields `w`, `pzz`, and `pxz`. The fields that vanish
in this flow contain round-off noise of order `1e-15` with `FUSED2` and of
order `1e-17` with the general kernels.

## Performance

Mean fused-kernel time per step and throughput of the whole time step on one
NVIDIA A30, for `256^3` lattices, averaged over two runs. The single-component
cases are the Taylor--Green vortex of
[`TAYLOR_GREEN_256_PERFORMANCE.md`](TAYLOR_GREEN_256_PERFORMANCE.md) with
`nsteps = 500`; the `BOUNCE_BACK` case is a force-driven channel, and the
two-component case is the static droplet of the default configuration.

| Case | Fused kernel, general (ms) | Fused kernel, `FUSED2` (ms) | Speed-up | GLUPS, general | GLUPS, `FUSED2` |
|---|---:|---:|---:|---:|---:|
| Single-component D3Q19, FP64 | 22.2 | 9.8 | 2.28 | 0.599 | 1.081 |
| Single-component D3Q27, FP64 | 29.7 | 15.5 | 1.91 | 0.472 | 0.787 |
| Single-component D3Q19, FP32 | 11.8 | 6.0 | 1.97 | 1.074 | 1.703 |
| Single-component D3Q27, FP32 | 15.1 | 9.1 | 1.65 | 0.886 | 1.293 |
| Two-component D3Q27, FP64 | 35.4 | 15.7 | 2.25 | 0.297 | 0.457 |
| Single-component D3Q27, FP64, `BOUNCE_BACK` | 31.3 | 17.8 | 1.76 | 0.471 | 0.759 |

The scaling tables in the README were measured with the general kernels.
