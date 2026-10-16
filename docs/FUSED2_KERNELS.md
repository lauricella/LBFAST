# Second-order fused kernels (`FUSED2`)

## Purpose

For D3Q15, D3Q19, and D3Q27 without `HIGHORDER`, LBFAST performs the fused
collide-and-stream step with the kernels `fused_LB_kernel_2` and
`fused_LB_kernel_bb_2`, selected by the `FUSED2` macro. They give the same
results as the general kernels `fused_LB_kernel` and `fused_LB_kernel_bb`, to
round-off, and the fused kernel runs 2.3 to 3.3 times faster on an NVIDIA A30.
This document describes their formulation, implementation, verification, and
performance.

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

The kernels have the same arguments and per-node prologue as the general
kernels, but a different thread organization. Each block has one thread per
interior node of an `8 x 8 x 8` tile, that is 512 threads. Every thread reads
the moments of its node and computes the ten post-collision coefficients. The
first 488 threads also compute the coefficients of one node of the halo around
the tile, so that each of the 1000 nodes of the tile and its halo is processed
once. All ten coefficients of these nodes are stored in dynamic shared memory,
80 KB in double precision and 40 KB in single precision. After a single
barrier, every thread executes the ten gather stages and writes the ten output
moments of its node.

The general kernels use blocks of `10 x 10 x 10` threads, one per node of the
tile and its halo. Every thread, including the 488 halo threads, reconstructs
all the populations, and each direction needs two barriers. With one block of
1000 threads per multiprocessor, the multiprocessor is idle at every barrier.
The second-order kernels instead fit two blocks per multiprocessor in double
precision and four in single precision on the A30, as reported by
`cudaOccupancyMaxActiveBlocksPerMultiprocessor`, so that the global-memory
reads of one block overlap the gather of another.

More than 48 KB of shared memory per block must be requested explicitly. At the
first launch, `fused_LB_cuda` calls `cudaFuncSetAttribute` to allow the
required dynamic shared-memory size and to prefer the largest shared-memory
carveout. If the GPU cannot provide the memory, the run stops with a message
that suggests `noFUSED2`. The A30, A100, and H100 provide enough shared memory
for two blocks per multiprocessor in double precision.

The ten stages are written out in `fused2_stages.h` for D3Q15, D3Q19, and
D3Q27, with the lattice components and weights as literals. A loop over the
directions would not be equivalent in cost. nvfortran unrolls such a loop, but
it does not replace elements of parameter arrays by their values, and
multiplications by zero components cannot be removed under IEEE arithmetic. The
generator reads the lattice vectors, their order, and the weight names from
`vars_module.f90`. Run it from the repository root after any change to the
D3Q15, D3Q19, or D3Q27 lattice definitions:

```bash
python3 gen_fused2_stages.py
```

Compiled for compute capability 8.0 with the `maxregcount:64` limit, the
kernels use the following resources:

| Kernel | Threads per block | Barriers | Registers | Spills | Shared memory per block |
|---|---:|---:|---:|---:|---:|
| `fused_LB_kernel`, D3Q19 or D3Q27, FP64 | 1000 | 36 or 52 | 64 | 20--28 bytes | 8 KB |
| `fused_LB_kernel_2`, D3Q15, FP64 | 512 | 1 | 49 | none | 80 KB |
| `fused_LB_kernel_2`, D3Q19 or D3Q27, FP64 | 512 | 1 | 51 | none | 80 KB |
| `fused_LB_kernel_2`, D3Q15, D3Q19, or D3Q27, FP32 | 512 | 1 | 32 | none | 40 KB |
| `fused_LB_kernel_2`, two-component D3Q27, FP64 | 512 | 1 | 60 | none | 80 KB |
| `fused_LB_kernel_bb_2`, D3Q27, FP64 | 512 | 1 | 64 | 4 bytes | 80 KB |

An earlier implementation of these kernels used blocks of `10 x 10 x 10`
threads, two shared buffers holding one coefficient at a time, and ten
barriers. It gave bitwise identical results and was 1.3 to 1.8 times slower
than the current one.

## Selection

`FUSED2` is defined in the default `defines.h`. It takes effect for
`LATTICE 15`, `LATTICE 19`, and `LATTICE 27` without `HIGHORDER`. For the
high-order D3Q27 reconstruction, the general kernels are used, and the build
prints the message `FUSED2: not applicable with HIGHORDER`. Define `noFUSED2`
to use the general kernels for every lattice. In the executables produced by
`compile.sh`, `FUSED2` therefore applies to the `_15`, `_19`, and `_27`
builds.

## Verification

The final moments of the `FUSED2` kernels were compared with those of the
general kernels, built with `noFUSED2`, in double precision unless stated
otherwise. Only round-off differences were found:

| Case | Lattice | Steps | Largest relative difference |
|---|---|---:|---:|
| Taylor--Green `128^3` | D3Q15 | 2000 | `1.3e-13` |
| Taylor--Green `128^3`, two MPI ranks | D3Q15 | 2000 | `1.4e-13` |
| Taylor--Green `128^3` | D3Q19 | 2000 | `1.3e-13` |
| Taylor--Green `128^3`, two MPI ranks | D3Q19 | 2000 | `1.3e-13` |
| Taylor--Green `128^3` | D3Q27 | 2000 | `1.1e-13` |
| Taylor--Green `128^3`, `SMAGORINSKI` | D3Q19 | 2000 | `1.1e-13` |
| Lamb droplet, two-component `80^3` | D3Q15 | 2000 | `5.9e-13` |
| Lamb droplet, two-component `80^3` | D3Q19 | 2000 | `6.2e-13` |
| Lamb droplet, two-component `80^3` | D3Q27 | 2000 | `7.8e-13` |
| Force-driven Poiseuille `32 x 8 x 32` | D3Q19 | 10000 | `8.1e-14` |
| Force-driven Poiseuille, `BOUNCE_BACK` | D3Q15, D3Q19, and D3Q27 | 10000 | `1.3e-13` |
| Taylor--Green `128^3`, single-precision storage (`MIXEDPRC`) | D3Q19 | 2000 | `9.4e-6` |
| Taylor--Green `128^3`, single precision | D3Q19 | 2000 | `6.8e-5` |
| Taylor--Green `128^3`, single precision | D3Q15 | 2000 | `1.2e-4` |

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
| Single-component D3Q15, FP64 | 18.4 | 7.3 | 2.5 | 0.694 | 1.286 |
| Single-component D3Q19, FP64 | 22.2 | 7.4 | 3.0 | 0.599 | 1.278 |
| Single-component D3Q27, FP64 | 29.7 | 9.1 | 3.3 | 0.472 | 1.128 |
| Single-component D3Q15, FP32 | 10.1 | 4.3 | 2.3 | 1.202 | 2.053 |
| Single-component D3Q19, FP32 | 11.8 | 4.7 | 2.5 | 1.074 | 1.972 |
| Single-component D3Q27, FP32 | 15.1 | 5.0 | 3.0 | 0.886 | 1.894 |
| Two-component D3Q27, FP64 | 35.4 | 10.6 | 3.3 | 0.297 | 0.527 |
| Single-component D3Q27, FP64, `BOUNCE_BACK` | 31.3 | 12.5 | 2.5 | 0.471 | 0.995 |

With D3Q19 in double precision, the second-order kernel reads and writes about
280 bytes per updated node, which corresponds to about 650 GB/s, or 70% of the
A30 peak memory bandwidth. With `FUSED2`, the fused kernel takes about 57% of the
single-component time step for D3Q19; the moment and halo kernels take most of
the rest.

The scaling tables in the README were measured with the general kernels.

## Experimental high-order kernels (`FUSEDH`)

For the high-order D3Q27 reconstruction (`HIGHORDER`), the macro `FUSEDH`,
disabled by default, selects the kernels `fused_LB_kernel_h` and
`fused_LB_kernel_bb_h` in `lb_cuda_fused_module.f90`. They use the thread
organization of the second-order kernels: one thread per interior node, the
halo nodes computed by the first 488 threads, and a single barrier.

The terms of orders three to six are nonlinear functions of the moments, so the
post-collision coefficients cannot be gathered linearly as in the second-order
kernels. Each thread instead stores 13 values for its node: `p*`, `u`,
`(1 - omega) Pi_neq`, and `F`. These take 104 KB of shared memory per block in
double precision and 52 KB in single precision. After the barrier, every thread
reconstructs the 27 incoming post-collision populations from the values of the
upstream nodes, with the same high-order expressions as `fused_LB_kernel`. With
`BOUNCE_BACK`, a population pulled from a solid node is reconstructed from the
values of the node itself for the opposite direction.

The shared memory limits double precision to one block per multiprocessor, and
the kernels are compiled with launch bounds that allow up to 128 registers. In
single precision, the launch bounds give two blocks per multiprocessor with 64
registers. The kernels have no register spills and require `8 x 8 x 8` tiles.

Compared with the general kernels in double precision, unless stated
otherwise, the final moments differ only by round-off:

| Case | Steps | Largest relative difference |
|---|---:|---:|
| Taylor--Green `128^3` | 2000 | `3.2e-14` |
| Taylor--Green `128^3`, two MPI ranks | 2000 | `3.2e-14` |
| Lamb droplet, two-component `80^3` | 2000 | `3.8e-14` |
| Force-driven Poiseuille, `BOUNCE_BACK` | 10000 | none in `w`, `pzz`, `pxz` |
| Taylor--Green `128^3`, single precision | 2000 | `2.0e-5` |

Mean fused-kernel time per step and whole-step throughput on one NVIDIA A30 for
the `256^3` cases of the performance section above, averaged over two runs:

| Case | Fused kernel, general (ms) | Fused kernel, `FUSEDH` (ms) | Speed-up | GLUPS, general | GLUPS, `FUSEDH` |
|---|---:|---:|---:|---:|---:|
| Single-component, FP64 | 82.7 | 46.0 | 1.80 | 0.190 | 0.324 |
| Single-component, FP32 | 32.1 | 22.5 | 1.43 | 0.467 | 0.635 |
| Two-component, FP64 | 108.5 | 46.5 | 2.33 | 0.129 | 0.248 |
| Single-component, FP64, `BOUNCE_BACK` | 114.4 | 47.3 | 2.42 | 0.141 | 0.325 |

The general high-order kernels remain the default until `FUSEDH` has been used
more widely.
