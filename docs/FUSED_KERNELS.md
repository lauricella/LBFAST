# Fused collide-and-stream kernels

## Purpose

At every time step, LBFAST performs the collision, the streaming, and the
computation of the new moments in a single GPU kernel, without storing the
populations in global memory. Two families of kernels implement this fused
collide-and-stream step:

| Lattice | Kernels |
|---|---|
| D3Q15, D3Q19, and D3Q27 (second order) | `fused_LB_kernel`, `fused_LB_kernel_bb` |
| High-order D3Q27 (`HIGHORDER`) | `fused_LB_kernel_h`, `fused_LB_kernel_bb_h` |

The `_bb` kernels are used when `BOUNCE_BACK` is defined. Both families
replace the kernels of LBFAST `v1.0.3` and give the same results to round-off.
On an NVIDIA A30, the fused step is 2.3 to 3.3 times faster for the
second-order lattices and 2.8 to 4.5 times faster for high-order D3Q27. This
document describes their formulation, implementation, verification, and
performance.

## Thread organization

The kernels are in `lb_cuda_fused_module.f90`, and `fused_LB_cuda` in
`lb_cuda_driver_module.f90` selects them from the lattice and from
`BOUNCE_BACK`. Only the selected kernel is compiled.

Each block has one thread per interior node of an `8 x 8 x 8` tile, that is
512 threads. Every thread reads the moments of its node and computes the values
that the neighbors need to reconstruct its post-collision populations. The
first 488 threads also process one node of the halo around the tile, so that
each of the 1000 nodes of the tile and its halo is processed once. These values
are stored in dynamic shared memory. After a barrier, every thread gathers
(pull) the contributions of all directions and writes the ten output moments of
its node. The second-order kernels need a single barrier; the high-order
kernels add a second pass for the coefficients of orders three to six and use
three barriers.

The kernels of `v1.0.3` used blocks of `10 x 10 x 10` threads, one per node of
the tile and its halo. Every thread, including the 488 halo threads,
reconstructed all the populations, each direction needed two barriers, and only
one block of 1000 threads fitted on a multiprocessor.

More than 48 KB of shared memory per block must be requested explicitly. At the
first launch, `fused_LB_cuda` calls `cudaFuncSetAttribute` to allow the
required dynamic shared-memory size and to prefer the largest shared-memory
carveout. If the GPU cannot provide the memory, the run stops with an error
message.

| Kernels | Values per node | Shared memory per block, FP64 / FP32 | Blocks per A30 multiprocessor, FP64 / FP32 |
|---|---:|---:|---:|
| Second order | 10 | 80 KB / 40 KB | 2 / 4 |
| High order | 19 | 152 KB / 76 KB | 1 / 2 |

With `MIXEDPRC`, the arithmetic and the shared buffer are in double precision,
so the double-precision requirements apply. The numbers of blocks were reported by
`cudaOccupancyMaxActiveBlocksPerMultiprocessor`. The A30, A100, and H100
provide enough shared memory for all kernels. A GPU with less than 152 KB of
shared memory per block, such as the V100 with 96 KB, cannot run the
high-order kernels in double precision.

## Second-order kernels

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
the relaxed non-equilibrium population, and half of the Guo forcing term.

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

In these kernels, the values stored per node are the ten coefficients `A0'`,
`A1`, and `A2`, and the gather consists of ten stages, one per coefficient.
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

## High-order kernels

The high-order D3Q27 reconstruction adds to the second-order populations the
terms of orders three to six,

```text
f*_l += sum_(k=10..26) Minv_mrt(k,l) h_k,
h_k = E_k(u) + N_k(u, (1 - omega) Pi_neq),
```

where `E_k` are the equilibrium coefficients, products of the velocity
components, and `N_k` the non-equilibrium coefficients obtained with the
recursive relation of Malaspinas; the forcing term is of second order only.
Since the coefficients `h_k` are nonlinear functions of the moments, they are
computed once per node and then gathered linearly in a second pass:

1. Every thread stores, for its node and its halo node, the ten second-order
   coefficients together with `u` and `(1 - omega) Pi_neq`: 19 values per node.
2. After a barrier, the second-order coefficients are gathered with the stages
   of `fused2_stages.h`, as in the second-order kernels.
3. After a second barrier, every thread computes the seventeen coefficients
   `h_k` of its nodes from `u` and `(1 - omega) Pi_neq` and stores them in the
   same shared buffer.
4. After a third barrier, the seventeen coefficients are gathered with the
   stages of `fusedh_stages.h`.

`fusedh_stages.h` is generated by `gen_fused2_stages.py`, which writes the
elements of `Minv_mrt` as exact fractions, `p_l H_k(c_l) / N_k` for the
tensor-product Hermite basis; all elements used coincide with the values of
`Minv_mrt` in double precision. By orthogonality, the sum over the directions
of `T_j(c_l) Minv_mrt(k,l)` vanishes for every output moment, so the value of
the node itself is subtracted from the pulled values without changing the exact
result. The high-order contributions are summed in separate accumulators and
added to the output moments at the end, which limits the round-off of the many
small terms in single precision.

With `BOUNCE_BACK`, a coefficient pulled from a solid node is replaced by the
value of the node itself, with the sign of the reflected direction, as in the
second-order kernels.

The shared memory limits double precision to one block per multiprocessor, and
the kernels are compiled with launch bounds that allow up to 128 registers. In
single precision, the launch bounds give two blocks per multiprocessor with at
most 64 registers. The launch bounds require `8 x 8 x 8` tiles; the build stops
with an error for other tile sizes.

An earlier version of the high-order kernels stored 13 values per node and
reconstructed every incoming population from the moments of the upstream node,
so the high-order coefficients of each node were recomputed for all 26
directions. It used about 5900 double-precision instructions per node, against
about 2600 for the current version, and was about twice as slow.

## Resources

Compiled for compute capability 8.0, the kernels use the following resources:

| Kernel | Threads per block | Barriers | Registers | Spills |
|---|---:|---:|---:|---:|
| `v1.0.3` kernels, D3Q19 or D3Q27, FP64 | 1000 | 36 or 52 | 64 | 20--28 bytes |
| `v1.0.3` kernels, high-order D3Q27, FP64 | 1000 | 52 | 64 | about 250 bytes |
| `fused_LB_kernel`, D3Q15, FP64 | 512 | 1 | 49 | none |
| `fused_LB_kernel`, D3Q19 or D3Q27, FP64 | 512 | 1 | 51 | none |
| `fused_LB_kernel`, FP32 | 512 | 1 | 32 | none |
| `fused_LB_kernel`, two-component D3Q27, FP64 | 512 | 1 | 60 | none |
| `fused_LB_kernel_bb`, D3Q19, FP64 | 512 | 1 | 64 | none |
| `fused_LB_kernel_h`, FP64 | 512 | 3 | 95 | none |
| `fused_LB_kernel_h`, two-component, FP64 | 512 | 3 | 92 | none |
| `fused_LB_kernel_h`, FP32 | 512 | 3 | 56 | none |
| `fused_LB_kernel_bb_h`, FP64 | 512 | 3 | 128 | none |

An earlier version of the second-order kernels used blocks of `10 x 10 x 10`
threads, two shared buffers holding one coefficient at a time, and ten
barriers. It gave bitwise identical results and was 1.3 to 1.8 times slower
than the current one.

## Verification

The final moments were compared with those of the `v1.0.3` kernels, in double
precision unless stated otherwise. Only round-off differences were found:

| Case | Lattice | Steps | Largest relative difference |
|---|---|---:|---:|
| Taylor--Green `128^3` | D3Q15 | 2000 | `1.3e-13` |
| Taylor--Green `128^3` | D3Q19 | 2000 | `1.3e-13` |
| Taylor--Green `128^3` | D3Q27 | 2000 | `1.1e-13` |
| Taylor--Green `128^3` | high-order D3Q27 | 2000 | `1.2e-13` |
| Taylor--Green `128^3`, two MPI ranks | D3Q15, D3Q19 | 2000 | `1.4e-13` |
| Taylor--Green `128^3`, two MPI ranks | high-order D3Q27 | 2000 | `1.2e-13` |
| Taylor--Green `128^3`, `U0 = 0.15` (Mach 0.26) | high-order D3Q27 | 500 | `1.6e-14` |
| Taylor--Green `128^3`, `SMAGORINSKI` | D3Q19 | 2000 | `1.1e-13` |
| Lamb droplet, two-component `80^3` | D3Q15, D3Q19, D3Q27 | 2000 | `7.8e-13` |
| Lamb droplet, two-component `80^3` | high-order D3Q27 | 2000 | `9.0e-14` |
| Force-driven Poiseuille `32 x 8 x 32` | D3Q19 | 10000 | `8.1e-14` |
| Force-driven Poiseuille, `BOUNCE_BACK` | D3Q15, D3Q19, D3Q27 | 10000 | `1.3e-13` |
| Force-driven Poiseuille, `BOUNCE_BACK` | high-order D3Q27 | 10000 | `9.9e-14` |
| Taylor--Green `128^3`, single-precision storage (`MIXEDPRC`) | D3Q19 | 2000 | `9.4e-6` |
| Taylor--Green `128^3`, single-precision storage (`MIXEDPRC`) | high-order D3Q27 | 2000 | `7.1e-6` |
| Taylor--Green `128^3`, single precision | D3Q15, D3Q19 | 2000 | `1.2e-4` |
| Taylor--Green `128^3`, single precision | high-order D3Q27 | 2000 | `1.3e-4` |

The relative difference is the largest absolute difference of each moment
field divided by the largest magnitude of that field. In the Poiseuille cases
it refers to the nonzero fields `w`, `pzz`, and `pxz`. The fields that vanish
in this flow contain round-off noise of order `1e-15` with the second-order
kernels and of order `1e-17` with the `v1.0.3` kernels.

The validation cases run at low Mach number, where the terms of orders three
to six are small. The Taylor--Green case with `U0 = 0.15` tests them directly:
there, the high-order and second-order D3Q27 reconstructions differ by `1e-2`
in relative terms, while the high-order kernels agree with the `v1.0.3` kernels
to `1.6e-14`.

With `MIXEDPRC`, the accumulation is in double precision, and the deviations
from the double-precision solution are those of the single-precision storage,
as with the `v1.0.3` kernels.

In single precision the gather sums more terms per output moment than the
`v1.0.3` kernels, so the round-off is somewhat larger. After 2000 steps of the
`128^3` Taylor--Green vortex, the largest deviations of the moment fields from
the double-precision solution are up to 2.6 times those of the `v1.0.3`
kernels, below `1e-4` in relative terms. The kinetic energy differs from the
`v1.0.3` single-precision run by `1e-5` to `1e-4` in relative terms, depending
on the lattice, much less than the deviation of both single-precision runs from
the double-precision result, about `1.2e-3`. Summing the contributions of each
second-order stage in partial sums before adding them to the output moments was
tried and did not reduce these deviations, so the larger round-off comes from
the coefficient formulation rather than from the number of accumulations.

The complete validation suite, which uses high-order D3Q27 for all its cases,
passes with the high-order kernels (`5 passed, 0 failed` on an A30). Compared
with the `v1.0.3` kernels, the Laplace, Lamb, and Poiseuille results are
identical at all printed digits, and the capillary-wave and Taylor--Green
summaries agree to `2e-11` and `1e-9` in relative terms. The reference CSV
files in `docs/` were regenerated from this run.

## Performance

Mean fused-kernel time per step and throughput of the whole time step on one
NVIDIA A30, for `256^3` lattices, averaged over two runs. The single-component
cases are the Taylor--Green vortex of
[`TAYLOR_GREEN_256_PERFORMANCE.md`](TAYLOR_GREEN_256_PERFORMANCE.md) with
`nsteps = 500`; the `BOUNCE_BACK` case is a force-driven channel, and the
two-component case is the static droplet of the default configuration.

| Case | Fused kernel, `v1.0.3` (ms) | Fused kernel, current (ms) | Speed-up | GLUPS, `v1.0.3` | GLUPS, current |
|---|---:|---:|---:|---:|---:|
| Single-component D3Q15, FP64 | 18.4 | 7.3 | 2.5 | 0.694 | 1.286 |
| Single-component D3Q19, FP64 | 22.2 | 7.4 | 3.0 | 0.599 | 1.278 |
| Single-component D3Q27, FP64 | 29.7 | 9.1 | 3.3 | 0.472 | 1.128 |
| Single-component D3Q15, FP32 | 10.1 | 4.3 | 2.3 | 1.202 | 2.053 |
| Single-component D3Q19, FP32 | 11.8 | 4.7 | 2.5 | 1.074 | 1.972 |
| Single-component D3Q27, FP32 | 15.1 | 5.0 | 3.0 | 0.886 | 1.894 |
| Two-component D3Q27, FP64 | 35.4 | 10.6 | 3.3 | 0.297 | 0.527 |
| Single-component D3Q27, FP64, `BOUNCE_BACK` | 31.3 | 12.5 | 2.5 | 0.471 | 0.995 |
| Single-component high-order D3Q27, FP64 | 82.7 | 23.6 | 3.5 | 0.190 | 0.572 |
| Single-component high-order D3Q27, FP32 | 32.1 | 11.3 | 2.8 | 0.467 | 1.108 |
| Two-component high-order D3Q27, FP64 | 108.5 | 24.3 | 4.5 | 0.129 | 0.369 |
| Single-component high-order D3Q27, FP64, `BOUNCE_BACK` | 114.4 | 26.9 | 4.3 | 0.141 | 0.538 |

With D3Q19 in double precision, the second-order kernel reads and writes about
280 bytes per updated node, which corresponds to about 650 GB/s, or 70% of the
A30 peak memory bandwidth. The fused kernel now takes about 57% of the
single-component time step for D3Q19; the moment and halo kernels take most of
the rest.

The scaling tables in the README were measured with the `v1.0.3` kernels.
