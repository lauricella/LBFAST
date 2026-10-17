# Taylor--Green vortex performance test on a 256-cubed lattice

## Purpose

This case measures the single-GPU throughput of the single-component solver.
It uses the Reynolds number of the large `taylorgreen.inp` case, `Re = 80`, on
a `256 x 256 x 256` lattice that fits on one GPU. The kernel timings reported
in [`RELEASE_NOTES_v1.0.3.md`](RELEASE_NOTES_v1.0.3.md) and in
[`FUSED_KERNELS.md`](FUSED_KERNELS.md) were measured with
this input.

The case is not part of the smoke or complete validation suites. It has no
automated driver and no pass/fail tolerance. For a quantitative viscosity
check, use [`TAYLOR_GREEN_128_VALIDATION.md`](TAYLOR_GREEN_128_VALIDATION.md).

## Taylor--Green inputs

All three inputs use the initial condition described in
[`TAYLOR_GREEN_128_VALIDATION.md`](TAYLOR_GREEN_128_VALIDATION.md), with
`U0 = uwall` and `nu = visc1`. The characteristic length is `L = N / (2 pi)`
and the Mach number is `U0 / cs`, with `cs = 1/sqrt(3)`.

| Input | Domain | `U0` | `nu` | `Re = U0 L / nu` | Mach | Time steps | Use |
|---|---|---:|---:|---:|---:|---:|---|
| `taylorgreen.inp` | `512^3` | `0.02` | `0.02037` | `80.0` | `0.035` | `200000` | Large reference case |
| `taylorgreen256.inp` | `256^3` | `0.04` | `0.02037` | `80.0` | `0.069` | `50000` | Single-GPU performance test |
| `taylorgreen128.inp` | `128^3` | `0.04` | `0.04074` | `20.0` | `0.069` | `2000` | Automated viscosity validation |

The `256^3` case keeps `Re = 80` by doubling `U0` with respect to the `512^3`
case, so its Mach number is also twice as large.

## Configuration

| Parameter | Value |
|---|---:|
| Domain | `256 x 256 x 256`, fully periodic |
| Time steps | `50000` |
| Energy sampling interval (`stamp_term`) | `50` steps |
| Initial velocity `U0` | `0.04` |
| Kinematic viscosity | `0.02037` |
| Characteristic length `L` | `40.7437` |
| Reynolds number | `80.0` |
| VTK and RAW field output | disabled |

With `50000` time steps, the kinetic energy decays to about `2.5%` of its
initial value. For a throughput measurement, a few hundred steps are enough;
the reference values below were obtained with `nsteps = 500`.

## Build and run

Edit `defines.h` as follows, and keep the other macros at their default
values:

| Macro | Setting |
|---|---|
| `TAYLORGREEN` | enabled (replace `noTAYLORGREEN`) |
| `TWOCOMPONENT`, `DENSRATIO`, `CSF` | disabled (`noTWOCOMPONENT`, `noDENSRATIO`, `noCSF`) |
| `LATTICE`, `HIGHORDER` | the lattice to be measured, for example `LATTICE 27` with `HIGHORDER` |
| `PRC`, `STRPRC` | `8` for double precision, `4` for single precision |

Then build and run from the repository root:

```bash
make clean
make nvfortran GPUCC=80
./main.x taylorgreen256.inp > tg256.log
```

With an MPI build (`make nvfortran-mpi GPUCC=80`), run the case on one rank
with

```bash
mpirun -np 1 ./main.x 1 1 1 taylorgreen256.inp > tg256.log
```

Restore the default configuration afterwards, for example with
`git checkout -- defines.h`, and rebuild before running other cases.

For timing, run one case at a time on an otherwise idle GPU. Concurrent runs
on the same node can distort the measured times.

## Output and metrics

| Output | Content |
|---|---|
| `glups:` in the log | Throughput of the whole time step, `lx*ly*lz*nsteps / (t_loop * 1e9)`. |
| `time elapsed:` in the log | Wall-clock time of the time loop. |
| `DEVICE memory occupied at the start` in the log | GPU memory allocated by the solver. |
| `time.dat` | Timing of each routine: time share, section, routine, mean time per call, number of calls, minimum, maximum, and standard deviation. The fused collide-and-stream kernel is the `fused` row. |
| `taylorgreen.dat` | Time step, `E/E0`, and `log(E/E0)` every `stamp_term` steps. |
| `taylorgreen_theory.dat` | Analytical decay `E/E0 = exp[-6 nu (2 pi / N)^2 t]` at the same steps. |
| `plot_taylorgreen.gp` | Gnuplot script comparing the two energy histories; its title reports the fitted viscosity `nu_fit`. |

## Reference results

The following values were obtained on one NVIDIA A30 24 GB GPU, using an MPI
build on one rank and `nsteps = 500`. The kernel time and the throughput are
the mean of two runs. The table reports the current fused kernels, which also
apply the moment update in single-component flows (see
[`FUSED_KERNELS.md`](FUSED_KERNELS.md)), and, for comparison, the kernels of
`v1.0.3`.

| Lattice | Precision | Fused kernels | GPU memory (GB) | Fused kernel (ms/step) | GLUPS | `nu_fit` |
|---|---|---|---:|---:|---:|---:|
| D3Q15 | FP64 | current | `4.08` | `8.0` | `1.738` | `0.020523` |
| D3Q15 | FP64 | `v1.0.3` | `4.08` | `18.4` | `0.694` | `0.020523` |
| D3Q19 | FP64 | current | `4.08` | `7.0` | `1.933` | `0.020527` |
| D3Q19 | FP64 | `v1.0.3` | `4.08` | `22.2` | `0.599` | `0.020527` |
| D3Q27 | FP64 | current | `4.08` | `9.4` | `1.509` | `0.020526` |
| D3Q27 | FP64 | `v1.0.3` | `4.08` | `29.8` | `0.472` | `0.020526` |
| High-order D3Q27 | FP64 | current | `4.08` | `23.5` | `0.667` | `0.020524` |
| High-order D3Q27 | FP64 | `v1.0.3` | `4.08` | `83.4` | `0.188` | `0.020524` |
| High-order D3Q27 | FP32 | current | `2.35` | `11.8` | `1.257` | `0.022233` |
| High-order D3Q27 | FP32 | `v1.0.3` | `2.35` | `32.1` | `0.466` | `0.022233` |

The GPU memory does not depend on the lattice or on the fused kernels, because
LBFAST stores the same ten moments per node in every case.

After `500` steps, the high-order D3Q27 FP64 run gives `E/E0 = 0.963387`,
against the analytical value `0.963857`. Over this short interval the fitted
viscosity in FP64 exceeds the input value by about `0.8%`, while in FP32 it
exceeds it by about `9%`. Use double precision when the decay rate matters.

At these throughputs, the complete `50000`-step run takes about `7` minutes
with D3Q19 FP64, about `9` minutes with D3Q27 FP64, and about `21` minutes
with high-order D3Q27 FP64 on an A30.
