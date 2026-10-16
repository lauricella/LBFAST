# Strong- and weak-scaling benchmarks

## Purpose

This document explains how the strong- and weak-scaling measurements reported
in the main README and in the LBFAST publication were produced, and how to
repeat them. It covers the input files, the executable matrix, the job
scripts, the decompositions, and the extraction of throughput and energy
figures from the run logs.

The published measurements were obtained on the Leonardo Booster system at
CINECA. Each node provides four NVIDIA A100 64 GB GPUs; the largest runs use
128 nodes, that is, 512 GPUs with one MPI rank per GPU.

## Benchmark cases

Two fully periodic `512^3` configurations are used, one for each physical
model.

| Model | Strong-scaling input | Weak-scaling input | Executables |
|---|---|---|---|
| One-component | `bench1.inp` | `bench1w.inp` | `main_1c_*.x` |
| Two-component | `bench.inp` | `benchw.inp` | `main_2c_*.x` |

The files `test512_1.inp`, `test512_1w.inp`, `test512.inp`, and
`test512_w.inp` are identical copies of the four inputs above and are used by
the older `run_benchmark.sh` template.

All four inputs run 100 time steps, disable VTK and RAW field output, and
print the probe values once at the end of the run. The weak-scaling inputs
differ from the strong-scaling ones only by `lweakscaling = .true.`.

### Initial conditions

The executables built by `compile.sh` keep every benchmark selector disabled
(`noTAYLORGREEN`, `noLAPLACE`, and so on), while `DOBENCHMARK` remains
defined. The runs therefore use the generic initial condition in
`init_conditions_module.f90`:

| Model | Initial state | Main physical parameters |
|---|---|---|
| Two-component | Static spherical droplet of radius `16` and width `8`, centred at `(112, 64, 64)`, with the Laplace pressure jump `2 sigma / R` preloaded inside the droplet. | `sigma = 0.01`, densities `5 / 0.5`, kinematic viscosities `0.2 / 0.025`, `tau_diff = 0.12`. |
| One-component | Fluid at rest with a localized pressure perturbation of amplitude `0.01` inside the same sphere. | Kinematic viscosity `0.1`. |

Both inputs also impose a uniform body acceleration `fx = 1.0e-6`.

The Taylor--Green velocity field is set only when `TAYLORGREEN` is defined
before compilation, and `compile.sh` does not define it. To run the
one-component scaling with the Taylor--Green initial condition, change
`noTAYLORGREEN` to `TAYLORGREEN` in `defines.h`, build the one-component
executables individually with `make nvfortran-nvml-mpi`, and set the initial
amplitude through `uwall` in the input. Do not run `compile.sh` with this
change, because `TAYLORGREEN` also disables the two-component model and the
`main_2c_*` executables would then contain the one-component solver.

## Building the executables

Build the complete executable matrix from the repository root with

```bash
./compile.sh
```

The script and its requirements are described in
[`COMPILE_SCRIPT.md`](COMPILE_SCRIPT.md). It uses the
`nvfortran-nvml-mpi` Make target for compute capability 8.0, so every
executable supports MPI and NVML energy measurement.

Each executable name encodes the model, the lattice, and the precision:

| Name component | Meaning |
|---|---|
| `main_1c_`, `main_2c_` | One-component or two-component solver. |
| `15`, `19`, `27`, `27high` | D3Q15, D3Q19, D3Q27, or high-order D3Q27 lattice. |
| `.x` | Single precision for arithmetic and storage. |
| `_d.x` | Double precision for arithmetic and storage. |
| `_sd.x` | Double-precision arithmetic with single-precision storage. |

The published tables report the single-precision (`.x`) and
double-precision (`_d.x`) executables for D3Q19, D3Q27, and high-order
D3Q27. The job scripts also run the D3Q15 and `_sd.x` variants, which are not
included in the tables.

The published measurements predate the current fused collide-and-stream
kernels (see [`FUSED_KERNELS.md`](FUSED_KERNELS.md)), which are two to three
times faster for the second-order lattices and about twice as fast for
high-order D3Q27. Executables built from the current sources are therefore
faster than the published figures. The code used for the publication
corresponds to the kernels of tag `v1.0.3` and earlier.

Rebuild the executables after any source change. The binaries must be present
in the directory from which the job is submitted, together with the four
`bench*.inp` files.

## Domain size and decomposition

**Strong scaling.** The global domain is the `512^3` box given by `lx`, `ly`,
and `lz` in the input. It is divided among `PX x PY x PZ` ranks.

**Weak scaling.** When `lweakscaling` is true, each rank receives a cubic
subdomain of side `MYSIDE`, and the global domain becomes
`MYSIDE*PX x MYSIDE*PY x MYSIDE*PZ`. The values of `lx`, `ly`, and `lz` in
the input are overridden. `MYSIDE` is a compile-time macro with default value
`512`, defined in `LBFAST.f90` inside the `DOBENCHMARK` block. Keep
`DOBENCHMARK` defined, as in the default `defines.h`. To use a different
subdomain, add for example `#define MYSIDE 256` to `defines.h` and rebuild.

In both modes every local extent must be a multiple of the `8 x 8 x 8` GPU
tile, and `PX*PY*PZ` must equal the number of MPI ranks.

The job scripts run the following decompositions:

| GPUs | Nodes | Weak scaling | Strong scaling |
|---:|---:|---|---|
| 1 | 1 | `1x1x1` | `1x1x1` |
| 2 | 1 | `1x1x2` | `1x1x2` |
| 4 | 1 | `1x1x4`, `1x2x2` | `1x1x4`, `1x2x2` |
| 8 | 2 | `1x1x8`, `1x2x4`, `2x2x2` | `1x2x4`, `2x2x2` |
| 16 | 4 | `1x1x16`, `1x4x4`, `2x2x4` | `1x4x4`, `2x2x4` |
| 32 | 8 | `1x1x32`, `1x4x8`, `2x4x4` | `1x4x8`, `2x4x4` |
| 64 | 16 | `1x1x64`, `1x8x8`, `2x4x8`, `4x4x4` | `1x8x8`, `2x4x8`, `4x4x4` |
| 128 | 32 | `1x1x128`, `1x8x16`, `4x4x8` | `1x8x16`, `4x4x8` |
| 256 | 64 | `1x1x256`, `1x16x16`, `4x8x8` | `1x16x16`, `4x8x8` |
| 512 | 128 | `1x1x512`, `1x16x32`, `4x8x16`, `8x8x8` | `1x16x32`, `4x8x16`, `8x8x8` |

The published tables list a subset of these decompositions.

## Job scripts

The scripts `run1.sh`, `run2.sh`, `run4.sh`, `run8.sh`, `run16.sh`,
`run32.sh`, `run64.sh`, and `run128.sh` request 1 to 128 Leonardo nodes. They
share the same body and differ only in the Slurm node count, the job name,
the quality of service, and the variable `mynodes`. The body contains one
block per node count, and only the block matching `mynodes` is executed.
Within that block, every executable is run in weak and strong mode for each
decomposition listed above.

The scripts assume the following environment:

| Setting | Value |
|---|---|
| Modules | `nvhpc/24.3`, `openmpi/4.1.6--nvhpc--24.3` |
| Ranks and GPUs | 4 ranks per node, one GPU per rank, exclusive nodes |
| Launcher | `mpirun --bind-to core --map-by ppr:4:node:pe=8` |
| MPI I/O | `OMPI_MCA_io=ompio`, `OMPI_MCA_sharedfp=individual` |
| Wall time | 2 hours per job |
| Account and partition | `IscrB_MIPLAST`, `boost_usr_prod`; `run128.sh` also uses `boost_qos_bprod` |

Replace the account, partition, and quality of service with those of your own
allocation. Submit one job per node count from the directory containing the
executables and inputs, for example

```bash
sbatch run16.sh
```

On a system without Slurm, a single run can be launched directly. For
example, a four-GPU weak-scaling run with the one-component D3Q27 executable
is

```bash
mpirun -np 4 ./main_1c_27.x 1 1 4 bench1w.inp > mio1_27_w_4_1_1_4.dat
```

The device memory must accommodate one `512^3` subdomain per GPU in weak
scaling, and the whole box divided by the number of ranks in strong scaling.
On GPUs with less memory than the A100 64 GB, reduce `MYSIDE` for weak
scaling or `lx`, `ly`, and `lz` for strong scaling.

## Output files and metrics

Each run writes its log to a file named

```text
mio<model>_<lattice>_<mode>_<N>_<PX>_<PY>_<PZ><precision>.dat
```

where `<model>` is `1` or `2`, `<lattice>` is `15`, `19`, `27`, or `27h`,
`<mode>` is `w` for weak or `s` for strong scaling, `<N>` is the number of
GPUs, and `<precision>` is empty for single, `_d` for double, and `_sd` for
mixed precision.

At the end of the run, rank 0 prints the following quantities:

| Log entry | Meaning |
|---|---|
| `glups:` | Giga lattice updates per second, `lx*ly*lz*nsteps / (t_loop * 1e9)`. |
| `time elapsed:` | `t_loop`, the wall-clock time of the time loop only, excluding setup and final I/O. |
| `J_per_LUP:` | Energy per lattice update in joules. |
| `E_tot:` | Total energy of the time loop in joules, summed over all GPUs. |
| `Pavg:` | Average power in watts, summed over all GPUs. |
| `eff_glups_per_w:` | Throughput per watt. |

The energy entries are available only in NVML builds. Each rank reads the
cumulative energy counter of its own GPU through
`nvmlDeviceGetTotalEnergyConsumption` at the start and at the end of the time
loop, and the differences are summed over all ranks. The published energy
metric is

```text
nJ/LUP = 1e9 * J_per_LUP.
```

The parallel efficiency is computed from the throughput as

```text
PE(N) = GLUPS(N) / (N * GLUPS(1)).
```

A quick tabulation of one series can be obtained with

```bash
for f in mio1_27_w_*.dat; do
  printf '%s ' "$f"; awk '/glups:/ {print $2}' "$f"
done
```

## Correspondence with the published results

| Published result | Model | Mode | Data |
|---|---|---|---|
| README Table 1 | One-component | Strong | `glups:` from `mio1_*_s_*.dat` |
| README Table 2 | One-component | Weak | `glups:` from `mio1_*_w_*.dat` |
| README Table 3 | Two-component | Strong | `glups:` from `mio2_*_s_*.dat` |
| README Table 4 | Two-component | Weak | `glups:` from `mio2_*_w_*.dat` |
| Publication, weak-scaling efficiency figure | One-component D3Q27 | Weak | `PE(N)` from `mio1_27_w_*.dat` |
| Publication, energy figure | One-component | Weak, 3D decompositions | `J_per_LUP:` from `mio1_*_w_*.dat` |

README Tables 1 and 3 together correspond to Table 1 of the publication, and
README Tables 2 and 4 to its Table 2.
