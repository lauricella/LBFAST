# Full quantitative validation suite

## Purpose

`tests/run_full_validation.py` runs every complete physical benchmark described
in the LBFAST documentation and converts the numerical comparisons into a
single pass/fail report. Unlike the short regression suite, it uses the
tracked production validation inputs and their full run lengths. A passing
case therefore requires quantitative agreement, not only successful execution
and finite output.

## Included benchmarks

| Case | Driver and full configuration | Quantitative acceptance criteria |
|---|---|---|
| CSF Laplace law | Three-radius `run_laplace_csf_sweep.py` sweep, 10000 steps per radius | Three radii; maximum pointwise error <= 3%; through-origin slope error <= 3%. |
| Planar capillary wave | `analyse_capillary_wave.py`, 80000 steps | Frequency error <= 12%; mass error <= `1e-6`; CSF projection error <= 2%; kinematic RMS error <= 2%. |
| Lamb oscillation | `run_lamb_validation.py`, `80^3`, 20000 steps | At least two peaks and period error <= 10%. |
| Taylor--Green vortex | `run_taylorgreen_validation.py`, `128^3`, 2000 steps | Viscosity error <= 3%; log-energy fit RMS <= `2e-3`. |
| Forced Poiseuille flow | `run_poiseuille_force_validation.py`, 10000 steps | Relative L2 and centerline errors <= 1%; independent analytical profiles agree within `1e-12`. |

These limits surround the current documented results while remaining strict
enough to detect material numerical regressions.

## Running all validations

From the repository root:

```bash
python3 tests/run_full_validation.py
```

The default build uses the `nvfortran` Make target and NVIDIA compute
capability 8.0. Both are configurable:

```bash
python3 tests/run_full_validation.py --gpu-cc 90
python3 tests/run_full_validation.py --make-target nvfortran
```

## Multi-GPU execution

The complete suite can run on two or four MPI ranks with one command:

```bash
python3 tests/run_full_validation.py --mpi-procs 2
python3 tests/run_full_validation.py --mpi-procs 4
```

For more than one rank, the default Make target automatically changes to
`nvfortran-mpi`, and the default launcher is `mpirun`. LBFAST requires
`PX*PY*PZ=N`, and every local dimension must remain a multiple of the
`8 x 8 x 8` GPU tile. The suite therefore selects a compatible decomposition
for each case:

| Case | 2-rank decomposition | 4-rank decomposition |
|---|---:|---:|
| CSF Laplace | `1 x 1 x 2` | `1 x 1 x 4` |
| Planar capillary wave | `1 x 2 x 1` | `2 x 2 x 1` |
| Lamb oscillation | `1 x 1 x 2` | `1 x 2 x 2` |
| Taylor--Green vortex | `1 x 1 x 2` | `1 x 1 x 4` |
| Forced Poiseuille flow | `1 x 1 x 2` | `1 x 1 x 4` |

The capillary wave cannot be divided along `z` because its global thickness is
only eight nodes. The `80^3` Lamb domain cannot be divided by four along a
single direction because the resulting local extent would be 20 rather than a
multiple of eight.

On a Slurm allocation, select `srun` without changing the decomposition logic:

```bash
python3 tests/run_full_validation.py --mpi-procs 4 --launcher srun
```

The individual benchmark drivers expose the same `--mpi-procs`,
`--decomposition PX PY PZ`, and `--launcher` options. For example:

```bash
python3 run_lamb_validation.py \
  --make-target nvfortran-mpi \
  --mpi-procs 4 --decomposition 1 2 2
```

The full suite is intentionally expensive. In particular, it performs three
separate 10000-step Laplace simulations and the 80000-step capillary-wave run,
in addition to the other documented cases.

Individual cases can be selected during development:

```bash
python3 tests/run_full_validation.py --case taylor_green
python3 tests/run_full_validation.py \
  --case laplace_csf --case poiseuille_force
```

Use `--fail-fast` to stop after the first failing case. By default, the runner
continues so that the final report exposes every regression. `--keep-results`
preserves result directories belonging to unselected cases when rerunning a
subset; selected cases are always executed again.

## Outputs and source restoration

Each validation driver writes its normal input copy, simulation logs,
diagnostics, plots, and machine-readable summary under
`tests/full_results/<case-output>/`. The orchestrator additionally records the
driver console output in `driver.log` and writes

```text
tests/full_results/full_validation_results.csv
```

with the case name, status, elapsed time, primary error metric, tolerance, and
detail message. `tests/full_results/` is ignored by Git.

Every underlying driver changes `defines.h` to match its executable. The
orchestrator saves the complete initial file and restores it in a `finally`
block, including when a driver or numerical check fails. It then removes build
objects and `main.x`, preventing a stale executable from remaining beside the
restored base configuration.

## Relationship to the short suite

Use `tests/run_smoke_tests.py` for rapid compilation and initial-physics checks
while editing. Use `tests/run_full_validation.py` before a release or after a
change that can affect numerical accuracy. The short and full suites share the
same five physical cases but serve different purposes and use separate output
directories.
