# Building the LBFAST executable matrix with `compile.sh`

## Purpose

`compile.sh` builds the standard LBFAST executable matrix used for performance
and validation runs. Starting from one known `defines.h` configuration, it
performs clean builds for several lattice and precision combinations and gives
each executable a descriptive name.

The script targets NVIDIA GPUs, CUDA-aware MPI, and NVML. Its default compute
capability is CC 8.0, which covers NVIDIA A30 and A100 GPUs.

## Prerequisites

The NVIDIA HPC SDK environment must already be loaded. In particular,
`NVHPC_HOME` must identify the active HPC SDK installation. The script derives
the CUDA NVML include and stub-library directories from

```text
$NVHPC_HOME/Linux_x86_64/24.3/cuda/12.3/targets/x86_64-linux
```

The exact path reflects the HPC SDK version for which the script is currently
configured and may need adjustment on installations using another version.

The build target is invoked as

```bash
make nvfortran-nvml-mpi \
  GPUCC=80 \
  NVML_INC=... \
  NNVML_LIB=...
```

Consequently, `nvfortran`, `nvcc`, MPI, CUDA, and the NVML headers and library
must be available.

## Required initial `defines.h`

Before performing any build, the script checks that `defines.h` contains all
of the following baseline definitions:

```c
#define LATTICE 27
#define HIGHORDER
#define TWOCOMPONENT
#define DENSRATIO
#define PRC 8
#define noMIXEDPRC
#define STRPRC 8
```

The check accepts flexible whitespace but requires the macro values shown
above. If any definition is missing or inconsistent, the script stops before
modifying the file.

The tracked baseline also keeps individual validation selectors disabled, for
example `noTAYLORGREEN`, `noLAMBTEST`, `noLAPLACE`, and
`noCAPILLARYWAVE`. This ensures that the generated binaries are general
benchmark executables rather than executables specialized for one validation
case.

## Precision naming convention

The generated filenames use the following suffixes:

| Suffix | Arithmetic (`PRC`) | Stored fields (`STRPRC`) | Mode |
|---|---:|---:|---|
| `_d.x` | 8 bytes | 8 bytes | full double precision |
| `_sd.x` | 8 bytes | 4 bytes | double arithmetic, single storage |
| `.x` | 4 bytes | 4 bytes | full single precision |

`MIXEDPRC` is enabled only for the `_sd.x` executables. In the other two
modes, `noMIXEDPRC` makes the storage kind follow `PRC`.

## Lattice naming convention

The central part of each filename identifies the velocity stencil:

| Filename component | Lattice |
|---|---|
| `_15` | D3Q15 |
| `_19` | D3Q19 |
| `_27` | D3Q27 |
| `_27high` | high-order D3Q27 equilibrium |

The prefixes distinguish two- and one-component solvers:

```text
main_2c_...   two-component solver
main_1c_...   single-component solver
```

## Generated executables

The script produces 24 executables: twelve two-component builds followed by
twelve single-component builds.

For each component mode it generates:

```text
main_[1c|2c]_15_d.x
main_[1c|2c]_19_d.x
main_[1c|2c]_27_d.x
main_[1c|2c]_27high_d.x

main_[1c|2c]_15_sd.x
main_[1c|2c]_19_sd.x
main_[1c|2c]_27_sd.x
main_[1c|2c]_27high_sd.x

main_[1c|2c]_15.x
main_[1c|2c]_19.x
main_[1c|2c]_27.x
main_[1c|2c]_27high.x
```

The order in which these files are built differs from this grouped listing,
but the final matrix contains the same combinations.

## How the configuration is advanced

Every executable is built from a clean object-file state:

```bash
make clean
make nvfortran-nvml-mpi ...
mv main.x <descriptive executable name>
```

Between builds, narrowly scoped `sed` substitutions advance `defines.h`
through lattice and precision variants. After completing the two-component
matrix, the script changes

```c
#define TWOCOMPONENT
#define DENSRATIO
#define CSF
```

to

```c
#define noTWOCOMPONENT
#define noDENSRATIO
#define noCSF
```

for the single-component matrix. CSF is explicitly disabled because capillary
forcing belongs to the two-component solver and its kernels are not present in
a single-component build.

After the final executable has been produced, the script performs one last
`make clean` and restores `TWOCOMPONENT`, `DENSRATIO`, and `CSF`. Under normal
completion, `defines.h` therefore returns to the required baseline state.

## Running the script

From the repository root:

```bash
./compile.sh
```

The script uses `set -euo pipefail`, so it stops at the first failed command,
undefined shell variable, or failed pipeline. This prevents subsequent
executables from being labeled as successful after a build error.

Because the script edits `defines.h` progressively, an interrupted or failed
run can leave the file in an intermediate configuration. Before restarting,
inspect it with

```bash
git diff -- defines.h
```

and restore the required baseline definitions deliberately. Do not discard
unrelated local changes. Once the baseline is restored, rerunning `compile.sh`
will repeat the complete matrix from the beginning.

## Relationship to validation scripts

The dedicated validation drivers such as `run_lamb_validation.py` and
`run_taylorgreen_validation.py` build a single executable with benchmark-specific
macros. They are preferable when only one validation case is required.
`compile.sh` is intended for generating the complete standard executable
matrix.
