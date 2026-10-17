# Building LBFAST with CMake

## Purpose

`CMakeLists.txt` is an alternative to the `Makefile`. Both build systems
compile the same sources with the same flags, so they produce the same
executable. CMake adds

- **out-of-source builds**: each configuration lives in its own build
  directory, so several lattice and precision variants can be built, and
  kept, side by side without touching the sources;
- **compile-time macros without editing `defines.h`**: they can be given as
  CMake cache variables, so builds and test suites no longer interfere through
  a shared `defines.h`;
- **dependency tracking**: objects are rebuilt when `defines.h`, a stage
  header or a module they use changes, so `make clean` is not needed;
- **parallel compilation** of the modules that do not depend on each other.

The `Makefile` remains the reference build. Both read the source lists and
flags from `lbfast_build.mk`, and the smoke suite checks that they agree
(see [Keeping the two build systems in step](#keeping-the-two-build-systems-in-step)).

## Requirements

- CMake 3.20 or newer (3.21 for the presets);
- the NVIDIA HPC SDK compilers `nvfortran` and `nvc`. Without `FC` and `CC`
  in the environment or on the command line, CMake selects them by name;
- for MPI builds, an MPI library built for the HPC SDK (on the reference
  system the HPC-X library of the `nvhpc-hpcx` module), found through CMake's
  `FindMPI`;
- for NVML builds, `nvml.h` and `libnvidia-ml`. CMake looks for them in the
  CUDA toolkits shipped with the HPC SDK of the selected compiler, so no path
  has to be given.

On the reference system:

```bash
source load_iac.sh
```

## Quick start

From the repository root:

```bash
cmake -S . -B build/default -DGPUCC=80
cmake --build build/default -j 8
```

or, equivalently, with the `default` preset:

```bash
cmake --preset default
cmake --build --preset default -j 8
```

The executable is `build/default/main.x`. It is run as described in the main
README, for example

```bash
mpirun -np 4 ./build/default/main.x 1 1 4 input.inp
```

An in-source build (`cmake .` in the repository root) is refused, because the
generated `Makefile` would replace the LBFAST one.

## Options

| Option | Default | Meaning | `make` equivalent |
|---|---|---|---|
| `CMAKE_BUILD_TYPE` | `Release` | `Release` or `Debug` flags of `lbfast_build.mk`. | `nvfortran*` / `nvfortran*-debug` |
| `GPUCC` | `all` | GPU compute capability, passed as `-gpu=cc<GPUCC>`. | `GPUCC=...` |
| `LBFAST_MPI` | `ON` | MPI build (`-DMPI`, MPI library). | `nvfortran-mpi` |
| `LBFAST_NVML` | `OFF` | NVML power monitoring (`-D_NVML`, `nvml_wrapper.c`, `libnvidia-ml`). | `nvfortran-nvml` |
| `LBFAST_USE_DEFINES_H` | `ON` | Take the macros from `defines.h`; `OFF` takes them from cache variables. | always `defines.h` |
| `LBFAST_EXTRA_DEFINES` | empty | Further macros for every source, as a list, e.g. `"BOUNCE_BACK;GETPOWER"`. | edit `defines.h` |
| `LBFAST_EXE_NAME` | `main.x` | Name of the executable. | `main.x` |
| `LBFAST_OUTPUT_DIR` | build directory | Directory of the executable. | repository root |
| `NVML_INCLUDE_DIR`, `NVML_LIBRARY` | found | Location of `nvml.h` and `libnvidia-ml`, to override the search. | `NVML_INC=...`, `NVML_LIB=...` |

The `make` targets therefore correspond to

| `make` target | CMake options |
|---|---|
| `nvfortran` | `-DLBFAST_MPI=OFF` |
| `nvfortran-mpi` | (defaults) |
| `nvfortran-nvml` | `-DLBFAST_MPI=OFF -DLBFAST_NVML=ON` |
| `nvfortran-nvml-mpi` | `-DLBFAST_NVML=ON` |
| `nvfortran-debug` | `-DLBFAST_MPI=OFF -DCMAKE_BUILD_TYPE=Debug` |
| `nvfortran-mpi-debug` | `-DCMAKE_BUILD_TYPE=Debug` |

An option of an existing build directory is changed by configuring it again,
e.g. `cmake build/default -DGPUCC=90`.

## Compile-time macros

### From `defines.h` (default)

With `LBFAST_USE_DEFINES_H=ON` the sources read `defines.h` exactly as in a
`make` build. After an edit of `defines.h`, `cmake --build` recompiles every
object.

### From cache variables

With `LBFAST_USE_DEFINES_H=OFF` every macro of the user block of `defines.h`
becomes a cache variable of the same name:

- `#define NAME` becomes an option `NAME=ON`, and `#define noNAME` an option
  `NAME=OFF`;
- `#define NAME value` becomes a string variable `NAME=value`.

For example, a one-component D3Q19 double-precision build:

```bash
cmake -S . -B build/q19 -DGPUCC=80 -DLBFAST_USE_DEFINES_H=OFF \
      -DLATTICE=19 -DHIGHORDER=OFF -DTWOCOMPONENT=OFF -DDENSRATIO=OFF -DCSF=OFF
cmake --build build/q19 -j 8
```

CMake writes the chosen values, in the layout of the user block, to
`lbfast_config.h` in the build directory. The user block of `defines.h` is
enclosed in `#ifndef LBFAST_CMAKE_CONFIG ... #else #include "lbfast_config.h"
#endif`, and these builds define `LBFAST_CMAKE_CONFIG`. The derived macros
that follow the block (for instance `LAPLACE` selecting `TWOCOMPONENT`, or
`FUSED_MOMENTS`) are applied in both modes.

The macros not given on the command line take their value from `defines.h`
when the build directory is configured for the first time. Later edits of
`defines.h` do not change an existing cache: change the cache variable
(`cmake build/q19 -DPRC=4`, or `ccmake build/q19`), or start again with
`cmake --fresh` or a new build directory.

Passing a macro such as `-DLATTICE=19` while `LBFAST_USE_DEFINES_H=ON` is an
error, since it would otherwise be ignored silently.

Macros that are not in the user block, such as `BOUNCE_BACK`, `MULTIHIT`,
`GETPOWER` or `WALLCLOCK`, are given in both modes with
`LBFAST_EXTRA_DEFINES`.

## The executable matrix

`CMakePresets.json` contains one preset for each executable built by
`compile.sh`, with the same name and the same settings
(`nvfortran-nvml-mpi`, `GPUCC=80`):

| Preset part | Macros |
|---|---|
| `1c_` / `2c_` | `TWOCOMPONENT`, `DENSRATIO`, `CSF` off / on |
| `15`, `19`, `27` | `LATTICE`, `HIGHORDER=OFF` |
| `27high` | `LATTICE=27`, `HIGHORDER=ON` |
| `_d` | `PRC=8`, `STRPRC=8`, `MIXEDPRC=OFF` |
| `_sd` | `PRC=8`, `STRPRC=4`, `MIXEDPRC=ON` |
| no suffix | `PRC=4`, `STRPRC=4`, `MIXEDPRC=OFF` |

As with `compile.sh`, the other macros come from `defines.h`. A single
variant is built with

```bash
cmake --preset 2c_27high_d
cmake --build --preset 2c_27high_d -j 8
```

and the whole matrix with

```bash
./compile_cmake.sh                # all 24 presets, 4 at a time, 4 jobs each
./compile_cmake.sh -p 6 -j 4      # 6 at a time
./compile_cmake.sh 2c_19_d 1c_19_d -- -DGPUCC=90
```

Each preset is configured in `build/<preset>` and its executable written to
`build/bin/main_<preset>.x`, e.g. `build/bin/main_2c_27high_d.x`; the logs are
in `build/logs/`. Unlike `compile.sh`, `compile_cmake.sh` does not modify
`defines.h` and can be interrupted and resumed: a second run only rebuilds
what changed.

Build times on the reference node (64 cores, NVIDIA HPC SDK 24.3, `GPUCC=80`):

| Build | `make` | CMake |
|---|---:|---:|
| One configuration (`nvfortran-mpi`, D3Q27 high-order two-component) | 63 s (serial) | 38 s (`-j 8`) |
| The 24-executable matrix | 1463 s (`compile.sh`) | 232 s (`compile_cmake.sh`, 4 variants x 4 jobs) |

The executables are the same: for each of the 24 variants the `compile.sh`
and `compile_cmake.sh` executables write byte-identical restart, RAW and
diagnostic files on a `64^3` run, and run at the same speed (A30, `256^3`:
D3Q19 FP64 2.17 GLUPS with both, two-component D3Q27 high-order FP64
0.396-0.399 GLUPS with both).

## Keeping the two build systems in step

**Shared lists and flags.** `lbfast_build.mk` defines the source lists, the
headers, and the release and debug flags. The `Makefile` includes it, and
`CMakeLists.txt` parses it. A new source file, or a change of optimisation or
GPU flags, is made there once. The file must keep a simple format: one
`NAME = value` assignment per variable, optionally continued with a trailing
backslash, comments on lines of their own, and no make functions or variable
references.

**Flag groups.** As in the `Makefile`, every Fortran source gets the common
flags; the CUDA Fortran sources, recognised by their `lb_cuda_*.f90` name, add
`-cuda`, and the other sources add `-DMPI` in MPI builds.

**Module order.** `LBFAST_F_SOURCES` must list the Fortran sources in a valid
build order, since the `Makefile` compiles them in that order. CMake reads the
`module` and `use` statements of every source and builds each source as an
object library that depends on the sources whose modules it uses; sources
that do not depend on each other are compiled in parallel. The search
ignores `#ifdef`, so it may add an unneeded dependency but never misses one.
If a source uses a module defined in a source listed after it, CMake stops
with an error naming both files. (CMake's own Fortran dependency scanner is
not used: it misses `use` statements whose continuation lines contain
preprocessor directives, which are common in LBFAST, and ignores per-source
macros such as `-DMPI`.)

**Automatic check.** `tests/check_build_systems.py` builds the current
`defines.h` configuration with `make` and with CMake and checks that

1. both compile the same source files;
2. every source and the link step get the same compiler and flags, apart from
   include paths and module, output and library options;
3. with `--input`, both executables write byte-identical output, including the
   final RAW fields and the double-precision restart file, which the check
   enables in its copy of the input;
4. a configuration with `LBFAST_USE_DEFINES_H=OFF` and no macro given writes a
   `lbfast_config.h` identical to the user block of `defines.h`.

```bash
python3 tests/check_build_systems.py --make-target nvfortran-mpi
```

The smoke suite runs it as the `build_systems` case on the Laplace
configuration (see [`SMOKE_TEST_SUITE.md`](SMOKE_TEST_SUITE.md)).

## Differences from the `Makefile`

- `make nvfortran` compiles the two C helpers `get_mem.c` and `get_ram.c`
  with `nvcc`, CMake with `nvc` (as `make nvfortran-mpi` does through
  `mpicc`). They only query the host memory.
- MPI builds call `nvfortran` with the include and library flags of the MPI
  library, found by `FindMPI`, instead of the `mpif90` wrapper; the compiler
  and its flags are the same. `-DCMAKE_Fortran_COMPILER=mpif90
  -DCMAKE_C_COMPILER=mpicc` selects the wrappers instead.
- The `Makefile` takes the NVML location from `NVML_INC` and `NVML_LIB`, in
  the environment or on the command line (`compile.sh` gets them from the
  environment). If `NVML_INC` is empty, `-I$(NVML_INC) -D_NVML` becomes
  `-I -D_NVML` and `nvml_wrapper.c` is compiled without `_NVML`. CMake finds
  NVML in the CUDA toolkit of the HPC SDK of the compiler
  (`$NVHPC_ROOT/cuda/<version>/targets/x86_64-linux`, `include/nvml.h` and
  `lib/stubs/libnvidia-ml.so`). The library there is a stub for linking only:
  both build systems link it with `-L<dir> -lnvidia-ml`, without adding the
  `stubs` directory to the run-time path, so that the executable loads the
  NVML library and the CUDA driver installed with the GPU driver.
- The `Makefile` has no dependencies between modules: after editing a module
  used by other sources, run `make clean`. CMake rebuilds the dependent
  objects itself.

## Troubleshooting

| Message | Cause and remedy |
|---|---|
| `LBFAST needs nvfortran from the NVIDIA HPC SDK` | Another Fortran compiler was found, typically because the HPC SDK module was not loaded. Load it and configure a new build directory. |
| `Could NOT find MPI` | No MPI library for the HPC SDK in the environment. Load the `nvhpc-hpcx` module or configure with `-DLBFAST_MPI=OFF`. |
| `NVML not found` | Set `NVML_INCLUDE_DIR` and `NVML_LIBRARY`, or configure with `-DLBFAST_NVML=OFF`. |
| `Macros set as cache variables (...) are ignored` | A macro was given while `LBFAST_USE_DEFINES_H=ON`. Add `-DLBFAST_USE_DEFINES_H=OFF`, or use a new build directory to go back to `defines.h`. |
| `... uses module ..., defined in ..., which comes after it` | Fix the order of `LBFAST_F_SOURCES` in `lbfast_build.mk`. |
