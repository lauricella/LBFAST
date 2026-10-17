# LBFAST v1.1.1: Changes Since v1.1.0

This document summarizes the changes included in LBFAST `v1.1.1`, released on
18 October 2026, relative to the `v1.1.0` release.

## Release scope

LBFAST `v1.1.1` adds a CMake build alongside the `Makefile`, and lets the
`Makefile` and `compile.sh` find the NVML library in the NVIDIA HPC SDK by
themselves. The solver, its numerical results, the input files, and the output
formats are unchanged: the `Makefile` issues the same compile and link commands
as in `v1.1.0`, and the executables built with `make` and with CMake write
byte-identical results. The comparison baseline is tag `v1.1.0` at commit
`be3acc4`.

## CMake build

`CMakeLists.txt` builds *LBFAST* in a separate build directory with the same
sources and flags as the `Makefile`
([`CMAKE.md`](CMAKE.md); procedure in the README, "Compiling *LBFAST*"):

```bash
source load_iac.sh
cmake -S . -B build/default -DGPUCC=80
cmake --build build/default -j 8
```

- The options `LBFAST_MPI`, `LBFAST_NVML`, `CMAKE_BUILD_TYPE` (`Release` or
  `Debug`), and `GPUCC` correspond to the `make` targets and to `GPUCC=...`.
- The compile-time macros come from `defines.h`, as with `make`, or, with
  `-DLBFAST_USE_DEFINES_H=OFF`, from CMake cache variables named after the
  macros (`-DLATTICE=19 -DTWOCOMPONENT=OFF ...`). Their defaults are the
  values in `defines.h`, which is left untouched, so several configurations
  can be built side by side.
- `CMakePresets.json` contains the 24 executables of `compile.sh`, with the
  same names, and `compile_cmake.sh` builds them in parallel into
  `build/bin/`.
- Objects are rebuilt when `defines.h`, a stage header, or a module they use
  changes, and modules that do not depend on each other compile in parallel.
  CMake reads the `module` and `use` statements of the sources itself, since
  its own Fortran dependency scanner misses `use` statements whose
  continuation lines contain preprocessor directives.

On the reference node (64 cores, NVIDIA HPC SDK 24.3, `GPUCC=80`):

| Build | `make` | CMake |
|---|---:|---:|
| One configuration (`nvfortran-mpi`, two-component high-order D3Q27) | 63 s (serial) | 38 s (`-j 8`) |
| The 24-executable matrix | 1463 s (`compile.sh`) | 232 s (`compile_cmake.sh`) |

## Shared build description

The source lists, the headers, and the release and debug flags are defined
once in `lbfast_build.mk`. The `Makefile` includes it and `CMakeLists.txt`
parses it, so a new source file or a change of flags is made in one place.
The user block of `defines.h` is enclosed in
`#ifndef LBFAST_CMAKE_CONFIG ... #endif`, so that CMake can replace it with
the header it generates; `make` builds are not affected.

`tests/check_build_systems.py` builds the current configuration with both
systems and checks that they compile the same sources with the same compiler
and flags, use the same NVML files, write byte-identical output (including the
final RAW fields and the double-precision restart file), and that the
generated macro header matches `defines.h`. The smoke suite runs it as a
sixth case, `build_systems`.

## Makefile and `compile.sh`

- The NVML targets take `nvml.h` and the `libnvidia-ml` stub from the default
  CUDA toolkit of the HPC SDK, `$NVHPC_ROOT/cuda/include` and
  `$NVHPC_ROOT/cuda/lib64/stubs`, with `NVHPC_ROOT` set by the SDK module or
  derived from the location of `nvfortran`. `NVML_INC` and `NVML_LIB` still
  override them. The build now stops with an error if `nvml.h` or
  `libnvidia-ml.so` is missing; in `v1.1.0` an empty `NVML_INC` turned
  `-I$(NVML_INC) -D_NVML` into `-I -D_NVML`, and `nvml_wrapper.c` was compiled
  without `_NVML`.
- `compile.sh` no longer needs `NVHPC_HOME`, nor contains the SDK and CUDA
  versions; the misspelled `NNVML_LIB` argument is gone. It builds the same 24
  configurations with `make nvfortran-nvml-mpi GPUCC=80`.
- Objects are rebuilt when `defines.h` or a stage header changes.
- `nvfortran-mpi-debug` uses the same flags as `nvfortran-debug`, which adds
  `lineinfo`; the unused variable `FCUDAFLAGSRID` was removed. The commands of
  all other targets are unchanged.

## Verification

| Check | Result |
|---|---|
| Smoke suite | `6 passed, 0 failed`, including `build_systems` |
| Complete suite, 1 rank | `5 passed, 0 failed` |
| Complete suite, 2 ranks (`1x1x2`; `1x2x1` for the capillary wave) | `5 passed, 0 failed` |
| `compile.sh` executable matrix | 24 of 24 executables built, no errors, without `NVHPC_HOME` |
| `compile_cmake.sh` executable matrix | 24 of 24 executables built, no errors |
| `make -n` before and after the `Makefile` changes | same commands, apart from `lineinfo` in `nvfortran-mpi-debug`, the NVML directories, and the compile order of the NVML targets |
| `check_build_systems.py` for the six `make` targets | same sources, flags, and NVML files |
| `compile.sh` and `compile_cmake.sh` executables, `64^3`, 20 steps | byte-identical output for all 24 variants |
| Speed of the two builds (A30, `256^3`) | D3Q19 FP64 2.17 GLUPS, two-component high-order D3Q27 FP64 0.396-0.399 GLUPS, with both |

The summaries of the complete suite, at one and at two ranks, are identical
to those of `v1.1.0` at all printed digits.

## Upgrade notes

Users updating from `v1.1.0` should review the following points:

1. Results, input files, `defines.h` macros, and output and restart formats
   are unchanged; executables built with `make` from `v1.1.0` remain valid.
2. `NVHPC_HOME` is no longer needed by `compile.sh`. If `NVML_INC` or
   `NVML_LIB` were set by hand, they are still used.
3. In `defines.h`, keep the compile-time choices inside the block between
   `#ifndef LBFAST_CMAKE_CONFIG` and `#else`: CMake turns only the macros of
   that block into cache variables.
4. CMake must be run from a separate build directory, for example `build/`,
   which Git ignores; `cmake .` in the main directory is refused, because it
   would replace the `Makefile`.

## Release checklist

Before publishing the `v1.1.1` tag:

- run the smoke suite from a clean build;
- run the complete suite on at least one supported NVIDIA GPU;
- run at least one two- or four-rank MPI validation;
- build the `compile.sh` and `compile_cmake.sh` executable matrices;
- confirm that `defines.h` is restored to the documented default;
- confirm that this document matches the final release commit;
- create the `v1.1.1` tag from the validated `main` commit;
- archive the release and record its DOI.
