# LBFAST v1.1.2: Changes Since v1.1.1

This document summarizes the changes included in LBFAST `v1.1.2`, released on
18 October 2026, relative to the `v1.1.1` release.

## Release scope

LBFAST `v1.1.2` corrects the documentation of `v1.1.1`. The source code, the
build files, the scripts, and the numerical results are unchanged. The
comparison baseline is tag `v1.1.1` at commit `17df40b`.

## Documentation

- The build procedures for `make` and CMake in the README, `docs/CMAKE.md`,
  and the `v1.1.1` release notes asked to run a script that loads the modules
  of one particular cluster. The instructions now ask to load the NVIDIA HPC
  SDK environment, including its MPI library, as required on each system,
  typically with the module provided by the site (for example
  `module load nvhpc-hpcx`; module names differ between systems).
- The README and `CITATION.cff` cite the archived software release.

## Verification

Only documentation files changed since `v1.1.1` (`git diff --stat v1.1.1
v1.1.2`: `README.md`, `CITATION.cff`, `docs/CMAKE.md`, and the release notes),
so the verification of `v1.1.1` applies unchanged: smoke suite `6 passed,
0 failed`, complete suite `5 passed, 0 failed` at one and two ranks with
summaries identical to `v1.1.0`, and 24 of 24 executables built by both
`compile.sh` and `compile_cmake.sh`. See
[`RELEASE_NOTES_v1.1.1.md`](RELEASE_NOTES_v1.1.1.md).

## Upgrade notes

No action is needed: executables, inputs, and outputs of `v1.1.1` remain valid.
