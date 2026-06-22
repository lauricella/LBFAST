# *LBFAST*

*LBFAST* is a specific-purpose research software for GPU-accelerated lattice
Boltzmann simulations based on a lightweight moment-represented formulation.

The code was developed by:

- Marco Lauricella, IAC-CNR, Rome, Italy
- Andrea Montessori, Roma Tre University, Rome, Italy
- Giorgio Amati, CINECA, Rome, Italy
- Filippo Spiga, NVIDIA Development UK Ltd, United Kingdom
- Adriano Tiribocchi, IAC-CNR, Rome, Italy
- Massimo Bernaschi, IAC-CNR, Rome, Italy
- Sauro Succi, IIT, Rome, Italy 

The software development process has received support from the Italian
Government through the PRIN grant MOBIOS, ID: 2022N4ZNH3,
CUP: F53C24001000006, from CINECA through the ISCRA-B project MIPLAST
HP10BZY7BK, from GNFM-INdAM, and from the European Research Council
through the ERC Proof of Concept Grant No. 101187935, *LBFAST*.

This is an experimental code. The authors accept no responsibility for the
performance of the code or for the correctness of the results.

The code is licensed under the **Non-Commercial Research License – Version 1.0**.
Use, copy, and modification are permitted for research, educational, and other
non-commercial purposes only. Commercial use requires a separate written
agreement with the copyright holders.

## Structure

*LBFAST* is supplied as a main UNIX directory with subdirectories.

All source files are contained in the `source` sub-directory. The `test_cases`
sub-directory contains example input files that can help the user prepare new
simulations. The `execute` sub-directory is intended as the working directory
from which jobs are submitted and output data are collected.

## Compiling *LBFAST*

**Important note:** *LBFAST* requires CUDA Fortran and must be compiled with the
NVIDIA HPC SDK, using the `nvfortran` compiler. The code also uses OpenACC
directives and CUDA-enabled compilation flags. For MPI builds, the MPI compiler
wrappers should be compatible with the NVIDIA HPC SDK and, on multi-GPU systems,
CUDA-aware MPI is recommended.

The `source` sub-directory contains a UNIX `Makefile` that builds the executable
version of the code in single-process, MPI, debug, and NVML-enabled variants.

The GPU architecture can be selected through the `GPUCC` variable. For example,
for NVIDIA A100 GPUs, one can compile the MPI version with:

```bash
make nvfortran-mpi GPUCC=80
```

A list of available targets can be obtained with:

```bash
make help
```

or simply:

```bash
make
```

Available targets are:

| Target | Description |
|---|---|
| `nvfortran` | Compile a single-process GPU version using `nvfortran` and `nvcc`. |
| `nvfortran-debug` | Compile a single-process GPU debug version with debugging and runtime-checking flags. |
| `nvfortran-mpi` | Compile a parallel MPI GPU version using `mpif90` and `mpicc`; this enables the `-DMPI` preprocessor flag. |
| `nvfortran-mpi-debug` | Compile a parallel MPI GPU debug version with debugging and runtime-checking flags. |
| `nvfortran-nvml` | Compile a single-process GPU version with NVIDIA Management Library support enabled through the `_NVML` preprocessor flag. |
| `nvfortran-nvml-mpi` | Compile a parallel MPI GPU version with both `_NVML` and `-DMPI` enabled. |
| `clean` | Remove object files, module files, and intermediate preprocessed files. |
| `clean-all` | Remove object files, module files, executable files, data files, and intermediate preprocessed files. |

The executable produced by the `Makefile` is named `main.x`.

## Scalability

## One-component model

## Performance in Strong Scaling

## Table 1: Performance in Strong Scaling

**Table 1**: MPI decomposition along the \(x, y, z\) axes and GLUPS versus the number of computing GPU devices \(N_{\mathrm{procs}}\) for a one-component fluid simulation in **strong scaling**.  
Results refer to the Taylor-Green benchmark on a fixed cubic box of side 512 lattice points.

| GPU | Nprocs | MPI decomp. | D3Q19 SP | D3Q27 SP | D3Q27h SP | D3Q19 DP | D3Q27 DP | D3Q27h DP |
|---|---:|---|---:|---:|---:|---:|---:|---:|
| 🔵 NVIDIA A100 64GB | 1   | 1 × 1 × 1   | 2.3  | 1.8  | 1.0  | 1.4  | 1.0  | 0.4  |
| 🔵 NVIDIA A100 64GB | 2   | 1 × 1 × 2   | 4.4  | 3.5  | 1.9  | 2.6  | 2.0  | 0.8  |
| 🔵 NVIDIA A100 64GB | 4   | 1 × 1 × 4   | 8.2  | 6.6  | 3.6  | 5.0  | 3.9  | 1.6  |
| 🔵 NVIDIA A100 64GB | 4   | 1 × 2 × 2   | 8.2  | 6.6  | 3.6  | 4.9  | 3.9  | 1.6  |
| 🔵 NVIDIA A100 64GB | 8   | 1 × 2 × 4   | 14.1 | 11.6 | 6.8  | 8.7  | 7.0  | 3.1  |
| 🔵 NVIDIA A100 64GB | 8   | 2 × 2 × 2   | 13.8 | 11.4 | 6.7  | 8.6  | 6.9  | 3.0  |
| 🔵 NVIDIA A100 64GB | 16  | 1 × 4 × 4   | 22.3 | 19.1 | 12.0 | 14.7 | 12.1 | 5.8  |
| 🔵 NVIDIA A100 64GB | 16  | 2 × 2 × 4   | 22.2 | 18.8 | 11.9 | 14.5 | 12.0 | 5.7  |
| 🔵 NVIDIA A100 64GB | 32  | 1 × 4 × 8   | 31.5 | 28.2 | 19.6 | 21.4 | 18.5 | 10.0 |
| 🔵 NVIDIA A100 64GB | 32  | 2 × 4 × 4   | 32.4 | 28.5 | 19.8 | 22.3 | 19.2 | 10.2 |
| 🔵 NVIDIA A100 64GB | 64  | 1 × 8 × 8   | 41.7 | 39.0 | 29.0 | 29.9 | 26.4 | 16.8 |
| 🔵 NVIDIA A100 64GB | 64  | 4 × 4 × 4   | 39.0 | 40.1 | 27.7 | 30.2 | 26.3 | 17.2 |
| 🔵 NVIDIA A100 64GB | 128 | 1 × 8 × 16  | 49.5 | 47.6 | 41.5 | 38.5 | 35.4 | 26.4 |
| 🔵 NVIDIA A100 64GB | 128 | 4 × 4 × 8   | 41.4 | 42.2 | 40.4 | 37.9 | 33.8 | 26.5 |
| 🔵 NVIDIA A100 64GB | 256 | 1 × 16 × 16 | 44.8 | 58.5 | 47.6 | 46.1 | 43.9 | 34.3 |
| 🔵 NVIDIA A100 64GB | 256 | 4 × 8 × 8   | 44.3 | 48.5 | 45.1 | 45.1 | 40.7 | 33.8 |
| 🔵 NVIDIA A100 64GB | 512 | 1 × 16 × 32 | 62.3 | 55.9 | 55.1 | 56.1 | 53.9 | 46.9 |
| 🔵 NVIDIA A100 64GB | 512 | 8 × 8 × 8   | 48.1 | 46.0 | 44.5 | 48.3 | 47.4 | 42.7 |

## Performance in Weak Scaling

## Table 2: Performance in Weak Scaling

**Table 2**: MPI decomposition along the \(x, y, z\) axes and GLUPS versus the number of computing GPU devices \(N_{\mathrm{procs}}\) for a one-component fluid simulation in **weak scaling**.  
Results refer to the Taylor-Green benchmark with a fixed cubic sub-domain of side 512 lattice points assigned to each GPU device.

| GPU | Nprocs | MPI decomp. | D3Q19 SP | D3Q27 SP | D3Q27h SP | D3Q19 DP | D3Q27 DP | D3Q27h DP |
|---|---:|---|---:|---:|---:|---:|---:|---:|
| 🔵 NVIDIA A100 64GB | 1   | 1 × 1 × 1     | 2.3    | 1.8   | 1.0   | 1.3   | 1.0   | 0.4   |
| 🔵 NVIDIA A100 64GB | 2   | 1 × 1 × 2     | 4.5    | 3.6   | 1.9   | 2.6   | 2.1   | 0.8   |
| 🔵 NVIDIA A100 64GB | 4   | 1 × 1 × 4     | 9.0    | 7.1   | 3.8   | 5.3   | 4.1   | 1.7   |
| 🔵 NVIDIA A100 64GB | 4   | 1 × 2 × 2     | 9.0    | 7.0   | 3.8   | 5.2   | 4.1   | 1.6   |
| 🔵 NVIDIA A100 64GB | 8   | 1 × 1 × 8     | 17.7   | 14.0  | 7.5   | 10.3  | 8.1   | 3.3   |
| 🔵 NVIDIA A100 64GB | 8   | 1 × 2 × 4     | 17.4   | 13.8  | 7.4   | 10.2  | 8.0   | 3.3   |
| 🔵 NVIDIA A100 64GB | 8   | 2 × 2 × 2     | 17.4   | 13.8  | 7.4   | 10.1  | 8.0   | 3.3   |
| 🔵 NVIDIA A100 64GB | 16  | 1 × 1 × 16    | 35.3   | 27.8  | 14.9  | 20.7  | 16.2  | 6.6   |
| 🔵 NVIDIA A100 64GB | 16  | 1 × 4 × 4     | 34.7   | 27.6  | 14.8  | 20.2  | 15.9  | 6.5   |
| 🔵 NVIDIA A100 64GB | 16  | 2 × 2 × 4     | 33.6   | 26.7  | 14.6  | 19.4  | 15.4  | 6.4   |
| 🔵 NVIDIA A100 64GB | 32  | 1 × 1 × 32    | 70.9   | 55.9  | 29.9  | 41.2  | 32.4  | 13.1  |
| 🔵 NVIDIA A100 64GB | 32  | 1 × 4 × 8     | 67.9   | 54.2  | 29.5  | 39.3  | 31.2  | 12.7  |
| 🔵 NVIDIA A100 64GB | 32  | 2 × 4 × 4     | 66.7   | 53.3  | 29.2  | 38.4  | 30.6  | 12.8  |
| 🔵 NVIDIA A100 64GB | 64  | 1 × 1 × 64    | 140.7  | 111.7 | 59.9  | 82.3  | 64.5  | 26.0  |
| 🔵 NVIDIA A100 64GB | 64  | 1 × 8 × 8     | 134.8  | 107.7 | 58.7  | 78.4  | 62.3  | 25.7  |
| 🔵 NVIDIA A100 64GB | 64  | 4 × 4 × 4     | 132.2  | 106.2 | 58.3  | 76.7  | 60.8  | 25.3  |
| 🔵 NVIDIA A100 64GB | 128 | 1 × 1 × 128   | 278.8  | 223.2 | 119.6 | 164.8 | 128.9 | 52.1  |
| 🔵 NVIDIA A100 64GB | 128 | 1 × 8 × 16    | 268.0  | 215.9 | 117.7 | 157.2 | 124.4 | 51.5  |
| 🔵 NVIDIA A100 64GB | 128 | 4 × 4 × 8     | 256.8  | 208.7 | 115.2 | 150.2 | 119.9 | 50.7  |
| 🔵 NVIDIA A100 64GB | 256 | 1 × 1 × 256   | 564.3  | 447.5 | 239.5 | 329.8 | 258.7 | 104.4 |
| 🔵 NVIDIA A100 64GB | 256 | 1 × 16 × 16   | 536.0  | 428.7 | 234.7 | 311.1 | 247.5 | 102.4 |
| 🔵 NVIDIA A100 64GB | 256 | 4 × 8 × 8     | 512.6  | 413.1 | 230.3 | 296.1 | 237.5 | 101.1 |
| 🔵 NVIDIA A100 64GB | 512 | 1 × 1 × 512   | 1131.4 | 890.2 | 479.2 | 658.9 | 517.3 | 209.4 |
| 🔵 NVIDIA A100 64GB | 512 | 1 × 16 × 32   | 1001.5 | 797.7 | 430.5 | 582.2 | 457.3 | 191.5 |
| 🔵 NVIDIA A100 64GB | 512 | 8 × 8 × 8     | 960.5  | 770.7 | 421.1 | 556.2 | 440.6 | 187.1 |

## Two-component model

## Performance in Strong Scaling

## Table 3: Performance in Strong Scaling

**Table 3**: MPI decomposition along the \(x, y, z\) axes and GLUPS versus the number of computing GPU devices \(N_{\mathrm{procs}}\) for a two-component fluid simulation in **strong scaling**.  
Results refer to the Laplace benchmark on a fixed cubic box of side 512 lattice points.

| GPU | Nprocs | MPI decomp. | D3Q19 SP | D3Q27 SP | D3Q27h SP | D3Q19 DP | D3Q27 DP | D3Q27h DP |
|---|---:|---|---:|---:|---:|---:|---:|---:|
| 🔵 NVIDIA A100 64GB | 1   | 1 × 1 × 1   | 1.6  | 1.4  | 0.8  | 0.8  | 0.7  | 0.2  |
| 🔵 NVIDIA A100 64GB | 2   | 1 × 1 × 2   | 3.0  | 2.6  | 1.6  | 1.5  | 1.2  | 0.5  |
| 🔵 NVIDIA A100 64GB | 4   | 1 × 1 × 4   | 5.7  | 4.8  | 3.0  | 2.9  | 2.5  | 0.9  |
| 🔵 NVIDIA A100 64GB | 4   | 1 × 2 × 2   | 5.6  | 4.8  | 3.0  | 2.9  | 2.4  | 0.9  |
| 🔵 NVIDIA A100 64GB | 8   | 1 × 2 × 4   | 9.3  | 8.2  | 5.4  | 5.1  | 4.3  | 1.8  |
| 🔵 NVIDIA A100 64GB | 8   | 2 × 2 × 2   | 9.1  | 8.0  | 5.3  | 5.0  | 4.3  | 1.7  |
| 🔵 NVIDIA A100 64GB | 16  | 1 × 4 × 4   | 14.1 | 12.6 | 9.1  | 8.4  | 7.4  | 3.3  |
| 🔵 NVIDIA A100 64GB | 16  | 2 × 2 × 4   | 13.5 | 12.4 | 8.9  | 8.4  | 7.2  | 3.2  |
| 🔵 NVIDIA A100 64GB | 32  | 1 × 4 × 8   | 18.7 | 16.9 | 13.3 | 11.8 | 10.8 | 5.6  |
| 🔵 NVIDIA A100 64GB | 32  | 2 × 4 × 4   | 18.6 | 17.2 | 13.4 | 12.3 | 11.1 | 5.7  |
| 🔵 NVIDIA A100 64GB | 64  | 1 × 8 × 8   | 21.9 | 21.6 | 14.9 | 15.7 | 14.5 | 9.1  |
| 🔵 NVIDIA A100 64GB | 64  | 4 × 4 × 4   | 19.7 | 20.5 | 17.5 | 14.4 | 14.5 | 9.2  |
| 🔵 NVIDIA A100 64GB | 128 | 1 × 8 × 16  | 24.4 | 23.6 | 23.4 | 19.3 | 16.5 | 13.2 |
| 🔵 NVIDIA A100 64GB | 128 | 4 × 4 × 8   | 19.8 | 22.0 | 21.6 | 17.8 | 14.8 | 13.1 |
| 🔵 NVIDIA A100 64GB | 256 | 1 × 16 × 16 | 21.6 | 26.1 | 24.7 | 19.5 | 21.3 | 17.0 |
| 🔵 NVIDIA A100 64GB | 256 | 4 × 8 × 8   | 19.1 | 23.9 | 19.0 | 18.1 | 18.0 | 16.8 |
| 🔵 NVIDIA A100 64GB | 512 | 1 × 16 × 32 | 27.3 | 24.8 | 24.5 | 23.4 | 23.2 | 20.7 |
| 🔵 NVIDIA A100 64GB | 512 | 8 × 8 × 8   | 21.8 | 19.4 | 21.5 | 20.4 | 19.0 | 20.5 |

## Performance in Weak Scaling

## Table 4: Performance in Weak Scaling

**Table 4**: MPI decomposition along the \(x, y, z\) axes and GLUPS versus the number of computing GPU devices \(N_{\mathrm{procs}}\) for a two-component fluid simulation in **weak scaling**.  
Results refer to the Laplace benchmark with a fixed cubic sub-domain of side 512 lattice points assigned to each GPU device.

| GPU | Nprocs | MPI decomp. | D3Q19 SP | D3Q27 SP | D3Q27h SP | D3Q19 DP | D3Q27 DP | D3Q27h DP |
|---|---:|---|---:|---:|---:|---:|---:|---:|
| 🔵 NVIDIA A100 64GB | 1   | 1 × 1 × 1     | 1.6   | 1.4   | 0.8   | 0.8   | 0.7   | 0.2   |
| 🔵 NVIDIA A100 64GB | 2   | 1 × 1 × 2     | 3.2   | 2.7   | 1.6   | 1.6   | 1.3   | 0.5   |
| 🔵 NVIDIA A100 64GB | 4   | 1 × 1 × 4     | 6.4   | 5.3   | 3.2   | 3.1   | 2.6   | 0.9   |
| 🔵 NVIDIA A100 64GB | 4   | 1 × 2 × 2     | 6.3   | 5.3   | 3.2   | 3.1   | 2.6   | 0.9   |
| 🔵 NVIDIA A100 64GB | 8   | 1 × 1 × 8     | 12.5  | 10.5  | 6.3   | 6.2   | 5.1   | 1.9   |
| 🔵 NVIDIA A100 64GB | 8   | 1 × 2 × 4     | 12.2  | 10.2  | 6.3   | 6.0   | 5.0   | 1.8   |
| 🔵 NVIDIA A100 64GB | 8   | 2 × 2 × 2     | 12.1  | 10.2  | 6.3   | 6.0   | 5.0   | 1.8   |
| 🔵 NVIDIA A100 64GB | 16  | 1 × 1 × 16    | 24.9  | 20.8  | 12.7  | 12.3  | 10.2  | 3.7   |
| 🔵 NVIDIA A100 64GB | 16  | 1 × 4 × 4     | 24.3  | 20.3  | 12.5  | 12.0  | 10.0  | 3.6   |
| 🔵 NVIDIA A100 64GB | 16  | 2 × 2 × 4     | 23.3  | 19.8  | 12.3  | 11.5  | 9.6   | 3.6   |
| 🔵 NVIDIA A100 64GB | 32  | 1 × 1 × 32    | 49.6  | 41.6  | 25.2  | 24.6  | 20.3  | 7.4   |
| 🔵 NVIDIA A100 64GB | 32  | 1 × 4 × 8     | 47.3  | 39.9  | 24.7  | 23.4  | 19.5  | 7.1   |
| 🔵 NVIDIA A100 64GB | 32  | 2 × 4 × 4     | 46.1  | 39.1  | 24.3  | 22.9  | 19.2  | 7.1   |
| 🔵 NVIDIA A100 64GB | 64  | 1 × 1 × 64    | 98.3  | 83.6  | 50.6  | 48.9  | 40.8  | 14.8  |
| 🔵 NVIDIA A100 64GB | 64  | 1 × 8 × 8     | 93.8  | 79.6  | 49.3  | 46.7  | 39.1  | 14.3  |
| 🔵 NVIDIA A100 64GB | 64  | 4 × 4 × 4     | 91.2  | 78.2  | 48.6  | 45.7  | 37.9  | 14.1  |
| 🔵 NVIDIA A100 64GB | 128 | 1 × 1 × 128   | 196.7 | 165.9 | 100.8 | 97.8  | 81.0  | 29.7  |
| 🔵 NVIDIA A100 64GB | 128 | 1 × 8 × 16    | 187.4 | 157.9 | 98.3  | 93.0  | 77.2  | 28.5  |
| 🔵 NVIDIA A100 64GB | 128 | 4 × 4 × 8     | 179.2 | 152.4 | 95.9  | 88.8  | 74.8  | 28.3  |
| 🔵 NVIDIA A100 64GB | 256 | 1 × 1 × 256   | 393.8 | 332.2 | 201.9 | 195.1 | 162.0 | 59.3  |
| 🔵 NVIDIA A100 64GB | 256 | 1 × 16 × 16   | 370.9 | 317.6 | 196.2 | 184.9 | 154.4 | 56.9  |
| 🔵 NVIDIA A100 64GB | 256 | 4 × 8 × 8     | 351.7 | 303.1 | 191.2 | 176.7 | 148.7 | 57.5  |
| 🔵 NVIDIA A100 64GB | 512 | 1 × 1 × 512   | 767.8 | 641.1 | 386.3 | 379.3 | 312.6 | 115.6 |
| 🔵 NVIDIA A100 64GB | 512 | 1 × 16 × 32   | 709.3 | 596.2 | 362.7 | 351.5 | 290.8 | 111.2 |
| 🔵 NVIDIA A100 64GB | 512 | 8 × 8 × 8     | 674.2 | 566.5 | 353.4 | 334.3 | 277.2 | 108.2 |


