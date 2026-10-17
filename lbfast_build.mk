# Source lists and compiler flags shared by Makefile and CMakeLists.txt.
#
# Add new source files and change optimisation or GPU flags here, so that
# both build systems stay in step. CMakeLists.txt parses this file itself,
# so keep the format simple: one "NAME = value" assignment per variable,
# optionally continued with a trailing backslash, comments only on lines of
# their own, and no make functions or references to other variables.

# C helpers, built in every configuration.
LBFAST_C_SOURCES = get_mem.c get_ram.c

# NVML power monitoring, built only with -D_NVML (make nvfortran-nvml*,
# cmake -DLBFAST_NVML=ON).
LBFAST_NVML_C_SOURCES = nvml_wrapper.c
LBFAST_NVML_F_SOURCES = nvml_interface_module.f90

# Fortran sources in module dependency order: the Makefile has no
# inter-module dependencies and compiles them in this order. Files named
# lb_cuda_*.f90 are CUDA Fortran and get LBFAST_CUDAFLAGS on top of the
# common flags; the others get -DMPI in MPI builds.
LBFAST_F_SOURCES = vars_module.f90 mpi_module.f90 profiling_m.f90 \
	lb_cuda_vars_module.f90 lb_cuda_auxfields_module.f90 \
	lb_cuda_repulsive_module.f90 lb_cuda_moments_module.f90 \
	lb_cuda_fused_module.f90 lb_cuda_update_phi_module.f90 \
	lb_cuda_boundary_module.f90 lb_cuda_driver_module.f90 \
	boundary_cds_module.f90 statistics.f90 init_conditions_module.f90 \
	print_module.f90 allocate_module.f90 integrator_module.f90 LBFAST.f90

# Headers included by the sources (defines.h also by nvml_wrapper.c): every
# Fortran object, and the NVML wrapper, is rebuilt when one of them changes.
LBFAST_HEADERS = defines.h fused2_stages.h fusedh_stages.h

# Flags common to all Fortran sources, and the extra flag of lb_cuda_*.f90.
LBFAST_FFLAGS = -cpp -acc
LBFAST_CUDAFLAGS = -cuda
LBFAST_INFO = -Minfo=accel

# Link flags; the GPU sub-options are appended to -gpu=cc<GPUCC>.
LBFAST_LDFLAGS = -acc -cuda
LBFAST_LINK_GPU_OPTS = lineinfo

# Release builds (make nvfortran, nvfortran-mpi, nvfortran-nvml*;
# cmake -DCMAKE_BUILD_TYPE=Release).
LBFAST_RELEASE_OPT = -O3
LBFAST_RELEASE_LINK_OPT = -O3
LBFAST_RELEASE_GPU_OPTS = lineinfo,ptxinfo,maxregcount:64

# Debug builds (make nvfortran-debug, nvfortran-mpi-debug;
# cmake -DCMAKE_BUILD_TYPE=Debug).
LBFAST_DEBUG_OPT = -O0 -g
LBFAST_DEBUG_LINK_OPT = -O0
LBFAST_DEBUG_GPU_OPTS = debug,keep,ptxinfo,lineinfo,maxregcount:64
LBFAST_DEBUG_CHECKS = -Mchkptr -Mchkstk -traceback
