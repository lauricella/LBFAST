#
# makefile
#
# Source lists and compiler flags are defined in lbfast_build.mk, which is
# shared with CMakeLists.txt (see docs/CMAKE.md). Add new source files and
# change optimisation or GPU flags there, so that both build systems stay
# in step.
#

include lbfast_build.mk

BINROOT = $(CURDIR)
FC=undefined
CC=undefined
CFLAGSNVML=undefined
NVML_INC ?= $(NVHPC_HOME)
NVML_LIB ?= $(NVHPC_HOME)
FCUDAFLAGS=undefined
FFLAGS=undefined
LDFLAGS=undefined
TYPE=undefined
EX=main.x
EXP=main.x
EXE = $(BINROOT)/$(EX)
EXEP = $(BINROOT)/$(EXP)
SHELL=/bin/sh

# default (change calling: make GPUCC=80)
GPUCC ?= all

# Flag sets composed from lbfast_build.mk
FREL  = $(LBFAST_RELEASE_OPT) $(LBFAST_FFLAGS)
GREL  = -gpu=cc$(GPUCC),$(LBFAST_RELEASE_GPU_OPTS) $(LBFAST_INFO)
LDREL = $(LBFAST_RELEASE_LINK_OPT) $(LBFAST_LDFLAGS) -gpu=cc$(GPUCC),$(LBFAST_LINK_GPU_OPTS) $(LBFAST_INFO)
FDBG  = $(LBFAST_DEBUG_OPT) $(LBFAST_FFLAGS)
GDBG  = -gpu=cc$(GPUCC),$(LBFAST_DEBUG_GPU_OPTS) $(LBFAST_INFO) $(LBFAST_DEBUG_CHECKS)
LDDBG = $(LBFAST_DEBUG_LINK_OPT) $(LBFAST_LDFLAGS) -gpu=cc$(GPUCC),$(LBFAST_LINK_GPU_OPTS) $(LBFAST_INFO) $(LBFAST_DEBUG_CHECKS)

# Object lists from lbfast_build.mk; lb_cuda_*.o are compiled with FCUDAFLAGS
C_OBJS      = $(LBFAST_C_SOURCES:.c=.o)
NVML_C_OBJS = $(LBFAST_NVML_C_SOURCES:.c=.o)
NVML_F_OBJS = $(LBFAST_NVML_F_SOURCES:.f90=.o)
F_OBJS      = $(LBFAST_F_SOURCES:.f90=.o)
CUDA_OBJS   = $(filter lb_cuda_%.o,$(F_OBJS))
HOST_OBJS   = $(filter-out lb_cuda_%.o,$(F_OBJS))

def:	all

all:
	@echo "Error - you should specify a target machine!"
	@echo "Possible choices are:              "
	@echo "                                   "
	@echo "nvfortran                          "
	@echo "nvfortran-debug                    "
	@echo "nvfortran-mpi                      "
	@echo "nvfortran-mpi-debug                "
	@echo "nvfortran-nvml                     "
	@echo "nvfortran-nvml-mpi                 "
	@echo "                                   "
	@echo "Possible choices for debugging are:"
	@echo "                                   "
	@echo "                                   "
	@echo "Please examine Makefile for further details "

help: all
	
nvfortran:
	$(MAKE) CC=nvcc \
	FC=nvfortran \
	CFLAGS="$(LBFAST_RELEASE_OPT) -I$(CURDIR) " \
	FCUDAFLAGS="$(FREL) $(LBFAST_CUDAFLAGS) $(GREL) -I$(CURDIR) " \
	FFLAGS="$(FREL) $(GREL) -I$(CURDIR) " \
	LDFLAGS="$(LDREL) -I$(CURDIR) -o" \
	TYPE=seq \
	EX=$(EXP) BINROOT=$(BINROOT) seq
	
nvfortran-debug:
	$(MAKE) CC=nvcc \
	FC=nvfortran \
	CFLAGS="$(LBFAST_DEBUG_OPT) -I$(CURDIR) " \
	FCUDAFLAGS="$(FDBG) $(LBFAST_CUDAFLAGS) $(GDBG) -I$(CURDIR) " \
	FFLAGS="$(FDBG) $(GDBG) -I$(CURDIR) " \
	LDFLAGS="$(LDDBG) -I$(CURDIR) -o" \
	TYPE=seq \
	EX=$(EX) BINROOT=$(BINROOT) seq

nvfortran-mpi:
	$(MAKE) CC=mpicc \
	FC=mpif90 \
	CFLAGS="$(LBFAST_RELEASE_OPT) -I$(CURDIR) " \
	FCUDAFLAGS="$(FREL) $(LBFAST_CUDAFLAGS) $(GREL) -I$(CURDIR) " \
	FFLAGS="$(FREL) $(GREL) -DMPI -I$(CURDIR) " \
	LDFLAGS="$(LDREL) -DMPI -I$(CURDIR) -o" \
	TYPE=seq \
	EX=$(EXP) BINROOT=$(BINROOT) seq

nvfortran-mpi-debug:
	$(MAKE) CC=mpicc \
	FC=mpif90 \
	CFLAGS="$(LBFAST_DEBUG_OPT) -I$(CURDIR) " \
	FCUDAFLAGS="$(FDBG) $(LBFAST_CUDAFLAGS) $(GDBG) -I$(CURDIR) " \
	FFLAGS="$(FDBG) $(GDBG) -DMPI -I$(CURDIR) " \
	LDFLAGS="$(LDDBG) -DMPI -I$(CURDIR) -o" \
	TYPE=seq \
	EX=$(EXP) BINROOT=$(BINROOT) seq

nvfortran-nvml:
	$(MAKE) CC=nvcc \
	FC=nvfortran \
	CFLAGS="$(LBFAST_RELEASE_OPT) -I$(CURDIR) " \
	CFLAGSNVML="$(LBFAST_RELEASE_OPT) -I$(NVML_INC) -D_NVML -I$(CURDIR) " \
	FCUDAFLAGS="$(FREL) $(LBFAST_CUDAFLAGS) $(GREL) -D_NVML -I$(CURDIR) " \
	FFLAGS="$(FREL) $(GREL) -D_NVML -I$(CURDIR) " \
	LDFLAGS="$(LDREL) -I$(CURDIR) -o" \
	TYPE=seqnvml \
	EX=$(EXP) BINROOT=$(BINROOT) seqnvml

nvfortran-nvml-mpi:
	$(MAKE) CC=mpicc \
	FC=mpif90 \
	CFLAGS="$(LBFAST_RELEASE_OPT) -I$(CURDIR) " \
	CFLAGSNVML="$(LBFAST_RELEASE_OPT) -I$(NVML_INC) -D_NVML -I$(CURDIR) " \
	FCUDAFLAGS="$(FREL) $(LBFAST_CUDAFLAGS) $(GREL) -D_NVML -I$(CURDIR) " \
	FFLAGS="$(FREL) $(GREL) -D_NVML -DMPI -I$(CURDIR) " \
	LDFLAGS="$(LDREL) -DMPI -I$(CURDIR) -o" \
	TYPE=seqnvml \
	EX=$(EXP) BINROOT=$(BINROOT) seqnvml

seq: $(C_OBJS) $(F_OBJS)
	$(FC) $(LDFLAGS) $(EX) $(C_OBJS) $(F_OBJS)

seqnvml: $(C_OBJS) $(NVML_C_OBJS) $(NVML_F_OBJS) $(F_OBJS)
	$(FC) $(LDFLAGS) $(EX) $(C_OBJS) $(NVML_C_OBJS) $(NVML_F_OBJS) $(F_OBJS) \
	-L$(NVML_LIB) -lnvidia-ml

$(C_OBJS): %.o: %.c
	$(CC) $(CFLAGS) -c $<

$(NVML_C_OBJS): %.o: %.c $(LBFAST_HEADERS)
	$(CC) $(CFLAGSNVML) -c $<

$(HOST_OBJS) $(NVML_F_OBJS): %.o: %.f90 $(LBFAST_HEADERS)
	$(FC) $(FFLAGS) -c $<

$(CUDA_OBJS): %.o: %.f90 $(LBFAST_HEADERS)
	$(FC) $(FCUDAFLAGS) -c $<

clean-all:
	rm -rf *.mod *.o *.x *.dat *.i

clean:
	rm -rf *.mod *.o *.i
