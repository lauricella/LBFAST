#
# makefile
#

BINROOT = $(CURDIR)
FC=undefined
CC=undefined
CFLAGSNVML=undefined
#NVML_INC = $(NVHPC_HOME)/Linux_x86_64/24.3/cuda/12.3/targets/x86_64-linux/include
#NVML_LIB = $(NVHPC_HOME)/Linux_x86_64/24.3/cuda/12.3/targets/x86_64-linux/lib/stubs
NVML_INC = $(CURDIR)
NVML_LIB = $(CURDIR)
FCUDAFLAGS=undefined
FCUDAFLAGSRID=undefined
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
	CFLAGS="-O3 -I$(CURDIR) " \
	FCUDAFLAGS="-O3 -cpp -acc -cuda -gpu=cc$(GPUCC),lineinfo,ptxinfo,maxregcount:64 -Minfo=accel -I$(CURDIR) " \
	FCUDAFLAGSRID="-O3 -cpp -acc -cuda -gpu=cc$(GPUCC),lineinfo,ptxinfo,maxregcount:64 -Minfo=accel -I$(CURDIR) " \
	FFLAGS="-O3 -cpp -acc -gpu=cc$(GPUCC),lineinfo,ptxinfo,maxregcount:64 -Minfo=accel -I$(CURDIR) " \
	LDFLAGS="-O3 -acc -cuda -gpu=cc$(GPUCC),lineinfo -Minfo=accel -I$(CURDIR) -o" \
	TYPE=seq \
	EX=$(EXP) BINROOT=$(BINROOT) seq
	
nvfortran-debug:
	$(MAKE) CC=nvcc \
	FC=nvfortran \
	CFLAGS="-O0 -g -I$(CURDIR) " \
	FCUDAFLAGS="-O0 -g -cpp -acc -cuda -gpu=cc$(GPUCC),debug,keep,ptxinfo,lineinfo,maxregcount:64 -Minfo=accel -Mchkptr -Mchkstk -traceback -I$(CURDIR) " \
	FCUDAFLAGSRID="-O0 -g -cpp -acc -cuda -gpu=cc$(GPUCC),lineinfo,ptxinfo,maxregcount:64 -Minfo=accel -Mchkptr -Mchkstk -traceback -I$(CURDIR) " \
	FFLAGS="-O0 -g -cpp -acc -gpu=cc$(GPUCC),debug,keep,ptxinfo,lineinfo,maxregcount:64 -Minfo=accel -Mchkptr -Mchkstk -traceback -I$(CURDIR) " \
	LDFLAGS="-O0 -acc -cuda -gpu=cc$(GPUCC),lineinfo -Minfo=accel -Mchkptr -Mchkstk -traceback -I$(CURDIR) -o" \
	TYPE=seq \
	EX=$(EX) BINROOT=$(BINROOT) seq

nvfortran-mpi:
	$(MAKE) CC=mpicc \
	FC=mpif90 \
	CFLAGS="-O3 -I$(CURDIR) " \
	FCUDAFLAGS="-O3 -cpp -acc -cuda -gpu=cc$(GPUCC),lineinfo,ptxinfo,maxregcount:64 -Minfo=accel -I$(CURDIR) " \
	FCUDAFLAGSRID="-O3 -cpp -acc -cuda -gpu=cc$(GPUCC),lineinfo,ptxinfo,maxregcount:64 -Minfo=accel -I$(CURDIR) " \
	FFLAGS="-O3 -cpp -acc -gpu=cc$(GPUCC),lineinfo,ptxinfo,maxregcount:64 -Minfo=accel -DMPI -I$(CURDIR) " \
	LDFLAGS="-O3 -acc -cuda -gpu=cc$(GPUCC),lineinfo -Minfo=accel -DMPI -I$(CURDIR) -o" \
	TYPE=seq \
	EX=$(EXP) BINROOT=$(BINROOT) seq

nvfortran-mpi-debug:
	$(MAKE) CC=mpicc \
	FC=mpif90 \
	CFLAGS="-O0 -g -I$(CURDIR) " \
	FCUDAFLAGS="-O0 -g -cpp -acc -cuda -gpu=cc$(GPUCC),debug,keep,ptxinfo,maxregcount:64 -Minfo=accel -Mchkptr -Mchkstk -traceback -I$(CURDIR) " \
	FCUDAFLAGSRID="-O0 -g -cpp -acc -cuda -gpu=cc$(GPUCC),lineinfo,ptxinfo,maxregcount:64 -Minfo=accel -Mchkptr -Mchkstk -traceback -I$(CURDIR) " \
	FFLAGS="-O0 -g -cpp -acc -gpu=cc$(GPUCC),debug,keep,ptxinfo,maxregcount:64 -Minfo=accel -Mchkptr -Mchkstk -traceback -DMPI -I$(CURDIR) " \
	LDFLAGS="-O0 -acc -cuda -gpu=cc$(GPUCC) -Minfo=accel -Mchkptr -Mchkstk -traceback -DMPI -I$(CURDIR) -o" \
	TYPE=seq \
	EX=$(EXP) BINROOT=$(BINROOT) seq

nvfortran-nvml:
	$(MAKE) CC=nvcc \
	FC=nvfortran \
	CFLAGS="-O3 -I$(CURDIR) " \
	CFLAGSNVML="-O3 -I$(NVML_INC) -D_NVML -I$(CURDIR) " \
	FCUDAFLAGS="-O3 -cpp -acc -cuda -gpu=cc$(GPUCC),lineinfo,ptxinfo,maxregcount:64 -Minfo=accel -D_NVML -I$(CURDIR) " \
	FCUDAFLAGSRID="-O3 -cpp -acc -cuda -gpu=cc$(GPUCC),lineinfo,ptxinfo,maxregcount:64 -Minfo=accel -D_NVML -I$(CURDIR) " \
	FFLAGS="-O3 -cpp -acc -gpu=cc$(GPUCC),lineinfo,ptxinfo,maxregcount:64 -Minfo=accel -D_NVML -I$(CURDIR) " \
	LDFLAGS="-O3 -acc -cuda -gpu=cc$(GPUCC),lineinfo -Minfo=accel -I$(CURDIR) -o" \
	TYPE=seqnvml \
	EX=$(EXP) BINROOT=$(BINROOT) seqnvml

nvfortran-nvml-mpi:
	$(MAKE) CC=mpicc \
	FC=mpif90 \
	CFLAGS="-O3 -I$(CURDIR) " \
	CFLAGSNVML="-O3 -I$(NVML_INC)  -I$(CURDIR) " \
	FCUDAFLAGS="-O3 -cpp -acc -cuda -gpu=cc$(GPUCC),lineinfo,ptxinfo,maxregcount:64 -Minfo=accel -D_NVML -I$(CURDIR) " \
	FCUDAFLAGSRID="-O3 -cpp -acc -cuda -gpu=cc$(GPUCC),lineinfo,ptxinfo,maxregcount:64 -Minfo=accel -D_NVML -I$(CURDIR) " \
	FFLAGS="-O3 -cpp -acc -gpu=cc$(GPUCC),lineinfo,ptxinfo,maxregcount:64 -Minfo=accel -D_NVML -DMPI -I$(CURDIR) " \
	LDFLAGS="-O3 -acc -cuda -gpu=cc$(GPUCC),lineinfo -Minfo=accel -DMPI -I$(CURDIR) -o" \
	TYPE=seqnvml \
	EX=$(EXP) BINROOT=$(BINROOT) seqnvml

seq:get_mem.o get_ram.o vars_module.o \
	mpi_module.o profiling_m.o lb_cuda_vars_module.o \
	lb_cuda_auxfields_module.o lb_cuda_repulsive_module.o lb_cuda_moments_module.o \
	lb_cuda_fused_module.o lb_cuda_update_phi_module.o lb_cuda_boundary_module.o \
	lb_cuda_driver_module.o boundary_cds_module.o init_conditions_module.o \
	statistics.o print_module.o allocate_module.o integrator_module.o \
	LBFAST.o
	$(FC) $(LDFLAGS) $(EX) get_mem.o get_ram.o \
	vars_module.o mpi_module.o \
	profiling_m.o lb_cuda_vars_module.o lb_cuda_auxfields_module.o \
	lb_cuda_repulsive_module.o lb_cuda_moments_module.o \
	lb_cuda_fused_module.o lb_cuda_update_phi_module.o lb_cuda_boundary_module.o \
	lb_cuda_driver_module.o boundary_cds_module.o init_conditions_module.o \
	statistics.o print_module.o allocate_module.o integrator_module.o \
	LBFAST.o

seqnvml:get_mem.o get_ram.o vars_module.o \
	nvml_wrapper.o nvml_interface_module.o \
	mpi_module.o profiling_m.o lb_cuda_vars_module.o \
	lb_cuda_auxfields_module.o lb_cuda_repulsive_module.o lb_cuda_moments_module.o \
	lb_cuda_fused_module.o lb_cuda_update_phi_module.o lb_cuda_boundary_module.o \
	lb_cuda_driver_module.o boundary_cds_module.o init_conditions_module.o \
	statistics.o print_module.o allocate_module.o integrator_module.o \
	LBFAST.o
	$(FC) $(LDFLAGS) $(EX) get_mem.o get_ram.o nvml_wrapper.o nvml_interface_module.o \
	vars_module.o mpi_module.o \
	profiling_m.o lb_cuda_vars_module.o lb_cuda_auxfields_module.o \
	lb_cuda_repulsive_module.o lb_cuda_moments_module.o \
	lb_cuda_fused_module.o lb_cuda_update_phi_module.o lb_cuda_boundary_module.o \
	lb_cuda_driver_module.o boundary_cds_module.o init_conditions_module.o \
	statistics.o print_module.o allocate_module.o integrator_module.o \
	LBFAST.o -L$(NVML_LIB) -lnvidia-ml
#	mv $(EXP) $(EXEP)

get_mem.o:get_mem.c
	$(CC) $(CFLAGS) -c get_mem.c

get_ram.o:get_ram.c
	$(CC) $(CFLAGS) -c get_ram.c

nvml_wrapper.o: nvml_wrapper.c
	$(CC) $(CFLAGSNVML) -c nvml_wrapper.c

nvml_interface_module.o: nvml_interface_module.f90
	$(FC) $(FFLAGS) -c nvml_interface_module.f90

vars_module.o:vars_module.f90
	$(FC) $(FFLAGS) -c vars_module.f90
	
mpi_module.o:mpi_module.f90
	$(FC) $(FFLAGS) -c mpi_module.f90

profiling_m.o: profiling_m.f90
	$(FC) $(FFLAGS) -c profiling_m.f90

lb_cuda_vars_module.o:lb_cuda_vars_module.f90
	$(FC) $(FCUDAFLAGS) -c lb_cuda_vars_module.f90

lb_cuda_auxfields_module.o:lb_cuda_auxfields_module.f90
	$(FC) $(FCUDAFLAGS) -c lb_cuda_auxfields_module.f90

lb_cuda_repulsive_module.o:lb_cuda_repulsive_module.f90
	$(FC) $(FCUDAFLAGS) -c lb_cuda_repulsive_module.f90

lb_cuda_moments_module.o:lb_cuda_moments_module.f90
	$(FC) $(FCUDAFLAGS) -c lb_cuda_moments_module.f90

lb_cuda_fused_module.o:lb_cuda_fused_module.f90
	$(FC) $(FCUDAFLAGS) -c lb_cuda_fused_module.f90

lb_cuda_update_phi_module.o:lb_cuda_update_phi_module.f90
	$(FC) $(FCUDAFLAGS) -c lb_cuda_update_phi_module.f90

lb_cuda_boundary_module.o:lb_cuda_boundary_module.f90
	$(FC) $(FCUDAFLAGS) -c lb_cuda_boundary_module.f90

lb_cuda_driver_module.o:lb_cuda_driver_module.f90
	$(FC) $(FCUDAFLAGS) -c lb_cuda_driver_module.f90

boundary_cds_module.o:boundary_cds_module.f90
	$(FC) $(FFLAGS) -c boundary_cds_module.f90

init_conditions_module.o:init_conditions_module.f90
	$(FC) $(FFLAGS) -c init_conditions_module.f90

statistics.o:statistics.f90
	$(FC) $(FFLAGS) -c statistics.f90

print_module.o:print_module.f90
	$(FC) $(FFLAGS) -c print_module.f90

allocate_module.o:allocate_module.f90
	$(FC) $(FFLAGS) -c allocate_module.f90

integrator_module.o:integrator_module.f90
	$(FC) $(FFLAGS) -c integrator_module.f90

LBFAST.o:LBFAST.f90
	$(FC) $(FFLAGS) -c LBFAST.f90

clean-all:
	rm -rf *.mod *.o *.x *.dat *.i

clean:
	rm -rf *.mod *.o *.i


