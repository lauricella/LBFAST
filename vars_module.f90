#include "defines.h"

#ifdef __NVCOMPILER
#warning "Compiling with NVIDIA!"
#endif

#ifdef TWOCOMPONENT
#warning "TWOCOMPONENT: the code is compiled for a two component system" 
#endif


#ifdef DENSRATIO
#warning "DENSRATIO: the 2c is compiled with the density ratio model" 
#endif

#ifdef DROPLET
#warning "DROPLET: the 2c is compiled with a single droplet initial configuration" 
#endif

#ifdef MONOD
#warning "MONOD_like growth of biofilms is activated" 
#endif

#ifdef DOXDMF
#warning "DOXDMF the xdmf file writing is activated" 
#endif

#ifdef MULTIHIT
#warning "HOMOGENEOUS ISOTROPIC TURBULENCE"
#endif

#ifdef CSF
#warning "CSF: geometric continuum-surface-force capillary model activated"
#endif

#ifdef CRAY
#warning "CRAY: CRAY activated for cray compilers"
#endif

#ifdef ACCNOKERNELS
#warning "ACCNOKERNELS: ACCNOKERNELS activated"
#endif

#ifdef PRINTHALF
#warning "PRINTHALF: printing output in half-precision FP16 raw files"
#endif

#ifdef EXCHANGEVEL
#warning "EXCHANGEVEL: mpi exchange velocity performed for boundary condition"
#endif

#ifdef AVOIDMPIIO
#warning "AVOIDMPIIO: Collective MPI-IO writing is disabled"
#endif

#ifdef INTERNAL_OBSTACLES
#warning "INTERNAL_OBSTACLES: actiaved"
#endif

#ifdef BOUNCE_BACK
#warning "BOUNCE_BACK: actiaved"
#endif

#ifdef BUOYANCY_FORCING
#warning "BUOYANCY_FORCING: actiaved"
#endif

#ifdef IMPOSED_PRESSURE_GRADIENT
#warning "IMPOSED pressure grad: actiaved"
#endif

#ifdef REPULSIVE_FLUX
#warning "Repulsive interfacial flux: actiaved"
#ifndef TWOCOMPONENT
#error "TWOCOMPONENT must be actiaved if REPULSIVE_FLUX is activated"
#endif
#endif

#ifdef INTERFACE_INCOMP
#warning "reinforce incompressibility at interface: actiaved"
#endif

#ifdef SMAGORINSKI
#warning "SMAGORINSKI LES: activated"
#endif

#ifdef ASYNCMPI
#warning "ASYNCMPI: activated"
#endif

#ifdef _NVML
#warning "_NVML: activated"
#endif

#ifdef GETPOWER
#warning "GETPOWER: activated"
#endif

#ifdef MIXEDPRC
#warning "MIXEDPRC: mixed precision activated"
#endif

#ifndef LATTICE
#error "LATTICE not defined. Use -DLATTICE=27 or #define LATTICE 27 (or 19 or 15)"
#endif

#if defined(HIGHORDER) && (!defined(LATTICE) || (LATTICE != 27))
#error "HIGHORDER: if the macro HIGHORDER is defined you MUST USE -DLATTICE=27 or #define LATTICE 27"
#endif

#ifdef DOBENCHMARK
#warning "DOBENCHMARK: activated"
#endif

#ifdef VELUNIFORMV
#warning "VELUNIFORMV: activated"
#endif

#ifdef POISEUILLE
#warning "POISEUILLE: activated"
#endif

#ifdef TWOPOISEUILLE
#warning "TWOPOISEUILLE: activated"
#endif

#ifdef TAYLORGREEN
#warning "TAYLORGREEN: activated"
#endif

#ifdef LAMBTEST
#warning "LAMBTEST: activated"
#endif

#ifdef LAPLACE
#warning "LAPLACE: activated"
#endif

#ifdef CAPILLARYWAVE
#warning "CAPILLARYWAVE: activated"
#endif

#ifdef PRINTPHI
#warning "PRINTPHI: activated"
#endif

#ifdef USEGNUPLOT
#warning "USEGNUPLOT: activated"
#endif


!     Module vars: global parameters and variables of LBFAST. It defines the precision kinds, the
!     GPU tile sizes, the lattice velocity sets and weights, the Hermite matrices of the high-order
!     D3Q27 model, the field arrays stored tile by tile with the numbers of their components, and
!     utility functions (random numbers on host and device, number formatting, cut-off functions).
module vars
#ifdef _OPENACC
   use openacc
#endif
   use iso_c_binding, only: c_long_long,c_int,c_float,c_double
   implicit none

   ! Precision kinds, selected in defines.h. db (PRC = 4 or 8) is the compute precision: the
   ! arithmetic and the shared-memory buffers of the kernels use it. strdb is the storage precision
   ! of the field arrays (hfields, phifields, forces, auxfields, locauxfields): with MIXEDPRC it is
   ! set by STRPRC (2, 4 or 8), otherwise strdb = db. The kernels convert the stored values to db
   ! when they read them and round the results to strdb when they write them.
#if PRC==4
#warning "PRC 4: single precision activated"
   integer, parameter :: db=4 !kind(1.0)
#elif PRC==8
#warning "PRC 8: double precision activated"
   integer, parameter :: db=8 !kind(1.0)
#else
   //#error "ERROR in specifying PRC"
#endif

#ifdef MIXEDPRC
#if STRPRC==2
#warning "STRPRC 2: half precision activated for storaging"
   integer, parameter :: strdb=2 !kind(1.0)
#elif STRPRC==4
#warning "STRPRC 4: single precision activated for storaging"
   integer, parameter :: strdb=4 !kind(1.0)
#elif STRPRC==8
#warning "STRPRC 8: double precision activated for storaging"
   integer, parameter :: strdb=8 !kind(1.0)
#else
   //#error "ERROR in specifying STRPRC"
#endif
#else
   integer, parameter :: strdb=db !kind(1.0)
#endif

   ! isf: integer kind of the isfluid array; printdb: real kind of the output arrays rhoprint,
   ! pressprint and velprint (FP16 with PRINTHALF, FP32 otherwise)
   integer, parameter :: isf=1 !kind(1.0)
   
#ifdef PRINTHALF
   integer, parameter :: printdb=2
#else
   integer, parameter :: printdb=4
#endif
   integer(kind=8) :: acc_device_radeon=5

   real(kind=db), parameter :: ZERO=real(0.d0,kind=db)
   real(kind=db), parameter :: HALF=real(0.5d0,kind=db)
   real(kind=db), parameter :: ONE=real(1.d0,kind=db)
   real(kind=db), parameter :: TWO=real(2.d0,kind=db)
   real(kind=db), parameter :: THREE=real(3.d0,kind=db)
   real(kind=db), parameter :: FOUR=real(4.d0,kind=db)
   real(kind=db), parameter :: FIVE=real(5.d0,kind=db)
   real(kind=db), parameter :: SIX=real(6.d0,kind=db)
   real(kind=db), parameter :: SEVEN=real(7.d0,kind=db)
   real(kind=db), parameter :: EIGHT=real(8.d0,kind=db)
   real(kind=db), parameter :: NINE=real(9.d0,kind=db)
   real(kind=db), parameter :: TEN=real(10.d0,kind=db)
   real(kind=db), parameter :: TWELVE=real(12.d0,kind=db)
   real(kind=db), parameter :: FOURTEEN=real(14.d0,kind=db)
   real(kind=db), parameter :: SIXTEEN=real(16.d0,kind=db)
   real(kind=db), parameter :: EIGHTEEN=real(18.d0,kind=db)
   real(kind=db), parameter :: TWENTYFOUR=real(24.d0,kind=db)
   real(kind=db), parameter :: TWENTYSEVEN=real(27.d0,kind=db)
   real(kind=db), parameter :: THERTYSIX=real(36.d0,kind=db)
   real(kind=db), parameter :: SEVENTYTWO=real(72.d0,kind=db)
   real(kind=db), parameter :: ONEHUNDREDEIGHT=real(108.d0,kind=db)
   real(kind=db), parameter :: TWOHUNDREDSIXTEEN=real(216.d0,kind=db)
   
   real(kind=db), parameter :: ZEROSTR=real(0.d0,kind=strdb)

   !aritra 
   ! !!!for WENO5
   ! (the WENO5 coefficients below are not referenced elsewhere in the source)
   real(kind=db), parameter ::  c11    =  2.0_db/6.0_db, &
                           c21    = -7.0_db/6.0_db, & 
                           c31    = 11.0_db/6.0_db, & 
                           c12    = -1.0_db/6.0_db, & 
                           c22    =  5.0_db/6.0_db, & 
                           c32    =  2.0_db/6.0_db, & 
                           c13    =  2.0_db/6.0_db, & 
                           c23    =  5.0_db/6.0_db, & 
                           c33    = -1.0_db/6.0_db 

   integer :: i,j,k,ll,l
   integer :: gi,gj,gk
   ! iprobe, jprobe, kprobe: global coordinates of the probe node; nx, ny, nz: size of the subdomain
   ! of this MPI rank; lx, ly, lz: size of the global lattice; nxskip, nyskip, nzskip: size of the
   ! output arrays, sampled every stepskip nodes
   integer :: iprobe=0,jprobe=0,kprobe=0
   integer :: nx,ny,nz,step,step_flip,stamp,stamp2D,nsteps,ngpus
   integer :: stamp_term=huge(1)
   integer :: nxskip,nyskip,nzskip,lxskip,lyskip,lzskip
   integer :: stepskip=1
   integer :: lx,ly,lz
   integer :: iframe,iframe2D
   integer :: physic_type
   integer :: openbc=0
   
   integer, dimension(2) :: openbc_type_x=0
   integer, dimension(2) :: openbc_type_y=0
   integer, dimension(2) :: openbc_type_z=0
   
   real(kind=db), dimension(2) :: openbc_press_x=ZERO
   real(kind=db), dimension(2) :: openbc_press_y=ZERO
   real(kind=db), dimension(2) :: openbc_press_z=ZERO
   real(kind=db), dimension(2) :: openbc_u_x=ZERO
   real(kind=db), dimension(2) :: openbc_u_y=ZERO
   real(kind=db), dimension(2) :: openbc_u_z=ZERO
   real(kind=db), dimension(2) :: openbc_v_x=ZERO
   real(kind=db), dimension(2) :: openbc_v_y=ZERO
   real(kind=db), dimension(2) :: openbc_v_z=ZERO
   real(kind=db), dimension(2) :: openbc_w_x=ZERO
   real(kind=db), dimension(2) :: openbc_w_y=ZERO
   real(kind=db), dimension(2) :: openbc_w_z=ZERO
   
   logical, save :: lprint,lvtk,lasync,lraw,lrestart
   logical, save :: lreadisfluid=.false.
   logical, save :: lreadinit=.false.
   logical, save :: lwriterestart=.false.
   logical, parameter :: lreadinput=.true.
   logical, save :: lweakscaling=.false.
   
   integer, save :: initseed=317
   
   integer :: narg,inumchar
   logical :: mydiagnostic
   integer :: tdiagnostic
   real(kind=db) :: smemory,sram

   integer :: nplanes
   integer, allocatable, dimension(:) :: ndir,npoint,skip_npoint

#ifdef _OPENACC
   integer :: devNum
   integer(acc_device_kind) :: devType
#endif

   ! Tile sizes. The subdomain of each MPI rank is divided into tiles of TILE_DIMx x TILE_DIMy x
   ! TILE_DIMz nodes (GPUTILEX, GPUTILEY, GPUTILEZ, set to 8 in defines.h; 4 if they are not
   ! defined). A thread block of the
   ! fused and moment kernels has one thread per node of a tile, and nx, ny, nz must be multiples
   ! of the tile sizes. TILE_DIM is the side of the two-dimensional thread block dimBlock2 set in
   ! setup_cuda.
#ifdef GPUTILEX
#warning "GPUTILEX: the TILE_DIMx value was defined in defines.h" 
   integer, parameter :: TILE_DIMx = GPUTILEX
#else
   integer, parameter :: TILE_DIMx = 4
#endif

#ifdef GPUTILEY
#warning "GPUTILEY: the TILE_DIMy value was defined in defines.h" 
   integer, parameter :: TILE_DIMy = GPUTILEY
#else
   integer, parameter :: TILE_DIMy = 4
#endif

#ifdef GPUTILEZ
#warning "GPUTILEZ: the TILE_DIMz value was defined in defines.h" 
   integer, parameter :: TILE_DIMz = GPUTILEZ
#else
   integer, parameter :: TILE_DIMz = 4
#endif
   integer, parameter :: TILE_DIM=16
   
   ! Tiles along x, y and z of the subdomain, including one layer of halo tiles on each side
   ! (nxblock = nx/TILE_DIMx + 2, set in setup_cuda); nxyblock = nxblock*nyblock, and nblocks is the
   ! total number of tiles. The fields are stored as field(TILE_DIMx,TILE_DIMy,TILE_DIMz,ncomp,nblocks):
   ! node i of the subdomain (1-TILE_DIMx <= i <= nx+TILE_DIMx) lies in the tile with x index
   ! xblock = (i+2*TILE_DIMx-1)/TILE_DIMx at the local index ii = i-(xblock-2)*TILE_DIMx, likewise
   ! for j and k, and the tile number is myblock = (xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1.
   integer, save :: nxblock,nyblock,nzblock,nxyblock,nblocks

   real(kind=db),parameter :: pi_greek=real(3.1415926535897932384626433832795028841971d0,kind=db)

   real(kind=db)  :: ts1,ts2
   real(kind=db), dimension(3) :: center=ZERO
   real(kind=db) :: uu,udotc,omega,omega_diff,radius,width
   real(kind=db) :: tau_diff,fx,fy,fz,temp,tau1,tau2
   real(kind=db) :: sharp_c,sigma,beta,kapp
   real(kind=db) :: visc1,visc2
   real(kind=db) :: wettab_r,wettab_b,pc_rate
   real(kind=db), parameter :: eps1=sqrt(2.0_db)*6.0_db
#ifdef IMPOSED_PRESSURE_GRADIENT
   real(kind=db) :: rhoIN,rhoOUT
#endif
   ! squared lattice speed of sound, cs^2 = 1/3 for the three lattices, and its inverse
   real(kind=db), parameter :: cssq=real(1.d0/3.d0,kind=db)
   real(kind=db), parameter :: invcssq=real(3.d0,kind=db)

   ! Lattice velocity sets, selected by LATTICE (15, 19 or 27) in defines.h. nlinks is the number of
   ! moving directions; ex, ey, ez are the integer components of the velocities c_l, l = 0..nlinks,
   ! with l = 0 the rest direction. The moving directions are ordered in opposite pairs,
   ! c_(l+1) = -c_l for odd l, and opp(l) is the index of the direction opposite to l.
   ! gen_fused2_stages.py relies on this order and reads ex, ey, ez and the weight names of p from
   ! this file: run it again after any change of the lattice definitions.
   ! p(l) are the lattice weights: sum_l p(l) = 1 and sum_l p(l) c_la c_lb = cs^2 delta_ab.
#if LATTICE == 27
#warning "LATTICE 27: the lattice D3Q27 is utilized"
   integer, parameter :: nlinks=26
   !lattice vectors
   integer, dimension(0:nlinks), parameter :: &
   !          0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26
      ex=   (/0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 0, 0, 0, 0, 1,-1,-1, 1, 1,-1, 1,-1,-1, 1, 1,-1/)
   integer, dimension(0:nlinks), parameter:: &
      ey=  (/ 0, 0, 0, 1,-1, 0, 0, 1,-1,-1, 1, 1,-1, 1,-1, 0, 0, 0, 0, 1,-1,-1, 1,-1, 1,-1, 1/)
   integer, dimension(0:nlinks), parameter:: &
      ez=  (/ 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 1,-1,-1, 1, 1,-1, 1,-1, 1,-1, 1,-1, 1,-1,-1, 1/)
   integer, dimension(0:nlinks), parameter:: &
      opp =(/ 0, 2, 1, 4, 3, 6, 5, 8, 7,10, 9,12,11,14,13,16,15,18,17,20,19,22,21,24,23,26,25/)

   ! D3Q27 weights: rest direction 0 (p0), axis directions 1-6 (p1), face diagonals 7-18 (p2),
   ! body diagonals 19-26 (p3)
   real(kind=db), parameter :: p0=real(8.d0/27.d0,kind=db)
   real(kind=db), parameter :: p1=real(2.d0/27.d0,kind=db)
   real(kind=db), parameter :: p2=real(1.d0/54.d0,kind=db)
   real(kind=db), parameter :: p3=real(1.d0/216.d0,kind=db)
   real(kind=db), dimension(0:nlinks), parameter :: &
      p=(/p0,p1,p1,p1,p1,p1,p1,p2,p2,p2,p2,p2,p2,p2,p2,p2,p2,p2,p2,p3,p3,p3,p3,p3,p3,p3,p3/)
      
#ifdef HIGHORDER
#warning "HIGHORDER: the lattice D3Q27 is utilized with high-order"
   ! Tensor-product Hermite basis of D3Q27 used by the high-order (recursive regularized) model.
   ! Basis function k = 0..26 is a product of one factor per Cartesian axis: 1, H1(c_a) = c_a or
   ! H2(c_a) = c_a^2 - cs^2. In the labels below, a single letter is a factor H1 along that axis
   ! and a repeated letter a factor H2:
   !   0: 1;  1-3: x, y, z;  4-6: xx, yy, zz;  7-9: xy, xz, yz;  10-15: xxy, xxz, xyy, xzz, yzz, yyz;
   !   16: xyz;  17-19: xxyy, xxzz, yyzz;  20-22: xyzz, xyyz, xxyz;  23-25: xxyzz, xxyyz, xyyzz;
   !   26: xxyyzz.
   ! Hermite_norm(k) = N_k = sum_l p(l) H_k(c_l)^2 is the norm of basis function k, and
   ! Hermite_prefact(k) = 1/N_k.
   real(kind=db), dimension(0:nlinks), parameter :: &   !eq. 37 Malaspinas 2015 !equivalent to b_k^-1 in 2.60 dissertation shiller
      Hermite_prefact=(/ONE,ONE/cssq,ONE/cssq,ONE/cssq, &
       ONE/(TWO*cssq**TWO),ONE/(TWO*cssq**TWO),ONE/(TWO*cssq**TWO), &
       TWO/(TWO*cssq**TWO),TWO/(TWO*cssq**TWO),TWO/(TWO*cssq**TWO), &
       ONE/(TWO*cssq**THREE),ONE/(TWO*cssq**THREE),ONE/(TWO*cssq**THREE), &
       ONE/(TWO*cssq**THREE),ONE/(TWO*cssq**THREE),ONE/(TWO*cssq**THREE), &
       TWO/(TWO*cssq**THREE), &
       ONE/(FOUR*cssq**FOUR),ONE/(FOUR*cssq**FOUR),ONE/(FOUR*cssq**FOUR), &
       TWO/(FOUR*cssq**FOUR),TWO/(FOUR*cssq**FOUR),TWO/(FOUR*cssq**FOUR), &
       ONE/(FOUR*cssq**FIVE),ONE/(FOUR*cssq**FIVE),ONE/(FOUR*cssq**FIVE), &
       ONE/(EIGHT*cssq**SIX)/)
       
   real(kind=db), dimension(0:nlinks), parameter :: &  
      Hermite_norm=ONE/Hermite_prefact   !equivalent to b_k in 2.60 dissertation shiller
      
   ! sqrtHermite_norm is not referenced elsewhere in the source
   real(kind=db), dimension(0:nlinks), parameter :: &  
      sqrtHermite_norm=sqrt(Hermite_norm)  !equivalent to sqrt(b_k) in 2.60 dissertation shiller   
   
   ! Minv_mrt(k,l) = p(l) H_k(c_l)/N_k (first index: basis function k, second index: direction l).
   ! A population is reconstructed from its Hermite coefficients a_k = sum_l H_k(c_l) f_l as
   ! f_l = sum_k Minv_mrt(k,l) a_k. The constructor lists the values direction by direction, and the
   ! column-major reshape stores each group of 27 values along the first index k.
   ! Minv_mrt is used in boundary_cds_module.f90; the high-order fused kernels use the same values,
   ! written as exact fractions in fusedh_stages.h by gen_fused2_stages.py.
   real(kind=db), parameter :: Minv_mrt(0:nlinks,0:nlinks) = reshape([ & !first index i-th pop, second index k-th basis !it is trasposed to be column major compliant!
    [8.0_db/27.0_db, 0.0_db, 0.0_db, 0.0_db, -4.0_db/9.0_db, &
	-4.0_db/9.0_db, -4.0_db/9.0_db, 0.0_db, 0.0_db, 0.0_db, 0.0_db, &
	0.0_db, 0.0_db, 0.0_db, 0.0_db, 0.0_db, 0.0_db, 2.0_db/3.0_db, &
	2.0_db/3.0_db, 2.0_db/3.0_db, 0.0_db, 0.0_db, 0.0_db, 0.0_db, 0.0_db, &
	0.0_db, -1.0_db], [2.0_db/27.0_db, 2.0_db/9.0_db, 0.0_db, 0.0_db, &
	2.0_db/9.0_db, -1.0_db/9.0_db, -1.0_db/9.0_db, 0.0_db, 0.0_db, &
	0.0_db, 0.0_db, 0.0_db, -1.0_db/3.0_db, -1.0_db/3.0_db, 0.0_db, &
	0.0_db, 0.0_db, -1.0_db/3.0_db, -1.0_db/3.0_db, 1.0_db/6.0_db, &
	0.0_db, 0.0_db, 0.0_db, 0.0_db, 0.0_db, 1.0_db/2.0_db, &
	1.0_db/2.0_db], [2.0_db/27.0_db, -2.0_db/9.0_db, 0.0_db, 0.0_db, &
	2.0_db/9.0_db, -1.0_db/9.0_db, -1.0_db/9.0_db, 0.0_db, 0.0_db, &
	0.0_db, 0.0_db, 0.0_db, 1.0_db/3.0_db, 1.0_db/3.0_db, 0.0_db, 0.0_db, &
	0.0_db, -1.0_db/3.0_db, -1.0_db/3.0_db, 1.0_db/6.0_db, 0.0_db, &
	0.0_db, 0.0_db, 0.0_db, 0.0_db, -1.0_db/2.0_db, 1.0_db/2.0_db], &
	[2.0_db/27.0_db, 0.0_db, 2.0_db/9.0_db, 0.0_db, -1.0_db/9.0_db, &
	2.0_db/9.0_db, -1.0_db/9.0_db, 0.0_db, 0.0_db, 0.0_db, &
	-1.0_db/3.0_db, 0.0_db, 0.0_db, 0.0_db, -1.0_db/3.0_db, 0.0_db, &
	0.0_db, -1.0_db/3.0_db, 1.0_db/6.0_db, -1.0_db/3.0_db, 0.0_db, &
	0.0_db, 0.0_db, 1.0_db/2.0_db, 0.0_db, 0.0_db, 1.0_db/2.0_db], &
	[2.0_db/27.0_db, 0.0_db, -2.0_db/9.0_db, 0.0_db, -1.0_db/9.0_db, &
	2.0_db/9.0_db, -1.0_db/9.0_db, 0.0_db, 0.0_db, 0.0_db, 1.0_db/3.0_db, &
	0.0_db, 0.0_db, 0.0_db, 1.0_db/3.0_db, 0.0_db, 0.0_db, &
	-1.0_db/3.0_db, 1.0_db/6.0_db, -1.0_db/3.0_db, 0.0_db, 0.0_db, &
	0.0_db, -1.0_db/2.0_db, 0.0_db, 0.0_db, 1.0_db/2.0_db], &
	[2.0_db/27.0_db, 0.0_db, 0.0_db, 2.0_db/9.0_db, -1.0_db/9.0_db, &
	-1.0_db/9.0_db, 2.0_db/9.0_db, 0.0_db, 0.0_db, 0.0_db, 0.0_db, &
	-1.0_db/3.0_db, 0.0_db, 0.0_db, 0.0_db, -1.0_db/3.0_db, 0.0_db, &
	1.0_db/6.0_db, -1.0_db/3.0_db, -1.0_db/3.0_db, 0.0_db, 0.0_db, &
	0.0_db, 0.0_db, 1.0_db/2.0_db, 0.0_db, 1.0_db/2.0_db], &
	[2.0_db/27.0_db, 0.0_db, 0.0_db, -2.0_db/9.0_db, -1.0_db/9.0_db, &
	-1.0_db/9.0_db, 2.0_db/9.0_db, 0.0_db, 0.0_db, 0.0_db, 0.0_db, &
	1.0_db/3.0_db, 0.0_db, 0.0_db, 0.0_db, 1.0_db/3.0_db, 0.0_db, &
	1.0_db/6.0_db, -1.0_db/3.0_db, -1.0_db/3.0_db, 0.0_db, 0.0_db, &
	0.0_db, 0.0_db, -1.0_db/2.0_db, 0.0_db, 1.0_db/2.0_db], &
	[1.0_db/54.0_db, 1.0_db/18.0_db, 1.0_db/18.0_db, 0.0_db, &
	1.0_db/18.0_db, 1.0_db/18.0_db, -1.0_db/36.0_db, 1.0_db/6.0_db, &
	0.0_db, 0.0_db, 1.0_db/6.0_db, 0.0_db, 1.0_db/6.0_db, &
	-1.0_db/12.0_db, -1.0_db/12.0_db, 0.0_db, 0.0_db, 1.0_db/6.0_db, &
	-1.0_db/12.0_db, -1.0_db/12.0_db, -1.0_db/4.0_db, 0.0_db, 0.0_db, &
	-1.0_db/4.0_db, 0.0_db, -1.0_db/4.0_db, -1.0_db/4.0_db], &
	[1.0_db/54.0_db, -1.0_db/18.0_db, -1.0_db/18.0_db, 0.0_db, &
	1.0_db/18.0_db, 1.0_db/18.0_db, -1.0_db/36.0_db, 1.0_db/6.0_db, &
	0.0_db, 0.0_db, -1.0_db/6.0_db, 0.0_db, -1.0_db/6.0_db, &
	1.0_db/12.0_db, 1.0_db/12.0_db, 0.0_db, 0.0_db, 1.0_db/6.0_db, &
	-1.0_db/12.0_db, -1.0_db/12.0_db, -1.0_db/4.0_db, 0.0_db, 0.0_db, &
	1.0_db/4.0_db, 0.0_db, 1.0_db/4.0_db, -1.0_db/4.0_db], &
	[1.0_db/54.0_db, 1.0_db/18.0_db, -1.0_db/18.0_db, 0.0_db, &
	1.0_db/18.0_db, 1.0_db/18.0_db, -1.0_db/36.0_db, -1.0_db/6.0_db, &
	0.0_db, 0.0_db, -1.0_db/6.0_db, 0.0_db, 1.0_db/6.0_db, &
	-1.0_db/12.0_db, 1.0_db/12.0_db, 0.0_db, 0.0_db, 1.0_db/6.0_db, &
	-1.0_db/12.0_db, -1.0_db/12.0_db, 1.0_db/4.0_db, 0.0_db, 0.0_db, &
	1.0_db/4.0_db, 0.0_db, -1.0_db/4.0_db, -1.0_db/4.0_db], &
	[1.0_db/54.0_db, -1.0_db/18.0_db, 1.0_db/18.0_db, 0.0_db, &
	1.0_db/18.0_db, 1.0_db/18.0_db, -1.0_db/36.0_db, -1.0_db/6.0_db, &
	0.0_db, 0.0_db, 1.0_db/6.0_db, 0.0_db, -1.0_db/6.0_db, &
	1.0_db/12.0_db, -1.0_db/12.0_db, 0.0_db, 0.0_db, 1.0_db/6.0_db, &
	-1.0_db/12.0_db, -1.0_db/12.0_db, 1.0_db/4.0_db, 0.0_db, 0.0_db, &
	-1.0_db/4.0_db, 0.0_db, 1.0_db/4.0_db, -1.0_db/4.0_db], &
	[1.0_db/54.0_db, 0.0_db, 1.0_db/18.0_db, 1.0_db/18.0_db, &
	-1.0_db/36.0_db, 1.0_db/18.0_db, 1.0_db/18.0_db, 0.0_db, 0.0_db, &
	1.0_db/6.0_db, -1.0_db/12.0_db, -1.0_db/12.0_db, 0.0_db, 0.0_db, &
	1.0_db/6.0_db, 1.0_db/6.0_db, 0.0_db, -1.0_db/12.0_db, &
	-1.0_db/12.0_db, 1.0_db/6.0_db, 0.0_db, 0.0_db, -1.0_db/4.0_db, &
	-1.0_db/4.0_db, -1.0_db/4.0_db, 0.0_db, -1.0_db/4.0_db], &
	[1.0_db/54.0_db, 0.0_db, -1.0_db/18.0_db, -1.0_db/18.0_db, &
	-1.0_db/36.0_db, 1.0_db/18.0_db, 1.0_db/18.0_db, 0.0_db, 0.0_db, &
	1.0_db/6.0_db, 1.0_db/12.0_db, 1.0_db/12.0_db, 0.0_db, 0.0_db, &
	-1.0_db/6.0_db, -1.0_db/6.0_db, 0.0_db, -1.0_db/12.0_db, &
	-1.0_db/12.0_db, 1.0_db/6.0_db, 0.0_db, 0.0_db, -1.0_db/4.0_db, &
	1.0_db/4.0_db, 1.0_db/4.0_db, 0.0_db, -1.0_db/4.0_db], &
	[1.0_db/54.0_db, 0.0_db, 1.0_db/18.0_db, -1.0_db/18.0_db, &
	-1.0_db/36.0_db, 1.0_db/18.0_db, 1.0_db/18.0_db, 0.0_db, 0.0_db, &
	-1.0_db/6.0_db, -1.0_db/12.0_db, 1.0_db/12.0_db, 0.0_db, 0.0_db, &
	1.0_db/6.0_db, -1.0_db/6.0_db, 0.0_db, -1.0_db/12.0_db, &
	-1.0_db/12.0_db, 1.0_db/6.0_db, 0.0_db, 0.0_db, 1.0_db/4.0_db, &
	-1.0_db/4.0_db, 1.0_db/4.0_db, 0.0_db, -1.0_db/4.0_db], &
	[1.0_db/54.0_db, 0.0_db, -1.0_db/18.0_db, 1.0_db/18.0_db, &
	-1.0_db/36.0_db, 1.0_db/18.0_db, 1.0_db/18.0_db, 0.0_db, 0.0_db, &
	-1.0_db/6.0_db, 1.0_db/12.0_db, -1.0_db/12.0_db, 0.0_db, 0.0_db, &
	-1.0_db/6.0_db, 1.0_db/6.0_db, 0.0_db, -1.0_db/12.0_db, &
	-1.0_db/12.0_db, 1.0_db/6.0_db, 0.0_db, 0.0_db, 1.0_db/4.0_db, &
	1.0_db/4.0_db, -1.0_db/4.0_db, 0.0_db, -1.0_db/4.0_db], &
	[1.0_db/54.0_db, 1.0_db/18.0_db, 0.0_db, 1.0_db/18.0_db, &
	1.0_db/18.0_db, -1.0_db/36.0_db, 1.0_db/18.0_db, 0.0_db, &
	1.0_db/6.0_db, 0.0_db, 0.0_db, 1.0_db/6.0_db, -1.0_db/12.0_db, &
	1.0_db/6.0_db, 0.0_db, -1.0_db/12.0_db, 0.0_db, -1.0_db/12.0_db, &
	1.0_db/6.0_db, -1.0_db/12.0_db, 0.0_db, -1.0_db/4.0_db, 0.0_db, &
	0.0_db, -1.0_db/4.0_db, -1.0_db/4.0_db, -1.0_db/4.0_db], &
	[1.0_db/54.0_db, -1.0_db/18.0_db, 0.0_db, -1.0_db/18.0_db, &
	1.0_db/18.0_db, -1.0_db/36.0_db, 1.0_db/18.0_db, 0.0_db, &
	1.0_db/6.0_db, 0.0_db, 0.0_db, -1.0_db/6.0_db, 1.0_db/12.0_db, &
	-1.0_db/6.0_db, 0.0_db, 1.0_db/12.0_db, 0.0_db, -1.0_db/12.0_db, &
	1.0_db/6.0_db, -1.0_db/12.0_db, 0.0_db, -1.0_db/4.0_db, 0.0_db, &
	0.0_db, 1.0_db/4.0_db, 1.0_db/4.0_db, -1.0_db/4.0_db], &
	[1.0_db/54.0_db, -1.0_db/18.0_db, 0.0_db, 1.0_db/18.0_db, &
	1.0_db/18.0_db, -1.0_db/36.0_db, 1.0_db/18.0_db, 0.0_db, &
	-1.0_db/6.0_db, 0.0_db, 0.0_db, 1.0_db/6.0_db, 1.0_db/12.0_db, &
	-1.0_db/6.0_db, 0.0_db, -1.0_db/12.0_db, 0.0_db, -1.0_db/12.0_db, &
	1.0_db/6.0_db, -1.0_db/12.0_db, 0.0_db, 1.0_db/4.0_db, 0.0_db, &
	0.0_db, -1.0_db/4.0_db, 1.0_db/4.0_db, -1.0_db/4.0_db], &
	[1.0_db/54.0_db, 1.0_db/18.0_db, 0.0_db, -1.0_db/18.0_db, &
	1.0_db/18.0_db, -1.0_db/36.0_db, 1.0_db/18.0_db, 0.0_db, &
	-1.0_db/6.0_db, 0.0_db, 0.0_db, -1.0_db/6.0_db, -1.0_db/12.0_db, &
	1.0_db/6.0_db, 0.0_db, 1.0_db/12.0_db, 0.0_db, -1.0_db/12.0_db, &
	1.0_db/6.0_db, -1.0_db/12.0_db, 0.0_db, 1.0_db/4.0_db, 0.0_db, &
	0.0_db, 1.0_db/4.0_db, -1.0_db/4.0_db, -1.0_db/4.0_db], &
	[1.0_db/216.0_db, 1.0_db/72.0_db, 1.0_db/72.0_db, 1.0_db/72.0_db, &
	1.0_db/72.0_db, 1.0_db/72.0_db, 1.0_db/72.0_db, 1.0_db/24.0_db, &
	1.0_db/24.0_db, 1.0_db/24.0_db, 1.0_db/24.0_db, 1.0_db/24.0_db, &
	1.0_db/24.0_db, 1.0_db/24.0_db, 1.0_db/24.0_db, 1.0_db/24.0_db, &
	1.0_db/8.0_db, 1.0_db/24.0_db, 1.0_db/24.0_db, 1.0_db/24.0_db, &
	1.0_db/8.0_db, 1.0_db/8.0_db, 1.0_db/8.0_db, 1.0_db/8.0_db, &
	1.0_db/8.0_db, 1.0_db/8.0_db, 1.0_db/8.0_db], [1.0_db/216.0_db, &
	-1.0_db/72.0_db, -1.0_db/72.0_db, -1.0_db/72.0_db, 1.0_db/72.0_db, &
	1.0_db/72.0_db, 1.0_db/72.0_db, 1.0_db/24.0_db, 1.0_db/24.0_db, &
	1.0_db/24.0_db, -1.0_db/24.0_db, -1.0_db/24.0_db, -1.0_db/24.0_db, &
	-1.0_db/24.0_db, -1.0_db/24.0_db, -1.0_db/24.0_db, -1.0_db/8.0_db, &
	1.0_db/24.0_db, 1.0_db/24.0_db, 1.0_db/24.0_db, 1.0_db/8.0_db, &
	1.0_db/8.0_db, 1.0_db/8.0_db, -1.0_db/8.0_db, -1.0_db/8.0_db, &
	-1.0_db/8.0_db, 1.0_db/8.0_db], [1.0_db/216.0_db, 1.0_db/72.0_db, &
	-1.0_db/72.0_db, 1.0_db/72.0_db, 1.0_db/72.0_db, 1.0_db/72.0_db, &
	1.0_db/72.0_db, -1.0_db/24.0_db, 1.0_db/24.0_db, -1.0_db/24.0_db, &
	-1.0_db/24.0_db, 1.0_db/24.0_db, 1.0_db/24.0_db, 1.0_db/24.0_db, &
	-1.0_db/24.0_db, 1.0_db/24.0_db, -1.0_db/8.0_db, 1.0_db/24.0_db, &
	1.0_db/24.0_db, 1.0_db/24.0_db, -1.0_db/8.0_db, 1.0_db/8.0_db, &
	-1.0_db/8.0_db, -1.0_db/8.0_db, 1.0_db/8.0_db, 1.0_db/8.0_db, &
	1.0_db/8.0_db], [1.0_db/216.0_db, -1.0_db/72.0_db, 1.0_db/72.0_db, &
	-1.0_db/72.0_db, 1.0_db/72.0_db, 1.0_db/72.0_db, 1.0_db/72.0_db, &
	-1.0_db/24.0_db, 1.0_db/24.0_db, -1.0_db/24.0_db, 1.0_db/24.0_db, &
	-1.0_db/24.0_db, -1.0_db/24.0_db, -1.0_db/24.0_db, 1.0_db/24.0_db, &
	-1.0_db/24.0_db, 1.0_db/8.0_db, 1.0_db/24.0_db, 1.0_db/24.0_db, &
	1.0_db/24.0_db, -1.0_db/8.0_db, 1.0_db/8.0_db, -1.0_db/8.0_db, &
	1.0_db/8.0_db, -1.0_db/8.0_db, -1.0_db/8.0_db, 1.0_db/8.0_db], &
	[1.0_db/216.0_db, -1.0_db/72.0_db, -1.0_db/72.0_db, 1.0_db/72.0_db, &
	1.0_db/72.0_db, 1.0_db/72.0_db, 1.0_db/72.0_db, 1.0_db/24.0_db, &
	-1.0_db/24.0_db, -1.0_db/24.0_db, -1.0_db/24.0_db, 1.0_db/24.0_db, &
	-1.0_db/24.0_db, -1.0_db/24.0_db, -1.0_db/24.0_db, 1.0_db/24.0_db, &
	1.0_db/8.0_db, 1.0_db/24.0_db, 1.0_db/24.0_db, 1.0_db/24.0_db, &
	1.0_db/8.0_db, -1.0_db/8.0_db, -1.0_db/8.0_db, -1.0_db/8.0_db, &
	1.0_db/8.0_db, -1.0_db/8.0_db, 1.0_db/8.0_db], [1.0_db/216.0_db, &
	1.0_db/72.0_db, 1.0_db/72.0_db, -1.0_db/72.0_db, 1.0_db/72.0_db, &
	1.0_db/72.0_db, 1.0_db/72.0_db, 1.0_db/24.0_db, -1.0_db/24.0_db, &
	-1.0_db/24.0_db, 1.0_db/24.0_db, -1.0_db/24.0_db, 1.0_db/24.0_db, &
	1.0_db/24.0_db, 1.0_db/24.0_db, -1.0_db/24.0_db, -1.0_db/8.0_db, &
	1.0_db/24.0_db, 1.0_db/24.0_db, 1.0_db/24.0_db, 1.0_db/8.0_db, &
	-1.0_db/8.0_db, -1.0_db/8.0_db, 1.0_db/8.0_db, -1.0_db/8.0_db, &
	1.0_db/8.0_db, 1.0_db/8.0_db], [1.0_db/216.0_db, 1.0_db/72.0_db, &
	-1.0_db/72.0_db, -1.0_db/72.0_db, 1.0_db/72.0_db, 1.0_db/72.0_db, &
	1.0_db/72.0_db, -1.0_db/24.0_db, -1.0_db/24.0_db, 1.0_db/24.0_db, &
	-1.0_db/24.0_db, -1.0_db/24.0_db, 1.0_db/24.0_db, 1.0_db/24.0_db, &
	-1.0_db/24.0_db, -1.0_db/24.0_db, 1.0_db/8.0_db, 1.0_db/24.0_db, &
	1.0_db/24.0_db, 1.0_db/24.0_db, -1.0_db/8.0_db, -1.0_db/8.0_db, &
	1.0_db/8.0_db, -1.0_db/8.0_db, -1.0_db/8.0_db, 1.0_db/8.0_db, &
	1.0_db/8.0_db], [1.0_db/216.0_db, -1.0_db/72.0_db, 1.0_db/72.0_db, &
	1.0_db/72.0_db, 1.0_db/72.0_db, 1.0_db/72.0_db, 1.0_db/72.0_db, &
	-1.0_db/24.0_db, -1.0_db/24.0_db, 1.0_db/24.0_db, 1.0_db/24.0_db, &
	1.0_db/24.0_db, -1.0_db/24.0_db, -1.0_db/24.0_db, 1.0_db/24.0_db, &
	1.0_db/24.0_db, -1.0_db/8.0_db, 1.0_db/24.0_db, 1.0_db/24.0_db, &
	1.0_db/24.0_db, -1.0_db/8.0_db, -1.0_db/8.0_db, 1.0_db/8.0_db, &
	1.0_db/8.0_db, 1.0_db/8.0_db, -1.0_db/8.0_db, 1.0_db/8.0_db]  &
    ], shape=[nlinks+1,nlinks+1])
   
   ! M_mrt(l,k) = H_k(c_l) (first index: direction l, second index: basis function k) is the inverse
   ! of Minv_mrt: a_k = sum_l M_mrt(l,k) f_l and sum_l Minv_mrt(k,l) M_mrt(l,m) = delta_km.
   ! M_mrt is not referenced elsewhere in the source.
   real(kind=db), parameter :: M_mrt(0:nlinks,0:nlinks) = reshape([ &
    [1.0_db, 1.0_db, 1.0_db, 1.0_db, 1.0_db, 1.0_db, 1.0_db, 1.0_db, &
	1.0_db, 1.0_db, 1.0_db, 1.0_db, 1.0_db, 1.0_db, 1.0_db, 1.0_db, &
	1.0_db, 1.0_db, 1.0_db, 1.0_db, 1.0_db, 1.0_db, 1.0_db, 1.0_db, &
	1.0_db, 1.0_db, 1.0_db], [0.0_db, 1.0_db, -1.0_db, 0.0_db, 0.0_db, &
	0.0_db, 0.0_db, 1.0_db, -1.0_db, 1.0_db, -1.0_db, 0.0_db, 0.0_db, &
	0.0_db, 0.0_db, 1.0_db, -1.0_db, -1.0_db, 1.0_db, 1.0_db, -1.0_db, &
	1.0_db, -1.0_db, -1.0_db, 1.0_db, 1.0_db, -1.0_db], [0.0_db, 0.0_db, &
	0.0_db, 1.0_db, -1.0_db, 0.0_db, 0.0_db, 1.0_db, -1.0_db, -1.0_db, &
	1.0_db, 1.0_db, -1.0_db, 1.0_db, -1.0_db, 0.0_db, 0.0_db, 0.0_db, &
	0.0_db, 1.0_db, -1.0_db, -1.0_db, 1.0_db, -1.0_db, 1.0_db, -1.0_db, &
	1.0_db], [0.0_db, 0.0_db, 0.0_db, 0.0_db, 0.0_db, 1.0_db, -1.0_db, &
	0.0_db, 0.0_db, 0.0_db, 0.0_db, 1.0_db, -1.0_db, -1.0_db, 1.0_db, &
	1.0_db, -1.0_db, 1.0_db, -1.0_db, 1.0_db, -1.0_db, 1.0_db, -1.0_db, &
	1.0_db, -1.0_db, -1.0_db, 1.0_db], [-1.0_db/3.0_db, 2.0_db/3.0_db, &
	2.0_db/3.0_db, -1.0_db/3.0_db, -1.0_db/3.0_db, -1.0_db/3.0_db, &
	-1.0_db/3.0_db, 2.0_db/3.0_db, 2.0_db/3.0_db, 2.0_db/3.0_db, &
	2.0_db/3.0_db, -1.0_db/3.0_db, -1.0_db/3.0_db, -1.0_db/3.0_db, &
	-1.0_db/3.0_db, 2.0_db/3.0_db, 2.0_db/3.0_db, 2.0_db/3.0_db, &
	2.0_db/3.0_db, 2.0_db/3.0_db, 2.0_db/3.0_db, 2.0_db/3.0_db, &
	2.0_db/3.0_db, 2.0_db/3.0_db, 2.0_db/3.0_db, 2.0_db/3.0_db, &
	2.0_db/3.0_db], [-1.0_db/3.0_db, -1.0_db/3.0_db, -1.0_db/3.0_db, &
	2.0_db/3.0_db, 2.0_db/3.0_db, -1.0_db/3.0_db, -1.0_db/3.0_db, &
	2.0_db/3.0_db, 2.0_db/3.0_db, 2.0_db/3.0_db, 2.0_db/3.0_db, &
	2.0_db/3.0_db, 2.0_db/3.0_db, 2.0_db/3.0_db, 2.0_db/3.0_db, &
	-1.0_db/3.0_db, -1.0_db/3.0_db, -1.0_db/3.0_db, -1.0_db/3.0_db, &
	2.0_db/3.0_db, 2.0_db/3.0_db, 2.0_db/3.0_db, 2.0_db/3.0_db, &
	2.0_db/3.0_db, 2.0_db/3.0_db, 2.0_db/3.0_db, 2.0_db/3.0_db], &
	[-1.0_db/3.0_db, -1.0_db/3.0_db, -1.0_db/3.0_db, -1.0_db/3.0_db, &
	-1.0_db/3.0_db, 2.0_db/3.0_db, 2.0_db/3.0_db, -1.0_db/3.0_db, &
	-1.0_db/3.0_db, -1.0_db/3.0_db, -1.0_db/3.0_db, 2.0_db/3.0_db, &
	2.0_db/3.0_db, 2.0_db/3.0_db, 2.0_db/3.0_db, 2.0_db/3.0_db, &
	2.0_db/3.0_db, 2.0_db/3.0_db, 2.0_db/3.0_db, 2.0_db/3.0_db, &
	2.0_db/3.0_db, 2.0_db/3.0_db, 2.0_db/3.0_db, 2.0_db/3.0_db, &
	2.0_db/3.0_db, 2.0_db/3.0_db, 2.0_db/3.0_db], [0.0_db, 0.0_db, &
	0.0_db, 0.0_db, 0.0_db, 0.0_db, 0.0_db, 1.0_db, 1.0_db, -1.0_db, &
	-1.0_db, 0.0_db, 0.0_db, 0.0_db, 0.0_db, 0.0_db, 0.0_db, 0.0_db, &
	0.0_db, 1.0_db, 1.0_db, -1.0_db, -1.0_db, 1.0_db, 1.0_db, -1.0_db, &
	-1.0_db], [0.0_db, 0.0_db, 0.0_db, 0.0_db, 0.0_db, 0.0_db, 0.0_db, &
	0.0_db, 0.0_db, 0.0_db, 0.0_db, 0.0_db, 0.0_db, 0.0_db, 0.0_db, &
	1.0_db, 1.0_db, -1.0_db, -1.0_db, 1.0_db, 1.0_db, 1.0_db, 1.0_db, &
	-1.0_db, -1.0_db, -1.0_db, -1.0_db], [0.0_db, 0.0_db, 0.0_db, 0.0_db, &
	0.0_db, 0.0_db, 0.0_db, 0.0_db, 0.0_db, 0.0_db, 0.0_db, 1.0_db, &
	1.0_db, -1.0_db, -1.0_db, 0.0_db, 0.0_db, 0.0_db, 0.0_db, 1.0_db, &
	1.0_db, -1.0_db, -1.0_db, -1.0_db, -1.0_db, 1.0_db, 1.0_db], [0.0_db, &
	0.0_db, 0.0_db, -1.0_db/3.0_db, 1.0_db/3.0_db, 0.0_db, 0.0_db, &
	2.0_db/3.0_db, -2.0_db/3.0_db, -2.0_db/3.0_db, 2.0_db/3.0_db, &
	-1.0_db/3.0_db, 1.0_db/3.0_db, -1.0_db/3.0_db, 1.0_db/3.0_db, 0.0_db, &
	0.0_db, 0.0_db, 0.0_db, 2.0_db/3.0_db, -2.0_db/3.0_db, &
	-2.0_db/3.0_db, 2.0_db/3.0_db, -2.0_db/3.0_db, 2.0_db/3.0_db, &
	-2.0_db/3.0_db, 2.0_db/3.0_db], [0.0_db, 0.0_db, 0.0_db, 0.0_db, &
	0.0_db, -1.0_db/3.0_db, 1.0_db/3.0_db, 0.0_db, 0.0_db, 0.0_db, &
	0.0_db, -1.0_db/3.0_db, 1.0_db/3.0_db, 1.0_db/3.0_db, -1.0_db/3.0_db, &
	2.0_db/3.0_db, -2.0_db/3.0_db, 2.0_db/3.0_db, -2.0_db/3.0_db, &
	2.0_db/3.0_db, -2.0_db/3.0_db, 2.0_db/3.0_db, -2.0_db/3.0_db, &
	2.0_db/3.0_db, -2.0_db/3.0_db, -2.0_db/3.0_db, 2.0_db/3.0_db], &
	[0.0_db, -1.0_db/3.0_db, 1.0_db/3.0_db, 0.0_db, 0.0_db, 0.0_db, &
	0.0_db, 2.0_db/3.0_db, -2.0_db/3.0_db, 2.0_db/3.0_db, -2.0_db/3.0_db, &
	0.0_db, 0.0_db, 0.0_db, 0.0_db, -1.0_db/3.0_db, 1.0_db/3.0_db, &
	1.0_db/3.0_db, -1.0_db/3.0_db, 2.0_db/3.0_db, -2.0_db/3.0_db, &
	2.0_db/3.0_db, -2.0_db/3.0_db, -2.0_db/3.0_db, 2.0_db/3.0_db, &
	2.0_db/3.0_db, -2.0_db/3.0_db], [0.0_db, -1.0_db/3.0_db, &
	1.0_db/3.0_db, 0.0_db, 0.0_db, 0.0_db, 0.0_db, -1.0_db/3.0_db, &
	1.0_db/3.0_db, -1.0_db/3.0_db, 1.0_db/3.0_db, 0.0_db, 0.0_db, 0.0_db, &
	0.0_db, 2.0_db/3.0_db, -2.0_db/3.0_db, -2.0_db/3.0_db, 2.0_db/3.0_db, &
	2.0_db/3.0_db, -2.0_db/3.0_db, 2.0_db/3.0_db, -2.0_db/3.0_db, &
	-2.0_db/3.0_db, 2.0_db/3.0_db, 2.0_db/3.0_db, -2.0_db/3.0_db], &
	[0.0_db, 0.0_db, 0.0_db, -1.0_db/3.0_db, 1.0_db/3.0_db, 0.0_db, &
	0.0_db, -1.0_db/3.0_db, 1.0_db/3.0_db, 1.0_db/3.0_db, -1.0_db/3.0_db, &
	2.0_db/3.0_db, -2.0_db/3.0_db, 2.0_db/3.0_db, -2.0_db/3.0_db, 0.0_db, &
	0.0_db, 0.0_db, 0.0_db, 2.0_db/3.0_db, -2.0_db/3.0_db, &
	-2.0_db/3.0_db, 2.0_db/3.0_db, -2.0_db/3.0_db, 2.0_db/3.0_db, &
	-2.0_db/3.0_db, 2.0_db/3.0_db], [0.0_db, 0.0_db, 0.0_db, 0.0_db, &
	0.0_db, -1.0_db/3.0_db, 1.0_db/3.0_db, 0.0_db, 0.0_db, 0.0_db, &
	0.0_db, 2.0_db/3.0_db, -2.0_db/3.0_db, -2.0_db/3.0_db, 2.0_db/3.0_db, &
	-1.0_db/3.0_db, 1.0_db/3.0_db, -1.0_db/3.0_db, 1.0_db/3.0_db, &
	2.0_db/3.0_db, -2.0_db/3.0_db, 2.0_db/3.0_db, -2.0_db/3.0_db, &
	2.0_db/3.0_db, -2.0_db/3.0_db, -2.0_db/3.0_db, 2.0_db/3.0_db], &
	[0.0_db, 0.0_db, 0.0_db, 0.0_db, 0.0_db, 0.0_db, 0.0_db, 0.0_db, &
	0.0_db, 0.0_db, 0.0_db, 0.0_db, 0.0_db, 0.0_db, 0.0_db, 0.0_db, &
	0.0_db, 0.0_db, 0.0_db, 1.0_db, -1.0_db, -1.0_db, 1.0_db, 1.0_db, &
	-1.0_db, 1.0_db, -1.0_db], [1.0_db/9.0_db, -2.0_db/9.0_db, &
	-2.0_db/9.0_db, -2.0_db/9.0_db, -2.0_db/9.0_db, 1.0_db/9.0_db, &
	1.0_db/9.0_db, 4.0_db/9.0_db, 4.0_db/9.0_db, 4.0_db/9.0_db, &
	4.0_db/9.0_db, -2.0_db/9.0_db, -2.0_db/9.0_db, -2.0_db/9.0_db, &
	-2.0_db/9.0_db, -2.0_db/9.0_db, -2.0_db/9.0_db, -2.0_db/9.0_db, &
	-2.0_db/9.0_db, 4.0_db/9.0_db, 4.0_db/9.0_db, 4.0_db/9.0_db, &
	4.0_db/9.0_db, 4.0_db/9.0_db, 4.0_db/9.0_db, 4.0_db/9.0_db, &
	4.0_db/9.0_db], [1.0_db/9.0_db, -2.0_db/9.0_db, -2.0_db/9.0_db, &
	1.0_db/9.0_db, 1.0_db/9.0_db, -2.0_db/9.0_db, -2.0_db/9.0_db, &
	-2.0_db/9.0_db, -2.0_db/9.0_db, -2.0_db/9.0_db, -2.0_db/9.0_db, &
	-2.0_db/9.0_db, -2.0_db/9.0_db, -2.0_db/9.0_db, -2.0_db/9.0_db, &
	4.0_db/9.0_db, 4.0_db/9.0_db, 4.0_db/9.0_db, 4.0_db/9.0_db, &
	4.0_db/9.0_db, 4.0_db/9.0_db, 4.0_db/9.0_db, 4.0_db/9.0_db, &
	4.0_db/9.0_db, 4.0_db/9.0_db, 4.0_db/9.0_db, 4.0_db/9.0_db], &
	[1.0_db/9.0_db, 1.0_db/9.0_db, 1.0_db/9.0_db, -2.0_db/9.0_db, &
	-2.0_db/9.0_db, -2.0_db/9.0_db, -2.0_db/9.0_db, -2.0_db/9.0_db, &
	-2.0_db/9.0_db, -2.0_db/9.0_db, -2.0_db/9.0_db, 4.0_db/9.0_db, &
	4.0_db/9.0_db, 4.0_db/9.0_db, 4.0_db/9.0_db, -2.0_db/9.0_db, &
	-2.0_db/9.0_db, -2.0_db/9.0_db, -2.0_db/9.0_db, 4.0_db/9.0_db, &
	4.0_db/9.0_db, 4.0_db/9.0_db, 4.0_db/9.0_db, 4.0_db/9.0_db, &
	4.0_db/9.0_db, 4.0_db/9.0_db, 4.0_db/9.0_db], [0.0_db, 0.0_db, &
	0.0_db, 0.0_db, 0.0_db, 0.0_db, 0.0_db, -1.0_db/3.0_db, &
	-1.0_db/3.0_db, 1.0_db/3.0_db, 1.0_db/3.0_db, 0.0_db, 0.0_db, 0.0_db, &
	0.0_db, 0.0_db, 0.0_db, 0.0_db, 0.0_db, 2.0_db/3.0_db, 2.0_db/3.0_db, &
	-2.0_db/3.0_db, -2.0_db/3.0_db, 2.0_db/3.0_db, 2.0_db/3.0_db, &
	-2.0_db/3.0_db, -2.0_db/3.0_db], [0.0_db, 0.0_db, 0.0_db, 0.0_db, &
	0.0_db, 0.0_db, 0.0_db, 0.0_db, 0.0_db, 0.0_db, 0.0_db, 0.0_db, &
	0.0_db, 0.0_db, 0.0_db, -1.0_db/3.0_db, -1.0_db/3.0_db, &
	1.0_db/3.0_db, 1.0_db/3.0_db, 2.0_db/3.0_db, 2.0_db/3.0_db, &
	2.0_db/3.0_db, 2.0_db/3.0_db, -2.0_db/3.0_db, -2.0_db/3.0_db, &
	-2.0_db/3.0_db, -2.0_db/3.0_db], [0.0_db, 0.0_db, 0.0_db, 0.0_db, &
	0.0_db, 0.0_db, 0.0_db, 0.0_db, 0.0_db, 0.0_db, 0.0_db, &
	-1.0_db/3.0_db, -1.0_db/3.0_db, 1.0_db/3.0_db, 1.0_db/3.0_db, 0.0_db, &
	0.0_db, 0.0_db, 0.0_db, 2.0_db/3.0_db, 2.0_db/3.0_db, -2.0_db/3.0_db, &
	-2.0_db/3.0_db, -2.0_db/3.0_db, -2.0_db/3.0_db, 2.0_db/3.0_db, &
	2.0_db/3.0_db], [0.0_db, 0.0_db, 0.0_db, 1.0_db/9.0_db, &
	-1.0_db/9.0_db, 0.0_db, 0.0_db, -2.0_db/9.0_db, 2.0_db/9.0_db, &
	2.0_db/9.0_db, -2.0_db/9.0_db, -2.0_db/9.0_db, 2.0_db/9.0_db, &
	-2.0_db/9.0_db, 2.0_db/9.0_db, 0.0_db, 0.0_db, 0.0_db, 0.0_db, &
	4.0_db/9.0_db, -4.0_db/9.0_db, -4.0_db/9.0_db, 4.0_db/9.0_db, &
	-4.0_db/9.0_db, 4.0_db/9.0_db, -4.0_db/9.0_db, 4.0_db/9.0_db], &
	[0.0_db, 0.0_db, 0.0_db, 0.0_db, 0.0_db, 1.0_db/9.0_db, &
	-1.0_db/9.0_db, 0.0_db, 0.0_db, 0.0_db, 0.0_db, -2.0_db/9.0_db, &
	2.0_db/9.0_db, 2.0_db/9.0_db, -2.0_db/9.0_db, -2.0_db/9.0_db, &
	2.0_db/9.0_db, -2.0_db/9.0_db, 2.0_db/9.0_db, 4.0_db/9.0_db, &
	-4.0_db/9.0_db, 4.0_db/9.0_db, -4.0_db/9.0_db, 4.0_db/9.0_db, &
	-4.0_db/9.0_db, -4.0_db/9.0_db, 4.0_db/9.0_db], [0.0_db, &
	1.0_db/9.0_db, -1.0_db/9.0_db, 0.0_db, 0.0_db, 0.0_db, 0.0_db, &
	-2.0_db/9.0_db, 2.0_db/9.0_db, -2.0_db/9.0_db, 2.0_db/9.0_db, 0.0_db, &
	0.0_db, 0.0_db, 0.0_db, -2.0_db/9.0_db, 2.0_db/9.0_db, 2.0_db/9.0_db, &
	-2.0_db/9.0_db, 4.0_db/9.0_db, -4.0_db/9.0_db, 4.0_db/9.0_db, &
	-4.0_db/9.0_db, -4.0_db/9.0_db, 4.0_db/9.0_db, 4.0_db/9.0_db, &
	-4.0_db/9.0_db], [-1.0_db/27.0_db, 2.0_db/27.0_db, 2.0_db/27.0_db, &
	2.0_db/27.0_db, 2.0_db/27.0_db, 2.0_db/27.0_db, 2.0_db/27.0_db, &
	-4.0_db/27.0_db, -4.0_db/27.0_db, -4.0_db/27.0_db, -4.0_db/27.0_db, &
	-4.0_db/27.0_db, -4.0_db/27.0_db, -4.0_db/27.0_db, -4.0_db/27.0_db, &
	-4.0_db/27.0_db, -4.0_db/27.0_db, -4.0_db/27.0_db, -4.0_db/27.0_db, &
	8.0_db/27.0_db, 8.0_db/27.0_db, 8.0_db/27.0_db, 8.0_db/27.0_db, &
	8.0_db/27.0_db, 8.0_db/27.0_db, 8.0_db/27.0_db, 8.0_db/27.0_db] &
    ], shape=[nlinks+1,nlinks+1])
#endif
      
#elif LATTICE == 19
#warning "LATTICE 19: the lattice D3Q19 is utilized"
   integer, parameter :: nlinks=18
   !lattice vectors
   integer, dimension(0:nlinks), parameter :: &
   !          0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18
      ex=   (/0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 0, 0, 0, 0, 1,-1,-1, 1/)
   integer, dimension(0:nlinks), parameter:: &
      ey=  (/ 0, 0, 0, 1,-1, 0, 0, 1,-1,-1, 1, 1,-1, 1,-1, 0, 0, 0, 0/)
   integer, dimension(0:nlinks), parameter:: &
      ez=  (/ 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 1,-1,-1, 1, 1,-1, 1,-1/)
   integer, dimension(0:nlinks), parameter:: &
      opp =(/ 0, 2, 1, 4, 3, 6, 5, 8, 7,10, 9,12,11,14,13,16,15,18,17/)

   ! D3Q19 weights: rest direction 0 (p0), axis directions 1-6 (p1), face diagonals 7-18 (p2)
   real(kind=db), parameter :: p0 = real(1.d0/3.d0 , kind=db)
   real(kind=db), parameter :: p1 = real(1.d0/18.d0, kind=db)  
   real(kind=db), parameter :: p2 = real(1.d0/36.d0, kind=db) 
   real(kind=db), dimension(0:nlinks), parameter :: &
      p=(/p0,p1,p1,p1,p1,p1,p1,p2,p2,p2,p2,p2,p2,p2,p2,p2,p2,p2,p2/)
#elif LATTICE == 15
#warning "LATTICE 15: the lattice D3Q15 is utilized"
   integer, parameter :: nlinks=14
      !lattice vectors
   integer, dimension(0:nlinks), parameter :: &
   !          0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14
      ex=   (/0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1,-1, 1, 1,-1/)
   integer, dimension(0:nlinks), parameter:: &
      ey=  (/ 0, 0, 0, 1,-1, 0, 0, 1,-1,-1, 1,-1, 1,-1, 1/)
   integer, dimension(0:nlinks), parameter:: &
      ez=  (/ 0, 0, 0, 0, 0, 1,-1, 1,-1, 1,-1, 1,-1,-1, 1/)
   integer, dimension(0:nlinks), parameter:: &
      opp =(/ 0, 2, 1, 4, 3, 6, 5, 8, 7,10, 9,12,11,14,13/)
      
   ! D3Q15 weights: rest direction 0 (p0), axis directions 1-6 (p1), body diagonals 7-14 (p2)
   real(kind=db), parameter :: p0 = real(2.d0/9.d0 , kind=db)
   real(kind=db), parameter :: p1 = real(1.d0/9.d0 , kind=db)
   real(kind=db), parameter :: p2 = real(1.d0/72.d0, kind=db)
   real(kind=db), dimension(0:nlinks), parameter :: &
      p=(/p0,p1,p1,p1,p1,p1,p1,p2,p2,p2,p2,p2,p2,p2,p2/)
#else
#error "LATTICE not supported"
#endif
   
   ! D3Q27 weights, defined for every lattice: lb_cuda_auxfields_module.f90 uses them in the
   ! isotropic finite-difference stencils of the phase field (gradient, Laplacian)
   real(kind=db), parameter :: p0d3q27=real(8.d0/27.d0,kind=db)
   real(kind=db), parameter :: p1d3q27=real(2.d0/27.d0,kind=db)
   real(kind=db), parameter :: p2d3q27=real(1.d0/54.d0,kind=db)
   real(kind=db), parameter :: p3d3q27=real(1.d0/216.d0,kind=db)

   real(kind=db), dimension(0:nlinks), parameter :: dex=real(ex,kind=db)
   real(kind=db), dimension(0:nlinks), parameter :: dey=real(ey,kind=db)
   real(kind=db), dimension(0:nlinks), parameter :: dez=real(ez,kind=db)


   integer(kind=isf), allocatable,dimension(:,:,:)   :: isfluid
   
   integer :: flip,flop

   real(kind=4), allocatable, dimension(:,:,:) :: arr_3d

   ! Field arrays, stored tile by tile as (TILE_DIMx,TILE_DIMy,TILE_DIMz,ncomp,nblocks) in the storage
   ! kind strdb. The flip and flop arrays are the two time levels of the moments and of the phase
   ! field, used alternately as input and output of the time step.
   real(kind=strdb), allocatable, dimension(:,:,:,:,:) :: hfields_flip,hfields_flop  !allocate hydro fields flip and flop
   real(kind=strdb), allocatable, dimension(:,:,:,:,:) :: phifields_flip,phifields_flop  !allocate phi fields flip and flop
   real(kind=strdb), allocatable, dimension(:,:,:,:,:) :: forces
   real(kind=strdb), allocatable, dimension(:,:,:,:,:) :: auxfields !allocate aux fields
   real(kind=strdb), allocatable, dimension(:,:,:,:,:) :: locauxfields !allocate aux fields
   
   
   ! Numbers of components of the field arrays:
   ! hfields: 1 press, the zeroth moment sum_l f_l (pressure-like variable p*); 2-4 velocity u, v, w
   !   (first moment); 5-10 second-order Hermite moments pxx, pyy, pzz, pxy, pxz, pyz, that is
   !   sum_l (c_la c_lb - cs^2 delta_ab) f_l;
   ! forces: x, y, z components of the force per unit volume used in the collision;
   ! phifields: phase field phi (allocated only with TWOCOMPONENT);
   ! auxfields: 1-3 unit normal n = grad phi/|grad phi|, 4 |grad phi|, 5-7 phi (1-phi) n;
   ! locauxfields: 1 Laplacian of phi, 2 divergence of phi (1-phi) n; with REPULSIVE_FLUX also
   !   3-5 indices of the paired node and 6-8 repulsive flux vector (lb_cuda_repulsive_module.f90).
   integer, parameter :: nhfields=10
   integer, parameter :: nforces=3
   integer, parameter :: nphifields=1
#ifdef TWOCOMPONENT
   integer, parameter :: nauxfields=7    ! 3 norm unit vec ! 1 modgrad ! 3 arr_ 
#ifdef REPULSIVE_FLUX
   integer, parameter :: nlocauxfields=8 !1 lap_phi !1 div_thetan ! 3 pair_i ! 3 Jx vector 
#else
   integer, parameter :: nlocauxfields=2 !1 lap_phi !1 div_thetan
#endif
   
#else   
   integer, parameter :: nauxfields=0    ! 
   integer, parameter :: nlocauxfields=0 !
#endif   
   
   ! total numbers of elements of the field arrays (set in allocate_struct)
   integer, save :: ntothfields
   integer, save :: ntotforces
   integer, save :: ntotphifields
   integer, save :: ntotauxfields
   integer, save :: ntotlocauxfields
   
#ifdef MULTIHIT
   real(kind=db), allocatable, dimension(:,:,:) ::ABCx,ABCy,ABCz
#endif
#ifdef ELASTIC_FORCE
   real(kind=db), allocatable, dimension(:,:,:) :: u_ref,v_ref,w_ref
   real(kind=db) :: lambda_rel,k_elastic
#endif 
#ifdef MONOD
	real(kind=db) :: mu_max,Ks
#endif

#ifdef REPULSIVE_FLUX
	integer(kind=isf), allocatable,dimension(:,:,:) :: rep_mask
	real(kind=db) :: q_th,win,cosOppT,pwr,A_rep	
#endif

   real(kind=db) :: global_phi_sum=ZERO,dphi=ZERO,corr=ZERO,global_phi_sum_ini=ZERO,global_phi_sum_new=ZERO
   real(kind=db) :: global_phi_change=ZERO,global_phi_change_new=ZERO
   integer :: global_count=0, global_count_new=0
   ! global correction of phi with a Lagrange multiplier at internal obstacles (fluid namelist,
   ! bcs_mesoscopic_phifields): llagrange switches it on, lagrange_every is the step interval
   logical :: llagrange=.true.
   integer :: lagrange_every=100
   real(kind=db) :: mymemory,totmemory
   real(kind=db) :: uwall
   real(kind=db) :: rrx,rry,rrz
   real(kind=db) :: rho_r=ONE,rho_b=ONE
   real(kind=db) :: invrho_r=ONE,invrho_b=ONE
   real(kind=db) :: time_init,time_actual,time_limit=-ONE
   real(kind=db) :: time_actual_old
   integer :: every_time_check=1000
   integer(c_int) :: mydev_c
#ifdef _NVML  
   integer :: num_p_w=0
   integer(c_int) :: p_mw
   real(kind=db) :: p_w=ZERO
   real(kind=db) :: tot_energy=ZERO
   real(kind=db) :: step_energy
   integer(c_long_long) :: energy_1,energy_2
#endif
   
   !****************************print vars**************************************!

   integer, parameter :: mxln=256
   character(len=mxln) :: inipFile
   character(len=8), allocatable, dimension(:) :: namevarvtk
   character(len=500), allocatable, dimension(:) :: headervtk
   character(len=30), allocatable, dimension(:) :: footervtk
   integer, allocatable, dimension(:) :: ndimvtk
   integer, allocatable, dimension(:) :: vtkoffset
   integer, allocatable, dimension(:) :: ndatavtk
   integer, allocatable, dimension(:) :: nheadervtk
   integer :: nfilevtk
   integer, allocatable, dimension(:) :: varlistvtk
   character :: delimiter
   character(len=*), parameter :: filenamevtk='out'

   real(kind=printdb), allocatable, dimension(:,:,:) :: rhoprint
   real(kind=printdb), allocatable, dimension(:,:,:) :: pressprint
   real(kind=printdb), allocatable, dimension(:,:,:,:) :: velprint
   logical :: lelittle
   character(len=mxln), save :: dir_out
   character(len=mxln) :: extentvtk
   character(len=mxln) :: sevt1,sevt2,sevt3,arg,directive
   character(len=1), allocatable, dimension(:) :: head1,head2,head3

contains

   ! Standard normal deviate computed as a pure function of the integers i, j, k, l, without stored
   ! state, so that it can be called from device code. Two uniform deviates are built from the seeds
   ! with the initialization sequence of the Marsaglia-Zaman (RANMAR) generator, the second one with
   ! l shifted by nint(50 u1), and combined with the Box-Muller transformation.
   pure function gauss_noseeded(i,j,k,l)
      !$acc routine seq
      !the directive above tells OpenACC to create a device copy; gang means that the routine can be called independently by several threads/vectors/workers/gangs
      implicit none
      integer, intent(in) :: i,j,k,l
      integer :: kk,ll
      real(kind=db) :: gauss_noseeded
      real(kind=db) :: dtemp1,dtemp2
      real(kind=db), parameter :: mylimit=real(1.d-30,kind=db)
      real(kind=db), parameter :: FIFTY = real(50.d0,kind=db)
      real(kind=db), parameter :: ONE = real(1.d0,kind=db)
      real(kind=db), parameter :: TWO = real(2.d0,kind=db)
      real(kind=db), parameter :: Pi = real(3.1415926535897932384626433832795028841971d0,kind=db)
	  !****************************************************
	  integer :: isub,jsub,ksub,lsub,msub
      integer ::ii,jj
      real(4) :: s,t,u33,u97,csub,uni
      real(4), parameter :: c =  362436.0/16777216.0
      real(4), parameter :: cd= 7654321.0/16777216.0
      real(4), parameter :: cm=16777213.0/16777216.0
	  
	  ! initial values of i,j,k must be in range 1 to 178 (not all 1)
      ! initial value of l must be in range 0 to 168.
      isub=mod(i,178)
      isub=isub+1

      jsub=mod(j,178)
      jsub=jsub+1

      ksub=mod(k,178)
      ksub=ksub+1

      lsub=mod(l,169)

      ! initialization on fly
      ii=97
      s=0.0
      t=0.5
      do jj=1,24
         msub=mod(mod(isub*jsub,179)*ksub,179)
         isub=jsub
         jsub=ksub
         ksub=msub
         lsub=mod(53*lsub+1,169)
         if(mod(lsub*msub,64).ge.32)s=s+t
         t=0.5*t
      enddo
      u97=s

      ii=33
      s=0.0
      t=0.5
      do jj=1,24
         msub=mod(mod(isub*jsub,179)*ksub,179)
         isub=jsub
         jsub=ksub
         ksub=msub
         lsub=mod(53*lsub+1,169)
         if(mod(lsub*msub,64).ge.32)s=s+t
         t=0.5*t
      enddo
      u33=s
      uni=u97-u33
      if (uni.lt.0.0) uni = uni + 1.0
      csub = c-cd
      if (csub.lt.0.0) csub = csub + cm
      uni = uni-csub
      if (uni.lt.0.0) uni = uni+1.0
      dtemp1 = real(uni,kind=db)
	  kk=nint(dtemp1*FIFTY)
	  !*********************************************
      ll=l+kk
      
      isub=mod(i,178)
      isub=isub+1

      jsub=mod(j,178)
      jsub=jsub+1

      ksub=mod(k,178)
      ksub=ksub+1

      lsub=mod(ll,169)

      ! initialization on fly
      ii=97
      s=0.0
      t=0.5
      do jj=1,24
         msub=mod(mod(isub*jsub,179)*ksub,179)
         isub=jsub
         jsub=ksub
         ksub=msub
         lsub=mod(53*lsub+1,169)
         if(mod(lsub*msub,64).ge.32)s=s+t
         t=0.5*t
      enddo
      u97=s

      ii=33
      s=0.0
      t=0.5
      do jj=1,24
         msub=mod(mod(isub*jsub,179)*ksub,179)
         isub=jsub
         jsub=ksub
         ksub=msub
         lsub=mod(53*lsub+1,169)
         if(mod(lsub*msub,64).ge.32)s=s+t
         t=0.5*t
      enddo
      u33=s
      uni=u97-u33
      if (uni.lt.0.0) uni = uni + 1.0
      csub = c-cd
      if (csub.lt.0.0) csub = csub + cm
      uni = uni-csub
      if (uni.lt.0.0) uni = uni+1.0
      dtemp2 = real(uni,kind=db)

      dtemp1=dtemp1*(ONE-mylimit)+mylimit

      ! Box-Muller transformation
      gauss_noseeded=sqrt(- TWO *log(dtemp1))*cos(TWO*pi*dtemp2)
   end function gauss_noseeded

   ! Uniform deviate in [0,1) computed as a pure function of i, j, k, l (see gauss_noseeded).
   pure function rand_noseeded(i,j,k,l)
      !$acc routine seq !!!not gang
      !the directive above tells OpenACC to create a device copy; gang means that the routine can be called independently by several threads/vectors/workers/gangs
      implicit none
      integer, intent(in) :: i,j,k,l
      integer :: isub,jsub,ksub,lsub,msub
      integer ::ii,jj
      real(4) :: s,t,u33,u97,csub,uni
      real(kind=db) :: rand_noseeded

      real(4), parameter :: c =  362436.0/16777216.0
      real(4), parameter :: cd= 7654321.0/16777216.0
      real(4), parameter :: cm=16777213.0/16777216.0


      ! initial values of i,j,k must be in range 1 to 178 (not all 1)
      ! initial value of l must be in range 0 to 168.
      isub=mod(i,178)
      isub=isub+1

      jsub=mod(j,178)
      jsub=jsub+1

      ksub=mod(k,178)
      ksub=ksub+1

      lsub=mod(l,169)

      ! initialization on fly
      ii=97
      s=0.0
      t=0.5
      do jj=1,24
         msub=mod(mod(isub*jsub,179)*ksub,179)
         isub=jsub
         jsub=ksub
         ksub=msub
         lsub=mod(53*lsub+1,169)
         if(mod(lsub*msub,64).ge.32)s=s+t
         t=0.5*t
      enddo
      u97=s

      ii=33
      s=0.0
      t=0.5
      do jj=1,24
         msub=mod(mod(isub*jsub,179)*ksub,179)
         isub=jsub
         jsub=ksub
         ksub=msub
         lsub=mod(53*lsub+1,169)
         if(mod(lsub*msub,64).ge.32)s=s+t
         t=0.5*t
      enddo
      u33=s
      uni=u97-u33
      if (uni.lt.0.0) uni = uni + 1.0
      csub = c-cd
      if (csub.lt.0.0) csub = csub + cm
      uni = uni-csub
      if (uni.lt.0.0) uni = uni+1.0
      rand_noseeded = real(uni,kind=db)
   end function rand_noseeded
   
   ! Seeds the intrinsic random_number generator with myseed, or with the system clock if myseed is
   ! absent, multiplied by myranksub+1 so that the MPI ranks get different sequences.
   subroutine init_random_seed_CPU(myseed,myranksub)
  
   implicit none
  
   integer,intent(in),optional :: myseed,myranksub
   integer :: i, n, clock
   integer :: idrank=0
  
   integer, allocatable :: seed(:)
   
   if(present(myranksub))idrank=myranksub
          
   call random_seed(size = n)
  
   allocate(seed(n))
  
   if(present(myseed))then
!    If the seed is given in input
     seed = myseed*(idrank+1) + 37 * (/ (i - 1, i = 1, n) /)
    
   else
!    If the seed is not given in input it is generated by the clock
     call system_clock(count=clock)
         
     seed = clock*(idrank+1) + 37 * (/ (i - 1, i = 1, n) /)
    
   endif
  
   call random_seed(put = seed)
       
   deallocate(seed)
  
   return
 
  end subroutine init_random_seed_CPU
  
  ! Fills gauss(1:nelement) with standard normal deviates of the intrinsic generator (Box-Muller
  ! transformation); myseed is not used.
  subroutine gauss_CPU(nelement,gauss,myseed)
  
  implicit none
  
   integer, intent(in) :: nelement
   real(kind=db), dimension(nelement) :: gauss
   integer,intent(in),optional :: myseed
  
   real(kind=db), parameter :: mylimit=real(1.d-30,kind=db)
   real(kind=db), parameter :: FIFTY = real(50.d0,kind=db)
   real(kind=db), parameter :: ONE = real(1.d0,kind=db)
   real(kind=db), parameter :: TWO = real(2.d0,kind=db)
   real(kind=db), parameter :: Pi = real(3.1415926535897932384626433832795028841971d0,kind=db)
  
   real(kind=db), dimension(nelement) :: dtemp1,dtemp2
   integer :: seedsub
  
   call random_number(dtemp1)
   call random_number(dtemp2)
  
   dtemp1=dtemp1*(ONE-mylimit)+mylimit
  
   ! Box-Muller transformation
   gauss(1:nelement)=sqrt(- TWO *log(dtemp1(1:nelement)))*cos(TWO*pi*dtemp2(1:nelement))
  
  
   end subroutine gauss_CPU
   
  ! One standard normal deviate of the intrinsic generator (Box-Muller transformation).
  function randgauss_CPU()
  
   implicit none
   
   real(kind=db) :: randgauss_CPU
  
   real(kind=db), parameter :: mylimit=real(1.d-30,kind=db)
   real(kind=db), parameter :: FIFTY = real(50.d0,kind=db)
   real(kind=db), parameter :: ONE = real(1.d0,kind=db)
   real(kind=db), parameter :: TWO = real(2.d0,kind=db)
   real(kind=db), parameter :: Pi = real(3.1415926535897932384626433832795028841971d0,kind=db)
  
   real(kind=db) :: dtemp1,dtemp2
  
   call random_number(dtemp1)
   call random_number(dtemp2)
  
   dtemp1=dtemp1*(ONE-mylimit)+mylimit
  
   ! Box-Muller transformation
   randgauss_CPU=sqrt(- TWO *log(dtemp1))*cos(TWO*pi*dtemp2)
  
  
   end function randgauss_CPU

   ! Copies the first nstring characters of mystring into the array mychar, allocated here.
   subroutine string_char(mychar,nstring,mystring)

      implicit none

      integer :: i
      character(1), allocatable, dimension(:) :: mychar
      integer, intent(in) :: nstring
      character(len=*), intent(in) :: mystring

      allocate(mychar(nstring))

      do i=1,nstring
         mychar(i)=mystring(i:i)
      enddo

   end subroutine string_char

   function space_fmtnumb(inum)

      !***********************************************************************
      !
      !     LBsoft function for returning the string of six characters
      !     with integer digits and leading spaces to the left
      !     originally written in JETSPIN by M. Lauricella et al.
      !
      !     licensed under Open Software License v. 3.0 (OSL-3.0)
      !     author: M. Lauricella
      !     last modification October 2019
      !
      !***********************************************************************

      implicit none

      integer,intent(in) :: inum
      character(len=6) :: space_fmtnumb
      integer :: numdigit,irest
      real(kind=8) :: tmp
      character(len=22) :: cnumberlabel

      numdigit=dimenumb(inum)
      irest=6-numdigit
      if(irest>0)then
         write(cnumberlabel,"(a,i8,a,i8,a)")"(a",irest,",i",numdigit,")"
         write(space_fmtnumb,fmt=cnumberlabel)repeat(' ',irest),inum
      else
         write(cnumberlabel,"(a,i8,a)")"(i",numdigit,")"
         write(space_fmtnumb,fmt=cnumberlabel)inum
      endif

      return

   end function space_fmtnumb

   function space_fmtnumb12(inum)

      !***********************************************************************
      !
      !     LBsoft function for returning the string of six characters
      !     with integer digits and leading TWELVE spaces to the left
      !     originally written in JETSPIN by M. Lauricella et al.
      !
      !     licensed under Open Software License v. 3.0 (OSL-3.0)
      !     author: M. Lauricella
      !     last modification October 2019
      !
      !***********************************************************************

      implicit none

      integer,intent(in) :: inum
      character(len=12) :: space_fmtnumb12
      integer :: numdigit,irest
      real(kind=8) :: tmp
      character(len=22) :: cnumberlabel

      numdigit=dimenumb(inum)
      irest=12-numdigit
      if(irest>0)then
         write(cnumberlabel,"(a,i8,a,i8,a)")"(a",irest,",i",numdigit,")"
         write(space_fmtnumb12,fmt=cnumberlabel)repeat(' ',irest),inum
      else
         write(cnumberlabel,"(a,i8,a)")"(i",numdigit,")"
         write(space_fmtnumb12,fmt=cnumberlabel)inum
      endif

      return

   end function space_fmtnumb12

   function dimenumb(inum)

      !***********************************************************************
      !
      !     LBsoft function for returning the number of digits
      !     of an integer number
      !     originally written in JETSPIN by M. Lauricella et al.
      !
      !     licensed under the 3-Clause BSD License (BSD-3-Clause)
      !     author: M. Lauricella
      !     last modification July 2018
      !
      !***********************************************************************

      implicit none

      integer,intent(in) :: inum
      integer :: dimenumb
      integer :: i
      real(kind=db) :: tmp

      i=1
      tmp=real(inum,kind=db)
      do
         if(tmp< 10.0_db )exit
         i=i+1
         tmp=tmp/ 10.0_db
      enddo

      dimenumb=i

      return

   end function dimenumb

   pure function write_fmtnumb(inum)

      !***********************************************************************
      !
      !     LBsoft function for returning the string of six characters
      !     with integer digits and leading zeros to the left
      !     originally written in JETSPIN by M. Lauricella et al.
      !
      !     licensed under the 3-Clause BSD License (BSD-3-Clause)
      !     author: M. Lauricella
      !     last modification July 2018
      !
      !***********************************************************************

      implicit none

      integer,intent(in) :: inum
      character(len=6) :: write_fmtnumb
      integer :: tmp
      
      tmp = max(0, min(999999, inum))  
      write(write_fmtnumb,'(I6.6)') tmp
      
      return
   end function write_fmtnumb

   function write_fmtnumb2(inum)

      !***********************************************************************
      !
      !     LBsoft function for returning the string of six characters
      !     with integer digits and leading zeros to the left
      !     originally written in JETSPIN by M. Lauricella et al.
      !
      !     licensed under the 3-Clause BSD License (BSD-3-Clause)
      !     author: M. Lauricella
      !     last modification July 2018
      !
      !***********************************************************************

      implicit none

      integer,intent(in) :: inum
      character(len=2) :: write_fmtnumb2
      integer :: tmp
      
      tmp = max(0, min(99, inum))  
      write(write_fmtnumb2,'(I2.2)') tmp

      return
   end function write_fmtnumb2
   
   ! Smooth cut-off: 1 for r <= inner_cut, 0 for r > outer_cut, and a cosine ramp in between.
   pure function fcut(r,inner_cut,outer_cut)
      !$acc routine seq !!!not gang
      !the directive above tells OpenACC to create a device copy; gang means that the routine can be called independently by several threads/vectors/workers/gangs
!***********************************************************************
!
!     LBsoft function for fading an observable (r) within a given
!     interval
!
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification January 2018
!
!***********************************************************************

  implicit none

  real(kind=db), intent(in) :: r,inner_cut,outer_cut
  real(kind=db) :: fcut
  
  real(kind=db),parameter :: Pi=3.14159265359793234626433

  if ( r <= inner_cut ) then
    fcut = 1.0_db
  elseif ( r > outer_cut ) then
      fcut = 0.0_db
  else
      fcut = 0.5_db*cos((r-inner_cut)*Pi/(outer_cut-inner_cut))+0.5_db
  endif

  return

 end function fcut
 
 ! Smooth step 0.5 (1 - tanh(2 (r - r_cut)/width_cut)): 1 for r << r_cut, 1/2 at r_cut, 0 for
 ! r >> r_cut.
 pure function fcut_tanh(r,r_cut,width_cut)
      !$acc routine seq !!!not gang
      !the directive above tells OpenACC to create a device copy; gang means that the routine can be called independently by several threads/vectors/workers/gangs
!***********************************************************************
!
!     LBsoft function for fading an observable (r) within a given
!     interval
!
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification January 2018
!
!***********************************************************************

  implicit none

  real(kind=db), intent(in) :: r,r_cut,width_cut
  real(kind=db) :: fcut_tanh
  
  fcut_tanh = 0.5_db * ( 1.0_db - tanh( (r-r_cut) / (0.5_db*width_cut) )) 

  return

 end function fcut_tanh

! Remainder of n/m in [0, m) for m > 0, as the intrinsic modulo; callable from device code.
function my_mod(n, m) result(res)
 !$acc routine seq
    implicit none
    integer, intent(in) :: n, m
    integer :: res

    res = n - (n / m) * m
    if (res < 0) then
        res = res + m
    endif
end function my_mod


elemental function coordblock(idblock,nxblock,nxyblock)
  !$acc routine seq
  !return block coordinate from id block (idblock start from 1 so we apply minus 1)
     implicit none
     integer, intent(in) :: idblock,nxblock,nxyblock
     integer, dimension(3) :: coordblock
   
     coordblock(3)=(idblock-1)/nxyblock +1
     coordblock(2)=((idblock-1)-(coordblock(3)-1)*nxyblock)/nxblock +1
     coordblock(1)=(idblock-1)-(coordblock(3)-1)*nxyblock-(coordblock(2)-1)*nxblock +1
      
end function coordblock

endmodule
