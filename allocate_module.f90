#include "defines.h"
module allocate_arrays
   use vars
   use mpi_template, only: nbuff,myrank,dostop,doerror
   implicit none

contains

   !*******************************************************************!
   ! allocate_struct: allocates the field arrays. setup_cuda must be called first, since it sets
   ! nxblock, nyblock, nzblock, nxyblock and nblocks.
   !
   ! Tiled layout: the field arrays are dimensioned (TILE_DIMx,TILE_DIMy,TILE_DIMz,nfields,nblocks).
   ! The local subdomain is split into tiles of TILE_DIMx x TILE_DIMy x TILE_DIMz nodes and is
   ! surrounded by one layer of halo tiles, so nxblock = nx/TILE_DIMx + 2 (and so on) and
   ! nblocks = nxblock*nyblock*nzblock. The node (i,j,k), with i = 1-TILE_DIMx,...,nx+TILE_DIMx
   ! (and so on), is stored at (ii,jj,kk,field,myblock) with
   !    xblock = (i+2*TILE_DIMx-1)/TILE_DIMx,   ii = i-(xblock-2)*TILE_DIMx   (same for y and z),
   !    myblock = xblock + (yblock-1)*nxblock + (zblock-1)*nxyblock,
   ! so the tiles of the subdomain have xblock = 2,...,nxblock-1.
   !
   ! Fields (index nfields):
   !   hfields_flip, hfields_flop (nhfields = 10): 1 = p* (zeroth moment), 2-4 = u, v, w,
   !     5-10 = second-order Hermite moments pxx, pyy, pzz, pxy, pxz, pyz, that is
   !     sum_l (c_a c_b - cs^2 delta_ab) f_l, equal to u_a u_b at equilibrium. The time steps read
   !     one buffer and write the other (see integrator_module.f90).
   !   phifields_flip, phifields_flop (nphifields = 1, TWOCOMPONENT only): phase field phi.
   !   auxfields (nauxfields = 7 with TWOCOMPONENT, else 0): 1-3 = unit normal n, 4 = |grad phi|,
   !     5-7 = phi (1-phi) n; exchanged between ranks.
   !   locauxfields (TWOCOMPONENT): 1 = Laplacian of phi, 2 = div(phi (1-phi) n), 3-8 used by
   !     REPULSIVE_FLUX; not exchanged.
   !   forces (nforces = 3): force on each node, computed by moments_LB_kernel for the fused kernel.
   subroutine allocate_struct
      implicit none
      
      integer :: mydim
      
      ! arr_3d: single-precision copy of one scalar field of the subdomain, used to read and write
      ! RAW files.
      !i used this also in single component to write and read the restart file
      allocate(arr_3d(1:nx,1:ny,1:nz))
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! ntot*: total number of elements of each array, passed to the kernels.
      ntothfields=TILE_DIMx*TILE_DIMy*TILE_DIMz*nhfields*nblocks
      allocate(hfields_flip(TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields,nblocks))
      allocate(hfields_flop(TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields,nblocks))
      
      ntotphifields=TILE_DIMx*TILE_DIMy*TILE_DIMz*nphifields*nblocks
#ifdef TWOCOMPONENT       
      allocate(phifields_flip(TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields,nblocks))
      allocate(phifields_flop(TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields,nblocks))
#endif      
      ntotauxfields=TILE_DIMx*TILE_DIMy*TILE_DIMz*nauxfields*nblocks
      allocate(auxfields(TILE_DIMx,TILE_DIMy,TILE_DIMz,nauxfields,nblocks))
      
      ntotlocauxfields=TILE_DIMx*TILE_DIMy*TILE_DIMz*nlocauxfields*nblocks
      allocate(locauxfields(TILE_DIMx,TILE_DIMy,TILE_DIMz,nlocauxfields,nblocks))
      
      ntotforces=TILE_DIMx*TILE_DIMy*TILE_DIMz*nforces*nblocks
      allocate(forces(TILE_DIMx,TILE_DIMy,TILE_DIMz,nforces,nblocks))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      
      
#ifdef TWOCOMPONENT
	  
#ifdef REPULSIVE_FLUX
	allocate(rep_mask(1:nx,1:ny,1:nz))
#endif

#endif

#ifdef MULTIHIT
	  allocate(ABCx(1:nx,1:ny,1:nz),ABCy(1:nx,1:ny,1:nz),ABCz(1:nx,1:ny,1:nz))
#endif

#ifdef ELASTIC_FORCE
	 allocate(u_ref(1:nx,1:ny,1:nz),v_ref(1:nx,1:ny,1:nz),w_ref(1:nx,1:ny,1:nz))
#endif
      !allocate(isfluid(0:nx+1,0:ny+1,0:nz+1))
      ! isfluid has nbuff halo layers (3 with REPULSIVE_FLUX, else 1): 1 = fluid, 0 = solid,
      ! -1 = fluid node next to a non-fluid node.
	  allocate(isfluid(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff))
	  
	  
      ! Output arrays on the subsampled grid (nxskip = nx/stepskip, and so on).
      if(lprint)then
         allocate(rhoprint(1:nxskip,1:nyskip,1:nzskip))
         allocate(velprint(1:3,1:nxskip,1:nyskip,1:nzskip))
         rhoprint(1:nxskip,1:nyskip,1:nzskip)=ZERO
         velprint(1:3,1:nxskip,1:nyskip,1:nzskip)=ZERO
#if defined(TWOCOMPONENT) && defined(WRITEPRESS)
         allocate(pressprint(1:nxskip,1:nyskip,1:nzskip))
         pressprint(1:nxskip,1:nyskip,1:nzskip)=ZERO
#endif
      endif
      
      
     
   endsubroutine
endmodule
