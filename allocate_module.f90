#include "defines.h"
module allocate_arrays
   use vars
   use mpi_template, only: nbuff,myrank,dostop,doerror
   implicit none

contains

   !*******************************************************************!
   subroutine allocate_struct
      implicit none
      
      integer :: mydim
      
      !i used this also in single component to write and read the restart file
      allocate(arr_3d(1:nx,1:ny,1:nz))
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ntothfields=TILE_DIMx*TILE_DIMy*TILE_DIMz*nhfields*nblocks
      allocate(hfields_flip(ntothfields))
      allocate(hfields_flop(ntothfields))
      
      ntotphifields=TILE_DIMx*TILE_DIMy*TILE_DIMz*nphifields*nblocks
#ifdef TWOCOMPONENT       
      allocate(phifields_flip(ntotphifields))
      allocate(phifields_flop(ntotphifields))
#endif      
      ntotauxfields=TILE_DIMx*TILE_DIMy*TILE_DIMz*nauxfields*nblocks
      allocate(auxfields(ntotauxfields))
      
      ntotlocauxfields=TILE_DIMx*TILE_DIMy*TILE_DIMz*nlocauxfields*nblocks
      allocate(locauxfields(ntotlocauxfields))
      
      ntotforces=TILE_DIMx*TILE_DIMy*TILE_DIMz*nforces*nblocks
      allocate(forces(ntotforces))
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
	  allocate(isfluid(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff))
	  
	  
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
