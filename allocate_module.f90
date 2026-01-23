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
      
      allocate(press_flip(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff))
      allocate(u_flip(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff))
      allocate(v_flip(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff))
      allocate(w_flip(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff))
      allocate(pxx_flip(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff))
      allocate(pyy_flip(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff))
      allocate(pzz_flip(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff))
      allocate(pxy_flip(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff))
      allocate(pxz_flip(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff))
      allocate(pyz_flip(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff))
      
      allocate(press_flop(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff))
      allocate(u_flop(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff))
      allocate(v_flop(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff))
      allocate(w_flop(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff))
      allocate(pxx_flop(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff))
      allocate(pyy_flop(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff))
      allocate(pzz_flop(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff))
      allocate(pxy_flop(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff))
      allocate(pxz_flop(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff))
      allocate(pyz_flop(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff))
      
      press_flip(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff)=ZERO
      u_flip(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff)=ZERO
      v_flip(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff)=ZERO
      w_flip(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff)=ZERO
      pxx_flip(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff)=ZERO
      pyy_flip(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff)=ZERO
      pzz_flip(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff)=ZERO
      pxy_flip(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff)=ZERO
      pxz_flip(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff)=ZERO
      pyz_flip(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff)=ZERO
      
      press_flop(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff)=ZERO
      u_flop(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff)=ZERO
      v_flop(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff)=ZERO
      w_flop(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff)=ZERO
      pxx_flop(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff)=ZERO
      pyy_flop(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff)=ZERO
      pzz_flop(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff)=ZERO
      pxy_flop(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff)=ZERO
      pxz_flop(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff)=ZERO
      pyz_flop(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff)=ZERO
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ntothfields=TILE_DIMx*TILE_DIMy*TILE_DIMz*nhfields*nblocks

      
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
