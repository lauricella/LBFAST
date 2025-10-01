#include "defines.h"
module allocate_arrays
   use vars
   use mpi_template, only: nbuff,myrank,dostop,doerror
   implicit none

contains

   !*******************************************************************!
   subroutine allocate_struct
      implicit none
      allocate(f(0:nx+1,0:ny+1,0:nz+1,0:nlinks))
      allocate(rho(0:nx+1,0:ny+1,0:nz+1))
      allocate(u(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff))
      allocate(v(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff))
      allocate(w(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff))
      allocate(pxx(0:nx+1,0:ny+1,0:nz+1),pxy(0:nx+1,0:ny+1,0:nz+1),pxz(0:nx+1,0:ny+1,0:nz+1),pyy(0:nx+1,0:ny+1,0:nz+1))
      allocate(pyz(0:nx+1,0:ny+1,0:nz+1),pzz(0:nx+1,0:ny+1,0:nz+1))
      !i used this also in single component to write and read the restart file
      allocate(lap_phi(1:nx,1:ny,1:nz))
      allocate(fux(1:nx,1:ny,1:nz))
      allocate(fvy(1:nx,1:ny,1:nz))
      allocate(fwz(1:nx,1:ny,1:nz))
#ifdef TWOCOMPONENT

      allocate(selphi(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff,2))
	  
      selphi=ZERO

	  allocate(modgrad(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff))
	  
	  allocate(arr_x(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff), arr_y(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff),arr_z(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff)) 
	  allocate(normx(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff),normy(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff),normz(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff))
#ifdef REPULSIVE_FLUX
	allocate(rep_mask(1:nx,1:ny,1:nz))
	allocate(pair_i(1:nx,1:ny,1:nz),pair_j(1:nx,1:ny,1:nz),pair_k(1:nx,1:ny,1:nz))
	allocate(Jx(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff),Jy(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff),Jz(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff))
#endif

#endif

#ifdef DENSRATIO    
      allocate(rhophi(0:nx+1,0:ny+1,0:nz+1))
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
         rhoprint(1:nxskip,1:nyskip,1:nzskip)=0.0
         velprint(1:3,1:nxskip,1:nyskip,1:nzskip)=0.0
#if defined(TWOCOMPONENT) && defined(WRITEPRESS)
         allocate(pressprint(1:nxskip,1:nyskip,1:nzskip))
         pressprint(1:nxskip,1:nyskip,1:nzskip)=0.0
#endif
      endif
      
      
     
   endsubroutine
endmodule
