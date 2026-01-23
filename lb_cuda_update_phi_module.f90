#include "defines.h"
#if !defined(_OPENACC)  
#error "To use this module the macros _OPENACC should be defined."
#endif

module lb_cuda_update_phi

   use vars
   use iso_c_binding
   use cudafor
   use mpi_template, only: coords,dostop,doerror,mydev,myrank,nprocs,nbuff,nbuffbvec
   use lb_cuda_vars

   implicit none
   

contains
   
   attributes(global) subroutine update_phifields_kernel(flop,nx,ny,nz,coords,isfluid &
#ifdef TWOCOMPONENT    
    ,visc2,rho_r,rho_b,invrho_r,invrho_b,sharp_c,beta,kapp,tau_diff,sigma &
#ifdef MONOD	
    ,mu_max,Ks &
#endif
#endif 
    ,ntothfields,ntotphifields,ntotauxfields,ntotlocauxfields,ntotforces &
       ,hfields_s,phifields_in,phifields_out,auxfields_s,locauxfields_s,forces_s)

      implicit none
      integer :: flop,nx,ny,nz
      integer, dimension(3) :: coords
      integer(kind=isf), dimension(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff) :: isfluid
     
 
#ifdef TWOCOMPONENT 
      real(kind=db) :: visc2,rho_r,rho_b,invrho_r,invrho_b,sharp_c,beta, &
       kapp,tau_diff,sigma,modgrad 
#ifdef MONOD
      real(kind=db) :: mu_max,Ks
#endif
#endif 
      integer :: ntothfields,ntotphifields,ntotauxfields,ntotlocauxfields,ntotforces
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields,nblocks_d) :: hfields_s
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields,nblocks_d) :: phifields_in,phifields_out
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nauxfields,nblocks_d) :: auxfields_s
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nlocauxfields,nblocks_d) :: locauxfields_s
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nforces,nblocks_d) :: forces_s
      
      integer :: i,j,k
      !integer :: gi,gj,gk
      integer :: myblock,ii,jj,kk
      real(kind=db) :: gradfix,gradfiy,gradfiz,phi_loc,phi_out,mytemp
      real(kind=db) :: loc_u,loc_v,loc_w,lap_phi_loc
#ifdef MONOD
	  real(kind=db) :: S_mono
#endif
      
      i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z
      
      if(abs(isfluid(i,j,k)) .ne. 1)return
      
      !gi=nx*coords(1)+i
      !gj=ny*coords(2)+j
      !gk=nz*coords(3)+k
      
      myblock=blockIdx%x+blockIdx%y*nxblock_d+blockIdx%z*nxyblock_d+1
      
      ii=threadIdx%x
      jj=threadIdx%y
      kk=threadIdx%z
	   
#ifdef TWOCOMPONENT
				  
	  mytemp= -sharp_c*locauxfields_s(ii,jj,kk,2,myblock)
	  			  
	  !reuse gradrhox,gradrhoy,gradrhoz as local velocity (reusing variables is saving register memory)
	  !reuse gradfix,gradfiy,gradfiz
      modgrad=auxfields_s(ii,jj,kk,4,myblock) !modgrad
	  gradfix=auxfields_s(ii,jj,kk,1,myblock)*modgrad !normx*modgrad
	  gradfiy=auxfields_s(ii,jj,kk,2,myblock)*modgrad !normy*modgrad
	  gradfiz=auxfields_s(ii,jj,kk,3,myblock)*modgrad !normz*modgrad
                  
      phi_loc=phifields_in(ii,jj,kk,1,myblock)
      lap_phi_loc=locauxfields_s(ii,jj,kk,1,myblock)
                  
      loc_u=hfields_s(ii,jj,kk,2,myblock) !velocity
      loc_v=hfields_s(ii,jj,kk,3,myblock)
      loc_w=hfields_s(ii,jj,kk,4,myblock)
                  
      phi_out = phi_loc &
        - loc_u*0.5_db*(gradfix) - loc_v*0.5_db*(gradfiy) &
        - loc_w*0.5_db*(gradfiz) + tau_diff*lap_phi_loc + mytemp 
#endif	

#ifdef MONOD
      S_mono = mu_max * phi_loc)/(Ks + phi_loc) * phi_loc * (1.0_db - phi_loc)
      phi_out=phi_out + S_mono
	  !phi_out = min(1.0_db, max(0.0_db, phi_out))		 
#endif
      phifields_out(ii,jj,kk,1,myblock)=phi_out
      
      return
      
   endsubroutine update_phifields_kernel

endmodule lb_cuda_update_phi
