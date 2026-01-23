#include "defines.h"
#if !defined(_OPENACC)  
#error "To use this module the macros _OPENACC should be defined."
#endif

module lb_cuda_auxfields

   use vars
   use iso_c_binding
   use cudafor
   use mpi_template, only: coords,dostop,doerror,mydev,myrank,nprocs,nbuff,nbuffbvec
   use lb_cuda_vars

   implicit none
   

contains

  attributes(global) subroutine compute_norm_interface_kernel(flop,nx,ny,nz,coords,isfluid, &
  rho_r,rho_b,ntotphifields,ntotauxfields,ntotlocauxfields,phifields_s,auxfields_s,locauxfields_s)

      implicit none
      integer :: flop,nx,ny,nz
      integer, dimension(3) :: coords
      integer(kind=isf), dimension(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff) :: isfluid
      real(kind=db) :: rho_r,rho_b
      integer :: ntotphifields,ntotauxfields,ntotlocauxfields
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields,nblocks_d) :: phifields_s
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nauxfields,nblocks_d) :: auxfields_s
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nlocauxfields,nblocks_d) :: locauxfields_s
      
      real(kind=db), shared :: myphi(0:TILE_DIMx+1,0:TILE_DIMy+1,0:TILE_DIMz+1)
      real(kind=db):: grad_fix,grad_fiy,grad_fiz,mod_grad
      
      integer :: i,j,k,gi,gj,gk,myblock,idblock
      integer :: ii,jj,kk
      integer :: li,lj,lk
      integer :: xblock,yblock,zblock
      
      li = threadIdx%x-1
      lj = threadIdx%y-1
      lk = threadIdx%z-1

      i = (blockIdx%x-1) * TILE_DIMx + li
      j = (blockIdx%y-1) * TILE_DIMy + lj
      k = (blockIdx%z-1) * TILE_DIMz + lk
      
      gi=nx*coords(1)+i
      gj=ny*coords(2)+j
      gk=nz*coords(3)+k
      
      xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
	  yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
	  zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
      
      myblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
      ii=i-xblock*TILE_DIMx+2*TILE_DIMx
      jj=j-yblock*TILE_DIMy+2*TILE_DIMy
      kk=k-zblock*TILE_DIMz+2*TILE_DIMz
      
      myphi(li,lj,lk)=phifields_s(ii,jj,kk,1,myblock)
      
      call syncthreads
      
      if(abs(isfluid(i,j,k)) .ne. 1)return
      
      !idblock is the index of the block of internal nodes without the surrounding halo
	  idblock=blockIdx%x+blockIdx%y*nxblock_d+blockIdx%z*nxyblock_d+1
	  !If my block index does not match the index of the internal-node block, it means my thread is on the outer halo and must exit.
	  if(myblock .ne. idblock)return

	  grad_fix=3.0_db*(p1*(myphi(li+1,lj,lk)-myphi(li-1,lj,lk)) + &
		 p2*( (myphi(li+1,lj+1,lk)-myphi(li-1,lj-1,lk))+ &
		 (myphi(li+1,lj-1,lk)-myphi(li-1,lj+1,lk))+ &
		 (myphi(li+1,lj,lk+1)-myphi(li-1,lj,lk-1))+ &
		 (myphi(li+1,lj,lk-1)-myphi(li-1,lj,lk+1)) )  + &
		 p3*((myphi(li+1,lj+1,lk+1)-myphi(li-1,lj-1,lk-1))+ &
		 (myphi(li+1,lj-1,lk-1)-myphi(li-1,lj+1,lk+1))+ &
		 (myphi(li+1,lj-1,lk+1)-myphi(li-1,lj+1,lk-1))+ &
		 (myphi(li+1,lj+1,lk-1)-myphi(li-1,lj-1,lk+1))))

	  grad_fiy=3.0_db*(p1*(myphi(li,lj+1,lk)-myphi(li,lj-1,lk)) + &
		 p2*((myphi(li+1,lj+1,lk)-myphi(li-1,lj-1,lk))+ &
		 (myphi(li-1,lj+1,lk)-myphi(li+1,lj-1,lk))+ &
		 (myphi(li,lj+1,lk+1)-myphi(li,lj-1,lk-1))+ &
		 (myphi(li,lj+1,lk-1)-myphi(li,lj-1,lk+1)) ) + &
		 p3*((myphi(li+1,lj+1,lk+1)-myphi(li-1,lj-1,lk-1))+ &
		 (myphi(li-1,lj+1,lk-1)-myphi(li+1,lj-1,lk+1))+ &
		 (myphi(li+1,lj+1,lk-1)-myphi(li-1,lj-1,lk+1))+ &
		 (myphi(li-1,lj+1,lk+1)-myphi(li+1,lj-1,lk-1))))

	  grad_fiz=3.0_db*(p1*(myphi(li,lj,lk+1)-myphi(li,lj,lk-1)) + &
		 p2*((myphi(li+1,lj,lk+1)-myphi(li-1,lj,lk-1))+ &
		 (myphi(li-1,lj,lk+1)-myphi(li+1,lj,lk-1))+ &
		 (myphi(li,lj+1,lk+1)-myphi(li,lj-1,lk-1))+ &
		 (myphi(li,lj-1,lk+1)-myphi(li,lj+1,lk-1)) ) + &
		 p3*((myphi(li+1,lj+1,lk+1)-myphi(li-1,lj-1,lk-1)) &
		 +(myphi(li-1,lj-1,lk+1)-myphi(li+1,lj+1,lk-1))+ &
		 (myphi(li+1,lj-1,lk+1)-myphi(li-1,lj+1,lk-1))+ &
		 (myphi(li-1,lj+1,lk+1)-myphi(li+1,lj-1,lk-1))))
      
	  mod_grad= sqrt(grad_fix**TWO + grad_fiy**TWO + grad_fiz**TWO)

	  auxfields_s(ii,jj,kk,1,myblock)= &
	   grad_fix/(mod_grad+1.0e-9_db)
	  auxfields_s(ii,jj,kk,2,myblock)= &
	   grad_fiy/(mod_grad+1.0e-9_db)
	  auxfields_s(ii,jj,kk,3,myblock)= &
	   grad_fiz/(mod_grad+1.0e-9_db)
	  
	  auxfields_s(ii,jj,kk,4,myblock)=mod_grad 

	  auxfields_s(ii,jj,kk,5,myblock)= &
	   myphi(li,lj,lk)*(1.0_db-myphi(li,lj,lk))*(grad_fix/(mod_grad+1.0e-9_db))
	  auxfields_s(ii,jj,kk,6,myblock)= &
	   myphi(li,lj,lk)*(1.0_db-myphi(li,lj,lk))*(grad_fiy/(mod_grad+1.0e-9_db))
	  auxfields_s(ii,jj,kk,7,myblock)= &
	   myphi(li,lj,lk)*(1.0_db-myphi(li,lj,lk))*(grad_fiz/(mod_grad+1.0e-9_db))
	   
      !lap_phi here
      locauxfields_s(ii,jj,kk,1,myblock)= &
                   (2.0_db*invcssq)*(myphi(li,lj,lk)*(p0-1.0_db) + &
                   ( p1*(myphi(li+1,lj,lk)+myphi(li-1,lj,lk) + &
                   myphi(li,lj+1,lk)+myphi(li,lj-1,lk) + &
                   myphi(li,lj,lk+1)+myphi(li,lj,lk-1)) + &
                   p2*( (myphi(li+1,lj+1,lk)+myphi(li-1,lj-1,lk))+ &
                   (myphi(li+1,lj-1,lk)+myphi(li-1,lj+1,lk))+ &
                   (myphi(li+1,lj,lk+1)+myphi(li-1,lj,lk-1))+ &
                   (myphi(li+1,lj,lk-1)+myphi(li-1,lj,lk+1)) + &
                   (myphi(li,lj+1,lk+1)+myphi(li,lj-1,lk-1))+ &
                   (myphi(li,lj+1,lk-1)+myphi(li,lj-1,lk+1)) )  + &
                   p3*((myphi(li+1,lj+1,lk+1)+myphi(li-1,lj-1,lk-1))+ &
                   (myphi(li+1,lj-1,lk-1)+myphi(li-1,lj+1,lk+1))+ &
                   (myphi(li+1,lj-1,lk+1)+myphi(li-1,lj+1,lk-1))+ &
                   (myphi(li+1,lj+1,lk-1)+myphi(li-1,lj-1,lk+1)))))
      
      return
      
   endsubroutine compute_norm_interface_kernel
   
   attributes(global) subroutine compute_div_theta_n_kernel(flop,nx,ny,nz,coords,isfluid, &
     rho_r,rho_b,ntotphifields,ntotauxfields,ntotlocauxfields,phifields_s,auxfields_s,locauxfields_s)

      implicit none
      integer :: flop,nx,ny,nz
      integer, dimension(3) :: coords
      integer(kind=isf), dimension(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff) :: isfluid
      real(kind=db) :: rho_r,rho_b
      integer :: ntotphifields,ntotauxfields,ntotlocauxfields
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields,nblocks_d) :: phifields_s
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nauxfields,nblocks_d) :: auxfields_s
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nlocauxfields,nblocks_d) :: locauxfields_s
      
      real(kind=db), shared :: myarrx(0:TILE_DIMx+1,0:TILE_DIMy+1,0:TILE_DIMz+1)
      real(kind=db), shared :: myarry(0:TILE_DIMx+1,0:TILE_DIMy+1,0:TILE_DIMz+1)
      real(kind=db), shared :: myarrz(0:TILE_DIMx+1,0:TILE_DIMy+1,0:TILE_DIMz+1)
      
      integer :: i,j,k,gi,gj,gk,myblock,idblock
      integer :: ii,jj,kk
      integer :: li,lj,lk
      integer :: xblock,yblock,zblock
      
      li = threadIdx%x-1
      lj = threadIdx%y-1
      lk = threadIdx%z-1
      
      i = (blockIdx%x-1) * TILE_DIMx + li
      j = (blockIdx%y-1) * TILE_DIMy + lj
      k = (blockIdx%z-1) * TILE_DIMz + lk
      
      gi=nx*coords(1)+i
      gj=ny*coords(2)+j
      gk=nz*coords(3)+k
      
      xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
	  yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
	  zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
      
      myblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
      ii=i-xblock*TILE_DIMx+2*TILE_DIMx
      jj=j-yblock*TILE_DIMy+2*TILE_DIMy
      kk=k-zblock*TILE_DIMz+2*TILE_DIMz
      
      myarrx(li,lj,lk)=auxfields_s(ii,jj,kk,5,myblock)
      myarry(li,lj,lk)=auxfields_s(ii,jj,kk,6,myblock)
      myarrz(li,lj,lk)=auxfields_s(ii,jj,kk,7,myblock)
      
      call syncthreads
      
      if(abs(isfluid(i,j,k)) .ne. 1)return
      
      !idblock is the index of the block of internal nodes without the surrounding halo
	  idblock=blockIdx%x+blockIdx%y*nxblock_d+blockIdx%z*nxyblock_d+1
	  !If my block index does not match the index of the internal-node block, it means my thread is on the outer halo and must exit.
	  if(myblock .ne. idblock)return


	   
      !div_thetan here
      locauxfields_s(ii,jj,kk,2,myblock)= &
       (( p1*(myarrx(li+1,lj,lk)-myarrx(li-1,lj,lk)) + &
       p2*((myarrx(li+1,lj+1,lk)-myarrx(li-1,lj-1,lk))+(myarrx(li+1,lj-1,lk)-myarrx(li-1,lj+1,lk))+ &
       (myarrx(li+1,lj,lk+1)-myarrx(li-1,lj,lk-1))+(myarrx(li+1,lj,lk-1)-myarrx(li-1,lj,lk+1)))+ &
       p3*((myarrx(li+1,lj+1,lk+1)-myarrx(li-1,lj-1,lk-1))+(myarrx(li+1,lj-1,lk-1)-myarrx(li-1,lj+1,lk+1))+ &
       (myarrx(li+1,lj-1,lk+1)-myarrx(li-1,lj+1,lk-1))+(myarrx(li+1,lj+1,lk-1)-myarrx(li-1,lj-1,lk+1))))+ &
       (p1*(myarry(li,lj+1,lk)-myarry(li,lj-1,lk)) + &
       p2*((myarry(li+1,lj+1,lk)-myarry(li-1,lj-1,lk))+(myarry(li-1,lj+1,lk)-myarry(li+1,lj-1,lk))+ &
       (myarry(li,lj+1,lk+1)-myarry(li,lj-1,lk-1))+(myarry(li,lj+1,lk-1)-myarry(li,lj-1,lk+1)))+ &
       p3*((myarry(li+1,lj+1,lk+1)-myarry(li-1,lj-1,lk-1))+(myarry(li-1,lj+1,lk-1)-myarry(li+1,lj-1,lk+1))+ &
       (myarry(li+1,lj+1,lk-1)-myarry(li-1,lj-1,lk+1))+(myarry(li-1,lj+1,lk+1)-myarry(li+1,lj-1,lk-1))))+ &
       (p1*(myarrz(li,lj,lk+1)-myarrz(li,lj,lk-1)) + &
       p2*((myarrz(li+1,lj,lk+1)-myarrz(li-1,lj,lk-1))+(myarrz(li-1,lj,lk+1)-myarrz(li+1,lj,lk-1))+ &
       (myarrz(li,lj+1,lk+1)-myarrz(li,lj-1,lk-1))+(myarrz(li,lj-1,lk+1)-myarrz(li,lj+1,lk-1)))+ &
       p3*((myarrz(li+1,lj+1,lk+1)-myarrz(li-1,lj-1,lk-1))+(myarrz(li-1,lj-1,lk+1)-myarrz(li+1,lj+1,lk-1))+ &
       (myarrz(li+1,lj-1,lk+1)-myarrz(li-1,lj+1,lk-1))+(myarrz(li-1,lj+1,lk+1)-myarrz(li+1,lj-1,lk-1)))))*invcssq
      
      return
      
   endsubroutine compute_div_theta_n_kernel

endmodule lb_cuda_auxfields
