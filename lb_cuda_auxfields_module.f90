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
#ifdef TWOCOMPONENT
  attributes(global) subroutine compute_norm_interface_kernel(step, &
   iprobe,jprobe,kprobe,flop,nx,ny,nz,coords,isfluid,rho_r,rho_b, &
   ntotphifields,ntotauxfields,ntotlocauxfields,phifields_s,auxfields_s,locauxfields_s)

      implicit none
      integer :: step,iprobe,jprobe,kprobe,flop,nx,ny,nz
      integer, dimension(3) :: coords
      integer(kind=isf), dimension(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff) :: isfluid
      real(kind=db) :: rho_r,rho_b
      integer :: ntotphifields,ntotauxfields,ntotlocauxfields
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields,nblocks_d) :: phifields_s
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nauxfields,nblocks_d) :: auxfields_s
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nlocauxfields,nblocks_d) :: locauxfields_s
      
      real(kind=db), shared :: myphi(0:TILE_DIMx+1,0:TILE_DIMy+1,0:TILE_DIMz+1)
      real(kind=db):: grad_fix,grad_fiy,grad_fiz,mod_grad
      
      integer :: i,j,k,gi,gj,gk,myblock,intblock
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
      
      !intblock is the index of the block of internal nodes without the surrounding halo
	  intblock=blockIdx%x+blockIdx%y*nxblock_d+blockIdx%z*nxyblock_d+1
	  !If my block index does not match the index of the internal-node block, it means my thread is on the outer halo and must exit.
	  if(myblock .ne. intblock)return

	  grad_fix=3.0_db*(p1d3q27*(myphi(li+1,lj,lk)-myphi(li-1,lj,lk)) + &
		 p2d3q27*( (myphi(li+1,lj+1,lk)-myphi(li-1,lj-1,lk))+ &
		 (myphi(li+1,lj-1,lk)-myphi(li-1,lj+1,lk))+ &
		 (myphi(li+1,lj,lk+1)-myphi(li-1,lj,lk-1))+ &
		 (myphi(li+1,lj,lk-1)-myphi(li-1,lj,lk+1)) )  + &
		 p3d3q27*((myphi(li+1,lj+1,lk+1)-myphi(li-1,lj-1,lk-1))+ &
		 (myphi(li+1,lj-1,lk-1)-myphi(li-1,lj+1,lk+1))+ &
		 (myphi(li+1,lj-1,lk+1)-myphi(li-1,lj+1,lk-1))+ &
		 (myphi(li+1,lj+1,lk-1)-myphi(li-1,lj-1,lk+1))))

	  grad_fiy=3.0_db*(p1d3q27*(myphi(li,lj+1,lk)-myphi(li,lj-1,lk)) + &
		 p2d3q27*((myphi(li+1,lj+1,lk)-myphi(li-1,lj-1,lk))+ &
		 (myphi(li-1,lj+1,lk)-myphi(li+1,lj-1,lk))+ &
		 (myphi(li,lj+1,lk+1)-myphi(li,lj-1,lk-1))+ &
		 (myphi(li,lj+1,lk-1)-myphi(li,lj-1,lk+1)) ) + &
		 p3d3q27*((myphi(li+1,lj+1,lk+1)-myphi(li-1,lj-1,lk-1))+ &
		 (myphi(li-1,lj+1,lk-1)-myphi(li+1,lj-1,lk+1))+ &
		 (myphi(li+1,lj+1,lk-1)-myphi(li-1,lj-1,lk+1))+ &
		 (myphi(li-1,lj+1,lk+1)-myphi(li+1,lj-1,lk-1))))

	  grad_fiz=3.0_db*(p1d3q27*(myphi(li,lj,lk+1)-myphi(li,lj,lk-1)) + &
		 p2d3q27*((myphi(li+1,lj,lk+1)-myphi(li-1,lj,lk-1))+ &
		 (myphi(li-1,lj,lk+1)-myphi(li+1,lj,lk-1))+ &
		 (myphi(li,lj+1,lk+1)-myphi(li,lj-1,lk-1))+ &
		 (myphi(li,lj-1,lk+1)-myphi(li,lj+1,lk-1)) ) + &
		 p3d3q27*((myphi(li+1,lj+1,lk+1)-myphi(li-1,lj-1,lk-1)) &
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
                   (2.0_db*invcssq)*(myphi(li,lj,lk)*(p0d3q27-1.0_db) + &
                   ( p1d3q27*(myphi(li+1,lj,lk)+myphi(li-1,lj,lk) + &
                   myphi(li,lj+1,lk)+myphi(li,lj-1,lk) + &
                   myphi(li,lj,lk+1)+myphi(li,lj,lk-1)) + &
                   p2d3q27*( (myphi(li+1,lj+1,lk)+myphi(li-1,lj-1,lk))+ &
                   (myphi(li+1,lj-1,lk)+myphi(li-1,lj+1,lk))+ &
                   (myphi(li+1,lj,lk+1)+myphi(li-1,lj,lk-1))+ &
                   (myphi(li+1,lj,lk-1)+myphi(li-1,lj,lk+1)) + &
                   (myphi(li,lj+1,lk+1)+myphi(li,lj-1,lk-1))+ &
                   (myphi(li,lj+1,lk-1)+myphi(li,lj-1,lk+1)) )  + &
                   p3d3q27*((myphi(li+1,lj+1,lk+1)+myphi(li-1,lj-1,lk-1))+ &
                   (myphi(li+1,lj-1,lk-1)+myphi(li-1,lj+1,lk+1))+ &
                   (myphi(li+1,lj-1,lk+1)+myphi(li-1,lj+1,lk-1))+ &
                   (myphi(li+1,lj+1,lk-1)+myphi(li-1,lj-1,lk+1)))))
      
      return
      
   endsubroutine compute_norm_interface_kernel
   
     attributes(global) subroutine compute_norm_interface_kernel_int(step, &
      iprobe,jprobe,kprobe,flop,nx,ny,nz,coords,isfluid,rho_r,rho_b, &
      ntotphifields,ntotauxfields,ntotlocauxfields,phifields_s,auxfields_s,locauxfields_s)

      implicit none
      integer :: step,iprobe,jprobe,kprobe,flop,nx,ny,nz
      integer, dimension(3) :: coords
      integer(kind=isf), dimension(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff) :: isfluid
      real(kind=db) :: rho_r,rho_b
      integer :: ntotphifields,ntotauxfields,ntotlocauxfields
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields,nblocks_d) :: phifields_s
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nauxfields,nblocks_d) :: auxfields_s
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nlocauxfields,nblocks_d) :: locauxfields_s
      
      real(kind=db), shared :: myphi(0:TILE_DIMx+1,0:TILE_DIMy+1,0:TILE_DIMz+1)
      real(kind=db):: grad_fix,grad_fiy,grad_fiz,mod_grad
      
      integer :: i,j,k,myblock,intblock
      integer :: ii,jj,kk
      integer :: li,lj,lk
      integer :: xblock,yblock,zblock
!      integer :: gi,gj,gk
      
      li = threadIdx%x-1
      lj = threadIdx%y-1
      lk = threadIdx%z-1
      
      i = (blockIdx%x) * TILE_DIMx + li
      j = (blockIdx%y) * TILE_DIMy + lj
      k = (blockIdx%z) * TILE_DIMz + lk
      
!      gi=nx*coords(1)+i
!      gj=ny*coords(2)+j
!      gk=nz*coords(3)+k
      
      xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
	  yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
	  zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
      
      myblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
      ii=i-xblock*TILE_DIMx+2*TILE_DIMx
      jj=j-yblock*TILE_DIMy+2*TILE_DIMy
      kk=k-zblock*TILE_DIMz+2*TILE_DIMz

      !i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      !j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      !k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z
      
!      gi=nx*coords(1)+i
!      gj=ny*coords(2)+j
!      gk=nz*coords(3)+k
      
      intblock=(blockIdx%x+1)+(blockIdx%y+1)*nxblock_d+(blockIdx%z+1)*nxyblock_d+1 !internal-node block
      
      myphi(li,lj,lk)=phifields_s(ii,jj,kk,1,myblock)
      
      call syncthreads
      
      if(abs(isfluid(i,j,k)) .ne. 1)return
      
	  !If my block index does not match the index of the internal-node block, it means my thread is on the outer halo and must exit.
	  if(myblock .ne. intblock)return

	  grad_fix=3.0_db*(p1d3q27*(myphi(li+1,lj,lk)-myphi(li-1,lj,lk)) + &
		 p2d3q27*( (myphi(li+1,lj+1,lk)-myphi(li-1,lj-1,lk))+ &
		 (myphi(li+1,lj-1,lk)-myphi(li-1,lj+1,lk))+ &
		 (myphi(li+1,lj,lk+1)-myphi(li-1,lj,lk-1))+ &
		 (myphi(li+1,lj,lk-1)-myphi(li-1,lj,lk+1)) )  + &
		 p3d3q27*((myphi(li+1,lj+1,lk+1)-myphi(li-1,lj-1,lk-1))+ &
		 (myphi(li+1,lj-1,lk-1)-myphi(li-1,lj+1,lk+1))+ &
		 (myphi(li+1,lj-1,lk+1)-myphi(li-1,lj+1,lk-1))+ &
		 (myphi(li+1,lj+1,lk-1)-myphi(li-1,lj-1,lk+1))))

	  grad_fiy=3.0_db*(p1d3q27*(myphi(li,lj+1,lk)-myphi(li,lj-1,lk)) + &
		 p2d3q27*((myphi(li+1,lj+1,lk)-myphi(li-1,lj-1,lk))+ &
		 (myphi(li-1,lj+1,lk)-myphi(li+1,lj-1,lk))+ &
		 (myphi(li,lj+1,lk+1)-myphi(li,lj-1,lk-1))+ &
		 (myphi(li,lj+1,lk-1)-myphi(li,lj-1,lk+1)) ) + &
		 p3d3q27*((myphi(li+1,lj+1,lk+1)-myphi(li-1,lj-1,lk-1))+ &
		 (myphi(li-1,lj+1,lk-1)-myphi(li+1,lj-1,lk+1))+ &
		 (myphi(li+1,lj+1,lk-1)-myphi(li-1,lj-1,lk+1))+ &
		 (myphi(li-1,lj+1,lk+1)-myphi(li+1,lj-1,lk-1))))

	  grad_fiz=3.0_db*(p1d3q27*(myphi(li,lj,lk+1)-myphi(li,lj,lk-1)) + &
		 p2d3q27*((myphi(li+1,lj,lk+1)-myphi(li-1,lj,lk-1))+ &
		 (myphi(li-1,lj,lk+1)-myphi(li+1,lj,lk-1))+ &
		 (myphi(li,lj+1,lk+1)-myphi(li,lj-1,lk-1))+ &
		 (myphi(li,lj-1,lk+1)-myphi(li,lj+1,lk-1)) ) + &
		 p3d3q27*((myphi(li+1,lj+1,lk+1)-myphi(li-1,lj-1,lk-1)) &
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
                   (2.0_db*invcssq)*(myphi(li,lj,lk)*(p0d3q27-1.0_db) + &
                   ( p1d3q27*(myphi(li+1,lj,lk)+myphi(li-1,lj,lk) + &
                   myphi(li,lj+1,lk)+myphi(li,lj-1,lk) + &
                   myphi(li,lj,lk+1)+myphi(li,lj,lk-1)) + &
                   p2d3q27*( (myphi(li+1,lj+1,lk)+myphi(li-1,lj-1,lk))+ &
                   (myphi(li+1,lj-1,lk)+myphi(li-1,lj+1,lk))+ &
                   (myphi(li+1,lj,lk+1)+myphi(li-1,lj,lk-1))+ &
                   (myphi(li+1,lj,lk-1)+myphi(li-1,lj,lk+1)) + &
                   (myphi(li,lj+1,lk+1)+myphi(li,lj-1,lk-1))+ &
                   (myphi(li,lj+1,lk-1)+myphi(li,lj-1,lk+1)) )  + &
                   p3d3q27*((myphi(li+1,lj+1,lk+1)+myphi(li-1,lj-1,lk-1))+ &
                   (myphi(li+1,lj-1,lk-1)+myphi(li-1,lj+1,lk+1))+ &
                   (myphi(li+1,lj-1,lk+1)+myphi(li-1,lj+1,lk-1))+ &
                   (myphi(li+1,lj+1,lk-1)+myphi(li-1,lj-1,lk+1)))))
      
      return
      
   endsubroutine compute_norm_interface_kernel_int
   
   attributes(global) subroutine compute_norm_interface_kernel_xplus(step, &
      iprobe,jprobe,kprobe,flop,nx,ny,nz,coords,isfluid,rho_r,rho_b, &
      ntotphifields,ntotauxfields,ntotlocauxfields,phifields_s,auxfields_s,locauxfields_s)

      implicit none
      integer :: step,iprobe,jprobe,kprobe,flop,nx,ny,nz
      integer, dimension(3) :: coords
      integer(kind=isf), dimension(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff) :: isfluid
      real(kind=db) :: rho_r,rho_b
      integer :: ntotphifields,ntotauxfields,ntotlocauxfields
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields,nblocks_d) :: phifields_s
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nauxfields,nblocks_d) :: auxfields_s
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nlocauxfields,nblocks_d) :: locauxfields_s
      
      real(kind=db), shared :: myphi(0:TILE_DIMx+1,0:TILE_DIMy+1,0:TILE_DIMz+1)
      real(kind=db):: grad_fix,grad_fiy,grad_fiz,mod_grad
      
      integer :: i,j,k,myblock,intblock
      integer :: ii,jj,kk
      integer :: li,lj,lk
      integer :: xblock,yblock,zblock
!      integer :: gi,gj,gk
      
      
      li = threadIdx%x-1
      lj = threadIdx%y-1
      lk = threadIdx%z-1
      
      i = ((nxblock_d-2)-1) * TILE_DIMx + li
      j = (blockIdx%y) * TILE_DIMy + lj 
      k = (blockIdx%z) * TILE_DIMz + lk 

!      gj=ny*coords(2)+j
!      gk=nz*coords(3)+k
      
	  yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
	  zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
      xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
      
      
	  !gi=nx*coords(1)+i
	  
      
      myblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
      ii=i-xblock*TILE_DIMx+2*TILE_DIMx
      jj=j-yblock*TILE_DIMy+2*TILE_DIMy
      kk=k-zblock*TILE_DIMz+2*TILE_DIMz

      !i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      !j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      !k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z
      
!      gi=nx*coords(1)+i
!      gj=ny*coords(2)+j
!      gk=nz*coords(3)+k
      
      !if(gi==iprobe)write(*,*)'cazzo',xblock,blockIdx%y,blockIdx%z
      
      
      intblock=(nxblock_d-2)+(blockIdx%y+1)*nxblock_d+(blockIdx%z+1)*nxyblock_d+1 !internal-node block
      
      myphi(li,lj,lk)=phifields_s(ii,jj,kk,1,myblock)
      
      call syncthreads
      
      if(abs(isfluid(i,j,k)) .ne. 1)return
      
	  !If my block index does not match the index of the internal-node block, it means my thread is on the outer halo and must exit.
	  if(myblock .ne. intblock)return

	  grad_fix=3.0_db*(p1d3q27*(myphi(li+1,lj,lk)-myphi(li-1,lj,lk)) + &
		 p2d3q27*( (myphi(li+1,lj+1,lk)-myphi(li-1,lj-1,lk))+ &
		 (myphi(li+1,lj-1,lk)-myphi(li-1,lj+1,lk))+ &
		 (myphi(li+1,lj,lk+1)-myphi(li-1,lj,lk-1))+ &
		 (myphi(li+1,lj,lk-1)-myphi(li-1,lj,lk+1)) )  + &
		 p3d3q27*((myphi(li+1,lj+1,lk+1)-myphi(li-1,lj-1,lk-1))+ &
		 (myphi(li+1,lj-1,lk-1)-myphi(li-1,lj+1,lk+1))+ &
		 (myphi(li+1,lj-1,lk+1)-myphi(li-1,lj+1,lk-1))+ &
		 (myphi(li+1,lj+1,lk-1)-myphi(li-1,lj-1,lk+1))))

	  grad_fiy=3.0_db*(p1d3q27*(myphi(li,lj+1,lk)-myphi(li,lj-1,lk)) + &
		 p2d3q27*((myphi(li+1,lj+1,lk)-myphi(li-1,lj-1,lk))+ &
		 (myphi(li-1,lj+1,lk)-myphi(li+1,lj-1,lk))+ &
		 (myphi(li,lj+1,lk+1)-myphi(li,lj-1,lk-1))+ &
		 (myphi(li,lj+1,lk-1)-myphi(li,lj-1,lk+1)) ) + &
		 p3d3q27*((myphi(li+1,lj+1,lk+1)-myphi(li-1,lj-1,lk-1))+ &
		 (myphi(li-1,lj+1,lk-1)-myphi(li+1,lj-1,lk+1))+ &
		 (myphi(li+1,lj+1,lk-1)-myphi(li-1,lj-1,lk+1))+ &
		 (myphi(li-1,lj+1,lk+1)-myphi(li+1,lj-1,lk-1))))

	  grad_fiz=3.0_db*(p1d3q27*(myphi(li,lj,lk+1)-myphi(li,lj,lk-1)) + &
		 p2d3q27*((myphi(li+1,lj,lk+1)-myphi(li-1,lj,lk-1))+ &
		 (myphi(li-1,lj,lk+1)-myphi(li+1,lj,lk-1))+ &
		 (myphi(li,lj+1,lk+1)-myphi(li,lj-1,lk-1))+ &
		 (myphi(li,lj-1,lk+1)-myphi(li,lj+1,lk-1)) ) + &
		 p3d3q27*((myphi(li+1,lj+1,lk+1)-myphi(li-1,lj-1,lk-1)) &
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
                   (2.0_db*invcssq)*(myphi(li,lj,lk)*(p0d3q27-1.0_db) + &
                   ( p1d3q27*(myphi(li+1,lj,lk)+myphi(li-1,lj,lk) + &
                   myphi(li,lj+1,lk)+myphi(li,lj-1,lk) + &
                   myphi(li,lj,lk+1)+myphi(li,lj,lk-1)) + &
                   p2d3q27*( (myphi(li+1,lj+1,lk)+myphi(li-1,lj-1,lk))+ &
                   (myphi(li+1,lj-1,lk)+myphi(li-1,lj+1,lk))+ &
                   (myphi(li+1,lj,lk+1)+myphi(li-1,lj,lk-1))+ &
                   (myphi(li+1,lj,lk-1)+myphi(li-1,lj,lk+1)) + &
                   (myphi(li,lj+1,lk+1)+myphi(li,lj-1,lk-1))+ &
                   (myphi(li,lj+1,lk-1)+myphi(li,lj-1,lk+1)) )  + &
                   p3d3q27*((myphi(li+1,lj+1,lk+1)+myphi(li-1,lj-1,lk-1))+ &
                   (myphi(li+1,lj-1,lk-1)+myphi(li-1,lj+1,lk+1))+ &
                   (myphi(li+1,lj-1,lk+1)+myphi(li-1,lj+1,lk-1))+ &
                   (myphi(li+1,lj+1,lk-1)+myphi(li-1,lj-1,lk+1)))))
      
      return
      
   endsubroutine compute_norm_interface_kernel_xplus
   
   attributes(global) subroutine compute_norm_interface_kernel_xminus(step, &
      iprobe,jprobe,kprobe,flop,nx,ny,nz,coords,isfluid,rho_r,rho_b, &
      ntotphifields,ntotauxfields,ntotlocauxfields,phifields_s,auxfields_s,locauxfields_s)

      implicit none
      integer :: step,iprobe,jprobe,kprobe,flop,nx,ny,nz
      integer, dimension(3) :: coords
      integer(kind=isf), dimension(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff) :: isfluid
      real(kind=db) :: rho_r,rho_b
      integer :: ntotphifields,ntotauxfields,ntotlocauxfields
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields,nblocks_d) :: phifields_s
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nauxfields,nblocks_d) :: auxfields_s
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nlocauxfields,nblocks_d) :: locauxfields_s
      
      real(kind=db), shared :: myphi(0:TILE_DIMx+1,0:TILE_DIMy+1,0:TILE_DIMz+1)
      real(kind=db):: grad_fix,grad_fiy,grad_fiz,mod_grad
      
      integer :: i,j,k,myblock,intblock
      integer :: ii,jj,kk
      integer :: li,lj,lk
      integer :: xblock,yblock,zblock
!      integer :: gi,gj,gk

      
      li = threadIdx%x-1
      lj = threadIdx%y-1
      lk = threadIdx%z-1
      
      i = li
      j = (blockIdx%y) * TILE_DIMy + lj 
      k = (blockIdx%z) * TILE_DIMz + lk 

!      gj=ny*coords(2)+j
!      gk=nz*coords(3)+k
      
	  yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
	  zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
      xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
      
      
	  !gi=nx*coords(1)+i
	  
      
      myblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
      ii=i-xblock*TILE_DIMx+2*TILE_DIMx
      jj=j-yblock*TILE_DIMy+2*TILE_DIMy
      kk=k-zblock*TILE_DIMz+2*TILE_DIMz

      !i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      !j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      !k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z
      
!      gi=nx*coords(1)+i
!      gj=ny*coords(2)+j
!      gk=nz*coords(3)+k
      
      !if(gi==iprobe)write(*,*)'cazzo',xblock,blockIdx%y,blockIdx%z
      
      
      intblock=1+(blockIdx%y+1)*nxblock_d+(blockIdx%z+1)*nxyblock_d+1 !internal-node block
      
      myphi(li,lj,lk)=phifields_s(ii,jj,kk,1,myblock)
      
      call syncthreads
      
      if(abs(isfluid(i,j,k)) .ne. 1)return
      
	  !If my block index does not match the index of the internal-node block, it means my thread is on the outer halo and must exit.
	  if(myblock .ne. intblock)return

	  grad_fix=3.0_db*(p1d3q27*(myphi(li+1,lj,lk)-myphi(li-1,lj,lk)) + &
		 p2d3q27*( (myphi(li+1,lj+1,lk)-myphi(li-1,lj-1,lk))+ &
		 (myphi(li+1,lj-1,lk)-myphi(li-1,lj+1,lk))+ &
		 (myphi(li+1,lj,lk+1)-myphi(li-1,lj,lk-1))+ &
		 (myphi(li+1,lj,lk-1)-myphi(li-1,lj,lk+1)) )  + &
		 p3d3q27*((myphi(li+1,lj+1,lk+1)-myphi(li-1,lj-1,lk-1))+ &
		 (myphi(li+1,lj-1,lk-1)-myphi(li-1,lj+1,lk+1))+ &
		 (myphi(li+1,lj-1,lk+1)-myphi(li-1,lj+1,lk-1))+ &
		 (myphi(li+1,lj+1,lk-1)-myphi(li-1,lj-1,lk+1))))

	  grad_fiy=3.0_db*(p1d3q27*(myphi(li,lj+1,lk)-myphi(li,lj-1,lk)) + &
		 p2d3q27*((myphi(li+1,lj+1,lk)-myphi(li-1,lj-1,lk))+ &
		 (myphi(li-1,lj+1,lk)-myphi(li+1,lj-1,lk))+ &
		 (myphi(li,lj+1,lk+1)-myphi(li,lj-1,lk-1))+ &
		 (myphi(li,lj+1,lk-1)-myphi(li,lj-1,lk+1)) ) + &
		 p3d3q27*((myphi(li+1,lj+1,lk+1)-myphi(li-1,lj-1,lk-1))+ &
		 (myphi(li-1,lj+1,lk-1)-myphi(li+1,lj-1,lk+1))+ &
		 (myphi(li+1,lj+1,lk-1)-myphi(li-1,lj-1,lk+1))+ &
		 (myphi(li-1,lj+1,lk+1)-myphi(li+1,lj-1,lk-1))))

	  grad_fiz=3.0_db*(p1d3q27*(myphi(li,lj,lk+1)-myphi(li,lj,lk-1)) + &
		 p2d3q27*((myphi(li+1,lj,lk+1)-myphi(li-1,lj,lk-1))+ &
		 (myphi(li-1,lj,lk+1)-myphi(li+1,lj,lk-1))+ &
		 (myphi(li,lj+1,lk+1)-myphi(li,lj-1,lk-1))+ &
		 (myphi(li,lj-1,lk+1)-myphi(li,lj+1,lk-1)) ) + &
		 p3d3q27*((myphi(li+1,lj+1,lk+1)-myphi(li-1,lj-1,lk-1)) &
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
                   (2.0_db*invcssq)*(myphi(li,lj,lk)*(p0d3q27-1.0_db) + &
                   ( p1d3q27*(myphi(li+1,lj,lk)+myphi(li-1,lj,lk) + &
                   myphi(li,lj+1,lk)+myphi(li,lj-1,lk) + &
                   myphi(li,lj,lk+1)+myphi(li,lj,lk-1)) + &
                   p2d3q27*( (myphi(li+1,lj+1,lk)+myphi(li-1,lj-1,lk))+ &
                   (myphi(li+1,lj-1,lk)+myphi(li-1,lj+1,lk))+ &
                   (myphi(li+1,lj,lk+1)+myphi(li-1,lj,lk-1))+ &
                   (myphi(li+1,lj,lk-1)+myphi(li-1,lj,lk+1)) + &
                   (myphi(li,lj+1,lk+1)+myphi(li,lj-1,lk-1))+ &
                   (myphi(li,lj+1,lk-1)+myphi(li,lj-1,lk+1)) )  + &
                   p3d3q27*((myphi(li+1,lj+1,lk+1)+myphi(li-1,lj-1,lk-1))+ &
                   (myphi(li+1,lj-1,lk-1)+myphi(li-1,lj+1,lk+1))+ &
                   (myphi(li+1,lj-1,lk+1)+myphi(li-1,lj+1,lk-1))+ &
                   (myphi(li+1,lj+1,lk-1)+myphi(li-1,lj-1,lk+1)))))
      
      return
      
   endsubroutine compute_norm_interface_kernel_xminus
   
   attributes(global) subroutine compute_norm_interface_kernel_yplus(step, &
      iprobe,jprobe,kprobe,flop,nx,ny,nz,coords,isfluid,rho_r,rho_b, &
      ntotphifields,ntotauxfields,ntotlocauxfields,phifields_s,auxfields_s,locauxfields_s)

      implicit none
      integer :: step,iprobe,jprobe,kprobe,flop,nx,ny,nz
      integer, dimension(3) :: coords
      integer(kind=isf), dimension(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff) :: isfluid
      real(kind=db) :: rho_r,rho_b
      integer :: ntotphifields,ntotauxfields,ntotlocauxfields
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields,nblocks_d) :: phifields_s
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nauxfields,nblocks_d) :: auxfields_s
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nlocauxfields,nblocks_d) :: locauxfields_s
      
      real(kind=db), shared :: myphi(0:TILE_DIMx+1,0:TILE_DIMy+1,0:TILE_DIMz+1)
      real(kind=db):: grad_fix,grad_fiy,grad_fiz,mod_grad
      
      integer :: i,j,k,myblock,intblock
      integer :: ii,jj,kk
      integer :: li,lj,lk
      integer :: xblock,yblock,zblock
!      integer :: gi,gj,gk

      
      li = threadIdx%x-1
      lj = threadIdx%y-1
      lk = threadIdx%z-1

      i = (blockIdx%x-1) * TILE_DIMx + li 
      j = ((nyblock_d-2)-1) * TILE_DIMy + lj
      k = (blockIdx%z) * TILE_DIMz + lk
      
      !gi=nx*coords(1)+i
      !gk=nz*coords(3)+k
      
      xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
      yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
	  zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
      
      
      
      myblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
      ii=i-xblock*TILE_DIMx+2*TILE_DIMx
      jj=j-yblock*TILE_DIMy+2*TILE_DIMy
      kk=k-zblock*TILE_DIMz+2*TILE_DIMz

      !i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      !j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      !k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z
      
!      gi=nx*coords(1)+i
!      gj=ny*coords(2)+j
!      gk=nz*coords(3)+k
      
      intblock=blockIdx%x+(nyblock_d-2)*nxblock_d+(blockIdx%z+1)*nxyblock_d+1 !internal-node block
      
      myphi(li,lj,lk)=phifields_s(ii,jj,kk,1,myblock)
      
      call syncthreads
      
      if(abs(isfluid(i,j,k)) .ne. 1)return
      
	  !If my block index does not match the index of the internal-node block, it means my thread is on the outer halo and must exit.
	  if(myblock .ne. intblock)return

	  grad_fix=3.0_db*(p1d3q27*(myphi(li+1,lj,lk)-myphi(li-1,lj,lk)) + &
		 p2d3q27*( (myphi(li+1,lj+1,lk)-myphi(li-1,lj-1,lk))+ &
		 (myphi(li+1,lj-1,lk)-myphi(li-1,lj+1,lk))+ &
		 (myphi(li+1,lj,lk+1)-myphi(li-1,lj,lk-1))+ &
		 (myphi(li+1,lj,lk-1)-myphi(li-1,lj,lk+1)) )  + &
		 p3d3q27*((myphi(li+1,lj+1,lk+1)-myphi(li-1,lj-1,lk-1))+ &
		 (myphi(li+1,lj-1,lk-1)-myphi(li-1,lj+1,lk+1))+ &
		 (myphi(li+1,lj-1,lk+1)-myphi(li-1,lj+1,lk-1))+ &
		 (myphi(li+1,lj+1,lk-1)-myphi(li-1,lj-1,lk+1))))

	  grad_fiy=3.0_db*(p1d3q27*(myphi(li,lj+1,lk)-myphi(li,lj-1,lk)) + &
		 p2d3q27*((myphi(li+1,lj+1,lk)-myphi(li-1,lj-1,lk))+ &
		 (myphi(li-1,lj+1,lk)-myphi(li+1,lj-1,lk))+ &
		 (myphi(li,lj+1,lk+1)-myphi(li,lj-1,lk-1))+ &
		 (myphi(li,lj+1,lk-1)-myphi(li,lj-1,lk+1)) ) + &
		 p3d3q27*((myphi(li+1,lj+1,lk+1)-myphi(li-1,lj-1,lk-1))+ &
		 (myphi(li-1,lj+1,lk-1)-myphi(li+1,lj-1,lk+1))+ &
		 (myphi(li+1,lj+1,lk-1)-myphi(li-1,lj-1,lk+1))+ &
		 (myphi(li-1,lj+1,lk+1)-myphi(li+1,lj-1,lk-1))))

	  grad_fiz=3.0_db*(p1d3q27*(myphi(li,lj,lk+1)-myphi(li,lj,lk-1)) + &
		 p2d3q27*((myphi(li+1,lj,lk+1)-myphi(li-1,lj,lk-1))+ &
		 (myphi(li-1,lj,lk+1)-myphi(li+1,lj,lk-1))+ &
		 (myphi(li,lj+1,lk+1)-myphi(li,lj-1,lk-1))+ &
		 (myphi(li,lj-1,lk+1)-myphi(li,lj+1,lk-1)) ) + &
		 p3d3q27*((myphi(li+1,lj+1,lk+1)-myphi(li-1,lj-1,lk-1)) &
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
                   (2.0_db*invcssq)*(myphi(li,lj,lk)*(p0d3q27-1.0_db) + &
                   ( p1d3q27*(myphi(li+1,lj,lk)+myphi(li-1,lj,lk) + &
                   myphi(li,lj+1,lk)+myphi(li,lj-1,lk) + &
                   myphi(li,lj,lk+1)+myphi(li,lj,lk-1)) + &
                   p2d3q27*( (myphi(li+1,lj+1,lk)+myphi(li-1,lj-1,lk))+ &
                   (myphi(li+1,lj-1,lk)+myphi(li-1,lj+1,lk))+ &
                   (myphi(li+1,lj,lk+1)+myphi(li-1,lj,lk-1))+ &
                   (myphi(li+1,lj,lk-1)+myphi(li-1,lj,lk+1)) + &
                   (myphi(li,lj+1,lk+1)+myphi(li,lj-1,lk-1))+ &
                   (myphi(li,lj+1,lk-1)+myphi(li,lj-1,lk+1)) )  + &
                   p3d3q27*((myphi(li+1,lj+1,lk+1)+myphi(li-1,lj-1,lk-1))+ &
                   (myphi(li+1,lj-1,lk-1)+myphi(li-1,lj+1,lk+1))+ &
                   (myphi(li+1,lj-1,lk+1)+myphi(li-1,lj+1,lk-1))+ &
                   (myphi(li+1,lj+1,lk-1)+myphi(li-1,lj-1,lk+1)))))
      
      return
      
   endsubroutine compute_norm_interface_kernel_yplus
   
   attributes(global) subroutine compute_norm_interface_kernel_yminus(step, &
      iprobe,jprobe,kprobe,flop,nx,ny,nz,coords,isfluid,rho_r,rho_b, &
      ntotphifields,ntotauxfields,ntotlocauxfields,phifields_s,auxfields_s,locauxfields_s)

      implicit none
      integer :: step,iprobe,jprobe,kprobe,flop,nx,ny,nz
      integer, dimension(3) :: coords
      integer(kind=isf), dimension(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff) :: isfluid
      real(kind=db) :: rho_r,rho_b
      integer :: ntotphifields,ntotauxfields,ntotlocauxfields
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields,nblocks_d) :: phifields_s
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nauxfields,nblocks_d) :: auxfields_s
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nlocauxfields,nblocks_d) :: locauxfields_s
      
      real(kind=db), shared :: myphi(0:TILE_DIMx+1,0:TILE_DIMy+1,0:TILE_DIMz+1)
      real(kind=db):: grad_fix,grad_fiy,grad_fiz,mod_grad
      
      integer :: i,j,k,myblock,intblock
      integer :: ii,jj,kk
      integer :: li,lj,lk
      integer :: xblock,yblock,zblock
!      integer :: gi,gj,gk

      
      li = threadIdx%x-1
      lj = threadIdx%y-1
      lk = threadIdx%z-1

      i = (blockIdx%x-1) * TILE_DIMx + li 
      j = lj
      k = (blockIdx%z) * TILE_DIMz + lk
      
      !gi=nx*coords(1)+i
      !gk=nz*coords(3)+k
      
      xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
      yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
	  zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
      
      
      
      myblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
      ii=i-xblock*TILE_DIMx+2*TILE_DIMx
      jj=j-yblock*TILE_DIMy+2*TILE_DIMy
      kk=k-zblock*TILE_DIMz+2*TILE_DIMz

      !i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      !j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      !k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z
      
!      gi=nx*coords(1)+i
!      gj=ny*coords(2)+j
!      gk=nz*coords(3)+k
      
      intblock=blockIdx%x+1*nxblock_d+(blockIdx%z+1)*nxyblock_d+1 !internal-node block
      
      myphi(li,lj,lk)=phifields_s(ii,jj,kk,1,myblock)
      
      call syncthreads
      
      if(abs(isfluid(i,j,k)) .ne. 1)return
      
	  !If my block index does not match the index of the internal-node block, it means my thread is on the outer halo and must exit.
	  if(myblock .ne. intblock)return

	  grad_fix=3.0_db*(p1d3q27*(myphi(li+1,lj,lk)-myphi(li-1,lj,lk)) + &
		 p2d3q27*( (myphi(li+1,lj+1,lk)-myphi(li-1,lj-1,lk))+ &
		 (myphi(li+1,lj-1,lk)-myphi(li-1,lj+1,lk))+ &
		 (myphi(li+1,lj,lk+1)-myphi(li-1,lj,lk-1))+ &
		 (myphi(li+1,lj,lk-1)-myphi(li-1,lj,lk+1)) )  + &
		 p3d3q27*((myphi(li+1,lj+1,lk+1)-myphi(li-1,lj-1,lk-1))+ &
		 (myphi(li+1,lj-1,lk-1)-myphi(li-1,lj+1,lk+1))+ &
		 (myphi(li+1,lj-1,lk+1)-myphi(li-1,lj+1,lk-1))+ &
		 (myphi(li+1,lj+1,lk-1)-myphi(li-1,lj-1,lk+1))))

	  grad_fiy=3.0_db*(p1d3q27*(myphi(li,lj+1,lk)-myphi(li,lj-1,lk)) + &
		 p2d3q27*((myphi(li+1,lj+1,lk)-myphi(li-1,lj-1,lk))+ &
		 (myphi(li-1,lj+1,lk)-myphi(li+1,lj-1,lk))+ &
		 (myphi(li,lj+1,lk+1)-myphi(li,lj-1,lk-1))+ &
		 (myphi(li,lj+1,lk-1)-myphi(li,lj-1,lk+1)) ) + &
		 p3d3q27*((myphi(li+1,lj+1,lk+1)-myphi(li-1,lj-1,lk-1))+ &
		 (myphi(li-1,lj+1,lk-1)-myphi(li+1,lj-1,lk+1))+ &
		 (myphi(li+1,lj+1,lk-1)-myphi(li-1,lj-1,lk+1))+ &
		 (myphi(li-1,lj+1,lk+1)-myphi(li+1,lj-1,lk-1))))

	  grad_fiz=3.0_db*(p1d3q27*(myphi(li,lj,lk+1)-myphi(li,lj,lk-1)) + &
		 p2d3q27*((myphi(li+1,lj,lk+1)-myphi(li-1,lj,lk-1))+ &
		 (myphi(li-1,lj,lk+1)-myphi(li+1,lj,lk-1))+ &
		 (myphi(li,lj+1,lk+1)-myphi(li,lj-1,lk-1))+ &
		 (myphi(li,lj-1,lk+1)-myphi(li,lj+1,lk-1)) ) + &
		 p3d3q27*((myphi(li+1,lj+1,lk+1)-myphi(li-1,lj-1,lk-1)) &
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
                   (2.0_db*invcssq)*(myphi(li,lj,lk)*(p0d3q27-1.0_db) + &
                   ( p1d3q27*(myphi(li+1,lj,lk)+myphi(li-1,lj,lk) + &
                   myphi(li,lj+1,lk)+myphi(li,lj-1,lk) + &
                   myphi(li,lj,lk+1)+myphi(li,lj,lk-1)) + &
                   p2d3q27*( (myphi(li+1,lj+1,lk)+myphi(li-1,lj-1,lk))+ &
                   (myphi(li+1,lj-1,lk)+myphi(li-1,lj+1,lk))+ &
                   (myphi(li+1,lj,lk+1)+myphi(li-1,lj,lk-1))+ &
                   (myphi(li+1,lj,lk-1)+myphi(li-1,lj,lk+1)) + &
                   (myphi(li,lj+1,lk+1)+myphi(li,lj-1,lk-1))+ &
                   (myphi(li,lj+1,lk-1)+myphi(li,lj-1,lk+1)) )  + &
                   p3d3q27*((myphi(li+1,lj+1,lk+1)+myphi(li-1,lj-1,lk-1))+ &
                   (myphi(li+1,lj-1,lk-1)+myphi(li-1,lj+1,lk+1))+ &
                   (myphi(li+1,lj-1,lk+1)+myphi(li-1,lj+1,lk-1))+ &
                   (myphi(li+1,lj+1,lk-1)+myphi(li-1,lj-1,lk+1)))))
      
      return
      
   endsubroutine compute_norm_interface_kernel_yminus
   
   attributes(global) subroutine compute_norm_interface_kernel_zplus(step, &
      iprobe,jprobe,kprobe,flop,nx,ny,nz,coords,isfluid,rho_r,rho_b, &
      ntotphifields,ntotauxfields,ntotlocauxfields,phifields_s,auxfields_s,locauxfields_s)

      implicit none
      integer :: step,iprobe,jprobe,kprobe,flop,nx,ny,nz
      integer, dimension(3) :: coords
      integer(kind=isf), dimension(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff) :: isfluid
      real(kind=db) :: rho_r,rho_b
      integer :: ntotphifields,ntotauxfields,ntotlocauxfields
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields,nblocks_d) :: phifields_s
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nauxfields,nblocks_d) :: auxfields_s
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nlocauxfields,nblocks_d) :: locauxfields_s
      
      real(kind=db), shared :: myphi(0:TILE_DIMx+1,0:TILE_DIMy+1,0:TILE_DIMz+1)
      real(kind=db):: grad_fix,grad_fiy,grad_fiz,mod_grad
      
      integer :: i,j,k,myblock,intblock
      integer :: ii,jj,kk
      integer :: li,lj,lk
      integer :: xblock,yblock,zblock
!      integer :: gi,gj,gk
      
            
      li = threadIdx%x-1
      lj = threadIdx%y-1
      lk = threadIdx%z-1

      i = (blockIdx%x-1) * TILE_DIMx + li
      j = (blockIdx%y-1) * TILE_DIMy + lj
      k = ((nzblock_d-2)-1) * TILE_DIMz + lk
      
      !gi=nx*coords(1)+i
      !gj=ny*coords(2)+j
      !gk=nz*coords(3)+k
      
      xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
	  yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
	  zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
      
      myblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
      ii=i-xblock*TILE_DIMx+2*TILE_DIMx
      jj=j-yblock*TILE_DIMy+2*TILE_DIMy
      kk=k-zblock*TILE_DIMz+2*TILE_DIMz

      !i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      !j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      !k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z
      
!      gi=nx*coords(1)+i
!      gj=ny*coords(2)+j
!      gk=nz*coords(3)+k
      
      intblock=blockIdx%x+blockIdx%y*nxblock_d+(nzblock_d-2)*nxyblock_d+1 !internal-node block
      
      myphi(li,lj,lk)=phifields_s(ii,jj,kk,1,myblock)
      
      call syncthreads
      
      if(abs(isfluid(i,j,k)) .ne. 1)return
      
	  !If my block index does not match the index of the internal-node block, it means my thread is on the outer halo and must exit.
	  if(myblock .ne. intblock)return

	  grad_fix=3.0_db*(p1d3q27*(myphi(li+1,lj,lk)-myphi(li-1,lj,lk)) + &
		 p2d3q27*( (myphi(li+1,lj+1,lk)-myphi(li-1,lj-1,lk))+ &
		 (myphi(li+1,lj-1,lk)-myphi(li-1,lj+1,lk))+ &
		 (myphi(li+1,lj,lk+1)-myphi(li-1,lj,lk-1))+ &
		 (myphi(li+1,lj,lk-1)-myphi(li-1,lj,lk+1)) )  + &
		 p3d3q27*((myphi(li+1,lj+1,lk+1)-myphi(li-1,lj-1,lk-1))+ &
		 (myphi(li+1,lj-1,lk-1)-myphi(li-1,lj+1,lk+1))+ &
		 (myphi(li+1,lj-1,lk+1)-myphi(li-1,lj+1,lk-1))+ &
		 (myphi(li+1,lj+1,lk-1)-myphi(li-1,lj-1,lk+1))))

	  grad_fiy=3.0_db*(p1d3q27*(myphi(li,lj+1,lk)-myphi(li,lj-1,lk)) + &
		 p2d3q27*((myphi(li+1,lj+1,lk)-myphi(li-1,lj-1,lk))+ &
		 (myphi(li-1,lj+1,lk)-myphi(li+1,lj-1,lk))+ &
		 (myphi(li,lj+1,lk+1)-myphi(li,lj-1,lk-1))+ &
		 (myphi(li,lj+1,lk-1)-myphi(li,lj-1,lk+1)) ) + &
		 p3d3q27*((myphi(li+1,lj+1,lk+1)-myphi(li-1,lj-1,lk-1))+ &
		 (myphi(li-1,lj+1,lk-1)-myphi(li+1,lj-1,lk+1))+ &
		 (myphi(li+1,lj+1,lk-1)-myphi(li-1,lj-1,lk+1))+ &
		 (myphi(li-1,lj+1,lk+1)-myphi(li+1,lj-1,lk-1))))

	  grad_fiz=3.0_db*(p1d3q27*(myphi(li,lj,lk+1)-myphi(li,lj,lk-1)) + &
		 p2d3q27*((myphi(li+1,lj,lk+1)-myphi(li-1,lj,lk-1))+ &
		 (myphi(li-1,lj,lk+1)-myphi(li+1,lj,lk-1))+ &
		 (myphi(li,lj+1,lk+1)-myphi(li,lj-1,lk-1))+ &
		 (myphi(li,lj-1,lk+1)-myphi(li,lj+1,lk-1)) ) + &
		 p3d3q27*((myphi(li+1,lj+1,lk+1)-myphi(li-1,lj-1,lk-1)) &
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
                   (2.0_db*invcssq)*(myphi(li,lj,lk)*(p0d3q27-1.0_db) + &
                   ( p1d3q27*(myphi(li+1,lj,lk)+myphi(li-1,lj,lk) + &
                   myphi(li,lj+1,lk)+myphi(li,lj-1,lk) + &
                   myphi(li,lj,lk+1)+myphi(li,lj,lk-1)) + &
                   p2d3q27*( (myphi(li+1,lj+1,lk)+myphi(li-1,lj-1,lk))+ &
                   (myphi(li+1,lj-1,lk)+myphi(li-1,lj+1,lk))+ &
                   (myphi(li+1,lj,lk+1)+myphi(li-1,lj,lk-1))+ &
                   (myphi(li+1,lj,lk-1)+myphi(li-1,lj,lk+1)) + &
                   (myphi(li,lj+1,lk+1)+myphi(li,lj-1,lk-1))+ &
                   (myphi(li,lj+1,lk-1)+myphi(li,lj-1,lk+1)) )  + &
                   p3d3q27*((myphi(li+1,lj+1,lk+1)+myphi(li-1,lj-1,lk-1))+ &
                   (myphi(li+1,lj-1,lk-1)+myphi(li-1,lj+1,lk+1))+ &
                   (myphi(li+1,lj-1,lk+1)+myphi(li-1,lj+1,lk-1))+ &
                   (myphi(li+1,lj+1,lk-1)+myphi(li-1,lj-1,lk+1)))))
      
      return
      
   endsubroutine compute_norm_interface_kernel_zplus
   
   attributes(global) subroutine compute_norm_interface_kernel_zminus(step, &
      iprobe,jprobe,kprobe,flop,nx,ny,nz,coords,isfluid,rho_r,rho_b, &
      ntotphifields,ntotauxfields,ntotlocauxfields,phifields_s,auxfields_s,locauxfields_s)

      implicit none
      integer :: step,iprobe,jprobe,kprobe,flop,nx,ny,nz
      integer, dimension(3) :: coords
      integer(kind=isf), dimension(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff) :: isfluid
      real(kind=db) :: rho_r,rho_b
      integer :: ntotphifields,ntotauxfields,ntotlocauxfields
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields,nblocks_d) :: phifields_s
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nauxfields,nblocks_d) :: auxfields_s
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nlocauxfields,nblocks_d) :: locauxfields_s
      
      real(kind=db), shared :: myphi(0:TILE_DIMx+1,0:TILE_DIMy+1,0:TILE_DIMz+1)
      real(kind=db):: grad_fix,grad_fiy,grad_fiz,mod_grad
      
      integer :: i,j,k,myblock,intblock
      integer :: ii,jj,kk
      integer :: li,lj,lk
      integer :: xblock,yblock,zblock
!      integer :: gi,gj,gk

      
      li = threadIdx%x-1
      lj = threadIdx%y-1
      lk = threadIdx%z-1

      i = (blockIdx%x-1) * TILE_DIMx + li
      j = (blockIdx%y-1) * TILE_DIMy + lj
      k = lk
      
      !gi=nx*coords(1)+i
      !gj=ny*coords(2)+j
      !gk=nz*coords(3)+k
      
      xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
	  yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
	  zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
      
      myblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
      ii=i-xblock*TILE_DIMx+2*TILE_DIMx
      jj=j-yblock*TILE_DIMy+2*TILE_DIMy
      kk=k-zblock*TILE_DIMz+2*TILE_DIMz

      !i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      !j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      !k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z
      
!      gi=nx*coords(1)+i
!      gj=ny*coords(2)+j
!      gk=nz*coords(3)+k
      
      intblock=blockIdx%x+blockIdx%y*nxblock_d+1*nxyblock_d+1 !internal-node block
      
      myphi(li,lj,lk)=phifields_s(ii,jj,kk,1,myblock)
      
      call syncthreads
      
      if(abs(isfluid(i,j,k)) .ne. 1)return
      
	  !If my block index does not match the index of the internal-node block, it means my thread is on the outer halo and must exit.
	  if(myblock .ne. intblock)return

	  grad_fix=3.0_db*(p1d3q27*(myphi(li+1,lj,lk)-myphi(li-1,lj,lk)) + &
		 p2d3q27*( (myphi(li+1,lj+1,lk)-myphi(li-1,lj-1,lk))+ &
		 (myphi(li+1,lj-1,lk)-myphi(li-1,lj+1,lk))+ &
		 (myphi(li+1,lj,lk+1)-myphi(li-1,lj,lk-1))+ &
		 (myphi(li+1,lj,lk-1)-myphi(li-1,lj,lk+1)) )  + &
		 p3d3q27*((myphi(li+1,lj+1,lk+1)-myphi(li-1,lj-1,lk-1))+ &
		 (myphi(li+1,lj-1,lk-1)-myphi(li-1,lj+1,lk+1))+ &
		 (myphi(li+1,lj-1,lk+1)-myphi(li-1,lj+1,lk-1))+ &
		 (myphi(li+1,lj+1,lk-1)-myphi(li-1,lj-1,lk+1))))

	  grad_fiy=3.0_db*(p1d3q27*(myphi(li,lj+1,lk)-myphi(li,lj-1,lk)) + &
		 p2d3q27*((myphi(li+1,lj+1,lk)-myphi(li-1,lj-1,lk))+ &
		 (myphi(li-1,lj+1,lk)-myphi(li+1,lj-1,lk))+ &
		 (myphi(li,lj+1,lk+1)-myphi(li,lj-1,lk-1))+ &
		 (myphi(li,lj+1,lk-1)-myphi(li,lj-1,lk+1)) ) + &
		 p3d3q27*((myphi(li+1,lj+1,lk+1)-myphi(li-1,lj-1,lk-1))+ &
		 (myphi(li-1,lj+1,lk-1)-myphi(li+1,lj-1,lk+1))+ &
		 (myphi(li+1,lj+1,lk-1)-myphi(li-1,lj-1,lk+1))+ &
		 (myphi(li-1,lj+1,lk+1)-myphi(li+1,lj-1,lk-1))))

	  grad_fiz=3.0_db*(p1d3q27*(myphi(li,lj,lk+1)-myphi(li,lj,lk-1)) + &
		 p2d3q27*((myphi(li+1,lj,lk+1)-myphi(li-1,lj,lk-1))+ &
		 (myphi(li-1,lj,lk+1)-myphi(li+1,lj,lk-1))+ &
		 (myphi(li,lj+1,lk+1)-myphi(li,lj-1,lk-1))+ &
		 (myphi(li,lj-1,lk+1)-myphi(li,lj+1,lk-1)) ) + &
		 p3d3q27*((myphi(li+1,lj+1,lk+1)-myphi(li-1,lj-1,lk-1)) &
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
                   (2.0_db*invcssq)*(myphi(li,lj,lk)*(p0d3q27-1.0_db) + &
                   ( p1d3q27*(myphi(li+1,lj,lk)+myphi(li-1,lj,lk) + &
                   myphi(li,lj+1,lk)+myphi(li,lj-1,lk) + &
                   myphi(li,lj,lk+1)+myphi(li,lj,lk-1)) + &
                   p2d3q27*( (myphi(li+1,lj+1,lk)+myphi(li-1,lj-1,lk))+ &
                   (myphi(li+1,lj-1,lk)+myphi(li-1,lj+1,lk))+ &
                   (myphi(li+1,lj,lk+1)+myphi(li-1,lj,lk-1))+ &
                   (myphi(li+1,lj,lk-1)+myphi(li-1,lj,lk+1)) + &
                   (myphi(li,lj+1,lk+1)+myphi(li,lj-1,lk-1))+ &
                   (myphi(li,lj+1,lk-1)+myphi(li,lj-1,lk+1)) )  + &
                   p3d3q27*((myphi(li+1,lj+1,lk+1)+myphi(li-1,lj-1,lk-1))+ &
                   (myphi(li+1,lj-1,lk-1)+myphi(li-1,lj+1,lk+1))+ &
                   (myphi(li+1,lj-1,lk+1)+myphi(li-1,lj+1,lk-1))+ &
                   (myphi(li+1,lj+1,lk-1)+myphi(li-1,lj-1,lk+1)))))
      
      return
      
   endsubroutine compute_norm_interface_kernel_zminus
  
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
       (( p1d3q27*(myarrx(li+1,lj,lk)-myarrx(li-1,lj,lk)) + &
       p2d3q27*((myarrx(li+1,lj+1,lk)-myarrx(li-1,lj-1,lk))+(myarrx(li+1,lj-1,lk)-myarrx(li-1,lj+1,lk))+ &
       (myarrx(li+1,lj,lk+1)-myarrx(li-1,lj,lk-1))+(myarrx(li+1,lj,lk-1)-myarrx(li-1,lj,lk+1)))+ &
       p3d3q27*((myarrx(li+1,lj+1,lk+1)-myarrx(li-1,lj-1,lk-1))+(myarrx(li+1,lj-1,lk-1)-myarrx(li-1,lj+1,lk+1))+ &
       (myarrx(li+1,lj-1,lk+1)-myarrx(li-1,lj+1,lk-1))+(myarrx(li+1,lj+1,lk-1)-myarrx(li-1,lj-1,lk+1))))+ &
       (p1d3q27*(myarry(li,lj+1,lk)-myarry(li,lj-1,lk)) + &
       p2d3q27*((myarry(li+1,lj+1,lk)-myarry(li-1,lj-1,lk))+(myarry(li-1,lj+1,lk)-myarry(li+1,lj-1,lk))+ &
       (myarry(li,lj+1,lk+1)-myarry(li,lj-1,lk-1))+(myarry(li,lj+1,lk-1)-myarry(li,lj-1,lk+1)))+ &
       p3d3q27*((myarry(li+1,lj+1,lk+1)-myarry(li-1,lj-1,lk-1))+(myarry(li-1,lj+1,lk-1)-myarry(li+1,lj-1,lk+1))+ &
       (myarry(li+1,lj+1,lk-1)-myarry(li-1,lj-1,lk+1))+(myarry(li-1,lj+1,lk+1)-myarry(li+1,lj-1,lk-1))))+ &
       (p1d3q27*(myarrz(li,lj,lk+1)-myarrz(li,lj,lk-1)) + &
       p2d3q27*((myarrz(li+1,lj,lk+1)-myarrz(li-1,lj,lk-1))+(myarrz(li-1,lj,lk+1)-myarrz(li+1,lj,lk-1))+ &
       (myarrz(li,lj+1,lk+1)-myarrz(li,lj-1,lk-1))+(myarrz(li,lj-1,lk+1)-myarrz(li,lj+1,lk-1)))+ &
       p3d3q27*((myarrz(li+1,lj+1,lk+1)-myarrz(li-1,lj-1,lk-1))+(myarrz(li-1,lj-1,lk+1)-myarrz(li+1,lj+1,lk-1))+ &
       (myarrz(li+1,lj-1,lk+1)-myarrz(li-1,lj+1,lk-1))+(myarrz(li-1,lj+1,lk+1)-myarrz(li+1,lj-1,lk-1)))))*invcssq
      
      return
      
   endsubroutine compute_div_theta_n_kernel
#endif 
endmodule lb_cuda_auxfields
