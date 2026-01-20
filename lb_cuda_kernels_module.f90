#include "defines.h"
#if !defined(_OPENACC)  
#error "To use this module the macros _OPENACC should be defined."
#endif


module lb_cuda_kernels

   use vars
   use iso_c_binding
   use cudafor
   use mpi_template, only: coords,dostop,doerror,mydev,myrank,nprocs,nbuff,nbuffbvec
   use lb_cuda_vars
   use lb_cuda_fused, only: fused_LB_kernel

   implicit none
   

contains

   subroutine setup_cuda
      implicit none
      
      istat = cudaSetDevice(mydev)
      if (istat/=0) then
        if(myrank==0)write(6,*) 'status after cudaSetDevice:', cudaGetErrorString(istat)
        call dostop
      endif
      
      istat = cudaGetDeviceProperties(prop, mydev)
	
      mshared = prop%sharedMemPerBlock
    
	  if(myrank==0)call printDeviceProperties(prop,6, mydev)
        

      if(TILE_DIMx==0 .or. TILE_DIMy==0 .or. TILE_DIMz==0 .or. TILE_DIM==0)then
        call doerror(6,'one of TILE_DIM values is zero')
      endif

      ! --- Check that the product of the TILE_DIMs does not exceed 1024 ---
      if (TILE_DIMx * TILE_DIMy * TILE_DIMz > 1024) then
        if(myrank==0) then
          write(6,*) 'ERROR: TILE_DIMx*TILE_DIMy*TILE_DIMz =',TILE_DIMx*TILE_DIMy*TILE_DIMz,' > 1024'
          write(6,*) 'Decrease at least one of the TILE_DIM values!'
        endif
        call doerror(6,'TILE_DIM product exceeds 1024')
      endif

      dimGrid  = dim3((nx+TILE_DIMx-1)/TILE_DIMx,(ny+TILE_DIMy-1)/TILE_DIMy,(nz+TILE_DIMz-1)/TILE_DIMz)
      dimBlock = dim3(TILE_DIMx, TILE_DIMy, TILE_DIMz)
      lx_d=lx
      ly_d=ly
      lz_d=lz
      nx_d=nx
      ny_d=ny
      nz_d=nz
      TILE_DIMx_d=TILE_DIMx
      TILE_DIMy_d=TILE_DIMy
      TILE_DIMz_d=TILE_DIMz
      TILE_DIM_d=TILE_DIM
      if (mod(nx, TILE_DIMx)/= 0) then
        if(myrank==0)write(6,*) 'nx must be a multiple of TILE_DIM'
        call dostop
      end if
      if (mod(ny, TILE_DIMy) /= 0) then
        if(myrank==0)write(6,*) 'ny must be a multiple of TILE_DIMy'
        call dostop
      end if
      if (mod(nz, TILE_DIMz) /= 0) then
        if(myrank==0)write(6,*) 'nz must be a multiple of TILE_DIMz'
        call dostop
      end if
      
      dimGridInt = dim3((nx+TILE_DIMx-1)/TILE_DIMx -2,(ny+TILE_DIMy-1)/TILE_DIMy -2,(nz+TILE_DIMz-1)/TILE_DIMz -2)
      dimGridhalo  = dim3((nx+TILE_DIMx-1)/TILE_DIMx +2,(ny+TILE_DIMy-1)/TILE_DIMy +2,(nz+TILE_DIMz-1)/TILE_DIMz +2)
      dimBlockhalo = dim3(TILE_DIMx, TILE_DIMy, TILE_DIMz)
      
      dimBlockshared = dim3(TILE_DIMx +2, TILE_DIMy +2, TILE_DIMz +2)
      
      dimGridx  = dim3(1,(ny+TILE_DIM-1)/TILE_DIM -2, (nz+TILE_DIM-1)/TILE_DIM -2) !only yz faces
      dimGridy  = dim3((nx+TILE_DIM-1)/TILE_DIM -2, 1, (nz+TILE_DIM-1)/TILE_DIM)  !xz faces also doing edge xy
      dimGridz  = dim3((nx+TILE_DIM-1)/TILE_DIM, (ny+TILE_DIM-1)/TILE_DIM, 1)    !xy faces also doing edges xz yz and corners
      
      
      dimBlock2 = dim3(TILE_DIM, TILE_DIM, 1)
      !plus 2 for the halo forward and backward
      nxblock=nx/TILE_DIMx +2
      nyblock=ny/TILE_DIMy +2
      nzblock=nz/TILE_DIMz +2
      
      nxyblock=nxblock*nyblock
      nblocks=nxblock*nyblock*nzblock
      
      nxblock_d=nxblock
      nyblock_d=nyblock
      nzblock_d=nzblock
      
      nxyblock_d=nxyblock
      nblocks_d=nblocks

#if 1
      if(myrank==0)then
        write(6,*)'nx,ny,nz',nx,ny,nz
        write(6,*)'TILE_DIMx,TILE_DIMy,TILE_DIMz',TILE_DIMx,TILE_DIMy,TILE_DIMz
        write(6,*)'nxblock,nyblock,nzblock',nxblock,nyblock,nzblock
        write(6,*)'nblocks',nblocks
      endif
#endif
      
   
   endsubroutine setup_cuda
   
   subroutine test_LB_cuda

   implicit none
      
!  if(myrank==0)write(6,*)'step ',step, __LINE__ , __FILE__
   !$acc wait
   istat = cudaDeviceSynchronize
       
   !$acc host_data use_device(myrank,nx,ny,nz,coords)
   !call test_LB_kernel_shared_x<<<dimGridx, dimBlockshared>>>(myrank,nx,ny,nz,coords)
   !call test_LB_kernel_shared_y<<<dimGridy, dimBlockshared>>>(myrank,nx,ny,nz,coords)
   !call test_LB_kernel_shared_z<<<dimGridz, dimBlockshared>>>(myrank,nx,ny,nz,coords)
   call  test_LB_kernel_shared_internal<<<dimGridInt, dimBlockshared>>>(myrank,nx,ny,nz,coords)
   !call test_LB_kernel_shared<<<dimGrid, dimBlockshared>>>(myrank,nx,ny,nz,coords)
   !call test_LB_kernel_halo<<<dimGridhalo,dimBlockhalo>>>(myrank,nx,ny,nz,coords)
   !$acc end host_data
   
   end subroutine test_LB_cuda
   
    
   subroutine printDeviceProperties(prop,iu,num)
      
          use cudafor
          type(cudadeviceprop) :: prop
          integer,intent(in) :: iu,num 
          real(kind=db), parameter :: convG=1024.0_db*1024.0_db*1024.0_db
          
          write(iu,907)"                                                                               "
          write(iu,907)"*****************************GPU FEATURE MONITOR*******************************"
          write(iu,907)"                                                                               "
          
          write (iu,900) "Device Number: "      ,num
          write (iu,901) "Device Name: "        ,trim(prop%name)
          write (iu,903) "Total Global Memory: ",real(prop%totalGlobalMem)/convG," Gbytes"
          write (iu,902) "sharedMemPerBlock: "  ,prop%sharedMemPerBlock," bytes"
          write (iu,900) "regsPerBlock: "       ,prop%regsPerBlock
          write (iu,900) "warpSize: "           ,prop%warpSize
          write (iu,900) "maxThreadsPerBlock: " ,prop%maxThreadsPerBlock
          write (iu,904) "maxThreadsDim: "      ,prop%maxThreadsDim
          write (iu,904) "maxGridSize: "        ,prop%maxGridSize
          write (iu,903) "ClockRate: "          ,real(prop%clockRate)/1e6," GHz"
          write (iu,902) "Total Const Memory: " ,prop%totalConstMem," bytes"
          write (iu,905) "Compute Capability Revision: ",prop%major,prop%minor
          write (iu,902) "TextureAlignment: "   ,prop%textureAlignment," bytes"
          write (iu,906) "deviceOverlap: "      ,prop%deviceOverlap
          write (iu,900) "multiProcessorCount: ",prop%multiProcessorCount
          write (iu,906) "integrated: "         ,prop%integrated
          write (iu,906) "canMapHostMemory: "   ,prop%canMapHostMemory
          write (iu,906) "ECCEnabled: "         ,prop%ECCEnabled
          write (iu,906) "UnifiedAddressing: "  ,prop%unifiedAddressing
          write (iu,900) "L2 Cache Size: "      ,prop%l2CacheSize
          write (iu,900) "maxThreadsPerSMP: "   ,prop%maxThreadsPerMultiProcessor
          
          write(iu,907)"                                                                               "
          write(iu,907)"*******************************************************************************"
          write(iu,907)"                                                                               "
          
          900 format (a,i0)
          901 format (a,a)
          902 format (a,i0,a)
          903 format (a,f16.8,a)
          904 format (a,2(i0,1x,'x',1x),i0)
          905 format (a,i0,'.',i0)
          906 format (a,l0)
          907 format (a)
          
          return
      
  end subroutine printDeviceProperties
   



   !****************************************************************************!


 attributes(global) subroutine test_LB_kernel(myrank,nx,ny,nz,coords)
      implicit none
      
      integer :: myrank,nx,ny,nz
      integer, dimension(3) :: coords
      
      integer :: i,j,k,gi,gj,gk,myblock

      i = (blockIdx%x-1) * TILE_DIMx_d + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy_d + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz_d + threadIdx%z
      
      gi=nx*coords(1)+i
      gj=ny*coords(2)+j
      gk=nz*coords(3)+k
      
      myblock=blockIdx%x+blockIdx%y*nxblock_d+blockIdx%z*nxyblock_d+1
      
      !if(myblock<344)then
      if(gi==1 .and. gj==1 .and. gk==31)then
        write(*,*)'i',gi,gj,gk,myblock,myrank
      endif
     !if(gi==1) write(*,*)'ciao',i,j,k,phi(i,j,k)
     
     !phi_old(i,j,k)=phi(i,j,k)
     
      
 end subroutine test_LB_kernel
 
 attributes(global) subroutine test_LB_kernel_shared(myrank,nx,ny,nz,coords)
      implicit none
      
      integer :: myrank,nx,ny,nz
      integer, dimension(3) :: coords
      
      integer :: i,j,k,gi,gj,gk,myblock
      integer :: li,lj,lk
      integer :: xblock,yblock,zblock
      
      li = threadIdx%x-1
      lj = threadIdx%y-1
      lk = threadIdx%z-1

      i = (blockIdx%x-1) * TILE_DIMx_d + li
      j = (blockIdx%y-1) * TILE_DIMy_d + lj
      k = (blockIdx%z-1) * TILE_DIMz_d + lk
      
      gi=nx*coords(1)+i
      gj=ny*coords(2)+j
      gk=nz*coords(3)+k
      
      xblock=(i+2*TILE_DIMx_d-1)/TILE_DIMx_d
	  yblock=(j+2*TILE_DIMy_d-1)/TILE_DIMy_d
	  zblock=(k+2*TILE_DIMz_d-1)/TILE_DIMz_d
      
      myblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
      
      !if(myblock<344)then
      if(gi==1 .and. gj==1 .and. gk==36)then
        write(*,*)'i',gi,gj,gk,myrank,lk
      endif
     !if(gi==1) write(*,*)'ciao',i,j,k,phi(i,j,k)
     
     !phi_old(i,j,k)=phi(i,j,k)
     
      
 end subroutine test_LB_kernel_shared 
 
 attributes(global) subroutine test_LB_kernel_shared_internal(myrank,nx,ny,nz,coords)
      implicit none
      
      integer :: myrank,nx,ny,nz
      integer, dimension(3) :: coords
      
      integer :: i,j,k,gi,gj,gk,myblock
      integer :: li,lj,lk
      integer :: xblock,yblock,zblock
      
      li = threadIdx%x-1
      lj = threadIdx%y-1
      lk = threadIdx%z-1

      i = (blockIdx%x) * TILE_DIMx_d + li
      j = (blockIdx%y) * TILE_DIMy_d + lj
      k = (blockIdx%z) * TILE_DIMz_d + lk
      
      gi=nx*coords(1)+i
      gj=ny*coords(2)+j
      gk=nz*coords(3)+k
      
      xblock=(i+2*TILE_DIMx_d-1)/TILE_DIMx_d
	  yblock=(j+2*TILE_DIMy_d-1)/TILE_DIMy_d
	  zblock=(k+2*TILE_DIMz_d-1)/TILE_DIMz_d
      
      myblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
      
      
      !if(myblock<344)then
      if(gi==17 .and. gj==17 .and. gk==17)then
        write(*,*)'i',li,lj,lk,myrank
      endif
     !if(gi==1) write(*,*)'ciao',i,j,k,phi(i,j,k)
     
     !phi_old(i,j,k)=phi(i,j,k)
     
      
 end subroutine test_LB_kernel_shared_internal 
 
 attributes(global) subroutine test_LB_kernel_shared_z(myrank,nx,ny,nz,coords)
      implicit none
      
      integer :: myrank,nx,ny,nz
      integer, dimension(3) :: coords
      
      integer :: i,j,k,gi,gj,gk,myblock
      integer :: li,lj,lk
      integer :: xblock,yblock,zblock
      
      li = threadIdx%x-1
      lj = threadIdx%y-1
      lk = threadIdx%z-1

      i = (blockIdx%x-1) * TILE_DIMx_d + li
      j = (blockIdx%y-1) * TILE_DIMy_d + lj
      
      gi=nx*coords(1)+i
      gj=ny*coords(2)+j
      
      xblock=(i+2*TILE_DIMx_d-1)/TILE_DIMx_d
	  yblock=(j+2*TILE_DIMy_d-1)/TILE_DIMy_d
	  
	  
	  k = lk
	  gk=nz*coords(3)+k
	  zblock=(k+2*TILE_DIMz_d-1)/TILE_DIMz_d
      
      myblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
      
      !if(myblock<344)then
      if(gi==1 .and. gj==1 .and. gk==38)then
        write(*,*)'asd',gi,gj,gk,myrank,lk
      endif
     !if(gi==1) write(*,*)'ciao',i,j,k,phi(i,j,k)
     
     !phi_old(i,j,k)=phi(i,j,k)
     
      k = ((nzblock_d-2)-1) * TILE_DIMz_d + lk
	  gk=nz*coords(3)+k
	  zblock=(k+2*TILE_DIMz_d-1)/TILE_DIMz_d
      
      myblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
      
      !if(myblock<344)then
      if(gi==1 .and. gj==1 .and. gk==36)then
        write(*,*)'asf',gi,gj,gk,myrank,lk
      endif
      
 end subroutine test_LB_kernel_shared_z 
 
 attributes(global) subroutine test_LB_kernel_shared_y(myrank,nx,ny,nz,coords)
      implicit none
      
      integer :: myrank,nx,ny,nz
      integer, dimension(3) :: coords
      
      integer :: i,j,k,gi,gj,gk,myblock
      integer :: li,lj,lk
      integer :: xblock,yblock,zblock
      
      li = threadIdx%x-1
      lj = threadIdx%y-1
      lk = threadIdx%z-1

      i = (blockIdx%x-1) * TILE_DIMx_d + li + TILE_DIMx_d
      k = (blockIdx%z-1) * TILE_DIMz_d + lk
      
      gi=nx*coords(1)+i
      gk=nz*coords(3)+k
      
      xblock=(i+2*TILE_DIMx_d-1)/TILE_DIMx_d
	  zblock=(k+2*TILE_DIMz_d-1)/TILE_DIMz_d
	  
	  
      j = lj
      gj=ny*coords(2)+j
      yblock=(j+2*TILE_DIMy_d-1)/TILE_DIMy_d
      
      myblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
      
      !if(myblock<344)then
      if(gi==8 .and. gj==6 .and. gk==1)then
        write(*,*)'i',gi,gj,gk,myrank,li
      endif
     !if(gi==1) write(*,*)'ciao',i,j,k,phi(i,j,k)
     
     !phi_old(i,j,k)=phi(i,j,k)
 
      j = ((nyblock_d-2)-1) * TILE_DIMy_d + lj
      gj=ny*coords(2)+j
      yblock=(j+2*TILE_DIMy_d-1)/TILE_DIMy_d
      
      myblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
      
      !if(myblock<344)then
      if(gi==5 .and. gj==5 .and. gk==1)then
        write(*,*)'ii',gi,gj,gk,myrank,lj
      endif     
      
 end subroutine test_LB_kernel_shared_y
 
 attributes(global) subroutine test_LB_kernel_shared_x(myrank,nx,ny,nz,coords)
      implicit none
      
      integer :: myrank,nx,ny,nz
      integer, dimension(3) :: coords
      
      integer :: i,j,k,gi,gj,gk,myblock
      integer :: li,lj,lk
      integer :: xblock,yblock,zblock
      
      li = threadIdx%x-1
      lj = threadIdx%y-1
      lk = threadIdx%z-1

      
      j = (blockIdx%y-1) * TILE_DIMy_d + lj + TILE_DIMy_d
      k = (blockIdx%z-1) * TILE_DIMz_d + lk + TILE_DIMz_d
      !if(myrank==0)write(*,*)'ciao',j
      
     
      gj=ny*coords(2)+j
      gk=nz*coords(3)+k
      
      
	  yblock=(j+2*TILE_DIMy_d-1)/TILE_DIMy_d
	  zblock=(k+2*TILE_DIMz_d-1)/TILE_DIMz_d
	  
	  i = li
	  gi=nx*coords(1)+i
	  xblock=(i+2*TILE_DIMx_d-1)/TILE_DIMx_d
      
      myblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
      
      !if(myblock<344)then
      if(gi==1 .and. gj==8 .and. gk==8)then
        write(*,*)'i',gi,gj,gk,lj,lk
      endif
     !if(gi==1) write(*,*)'ciao',i,j,k,phi(i,j,k)
     
     !phi_old(i,j,k)=phi(i,j,k)
 
	  i = ((nxblock_d-2)-1) * TILE_DIMx_d + li
	  gi=nx*coords(1)+i
	  xblock=(i+2*TILE_DIMx_d-1)/TILE_DIMx_d
      
      myblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
      
      !if(myblock<344)then
      if(gi==1 .and. gj==1 .and. gk==36)then
        write(*,*)'i',gi,gj,gk,myrank,lk
      endif     
      
 end subroutine test_LB_kernel_shared_x 
 
 attributes(global) subroutine test_LB_kernel_halo(myrank,nx,ny,nz,coords)
      implicit none
      
      integer :: myrank,nx,ny,nz
      integer, dimension(3) :: coords
      
      integer :: i,j,k,gi,gj,gk,myblock
      integer :: xblock,yblock,zblock

      i = (blockIdx%x-2) * TILE_DIMx_d + threadIdx%x
      j = (blockIdx%y-2) * TILE_DIMy_d + threadIdx%y
      k = (blockIdx%z-2) * TILE_DIMz_d + threadIdx%z
      
      gi=nx*coords(1)+i
      gj=ny*coords(2)+j
      gk=nz*coords(3)+k
      
      xblock=(i+2*TILE_DIMx_d-1)/TILE_DIMx_d
	  yblock=(j+2*TILE_DIMy_d-1)/TILE_DIMy_d
	  zblock=(k+2*TILE_DIMz_d-1)/TILE_DIMz_d
      
      myblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
      
      !if(myblock==20 .and. myrank==1)then
      !if(blockIdx%x==1 .and. blockIdx%y==1 .and. blockIdx%z==1)then
      if(gi==1 .and. gj==1 .and. gk==31)then
        write(*,*)'e',gi,gj,gk,myblock,myrank
      endif
     !if(gi==1) write(*,*)'ciao',i,j,k,phi(i,j,k)
     
     !phi_old(i,j,k)=phi(i,j,k)
     
      
 end subroutine test_LB_kernel_halo
 
 subroutine compute_densityratio_cuda(phifields_s)

      implicit none
      real(kind=db), allocatable, dimension(:) :: phifields_s
      
!      if(myrank==0)write(6,*)'step ',step, __LINE__ , __FILE__

#ifdef DENSRATIO
      !$acc wait
      istat = cudaDeviceSynchronize

!$acc host_data use_device(flop,nx,ny,nz,coords,isfluid &
       !$acc& ,rho_r,rho_b,ntotphifields,ntotlocauxfields,phifields_s,locauxfields)
       call compute_densityratio_kernel<<<dimGrid, dimBlock>>>(flop,nx,ny,nz,coords,isfluid, &
        rho_r,rho_b,ntotphifields,ntotlocauxfields,phifields_s,locauxfields)
!$acc end host_data
      istat = cudaDeviceSynchronize
      istat = cudaGetLastError()
      if (istat/= cudaSuccess) then
        if(myrank==0)write(6,*) 'status after at ', __LINE__ ,' of file ', __FILE__ ,' :'
        if(myrank==0)write(6,*) cudaGetErrorString(istat)
        call doerror(6,'ERROR in compute_densityratio_cuda')
      endif
      !$acc wait        
#endif

   endsubroutine compute_densityratio_cuda
 
 attributes(global) subroutine compute_densityratio_kernel(flop,nx,ny,nz,coords,isfluid, &
  rho_r,rho_b,ntotphifields,ntotlocauxfields,phifields_s,locauxfields_s)
      implicit none
      
      integer :: flop,nx,ny,nz
      integer, dimension(3) :: coords
      integer(kind=isf), dimension(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff) :: isfluid
      real(kind=db) :: rho_r,rho_b
      integer :: ntotphifields,ntotlocauxfields
      real(kind=db), dimension(ntotphifields) :: phifields_s
      real(kind=db), dimension(ntotlocauxfields) :: locauxfields_s
      
      integer :: i,j,k,gi,gj,gk,myblock,ii,jj,kk
      real(kind=db) :: phitemp
      
      i = (blockIdx%x-1) * TILE_DIMx_d + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy_d + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz_d + threadIdx%z
      
      if(abs(isfluid(i,j,k)) .ne. 1)return
      
      gi=nx*coords(1)+i
      gj=ny*coords(2)+j
      gk=nz*coords(3)+k
      
      myblock=blockIdx%x+blockIdx%y*nxblock_d+blockIdx%z*nxyblock_d+1
      
      ii=threadIdx%x
      jj=threadIdx%y
      kk=threadIdx%z
      
      phitemp=phifields_s(idx5d(ii,jj,kk,1,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nphifields))
      
      return
     
      
 end subroutine compute_densityratio_kernel 
 
 subroutine compute_norm_interface_cuda(phifields_s)

      implicit none
      real(kind=db), allocatable, dimension(:) :: phifields_s
      
      integer(8) :: mymshared
      
      mymshared = 1 * (TILE_DIMx + 2) * (TILE_DIMy + 2) * (TILE_DIMz + 2) * db
      
#ifdef TWOCOMPONENT
      !$acc wait
      istat = cudaDeviceSynchronize

!$acc host_data use_device(flop,nx,ny,nz,coords,isfluid &
       !$acc& ,rho_r,rho_b,ntotphifields,ntotauxfields,ntotlocauxfields,phifields_s,auxfields,locauxfields)
       call compute_norm_interface_kernel<<<dimGrid, dimBlockshared, mymshared>>>(flop,nx,ny,nz,coords,isfluid, &
        rho_r,rho_b,ntotphifields,ntotauxfields,ntotlocauxfields,phifields_s,auxfields,locauxfields)
!$acc end host_data
      istat = cudaDeviceSynchronize
      istat = cudaGetLastError()
      if (istat/= cudaSuccess) then
        if(myrank==0)write(6,*) 'status after at ', __LINE__ ,' of file ', __FILE__ ,' :'
        if(myrank==0)write(6,*) cudaGetErrorString(istat)
        call doerror(6,'ERROR in compute_norm_interface_cuda')
      endif
      !$acc wait        
#endif
      
      return
      
 endsubroutine compute_norm_interface_cuda
 
 attributes(global) subroutine compute_norm_interface_kernel(flop,nx,ny,nz,coords,isfluid, &
  rho_r,rho_b,ntotphifields,ntotauxfields,ntotlocauxfields,phifields_s,auxfields_s,locauxfields_s)

      implicit none
      integer :: flop,nx,ny,nz
      integer, dimension(3) :: coords
      integer(kind=isf), dimension(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff) :: isfluid
      real(kind=db) :: rho_r,rho_b
      integer :: ntotphifields,ntotauxfields,ntotlocauxfields
      real(kind=db), dimension(ntotphifields) :: phifields_s
      real(kind=db), dimension(ntotauxfields) :: auxfields_s
      real(kind=db), dimension(ntotlocauxfields) :: locauxfields_s
      
      real(kind=db), shared :: myphi(0:TILE_DIMx_d+1,0:TILE_DIMy_d+1,0:TILE_DIMz_d+1)
      real(kind=db):: grad_fix,grad_fiy,grad_fiz,mod_grad
      
      integer :: i,j,k,gi,gj,gk,myblock,idblock
      integer :: ii,jj,kk
      integer :: li,lj,lk
      integer :: xblock,yblock,zblock
      
      li = threadIdx%x-1
      lj = threadIdx%y-1
      lk = threadIdx%z-1

      i = (blockIdx%x-1) * TILE_DIMx_d + li
      j = (blockIdx%y-1) * TILE_DIMy_d + lj
      k = (blockIdx%z-1) * TILE_DIMz_d + lk
      
      gi=nx*coords(1)+i
      gj=ny*coords(2)+j
      gk=nz*coords(3)+k
      
      xblock=(i+2*TILE_DIMx_d-1)/TILE_DIMx_d
	  yblock=(j+2*TILE_DIMy_d-1)/TILE_DIMy_d
	  zblock=(k+2*TILE_DIMz_d-1)/TILE_DIMz_d
      
      myblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
      ii=i-xblock*TILE_DIMx_d+2*TILE_DIMx_d
      jj=j-yblock*TILE_DIMy_d+2*TILE_DIMy_d
      kk=k-zblock*TILE_DIMz_d+2*TILE_DIMz_d
      
      myphi(li,lj,lk)=phifields_s(idx5d(ii,jj,kk,1,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nphifields))
      
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

	  auxfields_s(idx5d(ii,jj,kk,1,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nauxfields))= &
	   grad_fix/(mod_grad+1.0e-9_db)
	  auxfields_s(idx5d(ii,jj,kk,2,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nauxfields))= &
	   grad_fiy/(mod_grad+1.0e-9_db)
	  auxfields_s(idx5d(ii,jj,kk,3,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nauxfields))= &
	   grad_fiz/(mod_grad+1.0e-9_db)
	  
	  auxfields_s(idx5d(ii,jj,kk,4,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nauxfields))=mod_grad

	  auxfields_s(idx5d(ii,jj,kk,5,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nauxfields))= &
	   myphi(li,lj,lk)*(1.0_db-myphi(li,lj,lk))*(grad_fix/(mod_grad+1.0e-9_db))
	  auxfields_s(idx5d(ii,jj,kk,6,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nauxfields))= &
	   myphi(li,lj,lk)*(1.0_db-myphi(li,lj,lk))*(grad_fiy/(mod_grad+1.0e-9_db))
	  auxfields_s(idx5d(ii,jj,kk,7,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nauxfields))= &
	   myphi(li,lj,lk)*(1.0_db-myphi(li,lj,lk))*(grad_fiz/(mod_grad+1.0e-9_db))
	   
      !lap_phi here
      locauxfields_s(idx5d(ii,jj,kk,1,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nlocauxfields))= &
                   (2.0_db/cssq)*(myphi(li,lj,lk)*(p0-1.0_db) + &
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
   
   subroutine compute_div_thetan(phifields_s)

      implicit none
      real(kind=db), allocatable, dimension(:) :: phifields_s
      integer(8) :: mymshared
      
      mymshared = 3 * (TILE_DIMx + 2) * (TILE_DIMy + 2) * (TILE_DIMz + 2) * db
      
#ifdef TWOCOMPONENT
      !$acc wait
      istat = cudaDeviceSynchronize

!$acc host_data use_device(flop,nx,ny,nz,coords,isfluid &
       !$acc& ,rho_r,rho_b,ntotphifields,ntotauxfields,ntotlocauxfields,phifields_s,auxfields,locauxfields)
       call compute_div_thetan_kernel<<<dimGrid, dimBlockshared, mymshared>>>(flop,nx,ny,nz,coords,isfluid, &
        rho_r,rho_b,ntotphifields,ntotauxfields,ntotlocauxfields,phifields_s,auxfields,locauxfields)
!$acc end host_data
      istat = cudaDeviceSynchronize
      istat = cudaGetLastError()
      if (istat/= cudaSuccess) then
        if(myrank==0)write(6,*) 'status after at ', __LINE__ ,' of file ', __FILE__ ,' :'
        if(myrank==0)write(6,*) cudaGetErrorString(istat)
        call doerror(6,'ERROR in compute_div_thetan')
      endif
      !$acc wait        
#endif
      
      return
      
   endsubroutine compute_div_thetan 
   
   attributes(global) subroutine compute_div_thetan_kernel(flop,nx,ny,nz,coords,isfluid, &
     rho_r,rho_b,ntotphifields,ntotauxfields,ntotlocauxfields,phifields_s,auxfields_s,locauxfields_s)

      implicit none
      integer :: flop,nx,ny,nz
      integer, dimension(3) :: coords
      integer(kind=isf), dimension(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff) :: isfluid
      real(kind=db) :: rho_r,rho_b
      integer :: ntotphifields,ntotauxfields,ntotlocauxfields
      real(kind=db), dimension(ntotphifields) :: phifields_s
      real(kind=db), dimension(ntotauxfields) :: auxfields_s
      real(kind=db), dimension(ntotlocauxfields) :: locauxfields_s
      
      real(kind=db), shared :: myarrx(0:TILE_DIMx_d+1,0:TILE_DIMy_d+1,0:TILE_DIMz_d+1)
      real(kind=db), shared :: myarry(0:TILE_DIMx_d+1,0:TILE_DIMy_d+1,0:TILE_DIMz_d+1)
      real(kind=db), shared :: myarrz(0:TILE_DIMx_d+1,0:TILE_DIMy_d+1,0:TILE_DIMz_d+1)
      
      integer :: i,j,k,gi,gj,gk,myblock,idblock
      integer :: ii,jj,kk
      integer :: li,lj,lk
      integer :: xblock,yblock,zblock
      
      li = threadIdx%x-1
      lj = threadIdx%y-1
      lk = threadIdx%z-1
      
      i = (blockIdx%x-1) * TILE_DIMx_d + li
      j = (blockIdx%y-1) * TILE_DIMy_d + lj
      k = (blockIdx%z-1) * TILE_DIMz_d + lk
      
      gi=nx*coords(1)+i
      gj=ny*coords(2)+j
      gk=nz*coords(3)+k
      
      xblock=(i+2*TILE_DIMx_d-1)/TILE_DIMx_d
	  yblock=(j+2*TILE_DIMy_d-1)/TILE_DIMy_d
	  zblock=(k+2*TILE_DIMz_d-1)/TILE_DIMz_d
      
      myblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
      ii=i-xblock*TILE_DIMx_d+2*TILE_DIMx_d
      jj=j-yblock*TILE_DIMy_d+2*TILE_DIMy_d
      kk=k-zblock*TILE_DIMz_d+2*TILE_DIMz_d
      
      myarrx(li,lj,lk)=auxfields_s(idx5d(ii,jj,kk,5,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nauxfields))
      myarry(li,lj,lk)=auxfields_s(idx5d(ii,jj,kk,6,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nauxfields))
      myarrz(li,lj,lk)=auxfields_s(idx5d(ii,jj,kk,7,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nauxfields))
      
      call syncthreads
      
      if(abs(isfluid(i,j,k)) .ne. 1)return
      
      !idblock is the index of the block of internal nodes without the surrounding halo
	  idblock=blockIdx%x+blockIdx%y*nxblock_d+blockIdx%z*nxyblock_d+1
	  !If my block index does not match the index of the internal-node block, it means my thread is on the outer halo and must exit.
	  if(myblock .ne. idblock)return


	   
      !div_thetan here
      locauxfields_s(idx5d(ii,jj,kk,2,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nlocauxfields))= &
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
      
   endsubroutine compute_div_thetan_kernel
      
   subroutine thinfilm_scan_mark_cuda(phifields_s)

      implicit none
      real(kind=db), allocatable, dimension(:) :: phifields_s
 
   
#if defined(TWOCOMPONENT) && defined(REPULSIVE_FLUX)
      !$acc wait
      istat = cudaDeviceSynchronize


!$acc host_data use_device(flop,nx,ny,nz,coords,q_th,win,cosOppT,pwr,A_rep,isfluid,rep_mask &
       !$acc& ,ntotphifields,ntotauxfields,ntotlocauxfields,phifields_s,auxfields,locauxfields)
       call thinfilm_scan_mark_kernel<<<dimGrid, dimBlock>>>(flop,nx,ny,nz,coords,q_th,win,cosOppT,pwr,A_rep,isfluid, &
        rep_mask,ntotphifields,ntotauxfields,ntotlocauxfields,phifields_s,auxfields,locauxfields)
!$acc end host_data
      istat = cudaDeviceSynchronize
      istat = cudaGetLastError()
      if (istat/= cudaSuccess) then
        if(myrank==0)write(6,*) 'status after at ', __LINE__ ,' of file ', __FILE__ ,' :'
        if(myrank==0)write(6,*) cudaGetErrorString(istat)
        call doerror(6,'ERROR in thinfilm_scan_mark_cuda')
      endif
      !$acc wait        
#endif

   
      return
      
   endsubroutine thinfilm_scan_mark_cuda
   
   attributes(global) subroutine thinfilm_scan_mark_kernel(flop,nx,ny,nz,coords,q_th,win,cosOppT,pwr,A_rep,isfluid, &
    rep_mask,ntotphifields,ntotauxfields,ntotlocauxfields,phifields_s,auxfields_s,locauxfields_s)
      implicit none
      
      integer :: flop,nx,ny,nz
      integer, dimension(3) :: coords
      real(kind=db) :: q_th,win,cosOppT,pwr,A_rep	
      integer(kind=isf), dimension(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff) :: isfluid
      integer(kind=isf), dimension(1:nx,1:ny,1:nz) :: rep_mask
      integer :: ntotphifields,ntotauxfields,ntotlocauxfields
      real(kind=db), dimension(ntotphifields) :: phifields_s
      real(kind=db), dimension(ntotauxfields) :: auxfields_s
      real(kind=db), dimension(ntotlocauxfields) :: locauxfields_s
      
      integer :: i,j,k,gi,gj,gk,myblock,ii,jj,kk,iii,jjj,kkk
      
      integer :: di,dj,dk
	  integer :: diii,djjj,dkkk
	  real(kind=db) :: nix,niy,niz, dotn, qloc, qneig, face
	  real(kind=db) :: best_r2, r2, best_face
	  integer :: iii_best, jjj_best, kkk_best
	  logical :: found
	  real(kind=db), parameter :: eps = 1.0e-12_db
	  integer :: oii,ojj,okk
      integer :: oxblock,oyblock,ozblock,omyblock
      
      i = (blockIdx%x-1) * TILE_DIMx_d + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy_d + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz_d + threadIdx%z
      
      if(abs(isfluid(i,j,k)) .ne. 1)return
      
      gi=nx*coords(1)+i
      gj=ny*coords(2)+j
      gk=nz*coords(3)+k
      
      myblock=blockIdx%x+blockIdx%y*nxblock_d+blockIdx%z*nxyblock_d+1
      
      ii=threadIdx%x
      jj=threadIdx%y
      kk=threadIdx%z
      
      rep_mask(i,j,k) = 0
      locauxfields_s(idx5d(ii,jj,kk,3,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nlocauxfields)) = 0
      locauxfields_s(idx5d(ii,jj,kk,4,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nlocauxfields)) = 0
      locauxfields_s(idx5d(ii,jj,kk,5,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nlocauxfields)) = 0

	  ! gate: interfacial cell (use clamped phi for q)
      qloc = phifields_s(idx5d(ii,jj,kk,1,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nphifields))
      qloc = min(max(qloc,0.0_db),1.0_db)
      qloc = qloc*(1.0_db - qloc)
      if (qloc < q_th) return

      nix = auxfields_s(idx5d(ii,jj,kk,1,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nauxfields)) 
      niy = auxfields_s(idx5d(ii,jj,kk,2,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nauxfields))
      niz = auxfields_s(idx5d(ii,jj,kk,3,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nauxfields))

      best_r2   = HUGE(1.0_db)
      best_face = -1.0_db
      found     = .false.

      do di = -win, win
        do dj = -win, win
          do dk = -win, win
                  if (di==0 .and. dj==0 .and. dk==0) cycle

				  ! ---- 
				  iii = i + di
				  jjj = j + dj
				  kkk = k + dk
				  
				  if(abs(isfluid(iii,jjj,kkk)) .ne. 1)cycle

				  ! ---- minimum-image index differences
				  diii = iii - i
				 
				  djjj = jjj - j
				  
				  dkkk = kkk - k
				  
                  
				  r2 = real(diii,db)*real(diii,db) + real(djjj,db)*real(djjj,db) + real(dkkk,db)*real(dkkk,db)
				  if (r2 < eps) cycle

				  ! ---- neighbor interfacial gate (clamped) + similarity
				  
				  oxblock=(iii+2*TILE_DIMx_d-1)/TILE_DIMx_d   
                  oyblock=(jjj+2*TILE_DIMy_d-1)/TILE_DIMy_d     
                  ozblock=(kkk+2*TILE_DIMz_d-1)/TILE_DIMz_d 
                  omyblock=(oxblock-1)+(oyblock-1)*nxblock_d+(ozblock-1)*nxyblock_d+1
                  oii=iii-oxblock*TILE_DIMx_d+2*TILE_DIMx_d
                  ojj=jjj-oyblock*TILE_DIMy_d+2*TILE_DIMy_d
                  okk=kkk-ozblock*TILE_DIMz_d+2*TILE_DIMz_d
				  
				  qneig = phifields_s(idx5d(oii,ojj,okk,1,omyblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nphifields))
				  qneig = min(max(qneig,0.0_db),1.0_db)
				  qneig = qneig*(1.0_db - qneig)
				  if ( (qneig < q_th) .or. (abs(qneig - qloc) > 0.1_db*max(qloc,1.0e-12_db)) ) cycle

				  ! ---- facing condition (opposite normals): dotn <= cosOppT
				  dotn = nix*auxfields_s(idx5d(oii,ojj,okk,1,omyblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nauxfields)) &
				   + niy*auxfields_s(idx5d(oii,ojj,okk,2,omyblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nauxfields)) &
				   + niz*auxfields_s(idx5d(oii,ojj,okk,3,omyblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nauxfields)) 
				  if (dotn > cosOppT) cycle
				  face = 0.5_db*(1.0_db - dotn)   ! in [0,1]

				  ! ---- pick nearest; tie-break by larger 'face'
				  if (r2 < best_r2 - 1.0e-14_db) then
					best_r2 = r2; best_face = face
					iii_best = iii; jjj_best = jjj; kkk_best = kkk
					found   = .true.
				  else if (abs(r2 - best_r2) <= 1.0e-14_db) then
					if (face > best_face) then
					  best_face = face
					  iii_best = iii; jjj_best = jjj; kkk_best = kkk
					  found   = .true.
					end if
				  end if

				end do
			  end do
			end do

			if (found) then
			  locauxfields_s(idx5d(ii,jj,kk,3,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nlocauxfields)) = iii_best
			  locauxfields_s(idx5d(ii,jj,kk,4,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nlocauxfields)) = jjj_best
			  locauxfields_s(idx5d(ii,jj,kk,5,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nlocauxfields)) = kkk_best
			  rep_mask(i,j,k) = 1
			end if
   
      return
      
   endsubroutine thinfilm_scan_mark_kernel
   
   subroutine repulsive_flux_normal_cuda(phifields_s)

      implicit none
      real(kind=db), allocatable, dimension(:) :: phifields_s
 
   
#if defined(TWOCOMPONENT) && defined(REPULSIVE_FLUX)
      !$acc wait
      istat = cudaDeviceSynchronize


!$acc host_data use_device(flop,nx,ny,nz,coords,width,q_th,win,cosOppT,pwr,A_rep,isfluid,rep_mask &
       !$acc& ,ntotphifields,ntotauxfields,ntotlocauxfields,phifields_s,auxfields,locauxfields)
       call repulsive_flux_normal_kernel<<<dimGrid, dimBlock>>>(flop,nx,ny,nz,coords,width,q_th,win,cosOppT,pwr,A_rep,isfluid, &
        rep_mask,ntotphifields,ntotauxfields,ntotlocauxfields,phifields_s,auxfields,locauxfields)
!$acc end host_data
      istat = cudaDeviceSynchronize
      istat = cudaGetLastError()
      if (istat/= cudaSuccess) then
        if(myrank==0)write(6,*) 'status after at ', __LINE__ ,' of file ', __FILE__ ,' :'
        if(myrank==0)write(6,*) cudaGetErrorString(istat)
        call doerror(6,'ERROR in repulsive_flux_normal_cuda')
      endif
      !$acc wait        
#endif

   
      return
      
   endsubroutine repulsive_flux_normal_cuda
   
   attributes(global) subroutine repulsive_flux_normal_kernel(flop,nx,ny,nz,coords,width,q_th,win,cosOppT,pwr,A_rep,isfluid, &
    rep_mask,ntotphifields,ntotauxfields,ntotlocauxfields,phifields_s,auxfields_s,locauxfields_s)
      implicit none
      
      integer :: flop,nx,ny,nz
      integer, dimension(3) :: coords
      real(kind=db) :: width,q_th,win,cosOppT,pwr,A_rep	
      integer(kind=isf), dimension(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff) :: isfluid
      integer(kind=isf), dimension(1:nx,1:ny,1:nz) :: rep_mask
      integer :: ntotphifields,ntotauxfields,ntotlocauxfields
      real(kind=db), dimension(ntotphifields) :: phifields_s
      real(kind=db), dimension(ntotauxfields) :: auxfields_s
      real(kind=db), dimension(ntotlocauxfields) :: locauxfields_s
      
      integer :: i,j,k,gi,gj,gk,myblock,ii,jj,kk,iii,jjj,kkk
      
      real(kind=db) :: q1,q2,qpair,qcl,loc_phi,loc_phi2
	  real(kind=db) :: nx1,ny1,nz1, nx2,ny2,nz2
	  real(kind=db) :: dx,dy,dz, r, rinv, face, arg_arcosh, ach, Wfilm, wdth
	  real(kind=db) :: nsx,nsy,nsz, nsmag,alpha,cap,scales
	  real(kind=db), parameter :: eps = 1.0e-9_db
      
	  integer :: oii,ojj,okk
      integer :: oxblock,oyblock,ozblock,omyblock
      
      i = (blockIdx%x-1) * TILE_DIMx_d + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy_d + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz_d + threadIdx%z
      
      if(abs(isfluid(i,j,k)) .ne. 1)return
      
      gi=nx*coords(1)+i
      gj=ny*coords(2)+j
      gk=nz*coords(3)+k
      
      myblock=blockIdx%x+blockIdx%y*nxblock_d+blockIdx%z*nxyblock_d+1
      
      ii=threadIdx%x
      jj=threadIdx%y
      kk=threadIdx%z

	  locauxfields_s(idx5d(ii,jj,kk,6,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nlocauxfields))=0.0_db
	  locauxfields_s(idx5d(ii,jj,kk,7,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nlocauxfields))=0.0_db
	  locauxfields_s(idx5d(ii,jj,kk,8,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nlocauxfields))=0.0_db
	
	  if (rep_mask(i,j,k) .ne. 1) return
      
      loc_phi=phifields_s(idx5d(ii,jj,kk,1,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nphifields))
      
	  q1 = loc_phi*(1.0_db - loc_phi)
	
	  if (q1 <= eps) return

	  iii = int(locauxfields_s(idx5d(ii,jj,kk,3,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nlocauxfields)))
	  jjj = int(locauxfields_s(idx5d(ii,jj,kk,4,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nlocauxfields)))
	  kkk = int(locauxfields_s(idx5d(ii,jj,kk,5,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nlocauxfields)))

	  !line-of-centers
	  dx = real(iii - i,db)
	  dy = real(jjj - j,db)
	  dz = real(kkk - k,db)
	  r  = sqrt(dx*dx + dy*dy + dz*dz)
	  if (r <= eps) return
	  rinv = 1.0_db / r
	  dx = dx*rinv; dy = dy*rinv; dz = dz*rinv      ! u

	  !normals
	  nx1 = auxfields_s(idx5d(ii,jj,kk,1,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nauxfields)) 
      ny1 = auxfields_s(idx5d(ii,jj,kk,2,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nauxfields))
      nz1 = auxfields_s(idx5d(ii,jj,kk,3,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nauxfields))
	  
	  
	  oxblock=(iii+2*TILE_DIMx_d-1)/TILE_DIMx_d   
      oyblock=(jjj+2*TILE_DIMy_d-1)/TILE_DIMy_d     
      ozblock=(kkk+2*TILE_DIMz_d-1)/TILE_DIMz_d 
      omyblock=(oxblock-1)+(oyblock-1)*nxblock_d+(ozblock-1)*nxyblock_d+1
      oii=iii-oxblock*TILE_DIMx_d+2*TILE_DIMx_d
      ojj=jjj-oyblock*TILE_DIMy_d+2*TILE_DIMy_d
      okk=kkk-ozblock*TILE_DIMz_d+2*TILE_DIMz_d
                  
	  nx2 = auxfields_s(idx5d(oii,ojj,okk,1,omyblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nauxfields))
	  ny2 = auxfields_s(idx5d(oii,ojj,okk,2,omyblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nauxfields))
	  nz2 = auxfields_s(idx5d(oii,ojj,okk,3,omyblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nauxfields))

	  !facing factor in [0,1]
	  face = max( 0.0_db, -(nx1*nx2 + ny1*ny2 + nz1*nz2) )

	  if (face <= eps) return

	  !symmetric normal: bisector n1 - n2 (for facing sheets)
	  nsx = nx1 - nx2
	  nsy = ny1 - ny2
	  nsz = nz1 - nz2
	  nsmag = sqrt(nsx*nsx + nsy*nsy + nsz*nsz)
	  if (nsmag <= eps) then
	    !fallback to line-of-centers
	    nsx = dx; nsy = dy; nsz = dz
	  else
	    nsx = nsx / nsmag
	    nsy = nsy / nsmag
	    nsz = nsz / nsmag
	  end if

	  !orient so u·nsym >= 0  (partner will flip)
	  if (dx*nsx + dy*nsy + dz*nsz < 0.0_db) then
	    nsx = -nsx; nsy = -nsy; nsz = -nsz
	  end if

	  !symmetric magnitude from qpair
	  loc_phi2=phifields_s(idx5d(oii,ojj,okk,1,omyblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nphifields))
	  q2 = loc_phi2*(1.0_db - loc_phi2)
	  qpair = 0.5_db*(q1 + q2)
	  qcl   = min( max(qpair, eps), 0.25_db - eps )

	  arg_arcosh = 1.0_db / ( 2.0_db*sqrt(qcl) )
	  if (arg_arcosh <= 1.0_db) return

	  ach   = log( arg_arcosh + sqrt(arg_arcosh*arg_arcosh - 1.0_db) )
	  Wfilm = width * ach
	  if (Wfilm <= 0.0_db) return

	  wdth  = 1.0_db /(1.0 + wfilm**4.0) ! ( 1.0_db + (1.0_db/Wfilm)**4 )

	  !final purely-normal, symmetric repulsive flux
	  dx = A_rep * wdth * qcl * face * nsx
	  dy = A_rep * wdth * qcl * face * nsy
	  dz = A_rep * wdth * qcl * face * nsz
	
	  alpha = 1.5_db
	  cap   = alpha * (abs(dx)+abs(dy)+abs(dz)) 
	  scales = min(1.0_db, loc_phi / max(cap, 1.0e-9_db))
	  
	  locauxfields_s(idx5d(ii,jj,kk,6,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nlocauxfields)) = dx * scales
	  locauxfields_s(idx5d(ii,jj,kk,7,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nlocauxfields)) = dy * scales
	  locauxfields_s(idx5d(ii,jj,kk,8,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nlocauxfields)) = dz * scales
       
 end subroutine repulsive_flux_normal_kernel
 
 subroutine moments_LB_cuda(hfields_old,hfields_s &
#ifdef TWOCOMPONENT    
  ,phifields_s &
#endif 
 )

      implicit none
      
      real(kind=db), allocatable, dimension(:) :: hfields_old,hfields_s
#ifdef TWOCOMPONENT          
      real(kind=db), allocatable, dimension(:) :: phifields_s
#endif 

!      if(myrank==0)write(6,*)'step ',step, __LINE__ , __FILE__
      !$acc wait
      istat = cudaDeviceSynchronize

!$acc host_data use_device(step,iprobe,jprobe,kprobe,flop,nx,ny,nz,coords,isfluid & 
#ifdef MULTIHIT
       !$acc& ,ABCx,ABCy,ABCz &
#endif 
#ifdef WETTABILITY
       !$acc& ,wettab_r,wettab_b &
#endif  
#ifdef TWOCOMPONENT 
       !$acc& ,visc1,visc2,rho_r,rho_b,invrho_r,invrho_b,beta,kapp,sigma,sharp_c,tau_diff,phifields_s &
#endif   
#if defined(ELASTIC_FORCE)
       !$acc& ,lambda_rel,k_elastic,u_ref,v_ref,w_ref &
#endif
       !$acc& ,fx,fy,fz,ntothfields,ntotphifields,ntotauxfields,ntotlocauxfields,ntotforces &
       !$acc& ,hfields_old,hfields_s,auxfields,locauxfields,forces)
       call moments_LB_kernel<<<dimGrid, dimBlock>>>(step,iprobe,jprobe,kprobe,flop,nx,ny,nz,coords,isfluid &   
#ifdef MULTIHIT
	   ,ABCx,ABCy,ABCz &
#endif 
#ifdef WETTABILITY
       ,wettab_r,wettab_b &
#endif  
#ifdef TWOCOMPONENT 
       ,visc1,visc2,rho_r,rho_b,invrho_r,invrho_b,beta,kapp,sigma,sharp_c,tau_diff,phifields_s &       
#endif  
#if defined(ELASTIC_FORCE)
       ,lambda_rel,k_elastic,u_ref,v_ref,w_ref &
#endif
	   ,fx,fy,fz,ntothfields,ntotphifields,ntotauxfields,ntotlocauxfields,ntotforces &
	   ,hfields_old,hfields_s,auxfields,locauxfields,forces)
!$acc end host_data
      istat = cudaDeviceSynchronize
      istat = cudaGetLastError()
      if (istat/= cudaSuccess) then
        if(myrank==0)write(6,*) 'status after at ', __LINE__ ,' of file ', __FILE__ ,' :'
        if(myrank==0)write(6,*) cudaGetErrorString(istat)
        call doerror(6,'ERROR in moments_LB_cuda')
      endif
      !$acc wait        
      
   endsubroutine moments_LB_cuda
 
 attributes(global) subroutine moments_LB_kernel(step,iprobe,jprobe,kprobe,flop,nx,ny,nz,coords,isfluid &   
#ifdef MULTIHIT
	   ,ABCx,ABCy,ABCz &
#endif 
#ifdef WETTABILITY
       ,wettab_r,wettab_b &
#endif  
#ifdef TWOCOMPONENT 
       ,visc1,visc2,rho_r,rho_b,invrho_r,invrho_b,beta,kapp,sigma,sharp_c,tau_diff,phifields_s &
#endif   
#if defined(ELASTIC_FORCE)
       ,lambda_rel,k_elastic,u_ref,v_ref,w_ref &
#endif
       ,fx,fy,fz,ntothfields,ntotphifields,ntotauxfields,ntotlocauxfields,ntotforces &
       ,hfields_old,hfields_s,auxfields_s,locauxfields_s,forces_s)
 

      implicit none
      
      integer :: step,iprobe,jprobe,kprobe,flop,nx,ny,nz
      
      integer, dimension(3) :: coords
      integer(kind=isf), dimension(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff) :: isfluid
      integer :: ntothfields,ntotphifields,ntotauxfields,ntotlocauxfields,ntotforces
#ifdef MULTIHIT
	  real(kind=db), dimension(1:nx,1:ny,1:nz) :: ABCx,ABCy,ABCz
#endif
#ifdef WETTABILITY  
      real(kind=db) :: wettab_r,wettab_b  
#endif  
#ifdef TWOCOMPONENT 
      real(kind=db) :: visc1,visc2,rho_r,rho_b,invrho_r,invrho_b,beta,kapp,sigma,sharp_c,tau_diff    
      real(kind=db), dimension(ntotphifields) :: phifields_s 
#endif           
#if defined(ELASTIC_FORCE)
      real(kind=db) :: lambda_rel,k_elastic
      real(kind=db), dimension(1:nx,1:ny,1:nz) :: u_ref,v_ref,w_ref
#endif
      real(kind=db) :: fx,fy,fz

      real(kind=db), dimension(ntothfields) :: hfields_old,hfields_s
      real(kind=db), dimension(ntotauxfields) :: auxfields_s
      real(kind=db), dimension(ntotlocauxfields) :: locauxfields_s
      real(kind=db), dimension(ntotforces) :: forces_s
      

      real(kind=db) :: pxx,pyy,pzz,pxy,pxz,pyz
      real(kind=db) :: mytemp,forcex,forcey,forcez,rhophi_loc,press_loc,u_loc,v_loc,w_loc
#ifdef EXPLICITEQ
      integer :: l
      real(kind=db) :: udotc,uu,F_discr,feq,fneq1,fpost
#endif
#ifdef TWOCOMPONENT
	  real(kind=db) ::gradfix,gradfiy,gradfiz,wet_loc,phi_loc,lap_phi_loc
#endif
      
#ifdef DENSRATIO
      real(kind=db) :: gradrhox,gradrhoy,gradrhoz,omega_loc, visc_loc,tau_loc
#endif
      
      integer :: i,j,k
      integer :: gi,gj,gk
      integer :: gif,gjf,gkf
      integer :: myblock,ii,jj,kk
      
      i = (blockIdx%x-1) * TILE_DIMx_d + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy_d + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz_d + threadIdx%z
      
      if(abs(isfluid(i,j,k)) .ne. 1)return
      
      gi=nx*coords(1)+i
      gj=ny*coords(2)+j
      gk=nz*coords(3)+k
      
      myblock=blockIdx%x+blockIdx%y*nxblock_d+blockIdx%z*nxyblock_d+1
      
      ii=threadIdx%x
      jj=threadIdx%y
      kk=threadIdx%z
      

				 
		  press_loc=hfields_s(idx5d(ii,jj,kk,1,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nhfields))
				 
#ifdef TWOCOMPONENT					 
		  phi_loc=phifields_s(idx5d(ii,jj,kk,1,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nphifields))
		  lap_phi_loc=locauxfields_s(idx5d(ii,jj,kk,1,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nlocauxfields))
#endif
#ifdef DENSRATIO
		  rhophi_loc = rho_r*phi_loc+(ONE-phi_loc)*rho_b 
#ifdef BUOYANCY_FORCING   
		  forcex=(rhophi_loc-(rho_r+rho_b)*HALF)*fx !0.0_db
		  forcey=(rhophi_loc-(rho_r+rho_b)*HALF)*fy !0.0_db
		  forcez=(rhophi_loc-(rho_r+rho_b)*HALF)*fz !0.0_db	
#else
		  forcex=rhophi_loc*fx !(rhophi_loc-(rho_r+rho_b)*HALF)*fx !0.0_db
		  forcey=rhophi_loc*fy !(rhophi_loc-(rho_r+rho_b)*HALF)*fy !0.0_db
		  forcez=rhophi_loc*fz !(rhophi_loc-(rho_r+rho_b)*HALF)*fz !0.0_db
#endif
#else
		  rhophi_loc = 1.0_db !press_loc
		  forcex=fx !0.0_db
		  forcey=fy !0.0_db
		  forcez=fz !0.0_db
#endif	
			 
#ifdef TWOCOMPONENT		
		   
		  !jaqmin 
		  mytemp=auxfields_s(idx5d(ii,jj,kk,4,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nauxfields)) !modgrad
		  gradfix=auxfields_s(idx5d(ii,jj,kk,1,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nauxfields))*mytemp !normx*modgrad
		  gradfiy=auxfields_s(idx5d(ii,jj,kk,2,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nauxfields))*mytemp !normy*modgrad
		  gradfiz=auxfields_s(idx5d(ii,jj,kk,3,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nauxfields))*mytemp !normz*modgrad
		  forcex = forcex + &
                   (4.0_db*beta*phi_loc*(phi_loc-1.0_db)*(phi_loc-0.5_db) - kapp*lap_phi_loc)*gradfix
		  forcey = forcey + &
                   (4.0_db*beta*phi_loc*(phi_loc-1.0_db)*(phi_loc-0.5_db) - kapp*lap_phi_loc)*gradfiy
		  forcez = forcez + &
                   (4.0_db*beta*phi_loc*(phi_loc-1.0_db)*(phi_loc-0.5_db) - kapp*lap_phi_loc)*gradfiz
				   				   

#ifdef REPULSIVE_FLUX
		  mytemp=locauxfields_s(idx5d(ii,jj,kk,6,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nlocauxfields))*rhophi_loc 
		  if(abs(mytemp)>1.0d-3) mytemp=1.0d-3*sign(1.0,mytemp)!mytemp*0.1_db
		  forcex=forcex + mytemp*rhophi_loc
				  
		  mytemp=locauxfields_s(idx5d(ii,jj,kk,7,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nlocauxfields))*rhophi_loc 
		  if(abs(mytemp)>1.0d-3) mytemp=1.0d-3*sign(1.0,mytemp)!mytemp*0.1_db
		  forcey=forcey + mytemp*rhophi_loc
				  
		  mytemp=locauxfields_s(idx5d(ii,jj,kk,8,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nlocauxfields))*rhophi_loc 
		  if(abs(mytemp)>1.0d-3) mytemp=1.0d-3*sign(1.0,mytemp)!mytemp*0.1_db
		  forcez=forcez + mytemp*rhophi_loc
#endif
#endif


#if defined(MULTIHIT)   
                  forcex=forcex + rhophi_loc*ABCx(i,j,k) !+ AAA*sin(k_zero*gk) + AAA*sin(k_zero*gj)  
		  forcey=forcey + rhophi_loc*ABCy(i,j,k) !+ AAA*sin(k_zero*gi) + AAA*sin(k_zero*gk)
		  forcez=forcez + rhophi_loc*ABCz(i,j,k) !+ AAA*sin(k_zero*gj) + AAA*sin(k_zero*gi) 	
#endif
                  
#ifdef DENSRATIO				  
		  ! pressure and viscous forces
				  
		  gradrhox=(rho_r-rho_b)*gradfix
		  gradrhoy=(rho_r-rho_b)*gradfiy
		  gradrhoz=(rho_r-rho_b)*gradfiz
				  
                  forcex=forcex - &
                   press_loc*cssq*gradrhox   !+ (rhophi_loc-(rho_r+rho_b)*HALF)*fx
                  forcey=forcey - &
                   press_loc*cssq*gradrhoy   !+ (rhophi_loc-(rho_r+rho_b)*HALF)*fy
		  forcez=forcez - &
		   press_loc*cssq*gradrhoz   !+ (rhophi_loc-(rho_r+rho_b)*HALF)*fz
		  !! from this point I compute the force terms that depend on the velocity
		  !! these terms should be not included in force arrays since they must be computed with the updated velocity
		  !! at the end of this subroutine
#endif
                  
                  gif=idx5d(ii,jj,kk,1,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nforces)
                  gjf=idx5d(ii,jj,kk,2,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nforces)
                  gkf=idx5d(ii,jj,kk,3,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nforces)
                  
		  forces_s(gif)=forcex
		  forces_s(gjf)=forcey
		  forces_s(gkf)=forcez
				  
		  u_loc=hfields_old(idx5d(ii,jj,kk,2,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nhfields)) !velocity
                  v_loc=hfields_old(idx5d(ii,jj,kk,3,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nhfields))
                  w_loc=hfields_old(idx5d(ii,jj,kk,4,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nhfields))
                  
                 
				  
#if defined(ELASTIC_FORCE)
                  forcex=forcex + &
                   rhophi_loc*(u_loc - u_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+   
		  forcey=forcey + &
		   rhophi_loc*(v_loc - v_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+ 
		  forcez=forcez + &
		   rhophi_loc*(w_loc - w_ref(i,j,k))*k_elastic*lambda_rel*phi_loc + rhophi_loc*fz !+  	
#endif

#ifdef DENSRATIO 
			  
                  pxx=hfields_s(idx5d(ii,jj,kk,5,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nhfields))
                  pyy=hfields_s(idx5d(ii,jj,kk,6,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nhfields))
                  pzz=hfields_s(idx5d(ii,jj,kk,7,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nhfields))
                  pxy=hfields_s(idx5d(ii,jj,kk,8,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nhfields))
                  pxz=hfields_s(idx5d(ii,jj,kk,9,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nhfields))
                  pyz=hfields_s(idx5d(ii,jj,kk,10,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nhfields))
                  
                  !1-2
                  !*1
                  ! 2nd order
#ifdef EXPLICITEQ                  
		  uu=HALF*(u_loc*u_loc+v_loc*v_loc+w_loc*w_loc)*invcssq
				  
                  do l=1,nlinks
                     udotc=(u_loc*dex(l) + v_loc*dey(l)+ w_loc*dez(l))*invcssq
		     feq=p(l)*(press_loc + (udotc+0.5_db*udotc*udotc - uu))
                     !fneq1=(HALF/(cssq*cssq))*( (dex(l)*dex(l)-cssq)*pxx &
		     ! + (dey(l)*dey(l)-cssq)*pyy + (dez(l)*dez(l)-cssq)*pzz &
	             ! + TWO*(dex(l)*dey(l))*pxy + TWO*(dex(l)*dez(l))*pxz &
		     ! + TWO*(dey(l)*dez(l))*pyz)
                    fpost=feq!+fneq1
                    pxx=pxx - fpost*(dex(l)*dex(l))
                    pyy=pyy - fpost*(dey(l)*dey(l))
                    pzz=pzz - fpost*(dez(l)*dez(l))
                    pxy=pxy - fpost*(dex(l)*dey(l))
                    pxz=pxz - fpost*(dex(l)*dez(l))
                    pyz=pyz - fpost*(dey(l)*dez(l))
                  enddo
#else
                  pxx=pxx - cssq*press_loc - u_loc*u_loc
                  pyy=pyy - cssq*press_loc - v_loc*v_loc
                  pzz=pzz - cssq*press_loc - w_loc*w_loc
                  pxy=pxy - u_loc*v_loc
                  pxz=pxz - u_loc*w_loc
                  pyz=pyz - v_loc*w_loc

#endif
		  visc_loc=(rho_r*visc1*phi_loc+(1.0_db-phi_loc)*visc2*rho_b)/rhophi_loc
				  
                  tau_loc=(visc_loc/cssq + HALF) !è una tau
				  
		  forcex=forcex - (visc_loc/(tau_loc*cssq))*(pxx*gradrhox + pxy*gradrhoy + pxz*gradrhoz)
		  forcey=forcey - (visc_loc/(tau_loc*cssq))*(pyy*gradrhoy + pxy*gradrhox + pyz*gradrhoz)
		  forcez=forcez - (visc_loc/(tau_loc*cssq))*(pzz*gradrhoz + pxz*gradrhox + pyz*gradrhoy)
#endif	
                  !I compute the new velocities
		  u_loc=hfields_s(idx5d(ii,jj,kk,2,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nhfields)) !velocity
                  v_loc=hfields_s(idx5d(ii,jj,kk,3,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nhfields))
                  w_loc=hfields_s(idx5d(ii,jj,kk,4,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nhfields))            
					 
#if defined(INTERFACE_INCOMP) && defined(DENSRATIO)
		  mytemp= -sharp_c*locauxfields_s(idx5d(ii,jj,kk,2,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nlocauxfields))
		  u_loc = u_loc + 0.5_db*forcex/(rhophi_loc)
                  v_loc = v_loc + 0.5_db*forcey/(rhophi_loc)
                  w_loc = w_loc + 0.5_db*forcez/(rhophi_loc)
				  
		  u_loc = u_loc/ &
		   (1.0_db - 0.5_db*(rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)/rhophi_loc )
				  
		  v_loc = v_loc/ &
		   (1.0_db - 0.5_db*(rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)/rhophi_loc ) 
				  
		  w_loc = w_loc/ &
	           (1.0_db - 0.5_db*(rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)/rhophi_loc ) 

		  forces_s(gif)= forces_s(gif) - &
		   (rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)*u_loc
		  forces_s(gjf)= forces_s(gjf) - &
		   (rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)*v_loc
		  forces_s(gkf)= forces_s(gkf) - &
		   (rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)*w_loc
					 
		  forcex=forcex - &
		   (rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)*u_loc
					 
		  forcey=forcey - &
		   (rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)*v_loc
					 
		  forcez=forcez - &
		   (rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)*w_loc
                  
#else
                  u_loc = u_loc + 0.5_db*forcex/rhophi_loc
                  v_loc = v_loc + 0.5_db*forcey/rhophi_loc
                  w_loc = w_loc + 0.5_db*forcez/rhophi_loc
#endif
                  

                  hfields_s(idx5d(ii,jj,kk,2,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nhfields))=u_loc   !put the new velocity in hfields_s
                  hfields_s(idx5d(ii,jj,kk,3,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nhfields))=v_loc
                  hfields_s(idx5d(ii,jj,kk,4,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nhfields))=w_loc
                  
                                    

!regularized 
		  pxx=hfields_s(idx5d(ii,jj,kk,5,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nhfields))
                  pyy=hfields_s(idx5d(ii,jj,kk,6,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nhfields))
                  pzz=hfields_s(idx5d(ii,jj,kk,7,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nhfields))
                  pxy=hfields_s(idx5d(ii,jj,kk,8,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nhfields))
                  pxz=hfields_s(idx5d(ii,jj,kk,9,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nhfields))
                  pyz=hfields_s(idx5d(ii,jj,kk,10,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nhfields))
                  
                  
#ifdef EXPLICITEQ 
                  do l=1,nlinks
                     
		     F_discr= &
		     p(l)*( (dex(l)-u_loc)*(forcex)+(dey(l)-v_loc)*(forcey)+(dez(l)-w_loc)*(forcez) + &
		      (1.0_db/(cssq))*( (u_loc*dex(l)+v_loc*dey(l)+w_loc*dez(l))*&
		      ( (forcex)*dex(l) + (forcey)*dey(l) + (forcez)*dez(l) ) ) )/(cssq*rhophi_loc) 
					 
                     udotc = 0.5_db*F_discr 
                     
                     pxx=pxx + udotc*dex(l)*dex(l)
                     pyy=pyy + udotc*dey(l)*dey(l)
                     pzz=pzz + udotc*dez(l)*dez(l)
                     pxy=pxy + udotc*dex(l)*dey(l)
                     pxz=pxz + udotc*dex(l)*dez(l)
                     pyz=pyz + udotc*dey(l)*dez(l)
                  enddo
#else
                  
                  pxx=pxx + forcex*u_loc/rhophi_loc
                  pyy=pyy + forcey*v_loc/rhophi_loc
                  pzz=pzz + forcez*w_loc/rhophi_loc
                  pxy=pxy + HALF*(forcey*u_loc+forcex*v_loc)/rhophi_loc
                  pxz=pxz + HALF*(forcez*u_loc+forcex*w_loc)/rhophi_loc
                  pyz=pyz + HALF*(forcez*v_loc+forcey*w_loc)/rhophi_loc

#endif
                  hfields_s(idx5d(ii,jj,kk,5,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nhfields))=pxx
                  hfields_s(idx5d(ii,jj,kk,6,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nhfields))=pyy
                  hfields_s(idx5d(ii,jj,kk,7,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nhfields))=pzz
                  hfields_s(idx5d(ii,jj,kk,8,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nhfields))=pxy
                  hfields_s(idx5d(ii,jj,kk,9,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nhfields))=pxz
                  hfields_s(idx5d(ii,jj,kk,10,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nhfields))=pyz
                  
                  
#if defined(ELASTIC_FORCE)
		  u_ref(i,j,k) = u_ref(i,j,k) + &
		   lambda_rel*(u_loc - u_ref(i,j,k))
		  v_ref(i,j,k) = v_ref(i,j,k) + &
		   lambda_rel*(v_loc - v_ref(i,j,k))
		  w_ref(i,j,k) = w_ref(i,j,k) + &
		   lambda_rel*(w_loc - w_ref(i,j,k))
				  
                  forces_s(gif)= forces_s(gif) + &
                   rhophi_loc*(u_loc - u_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+ 
		  forces_s(gjf)= forces_s(gjf) +&
		   rhophi_loc*(v_loc - v_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+ 
		  forces_s(gkf)= forces_s(gkf) + &
		   rhophi_loc*(w_loc - w_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+rhophi_loc*fz  
#endif               
                  
#if defined(DENSRATIO)	
                  pxx=pxx - cssq*press_loc - u_loc*u_loc 
                  pyy=pyy - cssq*press_loc - v_loc*v_loc
                  pzz=pzz - cssq*press_loc - w_loc*w_loc 
                  pxy=pxy - u_loc*v_loc
                  pxz=pxz - u_loc*w_loc
                  pyz=pyz - v_loc*w_loc
		  
		  forces_s(gif)= forces_s(gif) - &
		   (visc_loc/(tau_loc*cssq))*(pxx*gradrhox + pxy*gradrhoy + pxz*gradrhoz)
		  forces_s(gjf)= forces_s(gjf) - &
		   (visc_loc/(tau_loc*cssq))*(pyy*gradrhoy + pxy*gradrhox + pyz*gradrhoz)
		  forces_s(gkf)= forces_s(gkf) - &
		   (visc_loc/(tau_loc*cssq))*(pzz*gradrhoz + pxz*gradrhox + pyz*gradrhoy)
#endif               


   endsubroutine moments_LB_kernel  

   subroutine fused_LB_cuda(hfields_in,hfields_out &
#ifdef TWOCOMPONENT    
   ,phifields_s &
#endif   
   )

      implicit none
      
      real(kind=db), allocatable, dimension(:) :: hfields_in,hfields_out
#ifdef TWOCOMPONENT       
      real(kind=db), allocatable, dimension(:) :: phifields_s
#endif
      integer(8) :: mymshared
      
      mymshared = 4 * (TILE_DIMx + 2) * (TILE_DIMy + 2) * (TILE_DIMz + 2) * db
      
!      if(myrank==0)write(6,*)'step ',step, __LINE__ , __FILE__
      !$acc wait
      istat = cudaDeviceSynchronize

!$acc host_data use_device(step,flip,flop,nx,ny,nz,coords,isfluid & 
#ifdef MULTIHIT
	   !$acc& ,ABCx,ABCy,ABCz &
#endif 
#ifdef WETTABILITY
       !$acc& ,wettab_r,wettab_b &
#endif  
#ifdef TWOCOMPONENT 
       !$acc& ,visc2,rho_r,rho_b,invrho_r,invrho_b,sharp_c,beta,kapp,tau_diff,sigma,phifields_s &
#ifdef MONOD
	   !$acc& ,mu_max,Ks &
#endif
#endif   
       !$acc& ,visc1,omega,fx,fy,fz,ntothfields,ntotphifields,ntotauxfields,ntotlocauxfields,ntotforces &
       !$acc& ,hfields_in,hfields_out,auxfields,locauxfields,forces)
      call fused_LB_kernel<<<dimGrid,dimBlockshared,mymshared>>>(step,flip,flop,nx,ny,nz,coords,isfluid &    
#ifdef MULTIHIT
	   ,ABCx,ABCy,ABCz &
#endif 
#ifdef WETTABILITY
       ,wettab_r,wettab_b &
#endif  
#ifdef TWOCOMPONENT 
       ,visc2,rho_r,rho_b,invrho_r,invrho_b,sharp_c,beta,kapp,tau_diff,sigma,phifields_s &
#ifdef MONOD
	   ,mu_max,Ks &
#endif
#endif   
       ,visc1,omega,fx,fy,fz,ntothfields,ntotphifields,ntotauxfields,ntotlocauxfields,ntotforces &
	   ,hfields_in,hfields_out,auxfields,locauxfields,forces)
!$acc end host_data
      
      istat = cudaGetLastError()         ! oppure cudaPeekAtLastError
      if (istat /= cudaSuccess) then
        if(myrank==0) then
          write(6,*) 'launch error at ', __LINE__, ' file ', __FILE__
          write(6,*) cudaGetErrorString(istat)
        endif
        call doerror(6,'ERROR in fused_LB_cuda (launch)')
      endif

      istat = cudaDeviceSynchronize()
      if (istat /= cudaSuccess) then
       if(myrank==0) then
         write(6,*) 'sync error at ', __LINE__, ' file ', __FILE__
         write(6,*) cudaGetErrorString(istat)
       endif
       call doerror(6,'ERROR in fused_LB_cuda (sync)')
      endif
      
      !$acc wait

      
   end subroutine fused_LB_cuda

   
   subroutine update_phifields(hfields_s,phifields_in,phifields_out)

      implicit none
      real(kind=db), allocatable, dimension(:) :: hfields_s,phifields_in,phifields_out
 
      
#ifdef TWOCOMPONENT
      !$acc wait
      istat = cudaDeviceSynchronize

!$acc host_data use_device(flop,nx,ny,nz,coords,isfluid &
#ifdef TWOCOMPONENT    
       !$acc& ,visc1,visc2,rho_r,rho_b,invrho_r,invrho_b,sharp_c,beta,kapp,tau_diff,sigma &
#ifdef MONOD
	   !$acc& ,mu_max,Ks &
#endif
#endif
       !$acc& ,ntothfields,ntotphifields,ntotauxfields,ntotlocauxfields,ntotforces &
	   !$acc& ,hfields_s,phifields_in,phifields_out,auxfields,locauxfields,forces)
       call update_phifields_kernel<<<dimGrid, dimBlock>>>(flop,nx,ny,nz,coords,isfluid &
#ifdef TWOCOMPONENT           
       ,visc2,rho_r,rho_b,invrho_r,invrho_b,sharp_c,beta,kapp,tau_diff,sigma &
#ifdef MONOD
	   ,mu_max,Ks &
#endif   
#endif   
       ,ntothfields,ntotphifields,ntotauxfields,ntotlocauxfields,ntotforces &
	   ,hfields_s,phifields_in,phifields_out,auxfields,locauxfields,forces)
!$acc end host_data
      istat = cudaDeviceSynchronize
      istat = cudaGetLastError()
      if (istat/= cudaSuccess) then
        if(myrank==0)write(6,*) 'status after at ', __LINE__ ,' of file ', __FILE__ ,' :'
        if(myrank==0)write(6,*) cudaGetErrorString(istat)
        call doerror(6,'ERROR in moments_LB_cuda')
      endif
      !$acc wait        
#endif
      
      return
      
   endsubroutine update_phifields 
   
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
      real(kind=db), dimension(ntothfields) :: hfields_s
      real(kind=db), dimension(ntotphifields) :: phifields_in,phifields_out
      real(kind=db), dimension(ntotauxfields) :: auxfields_s
      real(kind=db), dimension(ntotlocauxfields) :: locauxfields_s
      real(kind=db), dimension(ntotforces) :: forces_s
      
      integer :: i,j,k
      !integer :: gi,gj,gk
      integer :: myblock,ii,jj,kk
      real(kind=db) :: gradfix,gradfiy,gradfiz,phi_loc,phi_out,mytemp
      real(kind=db) :: loc_u,loc_v,loc_w,lap_phi_loc
#ifdef MONOD
	  real(kind=db) :: S_mono
#endif
      
      i = (blockIdx%x-1) * TILE_DIMx_d + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy_d + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz_d + threadIdx%z
      
      if(abs(isfluid(i,j,k)) .ne. 1)return
      
      !gi=nx*coords(1)+i
      !gj=ny*coords(2)+j
      !gk=nz*coords(3)+k
      
      myblock=blockIdx%x+blockIdx%y*nxblock_d+blockIdx%z*nxyblock_d+1
      
      ii=threadIdx%x
      jj=threadIdx%y
      kk=threadIdx%z
	   
#ifdef TWOCOMPONENT
				  
	  mytemp= -sharp_c*locauxfields_s(idx5d(ii,jj,kk,2,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nlocauxfields))
	  			  
	  !reuse gradrhox,gradrhoy,gradrhoz as local velocity (reusing variables is saving register memory)
	  !reuse gradfix,gradfiy,gradfiz
      modgrad=auxfields_s(idx5d(ii,jj,kk,4,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nauxfields)) !modgrad
	  gradfix=auxfields_s(idx5d(ii,jj,kk,1,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nauxfields))*modgrad !normx*modgrad
	  gradfiy=auxfields_s(idx5d(ii,jj,kk,2,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nauxfields))*modgrad !normy*modgrad
	  gradfiz=auxfields_s(idx5d(ii,jj,kk,3,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nauxfields))*modgrad !normz*modgrad
                  
      phi_loc=phifields_in(idx5d(ii,jj,kk,1,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nphifields))
      lap_phi_loc=locauxfields_s(idx5d(ii,jj,kk,1,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nlocauxfields))
                  
      loc_u=hfields_s(idx5d(ii,jj,kk,2,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nhfields)) !velocity
      loc_v=hfields_s(idx5d(ii,jj,kk,3,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nhfields))
      loc_w=hfields_s(idx5d(ii,jj,kk,4,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nhfields))
                  
      phi_out = phi_loc &
        - loc_u*0.5_db*(gradfix) - loc_v*0.5_db*(gradfiy) &
        - loc_w*0.5_db*(gradfiz) + tau_diff*lap_phi_loc + mytemp 
#endif	

#ifdef MONOD
      S_mono = mu_max * phi_loc)/(Ks + phi_loc) * phi_loc * (1.0_db - phi_loc)
      phi_out=phi_out + S_mono
	  !phi_out = min(1.0_db, max(0.0_db, phi_out))		 
#endif
      phifields_out(idx5d(ii,jj,kk,1,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nphifields))=phi_out
      
      return
      
   endsubroutine update_phifields_kernel
   
   subroutine LB_int_boundary_cuda(hfields_in,hfields_out &
#ifdef TWOCOMPONENT	       
     ,phifields_s &
#endif     
     )

      implicit none
      real(kind=db), allocatable, dimension(:) :: hfields_in,hfields_out
#ifdef TWOCOMPONENT	       
      real(kind=db), allocatable, dimension(:) :: phifields_s
#endif     
      
!      if(myrank==0)write(6,*)'step ',step, __LINE__ , __FILE__
      !$acc wait
      istat = cudaDeviceSynchronize

!$acc host_data use_device(flip,flop,nx,ny,nz,coords,isfluid & 
#ifdef MULTIHIT
	   !$acc& ,ABCx,ABCy,ABCz &
#endif 
#ifdef WETTABILITY
       !$acc& ,wettab_r,wettab_b &
#endif  
#ifdef TWOCOMPONENT 
       !$acc& ,visc2,rho_r,rho_b,invrho_r,invrho_b,sharp_c,beta,kapp,tau_diff,sigma,phifields_s &
#ifdef MONOD
	   !$acc& ,mu_max,Ks &
#endif
#endif   
       !$acc& ,visc1,omega,fx,fy,fz,ntothfields,ntotphifields,ntotauxfields,ntotlocauxfields,ntotforces &
       !$acc& ,hfields_in,hfields_out,auxfields,locauxfields,forces)
      call LB_int_boundary_kernel<<<dimGrid,dimBlock>>>(flip,flop,nx,ny,nz,coords,isfluid &    
#ifdef MULTIHIT
	   ,ABCx,ABCy,ABCz &
#endif 
#ifdef WETTABILITY
       ,wettab_r,wettab_b &
#endif  
#ifdef TWOCOMPONENT 
       ,visc2,rho_r,rho_b,invrho_r,invrho_b,sharp_c,beta,kapp,tau_diff,sigma,phifields_s &
#ifdef MONOD
	   ,mu_max,Ks &
#endif
#endif   
       ,visc1,omega,fx,fy,fz,ntothfields,ntotphifields,ntotauxfields,ntotlocauxfields,ntotforces &
	   ,hfields_in,hfields_out,auxfields,locauxfields,forces)
!$acc end host_data
      istat = cudaDeviceSynchronize
      istat = cudaGetLastError()
      if (istat/= cudaSuccess) then
        if(myrank==0)write(6,*) 'status after at ', __LINE__ ,' of file ', __FILE__ ,' :'
        if(myrank==0)write(6,*) cudaGetErrorString(istat)
        call doerror(6,'ERROR in LB_int_boundary_cuda')
      endif
      !$acc wait

      
   end subroutine LB_int_boundary_cuda

   attributes(global) subroutine LB_int_boundary_kernel(flip,flop,nx,ny,nz,coords,isfluid &  
#ifdef MULTIHIT
       ,ABCx,ABCy,ABCz &
#endif 
#ifdef WETTABILITY
       ,wettab_r,wettab_b &
#endif  
#ifdef TWOCOMPONENT 
       ,visc2,rho_r,rho_b,invrho_r,invrho_b,sharp_c,beta,kapp,tau_diff,sigma,phifields_s &
#ifdef MONOD	
       ,mu_max,Ks &
#endif
#endif   
       ,visc1,omega,fx,fy,fz,ntothfields,ntotphifields,ntotauxfields,ntotlocauxfields,ntotforces &
       ,hfields_in,hfields_out,auxfields_s,locauxfields_s,forces_s)

      implicit none
      
      integer :: flip,flop,nx,ny,nz
      integer, dimension(3) :: coords
      integer(kind=isf), dimension(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff) :: isfluid
      integer :: ntothfields,ntotphifields,ntotauxfields,ntotlocauxfields,ntotforces
#ifdef MULTIHIT
	  real(kind=db), dimension(1:nx,1:ny,1:nz) :: ABCx,ABCy,ABCz
#endif
#ifdef WETTABILITY  
      real(kind=db) :: wettab_r,wettab_b  
#endif  
#ifdef TWOCOMPONENT 
      real(kind=db) :: visc2,rho_r,rho_b,invrho_r,invrho_b,sharp_c,beta,kapp,tau_diff,sigma  
      real(kind=db), dimension(ntotphifields) :: phifields_s
#ifdef MONOD
      real(kind=db) :: mu_max,Ks
#endif
#endif           
      real(kind=db) :: visc1,omega,fx,fy,fz
      
      real(kind=db), dimension(ntothfields) :: hfields_in,hfields_out

      real(kind=db), dimension(ntotauxfields) :: auxfields_s
      real(kind=db), dimension(ntotlocauxfields) :: locauxfields_s
      real(kind=db), dimension(ntotforces) :: forces_s
      

      real(kind=db) :: press,u,v,w,pxx,pyy,pzz,pxy,pxz,pyz
      real(kind=db) :: opress,ou,ov,ow,opxx,opyy,opzz,opxy,opxz,opyz
      real(kind=db) :: mytemp,rhophi_loc,press_loc
      !real(kind=db) :: forcex,forcey,forcez,F_discr
  
      real(kind=db) :: fneq1,feq,fpost,uu,udotc
#ifdef TWOCOMPONENT
      real(kind=db) :: wet_loc
#endif
      real(kind=db) :: omega_loc,phi_loc,visc_loc
#ifdef SMAGORINSKI
	  real(kind=db) :: QQ
#endif

      integer :: i,j,k,l,lopp
      integer :: gi,gj,gk
      integer :: myblock,ii,jj,kk
      integer :: iii,jjj,kkk
      
      i = (blockIdx%x-1) * TILE_DIMx_d + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy_d + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz_d + threadIdx%z
      
      if(isfluid(i,j,k) .ne. -1)return
      
      gi=nx*coords(1)+i
      gj=ny*coords(2)+j
      gk=nz*coords(3)+k
      
      myblock=blockIdx%x+blockIdx%y*nxblock_d+blockIdx%z*nxyblock_d+1
      
      ii=threadIdx%x
      jj=threadIdx%y
      kk=threadIdx%z
       
#ifdef TWOCOMPONENT	       
	  phi_loc=phifields_s(idx5d(ii,jj,kk,1,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nphifields))
#endif	
#ifdef DENSRATIO
	  rhophi_loc = rho_r*phi_loc+(ONE-phi_loc)*rho_b 
#else
	  rhophi_loc = 1.0_db !press_loc
#endif	

!	  forcex=force_s(idx5d(ii,jj,kk,1,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nforces))
!	  forcey=force_s(idx5d(ii,jj,kk,2,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nforces))
!	  forcez=force_s(idx5d(ii,jj,kk,3,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nforces))


	  press=hfields_in(idx5d(ii,jj,kk,1,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nhfields))
	  u=hfields_in(idx5d(ii,jj,kk,2,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nhfields)) 
	  v=hfields_in(idx5d(ii,jj,kk,3,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nhfields))
	  w=hfields_in(idx5d(ii,jj,kk,4,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nhfields))
	  pxx=hfields_in(idx5d(ii,jj,kk,5,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nhfields))
	  pyy=hfields_in(idx5d(ii,jj,kk,6,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nhfields))
	  pzz=hfields_in(idx5d(ii,jj,kk,7,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nhfields))
	  pxy=hfields_in(idx5d(ii,jj,kk,8,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nhfields))
	  pxz=hfields_in(idx5d(ii,jj,kk,9,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nhfields))
	  pyz=hfields_in(idx5d(ii,jj,kk,10,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nhfields))
	  
	  opress=hfields_out(idx5d(ii,jj,kk,1,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nhfields))
	  ou=hfields_out(idx5d(ii,jj,kk,2,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nhfields))
	  ov=hfields_out(idx5d(ii,jj,kk,3,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nhfields))
	  ow=hfields_out(idx5d(ii,jj,kk,4,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nhfields))
	  opxx=hfields_out(idx5d(ii,jj,kk,5,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nhfields))
	  opyy=hfields_out(idx5d(ii,jj,kk,6,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nhfields))
	  opzz=hfields_out(idx5d(ii,jj,kk,7,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nhfields))
	  opxy=hfields_out(idx5d(ii,jj,kk,8,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nhfields))
	  opxz=hfields_out(idx5d(ii,jj,kk,9,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nhfields))
	  opyz=hfields_out(idx5d(ii,jj,kk,10,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nhfields))
	  
#ifdef EXPLICITEQ 
	  uu=HALF*(u*u+v*v+w*w)*invcssq
	  do l=1,nlinks
		 udotc=(u*dex(l) + v*dey(l)+ w*dez(l))*invcssq
		 feq=p(l)*(press + (udotc+0.5_db*udotc*udotc - uu))
                 !fneq1=(HALF/(cssq*cssq))*( (dex(l)*dex(l)-cssq)*pxx &
		 !     + (dey(l)*dey(l)-cssq)*pyy + (dez(l)*dez(l)-cssq)*pzz &
	         !     + TWO*(dex(l)*dey(l))*pxy + TWO*(dex(l)*dez(l))*pxz &
		 !     + TWO*(dey(l)*dez(l))*pyz)
                 fpost=feq!+fneq1
                 pxx=pxx - fpost*(dex(l)*dex(l))
                 pyy=pyy - fpost*(dey(l)*dey(l))
                 pzz=pzz - fpost*(dez(l)*dez(l))
                 pxy=pxy - fpost*(dex(l)*dey(l))
                 pxz=pxz - fpost*(dex(l)*dez(l))
                 pyz=pyz - fpost*(dey(l)*dez(l))
	  enddo
#else
	  pxx=pxx - cssq*press - u*u 
	  pyy=pyy - cssq*press - v*v 
	  pzz=pzz - cssq*press - w*w 
	  pxy=pxy - u*v
	  pxz=pxz - u*w
	  pyz=pyz - v*w
#endif

#ifdef TWOCOMPONENT
	  !visc_loc it is used to store the local viscosity
	  visc_loc=(rho_r*visc1*phi_loc+(ONE-phi_loc)*visc2*rho_b)/rhophi_loc
#else
#ifdef SMAGORINSKI
	  visc_loc=visc1
#endif
#endif

#ifdef SMAGORINSKI
	  QQ=pxx*pxx + pyy*pyy + pzz*pzz + &
	   TWO*(pxy*pxy + pxz*pxz + pyz*pyz)  !QQ i sused to store the double contraction of flux tensor (Frobenius norm) 
	  !!!smago
	  omega_loc= 0.5_db + (1.0_db/6.0_db)*(3.0_db*visc_loc + &   !visc_loc it is used to store the local viscosity
	   sqrt((3.0*visc_loc)**2.0 + 0.053*18.0*sqrt(2.0*QQ)/rhophi_loc)) !it is tau
	  omega_loc=1.0_db/omega_loc !it is omega

#else
#ifdef TWOCOMPONENT
	  omega_loc=(visc_loc/cssq + 0.5_db) !it is tau   !visc_loc it is used to store the local viscosity
	  omega_loc=1.0_db/omega_loc !it is omega
#else
	  omega_loc=omega
#endif
#endif

      

      do l=1,nlinks
		  lopp=opp(l)
		  iii=i+ex(lopp)
		  jjj=j+ey(lopp)
		  kkk=k+ez(lopp)
		  if(isfluid(iii,jjj,kkk).ne.0) cycle 
		  feq=p(l)*(press)
		  fneq1=(HALF/(cssq*cssq))*( (dex(l)*dex(l)-cssq)*pxx &
		   + (dey(l)*dey(l)-cssq)*pyy + (dez(l)*dez(l)-cssq)*pzz &
	       + TWO*(dex(l)*dey(l))*pxy + TWO*(dex(l)*dez(l))*pxz &
		   + TWO*(dey(l)*dez(l))*pyz)
          ! F_discr = p(l)*(dex(l)*forcex &
           ! + dey(l)*forcey &
           ! + dez(l)*forcez)/cssq
		  fpost=feq + (ONE-omega_loc)*p(l)*fneq1 !+ HALF*(F_discr)	
		  opress=opress + fpost
		  ou=ou + fpost*dex(l)
		  ov=ov + fpost*dey(l)
		  ow=ow + fpost*dez(l)
		  opxx=opxx + fpost*dex(l)*dex(l)
          opyy=opyy + fpost*dey(l)*dey(l)
          opzz=opzz + fpost*dez(l)*dez(l)
          opxy=opxy + fpost*dex(l)*dey(l)
          opxz=opxz + fpost*dex(l)*dez(l)
          opyz=opyz + fpost*dey(l)*dez(l)		  
      enddo

		 
	  hfields_out(idx5d(ii,jj,kk,1,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nhfields))=opress
	  hfields_out(idx5d(ii,jj,kk,2,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nhfields))=ou
	  hfields_out(idx5d(ii,jj,kk,3,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nhfields))=ov
	  hfields_out(idx5d(ii,jj,kk,4,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nhfields))=ow
	  hfields_out(idx5d(ii,jj,kk,5,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nhfields))=opxx
	  hfields_out(idx5d(ii,jj,kk,6,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nhfields))=opyy
	  hfields_out(idx5d(ii,jj,kk,7,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nhfields))=opzz
	  hfields_out(idx5d(ii,jj,kk,8,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nhfields))=opxy
	  hfields_out(idx5d(ii,jj,kk,9,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nhfields))=opxz
	  hfields_out(idx5d(ii,jj,kk,10,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nhfields))=opyz
	              
     return

   endsubroutine LB_int_boundary_kernel   
   
   subroutine PHI_int_boundary_cuda(hfields_s,phifields_s)

      implicit none
      real(kind=db), allocatable, dimension(:) :: hfields_s,phifields_s
 
      
#ifdef TWOCOMPONENT
      !$acc wait
      istat = cudaDeviceSynchronize

!$acc host_data use_device(flop,nx,ny,nz,coords,isfluid &
       !$acc& ,visc1,visc2,rho_r,rho_b,invrho_r,invrho_b,sharp_c,beta,kapp,tau_diff,sigma &
#ifdef WETTABILITY
       !$acc& ,wettab_r,wettab_b &
#endif         
#ifdef MONOD
	   !$acc& ,mu_max,Ks &
#endif
       !$acc& ,ntothfields,ntotphifields,ntotauxfields,ntotlocauxfields,ntotforces &
	   !$acc& ,hfields_s,phifields_s,auxfields,locauxfields,forces)
       call PHI_int_boundary_kernel<<<dimGrid, dimBlock>>>(flop,nx,ny,nz,coords,isfluid &
       ,visc2,rho_r,rho_b,invrho_r,invrho_b,sharp_c,beta,kapp,tau_diff,sigma &
#ifdef WETTABILITY
       ,wettab_r,wettab_b &
#endif         
#ifdef MONOD
	   ,mu_max,Ks &
#endif      
       ,ntothfields,ntotphifields,ntotauxfields,ntotlocauxfields,ntotforces &
	   ,hfields_s,phifields_s,auxfields,locauxfields,forces)
!$acc end host_data
      istat = cudaDeviceSynchronize
      istat = cudaGetLastError()
      if (istat/= cudaSuccess) then
        if(myrank==0)write(6,*) 'status after at ', __LINE__ ,' of file ', __FILE__ ,' :'
        if(myrank==0)write(6,*) cudaGetErrorString(istat)
        call doerror(6,'ERROR in PHI_int_boundary_cuda')
      endif
      !$acc wait        
#endif
      
      return
      
   endsubroutine PHI_int_boundary_cuda 
   
   attributes(global) subroutine PHI_int_boundary_kernel(flop,nx,ny,nz,coords,isfluid &
    ,visc2,rho_r,rho_b,invrho_r,invrho_b,sharp_c,beta,kapp,tau_diff,sigma &
#ifdef WETTABILITY  
    ,wettab_r,wettab_b
#endif       
#ifdef MONOD	
    ,mu_max,Ks &
#endif
    ,ntothfields,ntotphifields,ntotauxfields,ntotlocauxfields,ntotforces &
       ,hfields_s,phifields_s,auxfields_s,locauxfields_s,forces_s)

      implicit none
      integer :: flop,nx,ny,nz
      integer, dimension(3) :: coords
      integer(kind=isf), dimension(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff) :: isfluid
     
#ifdef WETTABILITY  
      real(kind=db) :: wettab_r,wettab_b  
#endif   
      real(kind=db) :: visc2,rho_r,rho_b,invrho_r,invrho_b,sharp_c,beta,kapp,tau_diff,sigma  
#ifdef MONOD
      real(kind=db) :: mu_max,Ks
#endif
      integer :: ntothfields,ntotphifields,ntotauxfields,ntotlocauxfields,ntotforces
      real(kind=db), dimension(ntothfields) :: hfields_s
      real(kind=db), dimension(ntotphifields) :: phifields_s
      real(kind=db), dimension(ntotauxfields) :: auxfields_s
      real(kind=db), dimension(ntotlocauxfields) :: locauxfields_s
      real(kind=db), dimension(ntotforces) :: forces_s
      
      integer :: i,j,k,l
      integer :: gi,gj,gk
      integer :: myblock,ii,jj,kk
      integer :: iii,jjj,kkk
      real(kind=db) :: gradfix,gradfiy,gradfiz,grad_parallel,modgrad,phi_fluid,phi_ghost
      real(kind=db) :: loc_u,loc_v,loc_w,loc_phi,theta_rad,cot_theta,dphi_dz
      integer :: oii,ojj,okk
      integer :: oxblock,oyblock,ozblock,omyblock
      real(kind=db) :: wettab_r_sub=90.0_db
    
      i = (blockIdx%x-1) * TILE_DIMx_d + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy_d + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz_d + threadIdx%z
      
      if(isfluid(i,j,k) .ne. 0)return
      
      gi=nx*coords(1)+i
      gj=ny*coords(2)+j
      gk=nz*coords(3)+k
      
      myblock=blockIdx%x+blockIdx%y*nxblock_d+blockIdx%z*nxyblock_d+1
      
      ii=threadIdx%x
      jj=threadIdx%y
      kk=threadIdx%z
	   
#ifdef WETTABILITY       
      wettab_r_sub=wettab_r
#endif                    
      
	  do l = 1, nlinks
		iii = i + ex(l)
		jjj = j + ey(l)
		kkk = k + ez(l)

		if (isfluid(ii,jj,kk) .ne. -1) cycle  ! only fluid neighbor
		
		oxblock=(iii+2*TILE_DIMx_d-1)/TILE_DIMx_d   
		oyblock=(jjj+2*TILE_DIMy_d-1)/TILE_DIMy_d     
		ozblock=(kkk+2*TILE_DIMz_d-1)/TILE_DIMz_d 
		omyblock=(oxblock-1)+(oyblock-1)*nxblock_d+(ozblock-1)*nxyblock_d+1
		oii=iii-oxblock*TILE_DIMx_d+2*TILE_DIMx_d
		ojj=jjj-oyblock*TILE_DIMy_d+2*TILE_DIMy_d
		okk=kkk-ozblock*TILE_DIMz_d+2*TILE_DIMz_d

		! Found fluid neighbor: enforce contact angle via ghost node extrapolation
		phi_fluid = phifields_s(idx5d(oii,ojj,okk,1,omyblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nphifields))
		
		
		! Estimate gradient parallel to wall
		modgrad=auxfields_s(idx5d(oii,ojj,okk,4,omyblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nauxfields)) !modgrad
		gradfix=auxfields_s(idx5d(oii,ojj,okk,1,omyblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nauxfields))*modgrad !normx*modgrad
		gradfiy=auxfields_s(idx5d(oii,ojj,okk,2,omyblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nauxfields))*modgrad !normy*modgrad
		gradfiz=auxfields_s(idx5d(oii,ojj,okk,3,omyblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nauxfields))*modgrad !normz*modgrad
        
        grad_parallel=ZERO
		if(l.eq.1 .or. l.eq.2)then
			grad_parallel = sqrt(gradfiy**2 + gradfiz**2)
		elseif(l.eq.3 .or. l.eq.4)then
			grad_parallel = sqrt(gradfix**2 + gradfiz**2)
		elseif(l.eq.5 .or. l.eq.6)then
			grad_parallel = sqrt(gradfix**2 + gradfiy**2)
		endif
		
		! Contact angle correction
		theta_rad = (180.0_db-wettab_r_sub) * pi_greek / 180.0_db
		cot_theta = 1.0_db / tan(theta_rad)

		dphi_dz = - grad_parallel * cot_theta 

		  

		phi_ghost = phi_fluid + dphi_dz  ! extrapolate from fluid node

		! Clamp to [0,1]
		loc_phi = max(0.0_db, min(1.0_db, phi_ghost))
		
		exit
	  end do
      
      phifields_s(idx5d(ii,jj,kk,1,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nphifields))=loc_phi
       
      return
      
   endsubroutine PHI_int_boundary_kernel
   
   subroutine phi_sum_count_cuda(hfields_s,phifields_s)

      implicit none
      real(kind=db), allocatable, dimension(:) :: hfields_s,phifields_s
 
      
#ifdef TWOCOMPONENT
      !$acc wait
      istat = cudaDeviceSynchronize

!$acc host_data use_device(flop,nx,ny,nz,coords,isfluid &
       !$acc& ,visc1,visc2,rho_r,rho_b,invrho_r,invrho_b,sharp_c,beta,kapp,tau_diff,sigma &
       !$acc& ,ntothfields,ntotphifields,ntotauxfields,ntotlocauxfields &
	   !$acc& ,hfields_s,phifields_s,auxfields,locauxfields,global_phi_sum,global_count)
       call phi_sum_count_kernel<<<dimGrid, dimBlock>>>(flop,nx,ny,nz,coords,isfluid &
       ,visc2,rho_r,rho_b,invrho_r,invrho_b,sharp_c,beta,kapp,tau_diff,sigma &
       ,ntothfields,ntotphifields,ntotauxfields,ntotlocauxfields &
	   ,hfields_s,phifields_s,auxfields,locauxfields,global_phi_sum,global_count)
!$acc end host_data
      istat = cudaDeviceSynchronize
      istat = cudaGetLastError()
      if (istat/= cudaSuccess) then
        if(myrank==0)write(6,*) 'status after at ', __LINE__ ,' of file ', __FILE__ ,' :'
        if(myrank==0)write(6,*) cudaGetErrorString(istat)
        call doerror(6,'ERROR in phi_sum_count_cuda')
      endif
      !$acc wait        
#endif
      
      return
      
   endsubroutine phi_sum_count_cuda 
   
   attributes(global) subroutine phi_sum_count_kernel(flop,nx,ny,nz,coords,isfluid &
    ,visc2,rho_r,rho_b,invrho_r,invrho_b,sharp_c,beta,kapp,tau_diff,sigma &
    ,ntothfields,ntotphifields,ntotauxfields,ntotlocauxfields, &
       hfields_s,phifields_s,auxfields_s,locauxfields_s,loc_phi_sum,cnt)

      implicit none
      integer :: flop,nx,ny,nz
      integer, dimension(3) :: coords
      integer(kind=isf), dimension(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff) :: isfluid
      
      real(kind=db) :: visc2,rho_r,rho_b,invrho_r,invrho_b,sharp_c,beta,kapp,tau_diff,sigma  

      integer :: ntothfields,ntotphifields,ntotauxfields,ntotlocauxfields
      real(kind=db), dimension(ntothfields) :: hfields_s
      real(kind=db), dimension(ntotphifields) :: phifields_s
      real(kind=db), dimension(ntotauxfields) :: auxfields_s
      real(kind=db), dimension(ntotlocauxfields) :: locauxfields_s
      real(kind=db) :: loc_phi_sum
      integer :: cnt
      
      integer :: i,j,k,l
      integer :: gi,gj,gk
      integer :: myblock,ii,jj,kk
      real(kind=db) :: gradfix,gradfiy,gradfiz,grad_parallel,modgrad,phi_fluid
      real(kind=db) :: loc_u,loc_v,loc_w,loc_phi,theta_rad,cot_theta,dphi_dz
      real(kind=db) :: dummy
      integer :: dummy_i
    
      i = (blockIdx%x-1) * TILE_DIMx_d + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy_d + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz_d + threadIdx%z
      
      if(abs(isfluid(i,j,k)) /= 1)return
      
      gi=nx*coords(1)+i
      gj=ny*coords(2)+j
      gk=nz*coords(3)+k
      
      myblock=blockIdx%x+blockIdx%y*nxblock_d+blockIdx%z*nxyblock_d+1
      
      ii=threadIdx%x
      jj=threadIdx%y
      kk=threadIdx%z

      loc_phi = phifields_s(idx5d(ii,jj,kk,1,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nphifields))

      dummy=atomicAdd(loc_phi_sum, loc_phi)

      if(loc_phi > 0.5d0 .and. loc_phi < 0.9d0)then
        dummy_i=atomicAdd(cnt, 1)
      endif
    
      return
    
  end subroutine phi_sum_count_kernel
  
  subroutine apply_lagrangian_phi_cuda(hfields_s,phifields_s)

      implicit none
      real(kind=db), allocatable, dimension(:) :: hfields_s,phifields_s
 
      
#ifdef TWOCOMPONENT
      !$acc wait
      istat = cudaDeviceSynchronize

!$acc host_data use_device(flop,nx,ny,nz,coords,isfluid &
       !$acc& ,visc1,visc2,rho_r,rho_b,invrho_r,invrho_b,sharp_c,beta,kapp,tau_diff,sigma &
       !$acc& ,ntothfields,ntotphifields,ntotauxfields,ntotlocauxfields &
	   !$acc& ,hfields_s,phifields_s,auxfields,locauxfields,corr)
       call apply_lagrangian_phi_kernel<<<dimGrid, dimBlock>>>(flop,nx,ny,nz,coords,isfluid &
       ,visc2,rho_r,rho_b,invrho_r,invrho_b,sharp_c,beta,kapp,tau_diff,sigma &
       ,ntothfields,ntotphifields,ntotauxfields,ntotlocauxfields &
	   ,hfields_s,phifields_s,auxfields,locauxfields,corr)
!$acc end host_data
      istat = cudaDeviceSynchronize
      istat = cudaGetLastError()
      if (istat/= cudaSuccess) then
        if(myrank==0)write(6,*) 'status after at ', __LINE__ ,' of file ', __FILE__ ,' :'
        if(myrank==0)write(6,*) cudaGetErrorString(istat)
        call doerror(6,'ERROR in phi_sum_count_cuda')
      endif
      !$acc wait        
#endif
      
      return
      
   endsubroutine apply_lagrangian_phi_cuda 
   
      attributes(global) subroutine apply_lagrangian_phi_kernel(flop,nx,ny,nz,coords,isfluid &
    ,visc2,rho_r,rho_b,invrho_r,invrho_b,sharp_c,beta,kapp,tau_diff,sigma &
    ,ntothfields,ntotphifields,ntotauxfields,ntotlocauxfields, &
       hfields_s,phifields_s,auxfields_s,locauxfields_s,loc_corr)

      implicit none
      integer :: flop,nx,ny,nz
      integer, dimension(3) :: coords
      integer(kind=isf), dimension(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff) :: isfluid
      
      real(kind=db) :: visc2,rho_r,rho_b,invrho_r,invrho_b,sharp_c,beta,kapp,tau_diff,sigma  

      integer :: ntothfields,ntotphifields,ntotauxfields,ntotlocauxfields
      real(kind=db), dimension(ntothfields) :: hfields_s
      real(kind=db), dimension(ntotphifields) :: phifields_s
      real(kind=db), dimension(ntotauxfields) :: auxfields_s
      real(kind=db), dimension(ntotlocauxfields) :: locauxfields_s
      real(kind=db) :: loc_corr
      
      integer :: i,j,k,l
      integer :: gi,gj,gk
      integer :: myblock,ii,jj,kk
      real(kind=db) :: gradfix,gradfiy,gradfiz,grad_parallel,modgrad,phi_fluid
      real(kind=db) :: loc_u,loc_v,loc_w,loc_phi,theta_rad,cot_theta,dphi_dz
    
      i = (blockIdx%x-1) * TILE_DIMx_d + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy_d + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz_d + threadIdx%z
      
      if(abs(isfluid(i,j,k)) /= 1)return
      
      gi=nx*coords(1)+i
      gj=ny*coords(2)+j
      gk=nz*coords(3)+k
      
      myblock=blockIdx%x+blockIdx%y*nxblock_d+blockIdx%z*nxyblock_d+1
      
      ii=threadIdx%x
      jj=threadIdx%y
      kk=threadIdx%z

      loc_phi = phifields_s(idx5d(ii,jj,kk,1,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nphifields))

      if(loc_phi > 0.5d0 .and. loc_phi < 0.9d0)then
        phifields_s(idx5d(ii,jj,kk,1,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nphifields))=loc_phi + loc_corr
      endif
    
      return
    
  end subroutine apply_lagrangian_phi_kernel

endmodule
