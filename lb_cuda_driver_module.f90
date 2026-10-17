#include "defines.h"
#if !defined(_OPENACC)  
#error "To use this module the macros _OPENACC should be defined."
#endif


!------------------------------------------------------------------------------------------------
! Host wrappers that launch the CUDA Fortran kernels of the time step, and setup_cuda, which sets
! the launch configuration.
!
! Tiles and grids: the local subdomain is split into tiles of TILE_DIMx x TILE_DIMy x TILE_DIMz
! nodes (8x8x8 by default), and the field arrays hold one extra layer of halo tiles on each side
! (see allocate_module.f90). dimBlock has one thread per node of a tile and dimGrid one block per
! tile; the CUDA block (bx,by,bz) of dimGrid works on the array block (bx+1,by+1,bz+1). For the
! overlap of computation and communication (ASYNCMPI) the tiles are split into disjoint sets:
!   dimGridInt  interior tiles, which do not touch a face of the subdomain (_int wrappers);
!   dimGridz    tiles of the two z faces, including edges and corners;
!   dimGridy    tiles of the two y faces, without those of the z faces;
!   dimGridx    tiles of the two x faces, without those of the y and z faces.
! The _ext wrappers launch the minus and plus kernels of the z, y and x faces; together with the
! _int kernel every tile is processed once. dimBlockshared, (TILE_DIMx+2) x (TILE_DIMy+2) x
! (TILE_DIMz+2) = 10x10x10 threads for a tile and its one-node halo, and dimGridhalo/dimBlockhalo
! are used only by the test kernels of this module.
!
! The wrappers wait for the OpenACC queues and synchronize the device before the launch, and pass
! the device addresses of the arrays through the OpenACC host_data use_device construct. Kernel
! launches are asynchronous: after the launch, cudaGetLastError reports launch errors (for
! instance an invalid configuration or too much shared memory), and cudaDeviceSynchronize waits
! for the kernel to finish (most wrappers also check its status, which reports errors raised
! during the execution). The results are therefore complete before the host continues with
! OpenACC loops, MPI calls on device buffers or output, and an error is reported by the wrapper
! that launched the failing kernel.
!------------------------------------------------------------------------------------------------
module lb_cuda_driver

   use vars
   use iso_c_binding
   use cudafor
#ifdef _NVML
   use nvml_interface
#endif
   use mpi_template, only: coords,dostop,doerror,mydev,myrank,nprocs,nbuff,nbuffbvec
   use lb_cuda_vars
#ifdef TWOCOMPONENT
   use lb_cuda_auxfields, only: compute_norm_interface_kernel,compute_div_theta_n_kernel, &
#ifdef CSF
    compute_csf_force_kernel, &
#endif
    compute_norm_interface_kernel_int,compute_norm_interface_kernel_ext, &
    compute_norm_interface_kernel_xminus,compute_norm_interface_kernel_xplus, &
    compute_norm_interface_kernel_yminus,compute_norm_interface_kernel_yplus, &
    compute_norm_interface_kernel_zminus,compute_norm_interface_kernel_zplus
   use lb_cuda_repulsive, only: thinfilm_scan_mark_kernel,repulsive_flux_normal_kernel
#endif
   use lb_cuda_moments, only: moments_LB_kernel,moments_LB_kernel_int, &
    moments_LB_kernel_xminus,moments_LB_kernel_xplus, &
    moments_LB_kernel_yminus,moments_LB_kernel_yplus, &
    moments_LB_kernel_zminus,moments_LB_kernel_zplus
#if defined(HIGHORDER) && (LATTICE == 27)
   use lb_cuda_fused, only: fused_LB_kernel_h,fused_LB_kernel_bb_h
#else
   use lb_cuda_fused, only: fused_LB_kernel,fused_LB_kernel_bb
#endif
#ifdef TWOCOMPONENT
   use lb_cuda_update_phi, only: update_phifields_kernel,update_phifields_kernel_int, &
    update_phifields_kernel_xminus,update_phifields_kernel_xplus, &
    update_phifields_kernel_yminus,update_phifields_kernel_yplus, &
    update_phifields_kernel_zminus,update_phifields_kernel_zplus
#endif
   use lb_cuda_boundary, only: LB_int_boundary_kernel,PHI_int_boundary_kernel, &
    phi_sum_count_kernel,apply_lagrangian_phi_kernel

   implicit none
   

contains

   ! setup_cuda: selects the GPU of the rank, checks the tile sizes (nonzero, at most 1024 threads
   ! per tile, nx, ny, nz multiples of the tile sizes), copies the sizes to constant memory, sets
   ! the launch configurations, and computes the number of array blocks along each direction,
   ! nxblock = nx/TILE_DIMx + 2 (one halo tile on each side), and nblocks.
   subroutine setup_cuda
      implicit none
      
      integer :: nxb,nyb,nzb
      
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
      
      lx_d=lx
      ly_d=ly
      lz_d=lz
      nx_d=nx
      ny_d=ny
      nz_d=nz
      
      if (mod(nx, TILE_DIMx)/= 0) then
        if(myrank==0)write(6,'(a,i4)') 'nx must be a multiple of TILE_DIMx=',TILE_DIMx
        call dostop
      end if
      if (mod(ny, TILE_DIMy) /= 0) then
        if(myrank==0)write(6,'(a,i4)') 'ny must be a multiple of TILE_DIMy=',TILE_DIMy
        call dostop
      end if
      if (mod(nz, TILE_DIMz) /= 0) then
        if(myrank==0)write(6,'(a,i4)') 'nz must be a multiple of TILE_DIMz=',TILE_DIMz
        call dostop
      end if
      
      nxb = nx / TILE_DIMx
      nyb = ny / TILE_DIMy
      nzb = nz / TILE_DIMz
      
      dimGrid  = dim3(nxb,nyb,nzb)
      dimBlock = dim3(TILE_DIMx, TILE_DIMy, TILE_DIMz)
      
      
      dimGridInt = dim3(nxb -2,nyb -2,nzb -2)
      ldodimGridInt=.true.
      if(nxb-2<=0)ldodimGridInt=.false.
      if(nyb-2<=0)ldodimGridInt=.false.
      if(nzb-2<=0)ldodimGridInt=.false.
      

      
      dimGridhalo  = dim3(nxb +2,nyb +2,nzb +2)
      dimBlockhalo = dim3(TILE_DIMx, TILE_DIMy, TILE_DIMz)
      
      dimBlockshared = dim3(TILE_DIMx +2, TILE_DIMy +2, TILE_DIMz +2)
      
      dimGridz  = dim3(nxb,nyb,1)    !xy faces also doing edges xz yz and corners
      
      dimGridy  = dim3(nxb,1,nzb-2)  !xz faces also doing edge xy
      ldodimGridy=.true.
      if(nzb-2<=0)ldodimGridy=.false.
      
      dimGridx  = dim3(1,nyb-2,nzb-2) !only yz faces
      ldodimGridx=.true.
      if(nyb-2<=0)ldodimGridx=.false.
      if(nzb-2<=0)ldodimGridx=.false.
      
      ! With a single tile along a direction, the minus face kernel covers it and the plus face
      ! kernel is not launched.
      ldodimGridx_twice=.true.
      ldodimGridy_twice=.true.
      ldodimGridz_twice=.true.
      
      if(nxb<2)ldodimGridx_twice=.false.
      if(nyb<2)ldodimGridy_twice=.false.
      if(nzb<2)ldodimGridz_twice=.false.
      
      
      
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
#ifdef ASYNCMPI
        write(6,*)'ldodimGridInt',ldodimGridInt
        write(6,*)'ldodimGridx',ldodimGridx
        write(6,*)'ldodimGridy',ldodimGridy
#endif
      endif
#endif
      
   
   endsubroutine setup_cuda
   
   ! Debug helper: launches one of the test kernels below, which print the indices of selected
   ! nodes to check the mapping between threads, nodes and tiles. It is not called during a run.
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

      i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z
      
      gi=nx*coords(1)+i
      gj=ny*coords(2)+j
      gk=nz*coords(3)+k
      
      myblock=blockIdx%x+blockIdx%y*nxblock_d+blockIdx%z*nxyblock_d+1
      
      !if(myblock<344)then
      if(gi==1 .and. gj==1 .and. gk==31)then
        write(*,*)'i',gi,gj,gk,myblock,myrank
      endif
     !if(gi==1) write(*,*)'debug',i,j,k,phi(i,j,k)
     
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
      
      !if(myblock<344)then
      if(gi==1 .and. gj==1 .and. gk==36)then
        write(*,*)'i',gi,gj,gk,myrank,lk
      endif
     !if(gi==1) write(*,*)'debug',i,j,k,phi(i,j,k)
     
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

      i = (blockIdx%x) * TILE_DIMx + li
      j = (blockIdx%y) * TILE_DIMy + lj
      k = (blockIdx%z) * TILE_DIMz + lk
      
      gi=nx*coords(1)+i
      gj=ny*coords(2)+j
      gk=nz*coords(3)+k
      
      xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
	  yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
	  zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
      
      myblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
      
      
      !if(myblock<344)then
      if(gi==17 .and. gj==17 .and. gk==17)then
        write(*,*)'i',li,lj,lk,myrank
      endif
     !if(gi==1) write(*,*)'debug',i,j,k,phi(i,j,k)
     
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

      i = (blockIdx%x-1) * TILE_DIMx + li
      j = (blockIdx%y-1) * TILE_DIMy + lj
      
      gi=nx*coords(1)+i
      gj=ny*coords(2)+j
      
      xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
	  yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
	  
	  
	  k = lk
	  gk=nz*coords(3)+k
	  zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
      
      myblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
      
      !if(myblock<344)then
      if(gi==1 .and. gj==1 .and. gk==38)then
        write(*,*)'asd',gi,gj,gk,myrank,lk
      endif
     !if(gi==1) write(*,*)'debug',i,j,k,phi(i,j,k)
     
     !phi_old(i,j,k)=phi(i,j,k)
     
      k = ((nzblock_d-2)-1) * TILE_DIMz + lk
	  gk=nz*coords(3)+k
	  zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
      
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

      i = (blockIdx%x-1) * TILE_DIMx + li + TILE_DIMx
      k = (blockIdx%z-1) * TILE_DIMz + lk
      
      gi=nx*coords(1)+i
      gk=nz*coords(3)+k
      
      xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
	  zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
	  
	  
      j = lj
      gj=ny*coords(2)+j
      yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
      
      myblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
      
      !if(myblock<344)then
      if(gi==8 .and. gj==6 .and. gk==1)then
        write(*,*)'i',gi,gj,gk,myrank,li
      endif
     !if(gi==1) write(*,*)'debug',i,j,k,phi(i,j,k)
     
     !phi_old(i,j,k)=phi(i,j,k)
 
      j = ((nyblock_d-2)-1) * TILE_DIMy + lj
      gj=ny*coords(2)+j
      yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
      
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

      
      j = (blockIdx%y-1) * TILE_DIMy + lj + TILE_DIMy
      k = (blockIdx%z-1) * TILE_DIMz + lk + TILE_DIMz
      !if(myrank==0)write(*,*)'debug',j
      
     
      gj=ny*coords(2)+j
      gk=nz*coords(3)+k
      
      
	  yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
	  zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
	  
	  i = li
	  gi=nx*coords(1)+i
	  xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
      
      myblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
      
      !if(myblock<344)then
      if(gi==1 .and. gj==8 .and. gk==8)then
        write(*,*)'i',gi,gj,gk,lj,lk
      endif
     !if(gi==1) write(*,*)'debug',i,j,k,phi(i,j,k)
     
     !phi_old(i,j,k)=phi(i,j,k)
 
	  i = ((nxblock_d-2)-1) * TILE_DIMx + li
	  gi=nx*coords(1)+i
	  xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
      
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

      i = (blockIdx%x-2) * TILE_DIMx + threadIdx%x
      j = (blockIdx%y-2) * TILE_DIMy + threadIdx%y
      k = (blockIdx%z-2) * TILE_DIMz + threadIdx%z
      
      gi=nx*coords(1)+i
      gj=ny*coords(2)+j
      gk=nz*coords(3)+k
      
      xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
	  yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
	  zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
      
      myblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
      
      !if(myblock==20 .and. myrank==1)then
      !if(blockIdx%x==1 .and. blockIdx%y==1 .and. blockIdx%z==1)then
      if(gi==1 .and. gj==1 .and. gk==31)then
        write(*,*)'e',gi,gj,gk,myblock,myrank
      endif
     !if(gi==1) write(*,*)'debug',i,j,k,phi(i,j,k)
     
     !phi_old(i,j,k)=phi(i,j,k)
     
      
 end subroutine test_LB_kernel_halo
 
 ! Interface geometry (TWOCOMPONENT; otherwise the wrappers do nothing). The kernels read phi from
 ! phifields_s, including the halo, and store the unit normal n, |grad phi| and phi (1-phi) n in
 ! auxfields (fields 1-3, 4 and 5-7) and the Laplacian of phi in locauxfields (field 1).
 ! compute_norm_interface_cuda processes all tiles; the _int and _ext variants process the
 ! interior and the face tiles, so that the halo exchange of auxfields can overlap the interior
 ! computation (ASYNCMPI).
 subroutine compute_norm_interface_cuda(phifields_s)

      implicit none
      real(kind=strdb), allocatable, dimension(:,:,:,:,:) :: phifields_s
      
#ifdef TWOCOMPONENT
      !$acc wait
      istat = cudaDeviceSynchronize

!$acc host_data use_device(step &
       !$acc& ,iprobe,jprobe,kprobe,flop,nx,ny,nz,coords,isfluid,rho_r,rho_b &
       !$acc& ,ntotphifields,ntotauxfields,ntotlocauxfields,phifields_s,auxfields,locauxfields)
       call compute_norm_interface_kernel<<<dimGrid, dimBlock>>>(step, &
        iprobe,jprobe,kprobe,flop,nx,ny,nz,coords,isfluid,rho_r,rho_b, &
        ntotphifields,ntotauxfields,ntotlocauxfields,phifields_s,auxfields,locauxfields)
       !$acc end host_data
         
      
       istat = cudaGetLastError()
       if (istat/= cudaSuccess) then
         if(myrank==0) then
           write(6,*) 'status after at ', __LINE__ ,' of file ', __FILE__ ,' :'
           write(6,*) cudaGetErrorString(istat)
         endif
         call doerror(6,'ERROR in compute_norm_interface_cuda (launch)')
       endif
      
       istat = cudaDeviceSynchronize()
       if (istat /= cudaSuccess) then
        if(myrank==0) then
          write(6,*) 'sync error at ', __LINE__, ' file ', __FILE__
          write(6,*) cudaGetErrorString(istat)
        endif
        call doerror(6,'ERROR in compute_norm_interface_cuda (sync)')
       endif
       !$acc wait 
#endif
      
      return
      
 endsubroutine compute_norm_interface_cuda
 
 subroutine compute_norm_interface_cuda_int(phifields_s)

      implicit none
      real(kind=strdb), allocatable, dimension(:,:,:,:,:) :: phifields_s
      
#ifdef TWOCOMPONENT
      !$acc wait
      istat = cudaDeviceSynchronize
      if(ldodimGridInt)then
!$acc host_data use_device(step &
       !$acc& ,iprobe,jprobe,kprobe,flop,nx,ny,nz,coords,isfluid,rho_r,rho_b &
       !$acc& ,ntotphifields,ntotauxfields,ntotlocauxfields,phifields_s,auxfields,locauxfields)
	   call compute_norm_interface_kernel_int<<<dimGridInt, dimBlock>>>(step, &
        iprobe,jprobe,kprobe,flop,nx,ny,nz,coords,isfluid,rho_r,rho_b, &
        ntotphifields,ntotauxfields,ntotlocauxfields,phifields_s,auxfields,locauxfields)
       !$acc end host_data
        
       istat = cudaGetLastError()
       if (istat/= cudaSuccess) then
         if(myrank==0) then
           write(6,*) 'status after at ', __LINE__ ,' of file ', __FILE__ ,' :'
           write(6,*) cudaGetErrorString(istat)
         endif
         call doerror(6,'ERROR in compute_norm_interface_cuda_int (launch)')
       endif
      
       istat = cudaDeviceSynchronize()
       if (istat /= cudaSuccess) then
        if(myrank==0) then
          write(6,*) 'sync error at ', __LINE__, ' file ', __FILE__
          write(6,*) cudaGetErrorString(istat)
        endif
        call doerror(6,'ERROR in compute_norm_interface_cuda_int (sync)')
       endif
       !$acc wait  
       endif
#endif
      
      return
      
 endsubroutine compute_norm_interface_cuda_int
 
 subroutine compute_norm_interface_cuda_ext(phifields_s)

      implicit none
      real(kind=strdb), allocatable, dimension(:,:,:,:,:) :: phifields_s
      
#ifdef TWOCOMPONENT
      !$acc wait
      istat = cudaDeviceSynchronize

       !$acc host_data use_device(step &
       !$acc& ,iprobe,jprobe,kprobe,flop,nx,ny,nz,coords,isfluid,rho_r,rho_b &
       !$acc& ,ntotphifields,ntotauxfields,ntotlocauxfields,phifields_s,auxfields,locauxfields)
	   call compute_norm_interface_kernel_zminus<<<dimGridz, dimBlock>>>(step, &
        iprobe,jprobe,kprobe,flop,nx,ny,nz,coords,isfluid,rho_r,rho_b, &
        ntotphifields,ntotauxfields,ntotlocauxfields,phifields_s,auxfields,locauxfields)
       !$acc end host_data
          
        
       istat = cudaGetLastError()
       if(istat/= cudaSuccess) then
         if(myrank==0) then
           write(6,*) 'status after at ', __LINE__ ,' of file ', __FILE__ ,' :'
           write(6,*) cudaGetErrorString(istat)
         endif
         call doerror(6,'ERROR in compute_norm_interface_kernel_zminus (launch)')
       endif
      
       istat = cudaDeviceSynchronize()
       if(istat /= cudaSuccess) then
        if(myrank==0) then
          write(6,*) 'sync error at ', __LINE__, ' file ', __FILE__
          write(6,*) cudaGetErrorString(istat)
        endif
        call doerror(6,'ERROR in compute_norm_interface_kernel_zminus (sync)')
       endif
       !$acc wait
       if(ldodimGridz_twice)then
       !$acc host_data use_device(step &
       !$acc& ,iprobe,jprobe,kprobe,flop,nx,ny,nz,coords,isfluid,rho_r,rho_b &
       !$acc& ,ntotphifields,ntotauxfields,ntotlocauxfields,phifields_s,auxfields,locauxfields)
       call compute_norm_interface_kernel_zplus<<<dimGridz, dimBlock>>>(step, &
        iprobe,jprobe,kprobe,flop,nx,ny,nz,coords,isfluid,rho_r,rho_b, &
        ntotphifields,ntotauxfields,ntotlocauxfields,phifields_s,auxfields,locauxfields)
       !$acc end host_data
       
          
       istat = cudaGetLastError()
       if(istat/= cudaSuccess) then
         if(myrank==0) then
           write(6,*) 'status after at ', __LINE__ ,' of file ', __FILE__ ,' :'
           write(6,*) cudaGetErrorString(istat)
         endif
         call doerror(6,'ERROR in compute_norm_interface_kernel_zplus (launch)')
       endif
      
       istat = cudaDeviceSynchronize()
       if(istat /= cudaSuccess) then
        if(myrank==0) then
          write(6,*) 'sync error at ', __LINE__, ' file ', __FILE__
          write(6,*) cudaGetErrorString(istat)
        endif
        call doerror(6,'ERROR in compute_norm_interface_kernel_zplus (sync)')
       endif
       !$acc wait 
       endif
       
	   if(ldodimGridy)then
         !$acc host_data use_device(step &
         !$acc& ,iprobe,jprobe,kprobe,flop,nx,ny,nz,coords,isfluid,rho_r,rho_b &
         !$acc& ,ntotphifields,ntotauxfields,ntotlocauxfields,phifields_s,auxfields,locauxfields)	    
	     call compute_norm_interface_kernel_yminus<<<dimGridy, dimBlock>>>(step, &
          iprobe,jprobe,kprobe,flop,nx,ny,nz,coords,isfluid,rho_r,rho_b, &
          ntotphifields,ntotauxfields,ntotlocauxfields,phifields_s,auxfields,locauxfields)
         !$acc end host_data
             

         istat = cudaGetLastError()
         if(istat/= cudaSuccess) then
           if(myrank==0) then
             write(6,*) 'status after at ', __LINE__ ,' of file ', __FILE__ ,' :'
             write(6,*) cudaGetErrorString(istat)
           endif
           call doerror(6,'ERROR in compute_norm_interface_kernel_yminus (launch)')
         endif
      
         istat = cudaDeviceSynchronize()
         if(istat /= cudaSuccess) then
          if(myrank==0) then
            write(6,*) 'sync error at ', __LINE__, ' file ', __FILE__
            write(6,*) cudaGetErrorString(istat)
          endif
          call doerror(6,'ERROR in compute_norm_interface_kernel_yminus (sync)')
         endif
         !$acc wait
       endif
              
       if(ldodimGridy .and. ldodimGridy_twice)then
         !$acc host_data use_device(step &
         !$acc& ,iprobe,jprobe,kprobe,flop,nx,ny,nz,coords,isfluid,rho_r,rho_b &
         !$acc& ,ntotphifields,ntotauxfields,ntotlocauxfields,phifields_s,auxfields,locauxfields)
         call compute_norm_interface_kernel_yplus<<<dimGridy, dimBlock>>>(step, &
          iprobe,jprobe,kprobe,flop,nx,ny,nz,coords,isfluid,rho_r,rho_b, &
          ntotphifields,ntotauxfields,ntotlocauxfields,phifields_s,auxfields,locauxfields)
         !$acc end host_data
            

         istat = cudaGetLastError()
         if(istat/= cudaSuccess) then
           if(myrank==0) then
             write(6,*) 'status after at ', __LINE__ ,' of file ', __FILE__ ,' :'
             write(6,*) cudaGetErrorString(istat)
           endif
           call doerror(6,'ERROR in compute_norm_interface_kernel_yplus (launch)')
         endif
         
         istat = cudaDeviceSynchronize()
         if(istat /= cudaSuccess) then
          if(myrank==0) then
            write(6,*) 'sync error at ', __LINE__, ' file ', __FILE__
            write(6,*) cudaGetErrorString(istat)
          endif
          call doerror(6,'ERROR in compute_norm_interface_kernel_yplus (sync)')
         endif
         !$acc wait 
       endif
              
	   if(ldodimGridx)then
         !$acc host_data use_device(step &
         !$acc& ,iprobe,jprobe,kprobe,flop,nx,ny,nz,coords,isfluid,rho_r,rho_b &
         !$acc& ,ntotphifields,ntotauxfields,ntotlocauxfields,phifields_s,auxfields,locauxfields)
	     call compute_norm_interface_kernel_xminus<<<dimGridx, dimBlock>>>(step, &
          iprobe,jprobe,kprobe,flop,nx,ny,nz,coords,isfluid,rho_r,rho_b, &
          ntotphifields,ntotauxfields,ntotlocauxfields,phifields_s,auxfields,locauxfields)
         !$acc end host_data
            

         istat = cudaGetLastError()
         if(istat/= cudaSuccess) then
           if(myrank==0) then
             write(6,*) 'status after at ', __LINE__ ,' of file ', __FILE__ ,' :'
             write(6,*) cudaGetErrorString(istat)
           endif
           call doerror(6,'ERROR in compute_norm_interface_kernel_xminus (launch)')
         endif
      
         istat = cudaDeviceSynchronize()
         if(istat /= cudaSuccess) then
          if(myrank==0) then
            write(6,*) 'sync error at ', __LINE__, ' file ', __FILE__
            write(6,*) cudaGetErrorString(istat)
          endif
          call doerror(6,'ERROR in compute_norm_interface_kernel_xminus (sync)')
         endif
         !$acc wait 
       endif
              
       if(ldodimGridx .and. ldodimGridx_twice)then
         !$acc host_data use_device(step &
         !$acc& ,iprobe,jprobe,kprobe,flop,nx,ny,nz,coords,isfluid,rho_r,rho_b &
         !$acc& ,ntotphifields,ntotauxfields,ntotlocauxfields,phifields_s,auxfields,locauxfields)
         call compute_norm_interface_kernel_xplus<<<dimGridx, dimBlock>>>(step, &
          iprobe,jprobe,kprobe,flop,nx,ny,nz,coords,isfluid,rho_r,rho_b, &
          ntotphifields,ntotauxfields,ntotlocauxfields,phifields_s,auxfields,locauxfields)
         !$acc end host_data
             

         istat = cudaGetLastError()
         if(istat/= cudaSuccess) then
           if(myrank==0) then
             write(6,*) 'status after at ', __LINE__ ,' of file ', __FILE__ ,' :'
             write(6,*) cudaGetErrorString(istat)
           endif
           call doerror(6,'ERROR in compute_norm_interface_kernel_xplus (launch)')
         endif
      
         istat = cudaDeviceSynchronize()
         if(istat /= cudaSuccess) then
          if(myrank==0) then
            write(6,*) 'sync error at ', __LINE__, ' file ', __FILE__
            write(6,*) cudaGetErrorString(istat)
          endif
          call doerror(6,'ERROR in compute_norm_interface_kernel_xplus (sync)')
         endif
         !$acc wait
       endif       
    
#endif
      
      return
      
 endsubroutine compute_norm_interface_cuda_ext
   
 ! div(phi (1-phi) n), computed from auxfields 5-7 and stored in locauxfields field 2
 ! (TWOCOMPONENT only).
 subroutine compute_div_theta_n(phifields_s)

      implicit none
      real(kind=strdb), allocatable, dimension(:,:,:,:,:) :: phifields_s

      
#ifdef TWOCOMPONENT
      !$acc wait
      istat = cudaDeviceSynchronize

!$acc host_data use_device(flop,nx,ny,nz,coords,isfluid &
       !$acc& ,rho_r,rho_b,ntotphifields,ntotauxfields,ntotlocauxfields,phifields_s,auxfields,locauxfields)
       call compute_div_theta_n_kernel<<<dimGrid, dimBlock>>>(flop,nx,ny,nz,coords,isfluid, &
        rho_r,rho_b,ntotphifields,ntotauxfields,ntotlocauxfields,phifields_s,auxfields,locauxfields)
!$acc end host_data
      istat = cudaDeviceSynchronize
      istat = cudaGetLastError()
      if (istat/= cudaSuccess) then
        if(myrank==0)write(6,*) 'status after at ', __LINE__ ,' of file ', __FILE__ ,' :'
        if(myrank==0)write(6,*) cudaGetErrorString(istat)
        call doerror(6,'ERROR in compute_div_theta_n')
      endif
      !$acc wait        
#endif
      
      return
      
   endsubroutine compute_div_theta_n 

   ! CSF surface tension: curvature from the normals in auxfields, and force
   ! sigma * curvature * n * |grad phi| stored in forces; moments_LB_kernel adds it to the force.
#ifdef CSF
   subroutine compute_csf_force

      implicit none
      !$acc wait
      istat=cudaDeviceSynchronize()
      !$acc host_data use_device(nx,ny,nz,isfluid,sigma,ntotauxfields,ntotforces,auxfields,forces)
      call compute_csf_force_kernel<<<dimGrid,dimBlock>>>(nx,ny,nz,isfluid,sigma, &
       ntotauxfields,ntotforces,auxfields,forces)
      !$acc end host_data
      istat=cudaDeviceSynchronize()
      istat=cudaGetLastError()
      if(istat /= cudaSuccess)then
         if(myrank==0)write(6,*)cudaGetErrorString(istat)
         call doerror(6,'ERROR in compute_csf_force')
      endif
      !$acc wait
   endsubroutine compute_csf_force
#endif
      
   ! REPULSIVE_FLUX (TWOCOMPONENT only): thinfilm_scan_mark_kernel and repulsive_flux_normal_kernel,
   ! see lb_cuda_repulsive_module.f90. Without REPULSIVE_FLUX these wrappers do nothing.
   subroutine thinfilm_scan_mark_cuda(phifields_s)

      implicit none
      real(kind=strdb), allocatable, dimension(:,:,:,:,:) :: phifields_s
 
   
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
   
   subroutine repulsive_flux_normal_cuda(phifields_s)

      implicit none
      real(kind=strdb), allocatable, dimension(:,:,:,:,:) :: phifields_s
 
   
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
 
 ! Moment update before the collision (not used with FUSED_MOMENTS). moments_LB_kernel computes
 ! the force on every fluid node (body force and the model-dependent terms of TWOCOMPONENT,
 ! MULTIHIT and ELASTIC_FORCE), adds half the force per unit density to the velocity and
 ! (u F + F u)/(2 rho) to the second-order moments of hfields_s in place, and stores in forces the
 ! force used by the fused kernel. hfields_old is the other buffer: only its velocity, the
 ! velocity of the previous time step, is read (DENSRATIO viscous term, ELASTIC_FORCE).
 ! moments_LB_cuda processes all tiles; moments_LB_cuda_int and moments_LB_cuda_ext launch the
 ! interior and face variants used with ASYNCMPI (see the module header).
 subroutine moments_LB_cuda(hfields_old,hfields_s &
#ifdef TWOCOMPONENT    
  ,phifields_s &
#endif 
 )

      implicit none
      
      real(kind=strdb), allocatable, dimension(:,:,:,:,:) :: hfields_old,hfields_s
#ifdef TWOCOMPONENT          
      real(kind=strdb), allocatable, dimension(:,:,:,:,:) :: phifields_s
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
	   ,hfields_old,hfields_s &
#ifdef TWOCOMPONENT 	   
	   ,auxfields,locauxfields &
#endif
	   ,forces)
!$acc end host_data
      
      istat = cudaGetLastError()
      if (istat/= cudaSuccess) then
        if(myrank==0) then
          write(6,*) 'status after at ', __LINE__ ,' of file ', __FILE__ ,' :'
          write(6,*) cudaGetErrorString(istat)
        endif
        call doerror(6,'ERROR in moments_LB_cuda_int (launch)')
      endif
      
      istat = cudaDeviceSynchronize()
      if (istat /= cudaSuccess) then
       if(myrank==0) then
         write(6,*) 'sync error at ', __LINE__, ' file ', __FILE__
         write(6,*) cudaGetErrorString(istat)
       endif
       call doerror(6,'ERROR in moments_LB_cuda_int (sync)')
      endif
      
      !$acc wait        
      
   endsubroutine moments_LB_cuda
   
   subroutine moments_LB_cuda_int(hfields_old,hfields_s &
#ifdef TWOCOMPONENT    
  ,phifields_s &
#endif 
 )

      implicit none
      
      real(kind=strdb), allocatable, dimension(:,:,:,:,:) :: hfields_old,hfields_s
#ifdef TWOCOMPONENT          
      real(kind=strdb), allocatable, dimension(:,:,:,:,:) :: phifields_s
#endif 

!      if(myrank==0)write(6,*)'step ',step, __LINE__ , __FILE__
      !$acc wait
      istat = cudaDeviceSynchronize
      if(ldodimGridInt)then
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
       
       call moments_LB_kernel_int<<<dimGridInt, dimBlock>>>(step &
        ,iprobe,jprobe,kprobe,flop,nx,ny,nz,coords,isfluid &   
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
	    ,hfields_old,hfields_s &
#ifdef TWOCOMPONENT 	   
	   ,auxfields,locauxfields &
#endif
	   ,forces)
	   !$acc end host_data
      endif
      
            
      istat = cudaGetLastError()
      if (istat/= cudaSuccess) then
        if(myrank==0) then
          write(6,*) 'status after at ', __LINE__ ,' of file ', __FILE__ ,' :'
          write(6,*) cudaGetErrorString(istat)
        endif
        call doerror(6,'ERROR in moments_LB_cuda_int (launch)')
      endif
      
      istat = cudaDeviceSynchronize()
      if (istat /= cudaSuccess) then
       if(myrank==0) then
         write(6,*) 'sync error at ', __LINE__, ' file ', __FILE__
         write(6,*) cudaGetErrorString(istat)
       endif
       call doerror(6,'ERROR in moments_LB_cuda_int (sync)')
      endif
      !$acc wait        
      
   endsubroutine moments_LB_cuda_int
   
   subroutine moments_LB_cuda_ext(hfields_old,hfields_s &
#ifdef TWOCOMPONENT    
  ,phifields_s &
#endif 
 )

      implicit none
      
      real(kind=strdb), allocatable, dimension(:,:,:,:,:) :: hfields_old,hfields_s
#ifdef TWOCOMPONENT          
      real(kind=strdb), allocatable, dimension(:,:,:,:,:) :: phifields_s
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
       
       call moments_LB_kernel_zminus<<<dimGridz, dimBlock>>>(step,iprobe,jprobe,kprobe,flop,nx,ny,nz,coords,isfluid &   
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
	   ,hfields_old,hfields_s &
#ifdef TWOCOMPONENT 	   
	   ,auxfields,locauxfields &
#endif
	   ,forces)
      !$acc end host_data
      
      istat = cudaGetLastError()
      if (istat/= cudaSuccess) then
        if(myrank==0) then
          write(6,*) 'status after at ', __LINE__ ,' of file ', __FILE__ ,' :'
          write(6,*) cudaGetErrorString(istat)
        endif
        call doerror(6,'ERROR in moments_LB_cuda_zminus (launch)')
      endif
      
      istat = cudaDeviceSynchronize()
      if (istat /= cudaSuccess) then
       if(myrank==0) then
         write(6,*) 'sync error at ', __LINE__, ' file ', __FILE__
         write(6,*) cudaGetErrorString(istat)
       endif
       call doerror(6,'ERROR in moments_LB_cuda_zminus (sync)')
      endif
      !$acc wait  
      
      if(ldodimGridz_twice)then
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
      call moments_LB_kernel_zplus<<<dimGridz, dimBlock>>>(step &
       ,iprobe,jprobe,kprobe,flop,nx,ny,nz,coords,isfluid &   
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
	   ,hfields_old,hfields_s &
#ifdef TWOCOMPONENT 	   
	   ,auxfields,locauxfields &
#endif
	   ,forces)
      !$acc end host_data
      
      istat = cudaGetLastError()
      if (istat/= cudaSuccess) then
        if(myrank==0) then
          write(6,*) 'status after at ', __LINE__ ,' of file ', __FILE__ ,' :'
          write(6,*) cudaGetErrorString(istat)
        endif
        call doerror(6,'ERROR in moments_LB_cuda_zplus (launch)')
      endif
      
      istat = cudaDeviceSynchronize()
      if (istat /= cudaSuccess) then
       if(myrank==0) then
         write(6,*) 'sync error at ', __LINE__, ' file ', __FILE__
         write(6,*) cudaGetErrorString(istat)
       endif
       call doerror(6,'ERROR in moments_LB_cuda_zplus (sync)')
      endif
      !$acc wait  
      endif
      
	  if(ldodimGridy)then
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
	    call moments_LB_kernel_yminus<<<dimGridy, dimBlock>>>(step &
	   ,iprobe,jprobe,kprobe,flop,nx,ny,nz,coords,isfluid &   
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
	   ,hfields_old,hfields_s &
#ifdef TWOCOMPONENT 	   
	   ,auxfields,locauxfields &
#endif
	   ,forces)
      !$acc end host_data
      
      istat = cudaGetLastError()
      if (istat/= cudaSuccess) then
        if(myrank==0) then
          write(6,*) 'status after at ', __LINE__ ,' of file ', __FILE__ ,' :'
          write(6,*) cudaGetErrorString(istat)
        endif
        call doerror(6,'ERROR in moments_LB_cuda_yminus (launch)')
      endif
      
      istat = cudaDeviceSynchronize()
      if (istat /= cudaSuccess) then
       if(myrank==0) then
         write(6,*) 'sync error at ', __LINE__, ' file ', __FILE__
         write(6,*) cudaGetErrorString(istat)
       endif
       call doerror(6,'ERROR in moments_LB_cuda_yminus (sync)')
      endif
      !$acc wait  
      endif
      
      if(ldodimGridy .and. ldodimGridy_twice)then
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
        call moments_LB_kernel_yplus<<<dimGridy, dimBlock>>>(step &
       ,iprobe,jprobe,kprobe,flop,nx,ny,nz,coords,isfluid &   
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
	   ,hfields_old,hfields_s &
#ifdef TWOCOMPONENT 	   
	   ,auxfields,locauxfields &
#endif
	   ,forces)
      !$acc end host_data
            
      istat = cudaGetLastError()
      if (istat/= cudaSuccess) then
        if(myrank==0) then
          write(6,*) 'status after at ', __LINE__ ,' of file ', __FILE__ ,' :'
          write(6,*) cudaGetErrorString(istat)
        endif
        call doerror(6,'ERROR in moments_LB_cuda_yplus (launch)')
      endif
      
      istat = cudaDeviceSynchronize()
      if (istat /= cudaSuccess) then
       if(myrank==0) then
         write(6,*) 'sync error at ', __LINE__, ' file ', __FILE__
         write(6,*) cudaGetErrorString(istat)
       endif
       call doerror(6,'ERROR in moments_LB_cuda_yplus (sync)')
      endif
      !$acc wait  
      endif
      
      if(ldodimGridx)then
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
        call moments_LB_kernel_xminus<<<dimGridx, dimBlock>>>(step &
	   ,iprobe,jprobe,kprobe,flop,nx,ny,nz,coords,isfluid &   
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
	   ,hfields_old,hfields_s &
#ifdef TWOCOMPONENT 	   
	   ,auxfields,locauxfields &
#endif
	   ,forces)
      !$acc end host_data
      
      istat = cudaGetLastError()
      if (istat/= cudaSuccess) then
        if(myrank==0) then
          write(6,*) 'status after at ', __LINE__ ,' of file ', __FILE__ ,' :'
          write(6,*) cudaGetErrorString(istat)
        endif
        call doerror(6,'ERROR in moments_LB_cuda_xminus (launch)')
      endif
      
      istat = cudaDeviceSynchronize()
      if (istat /= cudaSuccess) then
       if(myrank==0) then
         write(6,*) 'sync error at ', __LINE__, ' file ', __FILE__
         write(6,*) cudaGetErrorString(istat)
       endif
       call doerror(6,'ERROR in moments_LB_cuda_xminus (sync)')
      endif
      !$acc wait  
      endif
      
      if(ldodimGridx .and. ldodimGridx_twice)then
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
        call moments_LB_kernel_xplus<<<dimGridx, dimBlock>>>(step &
       ,iprobe,jprobe,kprobe,flop,nx,ny,nz,coords,isfluid &   
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
	   ,hfields_old,hfields_s &
#ifdef TWOCOMPONENT 	   
	   ,auxfields,locauxfields &
#endif
	   ,forces)
      !$acc end host_data
      
      istat = cudaGetLastError()
      if (istat/= cudaSuccess) then
        if(myrank==0) then
          write(6,*) 'status after at ', __LINE__ ,' of file ', __FILE__ ,' :'
          write(6,*) cudaGetErrorString(istat)
        endif
        call doerror(6,'ERROR in moments_LB_cuda_xplus (launch)')
      endif
      
      istat = cudaDeviceSynchronize()
      if (istat /= cudaSuccess) then
       if(myrank==0) then
         write(6,*) 'sync error at ', __LINE__, ' file ', __FILE__
         write(6,*) cudaGetErrorString(istat)
       endif
       call doerror(6,'ERROR in moments_LB_cuda_xplus (sync)')
      endif
      !$acc wait  
      endif

        
      
   endsubroutine moments_LB_cuda_ext

   ! Fused collide-and-stream step (docs/FUSED_KERNELS.md): the kernel reads the moments of
   ! hfields_in on each tile and its one-node halo, and writes the post-streaming moments of the
   ! tile nodes into hfields_out. The kernel is chosen at compile time: fused_LB_kernel_h or
   ! fused_LB_kernel_bb_h for high-order D3Q27, fused_LB_kernel or fused_LB_kernel_bb otherwise,
   ! the _bb variants with BOUNCE_BACK. It runs one block of TILE_DIMx*TILE_DIMy*TILE_DIMz threads
   ! per tile (dimGrid, dimBlock) with nbytesf bytes of dynamic shared memory. More than 48 KB per
   ! block must be requested explicitly, so at the first call the limit is raised to nbytesf with
   ! cudaFuncSetAttribute and the largest shared-memory carveout is preferred; the run stops if
   ! the GPU cannot provide the memory. With FUSED_MOMENTS the kernel also applies the moment
   ! update of moments_LB_kernel.
   subroutine fused_LB_cuda(hfields_in,hfields_out &
#ifdef TWOCOMPONENT    
   ,phifields_s &
#endif   
   )

      implicit none
      
      real(kind=strdb), allocatable, dimension(:,:,:,:,:) :: hfields_in,hfields_out
#ifdef TWOCOMPONENT       
      real(kind=strdb), allocatable, dimension(:,:,:,:,:) :: phifields_s
#endif
      ! nbytesf is in bytes: the kind db (4 or 8) is used as the size of a real of the arithmetic
      ! precision.
#if defined(HIGHORDER) && (LATTICE == 27)
      ! ten second-order coefficients, u and (1-omega) Pi_neq of the tile and its halo in dynamic
      ! shared memory, later reused for the seventeen high-order coefficients
      integer, parameter :: nbytesf=19*(TILE_DIMx+2)*(TILE_DIMy+2)*(TILE_DIMz+2)*db
#else
      ! the ten post-collision Hermite coefficients of the tile and its halo in dynamic shared memory
      integer, parameter :: nbytesf=10*(TILE_DIMx+2)*(TILE_DIMy+2)*(TILE_DIMz+2)*db
#endif
      logical, save :: lfused_set=.false.
      
!      if(myrank==0)write(6,*)'step ',step, __LINE__ , __FILE__
      !$acc wait
      istat = cudaDeviceSynchronize
      if(.not. lfused_set)then
        ! allow more than 48 KB of dynamic shared memory and prefer the largest shared carveout
#if defined(HIGHORDER) && (LATTICE == 27) && defined(BOUNCE_BACK)
        istat = cudaFuncSetAttribute(fused_LB_kernel_bb_h,cudaFuncAttributeMaxDynamicSharedMemorySize,nbytesf)
#elif defined(HIGHORDER) && (LATTICE == 27)
        istat = cudaFuncSetAttribute(fused_LB_kernel_h,cudaFuncAttributeMaxDynamicSharedMemorySize,nbytesf)
#elif defined(BOUNCE_BACK)
        istat = cudaFuncSetAttribute(fused_LB_kernel_bb,cudaFuncAttributeMaxDynamicSharedMemorySize,nbytesf)
#else
        istat = cudaFuncSetAttribute(fused_LB_kernel,cudaFuncAttributeMaxDynamicSharedMemorySize,nbytesf)
#endif
        if(istat /= cudaSuccess)call doerror(6,'the fused collide-and-stream kernel needs more shared '// &
         'memory per block than this GPU provides')
#if defined(HIGHORDER) && (LATTICE == 27) && defined(BOUNCE_BACK)
        istat = cudaFuncSetAttribute(fused_LB_kernel_bb_h,cudaFuncAttributePreferredSharedMemoryCarveout,100)
#elif defined(HIGHORDER) && (LATTICE == 27)
        istat = cudaFuncSetAttribute(fused_LB_kernel_h,cudaFuncAttributePreferredSharedMemoryCarveout,100)
#elif defined(BOUNCE_BACK)
        istat = cudaFuncSetAttribute(fused_LB_kernel_bb,cudaFuncAttributePreferredSharedMemoryCarveout,100)
#else
        istat = cudaFuncSetAttribute(fused_LB_kernel,cudaFuncAttributePreferredSharedMemoryCarveout,100)
#endif
        lfused_set=.true.
      endif

!$acc host_data use_device(step,iprobe,jprobe,kprobe,flip,flop,nx,ny,nz,coords,isfluid & 
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
#if defined(HIGHORDER) && (LATTICE == 27) && defined(BOUNCE_BACK)
      call fused_LB_kernel_bb_h<<<dimGrid,dimBlock,nbytesf>>>(step &
#elif defined(HIGHORDER) && (LATTICE == 27)
      call fused_LB_kernel_h<<<dimGrid,dimBlock,nbytesf>>>(step &
#elif defined(BOUNCE_BACK)
      call fused_LB_kernel_bb<<<dimGrid,dimBlock,nbytesf>>>(step &
#else
      call fused_LB_kernel<<<dimGrid,dimBlock,nbytesf>>>(step &
#endif
       ,iprobe,jprobe,kprobe,flip,flop,nx,ny,nz,coords,isfluid &  
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
	   ,hfields_in,hfields_out, &
#ifdef TWOCOMPONENT 
	   auxfields,locauxfields, &
#endif   
	   forces)

      !$acc end host_data
      ! The GPU power is sampled while the kernel runs (before the synchronization below), once every
      ! stamp_term steps; integrator integrates it over time.
#if defined(_NVML) && defined(GETPOWER)
      if(mod(step-1,stamp_term).eq.0)then 
        p_mw = get_gpu_power_index(mydev_c)
      endif
#endif

      
      istat = cudaGetLastError()         ! or cudaPeekAtLastError
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
   
   ! Phase-field update (TWOCOMPONENT only): explicit step of the phase-field equation,
   !   phi_out = phi_in - u.grad(phi) + tau_diff lap(phi) - sharp_c div(phi (1-phi) n),
   ! with the velocity of hfields_s and the fields in auxfields and locauxfields (MONOD adds a
   ! source term). update_phifields processes all tiles; update_phifields_int and
   ! update_phifields_ext process the interior and the face tiles (ASYNCMPI).
   subroutine update_phifields(hfields_s,phifields_in,phifields_out)

      implicit none
      real(kind=strdb), allocatable, dimension(:,:,:,:,:) :: hfields_s,phifields_in,phifields_out
 
      
#ifdef TWOCOMPONENT
      !$acc wait
      istat = cudaDeviceSynchronize

!$acc host_data use_device(step,iprobe,jprobe,kprobe,flop,nx,ny,nz,coords,isfluid &
#ifdef TWOCOMPONENT    
       !$acc& ,visc1,visc2,rho_r,rho_b,invrho_r,invrho_b,sharp_c,beta,kapp,tau_diff,sigma &
#ifdef MONOD
	   !$acc& ,mu_max,Ks &
#endif
#endif
       !$acc& ,ntothfields,ntotphifields,ntotauxfields,ntotlocauxfields,ntotforces &
	   !$acc& ,hfields_s,phifields_in,phifields_out,auxfields,locauxfields,forces)
       call update_phifields_kernel<<<dimGrid, dimBlock>>>(step &
        ,iprobe,jprobe,kprobe,flop,nx,ny,nz,coords,isfluid &
#ifdef TWOCOMPONENT           
        ,visc2,rho_r,rho_b,invrho_r,invrho_b,sharp_c,beta,kapp,tau_diff,sigma &
#ifdef MONOD
	    ,mu_max,Ks &
#endif   
#endif   
        ,ntothfields,ntotphifields,ntotauxfields,ntotlocauxfields,ntotforces &
	    ,hfields_s,phifields_in,phifields_out,auxfields,locauxfields,forces)
      !$acc end host_data
      
      istat = cudaGetLastError()
      if (istat/= cudaSuccess) then
        if(myrank==0) then
          write(6,*) 'status after at ', __LINE__ ,' of file ', __FILE__ ,' :'
          write(6,*) cudaGetErrorString(istat)
        endif
        call doerror(6,'ERROR in update_phifields (launch)')
      endif
      
      istat = cudaDeviceSynchronize()
      if (istat /= cudaSuccess) then
       if(myrank==0) then
         write(6,*) 'sync error at ', __LINE__, ' file ', __FILE__
         write(6,*) cudaGetErrorString(istat)
       endif
       call doerror(6,'ERROR in update_phifields (sync)')
      endif

      !$acc wait        
#endif
      
      return
      
   endsubroutine update_phifields
   
   subroutine update_phifields_int(hfields_s,phifields_in,phifields_out)

      implicit none
      real(kind=strdb), allocatable, dimension(:,:,:,:,:) :: hfields_s,phifields_in,phifields_out
 
      
#ifdef TWOCOMPONENT
      !$acc wait
      istat = cudaDeviceSynchronize
      if(ldodimGridInt)then
!$acc host_data use_device(step,iprobe,jprobe,kprobe,flop,nx,ny,nz,coords,isfluid &
#ifdef TWOCOMPONENT    
       !$acc& ,visc1,visc2,rho_r,rho_b,invrho_r,invrho_b,sharp_c,beta,kapp,tau_diff,sigma &
#ifdef MONOD
	   !$acc& ,mu_max,Ks &
#endif
#endif
       !$acc& ,ntothfields,ntotphifields,ntotauxfields,ntotlocauxfields,ntotforces &
	   !$acc& ,hfields_s,phifields_in,phifields_out,auxfields,locauxfields,forces)
       call update_phifields_kernel_int<<<dimGridInt, dimBlock>>>(step &
        ,iprobe,jprobe,kprobe,flop &
        ,nx,ny,nz,coords,isfluid &
#ifdef TWOCOMPONENT           
        ,visc2,rho_r,rho_b,invrho_r,invrho_b,sharp_c,beta,kapp,tau_diff,sigma &
#ifdef MONOD
	    ,mu_max,Ks &
#endif   
#endif   
        ,ntothfields,ntotphifields,ntotauxfields,ntotlocauxfields,ntotforces &
	    ,hfields_s,phifields_in,phifields_out,auxfields,locauxfields,forces)
      !$acc end host_data
      
      istat = cudaGetLastError()
      if (istat/= cudaSuccess) then
        if(myrank==0) then
          write(6,*) 'status after at ', __LINE__ ,' of file ', __FILE__ ,' :'
          write(6,*) cudaGetErrorString(istat)
        endif
        call doerror(6,'ERROR in update_phifields_int (launch)')
      endif
      
      istat = cudaDeviceSynchronize()
      if (istat /= cudaSuccess) then
       if(myrank==0) then
         write(6,*) 'sync error at ', __LINE__, ' file ', __FILE__
         write(6,*) cudaGetErrorString(istat)
       endif
       call doerror(6,'ERROR in update_phifields_int (sync)')
      endif
      !$acc wait   
      endif     
#endif
      
      return
      
   endsubroutine update_phifields_int
   
   subroutine update_phifields_ext(hfields_s,phifields_in,phifields_out)

      implicit none
      real(kind=strdb), allocatable, dimension(:,:,:,:,:) :: hfields_s,phifields_in,phifields_out
 
      
#ifdef TWOCOMPONENT
      !$acc wait
      istat = cudaDeviceSynchronize

      !$acc host_data use_device(step,iprobe,jprobe,kprobe,flop,nx,ny,nz,coords,isfluid &
#ifdef TWOCOMPONENT    
       !$acc& ,visc1,visc2,rho_r,rho_b,invrho_r,invrho_b,sharp_c,beta,kapp,tau_diff,sigma &
#ifdef MONOD
	   !$acc& ,mu_max,Ks &
#endif
#endif
       !$acc& ,ntothfields,ntotphifields,ntotauxfields,ntotlocauxfields,ntotforces &
	   !$acc& ,hfields_s,phifields_in,phifields_out,auxfields,locauxfields,forces)
	   
       call update_phifields_kernel_zminus<<<dimGridz, dimBlock>>>(step &
        ,iprobe,jprobe,kprobe,flop &
        ,nx,ny,nz,coords,isfluid &
#ifdef TWOCOMPONENT           
        ,visc2,rho_r,rho_b,invrho_r,invrho_b,sharp_c,beta,kapp,tau_diff,sigma &
#ifdef MONOD
	    ,mu_max,Ks &
#endif   
#endif   
        ,ntothfields,ntotphifields,ntotauxfields,ntotlocauxfields,ntotforces &
	    ,hfields_s,phifields_in,phifields_out,auxfields,locauxfields,forces)
      !$acc end host_data
      istat = cudaGetLastError()
      if (istat/= cudaSuccess) then
        if(myrank==0) then
          write(6,*) 'status after at ', __LINE__ ,' of file ', __FILE__ ,' :'
          write(6,*) cudaGetErrorString(istat)
        endif
        call doerror(6,'ERROR in update_phifields_kernel_zminus (launch)')
      endif
      
      istat = cudaDeviceSynchronize()
      if (istat /= cudaSuccess) then
       if(myrank==0) then
         write(6,*) 'sync error at ', __LINE__, ' file ', __FILE__
         write(6,*) cudaGetErrorString(istat)
       endif
       call doerror(6,'ERROR in update_phifields_kernel_zminus (sync)')
      endif
      !$acc wait  
      if(ldodimGridz_twice)then
      !$acc host_data use_device(step,iprobe,jprobe,kprobe,flop,nx,ny,nz,coords,isfluid &
#ifdef TWOCOMPONENT    
       !$acc& ,visc1,visc2,rho_r,rho_b,invrho_r,invrho_b,sharp_c,beta,kapp,tau_diff,sigma &
#ifdef MONOD
	   !$acc& ,mu_max,Ks &
#endif
#endif
       !$acc& ,ntothfields,ntotphifields,ntotauxfields,ntotlocauxfields,ntotforces &
	   !$acc& ,hfields_s,phifields_in,phifields_out,auxfields,locauxfields,forces)
      call update_phifields_kernel_zplus<<<dimGridz, dimBlock>>>(step &
       ,iprobe,jprobe,kprobe,flop &
       ,nx,ny,nz,coords,isfluid &
#ifdef TWOCOMPONENT           
       ,visc2,rho_r,rho_b,invrho_r,invrho_b,sharp_c,beta,kapp,tau_diff,sigma &
#ifdef MONOD
	   ,mu_max,Ks &
#endif   
#endif   
       ,ntothfields,ntotphifields,ntotauxfields,ntotlocauxfields,ntotforces &
	   ,hfields_s,phifields_in,phifields_out,auxfields,locauxfields,forces)
      !$acc end host_data
      istat = cudaGetLastError()
      if (istat/= cudaSuccess) then
        if(myrank==0) then
          write(6,*) 'status after at ', __LINE__ ,' of file ', __FILE__ ,' :'
          write(6,*) cudaGetErrorString(istat)
        endif
        call doerror(6,'ERROR in update_phifields_kernel_zplus (launch)')
      endif
      
      istat = cudaDeviceSynchronize()
      if (istat /= cudaSuccess) then
       if(myrank==0) then
         write(6,*) 'sync error at ', __LINE__, ' file ', __FILE__
         write(6,*) cudaGetErrorString(istat)
       endif
       call doerror(6,'ERROR in update_phifields_kernel_zplus (sync)')
      endif
      !$acc wait  
      endif
      
       if(ldodimGridy)then
      !$acc host_data use_device(step,iprobe,jprobe,kprobe,flop,nx,ny,nz,coords,isfluid &
#ifdef TWOCOMPONENT    
       !$acc& ,visc1,visc2,rho_r,rho_b,invrho_r,invrho_b,sharp_c,beta,kapp,tau_diff,sigma &
#ifdef MONOD
	   !$acc& ,mu_max,Ks &
#endif
#endif
       !$acc& ,ntothfields,ntotphifields,ntotauxfields,ntotlocauxfields,ntotforces &
	   !$acc& ,hfields_s,phifields_in,phifields_out,auxfields,locauxfields,forces)
         call update_phifields_kernel_yminus<<<dimGridy, dimBlock>>>(step &
        ,iprobe,jprobe,kprobe,flop &
        ,nx,ny,nz,coords,isfluid &
#ifdef TWOCOMPONENT           
        ,visc2,rho_r,rho_b,invrho_r,invrho_b,sharp_c,beta,kapp,tau_diff,sigma &
#ifdef MONOD
	    ,mu_max,Ks &
#endif   
#endif   
        ,ntothfields,ntotphifields,ntotauxfields,ntotlocauxfields,ntotforces &
	    ,hfields_s,phifields_in,phifields_out,auxfields,locauxfields,forces)
      !$acc end host_data
      istat = cudaGetLastError()
      if (istat/= cudaSuccess) then
        if(myrank==0) then
          write(6,*) 'status after at ', __LINE__ ,' of file ', __FILE__ ,' :'
          write(6,*) cudaGetErrorString(istat)
        endif
        call doerror(6,'ERROR in update_phifields_kernel_yminus (launch)')
      endif
      
      istat = cudaDeviceSynchronize()
      if (istat /= cudaSuccess) then
       if(myrank==0) then
         write(6,*) 'sync error at ', __LINE__, ' file ', __FILE__
         write(6,*) cudaGetErrorString(istat)
       endif
       call doerror(6,'ERROR in update_phifields_kernel_yminus (sync)')
      endif
      !$acc wait  
      endif
      
      if(ldodimGridy .and. ldodimGridy_twice)then
      !$acc host_data use_device(step,iprobe,jprobe,kprobe,flop,nx,ny,nz,coords,isfluid &
#ifdef TWOCOMPONENT    
       !$acc& ,visc1,visc2,rho_r,rho_b,invrho_r,invrho_b,sharp_c,beta,kapp,tau_diff,sigma &
#ifdef MONOD
	   !$acc& ,mu_max,Ks &
#endif
#endif
       !$acc& ,ntothfields,ntotphifields,ntotauxfields,ntotlocauxfields,ntotforces &
	   !$acc& ,hfields_s,phifields_in,phifields_out,auxfields,locauxfields,forces)
        call update_phifields_kernel_yplus<<<dimGridy, dimBlock>>>(step &
       ,iprobe,jprobe,kprobe,flop &
       ,nx,ny,nz,coords,isfluid &
#ifdef TWOCOMPONENT           
       ,visc2,rho_r,rho_b,invrho_r,invrho_b,sharp_c,beta,kapp,tau_diff,sigma &
#ifdef MONOD
	   ,mu_max,Ks &
#endif   
#endif   
       ,ntothfields,ntotphifields,ntotauxfields,ntotlocauxfields,ntotforces &
	   ,hfields_s,phifields_in,phifields_out,auxfields,locauxfields,forces)
      !$acc end host_data
      istat = cudaGetLastError()
      if (istat/= cudaSuccess) then
        if(myrank==0) then
          write(6,*) 'status after at ', __LINE__ ,' of file ', __FILE__ ,' :'
          write(6,*) cudaGetErrorString(istat)
        endif
        call doerror(6,'ERROR in update_phifields_kernel_yplus (launch)')
      endif
      
      istat = cudaDeviceSynchronize()
      if (istat /= cudaSuccess) then
       if(myrank==0) then
         write(6,*) 'sync error at ', __LINE__, ' file ', __FILE__
         write(6,*) cudaGetErrorString(istat)
       endif
       call doerror(6,'ERROR in update_phifields_kernel_yplus (sync)')
      endif
      !$acc wait  
      endif
      
       if(ldodimGridx)then
      !$acc host_data use_device(step,iprobe,jprobe,kprobe,flop,nx,ny,nz,coords,isfluid &
#ifdef TWOCOMPONENT    
       !$acc& ,visc1,visc2,rho_r,rho_b,invrho_r,invrho_b,sharp_c,beta,kapp,tau_diff,sigma &
#ifdef MONOD
	   !$acc& ,mu_max,Ks &
#endif
#endif
       !$acc& ,ntothfields,ntotphifields,ntotauxfields,ntotlocauxfields,ntotforces &
	   !$acc& ,hfields_s,phifields_in,phifields_out,auxfields,locauxfields,forces)
         call update_phifields_kernel_xminus<<<dimGridx, dimBlock>>>(step &
        ,iprobe,jprobe,kprobe,flop &
        ,nx,ny,nz,coords,isfluid &
#ifdef TWOCOMPONENT           
        ,visc2,rho_r,rho_b,invrho_r,invrho_b,sharp_c,beta,kapp,tau_diff,sigma &
#ifdef MONOD
	    ,mu_max,Ks &
#endif   
#endif   
        ,ntothfields,ntotphifields,ntotauxfields,ntotlocauxfields,ntotforces &
	    ,hfields_s,phifields_in,phifields_out,auxfields,locauxfields,forces)
      !$acc end host_data
      istat = cudaGetLastError()
      if (istat/= cudaSuccess) then
        if(myrank==0) then
          write(6,*) 'status after at ', __LINE__ ,' of file ', __FILE__ ,' :'
          write(6,*) cudaGetErrorString(istat)
        endif
        call doerror(6,'ERROR in update_phifields_kernel_xminus (launch)')
      endif
      
      istat = cudaDeviceSynchronize()
      if (istat /= cudaSuccess) then
       if(myrank==0) then
         write(6,*) 'sync error at ', __LINE__, ' file ', __FILE__
         write(6,*) cudaGetErrorString(istat)
       endif
       call doerror(6,'ERROR in update_phifields_kernel_xminus (sync)')
      endif
      !$acc wait  
      endif
      
      if(ldodimGridx .and. ldodimGridx_twice)then
      !$acc host_data use_device(step,iprobe,jprobe,kprobe,flop,nx,ny,nz,coords,isfluid &
#ifdef TWOCOMPONENT    
       !$acc& ,visc1,visc2,rho_r,rho_b,invrho_r,invrho_b,sharp_c,beta,kapp,tau_diff,sigma &
#ifdef MONOD
	   !$acc& ,mu_max,Ks &
#endif
#endif
       !$acc& ,ntothfields,ntotphifields,ntotauxfields,ntotlocauxfields,ntotforces &
	   !$acc& ,hfields_s,phifields_in,phifields_out,auxfields,locauxfields,forces)
        call update_phifields_kernel_xplus<<<dimGridx, dimBlock>>>(step &
       ,iprobe,jprobe,kprobe,flop &
       ,nx,ny,nz,coords,isfluid &
#ifdef TWOCOMPONENT           
       ,visc2,rho_r,rho_b,invrho_r,invrho_b,sharp_c,beta,kapp,tau_diff,sigma &
#ifdef MONOD
	   ,mu_max,Ks &
#endif   
#endif   
       ,ntothfields,ntotphifields,ntotauxfields,ntotlocauxfields,ntotforces &
	   ,hfields_s,phifields_in,phifields_out,auxfields,locauxfields,forces)
      !$acc end host_data
      istat = cudaGetLastError()
      if (istat/= cudaSuccess) then
        if(myrank==0) then
          write(6,*) 'status after at ', __LINE__ ,' of file ', __FILE__ ,' :'
          write(6,*) cudaGetErrorString(istat)
        endif
        call doerror(6,'ERROR in update_phifields_kernel_xplus (launch)')
      endif
      
      istat = cudaDeviceSynchronize()
      if (istat /= cudaSuccess) then
       if(myrank==0) then
         write(6,*) 'sync error at ', __LINE__, ' file ', __FILE__
         write(6,*) cudaGetErrorString(istat)
       endif
       call doerror(6,'ERROR in update_phifields_kernel_xplus (sync)')
      endif
      !$acc wait  
      endif
      

      
#endif
      
      return
      
   endsubroutine update_phifields_ext
   
   ! Boundary condition at internal obstacles, called by bcs_mesoscopic_hfields when BOUNCE_BACK is
   ! not defined: corrects the moments in hfields_out of the fluid nodes next to a solid node
   ! (isfluid = -1), using also the moments in hfields_in.
   subroutine LB_int_boundary_cuda(hfields_in,hfields_out &
#ifdef TWOCOMPONENT	       
     ,phifields_s &
#endif     
     )

      implicit none
      real(kind=strdb), allocatable, dimension(:,:,:,:,:) :: hfields_in,hfields_out
#ifdef TWOCOMPONENT	       
      real(kind=strdb), allocatable, dimension(:,:,:,:,:) :: phifields_s
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
	   ,hfields_in,hfields_out, &
#ifdef TWOCOMPONENT 
	   auxfields,locauxfields, &
#endif   
	   forces)
!$acc end host_data
      
      istat = cudaGetLastError()
      if (istat/= cudaSuccess) then
        if(myrank==0)write(6,*) 'status after at ', __LINE__ ,' of file ', __FILE__ ,' :'
        if(myrank==0)write(6,*) cudaGetErrorString(istat)
        call doerror(6,'ERROR in LB_int_boundary_cuda')
      endif
      istat = cudaDeviceSynchronize
      !$acc wait

      
   end subroutine LB_int_boundary_cuda
   
   ! Phase-field boundary condition at internal obstacles (called by bcs_mesoscopic_phifields).
   subroutine PHI_int_boundary_cuda(hfields_s,phifields_s)

      implicit none
      real(kind=strdb), allocatable, dimension(:,:,:,:,:) :: hfields_s,phifields_s
 
      
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
      
      istat = cudaGetLastError()
      if (istat/= cudaSuccess) then
        if(myrank==0)write(6,*) 'status after at ', __LINE__ ,' of file ', __FILE__ ,' :'
        if(myrank==0)write(6,*) cudaGetErrorString(istat)
        call doerror(6,'ERROR in PHI_int_boundary_cuda')
      endif
      istat = cudaDeviceSynchronize
      !$acc wait        
#endif
      
      return
      
   endsubroutine PHI_int_boundary_cuda 
   
   ! Global conservation of phi with internal obstacles (called by bcs_mesoscopic_phifields):
   ! phi_sum_count_cuda sums phi over the fluid nodes into global_phi_sum and counts in
   ! global_count the nodes with 0.5 < phi < 0.9; apply_lagrangian_phi_cuda adds the correction
   ! corr, computed on the host from these sums, to the same nodes.
   subroutine phi_sum_count_cuda(hfields_s,phifields_s)

      implicit none
      real(kind=strdb), allocatable, dimension(:,:,:,:,:) :: hfields_s,phifields_s
 
      
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
      
      istat = cudaGetLastError()
      if (istat/= cudaSuccess) then
        if(myrank==0)write(6,*) 'status after at ', __LINE__ ,' of file ', __FILE__ ,' :'
        if(myrank==0)write(6,*) cudaGetErrorString(istat)
        call doerror(6,'ERROR in phi_sum_count_cuda')
      endif
      istat = cudaDeviceSynchronize
      !$acc wait        
#endif
      
      return
      
   endsubroutine phi_sum_count_cuda 
  
  subroutine apply_lagrangian_phi_cuda(hfields_s,phifields_s)

      implicit none
      real(kind=strdb), allocatable, dimension(:,:,:,:,:) :: hfields_s,phifields_s
 
      
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
      
      istat = cudaGetLastError()
      if (istat/= cudaSuccess) then
        if(myrank==0)write(6,*) 'status after at ', __LINE__ ,' of file ', __FILE__ ,' :'
        if(myrank==0)write(6,*) cudaGetErrorString(istat)
        call doerror(6,'ERROR in phi_sum_count_cuda')
      endif
      istat = cudaDeviceSynchronize
      !$acc wait        
#endif
      
      return
      
   endsubroutine apply_lagrangian_phi_cuda 
   

endmodule lb_cuda_driver
