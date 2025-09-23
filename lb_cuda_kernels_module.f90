#include "defines.h"
#if !defined(_KERNELCUDA) && !defined(_OPENACC)  
#error "To use this module the macros _KERNELCUDA and _OPENACC should be defined."
#endif


module lb_cuda_kernels

   use vars
   use cudafor
   use mpi_template, only: coords,dostop,doerror,mydev,myrank,nprocs

   implicit none
   
   integer :: istat
   integer, save :: TILE_DIMx=64
   integer, save :: TILE_DIMy=4
   integer, save :: TILE_DIMz=1
   integer, save :: TILE_DIM=16
   integer, constant :: TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,TILE_DIM_d
   integer, constant :: nx_d,ny_d,nz_d
   integer, constant :: lx_d,ly_d,lz_d
   type (dim3) :: dimGrid,dimBlock

contains

   subroutine setup_cuda
      implicit none
      
      istat = cudaSetDevice(mydev)
      if (istat/=0) then
        if(myrank==0)write(6,*) 'status after cudaSetDevice:', cudaGetErrorString(istat)
        call dostop
      endif

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
   
   endsubroutine setup_cuda
   
   subroutine moments_LB_cuda

      implicit none
      
!      if(myrank==0)write(6,*)'step ',step, __LINE__ , __FILE__
      !$acc wait
      istat = cudaDeviceSynchronize

!$acc host_data use_device(flop,nx,ny,nz,coords,isfluid,f &
       !$acc& ,rho,u,v,w,pxx,pyy,pzz,pxy,pxz,pyz,fux,fvy,fwz &
#ifdef TWOCOMPONENT
       !$acc& ,selphi,modgrad,lap_phi,normx,normy,normz &
#ifdef DENSRATIO
       !$acc& ,rhophi &
#endif
#endif    
#ifdef MULTIHIT
       !$acc& ,ABCx,ABCy,ABCz &
#endif 
#ifdef WETTABILITY
       !$acc& ,wettab_r,wettab_b &
#endif  
#ifdef TWOCOMPONENT 
       !$acc& ,visc1,visc2,rho_r,rho_b,invrho_r,invrho_b,beta,kapp,sigma &
#ifdef PHASE_CHANGE 
       !$acc& ,pc_rate,src &
#endif
#endif   
#if defined(ELASTIC_FORCE)
       !$acc& ,lambda_rel,k_elastic,u_ref,v_ref,w_ref &
#endif
       !$acc& ,fx,fy,fz)
       call moments_LB_kernel<<<dimGrid, dimBlock>>>(flop,nx,ny,nz,coords,isfluid,f &
       ,rho,u,v,w,pxx,pyy,pzz,pxy,pxz,pyz,fux,fvy,fwz &
#ifdef TWOCOMPONENT
       ,selphi,modgrad,lap_phi,normx,normy,normz &
#ifdef DENSRATIO
       ,rhophi &
#endif
#endif    
#ifdef MULTIHIT
	   ,ABCx,ABCy,ABCz &
#endif 
#ifdef WETTABILITY
       ,wettab_r,wettab_b &
#endif  
#ifdef TWOCOMPONENT 
       ,visc1,visc2,rho_r,rho_b,invrho_r,invrho_b,beta,kapp,sigma &
#ifdef PHASE_CHANGE 
       ,pc_rate,src &
#endif       
#endif  
#if defined(ELASTIC_FORCE)
       ,lambda_rel,k_elastic,u_ref,v_ref,w_ref &
#endif
       ,fx,fy,fz)
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

   subroutine fused_LB_cuda

      implicit none
      
!      if(myrank==0)write(6,*)'step ',step, __LINE__ , __FILE__
      !$acc wait
      istat = cudaDeviceSynchronize

!$acc host_data use_device(flip,flop,nx,ny,nz,coords,isfluid,f &
       !$acc& ,rho,u,v,w,pxx,pyy,pzz,pxy,pxz,pyz,fux,fvy,fwz &
#ifdef TWOCOMPONENT
       !$acc& ,selphi,lap_phi,normx,normy,normz,arr_x,arr_y,arr_z &
#ifdef DENSRATIO
       !$acc& ,rhophi &
#endif
#endif    
#ifdef MULTIHIT
	   !$acc& ,ABCx,ABCy,ABCz &
#endif 
#ifdef WETTABILITY
       !$acc& ,wettab_r,wettab_b &
#endif  
#ifdef TWOCOMPONENT 
       !$acc& ,visc1,visc2,rho_r,rho_b,invrho_r,invrho_b,sharp_c,beta,kapp,tau_diff,sigma &
#ifdef PHASE_CHANGE 
       !$acc& ,src &
#endif   
#ifdef MONOD
	   !$acc& ,mu_max,Ks &
#endif
#endif   
       !$acc& ,omega,fx,fy,fz)
      call fused_LB_kernel<<<dimGrid, dimBlock>>>(flip,flop,nx,ny,nz,coords,isfluid,f &
       ,rho,u,v,w,pxx,pyy,pzz,pxy,pxz,pyz,fux,fvy,fwz &
#ifdef TWOCOMPONENT
       ,selphi,lap_phi,normx,normy,normz,arr_x,arr_y,arr_z &
#ifdef DENSRATIO
       ,rhophi &
#endif
#endif    
#ifdef MULTIHIT
	   ,ABCx,ABCy,ABCz &
#endif 
#ifdef WETTABILITY
       ,wettab_r,wettab_b &
#endif  
#ifdef TWOCOMPONENT 
       ,visc1,visc2,rho_r,rho_b,invrho_r,invrho_b,sharp_c,beta,kapp,tau_diff,sigma &
#ifdef PHASE_CHANGE 
       ,src &
#endif 
#ifdef MONOD
	   ,mu_max,Ks &
#endif
#endif   
       ,omega,fx,fy,fz)
!$acc end host_data
      istat = cudaDeviceSynchronize
      istat = cudaGetLastError()
      if (istat/= cudaSuccess) then
        if(myrank==0)write(6,*) 'status after at ', __LINE__ ,' of file ', __FILE__ ,' :'
        if(myrank==0)write(6,*) cudaGetErrorString(istat)
        call doerror(6,'ERROR in fused_LB_cuda')
      endif
      !$acc wait      

      
      
      

!!$acc host_data use_device(nx,ny,nz,coords,phi,phi_old)
!      call test_LB_kernel<<<dimGrid, dimBlock>>>(nx,ny,nz,coords,phi,phi_old)
!!$acc end host_data

      
   end subroutine fused_LB_cuda


   !****************************************************************************!


 attributes(global) subroutine test_LB_kernel(nx,ny,nz,coords,phi,phi_old)
      implicit none
      
      integer :: nx,ny,nz
      integer, dimension(3) :: coords
      real(kind=db), dimension(0:nx+1,0:ny+1,0:nz+1) :: phi,phi_old
      
      integer :: i,j,k,gi,gj,gk

      i = (blockIdx%x-1) * TILE_DIMx_d + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy_d + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz_d + threadIdx%z
      
      gi=nx*coords(1)+i
      gj=ny*coords(2)+j
      gk=nz*coords(3)+k
      
      if(gi==8 .and. gj==8 .and. gk==8)then
        write(*,*)'eccomi',phi(i,j,k)
      endif
     if(gi==1) write(*,*)'ciao',i,j,k,phi(i,j,k)
     
     phi_old(i,j,k)=phi(i,j,k)
     
      
 end subroutine test_LB_kernel
 
 attributes(global) subroutine moments_LB_kernel(flop,nx,ny,nz,coords,isfluid,f &
       ,rho,u,v,w,pxx,pyy,pzz,pxy,pxz,pyz,fux,fvy,fwz &
#ifdef TWOCOMPONENT
       ,selphi,modgrad,lap_phi,normx,normy,normz &
#ifdef DENSRATIO
       ,rhophi &
#endif
#endif    
#ifdef MULTIHIT
	   ,ABCx,ABCy,ABCz &
#endif 
#ifdef WETTABILITY
       ,wettab_r,wettab_b &
#endif  
#ifdef TWOCOMPONENT 
       ,visc1,visc2,rho_r,rho_b,invrho_r,invrho_b,beta,kapp,sigma &
#ifdef PHASE_CHANGE 
       ,pc_rate,src &
#endif
#endif   
#if defined(ELASTIC_FORCE)
       ,lambda_rel,k_elastic,u_ref,v_ref,w_ref &
#endif
       ,fx,fy,fz)
 

      implicit none
      
      integer :: flop,nx,ny,nz
      integer, dimension(3) :: coords
      integer(kind=isf), dimension(0:nx+1,0:ny+1,0:nz+1) :: isfluid
      real(kind=db), dimension(0:nx+1,0:ny+1,0:nz+1,0:nlinks) :: f
      real(kind=db), dimension(0:nx+1,0:ny+1,0:nz+1) :: rho,u,v,w,pxx,pyy,pzz,pxy,pxz,pyz
      real(kind=db), dimension(1:nx,1:ny,1:nz) :: fux,fvy,fwz
#ifdef TWOCOMPONENT
#ifdef WENO
      real(kind=db), dimension(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff,2) :: selphi
#else
      real(kind=db), dimension(0:nx+1,0:ny+1,0:nz+1,2) :: selphi
#endif
      real(kind=db), dimension(0:nx+1,0:ny+1,0:nz+1) :: modgrad
      real(kind=db), dimension(1:nx,1:ny,1:nz) :: lap_phi
      real(kind=db), dimension(0:nx+1,0:ny+1,0:nz+1) :: normx,normy,normz
#ifdef DENSRATIO
      real(kind=db), dimension(0:nx+1,0:ny+1,0:nz+1) :: rhophi
#endif
#endif
#ifdef MULTIHIT
	  real(kind=db), dimension(1:nx,1:ny,1:nz) :: ABCx,ABCy,ABCz
#endif
#ifdef WETTABILITY  
      real(kind=db) :: wettab_r,wettab_b  
#endif  
#ifdef TWOCOMPONENT 
      real(kind=db) :: visc1,visc2,rho_r,rho_b,invrho_r,invrho_b,beta,kapp,sigma
#ifdef PHASE_CHANGE 
      real(kind=db) :: pc_rate
      real(kind=db), dimension(1:nx,1:ny,1:nz) :: src
#endif      
#endif           
#if defined(ELASTIC_FORCE)
      real(kind=db) :: lambda_rel,k_elastic
      real(kind=db), dimension(1:nx,1:ny,1:nz) :: u_ref,v_ref,w_ref
#endif
      real(kind=db) :: fx,fy,fz


      real(kind=db) :: fneq1,feq,forcex,forcey,forcez,F_discr,rhophi_loc,uu,udotc
#ifdef TWOCOMPONENT
	  real(kind=db) ::gradfix,gradfiy,gradfiz,wet_loc
#ifdef CSF
      real(kind=db) :: curvature
#endif
#endif
      
#ifdef DENSRATIO
      real(kind=db) :: gradrhox,gradrhoy,gradrhoz,omega_loc, visc_loc,tau_loc
#endif
      
       integer :: i,j,k!,l
!      integer :: gi,gj,gk

      i = (blockIdx%x-1) * TILE_DIMx_d + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy_d + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz_d + threadIdx%z
      
!      gi=nx*coords(1)+i
!      gj=ny*coords(2)+j
!      gk=nz*coords(3)+k
      
!      if(gi==8 .and. gj==8 .and. gk==8)then
!        write(*,*)'eccomi',phi(i,j,k)
!      endif
      
               if (abs(isfluid(i,j,k)) /= 1) return
                  !
				  !pressure
				  rho(i,j,k) = f(i,j,k,0)+f(i,j,k,1)+f(i,j,k,2)+f(i,j,k,3)+f(i,j,k,4)+f(i,j,k,5) &
                     +f(i,j,k,6)+f(i,j,k,7)+f(i,j,k,8)+f(i,j,k,9)+f(i,j,k,10)+f(i,j,k,11) &
                     +f(i,j,k,12)+f(i,j,k,13)+f(i,j,k,14)+f(i,j,k,15)+f(i,j,k,16)+f(i,j,k,17) &
                     +f(i,j,k,18)+f(i,j,k,19)+f(i,j,k,20)+f(i,j,k,21)+f(i,j,k,22)+f(i,j,k,23)+f(i,j,k,24) &
                     +f(i,j,k,25) +f(i,j,k,26)


                  !gi=nx*coords(1)+i
                  !gj=ny*coords(2)+j
                  !gk=nz*coords(3)+k
                  !if(gi==32 .and. gj==32 .and. gk==17 )then
                  !  write(*,*)pzz(i,j,k)
                  !endif
                  !total flux tensor
                 pxx(i,j,k)=f(i,j,k,1)+f(i,j,k,2)+f(i,j,k,7)+f(i,j,k,8) &
                  +f(i,j,k,9)+f(i,j,k,10)+f(i,j,k,15)+f(i,j,k,16)+f(i,j,k,17)+f(i,j,k,18) &
                  +f(i,j,k,19)+f(i,j,k,20)+f(i,j,k,21)+f(i,j,k,22)+f(i,j,k,23)+f(i,j,k,24) &
                  +f(i,j,k,25)+f(i,j,k,26)
                 pyy(i,j,k)=f(i,j,k,3)+f(i,j,k,4)+f(i,j,k,7)+f(i,j,k,8) &
                  +f(i,j,k,9)+f(i,j,k,10)+f(i,j,k,11)+f(i,j,k,12)+f(i,j,k,13)+f(i,j,k,14) &
                  +f(i,j,k,19)+f(i,j,k,20)+f(i,j,k,21)+f(i,j,k,22)+f(i,j,k,23)+f(i,j,k,24) &
                  +f(i,j,k,25)+f(i,j,k,26)
                 pzz(i,j,k)=f(i,j,k,5)+f(i,j,k,6)+f(i,j,k,11)+f(i,j,k,12)+f(i,j,k,13)+f(i,j,k,14) &
                  +f(i,j,k,15)+f(i,j,k,16)+f(i,j,k,17)+f(i,j,k,18) &
                  +f(i,j,k,19)+f(i,j,k,20)+f(i,j,k,21)+f(i,j,k,22)+f(i,j,k,23)+f(i,j,k,24) &
                  +f(i,j,k,25)+f(i,j,k,26)
                 pxy(i,j,k)=(f(i,j,k,7)+f(i,j,k,8)+f(i,j,k,19)+f(i,j,k,20)+f(i,j,k,23)+f(i,j,k,24)) &
                  -(f(i,j,k,9)+f(i,j,k,10)+f(i,j,k,21)+f(i,j,k,22)+f(i,j,k,25)+f(i,j,k,26))
                 pxz(i,j,k)=(f(i,j,k,15)+f(i,j,k,16)+f(i,j,k,19)+f(i,j,k,20)+f(i,j,k,21)+f(i,j,k,22)) &
                  -(f(i,j,k,17)+f(i,j,k,18)+f(i,j,k,23)+f(i,j,k,24)+f(i,j,k,25)+f(i,j,k,26))
                 pyz(i,j,k)=(f(i,j,k,11)+f(i,j,k,12)+f(i,j,k,19)+f(i,j,k,20)+f(i,j,k,25)+f(i,j,k,26)) &
                  -(f(i,j,k,13)+f(i,j,k,14)+f(i,j,k,21)+f(i,j,k,22)+f(i,j,k,23)+f(i,j,k,24))

                          ! gi=nx*coords(1)+i
                          ! gj=ny*coords(2)+j
                          ! gk=nz*coords(3)+k
                          ! if(gi==8 .and. gj==8 .and. gk==18)then
                          !   write(*,*)'m',step,pxx(i,j,k),pyy(i,j,k),pzz(i,j,k)
                          ! endif



                 fux(i,j,k)=0.0_db
				 fvy(i,j,k)=0.0_db
				 fwz(i,j,k)=0.0_db
				 
#ifdef DENSRATIO
                  rhophi_loc = rhophi(i,j,k)
#else
                  rhophi_loc = 1.0_db !rho(i,j,k)
#endif	
#ifdef TWOCOMPONENT	
                  wet_loc = 1.0_db	
#ifdef WETTABILITY

                  if(isfluid(i,j,k).eq.-1) wet_loc = wettab_r*selphi(i,j,k,flop)+ &
                   (1.0_db-selphi(i,j,k,flop))*wettab_b

#endif
#endif				 
#ifdef TWOCOMPONENT		
#ifdef CSF
                   gradfix=normx(i,j,k)*modgrad(i,j,k)
                   gradfiy=normy(i,j,k)*modgrad(i,j,k)
                   gradfiz=normz(i,j,k)*modgrad(i,j,k)
                                  ! modgrad_x=3.0_db*(p1*(modgrad(i+1,j,k)-modgrad(i-1,j,k)) + &
                     ! p2*( (modgrad(i+1,j+1,k)-modgrad(i-1,j-1,k))+(modgrad(i+1,j-1,k)-modgrad(i-1,j+1,k))+(modgrad(i+1,j,k+1)-modgrad(i-1,j,k-1))+(modgrad(i+1,j,k-1)-modgrad(i-1,j,k+1)) )  + &
                     ! p3*((modgrad(i+1,j+1,k+1)-modgrad(i-1,j-1,k-1))+(modgrad(i+1,j-1,k-1)-modgrad(i-1,j+1,k+1))+(modgrad(i+1,j-1,k+1)-modgrad(i-1,j+1,k-1))+(modgrad(i+1,j+1,k-1)-modgrad(i-1,j-1,k+1))))
                                  ! modgrad_y=3.0_db*(p1*(modgrad(i,j+1,k)-modgrad(i,j-1,k)) + &
                     ! p2*((modgrad(i+1,j+1,k)-modgrad(i-1,j-1,k))+(modgrad(i-1,j+1,k)-modgrad(i+1,j-1,k))+(modgrad(i,j+1,k+1)-modgrad(i,j-1,k-1))+(modgrad(i,j+1,k-1)-modgrad(i,j-1,k+1)) ) + &
                     ! p3*((modgrad(i+1,j+1,k+1)-modgrad(i-1,j-1,k-1))+(modgrad(i-1,j+1,k-1)-modgrad(i+1,j-1,k+1))+(modgrad(i+1,j+1,k-1)-modgrad(i-1,j-1,k+1))+(modgrad(i-1,j+1,k+1)-modgrad(i+1,j-1,k-1))))
                                  ! modgrad_z=3.0_db*(p1*(modgrad(i,j,k+1)-modgrad(i,j,k-1)) + &
                     ! p2*((modgrad(i+1,j,k+1)-modgrad(i-1,j,k-1))+(modgrad(i-1,j,k+1)-modgrad(i+1,j,k-1))+(modgrad(i,j+1,k+1)-modgrad(i,j-1,k-1))+(modgrad(i,j-1,k+1)-modgrad(i,j+1,k-1)) ) + &
                     ! p3*((modgrad(i+1,j+1,k+1)-modgrad(i-1,j-1,k-1))+(modgrad(i-1,j-1,k+1)-modgrad(i+1,j+1,k-1))+(modgrad(i+1,j-1,k+1)-modgrad(i-1,j+1,k-1))+(modgrad(i-1,j+1,k+1)-modgrad(i+1,j-1,k-1))))
                                  ! curvature=(1.0_db/(modgrad(i,j,k)+1.0e-9))*(lap_phi(i,j,k) - (normx(i,j,k)*modgrad_x + normy(i,j,k)*modgrad_y + normz(i,j,k)*modgrad_z))
                   curvature=3.0_db*(p1*(normx(i+1,j,k)-normx(i-1,j,k)) + &
                    p2*( (normx(i+1,j+1,k)-normx(i-1,j-1,k))+(normx(i+1,j-1,k)-normx(i-1,j+1,k))+(normx(i+1,j,k+1)-normx(i-1,j,k-1))+(normx(i+1,j,k-1)-normx(i-1,j,k+1)) )  + &
                    p3*((normx(i+1,j+1,k+1)-normx(i-1,j-1,k-1))+(normx(i+1,j-1,k-1)-normx(i-1,j+1,k+1))+(normx(i+1,j-1,k+1)-normx(i-1,j+1,k-1))+(normx(i+1,j+1,k-1)-normx(i-1,j-1,k+1)))) + &
                    3.0_db*(p1*(normy(i,j+1,k)-normy(i,j-1,k)) + &
                    p2*((normy(i+1,j+1,k)-normy(i-1,j-1,k))+(normy(i-1,j+1,k)-normy(i+1,j-1,k))+(normy(i,j+1,k+1)-normy(i,j-1,k-1))+(normy(i,j+1,k-1)-normy(i,j-1,k+1)) ) + &
                    p3*((normy(i+1,j+1,k+1)-normy(i-1,j-1,k-1))+(normy(i-1,j+1,k-1)-normy(i+1,j-1,k+1))+(normy(i+1,j+1,k-1)-normy(i-1,j-1,k+1))+(normy(i-1,j+1,k+1)-normy(i+1,j-1,k-1)))) +&
                    3.0_db*(p1*(normz(i,j,k+1)-normz(i,j,k-1)) + &
                    p2*((normz(i+1,j,k+1)-normz(i-1,j,k-1))+(normz(i-1,j,k+1)-normz(i+1,j,k-1))+(normz(i,j+1,k+1)-normz(i,j-1,k-1))+(normz(i,j-1,k+1)-normz(i,j+1,k-1)) ) + &
                    p3*((normz(i+1,j+1,k+1)-normz(i-1,j-1,k-1))+(normz(i-1,j-1,k+1)-normz(i+1,j+1,k-1))+(normz(i+1,j-1,k+1)-normz(i-1,j+1,k-1))+(normz(i-1,j+1,k+1)-normz(i+1,j-1,k-1))))

                   fux(i,j,k)=-sigma*eps1*curvature*gradfix*modgrad(i,j,k)
				   fvy(i,j,k)=-sigma*eps1*curvature*gradfiy*modgrad(i,j,k)
				   fwz(i,j,k)=-sigma*eps1*curvature*gradfiz*modgrad(i,j,k)
#endif
#ifdef JACQMIN
				   gradfix=normx(i,j,k)*modgrad(i,j,k)
				   gradfiy=normy(i,j,k)*modgrad(i,j,k)
				   gradfiz=normz(i,j,k)*modgrad(i,j,k)

                   fux(i,j,k)=(4.0_db*(wet_loc*beta)*(selphi(i,j,k,flop))* &
                    (selphi(i,j,k,flop)-1.0_db)*(selphi(i,j,k,flop)-0.5_db) - &
                    (wet_loc*kapp)*lap_phi(i,j,k))*gradfix
				   fvy(i,j,k)=(4.0_db*(wet_loc*beta)*(selphi(i,j,k,flop))* &
				    (selphi(i,j,k,flop)-1.0_db)*(selphi(i,j,k,flop)-0.5_db) - &
				    (wet_loc*kapp)*lap_phi(i,j,k))*gradfiy
				   fwz(i,j,k)=(4.0_db*(wet_loc*beta)*(selphi(i,j,k,flop))* &
				    (selphi(i,j,k,flop)-1.0_db)*(selphi(i,j,k,flop)-0.5_db) - &
				    (wet_loc*kapp)*lap_phi(i,j,k))*gradfiz


#endif
#endif

#if defined(MULTIHIT)   
                  fux(i,j,k)=fux(i,j,k) + rhophi_loc*ABCx(i,j,k) !+ AAA*sin(k_zero*gk) + AAA*sin(k_zero*gj)  
				  fvy(i,j,k)=fvy(i,j,k) + rhophi_loc*ABCy(i,j,k) !+ AAA*sin(k_zero*gi) + AAA*sin(k_zero*gk)
				  fwz(i,j,k)=fwz(i,j,k) + rhophi_loc*ABCz(i,j,k) !+ AAA*sin(k_zero*gj) + AAA*sin(k_zero*gi) 	
#endif
#if defined(PLUG_FLOW)   
                  
				  fwz(i,j,k)=fwz(i,j,k) + selphi(i,j,k,flop)*fz ! if mnulticomponent/phase (rhophi_loc-rho_r or rho_b)*fz 	
#endif


#ifdef DENSRATIO				  
				  ! pressure and viscous forces
				  
				  gradrhox=(rho_r-rho_b)*gradfix
				  gradrhoy=(rho_r-rho_b)*gradfiy
				  gradrhoz=(rho_r-rho_b)*gradfiz

                  fux(i,j,k)=fux(i,j,k) - rho(i,j,k)*cssq*gradrhox   !+ (rhophi_loc-(rho_r+rho_b)*HALF)*fx
                  fvy(i,j,k)=fvy(i,j,k) - rho(i,j,k)*cssq*gradrhoy   !+ (rhophi_loc-(rho_r+rho_b)*HALF)*fy
				  fwz(i,j,k)=fwz(i,j,k) - rho(i,j,k)*cssq*gradrhoz   !+ (rhophi_loc-(rho_r+rho_b)*HALF)*fz
				  !! from this point I compute the force terms that depend on the velocity
				  !! these terms should be not included in force arrays since they must be computed with the updated velocity
				  !! at the end of this subroutine
#endif                                  
				  forcex=fux(i,j,k) 
                  forcey=fvy(i,j,k) 
				  forcez=fwz(i,j,k) 
				  
#if defined(ELASTIC_FORCE)

                  forcex=forcex + rhophi_loc*(u(i,j,k) - u_ref(i,j,k))*k_elastic*lambda_rel*selphi(i,j,k,flop) !+   
				  forcey=forcey + rhophi_loc*(v(i,j,k) - v_ref(i,j,k))*k_elastic*lambda_rel*selphi(i,j,k,flop) !+ 
				  forcez=forcez + rhophi_loc*(w(i,j,k) - w_ref(i,j,k))*k_elastic*lambda_rel*selphi(i,j,k,flop) + rhophi_loc*fz !+  	
#endif
#ifdef DENSRATIO 
#ifdef PHASE_CHANGE
!!!aritra
!!!The force term going to LB following Fakhari/Rahimian paper
!!!Phase-change modeling based on a novel conservative phase-field method - Eq. 41
				  src(i,j,k) = 6.0_db*selphi(i,j,k,flop)*(1.0_db - selphi(i,j,k,flop))*pc_rate
                  forcex=forcex+rhophi_loc*u(i,j,k)*src(i,j,k)*(invrho_b-invrho_r)
				  forcey=forcey+rhophi_loc*v(i,j,k)*src(i,j,k)*(invrho_b-invrho_r)
				  forcez=forcez+rhophi_loc*w(i,j,k)*src(i,j,k)*(invrho_b-invrho_r)
#endif				  

                  visc_loc=(rho_r*visc1*selphi(i,j,k,flop)+ &
                   (1.0_db-selphi(i,j,k,flop))*visc2*rho_b)/rhophi_loc

				  !visc_loc=(visc1*selphi(i,j,k,flop)+(1.0_db-selphi(i,j,k,flop))*visc2)
                  !visc_loc=(visc1*selphi(i,j,k,flop)+(1.0_db-selphi(i,j,k,flop))*visc2)
                  tau_loc=(visc_loc/cssq + 0.5_db) !è una tau
				  !omega_loc=1.0_db/(1.0_db/tau2 + (selphi(i,j,k,flop))*((1.0_db/tau1)-(1.0_db/tau2))) !tau
				  !visc_loc=(omega_loc-0.5_db)*cssq
				  !note that the non-equilibrium flux tensor is computer removing on the fly the equilibrium part with the old velocities
				  forcex=forcex - (visc_loc/(tau_loc*cssq))* &
				   ((pxx(i,j,k)-(u(i,j,k)*u(i,j,k)+rho(i,j,k)*cssq))*gradrhox + &
				   (pxy(i,j,k)-(u(i,j,k)*v(i,j,k)))*gradrhoy + (pxz(i,j,k)-(u(i,j,k)*w(i,j,k)))*gradrhoz)
				  forcey=forcey - (visc_loc/(tau_loc*cssq))* &
				   ((pyy(i,j,k)-(v(i,j,k)*v(i,j,k)+rho(i,j,k)*cssq))*gradrhoy + &
				   (pxy(i,j,k)-(u(i,j,k)*v(i,j,k)))*gradrhox + (pyz(i,j,k)-(v(i,j,k)*w(i,j,k)))*gradrhoz)
				  forcez=forcez - (visc_loc/(tau_loc*cssq))* &
				   ((pzz(i,j,k)-(w(i,j,k)*w(i,j,k)+rho(i,j,k)*cssq))*gradrhoz + &
				   (pxz(i,j,k)-(u(i,j,k)*w(i,j,k)))*gradrhox + (pyz(i,j,k)-(v(i,j,k)*w(i,j,k)))*gradrhoy)
#endif	
                  !I compute the new velocities
				  u(i,j,k) = ((f(i,j,k,1)+f(i,j,k,7)+f(i,j,k,9)+f(i,j,k,15)+f(i,j,k,18)+f(i,j,k,19)+f(i,j,k,21)+f(i,j,k,24)+f(i,j,k,25)) &
                     -(f(i,j,k,2)+f(i,j,k,8)+f(i,j,k,10)+f(i,j,k,16)+f(i,j,k,17)+f(i,j,k,20)+f(i,j,k,22)+f(i,j,k,23)+f(i,j,k,26)))

                  v(i,j,k) = ((f(i,j,k,3)+f(i,j,k,7)+f(i,j,k,10)+f(i,j,k,11)+f(i,j,k,13)+f(i,j,k,19)+f(i,j,k,22)+f(i,j,k,24)+f(i,j,k,26)) &
                     -(f(i,j,k,4)+f(i,j,k,8)+f(i,j,k,9)+f(i,j,k,12)+f(i,j,k,14)+f(i,j,k,20)+f(i,j,k,21)+f(i,j,k,23)+f(i,j,k,25)))

                  w(i,j,k) = ((f(i,j,k,5)+f(i,j,k,11)+f(i,j,k,14)+f(i,j,k,15)+f(i,j,k,17)+f(i,j,k,19)+f(i,j,k,21)+f(i,j,k,23)+f(i,j,k,26)) &
                     -(f(i,j,k,6)+f(i,j,k,12)+f(i,j,k,13)+f(i,j,k,16)+f(i,j,k,18)+f(i,j,k,20)+f(i,j,k,22)+f(i,j,k,24)+f(i,j,k,25)))
					 

                  u(i,j,k) = u(i,j,k) + 0.5_db*forcex/rhophi_loc
                  v(i,j,k) = v(i,j,k) + 0.5_db*forcey/rhophi_loc
                  w(i,j,k) = w(i,j,k) + 0.5_db*forcez/rhophi_loc
			  			  
                  !subtracting the eq flux tensor with new velocities and adding half force term added as macroscopic formulation (obtained from wolfram)
                  pxx(i,j,k)=pxx(i,j,k)-(u(i,j,k)*u(i,j,k)+rho(i,j,k)*cssq -0.5_db*(2.0_db*forcex*u(i,j,k)/rhophi_loc))
				  pyy(i,j,k)=pyy(i,j,k)-(v(i,j,k)*v(i,j,k)+rho(i,j,k)*cssq -0.5_db*(2.0_db*forcey*v(i,j,k)/rhophi_loc))
				  pzz(i,j,k)=pzz(i,j,k)-(w(i,j,k)*w(i,j,k)+rho(i,j,k)*cssq -0.5_db*(2.0_db*forcez*w(i,j,k)/rhophi_loc))
				  pxy(i,j,k)=pxy(i,j,k)-(u(i,j,k)*v(i,j,k) -0.5_db*((forcex*v(i,j,k)+forcey*u(i,j,k))/rhophi_loc))
				  pxz(i,j,k)=pxz(i,j,k)-(u(i,j,k)*w(i,j,k) -0.5_db*((forcex*w(i,j,k)+forcez*u(i,j,k))/rhophi_loc))
				  pyz(i,j,k)=pyz(i,j,k)-(v(i,j,k)*w(i,j,k) -0.5_db*((forcey*w(i,j,k)+forcez*v(i,j,k))/rhophi_loc))

                  
                  !!!!!! now I compute the force terms depending on the updated velocities and are stored in force arrays
#ifdef PHASE_CHANGE
                  fux(i,j,k)=fux(i,j,k)+rhophi_loc*u(i,j,k)*src(i,j,k)*(invrho_b-invrho_r)
				  fvy(i,j,k)=fvy(i,j,k)+rhophi_loc*v(i,j,k)*src(i,j,k)*(invrho_b-invrho_r)
				  fwz(i,j,k)=fwz(i,j,k)+rhophi_loc*w(i,j,k)*src(i,j,k)*(invrho_b-invrho_r)
#endif	
#if defined(ELASTIC_FORCE)
				  u_ref(i,j,k) = u_ref(i,j,k) + lambda_rel*(u(i,j,k) - u_ref(i,j,k))
				  v_ref(i,j,k) = v_ref(i,j,k) + lambda_rel*(v(i,j,k) - v_ref(i,j,k))
				  w_ref(i,j,k) = w_ref(i,j,k) + lambda_rel*(w(i,j,k) - w_ref(i,j,k))
                  fux(i,j,k)=fux(i,j,k) + rhophi_loc*(u(i,j,k) - u_ref(i,j,k))*k_elastic*lambda_rel*selphi(i,j,k,flop) !+ 
				  fvy(i,j,k)=fvy(i,j,k) + rhophi_loc*(v(i,j,k) - v_ref(i,j,k))*k_elastic*lambda_rel*selphi(i,j,k,flop) !+ 
				  fwz(i,j,k)=fwz(i,j,k) + rhophi_loc*(w(i,j,k) - w_ref(i,j,k))*k_elastic*lambda_rel*selphi(i,j,k,flop)+rhophi_loc*fz  
#endif
#if defined(DENSRATIO)			  
				  fux(i,j,k)=fux(i,j,k) - (visc_loc/(tau_loc*cssq))*(pxx(i,j,k)*gradrhox + &
				   pxy(i,j,k)*gradrhoy + pxz(i,j,k)*gradrhoz)
                  fvy(i,j,k)=fvy(i,j,k) - (visc_loc/(tau_loc*cssq))*(pyy(i,j,k)*gradrhoy + &
                   pxy(i,j,k)*gradrhox + pyz(i,j,k)*gradrhoz)
                  fwz(i,j,k)=fwz(i,j,k) - (visc_loc/(tau_loc*cssq))*(pzz(i,j,k)*gradrhoz + &
                   pxz(i,j,k)*gradrhox + pyz(i,j,k)*gradrhoy)
#endif           


   endsubroutine moments_LB_kernel  


 attributes(global) subroutine fused_LB_kernel(flip,flop,nx,ny,nz,coords,isfluid,f &
       ,rho,u,v,w,pxx,pyy,pzz,pxy,pxz,pyz,fux,fvy,fwz &
#ifdef TWOCOMPONENT
       ,selphi,lap_phi,normx,normy,normz,arr_x,arr_y,arr_z &
#ifdef DENSRATIO
       ,rhophi &
#endif
#endif    
#ifdef MULTIHIT
	   ,ABCx,ABCy,ABCz &
#endif 
#ifdef WETTABILITY
       ,wettab_r,wettab_b &
#endif  
#ifdef TWOCOMPONENT 
       ,visc1,visc2,rho_r,rho_b,invrho_r,invrho_b,sharp_c,beta,kapp,tau_diff,sigma &
#ifdef PHASE_CHANGE 
       ,src &
#endif 
#ifdef MONOD	
	   ,mu_max,Ks &
#endif
#endif   
       ,omega,fx,fy,fz)

      implicit none
      
      integer :: flip,flop,nx,ny,nz
      integer, dimension(3) :: coords
      integer(kind=isf), dimension(0:nx+1,0:ny+1,0:nz+1) :: isfluid
      real(kind=db), dimension(0:nx+1,0:ny+1,0:nz+1,0:nlinks) :: f
      real(kind=db), dimension(0:nx+1,0:ny+1,0:nz+1) :: rho,u,v,w,pxx,pyy,pzz,pxy,pxz,pyz
      real(kind=db), dimension(1:nx,1:ny,1:nz) :: fux,fvy,fwz
#ifdef TWOCOMPONENT
#ifdef WENO
      real(kind=db), dimension(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff,2) :: selphi
#else
      real(kind=db), dimension(0:nx+1,0:ny+1,0:nz+1,2) :: selphi
#endif
      real(kind=db), dimension(1:nx,1:ny,1:nz) :: lap_phi
      real(kind=db), dimension(0:nx+1,0:ny+1,0:nz+1) :: normx,normy,normz,arr_x,arr_y,arr_z
#ifdef DENSRATIO
      real(kind=db), dimension(0:nx+1,0:ny+1,0:nz+1) :: rhophi
#endif
#endif
#ifdef MULTIHIT
	  real(kind=db), dimension(1:nx,1:ny,1:nz) :: ABCx,ABCy,ABCz
#endif
#ifdef WETTABILITY  
      real(kind=db) :: wettab_r,wettab_b  
#endif  
#ifdef TWOCOMPONENT 
      real(kind=db) :: visc1,visc2,rho_r,rho_b,invrho_r,invrho_b,sharp_c,beta,kapp,tau_diff,sigma
#ifdef PHASE_CHANGE
      real(kind=db), dimension(1:nx,1:ny,1:nz) :: src
#endif   
#ifdef MONOD
      real(kind=db) :: mu_max,Ks
#endif
#endif           
      real(kind=db) :: omega,fx,fy,fz

  
      real(kind=db) :: F_discr,fneq1,feq,forcex,forcey,forcez
#ifdef TWOCOMPONENT
      real(kind=db) :: tau_loc,gradfix,gradfiy,gradfiz
      real(kind=db) :: gradrhox,gradrhoy,gradrhoz,wet_loc,visc_loc
#ifdef MONOD
	  real(kind=db) :: S_mono
#endif
#ifdef CSF
      real(kind=db) :: curvature
#endif
#endif
      real(kind=db) :: omega_loc, rhophi_loc
#ifdef WENO
      real(kind=db) :: fm2,fm1,f0,fp1,fp2,a
      real(kind=db) :: beta1,beta2,beta3
      real(kind=db) :: we1,we2,we3,sum_we
#endif
      integer :: i,j,k!,l
!      integer :: gi,gj,gk

      i = (blockIdx%x-1) * TILE_DIMx_d + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy_d + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz_d + threadIdx%z
      
!      gi=nx*coords(1)+i
!      gj=ny*coords(2)+j
!      gk=nz*coords(3)+k
      
!      if(gi==8 .and. gj==8 .and. gk==8)then
!        write(*,*)'eccomi',phi(i,j,k)
!      endif

               if (abs(isfluid(i,j,k)) /= 1) return
               
               
#ifdef DENSRATIO
                  rhophi_loc = rhophi(i,j,k)
#else
                  rhophi_loc = 1.0_db! rho(i,j,k)
#endif
	  
				  forcex=fux(i,j,k)
				  forcey=fvy(i,j,k)
				  forcez=fwz(i,j,k)

                           !gi=nx*coords(1)+i
                           !gj=ny*coords(2)+j
                           !gk=nz*coords(3)+k
                           !if(gi==8 .and. gj==8 .and. gk==18)then
                           !  write(*,*)'c',step,pxx(i,j,k),pyy(i,j,k),pzz(i,j,k) 
                           !endif

#ifdef TWOCOMPONENT

                  visc_loc=(rho_r*visc1*selphi(i,j,k,flip)+(1.0_db-selphi(i,j,k,flip))*visc2*rho_b)/rhophi_loc  

				  !visc_loc=(visc1*phi(i,j,k)+(1.0_db-phi(i,j,k))*visc2)
                  tau_loc=(visc_loc/cssq + 0.5_db) !è una tau
				  !omega_loc=1.0_db/(1.0_db/tau2 + (phi(i,j,k))*((1.0_db/tau1)-(1.0_db/tau2))) !tau
				  !visc_loc=(omega_loc-0.5_db)*cssq
                  omega_loc=1.0_db/tau_loc !è una omega
				  !visc_loc=(tau_loc-0.5_db)*cssq
				  !omega_loc=1.0_db/(1.0_db/tau2 + (phi(i,j,k))*((1.0_db/tau1)-(1.0_db/tau2))) !tau
				  !visc_loc=(omega_loc-0.5_db)*cssq
#else
				  omega_loc=omega
#endif


!!!!!!!!!!!!!!!!!!!!!!!!!!0
#ifdef SECOND_ORDER
			      feq=(4.0_db*(2.0_db*rho(i,j,k) - 3.0_db &
			       *(u(i,j,k)**2.0_db + v(i,j,k)**2.0_db + w(i,j,k)**2.0_db)))/27.0_db
!0
				  fneq1=(-3.0_db*(pxx(i,j,k) + pyy(i,j,k) + pzz(i,j,k)))/2.0_db
#endif
#ifdef FORTH_ORDER
!0
			      feq=(8.0_db*rho(i,j,k) + 6.0_db*(u(i,j,k)**2.0_db*( &
			       -2.0_db + 3.0_db*v(i,j,k)**2.0_db) + 3.0_db*(u(i,j,k)**2.0_db &
			       + v(i,j,k)**2.0_db)*w(i,j,k)**2.0_db - 2.0_db*(v(i,j,k)**2.0_db &
			       + w(i,j,k)**2.0_db)))/27.0_db
!0
				  fneq1=(3.0_db*(-2.0_db*pxx(i,j,k) - 2.0_db*pzz(i,j,k) &
				   + pyy(i,j,k)*(-2.0_db + 3.0_db*u(i,j,k)**2.0_db) &
				   + 6.0_db*pxx(i,j,k)*v(i,j,k)**2.0_db + 3.0_db*pzz(i,j,k)*v(i,j,k)**2.0_db &
				   + 18.0_db*pyz(i,j,k)*v(i,j,k)*w(i,j,k) + 6.0_db*(pxx(i,j,k) &
				   + pyy(i,j,k))*w(i,j,k)**2.0_db + 3.0_db*u(i,j,k)*(pzz(i,j,k)*u(i,j,k) &
				   + 6.0_db*pxy(i,j,k)*v(i,j,k) + 6.0_db*pxz(i,j,k)*w(i,j,k))))/4.0_db
#endif
#ifdef SIXTH_ORDER
!0
			      feq=(8.0_db*rho(i,j,k) - 3.0_db*(4.0_db*w(i,j,k)**2.0_db &
			       + v(i,j,k)**2.0_db*(4.0_db - 6.0_db*w(i,j,k)**2.0_db) &
			       + u(i,j,k)**2.0_db*(-2.0_db + 3.0_db*v(i,j,k)**2.0_db)*(-2.0_db &
			       + 3.0_db*w(i,j,k)**2.0_db)))/27.0_db
!0
#ifdef GHOSTONE   
				  fneq1=(-3.0_db*(pxx(i,j,k) + pyy(i,j,k) + pzz(i,j,k)))/2.0_db
#else 
				  fneq1=(3.0_db*(2.0_db*(pyy(i,j,k) + pzz(i,j,k))*(-2.0_db &
				   + 3.0_db*u(i,j,k)**2.0_db) + 36.0_db*pxy(i,j,k)*u(i,j,k)*v(i,j,k) &
				   + 3.0_db*pzz(i,j,k)*(2.0_db - 3.0_db*u(i,j,k)**2.0_db)*v(i,j,k)**2.0_db &
				   + 12.0_db*w(i,j,k)*(3.0_db*pyz(i,j,k)*v(i,j,k) + pyy(i,j,k)*w(i,j,k)) &
				   - 9.0_db*u(i,j,k)*w(i,j,k)*(6.0_db*pyz(i,j,k)*u(i,j,k)*v(i,j,k) &
				   + pxz(i,j,k)*(-4.0_db + 6.0_db*v(i,j,k)**2.0_db) &
				   + 3.0_db*pyy(i,j,k)*u(i,j,k)*w(i,j,k) &
				   + 14.0_db*pxy(i,j,k)*v(i,j,k)*w(i,j,k)) - 4.0_db*pxx(i,j,k)*(-1.0_db &
				   + 3.0_db*v(i,j,k)**2.0_db)*(-1.0_db + 3.0_db*w(i,j,k)**2.0_db)))/8.0_db
#endif
#endif
!0
#ifdef PHASE_CHANGE
                  feq=feq+p0*src(i,j,k)*(invrho_b-invrho_r)
#endif
				  F_discr=(-8.0_db*(forcex*u(i,j,k) + forcey*v(i,j,k) &
				   + forcez*w(i,j,k)))/(9.0_db*rhophi_loc)
                  f(i,j,k,0)=feq + (1.0_db-omega_loc)*fneq1*p0 + 0.5_db*(F_discr)
                  !gi=nx*coords(1)+i
                  !gj=ny*coords(2)+j
                  !gk=nz*coords(3)+k
                  !if(gi==32 .and. gj==32 .and. gk==17 )then
                   ! write(*,*)(pzz(i,j,k) )
                  !endif

!!!!!!!!!!!!!!!!!!!!!!!!!!1
#ifdef SECOND_ORDER
			      feq=(2.0_db*rho(i,j,k) + 6.0_db*u(i,j,k) &
			       *(1.0_db + u(i,j,k)) - 3.0_db*v(i,j,k)**2.0_db &
			       - 3.0_db*w(i,j,k)**2.0_db)/27.0_db
!1
				  fneq1=(3.0_db*(2.0_db*pxx(i,j,k) - pyy(i,j,k) - pzz(i,j,k)))/2.0_db
#endif
#ifdef FORTH_ORDER
!1
			      feq=(4.0_db*rho(i,j,k) + 3.0_db*(4.0_db*u(i,j,k)&
			       *(1.0_db + u(i,j,k)) - 2.0_db*(1.0_db + 3.0_db*u(i,j,k) &
			       *(1.0_db + u(i,j,k)))*v(i,j,k)**2.0_db - (2.0_db &
			       + 6.0_db*u(i,j,k)*(1.0_db + u(i,j,k)) - 3.0_db*v(i,j,k)**2.0_db) &
			       *w(i,j,k)**2.0_db))/54.0_db
!1
				  fneq1=(9.0_db*pzz(i,j,k)*v(i,j,k)**2.0_db &
				   + 18.0_db*w(i,j,k)*(3.0_db*pyz(i,j,k)*v(i,j,k) + pyy(i,j,k)*w(i,j,k)) &
				   - 12.0_db*pxx(i,j,k)*(-1.0_db + 3.0_db*v(i,j,k)**2.0_db + 3.0_db*w(i,j,k)**2.0_db) &
				   - 6.0_db*(pyy(i,j,k) + pzz(i,j,k) + 3.0_db*pyy(i,j,k)*u(i,j,k)*(1.0_db &
				   + u(i,j,k)) + 3.0_db*pzz(i,j,k)*u(i,j,k)*(1.0_db + u(i,j,k)) &
				   + 6.0_db*(1.0_db + 3.0_db*u(i,j,k))*(pxy(i,j,k)*v(i,j,k) &
				   + pxz(i,j,k)*w(i,j,k))))/4.0_db
#endif
#ifdef SIXTH_ORDER
!1
			      feq=(4.0_db*rho(i,j,k) + 3.0_db*(4.0_db*u(i,j,k)*(1.0_db &
			       + u(i,j,k)) - 2.0_db*(1.0_db + 3.0_db*u(i,j,k)*(1.0_db &
			       + u(i,j,k)))*v(i,j,k)**2.0_db + (1.0_db + 3.0_db*u(i,j,k)*(1.0_db &
			       + u(i,j,k)))*(-2.0_db + 3.0_db*v(i,j,k)**2.0_db)*w(i,j,k)**2.0_db))/54.0_db
!1
#ifdef GHOSTONE   
				  fneq1=(3.0_db*(2.0_db*pxx(i,j,k) - pyy(i,j,k) - pzz(i,j,k)))/2.0_db
#else 
				  fneq1=(12.0_db*pxx(i,j,k) - 6.0_db*(pyy(i,j,k) + pzz(i,j,k))*(1.0_db &
				   + 3.0_db*u(i,j,k)*(1.0_db + u(i,j,k))) - 36.0_db*pxy(i,j,k)*(1.0_db &
				   + 3.0_db*u(i,j,k))*v(i,j,k) + 9.0_db*(-4.0_db*pxx(i,j,k) + pzz(i,j,k) &
				   + 3.0_db*pzz(i,j,k)*u(i,j,k)*(1.0_db + u(i,j,k)))*v(i,j,k)**2.0_db &
				   + 9.0_db*(-4.0_db*pxz(i,j,k)*(1.0_db + 3.0_db*u(i,j,k)) &
				   + 6.0_db*pyz(i,j,k)*(1.0_db + 3.0_db*u(i,j,k)*(1.0_db + u(i,j,k)))*v(i,j,k) &
				   + 9.0_db*pxz(i,j,k)*(1.0_db + 2.0_db*u(i,j,k))*v(i,j,k)**2.0_db)*w(i,j,k) &
				   + 9.0_db*(pyy(i,j,k)*(2.0_db + 9.0_db*u(i,j,k)*(1.0_db + u(i,j,k))) &
				   + 6.0_db*pxy(i,j,k)*(3.0_db + 7.0_db*u(i,j,k))*v(i,j,k) &
				   + 4.0_db*pxx(i,j,k)*(-1.0_db + 3.0_db*v(i,j,k)**2.0_db))*w(i,j,k)**2.0_db)/4.0_db
#endif
#endif
!1
#ifdef PHASE_CHANGE
                  feq=feq+p1*src(i,j,k)*(invrho_b-invrho_r)
#endif
				  F_discr=(2.0_db*(forcex + 2.0_db*forcex*u(i,j,k) - forcey*v(i,j,k) &
				   - forcez*w(i,j,k)))/(9.0_db*rhophi_loc)
                  f(i+1,j,k,1)=feq + (1.0_db-omega_loc)*fneq1*p1 + 0.5_db*(F_discr)


!!!!!!!!!!!!!!!!!!!!!!!!!!2
#ifdef SECOND_ORDER
			      feq=(2.0_db*rho(i,j,k) + 6.0_db*(-1.0_db &
			       + u(i,j,k))*u(i,j,k) - 3.0_db*v(i,j,k)**2.0_db &
			       - 3.0_db*w(i,j,k)**2.0_db)/27.0_db
!2
				  fneq1=(3.0_db*(2.0_db*pxx(i,j,k) - pyy(i,j,k) - pzz(i,j,k)))/2.0_db
#endif
#ifdef FORTH_ORDER
!2
			      feq=(4.0_db*rho(i,j,k) + 3.0_db*(4.0_db*( &
			       -1.0_db + u(i,j,k))*u(i,j,k) - 2.0_db*(1.0_db + 3.0_db*( &
			       -1.0_db + u(i,j,k))*u(i,j,k))*v(i,j,k)**2.0_db + (-2.0_db &
			       - 6.0_db*(-1.0_db + u(i,j,k))*u(i,j,k) + 3.0_db*v(i,j,k)**2.0_db) &
			       *w(i,j,k)**2.0_db))/54.0_db
!2
				  fneq1=(3.0_db*(3.0_db*pzz(i,j,k)*v(i,j,k)**2.0_db & 
				   + 6.0_db*w(i,j,k)*(3.0_db*pyz(i,j,k)*v(i,j,k) &
				   + pyy(i,j,k)*w(i,j,k)) - 4.0_db*pxx(i,j,k)*(-1.0_db &
				   + 3.0_db*v(i,j,k)**2.0_db + 3.0_db*w(i,j,k)**2.0_db) &
				   - 2.0_db*(pyy(i,j,k) + pzz(i,j,k) + 3.0_db*pyy(i,j,k)*(-1.0_db &
				   + u(i,j,k))*u(i,j,k) + 3.0_db*pzz(i,j,k)*(-1.0_db + u(i,j,k))*u(i,j,k) &
				   + 6.0_db*(-1.0_db + 3.0_db*u(i,j,k))*(pxy(i,j,k)*v(i,j,k) &
				   + pxz(i,j,k)*w(i,j,k)))))/4.0_db
#endif
#ifdef SIXTH_ORDER
!2
			      feq=(4.0_db*rho(i,j,k) + 3.0_db*(4.0_db*(-1.0_db &
			       + u(i,j,k))*u(i,j,k) - 2.0_db*(1.0_db + 3.0_db*(-1.0_db &
			       + u(i,j,k))*u(i,j,k))*v(i,j,k)**2.0_db + (1.0_db &
			       + 3.0_db*(-1.0_db + u(i,j,k))*u(i,j,k))*(-2.0_db &
			       + 3.0_db*v(i,j,k)**2.0_db)*w(i,j,k)**2.0_db))/54.0_db
!2
#ifdef GHOSTONE   
				  fneq1=(3.0_db*(2.0_db*pxx(i,j,k) - pyy(i,j,k) - pzz(i,j,k)))/2.0_db
#else 
				  fneq1=(9.0_db*pzz(i,j,k)*(1.0_db + 3.0_db*(-1.0_db &
				   + u(i,j,k))*u(i,j,k))*v(i,j,k)**2.0_db - 6.0_db*(pyy(i,j,k) &
				   + pzz(i,j,k) + 3.0_db*pyy(i,j,k)*(-1.0_db + u(i,j,k))*u(i,j,k) &
				   + 3.0_db*pzz(i,j,k)*(-1.0_db + u(i,j,k))*u(i,j,k) &
				   + 6.0_db*pxy(i,j,k)*(-1.0_db + 3.0_db*u(i,j,k))*v(i,j,k) &
				   - 6.0_db*pxz(i,j,k)*w(i,j,k)) + 9.0_db*w(i,j,k)*( &
				   -12.0_db*pxz(i,j,k)*u(i,j,k) + 6.0_db*pyz(i,j,k)*(1.0_db &
				   + 3.0_db*(-1.0_db + u(i,j,k))*u(i,j,k))*v(i,j,k) &
				   + 9.0_db*pxz(i,j,k)*(-1.0_db + 2.0_db*u(i,j,k))*v(i,j,k)**2.0_db &
				   + (pyy(i,j,k)*(2.0_db + 9.0_db*(-1.0_db + u(i,j,k))*u(i,j,k)) &
				   + 6.0_db*pxy(i,j,k)*(-3.0_db + 7.0_db*u(i,j,k))*v(i,j,k))*w(i,j,k)) &
				   + 12.0_db*pxx(i,j,k)*(-1.0_db + 3.0_db*v(i,j,k)**2.0_db)*( &
				   -1.0_db + 3.0_db*w(i,j,k)**2.0_db))/4.0_db
#endif
#endif
!2
#ifdef PHASE_CHANGE
                  feq=feq+p1*src(i,j,k)*(invrho_b-invrho_r)
#endif
				  F_discr=(-2.0_db*(forcex - 2.0_db*forcex*u(i,j,k) + forcey*v(i,j,k) &
				   + forcez*w(i,j,k)))/(9.0_db*rhophi_loc)
                  f(i-1,j,k,2)=feq + (1.0_db-omega_loc)*fneq1*p1 + 0.5_db*(F_discr)


!!!!!!!!!!!!!!!!!!!!!!!!!!3
#ifdef SECOND_ORDER
			      feq=(2.0_db*rho(i,j,k) - 3.0_db*(u(i,j,k)**2.0_db &
			       - 2.0_db*v(i,j,k)*(1.0_db + v(i,j,k)) + w(i,j,k)**2.0_db))/27.0_db
!3
				  fneq1=(-3.0_db*(pxx(i,j,k) - 2.0_db*pyy(i,j,k) + pzz(i,j,k)))/2.0_db
#endif
#ifdef FORTH_ORDER
!3
			      feq=(4.0_db*rho(i,j,k) + 12.0_db*v(i,j,k)*(1.0_db &
			       + v(i,j,k)) - 6.0_db*u(i,j,k)**2.0_db*(1.0_db &
			       + 3.0_db*v(i,j,k)*(1.0_db + v(i,j,k))) + 3.0_db*( &
			       -2.0_db + 3.0_db*u(i,j,k)**2.0_db - 6.0_db*v(i,j,k)*(1.0_db &
			       + v(i,j,k)))*w(i,j,k)**2.0_db)/54.0_db
!3
				  fneq1=(3.0_db*(-2.0_db*pxx(i,j,k) + 4.0_db*pyy(i,j,k) &
				   + 3.0_db*(-2.0_db*pyy(i,j,k) + pzz(i,j,k))*u(i,j,k)**2.0_db &
				   - 2.0_db*(pzz(i,j,k) + 6.0_db*pxy(i,j,k)*u(i,j,k)) - 6.0_db*(pxx(i,j,k) &
				   + pzz(i,j,k) + 6.0_db*pxy(i,j,k)*u(i,j,k))*v(i,j,k) - 6.0_db*(2.0_db*pxx(i,j,k) &
				   + pzz(i,j,k))*v(i,j,k)**2.0_db + 6.0_db*(3.0_db*pxz(i,j,k)*u(i,j,k) &
				   - 2.0_db*pyz(i,j,k)*(1.0_db + 3.0_db*v(i,j,k)))*w(i,j,k) &
				   + 6.0_db*(pxx(i,j,k) - 2.0_db*pyy(i,j,k))*w(i,j,k)**2.0_db))/4.0_db
#endif
#ifdef SIXTH_ORDER
!3
			      feq=(4.0_db*rho(i,j,k) + 12.0_db*v(i,j,k)*(1.0_db &
			       + v(i,j,k)) - 6.0_db*u(i,j,k)**2.0_db*(1.0_db &
			       + 3.0_db*v(i,j,k)*(1.0_db + v(i,j,k))) + 3.0_db*(-2.0_db &
			       + 3.0_db*u(i,j,k)**2.0_db)*(1.0_db + 3.0_db*v(i,j,k)*(1.0_db &
			       + v(i,j,k)))*w(i,j,k)**2.0_db)/54.0_db
!3
#ifdef GHOSTONE   
				  fneq1=(-3.0_db*(pxx(i,j,k) - 2.0_db*pyy(i,j,k) + pzz(i,j,k)))/2.0_db
#else 
				  fneq1=(12.0_db*pyy(i,j,k) + 9.0_db*(-2.0_db*pyy(i,j,k) &
				   + pzz(i,j,k))*u(i,j,k)**2.0_db - 6.0_db*pxx(i,j,k)*(1.0_db &
				   + 3.0_db*v(i,j,k) + 6.0_db*v(i,j,k)**2.0_db) + 9.0_db*v(i,j,k)*( &
				   -12.0_db*pxy(i,j,k)*u(i,j,k) + pzz(i,j,k)*(-2.0_db &
				   + 3.0_db*u(i,j,k)**2.0_db)*(1.0_db + v(i,j,k))) &
				   + 27.0_db*(2.0_db*pxz(i,j,k)*u(i,j,k) - 4.0_db*pyz(i,j,k)*v(i,j,k) &
				   + 6.0_db*pxz(i,j,k)*u(i,j,k)*v(i,j,k)*(1.0_db + v(i,j,k)) &
				   + pyz(i,j,k)*u(i,j,k)**2.0_db*(3.0_db + 6.0_db*v(i,j,k)))*w(i,j,k) &
				   + 9.0_db*(pyy(i,j,k)*(-4.0_db + 9.0_db*u(i,j,k)**2.0_db) &
				   + 6.0_db*pxy(i,j,k)*u(i,j,k)*(3.0_db + 7.0_db*v(i,j,k)) &
				   + pxx(i,j,k)*(2.0_db + 3.0_db*v(i,j,k)*(3.0_db &
				   + 4.0_db*v(i,j,k))))*w(i,j,k)**2.0_db - 6.0_db*(pzz(i,j,k) &
				   + 6.0_db*pxy(i,j,k)*u(i,j,k) + 6.0_db*pyz(i,j,k)*w(i,j,k)))/4.0_db
#endif
#endif
!3
#ifdef PHASE_CHANGE
                  feq=feq+p1*src(i,j,k)*(invrho_b-invrho_r)
#endif
				  F_discr=(2.0_db*(forcey - forcex*u(i,j,k) + 2.0_db*forcey*v(i,j,k) &
				   - forcez*w(i,j,k)))/(9.0_db*rhophi_loc)
                  f(i,j+1,k,3)=feq + (1.0_db-omega_loc)*fneq1*p1 + 0.5_db*(F_discr)


!!!!!!!!!!!!!!!!!!!!!!!!!!4
#ifdef SECOND_ORDER
			      feq=(2.0_db*rho(i,j,k) - 3.0_db*(u(i,j,k)**2.0_db &
			      - 2.0_db*(-1.0_db + v(i,j,k))*v(i,j,k) + w(i,j,k)**2.0_db))/27.0_db
!4
				  fneq1=(-3.0_db*(pxx(i,j,k) - 2.0_db*pyy(i,j,k) + pzz(i,j,k)))/2.0_db
#endif
#ifdef FORTH_ORDER
!4
			      feq=(4.0_db*rho(i,j,k) + 12.0_db*(-1.0_db &
			       + v(i,j,k))*v(i,j,k) - 6.0_db*u(i,j,k)**2.0_db*(1.0_db &
			       + 3.0_db*(-1.0_db + v(i,j,k))*v(i,j,k)) + 3.0_db*( &
			       -2.0_db + 3.0_db*u(i,j,k)**2.0_db - 6.0_db*(-1.0_db &
			       + v(i,j,k))*v(i,j,k))*w(i,j,k)**2.0_db)/54.0_db
!4
				  fneq1=(3.0_db*(4.0_db*pyy(i,j,k) - 2.0_db*pzz(i,j,k) &
				   + 3.0_db*(-2.0_db*pyy(i,j,k) + pzz(i,j,k))*u(i,j,k)**2.0_db &
				   + 12.0_db*(pxy(i,j,k)*u(i,j,k) + pyz(i,j,k)*w(i,j,k)) &
				   + 2.0_db*pxx(i,j,k)*(-1.0_db + 3.0_db*v(i,j,k) &
				   - 6.0_db*v(i,j,k)**2.0_db + 3.0_db*w(i,j,k)**2.0_db) &
				   - 6.0_db*(6.0_db*pxy(i,j,k)*u(i,j,k)*v(i,j,k) + pzz(i,j,k)*( &
				   -1.0_db + v(i,j,k))*v(i,j,k) + w(i,j,k)*(-3.0_db*pxz(i,j,k)*u(i,j,k) &
				   + 6.0_db*pyz(i,j,k)*v(i,j,k) + 2.0_db*pyy(i,j,k)*w(i,j,k)))))/4.0_db
#endif
#ifdef SIXTH_ORDER
!4
			      feq=(4.0_db*rho(i,j,k) + 12.0_db*(-1.0_db &
			       + v(i,j,k))*v(i,j,k) - 6.0_db*u(i,j,k)**2.0_db*(1.0_db &
			       + 3.0_db*(-1.0_db + v(i,j,k))*v(i,j,k)) + 3.0_db*(-2.0_db &
			       + 3.0_db*u(i,j,k)**2.0_db)*(1.0_db + 3.0_db*(-1.0_db &
			       + v(i,j,k))*v(i,j,k))*w(i,j,k)**2.0_db)/54.0_db
!4
#ifdef GHOSTONE   
				  fneq1=(-3.0_db*(pxx(i,j,k) - 2.0_db*pyy(i,j,k) + pzz(i,j,k)))/2.0_db
#else 
				  fneq1=(3.0_db*(4.0_db*pyy(i,j,k) - 2.0_db*pzz(i,j,k) &
				  + 3.0_db*(-2.0_db*pyy(i,j,k) + pzz(i,j,k))*u(i,j,k)**2.0_db &
				  + 3.0_db*(-12.0_db*pxy(i,j,k)*u(i,j,k) + pzz(i,j,k)*(-2.0_db &
				  + 3.0_db*u(i,j,k)**2.0_db)*(-1.0_db + v(i,j,k)))*v(i,j,k) &
				  - 2.0_db*pxx(i,j,k)*(1.0_db - 3.0_db*v(i,j,k) + 6.0_db*v(i,j,k)**2.0_db) &
				  + 9.0_db*(2.0_db*pxz(i,j,k)*u(i,j,k) - 4.0_db*pyz(i,j,k)*v(i,j,k) &
				  + 6.0_db*pxz(i,j,k)*u(i,j,k)*(-1.0_db + v(i,j,k))*v(i,j,k) &
				  + pyz(i,j,k)*u(i,j,k)**2.0_db*(-3.0_db + 6.0_db*v(i,j,k)))*w(i,j,k) &
				  + 3.0_db*pxx(i,j,k)*(2.0_db + 3.0_db*v(i,j,k)*(-3.0_db &
				  + 4.0_db*v(i,j,k)))*w(i,j,k)**2.0_db + 3.0_db*(pyy(i,j,k)*( &
				  -4.0_db + 9.0_db*u(i,j,k)**2.0_db) + 6.0_db*pxy(i,j,k)*u(i,j,k)*( &
				  -3.0_db + 7.0_db*v(i,j,k)))*w(i,j,k)**2.0_db + 12.0_db*(pxy(i,j,k)*u(i,j,k) &
				  + pyz(i,j,k)*w(i,j,k))))/4.0_db
#endif
#endif
!4
#ifdef PHASE_CHANGE
                  feq=feq+p1*src(i,j,k)*(invrho_b-invrho_r)
#endif
				  F_discr=(-2.0_db*(forcey + forcex*u(i,j,k) - 2.0_db*forcey*v(i,j,k) &
				   + forcez*w(i,j,k)))/(9.0_db*rhophi_loc)
                  f(i,j-1,k,4)=feq + (1.0_db-omega_loc)*fneq1*p1 + 0.5_db*(F_discr)


!!!!!!!!!!!!!!!!!!!!!!!!!!5
#ifdef SECOND_ORDER
			      feq=(2.0_db*rho(i,j,k) - 3.0_db*(u(i,j,k)**2.0_db &
			       + v(i,j,k)**2.0_db - 2.0_db*w(i,j,k)*(1.0_db + w(i,j,k))))/27.0_db
!5
				  fneq1=(-3.0_db*(pxx(i,j,k) + pyy(i,j,k) - 2.0_db*pzz(i,j,k)))/2.0_db
#endif
#ifdef FORTH_ORDER
!5
			      feq=(4.0_db*rho(i,j,k) + 3.0_db*(4.0_db*w(i,j,k)*(1.0_db &
			       + w(i,j,k)) + u(i,j,k)**2.0_db*(-2.0_db + 3.0_db*v(i,j,k)**2.0_db &
			       - 6.0_db*w(i,j,k)*(1.0_db + w(i,j,k))) - 2.0_db*v(i,j,k)**2.0_db*(1.0_db &
			       + 3.0_db*w(i,j,k)*(1.0_db + w(i,j,k)))))/54.0_db
!5
				  fneq1=(3.0_db*(-2.0_db*pxx(i,j,k) + pyy(i,j,k)*(-2.0_db &
				   + 3.0_db*u(i,j,k)**2.0_db) + 6.0_db*pxx(i,j,k)*v(i,j,k)**2.0_db &
				   - 6.0_db*pzz(i,j,k)*v(i,j,k)**2.0_db - 6.0_db*u(i,j,k)*(pzz(i,j,k)*u(i,j,k) &
				   - 3.0_db*pxy(i,j,k)*v(i,j,k)) + 4.0_db*(pzz(i,j,k) - 3.0_db*(pxz(i,j,k)*u(i,j,k) &
				   + pyz(i,j,k)*v(i,j,k))) - 6.0_db*(pxx(i,j,k) + pyy(i,j,k) &
				   + 6.0_db*pxz(i,j,k)*u(i,j,k) + 6.0_db*pyz(i,j,k)*v(i,j,k))*w(i,j,k) &
				   - 12.0_db*(pxx(i,j,k) + pyy(i,j,k))*w(i,j,k)**2.0_db))/4.0_db
#endif
#ifdef SIXTH_ORDER
!5
			      feq=(4.0_db*rho(i,j,k) + 3.0_db*(4.0_db*w(i,j,k)*(1.0_db &
			       + w(i,j,k)) - 2.0_db*v(i,j,k)**2.0_db*(1.0_db + 3.0_db*w(i,j,k)*(1.0_db &
			       + w(i,j,k))) + u(i,j,k)**2.0_db*(-2.0_db + 3.0_db*v(i,j,k)**2.0_db)*(1.0_db &
			       + 3.0_db*w(i,j,k)*(1.0_db + w(i,j,k)))))/54.0_db
!5
#ifdef GHOSTONE   
				  fneq1=(-3.0_db*(pxx(i,j,k) + pyy(i,j,k) - 2.0_db*pzz(i,j,k)))/2.0_db
#else 
				  fneq1=(3.0_db*(-6.0_db*pzz(i,j,k)*u(i,j,k)**2.0_db + 3.0_db*pzz(i,j,k)*( &
				   -2.0_db + 3.0_db*u(i,j,k)**2.0_db)*v(i,j,k)**2.0_db + 4.0_db*(pzz(i,j,k) &
				   - 3.0_db*(pxz(i,j,k)*u(i,j,k) + pyz(i,j,k)*v(i,j,k))) + 3.0_db*pyy(i,j,k)*(u(i,j,k) &
				   + 3.0_db*u(i,j,k)*w(i,j,k))**2.0_db - 2.0_db*pyy(i,j,k)*(1.0_db + 3.0_db*w(i,j,k) &
				   + 6.0_db*w(i,j,k)**2.0_db) + pxx(i,j,k)*(-2.0_db*(1.0_db + 3.0_db*w(i,j,k) &
				   + 6.0_db*w(i,j,k)**2.0_db) + 3.0_db*v(i,j,k)**2.0_db*(2.0_db &
				   + 3.0_db*w(i,j,k)*(3.0_db + 4.0_db*w(i,j,k)))) &
				   + 18.0_db*(pyz(i,j,k)*v(i,j,k)*(u(i,j,k)**2.0_db &
				   - 2.0_db*w(i,j,k) + 3.0_db*u(i,j,k)**2.0_db*w(i,j,k)) &
				   + pxz(i,j,k)*u(i,j,k)*(v(i,j,k)**2.0_db - 2.0_db*w(i,j,k) &
				   + 3.0_db*v(i,j,k)**2.0_db*w(i,j,k)) + pxy(i,j,k)*u(i,j,k)*v(i,j,k)*(1.0_db &
				   + w(i,j,k)*(5.0_db + 7.0_db*w(i,j,k))))))/4.0_db
#endif
#endif
!5
#ifdef PHASE_CHANGE
                  feq=feq+p1*src(i,j,k)*(invrho_b-invrho_r)
#endif
				  F_discr=(2.0_db*(forcez - forcex*u(i,j,k) - forcey*v(i,j,k) &
				   + 2.0_db*forcez*w(i,j,k)))/(9.0_db*rhophi_loc)
                  f(i,j,k+1,5)= feq + (1.0_db-omega_loc)*fneq1*p1 + 0.5_db*(F_discr)


!!!!!!!!!!!!!!!!!!!!!!!!!!6
#ifdef SECOND_ORDER
			      feq=(2.0_db*rho(i,j,k) - 3.0_db*(u(i,j,k)**2.0_db &
			       + v(i,j,k)**2.0_db - 2.0_db*(-1.0_db + w(i,j,k))*w(i,j,k)))/27.0_db
!6
				  fneq1=(-3.0_db*(pxx(i,j,k) + pyy(i,j,k) - 2.0_db*pzz(i,j,k)))/2.0_db
#endif
#ifdef FORTH_ORDER
!6
			      feq=(4.0_db*rho(i,j,k) + 3.0_db*(4.0_db*(-1.0_db &
			       + w(i,j,k))*w(i,j,k) + v(i,j,k)**2.0_db*(-2.0_db - 6.0_db*(-1.0_db &
			       + w(i,j,k))*w(i,j,k)) + u(i,j,k)**2.0_db*(-2.0_db + 3.0_db*v(i,j,k)**2.0_db &
			       - 6.0_db*(-1.0_db + w(i,j,k))*w(i,j,k))))/54.0_db
!6
				  fneq1=(3.0_db*(-2.0_db*pxx(i,j,k) + pyy(i,j,k)*(-2.0_db &
				   + 3.0_db*u(i,j,k)**2.0_db) + 6.0_db*pxx(i,j,k)*v(i,j,k)**2.0_db &
				   - 6.0_db*pzz(i,j,k)*v(i,j,k)**2.0_db - 6.0_db*u(i,j,k)*(pzz(i,j,k)*u(i,j,k) &
				   - 3.0_db*pxy(i,j,k)*v(i,j,k)) + 4.0_db*(pzz(i,j,k) &
				   + 3.0_db*pxz(i,j,k)*u(i,j,k) + 3.0_db*pyz(i,j,k)*v(i,j,k)) &
				   + 6.0_db*(pxx(i,j,k) + pyy(i,j,k) - 6.0_db*(pxz(i,j,k)*u(i,j,k) &
				   + pyz(i,j,k)*v(i,j,k)))*w(i,j,k) - 12.0_db*(pxx(i,j,k) &
				   + pyy(i,j,k))*w(i,j,k)**2.0_db))/4.0_db
#endif
#ifdef SIXTH_ORDER
!6
			      feq=(4.0_db*rho(i,j,k) + 3.0_db*(4.0_db*(-1.0_db &
			       + w(i,j,k))*w(i,j,k) + v(i,j,k)**2.0_db*(-2.0_db - 6.0_db*(-1.0_db &
			       + w(i,j,k))*w(i,j,k)) + u(i,j,k)**2.0_db*(-2.0_db &
			       + 3.0_db*v(i,j,k)**2.0_db)*(1.0_db + 3.0_db*(-1.0_db &
			       + w(i,j,k))*w(i,j,k))))/54.0_db
!6
#ifdef GHOSTONE   
				  fneq1=(-3.0_db*(pxx(i,j,k) + pyy(i,j,k) - 2.0_db*pzz(i,j,k)))/2.0_db
#else 
				  fneq1=(3.0_db*(4.0_db*(pzz(i,j,k) + 3.0_db*pxz(i,j,k)*u(i,j,k) &
				   + 3.0_db*pyz(i,j,k)*v(i,j,k)) + pyy(i,j,k)*(-2.0_db &
				   + 3.0_db*u(i,j,k)**2.0_db*(1.0_db - 3.0_db*w(i,j,k))**2.0_db &
				   + 6.0_db*(1.0_db - 2.0_db*w(i,j,k))*w(i,j,k)) + pxx(i,j,k)*( &
				   -2.0_db + 6.0_db*(1.0_db - 2.0_db*w(i,j,k))*w(i,j,k) &
				   + 3.0_db*v(i,j,k)**2.0_db*(2.0_db + 3.0_db*w(i,j,k)*(-3.0_db &
				   + 4.0_db*w(i,j,k)))) + 3.0_db*(-2.0_db*pzz(i,j,k)*u(i,j,k)**2.0_db &
				   + pzz(i,j,k)*(-2.0_db + 3.0_db*u(i,j,k)**2.0_db)*v(i,j,k)**2.0_db &
				   - 6.0_db*u(i,j,k)*v(i,j,k)*(pyz(i,j,k)*u(i,j,k) + pxz(i,j,k)*v(i,j,k)) &
				   + 6.0_db*(pyz(i,j,k)*(-2.0_db + 3.0_db*u(i,j,k)**2.0_db)*v(i,j,k) &
				   + pxz(i,j,k)*u(i,j,k)*(-2.0_db + 3.0_db*v(i,j,k)**2.0_db))*w(i,j,k) &
				   + 6.0_db*pxy(i,j,k)*u(i,j,k)*v(i,j,k)*(1.0_db + w(i,j,k)*(-5.0_db &
				   + 7.0_db*w(i,j,k))))))/4.0_db
#endif
#endif
!6
#ifdef PHASE_CHANGE
                  feq=feq+p1*src(i,j,k)*(invrho_b-invrho_r)
#endif
				  F_discr=(-2.0_db*(forcez + forcex*u(i,j,k) + forcey*v(i,j,k) &
				   - 2.0_db*forcez*w(i,j,k)))/(9.0_db*rhophi_loc)
                  f(i,j,k-1,6)=feq + (1.0_db-omega_loc)*fneq1*p1 + 0.5_db*(F_discr)


!!!!!!!!!!!!!!!!!!!!!!!!!!7
#ifdef SECOND_ORDER
			      feq=(2.0_db*rho(i,j,k) + 6.0_db*(u(i,j,k) &
			       + u(i,j,k)**2.0_db + v(i,j,k) + 3.0_db*u(i,j,k)*v(i,j,k) &
			       + v(i,j,k)**2.0_db) - 3.0_db*w(i,j,k)**2.0_db)/108.0_db
!7
				  fneq1=3.0_db*(pxx(i,j,k) + 3.0_db*pxy(i,j,k) + pyy(i,j,k)) &
				   - (3.0_db*pzz(i,j,k))/2.0_db
#endif
#ifdef FORTH_ORDER
!7
			      feq=(2.0_db*rho(i,j,k) + 6.0_db*(u(i,j,k) &
			       + v(i,j,k) + v(i,j,k)**2.0_db + 3.0_db*u(i,j,k)*v(i,j,k)*(1.0_db &
			       + v(i,j,k)) + u(i,j,k)**2.0_db*(1.0_db + 3.0_db*v(i,j,k)*(1.0_db &
			       + v(i,j,k)))) - 3.0_db*(1.0_db + 3.0_db*v(i,j,k) &
			       + 3.0_db*(u(i,j,k) + u(i,j,k)**2.0_db + 3.0_db*u(i,j,k)*v(i,j,k) &
			       + v(i,j,k)**2.0_db))*w(i,j,k)**2.0_db)/108.0_db
!7
				  fneq1=(3.0_db*(6.0_db*pxy(i,j,k)*(1.0_db + 2.0_db*u(i,j,k)) &
				   + (2.0_db*pyy(i,j,k) - pzz(i,j,k))*(1.0_db + 3.0_db*u(i,j,k)*(1.0_db + u(i,j,k))) &
				   + 12.0_db*pxy(i,j,k)*(1.0_db + 3.0_db*u(i,j,k))*v(i,j,k) &
				   - 3.0_db*pzz(i,j,k)*v(i,j,k)*(1.0_db + 3.0_db*u(i,j,k) + v(i,j,k)) &
				   - 6.0_db*(pxz(i,j,k) + pyz(i,j,k))*w(i,j,k) &
				   - 9.0_db*(2.0_db*pxz(i,j,k)*u(i,j,k) + 3.0_db*pyz(i,j,k)*u(i,j,k) &
				   + 3.0_db*pxz(i,j,k)*v(i,j,k) + 2.0_db*pyz(i,j,k)*v(i,j,k))*w(i,j,k) &
				   - 6.0_db*(3.0_db*pxy(i,j,k) + pyy(i,j,k))*w(i,j,k)**2.0_db &
				   + 2.0_db*pxx(i,j,k)*(1.0_db + 3.0_db*v(i,j,k) &
				   + 6.0_db*v(i,j,k)**2.0_db - 3.0_db*w(i,j,k)**2.0_db)))/2.0_db
#endif
#ifdef SIXTH_ORDER
!7
			      feq=(2.0_db*rho(i,j,k) + 6.0_db*(u(i,j,k) &
			       + v(i,j,k) + v(i,j,k)**2.0_db + 3.0_db*u(i,j,k)*v(i,j,k)*(1.0_db &
			       + v(i,j,k)) + u(i,j,k)**2.0_db*(1.0_db + 3.0_db*v(i,j,k)*(1.0_db &
			       + v(i,j,k)))) - 3.0_db*(1.0_db + 3.0_db*u(i,j,k)*(1.0_db &
			       + u(i,j,k)))*(1.0_db + 3.0_db*v(i,j,k)*(1.0_db &
			       + v(i,j,k)))*w(i,j,k)**2.0_db)/108.0_db
!7
#ifdef GHOSTONE   
				  fneq1=3.0_db*(pxx(i,j,k) + 3.0_db*pxy(i,j,k) + pyy(i,j,k)) &
				   - (3.0_db*pzz(i,j,k))/2.0_db
#else 
				  fneq1=(3.0_db*(2.0_db*pxx(i,j,k)*(1.0_db + 3.0_db*v(i,j,k) &
				   + 6.0_db*v(i,j,k)**2.0_db) + (1.0_db + 3.0_db*u(i,j,k)*(1.0_db &
				   + u(i,j,k)))*(2.0_db*pyy(i,j,k) - pzz(i,j,k)*(1.0_db &
				   + 3.0_db*v(i,j,k)*(1.0_db + v(i,j,k)))) &
				   - 3.0_db*(pyz(i,j,k)*(2.0_db + 6.0_db*v(i,j,k) &
				   + 9.0_db*u(i,j,k)*(1.0_db + u(i,j,k))*(1.0_db &
				   + 2.0_db*v(i,j,k))) + pxz(i,j,k)*(2.0_db + 9.0_db*v(i,j,k)*(1.0_db &
				   + v(i,j,k)) + 6.0_db*u(i,j,k)*(1.0_db + 3.0_db*v(i,j,k)*(1.0_db &
				   + v(i,j,k)))))*w(i,j,k) - 3.0_db*pyy(i,j,k)*(2.0_db &
				   + 9.0_db*u(i,j,k)*(1.0_db + u(i,j,k)))*w(i,j,k)**2.0_db &
				   - 3.0_db*pxx(i,j,k)*(2.0_db + 3.0_db*v(i,j,k)*(3.0_db &
				   + 4.0_db*v(i,j,k)))*w(i,j,k)**2.0_db + 6.0_db*pxy(i,j,k)*(1.0_db &
				   + 2.0_db*u(i,j,k) + 2.0_db*v(i,j,k) + 6.0_db*u(i,j,k)*v(i,j,k) &
				   - 3.0_db*(1.0_db + 3.0_db*v(i,j,k) + u(i,j,k)*(3.0_db &
				   + 7.0_db*v(i,j,k)))*w(i,j,k)**2.0_db)))/2.0_db
#endif
#endif
!7
#ifdef PHASE_CHANGE
                  feq=feq+p2*src(i,j,k)*(invrho_b-invrho_r)
#endif
				  F_discr=(forcex + forcey + 2.0_db*forcex*u(i,j,k) &
				   + 3.0_db*forcey*u(i,j,k) + 3.0_db*forcex*v(i,j,k) &
				   + 2.0_db*forcey*v(i,j,k) - forcez*w(i,j,k))/(18.0_db*rhophi_loc)
                  f(i+1,j+1,k,7)=feq + (1.0_db-omega_loc)*fneq1*p2 + 0.5_db*(F_discr)


!!!!!!!!!!!!!!!!!!!!!!!!!!8
#ifdef SECOND_ORDER
			      feq=(2.0_db*rho(i,j,k) + 6.0_db*(u(i,j,k)**2.0_db &
			       + (-1.0_db + v(i,j,k))*v(i,j,k) + u(i,j,k)*(-1.0_db &
			       + 3.0_db*v(i,j,k))) - 3.0_db*w(i,j,k)**2.0_db)/108.0_db
!8
				  fneq1=3.0_db*(pxx(i,j,k) + 3.0_db*pxy(i,j,k) + pyy(i,j,k)) &
				   - (3.0_db*pzz(i,j,k))/2.0_db
#endif
#ifdef FORTH_ORDER
!8
			      feq=(2.0_db*rho(i,j,k) + 6.0_db*((-1.0_db &
			       + v(i,j,k))*v(i,j,k) + u(i,j,k)*(-1.0_db - 3.0_db*(-1.0_db &
			       + v(i,j,k))*v(i,j,k)) + u(i,j,k)**2.0_db*(1.0_db + 3.0_db*(-1.0_db &
			       + v(i,j,k))*v(i,j,k))) - 3.0_db*(1.0_db &
			       + 3.0_db*u(i,j,k)**2.0_db + 3.0_db*(-1.0_db + v(i,j,k))*v(i,j,k) &
			       + u(i,j,k)*(-3.0_db + 9.0_db*v(i,j,k)))*w(i,j,k)**2.0_db)/108.0_db
!8
				  fneq1=(3.0_db*(pxy(i,j,k)*(6.0_db - 12.0_db*u(i,j,k)) &
				   + (2.0_db*pyy(i,j,k) - pzz(i,j,k))*(1.0_db + 3.0_db*( &
				   -1.0_db + u(i,j,k))*u(i,j,k)) + 3.0_db*v(i,j,k)*(pzz(i,j,k) &
				   + 4.0_db*pxy(i,j,k)*(-1.0_db + 3.0_db*u(i,j,k)) &
				   - pzz(i,j,k)*(3.0_db*u(i,j,k) + v(i,j,k))) + 6.0_db*(pxz(i,j,k) &
				   + pyz(i,j,k))*w(i,j,k) - 9.0_db*(2.0_db*pxz(i,j,k)*u(i,j,k) &
				   + 3.0_db*pyz(i,j,k)*u(i,j,k) + 3.0_db*pxz(i,j,k)*v(i,j,k) &
				   + 2.0_db*pyz(i,j,k)*v(i,j,k))*w(i,j,k) - 6.0_db*(3.0_db*pxy(i,j,k) &
				   + pyy(i,j,k))*w(i,j,k)**2.0_db + 2.0_db*pxx(i,j,k)*(1.0_db &
				   - 3.0_db*v(i,j,k) + 6.0_db*v(i,j,k)**2.0_db &
				   - 3.0_db*w(i,j,k)**2.0_db)))/2.0_db
#endif
#ifdef SIXTH_ORDER
!8
			      feq=(2.0_db*rho(i,j,k) + 6.0_db*((-1.0_db &
			       + v(i,j,k))*v(i,j,k) + u(i,j,k)*(-1.0_db - 3.0_db*(-1.0_db &
			       + v(i,j,k))*v(i,j,k)) + u(i,j,k)**2.0_db*(1.0_db + 3.0_db*(-1.0_db &
			       + v(i,j,k))*v(i,j,k))) - 3.0_db*(1.0_db + 3.0_db*(-1.0_db &
			       + u(i,j,k))*u(i,j,k))*(1.0_db + 3.0_db*(-1.0_db &
			       + v(i,j,k))*v(i,j,k))*w(i,j,k)**2.0_db)/108.0_db
!8
#ifdef GHOSTONE   
				  fneq1=3.0_db*(pxx(i,j,k) + 3.0_db*pxy(i,j,k) + pyy(i,j,k)) &
				   - (3.0_db*pzz(i,j,k))/2.0_db
#else 
				  fneq1=(3.0_db*(2.0_db*pxx(i,j,k)*(1.0_db - 3.0_db*v(i,j,k) &
				   + 6.0_db*v(i,j,k)**2.0_db) + (1.0_db + 3.0_db*(-1.0_db &
				   + u(i,j,k))*u(i,j,k))*(2.0_db*pyy(i,j,k) + pzz(i,j,k)*( &
				   -1.0_db - 3.0_db*(-1.0_db + v(i,j,k))*v(i,j,k))) &
				   - 3.0_db*(pyz(i,j,k)*(-2.0_db + 6.0_db*v(i,j,k) &
				   + 9.0_db*(-1.0_db + u(i,j,k))*u(i,j,k)*(-1.0_db &
				   + 2.0_db*v(i,j,k))) + pxz(i,j,k)*(-2.0_db - 9.0_db*( &
				   -1.0_db + v(i,j,k))*v(i,j,k) + 6.0_db*u(i,j,k)*(1.0_db &
				   + 3.0_db*(-1.0_db + v(i,j,k))*v(i,j,k))))*w(i,j,k) &
				   - 3.0_db*pyy(i,j,k)*(2.0_db + 9.0_db*(-1.0_db &
				   + u(i,j,k))*u(i,j,k))*w(i,j,k)**2.0_db &
				   - 3.0_db*pxx(i,j,k)*(2.0_db + 3.0_db*v(i,j,k)*( &
				   -3.0_db + 4.0_db*v(i,j,k)))*w(i,j,k)**2.0_db &
				   + 6.0_db*pxy(i,j,k)*(1.0_db - 2.0_db*u(i,j,k) &
				   - 2.0_db*v(i,j,k) + 6.0_db*u(i,j,k)*v(i,j,k) + 3.0_db*( &
				   -1.0_db + u(i,j,k)*(3.0_db - 7.0_db*v(i,j,k)) &
				   + 3.0_db*v(i,j,k))*w(i,j,k)**2.0_db)))/2.0_db
#endif
#endif
!8
#ifdef PHASE_CHANGE
                  feq=feq+p2*src(i,j,k)*(invrho_b-invrho_r)
#endif
				  F_discr=(forcex + forcey - 2.0_db*forcex*u(i,j,k) &
				   - 3.0_db*forcey*u(i,j,k) - 3.0_db*forcex*v(i,j,k) &
				   - 2.0_db*forcey*v(i,j,k) + forcez*w(i,j,k))/(-18.0_db*rhophi_loc)
                  f(i-1,j-1,k,8)=feq + (1.0_db-omega_loc)*fneq1*p2 + 0.5_db*(F_discr)


!!!!!!!!!!!!!!!!!!!!!!!!!!9
#ifdef SECOND_ORDER
			      feq=(2.0_db*rho(i,j,k) + 6.0_db*(u(i,j,k) &
			       + u(i,j,k)**2.0_db - 3.0_db*u(i,j,k)*v(i,j,k) + (-1.0_db &
			       + v(i,j,k))*v(i,j,k)) - 3.0_db*w(i,j,k)**2.0_db)/108.0_db
!9
				  fneq1=3.0_db*(pxx(i,j,k) - 3.0_db*pxy(i,j,k) + pyy(i,j,k)) &
				   - (3.0_db*pzz(i,j,k))/2.0_db
#endif
#ifdef FORTH_ORDER
!9
			      feq=(2.0_db*rho(i,j,k) + 6.0_db*(u(i,j,k) &
			       + (-1.0_db + v(i,j,k))*v(i,j,k) + 3.0_db*u(i,j,k)*(-1.0_db &
			       + v(i,j,k))*v(i,j,k) + u(i,j,k)**2.0_db*(1.0_db + 3.0_db*(-1.0_db &
			       + v(i,j,k))*v(i,j,k))) - 3.0_db*(1.0_db - 3.0_db*v(i,j,k) &
			       + 3.0_db*(u(i,j,k) + u(i,j,k)**2.0_db - 3.0_db*u(i,j,k)*v(i,j,k) &
			       + v(i,j,k)**2.0_db))*w(i,j,k)**2.0_db)/108.0_db
!9
				  fneq1=(3.0_db*(-6.0_db*pxy(i,j,k)*(1.0_db + 2.0_db*u(i,j,k)) &
				   + (2.0_db*pyy(i,j,k) - pzz(i,j,k))*(1.0_db + 3.0_db*u(i,j,k)*(1.0_db &
				   + u(i,j,k))) + 3.0_db*v(i,j,k)*((4.0_db*pxy(i,j,k) &
				   + pzz(i,j,k))*(1.0_db + 3.0_db*u(i,j,k)) - pzz(i,j,k)*v(i,j,k)) &
				   + 6.0_db*(-pxz(i,j,k) + pyz(i,j,k))*w(i,j,k) + 9.0_db*( &
				   -2.0_db*pxz(i,j,k)*u(i,j,k) + 3.0_db*pyz(i,j,k)*u(i,j,k) &
				   + 3.0_db*pxz(i,j,k)*v(i,j,k) - 2.0_db*pyz(i,j,k)*v(i,j,k))*w(i,j,k) &
				   + 6.0_db*(3.0_db*pxy(i,j,k) - pyy(i,j,k))*w(i,j,k)**2.0_db &
				   + 2.0_db*pxx(i,j,k)*(1.0_db - 3.0_db*v(i,j,k) &
				   + 6.0_db*v(i,j,k)**2.0_db - 3.0_db*w(i,j,k)**2.0_db)))/2.0_db
#endif
#ifdef SIXTH_ORDER
!9
			      feq=(2.0_db*rho(i,j,k) + 6.0_db*(u(i,j,k) &
			       + (-1.0_db + v(i,j,k))*v(i,j,k) + 3.0_db*u(i,j,k)*(-1.0_db &
			       + v(i,j,k))*v(i,j,k) + u(i,j,k)**2.0_db*(1.0_db + 3.0_db*(-1.0_db &
			       + v(i,j,k))*v(i,j,k))) - 3.0_db*(1.0_db &
			       + 3.0_db*u(i,j,k)*(1.0_db + u(i,j,k)))*(1.0_db + 3.0_db*(-1.0_db &
			       + v(i,j,k))*v(i,j,k))*w(i,j,k)**2.0_db)/108.0_db
!9
#ifdef GHOSTONE   
				  fneq1=3.0_db*(pxx(i,j,k) - 3.0_db*pxy(i,j,k) + pyy(i,j,k)) &
				   - (3.0_db*pzz(i,j,k))/2.0_db
#else 
				  fneq1=(3.0_db*(2.0_db*pxx(i,j,k)*(1.0_db - 3.0_db*v(i,j,k) &
				   + 6.0_db*v(i,j,k)**2.0_db) + (1.0_db + 3.0_db*u(i,j,k)*(1.0_db &
				   + u(i,j,k)))*(2.0_db*pyy(i,j,k) + pzz(i,j,k)*(-1.0_db &
				   - 3.0_db*(-1.0_db + v(i,j,k))*v(i,j,k))) - 3.0_db*(pyz(i,j,k)*( &
				   -2.0_db + 6.0_db*v(i,j,k) + 9.0_db*u(i,j,k)*(1.0_db &
				   + u(i,j,k))*(-1.0_db + 2.0_db*v(i,j,k))) + pxz(i,j,k)*(2.0_db &
				   + 9.0_db*(-1.0_db + v(i,j,k))*v(i,j,k) + 6.0_db*u(i,j,k)*(1.0_db &
				   + 3.0_db*(-1.0_db + v(i,j,k))*v(i,j,k))))*w(i,j,k) &
				   - 3.0_db*pyy(i,j,k)*(2.0_db + 9.0_db*u(i,j,k)*(1.0_db &
				   + u(i,j,k)))*w(i,j,k)**2.0_db - 3.0_db*pxx(i,j,k)*(2.0_db &
				   + 3.0_db*v(i,j,k)*(-3.0_db + 4.0_db*v(i,j,k)))*w(i,j,k)**2.0_db &
				   + 6.0_db*pxy(i,j,k)*(-1.0_db - 2.0_db*u(i,j,k) + 2.0_db*v(i,j,k) &
				   + 6.0_db*u(i,j,k)*v(i,j,k) + 3.0_db*(1.0_db + u(i,j,k)*(3.0_db &
				   - 7.0_db*v(i,j,k)) - 3.0_db*v(i,j,k))*w(i,j,k)**2.0_db)))/2.0_db
#endif
#endif
!9
#ifdef PHASE_CHANGE
                  feq=feq+p2*src(i,j,k)*(invrho_b-invrho_r)
#endif
				  F_discr=(forcey + 3.0_db*forcey*u(i,j,k) - 2.0_db*forcey*v(i,j,k) &
				   + forcex*(-1.0_db - 2.0_db*u(i,j,k) + 3.0_db*v(i,j,k)) &
				   + forcez*w(i,j,k))/(-18.0_db*rhophi_loc)
                  f(i+1,j-1,k,9)=feq + (1.0_db-omega_loc)*fneq1*p2 + 0.5_db*(F_discr)


!!!!!!!!!!!!!!!!!!!!!!!!!!10
#ifdef SECOND_ORDER
			      feq=(2.0_db*rho(i,j,k) + 6.0_db*((-1.0_db &
			       + u(i,j,k))*u(i,j,k) + v(i,j,k) - 3.0_db*u(i,j,k)*v(i,j,k) &
			       + v(i,j,k)**2.0_db) - 3.0_db*w(i,j,k)**2.0_db)/108.0_db
!10
				  fneq1=3.0_db*(pxx(i,j,k) - 3.0_db*pxy(i,j,k) + pyy(i,j,k)) &
				   - (3.0_db*pzz(i,j,k))/2.0_db
#endif
#ifdef FORTH_ORDER
!10
			      feq=(2.0_db*rho(i,j,k) + 6.0_db*((-1.0_db &
			       + u(i,j,k))*u(i,j,k) + v(i,j,k) + 3.0_db*(-1.0_db &
			       + u(i,j,k))*u(i,j,k)*v(i,j,k) + (1.0_db + 3.0_db*(-1.0_db &
			       + u(i,j,k))*u(i,j,k))*v(i,j,k)**2.0_db) - 3.0_db*(1.0_db &
			       + 3.0_db*u(i,j,k)**2.0_db + 3.0_db*v(i,j,k)*(1.0_db + v(i,j,k)) &
			       - 3.0_db*u(i,j,k)*(1.0_db + 3.0_db*v(i,j,k)))*w(i,j,k)**2.0_db)/108.0_db
!10
				  fneq1=(3.0_db*(6.0_db*pxy(i,j,k)*(-1.0_db + 2.0_db*u(i,j,k)) &
				   + (2.0_db*pyy(i,j,k) - pzz(i,j,k))*(1.0_db + 3.0_db*( &
				   -1.0_db + u(i,j,k))*u(i,j,k)) + 3.0_db*v(i,j,k)*((4.0_db*pxy(i,j,k) &
				   + pzz(i,j,k))*(-1.0_db + 3.0_db*u(i,j,k)) - pzz(i,j,k)*v(i,j,k)) &
				   + 6.0_db*(pxz(i,j,k) - pyz(i,j,k))*w(i,j,k) + 9.0_db*( &
				   -2.0_db*pxz(i,j,k)*u(i,j,k) + 3.0_db*pyz(i,j,k)*u(i,j,k) &
				   + 3.0_db*pxz(i,j,k)*v(i,j,k) - 2.0_db*pyz(i,j,k)*v(i,j,k))*w(i,j,k) &
				   + 6.0_db*(3.0_db*pxy(i,j,k) - pyy(i,j,k))*w(i,j,k)**2.0_db &
				   + 2.0_db*pxx(i,j,k)*(1.0_db + 3.0_db*v(i,j,k) &
				   + 6.0_db*v(i,j,k)**2.0_db - 3.0_db*w(i,j,k)**2.0_db)))/2.0_db
#endif
#ifdef SIXTH_ORDER
!10
			      feq=(2.0_db*rho(i,j,k) + 6.0_db*((-1.0_db &
			       + u(i,j,k))*u(i,j,k) + v(i,j,k) + 3.0_db*(-1.0_db &
			       + u(i,j,k))*u(i,j,k)*v(i,j,k) + (1.0_db + 3.0_db*(-1.0_db &
			       + u(i,j,k))*u(i,j,k))*v(i,j,k)**2.0_db) - 3.0_db*(1.0_db &
			       + 3.0_db*(-1.0_db + u(i,j,k))*u(i,j,k))*(1.0_db + 3.0_db*v(i,j,k)*(1.0_db &
			       + v(i,j,k)))*w(i,j,k)**2.0_db)/108.0_db
!10
#ifdef GHOSTONE   
				  fneq1=3.0_db*(pxx(i,j,k) - 3.0_db*pxy(i,j,k) + pyy(i,j,k)) &
				   - (3.0_db*pzz(i,j,k))/2.0_db
#else 
				  fneq1=(3.0_db*(2.0_db*pxx(i,j,k)*(1.0_db + 3.0_db*v(i,j,k) &
				   + 6.0_db*v(i,j,k)**2.0_db) + (1.0_db + 3.0_db*(-1.0_db &
				   + u(i,j,k))*u(i,j,k))*(2.0_db*pyy(i,j,k) - pzz(i,j,k)*(1.0_db &
				   + 3.0_db*v(i,j,k)*(1.0_db + v(i,j,k)))) - 3.0_db*(pyz(i,j,k)*(2.0_db &
				   + 6.0_db*v(i,j,k) + 9.0_db*(-1.0_db + u(i,j,k))*u(i,j,k)*(1.0_db &
				   + 2.0_db*v(i,j,k))) + pxz(i,j,k)*(-2.0_db - 9.0_db*v(i,j,k) &
				   - 9.0_db*v(i,j,k)**2.0_db + 6.0_db*u(i,j,k)*(1.0_db &
				   + 3.0_db*v(i,j,k)*(1.0_db + v(i,j,k)))))*w(i,j,k) &
				   - 3.0_db*pyy(i,j,k)*(2.0_db + 9.0_db*(-1.0_db &
				   + u(i,j,k))*u(i,j,k))*w(i,j,k)**2.0_db &
				   - 3.0_db*pxx(i,j,k)*(2.0_db + 3.0_db*v(i,j,k)*(3.0_db &
				   + 4.0_db*v(i,j,k)))*w(i,j,k)**2.0_db + 6.0_db*pxy(i,j,k)*(&
				   -1.0_db + 2.0_db*u(i,j,k) - 2.0_db*v(i,j,k) + 6.0_db*u(i,j,k)*v(i,j,k) &
				   + 3.0_db*(1.0_db + 3.0_db*v(i,j,k) - u(i,j,k)*(3.0_db &
				   + 7.0_db*v(i,j,k)))*w(i,j,k)**2.0_db)))/2.0_db
#endif
#endif
!10
#ifdef PHASE_CHANGE
                  feq=feq+p2*src(i,j,k)*(invrho_b-invrho_r)
#endif
				  F_discr=(forcex - forcey - 2.0_db*forcex*u(i,j,k) &
				   + 3.0_db*forcey*u(i,j,k) + 3.0_db*forcex*v(i,j,k) &
				   - 2.0_db*forcey*v(i,j,k) + forcez*w(i,j,k))/(-18.0_db*rhophi_loc)
                  f(i-1,j+1,k,10)=feq + (1.0_db-omega_loc)*fneq1*p2 + 0.5_db*(F_discr)


!!!!!!!!!!!!!!!!!!!!!!!!!!11
#ifdef SECOND_ORDER
			      feq=(2.0_db*rho(i,j,k) - 3.0_db*u(i,j,k)**2.0_db &
			       + 6.0_db*(v(i,j,k) + v(i,j,k)**2.0_db + w(i,j,k) &
			       + 3.0_db*v(i,j,k)*w(i,j,k) + w(i,j,k)**2.0_db))/108.0_db
!11
				  fneq1=(-3.0_db*pxx(i,j,k))/2.0_db + 3.0_db*(pyy(i,j,k) &
				   + 3.0_db*pyz(i,j,k) + pzz(i,j,k))
#endif
#ifdef FORTH_ORDER
!11
			      feq=(2.0_db*rho(i,j,k) - 3.0_db*u(i,j,k)**2.0_db*(1.0_db &
			       + 3.0_db*w(i,j,k) + 3.0_db*(v(i,j,k) + v(i,j,k)**2.0_db &
			       + 3.0_db*v(i,j,k)*w(i,j,k) + w(i,j,k)**2.0_db)) + 6.0_db &
			       *(v(i,j,k) + w(i,j,k) + w(i,j,k)**2.0_db + 3.0_db*v(i,j,k)*w(i,j,k)*(1.0_db &
			       + w(i,j,k)) + v(i,j,k)**2.0_db*(1.0_db + 3.0_db*w(i,j,k)*(1.0_db &
			       + w(i,j,k)))))/108.0_db
!11
				  fneq1=(3.0_db*(2.0_db*(pyy(i,j,k) + 3.0_db*pyz(i,j,k) + pzz(i,j,k)) &
				   - 6.0_db*(pxy(i,j,k) + pxz(i,j,k))*u(i,j,k) - 3.0_db*(pyy(i,j,k) &
				   + 3.0_db*pyz(i,j,k) + pzz(i,j,k))*u(i,j,k)**2.0_db &
				   - pxx(i,j,k)*(1.0_db + 6.0_db*v(i,j,k)**2.0_db + 3.0_db*w(i,j,k) &
				   + 6.0_db*w(i,j,k)**2.0_db + 3.0_db*v(i,j,k)*(1.0_db + 6.0_db*w(i,j,k))) &
				   + 6.0_db*(-3.0_db*(pxy(i,j,k) + pxz(i,j,k))*u(i,j,k)*v(i,j,k) &
				   + pzz(i,j,k)*v(i,j,k)*(1.0_db + v(i,j,k)) + (pyy(i,j,k) &
				   - 3.0_db*(2.0_db*pxy(i,j,k) + pxz(i,j,k))*u(i,j,k))*w(i,j,k) &
				   + 2.0_db*pyy(i,j,k)*w(i,j,k)**2.0_db + 2.0_db*pyz(i,j,k)*(v(i,j,k) &
				   + w(i,j,k) + 3.0_db*v(i,j,k)*w(i,j,k)))))/2.0_db
#endif
#ifdef SIXTH_ORDER
!11
			      feq=(2.0_db*rho(i,j,k) - 3.0_db*u(i,j,k)**2.0_db*(1.0_db &
			       + 3.0_db*v(i,j,k)*(1.0_db + v(i,j,k)))*(1.0_db + 3.0_db*w(i,j,k)*(1.0_db &
			       + w(i,j,k))) + 6.0_db*(v(i,j,k) + w(i,j,k) + w(i,j,k)**2.0_db &
			       + 3.0_db*v(i,j,k)*w(i,j,k)*(1.0_db + w(i,j,k)) + v(i,j,k)**2.0_db*(1.0_db &
			       + 3.0_db*w(i,j,k)*(1.0_db + w(i,j,k)))))/108.0_db
!11
#ifdef GHOSTONE   
				  fneq1=(-3.0_db*pxx(i,j,k))/2.0_db + 3.0_db*(pyy(i,j,k) &
				   + 3.0_db*pyz(i,j,k) + pzz(i,j,k))
#else 
				  fneq1=(3.0_db*(2.0_db*(pyy(i,j,k) + 3.0_db*pyz(i,j,k) + pzz(i,j,k)) &
				   - 6.0_db*(pxy(i,j,k) + pxz(i,j,k))*u(i,j,k) &
				   - 3.0_db*(pyy(i,j,k) + 3.0_db*pyz(i,j,k) + pzz(i,j,k))*u(i,j,k)**2.0_db &
				   - 3.0_db*v(i,j,k)*(pzz(i,j,k)*(-2.0_db + 3.0_db*u(i,j,k)**2.0_db)*(1.0_db &
				   + v(i,j,k)) + 6.0_db*u(i,j,k)*(pxy(i,j,k) + pxz(i,j,k) + pxz(i,j,k)*v(i,j,k))) &
				   + 6.0_db*pyy(i,j,k)*w(i,j,k) - 18.0_db*u(i,j,k)*(pxz(i,j,k) &
				   + pyy(i,j,k)*u(i,j,k) + 3.0_db*pxz(i,j,k)*v(i,j,k)*(1.0_db + v(i,j,k)) &
				   + pxy(i,j,k)*(2.0_db + 5.0_db*v(i,j,k)))*w(i,j,k) &
				   - 3.0_db*(pyy(i,j,k)*(-4.0_db + 9.0_db*u(i,j,k)**2.0_db) &
				   + 6.0_db*pxy(i,j,k)*u(i,j,k)*(3.0_db + 7.0_db*v(i,j,k)))*w(i,j,k)**2.0_db &
				   - 3.0_db*pyz(i,j,k)*((-4.0_db + 9.0_db*u(i,j,k)**2.0_db)*w(i,j,k) &
				   + 2.0_db*(-2.0_db + 3.0_db*u(i,j,k)**2.0_db)*v(i,j,k)*(1.0_db &
				   + 3.0_db*w(i,j,k))) + pxx(i,j,k)*(-1.0_db - 3.0_db*w(i,j,k)*(1.0_db &
				   + 2.0_db*w(i,j,k)) - 3.0_db*v(i,j,k)*(1.0_db + 3.0_db*w(i,j,k))**2.0_db &
				   - 3.0_db*v(i,j,k)**2.0_db*(2.0_db + 3.0_db*w(i,j,k)*(3.0_db &
				   + 4.0_db*w(i,j,k))))))/2.0_db
#endif
#endif
!11
#ifdef PHASE_CHANGE
                  feq=feq+p2*src(i,j,k)*(invrho_b-invrho_r)
#endif
				  F_discr=(forcey + forcez - forcex*u(i,j,k) + 2.0_db*forcey*v(i,j,k) &
				   + 3.0_db*forcez*v(i,j,k) + 3.0_db*forcey*w(i,j,k) &
				   + 2.0_db*forcez*w(i,j,k))/(18.0_db*rhophi_loc)
                  f(i,j+1,k+1,11)=feq + (1.0_db-omega_loc)*fneq1*p2 + 0.5_db*(F_discr)


!!!!!!!!!!!!!!!!!!!!!!!!!!12
#ifdef SECOND_ORDER
			      feq=(2.0_db*rho(i,j,k) - 3.0_db*u(i,j,k)**2.0_db &
			       + 6.0_db*(v(i,j,k)**2.0_db + (-1.0_db &
			       + w(i,j,k))*w(i,j,k) + v(i,j,k)*(-1.0_db + 3.0_db*w(i,j,k))))/108.0_db
!12
				  fneq1=(-3.0_db*pxx(i,j,k))/2.0_db + 3.0_db*(pyy(i,j,k) &
				   + 3.0_db*pyz(i,j,k) + pzz(i,j,k))
#endif
#ifdef FORTH_ORDER
!12
			      feq=(2.0_db*rho(i,j,k) + 3.0_db*(2.0_db*( &
			       -1.0_db + w(i,j,k))*w(i,j,k) + v(i,j,k)*(-2.0_db &
			       - 6.0_db*(-1.0_db + w(i,j,k))*w(i,j,k)) + v(i,j,k)**2.0_db*(2.0_db &
			       + 6.0_db*(-1.0_db + w(i,j,k))*w(i,j,k)) - u(i,j,k)**2.0_db*(1.0_db &
			       + 3.0_db*v(i,j,k)**2.0_db + 3.0_db*(-1.0_db + w(i,j,k))*w(i,j,k) &
			       + v(i,j,k)*(-3.0_db + 9.0_db*w(i,j,k)))))/108.0_db
!12
				  fneq1=(3.0_db*(2.0_db*(pyy(i,j,k) + 3.0_db*pyz(i,j,k) &
				   + pzz(i,j,k)) + 6.0_db*(pxy(i,j,k) + pxz(i,j,k))*u(i,j,k) &
				   - 3.0_db*(pyy(i,j,k) + 3.0_db*pyz(i,j,k) + pzz(i,j,k))*u(i,j,k)**2.0_db &
				   - 6.0_db*v(i,j,k)*(2.0_db*pyz(i,j,k) + pzz(i,j,k) &
				   + 3.0_db*(pxy(i,j,k) + pxz(i,j,k))*u(i,j,k) &
				   - pzz(i,j,k)*v(i,j,k)) - 6.0_db*(pyy(i,j,k) + 3.0_db*(2.0_db*pxy(i,j,k) &
				   + pxz(i,j,k))*u(i,j,k) + pyz(i,j,k)*(2.0_db - 6.0_db*v(i,j,k)))*w(i,j,k) &
				   + 12.0_db*pyy(i,j,k)*w(i,j,k)**2.0_db - pxx(i,j,k)*(1.0_db &
				   + 6.0_db*v(i,j,k)**2.0_db - 3.0_db*w(i,j,k) + 6.0_db*w(i,j,k)**2.0_db &
				   + 3.0_db*v(i,j,k)*(-1.0_db + 6.0_db*w(i,j,k)))))/2.0_db
#endif
#ifdef SIXTH_ORDER
!12
			      feq=(2.0_db*rho(i,j,k) + 3.0_db*(2.0_db*(-1.0_db &
			       + w(i,j,k))*w(i,j,k) + v(i,j,k)*(-2.0_db - 6.0_db*(-1.0_db &
			       + w(i,j,k))*w(i,j,k)) - u(i,j,k)**2.0_db*(1.0_db + 3.0_db*(-1.0_db &
			       + v(i,j,k))*v(i,j,k))*(1.0_db + 3.0_db*(-1.0_db + w(i,j,k))*w(i,j,k)) &
			       + v(i,j,k)**2.0_db*(2.0_db + 6.0_db*(-1.0_db + w(i,j,k))*w(i,j,k))))/108.0_db
!12
#ifdef GHOSTONE   
				  fneq1=(-3.0_db*pxx(i,j,k))/2.0_db + 3.0_db*(pyy(i,j,k) &
				   + 3.0_db*pyz(i,j,k) + pzz(i,j,k))
#else 
				  fneq1=(-3.0_db*pxx(i,j,k) + 6.0_db*(pyy(i,j,k) + 3.0_db*pyz(i,j,k) &
				   + pzz(i,j,k)) + 18.0_db*(pxy(i,j,k) + pxz(i,j,k))*u(i,j,k) &
				   - 9.0_db*(pyy(i,j,k) + 3.0_db*pyz(i,j,k) + pzz(i,j,k))*u(i,j,k)**2.0_db &
				   + 9.0_db*(pxx(i,j,k) - 2.0_db*(2.0_db*pyz(i,j,k) + pzz(i,j,k)) &
				   - 6.0_db*(pxy(i,j,k) + pxz(i,j,k))*u(i,j,k) + 3.0_db*(2.0_db*pyz(i,j,k) &
				   + pzz(i,j,k))*u(i,j,k)**2.0_db)*v(i,j,k) + 9.0_db*(-2.0_db*pxx(i,j,k) &
				   + 6.0_db*pxz(i,j,k)*u(i,j,k) + pzz(i,j,k)*(2.0_db &
				   - 3.0_db*u(i,j,k)**2.0_db))*v(i,j,k)**2.0_db &
				   + 9.0_db*(-6.0_db*(2.0_db*pxy(i,j,k) + pxz(i,j,k))*u(i,j,k) &
				   + pyy(i,j,k)*(-2.0_db + 6.0_db*u(i,j,k)**2.0_db) + pyz(i,j,k)*( &
				   -4.0_db + 9.0_db*u(i,j,k)**2.0_db) + pxx(i,j,k)*(1.0_db &
				   - 3.0_db*v(i,j,k))**2.0_db + 6.0_db*(5.0_db*pxy(i,j,k)*u(i,j,k) &
				   + pyz(i,j,k)*(2.0_db - 3.0_db*u(i,j,k)**2.0_db) &
				   - 3.0_db*pxz(i,j,k)*u(i,j,k)*(-1.0_db + v(i,j,k)))*v(i,j,k))*w(i,j,k) &
				   - 9.0_db*(pyy(i,j,k)*(-4.0_db + 9.0_db*u(i,j,k)**2.0_db) &
				   + 6.0_db*pxy(i,j,k)*u(i,j,k)*(-3.0_db + 7.0_db*v(i,j,k)) &
				   + pxx(i,j,k)*(2.0_db + 3.0_db*v(i,j,k)*(-3.0_db &
				   + 4.0_db*v(i,j,k))))*w(i,j,k)**2.0_db)/2.0_db
#endif
#endif
!12
#ifdef PHASE_CHANGE
                  feq=feq+p2*src(i,j,k)*(invrho_b-invrho_r)
#endif
				  F_discr=(forcey + forcez + forcex*u(i,j,k) - 2.0_db*forcey*v(i,j,k) &
				   - 3.0_db*forcez*v(i,j,k) - 3.0_db*forcey*w(i,j,k) &
				   - 2.0_db*forcez*w(i,j,k))/(-18.0_db*rhophi_loc)
                  f(i,j-1,k-1,12)=feq + (1.0_db-omega_loc)*fneq1*p2 + 0.5_db*(F_discr)


!!!!!!!!!!!!!!!!!!!!!!!!!!13
#ifdef SECOND_ORDER
			      feq=(2.0_db*rho(i,j,k) - 3.0_db*u(i,j,k)**2.0_db &
			       + 6.0_db*(v(i,j,k) + v(i,j,k)**2.0_db &
			       - 3.0_db*v(i,j,k)*w(i,j,k) + (-1.0_db + w(i,j,k))*w(i,j,k)))/108.0_db
!13
				  fneq1=(-3.0_db*pxx(i,j,k))/2.0_db + 3.0_db*(pyy(i,j,k) &
				   - 3.0_db*pyz(i,j,k) + pzz(i,j,k))
#endif
#ifdef FORTH_ORDER
!13
			      feq=(2.0_db*rho(i,j,k) + 3.0_db*(2.0_db*v(i,j,k)*(1.0_db &
			       + v(i,j,k)) - 2.0_db*(1.0_db + 3.0_db*v(i,j,k)*(1.0_db &
			       + v(i,j,k)))*w(i,j,k) + 2.0_db*(1.0_db + 3.0_db*v(i,j,k)*(1.0_db &
			       + v(i,j,k)))*w(i,j,k)**2.0_db - u(i,j,k)**2.0_db*(1.0_db &
			       - 3.0_db*w(i,j,k) + 3.0_db*(v(i,j,k) + v(i,j,k)**2.0_db &
			       - 3.0_db*v(i,j,k)*w(i,j,k) + w(i,j,k)**2.0_db))))/108.0_db
!13
				  fneq1=(3.0_db*(2.0_db*(pyy(i,j,k) - 3.0_db*pyz(i,j,k) &
				   + pzz(i,j,k)) + 6.0_db*(-pxy(i,j,k) + pxz(i,j,k))*u(i,j,k) &
				   - 3.0_db*(pyy(i,j,k) - 3.0_db*pyz(i,j,k) + pzz(i,j,k))*u(i,j,k)**2.0_db &
				   + 6.0_db*v(i,j,k)*(-2.0_db*pyz(i,j,k) + pzz(i,j,k) &
				   - 3.0_db*pxy(i,j,k)*u(i,j,k) + 3.0_db*pxz(i,j,k)*u(i,j,k) &
				   + pzz(i,j,k)*v(i,j,k)) - pxx(i,j,k)*(1.0_db + 3.0_db*v(i,j,k) &
				   + 6.0_db*v(i,j,k)**2.0_db) + 3.0_db*(pxx(i,j,k) - 2.0_db*pyy(i,j,k) &
				   + 4.0_db*pyz(i,j,k) + 12.0_db*pxy(i,j,k)*u(i,j,k) - 6.0_db*pxz(i,j,k)*u(i,j,k) &
				   + 6.0_db*(pxx(i,j,k) + 2.0_db*pyz(i,j,k))*v(i,j,k))*w(i,j,k) &
				   - 6.0_db*(pxx(i,j,k) - 2.0_db*pyy(i,j,k))*w(i,j,k)**2.0_db))/2.0_db
#endif
#ifdef SIXTH_ORDER
!13
			      feq=(2.0_db*rho(i,j,k) - 3.0_db*u(i,j,k)**2.0_db*(1.0_db &
			       + 3.0_db*v(i,j,k)*(1.0_db + v(i,j,k)))*(1.0_db + 3.0_db*(-1.0_db &
			       + w(i,j,k))*w(i,j,k)) + 6.0_db*(v(i,j,k) + (-1.0_db &
			       + w(i,j,k))*w(i,j,k) + 3.0_db*v(i,j,k)*(-1.0_db + w(i,j,k))*w(i,j,k) &
			       + v(i,j,k)**2.0_db*(1.0_db + 3.0_db*(-1.0_db + w(i,j,k))*w(i,j,k))))/108.0_db
!13
#ifdef GHOSTONE   
				  fneq1=(-3.0_db*pxx(i,j,k))/2.0_db + 3.0_db*(pyy(i,j,k) &
				   - 3.0_db*pyz(i,j,k) + pzz(i,j,k))
#else 
				  fneq1=(-3.0_db*pxx(i,j,k) + 6.0_db*(pyy(i,j,k) - 3.0_db*pyz(i,j,k) &
				   + pzz(i,j,k)) + 18.0_db*(-pxy(i,j,k) + pxz(i,j,k))*u(i,j,k) &
				   - 9.0_db*(pyy(i,j,k) - 3.0_db*pyz(i,j,k) + pzz(i,j,k))*u(i,j,k)**2.0_db &
				   - 9.0_db*(pxx(i,j,k) + 4.0_db*pyz(i,j,k) + 3.0_db*(-2.0_db*pyz(i,j,k) &
				   + pzz(i,j,k))*u(i,j,k)**2.0_db - 2.0_db*(pzz(i,j,k) + 3.0_db*(-pxy(i,j,k) &
				   + pxz(i,j,k))*u(i,j,k)))*v(i,j,k) - 9.0_db*(2.0_db*pxx(i,j,k) &
				   - 6.0_db*pxz(i,j,k)*u(i,j,k) + pzz(i,j,k)*(-2.0_db &
				   + 3.0_db*u(i,j,k)**2.0_db))*v(i,j,k)**2.0_db + 9.0_db*(4.0_db*pyz(i,j,k) &
				   + 12.0_db*pxy(i,j,k)*u(i,j,k) - 6.0_db*pxz(i,j,k)*u(i,j,k) &
				   - 9.0_db*pyz(i,j,k)*u(i,j,k)**2.0_db + pyy(i,j,k)*(-2.0_db &
				   + 6.0_db*u(i,j,k)**2.0_db) + pxx(i,j,k)*(1.0_db + 3.0_db*v(i,j,k))**2.0_db &
				   + 6.0_db*v(i,j,k)*(5.0_db*pxy(i,j,k)*u(i,j,k) + pyz(i,j,k)*(2.0_db &
				   - 3.0_db*u(i,j,k)**2.0_db) - 3.0_db*pxz(i,j,k)*u(i,j,k)*(1.0_db &
				   + v(i,j,k))))*w(i,j,k) - 9.0_db*(pyy(i,j,k)*(-4.0_db &
				   + 9.0_db*u(i,j,k)**2.0_db) + 6.0_db*pxy(i,j,k)*u(i,j,k)*(3.0_db &
				   + 7.0_db*v(i,j,k)) + pxx(i,j,k)*(2.0_db + 3.0_db*v(i,j,k)*(3.0_db &
				   + 4.0_db*v(i,j,k))))*w(i,j,k)**2.0_db)/2.0_db
#endif
#endif
!13
#ifdef PHASE_CHANGE
                  feq=feq+p2*src(i,j,k)*(invrho_b-invrho_r)
#endif
				  F_discr=(forcez + forcex*u(i,j,k) + 3.0_db*forcez*v(i,j,k) &
				   - 2.0_db*forcez*w(i,j,k) + forcey*(-1.0_db - 2.0_db*v(i,j,k) &
				   + 3.0_db*w(i,j,k)))/(-18.0_db*rhophi_loc)
                  f(i,j+1,k-1,13)=feq + (1.0_db-omega_loc)*fneq1*p2 + 0.5_db*(F_discr)



!!!!!!!!!!!!!!!!!!!!!!!!!!14
#ifdef SECOND_ORDER
			      feq=(2.0_db*rho(i,j,k) - 3.0_db*u(i,j,k)**2.0_db &
			       + 6.0_db*((-1.0_db + v(i,j,k))*v(i,j,k) + w(i,j,k) &
			       - 3.0_db*v(i,j,k)*w(i,j,k) + w(i,j,k)**2.0_db))/108.0_db
!14
				  fneq1=(-3.0_db*pxx(i,j,k))/2.0_db + 3.0_db*(pyy(i,j,k) &
				   - 3.0_db*pyz(i,j,k) + pzz(i,j,k))
#endif
#ifdef FORTH_ORDER
!14
			      feq=(2.0_db*rho(i,j,k) + 6.0_db*((-1.0_db &
			       + v(i,j,k))*v(i,j,k) + w(i,j,k) + 3.0_db*(-1.0_db &
			       + v(i,j,k))*v(i,j,k)*w(i,j,k) + (1.0_db + 3.0_db*( &
			       -1.0_db + v(i,j,k))*v(i,j,k))*w(i,j,k)**2.0_db) &
			       - 3.0_db*u(i,j,k)**2.0_db*(1.0_db &
			       + 3.0_db*v(i,j,k)**2.0_db + 3.0_db*w(i,j,k)*(1.0_db &
			       + w(i,j,k)) - 3.0_db*v(i,j,k)*(1.0_db + 3.0_db*w(i,j,k))))/108.0_db
!14
				  fneq1=(3.0_db*(2.0_db*(pyy(i,j,k) - 3.0_db*pyz(i,j,k) &
				   + pzz(i,j,k)) + 6.0_db*(pxy(i,j,k) - pxz(i,j,k))*u(i,j,k) &
				   - 3.0_db*(pyy(i,j,k) - 3.0_db*pyz(i,j,k) + pzz(i,j,k))*u(i,j,k)**2.0_db &
				   + 6.0_db*((2.0_db*pyz(i,j,k) + 3.0_db*(-pxy(i,j,k) + pxz(i,j,k))*u(i,j,k) &
				   + pzz(i,j,k)*(-1.0_db + v(i,j,k)))*v(i,j,k) + (pyy(i,j,k) &
				   + 6.0_db*pxy(i,j,k)*u(i,j,k) - 3.0_db*pxz(i,j,k)*u(i,j,k) &
				   + pyz(i,j,k)*(-2.0_db + 6.0_db*v(i,j,k)))*w(i,j,k) &
				   + 2.0_db*pyy(i,j,k)*w(i,j,k)**2.0_db) - pxx(i,j,k)*(1.0_db &
				   + 6.0_db*v(i,j,k)**2.0_db + 3.0_db*w(i,j,k) + 6.0_db*w(i,j,k)**2.0_db &
				   - 3.0_db*v(i,j,k)*(1.0_db + 6.0_db*w(i,j,k)))))/2.0_db
#endif
#ifdef SIXTH_ORDER
!14
			      feq=(2.0_db*rho(i,j,k) + 6.0_db*((-1.0_db &
			       + v(i,j,k))*v(i,j,k) + w(i,j,k) + 3.0_db*(-1.0_db &
			       + v(i,j,k))*v(i,j,k)*w(i,j,k) + (1.0_db + 3.0_db*(-1.0_db &
			       + v(i,j,k))*v(i,j,k))*w(i,j,k)**2.0_db) &
			       - 3.0_db*u(i,j,k)**2.0_db*(1.0_db + 3.0_db*(-1.0_db &
			       + v(i,j,k))*v(i,j,k))*(1.0_db + 3.0_db*w(i,j,k)*(1.0_db &
			       + w(i,j,k))))/108.0_db
!14
#ifdef GHOSTONE   
				  fneq1=(-3.0_db*pxx(i,j,k))/2.0_db + 3.0_db*(pyy(i,j,k) &
				   - 3.0_db*pyz(i,j,k) + pzz(i,j,k))
#else 
				  fneq1=(3.0_db*(2.0_db*(pyy(i,j,k) - 3.0_db*pyz(i,j,k) + pzz(i,j,k)) &
				   + 6.0_db*(pxy(i,j,k) - pxz(i,j,k))*u(i,j,k) - 3.0_db*(pyy(i,j,k) &
				   - 3.0_db*pyz(i,j,k) + pzz(i,j,k))*u(i,j,k)**2.0_db &
				   - 3.0_db*(6.0_db*u(i,j,k)*(pxy(i,j,k) + pxz(i,j,k)*(-1.0_db + v(i,j,k))) &
				   + pzz(i,j,k)*(-2.0_db + 3.0_db*u(i,j,k)**2.0_db)*(-1.0_db + v(i,j,k)))*v(i,j,k) &
				   + 6.0_db*pyy(i,j,k)*w(i,j,k) - 18.0_db*u(i,j,k)*(pxz(i,j,k) + pyy(i,j,k)*u(i,j,k) &
				   + 3.0_db*pxz(i,j,k)*(-1.0_db + v(i,j,k))*v(i,j,k) + pxy(i,j,k)*(-2.0_db &
				   + 5.0_db*v(i,j,k)))*w(i,j,k) + 3.0_db*(pyy(i,j,k)*(4.0_db &
				   - 9.0_db*u(i,j,k)**2.0_db) + 6.0_db*pxy(i,j,k)*u(i,j,k)*(3.0_db &
				   - 7.0_db*v(i,j,k)))*w(i,j,k)**2.0_db - 3.0_db*pyz(i,j,k)*((4.0_db &
				   - 9.0_db*u(i,j,k)**2.0_db)*w(i,j,k) + 2.0_db*(-2.0_db &
				   + 3.0_db*u(i,j,k)**2.0_db)*v(i,j,k)*(1.0_db + 3.0_db*w(i,j,k))) &
				   + pxx(i,j,k)*(-1.0_db - 3.0_db*w(i,j,k)*(1.0_db + 2.0_db*w(i,j,k)) &
				   + 3.0_db*v(i,j,k)*(1.0_db + 3.0_db*w(i,j,k))**2.0_db &
				   - 3.0_db*v(i,j,k)**2.0_db*(2.0_db + 3.0_db*w(i,j,k)*(3.0_db &
				   + 4.0_db*w(i,j,k))))))/2.0_db
#endif
#endif
!14
#ifdef PHASE_CHANGE
                  feq=feq+p2*src(i,j,k)*(invrho_b-invrho_r)
#endif
				  F_discr=(forcey - forcez + forcex*u(i,j,k) - 2.0_db*forcey*v(i,j,k) &
				   + 3.0_db*forcez*v(i,j,k) + 3.0_db*forcey*w(i,j,k) &
				   - 2.0_db*forcez*w(i,j,k))/(-18.0_db*rhophi_loc)
                  f(i,j-1,k+1,14)=feq + (1.0_db-omega_loc)*fneq1*p2 + 0.5_db*(F_discr)


!!!!!!!!!!!!!!!!!!!!!!!!!!15
#ifdef SECOND_ORDER
			      feq=(2.0_db*rho(i,j,k) - 3.0_db*v(i,j,k)**2.0_db &
			       + 6.0_db*w(i,j,k) + 6.0_db*(u(i,j,k) &
			       + u(i,j,k)**2.0_db + 3.0_db*u(i,j,k)*w(i,j,k) + w(i,j,k)**2.0_db))/108.0_db
!15
				  fneq1=3.0_db*pxx(i,j,k) + 9.0_db*pxz(i,j,k) &
				   - (3.0_db*pyy(i,j,k))/2.0_db + 3.0_db*pzz(i,j,k)
#endif
#ifdef FORTH_ORDER
!15
			      feq=(2.0_db*rho(i,j,k) + 3.0_db*(2.0_db*w(i,j,k)*(1.0_db &
			       + w(i,j,k)) - v(i,j,k)**2.0_db*(1.0_db + 3.0_db*w(i,j,k)*(1.0_db &
			       + w(i,j,k))) + u(i,j,k)**2.0_db*(2.0_db - 3.0_db*v(i,j,k)**2.0_db &
			       + 6.0_db*w(i,j,k)*(1.0_db + w(i,j,k))) + u(i,j,k)*(2.0_db &
			       + 6.0_db*w(i,j,k)*(1.0_db + w(i,j,k)) &
			       - 3.0_db*v(i,j,k)**2.0_db*(1.0_db + 3.0_db*w(i,j,k)))))/108.0_db
!15
				  fneq1=(3.0_db*(2.0_db*pxx(i,j,k) + 6.0_db*pxz(i,j,k)*(1.0_db &
				   + 2.0_db*u(i,j,k)) - (pyy(i,j,k) - 2.0_db*pzz(i,j,k))*(1.0_db &
				   + 3.0_db*u(i,j,k)*(1.0_db + u(i,j,k))) - 6.0_db*(pxy(i,j,k) &
				   + pyz(i,j,k))*(1.0_db + 3.0_db*u(i,j,k))*v(i,j,k) &
				   - 3.0_db*(2.0_db*pxx(i,j,k) + 3.0_db*pxz(i,j,k) &
				   + pzz(i,j,k))*v(i,j,k)**2.0_db + 3.0_db*(2.0_db*pxx(i,j,k) &
				   + 4.0_db*pxz(i,j,k)*(1.0_db + 3.0_db*u(i,j,k)) - pyy(i,j,k)*(1.0_db &
				   + 6.0_db*u(i,j,k)) - 6.0_db*(2.0_db*pxy(i,j,k) + pyz(i,j,k))*v(i,j,k))*w(i,j,k) &
				   + 6.0_db*(2.0_db*pxx(i,j,k) - pyy(i,j,k))*w(i,j,k)**2.0_db))/2.0_db

#endif
#ifdef SIXTH_ORDER
!15
			      feq=(2.0_db*rho(i,j,k) + 3.0_db*(2.0_db*w(i,j,k)*(1.0_db &
			       + w(i,j,k)) - v(i,j,k)**2.0_db*(1.0_db + 3.0_db*w(i,j,k)*(1.0_db &
			       + w(i,j,k))) - u(i,j,k)*(-2.0_db + 3.0_db*v(i,j,k)**2.0_db)*(1.0_db &
			       + 3.0_db*w(i,j,k)*(1.0_db + w(i,j,k))) - u(i,j,k)**2.0_db*(-2.0_db &
			       + 3.0_db*v(i,j,k)**2.0_db)*(1.0_db + 3.0_db*w(i,j,k)*(1.0_db &
			       + w(i,j,k)))))/108.0_db
!15
#ifdef GHOSTONE   
				  fneq1=3.0_db*pxx(i,j,k) + 9.0_db*pxz(i,j,k) &
				   - (3.0_db*pyy(i,j,k))/2.0_db + 3.0_db*pzz(i,j,k)
#else 
				  fneq1=(3.0_db*(2.0_db*pxx(i,j,k) + 6.0_db*pxz(i,j,k)*(1.0_db &
				   + 2.0_db*u(i,j,k)) - (pyy(i,j,k) - 2.0_db*pzz(i,j,k))*(1.0_db &
				   + 3.0_db*u(i,j,k)*(1.0_db + u(i,j,k))) - 6.0_db*(pxy(i,j,k) &
				   + pyz(i,j,k) + 3.0_db*pxy(i,j,k)*u(i,j,k) + 3.0_db*pyz(i,j,k)*u(i,j,k)*(1.0_db &
				   + u(i,j,k)))*v(i,j,k) - 3.0_db*(2.0_db*pxx(i,j,k) + pzz(i,j,k) &
				   + 3.0_db*pzz(i,j,k)*u(i,j,k)*(1.0_db + u(i,j,k)) + pxz(i,j,k)*(3.0_db &
				   + 6.0_db*u(i,j,k)))*v(i,j,k)**2.0_db - 3.0_db*(-2.0_db*pxx(i,j,k) &
				   + pyy(i,j,k) + 6.0_db*pyy(i,j,k)*u(i,j,k)*(1.0_db + u(i,j,k)) &
				   - 4.0_db*pxz(i,j,k)*(1.0_db + 3.0_db*u(i,j,k)) + 6.0_db*(pyz(i,j,k) &
				   + 3.0_db*pyz(i,j,k)*u(i,j,k)*(1.0_db + u(i,j,k)) + pxy(i,j,k)*(2.0_db &
				   + 5.0_db*u(i,j,k)))*v(i,j,k) + 9.0_db*(pxx(i,j,k) + pxz(i,j,k) &
				   + 2.0_db*pxz(i,j,k)*u(i,j,k))*v(i,j,k)**2.0_db)*w(i,j,k) &
				   - 3.0_db*(pyy(i,j,k)*(2.0_db + 9.0_db*u(i,j,k)*(1.0_db + u(i,j,k))) &
				   + 6.0_db*pxy(i,j,k)*(3.0_db + 7.0_db*u(i,j,k))*v(i,j,k) &
				   + 4.0_db*pxx(i,j,k)*(-1.0_db + 3.0_db*v(i,j,k)**2.0_db))*w(i,j,k)**2.0_db))/2.0_db
#endif
#endif
!15
#ifdef PHASE_CHANGE
                  feq=feq+p2*src(i,j,k)*(invrho_b-invrho_r)
#endif
				  F_discr=(forcex + forcez + 2.0_db*forcex*u(i,j,k) &
				   + 3.0_db*forcez*u(i,j,k) - forcey*v(i,j,k) + 3.0_db*forcex*w(i,j,k) &
				   + 2.0_db*forcez*w(i,j,k))/(18.0_db*rhophi_loc)
                  f(i+1,j,k+1,15)=feq + (1.0_db-omega_loc)*fneq1*p2 + 0.5_db*(F_discr)



!!!!!!!!!!!!!!!!!!!!!!!!!!16
#ifdef SECOND_ORDER
			      feq=(2.0_db*rho(i,j,k) + 3.0_db*(2.0_db*u(i,j,k)**2.0_db &
			       - v(i,j,k)**2.0_db + 2.0_db*(-1.0_db + w(i,j,k))*w(i,j,k) &
			       + u(i,j,k)*(-2.0_db + 6.0_db*w(i,j,k))))/108.0_db
!16
				  fneq1=3.0_db*pxx(i,j,k) + 9.0_db*pxz(i,j,k) &
				   - (3.0_db*pyy(i,j,k))/2.0_db + 3.0_db*pzz(i,j,k)
#endif
#ifdef FORTH_ORDER
!16
			      feq=(2.0_db*rho(i,j,k) + 3.0_db*(2.0_db*(-1.0_db &
			       + w(i,j,k))*w(i,j,k) + u(i,j,k)*(-2.0_db + v(i,j,k)**2.0_db*(3.0_db &
			       - 9.0_db*w(i,j,k)) - 6.0_db*(-1.0_db + w(i,j,k))*w(i,j,k)) &
			       + v(i,j,k)**2.0_db*(-1.0_db - 3.0_db*(-1.0_db + w(i,j,k))*w(i,j,k)) &
			       + u(i,j,k)**2.0_db*(2.0_db - 3.0_db*v(i,j,k)**2.0_db + 6.0_db*( &
			       -1.0_db + w(i,j,k))*w(i,j,k))))/108.0_db
!16
				  fneq1=(3.0_db*(2.0_db*pxx(i,j,k) + pxz(i,j,k)*(6.0_db &
				   - 12.0_db*u(i,j,k)) - (pyy(i,j,k) - 2.0_db*pzz(i,j,k))*(1.0_db &
				   + 3.0_db*(-1.0_db + u(i,j,k))*u(i,j,k)) - 6.0_db*(pxy(i,j,k) &
				   + pyz(i,j,k))*(-1.0_db + 3.0_db*u(i,j,k))*v(i,j,k) &
				   - 3.0_db*(2.0_db*pxx(i,j,k) + 3.0_db*pxz(i,j,k) &
				   + pzz(i,j,k))*v(i,j,k)**2.0_db - 3.0_db*(2.0_db*pxx(i,j,k) &
				   + pxz(i,j,k)*(4.0_db - 12.0_db*u(i,j,k)) + pyy(i,j,k)*(-1.0_db &
				   + 6.0_db*u(i,j,k)) + 6.0_db*(2.0_db*pxy(i,j,k) + pyz(i,j,k))*v(i,j,k))*w(i,j,k) &
				   + 6.0_db*(2.0_db*pxx(i,j,k) - pyy(i,j,k))*w(i,j,k)**2.0_db))/2.0_db
#endif
#ifdef SIXTH_ORDER
!16
			      feq=(2.0_db*rho(i,j,k) + 3.0_db*(2.0_db*(-1.0_db &
			       + w(i,j,k))*w(i,j,k) + v(i,j,k)**2.0_db*(-1.0_db - 3.0_db*(-1.0_db &
			       + w(i,j,k))*w(i,j,k)) + u(i,j,k)*(-2.0_db + 3.0_db*v(i,j,k)**2.0_db)*(1.0_db &
			       + 3.0_db*(-1.0_db + w(i,j,k))*w(i,j,k)) - u(i,j,k)**2.0_db*(-2.0_db &
			       + 3.0_db*v(i,j,k)**2.0_db)*(1.0_db + 3.0_db*(-1.0_db &
			       + w(i,j,k))*w(i,j,k))))/108.0_db
!16
#ifdef GHOSTONE   
				  fneq1=3.0_db*pxx(i,j,k) + 9.0_db*pxz(i,j,k) &
				   - (3.0_db*pyy(i,j,k))/2.0_db + 3.0_db*pzz(i,j,k)
#else 
				  fneq1=(3.0_db*(2.0_db*pxx(i,j,k) + pxz(i,j,k)*(6.0_db &
				   - 12.0_db*u(i,j,k)) - (pyy(i,j,k) - 2.0_db*pzz(i,j,k))*(1.0_db &
				   + 3.0_db*(-1.0_db + u(i,j,k))*u(i,j,k)) + 6.0_db*(pxy(i,j,k) &
				   + pyz(i,j,k) - 3.0_db*pxy(i,j,k)*u(i,j,k) + 3.0_db*pyz(i,j,k)*( &
				   -1.0_db + u(i,j,k))*u(i,j,k))*v(i,j,k) - 3.0_db*(2.0_db*pxx(i,j,k) &
				   + pzz(i,j,k) + pxz(i,j,k)*(3.0_db - 6.0_db*u(i,j,k)) &
				   + 3.0_db*pzz(i,j,k)*(-1.0_db + u(i,j,k))*u(i,j,k))*v(i,j,k)**2.0_db &
				   + 3.0_db*(-2.0_db*pxx(i,j,k) + pyy(i,j,k) + 4.0_db*pxz(i,j,k)*( &
				   -1.0_db + 3.0_db*u(i,j,k)) - 6.0_db*(2.0_db*pxy(i,j,k) &
				   + pyz(i,j,k))*v(i,j,k) + 9.0_db*(pxx(i,j,k) + pxz(i,j,k) &
				   - 2.0_db*pxz(i,j,k)*u(i,j,k))*v(i,j,k)**2.0_db + 6.0_db*u(i,j,k)*(pyy(i,j,k)*( &
				   -1.0_db + u(i,j,k)) + 5.0_db*pxy(i,j,k)*v(i,j,k) &
				   - 3.0_db*pyz(i,j,k)*(-1.0_db + u(i,j,k))*v(i,j,k)))*w(i,j,k) &
				   - 3.0_db*(pyy(i,j,k)*(2.0_db + 9.0_db*(-1.0_db + u(i,j,k))*u(i,j,k)) &
				   + 6.0_db*pxy(i,j,k)*(-3.0_db + 7.0_db*u(i,j,k))*v(i,j,k) + 4.0_db*pxx(i,j,k)*( &
				   -1.0_db + 3.0_db*v(i,j,k)**2.0_db))*w(i,j,k)**2.0_db))/2.0_db
#endif
#endif
!16
#ifdef PHASE_CHANGE
                  feq=feq+p2*src(i,j,k)*(invrho_b-invrho_r)
#endif
				  F_discr=(forcex + forcez - 2.0_db*forcex*u(i,j,k) &
				   - 3.0_db*forcez*u(i,j,k) + forcey*v(i,j,k) - 3.0_db*forcex*w(i,j,k) &
				   - 2.0_db*forcez*w(i,j,k))/(-18.0_db*rhophi_loc)
                  f(i-1,j,k-1,16)=feq + (1.0_db-omega_loc)*fneq1*p2 + 0.5_db*(F_discr)



!!!!!!!!!!!!!!!!!!!!!!!!!!17
#ifdef SECOND_ORDER
			      feq=(2.0_db*rho(i,j,k) + 3.0_db*(2.0_db*u(i,j,k)**2.0_db &
			       - v(i,j,k)**2.0_db + 2.0_db*w(i,j,k)*(1.0_db + w(i,j,k)) &
			       - 2.0_db*u(i,j,k)*(1.0_db + 3.0_db*w(i,j,k))))/108.0_db
!17
				  fneq1=3.0_db*pxx(i,j,k) - 9.0_db*pxz(i,j,k) &
				   - (3.0_db*pyy(i,j,k))/2.0_db + 3.0_db*pzz(i,j,k)
#endif
#ifdef FORTH_ORDER
!17
			      feq=(2.0_db*rho(i,j,k) + 3.0_db*(2.0_db*w(i,j,k)*(1.0_db &
			       + w(i,j,k)) - v(i,j,k)**2.0_db*(1.0_db + 3.0_db*w(i,j,k)*(1.0_db &
			       + w(i,j,k))) + u(i,j,k)**2.0_db*(2.0_db - 3.0_db*v(i,j,k)**2.0_db &
			       + 6.0_db*w(i,j,k)*(1.0_db + w(i,j,k))) + u(i,j,k)*(-2.0_db &
			       - 6.0_db*w(i,j,k)*(1.0_db + w(i,j,k)) + v(i,j,k)**2.0_db*(3.0_db &
			       + 9.0_db*w(i,j,k)))))/108.0_db
!17
				  fneq1=(3.0_db*(6.0_db*pxz(i,j,k)*(-1.0_db + 2.0_db*u(i,j,k)) &
				   - (pyy(i,j,k) - 2.0_db*pzz(i,j,k))*(1.0_db + 3.0_db*(-1.0_db &
				   + u(i,j,k))*u(i,j,k)) - 6.0_db*(pxy(i,j,k) - pyz(i,j,k))*(-1.0_db &
				   + 3.0_db*u(i,j,k))*v(i,j,k) + 3.0_db*(3.0_db*pxz(i,j,k) &
				   - pzz(i,j,k))*v(i,j,k)**2.0_db + 3.0_db*(-4.0_db*pxz(i,j,k) &
				   - pyy(i,j,k) + 6.0_db*(2.0_db*pxz(i,j,k) + pyy(i,j,k))*u(i,j,k) &
				   + 12.0_db*pxy(i,j,k)*v(i,j,k) - 6.0_db*pyz(i,j,k)*v(i,j,k))*w(i,j,k) &
				   - 6.0_db*pyy(i,j,k)*w(i,j,k)**2.0_db + pxx(i,j,k)*(2.0_db &
				   - 6.0_db*v(i,j,k)**2.0_db + 6.0_db*w(i,j,k)*(1.0_db &
				   + 2.0_db*w(i,j,k)))))/2.0_db
#endif
#ifdef SIXTH_ORDER
!17
			      feq=(2.0_db*rho(i,j,k) + 3.0_db*(2.0_db*w(i,j,k)*(1.0_db &
			       + w(i,j,k)) - v(i,j,k)**2.0_db*(1.0_db + 3.0_db*w(i,j,k)*(1.0_db &
			       + w(i,j,k))) + u(i,j,k)*(-2.0_db + 3.0_db*v(i,j,k)**2.0_db)*(1.0_db &
			       + 3.0_db*w(i,j,k)*(1.0_db + w(i,j,k))) - u(i,j,k)**2.0_db*(-2.0_db &
			       + 3.0_db*v(i,j,k)**2.0_db)*(1.0_db + 3.0_db*w(i,j,k)*(1.0_db &
			       + w(i,j,k)))))/108.0_db
!17
#ifdef GHOSTONE   
				  fneq1=3.0_db*pxx(i,j,k) - 9.0_db*pxz(i,j,k) &
				   - (3.0_db*pyy(i,j,k))/2.0_db + 3.0_db*pzz(i,j,k)
#else 
				  fneq1=(3.0_db*(2.0_db*pxx(i,j,k) + 6.0_db*pxz(i,j,k)*(-1.0_db &
				   + 2.0_db*u(i,j,k)) - (pyy(i,j,k) - 2.0_db*pzz(i,j,k))*(1.0_db &
				   + 3.0_db*(-1.0_db + u(i,j,k))*u(i,j,k)) - 6.0_db*(pyz(i,j,k) &
				   + 3.0_db*pyz(i,j,k)*(-1.0_db + u(i,j,k))*u(i,j,k) + pxy(i,j,k)*(-1.0_db &
				   + 3.0_db*u(i,j,k)))*v(i,j,k) - 3.0_db*(2.0_db*pxx(i,j,k) + pzz(i,j,k) &
				   + 3.0_db*pzz(i,j,k)*(-1.0_db + u(i,j,k))*u(i,j,k) + pxz(i,j,k)*(-3.0_db &
				   + 6.0_db*u(i,j,k)))*v(i,j,k)**2.0_db - 3.0_db*(-2.0_db*pxx(i,j,k) &
				   + pyy(i,j,k) + pxz(i,j,k)*(4.0_db - 12.0_db*u(i,j,k)) &
				   + 6.0_db*pyy(i,j,k)*(-1.0_db + u(i,j,k))*u(i,j,k) &
				   + 6.0_db*(pyz(i,j,k) + 3.0_db*pyz(i,j,k)*(-1.0_db &
				   + u(i,j,k))*u(i,j,k) + pxy(i,j,k)*(-2.0_db + 5.0_db*u(i,j,k)))*v(i,j,k) &
				   + 9.0_db*(pxx(i,j,k) + pxz(i,j,k)*(-1.0_db &
				   + 2.0_db*u(i,j,k)))*v(i,j,k)**2.0_db)*w(i,j,k) &
				   - 3.0_db*(pyy(i,j,k)*(2.0_db + 9.0_db*(-1.0_db + u(i,j,k))*u(i,j,k)) &
				   + 6.0_db*pxy(i,j,k)*(-3.0_db + 7.0_db*u(i,j,k))*v(i,j,k) &
				   + 4.0_db*pxx(i,j,k)*(-1.0_db + 3.0_db*v(i,j,k)**2.0_db))*w(i,j,k)**2.0_db))/2.0_db
#endif
#endif
!17
#ifdef PHASE_CHANGE
                  feq=feq+p2*src(i,j,k)*(invrho_b-invrho_r)
#endif
				  F_discr=(forcex - forcez - 2.0_db*forcex*u(i,j,k) &
				   + 3.0_db*forcez*u(i,j,k) + forcey*v(i,j,k) + 3.0_db*forcex*w(i,j,k) &
				   - 2.0_db*forcez*w(i,j,k))/(-18.0_db*rhophi_loc)
                  f(i-1,j,k+1,17)=feq + (1.0_db-omega_loc)*fneq1*p2 + 0.5_db*(F_discr)



!!!!!!!!!!!!!!!!!!!!!!!!!!18
#ifdef SECOND_ORDER
			      feq=(2.0_db*rho(i,j,k) - 3.0_db*v(i,j,k)**2.0_db &
			       - 6.0_db*w(i,j,k) + 6.0_db*(u(i,j,k) &
			       + u(i,j,k)**2.0_db - 3.0_db*u(i,j,k)*w(i,j,k) + w(i,j,k)**2.0_db))/108.0_db
!18
				  fneq1=3.0_db*pxx(i,j,k) - 9.0_db*pxz(i,j,k) &
				   - (3.0_db*pyy(i,j,k))/2.0_db + 3.0_db*pzz(i,j,k)
#endif
#ifdef FORTH_ORDER
!18
			      feq=(2.0_db*rho(i,j,k) + 3.0_db*(2.0_db*(-1.0_db &
			       + w(i,j,k))*w(i,j,k) + v(i,j,k)**2.0_db*(-1.0_db - 3.0_db*(-1.0_db &
			       + w(i,j,k))*w(i,j,k)) + u(i,j,k)**2.0_db*(2.0_db - 3.0_db*v(i,j,k)**2.0_db &
			       + 6.0_db*(-1.0_db + w(i,j,k))*w(i,j,k)) + u(i,j,k)*(2.0_db + 6.0_db*(-1.0_db &
			       + w(i,j,k))*w(i,j,k) + v(i,j,k)**2.0_db*(-3.0_db + 9.0_db*w(i,j,k)))))/108.0_db
!18
				  fneq1=(3.0_db*(2.0_db*pxx(i,j,k) - 6.0_db*pxz(i,j,k)*(1.0_db &
				   + 2.0_db*u(i,j,k)) - (pyy(i,j,k) - 2.0_db*pzz(i,j,k))*(1.0_db &
				   + 3.0_db*u(i,j,k)*(1.0_db + u(i,j,k))) - 6.0_db*(pxy(i,j,k) &
				   - pyz(i,j,k))*(1.0_db + 3.0_db*u(i,j,k))*v(i,j,k) &
				   - 3.0_db*(2.0_db*pxx(i,j,k) - 3.0_db*pxz(i,j,k) &
				   + pzz(i,j,k))*v(i,j,k)**2.0_db + 3.0_db*(-2.0_db*pxx(i,j,k) &
				   + 4.0_db*pxz(i,j,k) + pyy(i,j,k) + 6.0_db*(2.0_db*pxz(i,j,k) &
				   + pyy(i,j,k))*u(i,j,k) + 12.0_db*pxy(i,j,k)*v(i,j,k) &
				   - 6.0_db*pyz(i,j,k)*v(i,j,k))*w(i,j,k) + 6.0_db*(2.0_db*pxx(i,j,k) &
				   - pyy(i,j,k))*w(i,j,k)**2.0_db))/2.0_db
#endif
#ifdef SIXTH_ORDER
!18
			      feq=(2.0_db*rho(i,j,k) + 3.0_db*(2.0_db*(-1.0_db &
			       + w(i,j,k))*w(i,j,k) + v(i,j,k)**2.0_db*(-1.0_db &
			       - 3.0_db*(-1.0_db + w(i,j,k))*w(i,j,k)) - u(i,j,k)*(-2.0_db &
			       + 3.0_db*v(i,j,k)**2.0_db)*(1.0_db + 3.0_db*(-1.0_db &
			       + w(i,j,k))*w(i,j,k)) - u(i,j,k)**2.0_db*(-2.0_db &
			       + 3.0_db*v(i,j,k)**2.0_db)*(1.0_db + 3.0_db*(-1.0_db &
			       + w(i,j,k))*w(i,j,k))))/108.0_db
!18
#ifdef GHOSTONE   
				  fneq1=3.0_db*pxx(i,j,k) - 9.0_db*pxz(i,j,k) &
				   - (3.0_db*pyy(i,j,k))/2.0_db + 3.0_db*pzz(i,j,k)
#else 
				  fneq1=(3.0_db*(2.0_db*pxx(i,j,k) - 6.0_db*pxz(i,j,k)*(1.0_db &
				   + 2.0_db*u(i,j,k)) - (pyy(i,j,k) - 2.0_db*pzz(i,j,k))*(1.0_db &
				   + 3.0_db*u(i,j,k)*(1.0_db + u(i,j,k))) + 6.0_db*(pyz(i,j,k) &
				   + 3.0_db*pyz(i,j,k)*u(i,j,k)*(1.0_db + u(i,j,k)) - pxy(i,j,k)*(1.0_db &
				   + 3.0_db*u(i,j,k)))*v(i,j,k) - 3.0_db*(2.0_db*pxx(i,j,k) + pzz(i,j,k) &
				   + 3.0_db*pzz(i,j,k)*u(i,j,k)*(1.0_db + u(i,j,k)) - 3.0_db*pxz(i,j,k)*(1.0_db &
				   + 2.0_db*u(i,j,k)))*v(i,j,k)**2.0_db + 3.0_db*(-2.0_db*pxx(i,j,k) &
				   + pyy(i,j,k) + 6.0_db*pyy(i,j,k)*u(i,j,k)*(1.0_db + u(i,j,k)) &
				   + 4.0_db*pxz(i,j,k)*(1.0_db + 3.0_db*u(i,j,k)) + 6.0_db*pxy(i,j,k)*(2.0_db &
				   + 5.0_db*u(i,j,k))*v(i,j,k) - 6.0_db*pyz(i,j,k)*(1.0_db &
				   + 3.0_db*u(i,j,k)*(1.0_db + u(i,j,k)))*v(i,j,k) &
				   + 9.0_db*(pxx(i,j,k) - pxz(i,j,k)*(1.0_db &
				   + 2.0_db*u(i,j,k)))*v(i,j,k)**2.0_db)*w(i,j,k) &
				   - 3.0_db*(pyy(i,j,k)*(2.0_db + 9.0_db*u(i,j,k)*(1.0_db &
				   + u(i,j,k))) + 6.0_db*pxy(i,j,k)*(3.0_db + 7.0_db*u(i,j,k))*v(i,j,k) &
				   + 4.0_db*pxx(i,j,k)*(-1.0_db + 3.0_db*v(i,j,k)**2.0_db))*w(i,j,k)**2.0_db))/2.0_db
#endif
#endif
!18
#ifdef PHASE_CHANGE
                  feq=feq+p2*src(i,j,k)*(invrho_b-invrho_r)
#endif
				  F_discr=(forcez + 3.0_db*forcez*u(i,j,k) + forcey*v(i,j,k) &
				   - 2.0_db*forcez*w(i,j,k) + forcex*(-1.0_db - 2.0_db*u(i,j,k) &
				   + 3.0_db*w(i,j,k)))/(-18.0_db*rhophi_loc)
                  f(i+1,j,k-1,18)=feq + (1.0_db-omega_loc)*fneq1*p2 + 0.5_db*(F_discr)



!!!!!!!!!!!!!!!!!!!!!!!!!!19
#ifdef SECOND_ORDER
			      feq=(rho(i,j,k) + 3.0_db*(u(i,j,k) + u(i,j,k)**2.0_db &
			       + v(i,j,k) + 3.0_db*u(i,j,k)*v(i,j,k) + v(i,j,k)**2.0_db &
			       + w(i,j,k) + 3.0_db*(u(i,j,k) + v(i,j,k))*w(i,j,k) &
			       + w(i,j,k)**2.0_db))/216.0_db
!19
				  fneq1=3.0_db*(pxx(i,j,k) + 3.0_db*pxy(i,j,k) &
				   + 3.0_db*pxz(i,j,k) + pyy(i,j,k) + 3.0_db*pyz(i,j,k) + pzz(i,j,k))
#endif
#ifdef FORTH_ORDER
!19
			      feq=(rho(i,j,k) + 3.0_db*(v(i,j,k) + v(i,j,k)**2.0_db &
			       + w(i,j,k) + 3.0_db*v(i,j,k)*(1.0_db + v(i,j,k))*w(i,j,k) &
			       + (1.0_db + 3.0_db*v(i,j,k)*(1.0_db + v(i,j,k)))*w(i,j,k)**2.0_db &
			       + u(i,j,k)**2.0_db*(1.0_db + 3.0_db*w(i,j,k) + 3.0_db*(v(i,j,k) &
			       + v(i,j,k)**2.0_db + 3.0_db*v(i,j,k)*w(i,j,k) + w(i,j,k)**2.0_db)) &
			       + u(i,j,k)*(1.0_db + 3.0_db*w(i,j,k)*(1.0_db + w(i,j,k)) &
			       + v(i,j,k)**2.0_db*(3.0_db + 9.0_db*w(i,j,k)) + v(i,j,k)*(3.0_db &
			       + 9.0_db*w(i,j,k)*(1.0_db + w(i,j,k))))))/216.0_db
!19
				  fneq1=3.0_db*(pxy(i,j,k)*(3.0_db + 6.0_db*u(i,j,k)) &
				   + pxz(i,j,k)*(3.0_db + 6.0_db*u(i,j,k)) + (pyy(i,j,k) &
				   + 3.0_db*pyz(i,j,k) + pzz(i,j,k))*(1.0_db + 3.0_db*u(i,j,k)*(1.0_db &
				   + u(i,j,k))) + 3.0_db*v(i,j,k)*((2.0_db*pyz(i,j,k) + pzz(i,j,k))*(1.0_db &
				   + 3.0_db*u(i,j,k)) + pxy(i,j,k)*(2.0_db + 6.0_db*u(i,j,k)) &
				   + pzz(i,j,k)*v(i,j,k) + 3.0_db*pxz(i,j,k)*(1.0_db &
				   + 2.0_db*u(i,j,k) + v(i,j,k))) + 3.0_db*(pyy(i,j,k) &
				   + 2.0_db*pyz(i,j,k) + 6.0_db*pyy(i,j,k)*u(i,j,k) &
				   + 9.0_db*pyz(i,j,k)*u(i,j,k) + 6.0_db*pyz(i,j,k)*v(i,j,k) &
				   + 3.0_db*pxy(i,j,k)*(1.0_db + 4.0_db*u(i,j,k) + 4.0_db*v(i,j,k)) &
				   + pxz(i,j,k)*(2.0_db + 6.0_db*u(i,j,k) + 9.0_db*v(i,j,k)))*w(i,j,k) &
				   + 6.0_db*(3.0_db*pxy(i,j,k) + pyy(i,j,k))*w(i,j,k)**2.0_db &
				   + pxx(i,j,k)*(1.0_db + 6.0_db*v(i,j,k)**2.0_db + 3.0_db*w(i,j,k) &
				   + 6.0_db*w(i,j,k)**2.0_db + 3.0_db*v(i,j,k)*(1.0_db + 6.0_db*w(i,j,k))))
#endif
#ifdef SIXTH_ORDER
!19
			      feq=(rho(i,j,k) + 3.0_db*(v(i,j,k) + v(i,j,k)**2.0_db &
			       + w(i,j,k) + 3.0_db*v(i,j,k)*(1.0_db + v(i,j,k))*w(i,j,k) + (1.0_db &
			       + 3.0_db*v(i,j,k)*(1.0_db + v(i,j,k)))*w(i,j,k)**2.0_db + u(i,j,k)*(1.0_db &
			       + 3.0_db*v(i,j,k)*(1.0_db + v(i,j,k)))*(1.0_db + 3.0_db*w(i,j,k)*(1.0_db &
			       + w(i,j,k))) + u(i,j,k)**2.0_db*(1.0_db + 3.0_db*v(i,j,k)*(1.0_db &
			       + v(i,j,k)))*(1.0_db + 3.0_db*w(i,j,k)*(1.0_db + w(i,j,k)))))/216.0_db
!19
#ifdef GHOSTONE   
				  fneq1=3.0_db*(pxx(i,j,k) + 3.0_db*pxy(i,j,k) &
				   + 3.0_db*pxz(i,j,k) + pyy(i,j,k) + 3.0_db*pyz(i,j,k) + pzz(i,j,k))
#else 
				  fneq1=3.0_db*(pxy(i,j,k)*(3.0_db + 6.0_db*u(i,j,k)) &
				   + pxz(i,j,k)*(3.0_db + 6.0_db*u(i,j,k)) + (pyy(i,j,k) &
				   + 3.0_db*pyz(i,j,k) + pzz(i,j,k))*(1.0_db + 3.0_db*u(i,j,k)*(1.0_db &
				   + u(i,j,k))) + 3.0_db*(v(i,j,k)*(pxy(i,j,k)*(2.0_db + 6.0_db*u(i,j,k)) &
				   + 3.0_db*pxz(i,j,k)*(1.0_db + 2.0_db*u(i,j,k))*(1.0_db + v(i,j,k)) &
				   + (1.0_db + 3.0_db*u(i,j,k)*(1.0_db + u(i,j,k)))*(2.0_db*pyz(i,j,k) &
				   + pzz(i,j,k) + pzz(i,j,k)*v(i,j,k))) + (pyy(i,j,k) &
				   + 2.0_db*pyz(i,j,k) + 3.0_db*(2.0_db*pyy(i,j,k) &
				   + 3.0_db*pyz(i,j,k))*u(i,j,k)*(1.0_db + u(i,j,k)) &
				   + 6.0_db*pyz(i,j,k)*(1.0_db + 3.0_db*u(i,j,k)*(1.0_db &
				   + u(i,j,k)))*v(i,j,k) + 3.0_db*pxy(i,j,k)*(1.0_db + 4.0_db*v(i,j,k) &
				   + 2.0_db*u(i,j,k)*(2.0_db + 5.0_db*v(i,j,k))) &
				   + pxz(i,j,k)*(2.0_db + 9.0_db*v(i,j,k)*(1.0_db + v(i,j,k)) &
				   + 6.0_db*u(i,j,k)*(1.0_db + 3.0_db*v(i,j,k)*(1.0_db &
				   + v(i,j,k)))))*w(i,j,k) + (pyy(i,j,k)*(2.0_db &
				   + 9.0_db*u(i,j,k)*(1.0_db + u(i,j,k))) + 6.0_db*pxy(i,j,k)*(1.0_db &
				   + 3.0_db*v(i,j,k) + u(i,j,k)*(3.0_db + 7.0_db*v(i,j,k))))*w(i,j,k)**2.0_db) &
				   + pxx(i,j,k)*(1.0_db + 3.0_db*w(i,j,k) + 6.0_db*w(i,j,k)**2.0_db &
				   + 3.0_db*v(i,j,k)*(1.0_db + 3.0_db*w(i,j,k))**2.0_db &
				   + 3.0_db*v(i,j,k)**2.0_db*(2.0_db + 3.0_db*w(i,j,k)*(3.0_db &
				   + 4.0_db*w(i,j,k)))))
#endif
#endif
!19
#ifdef PHASE_CHANGE
                  feq=feq+p3*src(i,j,k)*(invrho_b-invrho_r)
#endif
				  F_discr=(forcez + 3.0_db*forcez*(u(i,j,k) + v(i,j,k)) &
				   + 2.0_db*forcez*w(i,j,k) + forcey*(1.0_db + 3.0_db*u(i,j,k) &
				   + 2.0_db*v(i,j,k) + 3.0_db*w(i,j,k)) + forcex*(1.0_db + 2.0_db*u(i,j,k) &
				   + 3.0_db*v(i,j,k) + 3.0_db*w(i,j,k)))/(72.0_db*rhophi_loc)
                  f(i+1,j+1,k+1,19)=feq + (1.0_db-omega_loc)*fneq1*p3 + 0.5_db*(F_discr)




!!!!!!!!!!!!!!!!!!!!!!!!!!20
#ifdef SECOND_ORDER
			      feq=(rho(i,j,k) + 3.0_db*(u(i,j,k)**2.0_db + (-1.0_db &
			       + v(i,j,k))*v(i,j,k) + (-1.0_db + 3.0_db*v(i,j,k))*w(i,j,k) &
			       + w(i,j,k)**2.0_db + u(i,j,k)*(-1.0_db + 3.0_db*v(i,j,k) &
			       + 3.0_db*w(i,j,k))))/216.0_db
!20
				  fneq1=3.0_db*(pxx(i,j,k) + 3.0_db*pxy(i,j,k) &
				   + 3.0_db*pxz(i,j,k) + pyy(i,j,k) + 3.0_db*pyz(i,j,k) + pzz(i,j,k))
#endif
#ifdef FORTH_ORDER
!20
			      feq=(rho(i,j,k) + 3.0_db*((-1.0_db + v(i,j,k))*v(i,j,k) &
			       + (-1.0_db - 3.0_db*(-1.0_db + v(i,j,k))*v(i,j,k))*w(i,j,k) &
			       + (1.0_db + 3.0_db*(-1.0_db + v(i,j,k))*v(i,j,k))*w(i,j,k)**2.0_db &
			       + u(i,j,k)**2.0_db*(1.0_db + 3.0_db*v(i,j,k)**2.0_db + 3.0_db*(-1.0_db &
			       + w(i,j,k))*w(i,j,k) + v(i,j,k)*(-3.0_db + 9.0_db*w(i,j,k))) &
			       + u(i,j,k)*(-1.0_db - 3.0_db*(-1.0_db + w(i,j,k))*w(i,j,k) &
			       + v(i,j,k)**2.0_db*(-3.0_db + 9.0_db*w(i,j,k)) + v(i,j,k)*(3.0_db &
			       + 9.0_db*(-1.0_db + w(i,j,k))*w(i,j,k)))))/216.0_db
!20
				  fneq1=3.0_db*(pxy(i,j,k)*(3.0_db - 6.0_db*u(i,j,k)) &
				   + pxz(i,j,k)*(3.0_db - 6.0_db*u(i,j,k)) + (pyy(i,j,k) &
				   + 3.0_db*pyz(i,j,k) + pzz(i,j,k))*(1.0_db + 3.0_db*(-1.0_db &
				   + u(i,j,k))*u(i,j,k)) + 3.0_db*v(i,j,k)*((2.0_db*pyz(i,j,k) &
				   + pzz(i,j,k))*(-1.0_db + 3.0_db*u(i,j,k)) + pxy(i,j,k)*(-2.0_db &
				   + 6.0_db*u(i,j,k)) + pzz(i,j,k)*v(i,j,k) + 3.0_db*pxz(i,j,k)*(-1.0_db &
				   + 2.0_db*u(i,j,k) + v(i,j,k))) + 3.0_db*(-pyy(i,j,k) &
				   - 2.0_db*pyz(i,j,k) + 6.0_db*pyy(i,j,k)*u(i,j,k) &
				   + 9.0_db*pyz(i,j,k)*u(i,j,k) + 6.0_db*pyz(i,j,k)*v(i,j,k) &
				   + 3.0_db*pxy(i,j,k)*(-1.0_db + 4.0_db*u(i,j,k) + 4.0_db*v(i,j,k)) &
				   + pxz(i,j,k)*(-2.0_db + 6.0_db*u(i,j,k) + 9.0_db*v(i,j,k)))*w(i,j,k) &
				   + 6.0_db*(3.0_db*pxy(i,j,k) + pyy(i,j,k))*w(i,j,k)**2.0_db &
				   + pxx(i,j,k)*(1.0_db + 6.0_db*v(i,j,k)**2.0_db - 3.0_db*w(i,j,k) &
				   + 6.0_db*w(i,j,k)**2.0_db + 3.0_db*v(i,j,k)*(-1.0_db + 6.0_db*w(i,j,k))))
#endif
#ifdef SIXTH_ORDER
!20
			      feq=(rho(i,j,k) + 3.0_db*((-1.0_db + v(i,j,k))*v(i,j,k) &
			       + (-1.0_db - 3.0_db*(-1.0_db + v(i,j,k))*v(i,j,k))*w(i,j,k) &
			       + (1.0_db + 3.0_db*(-1.0_db + v(i,j,k))*v(i,j,k))*w(i,j,k)**2.0_db &
			       - u(i,j,k)*(1.0_db + 3.0_db*(-1.0_db + v(i,j,k))*v(i,j,k))*(1.0_db &
			       + 3.0_db*(-1.0_db + w(i,j,k))*w(i,j,k)) + u(i,j,k)**2.0_db*(1.0_db &
			       + 3.0_db*(-1.0_db + v(i,j,k))*v(i,j,k))*(1.0_db + 3.0_db*(-1.0_db &
			       + w(i,j,k))*w(i,j,k))))/216.0_db
!20
#ifdef GHOSTONE   
				  fneq1=3.0_db*(pxx(i,j,k) + 3.0_db*pxy(i,j,k) &
				   + 3.0_db*pxz(i,j,k) + pyy(i,j,k) + 3.0_db*pyz(i,j,k) + pzz(i,j,k))    
#else 
				  fneq1=3.0_db*(pxy(i,j,k)*(3.0_db - 6.0_db*u(i,j,k)) &
				   + pxz(i,j,k)*(3.0_db - 6.0_db*u(i,j,k)) + (pyy(i,j,k) &
				   + 3.0_db*pyz(i,j,k) + pzz(i,j,k))*(1.0_db + 3.0_db*(-1.0_db &
				   + u(i,j,k))*u(i,j,k)) - 3.0_db*v(i,j,k)*(pxy(i,j,k)*(2.0_db &
				   - 6.0_db*u(i,j,k)) + 3.0_db*pxz(i,j,k)*(-1.0_db &
				   + 2.0_db*u(i,j,k))*(-1.0_db + v(i,j,k)) + (1.0_db &
				   + 3.0_db*(-1.0_db + u(i,j,k))*u(i,j,k))*(2.0_db*pyz(i,j,k) &
				   + pzz(i,j,k) - pzz(i,j,k)*v(i,j,k))) - 3.0_db*(pyy(i,j,k) &
				   + 2.0_db*pyz(i,j,k) + 3.0_db*(2.0_db*pyy(i,j,k) &
				   + 3.0_db*pyz(i,j,k))*(-1.0_db + u(i,j,k))*u(i,j,k) &
				   - 6.0_db*pyz(i,j,k)*(1.0_db + 3.0_db*(-1.0_db &
				   + u(i,j,k))*u(i,j,k))*v(i,j,k) + 3.0_db*pxy(i,j,k)*(1.0_db &
				   - 4.0_db*v(i,j,k) + 2.0_db*u(i,j,k)*(-2.0_db + 5.0_db*v(i,j,k))) &
				   + pxz(i,j,k)*(2.0_db + 9.0_db*(-1.0_db + v(i,j,k))*v(i,j,k) &
				   - 6.0_db*u(i,j,k)*(1.0_db + 3.0_db*(-1.0_db + v(i,j,k))*v(i,j,k))))*w(i,j,k) &
				   + 3.0_db*(pyy(i,j,k)*(2.0_db + 9.0_db*(-1.0_db + u(i,j,k))*u(i,j,k)) &
				   + 6.0_db*pxy(i,j,k)*(1.0_db - 3.0_db*v(i,j,k) + u(i,j,k)*(-3.0_db &
				   + 7.0_db*v(i,j,k))))*w(i,j,k)**2.0_db + pxx(i,j,k)*(1.0_db &
				   - 3.0_db*v(i,j,k)*(1.0_db - 3.0_db*w(i,j,k))**2.0_db &
				   - 3.0_db*w(i,j,k) + 6.0_db*w(i,j,k)**2.0_db + 3.0_db*v(i,j,k)**2.0_db*(2.0_db &
				   + 3.0_db*w(i,j,k)*(-3.0_db + 4.0_db*w(i,j,k)))))
#endif
#endif
!20
#ifdef PHASE_CHANGE
                  feq=feq+p3*src(i,j,k)*(invrho_b-invrho_r)
#endif
				  F_discr=(forcez*(-1.0_db + 3.0_db*u(i,j,k) + 3.0_db*v(i,j,k) &
				   + 2.0_db*w(i,j,k)) + forcey*(-1.0_db + 3.0_db*u(i,j,k) &
				   + 2.0_db*v(i,j,k) + 3.0_db*w(i,j,k)) + forcex*(-1.0_db &
				   + 2.0_db*u(i,j,k) + 3.0_db*v(i,j,k) + 3.0_db*w(i,j,k)))/(72.0_db*rhophi_loc)
                  f(i-1,j-1,k-1,20)=feq + (1.0_db-omega_loc)*fneq1*p3 + 0.5_db*(F_discr)




!!!!!!!!!!!!!!!!!!!!!!!!!!21
#ifdef SECOND_ORDER
			      feq=(rho(i,j,k) + 3.0_db*(u(i,j,k)**2.0_db + (-1.0_db &
			       + v(i,j,k))*v(i,j,k) + w(i,j,k) - 3.0_db*v(i,j,k)*w(i,j,k) &
			       + w(i,j,k)**2.0_db + u(i,j,k)*(1.0_db - 3.0_db*v(i,j,k) &
			       + 3.0_db*w(i,j,k))))/216.0_db
!21
				  fneq1=3.0_db*(pxx(i,j,k) - 3.0_db*pxy(i,j,k) &
				   + 3.0_db*pxz(i,j,k) + pyy(i,j,k) - 3.0_db*pyz(i,j,k) + pzz(i,j,k))
#endif
#ifdef FORTH_ORDER
!21
			      feq=(rho(i,j,k) + 3.0_db*((-1.0_db + v(i,j,k))*v(i,j,k) &
			       + w(i,j,k) + 3.0_db*(-1.0_db + v(i,j,k))*v(i,j,k)*w(i,j,k) + (1.0_db &
			       + 3.0_db*(-1.0_db + v(i,j,k))*v(i,j,k))*w(i,j,k)**2.0_db &
			       + u(i,j,k)**2.0_db*(1.0_db + 3.0_db*v(i,j,k)**2.0_db &
			       + 3.0_db*w(i,j,k)*(1.0_db + w(i,j,k)) - 3.0_db*v(i,j,k)*(1.0_db &
			       + 3.0_db*w(i,j,k))) + u(i,j,k)*(1.0_db + 3.0_db*w(i,j,k)*(1.0_db &
			       + w(i,j,k)) + v(i,j,k)**2.0_db*(3.0_db + 9.0_db*w(i,j,k)) &
			       - 3.0_db*v(i,j,k)*(1.0_db + 3.0_db*w(i,j,k)*(1.0_db + w(i,j,k))))))/216.0_db
!21
				  fneq1=3.0_db*(-3.0_db*pxy(i,j,k)*(1.0_db + 2.0_db*u(i,j,k)) &
				   + pxz(i,j,k)*(3.0_db + 6.0_db*u(i,j,k)) + (pyy(i,j,k) &
				   - 3.0_db*pyz(i,j,k) + pzz(i,j,k))*(1.0_db + 3.0_db*u(i,j,k)*(1.0_db &
				   + u(i,j,k))) + 3.0_db*v(i,j,k)*((2.0_db*pyz(i,j,k) - pzz(i,j,k))*(1.0_db &
				   + 3.0_db*u(i,j,k)) + pxy(i,j,k)*(2.0_db + 6.0_db*u(i,j,k)) &
				   + pzz(i,j,k)*v(i,j,k) + 3.0_db*pxz(i,j,k)*(-1.0_db - 2.0_db*u(i,j,k) &
				   + v(i,j,k))) + 3.0_db*(pyy(i,j,k) - 2.0_db*pyz(i,j,k) + 6.0_db*pyy(i,j,k)*u(i,j,k) &
				   - 9.0_db*pyz(i,j,k)*u(i,j,k) + pxz(i,j,k)*(2.0_db + 6.0_db*u(i,j,k) &
				   - 9.0_db*v(i,j,k)) - 3.0_db*pxy(i,j,k)*(1.0_db + 4.0_db*u(i,j,k) &
				   - 4.0_db*v(i,j,k)) + 6.0_db*pyz(i,j,k)*v(i,j,k))*w(i,j,k) &
				   + 6.0_db*(-3.0_db*pxy(i,j,k) + pyy(i,j,k))*w(i,j,k)**2.0_db &
				   + pxx(i,j,k)*(1.0_db + 6.0_db*v(i,j,k)**2.0_db + 3.0_db*w(i,j,k) &
				   + 6.0_db*w(i,j,k)**2.0_db - 3.0_db*v(i,j,k)*(1.0_db + 6.0_db*w(i,j,k))))
#endif
#ifdef SIXTH_ORDER
!21
			      feq=(rho(i,j,k) + 3.0_db*((-1.0_db + v(i,j,k))*v(i,j,k) &
			       + w(i,j,k) + 3.0_db*(-1.0_db + v(i,j,k))*v(i,j,k)*w(i,j,k) &
			       + (1.0_db + 3.0_db*(-1.0_db + v(i,j,k))*v(i,j,k))*w(i,j,k)**2.0_db &
			       + u(i,j,k)*(1.0_db + 3.0_db*(-1.0_db + v(i,j,k))*v(i,j,k))*(1.0_db &
			       + 3.0_db*w(i,j,k)*(1.0_db + w(i,j,k))) + u(i,j,k)**2.0_db*(1.0_db &
			       + 3.0_db*(-1.0_db + v(i,j,k))*v(i,j,k))*(1.0_db &
			       + 3.0_db*w(i,j,k)*(1.0_db + w(i,j,k)))))/216.0_db
!21
#ifdef GHOSTONE   
				  fneq1=3.0_db*(pxx(i,j,k) - 3.0_db*pxy(i,j,k) &
				   + 3.0_db*pxz(i,j,k) + pyy(i,j,k) - 3.0_db*pyz(i,j,k) + pzz(i,j,k))
#else 
				  fneq1=3.0_db*(-3.0_db*pxy(i,j,k)*(1.0_db + 2.0_db*u(i,j,k)) &
				   + pxz(i,j,k)*(3.0_db + 6.0_db*u(i,j,k)) + (pyy(i,j,k) &
				   - 3.0_db*pyz(i,j,k) + pzz(i,j,k))*(1.0_db + 3.0_db*u(i,j,k)*(1.0_db &
				   + u(i,j,k))) + pxx(i,j,k)*(1.0_db + 3.0_db*w(i,j,k) &
				   + 6.0_db*w(i,j,k)**2.0_db - 3.0_db*v(i,j,k)*(1.0_db &
				   + 3.0_db*w(i,j,k))**2.0_db + 3.0_db*v(i,j,k)**2.0_db*(2.0_db &
				   + 3.0_db*w(i,j,k)*(3.0_db + 4.0_db*w(i,j,k)))) + 3.0_db*(((1.0_db &
				   + 3.0_db*u(i,j,k)*(1.0_db + u(i,j,k)))*(2.0_db*pyz(i,j,k) &
				   + pzz(i,j,k)*(-1.0_db + v(i,j,k))) + 3.0_db*pxz(i,j,k)*(1.0_db &
				   + 2.0_db*u(i,j,k))*(-1.0_db + v(i,j,k)))*v(i,j,k) + (pyy(i,j,k) &
				   + 6.0_db*pyy(i,j,k)*u(i,j,k)*(1.0_db + u(i,j,k)) + pyz(i,j,k)*(-2.0_db &
				   + 6.0_db*v(i,j,k) + 9.0_db*u(i,j,k)*(1.0_db + u(i,j,k))*(-1.0_db &
				   + 2.0_db*v(i,j,k))) + pxz(i,j,k)*(2.0_db + 9.0_db*(-1.0_db &
				   + v(i,j,k))*v(i,j,k) + 6.0_db*u(i,j,k)*(1.0_db + 3.0_db*(-1.0_db &
				   + v(i,j,k))*v(i,j,k))))*w(i,j,k) + pyy(i,j,k)*(2.0_db + 9.0_db*u(i,j,k)*(1.0_db &
				   + u(i,j,k)))*w(i,j,k)**2.0_db + pxy(i,j,k)*(-3.0_db*w(i,j,k)*(1.0_db &
				   + 2.0_db*w(i,j,k) + u(i,j,k)*(4.0_db + 6.0_db*w(i,j,k))) &
				   + 2.0_db*v(i,j,k)*((1.0_db + 3.0_db*w(i,j,k))**2.0_db + 3.0_db*u(i,j,k)*(1.0_db &
				   + w(i,j,k)*(5.0_db + 7.0_db*w(i,j,k)))))))
#endif
#endif
!21
#ifdef PHASE_CHANGE
                  feq=feq+p3*src(i,j,k)*(invrho_b-invrho_r)
#endif
				  F_discr=(forcez*(1.0_db + 3.0_db*u(i,j,k) - 3.0_db*v(i,j,k) &
				   + 2.0_db*w(i,j,k)) + forcex*(1.0_db + 2.0_db*u(i,j,k) &
				   - 3.0_db*v(i,j,k) + 3.0_db*w(i,j,k)) - forcey*(1.0_db &
				   + 3.0_db*u(i,j,k) - 2.0_db*v(i,j,k) + 3.0_db*w(i,j,k)))/(72.0_db*rhophi_loc)
                  f(i+1,j-1,k+1,21)=feq + (1.0_db-omega_loc)*fneq1*p3 + 0.5_db*(F_discr)



!!!!!!!!!!!!!!!!!!!!!!!!!!22
#ifdef SECOND_ORDER
			      feq=(rho(i,j,k) + 3.0_db*(u(i,j,k)**2.0_db + v(i,j,k) &
			       + v(i,j,k)**2.0_db - 3.0_db*v(i,j,k)*w(i,j,k) + (-1.0_db &
			       + w(i,j,k))*w(i,j,k) + u(i,j,k)*(-1.0_db - 3.0_db*v(i,j,k) &
			       + 3.0_db*w(i,j,k))))/216.0_db
!22
				  fneq1=3.0_db*(pxx(i,j,k) - 3.0_db*pxy(i,j,k) &
				   + 3.0_db*pxz(i,j,k) + pyy(i,j,k) - 3.0_db*pyz(i,j,k) + pzz(i,j,k))
#endif
#ifdef FORTH_ORDER
!22
			      feq=(rho(i,j,k) + 3.0_db*(v(i,j,k) + v(i,j,k)**2.0_db &
			       - w(i,j,k) - 3.0_db*v(i,j,k)*(1.0_db + v(i,j,k))*w(i,j,k) &
			       + (1.0_db + 3.0_db*v(i,j,k)*(1.0_db + v(i,j,k)))*w(i,j,k)**2.0_db &
			       + u(i,j,k)*(-1.0_db - 3.0_db*v(i,j,k)*(1.0_db + v(i,j,k)) &
			       + 3.0_db*w(i,j,k) + 9.0_db*v(i,j,k)*(1.0_db + v(i,j,k))*w(i,j,k) &
			       - 3.0_db*(1.0_db + 3.0_db*v(i,j,k))*w(i,j,k)**2.0_db) &
			       + u(i,j,k)**2.0_db*(1.0_db - 3.0_db*w(i,j,k) + 3.0_db*(v(i,j,k) &
			       + v(i,j,k)**2.0_db - 3.0_db*v(i,j,k)*w(i,j,k) + w(i,j,k)**2.0_db))))/216.0_db
!22
				  fneq1=3.0_db*(pxz(i,j,k)*(3.0_db - 6.0_db*u(i,j,k)) &
				   + pxy(i,j,k)*(-3.0_db + 6.0_db*u(i,j,k)) + (pyy(i,j,k) &
				   - 3.0_db*pyz(i,j,k) + pzz(i,j,k))*(1.0_db + 3.0_db*(-1.0_db &
				   + u(i,j,k))*u(i,j,k)) + 3.0_db*v(i,j,k)*(-2.0_db*pyz(i,j,k) &
				   + pzz(i,j,k) + 6.0_db*pyz(i,j,k)*u(i,j,k) - 3.0_db*pzz(i,j,k)*u(i,j,k) &
				   + pxy(i,j,k)*(-2.0_db + 6.0_db*u(i,j,k)) + pzz(i,j,k)*v(i,j,k) &
				   + 3.0_db*pxz(i,j,k)*(1.0_db - 2.0_db*u(i,j,k) + v(i,j,k))) &
				   - 3.0_db*(pyy(i,j,k) - 2.0_db*pyz(i,j,k) - 6.0_db*pyy(i,j,k)*u(i,j,k) &
				   + 9.0_db*pyz(i,j,k)*u(i,j,k) + 3.0_db*pxy(i,j,k)*(-1.0_db &
				   + 4.0_db*u(i,j,k) - 4.0_db*v(i,j,k)) - 6.0_db*pyz(i,j,k)*v(i,j,k) &
				   + pxz(i,j,k)*(2.0_db - 6.0_db*u(i,j,k) + 9.0_db*v(i,j,k)))*w(i,j,k) &
				   + 6.0_db*(-3.0_db*pxy(i,j,k) + pyy(i,j,k))*w(i,j,k)**2.0_db &
				   + pxx(i,j,k)*(1.0_db + 3.0_db*v(i,j,k) + 6.0_db*v(i,j,k)**2.0_db &
				   - 3.0_db*(1.0_db + 6.0_db*v(i,j,k))*w(i,j,k) + 6.0_db*w(i,j,k)**2.0_db))
#endif
#ifdef SIXTH_ORDER
!22
			      feq=(rho(i,j,k) + 3.0_db*(v(i,j,k) + v(i,j,k)**2.0_db &
			       - w(i,j,k) - 3.0_db*v(i,j,k)*(1.0_db + v(i,j,k))*w(i,j,k) + (1.0_db &
			       + 3.0_db*v(i,j,k)*(1.0_db + v(i,j,k)))*w(i,j,k)**2.0_db - u(i,j,k)*(1.0_db &
			       + 3.0_db*v(i,j,k)*(1.0_db + v(i,j,k)))*(1.0_db + 3.0_db*(-1.0_db &
			       + w(i,j,k))*w(i,j,k)) + u(i,j,k)**2.0_db*(1.0_db + 3.0_db*v(i,j,k)*(1.0_db &
			       + v(i,j,k)))*(1.0_db + 3.0_db*(-1.0_db + w(i,j,k))*w(i,j,k))))/216.0_db
!22
#ifdef GHOSTONE   
				  fneq1=3.0_db*(pxx(i,j,k) - 3.0_db*pxy(i,j,k) &
				   + 3.0_db*pxz(i,j,k) + pyy(i,j,k) - 3.0_db*pyz(i,j,k) + pzz(i,j,k))
#else 
				  fneq1=3.0_db*(pxz(i,j,k)*(3.0_db - 6.0_db*u(i,j,k)) &
				   + pxy(i,j,k)*(-3.0_db + 6.0_db*u(i,j,k)) + (pyy(i,j,k) &
				   - 3.0_db*pyz(i,j,k) + pzz(i,j,k))*(1.0_db + 3.0_db*(-1.0_db &
				   + u(i,j,k))*u(i,j,k)) + 3.0_db*v(i,j,k)*(pxy(i,j,k)*(-2.0_db &
				   + 6.0_db*u(i,j,k)) - 3.0_db*pxz(i,j,k)*(-1.0_db + 2.0_db*u(i,j,k))*(1.0_db &
				   + v(i,j,k)) - (1.0_db + 3.0_db*(-1.0_db + u(i,j,k))*u(i,j,k))*(2.0_db*pyz(i,j,k) &
				   - pzz(i,j,k)*(1.0_db + v(i,j,k)))) - 3.0_db*(pyy(i,j,k) &
				   - 2.0_db*pyz(i,j,k) + 3.0_db*(2.0_db*pyy(i,j,k) - 3.0_db*pyz(i,j,k))*(-1.0_db &
				   + u(i,j,k))*u(i,j,k) - 6.0_db*pyz(i,j,k)*(1.0_db + 3.0_db*(-1.0_db &
				   + u(i,j,k))*u(i,j,k))*v(i,j,k) + 3.0_db*pxy(i,j,k)*(-1.0_db &
				   - 4.0_db*v(i,j,k) + 2.0_db*u(i,j,k)*(2.0_db + 5.0_db*v(i,j,k))) &
				   + pxz(i,j,k)*(2.0_db + 9.0_db*v(i,j,k)*(1.0_db + v(i,j,k)) &
				   - 6.0_db*u(i,j,k)*(1.0_db + 3.0_db*v(i,j,k)*(1.0_db + v(i,j,k)))))*w(i,j,k) &
				   + 3.0_db*(pyy(i,j,k)*(2.0_db + 9.0_db*(-1.0_db + u(i,j,k))*u(i,j,k)) &
				   + 6.0_db*pxy(i,j,k)*(-1.0_db - 3.0_db*v(i,j,k) + u(i,j,k)*(3.0_db &
				   + 7.0_db*v(i,j,k))))*w(i,j,k)**2.0_db + pxx(i,j,k)*(1.0_db &
				   + 3.0_db*v(i,j,k) + 6.0_db*v(i,j,k)**2.0_db - 3.0_db*(1.0_db &
				   + 3.0_db*v(i,j,k))**2.0_db*w(i,j,k) + 3.0_db*(2.0_db &
				   + 3.0_db*v(i,j,k)*(3.0_db + 4.0_db*v(i,j,k)))*w(i,j,k)**2.0_db))
#endif
#endif
!22
#ifdef PHASE_CHANGE
                  feq=feq+p3*src(i,j,k)*(invrho_b-invrho_r)
#endif
				  F_discr=(forcey*(1.0_db - 3.0_db*u(i,j,k) + 2.0_db*v(i,j,k) &
				   - 3.0_db*w(i,j,k)) + forcez*(-1.0_db + 3.0_db*u(i,j,k) &
				   - 3.0_db*v(i,j,k) + 2.0_db*w(i,j,k)) + forcex*(-1.0_db &
				   + 2.0_db*u(i,j,k) - 3.0_db*v(i,j,k) + 3.0_db*w(i,j,k)))/(72.0_db*rhophi_loc)
                  f(i-1,j+1,k-1,22)=feq + (1.0_db-omega_loc)*fneq1*p3 + 0.5_db*(F_discr)




!!!!!!!!!!!!!!!!!!!!!!!!!!23
#ifdef SECOND_ORDER
			      feq=(rho(i,j,k) + 3.0_db*(-u(i,j,k) + u(i,j,k)**2.0_db &
			       - v(i,j,k) + 3.0_db*u(i,j,k)*v(i,j,k) + v(i,j,k)**2.0_db + w(i,j,k) &
			       - 3.0_db*(u(i,j,k) + v(i,j,k))*w(i,j,k) + w(i,j,k)**2.0_db))/216.0_db
!23
				  fneq1=3.0_db*(pxx(i,j,k) + 3.0_db*pxy(i,j,k) &
				   - 3.0_db*pxz(i,j,k) + pyy(i,j,k) - 3.0_db*pyz(i,j,k) + pzz(i,j,k))
#endif
#ifdef FORTH_ORDER
!23
			      feq=(rho(i,j,k) + 3.0_db*((-1.0_db + v(i,j,k))*v(i,j,k) &
			       + w(i,j,k) + 3.0_db*(-1.0_db + v(i,j,k))*v(i,j,k)*w(i,j,k) + (1.0_db &
			       + 3.0_db*(-1.0_db + v(i,j,k))*v(i,j,k))*w(i,j,k)**2.0_db &
			       + u(i,j,k)**2.0_db*(1.0_db + 3.0_db*v(i,j,k)**2.0_db &
			       + 3.0_db*w(i,j,k)*(1.0_db + w(i,j,k)) - 3.0_db*v(i,j,k)*(1.0_db &
			       + 3.0_db*w(i,j,k))) - u(i,j,k)*(1.0_db + 3.0_db*w(i,j,k)*(1.0_db &
			       + w(i,j,k)) + v(i,j,k)**2.0_db*(3.0_db + 9.0_db*w(i,j,k)) &
			       - 3.0_db*v(i,j,k)*(1.0_db + 3.0_db*w(i,j,k)*(1.0_db + w(i,j,k))))))/216.0_db
!23
				  fneq1=3.0_db*(pxy(i,j,k)*(3.0_db - 6.0_db*u(i,j,k)) &
				   + pxz(i,j,k)*(-3.0_db + 6.0_db*u(i,j,k)) + (pyy(i,j,k) &
				   - 3.0_db*pyz(i,j,k) + pzz(i,j,k))*(1.0_db + 3.0_db*(-1.0_db &
				   + u(i,j,k))*u(i,j,k)) - 3.0_db*v(i,j,k)*(-2.0_db*pyz(i,j,k) &
				   + pzz(i,j,k) + pxy(i,j,k)*(2.0_db - 6.0_db*u(i,j,k)) &
				   + 6.0_db*pyz(i,j,k)*u(i,j,k) - 3.0_db*pzz(i,j,k)*u(i,j,k) &
				   - pzz(i,j,k)*v(i,j,k) + 3.0_db*pxz(i,j,k)*(-1.0_db + 2.0_db*u(i,j,k) &
				   + v(i,j,k))) + 3.0_db*(pyy(i,j,k) - 2.0_db*pyz(i,j,k) &
				   - 6.0_db*pyy(i,j,k)*u(i,j,k) + 9.0_db*pyz(i,j,k)*u(i,j,k) &
				   + 6.0_db*pyz(i,j,k)*v(i,j,k) - 3.0_db*pxy(i,j,k)*(-1.0_db &
				   + 4.0_db*u(i,j,k) + 4.0_db*v(i,j,k)) + pxz(i,j,k)*(-2.0_db &
				   + 6.0_db*u(i,j,k) + 9.0_db*v(i,j,k)))*w(i,j,k) + 6.0_db*(3.0_db*pxy(i,j,k) &
				   + pyy(i,j,k))*w(i,j,k)**2.0_db + pxx(i,j,k)*(1.0_db + 6.0_db*v(i,j,k)**2.0_db &
				   + 3.0_db*w(i,j,k) + 6.0_db*w(i,j,k)**2.0_db - 3.0_db*v(i,j,k)*(1.0_db &
				   + 6.0_db*w(i,j,k))))
#endif
#ifdef SIXTH_ORDER
!23
			      feq=(rho(i,j,k) + 3.0_db*((-1.0_db + v(i,j,k))*v(i,j,k) &
			       + w(i,j,k) + 3.0_db*(-1.0_db + v(i,j,k))*v(i,j,k)*w(i,j,k) + (1.0_db &
			       + 3.0_db*(-1.0_db + v(i,j,k))*v(i,j,k))*w(i,j,k)**2.0_db - u(i,j,k)*(1.0_db &
			       + 3.0_db*(-1.0_db + v(i,j,k))*v(i,j,k))*(1.0_db + 3.0_db*w(i,j,k)*(1.0_db &
			       + w(i,j,k))) + u(i,j,k)**2.0_db*(1.0_db + 3.0_db*(-1.0_db + v(i,j,k))*v(i,j,k))*(1.0_db &
			       + 3.0_db*w(i,j,k)*(1.0_db + w(i,j,k)))))/216.0_db
!23
#ifdef GHOSTONE   
				  fneq1=3.0_db*(pxx(i,j,k) + 3.0_db*pxy(i,j,k) &
				   - 3.0_db*pxz(i,j,k) + pyy(i,j,k) - 3.0_db*pyz(i,j,k) + pzz(i,j,k))
#else 
				  fneq1=3.0_db*(pxy(i,j,k)*(3.0_db - 6.0_db*u(i,j,k)) &
				   + pxz(i,j,k)*(-3.0_db + 6.0_db*u(i,j,k)) + (pyy(i,j,k) &
				   - 3.0_db*pyz(i,j,k) + pzz(i,j,k))*(1.0_db + 3.0_db*(-1.0_db &
				   + u(i,j,k))*u(i,j,k)) + pxx(i,j,k)*(1.0_db + 3.0_db*w(i,j,k) &
				   + 6.0_db*w(i,j,k)**2.0_db - 3.0_db*v(i,j,k)*(1.0_db &
				   + 3.0_db*w(i,j,k))**2.0_db + 3.0_db*v(i,j,k)**2.0_db*(2.0_db &
				   + 3.0_db*w(i,j,k)*(3.0_db + 4.0_db*w(i,j,k)))) + 3.0_db*(((1.0_db &
				   + 3.0_db*(-1.0_db + u(i,j,k))*u(i,j,k))*(2.0_db*pyz(i,j,k) &
				   + pzz(i,j,k)*(-1.0_db + v(i,j,k))) + 3.0_db*pxz(i,j,k)*(-1.0_db &
				   + 2.0_db*u(i,j,k))*(-1.0_db + v(i,j,k)))*v(i,j,k) + (pyy(i,j,k) &
				   + 6.0_db*pyy(i,j,k)*(-1.0_db + u(i,j,k))*u(i,j,k) + pyz(i,j,k)*(-2.0_db &
				   + 6.0_db*v(i,j,k) + 9.0_db*(-1.0_db + u(i,j,k))*u(i,j,k)*(-1.0_db &
				   + 2.0_db*v(i,j,k))) + pxz(i,j,k)*(-2.0_db - 9.0_db*(-1.0_db &
				   + v(i,j,k))*v(i,j,k) + 6.0_db*u(i,j,k)*(1.0_db + 3.0_db*(-1.0_db &
				   + v(i,j,k))*v(i,j,k))))*w(i,j,k) + pyy(i,j,k)*(2.0_db + 9.0_db*(-1.0_db &
				   + u(i,j,k))*u(i,j,k))*w(i,j,k)**2.0_db + pxy(i,j,k)*(-2.0_db*v(i,j,k)*(1.0_db &
				   + 3.0_db*w(i,j,k))**2.0_db + 3.0_db*w(i,j,k)*(1.0_db + 2.0_db*w(i,j,k) &
				   - 2.0_db*u(i,j,k)*(2.0_db + 3.0_db*w(i,j,k))) + 6.0_db*u(i,j,k)*v(i,j,k)*(1.0_db &
				   + w(i,j,k)*(5.0_db + 7.0_db*w(i,j,k))))))
#endif
#endif
!23
#ifdef PHASE_CHANGE
                  feq=feq+p3*src(i,j,k)*(invrho_b-invrho_r)
#endif
				  F_discr=(forcez - 3.0_db*forcez*(u(i,j,k) + v(i,j,k)) + forcey*(-1.0_db &
				   + 3.0_db*u(i,j,k) + 2.0_db*v(i,j,k) - 3.0_db*w(i,j,k)) &
				   + forcex*(-1.0_db + 2.0_db*u(i,j,k) + 3.0_db*v(i,j,k) &
				   - 3.0_db*w(i,j,k)) + 2.0_db*forcez*w(i,j,k))/(72.0_db*rhophi_loc)
                  f(i-1,j-1,k+1,23)=feq + (1.0_db-omega_loc)*fneq1*p3 + 0.5_db*(F_discr)


!!!!!!!!!!!!!!!!!!!!!!!!!!24
#ifdef SECOND_ORDER
			      feq=(rho(i,j,k) + 3.0_db*(u(i,j,k) + u(i,j,k)**2.0_db &
			       + v(i,j,k) + 3.0_db*u(i,j,k)*v(i,j,k) + v(i,j,k)**2.0_db - w(i,j,k) &
			       - 3.0_db*(u(i,j,k) + v(i,j,k))*w(i,j,k) + w(i,j,k)**2.0_db))/216.0_db
!24
				  fneq1=3.0_db*(pxx(i,j,k) + 3.0_db*pxy(i,j,k) &
				   - 3.0_db*pxz(i,j,k) + pyy(i,j,k) - 3.0_db*pyz(i,j,k) + pzz(i,j,k))
#endif
#ifdef FORTH_ORDER
!24
			      feq=(rho(i,j,k) + 3.0_db*(v(i,j,k) + v(i,j,k)**2.0_db &
			       - w(i,j,k) - 3.0_db*v(i,j,k)*(1.0_db + v(i,j,k))*w(i,j,k) &
			       + (1.0_db + 3.0_db*v(i,j,k)*(1.0_db + v(i,j,k)))*w(i,j,k)**2.0_db &
			       + u(i,j,k)*(1.0_db + v(i,j,k)**2.0_db*(3.0_db - 9.0_db*w(i,j,k)) &
			       + 3.0_db*(-1.0_db + w(i,j,k))*w(i,j,k) + v(i,j,k)*(3.0_db &
			       + 9.0_db*(-1.0_db + w(i,j,k))*w(i,j,k))) + u(i,j,k)**2.0_db*(1.0_db &
			       - 3.0_db*w(i,j,k) + 3.0_db*(v(i,j,k) + v(i,j,k)**2.0_db &
			       - 3.0_db*v(i,j,k)*w(i,j,k) + w(i,j,k)**2.0_db))))/216.0_db
!24
				  fneq1=3.0_db*(-3.0_db*pxz(i,j,k)*(1.0_db + 2.0_db*u(i,j,k)) &
				   + pxy(i,j,k)*(3.0_db + 6.0_db*u(i,j,k)) + (pyy(i,j,k) &
				   - 3.0_db*pyz(i,j,k) + pzz(i,j,k))*(1.0_db + 3.0_db*u(i,j,k)*(1.0_db &
				   + u(i,j,k))) + 3.0_db*v(i,j,k)*(-2.0_db*pyz(i,j,k) + pzz(i,j,k) &
				   - 6.0_db*pyz(i,j,k)*u(i,j,k) + 3.0_db*pzz(i,j,k)*u(i,j,k) &
				   + pxy(i,j,k)*(2.0_db + 6.0_db*u(i,j,k)) + pzz(i,j,k)*v(i,j,k) &
				   - 3.0_db*pxz(i,j,k)*(1.0_db + 2.0_db*u(i,j,k) + v(i,j,k))) &
				   - 3.0_db*(pyy(i,j,k) - 2.0_db*pyz(i,j,k) + 6.0_db*pyy(i,j,k)*u(i,j,k) &
				   - 9.0_db*pyz(i,j,k)*u(i,j,k) - 6.0_db*pyz(i,j,k)*v(i,j,k) &
				   + 3.0_db*pxy(i,j,k)*(1.0_db + 4.0_db*u(i,j,k) + 4.0_db*v(i,j,k)) &
				   - pxz(i,j,k)*(2.0_db + 6.0_db*u(i,j,k) + 9.0_db*v(i,j,k)))*w(i,j,k) &
				   + 6.0_db*(3.0_db*pxy(i,j,k) + pyy(i,j,k))*w(i,j,k)**2.0_db &
				   + pxx(i,j,k)*(1.0_db + 3.0_db*v(i,j,k) + 6.0_db*v(i,j,k)**2.0_db &
				   - 3.0_db*(1.0_db + 6.0_db*v(i,j,k))*w(i,j,k) + 6.0_db*w(i,j,k)**2.0_db))
#endif
#ifdef SIXTH_ORDER
!24
			      feq=(rho(i,j,k) + 3.0_db*(v(i,j,k) + v(i,j,k)**2.0_db &
			       - w(i,j,k) - 3.0_db*v(i,j,k)*(1.0_db + v(i,j,k))*w(i,j,k) + (1.0_db &
			       + 3.0_db*v(i,j,k)*(1.0_db + v(i,j,k)))*w(i,j,k)**2.0_db + u(i,j,k)*(1.0_db &
			       + 3.0_db*v(i,j,k)*(1.0_db + v(i,j,k)))*(1.0_db + 3.0_db*(-1.0_db &
			       + w(i,j,k))*w(i,j,k)) + u(i,j,k)**2.0_db*(1.0_db + 3.0_db*v(i,j,k)*(1.0_db &
			       + v(i,j,k)))*(1.0_db + 3.0_db*(-1.0_db + w(i,j,k))*w(i,j,k))))/216.0_db
!24
#ifdef GHOSTONE   
				  fneq1=3.0_db*(pxx(i,j,k) + 3.0_db*pxy(i,j,k) &
				   - 3.0_db*pxz(i,j,k) + pyy(i,j,k) - 3.0_db*pyz(i,j,k) + pzz(i,j,k))
#else   
				  fneq1=3.0_db*(-3.0_db*pxz(i,j,k)*(1.0_db + 2.0_db*u(i,j,k)) &
				   + pxy(i,j,k)*(3.0_db + 6.0_db*u(i,j,k)) + (pyy(i,j,k) - 3.0_db*pyz(i,j,k) &
				   + pzz(i,j,k))*(1.0_db + 3.0_db*u(i,j,k)*(1.0_db + u(i,j,k))) &
				   + 3.0_db*v(i,j,k)*(pxy(i,j,k)*(2.0_db + 6.0_db*u(i,j,k)) &
				   - 3.0_db*pxz(i,j,k)*(1.0_db + 2.0_db*u(i,j,k))*(1.0_db + v(i,j,k)) &
				   - (1.0_db + 3.0_db*u(i,j,k)*(1.0_db + u(i,j,k)))*(2.0_db*pyz(i,j,k) &
				   - pzz(i,j,k)*(1.0_db + v(i,j,k)))) - 3.0_db*(pyy(i,j,k) &
				   - 2.0_db*pyz(i,j,k) + 3.0_db*(2.0_db*pyy(i,j,k) - 3.0_db*pyz(i,j,k))*u(i,j,k)*(1.0_db &
				   + u(i,j,k)) - 6.0_db*pyz(i,j,k)*(1.0_db + 3.0_db*u(i,j,k)*(1.0_db &
				   + u(i,j,k)))*v(i,j,k) + 3.0_db*pxy(i,j,k)*(1.0_db + 4.0_db*v(i,j,k) &
				   + 2.0_db*u(i,j,k)*(2.0_db + 5.0_db*v(i,j,k))) - pxz(i,j,k)*(2.0_db &
				   + 9.0_db*v(i,j,k)*(1.0_db + v(i,j,k)) + 6.0_db*u(i,j,k)*(1.0_db &
				   + 3.0_db*v(i,j,k)*(1.0_db + v(i,j,k)))))*w(i,j,k) + 3.0_db*(pyy(i,j,k)*(2.0_db &
				   + 9.0_db*u(i,j,k)*(1.0_db + u(i,j,k))) + 6.0_db*pxy(i,j,k)*(1.0_db &
				   + 3.0_db*v(i,j,k) + u(i,j,k)*(3.0_db + 7.0_db*v(i,j,k))))*w(i,j,k)**2.0_db &
				   + pxx(i,j,k)*(1.0_db + 3.0_db*v(i,j,k) + 6.0_db*v(i,j,k)**2.0_db - 3.0_db*(1.0_db &
				   + 3.0_db*v(i,j,k))**2.0_db*w(i,j,k) + 3.0_db*(2.0_db + 3.0_db*v(i,j,k)*(3.0_db &
				   + 4.0_db*v(i,j,k)))*w(i,j,k)**2.0_db))
#endif
#endif
!24
#ifdef PHASE_CHANGE
                  feq=feq+p3*src(i,j,k)*(invrho_b-invrho_r)
#endif
				  F_discr=(-forcez - 3.0_db*forcez*(u(i,j,k) + v(i,j,k)) + forcey*(1.0_db &
				   + 3.0_db*u(i,j,k) + 2.0_db*v(i,j,k) - 3.0_db*w(i,j,k)) + forcex*(1.0_db &
				   + 2.0_db*u(i,j,k) + 3.0_db*v(i,j,k) - 3.0_db*w(i,j,k)) &
				   + 2.0_db*forcez*w(i,j,k))/(72.0_db*rhophi_loc)
                  f(i+1,j+1,k-1,24)=feq + (1.0_db-omega_loc)*fneq1*p3 + 0.5_db*(F_discr)



!!!!!!!!!!!!!!!!!!!!!!!!!!25
#ifdef SECOND_ORDER
			      feq=(rho(i,j,k) + 3.0_db*(u(i,j,k)**2.0_db + (-1.0_db &
			       + v(i,j,k))*v(i,j,k) + u(i,j,k)*(1.0_db - 3.0_db*v(i,j,k) &
			       - 3.0_db*w(i,j,k)) + (-1.0_db + 3.0_db*v(i,j,k))*w(i,j,k) &
			       + w(i,j,k)**2.0_db))/216.0_db
!25
				  fneq1=3.0_db*(pxx(i,j,k) - 3.0_db*pxy(i,j,k) - 3.0_db*pxz(i,j,k) &
				   + pyy(i,j,k) + 3.0_db*pyz(i,j,k) + pzz(i,j,k))
#endif
#ifdef FORTH_ORDER
!25
			      feq=(rho(i,j,k) + 3.0_db*((-1.0_db + v(i,j,k))*v(i,j,k) &
			       + (-1.0_db - 3.0_db*(-1.0_db + v(i,j,k))*v(i,j,k))*w(i,j,k) + (1.0_db &
			       + 3.0_db*(-1.0_db + v(i,j,k))*v(i,j,k))*w(i,j,k)**2.0_db &
			       + u(i,j,k)**2.0_db*(1.0_db + 3.0_db*v(i,j,k)**2.0_db + 3.0_db*(-1.0_db &
			       + w(i,j,k))*w(i,j,k) + v(i,j,k)*(-3.0_db + 9.0_db*w(i,j,k))) &
			       + u(i,j,k)*(1.0_db + v(i,j,k)**2.0_db*(3.0_db - 9.0_db*w(i,j,k)) &
			       + 3.0_db*(-1.0_db + w(i,j,k))*w(i,j,k) + v(i,j,k)*(-3.0_db - 9.0_db*(-1.0_db &
			       + w(i,j,k))*w(i,j,k)))))/216.0_db
!25
				  fneq1=3.0_db*(-3.0_db*pxy(i,j,k)*(1.0_db + 2.0_db*u(i,j,k)) &
				   - 3.0_db*pxz(i,j,k)*(1.0_db + 2.0_db*u(i,j,k)) + (pyy(i,j,k) &
				   + 3.0_db*pyz(i,j,k) + pzz(i,j,k))*(1.0_db + 3.0_db*u(i,j,k)*(1.0_db &
				   + u(i,j,k))) - 3.0_db*v(i,j,k)*(-2.0_db*pxy(i,j,k)*(1.0_db + 3.0_db*u(i,j,k)) &
				   + (2.0_db*pyz(i,j,k) + pzz(i,j,k))*(1.0_db + 3.0_db*u(i,j,k)) &
				   - pzz(i,j,k)*v(i,j,k) + 3.0_db*pxz(i,j,k)*(-1.0_db - 2.0_db*u(i,j,k) &
				   + v(i,j,k))) - 3.0_db*(pyy(i,j,k) + 2.0_db*pyz(i,j,k) + 6.0_db*pyy(i,j,k)*u(i,j,k) &
				   + 9.0_db*pyz(i,j,k)*u(i,j,k) - 3.0_db*pxy(i,j,k)*(1.0_db + 4.0_db*u(i,j,k) &
				   - 4.0_db*v(i,j,k)) - 6.0_db*pyz(i,j,k)*v(i,j,k) + pxz(i,j,k)*(-2.0_db &
				   - 6.0_db*u(i,j,k) + 9.0_db*v(i,j,k)))*w(i,j,k) + 6.0_db*(-3.0_db*pxy(i,j,k) &
				   + pyy(i,j,k))*w(i,j,k)**2.0_db + pxx(i,j,k)*(1.0_db + 6.0_db*v(i,j,k)**2.0_db &
				   - 3.0_db*w(i,j,k) + 6.0_db*w(i,j,k)**2.0_db + 3.0_db*v(i,j,k)*(-1.0_db &
				   + 6.0_db*w(i,j,k))))
#endif
#ifdef SIXTH_ORDER
!25
			      feq=(rho(i,j,k) + 3.0_db*((-1.0_db + v(i,j,k))*v(i,j,k) &
			       + (-1.0_db - 3.0_db*(-1.0_db + v(i,j,k))*v(i,j,k))*w(i,j,k) + (1.0_db &
			       + 3.0_db*(-1.0_db + v(i,j,k))*v(i,j,k))*w(i,j,k)**2.0_db + u(i,j,k)*(1.0_db &
			       + 3.0_db*(-1.0_db + v(i,j,k))*v(i,j,k))*(1.0_db + 3.0_db*(-1.0_db &
			       + w(i,j,k))*w(i,j,k)) + u(i,j,k)**2.0_db*(1.0_db + 3.0_db*(-1.0_db &
			       + v(i,j,k))*v(i,j,k))*(1.0_db + 3.0_db*(-1.0_db + w(i,j,k))*w(i,j,k))))/216.0_db
!25
#ifdef GHOSTONE   
				  fneq1=3.0_db*(pxx(i,j,k) - 3.0_db*pxy(i,j,k) - 3.0_db*pxz(i,j,k) &
				   + pyy(i,j,k) + 3.0_db*pyz(i,j,k) + pzz(i,j,k))
#else 
				  fneq1=3.0_db*(-3.0_db*pxy(i,j,k)*(1.0_db + 2.0_db*u(i,j,k)) - 3.0_db*pxz(i,j,k)*(1.0_db &
				   + 2.0_db*u(i,j,k)) + (pyy(i,j,k) + 3.0_db*pyz(i,j,k) + pzz(i,j,k))*(1.0_db &
				   + 3.0_db*u(i,j,k)*(1.0_db + u(i,j,k))) - 3.0_db*v(i,j,k)*(-2.0_db*pxy(i,j,k)*(1.0_db &
				   + 3.0_db*u(i,j,k)) + 3.0_db*pxz(i,j,k)*(1.0_db + 2.0_db*u(i,j,k))*(-1.0_db &
				   + v(i,j,k)) + (1.0_db + 3.0_db*u(i,j,k)*(1.0_db + u(i,j,k)))*(2.0_db*pyz(i,j,k) &
				   + pzz(i,j,k) - pzz(i,j,k)*v(i,j,k))) - 3.0_db*(pyy(i,j,k) + 2.0_db*pyz(i,j,k) &
				   + 3.0_db*(2.0_db*pyy(i,j,k) + 3.0_db*pyz(i,j,k))*u(i,j,k)*(1.0_db + u(i,j,k)) &
				   - 6.0_db*pyz(i,j,k)*(1.0_db + 3.0_db*u(i,j,k)*(1.0_db + u(i,j,k)))*v(i,j,k) &
				   + 3.0_db*pxy(i,j,k)*(-1.0_db + 4.0_db*v(i,j,k) + 2.0_db*u(i,j,k)*(-2.0_db &
				   + 5.0_db*v(i,j,k))) - pxz(i,j,k)*(2.0_db + 9.0_db*(-1.0_db + v(i,j,k))*v(i,j,k) &
				   + 6.0_db*u(i,j,k)*(1.0_db + 3.0_db*(-1.0_db + v(i,j,k))*v(i,j,k))))*w(i,j,k) &
				   + 3.0_db*(pyy(i,j,k)*(2.0_db + 9.0_db*u(i,j,k)*(1.0_db + u(i,j,k))) &
				   + 6.0_db*pxy(i,j,k)*(-1.0_db + 3.0_db*v(i,j,k) + u(i,j,k)*(-3.0_db &
				   + 7.0_db*v(i,j,k))))*w(i,j,k)**2.0_db + pxx(i,j,k)*(1.0_db &
				   - 3.0_db*v(i,j,k)*(1.0_db - 3.0_db*w(i,j,k))**2.0_db - 3.0_db*w(i,j,k) &
				   + 6.0_db*w(i,j,k)**2.0_db + 3.0_db*v(i,j,k)**2.0_db*(2.0_db &
				   + 3.0_db*w(i,j,k)*(-3.0_db + 4.0_db*w(i,j,k)))))
#endif
#endif
!25
#ifdef PHASE_CHANGE
                  feq=feq+p3*src(i,j,k)*(invrho_b-invrho_r)
#endif
				  F_discr=(forcex*(1.0_db + 2.0_db*u(i,j,k) - 3.0_db*v(i,j,k) &
				   - 3.0_db*w(i,j,k)) + forcez*(-1.0_db - 3.0_db*u(i,j,k) &
				   + 3.0_db*v(i,j,k) + 2.0_db*w(i,j,k)) + forcey*(-1.0_db &
				   - 3.0_db*u(i,j,k) + 2.0_db*v(i,j,k) + 3.0_db*w(i,j,k)))/(72.0_db*rhophi_loc)
                  f(i+1,j-1,k-1,25)=feq + (1.0_db-omega_loc)*fneq1*p3 + 0.5_db*(F_discr)



!!!!!!!!!!!!!!!!!!!!!!!!!!26
#ifdef SECOND_ORDER
			      feq=(rho(i,j,k) + 3.0_db*(u(i,j,k)**2.0_db + v(i,j,k) &
			       + v(i,j,k)**2.0_db + w(i,j,k) + 3.0_db*v(i,j,k)*w(i,j,k) &
			       + w(i,j,k)**2.0_db - u(i,j,k)*(1.0_db + 3.0_db*v(i,j,k) &
			       + 3.0_db*w(i,j,k))))/216.0_db
!26
				  fneq1=3.0_db*(pxx(i,j,k) - 3.0_db*pxy(i,j,k) &
				   - 3.0_db*pxz(i,j,k) + pyy(i,j,k) + 3.0_db*pyz(i,j,k) + pzz(i,j,k))
#endif
#ifdef FORTH_ORDER
!26
			      feq=(rho(i,j,k) + 3.0_db*(v(i,j,k) + v(i,j,k)**2.0_db &
			       + w(i,j,k) + 3.0_db*v(i,j,k)*(1.0_db + v(i,j,k))*w(i,j,k) + (1.0_db &
			       + 3.0_db*v(i,j,k)*(1.0_db + v(i,j,k)))*w(i,j,k)**2.0_db &
			       + u(i,j,k)**2.0_db*(1.0_db + 3.0_db*w(i,j,k) + 3.0_db*(v(i,j,k) &
			       + v(i,j,k)**2.0_db + 3.0_db*v(i,j,k)*w(i,j,k) + w(i,j,k)**2.0_db)) &
			       - u(i,j,k)*(1.0_db + 3.0_db*w(i,j,k)*(1.0_db + w(i,j,k)) &
			       + v(i,j,k)**2.0_db*(3.0_db + 9.0_db*w(i,j,k)) + v(i,j,k)*(3.0_db &
			       + 9.0_db*w(i,j,k)*(1.0_db + w(i,j,k))))))/216.0_db
!26
				  fneq1=3.0_db*(pxy(i,j,k)*(-3.0_db + 6.0_db*u(i,j,k)) &
				   + pxz(i,j,k)*(-3.0_db + 6.0_db*u(i,j,k)) + (pyy(i,j,k) &
				   + 3.0_db*pyz(i,j,k) + pzz(i,j,k))*(1.0_db + 3.0_db*(-1.0_db&
				   + u(i,j,k))*u(i,j,k)) + 3.0_db*v(i,j,k)*(2.0_db*pyz(i,j,k) &
				   + pzz(i,j,k) + 6.0_db*pxz(i,j,k)*u(i,j,k) - 3.0_db*(2.0_db*pyz(i,j,k) &
				   + pzz(i,j,k))*u(i,j,k) + pxy(i,j,k)*(-2.0_db + 6.0_db*u(i,j,k)) &
				   + pzz(i,j,k)*v(i,j,k) - 3.0_db*pxz(i,j,k)*(1.0_db + v(i,j,k))) &
				   + 3.0_db*(pyy(i,j,k) + 2.0_db*pyz(i,j,k) - 6.0_db*pyy(i,j,k)*u(i,j,k) &
				   - 9.0_db*pyz(i,j,k)*u(i,j,k) + pxz(i,j,k)*(-2.0_db + 6.0_db*u(i,j,k) &
				   - 9.0_db*v(i,j,k)) + 3.0_db*pxy(i,j,k)*(-1.0_db + 4.0_db*u(i,j,k) &
				   - 4.0_db*v(i,j,k)) + 6.0_db*pyz(i,j,k)*v(i,j,k))*w(i,j,k) &
				   + 6.0_db*(-3.0_db*pxy(i,j,k) + pyy(i,j,k))*w(i,j,k)**2.0_db &
				   + pxx(i,j,k)*(1.0_db + 6.0_db*v(i,j,k)**2.0_db + 3.0_db*w(i,j,k) &
				   + 6.0_db*w(i,j,k)**2.0_db + 3.0_db*v(i,j,k)*(1.0_db + 6.0_db*w(i,j,k))))
#endif
#ifdef SIXTH_ORDER
!26
			      feq=(rho(i,j,k) + 3.0_db*(v(i,j,k) + v(i,j,k)**2.0_db &
			       + w(i,j,k) + 3.0_db*v(i,j,k)*(1.0_db + v(i,j,k))*w(i,j,k) + (1.0_db &
			       + 3.0_db*v(i,j,k)*(1.0_db + v(i,j,k)))*w(i,j,k)**2.0_db - u(i,j,k)*(1.0_db &
			       + 3.0_db*v(i,j,k)*(1.0_db + v(i,j,k)))*(1.0_db + 3.0_db*w(i,j,k)*(1.0_db &
			       + w(i,j,k))) + u(i,j,k)**2.0_db*(1.0_db + 3.0_db*v(i,j,k)*(1.0_db &
			       + v(i,j,k)))*(1.0_db + 3.0_db*w(i,j,k)*(1.0_db + w(i,j,k)))))/216.0_db
!26
#ifdef GHOSTONE   
				  fneq1=3.0_db*(pxx(i,j,k) - 3.0_db*pxy(i,j,k) &
				   - 3.0_db*pxz(i,j,k) + pyy(i,j,k) + 3.0_db*pyz(i,j,k) + pzz(i,j,k))  
#else 
				  fneq1=3.0_db*(pxy(i,j,k)*(-3.0_db + 6.0_db*u(i,j,k)) &
				   + pxz(i,j,k)*(-3.0_db + 6.0_db*u(i,j,k)) + (pyy(i,j,k) &
				   + 3.0_db*pyz(i,j,k) + pzz(i,j,k))*(1.0_db + 3.0_db*(-1.0_db &
				   + u(i,j,k))*u(i,j,k)) + 3.0_db*(v(i,j,k)*(pxy(i,j,k)*(-2.0_db &
				   + 6.0_db*u(i,j,k)) + 3.0_db*pxz(i,j,k)*(-1.0_db + 2.0_db*u(i,j,k))*(1.0_db &
				   + v(i,j,k)) + (1.0_db + 3.0_db*(-1.0_db + u(i,j,k))*u(i,j,k))*(2.0_db*pyz(i,j,k) &
				   + pzz(i,j,k) + pzz(i,j,k)*v(i,j,k))) + (-2.0_db*pxz(i,j,k) &
				   + pyy(i,j,k) + 2.0_db*pyz(i,j,k) + 3.0_db*((2.0_db*pyy(i,j,k) &
				   + 3.0_db*pyz(i,j,k))*(-1.0_db + u(i,j,k))*u(i,j,k) + 2.0_db*pyz(i,j,k)*v(i,j,k) &
				   + 6.0_db*pyz(i,j,k)*(-1.0_db + u(i,j,k))*u(i,j,k)*v(i,j,k) &
				   - 3.0_db*pxz(i,j,k)*v(i,j,k)*(1.0_db + v(i,j,k)) + pxz(i,j,k)*u(i,j,k)*(2.0_db &
				   + 6.0_db*v(i,j,k)*(1.0_db + v(i,j,k))) + pxy(i,j,k)*(-1.0_db &
				   - 4.0_db*v(i,j,k) + 2.0_db*u(i,j,k)*(2.0_db + 5.0_db*v(i,j,k)))))*w(i,j,k) &
				   + (pyy(i,j,k)*(2.0_db + 9.0_db*(-1.0_db + u(i,j,k))*u(i,j,k)) &
				   + 6.0_db*pxy(i,j,k)*(-1.0_db - 3.0_db*v(i,j,k) + u(i,j,k)*(3.0_db &
				   + 7.0_db*v(i,j,k))))*w(i,j,k)**2.0_db) + pxx(i,j,k)*(1.0_db &
				   + 3.0_db*w(i,j,k) + 6.0_db*w(i,j,k)**2.0_db + 3.0_db*v(i,j,k)*(1.0_db &
				   + 3.0_db*w(i,j,k))**2.0_db + 3.0_db*v(i,j,k)**2.0_db*(2.0_db &
				   + 3.0_db*w(i,j,k)*(3.0_db + 4.0_db*w(i,j,k)))))
#endif
#endif
!26
#ifdef PHASE_CHANGE
                  feq=feq+p3*src(i,j,k)*(invrho_b-invrho_r)
#endif
				  F_discr=(forcex*(-1.0_db + 2.0_db*u(i,j,k) - 3.0_db*v(i,j,k) &
				   - 3.0_db*w(i,j,k)) + forcez*(1.0_db - 3.0_db*u(i,j,k) &
				   + 3.0_db*v(i,j,k) + 2.0_db*w(i,j,k)) + forcey*(1.0_db &
				   - 3.0_db*u(i,j,k) + 2.0_db*v(i,j,k) + 3.0_db*w(i,j,k)))/(72.0_db*rhophi_loc)
                  f(i-1,j+1,k+1,26)=feq + (1.0_db-omega_loc)*fneq1*p3 + 0.5_db*(F_discr)



#ifdef TWOCOMPONENT
                feq= (( p1*(arr_x(i+1,j,k)-arr_x(i-1,j,k)) + &
                     p2*( (arr_x(i+1,j+1,k)-arr_x(i-1,j-1,k))+(arr_x(i+1,j-1,k)-arr_x(i-1,j+1,k))+(arr_x(i+1,j,k+1)-arr_x(i-1,j,k-1))+(arr_x(i+1,j,k-1)-arr_x(i-1,j,k+1)) )  + &
                     p3*((arr_x(i+1,j+1,k+1)-arr_x(i-1,j-1,k-1))+(arr_x(i+1,j-1,k-1)-arr_x(i-1,j+1,k+1))+(arr_x(i+1,j-1,k+1)-arr_x(i-1,j+1,k-1))+(arr_x(i+1,j+1,k-1)-arr_x(i-1,j-1,k+1))))+ &

                     (p1*(arr_y(i,j+1,k)-arr_y(i,j-1,k)) + &
                     p2*((arr_y(i+1,j+1,k)-arr_y(i-1,j-1,k))+(arr_y(i-1,j+1,k)-arr_y(i+1,j-1,k))+(arr_y(i,j+1,k+1)-arr_y(i,j-1,k-1))+(arr_y(i,j+1,k-1)-arr_y(i,j-1,k+1)) ) + &
                     p3*((arr_y(i+1,j+1,k+1)-arr_y(i-1,j-1,k-1))+(arr_y(i-1,j+1,k-1)-arr_y(i+1,j-1,k+1))+(arr_y(i+1,j+1,k-1)-arr_y(i-1,j-1,k+1))+(arr_y(i-1,j+1,k+1)-arr_y(i+1,j-1,k-1))))+ &

                     (p1*(arr_z(i,j,k+1)-arr_z(i,j,k-1)) + &
                     p2*((arr_z(i+1,j,k+1)-arr_z(i-1,j,k-1))+(arr_z(i-1,j,k+1)-arr_z(i+1,j,k-1))+(arr_z(i,j+1,k+1)-arr_z(i,j-1,k-1))+(arr_z(i,j-1,k+1)-arr_z(i,j+1,k-1)) ) + &
                     p3*((arr_z(i+1,j+1,k+1)-arr_z(i-1,j-1,k-1))+(arr_z(i-1,j-1,k+1)-arr_z(i+1,j+1,k-1))+(arr_z(i+1,j-1,k+1)-arr_z(i-1,j+1,k-1))+(arr_z(i-1,j+1,k+1)-arr_z(i+1,j-1,k-1)))))
				  
				  feq= -sharp_c*feq/cssq
				  
				  !reuse gradrhox,gradrhoy,gradrhoz as local velocity (reusing variables is saving register memory)
				  !reuse gradfix,gradfiy,gradfiz
#ifdef SECFD
                  gradrhox=u(i,j,k) + 0.5_db*forcex/rhophi_loc
				  gradrhoy=v(i,j,k) + 0.5_db*forcey/rhophi_loc
				  gradrhoz=w(i,j,k) + 0.5_db*forcez/rhophi_loc
#else				  
				  gradrhox=u(i,j,k)
				  gradrhoy=v(i,j,k)
				  gradrhoz=w(i,j,k)
#endif
#ifdef UPWIND
				  !upwind
                  selphi(i,j,k,flop) = selphi(i,j,k,flip) - min(0.0_db,gradrhox)*(selphi(i+1,j,k,flip)-selphi(i,j,k,flip)) - max(0.0_db,gradrhox)*(selphi(i,j,k,flip)-selphi(i-1,j,k,flip)) &
				   - min(0.0_db,gradrhoy)*(selphi(i,j+1,k,flip)-selphi(i,j,k,flip)) - max(0.0_db,gradrhoy)*(selphi(i,j,k,flip)-selphi(i,j-1,k,flip)) &
				   - min(0.0_db,gradrhoz)*(selphi(i,j,k+1,flip)-selphi(i,j,k,flip)) - max(0.0_db,gradrhoz)*(selphi(i,j,k,flip)-selphi(i,j,k-1,flip))+ tau_diff*(lap_phi(i,j,k)) + feq

#endif

!**************************************aritra**************************************

#ifdef WENO
    !WENO5 
	!gradrhox = u, gradrhoy = w, gradrhoz = v
	!assuming dxi = dyi = dzi = dt =1

	!
        !!!!!!x-dirction!!!!!!!!!!!!

		!! i + 1/2
         
		  
		 beta1 = (13.0_db / 12.0_db) * (selphi(i-2,j,k,flip) - 2.0_db*selphi(i-1,j,k,flip) +     selphi(i,j,k,flip))**2 + &
          	       (1.0_db / 4.0_db) * (selphi(i-2,j,k,flip) - 4.0_db*selphi(i-1,j,k,flip) + 3.0_db*selphi(i,j,k,flip))**2

         beta2 = (13.0_db / 12.0_db) * (selphi(i-1,j,k,flip) - 2.0_db*selphi(i,j,k,flip) + selphi(i+1,j,k,flip))**2 + &
          	       (1.0_db / 4.0_db) * (selphi(i+1,j,k,flip) -     selphi(i-1,j,k,flip))**2

         beta3 = (13.0_db / 12.0_db) * (selphi(i,j,k,flip) - 2.0_db*selphi(i+1,j,k,flip) + selphi(i+2,j,k,flip))**2 + &
               (1.0_db / 4.0_db) * (3.0_db*selphi(i,j,k,flip) - 4.0_db*selphi(i+1,j,k,flip) + selphi(i+2,j,k,flip))**2
 
	
          fm1 = 0.1_db/(beta1+1e-7)**2   !alpha1
          f0 = 0.6_db/(beta2+1e-7)**2    !alpha2
          fp1 = 0.3_db/(beta3+1e-7)**2   !alpha3

          sum_we = fm1+f0+fp1

          we1 = fm1/sum_we
          we2 = f0/sum_we
          we3 = fp1/sum_we

		  fp2 = we1*(c11*selphi(i-2,j,k,flip) + c21*selphi(i-1,j,k,flip) + c31*selphi(i,j,k,flip)) + &
		        we2*(c12*selphi(i-1,j,k,flip) + c22*selphi(i,j,k,flip)   + c32*selphi(i+1,j,k,flip)) + &
				we3*(c13*selphi(i,j,k,flip)   + c23*selphi(i+1,j,k,flip) + c33*selphi(i+2,j,k,flip))

		 
		  !! i - 1/2
         
		  beta1 = (13.0_db / 12.0_db) * (selphi(i-3,j,k,flip) - 2.0_db*selphi(i-2,j,k,flip) +     selphi(i-1,j,k,flip))**2 + &
          	       (1.0_db / 4.0_db) * (selphi(i-3,j,k,flip) - 4.0_db*selphi(i-2,j,k,flip) + 3.0_db*selphi(i-1,j,k,flip))**2

         beta2 = (13.0_db / 12.0_db) * (selphi(i-2,j,k,flip) - 2.0_db*selphi(i-1,j,k,flip) + selphi(i,j,k,flip))**2 + &
          	       (1.0_db / 4.0_db) * (selphi(i,j,k,flip) -     selphi(i-2,j,k,flip))**2

         beta3 = (13.0_db / 12.0_db) * (selphi(i-1,j,k,flip) - 2.0_db*selphi(i,j,k,flip) + selphi(i+1,j,k,flip))**2 + &
               (1.0_db / 4.0_db) * (3.0_db*selphi(i-1,j,k,flip) - 4.0_db*selphi(i,j,k,flip) + selphi(i+1,j,k,flip))**2
 
	
          fm1 = 0.1_db/(beta1+1e-7)**2   !alpha1
          f0 = 0.6_db/(beta2+1e-7)**2    !alpha2
          fp1 = 0.3_db/(beta3+1e-7)**2   !alpha3

          sum_we = fm1+f0+fp1

          we1 = fm1/sum_we
          we2 = f0/sum_we
          we3 = fp1/sum_we

		  fm2 = we1*(c11*selphi(i-3,j,k,flip) + c21*selphi(i-2,j,k,flip) + c31*selphi(i-1,j,k,flip)) + &
		        we2*(c12*selphi(i-2,j,k,flip) + c22*selphi(i-1,j,k,flip) + c32*selphi(i,j,k,flip)) + &
				we3*(c13*selphi(i-1,j,k,flip) + c23*selphi(i,j,k,flip)   + c33*selphi(i+1,j,k,flip))
		 
	
          !!the gradient in x-direction

		  gradfix = fp2 - fm2

		  

          !
         !!!!!!Y-dirction!!!!!!!!!!!!

		 !! j + 1/2
         
		  
		 beta1 = (13.0_db / 12.0_db) * (selphi(i,j-2,k,flip) - 2.0_db*selphi(i,j-1,k,flip) +     selphi(i,j,k,flip))**2 + &
          	       (1.0_db / 4.0_db) * (selphi(i,j-2,k,flip) - 4.0_db*selphi(i,j-1,k,flip) + 3.0_db*selphi(i,j,k,flip))**2

         beta2 = (13.0_db / 12.0_db) * (selphi(i,j-1,k,flip) - 2.0_db*selphi(i,j,k,flip) + selphi(i,j+1,k,flip))**2 + &
          	       (1.0_db / 4.0_db) * (selphi(i,j+1,k,flip) -     selphi(i,j-1,k,flip))**2

         beta3 = (13.0_db / 12.0_db) * (selphi(i,j,k,flip) - 2.0_db*selphi(i,j+1,k,flip) + selphi(i,j+2,k,flip))**2 + &
               (1.0_db / 4.0_db) * (3.0_db*selphi(i,j,k,flip) - 4.0_db*selphi(i,j+1,k,flip) + selphi(i,j+2,k,flip))**2
 
	
          fm1 = 0.1_db/(beta1+1e-7)**2   !alpha1
          f0 = 0.6_db/(beta2+1e-7)**2    !alpha2
          fp1 = 0.3_db/(beta3+1e-7)**2   !alpha3

          sum_we = fm1+f0+fp1

          we1 = fm1/sum_we
          we2 = f0/sum_we
          we3 = fp1/sum_we

		  fp2 = we1*(c11*selphi(i,j-2,k,flip) + c21*selphi(i,j-1,k,flip) + c31*selphi(i,j,k,flip)) + &
		        we2*(c12*selphi(i,j-1,k,flip) + c22*selphi(i,j,k,flip)   + c32*selphi(i,j+1,k,flip)) + &
				we3*(c13*selphi(i,j,k,flip)   + c23*selphi(i,j+1,k,flip) + c33*selphi(i,j+2,k,flip))

		 
		  !! j - 1/2
         
		  beta1 = (13.0_db / 12.0_db) * (selphi(i,j-3,k,flip) - 2.0_db*selphi(i,j-2,k,flip) +     selphi(i,j-1,k,flip))**2 + &
          	       (1.0_db / 4.0_db) * (selphi(i,j-3,k,flip) - 4.0_db*selphi(i,j-2,k,flip) + 3.0_db*selphi(i,j-1,k,flip))**2

         beta2 = (13.0_db / 12.0_db) * (selphi(i,j-2,k,flip) - 2.0_db*selphi(i,j-1,k,flip) + selphi(i,j,k,flip))**2 + &
          	       (1.0_db / 4.0_db) * (selphi(i,j,k,flip) -     selphi(i,j-2,k,flip))**2

         beta3 = (13.0_db / 12.0_db) * (selphi(i,j-1,k,flip) - 2.0_db*selphi(i,j,k,flip) + selphi(i,j+1,k,flip))**2 + &
               (1.0_db / 4.0_db) * (3.0_db*selphi(i,j-1,k,flip) - 4.0_db*selphi(i,j,k,flip) + selphi(i,j+1,k,flip))**2
 
	
          fm1 = 0.1_db/(beta1+1e-7)**2   !alpha1
          f0 = 0.6_db/(beta2+1e-7)**2    !alpha2
          fp1 = 0.3_db/(beta3+1e-7)**2   !alpha3

          sum_we = fm1+f0+fp1

          we1 = fm1/sum_we
          we2 = f0/sum_we
          we3 = fp1/sum_we

		  fm2 = we1*(c11*selphi(i,j-3,k,flip) + c21*selphi(i,j-2,k,flip) + c31*selphi(i,j-1,k,flip)) + &
		        we2*(c12*selphi(i,j-2,k,flip) + c22*selphi(i,j-1,k,flip) + c32*selphi(i,j,k,flip)) + &
				we3*(c13*selphi(i,j-1,k,flip) + c23*selphi(i,j,k,flip)   + c33*selphi(i,j+1,k,flip))
		 
	
          !!the gradient in y-direction

		  gradfiy = fp2 - fm2


		

         !!!!!!z-dirction!!!!!!!!!!!!

		 
         !! k + 1/2
         
		  
		 beta1 = (13.0_db / 12.0_db) * (selphi(i,j,k-2,flip) - 2.0_db*selphi(i,j,k-1,flip) +     selphi(i,j,k,flip))**2 + &
          	       (1.0_db / 4.0_db) * (selphi(i,j,k-2,flip) - 4.0_db*selphi(i,j,k-1,flip) + 3.0_db*selphi(i,j,k,flip))**2

         beta2 = (13.0_db / 12.0_db) * (selphi(i,j,k-1,flip) - 2.0_db*selphi(i,j,k,flip) + selphi(i,j,k+1,flip))**2 + &
          	       (1.0_db / 4.0_db) * (selphi(i,j,k+1,flip) -     selphi(i,j,k-1,flip))**2

         beta3 = (13.0_db / 12.0_db) * (selphi(i,j,k,flip) - 2.0_db*selphi(i,j,k+1,flip) + selphi(i,j,k+2,flip))**2 + &
               (1.0_db / 4.0_db) * (3.0_db*selphi(i,j,k,flip) - 4.0_db*selphi(i,j,k+1,flip) + selphi(i,j,k+2,flip))**2
 
	
          fm1 = 0.1_db/(beta1+1e-7)**2   !alpha1
          f0 = 0.6_db/(beta2+1e-7)**2    !alpha2
          fp1 = 0.3_db/(beta3+1e-7)**2   !alpha3

          sum_we = fm1+f0+fp1

          we1 = fm1/sum_we
          we2 = f0/sum_we
          we3 = fp1/sum_we

		  fp2 = we1*(c11*selphi(i,j,k-2,flip) + c21*selphi(i,j,k-1,flip) + c31*selphi(i,j,k,flip)) + &
		        we2*(c12*selphi(i,j,k-1,flip) + c22*selphi(i,j,k,flip)   + c32*selphi(i,j,k+1,flip)) + &
				we3*(c13*selphi(i,j,k,flip)   + c23*selphi(i,j,k+1,flip) + c33*selphi(i,j,k+2,flip))

		 
		  !! k - 1/2
         
		  beta1 = (13.0_db / 12.0_db) * (selphi(i,j,k-3,flip) - 2.0_db*selphi(i,j,k-2,flip) +     selphi(i,j,k-1,flip))**2 + &
          	       (1.0_db / 4.0_db) * (selphi(i,j,k-3,flip) - 4.0_db*selphi(i,j,k-2,flip) + 3.0_db*selphi(i,j,k-1,flip))**2

         beta2 = (13.0_db / 12.0_db) * (selphi(i,j,k-2,flip) - 2.0_db*selphi(i,j,k-1,flip) + selphi(i,j,k,flip))**2 + &
          	       (1.0_db / 4.0_db) * (selphi(i,j,k,flip) -     selphi(i,j,k-2,flip))**2

         beta3 = (13.0_db / 12.0_db) * (selphi(i,j,k-1,flip) - 2.0_db*selphi(i,j,k,flip) + selphi(i,j,k+1,flip))**2 + &
               (1.0_db / 4.0_db) * (3.0*selphi(i,j,k-1,flip) - 4.0_db*selphi(i,j,k,flip) + selphi(i,j,k+1,flip))**2
 
	
          fm1 = 0.1_db/(beta1+1e-7)**2   !alpha1
          f0 = 0.6_db/(beta2+1e-7)**2    !alpha2
          fp1 = 0.3_db/(beta3+1e-7)**2   !alpha3

          sum_we = fm1+f0+fp1

          we1 = fm1/sum_we
          we2 = f0/sum_we
          we3 = fp1/sum_we

		  fm2 = we1*(c11*selphi(i,j,k-3,flip) + c21*selphi(i,j,k-2,flip) + c31*selphi(i,j,k-1,flip)) + &
		        we2*(c12*selphi(i,j,k-2,flip) + c22*selphi(i,j,k-1,flip) + c32*selphi(i,j,k,flip)) + &
				we3*(c13*selphi(i,j,k-1,flip) + c23*selphi(i,j,k,flip)   + c33*selphi(i,j,k+1,flip))
		 
	
          !!the gradient in z-direction

		  gradfiz = fp2 - fm2
		  
          selphi(i,j,k,flop) = selphi(i,j,k,flip) &
           - gradrhox*gradfix - gradrhoy*gradfiy - gradrhoz*gradfiz + & 
		   tau_diff*(lap_phi(i,j,k)) + feq	
#endif	
!**************************************aritra end**************************************

#ifdef CENTRAL

                  gradfix=selphi(i+1,j,k,flip)-selphi(i-1,j,k,flip)
				  gradfiy=selphi(i,j+1,k,flip)-selphi(i,j-1,k,flip)
		          gradfiz=selphi(i,j,k+1,flip)-selphi(i,j,k-1,flip)
                  selphi(i,j,k,flop) = selphi(i,j,k,flip) &
                   - gradrhox*0.5_db*(gradfix) - gradrhoy*0.5_db*(gradfiy) &
                   - gradrhoz*0.5_db*(gradfiz) + tau_diff*(lap_phi(i,j,k)) + feq
#endif	

#ifdef MONOD
			    S_mono = mu_max * selphi(i,j,k,flip)/(Ks + selphi(i,j,k,flip)) * selphi(i,j,k,flip) * (1.0_db - selphi(i,j,k,flip))
				selphi(i,j,k,flop)=selphi(i,j,k,flop) + S_mono
				selphi(i,j,k,flop) = min(1.0_db, max(0.0_db, selphi(i,j,k,flop)));
				
#endif			 
#endif
               

   endsubroutine fused_LB_kernel



endmodule
