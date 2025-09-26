#include "defines.h"
#if !defined(_KERNELCUDA) && !defined(_OPENACC)  
#error "To use this module the macros _KERNELCUDA and _OPENACC should be defined."
#endif


module lb_cuda_kernels

   use vars
   use cudafor
   use mpi_template, only: coords,dostop,doerror,mydev,myrank,nprocs,nbuff,nbuffbvec

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
       !$acc& ,selphi,modgrad,arr_x,arr_y,arr_z,lap_phi,normx,normy,normz &
#ifdef REPULSIVE_FLUX
       !$acc& ,Jx,Jy,Jz &
#endif          
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
       !$acc& ,visc1,visc2,rho_r,rho_b,invrho_r,invrho_b,beta,kapp,sigma,sharp_c,tau_diff &
#endif   
#if defined(ELASTIC_FORCE)
       !$acc& ,lambda_rel,k_elastic,u_ref,v_ref,w_ref &
#endif
       !$acc& ,fx,fy,fz)
       call moments_LB_kernel<<<dimGrid, dimBlock>>>(flop,nx,ny,nz,coords,isfluid,f &
       ,rho,u,v,w,pxx,pyy,pzz,pxy,pxz,pyz,fux,fvy,fwz &
#ifdef TWOCOMPONENT
       ,selphi,modgrad,arr_x,arr_y,arr_z,lap_phi,normx,normy,normz &
#ifdef REPULSIVE_FLUX
       ,Jx,Jy,Jz &
#endif           
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
       ,visc1,visc2,rho_r,rho_b,invrho_r,invrho_b,beta,kapp,sigma,sharp_c,tau_diff &       
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
       !$acc& ,selphi,modgrad,lap_phi,normx,normy,normz,arr_x,arr_y,arr_z &
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
#ifdef MONOD
	   !$acc& ,mu_max,Ks &
#endif
#endif   
       !$acc& ,omega,fx,fy,fz)
      call fused_LB_kernel<<<dimGrid, dimBlock>>>(flip,flop,nx,ny,nz,coords,isfluid,f &
       ,rho,u,v,w,pxx,pyy,pzz,pxy,pxz,pyz,fux,fvy,fwz &
#ifdef TWOCOMPONENT
       ,selphi,modgrad,lap_phi,normx,normy,normz,arr_x,arr_y,arr_z &
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
       ,selphi,modgrad,arr_x,arr_y,arr_z,lap_phi,normx,normy,normz &
#ifdef REPULSIVE_FLUX
       ,Jx,Jy,Jz &
#endif        
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
       ,visc1,visc2,rho_r,rho_b,invrho_r,invrho_b,beta,kapp,sigma,sharp_c,tau_diff &
#endif   
#if defined(ELASTIC_FORCE)
       ,lambda_rel,k_elastic,u_ref,v_ref,w_ref &
#endif
       ,fx,fy,fz)
 

      implicit none
      
      integer :: flop,nx,ny,nz
      integer, dimension(3) :: coords
      integer(kind=isf), dimension(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff) :: isfluid
      real(kind=db), dimension(0:nx+1,0:ny+1,0:nz+1,0:nlinks) :: f
      real(kind=db), dimension(0:nx+1,0:ny+1,0:nz+1) :: rho
      real(kind=db), dimension(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff) :: u,v,w
      real(kind=db), dimension(0:nx+1,0:ny+1,0:nz+1) :: pxx,pyy,pzz,pxy,pxz,pyz
      real(kind=db), dimension(1:nx,1:ny,1:nz) :: fux,fvy,fwz
#ifdef TWOCOMPONENT
      real(kind=db), dimension(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff,2) :: selphi
      real(kind=db), dimension(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff) :: modgrad
      real(kind=db), dimension(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff) :: arr_x,arr_y,arr_z
      real(kind=db), dimension(1:nx,1:ny,1:nz) :: lap_phi
      real(kind=db), dimension(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff) :: normx,normy,normz
#ifdef REPULSIVE_FLUX
      real(kind=db), dimension(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff) :: Jx,Jy,Jz
#endif      
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
      real(kind=db) :: visc1,visc2,rho_r,rho_b,invrho_r,invrho_b,beta,kapp,sigma,sharp_c,tau_diff     
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
      
       integer :: i,j,k,l
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

                 fux(i,j,k)=0.0_db
				 fvy(i,j,k)=0.0_db
				 fwz(i,j,k)=0.0_db
				 forcex=0.0_db
				 forcey=0.0_db
				 forcez=0.0_db
				 
#ifdef DENSRATIO
                  rhophi_loc = rhophi(i,j,k)
#else
                  rhophi_loc = 1.0_db !rho(i,j,k)
#endif	
			 
#ifdef TWOCOMPONENT		
#ifdef CSF
                   !csf
				   gradfix=normx(i,j,k)*modgrad(i,j,k)
                   gradfiy=normy(i,j,k)*modgrad(i,j,k)
                   gradfiz=normz(i,j,k)*modgrad(i,j,k)
                   curvature=3.0_db*(p1*(normx(i+1,j,k)-normx(i-1,j,k)) + &
                    p2*( (normx(i+1,j+1,k)-normx(i-1,j-1,k))+(normx(i+1,j-1,k)-normx(i-1,j+1,k))+(normx(i+1,j,k+1)-normx(i-1,j,k-1))+(normx(i+1,j,k-1)-normx(i-1,j,k+1)) )  + &
                    p3*((normx(i+1,j+1,k+1)-normx(i-1,j-1,k-1))+(normx(i+1,j-1,k-1)-normx(i-1,j+1,k+1))+(normx(i+1,j-1,k+1)-normx(i-1,j+1,k-1))+(normx(i+1,j+1,k-1)-normx(i-1,j-1,k+1)))) + &
                    3.0_db*(p1*(normy(i,j+1,k)-normy(i,j-1,k)) + &
                    p2*((normy(i+1,j+1,k)-normy(i-1,j-1,k))+(normy(i-1,j+1,k)-normy(i+1,j-1,k))+(normy(i,j+1,k+1)-normy(i,j-1,k-1))+(normy(i,j+1,k-1)-normy(i,j-1,k+1)) ) + &
                    p3*((normy(i+1,j+1,k+1)-normy(i-1,j-1,k-1))+(normy(i-1,j+1,k-1)-normy(i+1,j-1,k+1))+(normy(i+1,j+1,k-1)-normy(i-1,j-1,k+1))+(normy(i-1,j+1,k+1)-normy(i+1,j-1,k-1)))) +&
                    3.0_db*(p1*(normz(i,j,k+1)-normz(i,j,k-1)) + &
                    p2*((normz(i+1,j,k+1)-normz(i-1,j,k-1))+(normz(i-1,j,k+1)-normz(i+1,j,k-1))+(normz(i,j+1,k+1)-normz(i,j-1,k-1))+(normz(i,j-1,k+1)-normz(i,j+1,k-1)) ) + &
                    p3*((normz(i+1,j+1,k+1)-normz(i-1,j-1,k-1))+(normz(i-1,j-1,k+1)-normz(i+1,j+1,k-1))+(normz(i+1,j-1,k+1)-normz(i-1,j+1,k-1))+(normz(i-1,j+1,k+1)-normz(i+1,j-1,k-1))))

                    fux(i,j,k)=-eps1*sigma*curvature*gradfix*modgrad(i,j,k)
				    fvy(i,j,k)=-eps1*sigma*curvature*gradfiy*modgrad(i,j,k) 
				    fwz(i,j,k)=-eps1*sigma*curvature*gradfiz*modgrad(i,j,k)
				   
#endif
#ifdef JAQMIN			   
				   !jaqmin 
				   gradfix=normx(i,j,k)*modgrad(i,j,k)
				   gradfiy=normy(i,j,k)*modgrad(i,j,k)
				   gradfiz=normz(i,j,k)*modgrad(i,j,k)
				   fux(i,j,k)=(4.0_db*beta*selphi(i,j,k,flop)*(selphi(i,j,k,flop)-1.0_db)*(selphi(i,j,k,flop)-0.5_db) - kapp*lap_phi(i,j,k))*gradfix
				   fvy(i,j,k)=(4.0_db*beta*selphi(i,j,k,flop)*(selphi(i,j,k,flop)-1.0_db)*(selphi(i,j,k,flop)-0.5_db) - kapp*lap_phi(i,j,k))*gradfiy
				   fwz(i,j,k)=(4.0_db*beta*selphi(i,j,k,flop)*(selphi(i,j,k,flop)-1.0_db)*(selphi(i,j,k,flop)-0.5_db) - kapp*lap_phi(i,j,k))*gradfiz
				   				   
#endif

#ifdef REPULSIVE_FLUX
				  Jx(i,j,k)=Jx(i,j,k)*rhophi_loc 
				  Jy(i,j,k)=Jy(i,j,k)*rhophi_loc 
				  Jz(i,j,k)=Jz(i,j,k)*rhophi_loc 
				  if(abs(Jx(i,j,k))>1.0d-3) Jx(i,j,k)=1.0d-3*sign(1.0,Jx(i,j,k))!Jx(i,j,k)*0.1_db
				  if(abs(Jy(i,j,k))>1.0d-3) Jy(i,j,k)=1.0d-3*sign(1.0,Jy(i,j,k))!Jy(i,j,k)*0.1_db
				  if(abs(Jz(i,j,k))>1.0d-3) Jz(i,j,k)=1.0d-3*sign(1.0,Jz(i,j,k))!Jz(i,j,k)*0.1_db
				  fux(i,j,k)=fux(i,j,k) + Jx(i,j,k)*rhophi_loc
				  fvy(i,j,k)=fvy(i,j,k) + Jy(i,j,k)*rhophi_loc
				  fwz(i,j,k)=fwz(i,j,k) + Jz(i,j,k)*rhophi_loc
#endif
#endif

#if defined(PLUG_FLOW)   
                  
				  fwz(i,j,k) = fwz(i,j,k) + rhophi_loc*fz !fwz(i,j,k)=fwz(i,j,k) + selphi(i,j,k,flop)*fz ! if mnulticomponent/phase (rhophi_loc-rho_r or rho_b)*fz 	
#endif
#if defined(MULTIHIT)   
                  fux(i,j,k)=fux(i,j,k) + rhophi_loc*ABCx(i,j,k) !+ AAA*sin(k_zero*gk) + AAA*sin(k_zero*gj)  
				  fvy(i,j,k)=fvy(i,j,k) + rhophi_loc*ABCy(i,j,k) !+ AAA*sin(k_zero*gi) + AAA*sin(k_zero*gk)
				  fwz(i,j,k)=fwz(i,j,k) + rhophi_loc*ABCz(i,j,k) !+ AAA*sin(k_zero*gj) + AAA*sin(k_zero*gi) 	
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
			  

                  pxx(i,j,k)=0.0_db
                  pyy(i,j,k)=0.0_db
                  pzz(i,j,k)=0.0_db
                  pxy(i,j,k)=0.0_db
                  pxz(i,j,k)=0.0_db
                  pyz(i,j,k)=0.0_db
                  !1-2
                  !*1
                  ! 2nd order
				  uu=0.5_db*(u(i,j,k)*u(i,j,k)+v(i,j,k)*v(i,j,k)+w(i,j,k)*w(i,j,k))/cssq
				  
                  do l=1,nlinks
                     udotc=(u(i,j,k)*dex(l) + v(i,j,k)*dey(l)+ w(i,j,k)*dez(l))/cssq
					 feq=p(l)*(rho(i,j,k) + (udotc+0.5_db*udotc*udotc - uu))
                     !
                     fneq1=f(i,j,k,l)-feq !-0.5_db*F_discr) 
                     pxx(i,j,k)=pxx(i,j,k)+ fneq1*dex(l)*dex(l)
                     pyy(i,j,k)=pyy(i,j,k)+ fneq1*dey(l)*dey(l)
                     pzz(i,j,k)=pzz(i,j,k)+ fneq1*dez(l)*dez(l)
                     pxy(i,j,k)=pxy(i,j,k)+ fneq1*dex(l)*dey(l)
                     pxz(i,j,k)=pxz(i,j,k)+ fneq1*dex(l)*dez(l)
                     pyz(i,j,k)=pyz(i,j,k)+ fneq1*dey(l)*dez(l)
                  enddo

				  visc_loc=(rho_r*visc1*selphi(i,j,k,flop)+(1.0_db-selphi(i,j,k,flop))*visc2*rho_b)/rhophi_loc
				  
                  tau_loc=(visc_loc/cssq + 0.5_db) !è una tau
				  
				  forcex=forcex - (visc_loc/(tau_loc*cssq))*(pxx(i,j,k)*gradrhox + pxy(i,j,k)*gradrhoy + pxz(i,j,k)*gradrhoz)
				  forcey=forcey - (visc_loc/(tau_loc*cssq))*(pyy(i,j,k)*gradrhoy + pxy(i,j,k)*gradrhox + pyz(i,j,k)*gradrhoz)
				  forcez=forcez - (visc_loc/(tau_loc*cssq))*(pzz(i,j,k)*gradrhoz + pxz(i,j,k)*gradrhox + pyz(i,j,k)*gradrhoy)
#endif	
                  !I compute the new velocities
				  u(i,j,k) = ((f(i,j,k,1)+f(i,j,k,7)+f(i,j,k,9)+f(i,j,k,15)+f(i,j,k,18)+f(i,j,k,19)+f(i,j,k,21)+f(i,j,k,24)+f(i,j,k,25)) &
                     -(f(i,j,k,2)+f(i,j,k,8)+f(i,j,k,10)+f(i,j,k,16)+f(i,j,k,17)+f(i,j,k,20)+f(i,j,k,22)+f(i,j,k,23)+f(i,j,k,26)))

                  v(i,j,k) = ((f(i,j,k,3)+f(i,j,k,7)+f(i,j,k,10)+f(i,j,k,11)+f(i,j,k,13)+f(i,j,k,19)+f(i,j,k,22)+f(i,j,k,24)+f(i,j,k,26)) &
                     -(f(i,j,k,4)+f(i,j,k,8)+f(i,j,k,9)+f(i,j,k,12)+f(i,j,k,14)+f(i,j,k,20)+f(i,j,k,21)+f(i,j,k,23)+f(i,j,k,25)))

                  w(i,j,k) = ((f(i,j,k,5)+f(i,j,k,11)+f(i,j,k,14)+f(i,j,k,15)+f(i,j,k,17)+f(i,j,k,19)+f(i,j,k,21)+f(i,j,k,23)+f(i,j,k,26)) &
                     -(f(i,j,k,6)+f(i,j,k,12)+f(i,j,k,13)+f(i,j,k,16)+f(i,j,k,18)+f(i,j,k,20)+f(i,j,k,22)+f(i,j,k,24)+f(i,j,k,25)))
					 
#ifdef INTERFACE_INCOMP
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
				  u(i,j,k) = u(i,j,k) + 0.5_db*forcex/(rhophi_loc)
                  v(i,j,k) = v(i,j,k) + 0.5_db*forcey/(rhophi_loc)
                  w(i,j,k) = w(i,j,k) + 0.5_db*forcez/(rhophi_loc)
				  
				  u(i,j,k) = u(i,j,k)/(1.0_db - 0.5_db*(rho_r-rho_b)*(tau_diff*lap_phi(i,j,k) + feq)/rhophi_loc )
				  
				  v(i,j,k) = v(i,j,k)/(1.0_db - 0.5_db*(rho_r-rho_b)*(tau_diff*lap_phi(i,j,k) + feq)/rhophi_loc ) 
				  
				  w(i,j,k) = w(i,j,k)/(1.0_db - 0.5_db*(rho_r-rho_b)*(tau_diff*lap_phi(i,j,k) + feq)/rhophi_loc ) 
#else
                  u(i,j,k) = u(i,j,k) + 0.5_db*forcex/rhophi_loc
                  v(i,j,k) = v(i,j,k) + 0.5_db*forcey/rhophi_loc
                  w(i,j,k) = w(i,j,k) + 0.5_db*forcez/rhophi_loc
#endif
#ifdef INTERFACE_INCOMP
			         
					 fux(i,j,k)=fux(i,j,k) - (rho_r-rho_b)*(tau_diff*lap_phi(i,j,k) + feq)*u(i,j,k)
					 fvy(i,j,k)=fvy(i,j,k) - (rho_r-rho_b)*(tau_diff*lap_phi(i,j,k) + feq)*v(i,j,k)
					 fwz(i,j,k)=fwz(i,j,k) - (rho_r-rho_b)*(tau_diff*lap_phi(i,j,k) + feq)*w(i,j,k)
					 
					 forcex=forcex - (rho_r-rho_b)*(tau_diff*lap_phi(i,j,k) + feq)*u(i,j,k)
					 
					 forcey=forcey - (rho_r-rho_b)*(tau_diff*lap_phi(i,j,k) + feq)*v(i,j,k)
					 
					 forcez=forcez - (rho_r-rho_b)*(tau_diff*lap_phi(i,j,k) + feq)*w(i,j,k)
#endif 


!regularized 
				  pxx(i,j,k)=0.0_db
                  pyy(i,j,k)=0.0_db
                  pzz(i,j,k)=0.0_db
                  pxy(i,j,k)=0.0_db
                  pxz(i,j,k)=0.0_db
                  pyz(i,j,k)=0.0_db
				  uu=HALF*(u(i,j,k)*u(i,j,k)+v(i,j,k)*v(i,j,k)+w(i,j,k)*w(i,j,k))/cssq
				  
                  do l=1,nlinks
                     udotc=(u(i,j,k)*dex(l) + v(i,j,k)*dey(l)+ w(i,j,k)*dez(l))/cssq
					 feq=p(l)*(rho(i,j,k) + (udotc+0.5_db*udotc*udotc - uu))
					 F_discr= p(l)*( (dex(l)-u(i,j,k))*(forcex) + (dey(l)-v(i,j,k))*(forcey) + (dez(l)-w(i,j,k))*(forcez) + &
										                  (1.0_db/(cssq))*( (u(i,j,k)*dex(l) + v(i,j,k)*dey(l) + w(i,j,k)*dez(l))*( (forcex)*dex(l) + (forcey)*dey(l) + &
																	(forcez)*dez(l) ) ) )/(cssq*rhophi_loc) 
                     fneq1=f(i,j,k,l)-(feq-0.5_db*F_discr) 
					 

                     pxx(i,j,k)=pxx(i,j,k)+ fneq1*dex(l)*dex(l)
                     pyy(i,j,k)=pyy(i,j,k)+ fneq1*dey(l)*dey(l) 
                     pzz(i,j,k)=pzz(i,j,k)+ fneq1*dez(l)*dez(l)
                     pxy(i,j,k)=pxy(i,j,k)+ fneq1*dex(l)*dey(l)
                     pxz(i,j,k)=pxz(i,j,k)+ fneq1*dex(l)*dez(l)
                     pyz(i,j,k)=pyz(i,j,k)+ fneq1*dey(l)*dez(l)
                  enddo
#if defined(ELASTIC_FORCE)
				  u_ref(i,j,k) = u_ref(i,j,k) + lambda_rel*(u(i,j,k) - u_ref(i,j,k))
				  v_ref(i,j,k) = v_ref(i,j,k) + lambda_rel*(v(i,j,k) - v_ref(i,j,k))
				  w_ref(i,j,k) = w_ref(i,j,k) + lambda_rel*(w(i,j,k) - w_ref(i,j,k))
                  fux(i,j,k)=fux(i,j,k) + rhophi_loc*(u(i,j,k) - u_ref(i,j,k))*k_elastic*lambda_rel*selphi(i,j,k,flop) !+ 
				  fvy(i,j,k)=fvy(i,j,k) + rhophi_loc*(v(i,j,k) - v_ref(i,j,k))*k_elastic*lambda_rel*selphi(i,j,k,flop) !+ 
				  fwz(i,j,k)=fwz(i,j,k) + rhophi_loc*(w(i,j,k) - w_ref(i,j,k))*k_elastic*lambda_rel*selphi(i,j,k,flop)+rhophi_loc*fz  
#endif                  
#if defined(DENSRATIO)			  
				  
				  fux(i,j,k)=fux(i,j,k) - (visc_loc/(tau_loc*cssq))*(pxx(i,j,k)*gradrhox + pxy(i,j,k)*gradrhoy + pxz(i,j,k)*gradrhoz)
				  fvy(i,j,k)=fvy(i,j,k) - (visc_loc/(tau_loc*cssq))*(pyy(i,j,k)*gradrhoy + pxy(i,j,k)*gradrhox + pyz(i,j,k)*gradrhoz)
				  fwz(i,j,k)=fwz(i,j,k) - (visc_loc/(tau_loc*cssq))*(pzz(i,j,k)*gradrhoz + pxz(i,j,k)*gradrhox + pyz(i,j,k)*gradrhoy)
#endif               
 


   endsubroutine moments_LB_kernel  

   attributes(global) subroutine fused_LB_kernel(flip,flop,nx,ny,nz,coords,isfluid,f &
       ,rho,u,v,w,pxx,pyy,pzz,pxy,pxz,pyz,fux,fvy,fwz &
#ifdef TWOCOMPONENT
       ,selphi,modgrad,lap_phi,normx,normy,normz,arr_x,arr_y,arr_z &
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
#ifdef MONOD	
	   ,mu_max,Ks &
#endif
#endif   
       ,omega,fx,fy,fz)

      implicit none
      
      integer :: flip,flop,nx,ny,nz
      integer, dimension(3) :: coords
      integer(kind=isf), dimension(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff) :: isfluid
      real(kind=db), dimension(0:nx+1,0:ny+1,0:nz+1,0:nlinks) :: f
      real(kind=db), dimension(0:nx+1,0:ny+1,0:nz+1) :: rho
      real(kind=db), dimension(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff) :: u,v,w
      real(kind=db), dimension(0:nx+1,0:ny+1,0:nz+1) :: pxx,pyy,pzz,pxy,pxz,pyz
      real(kind=db), dimension(1:nx,1:ny,1:nz) :: fux,fvy,fwz
#ifdef TWOCOMPONENT
      real(kind=db), dimension(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff,2) :: selphi
      real(kind=db), dimension(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff) :: modgrad
      real(kind=db), dimension(1:nx,1:ny,1:nz) :: lap_phi
      real(kind=db), dimension(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff) :: normx,normy,normz,arr_x,arr_y,arr_z
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



#ifdef TWOCOMPONENT

                  visc_loc=(rho_r*visc1*selphi(i,j,k,flip)+(1.0_db-selphi(i,j,k,flip))*visc2*rho_b)/rhophi_loc  

				  
                  tau_loc=(visc_loc/cssq + 0.5_db) !è una tau
				  
                  omega_loc=1.0_db/tau_loc !è una omega
				  
#else
				  omega_loc=omega
#endif


!!!!!!!!!!!!!!!!!!!!!!!!!!0
#ifdef SECOND_ORDER
			      feq=(4.0_db*(2.0_db*rho(i,j,k) - 3.0_db &
			       *(u(i,j,k)**2.0_db + v(i,j,k)**2.0_db + w(i,j,k)**2.0_db)))/27.0_db

#else
!0
			      feq=(8.0_db*rho(i,j,k) - 3.0_db*(4.0_db*w(i,j,k)**2.0_db &
			       + v(i,j,k)**2.0_db*(4.0_db - 6.0_db*w(i,j,k)**2.0_db) &
			       + u(i,j,k)**2.0_db*(-2.0_db + 3.0_db*v(i,j,k)**2.0_db)*(-2.0_db &
			       + 3.0_db*w(i,j,k)**2.0_db)))/27.0_db
!0
#endif 

				  fneq1=(-3.0_db*(pxx(i,j,k) + pyy(i,j,k) + pzz(i,j,k)))/2.0_db


				  F_discr=(-8.0_db*(forcex*u(i,j,k) + forcey*v(i,j,k) &
				   + forcez*w(i,j,k)))/(9.0_db*rhophi_loc)
                  f(i,j,k,0)=feq + (1.0_db-omega_loc)*fneq1*p0 + 0.5_db*(F_discr)
                 

!!!!!!!!!!!!!!!!!!!!!!!!!!1
#ifdef SECOND_ORDER
			      feq=(2.0_db*rho(i,j,k) + 6.0_db*u(i,j,k) &
			       *(1.0_db + u(i,j,k)) - 3.0_db*v(i,j,k)**2.0_db &
			       - 3.0_db*w(i,j,k)**2.0_db)/27.0_db
#else
!1
			      feq=(4.0_db*rho(i,j,k) + 3.0_db*(4.0_db*u(i,j,k)*(1.0_db &
			       + u(i,j,k)) - 2.0_db*(1.0_db + 3.0_db*u(i,j,k)*(1.0_db &
			       + u(i,j,k)))*v(i,j,k)**2.0_db + (1.0_db + 3.0_db*u(i,j,k)*(1.0_db &
			       + u(i,j,k)))*(-2.0_db + 3.0_db*v(i,j,k)**2.0_db)*w(i,j,k)**2.0_db))/54.0_db
!1
#endif
				  fneq1=(3.0_db*(2.0_db*pxx(i,j,k) - pyy(i,j,k) - pzz(i,j,k)))/2.0_db


				  F_discr=(2.0_db*(forcex + 2.0_db*forcex*u(i,j,k) - forcey*v(i,j,k) &
				   - forcez*w(i,j,k)))/(9.0_db*rhophi_loc)
                  f(i+1,j,k,1)=feq + (1.0_db-omega_loc)*fneq1*p1 + 0.5_db*(F_discr)


!!!!!!!!!!!!!!!!!!!!!!!!!!2
#ifdef SECOND_ORDER
			      feq=(2.0_db*rho(i,j,k) + 6.0_db*(-1.0_db &
			       + u(i,j,k))*u(i,j,k) - 3.0_db*v(i,j,k)**2.0_db &
			       - 3.0_db*w(i,j,k)**2.0_db)/27.0_db
#else
!2
			      feq=(4.0_db*rho(i,j,k) + 3.0_db*(4.0_db*(-1.0_db &
			       + u(i,j,k))*u(i,j,k) - 2.0_db*(1.0_db + 3.0_db*(-1.0_db &
			       + u(i,j,k))*u(i,j,k))*v(i,j,k)**2.0_db + (1.0_db &
			       + 3.0_db*(-1.0_db + u(i,j,k))*u(i,j,k))*(-2.0_db &
			       + 3.0_db*v(i,j,k)**2.0_db)*w(i,j,k)**2.0_db))/54.0_db
!2
#endif
				  fneq1=(3.0_db*(2.0_db*pxx(i,j,k) - pyy(i,j,k) - pzz(i,j,k)))/2.0_db


				  F_discr=(-2.0_db*(forcex - 2.0_db*forcex*u(i,j,k) + forcey*v(i,j,k) &
				   + forcez*w(i,j,k)))/(9.0_db*rhophi_loc)
                  f(i-1,j,k,2)=feq + (1.0_db-omega_loc)*fneq1*p1 + 0.5_db*(F_discr)


!!!!!!!!!!!!!!!!!!!!!!!!!!3
#ifdef SECOND_ORDER
			      feq=(2.0_db*rho(i,j,k) - 3.0_db*(u(i,j,k)**2.0_db &
			       - 2.0_db*v(i,j,k)*(1.0_db + v(i,j,k)) + w(i,j,k)**2.0_db))/27.0_db
#else
!3
			      feq=(4.0_db*rho(i,j,k) + 12.0_db*v(i,j,k)*(1.0_db &
			       + v(i,j,k)) - 6.0_db*u(i,j,k)**2.0_db*(1.0_db &
			       + 3.0_db*v(i,j,k)*(1.0_db + v(i,j,k))) + 3.0_db*(-2.0_db &
			       + 3.0_db*u(i,j,k)**2.0_db)*(1.0_db + 3.0_db*v(i,j,k)*(1.0_db &
			       + v(i,j,k)))*w(i,j,k)**2.0_db)/54.0_db
!3
#endif
				  fneq1=(-3.0_db*(pxx(i,j,k) - 2.0_db*pyy(i,j,k) + pzz(i,j,k)))/2.0_db


				  F_discr=(2.0_db*(forcey - forcex*u(i,j,k) + 2.0_db*forcey*v(i,j,k) &
				   - forcez*w(i,j,k)))/(9.0_db*rhophi_loc)
                  f(i,j+1,k,3)=feq + (1.0_db-omega_loc)*fneq1*p1 + 0.5_db*(F_discr)


!!!!!!!!!!!!!!!!!!!!!!!!!!4
#ifdef SECOND_ORDER
			      feq=(2.0_db*rho(i,j,k) - 3.0_db*(u(i,j,k)**2.0_db &
			      - 2.0_db*(-1.0_db + v(i,j,k))*v(i,j,k) + w(i,j,k)**2.0_db))/27.0_db
#else
!4
			      feq=(4.0_db*rho(i,j,k) + 12.0_db*(-1.0_db &
			       + v(i,j,k))*v(i,j,k) - 6.0_db*u(i,j,k)**2.0_db*(1.0_db &
			       + 3.0_db*(-1.0_db + v(i,j,k))*v(i,j,k)) + 3.0_db*(-2.0_db &
			       + 3.0_db*u(i,j,k)**2.0_db)*(1.0_db + 3.0_db*(-1.0_db &
			       + v(i,j,k))*v(i,j,k))*w(i,j,k)**2.0_db)/54.0_db
!4
#endif
				  fneq1=(-3.0_db*(pxx(i,j,k) - 2.0_db*pyy(i,j,k) + pzz(i,j,k)))/2.0_db


				  F_discr=(-2.0_db*(forcey + forcex*u(i,j,k) - 2.0_db*forcey*v(i,j,k) &
				   + forcez*w(i,j,k)))/(9.0_db*rhophi_loc)
                  f(i,j-1,k,4)=feq + (1.0_db-omega_loc)*fneq1*p1 + 0.5_db*(F_discr)


!!!!!!!!!!!!!!!!!!!!!!!!!!5
#ifdef SECOND_ORDER
			      feq=(2.0_db*rho(i,j,k) - 3.0_db*(u(i,j,k)**2.0_db &
			       + v(i,j,k)**2.0_db - 2.0_db*w(i,j,k)*(1.0_db + w(i,j,k))))/27.0_db
#else
!5
			      feq=(4.0_db*rho(i,j,k) + 3.0_db*(4.0_db*w(i,j,k)*(1.0_db &
			       + w(i,j,k)) - 2.0_db*v(i,j,k)**2.0_db*(1.0_db + 3.0_db*w(i,j,k)*(1.0_db &
			       + w(i,j,k))) + u(i,j,k)**2.0_db*(-2.0_db + 3.0_db*v(i,j,k)**2.0_db)*(1.0_db &
			       + 3.0_db*w(i,j,k)*(1.0_db + w(i,j,k)))))/54.0_db
!5
#endif
				  fneq1=(-3.0_db*(pxx(i,j,k) + pyy(i,j,k) - 2.0_db*pzz(i,j,k)))/2.0_db


				  F_discr=(2.0_db*(forcez - forcex*u(i,j,k) - forcey*v(i,j,k) &
				   + 2.0_db*forcez*w(i,j,k)))/(9.0_db*rhophi_loc)
                  f(i,j,k+1,5)= feq + (1.0_db-omega_loc)*fneq1*p1 + 0.5_db*(F_discr)


!!!!!!!!!!!!!!!!!!!!!!!!!!6
#ifdef SECOND_ORDER
			      feq=(2.0_db*rho(i,j,k) - 3.0_db*(u(i,j,k)**2.0_db &
			       + v(i,j,k)**2.0_db - 2.0_db*(-1.0_db + w(i,j,k))*w(i,j,k)))/27.0_db
#else
!6
			      feq=(4.0_db*rho(i,j,k) + 3.0_db*(4.0_db*(-1.0_db &
			       + w(i,j,k))*w(i,j,k) + v(i,j,k)**2.0_db*(-2.0_db - 6.0_db*(-1.0_db &
			       + w(i,j,k))*w(i,j,k)) + u(i,j,k)**2.0_db*(-2.0_db &
			       + 3.0_db*v(i,j,k)**2.0_db)*(1.0_db + 3.0_db*(-1.0_db &
			       + w(i,j,k))*w(i,j,k))))/54.0_db
!6
#endif
				  fneq1=(-3.0_db*(pxx(i,j,k) + pyy(i,j,k) - 2.0_db*pzz(i,j,k)))/2.0_db


				  F_discr=(-2.0_db*(forcez + forcex*u(i,j,k) + forcey*v(i,j,k) &
				   - 2.0_db*forcez*w(i,j,k)))/(9.0_db*rhophi_loc)
                  f(i,j,k-1,6)=feq + (1.0_db-omega_loc)*fneq1*p1 + 0.5_db*(F_discr)


!!!!!!!!!!!!!!!!!!!!!!!!!!7
#ifdef SECOND_ORDER
			      feq=(2.0_db*rho(i,j,k) + 6.0_db*(u(i,j,k) &
			       + u(i,j,k)**2.0_db + v(i,j,k) + 3.0_db*u(i,j,k)*v(i,j,k) &
			       + v(i,j,k)**2.0_db) - 3.0_db*w(i,j,k)**2.0_db)/108.0_db
#else
!7
			      feq=(2.0_db*rho(i,j,k) + 6.0_db*(u(i,j,k) &
			       + v(i,j,k) + v(i,j,k)**2.0_db + 3.0_db*u(i,j,k)*v(i,j,k)*(1.0_db &
			       + v(i,j,k)) + u(i,j,k)**2.0_db*(1.0_db + 3.0_db*v(i,j,k)*(1.0_db &
			       + v(i,j,k)))) - 3.0_db*(1.0_db + 3.0_db*u(i,j,k)*(1.0_db &
			       + u(i,j,k)))*(1.0_db + 3.0_db*v(i,j,k)*(1.0_db &
			       + v(i,j,k)))*w(i,j,k)**2.0_db)/108.0_db
!7
#endif
				  fneq1=3.0_db*(pxx(i,j,k) + 3.0_db*pxy(i,j,k) + pyy(i,j,k)) &
				   - (3.0_db*pzz(i,j,k))/2.0_db


				  F_discr=(forcex + forcey + 2.0_db*forcex*u(i,j,k) &
				   + 3.0_db*forcey*u(i,j,k) + 3.0_db*forcex*v(i,j,k) &
				   + 2.0_db*forcey*v(i,j,k) - forcez*w(i,j,k))/(18.0_db*rhophi_loc)
                  f(i+1,j+1,k,7)=feq + (1.0_db-omega_loc)*fneq1*p2 + 0.5_db*(F_discr)


!!!!!!!!!!!!!!!!!!!!!!!!!!8
#ifdef SECOND_ORDER
			      feq=(2.0_db*rho(i,j,k) + 6.0_db*(u(i,j,k)**2.0_db &
			       + (-1.0_db + v(i,j,k))*v(i,j,k) + u(i,j,k)*(-1.0_db &
			       + 3.0_db*v(i,j,k))) - 3.0_db*w(i,j,k)**2.0_db)/108.0_db
#else
!8
			      feq=(2.0_db*rho(i,j,k) + 6.0_db*((-1.0_db &
			       + v(i,j,k))*v(i,j,k) + u(i,j,k)*(-1.0_db - 3.0_db*(-1.0_db &
			       + v(i,j,k))*v(i,j,k)) + u(i,j,k)**2.0_db*(1.0_db + 3.0_db*(-1.0_db &
			       + v(i,j,k))*v(i,j,k))) - 3.0_db*(1.0_db + 3.0_db*(-1.0_db &
			       + u(i,j,k))*u(i,j,k))*(1.0_db + 3.0_db*(-1.0_db &
			       + v(i,j,k))*v(i,j,k))*w(i,j,k)**2.0_db)/108.0_db
!8
#endif 
				  fneq1=3.0_db*(pxx(i,j,k) + 3.0_db*pxy(i,j,k) + pyy(i,j,k)) &
				   - (3.0_db*pzz(i,j,k))/2.0_db


				  F_discr=(forcex + forcey - 2.0_db*forcex*u(i,j,k) &
				   - 3.0_db*forcey*u(i,j,k) - 3.0_db*forcex*v(i,j,k) &
				   - 2.0_db*forcey*v(i,j,k) + forcez*w(i,j,k))/(-18.0_db*rhophi_loc)
                  f(i-1,j-1,k,8)=feq + (1.0_db-omega_loc)*fneq1*p2 + 0.5_db*(F_discr)


!!!!!!!!!!!!!!!!!!!!!!!!!!9
#ifdef SECOND_ORDER
			      feq=(2.0_db*rho(i,j,k) + 6.0_db*(u(i,j,k) &
			       + u(i,j,k)**2.0_db - 3.0_db*u(i,j,k)*v(i,j,k) + (-1.0_db &
			       + v(i,j,k))*v(i,j,k)) - 3.0_db*w(i,j,k)**2.0_db)/108.0_db
#else
!9
			      feq=(2.0_db*rho(i,j,k) + 6.0_db*(u(i,j,k) &
			       + (-1.0_db + v(i,j,k))*v(i,j,k) + 3.0_db*u(i,j,k)*(-1.0_db &
			       + v(i,j,k))*v(i,j,k) + u(i,j,k)**2.0_db*(1.0_db + 3.0_db*(-1.0_db &
			       + v(i,j,k))*v(i,j,k))) - 3.0_db*(1.0_db &
			       + 3.0_db*u(i,j,k)*(1.0_db + u(i,j,k)))*(1.0_db + 3.0_db*(-1.0_db &
			       + v(i,j,k))*v(i,j,k))*w(i,j,k)**2.0_db)/108.0_db
!9
#endif
				  fneq1=3.0_db*(pxx(i,j,k) - 3.0_db*pxy(i,j,k) + pyy(i,j,k)) &
				   - (3.0_db*pzz(i,j,k))/2.0_db


				  F_discr=(forcey + 3.0_db*forcey*u(i,j,k) - 2.0_db*forcey*v(i,j,k) &
				   + forcex*(-1.0_db - 2.0_db*u(i,j,k) + 3.0_db*v(i,j,k)) &
				   + forcez*w(i,j,k))/(-18.0_db*rhophi_loc)
                  f(i+1,j-1,k,9)=feq + (1.0_db-omega_loc)*fneq1*p2 + 0.5_db*(F_discr)


!!!!!!!!!!!!!!!!!!!!!!!!!!10
#ifdef SECOND_ORDER
			      feq=(2.0_db*rho(i,j,k) + 6.0_db*((-1.0_db &
			       + u(i,j,k))*u(i,j,k) + v(i,j,k) - 3.0_db*u(i,j,k)*v(i,j,k) &
			       + v(i,j,k)**2.0_db) - 3.0_db*w(i,j,k)**2.0_db)/108.0_db
#else
!10
			      feq=(2.0_db*rho(i,j,k) + 6.0_db*((-1.0_db &
			       + u(i,j,k))*u(i,j,k) + v(i,j,k) + 3.0_db*(-1.0_db &
			       + u(i,j,k))*u(i,j,k)*v(i,j,k) + (1.0_db + 3.0_db*(-1.0_db &
			       + u(i,j,k))*u(i,j,k))*v(i,j,k)**2.0_db) - 3.0_db*(1.0_db &
			       + 3.0_db*(-1.0_db + u(i,j,k))*u(i,j,k))*(1.0_db + 3.0_db*v(i,j,k)*(1.0_db &
			       + v(i,j,k)))*w(i,j,k)**2.0_db)/108.0_db
!10
#endif
				  fneq1=3.0_db*(pxx(i,j,k) - 3.0_db*pxy(i,j,k) + pyy(i,j,k)) &
				   - (3.0_db*pzz(i,j,k))/2.0_db


				  F_discr=(forcex - forcey - 2.0_db*forcex*u(i,j,k) &
				   + 3.0_db*forcey*u(i,j,k) + 3.0_db*forcex*v(i,j,k) &
				   - 2.0_db*forcey*v(i,j,k) + forcez*w(i,j,k))/(-18.0_db*rhophi_loc)
                  f(i-1,j+1,k,10)=feq + (1.0_db-omega_loc)*fneq1*p2 + 0.5_db*(F_discr)


!!!!!!!!!!!!!!!!!!!!!!!!!!11
#ifdef SECOND_ORDER
			      feq=(2.0_db*rho(i,j,k) - 3.0_db*u(i,j,k)**2.0_db &
			       + 6.0_db*(v(i,j,k) + v(i,j,k)**2.0_db + w(i,j,k) &
			       + 3.0_db*v(i,j,k)*w(i,j,k) + w(i,j,k)**2.0_db))/108.0_db
#else
!11
			      feq=(2.0_db*rho(i,j,k) - 3.0_db*u(i,j,k)**2.0_db*(1.0_db &
			       + 3.0_db*v(i,j,k)*(1.0_db + v(i,j,k)))*(1.0_db + 3.0_db*w(i,j,k)*(1.0_db &
			       + w(i,j,k))) + 6.0_db*(v(i,j,k) + w(i,j,k) + w(i,j,k)**2.0_db &
			       + 3.0_db*v(i,j,k)*w(i,j,k)*(1.0_db + w(i,j,k)) + v(i,j,k)**2.0_db*(1.0_db &
			       + 3.0_db*w(i,j,k)*(1.0_db + w(i,j,k)))))/108.0_db
!11
#endif
				  fneq1=(-3.0_db*pxx(i,j,k))/2.0_db + 3.0_db*(pyy(i,j,k) &
				   + 3.0_db*pyz(i,j,k) + pzz(i,j,k))


				  F_discr=(forcey + forcez - forcex*u(i,j,k) + 2.0_db*forcey*v(i,j,k) &
				   + 3.0_db*forcez*v(i,j,k) + 3.0_db*forcey*w(i,j,k) &
				   + 2.0_db*forcez*w(i,j,k))/(18.0_db*rhophi_loc)
                  f(i,j+1,k+1,11)=feq + (1.0_db-omega_loc)*fneq1*p2 + 0.5_db*(F_discr)


!!!!!!!!!!!!!!!!!!!!!!!!!!12
#ifdef SECOND_ORDER
			      feq=(2.0_db*rho(i,j,k) - 3.0_db*u(i,j,k)**2.0_db &
			       + 6.0_db*(v(i,j,k)**2.0_db + (-1.0_db &
			       + w(i,j,k))*w(i,j,k) + v(i,j,k)*(-1.0_db + 3.0_db*w(i,j,k))))/108.0_db
#else
!12
			      feq=(2.0_db*rho(i,j,k) + 3.0_db*(2.0_db*(-1.0_db &
			       + w(i,j,k))*w(i,j,k) + v(i,j,k)*(-2.0_db - 6.0_db*(-1.0_db &
			       + w(i,j,k))*w(i,j,k)) - u(i,j,k)**2.0_db*(1.0_db + 3.0_db*(-1.0_db &
			       + v(i,j,k))*v(i,j,k))*(1.0_db + 3.0_db*(-1.0_db + w(i,j,k))*w(i,j,k)) &
			       + v(i,j,k)**2.0_db*(2.0_db + 6.0_db*(-1.0_db + w(i,j,k))*w(i,j,k))))/108.0_db
!12
#endif
				  fneq1=(-3.0_db*pxx(i,j,k))/2.0_db + 3.0_db*(pyy(i,j,k) &
				   + 3.0_db*pyz(i,j,k) + pzz(i,j,k))


				  F_discr=(forcey + forcez + forcex*u(i,j,k) - 2.0_db*forcey*v(i,j,k) &
				   - 3.0_db*forcez*v(i,j,k) - 3.0_db*forcey*w(i,j,k) &
				   - 2.0_db*forcez*w(i,j,k))/(-18.0_db*rhophi_loc)
                  f(i,j-1,k-1,12)=feq + (1.0_db-omega_loc)*fneq1*p2 + 0.5_db*(F_discr)


!!!!!!!!!!!!!!!!!!!!!!!!!!13
#ifdef SECOND_ORDER
			      feq=(2.0_db*rho(i,j,k) - 3.0_db*u(i,j,k)**2.0_db &
			       + 6.0_db*(v(i,j,k) + v(i,j,k)**2.0_db &
			       - 3.0_db*v(i,j,k)*w(i,j,k) + (-1.0_db + w(i,j,k))*w(i,j,k)))/108.0_db
#else
!13
			      feq=(2.0_db*rho(i,j,k) - 3.0_db*u(i,j,k)**2.0_db*(1.0_db &
			       + 3.0_db*v(i,j,k)*(1.0_db + v(i,j,k)))*(1.0_db + 3.0_db*(-1.0_db &
			       + w(i,j,k))*w(i,j,k)) + 6.0_db*(v(i,j,k) + (-1.0_db &
			       + w(i,j,k))*w(i,j,k) + 3.0_db*v(i,j,k)*(-1.0_db + w(i,j,k))*w(i,j,k) &
			       + v(i,j,k)**2.0_db*(1.0_db + 3.0_db*(-1.0_db + w(i,j,k))*w(i,j,k))))/108.0_db
!13
#endif
				  fneq1=(-3.0_db*pxx(i,j,k))/2.0_db + 3.0_db*(pyy(i,j,k) &
				   - 3.0_db*pyz(i,j,k) + pzz(i,j,k))


				  F_discr=(forcez + forcex*u(i,j,k) + 3.0_db*forcez*v(i,j,k) &
				   - 2.0_db*forcez*w(i,j,k) + forcey*(-1.0_db - 2.0_db*v(i,j,k) &
				   + 3.0_db*w(i,j,k)))/(-18.0_db*rhophi_loc)
                  f(i,j+1,k-1,13)=feq + (1.0_db-omega_loc)*fneq1*p2 + 0.5_db*(F_discr)



!!!!!!!!!!!!!!!!!!!!!!!!!!14
#ifdef SECOND_ORDER
			      feq=(2.0_db*rho(i,j,k) - 3.0_db*u(i,j,k)**2.0_db &
			       + 6.0_db*((-1.0_db + v(i,j,k))*v(i,j,k) + w(i,j,k) &
			       - 3.0_db*v(i,j,k)*w(i,j,k) + w(i,j,k)**2.0_db))/108.0_db
#else
!14
			      feq=(2.0_db*rho(i,j,k) + 6.0_db*((-1.0_db &
			       + v(i,j,k))*v(i,j,k) + w(i,j,k) + 3.0_db*(-1.0_db &
			       + v(i,j,k))*v(i,j,k)*w(i,j,k) + (1.0_db + 3.0_db*(-1.0_db &
			       + v(i,j,k))*v(i,j,k))*w(i,j,k)**2.0_db) &
			       - 3.0_db*u(i,j,k)**2.0_db*(1.0_db + 3.0_db*(-1.0_db &
			       + v(i,j,k))*v(i,j,k))*(1.0_db + 3.0_db*w(i,j,k)*(1.0_db &
			       + w(i,j,k))))/108.0_db
!14
#endif
				  fneq1=(-3.0_db*pxx(i,j,k))/2.0_db + 3.0_db*(pyy(i,j,k) &
				   - 3.0_db*pyz(i,j,k) + pzz(i,j,k))


				  F_discr=(forcey - forcez + forcex*u(i,j,k) - 2.0_db*forcey*v(i,j,k) &
				   + 3.0_db*forcez*v(i,j,k) + 3.0_db*forcey*w(i,j,k) &
				   - 2.0_db*forcez*w(i,j,k))/(-18.0_db*rhophi_loc)
                  f(i,j-1,k+1,14)=feq + (1.0_db-omega_loc)*fneq1*p2 + 0.5_db*(F_discr)


!!!!!!!!!!!!!!!!!!!!!!!!!!15
#ifdef SECOND_ORDER
			      feq=(2.0_db*rho(i,j,k) - 3.0_db*v(i,j,k)**2.0_db &
			       + 6.0_db*w(i,j,k) + 6.0_db*(u(i,j,k) &
			       + u(i,j,k)**2.0_db + 3.0_db*u(i,j,k)*w(i,j,k) + w(i,j,k)**2.0_db))/108.0_db
#else
!15
			      feq=(2.0_db*rho(i,j,k) + 3.0_db*(2.0_db*w(i,j,k)*(1.0_db &
			       + w(i,j,k)) - v(i,j,k)**2.0_db*(1.0_db + 3.0_db*w(i,j,k)*(1.0_db &
			       + w(i,j,k))) - u(i,j,k)*(-2.0_db + 3.0_db*v(i,j,k)**2.0_db)*(1.0_db &
			       + 3.0_db*w(i,j,k)*(1.0_db + w(i,j,k))) - u(i,j,k)**2.0_db*(-2.0_db &
			       + 3.0_db*v(i,j,k)**2.0_db)*(1.0_db + 3.0_db*w(i,j,k)*(1.0_db &
			       + w(i,j,k)))))/108.0_db
!15
#endif
				  fneq1=3.0_db*pxx(i,j,k) + 9.0_db*pxz(i,j,k) &
				   - (3.0_db*pyy(i,j,k))/2.0_db + 3.0_db*pzz(i,j,k)


				  F_discr=(forcex + forcez + 2.0_db*forcex*u(i,j,k) &
				   + 3.0_db*forcez*u(i,j,k) - forcey*v(i,j,k) + 3.0_db*forcex*w(i,j,k) &
				   + 2.0_db*forcez*w(i,j,k))/(18.0_db*rhophi_loc)
                  f(i+1,j,k+1,15)=feq + (1.0_db-omega_loc)*fneq1*p2 + 0.5_db*(F_discr)



!!!!!!!!!!!!!!!!!!!!!!!!!!16
#ifdef SECOND_ORDER
			      feq=(2.0_db*rho(i,j,k) + 3.0_db*(2.0_db*u(i,j,k)**2.0_db &
			       - v(i,j,k)**2.0_db + 2.0_db*(-1.0_db + w(i,j,k))*w(i,j,k) &
			       + u(i,j,k)*(-2.0_db + 6.0_db*w(i,j,k))))/108.0_db
#else
!16
			      feq=(2.0_db*rho(i,j,k) + 3.0_db*(2.0_db*(-1.0_db &
			       + w(i,j,k))*w(i,j,k) + v(i,j,k)**2.0_db*(-1.0_db - 3.0_db*(-1.0_db &
			       + w(i,j,k))*w(i,j,k)) + u(i,j,k)*(-2.0_db + 3.0_db*v(i,j,k)**2.0_db)*(1.0_db &
			       + 3.0_db*(-1.0_db + w(i,j,k))*w(i,j,k)) - u(i,j,k)**2.0_db*(-2.0_db &
			       + 3.0_db*v(i,j,k)**2.0_db)*(1.0_db + 3.0_db*(-1.0_db &
			       + w(i,j,k))*w(i,j,k))))/108.0_db
!16
#endif
				  fneq1=3.0_db*pxx(i,j,k) + 9.0_db*pxz(i,j,k) &
				   - (3.0_db*pyy(i,j,k))/2.0_db + 3.0_db*pzz(i,j,k)


				  F_discr=(forcex + forcez - 2.0_db*forcex*u(i,j,k) &
				   - 3.0_db*forcez*u(i,j,k) + forcey*v(i,j,k) - 3.0_db*forcex*w(i,j,k) &
				   - 2.0_db*forcez*w(i,j,k))/(-18.0_db*rhophi_loc)
                  f(i-1,j,k-1,16)=feq + (1.0_db-omega_loc)*fneq1*p2 + 0.5_db*(F_discr)



!!!!!!!!!!!!!!!!!!!!!!!!!!17
#ifdef SECOND_ORDER
			      feq=(2.0_db*rho(i,j,k) + 3.0_db*(2.0_db*u(i,j,k)**2.0_db &
			       - v(i,j,k)**2.0_db + 2.0_db*w(i,j,k)*(1.0_db + w(i,j,k)) &
			       - 2.0_db*u(i,j,k)*(1.0_db + 3.0_db*w(i,j,k))))/108.0_db
#else
!17
			      feq=(2.0_db*rho(i,j,k) + 3.0_db*(2.0_db*w(i,j,k)*(1.0_db &
			       + w(i,j,k)) - v(i,j,k)**2.0_db*(1.0_db + 3.0_db*w(i,j,k)*(1.0_db &
			       + w(i,j,k))) + u(i,j,k)*(-2.0_db + 3.0_db*v(i,j,k)**2.0_db)*(1.0_db &
			       + 3.0_db*w(i,j,k)*(1.0_db + w(i,j,k))) - u(i,j,k)**2.0_db*(-2.0_db &
			       + 3.0_db*v(i,j,k)**2.0_db)*(1.0_db + 3.0_db*w(i,j,k)*(1.0_db &
			       + w(i,j,k)))))/108.0_db
!17
#endif
				  fneq1=3.0_db*pxx(i,j,k) - 9.0_db*pxz(i,j,k) &
				   - (3.0_db*pyy(i,j,k))/2.0_db + 3.0_db*pzz(i,j,k)


				  F_discr=(forcex - forcez - 2.0_db*forcex*u(i,j,k) &
				   + 3.0_db*forcez*u(i,j,k) + forcey*v(i,j,k) + 3.0_db*forcex*w(i,j,k) &
				   - 2.0_db*forcez*w(i,j,k))/(-18.0_db*rhophi_loc)
                  f(i-1,j,k+1,17)=feq + (1.0_db-omega_loc)*fneq1*p2 + 0.5_db*(F_discr)



!!!!!!!!!!!!!!!!!!!!!!!!!!18
#ifdef SECOND_ORDER
			      feq=(2.0_db*rho(i,j,k) - 3.0_db*v(i,j,k)**2.0_db &
			       - 6.0_db*w(i,j,k) + 6.0_db*(u(i,j,k) &
			       + u(i,j,k)**2.0_db - 3.0_db*u(i,j,k)*w(i,j,k) + w(i,j,k)**2.0_db))/108.0_db
#else
!18
			      feq=(2.0_db*rho(i,j,k) + 3.0_db*(2.0_db*(-1.0_db &
			       + w(i,j,k))*w(i,j,k) + v(i,j,k)**2.0_db*(-1.0_db &
			       - 3.0_db*(-1.0_db + w(i,j,k))*w(i,j,k)) - u(i,j,k)*(-2.0_db &
			       + 3.0_db*v(i,j,k)**2.0_db)*(1.0_db + 3.0_db*(-1.0_db &
			       + w(i,j,k))*w(i,j,k)) - u(i,j,k)**2.0_db*(-2.0_db &
			       + 3.0_db*v(i,j,k)**2.0_db)*(1.0_db + 3.0_db*(-1.0_db &
			       + w(i,j,k))*w(i,j,k))))/108.0_db
!18
#endif
				  fneq1=3.0_db*pxx(i,j,k) - 9.0_db*pxz(i,j,k) &
				   - (3.0_db*pyy(i,j,k))/2.0_db + 3.0_db*pzz(i,j,k)


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
#else
!19
			      feq=(rho(i,j,k) + 3.0_db*(v(i,j,k) + v(i,j,k)**2.0_db &
			       + w(i,j,k) + 3.0_db*v(i,j,k)*(1.0_db + v(i,j,k))*w(i,j,k) + (1.0_db &
			       + 3.0_db*v(i,j,k)*(1.0_db + v(i,j,k)))*w(i,j,k)**2.0_db + u(i,j,k)*(1.0_db &
			       + 3.0_db*v(i,j,k)*(1.0_db + v(i,j,k)))*(1.0_db + 3.0_db*w(i,j,k)*(1.0_db &
			       + w(i,j,k))) + u(i,j,k)**2.0_db*(1.0_db + 3.0_db*v(i,j,k)*(1.0_db &
			       + v(i,j,k)))*(1.0_db + 3.0_db*w(i,j,k)*(1.0_db + w(i,j,k)))))/216.0_db
!19
#endif
				  fneq1=3.0_db*(pxx(i,j,k) + 3.0_db*pxy(i,j,k) &
				   + 3.0_db*pxz(i,j,k) + pyy(i,j,k) + 3.0_db*pyz(i,j,k) + pzz(i,j,k))


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
#else
!20
			      feq=(rho(i,j,k) + 3.0_db*((-1.0_db + v(i,j,k))*v(i,j,k) &
			       + (-1.0_db - 3.0_db*(-1.0_db + v(i,j,k))*v(i,j,k))*w(i,j,k) &
			       + (1.0_db + 3.0_db*(-1.0_db + v(i,j,k))*v(i,j,k))*w(i,j,k)**2.0_db &
			       - u(i,j,k)*(1.0_db + 3.0_db*(-1.0_db + v(i,j,k))*v(i,j,k))*(1.0_db &
			       + 3.0_db*(-1.0_db + w(i,j,k))*w(i,j,k)) + u(i,j,k)**2.0_db*(1.0_db &
			       + 3.0_db*(-1.0_db + v(i,j,k))*v(i,j,k))*(1.0_db + 3.0_db*(-1.0_db &
			       + w(i,j,k))*w(i,j,k))))/216.0_db
!20
#endif
				  fneq1=3.0_db*(pxx(i,j,k) + 3.0_db*pxy(i,j,k) &
				   + 3.0_db*pxz(i,j,k) + pyy(i,j,k) + 3.0_db*pyz(i,j,k) + pzz(i,j,k))


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
#else
!21
			      feq=(rho(i,j,k) + 3.0_db*((-1.0_db + v(i,j,k))*v(i,j,k) &
			       + w(i,j,k) + 3.0_db*(-1.0_db + v(i,j,k))*v(i,j,k)*w(i,j,k) &
			       + (1.0_db + 3.0_db*(-1.0_db + v(i,j,k))*v(i,j,k))*w(i,j,k)**2.0_db &
			       + u(i,j,k)*(1.0_db + 3.0_db*(-1.0_db + v(i,j,k))*v(i,j,k))*(1.0_db &
			       + 3.0_db*w(i,j,k)*(1.0_db + w(i,j,k))) + u(i,j,k)**2.0_db*(1.0_db &
			       + 3.0_db*(-1.0_db + v(i,j,k))*v(i,j,k))*(1.0_db &
			       + 3.0_db*w(i,j,k)*(1.0_db + w(i,j,k)))))/216.0_db
!21
#endif 
				  fneq1=3.0_db*(pxx(i,j,k) - 3.0_db*pxy(i,j,k) &
				   + 3.0_db*pxz(i,j,k) + pyy(i,j,k) - 3.0_db*pyz(i,j,k) + pzz(i,j,k))


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
#else
!22
			      feq=(rho(i,j,k) + 3.0_db*(v(i,j,k) + v(i,j,k)**2.0_db &
			       - w(i,j,k) - 3.0_db*v(i,j,k)*(1.0_db + v(i,j,k))*w(i,j,k) + (1.0_db &
			       + 3.0_db*v(i,j,k)*(1.0_db + v(i,j,k)))*w(i,j,k)**2.0_db - u(i,j,k)*(1.0_db &
			       + 3.0_db*v(i,j,k)*(1.0_db + v(i,j,k)))*(1.0_db + 3.0_db*(-1.0_db &
			       + w(i,j,k))*w(i,j,k)) + u(i,j,k)**2.0_db*(1.0_db + 3.0_db*v(i,j,k)*(1.0_db &
			       + v(i,j,k)))*(1.0_db + 3.0_db*(-1.0_db + w(i,j,k))*w(i,j,k))))/216.0_db
!22
#endif
				  fneq1=3.0_db*(pxx(i,j,k) - 3.0_db*pxy(i,j,k) &
				   + 3.0_db*pxz(i,j,k) + pyy(i,j,k) - 3.0_db*pyz(i,j,k) + pzz(i,j,k))


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
#else
!23
			      feq=(rho(i,j,k) + 3.0_db*((-1.0_db + v(i,j,k))*v(i,j,k) &
			       + w(i,j,k) + 3.0_db*(-1.0_db + v(i,j,k))*v(i,j,k)*w(i,j,k) + (1.0_db &
			       + 3.0_db*(-1.0_db + v(i,j,k))*v(i,j,k))*w(i,j,k)**2.0_db - u(i,j,k)*(1.0_db &
			       + 3.0_db*(-1.0_db + v(i,j,k))*v(i,j,k))*(1.0_db + 3.0_db*w(i,j,k)*(1.0_db &
			       + w(i,j,k))) + u(i,j,k)**2.0_db*(1.0_db + 3.0_db*(-1.0_db + v(i,j,k))*v(i,j,k))*(1.0_db &
			       + 3.0_db*w(i,j,k)*(1.0_db + w(i,j,k)))))/216.0_db
!23
#endif
				  fneq1=3.0_db*(pxx(i,j,k) + 3.0_db*pxy(i,j,k) &
				   - 3.0_db*pxz(i,j,k) + pyy(i,j,k) - 3.0_db*pyz(i,j,k) + pzz(i,j,k))


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
#else
!24
			      feq=(rho(i,j,k) + 3.0_db*(v(i,j,k) + v(i,j,k)**2.0_db &
			       - w(i,j,k) - 3.0_db*v(i,j,k)*(1.0_db + v(i,j,k))*w(i,j,k) + (1.0_db &
			       + 3.0_db*v(i,j,k)*(1.0_db + v(i,j,k)))*w(i,j,k)**2.0_db + u(i,j,k)*(1.0_db &
			       + 3.0_db*v(i,j,k)*(1.0_db + v(i,j,k)))*(1.0_db + 3.0_db*(-1.0_db &
			       + w(i,j,k))*w(i,j,k)) + u(i,j,k)**2.0_db*(1.0_db + 3.0_db*v(i,j,k)*(1.0_db &
			       + v(i,j,k)))*(1.0_db + 3.0_db*(-1.0_db + w(i,j,k))*w(i,j,k))))/216.0_db
!24
#endif
				  fneq1=3.0_db*(pxx(i,j,k) + 3.0_db*pxy(i,j,k) &
				   - 3.0_db*pxz(i,j,k) + pyy(i,j,k) - 3.0_db*pyz(i,j,k) + pzz(i,j,k))


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
#else
!25
			      feq=(rho(i,j,k) + 3.0_db*((-1.0_db + v(i,j,k))*v(i,j,k) &
			       + (-1.0_db - 3.0_db*(-1.0_db + v(i,j,k))*v(i,j,k))*w(i,j,k) + (1.0_db &
			       + 3.0_db*(-1.0_db + v(i,j,k))*v(i,j,k))*w(i,j,k)**2.0_db + u(i,j,k)*(1.0_db &
			       + 3.0_db*(-1.0_db + v(i,j,k))*v(i,j,k))*(1.0_db + 3.0_db*(-1.0_db &
			       + w(i,j,k))*w(i,j,k)) + u(i,j,k)**2.0_db*(1.0_db + 3.0_db*(-1.0_db &
			       + v(i,j,k))*v(i,j,k))*(1.0_db + 3.0_db*(-1.0_db + w(i,j,k))*w(i,j,k))))/216.0_db
!25
#endif
				  fneq1=3.0_db*(pxx(i,j,k) - 3.0_db*pxy(i,j,k) - 3.0_db*pxz(i,j,k) &
				   + pyy(i,j,k) + 3.0_db*pyz(i,j,k) + pzz(i,j,k))


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
#else
!26
			      feq=(rho(i,j,k) + 3.0_db*(v(i,j,k) + v(i,j,k)**2.0_db &
			       + w(i,j,k) + 3.0_db*v(i,j,k)*(1.0_db + v(i,j,k))*w(i,j,k) + (1.0_db &
			       + 3.0_db*v(i,j,k)*(1.0_db + v(i,j,k)))*w(i,j,k)**2.0_db - u(i,j,k)*(1.0_db &
			       + 3.0_db*v(i,j,k)*(1.0_db + v(i,j,k)))*(1.0_db + 3.0_db*w(i,j,k)*(1.0_db &
			       + w(i,j,k))) + u(i,j,k)**2.0_db*(1.0_db + 3.0_db*v(i,j,k)*(1.0_db &
			       + v(i,j,k)))*(1.0_db + 3.0_db*w(i,j,k)*(1.0_db + w(i,j,k)))))/216.0_db
!26
#endif
				  fneq1=3.0_db*(pxx(i,j,k) - 3.0_db*pxy(i,j,k) &
				   - 3.0_db*pxz(i,j,k) + pyy(i,j,k) + 3.0_db*pyz(i,j,k) + pzz(i,j,k))


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

                  gradfix=normx(i,j,k)*modgrad(i,j,k)
				  gradfiy=normy(i,j,k)*modgrad(i,j,k)
		          gradfiz=normz(i,j,k)*modgrad(i,j,k)
                  selphi(i,j,k,flop) = selphi(i,j,k,flip) &
                   - u(i,j,k)*0.5_db*(gradfix) - v(i,j,k)*0.5_db*(gradfiy) &
                   - w(i,j,k)*0.5_db*(gradfiz) + tau_diff*lap_phi(i,j,k) + feq 
#endif	

#ifdef MONOD
			    S_mono = mu_max * selphi(i,j,k,flip)/(Ks + selphi(i,j,k,flip)) * selphi(i,j,k,flip) * (1.0_db - selphi(i,j,k,flip))
				selphi(i,j,k,flop)=selphi(i,j,k,flop) + S_mono
				!selphi(i,j,k,flop) = min(1.0_db, max(0.0_db, selphi(i,j,k,flop)))		 
#endif
               

   endsubroutine fused_LB_kernel   
   



endmodule
