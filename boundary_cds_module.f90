#include "defines.h"
module bcs3D

   use vars
   use mpi_template , only: coords,pbc_x,pbc_y,pbc_z,myoffset,myrank, &
    sum_world_float,sum_world_int
!   !$if _OPENACC
!   use openacc
!   !$endif
   implicit none

contains
   !***************************************************
   subroutine bcs_mesoscopic_all

      implicit none

      integer :: subchords(3)
      integer :: ii,jj,kk,l,lopp
	  real(kind=db) :: feq, fneq1,utmp,vtmp,wtmp,rhophi_loc
#ifdef TWOCOMPONENT	  
	  real(kind=db) :: phitemp
#endif
#if defined(PHASE_CHANGE) || defined(INTERNAL_OBSTACLES)
	  real(kind=db) :: visc_loc,omega_loc,tau_loc
#endif
#if defined(INTERNAL_OBSTACLES)
	  real(kind=db) :: F_discr,wet_R,rhotemp,phiavg, wet_thresh_low, wet_thresh_high,grad_thresh,phi_adj,weight,correc
	  real(kind=db) :: phi_fluid,gradfix,gradfiy,grad_parallel,theta_rad,cot_theta,phi_ghost,dphi_dz,gradfiz
	  integer :: conter
	  logical :: found
#endif




#if defined(INTERNAL_OBSTACLES)

#ifdef ACCNOKERNELS
		!$acc parallel loop collapse(3) present(pxx,pyy,pzz,pxy,pxz,pyz,rho,u,v,w,fux,fvy,fwz &
#ifdef TWOCOMPONENT
		!$acc& ,selphi &
#endif
#ifdef DENSRATIO
		!$acc& ,rhophi &
#endif
		!$acc& ) private(i,j,k,l,F_discr,fneq1,feq,omega_loc,udotc,uu,gi,gj,gk &
#ifdef TWOCOMPONENT
		!$acc& ,tau_loc,rhophi_loc,visc_loc &
#endif
		!$acc& )
#else
		!$acc kernels present(pxx,pyy,pzz,pxy,pxz,pyz,rho,u,v,w,fux,fvy,fwz &
#ifdef TWOCOMPONENT
		!$acc& ,selphi &
#endif
#ifdef DENSRATIO
		!$acc& ,rhophi &
#endif
		!$acc& )
		!$acc loop collapse(3) private(i,j,k,l,F_discr,fneq1,feq,omega_loc,udotc,uu,gi,gj,gk &
#ifdef TWOCOMPONENT
		!$acc& ,tau_loc,rhophi_loc,visc_loc &
#endif
		!$acc& )
#endif
		do k=1,nz
			do j=1,ny
				do i=1,nx

					    if(isfluid(i,j,k).ne.-1) cycle

						utmp=0.0_db!u(i,j,k)
						vtmp=0.0_db!v(i,j,k)	
						wtmp=0.0_db
						
#ifdef DENSRATIO
						rhophi_loc = rhophi(i,j,k)
#else
				        rhophi_loc = ONE 
#endif

#ifdef TWOCOMPONENT

						visc_loc=(rho_r*visc1*selphi(i,j,k,flip)+(ONE-selphi(i,j,k,flip))*visc2*rho_b)/rhophi_loc  

						
						tau_loc=(visc_loc/cssq + HALF) !è una tau
						
						omega_loc=ONE/tau_loc !è una omega
						
#else
						omega_loc=omega
#endif				        
				        
						!$acc loop seq
						do l=1,nlinks
						  lopp=opp(l)
						  ii=i+ex(lopp)
						  jj=j+ey(lopp)
						  kk=k+ez(lopp)
						  if(isfluid(ii,jj,kk).ne.0) cycle 
						  !w(i,j,k)=w(i,j,k+ez(l))
						  uu=HALF*(utmp*utmp + vtmp*vtmp+ wtmp*wtmp)/cssq
						  udotc=(utmp*dex(l) + vtmp*dey(l)+ wtmp*dez(l))/cssq
						  feq=p(l)*(rho(i+ex(l),j+ey(l),k+ez(l)) + udotc+ HALF*udotc*udotc - uu)
						  fneq1=(HALF/(cssq*cssq))*( (dex(l)*dex(l)-cssq)*pxx(i,j,k) &
						   + (dey(l)*dey(l)-cssq)*pyy(i,j,k) + (dez(l)*dez(l)-cssq)*pzz(i,j,k) &
					       + TWO*(dex(l)*dey(l))*pxy(i,j,k) + TWO*(dex(l)*dez(l))*pxz(i,j,k) &
						   + TWO*(dey(l)*dez(l))*pyz(i,j,k))
				          ! F_discr = p(l)*(((dex(l) - utmp) + udotc * dex(l))*fux(i,j,k) &
				           ! + ((dey(l) - vtmp) + udotc * dey(l))*fvy(i,j,k) &
				           ! + ((dez(l) - wtmp) + udotc * dez(l))*fwz(i,j,k))/cssq
						  f(i,j,k,l)=feq + (ONE-omega_loc)*p(l)*fneq1 !+ HALF*(F_discr)
						  
						enddo

				enddo
			enddo
		enddo
#ifdef ACCNOKERNELS
      !$acc end parallel loop
#else
      !$acc end kernels
#endif

!*****************************************

	   if(pbc_x.eq.0)then
#ifdef ACCNOKERNELS
		!$acc parallel loop collapse(2) independent present(selphi &
		!$acc& ) private(i,j,k,l,gi,gj,gk)
#else
		!$acc kernels present(selphi &
		!$acc& )
		!$acc loop collapse(2) independent private(i,j,k,l,gi,gj,gk)
#endif	    
	      do k=1,nz
			do j=1,ny
			  gj=ny*coords(2)+j
			  gk=nz*coords(3)+k
			  if(gj>1 .and. gj<ly .and. gk>1 .and. gk<lz)then
	          i=1
	          gi=nx*coords(1)+i
	          if(gi.eq.1)then
#ifdef DENSRATIO
			    rhophi_loc = rhophi(i+1,j,k)
#else
		        rhophi_loc = ONE 
#endif

#ifdef TWOCOMPONENT

			    visc_loc=(rho_r*visc1*selphi(i+1,j,k,flip)+(ONE-selphi(i+1,j,k,flip))*visc2*rho_b)/rhophi_loc  

						
			    tau_loc=(visc_loc/cssq + HALF) !è una tau
						
			    omega_loc=ONE/tau_loc !è una omega
						
#else
			    omega_loc=omega
#endif					
				utmp=0.0_db!u(i+1,j,k)
				vtmp=0.0_db!v(i,j,k)	
				wtmp=0.0_db
				!$acc loop seq
				do l=1,nlinks
				  ! lopp=opp(l)
				   ii=i+ex(l)
				   jj=j+ey(l)
				   kk=k+ez(l)
				  
				  if(isfluid(ii,jj,kk).ne.-1) cycle 
				  rhotemp=rho(ii,jj,kk)
				  !w(i,j,k)=w(i,j,k+ez(l))
				  uu=HALF*(utmp*utmp + vtmp*vtmp+ wtmp*wtmp)/cssq
				  udotc=(utmp*dex(l) + vtmp*dey(l)+ wtmp*dez(l))/cssq
				  feq=p(l)*(rhotemp + udotc+ HALF*udotc*udotc - uu)
				  fneq1=(HALF/(cssq*cssq))*( (dex(l)*dex(l)-cssq)*pxx(ii,jj,kk) &
				   + (dey(l)*dey(l)-cssq)*pyy(ii,jj,kk) + (dez(l)*dez(l)-cssq)*pzz(ii,jj,kk) &
				   + TWO*(dex(l)*dey(l))*pxy(ii,jj,kk) + TWO*(dex(l)*dez(l))*pxz(ii,jj,kk) &
				   + TWO*(dey(l)*dez(l))*pyz(ii,jj,kk))
				  ! F_discr = p(l)*(((dex(l) - utmp) + udotc * dex(l))*fux(i,j,k) &
				   ! + ((dey(l) - vtmp) + udotc * dey(l))*fvy(i,j,k) &
				   ! + ((dez(l) - wtmp) + udotc * dez(l))*fwz(i,j,k))/cssq
				  f(ii,jj,kk,l)=feq + (ONE-omega_loc)*p(l)*fneq1 !+ HALF*(F_discr)
				  
				enddo
			  endif
			  i=nx
	          gi=nx*coords(1)+i
	          if(gi.eq.lx)then
#ifdef DENSRATIO
			    rhophi_loc = rhophi(i-1,j,k)
#else
		        rhophi_loc = ONE 
#endif

#ifdef TWOCOMPONENT

			    visc_loc=(rho_r*visc1*selphi(i-1,j,k,flip)+(ONE-selphi(i-1,j,k,flip))*visc2*rho_b)/rhophi_loc  

						
			    tau_loc=(visc_loc/cssq + HALF) !è una tau
						
			    omega_loc=ONE/tau_loc !è una omega
						
#else
			    omega_loc=omega
#endif	
				utmp=0.0_db!u(i-1,j,k)
				vtmp=0.0_db!v(i,j,k)	
				wtmp=0.0_db
				!$acc loop seq
				do l=1,nlinks
				  ! lopp=opp(l)
				   ii=i+ex(l)
				   jj=j+ey(l)
				   kk=k+ez(l)
				  
				  if(isfluid(ii,jj,kk).ne.-1) cycle 
				  rhotemp=rho(ii,jj,kk)
				  !w(i,j,k)=w(i,j,k+ez(l))
				  uu=HALF*(utmp*utmp + vtmp*vtmp+ wtmp*wtmp)/cssq
				  udotc=(utmp*dex(l) + vtmp*dey(l)+ wtmp*dez(l))/cssq
				  feq=p(l)*(rhotemp + udotc+ HALF*udotc*udotc - uu)
				  fneq1=(HALF/(cssq*cssq))*( (dex(l)*dex(l)-cssq)*pxx(ii,jj,kk) &
				   + (dey(l)*dey(l)-cssq)*pyy(ii,jj,kk) + (dez(l)*dez(l)-cssq)*pzz(ii,jj,kk) &
				   + TWO*(dex(l)*dey(l))*pxy(ii,jj,kk) + TWO*(dex(l)*dez(l))*pxz(ii,jj,kk) &
				   + TWO*(dey(l)*dez(l))*pyz(ii,jj,kk))
				  ! F_discr = p(l)*(((dex(l) - utmp) + udotc * dex(l))*fux(i,j,k) &
				   ! + ((dey(l) - vtmp) + udotc * dey(l))*fvy(i,j,k) &
				   ! + ((dez(l) - wtmp) + udotc * dez(l))*fwz(i,j,k))/cssq
				  f(ii,jj,kk,l)=feq + (ONE-omega_loc)*p(l)*fneq1 !+ HALF*(F_discr)
				  
				enddo
			  endif
			  endif
	        enddo
	      enddo
#ifdef ACCNOKERNELS
      !$acc end parallel loop
#else
      !$acc end kernels
#endif
	    endif	
	    
	    if(pbc_y.eq.0)then
#ifdef ACCNOKERNELS
		!$acc parallel loop collapse(2) independent present(selphi)
#else
		!$acc kernels present(selphi &
		!$acc& )
		!$acc loop collapse(2) independent private(i,j,k,l,gi,gj,gk)
#endif	    
	      do k=1,nz
		    do i=1,nx
			  gi=nx*coords(1)+i
			  gk=nz*coords(3)+k
			  if(gi>1 .and. gi<lx .and. gk>1 .and. gk<lz)then
	          j=1
	          gj=ny*coords(2)+j
	          if(gj.eq.1)then
#ifdef DENSRATIO
			    rhophi_loc = rhophi(i,j+1,k)
#else
		        rhophi_loc = ONE 
#endif

#ifdef TWOCOMPONENT

			    visc_loc=(rho_r*visc1*selphi(i,j+1,k,flip)+(ONE-selphi(i,j+1,k,flip))*visc2*rho_b)/rhophi_loc  

						
			    tau_loc=(visc_loc/cssq + HALF) !è una tau
						
			    omega_loc=ONE/tau_loc !è una omega
						
#else
			    omega_loc=omega
#endif	
				utmp=0.0_db
				vtmp=0.0_db!v(i,j+1,k)	
				wtmp=0.0_db
				!$acc loop seq
				do l=1,nlinks
				  ! lopp=opp(l)
				   ii=i+ex(l)
				   jj=j+ey(l)
				   kk=k+ez(l)
				  
				  if(isfluid(ii,jj,kk).ne.-1) cycle 
				  rhotemp=rho(ii,jj,kk)
				  !w(i,j,k)=w(i,j,k+ez(l))
				  uu=HALF*(utmp*utmp + vtmp*vtmp+ wtmp*wtmp)/cssq
				  udotc=(utmp*dex(l) + vtmp*dey(l)+ wtmp*dez(l))/cssq
				  feq=p(l)*(rhotemp + udotc+ HALF*udotc*udotc - uu)
				  fneq1=(HALF/(cssq*cssq))*( (dex(l)*dex(l)-cssq)*pxx(ii,jj,kk) &
				   + (dey(l)*dey(l)-cssq)*pyy(ii,jj,kk) + (dez(l)*dez(l)-cssq)*pzz(ii,jj,kk) &
				   + TWO*(dex(l)*dey(l))*pxy(ii,jj,kk) + TWO*(dex(l)*dez(l))*pxz(ii,jj,kk) &
				   + TWO*(dey(l)*dez(l))*pyz(ii,jj,kk))
				  ! F_discr = p(l)*(((dex(l) - utmp) + udotc * dex(l))*fux(i,j,k) &
				   ! + ((dey(l) - vtmp) + udotc * dey(l))*fvy(i,j,k) &
				   ! + ((dez(l) - wtmp) + udotc * dez(l))*fwz(i,j,k))/cssq
				  f(ii,jj,kk,l)=feq + (ONE-omega_loc)*p(l)*fneq1 !+ HALF*(F_discr)
				  
				enddo
			  endif
			  j=ny
	          gj=ny*coords(2)+j
	          if(gj.eq.ly)then
#ifdef DENSRATIO
			    rhophi_loc = rhophi(i,j-1,k)
#else
		        rhophi_loc = ONE 
#endif

#ifdef TWOCOMPONENT

			    visc_loc=(rho_r*visc1*selphi(i,j-1,k,flip)+(ONE-selphi(i,j-1,k,flip))*visc2*rho_b)/rhophi_loc  

						
			    tau_loc=(visc_loc/cssq + HALF) !è una tau
						
			    omega_loc=ONE/tau_loc !è una omega
						
#else
			    omega_loc=omega
#endif
				utmp=0.0_db
				vtmp=0.0_db!v(i,j-1,k)	
				wtmp=0.0_db 
				!$acc loop seq
				do l=1,nlinks
				  ! lopp=opp(l)
				   ii=i+ex(l)
				   jj=j+ey(l)
				   kk=k+ez(l)
				  
				  if(isfluid(ii,jj,kk).ne.-1) cycle 
				  rhotemp=rho(ii,jj,kk)
				  !w(i,j,k)=w(i,j,k+ez(l))
				  uu=HALF*(utmp*utmp + vtmp*vtmp+ wtmp*wtmp)/cssq
				  udotc=(utmp*dex(l) + vtmp*dey(l)+ wtmp*dez(l))/cssq
				  feq=p(l)*(rhotemp + udotc+ HALF*udotc*udotc - uu)
				  fneq1=(HALF/(cssq*cssq))*( (dex(l)*dex(l)-cssq)*pxx(ii,jj,kk) &
				   + (dey(l)*dey(l)-cssq)*pyy(ii,jj,kk) + (dez(l)*dez(l)-cssq)*pzz(ii,jj,kk) &
				   + TWO*(dex(l)*dey(l))*pxy(ii,jj,kk) + TWO*(dex(l)*dez(l))*pxz(ii,jj,kk) &
				   + TWO*(dey(l)*dez(l))*pyz(ii,jj,kk))
				  ! F_discr = p(l)*(((dex(l) - utmp) + udotc * dex(l))*fux(i,j,k) &
				   ! + ((dey(l) - vtmp) + udotc * dey(l))*fvy(i,j,k) &
				   ! + ((dez(l) - wtmp) + udotc * dez(l))*fwz(i,j,k))/cssq
				  f(ii,jj,kk,l)=feq + (ONE-omega_loc)*p(l)*fneq1 !+ HALF*(F_discr)
				  
				enddo
			  endif
			  endif
	        enddo
	      enddo
#ifdef ACCNOKERNELS
      !$acc end parallel loop
#else
      !$acc end kernels
#endif
	    endif	    

	    if(pbc_z.eq.0)then
#ifdef ACCNOKERNELS
		!$acc parallel loop collapse(2) independent present(selphi &
		!$acc& ) private(i,j,k,l,gi,gj,gk)
#else
		!$acc kernels present(selphi &
		!$acc& )
		!$acc loop collapse(2) independent private(i,j,k,l,gi,gj,gk)
#endif	    
	      do j=1,ny
		    do i=1,nx
			  gi=nx*coords(1)+i
			  gj=ny*coords(2)+j
			  if(gi>1 .and. gi<lx .and. gj>1 .and. gj<ly)then
	          k=1
	          gk=nz*coords(3)+k
	          if(gk.eq.1)then
#ifdef DENSRATIO
			    rhophi_loc = rhophi(i,j,k+1)
#else
		        rhophi_loc = ONE 
#endif

#ifdef TWOCOMPONENT

			    visc_loc=(rho_r*visc1*selphi(i,j,k+1,flip)+(ONE-selphi(i,j,k+1,flip))*visc2*rho_b)/rhophi_loc  

						
			    tau_loc=(visc_loc/cssq + HALF) !è una tau
						
			    omega_loc=ONE/tau_loc !è una omega
						
#else
			    omega_loc=omega
#endif
				utmp=0.0_db
				vtmp=ZERO	
				wtmp=ZERO!w(i,j,k+1)
				!$acc loop seq
				do l=1,nlinks
				  ! lopp=opp(l)
				   ii=i+ex(l)
				   jj=j+ey(l)
				   kk=k+ez(l)
				  
				  if(isfluid(ii,jj,kk).ne.-1) cycle 
#ifdef IMPOSED_PRESSURE_GRADIENT
				  rhotemp=rhoIN
#else
				  rhotemp=rho(ii,jj,kk)
#endif
				  !w(i,j,k)=w(i,j,k+ez(l))
				  uu=HALF*(utmp*utmp + vtmp*vtmp+ wtmp*wtmp)/cssq
				  udotc=(utmp*dex(l) + vtmp*dey(l)+ wtmp*dez(l))/cssq
				  feq=p(l)*(rhotemp + udotc+ HALF*udotc*udotc - uu)
				  fneq1=(HALF/(cssq*cssq))*( (dex(l)*dex(l)-cssq)*pxx(ii,jj,kk) &
				   + (dey(l)*dey(l)-cssq)*pyy(ii,jj,kk) + (dez(l)*dez(l)-cssq)*pzz(ii,jj,kk) &
				   + TWO*(dex(l)*dey(l))*pxy(ii,jj,kk) + TWO*(dex(l)*dez(l))*pxz(ii,jj,kk) &
				   + TWO*(dey(l)*dez(l))*pyz(ii,jj,kk))
				  ! F_discr = p(l)*(((dex(l) - utmp) + udotc * dex(l))*fux(i,j,k) &
				   ! + ((dey(l) - vtmp) + udotc * dey(l))*fvy(i,j,k) &
				   ! + ((dez(l) - wtmp) + udotc * dez(l))*fwz(i,j,k))/cssq
				  f(ii,jj,kk,l)=feq + (ONE-omega_loc)*p(l)*fneq1 !+ HALF*(F_discr)
				  
				enddo
			  endif
			  k=nz
	          gk=nz*coords(3)+k
	          if(gk.eq.lz)then
#ifdef DENSRATIO
			    rhophi_loc = rhophi(i,j,k-1)
#else
		        rhophi_loc = ONE 
#endif

#ifdef TWOCOMPONENT

			    visc_loc=(rho_r*visc1*selphi(i,j,k-1,flip)+(ONE-selphi(i,j,k-1,flip))*visc2*rho_b)/rhophi_loc  

						
			    tau_loc=(visc_loc/cssq + HALF) !è una tau
						
			    omega_loc=ONE/tau_loc !è una omega
						
#else
			    omega_loc=omega
#endif
				utmp=0.0_db
				vtmp=ZERO	
				wtmp=ZERO!w(i,j,k-1)
				!$acc loop seq
				do l=1,nlinks
				  ! lopp=opp(l)
				   ii=i+ex(l)
				   jj=j+ey(l)
				   kk=k+ez(l)
				  
				  if(isfluid(ii,jj,kk).ne.-1) cycle 
#ifdef IMPOSED_PRESSURE_GRADIENT
				  rhotemp=rhoOUT
#else
				  rhotemp=rho(ii,jj,kk)
#endif
				  !w(i,j,k)=w(i,j,k+ez(l))
				  uu=HALF*(utmp*utmp + vtmp*vtmp+ wtmp*wtmp)/cssq
				  udotc=(utmp*dex(l) + vtmp*dey(l)+ wtmp*dez(l))/cssq
				  feq=p(l)*(rhotemp + udotc+ HALF*udotc*udotc - uu)
				  fneq1=(HALF/(cssq*cssq))*( (dex(l)*dex(l)-cssq)*pxx(ii,jj,kk) &
				   + (dey(l)*dey(l)-cssq)*pyy(ii,jj,kk) + (dez(l)*dez(l)-cssq)*pzz(ii,jj,kk) &
				   + TWO*(dex(l)*dey(l))*pxy(ii,jj,kk) + TWO*(dex(l)*dez(l))*pxz(ii,jj,kk) &
				   + TWO*(dey(l)*dez(l))*pyz(ii,jj,kk))
				  ! F_discr = p(l)*(((dex(l) - utmp) + udotc * dex(l))*fux(i,j,k) &
				   ! + ((dey(l) - vtmp) + udotc * dey(l))*fvy(i,j,k) &
				   ! + ((dez(l) - wtmp) + udotc * dez(l))*fwz(i,j,k))/cssq
				  f(ii,jj,kk,l)=feq + (ONE-omega_loc)*p(l)*fneq1 !+ HALF*(F_discr)
				  
				enddo
			  endif
			  endif
	        enddo
	      enddo
#ifdef ACCNOKERNELS
      !$acc end parallel loop
#else
      !$acc end kernels
#endif
	    endif	

!*******************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  phase field bcs3D
#ifdef TWOCOMPONENT

#ifdef ACCNOKERNELS
		!$acc parallel loop collapse(3) independent present(selphi &
		!$acc& ) private(i,j,k,l,gi,gj,gk &
		!$acc& ,phitemp,conter,phi_fluid,gradfix,gradfiy,grad_parallel,theta_rad,cot_theta,phi_ghost,dphi_dz&
		!$acc& ,gradfiz)
#else
		!$acc kernels present(selphi &
		!$acc& )
		!$acc loop collapse(3) independent private(i,j,k,l,gi,gj,gk &
		!$acc& ,phitemp,conter,phi_fluid,gradfix,gradfiy,grad_parallel,theta_rad,cot_theta,phi_ghost,dphi_dz &
		!$acc& ,gradfiz)
#endif
		do k = 1, nz
		  do j = 1, ny
			do i = 1, nx
			  if (isfluid(i,j,k) .ne. 0) cycle ! Only solid nodes

			  found = .false.
			  do l = 1, 6
				ii = i + ex(l)
				jj = j + ey(l)
				kk = k + ez(l)

				if (isfluid(ii,jj,kk) .ne. -1) cycle  ! only fluid neighbor

				! Found fluid neighbor: enforce contact angle via ghost node extrapolation
				phi_fluid = selphi(ii,jj,kk,flop)

				! Estimate gradient parallel to wall
				gradfix=normx(ii,jj,kk)*modgrad(ii,jj,kk)
				gradfiy=normy(ii,jj,kk)*modgrad(ii,jj,kk)
				gradfiz=normz(ii,jj,kk)*modgrad(ii,jj,kk)
				if(l.eq.1 .or. l.eq.2)then
					grad_parallel = sqrt(gradfiy**2 + gradfiz**2)
				elseif(l.eq.3 .or. l.eq.4)then
					grad_parallel = sqrt(gradfix**2 + gradfiz**2)
				elseif(l.eq.5 .or. l.eq.6)then
					grad_parallel = sqrt(gradfix**2 + gradfiy**2)
				endif
				
				! Contact angle correction
				theta_rad = (180.0_db-wettab_r) * pi_greek / 180.0_db
				cot_theta = 1.0_db / tan(theta_rad)

				dphi_dz = - grad_parallel * cot_theta 

				  

				phi_ghost = phi_fluid + dphi_dz  ! extrapolate from fluid node

				! Clamp to [0,1]
				selphi(i,j,k,flop) = max(0.0_db, min(1.0_db, phi_ghost))

				found = .true.
				exit
			  end do

			  if (.not. found) cycle  ! no fluid neighbor → skip
			end do
		  end do
		end do	
#ifdef ACCNOKERNELS
      !$acc end parallel loop
#else
      !$acc end kernels
#endif
					
	   if(pbc_x.eq.0)then
#ifdef ACCNOKERNELS
		!$acc parallel loop collapse(2) independent present(selphi &
		!$acc& ) private(i,j,k,l,gi,gj,gk)
#else
		!$acc kernels present(selphi &
		!$acc& )
		!$acc loop collapse(2) independent private(i,j,k,l,gi,gj,gk)
#endif	    
	      do k=1,nz
			do j=1,ny
			  gj=ny*coords(2)+j
			  gk=nz*coords(3)+k
			  if(gj>1 .and. gj<ly .and. gk>1 .and. gk<lz)then
	          i=1
	          gi=nx*coords(1)+i
	          if(gi.eq.1)then
				selphi(i,j,k,flop)=selphi(i+1,j,k,flop)
			  endif
			  i=nx
	          gi=nx*coords(1)+i
	          if(gi.eq.lx)then
				selphi(i,j,k,flop)=selphi(i-1,j,k,flop)
			  endif
			  endif
	        enddo
	      enddo
#ifdef ACCNOKERNELS
      !$acc end parallel loop
#else
      !$acc end kernels
#endif
	    endif	
	    
	    if(pbc_y.eq.0)then
#ifdef ACCNOKERNELS
		!$acc parallel loop collapse(2) independent present(selphi)
#else
		!$acc kernels present(selphi &
		!$acc& )
		!$acc loop collapse(2) independent private(i,j,k,l,gi,gj,gk)
#endif	    
	      do k=1,nz
		    do i=1,nx
			  gi=nx*coords(1)+i
			  gk=nz*coords(3)+k
			  if(gi>1 .and. gi<lx .and. gk>1 .and. gk<lz)then
	          j=1
	          gj=ny*coords(2)+j
	          if(gj.eq.1)then
				selphi(i,j,k,flop)=selphi(i,j+1,k,flop)
			  endif
			  j=ny
	          gj=ny*coords(2)+j
	          if(gj.eq.ly)then
				selphi(i,j,k,flop)=selphi(i,j-1,k,flop) 
			  endif
			  endif
	        enddo
	      enddo
#ifdef ACCNOKERNELS
      !$acc end parallel loop
#else
      !$acc end kernels
#endif
	    endif	    

	    if(pbc_z.eq.0)then
#ifdef ACCNOKERNELS
		!$acc parallel loop collapse(2) independent present(selphi &
		!$acc& ) private(i,j,k,l,gi,gj,gk)
#else
		!$acc kernels present(selphi &
		!$acc& )
		!$acc loop collapse(2) independent private(i,j,k,l,gi,gj,gk)
#endif	    
	      do j=1,ny
		    do i=1,nx
			  gi=nx*coords(1)+i
			  gj=ny*coords(2)+j
			  if(gi>1 .and. gi<lx .and. gj>1 .and. gj<ly)then
	          k=1
	          gk=nz*coords(3)+k
	          if(gk.eq.1)then
				selphi(i,j,k,flop)=selphi(i,j,k+1,flop)
			  endif
			  k=nz
	          gk=nz*coords(3)+k
	          if(gk.eq.lz)then
				selphi(i,j,k,flop)=selphi(i,j,k-1,flop)
			  endif
			  endif
	        enddo
	      enddo
#ifdef ACCNOKERNELS
      !$acc end parallel loop
#else
      !$acc end kernels
#endif
	    endif	
	    			

!******************lagrange multiplier
#ifdef BCPHIFLUX

	    if(pbc_x.eq.0)then
#ifdef ACCNOKERNELS
		!$acc parallel loop collapse(2) reduction(+:global_phi_change) present(selphi &
		!$acc& ) private(i,j,k,l,gi,gj,gk &
		!$acc& ,phitemp,conter)
#else
		!$acc kernels present(selphi &
		!$acc& )
		!$acc loop collapse(2) reduction(+:global_phi_change) private(i,j,k,l,gi,gj,gk &
		!$acc& ,phitemp,conter)
#endif	    
	      do k=1,nz
			do j=1,ny
			  gj=ny*coords(2)+j
			  gk=nz*coords(3)+k
	          i=2
	          gi=nx*coords(1)+i
	          if(gi.eq.2)then
				if(abs(isfluid(i,j,k)).eq.1)global_phi_change = global_phi_change + u(i,j,k)*selphi(i,j,k,flop)
			  endif
			  i=nx-1
	          gi=nx*coords(1)+i
	          if(gi.eq.lx-1)then
				if(abs(isfluid(i,j,k)).eq.1)global_phi_change = global_phi_change - u(i,j,k)*selphi(i,j,k,flop)
			  endif
	        enddo
	      enddo
#ifdef ACCNOKERNELS
      !$acc end parallel loop
#else
      !$acc end kernels
#endif
	    endif	
	    
	    if(pbc_y.eq.0)then
#ifdef ACCNOKERNELS
		!$acc parallel loop collapse(2) reduction(+:global_phi_change) present(selphi &
		!$acc& ) private(i,j,k,l,gi,gj,gk &
		!$acc& ,phitemp,conter)
#else
		!$acc kernels present(selphi &
		!$acc& )
		!$acc loop collapse(2) reduction(+:global_phi_change) private(i,j,k,l,gi,gj,gk &
		!$acc& ,phitemp,conter)
#endif	    
	      do k=1,nz
		    do i=1,nx
			  gi=nx*coords(1)+i
			  gk=nz*coords(3)+k
	          j=2
	          gj=ny*coords(2)+j
	          if(gj.eq.2)then
				if(abs(isfluid(i,j,k)).eq.1)global_phi_change = global_phi_change + v(i,j,k)*selphi(i,j,k,flop)
			  endif
			  j=ny-1
	          gj=ny*coords(2)+j
	          if(gj.eq.ly-1)then
				if(abs(isfluid(i,j,k)).eq.1)global_phi_change = global_phi_change - v(i,j,k)*selphi(i,j,k,flop)
			  endif
	        enddo
	      enddo
#ifdef ACCNOKERNELS
      !$acc end parallel loop
#else
      !$acc end kernels
#endif
	    endif	    

	    if(pbc_z.eq.0)then
#ifdef ACCNOKERNELS
		!$acc parallel loop collapse(2) reduction(+:global_phi_change) present(selphi &
		!$acc& ) private(i,j,k,l,gi,gj,gk &
		!$acc& ,phitemp,conter)
#else
		!$acc kernels present(selphi &
		!$acc& )
		!$acc loop collapse(2) reduction(+:global_phi_change) private(i,j,k,l,gi,gj,gk &
		!$acc& ,phitemp,conter)
#endif	    
	      do j=1,ny
		    do i=1,nx
			  gi=nx*coords(1)+i
			  gj=ny*coords(2)+j
	          k=2
	          gk=nz*coords(3)+k
	          if(gk.eq.2)then
				if(abs(isfluid(i,j,k)).eq.1)global_phi_change = global_phi_change + w(i,j,k)*selphi(i,j,k,flop)
			  endif
			  k=nz-1
	          gk=nz*coords(3)+k
	          if(gk.eq.lz-1)then
				if(abs(isfluid(i,j,k)).eq.1)global_phi_change = global_phi_change - w(i,j,k)*selphi(i,j,k,flop)
			  endif
	        enddo
	      enddo
#ifdef ACCNOKERNELS
      !$acc end parallel loop
#else
      !$acc end kernels
#endif
	    endif	
      !$acc wait
	  !$acc update host(global_phi_change)
	  !$acc wait
	  call sum_world_float(global_phi_change)
#endif


#ifdef ACCNOKERNELS
		!$acc parallel loop collapse(3) reduction(+:global_phi_sum) reduction(+:global_count) present(selphi &
		!$acc& ) private(i,j,k,l,gi,gj,gk &
		!$acc& )
#else
		!$acc kernels present(selphi &
		!$acc& )
		!$acc loop collapse(3) reduction(+:global_phi_sum) reduction(+:global_count) private(i,j,k,l,gi,gj,gk &
		!$acc& )
#endif
		do k=1,nz
			do j=1,ny
				do i=1,nx
						
					if(abs(isfluid(i,j,k)).ne.1) cycle

					global_phi_sum = global_phi_sum + selphi(i,j,k,flop)
					if(selphi(i,j,k,flop)>0.5_db .and. selphi(i,j,k,flop)<0.9_db)then
						global_count = global_count + 1
					endif

				enddo
			enddo
		enddo
#ifdef ACCNOKERNELS
      !$acc end parallel loop
#else
      !$acc end kernels
#endif
		!$acc wait
		!$acc update host(global_phi_sum,global_count)
		!$acc wait
		call sum_world_float(global_phi_sum)
		call sum_world_int(global_count)
		
		dphi = (global_phi_sum_ini+global_phi_change) - global_phi_sum
		corr = dphi / real(global_count)
		
		global_phi_sum_new=global_phi_sum
		global_phi_change_new=global_phi_change
		global_count_new=global_count
		
		global_phi_change=ZERO
		global_phi_sum=ZERO
	    global_count=0
		!$acc wait
		!$acc update device(corr,global_phi_sum,global_count)
#ifdef ACCNOKERNELS
		!$acc parallel loop collapse(3) present(selphi &
		!$acc& ) private(i,j,k,l,gi,gj,gk &
		!$acc& )
#else
		!$acc kernels present(selphi &
		!$acc& )
		!$acc loop collapse(3) private(i,j,k,l,gi,gj,gk &
		!$acc& )
#endif
		do k=1,nz
			do j=1,ny
				do i=1,nx
						
					if(abs(isfluid(i,j,k)).ne.1) cycle
					if(selphi(i,j,k,flop)>0.5_db .and. selphi(i,j,k,flop)<0.9_db)then
						selphi(i,j,k,flop)=selphi(i,j,k,flop) + corr
					endif
					! if(selphi(i,j,k,flop)>1.0_db)then
						! selphi(i,j,k,flop)=1.0_db
					! endif
				enddo
			enddo
		enddo
#ifdef ACCNOKERNELS
      !$acc end parallel loop
#else
      !$acc end kernels
#endif
#endif
#endif


   endsubroutine


endmodule
