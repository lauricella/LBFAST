#include "defines.h"

module lb_kernels

   use vars
   use mpi_template, only: coords,myoffset,sum_world_float

   implicit none

contains

   !****************************************************************************!

   subroutine compute_densityratio

      implicit none
	  real(kind=db) :: filoc
#ifdef DENSRATIO
!
#ifdef ACCNOKERNELS
      !$acc parallel loop collapse(3) present(selphi,rhophi,isfluid)
#else
      !$acc kernels present(selphi, rhophi,isfluid)
      !$acc loop collapse(3) private(i,j,k)
#endif
      do k=1,nz
         do j=1,ny
            do i=1,nx
               if(abs(isfluid(i,j,k)).eq.1)then
  
				  rhophi(i,j,k)=rho_r*selphi(i,j,k,flop)+(1.0_db-selphi(i,j,k,flop))*rho_b
				 
               endif
            enddo
         enddo
      enddo
#ifdef ACCNOKERNELS
      !$acc end parallel loop
#else
      !$acc end kernels
#endif
!
#endif
   endsubroutine compute_densityratio
  !****************************************************************************! 
   subroutine moments_LB

      implicit none

      real(kind=db) :: fneq1,feq,forcex,forcey,forcez,F_discr,rhophi_loc
	  
#ifdef TWOCOMPONENT
	  real(kind=db) ::gradfix,gradfiy,gradfiz
#ifdef CSF
      real(kind=db) :: curvature
#endif
#endif
      
#ifdef DENSRATIO
      real(kind=db) :: gradrhox,gradrhoy,gradrhoz,omega_loc, visc_loc,tau_loc
#endif

#ifdef ACCNOKERNELS
          !$acc parallel loop collapse(3) present(pxx,pyy,pzz,pxy,pxz,pyz,rho,u,v,w,fux,fvy,fwz &
#ifdef TWOCOMPONENT
          !$acc& ,selphi,modgrad, lap_phi,normx,normy,normz &
#endif            
#ifdef DENSRATIO
          !$acc& ,rhophi &
#endif
#ifdef MULTIHIT
          !$acc& ,ABCx,ABCy,ABCz &
#endif
#ifdef ELASTIC_FORCE
          !$acc& ,u_ref,v_ref,w_ref &
#endif
          !$acc& ) private(i,j,k,l,forcex,forcey,forcez,F_discr,fneq1,feq,rhophi_loc &
#ifdef DENSRATIO
      !$acc& ,gradrhox,gradrhoy,gradrhoz,omega_loc,visc_loc,tau_loc &
#endif
#ifdef CSF
	  !$acc& , curvature &
#endif
#ifdef TWOCOMPONENT
          !$acc& ,gradfix,gradfiy,gradfiz &
#endif
      !$acc& )
#else
	  !$acc kernels present(pxx,pyy,pzz,pxy,pxz,pyz,rho,u,v,w,fux,fvy,fwz &
#ifdef TWOCOMPONENT
          !$acc& ,selphi,modgrad, lap_phi,normx,normy,normz &
#endif            
#ifdef DENSRATIO
          !$acc& ,rhophi &
#endif
#ifdef MULTIHIT
          !$acc& ,ABCx,ABCy,ABCz &
#endif
#ifdef ELASTIC_FORCE
          !$acc& ,u_ref,v_ref,w_ref &
#endif
          !$acc& )
	  !$acc loop collapse(3) private(i,j,k,l,forcex,forcey,forcez,F_discr,fneq1,feq,rhophi_loc &
#ifdef DENSRATIO
      !$acc& ,gradrhox,gradrhoy,gradrhoz,omega_loc,visc_loc,tau_loc &
#endif
#ifdef CSF
	  !$acc& , curvature &
#endif
#ifdef TWOCOMPONENT
          !$acc& ,gradfix,gradfiy,gradfiz &
#endif
      !$acc& )
#endif
      do k=1,nz
         do j=1,ny
            do i=1,nx
               if(abs(isfluid(i,j,k)).eq.1)then
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
				  !$acc loop seq
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
				  !$acc loop seq
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
			   endif
            enddo
         enddo
      enddo
#ifdef ACCNOKERNELS
      !$acc end parallel loop
#else
      !$acc end kernels
#endif
   endsubroutine moments_LB  
   


#ifdef TWOCOMPONENT  
   subroutine compute_laplacian_phi

      implicit none
#ifdef ACCNOKERNELS
      !$acc parallel loop collapse(3) present(isfluid,selphi,lap_phi)
#else
      !$acc kernels present(isfluid,selphi,lap_phi)
      !$acc loop collapse(3) private(i,j,k)
#endif
      do k=1,nz
         do j=1,ny
            do i=1,nx
               if(abs(isfluid(i,j,k)).eq.1)then
#ifdef LAPLACIAN_SEVEN_POINTS
                  lap_phi(i,j,k)=selphi(i+1,j,k,flop)+selphi(i-1,j,k,flop)&
				  +selphi(i,j+1,k,flop)+selphi(i,j-1,k,flop)&
				  +selphi(i,j,k+1,flop)+selphi(i,j,k-1,flop)-6.0_db*selphi(i,j,k,flop)
#else
                  lap_phi(i,j,k)=(2.0_db/cssq)*(selphi(i,j,k,flop)*(p0-1.0_db) + &
                   ( p1*(selphi(i+1,j,k,flop)+selphi(i-1,j,k,flop) + &
                   selphi(i,j+1,k,flop)+selphi(i,j-1,k,flop) + &
                   selphi(i,j,k+1,flop)+selphi(i,j,k-1,flop)) + &
                   p2*( (selphi(i+1,j+1,k,flop)+selphi(i-1,j-1,k,flop))+ &
                   (selphi(i+1,j-1,k,flop)+selphi(i-1,j+1,k,flop))+ &
                   (selphi(i+1,j,k+1,flop)+selphi(i-1,j,k-1,flop))+ &
                   (selphi(i+1,j,k-1,flop)+selphi(i-1,j,k+1,flop)) + &
                   (selphi(i,j+1,k+1,flop)+selphi(i,j-1,k-1,flop))+ &
                   (selphi(i,j+1,k-1,flop)+selphi(i,j-1,k+1,flop)) )  + &
                   p3*((selphi(i+1,j+1,k+1,flop)+selphi(i-1,j-1,k-1,flop))+ &
                   (selphi(i+1,j-1,k-1,flop)+selphi(i-1,j+1,k+1,flop))+ &
                   (selphi(i+1,j-1,k+1,flop)+selphi(i-1,j+1,k-1,flop))+ &
                   (selphi(i+1,j+1,k-1,flop)+selphi(i-1,j-1,k+1,flop)))))
#endif
   
               endif
            enddo
         enddo
      enddo	
#ifdef ACCNOKERNELS
      !$acc end parallel loop
#else
      !$acc end kernels
#endif
   endsubroutine compute_laplacian_phi
   !****************************************************************************!
   subroutine compute_norm_interface

      implicit none
      real(kind=db):: grad_fix,grad_fiy,grad_fiz,mod_grad
#ifdef ACCNOKERNELS
      !$acc parallel loop collapse(3) present(isfluid,selphi, modgrad, normx,normy,normz,arr_x,arr_y,arr_z) private(l,grad_fix,grad_fiy,grad_fiz,mod_grad)
#else
      !$acc kernels present(isfluid,selphi, modgrad, normx,normy,normz,arr_x,arr_y,arr_z)
      !$acc loop collapse(3) private(l,grad_fix,grad_fiy,grad_fiz,mod_grad)
#endif
      do k=1,nz
         do j=1,ny
            do i=1,nx
               if(abs(isfluid(i,j,k)).eq.1)then


                  grad_fix=3.0_db*(p1*(selphi(i+1,j,k,flop)-selphi(i-1,j,k,flop)) + &
                     p2*( (selphi(i+1,j+1,k,flop)-selphi(i-1,j-1,k,flop))+ &
                     (selphi(i+1,j-1,k,flop)-selphi(i-1,j+1,k,flop))+ &
                     (selphi(i+1,j,k+1,flop)-selphi(i-1,j,k-1,flop))+ &
                     (selphi(i+1,j,k-1,flop)-selphi(i-1,j,k+1,flop)) )  + &
                     p3*((selphi(i+1,j+1,k+1,flop)-selphi(i-1,j-1,k-1,flop))+ &
                     (selphi(i+1,j-1,k-1,flop)-selphi(i-1,j+1,k+1,flop))+ &
                     (selphi(i+1,j-1,k+1,flop)-selphi(i-1,j+1,k-1,flop))+ &
                     (selphi(i+1,j+1,k-1,flop)-selphi(i-1,j-1,k+1,flop))))

                  grad_fiy=3.0_db*(p1*(selphi(i,j+1,k,flop)-selphi(i,j-1,k,flop)) + &
                     p2*((selphi(i+1,j+1,k,flop)-selphi(i-1,j-1,k,flop))+ &
                     (selphi(i-1,j+1,k,flop)-selphi(i+1,j-1,k,flop))+ &
                     (selphi(i,j+1,k+1,flop)-selphi(i,j-1,k-1,flop))+ &
                     (selphi(i,j+1,k-1,flop)-selphi(i,j-1,k+1,flop)) ) + &
                     p3*((selphi(i+1,j+1,k+1,flop)-selphi(i-1,j-1,k-1,flop))+ &
                     (selphi(i-1,j+1,k-1,flop)-selphi(i+1,j-1,k+1,flop))+ &
                     (selphi(i+1,j+1,k-1,flop)-selphi(i-1,j-1,k+1,flop))+ &
                     (selphi(i-1,j+1,k+1,flop)-selphi(i+1,j-1,k-1,flop))))

                  grad_fiz=3.0_db*(p1*(selphi(i,j,k+1,flop)-selphi(i,j,k-1,flop)) + &
                     p2*((selphi(i+1,j,k+1,flop)-selphi(i-1,j,k-1,flop))+ &
                     (selphi(i-1,j,k+1,flop)-selphi(i+1,j,k-1,flop))+ &
                     (selphi(i,j+1,k+1,flop)-selphi(i,j-1,k-1,flop))+ &
                     (selphi(i,j-1,k+1,flop)-selphi(i,j+1,k-1,flop)) ) + &
                     p3*((selphi(i+1,j+1,k+1,flop)-selphi(i-1,j-1,k-1,flop)) &
                     +(selphi(i-1,j-1,k+1,flop)-selphi(i+1,j+1,k-1,flop))+ &
                     (selphi(i+1,j-1,k+1,flop)-selphi(i-1,j+1,k-1,flop))+ &
                     (selphi(i-1,j+1,k+1,flop)-selphi(i+1,j-1,k-1,flop))))

  
                  mod_grad= sqrt(grad_fix**TWO + grad_fiy**TWO + grad_fiz**TWO)

                  normx(i,j,k)=grad_fix/(mod_grad+1.0e-9)
                  normy(i,j,k)=grad_fiy/(mod_grad+1.0e-9)
                  normz(i,j,k)=grad_fiz/(mod_grad+1.0e-9)

                  arr_x(i,j,k)= selphi(i,j,k,flop)*(1.0_db-selphi(i,j,k,flop))*normx(i,j,k)
				  arr_y(i,j,k)= selphi(i,j,k,flop)*(1.0_db-selphi(i,j,k,flop))*normy(i,j,k)
				  arr_z(i,j,k)= selphi(i,j,k,flop)*(1.0_db-selphi(i,j,k,flop))*normz(i,j,k)
			  
				  !lap_phi here
				
                  modgrad(i,j,k)=mod_grad

               endif
            enddo
         enddo
      enddo
#ifdef ACCNOKERNELS
      !$acc end parallel loop
#else
      !$acc end kernels
#endif
   endsubroutine compute_norm_interface
#endif
   !****************************************************************************!
   subroutine fused_LB

      implicit none
      integer :: gi,gj,gk

      real(kind=db) :: F_discr,fneq1,feq,forcex,forcey,forcez
#ifdef TWOCOMPONENT
      real(kind=db) :: tau_loc
      real(kind=db) :: visc_loc,gradfix,gradfiy,gradfiz
#ifdef MONOD
	  real(kind=db) :: S_mono
#endif
#endif
      real(kind=db) :: omega_loc, rhophi_loc



#ifdef ACCNOKERNELS
      !$acc parallel loop independent collapse(3) present(f,rho,u,v,w,pxx,pyy,pzz,pxy,pxz,pyz,fux,fvy,fwz,isfluid &
#ifdef DENSRATIO
          !$acc& ,rhophi &
#endif
#ifdef TWOCOMPONENT
          !$acc& ,selphi,lap_phi,arr_x,arr_y,arr_z,normx,normy,normz &
#endif
#ifdef MULTIHIT
          !$acc& ,ABCx,ABCy,ABCz &
#endif
#ifdef ELASTIC_FORCE
          !$acc& ,u_ref,v_ref,w_ref &
#endif
          !$acc& ) private(l,i,j,k,F_discr,forcex,forcey,forcez,feq,fneq1,uu,omega_loc,udotc,rhophi_loc &
#ifdef TWOCOMPONENT
          !$acc& ,visc_loc,tau_loc,gradfix,gradfiy,gradfiz &
#ifdef MONOD
		  !$acc& , S_mono &
#endif
#endif
          !$acc& )
#else
      !$acc kernels present(f,rho,u,v,w,pxx,pyy,pzz,pxy,pxz,pyz,fux,fvy,fwz,isfluid &
#ifdef DENSRATIO
          !$acc& ,rhophi &
#endif
#ifdef TWOCOMPONENT
          !$acc& ,selphi,lap_phi,arr_x,arr_y,arr_z,normx,normy,normz,modgrad &
#endif
#ifdef MULTIHIT
          !$acc& ,ABCx,ABCy,ABCz &
#endif
#ifdef ELASTIC_FORCE
          !$acc& ,u_ref,v_ref,w_ref &
#endif
          !$acc& )
      !$acc loop independent collapse(3) private(l,i,j,k,F_discr,forcex,forcey,forcez,feq,fneq1,uu,omega_loc,udotc,rhophi_loc &
#ifdef TWOCOMPONENT
          !$acc& ,visc_loc,tau_loc,gradfix,gradfiy,gradfiz &
#ifdef MONOD
		  !$acc& , S_mono &
#endif
#endif
          !$acc& )
#endif
      do k=1,nz
         do j=1,ny
            do i=1,nx
               if(abs(isfluid(i,j,k)).eq.1)then
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
		 
               endif
            enddo
         enddo
      enddo
#ifdef ACCNOKERNELS
      !$acc end parallel loop
#else
      !$acc end kernels
#endif
   endsubroutine fused_LB
   !****************************************************************************!
#ifdef REPULSIVE_FLUX	
	subroutine thinfilm_scan_mark
	  implicit none

	  integer :: i,j,k, ii,jj,kk, di,dj,dk
	  integer :: dii,djj,dkk
	  real(kind=db) :: nix,niy,niz, dotn, qloc, qneig, face
	  real(kind=db) :: best_r2, r2, best_face
	  integer :: ii_best, jj_best, kk_best
	  logical :: found
	  real(kind=db), parameter :: eps = 1.0e-12_db
#ifdef ACCNOKERNELS
    !$acc parallel loop independent collapse(3) present(normx,normy,normz,rep_mask,selphi,pair_i,pair_j,pair_k) &
	!$acc& private(i,j,k,di,dj,dk,ii,jj,kk,nix,niy,niz,found,dotn,qloc,qneig, &
	!$acc& face,best_r2,r2,best_face,ii_best,jj_best,kk_best,dii,djj,dkk)
#else
	!$acc kernels present(normx,normy,normz,rep_mask,selphi,pair_i,pair_j,pair_k)
	!$acc loop collapse(3) private(i,j,k,di,dj,dk,ii,jj,kk,nix,niy,niz,found,dotn,qloc,qneig, &
	!$acc& face,best_r2,r2,best_face,ii_best,jj_best,kk_best,dii,djj,dkk)
#endif
	  do k=1,nz
		do j=1,ny
		  do i=1,nx
			rep_mask(i,j,k) = 0
			pair_i(i,j,k) = 0; pair_j(i,j,k) = 0; pair_k(i,j,k) = 0

			! gate: interfacial cell (use clamped phi for q)
			qloc = selphi(i,j,k,flop); qloc = min(max(qloc,0.0_db),1.0_db)
			qloc = qloc*(1.0_db - qloc)
			if (qloc < q_th) cycle

			nix = normx(i,j,k); niy = normy(i,j,k); niz = normz(i,j,k)

			best_r2   = HUGE(1.0_db)
			best_face = -1.0_db
			found     = .false.

			!$acc loop seq
			do di = -win, win
			  !$acc loop seq
			  do dj = -win, win
				!$acc loop seq
				do dk = -win, win
				  if (di==0 .and. dj==0 .and. dk==0) cycle

				  ! ---- 
				  ii = i + di
				  jj = j + dj
				  kk = k + dk

				  ! ---- minimum-image index differences
				  dii = ii - i
				 
				  djj = jj - j
				  
				  dkk = kk - k
				  

				  r2 = real(dii,db)*real(dii,db) + real(djj,db)*real(djj,db) + real(dkk,db)*real(dkk,db)
				  if (r2 < eps) cycle

				  ! ---- neighbor interfacial gate (clamped) + similarity
				  qneig = selphi(ii,jj,kk,flop); qneig = min(max(qneig,0.0_db),1.0_db)
				  qneig = qneig*(1.0_db - qneig)
				  if ( (qneig < q_th) .or. (abs(qneig - qloc) > 0.1_db*max(qloc,1.0e-12_db)) ) cycle

				  ! ---- facing condition (opposite normals): dotn <= cosOppT
				  dotn = nix*normx(ii,jj,kk) + niy*normy(ii,jj,kk) + niz*normz(ii,jj,kk)
				  if (dotn > cosOppT) cycle
				  face = 0.5_db*(1.0_db - dotn)   ! in [0,1]

				  ! ---- pick nearest; tie-break by larger 'face'
				  if (r2 < best_r2 - 1.0e-14_db) then
					best_r2 = r2; best_face = face
					ii_best = ii; jj_best = jj; kk_best = kk
					found   = .true.
				  else if (abs(r2 - best_r2) <= 1.0e-14_db) then
					if (face > best_face) then
					  best_face = face
					  ii_best = ii; jj_best = jj; kk_best = kk
					  found   = .true.
					end if
				  end if

				end do
			  end do
			end do

			if (found) then
			  pair_i(i,j,k) = ii_best
			  pair_j(i,j,k) = jj_best
			  pair_k(i,j,k) = kk_best
			  rep_mask(i,j,k) = 1
			end if

		  end do
		end do
	  end do
#ifdef ACCNOKERNELS
      !$acc end parallel loop
#else
      !$acc end kernels
#endif
	end subroutine thinfilm_scan_mark
			
	subroutine repulsive_flux_normal
	  implicit none
	  integer :: ii,jj,kk
	  real(kind=db) :: q1,q2,qpair,qcl
	  real(kind=db) :: nx1,ny1,nz1, nx2,ny2,nz2
	  real(kind=db) :: dx,dy,dz, r, rinv, face, arg_arcosh, ach, Wfilm, wdth
	  real(kind=db) :: nsx,nsy,nsz, nsmag,alpha,cap,scales
	  real(kind=db), parameter :: eps = 1.0e-9_db
 #ifdef ACCNOKERNELS
		 !$acc parallel loop collapse(3) present(normx,normy,normz,rep_mask,pair_i,pair_j,pair_k,selphi,Jx,Jy,Jz)&
		 !$acc& private(rinv,wfilm,wdth,ach,qcl,q2,q1,nsz,nsy,nsx,nsmag,kk,nz1,ny1,nx1,ny2,nx2,dz,dy,dx,arg_arcosh,r,qpair,jj,ii,face,nz2)
 #else
		 !$acc kernels present(normx,normy,normz,rep_mask,pair_i,pair_j,pair_k,selphi,Jx,Jy,Jz)
		 !$acc loop collapse(3) private(rinv,wfilm,wdth,ach,qcl,q2,q1,nsz,nsy,nsx,nsmag,kk,nz1,ny1,nx1,ny2,nx2,dz,dy,dx,arg_arcosh,r,qpair,jj,ii,face,nz2) 
 #endif
	  do k=1,nz
		do j=1,ny
		  do i=1,nx
			Jx(i,j,k)=0.0_db
			Jy(i,j,k)=0.0_db
			Jz(i,j,k)=0.0_db
			
			if (rep_mask(i,j,k) .ne. 1) cycle

			q1 = selphi(i,j,k,flop)*(1.0_db - selphi(i,j,k,flop))
			
			if (q1 <= eps) cycle

			ii = pair_i(i,j,k)
			jj = pair_j(i,j,k)
			kk = pair_k(i,j,k)

			! line-of-centers
			dx = real(ii - i,db)
			dy = real(jj - j,db)
			dz = real(kk - k,db)
			r  = sqrt(dx*dx + dy*dy + dz*dz)
			if (r <= eps) cycle
			rinv = 1.0_db / r
			dx = dx*rinv; dy = dy*rinv; dz = dz*rinv      ! u

			! normals
			nx1 = normx(i ,j ,k ); ny1 = normy(i ,j ,k ); nz1 = normz(i ,j ,k )
			nx2 = normx(ii,jj,kk); ny2 = normy(ii,jj,kk); nz2 = normz(ii,jj,kk)

			! facing factor in [0,1]
			face = max( 0.0_db, -(nx1*nx2 + ny1*ny2 + nz1*nz2) )

			if (face <= eps) cycle

			! symmetric normal: bisector n1 - n2 (for facing sheets)
			nsx = nx1 - nx2
			nsy = ny1 - ny2
			nsz = nz1 - nz2
			nsmag = sqrt(nsx*nsx + nsy*nsy + nsz*nsz)
			if (nsmag <= eps) then
			  ! fallback to line-of-centers
			  nsx = dx; nsy = dy; nsz = dz
			else
			  nsx = nsx / nsmag
			  nsy = nsy / nsmag
			  nsz = nsz / nsmag
			end if

			! orient so u·nsym >= 0  (partner will flip)
			if (dx*nsx + dy*nsy + dz*nsz < 0.0_db) then
			  nsx = -nsx; nsy = -nsy; nsz = -nsz
			end if

			! symmetric magnitude from qpair
			q2 = selphi(ii,jj,kk,flop)*(1.0_db - selphi(ii,jj,kk,flop))
			qpair = 0.5_db*(q1 + q2)
			qcl   = min( max(qpair, eps), 0.25_db - eps )

			arg_arcosh = 1.0_db / ( 2.0_db*sqrt(qcl) )
			if (arg_arcosh <= 1.0_db) cycle

			ach   = log( arg_arcosh + sqrt(arg_arcosh*arg_arcosh - 1.0_db) )
			Wfilm = width * ach
			if (Wfilm <= 0.0_db) cycle

			wdth  = 1.0_db /(1.0 + wfilm**4.0) ! ( 1.0_db + (1.0_db/Wfilm)**4 )

			! final purely-normal, symmetric repulsive flux
			Jx(i,j,k) = A_rep * wdth * qcl * face * nsx
			Jy(i,j,k) = A_rep * wdth * qcl * face * nsy
			Jz(i,j,k) = A_rep * wdth * qcl * face * nsz
			
			! alpha = 1.5_db
			! cap   = alpha * (abs(Jx(i,j,k))+abs(Jy(i,j,k))+abs(Jz(i,j,k))) 
			! scales = min(1.0_db, selphi(i,j,k,flop) / max(cap, 1.0e-9_db))
			! Jx(i,j,k) = Jx(i,j,k) * scales; Jy(i,j,k) = Jy(i,j,k) * scales; Jz(i,j,k) = Jz(i,j,k) * scales
		  end do
		end do
	  end do
#ifdef ACCNOKERNELS
	  !$acc end parallel loop
#else
	  !$acc end kernels
#endif
	end subroutine repulsive_flux_normal
	
	subroutine repulsive_flux_tangential
	  implicit none
	  
	  integer ::  ii,jj,kk
	  real(kind=db) :: qloc, qcl, Wfilm, wdth,arg_arcosh,ach,qj,qpair
	  real(kind=db) :: nxv, nyv, nzv, dx, dy, dz, dotnd,rinv,nxv2, nyv2, nzv2,bn,bx,by,bz,sx,sy,sz
	  real(kind=db) :: tx, ty, tz, tnorm,alpha,cap,scales

#ifdef ACCNOKERNELS
		!$acc parallel loop collapse(3) present(normx,normy,normz,rep_mask,pair_i,pair_j,pair_k,selphi,Jx,Jy,Jz) &
		!$acc& private(i,j,k,ii,jj,kk,nxv,nyv,nzv,dotnd,wdth,arg_arcosh,ach,qcl,qloc,Wfilm,tx,ty,tz,tnorm,&
		!$acc& rinv,nxv2, nyv2, nzv2,bn,bx,by,bz,sx,sy,sz)
#else
		!$acc kernels present(normx,normy,normz,rep_mask,pair_i,pair_j,pair_k,selphi,Jx,Jy,Jz)
		!$acc loop collapse(3) private(i,j,k,ii,jj,kk,nxv,nyv,nzv,dotnd,wdth,arg_arcosh,ach,qcl,qloc,Wfilm,&
		!$acc& tx,ty,tz,tnorm,rinv,nxv2, nyv2, nzv2,bn,bx,by,bz,sx,sy,sz)
#endif
	  do k=1,nz
		do j=1,ny
		  do i=1,nx
			Jx(i,j,k) = 0.0_db
			Jy(i,j,k) = 0.0_db
			Jz(i,j,k) = 0.0_db
			if (rep_mask(i,j,k).ne.1) cycle
			qloc = selphi(i,j,k,flop)*(1.0_db-selphi(i,j,k,flop))
			if (qloc <= 1.0e-9_db) cycle
			
			ii = pair_i(i,j,k)
			jj = pair_j(i,j,k)
			kk = pair_k(i,j,k)
			
			dx = real(ii - i,db)
			dy = real(jj - j,db)
			dz = real(kk - k,db)
			
			rinv = 1.0_db / max(sqrt(dx*dx + dy*dy + dz*dz),1.0e-9_db) 
			dx = dx*rinv; dy = dy*rinv; dz = dz*rinv   ! u
			! normals
			nxv = normx(i,j,k);  nyv = normy(i,j,k);  nzv = normz(i,j,k)
			nxv2 = normx(ii,jj,kk); nyv2 = normy(ii,jj,kk); nzv2 = normz(ii,jj,kk)
			
			! symmetric tangent
			! b = n1 × n2
			bx = nyv*nzv2 - nzv*nyv2
			by = nzv*nxv2 - nxv*nzv2
			bz = nxv*nyv2 - nyv*nxv2
			bn = sqrt(bx*bx + by*by + bz*bz)

			if (bn > 1.0e-9_db) then
			  ! t_raw = b × (n1+n2) 
			  sx = nxv + nxv2;  sy = nyv + nyv2;  sz = nzv + nzv2
			  tx = by*sz - bz*sy
			  ty = bz*sx - bx*sz
			  tz = bx*sy - by*sx
			else
			  ! fallback: project u onto tangent(n1)
			  dotnd = nxv*dx + nyv*dy + nzv*dz
			  tx = dx - dotnd*nxv
			  ty = dy - dotnd*nyv
			  tz = dz - dotnd*nzv
			end if

			tnorm = sqrt(tx*tx + ty*ty + tz*tz) + 1.0e-9_db
			tx = tx/tnorm; ty = ty/tnorm; tz = tz/tnorm

			! orient so u·t >= 0 (partner has u'=-u -> t'=-t)
			if (dx*tx + dy*ty + dz*tz < 0.0_db) then
			  tx = -tx; ty = -ty; tz = -tz
			end if

			! symmetric magnitude via qpair
			qj = selphi(ii,jj,kk,flop)*(1.0_db - selphi(ii,jj,kk,flop))
			qpair = 0.5_db*(qloc + qj)                          ! or min(qi,qj)
			qcl   = min( max(qpair, 1.0e-9_db), 0.25_db - 1.0e-9_db )

			arg_arcosh = 1.0_db / ( 2.0_db*sqrt(qcl) )
			if (arg_arcosh <= 1.0_db) cycle
			ach   = log( arg_arcosh + sqrt(arg_arcosh*arg_arcosh - 1.0_db) )
			Wfilm = width * ach
			if (Wfilm <= 0.0_db) cycle
			wdth  = 1/abs(Wfilm) !1.0_db / ( 1.0_db + (1.0_db/Wfilm)**4 )

			! tangential flux
			Jx(i,j,k) = A_rep * wdth * qcl * tx
			Jy(i,j,k) = A_rep * wdth * qcl * ty
			Jz(i,j,k) = A_rep * wdth * qcl * tz
			! alpha = 1.5_db
			! cap   = alpha * (abs(Jx(i,j,k))+abs(Jy(i,j,k))+abs(Jz(i,j,k))) 
			! scales = min(1.0_db, selphi(i,j,k,flop) / max(cap, 1.0e-9_db))
			! Jx(i,j,k) = Jx(i,j,k) * scales; Jy(i,j,k) = Jy(i,j,k) * scales; Jz(i,j,k) = Jz(i,j,k) * scales
			
		  end do
		end do
	  end do
#ifdef ACCNOKERNELS
	  !$acc end parallel loop
#else
	  !$acc end kernels
#endif
	end subroutine repulsive_flux_tangential

#endif



endmodule
! normals calculation
 ! 0  1   2  3   4   5   6   7    8   9   10  11   12  13   14  15   16   17   18  19  20  21  22  23  24  25  26
!ex=(/0, 1, -1, 0,  0,  0,  0,  1,  -1,  1,  -1,  0,   0,  0,   0,  1,  -1,  -1,   1,  1, -1,  1, -1, -1,  1,  1, -1/)
!ey=(/0, 0,  0, 1, -1,  0,  0,  1,  -1, -1,   1,  1,  -1,  1,  -1,  0,   0,   0,   0,  1  -1, -1,  1, -1,  1, -1,  1/)
!ez=(/0, 0,  0, 0,  0,  1, -1,  0,   0,  0,   0,  1,  -1, -1,   1,  1,  -1,   1,  -1,  1, -1,  1, -1,  1, -1, -1,  1/)

!

