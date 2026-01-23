#include "defines.h"
module bcs3D

   use vars
   use mpi_template , only: coords,pbc_x,pbc_y,pbc_z,myoffset,myrank, &
    sum_world_float,sum_world_int
   use lb_cuda_driver, only: PHI_int_boundary_cuda,LB_int_boundary_cuda, &
    phi_sum_count_cuda,apply_lagrangian_phi_cuda
!   !$if _OPENACC
!   use openacc
!   !$endif
   implicit none

contains
   !***************************************************
   subroutine bcs_mesoscopic_hfields( &
#ifdef TWOCOMPONENT	    
   phifields_s, &
#endif
   press_in,u_in,v_in,w_in,pxx_in,pyy_in,pzz_in,pxy_in,pxz_in,pyz_in &
   ,press_out,u_out,v_out,w_out,pxx_out,pyy_out,pzz_out,pxy_out,pxz_out,pyz_out)

     implicit none
     
#ifdef TWOCOMPONENT	       
     real(kind=db), allocatable, dimension(:) :: phifields_s
#endif
     real(kind=db), allocatable, dimension(:,:,:) :: &
       press_in,u_in,v_in,w_in,pxx_in,pyy_in,pzz_in,pxy_in,pxz_in,pyz_in &
      ,press_out,u_out,v_out,w_out,pxx_out,pyy_out,pzz_out,pxy_out,pxz_out,pyz_out
      
     integer :: subchords(3)
     integer :: ii,jj,kk,l,lopp,iii,jjj,kkk
 	 real(kind=db) :: feq, fneq1,presstmp,utmp,vtmp,wtmp,rhophi_loc,fpost
#ifdef TWOCOMPONENT	  
 	 real(kind=db) :: phitemp,phi_loc
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
     integer :: xblock,yblock,zblock,myblock
     real(kind=db) :: pxx,pyy,pzz,pxy,pxz,pyz
     real(kind=db) :: opress,ou,ov,ow,opxx,opyy,opzz,opxy,opxz,opyz



#if defined(INTERNAL_OBSTACLES)

     call LB_int_boundary_cuda( &
#ifdef TWOCOMPONENT	       
     phifields_s, &
#endif     
     press_in,u_in,v_in,w_in,pxx_in,pyy_in,pzz_in,pxy_in,pxz_in,pyz_in &
     ,press_out,u_out,v_out,w_out,pxx_out,pyy_out,pzz_out,pxy_out,pxz_out,pyz_out)
     
     if(openbc==0)return
!*****************************************

	 if(pbc_x.eq.0)then
	   gi=2
       subchords(1)=(gi-1)/nx
	   if(subchords(1)==coords(1))then
#ifdef ACCNOKERNELS
		 !$acc parallel loop collapse(2) independent present(press_in,u_in,v_in,w_in,pxx_in,pyy_in,pzz_in,pxy_in,pxz_in,pyz_in &
#ifdef TWOCOMPONENT	
         !$acc& ,phifields_s &
#endif		 
		 !$acc& ) private(i,j,k,l,gi,gj,gk,ii,jj,kk,xblock,yblock,zblock,myblock)
#else
		 !$acc kernels present(press_in,u_in,v_in,w_in,pxx_in,pyy_in,pzz_in,pxy_in,pxz_in,pyz_in &
#ifdef TWOCOMPONENT	
         !$acc& ,phifields_s &
#endif			 
		 !$acc& )
		 !$acc loop collapse(2) independent private(i,j,k,l,gi,gj,gk,ii,jj,kk &
		 !$acc& ,xblock,yblock,zblock,myblock)
#endif	    
	     do k=1,nz
		   do j=1,ny
		 	 gj=ny*coords(2)+j
			 gk=nz*coords(3)+k
			 if(gj>1 .and. gj<ly .and. gk>1 .and. gk<lz)then
	           i=2
	           gi=nx*coords(1)+i
	           
	           if(isfluid(i,j,k) .ne. -1)cycle
	           
	             xblock=(i+2*TILE_DIMx-1)/TILE_DIMx   
                 yblock=(j+2*TILE_DIMy-1)/TILE_DIMy     
                 zblock=(k+2*TILE_DIMz-1)/TILE_DIMz  
                 myblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
                 ii=i-xblock*TILE_DIMx+2*TILE_DIMx
                 jj=j-yblock*TILE_DIMy+2*TILE_DIMy
                 kk=k-zblock*TILE_DIMz+2*TILE_DIMz 
				 
				 pxx=pxx_in(i,j,k)
				 pyy=pyy_in(i,j,k)
				 pzz=pzz_in(i,j,k)
				 pxy=pxy_in(i,j,k)
				 pxz=pxz_in(i,j,k)
				 pyz=pyz_in(i,j,k)
				 
				 
				 presstmp=press
				 utmp=0.0_db !
				 vtmp=0.0_db !
				 wtmp=0.0_db !
				 
	             pxx=pxx - cssq*press_in(i,j,k) - u_in(i,j,k)*u_in(i,j,k) 
	             pyy=pyy - cssq*press_in(i,j,k) - v_in(i,j,k)*v_in(i,j,k) 
	             pzz=pzz - cssq*press_in(i,j,k) - w_in(i,j,k)*w_in(i,j,k) 
	             pxy=pxy - u_in(i,j,k)*v_in(i,j,k)
	             pxz=pxz - u_in(i,j,k)*w_in(i,j,k)
	             pyz=pyz - v_in(i,j,k)*w_in(i,j,k)

#ifdef TWOCOMPONENT	  
                 phi_loc=phifields_s(idx5(ii,jj,kk,1,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields))
#endif	           
#ifdef DENSRATIO
			     rhophi_loc = rho_r*phi_loc+(ONE-phi_loc)*rho_b 
#else
		         rhophi_loc = ONE 
#endif

#ifdef TWOCOMPONENT

			     visc_loc=(rho_r*visc1*phi_loc+(1.0_db-phi_loc)*visc2*rho_b)/rhophi_loc

						
			     tau_loc=(visc_loc*invcssq + HALF) !è una tau
						
			     omega_loc=ONE/tau_loc !è una omega
						
#else
			     omega_loc=omega
#endif		
				 uu=HALF*(utmp*utmp + vtmp*vtmp+ wtmp*wtmp)*invcssq
				 !$acc loop seq
				 do l=1,nlinks
		           lopp=opp(l)
		           iii=i+ex(lopp)
		           jjj=j+ey(lopp)
		           kkk=k+ez(lopp)
		           if(isfluid(iii,jjj,kkk).ne.0) cycle 
		           feq=p(l)*press_in(i,j,k)
		           fneq1=(HALF/(cssq*cssq))*( (dex(l)*dex(l)-cssq)*pxx &
		            + (dey(l)*dey(l)-cssq)*pyy + (dez(l)*dez(l)-cssq)*pzz &
	                + TWO*(dex(l)*dey(l))*pxy + TWO*(dex(l)*dez(l))*pxz &
		            + TWO*(dey(l)*dez(l))*pyz)
                   !F_discr = p(l)*(dex(l)*forcex &
                   ! + dey(l)*forcey &
                   ! + dez(l)*forcez)/(cssq*rhophi_loc)
		           fpost=feq + (ONE-omega_loc)*p(l)*fneq1 !+ HALF*(F_discr)	
		           press_out(i,j,k)=press_out(i,j,k) - fpost
		           u_out(i,j,k)=u_out(i,j,k) - fpost*dex(l)
		           v_out(i,j,k)=v_out(i,j,k) - fpost*dey(l)
		           w_out(i,j,k)=w_out(i,j,k) - fpost*dez(l)
		           pxx_out(i,j,k)=pxx_out(i,j,k) - fpost*dex(l)*dex(l)
                   pyy_out(i,j,k)=pyy_out(i,j,k) - fpost*dey(l)*dey(l)
                   pzz_out(i,j,k)=pzz_out(i,j,k) - fpost*dez(l)*dez(l)
                   pxy_out(i,j,k)=pxy_out(i,j,k) - fpost*dex(l)*dey(l)
                   pxz_out(i,j,k)=pxz_out(i,j,k) - fpost*dex(l)*dez(l)
                   pyz_out(i,j,k)=pyz_out(i,j,k) - fpost*dey(l)*dez(l)	
		           udotc=(utmp*dex(l) + vtmp*dey(l)+ wtmp*dez(l))*invcssq
		           feq=p(l)*(presstmp + udotc+ HALF*udotc*udotc - uu)
		           fneq1=(HALF/(cssq*cssq))*( (dex(l)*dex(l)-cssq)*pxx &
		            + (dey(l)*dey(l)-cssq)*pyy + (dez(l)*dez(l)-cssq)*pzz &
	                + TWO*(dex(l)*dey(l))*pxy + TWO*(dex(l)*dez(l))*pxz &
		            + TWO*(dey(l)*dez(l))*pyz)
                   ! F_discr = p(l)*(((dex(l) - utmp) + udotc * dex(l))*forcex &
                   ! + ((dey(l) - vtmp) + udotc * dey(l))*forcey &
                   ! + ((dez(l) - wtmp) + udotc * dez(l))*forcez)/(cssq*rhophi_loc)
		           fpost=feq + (ONE-omega_loc)*p(l)*fneq1 !+ HALF*(F_discr)	
		           press_out(i,j,k)=press_out(i,j,k) + fpost
		           u_out(i,j,k)=u_out(i,j,k) + fpost*dex(l)
	       	       v_out(i,j,k)=v_out(i,j,k) + fpost*dey(l)
		           w_out(i,j,k)=w_out(i,j,k) + fpost*dez(l)
		           pxx_out(i,j,k)=pxx_out(i,j,k) + fpost*dex(l)*dex(l)
                   pyy_out(i,j,k)=pyy_out(i,j,k) + fpost*dey(l)*dey(l)
                   pzz_out(i,j,k)=pzz_out(i,j,k) + fpost*dez(l)*dez(l)
                   pxy_out(i,j,k)=pxy_out(i,j,k) + fpost*dex(l)*dey(l)
                   pxz_out(i,j,k)=pxz_out(i,j,k) + fpost*dex(l)*dez(l)
                   pyz_out(i,j,k)=pyz_out(i,j,k) + fpost*dey(l)*dez(l)
                 enddo
                 
           
                 
			 endif
	       enddo
	     enddo
#ifdef ACCNOKERNELS
         !$acc end parallel loop
#else
         !$acc end kernels
#endif   
       endif

	   gi=lx-1
       subchords(1)=(gi-1)/nx
	   if(subchords(1)==coords(1))then
#ifdef ACCNOKERNELS
		 !$acc parallel loop collapse(2) independent present(press_in,u_in,v_in,w_in,pxx_in,pyy_in,pzz_in,pxy_in,pxz_in,pyz_in &
#ifdef TWOCOMPONENT	
         !$acc& ,phifields_s &
#endif			 
		 !$acc& ) private(i,j,k,l,gi,gj,gk)
#else
		 !$acc kernels present(press_in,u_in,v_in,w_in,pxx_in,pyy_in,pzz_in,pxy_in,pxz_in,pyz_in &
#ifdef TWOCOMPONENT	
         !$acc& ,phifields_s &
#endif			 
		 !$acc& )
		 !$acc loop collapse(2) independent private(i,j,k,l,gi,gj,gk)
#endif	
         do k=1,nz
		   do j=1,ny
		 	 gj=ny*coords(2)+j
			 gk=nz*coords(3)+k
			 if(gj>1 .and. gj<ly .and. gk>1 .and. gk<lz)then
	           i=nx-1
               gi=nx*coords(1)+i
	           if(isfluid(i,j,k) .ne. -1)cycle
	           
	             xblock=(i+2*TILE_DIMx-1)/TILE_DIMx   
                 yblock=(j+2*TILE_DIMy-1)/TILE_DIMy     
                 zblock=(k+2*TILE_DIMz-1)/TILE_DIMz  
                 myblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
                 ii=i-xblock*TILE_DIMx+2*TILE_DIMx
                 jj=j-yblock*TILE_DIMy+2*TILE_DIMy
                 kk=k-zblock*TILE_DIMz+2*TILE_DIMz 
				 
				 pxx=pxx_in(i,j,k)
				 pyy=pyy_in(i,j,k)
				 pzz=pzz_in(i,j,k)
				 pxy=pxy_in(i,j,k)
				 pxz=pxz_in(i,j,k)
				 pyz=pyz_in(i,j,k)
				 
				 
				 presstmp=press
				 utmp=0.0_db !
				 vtmp=0.0_db !
				 wtmp=0.0_db !
				 
	             pxx=pxx - cssq*press_in(i,j,k) - u_in(i,j,k)*u_in(i,j,k) 
	             pyy=pyy - cssq*press_in(i,j,k) - v_in(i,j,k)*v_in(i,j,k) 
	             pzz=pzz - cssq*press_in(i,j,k) - w_in(i,j,k)*w_in(i,j,k) 
	             pxy=pxy - u_in(i,j,k)*v_in(i,j,k)
	             pxz=pxz - u_in(i,j,k)*w_in(i,j,k)
	             pyz=pyz - v_in(i,j,k)*w_in(i,j,k)
#ifdef TWOCOMPONENT	 
                 phi_loc=phifields_s(idx5(ii,jj,kk,1,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields))
#endif	           
#ifdef DENSRATIO
			     rhophi_loc = rho_r*phi_loc+(ONE-phi_loc)*rho_b 
#else
		         rhophi_loc = ONE 
#endif

#ifdef TWOCOMPONENT

			     visc_loc=(rho_r*visc1*phi_loc+(1.0_db-phi_loc)*visc2*rho_b)/rhophi_loc

						
			     tau_loc=(visc_loc*invcssq + HALF) !è una tau
						
			     omega_loc=ONE/tau_loc !è una omega
						
#else
			     omega_loc=omega
#endif		
                 uu=HALF*(utmp*utmp + vtmp*vtmp+ wtmp*wtmp)*invcssq
				 !$acc loop seq
				 do l=1,nlinks
		           lopp=opp(l)
		           iii=i+ex(lopp)
		           jjj=j+ey(lopp)
		           kkk=k+ez(lopp)
		           if(isfluid(iii,jjj,kkk).ne.0) cycle 
		           feq=p(l)*press_in(i,j,k)
		           fneq1=(HALF/(cssq*cssq))*( (dex(l)*dex(l)-cssq)*pxx &
		            + (dey(l)*dey(l)-cssq)*pyy + (dez(l)*dez(l)-cssq)*pzz &
	                + TWO*(dex(l)*dey(l))*pxy + TWO*(dex(l)*dez(l))*pxz &
		            + TWO*(dey(l)*dez(l))*pyz)
                   !F_discr = p(l)*(dex(l)*forcex &
                   ! + dey(l)*forcey &
                   ! + dez(l)*forcez)/(cssq*rhophi_loc)
		           fpost=feq + (ONE-omega_loc)*p(l)*fneq1 !+ HALF*(F_discr)	
		           press_out(i,j,k)=press_out(i,j,k) - fpost
		           u_out(i,j,k)=u_out(i,j,k) - fpost*dex(l)
		           v_out(i,j,k)=v_out(i,j,k) - fpost*dey(l)
		           w_out(i,j,k)=w_out(i,j,k) - fpost*dez(l)
		           pxx_out(i,j,k)=pxx_out(i,j,k) - fpost*dex(l)*dex(l)
                   pyy_out(i,j,k)=pyy_out(i,j,k) - fpost*dey(l)*dey(l)
                   pzz_out(i,j,k)=pzz_out(i,j,k) - fpost*dez(l)*dez(l)
                   pxy_out(i,j,k)=pxy_out(i,j,k) - fpost*dex(l)*dey(l)
                   pxz_out(i,j,k)=pxz_out(i,j,k) - fpost*dex(l)*dez(l)
                   pyz_out(i,j,k)=pyz_out(i,j,k) - fpost*dey(l)*dez(l)	
		           udotc=(utmp*dex(l) + vtmp*dey(l)+ wtmp*dez(l))*invcssq
		           feq=p(l)*(presstmp + udotc+ HALF*udotc*udotc - uu)
		           fneq1=(HALF/(cssq*cssq))*( (dex(l)*dex(l)-cssq)*pxx &
		            + (dey(l)*dey(l)-cssq)*pyy + (dez(l)*dez(l)-cssq)*pzz &
	                + TWO*(dex(l)*dey(l))*pxy + TWO*(dex(l)*dez(l))*pxz &
		            + TWO*(dey(l)*dez(l))*pyz)
                   ! F_discr = p(l)*(((dex(l) - utmp) + udotc * dex(l))*forcex &
                   ! + ((dey(l) - vtmp) + udotc * dey(l))*forcey &
                   ! + ((dez(l) - wtmp) + udotc * dez(l))*forcez)/(cssq*rhophi_loc)
		           fpost=feq + (ONE-omega_loc)*p(l)*fneq1 !+ HALF*(F_discr)	
		           press_out(i,j,k)=press_out(i,j,k) + fpost
		           u_out(i,j,k)=u_out(i,j,k) + fpost*dex(l)
	       	       v_out(i,j,k)=v_out(i,j,k) + fpost*dey(l)
		           w_out(i,j,k)=w_out(i,j,k) + fpost*dez(l)
		           pxx_out(i,j,k)=pxx_out(i,j,k) + fpost*dex(l)*dex(l)
                   pyy_out(i,j,k)=pyy_out(i,j,k) + fpost*dey(l)*dey(l)
                   pzz_out(i,j,k)=pzz_out(i,j,k) + fpost*dez(l)*dez(l)
                   pxy_out(i,j,k)=pxy_out(i,j,k) + fpost*dex(l)*dey(l)
                   pxz_out(i,j,k)=pxz_out(i,j,k) + fpost*dex(l)*dez(l)
                   pyz_out(i,j,k)=pyz_out(i,j,k) + fpost*dey(l)*dez(l)
                 enddo
                 
              
                 
			 endif
	       enddo
	     enddo
#ifdef ACCNOKERNELS
         !$acc end parallel loop
#else
         !$acc end kernels
#endif
       endif
	 endif	
	   
	 if(pbc_y.eq.0)then
	   gj=2
       subchords(2)=(gj-1)/ny
	   if(subchords(2)==coords(2))then
#ifdef ACCNOKERNELS
	     !$acc parallel loop collapse(2) independent present(press_in,u_in,v_in,w_in,pxx_in,pyy_in,pzz_in,pxy_in,pxz_in,pyz_in &
#ifdef TWOCOMPONENT	
         !$acc& ,phifields_s &
#endif		
	     !$acc& )     
#else
	     !$acc kernels present(press_in,u_in,v_in,w_in,pxx_in,pyy_in,pzz_in,pxy_in,pxz_in,pyz_in &
#ifdef TWOCOMPONENT	
         !$acc& ,phifields_s &
#endif		     
	     !$acc& )
	     !$acc loop collapse(2) independent private(i,j,k,l,gi,gj,gk)
#endif	    
	     do k=1,nz
		   do i=1,nx
			 gi=nx*coords(1)+i
			 gk=nz*coords(3)+k
			 if(gi>1 .and. gi<lx .and. gk>1 .and. gk<lz)then
	           j=2
               gj=ny*coords(2)+j
               if(isfluid(i,j,k) .ne. -1)cycle
	           
	             xblock=(i+2*TILE_DIMx-1)/TILE_DIMx   
                 yblock=(j+2*TILE_DIMy-1)/TILE_DIMy     
                 zblock=(k+2*TILE_DIMz-1)/TILE_DIMz  
                 myblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
                 ii=i-xblock*TILE_DIMx+2*TILE_DIMx
                 jj=j-yblock*TILE_DIMy+2*TILE_DIMy
                 kk=k-zblock*TILE_DIMz+2*TILE_DIMz 
				 
				 pxx=pxx_in(i,j,k)
				 pyy=pyy_in(i,j,k)
				 pzz=pzz_in(i,j,k)
				 pxy=pxy_in(i,j,k)
				 pxz=pxz_in(i,j,k)
				 pyz=pyz_in(i,j,k)
				 
				 
				 presstmp=press
				 utmp=0.0_db !
				 vtmp=0.0_db !
				 wtmp=0.0_db !
				 
	             pxx=pxx - cssq*press_in(i,j,k) - u_in(i,j,k)*u_in(i,j,k) 
	             pyy=pyy - cssq*press_in(i,j,k) - v_in(i,j,k)*v_in(i,j,k) 
	             pzz=pzz - cssq*press_in(i,j,k) - w_in(i,j,k)*w_in(i,j,k) 
	             pxy=pxy - u_in(i,j,k)*v_in(i,j,k)
	             pxz=pxz - u_in(i,j,k)*w_in(i,j,k)
	             pyz=pyz - v_in(i,j,k)*w_in(i,j,k)

#ifdef TWOCOMPONENT	 
                 phi_loc=phifields_s(idx5(ii,jj,kk,1,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields))
#endif	           
#ifdef DENSRATIO
			     rhophi_loc = rho_r*phi_loc+(ONE-phi_loc)*rho_b 
#else
		         rhophi_loc = ONE 
#endif

#ifdef TWOCOMPONENT

			     visc_loc=(rho_r*visc1*phi_loc+(1.0_db-phi_loc)*visc2*rho_b)/rhophi_loc

						
			     tau_loc=(visc_loc*invcssq + HALF) !è una tau
						
			     omega_loc=ONE/tau_loc !è una omega
						
#else
			     omega_loc=omega
#endif		
                 uu=HALF*(utmp*utmp + vtmp*vtmp+ wtmp*wtmp)*invcssq
				 !$acc loop seq
				 do l=1,nlinks
		           lopp=opp(l)
		           iii=i+ex(lopp)
		           jjj=j+ey(lopp)
		           kkk=k+ez(lopp)
		           if(isfluid(iii,jjj,kkk).ne.0) cycle 
		           feq=p(l)*press_in(i,j,k)
		           fneq1=(HALF/(cssq*cssq))*( (dex(l)*dex(l)-cssq)*pxx &
		            + (dey(l)*dey(l)-cssq)*pyy + (dez(l)*dez(l)-cssq)*pzz &
	                + TWO*(dex(l)*dey(l))*pxy + TWO*(dex(l)*dez(l))*pxz &
		            + TWO*(dey(l)*dez(l))*pyz)
                   !F_discr = p(l)*(dex(l)*forcex &
                   ! + dey(l)*forcey &
                   ! + dez(l)*forcez)/(cssq*rhophi_loc)
		           fpost=feq + (ONE-omega_loc)*p(l)*fneq1 !+ HALF*(F_discr)	
		           press_out(i,j,k)=press_out(i,j,k) - fpost
		           u_out(i,j,k)=u_out(i,j,k) - fpost*dex(l)
		           v_out(i,j,k)=v_out(i,j,k) - fpost*dey(l)
		           w_out(i,j,k)=w_out(i,j,k) - fpost*dez(l)
		           pxx_out(i,j,k)=pxx_out(i,j,k) - fpost*dex(l)*dex(l)
                   pyy_out(i,j,k)=pyy_out(i,j,k) - fpost*dey(l)*dey(l)
                   pzz_out(i,j,k)=pzz_out(i,j,k) - fpost*dez(l)*dez(l)
                   pxy_out(i,j,k)=pxy_out(i,j,k) - fpost*dex(l)*dey(l)
                   pxz_out(i,j,k)=pxz_out(i,j,k) - fpost*dex(l)*dez(l)
                   pyz_out(i,j,k)=pyz_out(i,j,k) - fpost*dey(l)*dez(l)	
		           udotc=(utmp*dex(l) + vtmp*dey(l)+ wtmp*dez(l))*invcssq
		           feq=p(l)*(presstmp + udotc+ HALF*udotc*udotc - uu)
		           fneq1=(HALF/(cssq*cssq))*( (dex(l)*dex(l)-cssq)*pxx &
		            + (dey(l)*dey(l)-cssq)*pyy + (dez(l)*dez(l)-cssq)*pzz &
	                + TWO*(dex(l)*dey(l))*pxy + TWO*(dex(l)*dez(l))*pxz &
		            + TWO*(dey(l)*dez(l))*pyz)
                   ! F_discr = p(l)*(((dex(l) - utmp) + udotc * dex(l))*forcex &
                   ! + ((dey(l) - vtmp) + udotc * dey(l))*forcey &
                   ! + ((dez(l) - wtmp) + udotc * dez(l))*forcez)/(cssq*rhophi_loc)
		           fpost=feq + (ONE-omega_loc)*p(l)*fneq1 !+ HALF*(F_discr)	
		           press_out(i,j,k)=press_out(i,j,k) + fpost
		           u_out(i,j,k)=u_out(i,j,k) + fpost*dex(l)
	       	       v_out(i,j,k)=v_out(i,j,k) + fpost*dey(l)
		           w_out(i,j,k)=w_out(i,j,k) + fpost*dez(l)
		           pxx_out(i,j,k)=pxx_out(i,j,k) + fpost*dex(l)*dex(l)
                   pyy_out(i,j,k)=pyy_out(i,j,k) + fpost*dey(l)*dey(l)
                   pzz_out(i,j,k)=pzz_out(i,j,k) + fpost*dez(l)*dez(l)
                   pxy_out(i,j,k)=pxy_out(i,j,k) + fpost*dex(l)*dey(l)
                   pxz_out(i,j,k)=pxz_out(i,j,k) + fpost*dex(l)*dez(l)
                   pyz_out(i,j,k)=pyz_out(i,j,k) + fpost*dey(l)*dez(l)
                 enddo
                 
            
                 
			 endif
	       enddo
	     enddo
#ifdef ACCNOKERNELS
         !$acc end parallel loop
#else
         !$acc end kernels
#endif
       endif
       
	   gj=ly-1
       subchords(2)=(gj-1)/ny
	   if(subchords(2)==coords(2))then
#ifdef ACCNOKERNELS
	     !$acc parallel loop collapse(2) independent present(press_in,u_in,v_in,w_in,pxx_in,pyy_in,pzz_in,pxy_in,pxz_in,pyz_in &
#ifdef TWOCOMPONENT	
         !$acc& ,phifields_s &
#endif		     
         !$acc& )
#else
	     !$acc kernels present(press_in,u_in,v_in,w_in,pxx_in,pyy_in,pzz_in,pxy_in,pxz_in,pyz_in &
#ifdef TWOCOMPONENT	
         !$acc& ,phifields_s &
#endif		     
	     !$acc& )
	     !$acc loop collapse(2) independent private(i,j,k,l,gi,gj,gk)
#endif	    
	     do k=1,nz
		   do i=1,nx
			 gi=nx*coords(1)+i
			 gk=nz*coords(3)+k
			 if(gi>1 .and. gi<lx .and. gk>1 .and. gk<lz)then  
			   j=ny-1
			   gj=ny*coords(2)+j
               if(isfluid(i,j,k) .ne. -1)cycle
	           
	             xblock=(i+2*TILE_DIMx-1)/TILE_DIMx   
                 yblock=(j+2*TILE_DIMy-1)/TILE_DIMy     
                 zblock=(k+2*TILE_DIMz-1)/TILE_DIMz  
                 myblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
                 ii=i-xblock*TILE_DIMx+2*TILE_DIMx
                 jj=j-yblock*TILE_DIMy+2*TILE_DIMy
                 kk=k-zblock*TILE_DIMz+2*TILE_DIMz 
				 
				 pxx=pxx_in(i,j,k)
				 pyy=pyy_in(i,j,k)
				 pzz=pzz_in(i,j,k)
				 pxy=pxy_in(i,j,k)
				 pxz=pxz_in(i,j,k)
				 pyz=pyz_in(i,j,k)
				 
				 
				 presstmp=press
				 utmp=0.0_db !
				 vtmp=0.0_db !
				 wtmp=0.0_db !
				 
	             pxx=pxx - cssq*press_in(i,j,k) - u_in(i,j,k)*u_in(i,j,k) 
	             pyy=pyy - cssq*press_in(i,j,k) - v_in(i,j,k)*v_in(i,j,k) 
	             pzz=pzz - cssq*press_in(i,j,k) - w_in(i,j,k)*w_in(i,j,k) 
	             pxy=pxy - u_in(i,j,k)*v_in(i,j,k)
	             pxz=pxz - u_in(i,j,k)*w_in(i,j,k)
	             pyz=pyz - v_in(i,j,k)*w_in(i,j,k)
	             
#ifdef TWOCOMPONENT	 
                 phi_loc=phifields_s(idx5(ii,jj,kk,1,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields))
#endif	           
#ifdef DENSRATIO
			     rhophi_loc = rho_r*phi_loc+(ONE-phi_loc)*rho_b 
#else
		         rhophi_loc = ONE 
#endif

#ifdef TWOCOMPONENT

			     visc_loc=(rho_r*visc1*phi_loc+(1.0_db-phi_loc)*visc2*rho_b)/rhophi_loc

						
			     tau_loc=(visc_loc*invcssq + HALF) !è una tau
						
			     omega_loc=ONE/tau_loc !è una omega
						
#else
			     omega_loc=omega
#endif		
                 uu=HALF*(utmp*utmp + vtmp*vtmp+ wtmp*wtmp)*invcssq
				 !$acc loop seq
				 do l=1,nlinks
		           lopp=opp(l)
		           iii=i+ex(lopp)
		           jjj=j+ey(lopp)
		           kkk=k+ez(lopp)
		           if(isfluid(iii,jjj,kkk).ne.0) cycle 
		           feq=p(l)*press_in(i,j,k)
		           fneq1=(HALF/(cssq*cssq))*( (dex(l)*dex(l)-cssq)*pxx &
		            + (dey(l)*dey(l)-cssq)*pyy + (dez(l)*dez(l)-cssq)*pzz &
	                + TWO*(dex(l)*dey(l))*pxy + TWO*(dex(l)*dez(l))*pxz &
		            + TWO*(dey(l)*dez(l))*pyz)
                   !F_discr = p(l)*(dex(l)*forcex &
                   ! + dey(l)*forcey &
                   ! + dez(l)*forcez)/(cssq*rhophi_loc)
		           fpost=feq + (ONE-omega_loc)*p(l)*fneq1 !+ HALF*(F_discr)	
		           press_out(i,j,k)=press_out(i,j,k) - fpost
		           u_out(i,j,k)=u_out(i,j,k) - fpost*dex(l)
		           v_out(i,j,k)=v_out(i,j,k) - fpost*dey(l)
		           w_out(i,j,k)=w_out(i,j,k) - fpost*dez(l)
		           pxx_out(i,j,k)=pxx_out(i,j,k) - fpost*dex(l)*dex(l)
                   pyy_out(i,j,k)=pyy_out(i,j,k) - fpost*dey(l)*dey(l)
                   pzz_out(i,j,k)=pzz_out(i,j,k) - fpost*dez(l)*dez(l)
                   pxy_out(i,j,k)=pxy_out(i,j,k) - fpost*dex(l)*dey(l)
                   pxz_out(i,j,k)=pxz_out(i,j,k) - fpost*dex(l)*dez(l)
                   pyz_out(i,j,k)=pyz_out(i,j,k) - fpost*dey(l)*dez(l)	
		           udotc=(utmp*dex(l) + vtmp*dey(l)+ wtmp*dez(l))*invcssq
		           feq=p(l)*(presstmp + udotc+ HALF*udotc*udotc - uu)
		           fneq1=(HALF/(cssq*cssq))*( (dex(l)*dex(l)-cssq)*pxx &
		            + (dey(l)*dey(l)-cssq)*pyy + (dez(l)*dez(l)-cssq)*pzz &
	                + TWO*(dex(l)*dey(l))*pxy + TWO*(dex(l)*dez(l))*pxz &
		            + TWO*(dey(l)*dez(l))*pyz)
                   ! F_discr = p(l)*(((dex(l) - utmp) + udotc * dex(l))*forcex &
                   ! + ((dey(l) - vtmp) + udotc * dey(l))*forcey &
                   ! + ((dez(l) - wtmp) + udotc * dez(l))*forcez)/(cssq*rhophi_loc)
		           fpost=feq + (ONE-omega_loc)*p(l)*fneq1 !+ HALF*(F_discr)	
		           press_out(i,j,k)=press_out(i,j,k) + fpost
		           u_out(i,j,k)=u_out(i,j,k) + fpost*dex(l)
	       	       v_out(i,j,k)=v_out(i,j,k) + fpost*dey(l)
		           w_out(i,j,k)=w_out(i,j,k) + fpost*dez(l)
		           pxx_out(i,j,k)=pxx_out(i,j,k) + fpost*dex(l)*dex(l)
                   pyy_out(i,j,k)=pyy_out(i,j,k) + fpost*dey(l)*dey(l)
                   pzz_out(i,j,k)=pzz_out(i,j,k) + fpost*dez(l)*dez(l)
                   pxy_out(i,j,k)=pxy_out(i,j,k) + fpost*dex(l)*dey(l)
                   pxz_out(i,j,k)=pxz_out(i,j,k) + fpost*dex(l)*dez(l)
                   pyz_out(i,j,k)=pyz_out(i,j,k) + fpost*dey(l)*dez(l)
                 enddo
                 
           
                 
			 endif
	       enddo
	     enddo
#ifdef ACCNOKERNELS
         !$acc end parallel loop
#else
         !$acc end kernels
#endif
	   endif
	 endif

	 if(pbc_z.eq.0)then
	   gk=2
       subchords(3)=(gk-1)/nz
	   if(subchords(3)==coords(3))then
#ifdef ACCNOKERNELS
	     !$acc parallel loop collapse(2) independent present(press_in,u_in,v_in,w_in,pxx_in,pyy_in,pzz_in,pxy_in,pxz_in,pyz_in &
#ifdef TWOCOMPONENT	
         !$acc& ,phifields_s &
#endif		     
	     !$acc& ) private(i,j,k,l,gi,gj,gk)
#else
	     !$acc kernels present(press_in,u_in,v_in,w_in,pxx_in,pyy_in,pzz_in,pxy_in,pxz_in,pyz_in &
#ifdef TWOCOMPONENT	
         !$acc& ,phifields_s &
#endif		     
	     !$acc& )
	     !$acc loop collapse(2) independent private(i,j,k,l,gi,gj,gk)
#endif	    
	     do j=1,ny
		   do i=1,nx
			 gi=nx*coords(1)+i
			 gj=ny*coords(2)+j
			 if(gi>1 .and. gi<lx .and. gj>1 .and. gj<ly)then
	           k=2
               gk=nz*coords(3)+k
               if(isfluid(i,j,k) .ne. -1)cycle
	           
	             xblock=(i+2*TILE_DIMx-1)/TILE_DIMx   
                 yblock=(j+2*TILE_DIMy-1)/TILE_DIMy     
                 zblock=(k+2*TILE_DIMz-1)/TILE_DIMz  
                 myblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
                 ii=i-xblock*TILE_DIMx+2*TILE_DIMx
                 jj=j-yblock*TILE_DIMy+2*TILE_DIMy
                 kk=k-zblock*TILE_DIMz+2*TILE_DIMz 
				 
				 pxx=pxx_in(i,j,k)
				 pyy=pyy_in(i,j,k)
				 pzz=pzz_in(i,j,k)
				 pxy=pxy_in(i,j,k)
				 pxz=pxz_in(i,j,k)
				 pyz=pyz_in(i,j,k)
				 
				 
				 presstmp=press
				 utmp=0.0_db !
				 vtmp=0.0_db !
				 wtmp=0.0_db !
				 
	             pxx=pxx - cssq*press_in(i,j,k) - u_in(i,j,k)*u_in(i,j,k) 
	             pyy=pyy - cssq*press_in(i,j,k) - v_in(i,j,k)*v_in(i,j,k) 
	             pzz=pzz - cssq*press_in(i,j,k) - w_in(i,j,k)*w_in(i,j,k) 
	             pxy=pxy - u_in(i,j,k)*v_in(i,j,k)
	             pxz=pxz - u_in(i,j,k)*w_in(i,j,k)
	             pyz=pyz - v_in(i,j,k)*w_in(i,j,k)

#ifdef TWOCOMPONENT	 
                 phi_loc=phifields_s(idx5(ii,jj,kk,1,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields))
#endif	           
#ifdef DENSRATIO
			     rhophi_loc = rho_r*phi_loc+(ONE-phi_loc)*rho_b 
#else
		         rhophi_loc = ONE 
#endif

#ifdef TWOCOMPONENT

			     visc_loc=(rho_r*visc1*phi_loc+(1.0_db-phi_loc)*visc2*rho_b)/rhophi_loc

						
			     tau_loc=(visc_loc*invcssq + HALF) !è una tau
						
			     omega_loc=ONE/tau_loc !è una omega
						
#else
			     omega_loc=omega
#endif		
                 uu=HALF*(utmp*utmp + vtmp*vtmp+ wtmp*wtmp)*invcssq
				 !$acc loop seq
				 do l=1,nlinks
		           lopp=opp(l)
		           iii=i+ex(lopp)
		           jjj=j+ey(lopp)
		           kkk=k+ez(lopp)
		           if(isfluid(iii,jjj,kkk).ne.0) cycle 
		           feq=p(l)*press_in(i,j,k)
		           fneq1=(HALF/(cssq*cssq))*( (dex(l)*dex(l)-cssq)*pxx &
		            + (dey(l)*dey(l)-cssq)*pyy + (dez(l)*dez(l)-cssq)*pzz &
	                + TWO*(dex(l)*dey(l))*pxy + TWO*(dex(l)*dez(l))*pxz &
		            + TWO*(dey(l)*dez(l))*pyz)
                   !F_discr = p(l)*(dex(l)*forcex &
                   ! + dey(l)*forcey &
                   ! + dez(l)*forcez)/(cssq*rhophi_loc)
		           fpost=feq + (ONE-omega_loc)*p(l)*fneq1 !+ HALF*(F_discr)	
		           press_out(i,j,k)=press_out(i,j,k) - fpost
		           u_out(i,j,k)=u_out(i,j,k) - fpost*dex(l)
		           v_out(i,j,k)=v_out(i,j,k) - fpost*dey(l)
		           w_out(i,j,k)=w_out(i,j,k) - fpost*dez(l)
		           pxx_out(i,j,k)=pxx_out(i,j,k) - fpost*dex(l)*dex(l)
                   pyy_out(i,j,k)=pyy_out(i,j,k) - fpost*dey(l)*dey(l)
                   pzz_out(i,j,k)=pzz_out(i,j,k) - fpost*dez(l)*dez(l)
                   pxy_out(i,j,k)=pxy_out(i,j,k) - fpost*dex(l)*dey(l)
                   pxz_out(i,j,k)=pxz_out(i,j,k) - fpost*dex(l)*dez(l)
                   pyz_out(i,j,k)=pyz_out(i,j,k) - fpost*dey(l)*dez(l)	
		           udotc=(utmp*dex(l) + vtmp*dey(l)+ wtmp*dez(l))*invcssq
		           feq=p(l)*(presstmp + udotc+ HALF*udotc*udotc - uu)
		           fneq1=(HALF/(cssq*cssq))*( (dex(l)*dex(l)-cssq)*pxx &
		            + (dey(l)*dey(l)-cssq)*pyy + (dez(l)*dez(l)-cssq)*pzz &
	                + TWO*(dex(l)*dey(l))*pxy + TWO*(dex(l)*dez(l))*pxz &
		            + TWO*(dey(l)*dez(l))*pyz)
                   ! F_discr = p(l)*(((dex(l) - utmp) + udotc * dex(l))*forcex &
                   ! + ((dey(l) - vtmp) + udotc * dey(l))*forcey &
                   ! + ((dez(l) - wtmp) + udotc * dez(l))*forcez)/(cssq*rhophi_loc)
		           fpost=feq + (ONE-omega_loc)*p(l)*fneq1 !+ HALF*(F_discr)	
		           press_out(i,j,k)=press_out(i,j,k) + fpost
		           u_out(i,j,k)=u_out(i,j,k) + fpost*dex(l)
	       	       v_out(i,j,k)=v_out(i,j,k) + fpost*dey(l)
		           w_out(i,j,k)=w_out(i,j,k) + fpost*dez(l)
		           pxx_out(i,j,k)=pxx_out(i,j,k) + fpost*dex(l)*dex(l)
                   pyy_out(i,j,k)=pyy_out(i,j,k) + fpost*dey(l)*dey(l)
                   pzz_out(i,j,k)=pzz_out(i,j,k) + fpost*dez(l)*dez(l)
                   pxy_out(i,j,k)=pxy_out(i,j,k) + fpost*dex(l)*dey(l)
                   pxz_out(i,j,k)=pxz_out(i,j,k) + fpost*dex(l)*dez(l)
                   pyz_out(i,j,k)=pyz_out(i,j,k) + fpost*dey(l)*dez(l)
                 enddo
                 
    
                 
			 endif
	       enddo
	     enddo
#ifdef ACCNOKERNELS
         !$acc end parallel loop
#else
         !$acc end kernels
#endif
	   endif
	   
	   gk=lz-1
       subchords(3)=(gk-1)/nz
	   if(subchords(3)==coords(3))then			   
#ifdef ACCNOKERNELS
	     !$acc parallel loop collapse(2) independent present(press_in,u_in,v_in,w_in,pxx_in,pyy_in,pzz_in,pxy_in,pxz_in,pyz_in &
#ifdef TWOCOMPONENT	
         !$acc& ,phifields_s &
#endif		     
	     !$acc& ) private(i,j,k,l,gi,gj,gk)
#else
	     !$acc kernels present(press_in,u_in,v_in,w_in,pxx_in,pyy_in,pzz_in,pxy_in,pxz_in,pyz_in &
#ifdef TWOCOMPONENT	
         !$acc& ,phifields_s &
#endif		     
	     !$acc& )
	     !$acc loop collapse(2) independent private(i,j,k,l,gi,gj,gk)
#endif	    
	     do j=1,ny
		   do i=1,nx
			 gi=nx*coords(1)+i
			 gj=ny*coords(2)+j
			 if(gi>1 .and. gi<lx .and. gj>1 .and. gj<ly)then  
			   k=nz-1
               gk=nz*coords(3)+k
               if(isfluid(i,j,k) .ne. -1)cycle
	           
	             xblock=(i+2*TILE_DIMx-1)/TILE_DIMx   
                 yblock=(j+2*TILE_DIMy-1)/TILE_DIMy     
                 zblock=(k+2*TILE_DIMz-1)/TILE_DIMz  
                 myblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
                 ii=i-xblock*TILE_DIMx+2*TILE_DIMx
                 jj=j-yblock*TILE_DIMy+2*TILE_DIMy
                 kk=k-zblock*TILE_DIMz+2*TILE_DIMz 
				 
				 pxx=pxx_in(i,j,k)
				 pyy=pyy_in(i,j,k)
				 pzz=pzz_in(i,j,k)
				 pxy=pxy_in(i,j,k)
				 pxz=pxz_in(i,j,k)
				 pyz=pyz_in(i,j,k)
				 
				 
				 presstmp=press
				 utmp=0.0_db !
				 vtmp=0.0_db !
				 wtmp=0.0_db !
				 
	             pxx=pxx - cssq*press_in(i,j,k) - u_in(i,j,k)*u_in(i,j,k) 
	             pyy=pyy - cssq*press_in(i,j,k) - v_in(i,j,k)*v_in(i,j,k) 
	             pzz=pzz - cssq*press_in(i,j,k) - w_in(i,j,k)*w_in(i,j,k) 
	             pxy=pxy - u_in(i,j,k)*v_in(i,j,k)
	             pxz=pxz - u_in(i,j,k)*w_in(i,j,k)
	             pyz=pyz - v_in(i,j,k)*w_in(i,j,k)

#ifdef TWOCOMPONENT	 
                 phi_loc=phifields_s(idx5(ii,jj,kk,1,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields))
#endif	           
#ifdef DENSRATIO
			     rhophi_loc = rho_r*phi_loc+(ONE-phi_loc)*rho_b 
#else
		         rhophi_loc = ONE 
#endif

#ifdef TWOCOMPONENT

			     visc_loc=(rho_r*visc1*phi_loc+(1.0_db-phi_loc)*visc2*rho_b)/rhophi_loc

						
			     tau_loc=(visc_loc*invcssq + HALF) !è una tau
						
			     omega_loc=ONE/tau_loc !è una omega
						
#else
			     omega_loc=omega
#endif		

                 uu=HALF*(utmp*utmp + vtmp*vtmp+ wtmp*wtmp)*invcssq
				 !$acc loop seq
				 do l=1,nlinks
		           lopp=opp(l)
		           iii=i+ex(lopp)
		           jjj=j+ey(lopp)
		           kkk=k+ez(lopp)
		           if(isfluid(iii,jjj,kkk).ne.0) cycle 
		           feq=p(l)*press_in(i,j,k)
		           fneq1=(HALF/(cssq*cssq))*( (dex(l)*dex(l)-cssq)*pxx &
		            + (dey(l)*dey(l)-cssq)*pyy + (dez(l)*dez(l)-cssq)*pzz &
	                + TWO*(dex(l)*dey(l))*pxy + TWO*(dex(l)*dez(l))*pxz &
		            + TWO*(dey(l)*dez(l))*pyz)
                   !F_discr = p(l)*(dex(l)*forcex &
                   ! + dey(l)*forcey &
                   ! + dez(l)*forcez)/(cssq*rhophi_loc)
		           fpost=feq + (ONE-omega_loc)*p(l)*fneq1 !+ HALF*(F_discr)	
		           press_out(i,j,k)=press_out(i,j,k) - fpost
		           u_out(i,j,k)=u_out(i,j,k) - fpost*dex(l)
		           v_out(i,j,k)=v_out(i,j,k) - fpost*dey(l)
		           w_out(i,j,k)=w_out(i,j,k) - fpost*dez(l)
		           pxx_out(i,j,k)=pxx_out(i,j,k) - fpost*dex(l)*dex(l)
                   pyy_out(i,j,k)=pyy_out(i,j,k) - fpost*dey(l)*dey(l)
                   pzz_out(i,j,k)=pzz_out(i,j,k) - fpost*dez(l)*dez(l)
                   pxy_out(i,j,k)=pxy_out(i,j,k) - fpost*dex(l)*dey(l)
                   pxz_out(i,j,k)=pxz_out(i,j,k) - fpost*dex(l)*dez(l)
                   pyz_out(i,j,k)=pyz_out(i,j,k) - fpost*dey(l)*dez(l)	
		           udotc=(utmp*dex(l) + vtmp*dey(l)+ wtmp*dez(l))*invcssq
		           feq=p(l)*(presstmp + udotc+ HALF*udotc*udotc - uu)
		           fneq1=(HALF/(cssq*cssq))*( (dex(l)*dex(l)-cssq)*pxx &
		            + (dey(l)*dey(l)-cssq)*pyy + (dez(l)*dez(l)-cssq)*pzz &
	                + TWO*(dex(l)*dey(l))*pxy + TWO*(dex(l)*dez(l))*pxz &
		            + TWO*(dey(l)*dez(l))*pyz)
                   ! F_discr = p(l)*(((dex(l) - utmp) + udotc * dex(l))*forcex &
                   ! + ((dey(l) - vtmp) + udotc * dey(l))*forcey &
                   ! + ((dez(l) - wtmp) + udotc * dez(l))*forcez)/(cssq*rhophi_loc)
		           fpost=feq + (ONE-omega_loc)*p(l)*fneq1 !+ HALF*(F_discr)	
		           press_out(i,j,k)=press_out(i,j,k) + fpost
		           u_out(i,j,k)=u_out(i,j,k) + fpost*dex(l)
	       	       v_out(i,j,k)=v_out(i,j,k) + fpost*dey(l)
		           w_out(i,j,k)=w_out(i,j,k) + fpost*dez(l)
		           pxx_out(i,j,k)=pxx_out(i,j,k) + fpost*dex(l)*dex(l)
                   pyy_out(i,j,k)=pyy_out(i,j,k) + fpost*dey(l)*dey(l)
                   pzz_out(i,j,k)=pzz_out(i,j,k) + fpost*dez(l)*dez(l)
                   pxy_out(i,j,k)=pxy_out(i,j,k) + fpost*dex(l)*dey(l)
                   pxz_out(i,j,k)=pxz_out(i,j,k) + fpost*dex(l)*dez(l)
                   pyz_out(i,j,k)=pyz_out(i,j,k) + fpost*dey(l)*dez(l)
                 enddo
                 
     
             
			 endif
	       enddo
	     enddo
#ifdef ACCNOKERNELS
         !$acc end parallel loop
#else
         !$acc end kernels
#endif
	   endif
	 endif
#endif


   endsubroutine bcs_mesoscopic_hfields

   subroutine bcs_mesoscopic_phifields(phifields_s, &
    press_s,u_s,v_s,w_s,pxx_s,pyy_s,pzz_s,pxy_s,pxz_s,pyz_s)

     implicit none
     real(kind=db), allocatable, dimension(:) :: phifields_s
     
     real(kind=db), allocatable, dimension(:,:,:) :: &
       press_s,u_s,v_s,w_s,pxx_s,pyy_s,pzz_s,pxy_s,pxz_s,pyz_s
     
     integer :: subchords(3)
     integer :: ii,jj,kk,l,lopp
	 real(kind=db) :: feq, fneq1,utmp,vtmp,wtmp,rhophi_loc,phi_loc
#ifdef TWOCOMPONENT	  
	 real(kind=db) :: phitemp
#endif
#if defined(PHASE_CHANGE) || defined(INTERNAL_OBSTACLES)
	 real(kind=db) :: visc_loc,omega_loc,tau_loc
#endif
#if defined(INTERNAL_OBSTACLES)
	 real(kind=db) :: wet_R,rhotemp,phiavg, wet_thresh_low, wet_thresh_high,grad_thresh,phi_adj,weight,correc
	 real(kind=db) :: phi_fluid,gradfix,gradfiy,grad_parallel,theta_rad,cot_theta,phi_ghost,dphi_dz,gradfiz
	 integer :: conter
	 logical :: found
#endif
     integer :: xblock,yblock,zblock,myblock
     integer :: oxblock,oyblock,ozblock,omyblock
     real(kind=db) :: pxx,pyy,pzz,pxy,pxz,pyz
     real(kind=db) :: opress,ou,ov,ow,opxx,opyy,opzz,opxy,opxz,opyz
     integer :: oi,oj,ok
     integer :: iii,jjj,kkk
     integer :: oii,ojj,okk


#if defined(INTERNAL_OBSTACLES)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  phase field bcs3D
#ifdef TWOCOMPONENT

     call PHI_int_boundary_cuda(phifields_s, &
      press_s,u_s,v_s,w_s,pxx_s,pyy_s,pzz_s,pxy_s,pxz_s,pyz_s)
	 
	 if(openbc==1)then
	 
	 if(pbc_x.eq.0)then
	   gi=1
       subchords(1)=(gi-1)/nx
	   if(subchords(1)==coords(1))then
#ifdef ACCNOKERNELS
		 !$acc parallel loop collapse(2) independent present(phifields_s &
		 !$acc& ) private(i,j,k,l,gi,gj,gk)
#else
		 !$acc kernels present(phifields_s &
		 !$acc& )
		 !$acc loop collapse(2) independent private(i,j,k,l,gi,gj,gk)
#endif	    
	     do k=1,nz
		   do j=1,ny
			 gj=ny*coords(2)+j
			 gk=nz*coords(3)+k
			 if(gj>1 .and. gj<ly .and. gk>1 .and. gk<lz)then
	           i=1
	           
	           oi=i+1
               oj=j
               ok=k
	           
	           xblock=(i+2*TILE_DIMx-1)/TILE_DIMx   
               yblock=(j+2*TILE_DIMy-1)/TILE_DIMy     
               zblock=(k+2*TILE_DIMz-1)/TILE_DIMz  
               myblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
               ii=i-xblock*TILE_DIMx+2*TILE_DIMx
               jj=j-yblock*TILE_DIMy+2*TILE_DIMy
               kk=k-zblock*TILE_DIMz+2*TILE_DIMz 
                  
               oxblock=(oi+2*TILE_DIMx-1)/TILE_DIMx   
               oyblock=(oj+2*TILE_DIMy-1)/TILE_DIMy     
               ozblock=(ok+2*TILE_DIMz-1)/TILE_DIMz 
               omyblock=(oxblock-1)+(oyblock-1)*nxblock+(ozblock-1)*nxyblock+1
               oii=oi-oxblock*TILE_DIMx+2*TILE_DIMx
               ojj=oj-oyblock*TILE_DIMy+2*TILE_DIMy
               okk=ok-ozblock*TILE_DIMz+2*TILE_DIMz
	           
	           phifields_s(idx5(ii,jj,kk,1,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields))= &
	            phifields_s(idx5(oii,ojj,okk,1,omyblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields))
			 endif
	       enddo
	     enddo
#ifdef ACCNOKERNELS
         !$acc end parallel loop
#else
         !$acc end kernels
#endif
	   endif
			  
			  
	   gi=lx
       subchords(1)=(gi-1)/nx
	   if(subchords(1)==coords(1))then
#ifdef ACCNOKERNELS
		 !$acc parallel loop collapse(2) independent present(phifields_s &
		 !$acc& ) private(i,j,k,l,gi,gj,gk)
#else
		 !$acc kernels present(phifields_s &
		 !$acc& )
		 !$acc loop collapse(2) independent private(i,j,k,l,gi,gj,gk)
#endif	    
	     do k=1,nz
		   do j=1,ny
			 gj=ny*coords(2)+j
			 gk=nz*coords(3)+k
			 if(gj>1 .and. gj<ly .and. gk>1 .and. gk<lz)then
			   i=nx

	           oi=i-1
               oj=j
               ok=k
	           
	           xblock=(i+2*TILE_DIMx-1)/TILE_DIMx   
               yblock=(j+2*TILE_DIMy-1)/TILE_DIMy     
               zblock=(k+2*TILE_DIMz-1)/TILE_DIMz  
               myblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
               ii=i-xblock*TILE_DIMx+2*TILE_DIMx
               jj=j-yblock*TILE_DIMy+2*TILE_DIMy
               kk=k-zblock*TILE_DIMz+2*TILE_DIMz 
                  
               oxblock=(oi+2*TILE_DIMx-1)/TILE_DIMx   
               oyblock=(oj+2*TILE_DIMy-1)/TILE_DIMy     
               ozblock=(ok+2*TILE_DIMz-1)/TILE_DIMz 
               omyblock=(oxblock-1)+(oyblock-1)*nxblock+(ozblock-1)*nxyblock+1
               oii=oi-oxblock*TILE_DIMx+2*TILE_DIMx
               ojj=oj-oyblock*TILE_DIMy+2*TILE_DIMy
               okk=ok-ozblock*TILE_DIMz+2*TILE_DIMz
	           
	           phifields_s(idx5(ii,jj,kk,1,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields))= &
	            phifields_s(idx5(oii,ojj,okk,1,omyblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields))

			 endif
	       enddo
	     enddo
#ifdef ACCNOKERNELS
         !$acc end parallel loop
#else
         !$acc end kernels
#endif
	   endif
	 endif	
	    
	 if(pbc_y.eq.0)then
	   gj=1
       subchords(2)=(gj-1)/ny
	   if(subchords(2)==coords(2))then
#ifdef ACCNOKERNELS
		 !$acc parallel loop collapse(2) independent present(phifields_s)
#else
		 !$acc kernels present(phifields_s &
		 !$acc& )
		 !$acc loop collapse(2) independent private(i,j,k,l,gi,gj,gk)
#endif	    
	     do k=1,nz
		   do i=1,nx
		     gi=nx*coords(1)+i
			 gk=nz*coords(3)+k
			 if(gi>1 .and. gi<lx .and. gk>1 .and. gk<lz)then
	           j=1

	           oi=i
               oj=j+1
               ok=k
	           
	           xblock=(i+2*TILE_DIMx-1)/TILE_DIMx   
               yblock=(j+2*TILE_DIMy-1)/TILE_DIMy     
               zblock=(k+2*TILE_DIMz-1)/TILE_DIMz  
               myblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
               ii=i-xblock*TILE_DIMx+2*TILE_DIMx
               jj=j-yblock*TILE_DIMy+2*TILE_DIMy
               kk=k-zblock*TILE_DIMz+2*TILE_DIMz 
                  
               oxblock=(oi+2*TILE_DIMx-1)/TILE_DIMx   
               oyblock=(oj+2*TILE_DIMy-1)/TILE_DIMy     
               ozblock=(ok+2*TILE_DIMz-1)/TILE_DIMz 
               omyblock=(oxblock-1)+(oyblock-1)*nxblock+(ozblock-1)*nxyblock+1
               oii=oi-oxblock*TILE_DIMx+2*TILE_DIMx
               ojj=oj-oyblock*TILE_DIMy+2*TILE_DIMy
               okk=ok-ozblock*TILE_DIMz+2*TILE_DIMz
	           
	           phifields_s(idx5(ii,jj,kk,1,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields))= &
	            phifields_s(idx5(oii,ojj,okk,1,omyblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields))

			 endif
	       enddo
	     enddo
#ifdef ACCNOKERNELS
         !$acc end parallel loop
#else
         !$acc end kernels
#endif			   
	   endif
	    
	   gj=ly
       subchords(2)=(gj-1)/ny
	   if(subchords(2)==coords(2))then
#ifdef ACCNOKERNELS
		 !$acc parallel loop collapse(2) independent present(phifields_s)
#else
		 !$acc kernels present(phifields_s &
		 !$acc& )
		 !$acc loop collapse(2) independent private(i,j,k,l,gi,gj,gk)
#endif	  
	     do k=1,nz
		   do i=1,nx
		     gi=nx*coords(1)+i
			 gk=nz*coords(3)+k
			 if(gi>1 .and. gi<lx .and. gk>1 .and. gk<lz)then  
			   j=ny

	           oi=i
               oj=j-1
               ok=k
	           
	           xblock=(i+2*TILE_DIMx-1)/TILE_DIMx   
               yblock=(j+2*TILE_DIMy-1)/TILE_DIMy     
               zblock=(k+2*TILE_DIMz-1)/TILE_DIMz  
               myblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
               ii=i-xblock*TILE_DIMx+2*TILE_DIMx
               jj=j-yblock*TILE_DIMy+2*TILE_DIMy
               kk=k-zblock*TILE_DIMz+2*TILE_DIMz 
                  
               oxblock=(oi+2*TILE_DIMx-1)/TILE_DIMx   
               oyblock=(oj+2*TILE_DIMy-1)/TILE_DIMy     
               ozblock=(ok+2*TILE_DIMz-1)/TILE_DIMz 
               omyblock=(oxblock-1)+(oyblock-1)*nxblock+(ozblock-1)*nxyblock+1
               oii=oi-oxblock*TILE_DIMx+2*TILE_DIMx
               ojj=oj-oyblock*TILE_DIMy+2*TILE_DIMy
               okk=ok-ozblock*TILE_DIMz+2*TILE_DIMz
	           
	           phifields_s(idx5(ii,jj,kk,1,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields))= &
	            phifields_s(idx5(oii,ojj,okk,1,omyblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields))

			 endif
	       enddo
	     enddo
#ifdef ACCNOKERNELS
         !$acc end parallel loop
#else
         !$acc end kernels
#endif
	   endif
	 endif	    

	 if(pbc_z.eq.0)then
	   gk=1
       subchords(3)=(gk-1)/nz
	   if(subchords(3)==coords(3))then
#ifdef ACCNOKERNELS
		 !$acc parallel loop collapse(2) independent present(phifields_s &
		 !$acc& ) private(i,j,k,l,gi,gj,gk)
#else
		 !$acc kernels present(phifields_s &
		 !$acc& )
		 !$acc loop collapse(2) independent private(i,j,k,l,gi,gj,gk)
#endif	    
	     do j=1,ny
		   do i=1,nx
			 gi=nx*coords(1)+i
			 gj=ny*coords(2)+j
			 if(gi>1 .and. gi<lx .and. gj>1 .and. gj<ly)then
	           k=1

	           oi=i
               oj=j
               ok=k+1
	           
	           xblock=(i+2*TILE_DIMx-1)/TILE_DIMx   
               yblock=(j+2*TILE_DIMy-1)/TILE_DIMy     
               zblock=(k+2*TILE_DIMz-1)/TILE_DIMz  
               myblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
               ii=i-xblock*TILE_DIMx+2*TILE_DIMx
               jj=j-yblock*TILE_DIMy+2*TILE_DIMy
               kk=k-zblock*TILE_DIMz+2*TILE_DIMz 
                  
               oxblock=(oi+2*TILE_DIMx-1)/TILE_DIMx   
               oyblock=(oj+2*TILE_DIMy-1)/TILE_DIMy     
               ozblock=(ok+2*TILE_DIMz-1)/TILE_DIMz 
               omyblock=(oxblock-1)+(oyblock-1)*nxblock+(ozblock-1)*nxyblock+1
               oii=oi-oxblock*TILE_DIMx+2*TILE_DIMx
               ojj=oj-oyblock*TILE_DIMy+2*TILE_DIMy
               okk=ok-ozblock*TILE_DIMz+2*TILE_DIMz
	           
	           phifields_s(idx5(ii,jj,kk,1,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields))= &
	            phifields_s(idx5(oii,ojj,okk,1,omyblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields))

			 endif
	       enddo
	     enddo
#ifdef ACCNOKERNELS
         !$acc end parallel loop
#else
         !$acc end kernels
#endif
	   endif
			   
	   gk=lz
       subchords(3)=(gk-1)/nz
	   if(subchords(3)==coords(3))then
#ifdef ACCNOKERNELS
		 !$acc parallel loop collapse(2) independent present(phifields_s &
		 !$acc& ) private(i,j,k,l,gi,gj,gk)
#else
		 !$acc kernels present(phifields_s &
		 !$acc& )
		 !$acc loop collapse(2) independent private(i,j,k,l,gi,gj,gk)
#endif	    
	     do j=1,ny
		   do i=1,nx
			 gi=nx*coords(1)+i
			 gj=ny*coords(2)+j
			 if(gi>1 .and. gi<lx .and. gj>1 .and. gj<ly)then
			   k=nz

	           oi=i
               oj=j
               ok=k-1
	           
	           xblock=(i+2*TILE_DIMx-1)/TILE_DIMx   
               yblock=(j+2*TILE_DIMy-1)/TILE_DIMy     
               zblock=(k+2*TILE_DIMz-1)/TILE_DIMz  
               myblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
               ii=i-xblock*TILE_DIMx+2*TILE_DIMx
               jj=j-yblock*TILE_DIMy+2*TILE_DIMy
               kk=k-zblock*TILE_DIMz+2*TILE_DIMz 
                  
               oxblock=(oi+2*TILE_DIMx-1)/TILE_DIMx   
               oyblock=(oj+2*TILE_DIMy-1)/TILE_DIMy     
               ozblock=(ok+2*TILE_DIMz-1)/TILE_DIMz 
               omyblock=(oxblock-1)+(oyblock-1)*nxblock+(ozblock-1)*nxyblock+1
               oii=oi-oxblock*TILE_DIMx+2*TILE_DIMx
               ojj=oj-oyblock*TILE_DIMy+2*TILE_DIMy
               okk=ok-ozblock*TILE_DIMz+2*TILE_DIMz
	           
	           phifields_s(idx5(ii,jj,kk,1,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields))= &
	            phifields_s(idx5(oii,ojj,okk,1,omyblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields))

			 endif
	       enddo
	     enddo
#ifdef ACCNOKERNELS
         !$acc end parallel loop
#else
         !$acc end kernels
#endif
	   endif
	 endif	
	    			

!******************lagrange multiplier
#ifdef BCPHIFLUX

	 if(pbc_x.eq.0)then
	   gi=2
       subchords(1)=(gi-1)/nx
	   if(subchords(1)==coords(1))then
		 !$acc parallel loop collapse(2) reduction(+:global_phi_change) present(phifields_s &
		 !$acc& ,u_s,v_s,w_s) private(i,j,k,l,gi,gj,gk &
		 !$acc& ,phitemp,conter,xblock,yblock,zblock,myblock,ii,jj,kk,phi_loc,u)  
	     do k=1,nz
		   do j=1,ny
			 gj=ny*coords(2)+j
			 gk=nz*coords(3)+k
	         i=2
	         gi=nx*coords(1)+i
	         if(isfluid(i,j,k) .ne. -1)cycle
	         xblock=(i+2*TILE_DIMx-1)/TILE_DIMx   
             yblock=(j+2*TILE_DIMy-1)/TILE_DIMy     
             zblock=(k+2*TILE_DIMz-1)/TILE_DIMz  
             myblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
             ii=i-xblock*TILE_DIMx+2*TILE_DIMx
             jj=j-yblock*TILE_DIMy+2*TILE_DIMy
             kk=k-zblock*TILE_DIMz+2*TILE_DIMz 
	         phi_loc=phifields_s(idx5(ii,jj,kk,1,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields))
	         u=u_s(i,j,k)
	         
			 global_phi_change = global_phi_change + u*phi_loc
	       enddo
	     enddo
         !$acc end parallel loop
	   endif
	   gi=lx-1
       subchords(1)=(gi-1)/nx
	   if(subchords(1)==coords(1))then
		 !$acc parallel loop collapse(2) reduction(+:global_phi_change) present(phifields_s &
		 !$acc& ,u_s,v_s,w_s) private(i,j,k,l,gi,gj,gk &
		 !$acc& ,phitemp,conter,xblock,yblock,zblock,myblock,ii,jj,kk,phi_loc,u)   
	     do k=1,nz
		   do j=1,ny
			 gj=ny*coords(2)+j
			 gk=nz*coords(3)+k
			 i=nx-1
			 gi=nx*coords(1)+i
			 if(isfluid(i,j,k) .ne. -1)cycle
			 xblock=(i+2*TILE_DIMx-1)/TILE_DIMx   
             yblock=(j+2*TILE_DIMy-1)/TILE_DIMy     
             zblock=(k+2*TILE_DIMz-1)/TILE_DIMz  
             myblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
             ii=i-xblock*TILE_DIMx+2*TILE_DIMx
             jj=j-yblock*TILE_DIMy+2*TILE_DIMy
             kk=k-zblock*TILE_DIMz+2*TILE_DIMz 
	         phi_loc=phifields_s(idx5(ii,jj,kk,1,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields))
	         u=u_s(i,j,k)
			 
			 global_phi_change = global_phi_change - u*phi_loc
	       enddo
	     enddo
         !$acc end parallel loop
	   endif
	 endif	
	    
	 if(pbc_y.eq.0)then
	   gj=2
       subchords(2)=(gj-1)/ny
	   if(subchords(2)==coords(2))then
		 !$acc parallel loop collapse(2) reduction(+:global_phi_change) present(phifields_s &
		 !$acc& ,u_s,v_s,w_s) private(i,j,k,l,gi,gj,gk &
		 !$acc& ,phitemp,conter,xblock,yblock,zblock,myblock,ii,jj,kk,phi_loc,v)    
	     do k=1,nz
		   do i=1,nx
			 gi=nx*coords(1)+i
			 gk=nz*coords(3)+k
	         j=2
	         gj=ny*coords(2)+j
	         if(isfluid(i,j,k) .ne. -1)cycle
	         xblock=(i+2*TILE_DIMx-1)/TILE_DIMx   
             yblock=(j+2*TILE_DIMy-1)/TILE_DIMy     
             zblock=(k+2*TILE_DIMz-1)/TILE_DIMz  
             myblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
             ii=i-xblock*TILE_DIMx+2*TILE_DIMx
             jj=j-yblock*TILE_DIMy+2*TILE_DIMy
             kk=k-zblock*TILE_DIMz+2*TILE_DIMz 
	         phi_loc=phifields_s(idx5(ii,jj,kk,1,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields))
	         v=v(i,j,k)
	         
			 global_phi_change = global_phi_change + v*phi_loc
	       enddo
	     enddo
         !$acc end parallel loop
	   endif
			 
	   gj=ly-1
       subchords(2)=(gj-1)/ny
	   if(subchords(2)==coords(2))then
		 !$acc parallel loop collapse(2) reduction(+:global_phi_change) present(phifields_s &
		 !$acc& ,u_s,v_s,w_s) private(i,j,k,l,gi,gj,gk &
		 !$acc& ,phitemp,conter,xblock,yblock,zblock,myblock,ii,jj,kk,phi_loc,v)	    
	     do k=1,nz
		   do i=1,nx
			 gi=nx*coords(1)+i
			 gk=nz*coords(3)+k			 
			 j=ny-1
	         gj=ny*coords(2)+j
	         if(isfluid(i,j,k) .ne. -1)cycle
	         xblock=(i+2*TILE_DIMx-1)/TILE_DIMx   
             yblock=(j+2*TILE_DIMy-1)/TILE_DIMy     
             zblock=(k+2*TILE_DIMz-1)/TILE_DIMz  
             myblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
             ii=i-xblock*TILE_DIMx+2*TILE_DIMx
             jj=j-yblock*TILE_DIMy+2*TILE_DIMy
             kk=k-zblock*TILE_DIMz+2*TILE_DIMz 
	         phi_loc=phifields_s(idx5(ii,jj,kk,1,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields))
	         v=v_s(i,j,k)
	         
			 global_phi_change = global_phi_change - v*phi_loc
	       enddo
	     enddo
      !$acc end parallel loop
	   endif
	 endif

	 if(pbc_z.eq.0)then
	   gk=2
       subchords(3)=(gk-1)/nz
	   if(subchords(3)==coords(3))then
		 !$acc parallel loop collapse(2) reduction(+:global_phi_change) present(phifields_s &
		 !$acc& ,u_s,v_s,w_s) private(i,j,k,l,gi,gj,gk &
		 !$acc& ,phitemp,conter,xblock,yblock,zblock,myblock,ii,jj,kk,phi_loc,w)  
	     do j=1,ny
		   do i=1,nx
			 gi=nx*coords(1)+i
			 gj=ny*coords(2)+j
	         k=2
	         gk=nz*coords(3)+k
	         if(isfluid(i,j,k) .ne. -1)cycle
	         xblock=(i+2*TILE_DIMx-1)/TILE_DIMx   
             yblock=(j+2*TILE_DIMy-1)/TILE_DIMy     
             zblock=(k+2*TILE_DIMz-1)/TILE_DIMz  
             myblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
             ii=i-xblock*TILE_DIMx+2*TILE_DIMx
             jj=j-yblock*TILE_DIMy+2*TILE_DIMy
             kk=k-zblock*TILE_DIMz+2*TILE_DIMz 
	         phi_loc=phifields_s(idx5(ii,jj,kk,1,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields))
	         w=w_s(i,j,k)
	         
			 global_phi_change = global_phi_change + w*phi_loc
	       enddo
	     enddo
      !$acc end parallel loop

	   endif
	   
	   gk=lz-1
       subchords(3)=(gk-1)/nz
	   if(subchords(3)==coords(3))then

		 !$acc parallel loop collapse(2) reduction(+:global_phi_change) present(phifields_s &
		 !$acc& ,u_s,v_s,w_s) private(i,j,k,l,gi,gj,gk &
		 !$acc& ,phitemp,conter,xblock,yblock,zblock,myblock,ii,jj,kk,phi_loc,wtmp)
  
	     do j=1,ny
		   do i=1,nx
			 gi=nx*coords(1)+i
			 gj=ny*coords(2)+j
			 k=nz-1
	         gk=nz*coords(3)+k
	         if(isfluid(i,j,k) .ne. -1)cycle
	         xblock=(i+2*TILE_DIMx-1)/TILE_DIMx   
             yblock=(j+2*TILE_DIMy-1)/TILE_DIMy     
             zblock=(k+2*TILE_DIMz-1)/TILE_DIMz  
             myblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
             ii=i-xblock*TILE_DIMx+2*TILE_DIMx
             jj=j-yblock*TILE_DIMy+2*TILE_DIMy
             kk=k-zblock*TILE_DIMz+2*TILE_DIMz 
	         phi_loc=phifields_s(idx5(ii,jj,kk,1,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields))
	         wtmp=w_s(i,j,k)
	         
			 global_phi_change = global_phi_change - wtmp*phi_loc
	       enddo
	     enddo
      !$acc end parallel loop

	    endif	
	 endif
	 
     !$acc wait
	 !$acc update host(global_phi_change)
	 !$acc wait
	 call sum_world_float(global_phi_change)

        
        call phi_sum_count_cuda( &
         press_s,u_s,v_s,w_s,pxx_s,pyy_s,pzz_s,pxy_s,pxz_s,pyz_s)
        
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
		
		call apply_lagrangian_phi_cuda(phifields_s, &
         press_s,u_s,v_s,w_s,pxx_s,pyy_s,pzz_s,pxy_s,pxz_s,pyz_s)

	 endif
#endif

#endif
#endif


   endsubroutine bcs_mesoscopic_phifields
   
endmodule
