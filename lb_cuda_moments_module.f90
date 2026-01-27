#include "defines.h"
#if !defined(_OPENACC)  
#error "To use this module the macros _OPENACC should be defined."
#endif

module lb_cuda_moments

   use vars
   use iso_c_binding
   use cudafor
   use mpi_template, only: coords,dostop,doerror,mydev,myrank,nprocs,nbuff,nbuffbvec
   use lb_cuda_vars

   implicit none
   

contains
   
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
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields,nblocks_d) :: phifields_s 
#endif           
#if defined(ELASTIC_FORCE)
      real(kind=db) :: lambda_rel,k_elastic
      real(kind=db), dimension(1:nx,1:ny,1:nz) :: u_ref,v_ref,w_ref
#endif
      real(kind=db) :: fx,fy,fz

      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields,nblocks_d) :: hfields_old,hfields_s
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nauxfields,nblocks_d) :: auxfields_s
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nlocauxfields,nblocks_d) :: locauxfields_s
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nforces,nblocks_d) :: forces_s
      

      real(kind=db) :: pxx,pyy,pzz,pxy,pxz,pyz
      real(kind=db) :: mytemp,forcex,forcey,forcez,rhophi_loc,press_loc,u_loc,v_loc,w_loc

#ifdef TWOCOMPONENT
	  real(kind=db) ::gradfix,gradfiy,gradfiz,wet_loc,phi_loc,lap_phi_loc
#endif
      
#ifdef DENSRATIO
      real(kind=db) :: gradrhox,gradrhoy,gradrhoz,omega_loc, visc_loc,tau_loc
#endif
      
      integer :: i,j,k
      integer :: gi,gj,gk
      integer :: myblock,ii,jj,kk
      
      i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z
      
      if(abs(isfluid(i,j,k)) .ne. 1)return
      
      gi=nx*coords(1)+i
      gj=ny*coords(2)+j
      gk=nz*coords(3)+k
      
      myblock=blockIdx%x+blockIdx%y*nxblock_d+blockIdx%z*nxyblock_d+1
      
      ii=threadIdx%x
      jj=threadIdx%y
      kk=threadIdx%z
      

				 
		  press_loc=hfields_s(ii,jj,kk,1,myblock)
				 
#ifdef TWOCOMPONENT					 
		  phi_loc=phifields_s(ii,jj,kk,1,myblock)
		  lap_phi_loc=locauxfields_s(ii,jj,kk,1,myblock)
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
		  mytemp=auxfields_s(ii,jj,kk,4,myblock) !modgrad
		  gradfix=auxfields_s(ii,jj,kk,1,myblock)*mytemp !normx*modgrad
		  gradfiy=auxfields_s(ii,jj,kk,2,myblock)*mytemp !normy*modgrad
		  gradfiz=auxfields_s(ii,jj,kk,3,myblock)*mytemp !normz*modgrad
		  forcex = forcex + &
                   (4.0_db*beta*phi_loc*(phi_loc-1.0_db)*(phi_loc-0.5_db) - kapp*lap_phi_loc)*gradfix
		  forcey = forcey + &
                   (4.0_db*beta*phi_loc*(phi_loc-1.0_db)*(phi_loc-0.5_db) - kapp*lap_phi_loc)*gradfiy
		  forcez = forcez + &
                   (4.0_db*beta*phi_loc*(phi_loc-1.0_db)*(phi_loc-0.5_db) - kapp*lap_phi_loc)*gradfiz
				   				   

#ifdef REPULSIVE_FLUX
		  mytemp=locauxfields_s(ii,jj,kk,6,myblock)*rhophi_loc 
		  if(abs(mytemp)>1.0d-3) mytemp=1.0d-3*sign(1.0,mytemp)!mytemp*0.1_db
		  forcex=forcex + mytemp*rhophi_loc
				  
		  mytemp=locauxfields_s(ii,jj,kk,7,myblock)*rhophi_loc 
		  if(abs(mytemp)>1.0d-3) mytemp=1.0d-3*sign(1.0,mytemp)!mytemp*0.1_db
		  forcey=forcey + mytemp*rhophi_loc
				  
		  mytemp=locauxfields_s(ii,jj,kk,8,myblock)*rhophi_loc 
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
                  
              
                  
		  forces_s(ii,jj,kk,1,myblock)=forcex
		  forces_s(ii,jj,kk,2,myblock)=forcey
		  forces_s(ii,jj,kk,3,myblock)=forcez
				  
		  u_loc=hfields_old(ii,jj,kk,2,myblock) !velocity
                  v_loc=hfields_old(ii,jj,kk,3,myblock)
                  w_loc=hfields_old(ii,jj,kk,4,myblock)
                  
                 
				  
#if defined(ELASTIC_FORCE)
                  forcex=forcex + &
                   rhophi_loc*(u_loc - u_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+   
		  forcey=forcey + &
		   rhophi_loc*(v_loc - v_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+ 
		  forcez=forcez + &
		   rhophi_loc*(w_loc - w_ref(i,j,k))*k_elastic*lambda_rel*phi_loc + rhophi_loc*fz !+  	
#endif

#ifdef DENSRATIO 
			  
                  pxx=hfields_s(ii,jj,kk,5,myblock)
                  pyy=hfields_s(ii,jj,kk,6,myblock)
                  pzz=hfields_s(ii,jj,kk,7,myblock)
                  pxy=hfields_s(ii,jj,kk,8,myblock)
                  pxz=hfields_s(ii,jj,kk,9,myblock)
                  pyz=hfields_s(ii,jj,kk,10,myblock)
                  
                  !1-2
                  !*1
                  ! 2nd order
                  pxx=pxx - cssq*press_loc - u_loc*u_loc
                  pyy=pyy - cssq*press_loc - v_loc*v_loc
                  pzz=pzz - cssq*press_loc - w_loc*w_loc
                  pxy=pxy - u_loc*v_loc
                  pxz=pxz - u_loc*w_loc
                  pyz=pyz - v_loc*w_loc

		  visc_loc=(rho_r*visc1*phi_loc+(1.0_db-phi_loc)*visc2*rho_b)/rhophi_loc
				  
                  tau_loc=(visc_loc*invcssq + HALF) !è una tau
				  
		  forcex=forcex - (visc_loc/(tau_loc*cssq))*(pxx*gradrhox + pxy*gradrhoy + pxz*gradrhoz)
		  forcey=forcey - (visc_loc/(tau_loc*cssq))*(pyy*gradrhoy + pxy*gradrhox + pyz*gradrhoz)
		  forcez=forcez - (visc_loc/(tau_loc*cssq))*(pzz*gradrhoz + pxz*gradrhox + pyz*gradrhoy)
#endif	
                  !I compute the new velocities
		  u_loc=hfields_s(ii,jj,kk,2,myblock) !velocity
                  v_loc=hfields_s(ii,jj,kk,3,myblock)
                  w_loc=hfields_s(ii,jj,kk,4,myblock)            
					 
#if defined(INTERFACE_INCOMP) && defined(DENSRATIO)
		  mytemp= -sharp_c*locauxfields_s(ii,jj,kk,2,myblock)
		  u_loc = u_loc + 0.5_db*forcex/(rhophi_loc)
                  v_loc = v_loc + 0.5_db*forcey/(rhophi_loc)
                  w_loc = w_loc + 0.5_db*forcez/(rhophi_loc)
				  
		  u_loc = u_loc/ &
		   (1.0_db - 0.5_db*(rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)/rhophi_loc )
				  
		  v_loc = v_loc/ &
		   (1.0_db - 0.5_db*(rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)/rhophi_loc ) 
				  
		  w_loc = w_loc/ &
	           (1.0_db - 0.5_db*(rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)/rhophi_loc ) 

		  forces_s(ii,jj,kk,1,myblock)= forces_s(ii,jj,kk,1,myblock) - &
		   (rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)*u_loc
		  forces_s(ii,jj,kk,2,myblock)= forces_s(ii,jj,kk,2,myblock) - &
		   (rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)*v_loc
		  forces_s(ii,jj,kk,3,myblock)= forces_s(ii,jj,kk,3,myblock) - &
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
                  

                  hfields_s(ii,jj,kk,2,myblock)=u_loc   !put the new velocity in hfields_s
                  hfields_s(ii,jj,kk,3,myblock)=v_loc
                  hfields_s(ii,jj,kk,4,myblock)=w_loc
                  
                                    

!regularized 
		  pxx=hfields_s(ii,jj,kk,5,myblock)
                  pyy=hfields_s(ii,jj,kk,6,myblock)
                  pzz=hfields_s(ii,jj,kk,7,myblock)
                  pxy=hfields_s(ii,jj,kk,8,myblock)
                  pxz=hfields_s(ii,jj,kk,9,myblock)
                  pyz=hfields_s(ii,jj,kk,10,myblock)
                  
                  
                  pxx=pxx + forcex*u_loc/rhophi_loc
                  pyy=pyy + forcey*v_loc/rhophi_loc
                  pzz=pzz + forcez*w_loc/rhophi_loc
                  pxy=pxy + HALF*(forcey*u_loc+forcex*v_loc)/rhophi_loc
                  pxz=pxz + HALF*(forcez*u_loc+forcex*w_loc)/rhophi_loc
                  pyz=pyz + HALF*(forcez*v_loc+forcey*w_loc)/rhophi_loc


                  hfields_s(ii,jj,kk,5,myblock)=pxx
                  hfields_s(ii,jj,kk,6,myblock)=pyy
                  hfields_s(ii,jj,kk,7,myblock)=pzz
                  hfields_s(ii,jj,kk,8,myblock)=pxy
                  hfields_s(ii,jj,kk,9,myblock)=pxz
                  hfields_s(ii,jj,kk,10,myblock)=pyz
                  
                  
#if defined(ELASTIC_FORCE)
		  u_ref(i,j,k) = u_ref(i,j,k) + &
		   lambda_rel*(u_loc - u_ref(i,j,k))
		  v_ref(i,j,k) = v_ref(i,j,k) + &
		   lambda_rel*(v_loc - v_ref(i,j,k))
		  w_ref(i,j,k) = w_ref(i,j,k) + &
		   lambda_rel*(w_loc - w_ref(i,j,k))
				  
                  forces_s(ii,jj,kk,1,myblock)= forces_s(ii,jj,kk,1,myblock) + &
                   rhophi_loc*(u_loc - u_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+ 
		  forces_s(ii,jj,kk,2,myblock)= forces_s(ii,jj,kk,2,myblock) +&
		   rhophi_loc*(v_loc - v_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+ 
		  forces_s(ii,jj,kk,3,myblock)= forces_s(ii,jj,kk,3,myblock) + &
		   rhophi_loc*(w_loc - w_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+rhophi_loc*fz  
#endif               
                  
#if defined(DENSRATIO)	
                  pxx=pxx - cssq*press_loc - u_loc*u_loc 
                  pyy=pyy - cssq*press_loc - v_loc*v_loc
                  pzz=pzz - cssq*press_loc - w_loc*w_loc 
                  pxy=pxy - u_loc*v_loc
                  pxz=pxz - u_loc*w_loc
                  pyz=pyz - v_loc*w_loc
		  
		  forces_s(ii,jj,kk,1,myblock)= forces_s(ii,jj,kk,1,myblock) - &
		   (visc_loc/(tau_loc*cssq))*(pxx*gradrhox + pxy*gradrhoy + pxz*gradrhoz)
		  forces_s(ii,jj,kk,2,myblock)= forces_s(ii,jj,kk,2,myblock) - &
		   (visc_loc/(tau_loc*cssq))*(pyy*gradrhoy + pxy*gradrhox + pyz*gradrhoz)
		  forces_s(ii,jj,kk,3,myblock)= forces_s(ii,jj,kk,3,myblock) - &
		   (visc_loc/(tau_loc*cssq))*(pzz*gradrhoz + pxz*gradrhox + pyz*gradrhoy)
#endif               


   endsubroutine moments_LB_kernel  
   
   attributes(global) subroutine moments_LB_kernel_int(step,iprobe,jprobe,kprobe,flop,nx,ny,nz,coords,isfluid &   
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
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields,nblocks_d) :: phifields_s 
#endif           
#if defined(ELASTIC_FORCE)
      real(kind=db) :: lambda_rel,k_elastic
      real(kind=db), dimension(1:nx,1:ny,1:nz) :: u_ref,v_ref,w_ref
#endif
      real(kind=db) :: fx,fy,fz

      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields,nblocks_d) :: hfields_old,hfields_s
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nauxfields,nblocks_d) :: auxfields_s
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nlocauxfields,nblocks_d) :: locauxfields_s
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nforces,nblocks_d) :: forces_s
      

      real(kind=db) :: pxx,pyy,pzz,pxy,pxz,pyz
      real(kind=db) :: mytemp,forcex,forcey,forcez,rhophi_loc,press_loc,u_loc,v_loc,w_loc

#ifdef TWOCOMPONENT
	  real(kind=db) ::gradfix,gradfiy,gradfiz,wet_loc,phi_loc,lap_phi_loc
#endif
      
#ifdef DENSRATIO
      real(kind=db) :: gradrhox,gradrhoy,gradrhoz,omega_loc, visc_loc,tau_loc
#endif
      
      integer :: i,j,k
      integer :: gi,gj,gk
      integer :: myblock,ii,jj,kk
      
      ii=threadIdx%x
      jj=threadIdx%y
      kk=threadIdx%z
      
      i = (blockIdx%x) * TILE_DIMx + ii
      j = (blockIdx%y) * TILE_DIMy + jj
      k = (blockIdx%z) * TILE_DIMz + kk
      
      if(abs(isfluid(i,j,k)) .ne. 1)return
      
      !gi=nx*coords(1)+i
      !gj=ny*coords(2)+j
      !gk=nz*coords(3)+k
      
      myblock=(blockIdx%x+1)+(blockIdx%y+1)*nxblock_d+(blockIdx%z+1)*nxyblock_d+1
      
				 
		  press_loc=hfields_s(ii,jj,kk,1,myblock)
				 
#ifdef TWOCOMPONENT					 
		  phi_loc=phifields_s(ii,jj,kk,1,myblock)
		  lap_phi_loc=locauxfields_s(ii,jj,kk,1,myblock)
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
		  mytemp=auxfields_s(ii,jj,kk,4,myblock) !modgrad
		  gradfix=auxfields_s(ii,jj,kk,1,myblock)*mytemp !normx*modgrad
		  gradfiy=auxfields_s(ii,jj,kk,2,myblock)*mytemp !normy*modgrad
		  gradfiz=auxfields_s(ii,jj,kk,3,myblock)*mytemp !normz*modgrad
		  forcex = forcex + &
                   (4.0_db*beta*phi_loc*(phi_loc-1.0_db)*(phi_loc-0.5_db) - kapp*lap_phi_loc)*gradfix
		  forcey = forcey + &
                   (4.0_db*beta*phi_loc*(phi_loc-1.0_db)*(phi_loc-0.5_db) - kapp*lap_phi_loc)*gradfiy
		  forcez = forcez + &
                   (4.0_db*beta*phi_loc*(phi_loc-1.0_db)*(phi_loc-0.5_db) - kapp*lap_phi_loc)*gradfiz
				   				   

#ifdef REPULSIVE_FLUX
		  mytemp=locauxfields_s(ii,jj,kk,6,myblock)*rhophi_loc 
		  if(abs(mytemp)>1.0d-3) mytemp=1.0d-3*sign(1.0,mytemp)!mytemp*0.1_db
		  forcex=forcex + mytemp*rhophi_loc
				  
		  mytemp=locauxfields_s(ii,jj,kk,7,myblock)*rhophi_loc 
		  if(abs(mytemp)>1.0d-3) mytemp=1.0d-3*sign(1.0,mytemp)!mytemp*0.1_db
		  forcey=forcey + mytemp*rhophi_loc
				  
		  mytemp=locauxfields_s(ii,jj,kk,8,myblock)*rhophi_loc 
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
                  
              
                  
		  forces_s(ii,jj,kk,1,myblock)=forcex
		  forces_s(ii,jj,kk,2,myblock)=forcey
		  forces_s(ii,jj,kk,3,myblock)=forcez
				  
		  u_loc=hfields_old(ii,jj,kk,2,myblock) !velocity
                  v_loc=hfields_old(ii,jj,kk,3,myblock)
                  w_loc=hfields_old(ii,jj,kk,4,myblock)
                  
                 
				  
#if defined(ELASTIC_FORCE)
                  forcex=forcex + &
                   rhophi_loc*(u_loc - u_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+   
		  forcey=forcey + &
		   rhophi_loc*(v_loc - v_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+ 
		  forcez=forcez + &
		   rhophi_loc*(w_loc - w_ref(i,j,k))*k_elastic*lambda_rel*phi_loc + rhophi_loc*fz !+  	
#endif

#ifdef DENSRATIO 
			  
                  pxx=hfields_s(ii,jj,kk,5,myblock)
                  pyy=hfields_s(ii,jj,kk,6,myblock)
                  pzz=hfields_s(ii,jj,kk,7,myblock)
                  pxy=hfields_s(ii,jj,kk,8,myblock)
                  pxz=hfields_s(ii,jj,kk,9,myblock)
                  pyz=hfields_s(ii,jj,kk,10,myblock)
                  
                  !1-2
                  !*1
                  ! 2nd order
                  pxx=pxx - cssq*press_loc - u_loc*u_loc
                  pyy=pyy - cssq*press_loc - v_loc*v_loc
                  pzz=pzz - cssq*press_loc - w_loc*w_loc
                  pxy=pxy - u_loc*v_loc
                  pxz=pxz - u_loc*w_loc
                  pyz=pyz - v_loc*w_loc

		  visc_loc=(rho_r*visc1*phi_loc+(1.0_db-phi_loc)*visc2*rho_b)/rhophi_loc
				  
                  tau_loc=(visc_loc*invcssq + HALF) !è una tau
				  
		  forcex=forcex - (visc_loc/(tau_loc*cssq))*(pxx*gradrhox + pxy*gradrhoy + pxz*gradrhoz)
		  forcey=forcey - (visc_loc/(tau_loc*cssq))*(pyy*gradrhoy + pxy*gradrhox + pyz*gradrhoz)
		  forcez=forcez - (visc_loc/(tau_loc*cssq))*(pzz*gradrhoz + pxz*gradrhox + pyz*gradrhoy)
#endif	
                  !I compute the new velocities
		  u_loc=hfields_s(ii,jj,kk,2,myblock) !velocity
                  v_loc=hfields_s(ii,jj,kk,3,myblock)
                  w_loc=hfields_s(ii,jj,kk,4,myblock)            
					 
#if defined(INTERFACE_INCOMP) && defined(DENSRATIO)
		  mytemp= -sharp_c*locauxfields_s(ii,jj,kk,2,myblock)
		  u_loc = u_loc + 0.5_db*forcex/(rhophi_loc)
                  v_loc = v_loc + 0.5_db*forcey/(rhophi_loc)
                  w_loc = w_loc + 0.5_db*forcez/(rhophi_loc)
				  
		  u_loc = u_loc/ &
		   (1.0_db - 0.5_db*(rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)/rhophi_loc )
				  
		  v_loc = v_loc/ &
		   (1.0_db - 0.5_db*(rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)/rhophi_loc ) 
				  
		  w_loc = w_loc/ &
	           (1.0_db - 0.5_db*(rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)/rhophi_loc ) 

		  forces_s(ii,jj,kk,1,myblock)= forces_s(ii,jj,kk,1,myblock) - &
		   (rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)*u_loc
		  forces_s(ii,jj,kk,2,myblock)= forces_s(ii,jj,kk,2,myblock) - &
		   (rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)*v_loc
		  forces_s(ii,jj,kk,3,myblock)= forces_s(ii,jj,kk,3,myblock) - &
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
                  

                  hfields_s(ii,jj,kk,2,myblock)=u_loc   !put the new velocity in hfields_s
                  hfields_s(ii,jj,kk,3,myblock)=v_loc
                  hfields_s(ii,jj,kk,4,myblock)=w_loc
                  
                                    

!regularized 
		  pxx=hfields_s(ii,jj,kk,5,myblock)
                  pyy=hfields_s(ii,jj,kk,6,myblock)
                  pzz=hfields_s(ii,jj,kk,7,myblock)
                  pxy=hfields_s(ii,jj,kk,8,myblock)
                  pxz=hfields_s(ii,jj,kk,9,myblock)
                  pyz=hfields_s(ii,jj,kk,10,myblock)
                  
                  
                  pxx=pxx + forcex*u_loc/rhophi_loc
                  pyy=pyy + forcey*v_loc/rhophi_loc
                  pzz=pzz + forcez*w_loc/rhophi_loc
                  pxy=pxy + HALF*(forcey*u_loc+forcex*v_loc)/rhophi_loc
                  pxz=pxz + HALF*(forcez*u_loc+forcex*w_loc)/rhophi_loc
                  pyz=pyz + HALF*(forcez*v_loc+forcey*w_loc)/rhophi_loc


                  hfields_s(ii,jj,kk,5,myblock)=pxx
                  hfields_s(ii,jj,kk,6,myblock)=pyy
                  hfields_s(ii,jj,kk,7,myblock)=pzz
                  hfields_s(ii,jj,kk,8,myblock)=pxy
                  hfields_s(ii,jj,kk,9,myblock)=pxz
                  hfields_s(ii,jj,kk,10,myblock)=pyz
                  
                  
#if defined(ELASTIC_FORCE)
		  u_ref(i,j,k) = u_ref(i,j,k) + &
		   lambda_rel*(u_loc - u_ref(i,j,k))
		  v_ref(i,j,k) = v_ref(i,j,k) + &
		   lambda_rel*(v_loc - v_ref(i,j,k))
		  w_ref(i,j,k) = w_ref(i,j,k) + &
		   lambda_rel*(w_loc - w_ref(i,j,k))
				  
                  forces_s(ii,jj,kk,1,myblock)= forces_s(ii,jj,kk,1,myblock) + &
                   rhophi_loc*(u_loc - u_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+ 
		  forces_s(ii,jj,kk,2,myblock)= forces_s(ii,jj,kk,2,myblock) +&
		   rhophi_loc*(v_loc - v_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+ 
		  forces_s(ii,jj,kk,3,myblock)= forces_s(ii,jj,kk,3,myblock) + &
		   rhophi_loc*(w_loc - w_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+rhophi_loc*fz  
#endif               
                  
#if defined(DENSRATIO)	
                  pxx=pxx - cssq*press_loc - u_loc*u_loc 
                  pyy=pyy - cssq*press_loc - v_loc*v_loc
                  pzz=pzz - cssq*press_loc - w_loc*w_loc 
                  pxy=pxy - u_loc*v_loc
                  pxz=pxz - u_loc*w_loc
                  pyz=pyz - v_loc*w_loc
		  
		  forces_s(ii,jj,kk,1,myblock)= forces_s(ii,jj,kk,1,myblock) - &
		   (visc_loc/(tau_loc*cssq))*(pxx*gradrhox + pxy*gradrhoy + pxz*gradrhoz)
		  forces_s(ii,jj,kk,2,myblock)= forces_s(ii,jj,kk,2,myblock) - &
		   (visc_loc/(tau_loc*cssq))*(pyy*gradrhoy + pxy*gradrhox + pyz*gradrhoz)
		  forces_s(ii,jj,kk,3,myblock)= forces_s(ii,jj,kk,3,myblock) - &
		   (visc_loc/(tau_loc*cssq))*(pzz*gradrhoz + pxz*gradrhox + pyz*gradrhoy)
#endif               


   endsubroutine moments_LB_kernel_int
   
      attributes(global) subroutine moments_LB_kernel_xplus(step,iprobe,jprobe,kprobe,flop,nx,ny,nz,coords,isfluid &   
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
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields,nblocks_d) :: phifields_s 
#endif           
#if defined(ELASTIC_FORCE)
      real(kind=db) :: lambda_rel,k_elastic
      real(kind=db), dimension(1:nx,1:ny,1:nz) :: u_ref,v_ref,w_ref
#endif
      real(kind=db) :: fx,fy,fz

      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields,nblocks_d) :: hfields_old,hfields_s
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nauxfields,nblocks_d) :: auxfields_s
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nlocauxfields,nblocks_d) :: locauxfields_s
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nforces,nblocks_d) :: forces_s
      

      real(kind=db) :: pxx,pyy,pzz,pxy,pxz,pyz
      real(kind=db) :: mytemp,forcex,forcey,forcez,rhophi_loc,press_loc,u_loc,v_loc,w_loc

#ifdef TWOCOMPONENT
	  real(kind=db) ::gradfix,gradfiy,gradfiz,wet_loc,phi_loc,lap_phi_loc
#endif
      
#ifdef DENSRATIO
      real(kind=db) :: gradrhox,gradrhoy,gradrhoz,omega_loc, visc_loc,tau_loc
#endif
      
      integer :: i,j,k
      integer :: gi,gj,gk
      integer :: myblock,ii,jj,kk
      
      ii=threadIdx%x
      jj=threadIdx%y
      kk=threadIdx%z
      
      i = ((nxblock_d-2)-1) * TILE_DIMx + ii
      j = (blockIdx%y) * TILE_DIMy + jj 
      k = (blockIdx%z) * TILE_DIMz + kk 
      
      if(abs(isfluid(i,j,k)) .ne. 1)return
      
      !gi=nx*coords(1)+i
      !gj=ny*coords(2)+j
      !gk=nz*coords(3)+k
      
      myblock=(nxblock_d-2)+(blockIdx%y+1)*nxblock_d+(blockIdx%z+1)*nxyblock_d+1
      
				 
		  press_loc=hfields_s(ii,jj,kk,1,myblock)
				 
#ifdef TWOCOMPONENT					 
		  phi_loc=phifields_s(ii,jj,kk,1,myblock)
		  lap_phi_loc=locauxfields_s(ii,jj,kk,1,myblock)
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
		  mytemp=auxfields_s(ii,jj,kk,4,myblock) !modgrad
		  gradfix=auxfields_s(ii,jj,kk,1,myblock)*mytemp !normx*modgrad
		  gradfiy=auxfields_s(ii,jj,kk,2,myblock)*mytemp !normy*modgrad
		  gradfiz=auxfields_s(ii,jj,kk,3,myblock)*mytemp !normz*modgrad
		  forcex = forcex + &
                   (4.0_db*beta*phi_loc*(phi_loc-1.0_db)*(phi_loc-0.5_db) - kapp*lap_phi_loc)*gradfix
		  forcey = forcey + &
                   (4.0_db*beta*phi_loc*(phi_loc-1.0_db)*(phi_loc-0.5_db) - kapp*lap_phi_loc)*gradfiy
		  forcez = forcez + &
                   (4.0_db*beta*phi_loc*(phi_loc-1.0_db)*(phi_loc-0.5_db) - kapp*lap_phi_loc)*gradfiz
				   				   

#ifdef REPULSIVE_FLUX
		  mytemp=locauxfields_s(ii,jj,kk,6,myblock)*rhophi_loc 
		  if(abs(mytemp)>1.0d-3) mytemp=1.0d-3*sign(1.0,mytemp)!mytemp*0.1_db
		  forcex=forcex + mytemp*rhophi_loc
				  
		  mytemp=locauxfields_s(ii,jj,kk,7,myblock)*rhophi_loc 
		  if(abs(mytemp)>1.0d-3) mytemp=1.0d-3*sign(1.0,mytemp)!mytemp*0.1_db
		  forcey=forcey + mytemp*rhophi_loc
				  
		  mytemp=locauxfields_s(ii,jj,kk,8,myblock)*rhophi_loc 
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
                  
              
                  
		  forces_s(ii,jj,kk,1,myblock)=forcex
		  forces_s(ii,jj,kk,2,myblock)=forcey
		  forces_s(ii,jj,kk,3,myblock)=forcez
				  
		  u_loc=hfields_old(ii,jj,kk,2,myblock) !velocity
                  v_loc=hfields_old(ii,jj,kk,3,myblock)
                  w_loc=hfields_old(ii,jj,kk,4,myblock)
                  
                 
				  
#if defined(ELASTIC_FORCE)
                  forcex=forcex + &
                   rhophi_loc*(u_loc - u_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+   
		  forcey=forcey + &
		   rhophi_loc*(v_loc - v_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+ 
		  forcez=forcez + &
		   rhophi_loc*(w_loc - w_ref(i,j,k))*k_elastic*lambda_rel*phi_loc + rhophi_loc*fz !+  	
#endif

#ifdef DENSRATIO 
			  
                  pxx=hfields_s(ii,jj,kk,5,myblock)
                  pyy=hfields_s(ii,jj,kk,6,myblock)
                  pzz=hfields_s(ii,jj,kk,7,myblock)
                  pxy=hfields_s(ii,jj,kk,8,myblock)
                  pxz=hfields_s(ii,jj,kk,9,myblock)
                  pyz=hfields_s(ii,jj,kk,10,myblock)
                  
                  !1-2
                  !*1
                  ! 2nd order
                  pxx=pxx - cssq*press_loc - u_loc*u_loc
                  pyy=pyy - cssq*press_loc - v_loc*v_loc
                  pzz=pzz - cssq*press_loc - w_loc*w_loc
                  pxy=pxy - u_loc*v_loc
                  pxz=pxz - u_loc*w_loc
                  pyz=pyz - v_loc*w_loc

		  visc_loc=(rho_r*visc1*phi_loc+(1.0_db-phi_loc)*visc2*rho_b)/rhophi_loc
				  
                  tau_loc=(visc_loc*invcssq + HALF) !è una tau
				  
		  forcex=forcex - (visc_loc/(tau_loc*cssq))*(pxx*gradrhox + pxy*gradrhoy + pxz*gradrhoz)
		  forcey=forcey - (visc_loc/(tau_loc*cssq))*(pyy*gradrhoy + pxy*gradrhox + pyz*gradrhoz)
		  forcez=forcez - (visc_loc/(tau_loc*cssq))*(pzz*gradrhoz + pxz*gradrhox + pyz*gradrhoy)
#endif	
                  !I compute the new velocities
		  u_loc=hfields_s(ii,jj,kk,2,myblock) !velocity
                  v_loc=hfields_s(ii,jj,kk,3,myblock)
                  w_loc=hfields_s(ii,jj,kk,4,myblock)            
					 
#if defined(INTERFACE_INCOMP) && defined(DENSRATIO)
		  mytemp= -sharp_c*locauxfields_s(ii,jj,kk,2,myblock)
		  u_loc = u_loc + 0.5_db*forcex/(rhophi_loc)
                  v_loc = v_loc + 0.5_db*forcey/(rhophi_loc)
                  w_loc = w_loc + 0.5_db*forcez/(rhophi_loc)
				  
		  u_loc = u_loc/ &
		   (1.0_db - 0.5_db*(rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)/rhophi_loc )
				  
		  v_loc = v_loc/ &
		   (1.0_db - 0.5_db*(rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)/rhophi_loc ) 
				  
		  w_loc = w_loc/ &
	           (1.0_db - 0.5_db*(rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)/rhophi_loc ) 

		  forces_s(ii,jj,kk,1,myblock)= forces_s(ii,jj,kk,1,myblock) - &
		   (rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)*u_loc
		  forces_s(ii,jj,kk,2,myblock)= forces_s(ii,jj,kk,2,myblock) - &
		   (rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)*v_loc
		  forces_s(ii,jj,kk,3,myblock)= forces_s(ii,jj,kk,3,myblock) - &
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
                  

                  hfields_s(ii,jj,kk,2,myblock)=u_loc   !put the new velocity in hfields_s
                  hfields_s(ii,jj,kk,3,myblock)=v_loc
                  hfields_s(ii,jj,kk,4,myblock)=w_loc
                  
                                    

!regularized 
		  pxx=hfields_s(ii,jj,kk,5,myblock)
                  pyy=hfields_s(ii,jj,kk,6,myblock)
                  pzz=hfields_s(ii,jj,kk,7,myblock)
                  pxy=hfields_s(ii,jj,kk,8,myblock)
                  pxz=hfields_s(ii,jj,kk,9,myblock)
                  pyz=hfields_s(ii,jj,kk,10,myblock)
                  
                  
                  pxx=pxx + forcex*u_loc/rhophi_loc
                  pyy=pyy + forcey*v_loc/rhophi_loc
                  pzz=pzz + forcez*w_loc/rhophi_loc
                  pxy=pxy + HALF*(forcey*u_loc+forcex*v_loc)/rhophi_loc
                  pxz=pxz + HALF*(forcez*u_loc+forcex*w_loc)/rhophi_loc
                  pyz=pyz + HALF*(forcez*v_loc+forcey*w_loc)/rhophi_loc


                  hfields_s(ii,jj,kk,5,myblock)=pxx
                  hfields_s(ii,jj,kk,6,myblock)=pyy
                  hfields_s(ii,jj,kk,7,myblock)=pzz
                  hfields_s(ii,jj,kk,8,myblock)=pxy
                  hfields_s(ii,jj,kk,9,myblock)=pxz
                  hfields_s(ii,jj,kk,10,myblock)=pyz
                  
                  
#if defined(ELASTIC_FORCE)
		  u_ref(i,j,k) = u_ref(i,j,k) + &
		   lambda_rel*(u_loc - u_ref(i,j,k))
		  v_ref(i,j,k) = v_ref(i,j,k) + &
		   lambda_rel*(v_loc - v_ref(i,j,k))
		  w_ref(i,j,k) = w_ref(i,j,k) + &
		   lambda_rel*(w_loc - w_ref(i,j,k))
				  
                  forces_s(ii,jj,kk,1,myblock)= forces_s(ii,jj,kk,1,myblock) + &
                   rhophi_loc*(u_loc - u_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+ 
		  forces_s(ii,jj,kk,2,myblock)= forces_s(ii,jj,kk,2,myblock) +&
		   rhophi_loc*(v_loc - v_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+ 
		  forces_s(ii,jj,kk,3,myblock)= forces_s(ii,jj,kk,3,myblock) + &
		   rhophi_loc*(w_loc - w_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+rhophi_loc*fz  
#endif               
                  
#if defined(DENSRATIO)	
                  pxx=pxx - cssq*press_loc - u_loc*u_loc 
                  pyy=pyy - cssq*press_loc - v_loc*v_loc
                  pzz=pzz - cssq*press_loc - w_loc*w_loc 
                  pxy=pxy - u_loc*v_loc
                  pxz=pxz - u_loc*w_loc
                  pyz=pyz - v_loc*w_loc
		  
		  forces_s(ii,jj,kk,1,myblock)= forces_s(ii,jj,kk,1,myblock) - &
		   (visc_loc/(tau_loc*cssq))*(pxx*gradrhox + pxy*gradrhoy + pxz*gradrhoz)
		  forces_s(ii,jj,kk,2,myblock)= forces_s(ii,jj,kk,2,myblock) - &
		   (visc_loc/(tau_loc*cssq))*(pyy*gradrhoy + pxy*gradrhox + pyz*gradrhoz)
		  forces_s(ii,jj,kk,3,myblock)= forces_s(ii,jj,kk,3,myblock) - &
		   (visc_loc/(tau_loc*cssq))*(pzz*gradrhoz + pxz*gradrhox + pyz*gradrhoy)
#endif               


   endsubroutine moments_LB_kernel_xplus
   
   attributes(global) subroutine moments_LB_kernel_xminus(step,iprobe,jprobe,kprobe,flop,nx,ny,nz,coords,isfluid &   
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
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields,nblocks_d) :: phifields_s 
#endif           
#if defined(ELASTIC_FORCE)
      real(kind=db) :: lambda_rel,k_elastic
      real(kind=db), dimension(1:nx,1:ny,1:nz) :: u_ref,v_ref,w_ref
#endif
      real(kind=db) :: fx,fy,fz

      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields,nblocks_d) :: hfields_old,hfields_s
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nauxfields,nblocks_d) :: auxfields_s
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nlocauxfields,nblocks_d) :: locauxfields_s
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nforces,nblocks_d) :: forces_s
      

      real(kind=db) :: pxx,pyy,pzz,pxy,pxz,pyz
      real(kind=db) :: mytemp,forcex,forcey,forcez,rhophi_loc,press_loc,u_loc,v_loc,w_loc

#ifdef TWOCOMPONENT
	  real(kind=db) ::gradfix,gradfiy,gradfiz,wet_loc,phi_loc,lap_phi_loc
#endif
      
#ifdef DENSRATIO
      real(kind=db) :: gradrhox,gradrhoy,gradrhoz,omega_loc, visc_loc,tau_loc
#endif
      
      integer :: i,j,k
      integer :: gi,gj,gk
      integer :: myblock,ii,jj,kk
      
      ii=threadIdx%x
      jj=threadIdx%y
      kk=threadIdx%z
      
      i = ii
      j = (blockIdx%y) * TILE_DIMy + jj 
      k = (blockIdx%z) * TILE_DIMz + kk 
      
      if(abs(isfluid(i,j,k)) .ne. 1)return
      
      !gi=nx*coords(1)+i
      !gj=ny*coords(2)+j
      !gk=nz*coords(3)+k
      
      myblock=1+(blockIdx%y+1)*nxblock_d+(blockIdx%z+1)*nxyblock_d+1
      
				 
		  press_loc=hfields_s(ii,jj,kk,1,myblock)
				 
#ifdef TWOCOMPONENT					 
		  phi_loc=phifields_s(ii,jj,kk,1,myblock)
		  lap_phi_loc=locauxfields_s(ii,jj,kk,1,myblock)
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
		  mytemp=auxfields_s(ii,jj,kk,4,myblock) !modgrad
		  gradfix=auxfields_s(ii,jj,kk,1,myblock)*mytemp !normx*modgrad
		  gradfiy=auxfields_s(ii,jj,kk,2,myblock)*mytemp !normy*modgrad
		  gradfiz=auxfields_s(ii,jj,kk,3,myblock)*mytemp !normz*modgrad
		  forcex = forcex + &
                   (4.0_db*beta*phi_loc*(phi_loc-1.0_db)*(phi_loc-0.5_db) - kapp*lap_phi_loc)*gradfix
		  forcey = forcey + &
                   (4.0_db*beta*phi_loc*(phi_loc-1.0_db)*(phi_loc-0.5_db) - kapp*lap_phi_loc)*gradfiy
		  forcez = forcez + &
                   (4.0_db*beta*phi_loc*(phi_loc-1.0_db)*(phi_loc-0.5_db) - kapp*lap_phi_loc)*gradfiz
				   				   

#ifdef REPULSIVE_FLUX
		  mytemp=locauxfields_s(ii,jj,kk,6,myblock)*rhophi_loc 
		  if(abs(mytemp)>1.0d-3) mytemp=1.0d-3*sign(1.0,mytemp)!mytemp*0.1_db
		  forcex=forcex + mytemp*rhophi_loc
				  
		  mytemp=locauxfields_s(ii,jj,kk,7,myblock)*rhophi_loc 
		  if(abs(mytemp)>1.0d-3) mytemp=1.0d-3*sign(1.0,mytemp)!mytemp*0.1_db
		  forcey=forcey + mytemp*rhophi_loc
				  
		  mytemp=locauxfields_s(ii,jj,kk,8,myblock)*rhophi_loc 
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
                  
              
                  
		  forces_s(ii,jj,kk,1,myblock)=forcex
		  forces_s(ii,jj,kk,2,myblock)=forcey
		  forces_s(ii,jj,kk,3,myblock)=forcez
				  
		  u_loc=hfields_old(ii,jj,kk,2,myblock) !velocity
                  v_loc=hfields_old(ii,jj,kk,3,myblock)
                  w_loc=hfields_old(ii,jj,kk,4,myblock)
                  
                 
				  
#if defined(ELASTIC_FORCE)
                  forcex=forcex + &
                   rhophi_loc*(u_loc - u_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+   
		  forcey=forcey + &
		   rhophi_loc*(v_loc - v_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+ 
		  forcez=forcez + &
		   rhophi_loc*(w_loc - w_ref(i,j,k))*k_elastic*lambda_rel*phi_loc + rhophi_loc*fz !+  	
#endif

#ifdef DENSRATIO 
			  
                  pxx=hfields_s(ii,jj,kk,5,myblock)
                  pyy=hfields_s(ii,jj,kk,6,myblock)
                  pzz=hfields_s(ii,jj,kk,7,myblock)
                  pxy=hfields_s(ii,jj,kk,8,myblock)
                  pxz=hfields_s(ii,jj,kk,9,myblock)
                  pyz=hfields_s(ii,jj,kk,10,myblock)
                  
                  !1-2
                  !*1
                  ! 2nd order
                  pxx=pxx - cssq*press_loc - u_loc*u_loc
                  pyy=pyy - cssq*press_loc - v_loc*v_loc
                  pzz=pzz - cssq*press_loc - w_loc*w_loc
                  pxy=pxy - u_loc*v_loc
                  pxz=pxz - u_loc*w_loc
                  pyz=pyz - v_loc*w_loc

		  visc_loc=(rho_r*visc1*phi_loc+(1.0_db-phi_loc)*visc2*rho_b)/rhophi_loc
				  
                  tau_loc=(visc_loc*invcssq + HALF) !è una tau
				  
		  forcex=forcex - (visc_loc/(tau_loc*cssq))*(pxx*gradrhox + pxy*gradrhoy + pxz*gradrhoz)
		  forcey=forcey - (visc_loc/(tau_loc*cssq))*(pyy*gradrhoy + pxy*gradrhox + pyz*gradrhoz)
		  forcez=forcez - (visc_loc/(tau_loc*cssq))*(pzz*gradrhoz + pxz*gradrhox + pyz*gradrhoy)
#endif	
                  !I compute the new velocities
		  u_loc=hfields_s(ii,jj,kk,2,myblock) !velocity
                  v_loc=hfields_s(ii,jj,kk,3,myblock)
                  w_loc=hfields_s(ii,jj,kk,4,myblock)            
					 
#if defined(INTERFACE_INCOMP) && defined(DENSRATIO)
		  mytemp= -sharp_c*locauxfields_s(ii,jj,kk,2,myblock)
		  u_loc = u_loc + 0.5_db*forcex/(rhophi_loc)
                  v_loc = v_loc + 0.5_db*forcey/(rhophi_loc)
                  w_loc = w_loc + 0.5_db*forcez/(rhophi_loc)
				  
		  u_loc = u_loc/ &
		   (1.0_db - 0.5_db*(rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)/rhophi_loc )
				  
		  v_loc = v_loc/ &
		   (1.0_db - 0.5_db*(rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)/rhophi_loc ) 
				  
		  w_loc = w_loc/ &
	           (1.0_db - 0.5_db*(rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)/rhophi_loc ) 

		  forces_s(ii,jj,kk,1,myblock)= forces_s(ii,jj,kk,1,myblock) - &
		   (rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)*u_loc
		  forces_s(ii,jj,kk,2,myblock)= forces_s(ii,jj,kk,2,myblock) - &
		   (rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)*v_loc
		  forces_s(ii,jj,kk,3,myblock)= forces_s(ii,jj,kk,3,myblock) - &
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
                  

                  hfields_s(ii,jj,kk,2,myblock)=u_loc   !put the new velocity in hfields_s
                  hfields_s(ii,jj,kk,3,myblock)=v_loc
                  hfields_s(ii,jj,kk,4,myblock)=w_loc
                  
                                    

!regularized 
		  pxx=hfields_s(ii,jj,kk,5,myblock)
                  pyy=hfields_s(ii,jj,kk,6,myblock)
                  pzz=hfields_s(ii,jj,kk,7,myblock)
                  pxy=hfields_s(ii,jj,kk,8,myblock)
                  pxz=hfields_s(ii,jj,kk,9,myblock)
                  pyz=hfields_s(ii,jj,kk,10,myblock)
                  
                  
                  pxx=pxx + forcex*u_loc/rhophi_loc
                  pyy=pyy + forcey*v_loc/rhophi_loc
                  pzz=pzz + forcez*w_loc/rhophi_loc
                  pxy=pxy + HALF*(forcey*u_loc+forcex*v_loc)/rhophi_loc
                  pxz=pxz + HALF*(forcez*u_loc+forcex*w_loc)/rhophi_loc
                  pyz=pyz + HALF*(forcez*v_loc+forcey*w_loc)/rhophi_loc


                  hfields_s(ii,jj,kk,5,myblock)=pxx
                  hfields_s(ii,jj,kk,6,myblock)=pyy
                  hfields_s(ii,jj,kk,7,myblock)=pzz
                  hfields_s(ii,jj,kk,8,myblock)=pxy
                  hfields_s(ii,jj,kk,9,myblock)=pxz
                  hfields_s(ii,jj,kk,10,myblock)=pyz
                  
                  
#if defined(ELASTIC_FORCE)
		  u_ref(i,j,k) = u_ref(i,j,k) + &
		   lambda_rel*(u_loc - u_ref(i,j,k))
		  v_ref(i,j,k) = v_ref(i,j,k) + &
		   lambda_rel*(v_loc - v_ref(i,j,k))
		  w_ref(i,j,k) = w_ref(i,j,k) + &
		   lambda_rel*(w_loc - w_ref(i,j,k))
				  
                  forces_s(ii,jj,kk,1,myblock)= forces_s(ii,jj,kk,1,myblock) + &
                   rhophi_loc*(u_loc - u_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+ 
		  forces_s(ii,jj,kk,2,myblock)= forces_s(ii,jj,kk,2,myblock) +&
		   rhophi_loc*(v_loc - v_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+ 
		  forces_s(ii,jj,kk,3,myblock)= forces_s(ii,jj,kk,3,myblock) + &
		   rhophi_loc*(w_loc - w_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+rhophi_loc*fz  
#endif               
                  
#if defined(DENSRATIO)	
                  pxx=pxx - cssq*press_loc - u_loc*u_loc 
                  pyy=pyy - cssq*press_loc - v_loc*v_loc
                  pzz=pzz - cssq*press_loc - w_loc*w_loc 
                  pxy=pxy - u_loc*v_loc
                  pxz=pxz - u_loc*w_loc
                  pyz=pyz - v_loc*w_loc
		  
		  forces_s(ii,jj,kk,1,myblock)= forces_s(ii,jj,kk,1,myblock) - &
		   (visc_loc/(tau_loc*cssq))*(pxx*gradrhox + pxy*gradrhoy + pxz*gradrhoz)
		  forces_s(ii,jj,kk,2,myblock)= forces_s(ii,jj,kk,2,myblock) - &
		   (visc_loc/(tau_loc*cssq))*(pyy*gradrhoy + pxy*gradrhox + pyz*gradrhoz)
		  forces_s(ii,jj,kk,3,myblock)= forces_s(ii,jj,kk,3,myblock) - &
		   (visc_loc/(tau_loc*cssq))*(pzz*gradrhoz + pxz*gradrhox + pyz*gradrhoy)
#endif               


   endsubroutine moments_LB_kernel_xminus
   
         attributes(global) subroutine moments_LB_kernel_yplus(step,iprobe,jprobe,kprobe,flop,nx,ny,nz,coords,isfluid &   
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
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields,nblocks_d) :: phifields_s 
#endif           
#if defined(ELASTIC_FORCE)
      real(kind=db) :: lambda_rel,k_elastic
      real(kind=db), dimension(1:nx,1:ny,1:nz) :: u_ref,v_ref,w_ref
#endif
      real(kind=db) :: fx,fy,fz

      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields,nblocks_d) :: hfields_old,hfields_s
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nauxfields,nblocks_d) :: auxfields_s
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nlocauxfields,nblocks_d) :: locauxfields_s
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nforces,nblocks_d) :: forces_s
      

      real(kind=db) :: pxx,pyy,pzz,pxy,pxz,pyz
      real(kind=db) :: mytemp,forcex,forcey,forcez,rhophi_loc,press_loc,u_loc,v_loc,w_loc

#ifdef TWOCOMPONENT
	  real(kind=db) ::gradfix,gradfiy,gradfiz,wet_loc,phi_loc,lap_phi_loc
#endif
      
#ifdef DENSRATIO
      real(kind=db) :: gradrhox,gradrhoy,gradrhoz,omega_loc, visc_loc,tau_loc
#endif
      
      integer :: i,j,k
      integer :: gi,gj,gk
      integer :: myblock,ii,jj,kk
      
      ii=threadIdx%x
      jj=threadIdx%y
      kk=threadIdx%z
      
      i = (blockIdx%x-1) * TILE_DIMx + ii
      j = ((nyblock_d-2)-1) * TILE_DIMy + jj
      k = (blockIdx%z) * TILE_DIMz + kk
      
      if(abs(isfluid(i,j,k)) .ne. 1)return
      
      !gi=nx*coords(1)+i
      !gj=ny*coords(2)+j
      !gk=nz*coords(3)+k
      
      myblock=blockIdx%x+(nyblock_d-2)*nxblock_d+(blockIdx%z+1)*nxyblock_d+1
      
				 
		  press_loc=hfields_s(ii,jj,kk,1,myblock)
				 
#ifdef TWOCOMPONENT					 
		  phi_loc=phifields_s(ii,jj,kk,1,myblock)
		  lap_phi_loc=locauxfields_s(ii,jj,kk,1,myblock)
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
		  mytemp=auxfields_s(ii,jj,kk,4,myblock) !modgrad
		  gradfix=auxfields_s(ii,jj,kk,1,myblock)*mytemp !normx*modgrad
		  gradfiy=auxfields_s(ii,jj,kk,2,myblock)*mytemp !normy*modgrad
		  gradfiz=auxfields_s(ii,jj,kk,3,myblock)*mytemp !normz*modgrad
		  forcex = forcex + &
                   (4.0_db*beta*phi_loc*(phi_loc-1.0_db)*(phi_loc-0.5_db) - kapp*lap_phi_loc)*gradfix
		  forcey = forcey + &
                   (4.0_db*beta*phi_loc*(phi_loc-1.0_db)*(phi_loc-0.5_db) - kapp*lap_phi_loc)*gradfiy
		  forcez = forcez + &
                   (4.0_db*beta*phi_loc*(phi_loc-1.0_db)*(phi_loc-0.5_db) - kapp*lap_phi_loc)*gradfiz
				   				   

#ifdef REPULSIVE_FLUX
		  mytemp=locauxfields_s(ii,jj,kk,6,myblock)*rhophi_loc 
		  if(abs(mytemp)>1.0d-3) mytemp=1.0d-3*sign(1.0,mytemp)!mytemp*0.1_db
		  forcex=forcex + mytemp*rhophi_loc
				  
		  mytemp=locauxfields_s(ii,jj,kk,7,myblock)*rhophi_loc 
		  if(abs(mytemp)>1.0d-3) mytemp=1.0d-3*sign(1.0,mytemp)!mytemp*0.1_db
		  forcey=forcey + mytemp*rhophi_loc
				  
		  mytemp=locauxfields_s(ii,jj,kk,8,myblock)*rhophi_loc 
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
                  
              
                  
		  forces_s(ii,jj,kk,1,myblock)=forcex
		  forces_s(ii,jj,kk,2,myblock)=forcey
		  forces_s(ii,jj,kk,3,myblock)=forcez
				  
		  u_loc=hfields_old(ii,jj,kk,2,myblock) !velocity
                  v_loc=hfields_old(ii,jj,kk,3,myblock)
                  w_loc=hfields_old(ii,jj,kk,4,myblock)
                  
                 
				  
#if defined(ELASTIC_FORCE)
                  forcex=forcex + &
                   rhophi_loc*(u_loc - u_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+   
		  forcey=forcey + &
		   rhophi_loc*(v_loc - v_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+ 
		  forcez=forcez + &
		   rhophi_loc*(w_loc - w_ref(i,j,k))*k_elastic*lambda_rel*phi_loc + rhophi_loc*fz !+  	
#endif

#ifdef DENSRATIO 
			  
                  pxx=hfields_s(ii,jj,kk,5,myblock)
                  pyy=hfields_s(ii,jj,kk,6,myblock)
                  pzz=hfields_s(ii,jj,kk,7,myblock)
                  pxy=hfields_s(ii,jj,kk,8,myblock)
                  pxz=hfields_s(ii,jj,kk,9,myblock)
                  pyz=hfields_s(ii,jj,kk,10,myblock)
                  
                  !1-2
                  !*1
                  ! 2nd order
                  pxx=pxx - cssq*press_loc - u_loc*u_loc
                  pyy=pyy - cssq*press_loc - v_loc*v_loc
                  pzz=pzz - cssq*press_loc - w_loc*w_loc
                  pxy=pxy - u_loc*v_loc
                  pxz=pxz - u_loc*w_loc
                  pyz=pyz - v_loc*w_loc

		  visc_loc=(rho_r*visc1*phi_loc+(1.0_db-phi_loc)*visc2*rho_b)/rhophi_loc
				  
                  tau_loc=(visc_loc*invcssq + HALF) !è una tau
				  
		  forcex=forcex - (visc_loc/(tau_loc*cssq))*(pxx*gradrhox + pxy*gradrhoy + pxz*gradrhoz)
		  forcey=forcey - (visc_loc/(tau_loc*cssq))*(pyy*gradrhoy + pxy*gradrhox + pyz*gradrhoz)
		  forcez=forcez - (visc_loc/(tau_loc*cssq))*(pzz*gradrhoz + pxz*gradrhox + pyz*gradrhoy)
#endif	
                  !I compute the new velocities
		  u_loc=hfields_s(ii,jj,kk,2,myblock) !velocity
                  v_loc=hfields_s(ii,jj,kk,3,myblock)
                  w_loc=hfields_s(ii,jj,kk,4,myblock)            
					 
#if defined(INTERFACE_INCOMP) && defined(DENSRATIO)
		  mytemp= -sharp_c*locauxfields_s(ii,jj,kk,2,myblock)
		  u_loc = u_loc + 0.5_db*forcex/(rhophi_loc)
                  v_loc = v_loc + 0.5_db*forcey/(rhophi_loc)
                  w_loc = w_loc + 0.5_db*forcez/(rhophi_loc)
				  
		  u_loc = u_loc/ &
		   (1.0_db - 0.5_db*(rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)/rhophi_loc )
				  
		  v_loc = v_loc/ &
		   (1.0_db - 0.5_db*(rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)/rhophi_loc ) 
				  
		  w_loc = w_loc/ &
	           (1.0_db - 0.5_db*(rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)/rhophi_loc ) 

		  forces_s(ii,jj,kk,1,myblock)= forces_s(ii,jj,kk,1,myblock) - &
		   (rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)*u_loc
		  forces_s(ii,jj,kk,2,myblock)= forces_s(ii,jj,kk,2,myblock) - &
		   (rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)*v_loc
		  forces_s(ii,jj,kk,3,myblock)= forces_s(ii,jj,kk,3,myblock) - &
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
                  

                  hfields_s(ii,jj,kk,2,myblock)=u_loc   !put the new velocity in hfields_s
                  hfields_s(ii,jj,kk,3,myblock)=v_loc
                  hfields_s(ii,jj,kk,4,myblock)=w_loc
                  
                                    

!regularized 
		  pxx=hfields_s(ii,jj,kk,5,myblock)
                  pyy=hfields_s(ii,jj,kk,6,myblock)
                  pzz=hfields_s(ii,jj,kk,7,myblock)
                  pxy=hfields_s(ii,jj,kk,8,myblock)
                  pxz=hfields_s(ii,jj,kk,9,myblock)
                  pyz=hfields_s(ii,jj,kk,10,myblock)
                  
                  
                  pxx=pxx + forcex*u_loc/rhophi_loc
                  pyy=pyy + forcey*v_loc/rhophi_loc
                  pzz=pzz + forcez*w_loc/rhophi_loc
                  pxy=pxy + HALF*(forcey*u_loc+forcex*v_loc)/rhophi_loc
                  pxz=pxz + HALF*(forcez*u_loc+forcex*w_loc)/rhophi_loc
                  pyz=pyz + HALF*(forcez*v_loc+forcey*w_loc)/rhophi_loc


                  hfields_s(ii,jj,kk,5,myblock)=pxx
                  hfields_s(ii,jj,kk,6,myblock)=pyy
                  hfields_s(ii,jj,kk,7,myblock)=pzz
                  hfields_s(ii,jj,kk,8,myblock)=pxy
                  hfields_s(ii,jj,kk,9,myblock)=pxz
                  hfields_s(ii,jj,kk,10,myblock)=pyz
                  
                  
#if defined(ELASTIC_FORCE)
		  u_ref(i,j,k) = u_ref(i,j,k) + &
		   lambda_rel*(u_loc - u_ref(i,j,k))
		  v_ref(i,j,k) = v_ref(i,j,k) + &
		   lambda_rel*(v_loc - v_ref(i,j,k))
		  w_ref(i,j,k) = w_ref(i,j,k) + &
		   lambda_rel*(w_loc - w_ref(i,j,k))
				  
                  forces_s(ii,jj,kk,1,myblock)= forces_s(ii,jj,kk,1,myblock) + &
                   rhophi_loc*(u_loc - u_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+ 
		  forces_s(ii,jj,kk,2,myblock)= forces_s(ii,jj,kk,2,myblock) +&
		   rhophi_loc*(v_loc - v_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+ 
		  forces_s(ii,jj,kk,3,myblock)= forces_s(ii,jj,kk,3,myblock) + &
		   rhophi_loc*(w_loc - w_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+rhophi_loc*fz  
#endif               
                  
#if defined(DENSRATIO)	
                  pxx=pxx - cssq*press_loc - u_loc*u_loc 
                  pyy=pyy - cssq*press_loc - v_loc*v_loc
                  pzz=pzz - cssq*press_loc - w_loc*w_loc 
                  pxy=pxy - u_loc*v_loc
                  pxz=pxz - u_loc*w_loc
                  pyz=pyz - v_loc*w_loc
		  
		  forces_s(ii,jj,kk,1,myblock)= forces_s(ii,jj,kk,1,myblock) - &
		   (visc_loc/(tau_loc*cssq))*(pxx*gradrhox + pxy*gradrhoy + pxz*gradrhoz)
		  forces_s(ii,jj,kk,2,myblock)= forces_s(ii,jj,kk,2,myblock) - &
		   (visc_loc/(tau_loc*cssq))*(pyy*gradrhoy + pxy*gradrhox + pyz*gradrhoz)
		  forces_s(ii,jj,kk,3,myblock)= forces_s(ii,jj,kk,3,myblock) - &
		   (visc_loc/(tau_loc*cssq))*(pzz*gradrhoz + pxz*gradrhox + pyz*gradrhoy)
#endif               


   endsubroutine moments_LB_kernel_yplus
   
   attributes(global) subroutine moments_LB_kernel_yminus(step,iprobe,jprobe,kprobe,flop,nx,ny,nz,coords,isfluid &   
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
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields,nblocks_d) :: phifields_s 
#endif           
#if defined(ELASTIC_FORCE)
      real(kind=db) :: lambda_rel,k_elastic
      real(kind=db), dimension(1:nx,1:ny,1:nz) :: u_ref,v_ref,w_ref
#endif
      real(kind=db) :: fx,fy,fz

      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields,nblocks_d) :: hfields_old,hfields_s
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nauxfields,nblocks_d) :: auxfields_s
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nlocauxfields,nblocks_d) :: locauxfields_s
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nforces,nblocks_d) :: forces_s
      

      real(kind=db) :: pxx,pyy,pzz,pxy,pxz,pyz
      real(kind=db) :: mytemp,forcex,forcey,forcez,rhophi_loc,press_loc,u_loc,v_loc,w_loc

#ifdef TWOCOMPONENT
	  real(kind=db) ::gradfix,gradfiy,gradfiz,wet_loc,phi_loc,lap_phi_loc
#endif
      
#ifdef DENSRATIO
      real(kind=db) :: gradrhox,gradrhoy,gradrhoz,omega_loc, visc_loc,tau_loc
#endif
      
      integer :: i,j,k
      integer :: gi,gj,gk
      integer :: myblock,ii,jj,kk
      
      ii=threadIdx%x
      jj=threadIdx%y
      kk=threadIdx%z
      
      i = (blockIdx%x-1) * TILE_DIMx + ii 
      j = jj
      k = (blockIdx%z) * TILE_DIMz + kk
      
      if(abs(isfluid(i,j,k)) .ne. 1)return
      
      !gi=nx*coords(1)+i
      !gj=ny*coords(2)+j
      !gk=nz*coords(3)+k
      
      myblock=blockIdx%x+1*nxblock_d+(blockIdx%z+1)*nxyblock_d+1
      
				 
		  press_loc=hfields_s(ii,jj,kk,1,myblock)
				 
#ifdef TWOCOMPONENT					 
		  phi_loc=phifields_s(ii,jj,kk,1,myblock)
		  lap_phi_loc=locauxfields_s(ii,jj,kk,1,myblock)
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
		  mytemp=auxfields_s(ii,jj,kk,4,myblock) !modgrad
		  gradfix=auxfields_s(ii,jj,kk,1,myblock)*mytemp !normx*modgrad
		  gradfiy=auxfields_s(ii,jj,kk,2,myblock)*mytemp !normy*modgrad
		  gradfiz=auxfields_s(ii,jj,kk,3,myblock)*mytemp !normz*modgrad
		  forcex = forcex + &
                   (4.0_db*beta*phi_loc*(phi_loc-1.0_db)*(phi_loc-0.5_db) - kapp*lap_phi_loc)*gradfix
		  forcey = forcey + &
                   (4.0_db*beta*phi_loc*(phi_loc-1.0_db)*(phi_loc-0.5_db) - kapp*lap_phi_loc)*gradfiy
		  forcez = forcez + &
                   (4.0_db*beta*phi_loc*(phi_loc-1.0_db)*(phi_loc-0.5_db) - kapp*lap_phi_loc)*gradfiz
				   				   

#ifdef REPULSIVE_FLUX
		  mytemp=locauxfields_s(ii,jj,kk,6,myblock)*rhophi_loc 
		  if(abs(mytemp)>1.0d-3) mytemp=1.0d-3*sign(1.0,mytemp)!mytemp*0.1_db
		  forcex=forcex + mytemp*rhophi_loc
				  
		  mytemp=locauxfields_s(ii,jj,kk,7,myblock)*rhophi_loc 
		  if(abs(mytemp)>1.0d-3) mytemp=1.0d-3*sign(1.0,mytemp)!mytemp*0.1_db
		  forcey=forcey + mytemp*rhophi_loc
				  
		  mytemp=locauxfields_s(ii,jj,kk,8,myblock)*rhophi_loc 
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
                  
              
                  
		  forces_s(ii,jj,kk,1,myblock)=forcex
		  forces_s(ii,jj,kk,2,myblock)=forcey
		  forces_s(ii,jj,kk,3,myblock)=forcez
				  
		  u_loc=hfields_old(ii,jj,kk,2,myblock) !velocity
                  v_loc=hfields_old(ii,jj,kk,3,myblock)
                  w_loc=hfields_old(ii,jj,kk,4,myblock)
                  
                 
				  
#if defined(ELASTIC_FORCE)
                  forcex=forcex + &
                   rhophi_loc*(u_loc - u_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+   
		  forcey=forcey + &
		   rhophi_loc*(v_loc - v_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+ 
		  forcez=forcez + &
		   rhophi_loc*(w_loc - w_ref(i,j,k))*k_elastic*lambda_rel*phi_loc + rhophi_loc*fz !+  	
#endif

#ifdef DENSRATIO 
			  
                  pxx=hfields_s(ii,jj,kk,5,myblock)
                  pyy=hfields_s(ii,jj,kk,6,myblock)
                  pzz=hfields_s(ii,jj,kk,7,myblock)
                  pxy=hfields_s(ii,jj,kk,8,myblock)
                  pxz=hfields_s(ii,jj,kk,9,myblock)
                  pyz=hfields_s(ii,jj,kk,10,myblock)
                  
                  !1-2
                  !*1
                  ! 2nd order
                  pxx=pxx - cssq*press_loc - u_loc*u_loc
                  pyy=pyy - cssq*press_loc - v_loc*v_loc
                  pzz=pzz - cssq*press_loc - w_loc*w_loc
                  pxy=pxy - u_loc*v_loc
                  pxz=pxz - u_loc*w_loc
                  pyz=pyz - v_loc*w_loc

		  visc_loc=(rho_r*visc1*phi_loc+(1.0_db-phi_loc)*visc2*rho_b)/rhophi_loc
				  
                  tau_loc=(visc_loc*invcssq + HALF) !è una tau
				  
		  forcex=forcex - (visc_loc/(tau_loc*cssq))*(pxx*gradrhox + pxy*gradrhoy + pxz*gradrhoz)
		  forcey=forcey - (visc_loc/(tau_loc*cssq))*(pyy*gradrhoy + pxy*gradrhox + pyz*gradrhoz)
		  forcez=forcez - (visc_loc/(tau_loc*cssq))*(pzz*gradrhoz + pxz*gradrhox + pyz*gradrhoy)
#endif	
                  !I compute the new velocities
		  u_loc=hfields_s(ii,jj,kk,2,myblock) !velocity
                  v_loc=hfields_s(ii,jj,kk,3,myblock)
                  w_loc=hfields_s(ii,jj,kk,4,myblock)            
					 
#if defined(INTERFACE_INCOMP) && defined(DENSRATIO)
		  mytemp= -sharp_c*locauxfields_s(ii,jj,kk,2,myblock)
		  u_loc = u_loc + 0.5_db*forcex/(rhophi_loc)
                  v_loc = v_loc + 0.5_db*forcey/(rhophi_loc)
                  w_loc = w_loc + 0.5_db*forcez/(rhophi_loc)
				  
		  u_loc = u_loc/ &
		   (1.0_db - 0.5_db*(rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)/rhophi_loc )
				  
		  v_loc = v_loc/ &
		   (1.0_db - 0.5_db*(rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)/rhophi_loc ) 
				  
		  w_loc = w_loc/ &
	           (1.0_db - 0.5_db*(rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)/rhophi_loc ) 

		  forces_s(ii,jj,kk,1,myblock)= forces_s(ii,jj,kk,1,myblock) - &
		   (rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)*u_loc
		  forces_s(ii,jj,kk,2,myblock)= forces_s(ii,jj,kk,2,myblock) - &
		   (rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)*v_loc
		  forces_s(ii,jj,kk,3,myblock)= forces_s(ii,jj,kk,3,myblock) - &
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
                  

                  hfields_s(ii,jj,kk,2,myblock)=u_loc   !put the new velocity in hfields_s
                  hfields_s(ii,jj,kk,3,myblock)=v_loc
                  hfields_s(ii,jj,kk,4,myblock)=w_loc
                  
                                    

!regularized 
		  pxx=hfields_s(ii,jj,kk,5,myblock)
                  pyy=hfields_s(ii,jj,kk,6,myblock)
                  pzz=hfields_s(ii,jj,kk,7,myblock)
                  pxy=hfields_s(ii,jj,kk,8,myblock)
                  pxz=hfields_s(ii,jj,kk,9,myblock)
                  pyz=hfields_s(ii,jj,kk,10,myblock)
                  
                  
                  pxx=pxx + forcex*u_loc/rhophi_loc
                  pyy=pyy + forcey*v_loc/rhophi_loc
                  pzz=pzz + forcez*w_loc/rhophi_loc
                  pxy=pxy + HALF*(forcey*u_loc+forcex*v_loc)/rhophi_loc
                  pxz=pxz + HALF*(forcez*u_loc+forcex*w_loc)/rhophi_loc
                  pyz=pyz + HALF*(forcez*v_loc+forcey*w_loc)/rhophi_loc


                  hfields_s(ii,jj,kk,5,myblock)=pxx
                  hfields_s(ii,jj,kk,6,myblock)=pyy
                  hfields_s(ii,jj,kk,7,myblock)=pzz
                  hfields_s(ii,jj,kk,8,myblock)=pxy
                  hfields_s(ii,jj,kk,9,myblock)=pxz
                  hfields_s(ii,jj,kk,10,myblock)=pyz
                  
                  
#if defined(ELASTIC_FORCE)
		  u_ref(i,j,k) = u_ref(i,j,k) + &
		   lambda_rel*(u_loc - u_ref(i,j,k))
		  v_ref(i,j,k) = v_ref(i,j,k) + &
		   lambda_rel*(v_loc - v_ref(i,j,k))
		  w_ref(i,j,k) = w_ref(i,j,k) + &
		   lambda_rel*(w_loc - w_ref(i,j,k))
				  
                  forces_s(ii,jj,kk,1,myblock)= forces_s(ii,jj,kk,1,myblock) + &
                   rhophi_loc*(u_loc - u_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+ 
		  forces_s(ii,jj,kk,2,myblock)= forces_s(ii,jj,kk,2,myblock) +&
		   rhophi_loc*(v_loc - v_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+ 
		  forces_s(ii,jj,kk,3,myblock)= forces_s(ii,jj,kk,3,myblock) + &
		   rhophi_loc*(w_loc - w_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+rhophi_loc*fz  
#endif               
                  
#if defined(DENSRATIO)	
                  pxx=pxx - cssq*press_loc - u_loc*u_loc 
                  pyy=pyy - cssq*press_loc - v_loc*v_loc
                  pzz=pzz - cssq*press_loc - w_loc*w_loc 
                  pxy=pxy - u_loc*v_loc
                  pxz=pxz - u_loc*w_loc
                  pyz=pyz - v_loc*w_loc
		  
		  forces_s(ii,jj,kk,1,myblock)= forces_s(ii,jj,kk,1,myblock) - &
		   (visc_loc/(tau_loc*cssq))*(pxx*gradrhox + pxy*gradrhoy + pxz*gradrhoz)
		  forces_s(ii,jj,kk,2,myblock)= forces_s(ii,jj,kk,2,myblock) - &
		   (visc_loc/(tau_loc*cssq))*(pyy*gradrhoy + pxy*gradrhox + pyz*gradrhoz)
		  forces_s(ii,jj,kk,3,myblock)= forces_s(ii,jj,kk,3,myblock) - &
		   (visc_loc/(tau_loc*cssq))*(pzz*gradrhoz + pxz*gradrhox + pyz*gradrhoy)
#endif               


   endsubroutine moments_LB_kernel_yminus
   
         attributes(global) subroutine moments_LB_kernel_zplus(step,iprobe,jprobe,kprobe,flop,nx,ny,nz,coords,isfluid &   
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
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields,nblocks_d) :: phifields_s 
#endif           
#if defined(ELASTIC_FORCE)
      real(kind=db) :: lambda_rel,k_elastic
      real(kind=db), dimension(1:nx,1:ny,1:nz) :: u_ref,v_ref,w_ref
#endif
      real(kind=db) :: fx,fy,fz

      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields,nblocks_d) :: hfields_old,hfields_s
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nauxfields,nblocks_d) :: auxfields_s
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nlocauxfields,nblocks_d) :: locauxfields_s
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nforces,nblocks_d) :: forces_s
      

      real(kind=db) :: pxx,pyy,pzz,pxy,pxz,pyz
      real(kind=db) :: mytemp,forcex,forcey,forcez,rhophi_loc,press_loc,u_loc,v_loc,w_loc

#ifdef TWOCOMPONENT
	  real(kind=db) ::gradfix,gradfiy,gradfiz,wet_loc,phi_loc,lap_phi_loc
#endif
      
#ifdef DENSRATIO
      real(kind=db) :: gradrhox,gradrhoy,gradrhoz,omega_loc, visc_loc,tau_loc
#endif
      
      integer :: i,j,k
      integer :: gi,gj,gk
      integer :: myblock,ii,jj,kk
      
      ii=threadIdx%x
      jj=threadIdx%y
      kk=threadIdx%z
      
      i = (blockIdx%x-1) * TILE_DIMx + ii
      j = (blockIdx%y-1) * TILE_DIMy + jj
      k = ((nzblock_d-2)-1) * TILE_DIMz + kk
      
      if(abs(isfluid(i,j,k)) .ne. 1)return
      
      !gi=nx*coords(1)+i
      !gj=ny*coords(2)+j
      !gk=nz*coords(3)+k
      
      myblock=blockIdx%x+blockIdx%y*nxblock_d+(nzblock_d-2)*nxyblock_d+1
      
				 
		  press_loc=hfields_s(ii,jj,kk,1,myblock)
				 
#ifdef TWOCOMPONENT					 
		  phi_loc=phifields_s(ii,jj,kk,1,myblock)
		  lap_phi_loc=locauxfields_s(ii,jj,kk,1,myblock)
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
		  mytemp=auxfields_s(ii,jj,kk,4,myblock) !modgrad
		  gradfix=auxfields_s(ii,jj,kk,1,myblock)*mytemp !normx*modgrad
		  gradfiy=auxfields_s(ii,jj,kk,2,myblock)*mytemp !normy*modgrad
		  gradfiz=auxfields_s(ii,jj,kk,3,myblock)*mytemp !normz*modgrad
		  forcex = forcex + &
                   (4.0_db*beta*phi_loc*(phi_loc-1.0_db)*(phi_loc-0.5_db) - kapp*lap_phi_loc)*gradfix
		  forcey = forcey + &
                   (4.0_db*beta*phi_loc*(phi_loc-1.0_db)*(phi_loc-0.5_db) - kapp*lap_phi_loc)*gradfiy
		  forcez = forcez + &
                   (4.0_db*beta*phi_loc*(phi_loc-1.0_db)*(phi_loc-0.5_db) - kapp*lap_phi_loc)*gradfiz
				   				   

#ifdef REPULSIVE_FLUX
		  mytemp=locauxfields_s(ii,jj,kk,6,myblock)*rhophi_loc 
		  if(abs(mytemp)>1.0d-3) mytemp=1.0d-3*sign(1.0,mytemp)!mytemp*0.1_db
		  forcex=forcex + mytemp*rhophi_loc
				  
		  mytemp=locauxfields_s(ii,jj,kk,7,myblock)*rhophi_loc 
		  if(abs(mytemp)>1.0d-3) mytemp=1.0d-3*sign(1.0,mytemp)!mytemp*0.1_db
		  forcey=forcey + mytemp*rhophi_loc
				  
		  mytemp=locauxfields_s(ii,jj,kk,8,myblock)*rhophi_loc 
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
                  
              
                  
		  forces_s(ii,jj,kk,1,myblock)=forcex
		  forces_s(ii,jj,kk,2,myblock)=forcey
		  forces_s(ii,jj,kk,3,myblock)=forcez
				  
		  u_loc=hfields_old(ii,jj,kk,2,myblock) !velocity
                  v_loc=hfields_old(ii,jj,kk,3,myblock)
                  w_loc=hfields_old(ii,jj,kk,4,myblock)
                  
                 
				  
#if defined(ELASTIC_FORCE)
                  forcex=forcex + &
                   rhophi_loc*(u_loc - u_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+   
		  forcey=forcey + &
		   rhophi_loc*(v_loc - v_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+ 
		  forcez=forcez + &
		   rhophi_loc*(w_loc - w_ref(i,j,k))*k_elastic*lambda_rel*phi_loc + rhophi_loc*fz !+  	
#endif

#ifdef DENSRATIO 
			  
                  pxx=hfields_s(ii,jj,kk,5,myblock)
                  pyy=hfields_s(ii,jj,kk,6,myblock)
                  pzz=hfields_s(ii,jj,kk,7,myblock)
                  pxy=hfields_s(ii,jj,kk,8,myblock)
                  pxz=hfields_s(ii,jj,kk,9,myblock)
                  pyz=hfields_s(ii,jj,kk,10,myblock)
                  
                  !1-2
                  !*1
                  ! 2nd order
                  pxx=pxx - cssq*press_loc - u_loc*u_loc
                  pyy=pyy - cssq*press_loc - v_loc*v_loc
                  pzz=pzz - cssq*press_loc - w_loc*w_loc
                  pxy=pxy - u_loc*v_loc
                  pxz=pxz - u_loc*w_loc
                  pyz=pyz - v_loc*w_loc

		  visc_loc=(rho_r*visc1*phi_loc+(1.0_db-phi_loc)*visc2*rho_b)/rhophi_loc
				  
                  tau_loc=(visc_loc*invcssq + HALF) !è una tau
				  
		  forcex=forcex - (visc_loc/(tau_loc*cssq))*(pxx*gradrhox + pxy*gradrhoy + pxz*gradrhoz)
		  forcey=forcey - (visc_loc/(tau_loc*cssq))*(pyy*gradrhoy + pxy*gradrhox + pyz*gradrhoz)
		  forcez=forcez - (visc_loc/(tau_loc*cssq))*(pzz*gradrhoz + pxz*gradrhox + pyz*gradrhoy)
#endif	
                  !I compute the new velocities
		  u_loc=hfields_s(ii,jj,kk,2,myblock) !velocity
                  v_loc=hfields_s(ii,jj,kk,3,myblock)
                  w_loc=hfields_s(ii,jj,kk,4,myblock)            
					 
#if defined(INTERFACE_INCOMP) && defined(DENSRATIO)
		  mytemp= -sharp_c*locauxfields_s(ii,jj,kk,2,myblock)
		  u_loc = u_loc + 0.5_db*forcex/(rhophi_loc)
                  v_loc = v_loc + 0.5_db*forcey/(rhophi_loc)
                  w_loc = w_loc + 0.5_db*forcez/(rhophi_loc)
				  
		  u_loc = u_loc/ &
		   (1.0_db - 0.5_db*(rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)/rhophi_loc )
				  
		  v_loc = v_loc/ &
		   (1.0_db - 0.5_db*(rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)/rhophi_loc ) 
				  
		  w_loc = w_loc/ &
	           (1.0_db - 0.5_db*(rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)/rhophi_loc ) 

		  forces_s(ii,jj,kk,1,myblock)= forces_s(ii,jj,kk,1,myblock) - &
		   (rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)*u_loc
		  forces_s(ii,jj,kk,2,myblock)= forces_s(ii,jj,kk,2,myblock) - &
		   (rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)*v_loc
		  forces_s(ii,jj,kk,3,myblock)= forces_s(ii,jj,kk,3,myblock) - &
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
                  

                  hfields_s(ii,jj,kk,2,myblock)=u_loc   !put the new velocity in hfields_s
                  hfields_s(ii,jj,kk,3,myblock)=v_loc
                  hfields_s(ii,jj,kk,4,myblock)=w_loc
                  
                                    

!regularized 
		  pxx=hfields_s(ii,jj,kk,5,myblock)
                  pyy=hfields_s(ii,jj,kk,6,myblock)
                  pzz=hfields_s(ii,jj,kk,7,myblock)
                  pxy=hfields_s(ii,jj,kk,8,myblock)
                  pxz=hfields_s(ii,jj,kk,9,myblock)
                  pyz=hfields_s(ii,jj,kk,10,myblock)
                  
                  
                  pxx=pxx + forcex*u_loc/rhophi_loc
                  pyy=pyy + forcey*v_loc/rhophi_loc
                  pzz=pzz + forcez*w_loc/rhophi_loc
                  pxy=pxy + HALF*(forcey*u_loc+forcex*v_loc)/rhophi_loc
                  pxz=pxz + HALF*(forcez*u_loc+forcex*w_loc)/rhophi_loc
                  pyz=pyz + HALF*(forcez*v_loc+forcey*w_loc)/rhophi_loc


                  hfields_s(ii,jj,kk,5,myblock)=pxx
                  hfields_s(ii,jj,kk,6,myblock)=pyy
                  hfields_s(ii,jj,kk,7,myblock)=pzz
                  hfields_s(ii,jj,kk,8,myblock)=pxy
                  hfields_s(ii,jj,kk,9,myblock)=pxz
                  hfields_s(ii,jj,kk,10,myblock)=pyz
                  
                  
#if defined(ELASTIC_FORCE)
		  u_ref(i,j,k) = u_ref(i,j,k) + &
		   lambda_rel*(u_loc - u_ref(i,j,k))
		  v_ref(i,j,k) = v_ref(i,j,k) + &
		   lambda_rel*(v_loc - v_ref(i,j,k))
		  w_ref(i,j,k) = w_ref(i,j,k) + &
		   lambda_rel*(w_loc - w_ref(i,j,k))
				  
                  forces_s(ii,jj,kk,1,myblock)= forces_s(ii,jj,kk,1,myblock) + &
                   rhophi_loc*(u_loc - u_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+ 
		  forces_s(ii,jj,kk,2,myblock)= forces_s(ii,jj,kk,2,myblock) +&
		   rhophi_loc*(v_loc - v_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+ 
		  forces_s(ii,jj,kk,3,myblock)= forces_s(ii,jj,kk,3,myblock) + &
		   rhophi_loc*(w_loc - w_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+rhophi_loc*fz  
#endif               
                  
#if defined(DENSRATIO)	
                  pxx=pxx - cssq*press_loc - u_loc*u_loc 
                  pyy=pyy - cssq*press_loc - v_loc*v_loc
                  pzz=pzz - cssq*press_loc - w_loc*w_loc 
                  pxy=pxy - u_loc*v_loc
                  pxz=pxz - u_loc*w_loc
                  pyz=pyz - v_loc*w_loc
		  
		  forces_s(ii,jj,kk,1,myblock)= forces_s(ii,jj,kk,1,myblock) - &
		   (visc_loc/(tau_loc*cssq))*(pxx*gradrhox + pxy*gradrhoy + pxz*gradrhoz)
		  forces_s(ii,jj,kk,2,myblock)= forces_s(ii,jj,kk,2,myblock) - &
		   (visc_loc/(tau_loc*cssq))*(pyy*gradrhoy + pxy*gradrhox + pyz*gradrhoz)
		  forces_s(ii,jj,kk,3,myblock)= forces_s(ii,jj,kk,3,myblock) - &
		   (visc_loc/(tau_loc*cssq))*(pzz*gradrhoz + pxz*gradrhox + pyz*gradrhoy)
#endif               


   endsubroutine moments_LB_kernel_zplus
   
   attributes(global) subroutine moments_LB_kernel_zminus(step,iprobe,jprobe,kprobe,flop,nx,ny,nz,coords,isfluid &   
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
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields,nblocks_d) :: phifields_s 
#endif           
#if defined(ELASTIC_FORCE)
      real(kind=db) :: lambda_rel,k_elastic
      real(kind=db), dimension(1:nx,1:ny,1:nz) :: u_ref,v_ref,w_ref
#endif
      real(kind=db) :: fx,fy,fz

      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields,nblocks_d) :: hfields_old,hfields_s
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nauxfields,nblocks_d) :: auxfields_s
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nlocauxfields,nblocks_d) :: locauxfields_s
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nforces,nblocks_d) :: forces_s
      

      real(kind=db) :: pxx,pyy,pzz,pxy,pxz,pyz
      real(kind=db) :: mytemp,forcex,forcey,forcez,rhophi_loc,press_loc,u_loc,v_loc,w_loc

#ifdef TWOCOMPONENT
	  real(kind=db) ::gradfix,gradfiy,gradfiz,wet_loc,phi_loc,lap_phi_loc
#endif
      
#ifdef DENSRATIO
      real(kind=db) :: gradrhox,gradrhoy,gradrhoz,omega_loc, visc_loc,tau_loc
#endif
      
      integer :: i,j,k
      integer :: gi,gj,gk
      integer :: myblock,ii,jj,kk
      
      ii=threadIdx%x
      jj=threadIdx%y
      kk=threadIdx%z
      
      i = (blockIdx%x-1) * TILE_DIMx + ii
      j = (blockIdx%y-1) * TILE_DIMy + jj
      k = kk
      
      if(abs(isfluid(i,j,k)) .ne. 1)return
      
      !gi=nx*coords(1)+i
      !gj=ny*coords(2)+j
      !gk=nz*coords(3)+k
      
      myblock=blockIdx%x+blockIdx%y*nxblock_d+1*nxyblock_d+1
      
				 
		  press_loc=hfields_s(ii,jj,kk,1,myblock)
				 
#ifdef TWOCOMPONENT					 
		  phi_loc=phifields_s(ii,jj,kk,1,myblock)
		  lap_phi_loc=locauxfields_s(ii,jj,kk,1,myblock)
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
		  mytemp=auxfields_s(ii,jj,kk,4,myblock) !modgrad
		  gradfix=auxfields_s(ii,jj,kk,1,myblock)*mytemp !normx*modgrad
		  gradfiy=auxfields_s(ii,jj,kk,2,myblock)*mytemp !normy*modgrad
		  gradfiz=auxfields_s(ii,jj,kk,3,myblock)*mytemp !normz*modgrad
		  forcex = forcex + &
                   (4.0_db*beta*phi_loc*(phi_loc-1.0_db)*(phi_loc-0.5_db) - kapp*lap_phi_loc)*gradfix
		  forcey = forcey + &
                   (4.0_db*beta*phi_loc*(phi_loc-1.0_db)*(phi_loc-0.5_db) - kapp*lap_phi_loc)*gradfiy
		  forcez = forcez + &
                   (4.0_db*beta*phi_loc*(phi_loc-1.0_db)*(phi_loc-0.5_db) - kapp*lap_phi_loc)*gradfiz
				   				   

#ifdef REPULSIVE_FLUX
		  mytemp=locauxfields_s(ii,jj,kk,6,myblock)*rhophi_loc 
		  if(abs(mytemp)>1.0d-3) mytemp=1.0d-3*sign(1.0,mytemp)!mytemp*0.1_db
		  forcex=forcex + mytemp*rhophi_loc
				  
		  mytemp=locauxfields_s(ii,jj,kk,7,myblock)*rhophi_loc 
		  if(abs(mytemp)>1.0d-3) mytemp=1.0d-3*sign(1.0,mytemp)!mytemp*0.1_db
		  forcey=forcey + mytemp*rhophi_loc
				  
		  mytemp=locauxfields_s(ii,jj,kk,8,myblock)*rhophi_loc 
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
                  
              
                  
		  forces_s(ii,jj,kk,1,myblock)=forcex
		  forces_s(ii,jj,kk,2,myblock)=forcey
		  forces_s(ii,jj,kk,3,myblock)=forcez
				  
		  u_loc=hfields_old(ii,jj,kk,2,myblock) !velocity
                  v_loc=hfields_old(ii,jj,kk,3,myblock)
                  w_loc=hfields_old(ii,jj,kk,4,myblock)
                  
                 
				  
#if defined(ELASTIC_FORCE)
                  forcex=forcex + &
                   rhophi_loc*(u_loc - u_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+   
		  forcey=forcey + &
		   rhophi_loc*(v_loc - v_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+ 
		  forcez=forcez + &
		   rhophi_loc*(w_loc - w_ref(i,j,k))*k_elastic*lambda_rel*phi_loc + rhophi_loc*fz !+  	
#endif

#ifdef DENSRATIO 
			  
                  pxx=hfields_s(ii,jj,kk,5,myblock)
                  pyy=hfields_s(ii,jj,kk,6,myblock)
                  pzz=hfields_s(ii,jj,kk,7,myblock)
                  pxy=hfields_s(ii,jj,kk,8,myblock)
                  pxz=hfields_s(ii,jj,kk,9,myblock)
                  pyz=hfields_s(ii,jj,kk,10,myblock)
                  
                  !1-2
                  !*1
                  ! 2nd order
                  pxx=pxx - cssq*press_loc - u_loc*u_loc
                  pyy=pyy - cssq*press_loc - v_loc*v_loc
                  pzz=pzz - cssq*press_loc - w_loc*w_loc
                  pxy=pxy - u_loc*v_loc
                  pxz=pxz - u_loc*w_loc
                  pyz=pyz - v_loc*w_loc

		  visc_loc=(rho_r*visc1*phi_loc+(1.0_db-phi_loc)*visc2*rho_b)/rhophi_loc
				  
                  tau_loc=(visc_loc*invcssq + HALF) !è una tau
				  
		  forcex=forcex - (visc_loc/(tau_loc*cssq))*(pxx*gradrhox + pxy*gradrhoy + pxz*gradrhoz)
		  forcey=forcey - (visc_loc/(tau_loc*cssq))*(pyy*gradrhoy + pxy*gradrhox + pyz*gradrhoz)
		  forcez=forcez - (visc_loc/(tau_loc*cssq))*(pzz*gradrhoz + pxz*gradrhox + pyz*gradrhoy)
#endif	
                  !I compute the new velocities
		  u_loc=hfields_s(ii,jj,kk,2,myblock) !velocity
                  v_loc=hfields_s(ii,jj,kk,3,myblock)
                  w_loc=hfields_s(ii,jj,kk,4,myblock)            
					 
#if defined(INTERFACE_INCOMP) && defined(DENSRATIO)
		  mytemp= -sharp_c*locauxfields_s(ii,jj,kk,2,myblock)
		  u_loc = u_loc + 0.5_db*forcex/(rhophi_loc)
                  v_loc = v_loc + 0.5_db*forcey/(rhophi_loc)
                  w_loc = w_loc + 0.5_db*forcez/(rhophi_loc)
				  
		  u_loc = u_loc/ &
		   (1.0_db - 0.5_db*(rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)/rhophi_loc )
				  
		  v_loc = v_loc/ &
		   (1.0_db - 0.5_db*(rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)/rhophi_loc ) 
				  
		  w_loc = w_loc/ &
	           (1.0_db - 0.5_db*(rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)/rhophi_loc ) 

		  forces_s(ii,jj,kk,1,myblock)= forces_s(ii,jj,kk,1,myblock) - &
		   (rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)*u_loc
		  forces_s(ii,jj,kk,2,myblock)= forces_s(ii,jj,kk,2,myblock) - &
		   (rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)*v_loc
		  forces_s(ii,jj,kk,3,myblock)= forces_s(ii,jj,kk,3,myblock) - &
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
                  

                  hfields_s(ii,jj,kk,2,myblock)=u_loc   !put the new velocity in hfields_s
                  hfields_s(ii,jj,kk,3,myblock)=v_loc
                  hfields_s(ii,jj,kk,4,myblock)=w_loc
                  
                                    

!regularized 
		  pxx=hfields_s(ii,jj,kk,5,myblock)
                  pyy=hfields_s(ii,jj,kk,6,myblock)
                  pzz=hfields_s(ii,jj,kk,7,myblock)
                  pxy=hfields_s(ii,jj,kk,8,myblock)
                  pxz=hfields_s(ii,jj,kk,9,myblock)
                  pyz=hfields_s(ii,jj,kk,10,myblock)
                  
                  
                  pxx=pxx + forcex*u_loc/rhophi_loc
                  pyy=pyy + forcey*v_loc/rhophi_loc
                  pzz=pzz + forcez*w_loc/rhophi_loc
                  pxy=pxy + HALF*(forcey*u_loc+forcex*v_loc)/rhophi_loc
                  pxz=pxz + HALF*(forcez*u_loc+forcex*w_loc)/rhophi_loc
                  pyz=pyz + HALF*(forcez*v_loc+forcey*w_loc)/rhophi_loc


                  hfields_s(ii,jj,kk,5,myblock)=pxx
                  hfields_s(ii,jj,kk,6,myblock)=pyy
                  hfields_s(ii,jj,kk,7,myblock)=pzz
                  hfields_s(ii,jj,kk,8,myblock)=pxy
                  hfields_s(ii,jj,kk,9,myblock)=pxz
                  hfields_s(ii,jj,kk,10,myblock)=pyz
                  
                  
#if defined(ELASTIC_FORCE)
		  u_ref(i,j,k) = u_ref(i,j,k) + &
		   lambda_rel*(u_loc - u_ref(i,j,k))
		  v_ref(i,j,k) = v_ref(i,j,k) + &
		   lambda_rel*(v_loc - v_ref(i,j,k))
		  w_ref(i,j,k) = w_ref(i,j,k) + &
		   lambda_rel*(w_loc - w_ref(i,j,k))
				  
                  forces_s(ii,jj,kk,1,myblock)= forces_s(ii,jj,kk,1,myblock) + &
                   rhophi_loc*(u_loc - u_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+ 
		  forces_s(ii,jj,kk,2,myblock)= forces_s(ii,jj,kk,2,myblock) +&
		   rhophi_loc*(v_loc - v_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+ 
		  forces_s(ii,jj,kk,3,myblock)= forces_s(ii,jj,kk,3,myblock) + &
		   rhophi_loc*(w_loc - w_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+rhophi_loc*fz  
#endif               
                  
#if defined(DENSRATIO)	
                  pxx=pxx - cssq*press_loc - u_loc*u_loc 
                  pyy=pyy - cssq*press_loc - v_loc*v_loc
                  pzz=pzz - cssq*press_loc - w_loc*w_loc 
                  pxy=pxy - u_loc*v_loc
                  pxz=pxz - u_loc*w_loc
                  pyz=pyz - v_loc*w_loc
		  
		  forces_s(ii,jj,kk,1,myblock)= forces_s(ii,jj,kk,1,myblock) - &
		   (visc_loc/(tau_loc*cssq))*(pxx*gradrhox + pxy*gradrhoy + pxz*gradrhoz)
		  forces_s(ii,jj,kk,2,myblock)= forces_s(ii,jj,kk,2,myblock) - &
		   (visc_loc/(tau_loc*cssq))*(pyy*gradrhoy + pxy*gradrhox + pyz*gradrhoz)
		  forces_s(ii,jj,kk,3,myblock)= forces_s(ii,jj,kk,3,myblock) - &
		   (visc_loc/(tau_loc*cssq))*(pzz*gradrhoz + pxz*gradrhox + pyz*gradrhoy)
#endif               


   endsubroutine moments_LB_kernel_zminus
   
endmodule lb_cuda_moments
