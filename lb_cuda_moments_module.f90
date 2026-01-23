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
       ,auxfields_s,locauxfields_s,forces_s &
       ,press_old,u_old,v_old,w_old,pxx_old,pyy_old,pzz_old,pxy_old,pxz_old,pyz_old &
       ,press_s,u_s,v_s,w_s,pxx_s,pyy_s,pzz_s,pxy_s,pxz_s,pyz_s)
 

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

      
      real(kind=db), dimension(ntotauxfields) :: auxfields_s
      real(kind=db), dimension(ntotlocauxfields) :: locauxfields_s
      real(kind=db), dimension(ntotforces) :: forces_s
      
      real(kind=db), dimension(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff) :: &
       press_old,u_old,v_old,w_old,pxx_old,pyy_old,pzz_old,pxy_old,pxz_old,pyz_old, &
       press_s,u_s,v_s,w_s,pxx_s,pyy_s,pzz_s,pxy_s,pxz_s,pyz_s
      

      real(kind=db) :: pxx,pyy,pzz,pxy,pxz,pyz
      real(kind=db) :: mytemp,forcex,forcey,forcez,rhophi_loc

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
      

                 
                 
#ifdef TWOCOMPONENT                     
          phi_loc=phifields_s(idx5d(ii,jj,kk,1,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields))
          lap_phi_loc=locauxfields_s(idx5d(ii,jj,kk,1,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nlocauxfields))
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
          rhophi_loc = 1.0_db !press_s(i,j,k)
          forcex=fx !0.0_db
          forcey=fy !0.0_db
          forcez=fz !0.0_db
#endif    
             
#ifdef TWOCOMPONENT        
           
          !jaqmin 
          mytemp=auxfields_s(idx5d(ii,jj,kk,4,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nauxfields)) !modgrad
          gradfix=auxfields_s(idx5d(ii,jj,kk,1,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nauxfields))*mytemp !normx*modgrad
          gradfiy=auxfields_s(idx5d(ii,jj,kk,2,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nauxfields))*mytemp !normy*modgrad
          gradfiz=auxfields_s(idx5d(ii,jj,kk,3,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nauxfields))*mytemp !normz*modgrad
          forcex = forcex + &
                   (4.0_db*beta*phi_loc*(phi_loc-1.0_db)*(phi_loc-0.5_db) - kapp*lap_phi_loc)*gradfix
          forcey = forcey + &
                   (4.0_db*beta*phi_loc*(phi_loc-1.0_db)*(phi_loc-0.5_db) - kapp*lap_phi_loc)*gradfiy
          forcez = forcez + &
                   (4.0_db*beta*phi_loc*(phi_loc-1.0_db)*(phi_loc-0.5_db) - kapp*lap_phi_loc)*gradfiz
                                      

#ifdef REPULSIVE_FLUX
          mytemp=locauxfields_s(idx5d(ii,jj,kk,6,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nlocauxfields))*rhophi_loc 
          if(abs(mytemp)>1.0d-3) mytemp=1.0d-3*sign(1.0,mytemp)!mytemp*0.1_db
          forcex=forcex + mytemp*rhophi_loc
                  
          mytemp=locauxfields_s(idx5d(ii,jj,kk,7,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nlocauxfields))*rhophi_loc 
          if(abs(mytemp)>1.0d-3) mytemp=1.0d-3*sign(1.0,mytemp)!mytemp*0.1_db
          forcey=forcey + mytemp*rhophi_loc
                  
          mytemp=locauxfields_s(idx5d(ii,jj,kk,8,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nlocauxfields))*rhophi_loc 
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
           press_s(i,j,k)*cssq*gradrhox   !+ (rhophi_loc-(rho_r+rho_b)*HALF)*fx
          forcey=forcey - &
           press_s(i,j,k)*cssq*gradrhoy   !+ (rhophi_loc-(rho_r+rho_b)*HALF)*fy
          forcez=forcez - &
           press_s(i,j,k)*cssq*gradrhoz   !+ (rhophi_loc-(rho_r+rho_b)*HALF)*fz
          !! from this point I compute the force terms that depend on the velocity
          !! these terms should be not included in force arrays since they must be computed with the updated velocity
          !! at the end of this subroutine
#endif
                  
          gif=idx5d(ii,jj,kk,1,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nforces)
          gjf=idx5d(ii,jj,kk,2,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nforces)
          gkf=idx5d(ii,jj,kk,3,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nforces)
                  
          forces_s(gif)=forcex
          forces_s(gjf)=forcey
          forces_s(gkf)=forcez
                  
                  
                 
                  
#if defined(ELASTIC_FORCE)
          forcex=forcex + &
           rhophi_loc*(u_old(i,j,k) - u_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+   
          forcey=forcey + &
           rhophi_loc*(v_old(i,j,k) - v_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+ 
          forcez=forcez + &
           rhophi_loc*(w_old(i,j,k) - w_ref(i,j,k))*k_elastic*lambda_rel*phi_loc + rhophi_loc*fz !+      
#endif

#ifdef DENSRATIO 
              
          pxx=pxx_s(i,j,k)
          pyy=pyy_s(i,j,k)
          pzz=pzz_s(i,j,k)
          pxy=pxy_s(i,j,k)
          pxz=pxz_s(i,j,k)
          pyz=pyz_s(i,j,k)
                  
                  !1-2
                  !*1
                  ! 2nd order
          pxx=pxx - cssq*press_s(i,j,k) - u_old(i,j,k)*u_old(i,j,k)
          pyy=pyy - cssq*press_s(i,j,k) - v_old(i,j,k)*v_old(i,j,k)
          pzz=pzz - cssq*press_s(i,j,k) - w_old(i,j,k)*w_old(i,j,k)
          pxy=pxy - u_old(i,j,k)*v_old(i,j,k)
          pxz=pxz - u_old(i,j,k)*w_old(i,j,k)
          pyz=pyz - v_old(i,j,k)*w_old(i,j,k)

          visc_loc=(rho_r*visc1*phi_loc+(1.0_db-phi_loc)*visc2*rho_b)/rhophi_loc
                  
          tau_loc=(visc_loc*invcssq + HALF) !è una tau
                  
          forcex=forcex - (visc_loc/(tau_loc*cssq))*(pxx*gradrhox + pxy*gradrhoy + pxz*gradrhoz)
          forcey=forcey - (visc_loc/(tau_loc*cssq))*(pyy*gradrhoy + pxy*gradrhox + pyz*gradrhoz)
          forcez=forcez - (visc_loc/(tau_loc*cssq))*(pzz*gradrhoz + pxz*gradrhox + pyz*gradrhoy)
#endif    
          !I use the new velocities          
                     
#if defined(INTERFACE_INCOMP) && defined(DENSRATIO)
          mytemp= -sharp_c*locauxfields_s(idx5d(ii,jj,kk,2,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nlocauxfields))
          u_s(i,j,k) = u_s(i,j,k) + 0.5_db*forcex/(rhophi_loc)
          v_s(i,j,k) = v_s(i,j,k) + 0.5_db*forcey/(rhophi_loc)
          w_s(i,j,k) = w_s(i,j,k) + 0.5_db*forcez/(rhophi_loc)
                  
          u_s(i,j,k) = u_s(i,j,k)/ &
           (1.0_db - 0.5_db*(rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)/rhophi_loc )
                  
          v_s(i,j,k) = v_s(i,j,k)/ &
           (1.0_db - 0.5_db*(rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)/rhophi_loc ) 
                  
          w_s(i,j,k) = w_s(i,j,k)/ &
               (1.0_db - 0.5_db*(rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)/rhophi_loc ) 

          forces_s(gif)= forces_s(gif) - &
           (rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)*u_s(i,j,k)
          forces_s(gjf)= forces_s(gjf) - &
           (rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)*v_s(i,j,k)
          forces_s(gkf)= forces_s(gkf) - &
           (rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)*w_s(i,j,k)
                     
          forcex=forcex - &
           (rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)*u_s(i,j,k)
                     
          forcey=forcey - &
           (rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)*v_s(i,j,k)
                     
          forcez=forcez - &
           (rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)*w_s(i,j,k)
                  
#else
          u_s(i,j,k) = u_s(i,j,k) + 0.5_db*forcex/rhophi_loc
          v_s(i,j,k) = v_s(i,j,k) + 0.5_db*forcey/rhophi_loc
          w_s(i,j,k) = w_s(i,j,k) + 0.5_db*forcez/rhophi_loc
#endif
                 
                  
                                    

!regularized 
                  
                  
          pxx_s(i,j,k)=pxx_s(i,j,k) + forcex*u_s(i,j,k)/rhophi_loc
          pyy_s(i,j,k)=pyy_s(i,j,k) + forcey*v_s(i,j,k)/rhophi_loc
          pzz_s(i,j,k)=pzz_s(i,j,k) + forcez*w_s(i,j,k)/rhophi_loc
          pxy_s(i,j,k)=pxy_s(i,j,k) + HALF*(forcey*u_s(i,j,k)+forcex*v_s(i,j,k))/rhophi_loc
          pxz_s(i,j,k)=pxz_s(i,j,k) + HALF*(forcez*u_s(i,j,k)+forcex*w_s(i,j,k))/rhophi_loc
          pyz_s(i,j,k)=pyz_s(i,j,k) + HALF*(forcez*v_s(i,j,k)+forcey*w_s(i,j,k))/rhophi_loc


                  
                  
#if defined(ELASTIC_FORCE)
          u_ref(i,j,k) = u_ref(i,j,k) + &
           lambda_rel*(u_s(i,j,k) - u_ref(i,j,k))
          v_ref(i,j,k) = v_ref(i,j,k) + &
           lambda_rel*(v_s(i,j,k) - v_ref(i,j,k))
          w_ref(i,j,k) = w_ref(i,j,k) + &
           lambda_rel*(w_s(i,j,k) - w_ref(i,j,k))
                  
          forces_s(gif)= forces_s(gif) + &
           rhophi_loc*(u_s(i,j,k) - u_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+ 
          forces_s(gjf)= forces_s(gjf) +&
           rhophi_loc*(v_s(i,j,k) - v_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+ 
          forces_s(gkf)= forces_s(gkf) + &
           rhophi_loc*(w_s(i,j,k) - w_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+rhophi_loc*fz  
#endif               
                  
#if defined(DENSRATIO)    
          pxx=pxx_s(i,j,k) - cssq*press_s(i,j,k) - u_s(i,j,k)*u_s(i,j,k) 
          pyy=pyy_s(i,j,k) - cssq*press_s(i,j,k) - v_s(i,j,k)*v_s(i,j,k)
          pzz=pzz_s(i,j,k) - cssq*press_s(i,j,k) - w_s(i,j,k)*w_s(i,j,k) 
          pxy=pxy_s(i,j,k) - u_s(i,j,k)*v_s(i,j,k)
          pxz=pxz_s(i,j,k) - u_s(i,j,k)*w_s(i,j,k)
          pyz=pyz_s(i,j,k) - v_s(i,j,k)*w_s(i,j,k)
          
          forces_s(gif)= forces_s(gif) - &
           (visc_loc/(tau_loc*cssq))*(pxx*gradrhox + pxy*gradrhoy + pxz*gradrhoz)
          forces_s(gjf)= forces_s(gjf) - &
           (visc_loc/(tau_loc*cssq))*(pyy*gradrhoy + pxy*gradrhox + pyz*gradrhoz)
          forces_s(gkf)= forces_s(gkf) - &
           (visc_loc/(tau_loc*cssq))*(pzz*gradrhoz + pxz*gradrhox + pyz*gradrhoy)
#endif               


   endsubroutine moments_LB_kernel  
   
endmodule lb_cuda_moments
