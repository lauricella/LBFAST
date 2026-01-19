#include "defines.h"
#if !defined(_OPENACC)  
#error "To use this module the macros _OPENACC should be defined."
#endif

module lb_cuda_fused

   use vars
   use iso_c_binding
   use cudafor
   use mpi_template, only: coords,dostop,doerror,mydev,myrank,nprocs,nbuff,nbuffbvec
   use lb_cuda_vars

   implicit none
   

contains


   attributes(global) subroutine fused_LB_kernel(step,flip,flop,nx,ny,nz,coords,isfluid &  
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
       ,visc1,omega,fx,fy,fz,ntothfields,ntotphifields,ntotauxfields,ntotlocauxfields,ntotforces, &
       hfields_in,hfields_out,auxfields_s,locauxfields_s,forces_s)

      implicit none
      
      integer :: step,flip,flop,nx,ny,nz
      
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
      real(kind=db) :: mytemp,forcex,forcey,forcez,rhophi_loc,press_loc,uu,udotc
      
      real(kind=db), shared :: f1(0:TILE_DIMx_d+1,0:TILE_DIMy_d+1,0:TILE_DIMz_d+1)
      real(kind=db), shared :: f2(0:TILE_DIMx_d+1,0:TILE_DIMy_d+1,0:TILE_DIMz_d+1)
      real(kind=db), shared :: f3(0:TILE_DIMx_d+1,0:TILE_DIMy_d+1,0:TILE_DIMz_d+1)
      real(kind=db), shared :: f4(0:TILE_DIMx_d+1,0:TILE_DIMy_d+1,0:TILE_DIMz_d+1)
  
      real(kind=db) :: F_discr,fneq1,feq,fpost
#ifdef TWOCOMPONENT
      real(kind=db) :: wet_loc
#endif
      real(kind=db) :: omega_loc,phi_loc,visc_loc
#ifdef SMAGORINSKI
	  real(kind=db) :: QQ
#endif

      integer :: i,j,k,myblock,imio,l
      integer :: ii,jj,kk
      integer :: li,lj,lk
      integer :: lii,ljj,lkk
      integer :: xblock,yblock,zblock
      integer :: gi,gj,gk,intblock
      
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

      !i = (blockIdx%x-1) * TILE_DIMx_d + threadIdx%x
      !j = (blockIdx%y-1) * TILE_DIMy_d + threadIdx%y
      !k = (blockIdx%z-1) * TILE_DIMz_d + threadIdx%z
      
!      gi=nx*coords(1)+i
!      gj=ny*coords(2)+j
!      gk=nz*coords(3)+k
      
      intblock=blockIdx%x+blockIdx%y*nxblock_d+blockIdx%z*nxyblock_d+1 !internal-node block


               !if (abs(isfluid(i,j,k)) /= 1) return
       
     
#ifdef TWOCOMPONENT	  
                  phi_loc=phifields_s(idx5d(ii,jj,kk,1,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nphifields))
#endif                  
#ifdef DENSRATIO
                  
                  rhophi_loc = rho_r*phi_loc+(ONE-phi_loc)*rho_b 
#else
                  rhophi_loc = 1.0_db !press_loc
#endif	

				  forcex=forces_s(idx5d(ii,jj,kk,1,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nforces))
				  forcey=forces_s(idx5d(ii,jj,kk,2,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nforces))
				  forcez=forces_s(idx5d(ii,jj,kk,3,myblock,TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,nforces))
                  
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
                  
                  
#ifdef INTERNAL_OBSTACLES
                  if(isfluid(i,j,k) == 0)then
                    forcex=ZERO
                    forcey=ZERO
                    forcez=ZERO
                    press=ZERO
                    u=ZERO
                    v=ZERO
                    w=ZERO
                    pxx=ZERO
                    pyy=ZERO
                    pzz=ZERO
                    pxy=ZERO
                    pxz=ZERO
                    pyz=ZERO
                  endif
#endif
                  
#ifdef EXPLICITEQ 
				  uu=HALF*(u*u+v*v+w*w)*invcssq
				  
                  do lii=1,nlinks
                     udotc=(u*dex(lii) + v*dey(lii)+ w*dez(lii))*invcssq
		     feq=p(lii)*(press + (udotc+0.5_db*udotc*udotc - uu))
                     !fneq1=(HALF/(cssq*cssq))*( (dex(lii)*dex(lii)-cssq)*pxx &
		     ! + (dey(lii)*dey(lii)-cssq)*pyy + (dez(lii)*dez(lii)-cssq)*pzz &
	             ! + TWO*(dex(lii)*dey(lii))*pxy + TWO*(dex(lii)*dez(lii))*pxz &
		     ! + TWO*(dey(lii)*dez(lii))*pyz)
                     fpost=feq!+fneq1
                     pxx=pxx - fpost*(dex(lii)*dex(lii))
                     pyy=pyy - fpost*(dey(lii)*dey(lii))
                     pzz=pzz - fpost*(dez(lii)*dez(lii))
                     pxy=pxy - fpost*(dex(lii)*dey(lii))
                     pxz=pxz - fpost*(dex(lii)*dez(lii))
                     pyz=pyz - fpost*(dey(lii)*dez(lii))
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

!                  rho(i,j,k) = f(i,j,k,0)+f(i,j,k,1)+f(i,j,k,2)+f(i,j,k,3)+f(i,j,k,4)+f(i,j,k,5) &
!                     +f(i,j,k,6)+f(i,j,k,7)+f(i,j,k,8)+f(i,j,k,9)+f(i,j,k,10)+f(i,j,k,11) &
!                     +f(i,j,k,12)+f(i,j,k,13)+f(i,j,k,14)+f(i,j,k,15)+f(i,j,k,16)+f(i,j,k,17) &
!                     +f(i,j,k,18)+f(i,j,k,19)+f(i,j,k,20)+f(i,j,k,21)+f(i,j,k,22)+f(i,j,k,23)+f(i,j,k,24) &
!                     +f(i,j,k,25) +f(i,j,k,26)
                 
!                 u(i,j,k) = ((f(i,j,k,1)+f(i,j,k,7)+f(i,j,k,9)+f(i,j,k,15)+f(i,j,k,18)+f(i,j,k,19)+f(i,j,k,21)+f(i,j,k,24)+f(i,j,k,25)) &
!                     -(f(i,j,k,2)+f(i,j,k,8)+f(i,j,k,10)+f(i,j,k,16)+f(i,j,k,17)+f(i,j,k,20)+f(i,j,k,22)+f(i,j,k,23)+f(i,j,k,26)))

!                  v(i,j,k) = ((f(i,j,k,3)+f(i,j,k,7)+f(i,j,k,10)+f(i,j,k,11)+f(i,j,k,13)+f(i,j,k,19)+f(i,j,k,22)+f(i,j,k,24)+f(i,j,k,26)) &
!                     -(f(i,j,k,4)+f(i,j,k,8)+f(i,j,k,9)+f(i,j,k,12)+f(i,j,k,14)+f(i,j,k,20)+f(i,j,k,21)+f(i,j,k,23)+f(i,j,k,25)))

!                  w(i,j,k) = ((f(i,j,k,5)+f(i,j,k,11)+f(i,j,k,14)+f(i,j,k,15)+f(i,j,k,17)+f(i,j,k,19)+f(i,j,k,21)+f(i,j,k,23)+f(i,j,k,26)) &
!                     -(f(i,j,k,6)+f(i,j,k,12)+f(i,j,k,13)+f(i,j,k,16)+f(i,j,k,18)+f(i,j,k,20)+f(i,j,k,22)+f(i,j,k,24)+f(i,j,k,25)))
                 
!                 !total flux tensor
!                 pxx(i,j,k)=f(i,j,k,1)+f(i,j,k,2)+f(i,j,k,7)+f(i,j,k,8) &
!                  +f(i,j,k,9)+f(i,j,k,10)+f(i,j,k,15)+f(i,j,k,16)+f(i,j,k,17)+f(i,j,k,18) &
!                  +f(i,j,k,19)+f(i,j,k,20)+f(i,j,k,21)+f(i,j,k,22)+f(i,j,k,23)+f(i,j,k,24) &
!                  +f(i,j,k,25)+f(i,j,k,26)
!                 pyy(i,j,k)=f(i,j,k,3)+f(i,j,k,4)+f(i,j,k,7)+f(i,j,k,8) &
!                  +f(i,j,k,9)+f(i,j,k,10)+f(i,j,k,11)+f(i,j,k,12)+f(i,j,k,13)+f(i,j,k,14) &
!                  +f(i,j,k,19)+f(i,j,k,20)+f(i,j,k,21)+f(i,j,k,22)+f(i,j,k,23)+f(i,j,k,24) &
!                  +f(i,j,k,25)+f(i,j,k,26)
!                 pzz(i,j,k)=f(i,j,k,5)+f(i,j,k,6)+f(i,j,k,11)+f(i,j,k,12)+f(i,j,k,13)+f(i,j,k,14) &
!                  +f(i,j,k,15)+f(i,j,k,16)+f(i,j,k,17)+f(i,j,k,18) &
!                  +f(i,j,k,19)+f(i,j,k,20)+f(i,j,k,21)+f(i,j,k,22)+f(i,j,k,23)+f(i,j,k,24) &
!                  +f(i,j,k,25)+f(i,j,k,26)
!                 pxy(i,j,k)=(f(i,j,k,7)+f(i,j,k,8)+f(i,j,k,19)+f(i,j,k,20)+f(i,j,k,23)+f(i,j,k,24)) &
!                  -(f(i,j,k,9)+f(i,j,k,10)+f(i,j,k,21)+f(i,j,k,22)+f(i,j,k,25)+f(i,j,k,26))
!                 pxz(i,j,k)=(f(i,j,k,15)+f(i,j,k,16)+f(i,j,k,19)+f(i,j,k,20)+f(i,j,k,21)+f(i,j,k,22)) &
!                  -(f(i,j,k,17)+f(i,j,k,18)+f(i,j,k,23)+f(i,j,k,24)+f(i,j,k,25)+f(i,j,k,26))
!                 pyz(i,j,k)=(f(i,j,k,11)+f(i,j,k,12)+f(i,j,k,19)+f(i,j,k,20)+f(i,j,k,25)+f(i,j,k,26)) &
!                  -(f(i,j,k,13)+f(i,j,k,14)+f(i,j,k,21)+f(i,j,k,22)+f(i,j,k,23)+f(i,j,k,24))
                                     
                  opress=ZERO
                  ou=ZERO
                  ov=ZERO
                  ow=ZERO
                  opxx=ZERO
                  opyy=ZERO
                  opzz=ZERO
                  opxy=ZERO
                  opxz=ZERO
                  opyz=ZERO
!!!!!!!!!!!!!!!!!!!!!!!!!!0

#ifdef SECOND_ORDER
			      feq=(4.0_db*(2.0_db*press - 3.0_db &
			       *(u**2.0_db + v**2.0_db + w**2.0_db)))/27.0_db

#else
!0
			      feq=(8.0_db*press - 3.0_db*(4.0_db*w**2.0_db &
			       + v**2.0_db*(4.0_db - 6.0_db*w**2.0_db) &
			       + u**2.0_db*(-2.0_db + 3.0_db*v**2.0_db)*(-2.0_db &
			       + 3.0_db*w**2.0_db)))/27.0_db
!0
#endif 

				  fneq1=(-3.0_db*(pxx + pyy + pzz))/2.0_db


				  F_discr=(-8.0_db*(forcex*u + forcey*v &
				   + forcez*w))/(9.0_db*rhophi_loc)
                  opress=feq + (1.0_db-omega_loc)*fneq1*p0 + 0.5_db*(F_discr)

!!!!!!!!!!!!!!!!!!!!!!!!!!1
                  lii=li+1
                  ljj=lj
                  lkk=lk
                  lii=mod(lii+TILE_DIMx_d+2,(TILE_DIMx_d+2))
                  !ljj=mod(ljj+TILE_DIMy_d+2,(TILE_DIMy_d+2))
                  !lkk=mod(lkk+TILE_DIMz_d+2,(TILE_DIMz_d+2))
#ifdef SECOND_ORDER
			      feq=(2.0_db*press + 6.0_db*u &
			       *(1.0_db + u) - 3.0_db*v**2.0_db &
			       - 3.0_db*w**2.0_db)/27.0_db
#else
!1
			      feq=(4.0_db*press + 3.0_db*(4.0_db*u*(1.0_db &
			       + u) - 2.0_db*(1.0_db + 3.0_db*u*(1.0_db &
			       + u))*v**2.0_db + (1.0_db + 3.0_db*u*(1.0_db &
			       + u))*(-2.0_db + 3.0_db*v**2.0_db)*w**2.0_db))/54.0_db
!1
#endif
				  fneq1=(3.0_db*(2.0_db*pxx - pyy - pzz))/2.0_db


				  F_discr=(2.0_db*(forcex + 2.0_db*forcex*u - forcey*v &
				   - forcez*w))/(9.0_db*rhophi_loc)
                  f1(lii,ljj,lkk)=feq + (1.0_db-omega_loc)*fneq1*p1 + 0.5_db*(F_discr)


!!!!!!!!!!!!!!!!!!!!!!!!!!2
                  lii=li-1
                  ljj=lj
                  lkk=lk
                  lii=mod(lii+TILE_DIMx_d+2,(TILE_DIMx_d+2))
                  !ljj=mod(ljj+TILE_DIMy_d+2,(TILE_DIMy_d+2))
                  !lkk=mod(lkk+TILE_DIMz_d+2,(TILE_DIMz_d+2))
#ifdef SECOND_ORDER
			      feq=(2.0_db*press + 6.0_db*(-1.0_db &
			       + u)*u - 3.0_db*v**2.0_db &
			       - 3.0_db*w**2.0_db)/27.0_db
#else
!2
			      feq=(4.0_db*press + 3.0_db*(4.0_db*(-1.0_db &
			       + u)*u - 2.0_db*(1.0_db + 3.0_db*(-1.0_db &
			       + u)*u)*v**2.0_db + (1.0_db &
			       + 3.0_db*(-1.0_db + u)*u)*(-2.0_db &
			       + 3.0_db*v**2.0_db)*w**2.0_db))/54.0_db
!2
#endif
				  fneq1=(3.0_db*(2.0_db*pxx - pyy - pzz))/2.0_db


				  F_discr=(-2.0_db*(forcex - 2.0_db*forcex*u + forcey*v &
				   + forcez*w))/(9.0_db*rhophi_loc)
                  f2(lii,ljj,lkk)=feq + (1.0_db-omega_loc)*fneq1*p1 + 0.5_db*(F_discr)


!!!!!!!!!!!!!!!!!!!!!!!!!!3
                  lii=li
                  ljj=lj+1
                  lkk=lk
                  !lii=mod(lii+TILE_DIMx_d+2,(TILE_DIMx_d+2))
                  ljj=mod(ljj+TILE_DIMy_d+2,(TILE_DIMy_d+2))
                  !lkk=mod(lkk+TILE_DIMz_d+2,(TILE_DIMz_d+2))
#ifdef SECOND_ORDER
			      feq=(2.0_db*press - 3.0_db*(u**2.0_db &
			       - 2.0_db*v*(1.0_db + v) + w**2.0_db))/27.0_db
#else
!3
			      feq=(4.0_db*press + 12.0_db*v*(1.0_db &
			       + v) - 6.0_db*u**2.0_db*(1.0_db &
			       + 3.0_db*v*(1.0_db + v)) + 3.0_db*(-2.0_db &
			       + 3.0_db*u**2.0_db)*(1.0_db + 3.0_db*v*(1.0_db &
			       + v))*w**2.0_db)/54.0_db
!3
#endif
				  fneq1=(-3.0_db*(pxx - 2.0_db*pyy + pzz))/2.0_db


				  F_discr=(2.0_db*(forcey - forcex*u + 2.0_db*forcey*v &
				   - forcez*w))/(9.0_db*rhophi_loc)
                  f3(lii,ljj,lkk)=feq + (1.0_db-omega_loc)*fneq1*p1 + 0.5_db*(F_discr)


!!!!!!!!!!!!!!!!!!!!!!!!!!4
                  lii=li
                  ljj=lj-1
                  lkk=lk
                  !lii=mod(lii+TILE_DIMx_d+2,(TILE_DIMx_d+2))
                  ljj=mod(ljj+TILE_DIMy_d+2,(TILE_DIMy_d+2))
                  !lkk=mod(lkk+TILE_DIMz_d+2,(TILE_DIMz_d+2))
#ifdef SECOND_ORDER
			      feq=(2.0_db*press - 3.0_db*(u**2.0_db &
			      - 2.0_db*(-1.0_db + v)*v + w**2.0_db))/27.0_db
#else
!4
			      feq=(4.0_db*press + 12.0_db*(-1.0_db &
			       + v)*v - 6.0_db*u**2.0_db*(1.0_db &
			       + 3.0_db*(-1.0_db + v)*v) + 3.0_db*(-2.0_db &
			       + 3.0_db*u**2.0_db)*(1.0_db + 3.0_db*(-1.0_db &
			       + v)*v)*w**2.0_db)/54.0_db
!4
#endif
				  fneq1=(-3.0_db*(pxx - 2.0_db*pyy + pzz))/2.0_db


				  F_discr=(-2.0_db*(forcey + forcex*u - 2.0_db*forcey*v &
				   + forcez*w))/(9.0_db*rhophi_loc)
                  f4(lii,ljj,lkk)=feq + (1.0_db-omega_loc)*fneq1*p1 + 0.5_db*(F_discr)
                  
                  call syncthreads
                            
                  opress = opress + (f1(li,lj,lk)+f2(li,lj,lk)+f3(li,lj,lk)+f4(li,lj,lk))
                  ou = f1(li,lj,lk)-f2(li,lj,lk)
                  ov = f3(li,lj,lk)-f4(li,lj,lk)
                  opxx = f1(li,lj,lk)+f2(li,lj,lk)
                  opyy = f3(li,lj,lk)+f4(li,lj,lk)
                        
                  call syncthreads
                
!!!!!!!!!!!!!!!!!!!!!!!!!!5
                  lii=li
                  ljj=lj
                  lkk=lk+1
                  !lii=mod(lii+TILE_DIMx_d+2,(TILE_DIMx_d+2))
                  !ljj=mod(ljj+TILE_DIMy_d+2,(TILE_DIMy_d+2))
                  lkk=mod(lkk+TILE_DIMz_d+2,(TILE_DIMz_d+2))
#ifdef SECOND_ORDER
			      feq=(2.0_db*press - 3.0_db*(u**2.0_db &
			       + v**2.0_db - 2.0_db*w*(1.0_db + w)))/27.0_db
#else
!5
			      feq=(4.0_db*press + 3.0_db*(4.0_db*w*(1.0_db &
			       + w) - 2.0_db*v**2.0_db*(1.0_db + 3.0_db*w*(1.0_db &
			       + w)) + u**2.0_db*(-2.0_db + 3.0_db*v**2.0_db)*(1.0_db &
			       + 3.0_db*w*(1.0_db + w))))/54.0_db
!5
#endif
				  fneq1=(-3.0_db*(pxx + pyy - 2.0_db*pzz))/2.0_db


				  F_discr=(2.0_db*(forcez - forcex*u - forcey*v &
				   + 2.0_db*forcez*w))/(9.0_db*rhophi_loc)
                  f1(lii,ljj,lkk)= feq + (1.0_db-omega_loc)*fneq1*p1 + 0.5_db*(F_discr)


!!!!!!!!!!!!!!!!!!!!!!!!!!6
                  lii=li
                  ljj=lj
                  lkk=lk-1
                  !lii=mod(lii+TILE_DIMx_d+2,(TILE_DIMx_d+2))
                  !ljj=mod(ljj+TILE_DIMy_d+2,(TILE_DIMy_d+2))
                  lkk=mod(lkk+TILE_DIMz_d+2,(TILE_DIMz_d+2))
#ifdef SECOND_ORDER
			      feq=(2.0_db*press - 3.0_db*(u**2.0_db &
			       + v**2.0_db - 2.0_db*(-1.0_db + w)*w))/27.0_db
#else
!6
			      feq=(4.0_db*press + 3.0_db*(4.0_db*(-1.0_db &
			       + w)*w + v**2.0_db*(-2.0_db - 6.0_db*(-1.0_db &
			       + w)*w) + u**2.0_db*(-2.0_db &
			       + 3.0_db*v**2.0_db)*(1.0_db + 3.0_db*(-1.0_db &
			       + w)*w)))/54.0_db
!6
#endif
				  fneq1=(-3.0_db*(pxx + pyy - 2.0_db*pzz))/2.0_db


				  F_discr=(-2.0_db*(forcez + forcex*u + forcey*v &
				   - 2.0_db*forcez*w))/(9.0_db*rhophi_loc)
                  f2(lii,ljj,lkk)=feq + (1.0_db-omega_loc)*fneq1*p1 + 0.5_db*(F_discr)


!!!!!!!!!!!!!!!!!!!!!!!!!!7
                  lii=li+1
                  ljj=lj+1
                  lkk=lk
                  lii=mod(lii+TILE_DIMx_d+2,(TILE_DIMx_d+2))
                  ljj=mod(ljj+TILE_DIMy_d+2,(TILE_DIMy_d+2))
                  !lkk=mod(lkk+TILE_DIMz_d+2,(TILE_DIMz_d+2))
#ifdef SECOND_ORDER
			      feq=(2.0_db*press + 6.0_db*(u &
			       + u**2.0_db + v + 3.0_db*u*v &
			       + v**2.0_db) - 3.0_db*w**2.0_db)/108.0_db
#else
!7
			      feq=(2.0_db*press + 6.0_db*(u &
			       + v + v**2.0_db + 3.0_db*u*v*(1.0_db &
			       + v) + u**2.0_db*(1.0_db + 3.0_db*v*(1.0_db &
			       + v))) - 3.0_db*(1.0_db + 3.0_db*u*(1.0_db &
			       + u))*(1.0_db + 3.0_db*v*(1.0_db &
			       + v))*w**2.0_db)/108.0_db
!7
#endif
				  fneq1=3.0_db*(pxx + 3.0_db*pxy + pyy) &
				   - (3.0_db*pzz)/2.0_db


				  F_discr=(forcex + forcey + 2.0_db*forcex*u &
				   + 3.0_db*forcey*u + 3.0_db*forcex*v &
				   + 2.0_db*forcey*v - forcez*w)/(18.0_db*rhophi_loc)
                  f3(lii,ljj,lkk)=feq + (1.0_db-omega_loc)*fneq1*p2 + 0.5_db*(F_discr)


!!!!!!!!!!!!!!!!!!!!!!!!!!8
                  lii=li-1
                  ljj=lj-1
                  lkk=lk
                  lii=mod(lii+TILE_DIMx_d+2,(TILE_DIMx_d+2))
                  ljj=mod(ljj+TILE_DIMy_d+2,(TILE_DIMy_d+2))
                  !lkk=mod(lkk+TILE_DIMz_d+2,(TILE_DIMz_d+2))
#ifdef SECOND_ORDER
			      feq=(2.0_db*press + 6.0_db*(u**2.0_db &
			       + (-1.0_db + v)*v + u*(-1.0_db &
			       + 3.0_db*v)) - 3.0_db*w**2.0_db)/108.0_db
#else
!8
			      feq=(2.0_db*press + 6.0_db*((-1.0_db &
			       + v)*v + u*(-1.0_db - 3.0_db*(-1.0_db &
			       + v)*v) + u**2.0_db*(1.0_db + 3.0_db*(-1.0_db &
			       + v)*v)) - 3.0_db*(1.0_db + 3.0_db*(-1.0_db &
			       + u)*u)*(1.0_db + 3.0_db*(-1.0_db &
			       + v)*v)*w**2.0_db)/108.0_db
!8
#endif 
				  fneq1=3.0_db*(pxx + 3.0_db*pxy + pyy) &
				   - (3.0_db*pzz)/2.0_db


				  F_discr=(forcex + forcey - 2.0_db*forcex*u &
				   - 3.0_db*forcey*u - 3.0_db*forcex*v &
				   - 2.0_db*forcey*v + forcez*w)/(-18.0_db*rhophi_loc)
                  f4(lii,ljj,lkk)=feq + (1.0_db-omega_loc)*fneq1*p2 + 0.5_db*(F_discr)
 
                  call syncthreads
                       
                  opress = opress + (f1(li,lj,lk)+f2(li,lj,lk)+f3(li,lj,lk)+f4(li,lj,lk))
                  ou = ou + f3(li,lj,lk)-f4(li,lj,lk)
                  ov = ov + f3(li,lj,lk)-f4(li,lj,lk)
                  ow = f1(li,lj,lk)-f2(li,lj,lk)
                  opxx = opxx + f3(li,lj,lk)+f4(li,lj,lk)
                  opyy = opyy + f3(li,lj,lk)+f4(li,lj,lk)
                  opzz = f1(li,lj,lk)+f2(li,lj,lk)
                  opxy = f3(li,lj,lk)+f4(li,lj,lk)
!                   if(gi==1 .and. gj==2 .and. gk==16 .and. myblock==intblock)then
!                    write(*,*)'5',f1(li,lj,lk),step
!                    write(*,*)'6',f2(li,lj,lk),step
!                  endif

                  call syncthreads

!!!!!!!!!!!!!!!!!!!!!!!!!!9
                  lii=li+1
                  ljj=lj-1
                  lkk=lk
                  lii=mod(lii+TILE_DIMx_d+2,(TILE_DIMx_d+2))
                  ljj=mod(ljj+TILE_DIMy_d+2,(TILE_DIMy_d+2))
                  !lkk=mod(lkk+TILE_DIMz_d+2,(TILE_DIMz_d+2))
#ifdef SECOND_ORDER
			      feq=(2.0_db*press + 6.0_db*(u &
			       + u**2.0_db - 3.0_db*u*v + (-1.0_db &
			       + v)*v) - 3.0_db*w**2.0_db)/108.0_db
#else
!9
			      feq=(2.0_db*press + 6.0_db*(u &
			       + (-1.0_db + v)*v + 3.0_db*u*(-1.0_db &
			       + v)*v + u**2.0_db*(1.0_db + 3.0_db*(-1.0_db &
			       + v)*v)) - 3.0_db*(1.0_db &
			       + 3.0_db*u*(1.0_db + u))*(1.0_db + 3.0_db*(-1.0_db &
			       + v)*v)*w**2.0_db)/108.0_db
!9
#endif
				  fneq1=3.0_db*(pxx - 3.0_db*pxy + pyy) &
				   - (3.0_db*pzz)/2.0_db


				  F_discr=(forcey + 3.0_db*forcey*u - 2.0_db*forcey*v &
				   + forcex*(-1.0_db - 2.0_db*u + 3.0_db*v) &
				   + forcez*w)/(-18.0_db*rhophi_loc)
                  f1(lii,ljj,lkk)=feq + (1.0_db-omega_loc)*fneq1*p2 + 0.5_db*(F_discr)


!!!!!!!!!!!!!!!!!!!!!!!!!!10
                  lii=li-1
                  ljj=lj+1
                  lkk=lk
                  lii=mod(lii+TILE_DIMx_d+2,(TILE_DIMx_d+2))
                  ljj=mod(ljj+TILE_DIMy_d+2,(TILE_DIMy_d+2))
                  !lkk=mod(lkk+TILE_DIMz_d+2,(TILE_DIMz_d+2))
#ifdef SECOND_ORDER
			      feq=(2.0_db*press + 6.0_db*((-1.0_db &
			       + u)*u + v - 3.0_db*u*v &
			       + v**2.0_db) - 3.0_db*w**2.0_db)/108.0_db
#else
!10
			      feq=(2.0_db*press + 6.0_db*((-1.0_db &
			       + u)*u + v + 3.0_db*(-1.0_db &
			       + u)*u*v + (1.0_db + 3.0_db*(-1.0_db &
			       + u)*u)*v**2.0_db) - 3.0_db*(1.0_db &
			       + 3.0_db*(-1.0_db + u)*u)*(1.0_db + 3.0_db*v*(1.0_db &
			       + v))*w**2.0_db)/108.0_db
!10
#endif
				  fneq1=3.0_db*(pxx - 3.0_db*pxy + pyy) &
				   - (3.0_db*pzz)/2.0_db


				  F_discr=(forcex - forcey - 2.0_db*forcex*u &
				   + 3.0_db*forcey*u + 3.0_db*forcex*v &
				   - 2.0_db*forcey*v + forcez*w)/(-18.0_db*rhophi_loc)
                  f2(lii,ljj,lkk)=feq + (1.0_db-omega_loc)*fneq1*p2 + 0.5_db*(F_discr)


!!!!!!!!!!!!!!!!!!!!!!!!!!11
                  lii=li
                  ljj=lj+1
                  lkk=lk+1
                  !lii=mod(lii+TILE_DIMx_d+2,(TILE_DIMx_d+2))
                  ljj=mod(ljj+TILE_DIMy_d+2,(TILE_DIMy_d+2))
                  lkk=mod(lkk+TILE_DIMz_d+2,(TILE_DIMz_d+2))
#ifdef SECOND_ORDER
			      feq=(2.0_db*press - 3.0_db*u**2.0_db &
			       + 6.0_db*(v + v**2.0_db + w &
			       + 3.0_db*v*w + w**2.0_db))/108.0_db
#else
!11
			      feq=(2.0_db*press - 3.0_db*u**2.0_db*(1.0_db &
			       + 3.0_db*v*(1.0_db + v))*(1.0_db + 3.0_db*w*(1.0_db &
			       + w)) + 6.0_db*(v + w + w**2.0_db &
			       + 3.0_db*v*w*(1.0_db + w) + v**2.0_db*(1.0_db &
			       + 3.0_db*w*(1.0_db + w))))/108.0_db
!11
#endif
				  fneq1=(-3.0_db*pxx)/2.0_db + 3.0_db*(pyy &
				   + 3.0_db*pyz + pzz)


				  F_discr=(forcey + forcez - forcex*u + 2.0_db*forcey*v &
				   + 3.0_db*forcez*v + 3.0_db*forcey*w &
				   + 2.0_db*forcez*w)/(18.0_db*rhophi_loc)
                  f3(lii,ljj,lkk)=feq + (1.0_db-omega_loc)*fneq1*p2 + 0.5_db*(F_discr)


!!!!!!!!!!!!!!!!!!!!!!!!!!12
                  lii=li
                  ljj=lj-1
                  lkk=lk-1
                  !lii=mod(lii+TILE_DIMx_d+2,(TILE_DIMx_d+2))
                  ljj=mod(ljj+TILE_DIMy_d+2,(TILE_DIMy_d+2))
                  lkk=mod(lkk+TILE_DIMz_d+2,(TILE_DIMz_d+2))
#ifdef SECOND_ORDER
			      feq=(2.0_db*press - 3.0_db*u**2.0_db &
			       + 6.0_db*(v**2.0_db + (-1.0_db &
			       + w)*w + v*(-1.0_db + 3.0_db*w)))/108.0_db
#else
!12
			      feq=(2.0_db*press + 3.0_db*(2.0_db*(-1.0_db &
			       + w)*w + v*(-2.0_db - 6.0_db*(-1.0_db &
			       + w)*w) - u**2.0_db*(1.0_db + 3.0_db*(-1.0_db &
			       + v)*v)*(1.0_db + 3.0_db*(-1.0_db + w)*w) &
			       + v**2.0_db*(2.0_db + 6.0_db*(-1.0_db + w)*w)))/108.0_db
!12
#endif
				  fneq1=(-3.0_db*pxx)/2.0_db + 3.0_db*(pyy &
				   + 3.0_db*pyz + pzz)


				  F_discr=(forcey + forcez + forcex*u - 2.0_db*forcey*v &
				   - 3.0_db*forcez*v - 3.0_db*forcey*w &
				   - 2.0_db*forcez*w)/(-18.0_db*rhophi_loc)
                  f4(lii,ljj,lkk)=feq + (1.0_db-omega_loc)*fneq1*p2 + 0.5_db*(F_discr)
                 
                  call syncthreads
                  
                  opress = opress + (f1(li,lj,lk)+f2(li,lj,lk)+f3(li,lj,lk)+f4(li,lj,lk))
                  ou = ou + f1(li,lj,lk)-f2(li,lj,lk)
                  ov = ov - f1(li,lj,lk)+f2(li,lj,lk)+f3(li,lj,lk)-f4(li,lj,lk)
                  ow = ow + f3(li,lj,lk)-f4(li,lj,lk)
                  opxx = opxx + f1(li,lj,lk)+f2(li,lj,lk)
                  opyy = opyy + f1(li,lj,lk)+f2(li,lj,lk)+f3(li,lj,lk)+f4(li,lj,lk)
                  opzz = opzz + f3(li,lj,lk)+f4(li,lj,lk)
                  opxy = opxy - f1(li,lj,lk)-f2(li,lj,lk)
                  opyz = f3(li,lj,lk)+f4(li,lj,lk)
!                  if(gi==1 .and. gj==2 .and. gk==16 .and. myblock==intblock)then
!                    write(*,*)'11',f3(li,lj,lk),step
!                    write(*,*)'12',f4(li,lj,lk),step
!                  endif
                  call syncthreads

!!!!!!!!!!!!!!!!!!!!!!!!!!13
                  lii=li
                  ljj=lj+1
                  lkk=lk-1
                  !lii=mod(lii+TILE_DIMx_d+2,(TILE_DIMx_d+2))
                  ljj=mod(ljj+TILE_DIMy_d+2,(TILE_DIMy_d+2))
                  lkk=mod(lkk+TILE_DIMz_d+2,(TILE_DIMz_d+2))
#ifdef SECOND_ORDER
			      feq=(2.0_db*press - 3.0_db*u**2.0_db &
			       + 6.0_db*(v + v**2.0_db &
			       - 3.0_db*v*w + (-1.0_db + w)*w))/108.0_db
#else
!13
			      feq=(2.0_db*press - 3.0_db*u**2.0_db*(1.0_db &
			       + 3.0_db*v*(1.0_db + v))*(1.0_db + 3.0_db*(-1.0_db &
			       + w)*w) + 6.0_db*(v + (-1.0_db &
			       + w)*w + 3.0_db*v*(-1.0_db + w)*w &
			       + v**2.0_db*(1.0_db + 3.0_db*(-1.0_db + w)*w)))/108.0_db
!13
#endif
				  fneq1=(-3.0_db*pxx)/2.0_db + 3.0_db*(pyy &
				   - 3.0_db*pyz + pzz)


				  F_discr=(forcez + forcex*u + 3.0_db*forcez*v &
				   - 2.0_db*forcez*w + forcey*(-1.0_db - 2.0_db*v &
				   + 3.0_db*w))/(-18.0_db*rhophi_loc)
                  f1(lii,ljj,lkk)=feq + (1.0_db-omega_loc)*fneq1*p2 + 0.5_db*(F_discr)


!!!!!!!!!!!!!!!!!!!!!!!!!!14
                  lii=li
                  ljj=lj-1
                  lkk=lk+1
                  !lii=mod(lii+TILE_DIMx_d+2,(TILE_DIMx_d+2))
                  ljj=mod(ljj+TILE_DIMy_d+2,(TILE_DIMy_d+2))
                  lkk=mod(lkk+TILE_DIMz_d+2,(TILE_DIMz_d+2))
#ifdef SECOND_ORDER
			      feq=(2.0_db*press - 3.0_db*u**2.0_db &
			       + 6.0_db*((-1.0_db + v)*v + w &
			       - 3.0_db*v*w + w**2.0_db))/108.0_db
#else
!14
			      feq=(2.0_db*press + 6.0_db*((-1.0_db &
			       + v)*v + w + 3.0_db*(-1.0_db &
			       + v)*v*w + (1.0_db + 3.0_db*(-1.0_db &
			       + v)*v)*w**2.0_db) &
			       - 3.0_db*u**2.0_db*(1.0_db + 3.0_db*(-1.0_db &
			       + v)*v)*(1.0_db + 3.0_db*w*(1.0_db &
			       + w)))/108.0_db
!14
#endif
				  fneq1=(-3.0_db*pxx)/2.0_db + 3.0_db*(pyy &
				   - 3.0_db*pyz + pzz)


				  F_discr=(forcey - forcez + forcex*u - 2.0_db*forcey*v &
				   + 3.0_db*forcez*v + 3.0_db*forcey*w &
				   - 2.0_db*forcez*w)/(-18.0_db*rhophi_loc)
                  f2(lii,ljj,lkk)=feq + (1.0_db-omega_loc)*fneq1*p2 + 0.5_db*(F_discr)


!!!!!!!!!!!!!!!!!!!!!!!!!!15
                  lii=li+1
                  ljj=lj
                  lkk=lk+1
                  lii=mod(lii+TILE_DIMx_d+2,(TILE_DIMx_d+2))
                  !ljj=mod(ljj+TILE_DIMy_d+2,(TILE_DIMy_d+2))
                  lkk=mod(lkk+TILE_DIMz_d+2,(TILE_DIMz_d+2))
#ifdef SECOND_ORDER
			      feq=(2.0_db*press - 3.0_db*v**2.0_db &
			       + 6.0_db*w + 6.0_db*(u &
			       + u**2.0_db + 3.0_db*u*w + w**2.0_db))/108.0_db
#else
!15
			      feq=(2.0_db*press + 3.0_db*(2.0_db*w*(1.0_db &
			       + w) - v**2.0_db*(1.0_db + 3.0_db*w*(1.0_db &
			       + w)) - u*(-2.0_db + 3.0_db*v**2.0_db)*(1.0_db &
			       + 3.0_db*w*(1.0_db + w)) - u**2.0_db*(-2.0_db &
			       + 3.0_db*v**2.0_db)*(1.0_db + 3.0_db*w*(1.0_db &
			       + w))))/108.0_db
!15
#endif
				  fneq1=3.0_db*pxx + 9.0_db*pxz &
				   - (3.0_db*pyy)/2.0_db + 3.0_db*pzz


				  F_discr=(forcex + forcez + 2.0_db*forcex*u &
				   + 3.0_db*forcez*u - forcey*v + 3.0_db*forcex*w &
				   + 2.0_db*forcez*w)/(18.0_db*rhophi_loc)
                  f3(lii,ljj,lkk)=feq + (1.0_db-omega_loc)*fneq1*p2 + 0.5_db*(F_discr)
!                  if(gi==0 .and. gj==2 .and. gk==15 .and. intblock==236)then
!                    write(*,*)'15 check',step,forcez
!                  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!16
                  lii=li-1
                  ljj=lj
                  lkk=lk-1
                  lii=mod(lii+TILE_DIMx_d+2,(TILE_DIMx_d+2))
                  !ljj=mod(ljj+TILE_DIMy_d+2,(TILE_DIMy_d+2))
                  lkk=mod(lkk+TILE_DIMz_d+2,(TILE_DIMz_d+2))
#ifdef SECOND_ORDER
			      feq=(2.0_db*press + 3.0_db*(2.0_db*u**2.0_db &
			       - v**2.0_db + 2.0_db*(-1.0_db + w)*w &
			       + u*(-2.0_db + 6.0_db*w)))/108.0_db
#else
!16
			      feq=(2.0_db*press + 3.0_db*(2.0_db*(-1.0_db &
			       + w)*w + v**2.0_db*(-1.0_db - 3.0_db*(-1.0_db &
			       + w)*w) + u*(-2.0_db + 3.0_db*v**2.0_db)*(1.0_db &
			       + 3.0_db*(-1.0_db + w)*w) - u**2.0_db*(-2.0_db &
			       + 3.0_db*v**2.0_db)*(1.0_db + 3.0_db*(-1.0_db &
			       + w)*w)))/108.0_db
!16
#endif
				  fneq1=3.0_db*pxx + 9.0_db*pxz &
				   - (3.0_db*pyy)/2.0_db + 3.0_db*pzz


				  F_discr=(forcex + forcez - 2.0_db*forcex*u &
				   - 3.0_db*forcez*u + forcey*v - 3.0_db*forcex*w &
				   - 2.0_db*forcez*w)/(-18.0_db*rhophi_loc)
                  f4(lii,ljj,lkk)=feq + (1.0_db-omega_loc)*fneq1*p2 + 0.5_db*(F_discr)
        
                  call syncthreads
                  
                  opress = opress + (f1(li,lj,lk)+f2(li,lj,lk)+f3(li,lj,lk)+f4(li,lj,lk))
                  ou = ou + f3(li,lj,lk)-f4(li,lj,lk)
                  ov = ov + f1(li,lj,lk)-f2(li,lj,lk)
                  ow = ow - f1(li,lj,lk)+f2(li,lj,lk)+f3(li,lj,lk)-f4(li,lj,lk)
                  opxx = opxx + f3(li,lj,lk)+f4(li,lj,lk)
                  opyy = opyy + f1(li,lj,lk) +f2(li,lj,lk)
                  opzz = opzz + f1(li,lj,lk) +f2(li,lj,lk)+f3(li,lj,lk)+f4(li,lj,lk)
                  opxz = f3(li,lj,lk)+f4(li,lj,lk)
                  opyz = opyz - f1(li,lj,lk) -f2(li,lj,lk)
!                  if(gi==1 .and. gj==2 .and. gk==16 .and. myblock==intblock)then
!                    write(*,*)'13',f1(li,lj,lk),step
!                    write(*,*)'14',f2(li,lj,lk),step
!                    write(*,*)'15',f3(li,lj,lk),step
!                    write(*,*)'16',f4(li,lj,lk),step
!                  endif
                  call syncthreads

!!!!!!!!!!!!!!!!!!!!!!!!!!17
                  lii=li-1
                  ljj=lj
                  lkk=lk+1
                  lii=mod(lii+TILE_DIMx_d+2,(TILE_DIMx_d+2))
                  !ljj=mod(ljj+TILE_DIMy_d+2,(TILE_DIMy_d+2))
                  lkk=mod(lkk+TILE_DIMz_d+2,(TILE_DIMz_d+2))
#ifdef SECOND_ORDER
			      feq=(2.0_db*press + 3.0_db*(2.0_db*u**2.0_db &
			       - v**2.0_db + 2.0_db*w*(1.0_db + w) &
			       - 2.0_db*u*(1.0_db + 3.0_db*w)))/108.0_db
#else
!17
			      feq=(2.0_db*press + 3.0_db*(2.0_db*w*(1.0_db &
			       + w) - v**2.0_db*(1.0_db + 3.0_db*w*(1.0_db &
			       + w)) + u*(-2.0_db + 3.0_db*v**2.0_db)*(1.0_db &
			       + 3.0_db*w*(1.0_db + w)) - u**2.0_db*(-2.0_db &
			       + 3.0_db*v**2.0_db)*(1.0_db + 3.0_db*w*(1.0_db &
			       + w))))/108.0_db
!17
#endif
				  fneq1=3.0_db*pxx - 9.0_db*pxz &
				   - (3.0_db*pyy)/2.0_db + 3.0_db*pzz


				  F_discr=(forcex - forcez - 2.0_db*forcex*u &
				   + 3.0_db*forcez*u + forcey*v + 3.0_db*forcex*w &
				   - 2.0_db*forcez*w)/(-18.0_db*rhophi_loc)
                  f1(lii,ljj,lkk)=feq + (1.0_db-omega_loc)*fneq1*p2 + 0.5_db*(F_discr)


!!!!!!!!!!!!!!!!!!!!!!!!!!18
                  lii=li+1
                  ljj=lj
                  lkk=lk-1
                  lii=mod(lii+TILE_DIMx_d+2,(TILE_DIMx_d+2))
                  !ljj=mod(ljj+TILE_DIMy_d+2,(TILE_DIMy_d+2))
                  lkk=mod(lkk+TILE_DIMz_d+2,(TILE_DIMz_d+2))
#ifdef SECOND_ORDER
			      feq=(2.0_db*press - 3.0_db*v**2.0_db &
			       - 6.0_db*w + 6.0_db*(u &
			       + u**2.0_db - 3.0_db*u*w + w**2.0_db))/108.0_db
#else
!18
			      feq=(2.0_db*press + 3.0_db*(2.0_db*(-1.0_db &
			       + w)*w + v**2.0_db*(-1.0_db &
			       - 3.0_db*(-1.0_db + w)*w) - u*(-2.0_db &
			       + 3.0_db*v**2.0_db)*(1.0_db + 3.0_db*(-1.0_db &
			       + w)*w) - u**2.0_db*(-2.0_db &
			       + 3.0_db*v**2.0_db)*(1.0_db + 3.0_db*(-1.0_db &
			       + w)*w)))/108.0_db
!18
#endif
				  fneq1=3.0_db*pxx - 9.0_db*pxz &
				   - (3.0_db*pyy)/2.0_db + 3.0_db*pzz


				  F_discr=(forcez + 3.0_db*forcez*u + forcey*v &
				   - 2.0_db*forcez*w + forcex*(-1.0_db - 2.0_db*u &
				   + 3.0_db*w))/(-18.0_db*rhophi_loc)
                  f2(lii,ljj,lkk)=feq + (1.0_db-omega_loc)*fneq1*p2 + 0.5_db*(F_discr)


!!!!!!!!!!!!!!!!!!!!!!!!!!19
                  lii=li+1
                  ljj=lj+1
                  lkk=lk+1
                  lii=mod(lii+TILE_DIMx_d+2,(TILE_DIMx_d+2))
                  ljj=mod(ljj+TILE_DIMy_d+2,(TILE_DIMy_d+2))
                  lkk=mod(lkk+TILE_DIMz_d+2,(TILE_DIMz_d+2))
#ifdef SECOND_ORDER
			      feq=(press + 3.0_db*(u + u**2.0_db &
			       + v + 3.0_db*u*v + v**2.0_db &
			       + w + 3.0_db*(u + v)*w &
			       + w**2.0_db))/216.0_db
#else
!19
			      feq=(press + 3.0_db*(v + v**2.0_db &
			       + w + 3.0_db*v*(1.0_db + v)*w + (1.0_db &
			       + 3.0_db*v*(1.0_db + v))*w**2.0_db + u*(1.0_db &
			       + 3.0_db*v*(1.0_db + v))*(1.0_db + 3.0_db*w*(1.0_db &
			       + w)) + u**2.0_db*(1.0_db + 3.0_db*v*(1.0_db &
			       + v))*(1.0_db + 3.0_db*w*(1.0_db + w))))/216.0_db
!19
#endif
				  fneq1=3.0_db*(pxx + 3.0_db*pxy &
				   + 3.0_db*pxz + pyy + 3.0_db*pyz + pzz)


				  F_discr=(forcez + 3.0_db*forcez*(u + v) &
				   + 2.0_db*forcez*w + forcey*(1.0_db + 3.0_db*u &
				   + 2.0_db*v + 3.0_db*w) + forcex*(1.0_db + 2.0_db*u &
				   + 3.0_db*v + 3.0_db*w))/(72.0_db*rhophi_loc)
                  f3(lii,ljj,lkk)=feq + (1.0_db-omega_loc)*fneq1*p3 + 0.5_db*(F_discr)


!!!!!!!!!!!!!!!!!!!!!!!!!!20
                  lii=li-1
                  ljj=lj-1
                  lkk=lk-1
                  lii=mod(lii+TILE_DIMx_d+2,(TILE_DIMx_d+2))
                  ljj=mod(ljj+TILE_DIMy_d+2,(TILE_DIMy_d+2))
                  lkk=mod(lkk+TILE_DIMz_d+2,(TILE_DIMz_d+2))
#ifdef SECOND_ORDER
			      feq=(press + 3.0_db*(u**2.0_db + (-1.0_db &
			       + v)*v + (-1.0_db + 3.0_db*v)*w &
			       + w**2.0_db + u*(-1.0_db + 3.0_db*v &
			       + 3.0_db*w)))/216.0_db
#else
!20
			      feq=(press + 3.0_db*((-1.0_db + v)*v &
			       + (-1.0_db - 3.0_db*(-1.0_db + v)*v)*w &
			       + (1.0_db + 3.0_db*(-1.0_db + v)*v)*w**2.0_db &
			       - u*(1.0_db + 3.0_db*(-1.0_db + v)*v)*(1.0_db &
			       + 3.0_db*(-1.0_db + w)*w) + u**2.0_db*(1.0_db &
			       + 3.0_db*(-1.0_db + v)*v)*(1.0_db + 3.0_db*(-1.0_db &
			       + w)*w)))/216.0_db
!20
#endif
				  fneq1=3.0_db*(pxx + 3.0_db*pxy &
				   + 3.0_db*pxz + pyy + 3.0_db*pyz + pzz)


				  F_discr=(forcez*(-1.0_db + 3.0_db*u + 3.0_db*v &
				   + 2.0_db*w) + forcey*(-1.0_db + 3.0_db*u &
				   + 2.0_db*v + 3.0_db*w) + forcex*(-1.0_db &
				   + 2.0_db*u + 3.0_db*v + 3.0_db*w))/(72.0_db*rhophi_loc)
                  f4(lii,ljj,lkk)=feq + (1.0_db-omega_loc)*fneq1*p3 + 0.5_db*(F_discr)
     
                  call syncthreads
                  
                  opress = opress + (f1(li,lj,lk)+f2(li,lj,lk)+f3(li,lj,lk)+f4(li,lj,lk))
                  ou = ou - f1(li,lj,lk) +f2(li,lj,lk)+f3(li,lj,lk)-f4(li,lj,lk)
                  ov = ov + f3(li,lj,lk)-f4(li,lj,lk)
                  ow = ow + f1(li,lj,lk) -f2(li,lj,lk)+f3(li,lj,lk)-f4(li,lj,lk)
                  opxx = opxx + f1(li,lj,lk)+f2(li,lj,lk)+f3(li,lj,lk)+f4(li,lj,lk)
                  opyy = opyy + f3(li,lj,lk)+f4(li,lj,lk)
                  opzz = opzz + f1(li,lj,lk)+f2(li,lj,lk)+f3(li,lj,lk)+f4(li,lj,lk)
                  opxy = opxy + f3(li,lj,lk)+f4(li,lj,lk)
                  opxz = opxz - f1(li,lj,lk)-f2(li,lj,lk)+f3(li,lj,lk)+f4(li,lj,lk)
                  opyz = opyz + f3(li,lj,lk)+f4(li,lj,lk)
!                  if(gi==1 .and. gj==2 .and. gk==16 .and. myblock==intblock)then
!                    write(*,*)'17',f1(li,lj,lk),step
!                    write(*,*)'18',f2(li,lj,lk),step
!                    write(*,*)'19',f3(li,lj,lk),step
!                    write(*,*)'20',f4(li,lj,lk),step
!                  endif
                  call syncthreads

!!!!!!!!!!!!!!!!!!!!!!!!!!21
                  lii=li+1
                  ljj=lj-1
                  lkk=lk+1
                  lii=mod(lii+TILE_DIMx_d+2,(TILE_DIMx_d+2))
                  ljj=mod(ljj+TILE_DIMy_d+2,(TILE_DIMy_d+2))
                  lkk=mod(lkk+TILE_DIMz_d+2,(TILE_DIMz_d+2))
#ifdef SECOND_ORDER
			      feq=(press + 3.0_db*(u**2.0_db + (-1.0_db &
			       + v)*v + w - 3.0_db*v*w &
			       + w**2.0_db + u*(1.0_db - 3.0_db*v &
			       + 3.0_db*w)))/216.0_db
#else
!21
			      feq=(press + 3.0_db*((-1.0_db + v)*v &
			       + w + 3.0_db*(-1.0_db + v)*v*w &
			       + (1.0_db + 3.0_db*(-1.0_db + v)*v)*w**2.0_db &
			       + u*(1.0_db + 3.0_db*(-1.0_db + v)*v)*(1.0_db &
			       + 3.0_db*w*(1.0_db + w)) + u**2.0_db*(1.0_db &
			       + 3.0_db*(-1.0_db + v)*v)*(1.0_db &
			       + 3.0_db*w*(1.0_db + w))))/216.0_db
!21
#endif 
				  fneq1=3.0_db*(pxx - 3.0_db*pxy &
				   + 3.0_db*pxz + pyy - 3.0_db*pyz + pzz)


				  F_discr=(forcez*(1.0_db + 3.0_db*u - 3.0_db*v &
				   + 2.0_db*w) + forcex*(1.0_db + 2.0_db*u &
				   - 3.0_db*v + 3.0_db*w) - forcey*(1.0_db &
				   + 3.0_db*u - 2.0_db*v + 3.0_db*w))/(72.0_db*rhophi_loc)
                  f1(lii,ljj,lkk)=feq + (1.0_db-omega_loc)*fneq1*p3 + 0.5_db*(F_discr)



!!!!!!!!!!!!!!!!!!!!!!!!!!22
                  lii=li-1
                  ljj=lj+1
                  lkk=lk-1
                  lii=mod(lii+TILE_DIMx_d+2,(TILE_DIMx_d+2))
                  ljj=mod(ljj+TILE_DIMy_d+2,(TILE_DIMy_d+2))
                  lkk=mod(lkk+TILE_DIMz_d+2,(TILE_DIMz_d+2))
#ifdef SECOND_ORDER
			      feq=(press + 3.0_db*(u**2.0_db + v &
			       + v**2.0_db - 3.0_db*v*w + (-1.0_db &
			       + w)*w + u*(-1.0_db - 3.0_db*v &
			       + 3.0_db*w)))/216.0_db
#else
!22
			      feq=(press + 3.0_db*(v + v**2.0_db &
			       - w - 3.0_db*v*(1.0_db + v)*w + (1.0_db &
			       + 3.0_db*v*(1.0_db + v))*w**2.0_db - u*(1.0_db &
			       + 3.0_db*v*(1.0_db + v))*(1.0_db + 3.0_db*(-1.0_db &
			       + w)*w) + u**2.0_db*(1.0_db + 3.0_db*v*(1.0_db &
			       + v))*(1.0_db + 3.0_db*(-1.0_db + w)*w)))/216.0_db
!22
#endif
				  fneq1=3.0_db*(pxx - 3.0_db*pxy &
				   + 3.0_db*pxz + pyy - 3.0_db*pyz + pzz)


				  F_discr=(forcey*(1.0_db - 3.0_db*u + 2.0_db*v &
				   - 3.0_db*w) + forcez*(-1.0_db + 3.0_db*u &
				   - 3.0_db*v + 2.0_db*w) + forcex*(-1.0_db &
				   + 2.0_db*u - 3.0_db*v + 3.0_db*w))/(72.0_db*rhophi_loc)
                  f2(lii,ljj,lkk)=feq + (1.0_db-omega_loc)*fneq1*p3 + 0.5_db*(F_discr)




!!!!!!!!!!!!!!!!!!!!!!!!!!23
                  lii=li-1
                  ljj=lj-1
                  lkk=lk+1
                  lii=mod(lii+TILE_DIMx_d+2,(TILE_DIMx_d+2))
                  ljj=mod(ljj+TILE_DIMy_d+2,(TILE_DIMy_d+2))
                  lkk=mod(lkk+TILE_DIMz_d+2,(TILE_DIMz_d+2))
#ifdef SECOND_ORDER
			      feq=(press + 3.0_db*(-u + u**2.0_db &
			       - v + 3.0_db*u*v + v**2.0_db + w &
			       - 3.0_db*(u + v)*w + w**2.0_db))/216.0_db
#else
!23
			      feq=(press + 3.0_db*((-1.0_db + v)*v &
			       + w + 3.0_db*(-1.0_db + v)*v*w + (1.0_db &
			       + 3.0_db*(-1.0_db + v)*v)*w**2.0_db - u*(1.0_db &
			       + 3.0_db*(-1.0_db + v)*v)*(1.0_db + 3.0_db*w*(1.0_db &
			       + w)) + u**2.0_db*(1.0_db + 3.0_db*(-1.0_db + v)*v)*(1.0_db &
			       + 3.0_db*w*(1.0_db + w))))/216.0_db
!23
#endif
				  fneq1=3.0_db*(pxx + 3.0_db*pxy &
				   - 3.0_db*pxz + pyy - 3.0_db*pyz + pzz)


				  F_discr=(forcez - 3.0_db*forcez*(u + v) + forcey*(-1.0_db &
				   + 3.0_db*u + 2.0_db*v - 3.0_db*w) &
				   + forcex*(-1.0_db + 2.0_db*u + 3.0_db*v &
				   - 3.0_db*w) + 2.0_db*forcez*w)/(72.0_db*rhophi_loc)
                  f3(lii,ljj,lkk)=feq + (1.0_db-omega_loc)*fneq1*p3 + 0.5_db*(F_discr)


!!!!!!!!!!!!!!!!!!!!!!!!!!24
                  lii=li+1
                  ljj=lj+1
                  lkk=lk-1
                  lii=mod(lii+TILE_DIMx_d+2,(TILE_DIMx_d+2))
                  ljj=mod(ljj+TILE_DIMy_d+2,(TILE_DIMy_d+2))
                  lkk=mod(lkk+TILE_DIMz_d+2,(TILE_DIMz_d+2))
#ifdef SECOND_ORDER
			      feq=(press + 3.0_db*(u + u**2.0_db &
			       + v + 3.0_db*u*v + v**2.0_db - w &
			       - 3.0_db*(u + v)*w + w**2.0_db))/216.0_db
#else
!24
			      feq=(press + 3.0_db*(v + v**2.0_db &
			       - w - 3.0_db*v*(1.0_db + v)*w + (1.0_db &
			       + 3.0_db*v*(1.0_db + v))*w**2.0_db + u*(1.0_db &
			       + 3.0_db*v*(1.0_db + v))*(1.0_db + 3.0_db*(-1.0_db &
			       + w)*w) + u**2.0_db*(1.0_db + 3.0_db*v*(1.0_db &
			       + v))*(1.0_db + 3.0_db*(-1.0_db + w)*w)))/216.0_db
!24
#endif
				  fneq1=3.0_db*(pxx + 3.0_db*pxy &
				   - 3.0_db*pxz + pyy - 3.0_db*pyz + pzz)


				  F_discr=(-forcez - 3.0_db*forcez*(u + v) + forcey*(1.0_db &
				   + 3.0_db*u + 2.0_db*v - 3.0_db*w) + forcex*(1.0_db &
				   + 2.0_db*u + 3.0_db*v - 3.0_db*w) &
				   + 2.0_db*forcez*w)/(72.0_db*rhophi_loc)
                  f4(lii,ljj,lkk)=feq + (1.0_db-omega_loc)*fneq1*p3 + 0.5_db*(F_discr)
           
                  call syncthreads
                  
                  opress = opress + (f1(li,lj,lk)+f2(li,lj,lk)+f3(li,lj,lk)+f4(li,lj,lk))
                  ou = ou +f1(li,lj,lk)-f2(li,lj,lk)-f3(li,lj,lk)+f4(li,lj,lk)
                  ov = ov -f1(li,lj,lk)+f2(li,lj,lk)-f3(li,lj,lk)+f4(li,lj,lk)
                  ow = ow +f1(li,lj,lk)-f2(li,lj,lk)+f3(li,lj,lk)-f4(li,lj,lk)
                  opxx = opxx + f1(li,lj,lk)+f2(li,lj,lk)+f3(li,lj,lk)+f4(li,lj,lk)
                  opyy = opyy + f1(li,lj,lk)+f2(li,lj,lk)+f3(li,lj,lk)+f4(li,lj,lk)
                  opzz = opzz + f1(li,lj,lk)+f2(li,lj,lk)+f3(li,lj,lk)+f4(li,lj,lk)
                  opxy = opxy - f1(li,lj,lk)-f2(li,lj,lk)+f3(li,lj,lk)+f4(li,lj,lk)
                  opxz = opxz + f1(li,lj,lk)+f2(li,lj,lk)-f3(li,lj,lk)-f4(li,lj,lk)
                  opyz = opyz - f1(li,lj,lk)-f2(li,lj,lk)-f3(li,lj,lk)-f4(li,lj,lk)
!                  if(gi==1 .and. gj==2 .and. gk==16 .and. myblock==intblock)then
!                    write(*,*)'21',f1(li,lj,lk),step
!                    write(*,*)'22',f2(li,lj,lk),step
!                    write(*,*)'23',f3(li,lj,lk),step
!                    write(*,*)'24',f4(li,lj,lk),step
!                  endif
                  call syncthreads

!!!!!!!!!!!!!!!!!!!!!!!!!!25
                  lii=li+1
                  ljj=lj-1
                  lkk=lk-1
                  lii=mod(lii+TILE_DIMx_d+2,(TILE_DIMx_d+2))
                  ljj=mod(ljj+TILE_DIMy_d+2,(TILE_DIMy_d+2))
                  lkk=mod(lkk+TILE_DIMz_d+2,(TILE_DIMz_d+2))
#ifdef SECOND_ORDER
			      feq=(press + 3.0_db*(u**2.0_db + (-1.0_db &
			       + v)*v + u*(1.0_db - 3.0_db*v &
			       - 3.0_db*w) + (-1.0_db + 3.0_db*v)*w &
			       + w**2.0_db))/216.0_db
#else
!25
			      feq=(press + 3.0_db*((-1.0_db + v)*v &
			       + (-1.0_db - 3.0_db*(-1.0_db + v)*v)*w + (1.0_db &
			       + 3.0_db*(-1.0_db + v)*v)*w**2.0_db + u*(1.0_db &
			       + 3.0_db*(-1.0_db + v)*v)*(1.0_db + 3.0_db*(-1.0_db &
			       + w)*w) + u**2.0_db*(1.0_db + 3.0_db*(-1.0_db &
			       + v)*v)*(1.0_db + 3.0_db*(-1.0_db + w)*w)))/216.0_db
!25
#endif
				  fneq1=3.0_db*(pxx - 3.0_db*pxy - 3.0_db*pxz &
				   + pyy + 3.0_db*pyz + pzz)


				  F_discr=(forcex*(1.0_db + 2.0_db*u - 3.0_db*v &
				   - 3.0_db*w) + forcez*(-1.0_db - 3.0_db*u &
				   + 3.0_db*v + 2.0_db*w) + forcey*(-1.0_db &
				   - 3.0_db*u + 2.0_db*v + 3.0_db*w))/(72.0_db*rhophi_loc)
                  f1(lii,ljj,lkk)=feq + (1.0_db-omega_loc)*fneq1*p3 + 0.5_db*(F_discr)


!!!!!!!!!!!!!!!!!!!!!!!!!!26
                  lii=li-1
                  ljj=lj+1
                  lkk=lk+1
                  lii=mod(lii+TILE_DIMx_d+2,(TILE_DIMx_d+2))
                  ljj=mod(ljj+TILE_DIMy_d+2,(TILE_DIMy_d+2))
                  lkk=mod(lkk+TILE_DIMz_d+2,(TILE_DIMz_d+2))
#ifdef SECOND_ORDER
			      feq=(press + 3.0_db*(u**2.0_db + v &
			       + v**2.0_db + w + 3.0_db*v*w &
			       + w**2.0_db - u*(1.0_db + 3.0_db*v &
			       + 3.0_db*w)))/216.0_db
#else
!26
			      feq=(press + 3.0_db*(v + v**2.0_db &
			       + w + 3.0_db*v*(1.0_db + v)*w + (1.0_db &
			       + 3.0_db*v*(1.0_db + v))*w**2.0_db - u*(1.0_db &
			       + 3.0_db*v*(1.0_db + v))*(1.0_db + 3.0_db*w*(1.0_db &
			       + w)) + u**2.0_db*(1.0_db + 3.0_db*v*(1.0_db &
			       + v))*(1.0_db + 3.0_db*w*(1.0_db + w))))/216.0_db
!26
#endif
				  fneq1=3.0_db*(pxx - 3.0_db*pxy &
				   - 3.0_db*pxz + pyy + 3.0_db*pyz + pzz)


				  F_discr=(forcex*(-1.0_db + 2.0_db*u - 3.0_db*v &
				   - 3.0_db*w) + forcez*(1.0_db - 3.0_db*u &
				   + 3.0_db*v + 2.0_db*w) + forcey*(1.0_db &
				   - 3.0_db*u + 2.0_db*v + 3.0_db*w))/(72.0_db*rhophi_loc)
                  f2(lii,ljj,lkk)=feq + (1.0_db-omega_loc)*fneq1*p3 + 0.5_db*(F_discr)
                  
                  call syncthreads
                           
                  opress = opress + (f1(li,lj,lk)+f2(li,lj,lk))
                  ou = ou + f1(li,lj,lk)-f2(li,lj,lk)
                  ov = ov - f1(li,lj,lk)+f2(li,lj,lk)
                  ow = ow - f1(li,lj,lk)+f2(li,lj,lk)
                  opxx = opxx + f1(li,lj,lk)+f2(li,lj,lk)
                  opyy = opyy + f1(li,lj,lk)+f2(li,lj,lk)
                  opzz = opzz + f1(li,lj,lk)+f2(li,lj,lk)
                  opxy = opxy - f1(li,lj,lk)-f2(li,lj,lk)
                  opxz = opxz - f1(li,lj,lk)-f2(li,lj,lk) 
                  opyz = opyz + f1(li,lj,lk)+f2(li,lj,lk)
!                  if(gi==1 .and. gj==2 .and. gk==16 .and. myblock==intblock)then
!                    write(*,*)'25',f1(li,lj,lk),step
!                    write(*,*)'26',f2(li,lj,lk),step
!                  endif
                  !internal-node block is the index of the block of internal nodes without the surrounding halo
	              intblock=blockIdx%x+blockIdx%y*nxblock_d+blockIdx%z*nxyblock_d+1 !internal-node block
                  
                  !If my block index does not match the index of the internal-node block (lii), it means my thread is on the outer. I must exit
	              if(myblock .ne. intblock)return
	                 
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

   endsubroutine fused_LB_kernel   
   

endmodule
