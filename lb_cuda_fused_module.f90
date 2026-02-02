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
   
   
   attributes(global) subroutine fused_LB_kernel2(step,iprobe,jprobe,kprobe,flip,flop,nx,ny,nz,coords,isfluid &  
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
       hfields_in,hfields_out &
#ifdef TWOCOMPONENT 
       ,auxfields_s,locauxfields_s &
#endif   
       ,forces_s)

      implicit none
      
      integer :: step,iprobe,jprobe,kprobe,flip,flop,nx,ny,nz
      
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
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields,nblocks_d) :: phifields_s
#ifdef MONOD
      real(kind=db) :: mu_max,Ks
#endif
#endif           
      real(kind=db) :: visc1,omega,fx,fy,fz
      
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields,nblocks_d) :: hfields_in,hfields_out
#ifdef TWOCOMPONENT 
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nauxfields,nblocks_d) :: auxfields_s
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nlocauxfields,nblocks_d) :: locauxfields_s
#endif  
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nforces,nblocks_d) :: forces_s
      
      real(kind=db), shared :: f_front(0:TILE_DIMx+1,0:TILE_DIMy+1,0:TILE_DIMz+1)
      real(kind=db), shared :: f_rear(0:TILE_DIMx+1,0:TILE_DIMy+1,0:TILE_DIMz+1)
      
      real(kind=db) :: press,u,v,w,pxx,pyy,pzz,pxy,pxz,pyz
      real(kind=db) :: opress,ou,ov,ow,opxx,opyy,opzz,opxy,opxz,opyz,feq,fneq1,f_discr
      real(kind=db) :: mytemp,forcex,forcey,forcez,invrhophi_loc,uu,udotc
      
      real(kind=db) :: omega_loc,phi_loc,visc_loc
#ifdef SMAGORINSKI
	  real(kind=db) :: QQ
#endif

      integer :: i,j,k,myblock,l,intblock
      integer :: ii,jj,kk
      integer :: li,lj,lk
      integer :: lii,ljj,lkk
      integer :: xblock,yblock,zblock
      !integer :: gi,gj,gk

      
      li = threadIdx%x-1
      lj = threadIdx%y-1
      lk = threadIdx%z-1
      
      i = (blockIdx%x-1) * TILE_DIMx + li
      j = (blockIdx%y-1) * TILE_DIMy + lj
      k = (blockIdx%z-1) * TILE_DIMz + lk
      
!      gi=nx*coords(1)+i
!      gj=ny*coords(2)+j
!      gk=nz*coords(3)+k
      
      xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
	  yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
	  zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
      
      myblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
      ii=i-xblock*TILE_DIMx+2*TILE_DIMx
      jj=j-yblock*TILE_DIMy+2*TILE_DIMy
      kk=k-zblock*TILE_DIMz+2*TILE_DIMz

      !i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      !j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      !k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z
      
!      gi=nx*coords(1)+i
!      gj=ny*coords(2)+j
!      gk=nz*coords(3)+k
      
      intblock=blockIdx%x+blockIdx%y*nxblock_d+blockIdx%z*nxyblock_d+1 !internal-node block


               !if (abs(isfluid(i,j,k)) /= 1) return
       
     
#ifdef TWOCOMPONENT
                  phi_loc=real(phifields_s(ii,jj,kk,1,myblock),kind=db)
#endif                  
#ifdef DENSRATIO
                  
                  invrhophi_loc = 1.0_db/(rho_r*phi_loc+(ONE-phi_loc)*rho_b) 
#else
                  invrhophi_loc = 1.0_db !press_loc
#endif	

				  forcex=real(forces_s(ii,jj,kk,1,myblock),kind=db)*invrhophi_loc
				  forcey=real(forces_s(ii,jj,kk,2,myblock),kind=db)*invrhophi_loc
				  forcez=real(forces_s(ii,jj,kk,3,myblock),kind=db)*invrhophi_loc
                  
                  press=real(hfields_in(ii,jj,kk,1,myblock),kind=db)
                  u=real(hfields_in(ii,jj,kk,2,myblock),kind=db)
                  v=real(hfields_in(ii,jj,kk,3,myblock),kind=db)
                  w=real(hfields_in(ii,jj,kk,4,myblock),kind=db)
                  pxx=real(hfields_in(ii,jj,kk,5,myblock),kind=db)
                  pyy=real(hfields_in(ii,jj,kk,6,myblock),kind=db)
                  pzz=real(hfields_in(ii,jj,kk,7,myblock),kind=db)
                  pxy=real(hfields_in(ii,jj,kk,8,myblock),kind=db)
                  pxz=real(hfields_in(ii,jj,kk,9,myblock),kind=db)
                  pyz=real(hfields_in(ii,jj,kk,10,myblock),kind=db)
                 
                  
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
                  

!				  uu=HALF*(u*u+v*v+w*w)*invcssq
				  
!                  do lii=1,nlinks
!                     udotc=(u*dex(lii) + v*dey(lii)+ w*dez(lii))*invcssq
!		     feq=p(lii)*(press + (udotc+0.5_db*udotc*udotc - uu))
!                     !fneq1=p(l)*(HALF/(cssq*cssq))*( (dex(lii)*dex(lii)-cssq)*pxx &
!		     ! + (dey(lii)*dey(lii)-cssq)*pyy + (dez(lii)*dez(lii)-cssq)*pzz &
!	             ! + TWO*(dex(lii)*dey(lii))*pxy + TWO*(dex(lii)*dez(lii))*pxz &
!		     ! + TWO*(dey(lii)*dez(lii))*pyz)
!                     fpost=feq!+fneq1
!                     pxx=pxx - fpost*(dex(lii)*dex(lii))
!                     pyy=pyy - fpost*(dey(lii)*dey(lii))
!                     pzz=pzz - fpost*(dez(lii)*dez(lii))
!                     pxy=pxy - fpost*(dex(lii)*dey(lii))
!                     pxz=pxz - fpost*(dex(lii)*dez(lii))
!                     pyz=pyz - fpost*(dey(lii)*dez(lii))
!                  enddo

                  pxx=pxx - cssq*press - u*u 
                  pyy=pyy - cssq*press - v*v 
                  pzz=pzz - cssq*press - w*w 
                  pxy=pxy - u*v
                  pxz=pxz - u*w
                  pyz=pyz - v*w

#ifdef TWOCOMPONENT
                  !visc_loc it is used to store the local viscosity
                  visc_loc=(rho_r*visc1*phi_loc+(ONE-phi_loc)*visc2*rho_b)*invrhophi_loc
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
                   sqrt((3.0*visc_loc)**2.0 + 0.053*18.0*sqrt(2.0*QQ)*invrhophi_loc)) !it is tau
                  omega_loc=1.0_db/omega_loc !it is omega

#else
#ifdef TWOCOMPONENT
                  omega_loc=(visc_loc*invcssq + 0.5_db) !it is tau   !visc_loc it is used to store the local viscosity
                  omega_loc=1.0_db/omega_loc !it is omega
#else
                  omega_loc=omega
#endif
#endif      
      
                  !opress=ZERO
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
                  uu=HALF*(u*u+v*v+w*w)*invcssq

			      feq=p(0)*(press - uu)
				  fneq1=p(0)*(HALF*invcssq)*(-pxx-pyy-pzz)
#if defined(HIGHORDER) && (LATTICE == 27)
                  feq=feq + &
                    Minv_mrt(17,0)*u*u*v*v + &
                    Minv_mrt(18,0)*u*u*w*w + &
                    Minv_mrt(19,0)*v*v*w*w + &
                    Minv_mrt(26,0)*u*u*v*v*w*w
                  fneq1=fneq1 + &
                    Minv_mrt(17,0)*(pyy*u**TWO + TWO*v*(THREE*pxy*u + pxx*v)) + &
                    Minv_mrt(18,0)*(pzz*u**TWO + TWO*w*(THREE*pxz*u + pxx*w)) + &
                    Minv_mrt(19,0)*(pzz*v**TWO + TWO*w*(THREE*pyz*v + pyy*w)) + &
                    Minv_mrt(26,0)*(pzz*u**TWO*v**TWO + &
                    w*(SIX*u*v*(pyz*u + pxz*v) + THREE*pyy*u**TWO*w + TWO*v*(SEVEN*pxy*u + TWO*pxx*v)*w))
#endif
				  F_discr = p(0)*(- u*forcex - v*forcey - w*forcez)*invcssq
                  
                  opress=feq + (1.0_db-omega_loc)*fneq1 + HALF*(F_discr)
                  
                  do l=1,nlinks,2
                     udotc=(u*dex(l) + v*dey(l)+ w*dez(l))*invcssq
		             feq=p(l)*(press + (udotc+0.5_db*udotc*udotc - uu))
                     fneq1=p(l)*(HALF/(cssq*cssq))*( (dex(l)*dex(l)-cssq)*pxx &
		              + (dey(l)*dey(l)-cssq)*pyy + (dez(l)*dez(l)-cssq)*pzz &
	                  + TWO*(dex(l)*dey(l))*pxy + TWO*(dex(l)*dez(l))*pxz &
		              + TWO*(dey(l)*dez(l))*pyz)
#if defined(HIGHORDER) && (LATTICE == 27)
                     feq=feq+Minv_mrt(10,l)*u*u*v + &
                             Minv_mrt(11,l)*u*u*w + &
                             Minv_mrt(12,l)*u*v*v + &
                             Minv_mrt(13,l)*u*w*w + &
                             Minv_mrt(14,l)*v*w*w + &
                             Minv_mrt(15,l)*v*v*w + &
                             Minv_mrt(16,l)*u*v*w + &
                             Minv_mrt(17,l)*u*u*v*v + &
                             Minv_mrt(18,l)*u*u*w*w + &
                             Minv_mrt(19,l)*v*v*w*w + &
                             Minv_mrt(20,l)*u*v*w*w + &
                             Minv_mrt(21,l)*u*v*v*w + &
                             Minv_mrt(22,l)*u*u*v*w + &
                             Minv_mrt(23,l)*u*u*v*w*w + &
                             Minv_mrt(24,l)*u*u*v*v*w + &
                             Minv_mrt(25,l)*u*v*v*w*w + &
                             Minv_mrt(26,l)*u*u*v*v*w*w
                     fneq1=fneq1+Minv_mrt(10,l)*(TWO*pxy*u + pxx*v) + &
                             Minv_mrt(11,l)*(TWO*pxz*u + pxx*w) + &
                             Minv_mrt(12,l)*(pyy*u + TWO*pxy*v) + &
                             Minv_mrt(13,l)*(pzz*u + TWO*pxz*w)+ &
                             Minv_mrt(14,l)*(pzz*v + TWO*pyz*w) + &
                             Minv_mrt(15,l)*(TWO*pyz*v + pyy*w) + &
                             Minv_mrt(16,l)*(pyz*u + pxz*v + pxy*w) + &
                             Minv_mrt(17,l)*(pyy*u**TWO + TWO*v*(THREE*pxy*u + pxx*v)) + &
                             Minv_mrt(18,l)*(pzz*u**TWO + TWO*w*(THREE*pxz*u + pxx*w)) + &
                             Minv_mrt(19,l)*(pzz*v**TWO + TWO*w*(THREE*pyz*v + pyy*w)) + &
                             Minv_mrt(20,l)*(pzz*u*v + w*(THREE*pyz*u + THREE*pxz*v + TWO*pxy*w)) + &
                             Minv_mrt(21,l)*(TWO*pyz*u*v + pxz*v**TWO + TWO*pyy*u*w + FOUR*pxy*v*w) + &
                             Minv_mrt(22,l)*(pyz*u**TWO + TWO*pxz*u*v + FOUR*pxy*u*w + TWO*pxx*v*w) + &
                             Minv_mrt(23,l)*(pzz*u**TWO*v + THREE*w*(pyz*u**TWO + TWO*pxz*u*v + TWO*pxy*u*w + pxx*v*w)) + &
                             Minv_mrt(24,l)*(TWO*u*v*(pyz*u + pxz*v) + TWO*pyy*u**TWO*w + v*(TEN*pxy*u + THREE*pxx*v)*w) + &
                             Minv_mrt(25,l)*(pzz*u*v**TWO + THREE*w*(TWO*pyz*u*v + pxz*v**TWO + pyy*u*w + TWO*pxy*v*w)) + &
                             Minv_mrt(26,l)*(pzz*u**TWO*v**TWO + &
                              w*(SIX*u*v*(pyz*u + pxz*v) + THREE*pyy*u**TWO*w + TWO*v*(SEVEN*pxy*u + TWO*pxx*v)*w))
#endif
		             F_discr = p(l)*(((dex(l) - u) + udotc * dex(l))*forcex &
                      + ((dey(l) - v) + udotc * dey(l))*forcey &
                      + ((dez(l) - w) + udotc * dez(l))*forcez)*invcssq
                     lii=li+ex(l)
                     ljj=lj+ey(l)
                     lkk=lk+ez(l)
                     lii=mod(lii+TILE_DIMx+2,(TILE_DIMx+2))
                     ljj=mod(ljj+TILE_DIMy+2,(TILE_DIMy+2))
                     lkk=mod(lkk+TILE_DIMz+2,(TILE_DIMz+2)) 
		             f_front(lii,ljj,lkk)=feq + (ONE-omega_loc)*fneq1 + HALF*F_discr
		             
                     udotc=(u*dex(l+1) + v*dey(l+1)+ w*dez(l+1))*invcssq
		             feq=p(l+1)*(press + (udotc+0.5_db*udotc*udotc - uu))
                     fneq1=p(l+1)*(HALF/(cssq*cssq))*( (dex(l+1)*dex(l+1)-cssq)*pxx &
		              + (dey(l+1)*dey(l+1)-cssq)*pyy + (dez(l+1)*dez(l+1)-cssq)*pzz &
	                  + TWO*(dex(l+1)*dey(l+1))*pxy + TWO*(dex(l+1)*dez(l+1))*pxz &
		              + TWO*(dey(l+1)*dez(l+1))*pyz)
#if defined(HIGHORDER) && (LATTICE == 27)
                     feq=feq+Minv_mrt(10,l+1)*u*u*v + &
                             Minv_mrt(11,l+1)*u*u*w + &
                             Minv_mrt(12,l+1)*u*v*v + &
                             Minv_mrt(13,l+1)*u*w*w + &
                             Minv_mrt(14,l+1)*v*w*w + &
                             Minv_mrt(15,l+1)*v*v*w + &
                             Minv_mrt(16,l+1)*u*v*w + &
                             Minv_mrt(17,l+1)*u*u*v*v + &
                             Minv_mrt(18,l+1)*u*u*w*w + &
                             Minv_mrt(19,l+1)*v*v*w*w + &
                             Minv_mrt(20,l+1)*u*v*w*w + &
                             Minv_mrt(21,l+1)*u*v*v*w + &
                             Minv_mrt(22,l+1)*u*u*v*w + &
                             Minv_mrt(23,l+1)*u*u*v*w*w + &
                             Minv_mrt(24,l+1)*u*u*v*v*w + &
                             Minv_mrt(25,l+1)*u*v*v*w*w + &
                             Minv_mrt(26,l+1)*u*u*v*v*w*w
                     fneq1=fneq1+Minv_mrt(10,l+1)*(TWO*pxy*u + pxx*v) + &
                             Minv_mrt(11,l+1)*(TWO*pxz*u + pxx*w) + &
                             Minv_mrt(12,l+1)*(pyy*u + TWO*pxy*v) + &
                             Minv_mrt(13,l+1)*(pzz*u + TWO*pxz*w)+ &
                             Minv_mrt(14,l+1)*(pzz*v + TWO*pyz*w) + &
                             Minv_mrt(15,l+1)*(TWO*pyz*v + pyy*w) + &
                             Minv_mrt(16,l+1)*(pyz*u + pxz*v + pxy*w) + &
                             Minv_mrt(17,l+1)*(pyy*u**TWO + TWO*v*(THREE*pxy*u + pxx*v)) + &
                             Minv_mrt(18,l+1)*(pzz*u**TWO + TWO*w*(THREE*pxz*u + pxx*w)) + &
                             Minv_mrt(19,l+1)*(pzz*v**TWO + TWO*w*(THREE*pyz*v + pyy*w)) + &
                             Minv_mrt(20,l+1)*(pzz*u*v + w*(THREE*pyz*u + THREE*pxz*v + TWO*pxy*w)) + &
                             Minv_mrt(21,l+1)*(TWO*pyz*u*v + pxz*v**TWO + TWO*pyy*u*w + FOUR*pxy*v*w) + &
                             Minv_mrt(22,l+1)*(pyz*u**TWO + TWO*pxz*u*v + FOUR*pxy*u*w + TWO*pxx*v*w) + &
                             Minv_mrt(23,l+1)*(pzz*u**TWO*v + THREE*w*(pyz*u**TWO + TWO*pxz*u*v + TWO*pxy*u*w + pxx*v*w)) + &
                             Minv_mrt(24,l+1)*(TWO*u*v*(pyz*u + pxz*v) + TWO*pyy*u**TWO*w + v*(TEN*pxy*u + THREE*pxx*v)*w) + &
                             Minv_mrt(25,l+1)*(pzz*u*v**TWO + THREE*w*(TWO*pyz*u*v + pxz*v**TWO + pyy*u*w + TWO*pxy*v*w)) + &
                             Minv_mrt(26,l+1)*(pzz*u**TWO*v**TWO + &
                              w*(SIX*u*v*(pyz*u + pxz*v) + THREE*pyy*u**TWO*w + TWO*v*(SEVEN*pxy*u + TWO*pxx*v)*w))
#endif
		             F_discr = p(l+1)*(((dex(l+1) - u) + udotc * dex(l+1))*forcex &
                      + ((dey(l+1) - v) + udotc * dey(l+1))*forcey &
                      + ((dez(l+1) - w) + udotc * dez(l+1))*forcez)*invcssq
                     lii=li+ex(l+1)
                     ljj=lj+ey(l+1)
                     lkk=lk+ez(l+1)
                     lii=mod(lii+TILE_DIMx+2,(TILE_DIMx+2))
                     ljj=mod(ljj+TILE_DIMy+2,(TILE_DIMy+2))
                     lkk=mod(lkk+TILE_DIMz+2,(TILE_DIMz+2)) 
		             f_rear(lii,ljj,lkk)=feq + (ONE-omega_loc)*fneq1 + HALF*F_discr
		             
		             call syncthreads
		             opress=opress + f_front(li,lj,lk)
		             ou=ou + f_front(li,lj,lk)*dex(l)
                     ov=ov + f_front(li,lj,lk)*dey(l)
                     ow=ow + f_front(li,lj,lk)*dez(l)
                     opxx=opxx + f_front(li,lj,lk)*dex(l)*dex(l)
                     opyy=opyy + f_front(li,lj,lk)*dey(l)*dey(l)
                     opzz=opzz + f_front(li,lj,lk)*dez(l)*dez(l)
                     opxy=opxy + f_front(li,lj,lk)*dex(l)*dey(l)
                     opxz=opxz + f_front(li,lj,lk)*dex(l)*dez(l)
                     opyz=opyz + f_front(li,lj,lk)*dey(l)*dez(l)
                     
                     opress=opress + f_rear(li,lj,lk)
		             ou=ou + f_rear(li,lj,lk)*dex(l+1)
                     ov=ov + f_rear(li,lj,lk)*dey(l+1)
                     ow=ow + f_rear(li,lj,lk)*dez(l+1)
                     opxx=opxx + f_rear(li,lj,lk)*dex(l+1)*dex(l+1)
                     opyy=opyy + f_rear(li,lj,lk)*dey(l+1)*dey(l+1)
                     opzz=opzz + f_rear(li,lj,lk)*dez(l+1)*dez(l+1)
                     opxy=opxy + f_rear(li,lj,lk)*dex(l+1)*dey(l+1)
                     opxz=opxz + f_rear(li,lj,lk)*dex(l+1)*dez(l+1)
                     opyz=opyz + f_rear(li,lj,lk)*dey(l+1)*dez(l+1)
                     call syncthreads
                     
                  enddo
                  
                  
                  !If my block index does not match the index of the internal-node block (lii), it means my thread is on the outer. I must exit
	              if(myblock .ne. (blockIdx%x+blockIdx%y*nxblock_d+blockIdx%z*nxyblock_d+1) )return

	              hfields_out(ii,jj,kk,1,myblock)=real(opress,kind=strdb)
                  hfields_out(ii,jj,kk,2,myblock)=real(ou,kind=strdb)
                  hfields_out(ii,jj,kk,3,myblock)=real(ov,kind=strdb)
                  hfields_out(ii,jj,kk,4,myblock)=real(ow,kind=strdb)
                  hfields_out(ii,jj,kk,5,myblock)=real(opxx,kind=strdb)
                  hfields_out(ii,jj,kk,6,myblock)=real(opyy,kind=strdb)
                  hfields_out(ii,jj,kk,7,myblock)=real(opzz,kind=strdb)
                  hfields_out(ii,jj,kk,8,myblock)=real(opxy,kind=strdb)
                  hfields_out(ii,jj,kk,9,myblock)=real(opxz,kind=strdb)
                  hfields_out(ii,jj,kk,10,myblock)=real(opyz,kind=strdb)
              
    endsubroutine fused_LB_kernel2     

      attributes(global) subroutine fused_LB_kernel1(step,iprobe,jprobe,kprobe,flip,flop,nx,ny,nz,coords,isfluid &  
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
       hfields_in,hfields_out &
#ifdef TWOCOMPONENT 
       ,auxfields_s,locauxfields_s &
#endif   
       ,forces_s)

      implicit none
      
      integer :: step,iprobe,jprobe,kprobe,flip,flop,nx,ny,nz
      
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
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields,nblocks_d) :: phifields_s
#ifdef MONOD
      real(kind=db) :: mu_max,Ks
#endif
#endif           
      real(kind=db) :: visc1,omega,fx,fy,fz
      
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields,nblocks_d) :: hfields_in,hfields_out
#ifdef TWOCOMPONENT 
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nauxfields,nblocks_d) :: auxfields_s
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nlocauxfields,nblocks_d) :: locauxfields_s
#endif  
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nforces,nblocks_d) :: forces_s
      
      real(kind=db), shared :: ftemp(0:TILE_DIMx+1,0:TILE_DIMy+1,0:TILE_DIMz+1)
      
      real(kind=db) :: press,u,v,w,pxx,pyy,pzz,pxy,pxz,pyz
      real(kind=db) :: opress,ou,ov,ow,opxx,opyy,opzz,opxy,opxz,opyz,feq,fneq1,f_discr
      real(kind=db) :: mytemp,forcex,forcey,forcez,invrhophi_loc,uu,udotc
      
      real(kind=db) :: omega_loc,phi_loc,visc_loc
#ifdef SMAGORINSKI
	  real(kind=db) :: QQ
#endif

      integer :: i,j,k,myblock,l,intblock
      integer :: ii,jj,kk
      integer :: li,lj,lk
      integer :: lii,ljj,lkk
      integer :: xblock,yblock,zblock
!      integer :: gi,gj,gk

      
      li = threadIdx%x-1
      lj = threadIdx%y-1
      lk = threadIdx%z-1
      
      i = (blockIdx%x-1) * TILE_DIMx + li
      j = (blockIdx%y-1) * TILE_DIMy + lj
      k = (blockIdx%z-1) * TILE_DIMz + lk
      
!      gi=nx*coords(1)+i
!      gj=ny*coords(2)+j
!      gk=nz*coords(3)+k
      
      xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
	  yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
	  zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
      
      myblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
      ii=i-xblock*TILE_DIMx+2*TILE_DIMx
      jj=j-yblock*TILE_DIMy+2*TILE_DIMy
      kk=k-zblock*TILE_DIMz+2*TILE_DIMz

      !i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      !j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      !k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z
      
!      gi=nx*coords(1)+i
!      gj=ny*coords(2)+j
!      gk=nz*coords(3)+k
      
      intblock=blockIdx%x+blockIdx%y*nxblock_d+blockIdx%z*nxyblock_d+1 !internal-node block
      
      !if(li==1 .and. lj==1 .and. lk==1 .and. step==1 .and. blockIdx%x==1 .and. blockIdx%y==1)write(*,*)blockIdx%z,gk

               !if (abs(isfluid(i,j,k)) /= 1) return
       
     
#ifdef TWOCOMPONENT	  
                  phi_loc=real(phifields_s(ii,jj,kk,1,myblock),kind=db)
#endif                  
#ifdef DENSRATIO
                  
                  invrhophi_loc = 1.0_db/(rho_r*phi_loc+(ONE-phi_loc)*rho_b) 
#else
                  invrhophi_loc = 1.0_db !press_loc
#endif	

				  forcex=real(forces_s(ii,jj,kk,1,myblock),kind=db)*invrhophi_loc
				  forcey=real(forces_s(ii,jj,kk,2,myblock),kind=db)*invrhophi_loc
				  forcez=real(forces_s(ii,jj,kk,3,myblock),kind=db)*invrhophi_loc
                  
                  press=real(hfields_in(ii,jj,kk,1,myblock),kind=db)
                  u=real(hfields_in(ii,jj,kk,2,myblock),kind=db)
                  v=real(hfields_in(ii,jj,kk,3,myblock),kind=db)
                  w=real(hfields_in(ii,jj,kk,4,myblock),kind=db)
                  pxx=real(hfields_in(ii,jj,kk,5,myblock),kind=db)
                  pyy=real(hfields_in(ii,jj,kk,6,myblock),kind=db)
                  pzz=real(hfields_in(ii,jj,kk,7,myblock),kind=db)
                  pxy=real(hfields_in(ii,jj,kk,8,myblock),kind=db)
                  pxz=real(hfields_in(ii,jj,kk,9,myblock),kind=db)
                  pyz=real(hfields_in(ii,jj,kk,10,myblock),kind=db) 
                  
                  
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
                  

!				  uu=HALF*(u*u+v*v+w*w)*invcssq
				  
!                  do lii=1,nlinks
!                     udotc=(u*dex(lii) + v*dey(lii)+ w*dez(lii))*invcssq
!		     feq=p(lii)*(press + (udotc+0.5_db*udotc*udotc - uu))
!                     !fneq1=p(l)*(HALF/(cssq*cssq))*( (dex(lii)*dex(lii)-cssq)*pxx &
!		     ! + (dey(lii)*dey(lii)-cssq)*pyy + (dez(lii)*dez(lii)-cssq)*pzz &
!	             ! + TWO*(dex(lii)*dey(lii))*pxy + TWO*(dex(lii)*dez(lii))*pxz &
!		     ! + TWO*(dey(lii)*dez(lii))*pyz)
!                     fpost=feq!+fneq1
!                     pxx=pxx - fpost*(dex(lii)*dex(lii))
!                     pyy=pyy - fpost*(dey(lii)*dey(lii))
!                     pzz=pzz - fpost*(dez(lii)*dez(lii))
!                     pxy=pxy - fpost*(dex(lii)*dey(lii))
!                     pxz=pxz - fpost*(dex(lii)*dez(lii))
!                     pyz=pyz - fpost*(dey(lii)*dez(lii))
!                  enddo

                  pxx=pxx - cssq*press - u*u 
                  pyy=pyy - cssq*press - v*v 
                  pzz=pzz - cssq*press - w*w 
                  pxy=pxy - u*v
                  pxz=pxz - u*w
                  pyz=pyz - v*w

#ifdef TWOCOMPONENT
                  !visc_loc it is used to store the local viscosity
                  visc_loc=(rho_r*visc1*phi_loc+(ONE-phi_loc)*visc2*rho_b)*invrhophi_loc
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
                   sqrt((3.0*visc_loc)**2.0 + 0.053*18.0*sqrt(2.0*QQ)*invrhophi_loc)) !it is tau
                  omega_loc=1.0_db/omega_loc !it is omega

#else
#ifdef TWOCOMPONENT
                  omega_loc=(visc_loc*invcssq + 0.5_db) !it is tau   !visc_loc it is used to store the local viscosity
                  omega_loc=1.0_db/omega_loc !it is omega
#else
                  omega_loc=omega
#endif
#endif      
      
                  !opress=ZERO
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
                  uu=HALF*(u*u+v*v+w*w)*invcssq

			      feq=p(0)*(press - uu)
				  fneq1=p(0)*(HALF*invcssq)*(-pxx-pyy-pzz)
#if defined(HIGHORDER) && (LATTICE == 27)
                  feq=feq + &
                    Minv_mrt(17,0)*u*u*v*v + &
                    Minv_mrt(18,0)*u*u*w*w + &
                    Minv_mrt(19,0)*v*v*w*w + &
                    Minv_mrt(26,0)*u*u*v*v*w*w
                  fneq1=fneq1 + &
                    Minv_mrt(17,0)*(pyy*u**TWO + TWO*v*(THREE*pxy*u + pxx*v)) + &
                    Minv_mrt(18,0)*(pzz*u**TWO + TWO*w*(THREE*pxz*u + pxx*w)) + &
                    Minv_mrt(19,0)*(pzz*v**TWO + TWO*w*(THREE*pyz*v + pyy*w)) + &
                    Minv_mrt(26,0)*(pzz*u**TWO*v**TWO + &
                    w*(SIX*u*v*(pyz*u + pxz*v) + THREE*pyy*u**TWO*w + TWO*v*(SEVEN*pxy*u + TWO*pxx*v)*w))
#endif
				  F_discr = p(0)*(- u*forcex - v*forcey - w*forcez)*invcssq
                  
                  opress=feq + (1.0_db-omega_loc)*fneq1 + HALF*(F_discr)
                  
                  do l=1,nlinks
                     udotc=(u*dex(l) + v*dey(l)+ w*dez(l))*invcssq
		             feq=p(l)*(press + (udotc+0.5_db*udotc*udotc - uu))
                     fneq1=p(l)*(HALF/(cssq*cssq))*( (dex(l)*dex(l)-cssq)*pxx &
		              + (dey(l)*dey(l)-cssq)*pyy + (dez(l)*dez(l)-cssq)*pzz &
	                  + TWO*(dex(l)*dey(l))*pxy + TWO*(dex(l)*dez(l))*pxz &
		              + TWO*(dey(l)*dez(l))*pyz)
#if defined(HIGHORDER) && (LATTICE == 27)
                     feq=feq+Minv_mrt(10,l)*u*u*v + &
                             Minv_mrt(11,l)*u*u*w + &
                             Minv_mrt(12,l)*u*v*v + &
                             Minv_mrt(13,l)*u*w*w + &
                             Minv_mrt(14,l)*v*w*w + &
                             Minv_mrt(15,l)*v*v*w + &
                             Minv_mrt(16,l)*u*v*w + &
                             Minv_mrt(17,l)*u*u*v*v + &
                             Minv_mrt(18,l)*u*u*w*w + &
                             Minv_mrt(19,l)*v*v*w*w + &
                             Minv_mrt(20,l)*u*v*w*w + &
                             Minv_mrt(21,l)*u*v*v*w + &
                             Minv_mrt(22,l)*u*u*v*w + &
                             Minv_mrt(23,l)*u*u*v*w*w + &
                             Minv_mrt(24,l)*u*u*v*v*w + &
                             Minv_mrt(25,l)*u*v*v*w*w + &
                             Minv_mrt(26,l)*u*u*v*v*w*w
                     fneq1=fneq1+Minv_mrt(10,l)*(TWO*pxy*u + pxx*v) + &
                             Minv_mrt(11,l)*(TWO*pxz*u + pxx*w) + &
                             Minv_mrt(12,l)*(pyy*u + TWO*pxy*v) + &
                             Minv_mrt(13,l)*(pzz*u + TWO*pxz*w)+ &
                             Minv_mrt(14,l)*(pzz*v + TWO*pyz*w) + &
                             Minv_mrt(15,l)*(TWO*pyz*v + pyy*w) + &
                             Minv_mrt(16,l)*(pyz*u + pxz*v + pxy*w) + &
                             Minv_mrt(17,l)*(pyy*u**TWO + TWO*v*(THREE*pxy*u + pxx*v)) + &
                             Minv_mrt(18,l)*(pzz*u**TWO + TWO*w*(THREE*pxz*u + pxx*w)) + &
                             Minv_mrt(19,l)*(pzz*v**TWO + TWO*w*(THREE*pyz*v + pyy*w)) + &
                             Minv_mrt(20,l)*(pzz*u*v + w*(THREE*pyz*u + THREE*pxz*v + TWO*pxy*w)) + &
                             Minv_mrt(21,l)*(TWO*pyz*u*v + pxz*v**TWO + TWO*pyy*u*w + FOUR*pxy*v*w) + &
                             Minv_mrt(22,l)*(pyz*u**TWO + TWO*pxz*u*v + FOUR*pxy*u*w + TWO*pxx*v*w) + &
                             Minv_mrt(23,l)*(pzz*u**TWO*v + THREE*w*(pyz*u**TWO + TWO*pxz*u*v + TWO*pxy*u*w + pxx*v*w)) + &
                             Minv_mrt(24,l)*(TWO*u*v*(pyz*u + pxz*v) + TWO*pyy*u**TWO*w + v*(TEN*pxy*u + THREE*pxx*v)*w) + &
                             Minv_mrt(25,l)*(pzz*u*v**TWO + THREE*w*(TWO*pyz*u*v + pxz*v**TWO + pyy*u*w + TWO*pxy*v*w)) + &
                             Minv_mrt(26,l)*(pzz*u**TWO*v**TWO + &
                              w*(SIX*u*v*(pyz*u + pxz*v) + THREE*pyy*u**TWO*w + TWO*v*(SEVEN*pxy*u + TWO*pxx*v)*w))
#endif
		             F_discr = p(l)*(((dex(l) - u) + udotc * dex(l))*forcex &
                      + ((dey(l) - v) + udotc * dey(l))*forcey &
                      + ((dez(l) - w) + udotc * dez(l))*forcez)*invcssq
                     lii=li+ex(l)
                     ljj=lj+ey(l)
                     lkk=lk+ez(l)
                     lii=mod(lii+TILE_DIMx+2,(TILE_DIMx+2))
                     ljj=mod(ljj+TILE_DIMy+2,(TILE_DIMy+2))
                     lkk=mod(lkk+TILE_DIMz+2,(TILE_DIMz+2)) 
		             ftemp(lii,ljj,lkk)=feq + (ONE-omega_loc)*fneq1 + HALF*F_discr
		             !if(gi==32 .and. gj==32 .and. gk==2 .and. myblock==intblock)write(*,*)'a',l,ftemp(lii,ljj,lkk)
!		             if(gi==iprobe .and. gj==jprobe .and. gk==kprobe .and. myblock==intblock)write(*,*)l,feq
!		             if(gi==iprobe .and. gj==jprobe .and. gk==kprobe .and. myblock==intblock)write(*,*)l,(ONE-omega_loc)*fneq1
!		             if(gi==iprobe .and. gj==jprobe .and. gk==kprobe .and. myblock==intblock)write(*,*)l,HALF*F_discr
		             call syncthreads
		             
		             opress=opress + ftemp(li,lj,lk)
		             ou=ou + ftemp(li,lj,lk)*dex(l)
                     ov=ov + ftemp(li,lj,lk)*dey(l)
                     ow=ow + ftemp(li,lj,lk)*dez(l)
                     opxx=opxx + ftemp(li,lj,lk)*dex(l)*dex(l)
                     opyy=opyy + ftemp(li,lj,lk)*dey(l)*dey(l)
                     opzz=opzz + ftemp(li,lj,lk)*dez(l)*dez(l)
                     opxy=opxy + ftemp(li,lj,lk)*dex(l)*dey(l)
                     opxz=opxz + ftemp(li,lj,lk)*dex(l)*dez(l)
                     opyz=opyz + ftemp(li,lj,lk)*dey(l)*dez(l)
                     
                     call syncthreads
                     
                  enddo
                  
                  
                  !If my block index does not match the index of the internal-node block (lii), it means my thread is on the outer. I must exit
	              if(myblock .ne. intblock)return
	                 
	              hfields_out(ii,jj,kk,1,myblock)=real(opress,kind=strdb)
                  hfields_out(ii,jj,kk,2,myblock)=real(ou,kind=strdb)
                  hfields_out(ii,jj,kk,3,myblock)=real(ov,kind=strdb)
                  hfields_out(ii,jj,kk,4,myblock)=real(ow,kind=strdb)
                  hfields_out(ii,jj,kk,5,myblock)=real(opxx,kind=strdb)
                  hfields_out(ii,jj,kk,6,myblock)=real(opyy,kind=strdb)
                  hfields_out(ii,jj,kk,7,myblock)=real(opzz,kind=strdb)
                  hfields_out(ii,jj,kk,8,myblock)=real(opxy,kind=strdb)
                  hfields_out(ii,jj,kk,9,myblock)=real(opxz,kind=strdb)
                  hfields_out(ii,jj,kk,10,myblock)=real(opyz,kind=strdb)   
                   
    endsubroutine fused_LB_kernel1
    

endmodule
