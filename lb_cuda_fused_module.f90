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
   
   
   attributes(global) subroutine fused_LB_kernel_bb(step,iprobe,jprobe,kprobe,flip,flop,nx,ny,nz,coords,isfluid &  
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
      real(kind=db) :: ffront,frear
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

                  pxx=pxx - u*u 
                  pyy=pyy - v*v 
                  pzz=pzz - w*w 
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
                    Minv_mrt(17,0)*(pxx*v**TWO + FOUR*pxy*u*v + pyy*u**TWO) + &
                    Minv_mrt(18,0)*(pxx*w**TWO + FOUR*pxz*u*w + pzz*u**TWO) + &
                    Minv_mrt(19,0)*(pyy*w**TWO + FOUR*pyz*v*w + pzz*v**TWO) + &
                    Minv_mrt(26,0)*(pxx*v**TWO*w**TWO + FOUR*pxy*u*v*w**TWO + FOUR*pxz*u*v**TWO*w + &
                                     pyy*u**TWO*w**TWO + FOUR*pyz*u**TWO*v*w + pzz*u**TWO*v**TWO)
#endif
				  F_discr = p(0)*(- u*forcex - v*forcey - w*forcez)*invcssq
                  udotc=feq + (1.0_db-omega_loc)*fneq1 + HALF*(F_discr)
		          opress=opress + udotc
                  opxx=opxx + udotc*(-cssq)
                  opyy=opyy + udotc*(-cssq)
                  opzz=opzz + udotc*(-cssq)
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
                             Minv_mrt(17,l)*(pxx*v**TWO + FOUR*pxy*u*v + pyy*u**TWO) + &
                             Minv_mrt(18,l)*(pxx*w**TWO + FOUR*pxz*u*w + pzz*u**TWO) + &
                             Minv_mrt(19,l)*(pyy*w**TWO + FOUR*pyz*v*w + pzz*v**TWO) + &
                             Minv_mrt(20,l)*(pxy*w**TWO + TWO*pxz*v*w + TWO*pyz*u*w + pzz*u*v) + &
                             Minv_mrt(21,l)*(TWO*pxy*v*w + pxz*v**TWO + pyy*u*w + TWO*pyz*u*v) + &
                             Minv_mrt(22,l)*(pxx*v*w + TWO*pxy*u*w + TWO*pxz*u*v + pyz*u**TWO) + &
                             Minv_mrt(23,l)*(pxx*v*w**TWO + TWO*pxy*u*w**TWO + FOUR*pxz*u*v*w + TWO*pyz*u**TWO*w + pzz*u**TWO*v) + &
                             Minv_mrt(24,l)*(pxx*v**TWO*w + FOUR*pxy*u*v*w + TWO*pxz*u*v**TWO + pyy*u**TWO*w + TWO*pyz*u**TWO*v) + &
                             Minv_mrt(25,l)*(TWO*pxy*v*w**TWO + TWO*pxz*v**TWO*w + pyy*u*w**TWO + FOUR*pyz*u*v*w + pzz*u*v**TWO) + &
                             Minv_mrt(26,l)*(pxx*v**TWO*w**TWO + FOUR*pxy*u*v*w**TWO + FOUR*pxz*u*v**TWO*w + &
                                              pyy*u**TWO*w**TWO + FOUR*pyz*u**TWO*v*w + pzz*u**TWO*v**TWO)
#endif
		             F_discr = p(l)*(((dex(l) - u) + udotc * dex(l))*forcex &
                      + ((dey(l) - v) + udotc * dey(l))*forcey &
                      + ((dez(l) - w) + udotc * dez(l))*forcez)*invcssq
                     ! pull scheme: stored in place
                     f_front(li,lj,lk)=feq + (ONE-omega_loc)*fneq1 + HALF*F_discr
		             
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
                             Minv_mrt(17,l+1)*(pxx*v**TWO + FOUR*pxy*u*v + pyy*u**TWO) + &
                             Minv_mrt(18,l+1)*(pxx*w**TWO + FOUR*pxz*u*w + pzz*u**TWO) + &
                             Minv_mrt(19,l+1)*(pyy*w**TWO + FOUR*pyz*v*w + pzz*v**TWO) + &
                             Minv_mrt(20,l+1)*(pxy*w**TWO + TWO*pxz*v*w + TWO*pyz*u*w + pzz*u*v) + &
                             Minv_mrt(21,l+1)*(TWO*pxy*v*w + pxz*v**TWO + pyy*u*w + TWO*pyz*u*v) + &
                             Minv_mrt(22,l+1)*(pxx*v*w + TWO*pxy*u*w + TWO*pxz*u*v + pyz*u**TWO) + &
                             Minv_mrt(23,l+1)*(pxx*v*w**TWO + TWO*pxy*u*w**TWO + FOUR*pxz*u*v*w + TWO*pyz*u**TWO*w + pzz*u**TWO*v) + &
                             Minv_mrt(24,l+1)*(pxx*v**TWO*w + FOUR*pxy*u*v*w + TWO*pxz*u*v**TWO + pyy*u**TWO*w + TWO*pyz*u**TWO*v) + &
                             Minv_mrt(25,l+1)*(TWO*pxy*v*w**TWO + TWO*pxz*v**TWO*w + pyy*u*w**TWO + FOUR*pyz*u*v*w + pzz*u*v**TWO) + &
                             Minv_mrt(26,l+1)*(pxx*v**TWO*w**TWO + FOUR*pxy*u*v*w**TWO + FOUR*pxz*u*v**TWO*w + &
                                                pyy*u**TWO*w**TWO + FOUR*pyz*u**TWO*v*w + pzz*u**TWO*v**TWO)
#endif
		             F_discr = p(l+1)*(((dex(l+1) - u) + udotc * dex(l+1))*forcex &
                      + ((dey(l+1) - v) + udotc * dey(l+1))*forcey &
                      + ((dez(l+1) - w) + udotc * dez(l+1))*forcez)*invcssq
                     ! pull scheme: stored in place
                     f_rear(li,lj,lk)=feq + (ONE-omega_loc)*fneq1 + HALF*F_discr
		             
		             call syncthreads

                     ! gather (pull) on the interior nodes of the tile only; l and l+1 are opposite directions
                     if(myblock == intblock)then
                       ffront=f_front(li-ex(l),lj-ey(l),lk-ez(l))
                       frear=f_rear(li-ex(l+1),lj-ey(l+1),lk-ez(l+1))
#ifdef BOUNCE_BACK
                       ! a population pulled from a solid node is the reflected opposite population of this node
                       if(isfluid(i-ex(l),j-ey(l),k-ez(l)) == 0)ffront=f_rear(li,lj,lk)
                       if(isfluid(i-ex(l+1),j-ey(l+1),k-ez(l+1)) == 0)frear=f_front(li,lj,lk)
#endif
                       opress=opress + ffront
                       ou=ou + ffront*dex(l)
                       ov=ov + ffront*dey(l)
                       ow=ow + ffront*dez(l)
                       opxx=opxx + ffront*(dex(l)*dex(l)-cssq)
                       opyy=opyy + ffront*(dey(l)*dey(l)-cssq)
                       opzz=opzz + ffront*(dez(l)*dez(l)-cssq)
                       opxy=opxy + ffront*dex(l)*dey(l)
                       opxz=opxz + ffront*dex(l)*dez(l)
                       opyz=opyz + ffront*dey(l)*dez(l)
                       opress=opress + frear
                       ou=ou + frear*dex(l+1)
                       ov=ov + frear*dey(l+1)
                       ow=ow + frear*dez(l+1)
                       opxx=opxx + frear*(dex(l+1)*dex(l+1)-cssq)
                       opyy=opyy + frear*(dey(l+1)*dey(l+1)-cssq)
                       opzz=opzz + frear*(dez(l+1)*dez(l+1)-cssq)
                       opxy=opxy + frear*dex(l+1)*dey(l+1)
                       opxz=opxz + frear*dex(l+1)*dez(l+1)
                       opyz=opyz + frear*dey(l+1)*dez(l+1)
                     endif
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
              
    endsubroutine fused_LB_kernel_bb  

      attributes(global) subroutine fused_LB_kernel(step,iprobe,jprobe,kprobe,flip,flop,nx,ny,nz,coords,isfluid &  
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
      real(kind=db) :: fpull
      
      real(kind=db) :: omega_loc,phi_loc,visc_loc
#ifdef SMAGORINSKI
	  real(kind=db) :: QQ
#endif

      integer :: i,j,k,myblock,l,intblock
      integer :: ii,jj,kk
      integer :: li,lj,lk
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

                  pxx=pxx - u*u 
                  pyy=pyy - v*v 
                  pzz=pzz - w*w 
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
                    Minv_mrt(17,0)*(pxx*v**TWO + FOUR*pxy*u*v + pyy*u**TWO) + &
                    Minv_mrt(18,0)*(pxx*w**TWO + FOUR*pxz*u*w + pzz*u**TWO) + &
                    Minv_mrt(19,0)*(pyy*w**TWO + FOUR*pyz*v*w + pzz*v**TWO) + &
                    Minv_mrt(26,0)*(pxx*v**TWO*w**TWO + FOUR*pxy*u*v*w**TWO + FOUR*pxz*u*v**TWO*w + &
                                     pyy*u**TWO*w**TWO + FOUR*pyz*u**TWO*v*w + pzz*u**TWO*v**TWO)
#endif
				  F_discr = p(0)*(- u*forcex - v*forcey - w*forcez)*invcssq
                  udotc=feq + (1.0_db-omega_loc)*fneq1 + HALF*(F_discr)
		          opress=opress + udotc
                  opxx=opxx + udotc*(-cssq)
                  opyy=opyy + udotc*(-cssq)
                  opzz=opzz + udotc*(-cssq)
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
                             Minv_mrt(17,l)*(pxx*v**TWO + FOUR*pxy*u*v + pyy*u**TWO) + &
                             Minv_mrt(18,l)*(pxx*w**TWO + FOUR*pxz*u*w + pzz*u**TWO) + &
                             Minv_mrt(19,l)*(pyy*w**TWO + FOUR*pyz*v*w + pzz*v**TWO) + &
                             Minv_mrt(20,l)*(pxy*w**TWO + TWO*pxz*v*w + TWO*pyz*u*w + pzz*u*v) + &
                             Minv_mrt(21,l)*(TWO*pxy*v*w + pxz*v**TWO + pyy*u*w + TWO*pyz*u*v) + &
                             Minv_mrt(22,l)*(pxx*v*w + TWO*pxy*u*w + TWO*pxz*u*v + pyz*u**TWO) + &
                             Minv_mrt(23,l)*(pxx*v*w**TWO + TWO*pxy*u*w**TWO + FOUR*pxz*u*v*w + TWO*pyz*u**TWO*w + pzz*u**TWO*v) + &
                             Minv_mrt(24,l)*(pxx*v**TWO*w + FOUR*pxy*u*v*w + TWO*pxz*u*v**TWO + pyy*u**TWO*w + TWO*pyz*u**TWO*v) + &
                             Minv_mrt(25,l)*(TWO*pxy*v*w**TWO + TWO*pxz*v**TWO*w + pyy*u*w**TWO + FOUR*pyz*u*v*w + pzz*u*v**TWO) + &
                             Minv_mrt(26,l)*(pxx*v**TWO*w**TWO + FOUR*pxy*u*v*w**TWO + FOUR*pxz*u*v**TWO*w + &
                                              pyy*u**TWO*w**TWO + FOUR*pyz*u**TWO*v*w + pzz*u**TWO*v**TWO)
#endif
		             F_discr = p(l)*(((dex(l) - u) + udotc * dex(l))*forcex &
                      + ((dey(l) - v) + udotc * dey(l))*forcey &
                      + ((dez(l) - w) + udotc * dez(l))*forcez)*invcssq
                     ! pull scheme: the population of direction l is stored in place
                     ftemp(li,lj,lk)=feq + (ONE-omega_loc)*fneq1 + HALF*F_discr
		             !if(gi==32 .and. gj==32 .and. gk==2 .and. myblock==intblock)write(*,*)'a',l,ftemp(lii,ljj,lkk)
!		             if(gi==iprobe .and. gj==jprobe .and. gk==kprobe .and. myblock==intblock)write(*,*)l,feq
!		             if(gi==iprobe .and. gj==jprobe .and. gk==kprobe .and. myblock==intblock)write(*,*)l,(ONE-omega_loc)*fneq1
!		             if(gi==iprobe .and. gj==jprobe .and. gk==kprobe .and. myblock==intblock)write(*,*)l,HALF*F_discr
                     call syncthreads

                     ! gather (pull) on the interior nodes of the tile only
                     if(myblock == intblock)then
                       fpull=ftemp(li-ex(l),lj-ey(l),lk-ez(l))
                       opress=opress + fpull
                       ou=ou + fpull*dex(l)
                       ov=ov + fpull*dey(l)
                       ow=ow + fpull*dez(l)
                       opxx=opxx + fpull*(dex(l)*dex(l)-cssq)
                       opyy=opyy + fpull*(dey(l)*dey(l)-cssq)
                       opzz=opzz + fpull*(dez(l)*dez(l)-cssq)
                       opxy=opxy + fpull*dex(l)*dey(l)
                       opxz=opxz + fpull*dex(l)*dez(l)
                       opyz=opyz + fpull*dey(l)*dez(l)
                     endif

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
                   
    endsubroutine fused_LB_kernel

#ifdef FUSED2
!     Second-order fused collide-and-stream kernels for D3Q15, D3Q19 and D3Q27 without HIGHORDER.
!     The post-collision populations are a linear combination of ten Hermite coefficients.
!     Each block has one thread per interior node of the tile. Every thread computes the ten
!     post-collision coefficients of its node and, for the first threads, of one halo node;
!     all of them are stored in dynamic shared memory, and after a single barrier each thread
!     gathers (pull) every coefficient from all directions. The gather stages are in
!     fused2_stages.h, generated by gen_fused2_stages.py. Only the kernel selected by
!     BOUNCE_BACK is compiled.
#ifndef BOUNCE_BACK
      attributes(global) subroutine fused_LB_kernel_2(step,iprobe,jprobe,kprobe,flip,flop,nx,ny,nz,coords,isfluid &  
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
      
      ! all ten post-collision Hermite coefficients of the tile and its halo; the last
      ! dimension is set at launch time through the dynamic shared-memory size
      real(kind=db), shared :: abuf(0:TILE_DIMx+1,0:TILE_DIMy+1,0:TILE_DIMz+1,*)
      
      real(kind=db) :: press,u,v,w,pxx,pyy,pzz,pxy,pxz,pyz
      real(kind=db) :: o0,ox,oy,oz,qxx,qyy,qzz,oxy,oxz,oyz
      real(kind=db) :: forcex,forcey,forcez,invrhophi_loc
      real(kind=db) :: fa,fb,ge,go
      real(kind=db) :: omega_loc,phi_loc,visc_loc
#ifdef SMAGORINSKI
	  real(kind=db) :: QQ
#endif

      integer, parameter :: nhx=TILE_DIMx+2, nhy=TILE_DIMy+2, nhz=TILE_DIMz+2
      integer, parameter :: nthreads=TILE_DIMx*TILE_DIMy*TILE_DIMz
      ! halo nodes of a tile: two full planes, and a ring of 2*nhx+2*TILE_DIMy nodes per interior plane
      integer, parameter :: nring=2*nhx+2*TILE_DIMy, nhalo=nhx*nhy*nhz-nthreads
      integer :: i,j,k,myblock,l,intblock
      integer :: ii,jj,kk
      integer :: li,lj,lk,hi,hj,hk
      integer :: xblock,yblock,zblock
      integer :: ithread,inode,ih,ir

      ! one thread per interior node of the tile; li, lj, lk are its coordinates in the
      ! shared buffer, which includes the halo
      li = threadIdx%x
      lj = threadIdx%y
      lk = threadIdx%z
      ithread=(li-1)+TILE_DIMx*(lj-1)+TILE_DIMx*TILE_DIMy*(lk-1)
      intblock=blockIdx%x+blockIdx%y*nxblock_d+blockIdx%z*nxyblock_d+1 !internal-node block

      ! inode=0 is the node of this thread; inode>0 are the halo nodes ithread, ithread+nthreads, ...
      do inode=0,(nhalo-1)/nthreads+1
        if(inode == 0)then
          hi=li
          hj=lj
          hk=lk
        else
          ih=ithread+(inode-1)*nthreads
          if(ih >= nhalo)exit
          if(ih < nhx*nhy)then
            hk=0
            hj=ih/nhx
            hi=ih-hj*nhx
          elseif(ih < nhx*nhy+TILE_DIMz*nring)then
            ir=ih-nhx*nhy
            hk=ir/nring+1
            ir=ir-(hk-1)*nring
            if(ir < nhx)then
              hj=0
              hi=ir
            elseif(ir < nhx+2*TILE_DIMy)then
              ir=ir-nhx
              hj=ir/2+1
              hi=(nhx-1)*(ir-2*(ir/2))
            else
              hj=nhy-1
              hi=ir-nhx-2*TILE_DIMy
            endif
          else
            ir=ih-nhx*nhy-TILE_DIMz*nring
            hk=nhz-1
            hj=ir/nhx
            hi=ir-hj*nhx
          endif
        endif

        i = (blockIdx%x-1) * TILE_DIMx + hi
        j = (blockIdx%y-1) * TILE_DIMy + hj
        k = (blockIdx%z-1) * TILE_DIMz + hk
        xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
        yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
        zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
        myblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
        ii=i-xblock*TILE_DIMx+2*TILE_DIMx
        jj=j-yblock*TILE_DIMy+2*TILE_DIMy
        kk=k-zblock*TILE_DIMz+2*TILE_DIMz

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
                  


                    pxx=pxx - u*u 
                    pyy=pyy - v*v 
                    pzz=pzz - w*w 
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

                    ! post-collision Hermite coefficients up to second order: equilibrium,
                    ! relaxed non-equilibrium, and half of the Guo forcing term
                    pxx=u*u + (ONE-omega_loc)*pxx + u*forcex
                    pyy=v*v + (ONE-omega_loc)*pyy + v*forcey
                    pzz=w*w + (ONE-omega_loc)*pzz + w*forcez
                    pxy=u*v + (ONE-omega_loc)*pxy + HALF*(u*forcey + v*forcex)
                    pxz=u*w + (ONE-omega_loc)*pxz + HALF*(u*forcez + w*forcex)
                    pyz=v*w + (ONE-omega_loc)*pyz + HALF*(v*forcez + w*forcey)
                    ! trace of the second-order term moved into the zeroth-order coefficient
                    abuf(hi,hj,hk,1)=press - HALF*invcssq*(pxx + pyy + pzz)
                    abuf(hi,hj,hk,2)=u + HALF*forcex
                    abuf(hi,hj,hk,3)=v + HALF*forcey
                    abuf(hi,hj,hk,4)=w + HALF*forcez
                    abuf(hi,hj,hk,5)=pxx
                    abuf(hi,hj,hk,6)=pyy
                    abuf(hi,hj,hk,7)=pzz
                    abuf(hi,hj,hk,8)=pxy
                    abuf(hi,hj,hk,9)=pxz
                    abuf(hi,hj,hk,10)=pyz
      enddo

      call syncthreads

                  o0=ZERO
                  ox=ZERO
                  oy=ZERO
                  oz=ZERO
                  qxx=ZERO
                  qyy=ZERO
                  qzz=ZERO
                  oxy=ZERO
                  oxz=ZERO
                  oyz=ZERO

#include "fused2_stages.h"

                  ! the node of this thread is an interior node: its tile is intblock and ii=li, jj=lj, kk=lk
                  hfields_out(li,lj,lk,1,intblock)=real(o0,kind=strdb)
                  hfields_out(li,lj,lk,2,intblock)=real(ox,kind=strdb)
                  hfields_out(li,lj,lk,3,intblock)=real(oy,kind=strdb)
                  hfields_out(li,lj,lk,4,intblock)=real(oz,kind=strdb)
                  hfields_out(li,lj,lk,5,intblock)=real(qxx-cssq*o0,kind=strdb)
                  hfields_out(li,lj,lk,6,intblock)=real(qyy-cssq*o0,kind=strdb)
                  hfields_out(li,lj,lk,7,intblock)=real(qzz-cssq*o0,kind=strdb)
                  hfields_out(li,lj,lk,8,intblock)=real(oxy,kind=strdb)
                  hfields_out(li,lj,lk,9,intblock)=real(oxz,kind=strdb)
                  hfields_out(li,lj,lk,10,intblock)=real(oyz,kind=strdb)

    endsubroutine fused_LB_kernel_2
#else
      attributes(global) subroutine fused_LB_kernel_bb_2(step,iprobe,jprobe,kprobe,flip,flop,nx,ny,nz,coords,isfluid &  
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
      
      ! all ten post-collision Hermite coefficients of the tile and its halo; the last
      ! dimension is set at launch time through the dynamic shared-memory size
      real(kind=db), shared :: abuf(0:TILE_DIMx+1,0:TILE_DIMy+1,0:TILE_DIMz+1,*)
      
      real(kind=db) :: press,u,v,w,pxx,pyy,pzz,pxy,pxz,pyz
      real(kind=db) :: o0,ox,oy,oz,qxx,qyy,qzz,oxy,oxz,oyz
      real(kind=db) :: forcex,forcey,forcez,invrhophi_loc
      real(kind=db) :: fa,fb,ge,go
      real(kind=db) :: omega_loc,phi_loc,visc_loc
#ifdef SMAGORINSKI
	  real(kind=db) :: QQ
#endif

      integer, parameter :: nhx=TILE_DIMx+2, nhy=TILE_DIMy+2, nhz=TILE_DIMz+2
      integer, parameter :: nthreads=TILE_DIMx*TILE_DIMy*TILE_DIMz
      ! halo nodes of a tile: two full planes, and a ring of 2*nhx+2*TILE_DIMy nodes per interior plane
      integer, parameter :: nring=2*nhx+2*TILE_DIMy, nhalo=nhx*nhy*nhz-nthreads
      integer :: i,j,k,myblock,l,intblock
      integer :: ii,jj,kk
      integer :: li,lj,lk,hi,hj,hk
      integer :: xblock,yblock,zblock
      integer :: ithread,inode,ih,ir
#ifdef BOUNCE_BACK
      integer :: isolid
#endif

      ! one thread per interior node of the tile; li, lj, lk are its coordinates in the
      ! shared buffer, which includes the halo
      li = threadIdx%x
      lj = threadIdx%y
      lk = threadIdx%z
      ithread=(li-1)+TILE_DIMx*(lj-1)+TILE_DIMx*TILE_DIMy*(lk-1)
      intblock=blockIdx%x+blockIdx%y*nxblock_d+blockIdx%z*nxyblock_d+1 !internal-node block

      ! inode=0 is the node of this thread; inode>0 are the halo nodes ithread, ithread+nthreads, ...
      do inode=0,(nhalo-1)/nthreads+1
        if(inode == 0)then
          hi=li
          hj=lj
          hk=lk
        else
          ih=ithread+(inode-1)*nthreads
          if(ih >= nhalo)exit
          if(ih < nhx*nhy)then
            hk=0
            hj=ih/nhx
            hi=ih-hj*nhx
          elseif(ih < nhx*nhy+TILE_DIMz*nring)then
            ir=ih-nhx*nhy
            hk=ir/nring+1
            ir=ir-(hk-1)*nring
            if(ir < nhx)then
              hj=0
              hi=ir
            elseif(ir < nhx+2*TILE_DIMy)then
              ir=ir-nhx
              hj=ir/2+1
              hi=(nhx-1)*(ir-2*(ir/2))
            else
              hj=nhy-1
              hi=ir-nhx-2*TILE_DIMy
            endif
          else
            ir=ih-nhx*nhy-TILE_DIMz*nring
            hk=nhz-1
            hj=ir/nhx
            hi=ir-hj*nhx
          endif
        endif

        i = (blockIdx%x-1) * TILE_DIMx + hi
        j = (blockIdx%y-1) * TILE_DIMy + hj
        k = (blockIdx%z-1) * TILE_DIMz + hk
        xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
        yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
        zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
        myblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
        ii=i-xblock*TILE_DIMx+2*TILE_DIMx
        jj=j-yblock*TILE_DIMy+2*TILE_DIMy
        kk=k-zblock*TILE_DIMz+2*TILE_DIMz

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
                  


                    pxx=pxx - u*u 
                    pyy=pyy - v*v 
                    pzz=pzz - w*w 
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

                    ! post-collision Hermite coefficients up to second order: equilibrium,
                    ! relaxed non-equilibrium, and half of the Guo forcing term
                    pxx=u*u + (ONE-omega_loc)*pxx + u*forcex
                    pyy=v*v + (ONE-omega_loc)*pyy + v*forcey
                    pzz=w*w + (ONE-omega_loc)*pzz + w*forcez
                    pxy=u*v + (ONE-omega_loc)*pxy + HALF*(u*forcey + v*forcex)
                    pxz=u*w + (ONE-omega_loc)*pxz + HALF*(u*forcez + w*forcex)
                    pyz=v*w + (ONE-omega_loc)*pyz + HALF*(v*forcez + w*forcey)
                    ! trace of the second-order term moved into the zeroth-order coefficient
                    abuf(hi,hj,hk,1)=press - HALF*invcssq*(pxx + pyy + pzz)
                    abuf(hi,hj,hk,2)=u + HALF*forcex
                    abuf(hi,hj,hk,3)=v + HALF*forcey
                    abuf(hi,hj,hk,4)=w + HALF*forcez
                    abuf(hi,hj,hk,5)=pxx
                    abuf(hi,hj,hk,6)=pyy
                    abuf(hi,hj,hk,7)=pzz
                    abuf(hi,hj,hk,8)=pxy
                    abuf(hi,hj,hk,9)=pxz
                    abuf(hi,hj,hk,10)=pyz
      enddo

      call syncthreads

#ifdef BOUNCE_BACK
                  ! bit l of isolid is set when the upstream node along direction l is solid
                  i = (blockIdx%x-1) * TILE_DIMx + li
                  j = (blockIdx%y-1) * TILE_DIMy + lj
                  k = (blockIdx%z-1) * TILE_DIMz + lk
                  isolid=0
                  do l=1,nlinks
                    if(isfluid(i-ex(l),j-ey(l),k-ez(l)) == 0)isolid=ibset(isolid,l)
                  enddo
#endif

                  o0=ZERO
                  ox=ZERO
                  oy=ZERO
                  oz=ZERO
                  qxx=ZERO
                  qyy=ZERO
                  qzz=ZERO
                  oxy=ZERO
                  oxz=ZERO
                  oyz=ZERO

#include "fused2_stages.h"

                  ! the node of this thread is an interior node: its tile is intblock and ii=li, jj=lj, kk=lk
                  hfields_out(li,lj,lk,1,intblock)=real(o0,kind=strdb)
                  hfields_out(li,lj,lk,2,intblock)=real(ox,kind=strdb)
                  hfields_out(li,lj,lk,3,intblock)=real(oy,kind=strdb)
                  hfields_out(li,lj,lk,4,intblock)=real(oz,kind=strdb)
                  hfields_out(li,lj,lk,5,intblock)=real(qxx-cssq*o0,kind=strdb)
                  hfields_out(li,lj,lk,6,intblock)=real(qyy-cssq*o0,kind=strdb)
                  hfields_out(li,lj,lk,7,intblock)=real(qzz-cssq*o0,kind=strdb)
                  hfields_out(li,lj,lk,8,intblock)=real(oxy,kind=strdb)
                  hfields_out(li,lj,lk,9,intblock)=real(oxz,kind=strdb)
                  hfields_out(li,lj,lk,10,intblock)=real(oyz,kind=strdb)

    endsubroutine fused_LB_kernel_bb_2
#endif
#endif

endmodule
