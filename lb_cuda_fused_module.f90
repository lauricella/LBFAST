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
       hfields_in,hfields_out,auxfields_s,locauxfields_s,forces_s)

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
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields,nblocks_d) :: phifields_s
#ifdef MONOD
      real(kind=db) :: mu_max,Ks
#endif
#endif           
      real(kind=db) :: visc1,omega,fx,fy,fz
      
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields,nblocks_d) :: hfields_in,hfields_out
      
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nauxfields,nblocks_d) :: auxfields_s
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nlocauxfields,nblocks_d) :: locauxfields_s
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nforces,nblocks_d) :: forces_s
      
      real(kind=db), shared :: f1(0:TILE_DIMx+1,0:TILE_DIMy+1,0:TILE_DIMz+1)
      real(kind=db), shared :: f2(0:TILE_DIMx+1,0:TILE_DIMy+1,0:TILE_DIMz+1)
      
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
                  phi_loc=phifields_s(ii,jj,kk,1,myblock)
#endif                  
#ifdef DENSRATIO
                  
                  invrhophi_loc = 1.0_db/(rho_r*phi_loc+(ONE-phi_loc)*rho_b) 
#else
                  invrhophi_loc = 1.0_db !press_loc
#endif	

				  forcex=forces_s(ii,jj,kk,1,myblock)*invrhophi_loc
				  forcey=forces_s(ii,jj,kk,2,myblock)*invrhophi_loc
				  forcez=forces_s(ii,jj,kk,3,myblock)*invrhophi_loc
                  
                  press=hfields_in(ii,jj,kk,1,myblock)
                  u=hfields_in(ii,jj,kk,2,myblock) 
                  v=hfields_in(ii,jj,kk,3,myblock)
                  w=hfields_in(ii,jj,kk,4,myblock)
                  pxx=hfields_in(ii,jj,kk,5,myblock)
                  pyy=hfields_in(ii,jj,kk,6,myblock)
                  pzz=hfields_in(ii,jj,kk,7,myblock)
                  pxy=hfields_in(ii,jj,kk,8,myblock)
                  pxz=hfields_in(ii,jj,kk,9,myblock)
                  pyz=hfields_in(ii,jj,kk,10,myblock)
                  
                  
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
!                     !fneq1=(HALF/(cssq*cssq))*( (dex(lii)*dex(lii)-cssq)*pxx &
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
				  fneq1=(HALF*invcssq)*(-pxx-pyy-pzz)
				  F_discr = p(0)*(- u*forcex - v*forcey - w*forcez)*invcssq
                  
                  opress=feq + (1.0_db-omega_loc)*fneq1*p(0) + HALF*(F_discr)
                  
                  do l=1,nlinks,2
                     udotc=(u*dex(l) + v*dey(l)+ w*dez(l))*invcssq
		             feq=p(l)*(press + (udotc+0.5_db*udotc*udotc - uu))
                     fneq1=(HALF/(cssq*cssq))*( (dex(l)*dex(l)-cssq)*pxx &
		              + (dey(l)*dey(l)-cssq)*pyy + (dez(l)*dez(l)-cssq)*pzz &
	                  + TWO*(dex(l)*dey(l))*pxy + TWO*(dex(l)*dez(l))*pxz &
		              + TWO*(dey(l)*dez(l))*pyz)
		             F_discr = p(l)*(((dex(l) - u) + udotc * dex(l))*forcex &
                      + ((dey(l) - v) + udotc * dey(l))*forcey &
                      + ((dez(l) - w) + udotc * dez(l))*forcez)*invcssq
                     lii=li+ex(l)
                     ljj=lj+ey(l)
                     lkk=lk+ez(l)
                     lii=mod(lii+TILE_DIMx+2,(TILE_DIMx+2))
                     ljj=mod(ljj+TILE_DIMy+2,(TILE_DIMy+2))
                     lkk=mod(lkk+TILE_DIMz+2,(TILE_DIMz+2)) 
		             f1(lii,ljj,lkk)=feq + (ONE-omega_loc)*p(l)*fneq1 + HALF*F_discr
		             
                     udotc=(u*dex(l+1) + v*dey(l+1)+ w*dez(l+1))*invcssq
		             feq=p(l+1)*(press + (udotc+0.5_db*udotc*udotc - uu))
                     fneq1=(HALF/(cssq*cssq))*( (dex(l+1)*dex(l+1)-cssq)*pxx &
		              + (dey(l+1)*dey(l+1)-cssq)*pyy + (dez(l+1)*dez(l+1)-cssq)*pzz &
	                  + TWO*(dex(l+1)*dey(l+1))*pxy + TWO*(dex(l+1)*dez(l+1))*pxz &
		              + TWO*(dey(l+1)*dez(l+1))*pyz)
		             F_discr = p(l+1)*(((dex(l+1) - u) + udotc * dex(l+1))*forcex &
                      + ((dey(l+1) - v) + udotc * dey(l+1))*forcey &
                      + ((dez(l+1) - w) + udotc * dez(l+1))*forcez)*invcssq
                     lii=li+ex(l+1)
                     ljj=lj+ey(l+1)
                     lkk=lk+ez(l+1)
                     lii=mod(lii+TILE_DIMx+2,(TILE_DIMx+2))
                     ljj=mod(ljj+TILE_DIMy+2,(TILE_DIMy+2))
                     lkk=mod(lkk+TILE_DIMz+2,(TILE_DIMz+2)) 
		             f2(lii,ljj,lkk)=feq + (ONE-omega_loc)*p(l+1)*fneq1 + HALF*F_discr
		             
		             call syncthreads
		             opress=opress + f1(li,lj,lk)
		             ou=ou + f1(li,lj,lk)*dex(l)
                     ov=ov + f1(li,lj,lk)*dey(l)
                     ow=ow + f1(li,lj,lk)*dez(l)
                     opxx=opxx + f1(li,lj,lk)*dex(l)*dex(l)
                     opyy=opyy + f1(li,lj,lk)*dey(l)*dey(l)
                     opzz=opzz + f1(li,lj,lk)*dez(l)*dez(l)
                     opxy=opxy + f1(li,lj,lk)*dex(l)*dey(l)
                     opxz=opxz + f1(li,lj,lk)*dex(l)*dez(l)
                     opyz=opyz + f1(li,lj,lk)*dey(l)*dez(l)
                     
                     opress=opress + f2(li,lj,lk)
		             ou=ou + f2(li,lj,lk)*dex(l+1)
                     ov=ov + f2(li,lj,lk)*dey(l+1)
                     ow=ow + f2(li,lj,lk)*dez(l+1)
                     opxx=opxx + f2(li,lj,lk)*dex(l+1)*dex(l+1)
                     opyy=opyy + f2(li,lj,lk)*dey(l+1)*dey(l+1)
                     opzz=opzz + f2(li,lj,lk)*dez(l+1)*dez(l+1)
                     opxy=opxy + f2(li,lj,lk)*dex(l+1)*dey(l+1)
                     opxz=opxz + f2(li,lj,lk)*dex(l+1)*dez(l+1)
                     opyz=opyz + f2(li,lj,lk)*dey(l+1)*dez(l+1)
                     call syncthreads
                     
                  enddo
                  
                  
                  !If my block index does not match the index of the internal-node block (lii), it means my thread is on the outer. I must exit
	              if(myblock .ne. (blockIdx%x+blockIdx%y*nxblock_d+blockIdx%z*nxyblock_d+1) )return
	                 
	              hfields_out(ii,jj,kk,1,myblock)=opress
                  hfields_out(ii,jj,kk,2,myblock)=ou
                  hfields_out(ii,jj,kk,3,myblock)=ov
                  hfields_out(ii,jj,kk,4,myblock)=ow
                  hfields_out(ii,jj,kk,5,myblock)=opxx
                  hfields_out(ii,jj,kk,6,myblock)=opyy
                  hfields_out(ii,jj,kk,7,myblock)=opzz
                  hfields_out(ii,jj,kk,8,myblock)=opxy
                  hfields_out(ii,jj,kk,9,myblock)=opxz
                  hfields_out(ii,jj,kk,10,myblock)=opyz
                  
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
       hfields_in,hfields_out,auxfields_s,locauxfields_s,forces_s)

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
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields,nblocks_d) :: phifields_s
#ifdef MONOD
      real(kind=db) :: mu_max,Ks
#endif
#endif           
      real(kind=db) :: visc1,omega,fx,fy,fz
      
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields,nblocks_d) :: hfields_in,hfields_out
      
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nauxfields,nblocks_d) :: auxfields_s
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nlocauxfields,nblocks_d) :: locauxfields_s
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nforces,nblocks_d) :: forces_s
      
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


               !if (abs(isfluid(i,j,k)) /= 1) return
       
     
#ifdef TWOCOMPONENT	  
                  phi_loc=phifields_s(ii,jj,kk,1,myblock)
#endif                  
#ifdef DENSRATIO
                  
                  invrhophi_loc = 1.0_db/(rho_r*phi_loc+(ONE-phi_loc)*rho_b) 
#else
                  invrhophi_loc = 1.0_db !press_loc
#endif	

				  forcex=forces_s(ii,jj,kk,1,myblock)*invrhophi_loc
				  forcey=forces_s(ii,jj,kk,2,myblock)*invrhophi_loc
				  forcez=forces_s(ii,jj,kk,3,myblock)*invrhophi_loc
                  
                  press=hfields_in(ii,jj,kk,1,myblock)
                  u=hfields_in(ii,jj,kk,2,myblock) 
                  v=hfields_in(ii,jj,kk,3,myblock)
                  w=hfields_in(ii,jj,kk,4,myblock)
                  pxx=hfields_in(ii,jj,kk,5,myblock)
                  pyy=hfields_in(ii,jj,kk,6,myblock)
                  pzz=hfields_in(ii,jj,kk,7,myblock)
                  pxy=hfields_in(ii,jj,kk,8,myblock)
                  pxz=hfields_in(ii,jj,kk,9,myblock)
                  pyz=hfields_in(ii,jj,kk,10,myblock)
                  
                  
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
!                     !fneq1=(HALF/(cssq*cssq))*( (dex(lii)*dex(lii)-cssq)*pxx &
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
				  fneq1=(HALF*invcssq)*(-pxx-pyy-pzz)
				  F_discr = p(0)*(- u*forcex - v*forcey - w*forcez)*invcssq
                  
                  opress=feq + (1.0_db-omega_loc)*fneq1*p(0) + HALF*(F_discr)
                  
                  do l=1,nlinks
                     udotc=(u*dex(l) + v*dey(l)+ w*dez(l))*invcssq
		             feq=p(l)*(press + (udotc+0.5_db*udotc*udotc - uu))
                     fneq1=(HALF/(cssq*cssq))*( (dex(l)*dex(l)-cssq)*pxx &
		              + (dey(l)*dey(l)-cssq)*pyy + (dez(l)*dez(l)-cssq)*pzz &
	                  + TWO*(dex(l)*dey(l))*pxy + TWO*(dex(l)*dez(l))*pxz &
		              + TWO*(dey(l)*dez(l))*pyz)
		             F_discr = p(l)*(((dex(l) - u) + udotc * dex(l))*forcex &
                      + ((dey(l) - v) + udotc * dey(l))*forcey &
                      + ((dez(l) - w) + udotc * dez(l))*forcez)*invcssq
                     lii=li+ex(l)
                     ljj=lj+ey(l)
                     lkk=lk+ez(l)
                     lii=mod(lii+TILE_DIMx+2,(TILE_DIMx+2))
                     ljj=mod(ljj+TILE_DIMy+2,(TILE_DIMy+2))
                     lkk=mod(lkk+TILE_DIMz+2,(TILE_DIMz+2)) 
		             ftemp(lii,ljj,lkk)=feq + (ONE-omega_loc)*p(l)*fneq1 + HALF*F_discr
!		             if(gi==iprobe .and. gj==jprobe .and. gk==kprobe .and. myblock==intblock)write(*,*)l,ftemp(lii,ljj,lkk)
!		             if(gi==iprobe .and. gj==jprobe .and. gk==kprobe .and. myblock==intblock)write(*,*)l,feq
!		             if(gi==iprobe .and. gj==jprobe .and. gk==kprobe .and. myblock==intblock)write(*,*)l,(ONE-omega_loc)*p(l)*fneq1
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
	                 
	              hfields_out(ii,jj,kk,1,myblock)=opress
                  hfields_out(ii,jj,kk,2,myblock)=ou
                  hfields_out(ii,jj,kk,3,myblock)=ov
                  hfields_out(ii,jj,kk,4,myblock)=ow
                  hfields_out(ii,jj,kk,5,myblock)=opxx
                  hfields_out(ii,jj,kk,6,myblock)=opyy
                  hfields_out(ii,jj,kk,7,myblock)=opzz
                  hfields_out(ii,jj,kk,8,myblock)=opxy
                  hfields_out(ii,jj,kk,9,myblock)=opxz
                  hfields_out(ii,jj,kk,10,myblock)=opyz
                   
    endsubroutine fused_LB_kernel1
    
    
    attributes(global) subroutine fused_LB_kernel_int(step,iprobe,jprobe,kprobe,flip,flop,nx,ny,nz,coords,isfluid &  
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
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields,nblocks_d) :: phifields_s
#ifdef MONOD
      real(kind=db) :: mu_max,Ks
#endif
#endif           
      real(kind=db) :: visc1,omega,fx,fy,fz
      
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields,nblocks_d) :: hfields_in,hfields_out
      
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nauxfields,nblocks_d) :: auxfields_s
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nlocauxfields,nblocks_d) :: locauxfields_s
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nforces,nblocks_d) :: forces_s
      
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
      
      i = (blockIdx%x) * TILE_DIMx + li
      j = (blockIdx%y) * TILE_DIMy + lj
      k = (blockIdx%z) * TILE_DIMz + lk
      
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
      
      intblock=(blockIdx%x+1)+(blockIdx%y+1)*nxblock_d+(blockIdx%z+1)*nxyblock_d+1 !internal-node block


               !if (abs(isfluid(i,j,k)) /= 1) return
       
     
#ifdef TWOCOMPONENT	  
                  phi_loc=phifields_s(ii,jj,kk,1,myblock)
#endif                  
#ifdef DENSRATIO
                  
                  invrhophi_loc = 1.0_db/(rho_r*phi_loc+(ONE-phi_loc)*rho_b) 
#else
                  invrhophi_loc = 1.0_db !press_loc
#endif	

				  forcex=forces_s(ii,jj,kk,1,myblock)*invrhophi_loc
				  forcey=forces_s(ii,jj,kk,2,myblock)*invrhophi_loc
				  forcez=forces_s(ii,jj,kk,3,myblock)*invrhophi_loc
                  
                  press=hfields_in(ii,jj,kk,1,myblock)
                  u=hfields_in(ii,jj,kk,2,myblock) 
                  v=hfields_in(ii,jj,kk,3,myblock)
                  w=hfields_in(ii,jj,kk,4,myblock)
                  pxx=hfields_in(ii,jj,kk,5,myblock)
                  pyy=hfields_in(ii,jj,kk,6,myblock)
                  pzz=hfields_in(ii,jj,kk,7,myblock)
                  pxy=hfields_in(ii,jj,kk,8,myblock)
                  pxz=hfields_in(ii,jj,kk,9,myblock)
                  pyz=hfields_in(ii,jj,kk,10,myblock)
                  
                  
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
				  fneq1=(HALF*invcssq)*(-pxx-pyy-pzz)
				  F_discr = p(0)*(- u*forcex - v*forcey - w*forcez)*invcssq
                  
                  opress=feq + (1.0_db-omega_loc)*fneq1*p(0) + HALF*(F_discr)
                  
                  do l=1,nlinks
                     udotc=(u*dex(l) + v*dey(l)+ w*dez(l))*invcssq
		             feq=p(l)*(press + (udotc+0.5_db*udotc*udotc - uu))
                     fneq1=(HALF/(cssq*cssq))*( (dex(l)*dex(l)-cssq)*pxx &
		              + (dey(l)*dey(l)-cssq)*pyy + (dez(l)*dez(l)-cssq)*pzz &
	                  + TWO*(dex(l)*dey(l))*pxy + TWO*(dex(l)*dez(l))*pxz &
		              + TWO*(dey(l)*dez(l))*pyz)
		             F_discr = p(l)*(((dex(l) - u) + udotc * dex(l))*forcex &
                      + ((dey(l) - v) + udotc * dey(l))*forcey &
                      + ((dez(l) - w) + udotc * dez(l))*forcez)*invcssq
                     lii=li+ex(l)
                     ljj=lj+ey(l)
                     lkk=lk+ez(l)
                     lii=mod(lii+TILE_DIMx+2,(TILE_DIMx+2))
                     ljj=mod(ljj+TILE_DIMy+2,(TILE_DIMy+2))
                     lkk=mod(lkk+TILE_DIMz+2,(TILE_DIMz+2)) 
		             ftemp(lii,ljj,lkk)=feq + (ONE-omega_loc)*p(l)*fneq1 + HALF*F_discr
!		             if(gi==iprobe .and. gj==jprobe .and. gk==kprobe .and. myblock==intblock)write(*,*)l,ftemp(lii,ljj,lkk)
!		             if(gi==iprobe .and. gj==jprobe .and. gk==kprobe .and. myblock==intblock)write(*,*)l,feq
!		             if(gi==iprobe .and. gj==jprobe .and. gk==kprobe .and. myblock==intblock)write(*,*)l,(ONE-omega_loc)*p(l)*fneq1
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
	                 
	              hfields_out(ii,jj,kk,1,myblock)=opress
                  hfields_out(ii,jj,kk,2,myblock)=ou
                  hfields_out(ii,jj,kk,3,myblock)=ov
                  hfields_out(ii,jj,kk,4,myblock)=ow
                  hfields_out(ii,jj,kk,5,myblock)=opxx
                  hfields_out(ii,jj,kk,6,myblock)=opyy
                  hfields_out(ii,jj,kk,7,myblock)=opzz
                  hfields_out(ii,jj,kk,8,myblock)=opxy
                  hfields_out(ii,jj,kk,9,myblock)=opxz
                  hfields_out(ii,jj,kk,10,myblock)=opyz
                   
    endsubroutine fused_LB_kernel_int
    
    attributes(global) subroutine fused_LB_kernel_x(step,iprobe,jprobe,kprobe,flip,flop,nx,ny,nz,coords,isfluid &  
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
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields,nblocks_d) :: phifields_s
#ifdef MONOD
      real(kind=db) :: mu_max,Ks
#endif
#endif           
      real(kind=db) :: visc1,omega,fx,fy,fz
      
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields,nblocks_d) :: hfields_in,hfields_out
      
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nauxfields,nblocks_d) :: auxfields_s
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nlocauxfields,nblocks_d) :: locauxfields_s
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nforces,nblocks_d) :: forces_s
      
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
      integer :: intblockx
!      integer :: gi,gj,gk

      
      li = threadIdx%x-1
      lj = threadIdx%y-1
      lk = threadIdx%z-1
      
      j = (blockIdx%y-1) * TILE_DIMy + lj + TILE_DIMy
      k = (blockIdx%z-1) * TILE_DIMz + lk + TILE_DIMz
      
	  yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
	  zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
      
      i = li
	  !gi=nx*coords(1)+i
	  xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
      
      myblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1

      !i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      !j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      !k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z
      
!      gi=nx*coords(1)+i
!      gj=ny*coords(2)+j
!      gk=nz*coords(3)+k
      intblockx=(1+2*TILE_DIMx-1)/TILE_DIMx
      
      intblock=(intblockx-1)+blockIdx%y*nxblock_d+blockIdx%z*nxyblock_d+1 !internal-node block


               !if (abs(isfluid(i,j,k)) /= 1) return
       
     
#ifdef TWOCOMPONENT	  
                  phi_loc=phifields_s(ii,jj,kk,1,myblock)
#endif                  
#ifdef DENSRATIO
                  
                  invrhophi_loc = 1.0_db/(rho_r*phi_loc+(ONE-phi_loc)*rho_b) 
#else
                  invrhophi_loc = 1.0_db !press_loc
#endif	

				  forcex=forces_s(ii,jj,kk,1,myblock)*invrhophi_loc
				  forcey=forces_s(ii,jj,kk,2,myblock)*invrhophi_loc
				  forcez=forces_s(ii,jj,kk,3,myblock)*invrhophi_loc
                  
                  press=hfields_in(ii,jj,kk,1,myblock)
                  u=hfields_in(ii,jj,kk,2,myblock) 
                  v=hfields_in(ii,jj,kk,3,myblock)
                  w=hfields_in(ii,jj,kk,4,myblock)
                  pxx=hfields_in(ii,jj,kk,5,myblock)
                  pyy=hfields_in(ii,jj,kk,6,myblock)
                  pzz=hfields_in(ii,jj,kk,7,myblock)
                  pxy=hfields_in(ii,jj,kk,8,myblock)
                  pxz=hfields_in(ii,jj,kk,9,myblock)
                  pyz=hfields_in(ii,jj,kk,10,myblock)
                  
                  
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
				  fneq1=(HALF*invcssq)*(-pxx-pyy-pzz)
				  F_discr = p(0)*(- u*forcex - v*forcey - w*forcez)*invcssq
                  
                  opress=feq + (1.0_db-omega_loc)*fneq1*p(0) + HALF*(F_discr)
                  
                  do l=1,nlinks
                     udotc=(u*dex(l) + v*dey(l)+ w*dez(l))*invcssq
		             feq=p(l)*(press + (udotc+0.5_db*udotc*udotc - uu))
                     fneq1=(HALF/(cssq*cssq))*( (dex(l)*dex(l)-cssq)*pxx &
		              + (dey(l)*dey(l)-cssq)*pyy + (dez(l)*dez(l)-cssq)*pzz &
	                  + TWO*(dex(l)*dey(l))*pxy + TWO*(dex(l)*dez(l))*pxz &
		              + TWO*(dey(l)*dez(l))*pyz)
		             F_discr = p(l)*(((dex(l) - u) + udotc * dex(l))*forcex &
                      + ((dey(l) - v) + udotc * dey(l))*forcey &
                      + ((dez(l) - w) + udotc * dez(l))*forcez)*invcssq
                     lii=li+ex(l)
                     ljj=lj+ey(l)
                     lkk=lk+ez(l)
                     lii=mod(lii+TILE_DIMx+2,(TILE_DIMx+2))
                     ljj=mod(ljj+TILE_DIMy+2,(TILE_DIMy+2))
                     lkk=mod(lkk+TILE_DIMz+2,(TILE_DIMz+2)) 
		             ftemp(lii,ljj,lkk)=feq + (ONE-omega_loc)*p(l)*fneq1 + HALF*F_discr
!		             if(gi==iprobe .and. gj==jprobe .and. gk==kprobe .and. myblock==intblock)write(*,*)l,ftemp(lii,ljj,lkk)
!		             if(gi==iprobe .and. gj==jprobe .and. gk==kprobe .and. myblock==intblock)write(*,*)l,feq
!		             if(gi==iprobe .and. gj==jprobe .and. gk==kprobe .and. myblock==intblock)write(*,*)l,(ONE-omega_loc)*p(l)*fneq1
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
	              if(myblock .eq. intblock)then
	                 
	              hfields_out(ii,jj,kk,1,myblock)=opress
                  hfields_out(ii,jj,kk,2,myblock)=ou
                  hfields_out(ii,jj,kk,3,myblock)=ov
                  hfields_out(ii,jj,kk,4,myblock)=ow
                  hfields_out(ii,jj,kk,5,myblock)=opxx
                  hfields_out(ii,jj,kk,6,myblock)=opyy
                  hfields_out(ii,jj,kk,7,myblock)=opzz
                  hfields_out(ii,jj,kk,8,myblock)=opxy
                  hfields_out(ii,jj,kk,9,myblock)=opxz
                  hfields_out(ii,jj,kk,10,myblock)=opyz
                  
                  endif
                  
      i = ((nxblock_d-2)-1) * TILE_DIMx + li
	  !gi=nx*coords(1)+i
	  xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
      
      myblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1

      !i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      !j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      !k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z
      
!      gi=nx*coords(1)+i
!      gj=ny*coords(2)+j
!      gk=nz*coords(3)+k
      intblockx=(((nxblock_d-2)-1) * TILE_DIMx +1+2*TILE_DIMx-1)/TILE_DIMx
      
      intblock=(intblockx-1)+blockIdx%y*nxblock_d+blockIdx%z*nxyblock_d+1 !internal-node block
      
               !if (abs(isfluid(i,j,k)) /= 1) return
       
     
#ifdef TWOCOMPONENT	  
                  phi_loc=phifields_s(ii,jj,kk,1,myblock)
#endif                  
#ifdef DENSRATIO
                  
                  invrhophi_loc = 1.0_db/(rho_r*phi_loc+(ONE-phi_loc)*rho_b) 
#else
                  invrhophi_loc = 1.0_db !press_loc
#endif	

				  forcex=forces_s(ii,jj,kk,1,myblock)*invrhophi_loc
				  forcey=forces_s(ii,jj,kk,2,myblock)*invrhophi_loc
				  forcez=forces_s(ii,jj,kk,3,myblock)*invrhophi_loc
                  
                  press=hfields_in(ii,jj,kk,1,myblock)
                  u=hfields_in(ii,jj,kk,2,myblock) 
                  v=hfields_in(ii,jj,kk,3,myblock)
                  w=hfields_in(ii,jj,kk,4,myblock)
                  pxx=hfields_in(ii,jj,kk,5,myblock)
                  pyy=hfields_in(ii,jj,kk,6,myblock)
                  pzz=hfields_in(ii,jj,kk,7,myblock)
                  pxy=hfields_in(ii,jj,kk,8,myblock)
                  pxz=hfields_in(ii,jj,kk,9,myblock)
                  pyz=hfields_in(ii,jj,kk,10,myblock)
                  
                  
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
				  fneq1=(HALF*invcssq)*(-pxx-pyy-pzz)
				  F_discr = p(0)*(- u*forcex - v*forcey - w*forcez)*invcssq
                  
                  opress=feq + (1.0_db-omega_loc)*fneq1*p(0) + HALF*(F_discr)
                  
                  do l=1,nlinks
                     udotc=(u*dex(l) + v*dey(l)+ w*dez(l))*invcssq
		             feq=p(l)*(press + (udotc+0.5_db*udotc*udotc - uu))
                     fneq1=(HALF/(cssq*cssq))*( (dex(l)*dex(l)-cssq)*pxx &
		              + (dey(l)*dey(l)-cssq)*pyy + (dez(l)*dez(l)-cssq)*pzz &
	                  + TWO*(dex(l)*dey(l))*pxy + TWO*(dex(l)*dez(l))*pxz &
		              + TWO*(dey(l)*dez(l))*pyz)
		             F_discr = p(l)*(((dex(l) - u) + udotc * dex(l))*forcex &
                      + ((dey(l) - v) + udotc * dey(l))*forcey &
                      + ((dez(l) - w) + udotc * dez(l))*forcez)*invcssq
                     lii=li+ex(l)
                     ljj=lj+ey(l)
                     lkk=lk+ez(l)
                     lii=mod(lii+TILE_DIMx+2,(TILE_DIMx+2))
                     ljj=mod(ljj+TILE_DIMy+2,(TILE_DIMy+2))
                     lkk=mod(lkk+TILE_DIMz+2,(TILE_DIMz+2)) 
		             ftemp(lii,ljj,lkk)=feq + (ONE-omega_loc)*p(l)*fneq1 + HALF*F_discr
!		             if(gi==iprobe .and. gj==jprobe .and. gk==kprobe .and. myblock==intblock)write(*,*)l,ftemp(lii,ljj,lkk)
!		             if(gi==iprobe .and. gj==jprobe .and. gk==kprobe .and. myblock==intblock)write(*,*)l,feq
!		             if(gi==iprobe .and. gj==jprobe .and. gk==kprobe .and. myblock==intblock)write(*,*)l,(ONE-omega_loc)*p(l)*fneq1
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
	                 
	              hfields_out(ii,jj,kk,1,myblock)=opress
                  hfields_out(ii,jj,kk,2,myblock)=ou
                  hfields_out(ii,jj,kk,3,myblock)=ov
                  hfields_out(ii,jj,kk,4,myblock)=ow
                  hfields_out(ii,jj,kk,5,myblock)=opxx
                  hfields_out(ii,jj,kk,6,myblock)=opyy
                  hfields_out(ii,jj,kk,7,myblock)=opzz
                  hfields_out(ii,jj,kk,8,myblock)=opxy
                  hfields_out(ii,jj,kk,9,myblock)=opxz
                  hfields_out(ii,jj,kk,10,myblock)=opyz      
                   
    endsubroutine fused_LB_kernel_x

    attributes(global) subroutine fused_LB_kernel_y(step,iprobe,jprobe,kprobe,flip,flop,nx,ny,nz,coords,isfluid &  
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
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields,nblocks_d) :: phifields_s
#ifdef MONOD
      real(kind=db) :: mu_max,Ks
#endif
#endif           
      real(kind=db) :: visc1,omega,fx,fy,fz
      
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields,nblocks_d) :: hfields_in,hfields_out
      
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nauxfields,nblocks_d) :: auxfields_s
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nlocauxfields,nblocks_d) :: locauxfields_s
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nforces,nblocks_d) :: forces_s
      
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
      integer :: intblocky
!      integer :: gi,gj,gk

      
      li = threadIdx%x-1
      lj = threadIdx%y-1
      lk = threadIdx%z-1

      i = (blockIdx%x-1) * TILE_DIMx + li + TILE_DIMx
      k = (blockIdx%z-1) * TILE_DIMz + lk
      
      !gi=nx*coords(1)+i
      !gk=nz*coords(3)+k
      
      xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
	  zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
      
      j = lj
      !gj=ny*coords(2)+j
      yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
      
      myblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1

      !i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      !j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      !k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z
      
!      gi=nx*coords(1)+i
!      gj=ny*coords(2)+j
!      gk=nz*coords(3)+k
      intblocky=(1+2*TILE_DIMy-1)/TILE_DIMy
      
      intblock=blockIdx%x+(intblocky-1)*nxblock_d+blockIdx%z*nxyblock_d+1 !internal-node block


               !if (abs(isfluid(i,j,k)) /= 1) return
       
     
#ifdef TWOCOMPONENT	  
                  phi_loc=phifields_s(ii,jj,kk,1,myblock)
#endif                  
#ifdef DENSRATIO
                  
                  invrhophi_loc = 1.0_db/(rho_r*phi_loc+(ONE-phi_loc)*rho_b) 
#else
                  invrhophi_loc = 1.0_db !press_loc
#endif	

				  forcex=forces_s(ii,jj,kk,1,myblock)*invrhophi_loc
				  forcey=forces_s(ii,jj,kk,2,myblock)*invrhophi_loc
				  forcez=forces_s(ii,jj,kk,3,myblock)*invrhophi_loc
                  
                  press=hfields_in(ii,jj,kk,1,myblock)
                  u=hfields_in(ii,jj,kk,2,myblock) 
                  v=hfields_in(ii,jj,kk,3,myblock)
                  w=hfields_in(ii,jj,kk,4,myblock)
                  pxx=hfields_in(ii,jj,kk,5,myblock)
                  pyy=hfields_in(ii,jj,kk,6,myblock)
                  pzz=hfields_in(ii,jj,kk,7,myblock)
                  pxy=hfields_in(ii,jj,kk,8,myblock)
                  pxz=hfields_in(ii,jj,kk,9,myblock)
                  pyz=hfields_in(ii,jj,kk,10,myblock)
                  
                  
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
				  fneq1=(HALF*invcssq)*(-pxx-pyy-pzz)
				  F_discr = p(0)*(- u*forcex - v*forcey - w*forcez)*invcssq
                  
                  opress=feq + (1.0_db-omega_loc)*fneq1*p(0) + HALF*(F_discr)
                  
                  do l=1,nlinks
                     udotc=(u*dex(l) + v*dey(l)+ w*dez(l))*invcssq
		             feq=p(l)*(press + (udotc+0.5_db*udotc*udotc - uu))
                     fneq1=(HALF/(cssq*cssq))*( (dex(l)*dex(l)-cssq)*pxx &
		              + (dey(l)*dey(l)-cssq)*pyy + (dez(l)*dez(l)-cssq)*pzz &
	                  + TWO*(dex(l)*dey(l))*pxy + TWO*(dex(l)*dez(l))*pxz &
		              + TWO*(dey(l)*dez(l))*pyz)
		             F_discr = p(l)*(((dex(l) - u) + udotc * dex(l))*forcex &
                      + ((dey(l) - v) + udotc * dey(l))*forcey &
                      + ((dez(l) - w) + udotc * dez(l))*forcez)*invcssq
                     lii=li+ex(l)
                     ljj=lj+ey(l)
                     lkk=lk+ez(l)
                     lii=mod(lii+TILE_DIMx+2,(TILE_DIMx+2))
                     ljj=mod(ljj+TILE_DIMy+2,(TILE_DIMy+2))
                     lkk=mod(lkk+TILE_DIMz+2,(TILE_DIMz+2)) 
		             ftemp(lii,ljj,lkk)=feq + (ONE-omega_loc)*p(l)*fneq1 + HALF*F_discr
!		             if(gi==iprobe .and. gj==jprobe .and. gk==kprobe .and. myblock==intblock)write(*,*)l,ftemp(lii,ljj,lkk)
!		             if(gi==iprobe .and. gj==jprobe .and. gk==kprobe .and. myblock==intblock)write(*,*)l,feq
!		             if(gi==iprobe .and. gj==jprobe .and. gk==kprobe .and. myblock==intblock)write(*,*)l,(ONE-omega_loc)*p(l)*fneq1
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
	              if(myblock .eq. intblock)then
	                 
	              hfields_out(ii,jj,kk,1,myblock)=opress
                  hfields_out(ii,jj,kk,2,myblock)=ou
                  hfields_out(ii,jj,kk,3,myblock)=ov
                  hfields_out(ii,jj,kk,4,myblock)=ow
                  hfields_out(ii,jj,kk,5,myblock)=opxx
                  hfields_out(ii,jj,kk,6,myblock)=opyy
                  hfields_out(ii,jj,kk,7,myblock)=opzz
                  hfields_out(ii,jj,kk,8,myblock)=opxy
                  hfields_out(ii,jj,kk,9,myblock)=opxz
                  hfields_out(ii,jj,kk,10,myblock)=opyz
                  
                  endif
                  
      j = ((nyblock_d-2)-1) * TILE_DIMy + lj
      !gj=ny*coords(2)+j
      yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
      
      myblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1

      !i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      !j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      !k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z
      
!      gi=nx*coords(1)+i
!      gj=ny*coords(2)+j
!      gk=nz*coords(3)+k
      intblocky=(((nyblock_d-2)-1) * TILE_DIMy +1+2*TILE_DIMy-1)/TILE_DIMy
      
      intblock=blockIdx%x+(intblocky-1)*nxblock_d+blockIdx%z*nxyblock_d+1 !internal-node block
      
               !if (abs(isfluid(i,j,k)) /= 1) return
       
     
#ifdef TWOCOMPONENT	  
                  phi_loc=phifields_s(ii,jj,kk,1,myblock)
#endif                  
#ifdef DENSRATIO
                  
                  invrhophi_loc = 1.0_db/(rho_r*phi_loc+(ONE-phi_loc)*rho_b) 
#else
                  invrhophi_loc = 1.0_db !press_loc
#endif	

				  forcex=forces_s(ii,jj,kk,1,myblock)*invrhophi_loc
				  forcey=forces_s(ii,jj,kk,2,myblock)*invrhophi_loc
				  forcez=forces_s(ii,jj,kk,3,myblock)*invrhophi_loc
                  
                  press=hfields_in(ii,jj,kk,1,myblock)
                  u=hfields_in(ii,jj,kk,2,myblock) 
                  v=hfields_in(ii,jj,kk,3,myblock)
                  w=hfields_in(ii,jj,kk,4,myblock)
                  pxx=hfields_in(ii,jj,kk,5,myblock)
                  pyy=hfields_in(ii,jj,kk,6,myblock)
                  pzz=hfields_in(ii,jj,kk,7,myblock)
                  pxy=hfields_in(ii,jj,kk,8,myblock)
                  pxz=hfields_in(ii,jj,kk,9,myblock)
                  pyz=hfields_in(ii,jj,kk,10,myblock)
                  
                  
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
				  fneq1=(HALF*invcssq)*(-pxx-pyy-pzz)
				  F_discr = p(0)*(- u*forcex - v*forcey - w*forcez)*invcssq
                  
                  opress=feq + (1.0_db-omega_loc)*fneq1*p(0) + HALF*(F_discr)
                  
                  do l=1,nlinks
                     udotc=(u*dex(l) + v*dey(l)+ w*dez(l))*invcssq
		             feq=p(l)*(press + (udotc+0.5_db*udotc*udotc - uu))
                     fneq1=(HALF/(cssq*cssq))*( (dex(l)*dex(l)-cssq)*pxx &
		              + (dey(l)*dey(l)-cssq)*pyy + (dez(l)*dez(l)-cssq)*pzz &
	                  + TWO*(dex(l)*dey(l))*pxy + TWO*(dex(l)*dez(l))*pxz &
		              + TWO*(dey(l)*dez(l))*pyz)
		             F_discr = p(l)*(((dex(l) - u) + udotc * dex(l))*forcex &
                      + ((dey(l) - v) + udotc * dey(l))*forcey &
                      + ((dez(l) - w) + udotc * dez(l))*forcez)*invcssq
                     lii=li+ex(l)
                     ljj=lj+ey(l)
                     lkk=lk+ez(l)
                     lii=mod(lii+TILE_DIMx+2,(TILE_DIMx+2))
                     ljj=mod(ljj+TILE_DIMy+2,(TILE_DIMy+2))
                     lkk=mod(lkk+TILE_DIMz+2,(TILE_DIMz+2)) 
		             ftemp(lii,ljj,lkk)=feq + (ONE-omega_loc)*p(l)*fneq1 + HALF*F_discr
!		             if(gi==iprobe .and. gj==jprobe .and. gk==kprobe .and. myblock==intblock)write(*,*)l,ftemp(lii,ljj,lkk)
!		             if(gi==iprobe .and. gj==jprobe .and. gk==kprobe .and. myblock==intblock)write(*,*)l,feq
!		             if(gi==iprobe .and. gj==jprobe .and. gk==kprobe .and. myblock==intblock)write(*,*)l,(ONE-omega_loc)*p(l)*fneq1
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
	                 
	              hfields_out(ii,jj,kk,1,myblock)=opress
                  hfields_out(ii,jj,kk,2,myblock)=ou
                  hfields_out(ii,jj,kk,3,myblock)=ov
                  hfields_out(ii,jj,kk,4,myblock)=ow
                  hfields_out(ii,jj,kk,5,myblock)=opxx
                  hfields_out(ii,jj,kk,6,myblock)=opyy
                  hfields_out(ii,jj,kk,7,myblock)=opzz
                  hfields_out(ii,jj,kk,8,myblock)=opxy
                  hfields_out(ii,jj,kk,9,myblock)=opxz
                  hfields_out(ii,jj,kk,10,myblock)=opyz      
                   
    endsubroutine fused_LB_kernel_y

    attributes(global) subroutine fused_LB_kernel_z(step,iprobe,jprobe,kprobe,flip,flop,nx,ny,nz,coords,isfluid &  
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
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields,nblocks_d) :: phifields_s
#ifdef MONOD
      real(kind=db) :: mu_max,Ks
#endif
#endif           
      real(kind=db) :: visc1,omega,fx,fy,fz
      
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields,nblocks_d) :: hfields_in,hfields_out
      
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nauxfields,nblocks_d) :: auxfields_s
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nlocauxfields,nblocks_d) :: locauxfields_s
      real(kind=db), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nforces,nblocks_d) :: forces_s
      
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
      integer :: intblockz
!      integer :: gi,gj,gk

      
      li = threadIdx%x-1
      lj = threadIdx%y-1
      lk = threadIdx%z-1

      i = (blockIdx%x-1) * TILE_DIMx + li
      j = (blockIdx%y-1) * TILE_DIMy + lj
      
      !gi=nx*coords(1)+i
      !gj=ny*coords(2)+j
      
      xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
	  yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
	  
	  
	  k = lk
	  !gk=nz*coords(3)+k
	  zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
      
      myblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1

      !i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      !j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      !k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z
      
!      gi=nx*coords(1)+i
!      gj=ny*coords(2)+j
!      gk=nz*coords(3)+k
      intblockz=(1+2*TILE_DIMz-1)/TILE_DIMz
      
      intblock=blockIdx%x+threadIdx%y*nxblock_d+(intblockz-1)*nxyblock_d+1 !internal-node block


               !if (abs(isfluid(i,j,k)) /= 1) return
       
     
#ifdef TWOCOMPONENT	  
                  phi_loc=phifields_s(ii,jj,kk,1,myblock)
#endif                  
#ifdef DENSRATIO
                  
                  invrhophi_loc = 1.0_db/(rho_r*phi_loc+(ONE-phi_loc)*rho_b) 
#else
                  invrhophi_loc = 1.0_db !press_loc
#endif	

				  forcex=forces_s(ii,jj,kk,1,myblock)*invrhophi_loc
				  forcey=forces_s(ii,jj,kk,2,myblock)*invrhophi_loc
				  forcez=forces_s(ii,jj,kk,3,myblock)*invrhophi_loc
                  
                  press=hfields_in(ii,jj,kk,1,myblock)
                  u=hfields_in(ii,jj,kk,2,myblock) 
                  v=hfields_in(ii,jj,kk,3,myblock)
                  w=hfields_in(ii,jj,kk,4,myblock)
                  pxx=hfields_in(ii,jj,kk,5,myblock)
                  pyy=hfields_in(ii,jj,kk,6,myblock)
                  pzz=hfields_in(ii,jj,kk,7,myblock)
                  pxy=hfields_in(ii,jj,kk,8,myblock)
                  pxz=hfields_in(ii,jj,kk,9,myblock)
                  pyz=hfields_in(ii,jj,kk,10,myblock)
                  
                  
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
				  fneq1=(HALF*invcssq)*(-pxx-pyy-pzz)
				  F_discr = p(0)*(- u*forcex - v*forcey - w*forcez)*invcssq
                  
                  opress=feq + (1.0_db-omega_loc)*fneq1*p(0) + HALF*(F_discr)
                  
                  do l=1,nlinks
                     udotc=(u*dex(l) + v*dey(l)+ w*dez(l))*invcssq
		             feq=p(l)*(press + (udotc+0.5_db*udotc*udotc - uu))
                     fneq1=(HALF/(cssq*cssq))*( (dex(l)*dex(l)-cssq)*pxx &
		              + (dey(l)*dey(l)-cssq)*pyy + (dez(l)*dez(l)-cssq)*pzz &
	                  + TWO*(dex(l)*dey(l))*pxy + TWO*(dex(l)*dez(l))*pxz &
		              + TWO*(dey(l)*dez(l))*pyz)
		             F_discr = p(l)*(((dex(l) - u) + udotc * dex(l))*forcex &
                      + ((dey(l) - v) + udotc * dey(l))*forcey &
                      + ((dez(l) - w) + udotc * dez(l))*forcez)*invcssq
                     lii=li+ex(l)
                     ljj=lj+ey(l)
                     lkk=lk+ez(l)
                     lii=mod(lii+TILE_DIMx+2,(TILE_DIMx+2))
                     ljj=mod(ljj+TILE_DIMy+2,(TILE_DIMy+2))
                     lkk=mod(lkk+TILE_DIMz+2,(TILE_DIMz+2)) 
		             ftemp(lii,ljj,lkk)=feq + (ONE-omega_loc)*p(l)*fneq1 + HALF*F_discr
!		             if(gi==iprobe .and. gj==jprobe .and. gk==kprobe .and. myblock==intblock)write(*,*)l,ftemp(lii,ljj,lkk)
!		             if(gi==iprobe .and. gj==jprobe .and. gk==kprobe .and. myblock==intblock)write(*,*)l,feq
!		             if(gi==iprobe .and. gj==jprobe .and. gk==kprobe .and. myblock==intblock)write(*,*)l,(ONE-omega_loc)*p(l)*fneq1
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
	              if(myblock .eq. intblock)then
	                 
	              hfields_out(ii,jj,kk,1,myblock)=opress
                  hfields_out(ii,jj,kk,2,myblock)=ou
                  hfields_out(ii,jj,kk,3,myblock)=ov
                  hfields_out(ii,jj,kk,4,myblock)=ow
                  hfields_out(ii,jj,kk,5,myblock)=opxx
                  hfields_out(ii,jj,kk,6,myblock)=opyy
                  hfields_out(ii,jj,kk,7,myblock)=opzz
                  hfields_out(ii,jj,kk,8,myblock)=opxy
                  hfields_out(ii,jj,kk,9,myblock)=opxz
                  hfields_out(ii,jj,kk,10,myblock)=opyz
                  
                  endif
                  
      k = ((nzblock_d-2)-1) * TILE_DIMz + lk
	  !gk=nz*coords(3)+k
	  zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
      
      myblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1

      !i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      !j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      !k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z
      
!      gi=nx*coords(1)+i
!      gj=ny*coords(2)+j
!      gk=nz*coords(3)+k
      intblockz=(((nzblock_d-2)-1) * TILE_DIMz +1+2*TILE_DIMz-1)/TILE_DIMz
      
      intblock=blockIdx%x+threadIdx%y*nxblock_d+(intblockz-1)*nxyblock_d+1 !internal-node block
      
               !if (abs(isfluid(i,j,k)) /= 1) return
       
     
#ifdef TWOCOMPONENT	  
                  phi_loc=phifields_s(ii,jj,kk,1,myblock)
#endif                  
#ifdef DENSRATIO
                  
                  invrhophi_loc = 1.0_db/(rho_r*phi_loc+(ONE-phi_loc)*rho_b) 
#else
                  invrhophi_loc = 1.0_db !press_loc
#endif	

				  forcex=forces_s(ii,jj,kk,1,myblock)*invrhophi_loc
				  forcey=forces_s(ii,jj,kk,2,myblock)*invrhophi_loc
				  forcez=forces_s(ii,jj,kk,3,myblock)*invrhophi_loc
                  
                  press=hfields_in(ii,jj,kk,1,myblock)
                  u=hfields_in(ii,jj,kk,2,myblock) 
                  v=hfields_in(ii,jj,kk,3,myblock)
                  w=hfields_in(ii,jj,kk,4,myblock)
                  pxx=hfields_in(ii,jj,kk,5,myblock)
                  pyy=hfields_in(ii,jj,kk,6,myblock)
                  pzz=hfields_in(ii,jj,kk,7,myblock)
                  pxy=hfields_in(ii,jj,kk,8,myblock)
                  pxz=hfields_in(ii,jj,kk,9,myblock)
                  pyz=hfields_in(ii,jj,kk,10,myblock)
                  
                  
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
				  fneq1=(HALF*invcssq)*(-pxx-pyy-pzz)
				  F_discr = p(0)*(- u*forcex - v*forcey - w*forcez)*invcssq
                  
                  opress=feq + (1.0_db-omega_loc)*fneq1*p(0) + HALF*(F_discr)
                  
                  do l=1,nlinks
                     udotc=(u*dex(l) + v*dey(l)+ w*dez(l))*invcssq
		             feq=p(l)*(press + (udotc+0.5_db*udotc*udotc - uu))
                     fneq1=(HALF/(cssq*cssq))*( (dex(l)*dex(l)-cssq)*pxx &
		              + (dey(l)*dey(l)-cssq)*pyy + (dez(l)*dez(l)-cssq)*pzz &
	                  + TWO*(dex(l)*dey(l))*pxy + TWO*(dex(l)*dez(l))*pxz &
		              + TWO*(dey(l)*dez(l))*pyz)
		             F_discr = p(l)*(((dex(l) - u) + udotc * dex(l))*forcex &
                      + ((dey(l) - v) + udotc * dey(l))*forcey &
                      + ((dez(l) - w) + udotc * dez(l))*forcez)*invcssq
                     lii=li+ex(l)
                     ljj=lj+ey(l)
                     lkk=lk+ez(l)
                     lii=mod(lii+TILE_DIMx+2,(TILE_DIMx+2))
                     ljj=mod(ljj+TILE_DIMy+2,(TILE_DIMy+2))
                     lkk=mod(lkk+TILE_DIMz+2,(TILE_DIMz+2)) 
		             ftemp(lii,ljj,lkk)=feq + (ONE-omega_loc)*p(l)*fneq1 + HALF*F_discr
!		             if(gi==iprobe .and. gj==jprobe .and. gk==kprobe .and. myblock==intblock)write(*,*)l,ftemp(lii,ljj,lkk)
!		             if(gi==iprobe .and. gj==jprobe .and. gk==kprobe .and. myblock==intblock)write(*,*)l,feq
!		             if(gi==iprobe .and. gj==jprobe .and. gk==kprobe .and. myblock==intblock)write(*,*)l,(ONE-omega_loc)*p(l)*fneq1
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
	                 
	              hfields_out(ii,jj,kk,1,myblock)=opress
                  hfields_out(ii,jj,kk,2,myblock)=ou
                  hfields_out(ii,jj,kk,3,myblock)=ov
                  hfields_out(ii,jj,kk,4,myblock)=ow
                  hfields_out(ii,jj,kk,5,myblock)=opxx
                  hfields_out(ii,jj,kk,6,myblock)=opyy
                  hfields_out(ii,jj,kk,7,myblock)=opzz
                  hfields_out(ii,jj,kk,8,myblock)=opxy
                  hfields_out(ii,jj,kk,9,myblock)=opxz
                  hfields_out(ii,jj,kk,10,myblock)=opyz      
                   
    endsubroutine fused_LB_kernel_z

endmodule
