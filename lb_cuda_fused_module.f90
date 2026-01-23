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
       ,visc1,omega,fx,fy,fz,ntothfields,ntotphifields,ntotauxfields,ntotlocauxfields,ntotforces &
       ,auxfields_s,locauxfields_s,forces_s &
       ,press_in,u_in,v_in,w_in,pxx_in,pyy_in,pzz_in,pxy_in,pxz_in,pyz_in &
       ,press_out,u_out,v_out,w_out,pxx_out,pyy_out,pzz_out,pxy_out,pxz_out,pyz_out)

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
      real(kind=db), dimension(ntotphifields) :: phifields_s
#ifdef MONOD
      real(kind=db) :: mu_max,Ks
#endif
#endif           
      real(kind=db) :: visc1,omega,fx,fy,fz
      
      
      
      real(kind=db), dimension(ntotauxfields) :: auxfields_s
      real(kind=db), dimension(ntotlocauxfields) :: locauxfields_s
      real(kind=db), dimension(ntotforces) :: forces_s
      
      real(kind=db), dimension(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff) :: &
       press_in,u_in,v_in,w_in,pxx_in,pyy_in,pzz_in,pxy_in,pxz_in,pyz_in, &
       press_out,u_out,v_out,w_out,pxx_out,pyy_out,pzz_out,pxy_out,pxz_out,pyz_out
      
      real(kind=db), shared :: f1(0:TILE_DIMx+1,0:TILE_DIMy+1,0:TILE_DIMz+1)
      real(kind=db), shared :: f2(0:TILE_DIMx+1,0:TILE_DIMy+1,0:TILE_DIMz+1)
      
      real(kind=db) :: pxx,pyy,pzz,pxy,pxz,pyz
      real(kind=db) :: opress,ou,ov,ow,opxx,opyy,opzz,opxy,opxz,opyz,feq,fneq1,f_discr
      real(kind=db) :: mytemp,forcex,forcey,forcez,rhophi_loc,uu,udotc
      
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
                  phi_loc=phifields_s(idx5d(ii,jj,kk,1,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields))
#endif                  
#ifdef DENSRATIO
                  
                  rhophi_loc = rho_r*phi_loc+(ONE-phi_loc)*rho_b 
#else
                  rhophi_loc = 1.0_db !press_loc
#endif	

				  forcex=forces_s(idx5d(ii,jj,kk,1,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nforces))/rhophi_loc
				  forcey=forces_s(idx5d(ii,jj,kk,2,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nforces))/rhophi_loc
				  forcez=forces_s(idx5d(ii,jj,kk,3,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nforces))/rhophi_loc
                  

                  pxx=pxx_in(i,j,k)
                  pyy=pyy_in(i,j,k)
                  pzz=pzz_in(i,j,k)
                  pxy=pxy_in(i,j,k)
                  pxz=pxz_in(i,j,k)
                  pyz=pyz_in(i,j,k)
                  
                  
                  

!				  uu=HALF*(u_in(i,j,k)*u_in(i,j,k)+v_in(i,j,k)*v_in(i,j,k)+w_in(i,j,k)*w_in(i,j,k))*invcssq
				  
!                  do lii=1,nlinks
!                     udotc=(u_in(i,j,k)*dex(lii) + v_in(i,j,k)*dey(lii)+ w_in(i,j,k)*dez(lii))*invcssq
!		     feq=p(lii)*(press_in(i,j,k) + (udotc+0.5_db*udotc*udotc - uu))
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

                  pxx=pxx - cssq*press_in(i,j,k) - u_in(i,j,k)*u_in(i,j,k) 
                  pyy=pyy - cssq*press_in(i,j,k) - v_in(i,j,k)*v_in(i,j,k) 
                  pzz=pzz - cssq*press_in(i,j,k) - w_in(i,j,k)*w_in(i,j,k) 
                  pxy=pxy - u_in(i,j,k)*v_in(i,j,k)
                  pxz=pxz - u_in(i,j,k)*w_in(i,j,k)
                  pyz=pyz - v_in(i,j,k)*w_in(i,j,k)

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
                  omega_loc=(visc_loc*invcssq + 0.5_db) !it is tau   !visc_loc it is used to store the local viscosity
                  omega_loc=1.0_db/omega_loc !it is omega
#else
                  omega_loc=omega
#endif
#endif      
      
                  !press_out(i,j,k)=ZERO
                  u_out(i,j,k)=ZERO
                  v_out(i,j,k)=ZERO
                  w_out(i,j,k)=ZERO
                  pxx_out(i,j,k)=ZERO
                  pyy_out(i,j,k)=ZERO
                  pzz_out(i,j,k)=ZERO
                  pxy_out(i,j,k)=ZERO
                  pxz_out(i,j,k)=ZERO
                  pyz_out(i,j,k)=ZERO
!!!!!!!!!!!!!!!!!!!!!!!!!!0
                  uu=HALF*(u_in(i,j,k)*u_in(i,j,k)+v_in(i,j,k)*v_in(i,j,k)+w_in(i,j,k)*w_in(i,j,k))*invcssq

			      feq=p(0)*(press_in(i,j,k) - uu)
				  fneq1=(HALF*invcssq)*(-pxx-pyy-pzz)
				  F_discr = p(0)*(- u_in(i,j,k)*forcex - v_in(i,j,k)*forcey - w_in(i,j,k)*forcez)*invcssq
                  
                  press_out(i,j,k)=feq + (1.0_db-omega_loc)*fneq1*p(0) + HALF*(F_discr)
                  
                  do l=1,nlinks,2
                     udotc=(u_in(i,j,k)*dex(l) + v_in(i,j,k)*dey(l)+ w_in(i,j,k)*dez(l))*invcssq
		             feq=p(l)*(press_in(i,j,k) + (udotc+0.5_db*udotc*udotc - uu))
                     fneq1=(HALF/(cssq*cssq))*( (dex(l)*dex(l)-cssq)*pxx &
		              + (dey(l)*dey(l)-cssq)*pyy + (dez(l)*dez(l)-cssq)*pzz &
	                  + TWO*(dex(l)*dey(l))*pxy + TWO*(dex(l)*dez(l))*pxz &
		              + TWO*(dey(l)*dez(l))*pyz)
		             F_discr = p(l)*(((dex(l) - u_in(i,j,k)) + udotc * dex(l))*forcex &
                      + ((dey(l) - v_in(i,j,k)) + udotc * dey(l))*forcey &
                      + ((dez(l) - w_in(i,j,k)) + udotc * dez(l))*forcez)*invcssq
                     lii=li+ex(l)
                     ljj=lj+ey(l)
                     lkk=lk+ez(l)
                     lii=mod(lii+TILE_DIMx+2,(TILE_DIMx+2))
                     ljj=mod(ljj+TILE_DIMy+2,(TILE_DIMy+2))
                     lkk=mod(lkk+TILE_DIMz+2,(TILE_DIMz+2)) 
		             f1(lii,ljj,lkk)=feq + (ONE-omega_loc)*p(l)*fneq1 + HALF*F_discr
		             
                     udotc=(u_in(i,j,k)*dex(l+1) + v_in(i,j,k)*dey(l+1)+ w_in(i,j,k)*dez(l+1))*invcssq
		             feq=p(l+1)*(press_in(i,j,k) + (udotc+0.5_db*udotc*udotc - uu))
                     fneq1=(HALF/(cssq*cssq))*( (dex(l+1)*dex(l+1)-cssq)*pxx &
		              + (dey(l+1)*dey(l+1)-cssq)*pyy + (dez(l+1)*dez(l+1)-cssq)*pzz &
	                  + TWO*(dex(l+1)*dey(l+1))*pxy + TWO*(dex(l+1)*dez(l+1))*pxz &
		              + TWO*(dey(l+1)*dez(l+1))*pyz)
		             F_discr = p(l+1)*(((dex(l+1) - u_in(i,j,k)) + udotc * dex(l+1))*forcex &
                      + ((dey(l+1) - v_in(i,j,k)) + udotc * dey(l+1))*forcey &
                      + ((dez(l+1) - w_in(i,j,k)) + udotc * dez(l+1))*forcez)*invcssq
                     lii=li+ex(l+1)
                     ljj=lj+ey(l+1)
                     lkk=lk+ez(l+1)
                     lii=mod(lii+TILE_DIMx+2,(TILE_DIMx+2))
                     ljj=mod(ljj+TILE_DIMy+2,(TILE_DIMy+2))
                     lkk=mod(lkk+TILE_DIMz+2,(TILE_DIMz+2)) 
		             f2(lii,ljj,lkk)=feq + (ONE-omega_loc)*p(l+1)*fneq1 + HALF*F_discr
		             
		             call syncthreads
		             if(myblock .eq. intblock)then
		             press_out(i,j,k)=press_out(i,j,k) + f1(li,lj,lk)
		             u_out(i,j,k)=u_out(i,j,k) + f1(li,lj,lk)*dex(l)
                     v_out(i,j,k)=v_out(i,j,k) + f1(li,lj,lk)*dey(l)
                     w_out(i,j,k)=w_out(i,j,k) + f1(li,lj,lk)*dez(l)
                     pxx_out(i,j,k)=pxx_out(i,j,k) + f1(li,lj,lk)*dex(l)*dex(l)
                     pyy_out(i,j,k)=pyy_out(i,j,k) + f1(li,lj,lk)*dey(l)*dey(l)
                     pzz_out(i,j,k)=pzz_out(i,j,k) + f1(li,lj,lk)*dez(l)*dez(l)
                     pxy_out(i,j,k)=pxy_out(i,j,k) + f1(li,lj,lk)*dex(l)*dey(l)
                     pxz_out(i,j,k)=pxz_out(i,j,k) + f1(li,lj,lk)*dex(l)*dez(l)
                     pyz_out(i,j,k)=pyz_out(i,j,k) + f1(li,lj,lk)*dey(l)*dez(l)
                     
                     press_out(i,j,k)=press_out(i,j,k) + f2(li,lj,lk)
		             u_out(i,j,k)=u_out(i,j,k) + f2(li,lj,lk)*dex(l+1)
                     v_out(i,j,k)=v_out(i,j,k) + f2(li,lj,lk)*dey(l+1)
                     w_out(i,j,k)=w_out(i,j,k) + f2(li,lj,lk)*dez(l+1)
                     pxx_out(i,j,k)=pxx_out(i,j,k) + f2(li,lj,lk)*dex(l+1)*dex(l+1)
                     pyy_out(i,j,k)=pyy_out(i,j,k) + f2(li,lj,lk)*dey(l+1)*dey(l+1)
                     pzz_out(i,j,k)=pzz_out(i,j,k) + f2(li,lj,lk)*dez(l+1)*dez(l+1)
                     pxy_out(i,j,k)=pxy_out(i,j,k) + f2(li,lj,lk)*dex(l+1)*dey(l+1)
                     pxz_out(i,j,k)=pxz_out(i,j,k) + f2(li,lj,lk)*dex(l+1)*dez(l+1)
                     pyz_out(i,j,k)=pyz_out(i,j,k) + f2(li,lj,lk)*dey(l+1)*dez(l+1)
                     endif
                     call syncthreads
                     
                  enddo
                  
                  
                  !If my block index does not match the index of the internal-node block (lii), it means my thread is on the outer. I must exit
	              if(myblock .ne. (blockIdx%x+blockIdx%y*nxblock_d+blockIdx%z*nxyblock_d+1) )return
	                 

                  
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
       auxfields_s,locauxfields_s,forces_s &
       ,press_in,u_in,v_in,w_in,pxx_in,pyy_in,pzz_in,pxy_in,pxz_in,pyz_in &
       ,press_out,u_out,v_out,w_out,pxx_out,pyy_out,pzz_out,pxy_out,pxz_out,pyz_out)

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
      real(kind=db), dimension(ntotphifields) :: phifields_s
#ifdef MONOD
      real(kind=db) :: mu_max,Ks
#endif
#endif           
      real(kind=db) :: visc1,omega,fx,fy,fz
      
      
      
      real(kind=db), dimension(ntotauxfields) :: auxfields_s
      real(kind=db), dimension(ntotlocauxfields) :: locauxfields_s
      real(kind=db), dimension(ntotforces) :: forces_s
      
      real(kind=db), dimension(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff) :: &
       press_in,u_in,v_in,w_in,pxx_in,pyy_in,pzz_in,pxy_in,pxz_in,pyz_in, &
       press_out,u_out,v_out,w_out,pxx_out,pyy_out,pzz_out,pxy_out,pxz_out,pyz_out
      
      real(kind=db), shared :: f1(0:TILE_DIMx+1,0:TILE_DIMy+1,0:TILE_DIMz+1)
      
      real(kind=db) :: press,u,v,w,pxx,pyy,pzz,pxy,pxz,pyz
      real(kind=db) :: opress,ou,ov,ow,opxx,opyy,opzz,opxy,opxz,opyz,feq,fneq1,f_discr
      real(kind=db) :: mytemp,forcex,forcey,forcez,rhophi_loc,uu,udotc
      
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
                  phi_loc=phifields_s(idx5d(ii,jj,kk,1,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields))
#endif                  
#ifdef DENSRATIO
                  
                  rhophi_loc = rho_r*phi_loc+(ONE-phi_loc)*rho_b 
#else
                  rhophi_loc = 1.0_db !press_loc
#endif	

				  forcex=forces_s(idx5d(ii,jj,kk,1,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nforces))/rhophi_loc
				  forcey=forces_s(idx5d(ii,jj,kk,2,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nforces))/rhophi_loc
				  forcez=forces_s(idx5d(ii,jj,kk,3,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nforces))/rhophi_loc
                  
                  pxx=pxx_in(i,j,k)
                  pyy=pyy_in(i,j,k)
                  pzz=pzz_in(i,j,k)
                  pxy=pxy_in(i,j,k)
                  pxz=pxz_in(i,j,k)
                  pyz=pyz_in(i,j,k)
                  
                  

                  

!				  uu=HALF*(u_in(i,j,k)*u_in(i,j,k)+v_in(i,j,k)*v_in(i,j,k)+w_in(i,j,k)*w_in(i,j,k))*invcssq
				  
!                  do lii=1,nlinks
!                     udotc=(u_in(i,j,k)*dex(lii) + v_in(i,j,k)*dey(lii)+ w_in(i,j,k)*dez(lii))*invcssq
!		     feq=p(lii)*(press_in(i,j,k) + (udotc+0.5_db*udotc*udotc - uu))
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

                  pxx=pxx - cssq*press_in(i,j,k) - u_in(i,j,k)*u_in(i,j,k) 
                  pyy=pyy - cssq*press_in(i,j,k) - v_in(i,j,k)*v_in(i,j,k) 
                  pzz=pzz - cssq*press_in(i,j,k) - w_in(i,j,k)*w_in(i,j,k) 
                  pxy=pxy - u_in(i,j,k)*v_in(i,j,k)
                  pxz=pxz - u_in(i,j,k)*w_in(i,j,k)
                  pyz=pyz - v_in(i,j,k)*w_in(i,j,k)

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
                  omega_loc=(visc_loc*invcssq + 0.5_db) !it is tau   !visc_loc it is used to store the local viscosity
                  omega_loc=1.0_db/omega_loc !it is omega
#else
                  omega_loc=omega
#endif
#endif      
      
                  !press_out(i,j,k)=ZERO
                  u_out(i,j,k)=ZERO
                  v_out(i,j,k)=ZERO
                  w_out(i,j,k)=ZERO
                  pxx_out(i,j,k)=ZERO
                  pyy_out(i,j,k)=ZERO
                  pzz_out(i,j,k)=ZERO
                  pxy_out(i,j,k)=ZERO
                  pxz_out(i,j,k)=ZERO
                  pyz_out(i,j,k)=ZERO
!!!!!!!!!!!!!!!!!!!!!!!!!!0
                  uu=HALF*(u_in(i,j,k)*u_in(i,j,k)+v_in(i,j,k)*v_in(i,j,k)+w_in(i,j,k)*w_in(i,j,k))*invcssq

			      feq=p(0)*(press_in(i,j,k) - uu)
				  fneq1=(HALF*invcssq)*(-pxx-pyy-pzz)
				  F_discr = p(0)*(- u_in(i,j,k)*forcex - v_in(i,j,k)*forcey - w_in(i,j,k)*forcez)*invcssq
                  
                  press_out(i,j,k)=feq + (1.0_db-omega_loc)*fneq1*p(0) + HALF*(F_discr)
                  
                  do l=1,nlinks
                     udotc=(u_in(i,j,k)*dex(l) + v_in(i,j,k)*dey(l)+ w_in(i,j,k)*dez(l))*invcssq
		             feq=p(l)*(press_in(i,j,k) + (udotc+0.5_db*udotc*udotc - uu))
                     fneq1=(HALF/(cssq*cssq))*( (dex(l)*dex(l)-cssq)*pxx &
		              + (dey(l)*dey(l)-cssq)*pyy + (dez(l)*dez(l)-cssq)*pzz &
	                  + TWO*(dex(l)*dey(l))*pxy + TWO*(dex(l)*dez(l))*pxz &
		              + TWO*(dey(l)*dez(l))*pyz)
		             F_discr = p(l)*(((dex(l) - u_in(i,j,k)) + udotc * dex(l))*forcex &
                      + ((dey(l) - v_in(i,j,k)) + udotc * dey(l))*forcey &
                      + ((dez(l) - w_in(i,j,k)) + udotc * dez(l))*forcez)*invcssq
                     lii=li+ex(l)
                     ljj=lj+ey(l)
                     lkk=lk+ez(l)
                     lii=mod(lii+TILE_DIMx+2,(TILE_DIMx+2))
                     ljj=mod(ljj+TILE_DIMy+2,(TILE_DIMy+2))
                     lkk=mod(lkk+TILE_DIMz+2,(TILE_DIMz+2)) 
		             f1(lii,ljj,lkk)=feq + (ONE-omega_loc)*p(l)*fneq1 + HALF*F_discr
!		             if(gi==iprobe .and. gj==jprobe .and. gk==kprobe .and. myblock==intblock)write(*,*)l,f1(lii,ljj,lkk)
!		             if(gi==iprobe .and. gj==jprobe .and. gk==kprobe .and. myblock==intblock)write(*,*)l,feq
!		             if(gi==iprobe .and. gj==jprobe .and. gk==kprobe .and. myblock==intblock)write(*,*)l,(ONE-omega_loc)*p(l)*fneq1
!		             if(gi==iprobe .and. gj==jprobe .and. gk==kprobe .and. myblock==intblock)write(*,*)l,HALF*F_discr
		             call syncthreads
		             if(myblock .eq. intblock)then
		             press_out(i,j,k)=press_out(i,j,k) + f1(li,lj,lk)
		             u_out(i,j,k)=u_out(i,j,k) + f1(li,lj,lk)*dex(l)
                     v_out(i,j,k)=v_out(i,j,k) + f1(li,lj,lk)*dey(l)
                     w_out(i,j,k)=w_out(i,j,k) + f1(li,lj,lk)*dez(l)
                     pxx_out(i,j,k)=pxx_out(i,j,k) + f1(li,lj,lk)*dex(l)*dex(l)
                     pyy_out(i,j,k)=pyy_out(i,j,k) + f1(li,lj,lk)*dey(l)*dey(l)
                     pzz_out(i,j,k)=pzz_out(i,j,k) + f1(li,lj,lk)*dez(l)*dez(l)
                     pxy_out(i,j,k)=pxy_out(i,j,k) + f1(li,lj,lk)*dex(l)*dey(l)
                     pxz_out(i,j,k)=pxz_out(i,j,k) + f1(li,lj,lk)*dex(l)*dez(l)
                     pyz_out(i,j,k)=pyz_out(i,j,k) + f1(li,lj,lk)*dey(l)*dez(l)
                     endif
                     call syncthreads
                     
                  enddo
                  
                  
                  !If my block index does not match the index of the internal-node block (lii), it means my thread is on the outer. I must exit
	              if(myblock .ne. intblock)return
	                 

                   
    endsubroutine fused_LB_kernel1

endmodule
