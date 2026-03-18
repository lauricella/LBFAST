#include "defines.h"
#if !defined(_OPENACC)  
#error "To use this module the macros _OPENACC should be defined."
#endif

module lb_cuda_boundary

   use vars
   use iso_c_binding
   use cudafor
   use mpi_template, only: coords,dostop,doerror,mydev,myrank,nprocs,nbuff,nbuffbvec
   use lb_cuda_vars

   implicit none
   

contains
   
   attributes(global) subroutine LB_int_boundary_kernel(flip,flop,nx,ny,nz,coords,isfluid &  
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
       ,hfields_in,hfields_out &
#ifdef TWOCOMPONENT 
       ,auxfields_s,locauxfields_s &
#endif  
       ,forces_s)

      implicit none
      
      integer :: flip,flop,nx,ny,nz
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
      

      real(kind=db) :: press,u,v,w,pxx,pyy,pzz,pxy,pxz,pyz
      real(kind=db) :: opress,ou,ov,ow,opxx,opyy,opzz,opxy,opxz,opyz
      real(kind=db) :: mytemp,rhophi_loc,press_loc
      !real(kind=db) :: forcex,forcey,forcez,F_discr
  
      real(kind=db) :: fneq1,feq,fpost,uu,udotc
#ifdef TWOCOMPONENT
      real(kind=db) :: wet_loc
#endif
      real(kind=db) :: omega_loc,phi_loc,visc_loc
#ifdef SMAGORINSKI
	  real(kind=db) :: QQ
#endif

      integer :: i,j,k,l,lopp
      integer :: gi,gj,gk
      integer :: myblock,ii,jj,kk
      integer :: iii,jjj,kkk
      
      i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z
      
      if(isfluid(i,j,k) .ne. -1)return
      
      gi=nx*coords(1)+i
      gj=ny*coords(2)+j
      gk=nz*coords(3)+k
      
      myblock=blockIdx%x+blockIdx%y*nxblock_d+blockIdx%z*nxyblock_d+1
      
      ii=threadIdx%x
      jj=threadIdx%y
      kk=threadIdx%z
       
#ifdef TWOCOMPONENT	       
      phi_loc=real(phifields_s(ii,jj,kk,1,myblock),kind=db)
#endif	
#ifdef DENSRATIO
	  rhophi_loc = rho_r*phi_loc+(ONE-phi_loc)*rho_b 
#else
	  rhophi_loc = 1.0_db !press_loc
#endif	

!	  forcex=force_s(ii,jj,kk,1,myblock)
!	  forcey=force_s(ii,jj,kk,2,myblock)
!	  forcez=force_s(ii,jj,kk,3,myblock)


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
      
      opress=real(hfields_out(ii,jj,kk,1,myblock),kind=db)
      ou=real(hfields_out(ii,jj,kk,2,myblock),kind=db)
      ov=real(hfields_out(ii,jj,kk,3,myblock),kind=db)
      ow=real(hfields_out(ii,jj,kk,4,myblock),kind=db)
      opxx=real(hfields_out(ii,jj,kk,5,myblock),kind=db)
      opyy=real(hfields_out(ii,jj,kk,6,myblock),kind=db)
      opzz=real(hfields_out(ii,jj,kk,7,myblock),kind=db)
      opxy=real(hfields_out(ii,jj,kk,8,myblock),kind=db)
      opxz=real(hfields_out(ii,jj,kk,9,myblock),kind=db)
      opyz=real(hfields_out(ii,jj,kk,10,myblock),kind=db)
	  
	  pxx=pxx - u*u 
	  pyy=pyy - v*v 
	  pzz=pzz - w*w 
	  pxy=pxy - u*v
	  pxz=pxz - u*w
	  pyz=pyz - v*w

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

      

      do l=1,nlinks
		  lopp=opp(l)
		  iii=i+ex(lopp)
		  jjj=j+ey(lopp)
		  kkk=k+ez(lopp)
		  if(isfluid(iii,jjj,kkk).ne.0) cycle 
		  feq=p(l)*(press)
		  fneq1=p(l)*(HALF/(cssq*cssq))*( (dex(l)*dex(l)-cssq)*pxx &
		   + (dey(l)*dey(l)-cssq)*pyy + (dez(l)*dez(l)-cssq)*pzz &
	       + TWO*(dex(l)*dey(l))*pxy + TWO*(dex(l)*dez(l))*pxz &
		   + TWO*(dey(l)*dez(l))*pyz)
          ! F_discr = p(l)*(dex(l)*forcex &
           ! + dey(l)*forcey &
           ! + dez(l)*forcez)*invcssq
		  fpost=feq + (ONE-omega_loc)*fneq1 !+ HALF*(F_discr)	
		  opress=opress + fpost
		  ou=ou + fpost*dex(l)
		  ov=ov + fpost*dey(l)
		  ow=ow + fpost*dez(l)
		  opxx=opxx + fpost*(dex(l)*dex(l)-cssq)
          opyy=opyy + fpost*(dey(l)*dey(l)-cssq)
          opzz=opzz + fpost*(dez(l)*dez(l)-cssq)
          opxy=opxy + fpost*dex(l)*dey(l)
          opxz=opxz + fpost*dex(l)*dez(l)
          opyz=opyz + fpost*dey(l)*dez(l)		  
      enddo
      

		 
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
	              
     return

   endsubroutine LB_int_boundary_kernel   
#ifdef TWOCOMPONENT	    
   attributes(global) subroutine PHI_int_boundary_kernel(flop,nx,ny,nz,coords,isfluid &
    ,visc2,rho_r,rho_b,invrho_r,invrho_b,sharp_c,beta,kapp,tau_diff,sigma &
#ifdef WETTABILITY  
    ,wettab_r,wettab_b
#endif       
#ifdef MONOD	
    ,mu_max,Ks &
#endif
    ,ntothfields,ntotphifields,ntotauxfields,ntotlocauxfields,ntotforces &
       ,hfields_s,phifields_s,auxfields_s,locauxfields_s,forces_s)

      implicit none
      integer :: flop,nx,ny,nz
      integer, dimension(3) :: coords
      integer(kind=isf), dimension(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff) :: isfluid
     
#ifdef WETTABILITY  
      real(kind=db) :: wettab_r,wettab_b  
#endif   
      real(kind=db) :: visc2,rho_r,rho_b,invrho_r,invrho_b,sharp_c,beta,kapp,tau_diff,sigma  
#ifdef MONOD
      real(kind=db) :: mu_max,Ks
#endif
      integer :: ntothfields,ntotphifields,ntotauxfields,ntotlocauxfields,ntotforces
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields,nblocks_d) :: hfields_s
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields,nblocks_d) :: phifields_s
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nauxfields,nblocks_d) :: auxfields_s
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nlocauxfields,nblocks_d) :: locauxfields_s
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nforces,nblocks_d) :: forces_s
      
      integer :: i,j,k,l
      integer :: gi,gj,gk
      integer :: myblock,ii,jj,kk
      integer :: iii,jjj,kkk
      real(kind=db) :: gradfix,gradfiy,gradfiz,grad_parallel,modgrad,phi_fluid,phi_ghost
      real(kind=db) :: loc_u,loc_v,loc_w,loc_phi,theta_rad,cot_theta,dphi_dz
      integer :: oii,ojj,okk
      integer :: oxblock,oyblock,ozblock,omyblock
      real(kind=db) :: wettab_r_sub=90.0_db
    
      i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z
      
      if(isfluid(i,j,k) .ne. 0)return
      
      gi=nx*coords(1)+i
      gj=ny*coords(2)+j
      gk=nz*coords(3)+k
      
      myblock=blockIdx%x+blockIdx%y*nxblock_d+blockIdx%z*nxyblock_d+1
      
      ii=threadIdx%x
      jj=threadIdx%y
      kk=threadIdx%z
	   
#ifdef WETTABILITY       
      wettab_r_sub=wettab_r
#endif                    
      
	  do l = 1, nlinks
		iii = i + ex(l)
		jjj = j + ey(l)
		kkk = k + ez(l)

		if (isfluid(iii,jjj,kkk) .ne. -1) cycle  ! only fluid neighbor
		
		oxblock=(iii+2*TILE_DIMx-1)/TILE_DIMx   
		oyblock=(jjj+2*TILE_DIMy-1)/TILE_DIMy     
		ozblock=(kkk+2*TILE_DIMz-1)/TILE_DIMz 
		omyblock=(oxblock-1)+(oyblock-1)*nxblock_d+(ozblock-1)*nxyblock_d+1
		oii=iii-oxblock*TILE_DIMx+2*TILE_DIMx
		ojj=jjj-oyblock*TILE_DIMy+2*TILE_DIMy
		okk=kkk-ozblock*TILE_DIMz+2*TILE_DIMz

		! Found fluid neighbor: enforce contact angle via ghost node extrapolation
		phi_fluid = real(phifields_s(oii,ojj,okk,1,omyblock),kind=db)
		
		
		! Estimate gradient parallel to wall
		modgrad=real(auxfields_s(oii,ojj,okk,4,omyblock),kind=db) !modgrad
		gradfix=real(auxfields_s(oii,ojj,okk,1,omyblock),kind=db)*modgrad !normx*modgrad
		gradfiy=real(auxfields_s(oii,ojj,okk,2,omyblock),kind=db)*modgrad !normy*modgrad
		gradfiz=real(auxfields_s(oii,ojj,okk,3,omyblock),kind=db)*modgrad !normz*modgrad      
      
        grad_parallel=ZERO
		if(l.eq.1 .or. l.eq.2)then
			grad_parallel = sqrt(gradfiy**2 + gradfiz**2)
		elseif(l.eq.3 .or. l.eq.4)then
			grad_parallel = sqrt(gradfix**2 + gradfiz**2)
		elseif(l.eq.5 .or. l.eq.6)then
			grad_parallel = sqrt(gradfix**2 + gradfiy**2)
		endif
		
		! Contact angle correction
		theta_rad = (180.0_db-wettab_r_sub) * pi_greek / 180.0_db
		cot_theta = 1.0_db / tan(theta_rad)

		dphi_dz = - grad_parallel * cot_theta 

		  

		phi_ghost = phi_fluid + dphi_dz  ! extrapolate from fluid node

		! Clamp to [0,1]
		loc_phi = max(0.0_db, min(1.0_db, phi_ghost))
		
		exit
	  end do

      phifields_s(ii,jj,kk,1,myblock)=real(loc_phi,kind=strdb)
      
      return
      
   endsubroutine PHI_int_boundary_kernel
   
   attributes(global) subroutine phi_sum_count_kernel(flop,nx,ny,nz,coords,isfluid &
    ,visc2,rho_r,rho_b,invrho_r,invrho_b,sharp_c,beta,kapp,tau_diff,sigma &
    ,ntothfields,ntotphifields,ntotauxfields,ntotlocauxfields, &
       hfields_s,phifields_s,auxfields_s,locauxfields_s,loc_phi_sum,cnt)

      implicit none
      integer :: flop,nx,ny,nz
      integer, dimension(3) :: coords
      integer(kind=isf), dimension(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff) :: isfluid
      
      real(kind=db) :: visc2,rho_r,rho_b,invrho_r,invrho_b,sharp_c,beta,kapp,tau_diff,sigma  

      integer :: ntothfields,ntotphifields,ntotauxfields,ntotlocauxfields
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields,nblocks_d) :: hfields_s
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields,nblocks_d) :: phifields_s
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nauxfields,nblocks_d) :: auxfields_s
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nlocauxfields,nblocks_d) :: locauxfields_s
      real(kind=db) :: loc_phi_sum
      integer :: cnt
      
      integer :: i,j,k,l
      integer :: gi,gj,gk
      integer :: myblock,ii,jj,kk
      real(kind=db) :: gradfix,gradfiy,gradfiz,grad_parallel,modgrad,phi_fluid
      real(kind=db) :: loc_u,loc_v,loc_w,loc_phi,theta_rad,cot_theta,dphi_dz
      real(kind=db) :: dummy
      integer :: dummy_i
    
      i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z
      
      if(abs(isfluid(i,j,k)) /= 1)return
      
      gi=nx*coords(1)+i
      gj=ny*coords(2)+j
      gk=nz*coords(3)+k
      
      myblock=blockIdx%x+blockIdx%y*nxblock_d+blockIdx%z*nxyblock_d+1
      
      ii=threadIdx%x
      jj=threadIdx%y
      kk=threadIdx%z

      loc_phi=real(phifields_s(ii,jj,kk,1,myblock),kind=db)

      dummy=atomicAdd(loc_phi_sum, loc_phi)

      if(loc_phi > 0.5d0 .and. loc_phi < 0.9d0)then
        dummy_i=atomicAdd(cnt, 1)
      endif
    
      return
    
  end subroutine phi_sum_count_kernel

   attributes(global) subroutine apply_lagrangian_phi_kernel(flop,nx,ny,nz,coords,isfluid &
    ,visc2,rho_r,rho_b,invrho_r,invrho_b,sharp_c,beta,kapp,tau_diff,sigma &
    ,ntothfields,ntotphifields,ntotauxfields,ntotlocauxfields, &
       hfields_s,phifields_s,auxfields_s,locauxfields_s,loc_corr)

      implicit none
      integer :: flop,nx,ny,nz
      integer, dimension(3) :: coords
      integer(kind=isf), dimension(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff) :: isfluid
      
      real(kind=db) :: visc2,rho_r,rho_b,invrho_r,invrho_b,sharp_c,beta,kapp,tau_diff,sigma  

      integer :: ntothfields,ntotphifields,ntotauxfields,ntotlocauxfields
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields,nblocks_d) :: hfields_s
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields,nblocks_d) :: phifields_s
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nauxfields,nblocks_d) :: auxfields_s
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nlocauxfields,nblocks_d) :: locauxfields_s
      real(kind=db) :: loc_corr
      
      integer :: i,j,k,l
      integer :: gi,gj,gk
      integer :: myblock,ii,jj,kk
      real(kind=db) :: gradfix,gradfiy,gradfiz,grad_parallel,modgrad,phi_fluid
      real(kind=db) :: loc_u,loc_v,loc_w,loc_phi,theta_rad,cot_theta,dphi_dz
    
      i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z
      
      if(abs(isfluid(i,j,k)) /= 1)return
      
      gi=nx*coords(1)+i
      gj=ny*coords(2)+j
      gk=nz*coords(3)+k
      
      myblock=blockIdx%x+blockIdx%y*nxblock_d+blockIdx%z*nxyblock_d+1
      
      ii=threadIdx%x
      jj=threadIdx%y
      kk=threadIdx%z

      loc_phi=real(phifields_s(ii,jj,kk,1,myblock),kind=db)
      if(loc_phi > 0.5d0 .and. loc_phi < 0.9d0)then
        phifields_s(ii,jj,kk,1,myblock)=real(loc_phi + loc_corr ,kind=strdb)

      endif
    
      return
    
  end subroutine apply_lagrangian_phi_kernel
#endif
endmodule lb_cuda_boundary
