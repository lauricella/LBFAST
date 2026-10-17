#include "defines.h"
#if !defined(_OPENACC)  
#error "To use this module the macros _OPENACC should be defined."
#endif

! Phase-field update after the collision: one explicit time step of the conservative Allen-Cahn
! equation for phi, Eq. (22) of Lauricella et al., Phys. Fluids 37, 072111 (2025).
! update_phifields_kernel covers all the tiles of the subdomain. With ASYNCMPI, update_phifields_ext
! first launches the six face kernels (tiles that touch the subdomain boundary), then the halo
! exchange of phi is started and update_phifields_kernel_int updates the interior tiles before it
! is completed. The variants have the same body and cover the same tiles as the corresponding
! kernels of lb_cuda_moments_module.f90, where the tiled storage is described.
module lb_cuda_update_phi

   use vars
   use iso_c_binding
   use cudafor
   use mpi_template, only: coords,dostop,doerror,mydev,myrank,nprocs,nbuff,nbuffbvec
   use lb_cuda_vars

   implicit none
   

contains
   ! the phase field exists only in two-component flows
#ifdef TWOCOMPONENT   
   ! update_phifields_kernel: new phase field for all the fluid nodes of the subdomain (dimGrid: one
   ! block per tile, one thread per node).
   ! Reads phifields_in(1) = phi^n, hfields_s(2:4) (velocity stored by moments_LB_kernel in this time
   ! step, including the half-force shift), auxfields_s(1:4) and locauxfields_s(1:2), computed from
   ! phi^n earlier in the time step. Writes phifields_out(1) = phi^(n+1).
   ! tau_diff is the interface diffusivity and sharp_c = 4 tau_diff/width the coefficient of the
   ! sharpening flux. visc2, rho_r, rho_b, invrho_r, invrho_b, beta, kapp, sigma and forces_s are
   ! not used.
   attributes(global) subroutine update_phifields_kernel(step,iprobe,jprobe,kprobe,flop,nx,ny,nz,coords,isfluid &
#ifdef TWOCOMPONENT    
    ,visc2,rho_r,rho_b,invrho_r,invrho_b,sharp_c,beta,kapp,tau_diff,sigma &
#ifdef MONOD	
    ,mu_max,Ks &
#endif
#endif 
    ,ntothfields,ntotphifields,ntotauxfields,ntotlocauxfields,ntotforces &
       ,hfields_s,phifields_in,phifields_out,auxfields_s,locauxfields_s,forces_s)

      implicit none
      integer :: step,iprobe,jprobe,kprobe,flop,nx,ny,nz
      integer, dimension(3) :: coords
      integer(kind=isf), dimension(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff) :: isfluid
     
 
#ifdef TWOCOMPONENT 
      real(kind=db) :: visc2,rho_r,rho_b,invrho_r,invrho_b,sharp_c,beta, &
       kapp,tau_diff,sigma,modgrad 
#ifdef MONOD
      real(kind=db) :: mu_max,Ks
#endif
#endif 
      integer :: ntothfields,ntotphifields,ntotauxfields,ntotlocauxfields,ntotforces
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields,nblocks_d) :: hfields_s
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields,nblocks_d) :: phifields_in,phifields_out
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nauxfields,nblocks_d) :: auxfields_s
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nlocauxfields,nblocks_d) :: locauxfields_s
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nforces,nblocks_d) :: forces_s
      
      integer :: i,j,k
      !integer :: gi,gj,gk
      integer :: myblock,ii,jj,kk
      real(kind=db) :: gradfix,gradfiy,gradfiz,phi_loc,phi_out,mytemp
      real(kind=db) :: loc_u,loc_v,loc_w,lap_phi_loc
#ifdef MONOD
	  real(kind=db) :: S_mono
#endif
      
      ! one block per tile (tile number = blockIdx) and one thread per node of the tile
      i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z
      
      if(abs(isfluid(i,j,k)) .ne. 1)return
      
      !gi=nx*coords(1)+i
      !gj=ny*coords(2)+j
      !gk=nz*coords(3)+k
      
      ! storage tile: block coordinates equal the tile numbers, since coordinate 0 is the halo tile
      myblock=blockIdx%x+blockIdx%y*nxblock_d+blockIdx%z*nxyblock_d+1
      
      ii=threadIdx%x
      jj=threadIdx%y
      kk=threadIdx%z
	   
#ifdef TWOCOMPONENT		  
	  !reuse gradrhox,gradrhoy,gradrhoz as local velocity (reusing variables is saving register memory)
	  !reuse gradfix,gradfiy,gradfiz
          ! sharpening term -sharp_c div(theta), theta = phi(1-phi) n (locauxfields 2)
	  mytemp= -sharp_c*real(locauxfields_s(ii,jj,kk,2,myblock),kind=db)
	  
      ! grad(phi) = n |grad(phi)| (auxfields 1-3 and 4)
      modgrad=real(auxfields_s(ii,jj,kk,4,myblock),kind=db) !modgrad
	  gradfix=real(auxfields_s(ii,jj,kk,1,myblock),kind=db)*modgrad !normx*modgrad
	  gradfiy=real(auxfields_s(ii,jj,kk,2,myblock),kind=db)*modgrad !normy*modgrad
	  gradfiz=real(auxfields_s(ii,jj,kk,3,myblock),kind=db)*modgrad !normz*modgrad
                  
      phi_loc=real(phifields_in(ii,jj,kk,1,myblock),kind=db)
      lap_phi_loc=real(locauxfields_s(ii,jj,kk,1,myblock),kind=db)
                  
      ! velocity stored by moments_LB_kernel in this time step, including the half-force shift
      loc_u=real(hfields_s(ii,jj,kk,2,myblock),kind=db) !velocity
      loc_v=real(hfields_s(ii,jj,kk,3,myblock),kind=db)
      loc_w=real(hfields_s(ii,jj,kk,4,myblock),kind=db)
                
      ! phi^(n+1) = phi^n - u.grad(phi) + tau_diff lap(phi) - sharp_c div(theta), with time step 1
      phi_out = phi_loc &
        - loc_u*gradfix - loc_v*gradfiy &
        - loc_w*gradfiz + tau_diff*lap_phi_loc + mytemp
#endif	

      ! MONOD: source term S_mono, a growth term with the Monod parameters mu_max and Ks
#ifdef MONOD
      S_mono = mu_max * phi_loc/(Ks + phi_loc) * phi_loc * (1.0_db - phi_loc)
      phi_out=phi_out + S_mono
	  !phi_out = min(1.0_db, max(0.0_db, phi_out))		 
#endif
      phifields_out(ii,jj,kk,1,myblock)=real(phi_out,kind=strdb)
      return
      
   endsubroutine update_phifields_kernel
   
      ! update_phifields_kernel_int: as update_phifields_kernel, for the interior tiles (tiles
      ! 2..n/TILE_DIM-1 in each direction, dimGridInt); block b handles tile t = b+1, which has block
      ! coordinate t (coordinate 0 is the halo tile).
      attributes(global) subroutine update_phifields_kernel_int(step,iprobe,jprobe,kprobe,flop,nx,ny,nz,coords,isfluid &
#ifdef TWOCOMPONENT    
    ,visc2,rho_r,rho_b,invrho_r,invrho_b,sharp_c,beta,kapp,tau_diff,sigma &
#ifdef MONOD	
    ,mu_max,Ks &
#endif
#endif 
    ,ntothfields,ntotphifields,ntotauxfields,ntotlocauxfields,ntotforces &
       ,hfields_s,phifields_in,phifields_out,auxfields_s,locauxfields_s,forces_s)

      implicit none
      integer :: step,iprobe,jprobe,kprobe,flop,nx,ny,nz
      integer, dimension(3) :: coords
      integer(kind=isf), dimension(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff) :: isfluid
     
 
#ifdef TWOCOMPONENT 
      real(kind=db) :: visc2,rho_r,rho_b,invrho_r,invrho_b,sharp_c,beta, &
       kapp,tau_diff,sigma,modgrad 
#ifdef MONOD
      real(kind=db) :: mu_max,Ks
#endif
#endif 
      integer :: ntothfields,ntotphifields,ntotauxfields,ntotlocauxfields,ntotforces
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields,nblocks_d) :: hfields_s
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields,nblocks_d) :: phifields_in,phifields_out
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nauxfields,nblocks_d) :: auxfields_s
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nlocauxfields,nblocks_d) :: locauxfields_s
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nforces,nblocks_d) :: forces_s
      
      integer :: i,j,k
      !integer :: gi,gj,gk
      integer :: myblock,ii,jj,kk
      real(kind=db) :: gradfix,gradfiy,gradfiz,phi_loc,phi_out,mytemp
      real(kind=db) :: loc_u,loc_v,loc_w,lap_phi_loc
#ifdef MONOD
	  real(kind=db) :: S_mono
#endif
      
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
      
      
	   
#ifdef TWOCOMPONENT
				  
	  mytemp= -sharp_c*real(locauxfields_s(ii,jj,kk,2,myblock),kind=db)
	  
      modgrad=real(auxfields_s(ii,jj,kk,4,myblock),kind=db) !modgrad
	  gradfix=real(auxfields_s(ii,jj,kk,1,myblock),kind=db)*modgrad !normx*modgrad
	  gradfiy=real(auxfields_s(ii,jj,kk,2,myblock),kind=db)*modgrad !normy*modgrad
	  gradfiz=real(auxfields_s(ii,jj,kk,3,myblock),kind=db)*modgrad !normz*modgrad
                  
      phi_loc=real(phifields_in(ii,jj,kk,1,myblock),kind=db)
      lap_phi_loc=real(locauxfields_s(ii,jj,kk,1,myblock),kind=db)
                  
      loc_u=real(hfields_s(ii,jj,kk,2,myblock),kind=db) !velocity
      loc_v=real(hfields_s(ii,jj,kk,3,myblock),kind=db)
      loc_w=real(hfields_s(ii,jj,kk,4,myblock),kind=db)

                  
      phi_out = phi_loc &
        - loc_u*gradfix - loc_v*gradfiy &
        - loc_w*gradfiz + tau_diff*lap_phi_loc + mytemp
#endif	

#ifdef MONOD
      S_mono = mu_max * phi_loc/(Ks + phi_loc) * phi_loc * (1.0_db - phi_loc)
      phi_out=phi_out + S_mono
	  !phi_out = min(1.0_db, max(0.0_db, phi_out))		 
#endif
      phifields_out(ii,jj,kk,1,myblock)=real(phi_out,kind=strdb)
      
      return
      
   endsubroutine update_phifields_kernel_int
   
   ! update_phifields_kernel_xplus: as update_phifields_kernel, for the last tile layer along x
   ! (bx = nxblock_d-2), y and z tiles 2..n/TILE_DIM-1 (dimGridx); not launched when there is a
   ! single tile along x.
   attributes(global) subroutine update_phifields_kernel_xplus(step,iprobe,jprobe,kprobe,flop,nx,ny,nz,coords,isfluid &
#ifdef TWOCOMPONENT    
    ,visc2,rho_r,rho_b,invrho_r,invrho_b,sharp_c,beta,kapp,tau_diff,sigma &
#ifdef MONOD	
    ,mu_max,Ks &
#endif
#endif 
    ,ntothfields,ntotphifields,ntotauxfields,ntotlocauxfields,ntotforces &
       ,hfields_s,phifields_in,phifields_out,auxfields_s,locauxfields_s,forces_s)

      implicit none
      integer :: step,iprobe,jprobe,kprobe,flop,nx,ny,nz
      integer, dimension(3) :: coords
      integer(kind=isf), dimension(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff) :: isfluid
     
 
#ifdef TWOCOMPONENT 
      real(kind=db) :: visc2,rho_r,rho_b,invrho_r,invrho_b,sharp_c,beta, &
       kapp,tau_diff,sigma,modgrad 
#ifdef MONOD
      real(kind=db) :: mu_max,Ks
#endif
#endif 
      integer :: ntothfields,ntotphifields,ntotauxfields,ntotlocauxfields,ntotforces
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields,nblocks_d) :: hfields_s
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields,nblocks_d) :: phifields_in,phifields_out
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nauxfields,nblocks_d) :: auxfields_s
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nlocauxfields,nblocks_d) :: locauxfields_s
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nforces,nblocks_d) :: forces_s
      
      integer :: i,j,k
      !integer :: gi,gj,gk
      integer :: myblock,ii,jj,kk
      real(kind=db) :: gradfix,gradfiy,gradfiz,phi_loc,phi_out,mytemp
      real(kind=db) :: loc_u,loc_v,loc_w,lap_phi_loc
#ifdef MONOD
	  real(kind=db) :: S_mono
#endif
      
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
      
      
	   
#ifdef TWOCOMPONENT
				  
	  mytemp= -sharp_c*real(locauxfields_s(ii,jj,kk,2,myblock),kind=db)
	  
      modgrad=real(auxfields_s(ii,jj,kk,4,myblock),kind=db) !modgrad
	  gradfix=real(auxfields_s(ii,jj,kk,1,myblock),kind=db)*modgrad !normx*modgrad
	  gradfiy=real(auxfields_s(ii,jj,kk,2,myblock),kind=db)*modgrad !normy*modgrad
	  gradfiz=real(auxfields_s(ii,jj,kk,3,myblock),kind=db)*modgrad !normz*modgrad
                  
      phi_loc=real(phifields_in(ii,jj,kk,1,myblock),kind=db)
      lap_phi_loc=real(locauxfields_s(ii,jj,kk,1,myblock),kind=db)
                  
      loc_u=real(hfields_s(ii,jj,kk,2,myblock),kind=db) !velocity
      loc_v=real(hfields_s(ii,jj,kk,3,myblock),kind=db)
      loc_w=real(hfields_s(ii,jj,kk,4,myblock),kind=db)

                  
      phi_out = phi_loc &
        - loc_u*gradfix - loc_v*gradfiy &
        - loc_w*gradfiz + tau_diff*lap_phi_loc + mytemp
#endif	

#ifdef MONOD
      S_mono = mu_max * phi_loc/(Ks + phi_loc) * phi_loc * (1.0_db - phi_loc)
      phi_out=phi_out + S_mono
	  !phi_out = min(1.0_db, max(0.0_db, phi_out))		 
#endif
      phifields_out(ii,jj,kk,1,myblock)=real(phi_out,kind=strdb)
      
      return
      
   endsubroutine update_phifields_kernel_xplus
   
   ! update_phifields_kernel_xminus: as update_phifields_kernel, for the first tile layer along x
   ! (bx = 1), y and z tiles 2..n/TILE_DIM-1 (dimGridx).
   attributes(global) subroutine update_phifields_kernel_xminus(step,iprobe,jprobe,kprobe,flop,nx,ny,nz,coords,isfluid &
#ifdef TWOCOMPONENT    
    ,visc2,rho_r,rho_b,invrho_r,invrho_b,sharp_c,beta,kapp,tau_diff,sigma &
#ifdef MONOD	
    ,mu_max,Ks &
#endif
#endif 
    ,ntothfields,ntotphifields,ntotauxfields,ntotlocauxfields,ntotforces &
       ,hfields_s,phifields_in,phifields_out,auxfields_s,locauxfields_s,forces_s)

      implicit none
      integer :: step,iprobe,jprobe,kprobe,flop,nx,ny,nz
      integer, dimension(3) :: coords
      integer(kind=isf), dimension(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff) :: isfluid
     
 
#ifdef TWOCOMPONENT 
      real(kind=db) :: visc2,rho_r,rho_b,invrho_r,invrho_b,sharp_c,beta, &
       kapp,tau_diff,sigma,modgrad 
#ifdef MONOD
      real(kind=db) :: mu_max,Ks
#endif
#endif 
      integer :: ntothfields,ntotphifields,ntotauxfields,ntotlocauxfields,ntotforces
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields,nblocks_d) :: hfields_s
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields,nblocks_d) :: phifields_in,phifields_out
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nauxfields,nblocks_d) :: auxfields_s
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nlocauxfields,nblocks_d) :: locauxfields_s
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nforces,nblocks_d) :: forces_s
      
      integer :: i,j,k
      !integer :: gi,gj,gk
      integer :: myblock,ii,jj,kk
      real(kind=db) :: gradfix,gradfiy,gradfiz,phi_loc,phi_out,mytemp
      real(kind=db) :: loc_u,loc_v,loc_w,lap_phi_loc
#ifdef MONOD
	  real(kind=db) :: S_mono
#endif
      
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
      
	   
#ifdef TWOCOMPONENT
				  
	  mytemp= -sharp_c*real(locauxfields_s(ii,jj,kk,2,myblock),kind=db)
	  
      modgrad=real(auxfields_s(ii,jj,kk,4,myblock),kind=db) !modgrad
	  gradfix=real(auxfields_s(ii,jj,kk,1,myblock),kind=db)*modgrad !normx*modgrad
	  gradfiy=real(auxfields_s(ii,jj,kk,2,myblock),kind=db)*modgrad !normy*modgrad
	  gradfiz=real(auxfields_s(ii,jj,kk,3,myblock),kind=db)*modgrad !normz*modgrad
                  
      phi_loc=real(phifields_in(ii,jj,kk,1,myblock),kind=db)
      lap_phi_loc=real(locauxfields_s(ii,jj,kk,1,myblock),kind=db)
                  
      loc_u=real(hfields_s(ii,jj,kk,2,myblock),kind=db) !velocity
      loc_v=real(hfields_s(ii,jj,kk,3,myblock),kind=db)
      loc_w=real(hfields_s(ii,jj,kk,4,myblock),kind=db)

                  
      phi_out = phi_loc &
        - loc_u*gradfix - loc_v*gradfiy &
        - loc_w*gradfiz + tau_diff*lap_phi_loc + mytemp
#endif	

#ifdef MONOD
      S_mono = mu_max * phi_loc/(Ks + phi_loc) * phi_loc * (1.0_db - phi_loc)
      phi_out=phi_out + S_mono
	  !phi_out = min(1.0_db, max(0.0_db, phi_out))		 
#endif
      phifields_out(ii,jj,kk,1,myblock)=real(phi_out,kind=strdb)
      
      return
      
   endsubroutine update_phifields_kernel_xminus
   
   ! update_phifields_kernel_yplus: as update_phifields_kernel, for the last tile layer along y
   ! (by = nyblock_d-2), all x tiles, z tiles 2..nz/TILE_DIMz-1 (dimGridy); not launched when
   ! there is a single tile along y.
   attributes(global) subroutine update_phifields_kernel_yplus(step,iprobe,jprobe,kprobe,flop,nx,ny,nz,coords,isfluid &
#ifdef TWOCOMPONENT    
    ,visc2,rho_r,rho_b,invrho_r,invrho_b,sharp_c,beta,kapp,tau_diff,sigma &
#ifdef MONOD	
    ,mu_max,Ks &
#endif
#endif 
    ,ntothfields,ntotphifields,ntotauxfields,ntotlocauxfields,ntotforces &
       ,hfields_s,phifields_in,phifields_out,auxfields_s,locauxfields_s,forces_s)

      implicit none
      integer :: step,iprobe,jprobe,kprobe,flop,nx,ny,nz
      integer, dimension(3) :: coords
      integer(kind=isf), dimension(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff) :: isfluid
     
 
#ifdef TWOCOMPONENT 
      real(kind=db) :: visc2,rho_r,rho_b,invrho_r,invrho_b,sharp_c,beta, &
       kapp,tau_diff,sigma,modgrad 
#ifdef MONOD
      real(kind=db) :: mu_max,Ks
#endif
#endif 
      integer :: ntothfields,ntotphifields,ntotauxfields,ntotlocauxfields,ntotforces
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields,nblocks_d) :: hfields_s
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields,nblocks_d) :: phifields_in,phifields_out
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nauxfields,nblocks_d) :: auxfields_s
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nlocauxfields,nblocks_d) :: locauxfields_s
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nforces,nblocks_d) :: forces_s
      
      integer :: i,j,k
      !integer :: gi,gj,gk
      integer :: myblock,ii,jj,kk
      real(kind=db) :: gradfix,gradfiy,gradfiz,phi_loc,phi_out,mytemp
      real(kind=db) :: loc_u,loc_v,loc_w,lap_phi_loc
#ifdef MONOD
	  real(kind=db) :: S_mono
#endif
      
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
      
       
#ifdef TWOCOMPONENT
				  
	  mytemp= -sharp_c*real(locauxfields_s(ii,jj,kk,2,myblock),kind=db)
	  
      modgrad=real(auxfields_s(ii,jj,kk,4,myblock),kind=db) !modgrad
	  gradfix=real(auxfields_s(ii,jj,kk,1,myblock),kind=db)*modgrad !normx*modgrad
	  gradfiy=real(auxfields_s(ii,jj,kk,2,myblock),kind=db)*modgrad !normy*modgrad
	  gradfiz=real(auxfields_s(ii,jj,kk,3,myblock),kind=db)*modgrad !normz*modgrad
                  
      phi_loc=real(phifields_in(ii,jj,kk,1,myblock),kind=db)
      lap_phi_loc=real(locauxfields_s(ii,jj,kk,1,myblock),kind=db)
                  
      loc_u=real(hfields_s(ii,jj,kk,2,myblock),kind=db) !velocity
      loc_v=real(hfields_s(ii,jj,kk,3,myblock),kind=db)
      loc_w=real(hfields_s(ii,jj,kk,4,myblock),kind=db)

                  
      phi_out = phi_loc &
        - loc_u*gradfix - loc_v*gradfiy &
        - loc_w*gradfiz + tau_diff*lap_phi_loc + mytemp
#endif	

#ifdef MONOD
      S_mono = mu_max * phi_loc/(Ks + phi_loc) * phi_loc * (1.0_db - phi_loc)
      phi_out=phi_out + S_mono
	  !phi_out = min(1.0_db, max(0.0_db, phi_out))		 
#endif
      phifields_out(ii,jj,kk,1,myblock)=real(phi_out,kind=strdb)
      
      return
      
   endsubroutine update_phifields_kernel_yplus
   
   ! update_phifields_kernel_yminus: as update_phifields_kernel, for the first tile layer along y
   ! (by = 1), all x tiles, z tiles 2..nz/TILE_DIMz-1 (dimGridy).
   attributes(global) subroutine update_phifields_kernel_yminus(step,iprobe,jprobe,kprobe,flop,nx,ny,nz,coords,isfluid &
#ifdef TWOCOMPONENT    
    ,visc2,rho_r,rho_b,invrho_r,invrho_b,sharp_c,beta,kapp,tau_diff,sigma &
#ifdef MONOD	
    ,mu_max,Ks &
#endif
#endif 
    ,ntothfields,ntotphifields,ntotauxfields,ntotlocauxfields,ntotforces &
       ,hfields_s,phifields_in,phifields_out,auxfields_s,locauxfields_s,forces_s)

      implicit none
      integer :: step,iprobe,jprobe,kprobe,flop,nx,ny,nz
      integer, dimension(3) :: coords
      integer(kind=isf), dimension(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff) :: isfluid
     
 
#ifdef TWOCOMPONENT 
      real(kind=db) :: visc2,rho_r,rho_b,invrho_r,invrho_b,sharp_c,beta, &
       kapp,tau_diff,sigma,modgrad 
#ifdef MONOD
      real(kind=db) :: mu_max,Ks
#endif
#endif 
      integer :: ntothfields,ntotphifields,ntotauxfields,ntotlocauxfields,ntotforces
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields,nblocks_d) :: hfields_s
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields,nblocks_d) :: phifields_in,phifields_out
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nauxfields,nblocks_d) :: auxfields_s
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nlocauxfields,nblocks_d) :: locauxfields_s
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nforces,nblocks_d) :: forces_s
      
      integer :: i,j,k
      !integer :: gi,gj,gk
      integer :: myblock,ii,jj,kk
      real(kind=db) :: gradfix,gradfiy,gradfiz,phi_loc,phi_out,mytemp
      real(kind=db) :: loc_u,loc_v,loc_w,lap_phi_loc
#ifdef MONOD
	  real(kind=db) :: S_mono
#endif
      
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
      
	   
#ifdef TWOCOMPONENT
				  
	  mytemp= -sharp_c*real(locauxfields_s(ii,jj,kk,2,myblock),kind=db)
	  
      modgrad=real(auxfields_s(ii,jj,kk,4,myblock),kind=db) !modgrad
	  gradfix=real(auxfields_s(ii,jj,kk,1,myblock),kind=db)*modgrad !normx*modgrad
	  gradfiy=real(auxfields_s(ii,jj,kk,2,myblock),kind=db)*modgrad !normy*modgrad
	  gradfiz=real(auxfields_s(ii,jj,kk,3,myblock),kind=db)*modgrad !normz*modgrad
                  
      phi_loc=real(phifields_in(ii,jj,kk,1,myblock),kind=db)
      lap_phi_loc=real(locauxfields_s(ii,jj,kk,1,myblock),kind=db)
                  
      loc_u=real(hfields_s(ii,jj,kk,2,myblock),kind=db) !velocity
      loc_v=real(hfields_s(ii,jj,kk,3,myblock),kind=db)
      loc_w=real(hfields_s(ii,jj,kk,4,myblock),kind=db)

                  
      phi_out = phi_loc &
        - loc_u*gradfix - loc_v*gradfiy &
        - loc_w*gradfiz + tau_diff*lap_phi_loc + mytemp
#endif	

#ifdef MONOD
      S_mono = mu_max * phi_loc/(Ks + phi_loc) * phi_loc * (1.0_db - phi_loc)
      phi_out=phi_out + S_mono
	  !phi_out = min(1.0_db, max(0.0_db, phi_out))		 
#endif
      phifields_out(ii,jj,kk,1,myblock)=real(phi_out,kind=strdb)
      
      return
      
   endsubroutine update_phifields_kernel_yminus
   
   ! update_phifields_kernel_zplus: as update_phifields_kernel, for the last tile layer along z
   ! (bz = nzblock_d-2), all x and y tiles (dimGridz); not launched when there is a single tile
   ! along z.
   attributes(global) subroutine update_phifields_kernel_zplus(step,iprobe,jprobe,kprobe,flop,nx,ny,nz,coords,isfluid &
#ifdef TWOCOMPONENT    
    ,visc2,rho_r,rho_b,invrho_r,invrho_b,sharp_c,beta,kapp,tau_diff,sigma &
#ifdef MONOD	
    ,mu_max,Ks &
#endif
#endif 
    ,ntothfields,ntotphifields,ntotauxfields,ntotlocauxfields,ntotforces &
       ,hfields_s,phifields_in,phifields_out,auxfields_s,locauxfields_s,forces_s)

      implicit none
      integer :: step,iprobe,jprobe,kprobe,flop,nx,ny,nz
      integer, dimension(3) :: coords
      integer(kind=isf), dimension(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff) :: isfluid
     
 
#ifdef TWOCOMPONENT 
      real(kind=db) :: visc2,rho_r,rho_b,invrho_r,invrho_b,sharp_c,beta, &
       kapp,tau_diff,sigma,modgrad 
#ifdef MONOD
      real(kind=db) :: mu_max,Ks
#endif
#endif 
      integer :: ntothfields,ntotphifields,ntotauxfields,ntotlocauxfields,ntotforces
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields,nblocks_d) :: hfields_s
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields,nblocks_d) :: phifields_in,phifields_out
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nauxfields,nblocks_d) :: auxfields_s
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nlocauxfields,nblocks_d) :: locauxfields_s
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nforces,nblocks_d) :: forces_s
      
      integer :: i,j,k
      !integer :: gi,gj,gk
      integer :: myblock,ii,jj,kk
      real(kind=db) :: gradfix,gradfiy,gradfiz,phi_loc,phi_out,mytemp
      real(kind=db) :: loc_u,loc_v,loc_w,lap_phi_loc
#ifdef MONOD
	  real(kind=db) :: S_mono
#endif
      
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
      
       
#ifdef TWOCOMPONENT
				  
	  mytemp= -sharp_c*real(locauxfields_s(ii,jj,kk,2,myblock),kind=db)
	  
      modgrad=real(auxfields_s(ii,jj,kk,4,myblock),kind=db) !modgrad
	  gradfix=real(auxfields_s(ii,jj,kk,1,myblock),kind=db)*modgrad !normx*modgrad
	  gradfiy=real(auxfields_s(ii,jj,kk,2,myblock),kind=db)*modgrad !normy*modgrad
	  gradfiz=real(auxfields_s(ii,jj,kk,3,myblock),kind=db)*modgrad !normz*modgrad
                  
      phi_loc=real(phifields_in(ii,jj,kk,1,myblock),kind=db)
      lap_phi_loc=real(locauxfields_s(ii,jj,kk,1,myblock),kind=db)
                  
      loc_u=real(hfields_s(ii,jj,kk,2,myblock),kind=db) !velocity
      loc_v=real(hfields_s(ii,jj,kk,3,myblock),kind=db)
      loc_w=real(hfields_s(ii,jj,kk,4,myblock),kind=db)

                  
      phi_out = phi_loc &
        - loc_u*gradfix - loc_v*gradfiy &
        - loc_w*gradfiz + tau_diff*lap_phi_loc + mytemp
#endif	

#ifdef MONOD
      S_mono = mu_max * phi_loc/(Ks + phi_loc) * phi_loc * (1.0_db - phi_loc)
      phi_out=phi_out + S_mono
	  !phi_out = min(1.0_db, max(0.0_db, phi_out))		 
#endif
      phifields_out(ii,jj,kk,1,myblock)=real(phi_out,kind=strdb)
      
      return
      
   endsubroutine update_phifields_kernel_zplus
   
   ! update_phifields_kernel_zminus: as update_phifields_kernel, for the first tile layer along z
   ! (bz = 1), all x and y tiles (dimGridz).
   attributes(global) subroutine update_phifields_kernel_zminus(step,iprobe,jprobe,kprobe,flop,nx,ny,nz,coords,isfluid &
#ifdef TWOCOMPONENT    
    ,visc2,rho_r,rho_b,invrho_r,invrho_b,sharp_c,beta,kapp,tau_diff,sigma &
#ifdef MONOD	
    ,mu_max,Ks &
#endif
#endif 
    ,ntothfields,ntotphifields,ntotauxfields,ntotlocauxfields,ntotforces &
       ,hfields_s,phifields_in,phifields_out,auxfields_s,locauxfields_s,forces_s)

      implicit none
      integer :: step,iprobe,jprobe,kprobe,flop,nx,ny,nz
      integer, dimension(3) :: coords
      integer(kind=isf), dimension(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff) :: isfluid
     
 
#ifdef TWOCOMPONENT 
      real(kind=db) :: visc2,rho_r,rho_b,invrho_r,invrho_b,sharp_c,beta, &
       kapp,tau_diff,sigma,modgrad 
#ifdef MONOD
      real(kind=db) :: mu_max,Ks
#endif
#endif 
      integer :: ntothfields,ntotphifields,ntotauxfields,ntotlocauxfields,ntotforces
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields,nblocks_d) :: hfields_s
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields,nblocks_d) :: phifields_in,phifields_out
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nauxfields,nblocks_d) :: auxfields_s
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nlocauxfields,nblocks_d) :: locauxfields_s
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nforces,nblocks_d) :: forces_s
      
      integer :: i,j,k
      !integer :: gi,gj,gk
      integer :: myblock,ii,jj,kk
      real(kind=db) :: gradfix,gradfiy,gradfiz,phi_loc,phi_out,mytemp
      real(kind=db) :: loc_u,loc_v,loc_w,lap_phi_loc
#ifdef MONOD
	  real(kind=db) :: S_mono
#endif
      
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
      
	   
#ifdef TWOCOMPONENT
				  
	  mytemp= -sharp_c*real(locauxfields_s(ii,jj,kk,2,myblock),kind=db)
	  
      modgrad=real(auxfields_s(ii,jj,kk,4,myblock),kind=db) !modgrad
	  gradfix=real(auxfields_s(ii,jj,kk,1,myblock),kind=db)*modgrad !normx*modgrad
	  gradfiy=real(auxfields_s(ii,jj,kk,2,myblock),kind=db)*modgrad !normy*modgrad
	  gradfiz=real(auxfields_s(ii,jj,kk,3,myblock),kind=db)*modgrad !normz*modgrad
                  
      phi_loc=real(phifields_in(ii,jj,kk,1,myblock),kind=db)
      lap_phi_loc=real(locauxfields_s(ii,jj,kk,1,myblock),kind=db)
                  
      loc_u=real(hfields_s(ii,jj,kk,2,myblock),kind=db) !velocity
      loc_v=real(hfields_s(ii,jj,kk,3,myblock),kind=db)
      loc_w=real(hfields_s(ii,jj,kk,4,myblock),kind=db)

                  
      phi_out = phi_loc &
        - loc_u*gradfix - loc_v*gradfiy &
        - loc_w*gradfiz + tau_diff*lap_phi_loc + mytemp
#endif	

#ifdef MONOD
      S_mono = mu_max * phi_loc/(Ks + phi_loc) * phi_loc * (1.0_db - phi_loc)
      phi_out=phi_out + S_mono
	  !phi_out = min(1.0_db, max(0.0_db, phi_out))		 
#endif
      phifields_out(ii,jj,kk,1,myblock)=real(phi_out,kind=strdb)
      
      return
      
   endsubroutine update_phifields_kernel_zminus
#endif

endmodule lb_cuda_update_phi
