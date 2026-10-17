#include "defines.h"
#if !defined(_OPENACC)  
#error "To use this module the macros _OPENACC should be defined."
#endif

! Moment update before the collision. The kernels of this module compute the total force on every
! fluid node, add the first half of the forcing to the stored velocity and second-order moments, and
! store the force that the fused collide-and-stream kernels (lb_cuda_fused_module.f90) use for the
! second half. They are called when FUSED_MOMENTS is not defined, i.e. in two-component flows and
! with INTERNAL_OBSTACLES, MULTIHIT or ELASTIC_FORCE; otherwise the fused kernels apply the update
! for the uniform body force themselves (see docs/FUSED_KERNELS.md).
!
! Storage by tiles: field(ii,jj,kk,ifield,myblock), with (ii,jj,kk) the position of the node in a
! TILE_DIMx x TILE_DIMy x TILE_DIMz tile. Along x, tile t = 1..nx/TILE_DIMx contains the nodes
! i = (t-1)*TILE_DIMx+1..t*TILE_DIMx and has the zero-based block coordinate bx = t, while bx = 0
! and bx = nxblock_d-1 = nx/TILE_DIMx+1 are halo tiles; likewise along y and z, and
! myblock = bx + by*nxblock_d + bz*nxyblock_d + 1.
! Fields used in this module:
!   hfields      1 p* = p/(rho cs^2); 2-4 velocity u, v, w; 5-10 second-order Hermite moments
!                a2_ab = sum_l (c_la c_lb - cs^2 delta_ab) f_l, stored as pxx, pyy, pzz, pxy, pxz, pyz,
!                whose equilibrium value is u_a u_b
!   phifields    1 phase field phi: 1 in component r (density rho_r), 0 in component b (rho_b)
!   auxfields    1-3 unit normal n = grad(phi)/|grad(phi)|, 4 |grad(phi)|
!   locauxfields 1 lap(phi), 2 div(theta) with theta = phi(1-phi) n; 6-8 repulsive flux (REPULSIVE_FLUX)
!   forces       1-3 force per unit volume
!
! Variants: moments_LB_kernel covers all the tiles. With ASYNCMPI, moments_LB_cuda_ext first launches
! the six face kernels, which cover the tiles that touch the subdomain boundary, then the halo exchange
! of hfields and forces is started and moments_LB_kernel_int updates the interior tiles before it is
! completed. All variants have the same body and differ only in the tiles they cover.
module lb_cuda_moments

   use vars
   use iso_c_binding
   use cudafor
   use mpi_template, only: coords,dostop,doerror,mydev,myrank,nprocs,nbuff,nbuffbvec
   use lb_cuda_vars

   implicit none
   

contains
   
    ! moments_LB_kernel: total force and moment update for all the fluid nodes of the subdomain
    ! (dimGrid: one block per tile, one thread per node). The force is built in the order of
    ! Lauricella et al., Phys. Fluids 37, 072111 (2025), Secs. II C-D:
    !  1. velocity-independent part F0: body force, capillary force (CSF or chemical potential),
    !     REPULSIVE_FLUX and MULTIHIT forcing, and with DENSRATIO the pressure force -p* cs^2 grad(rho);
    !  2. predictor F = F0 + velocity-dependent terms (ELASTIC_FORCE; with DENSRATIO the viscous force)
    !     evaluated with the velocity of the previous time step, read from hfields_old;
    !  3. half-force update with F: u = u + F/(2 rho) (with INTERFACE_INCOMP also divided by
    !     1 - q/(2 rho), see below) and a2 = a2 + (u F + F u)/(2 rho), written to hfields_s;
    !  4. corrected force tforce = F0 + velocity-dependent terms evaluated with the updated u and a2,
    !     written to forces_s; the fused kernel uses it for the second half of the forcing.
    ! rho = rho_r phi + rho_b (1-phi) with DENSRATIO and rho = 1 otherwise.
    ! Reads hfields_s(1:10) (moments after streaming), hfields_old(2:4) (hfields_old is the other
    ! buffer of the flip/flop pair, holding the moments updated by this kernel at the previous step),
    ! phifields_s(1), auxfields_s(1:4), locauxfields_s(1:2) (and 6:8 with REPULSIVE_FLUX), forces_s
    ! (only with CSF). Writes hfields_s(2:10), forces_s(1:3), and with ELASTIC_FORCE u_ref, v_ref,
    ! w_ref. step, iprobe, jprobe, kprobe, flop, the ntot* sizes, wettab_r, wettab_b, invrho_r,
    ! invrho_b and sigma are not used.
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
       ,hfields_old,hfields_s &
#ifdef TWOCOMPONENT        
       ,auxfields_s,locauxfields_s &
#endif
       ,forces_s)
 

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
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields,nblocks_d) :: phifields_s 
#endif           
#if defined(ELASTIC_FORCE)
      real(kind=db) :: lambda_rel,k_elastic
      real(kind=db), dimension(1:nx,1:ny,1:nz) :: u_ref,v_ref,w_ref
#endif
      real(kind=db) :: fx,fy,fz

      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields,nblocks_d) :: hfields_old,hfields_s
#ifdef TWOCOMPONENT 
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nauxfields,nblocks_d) :: auxfields_s
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nlocauxfields,nblocks_d) :: locauxfields_s
#endif
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nforces,nblocks_d) :: forces_s
      

      real(kind=db) :: pxx,pyy,pzz,pxy,pxz,pyz
      real(kind=db) :: mytemp,forcex,forcey,forcez,rhophi_loc,press_loc,u_loc,v_loc,w_loc
      real(kind=db) :: tforcex,tforcey,tforcez

#ifdef TWOCOMPONENT
	  real(kind=db) ::gradfix,gradfiy,gradfiz,wet_loc,phi_loc,lap_phi_loc,mymask
#endif
      
#ifdef DENSRATIO
      real(kind=db) :: gradrhox,gradrhoy,gradrhoz,omega_loc, visc_loc,tau_loc
#endif
      
      integer :: i,j,k
      integer :: gi,gj,gk
      integer :: myblock,ii,jj,kk
      
      ! tile number = blockIdx in each direction; i, j, k are the node indices in the subdomain
      i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z
      
      ! only fluid nodes are updated (isfluid = 1, or -1 for fluid nodes next to a non-fluid node)
      if(abs(isfluid(i,j,k)) .ne. 1)return
      
      ! position in the global lattice (not used below)
      gi=nx*coords(1)+i
      gj=ny*coords(2)+j
      gk=nz*coords(3)+k
      
      ! storage tile: the block coordinates equal the tile numbers, because coordinate 0 is the halo
      ! layer of tiles; ii, jj, kk are the node position inside the tile
      myblock=blockIdx%x+blockIdx%y*nxblock_d+blockIdx%z*nxyblock_d+1
      
      ii=threadIdx%x
      jj=threadIdx%y
      kk=threadIdx%z
      

                  ! p* = p/(rho cs^2)
		  press_loc=real(hfields_s(ii,jj,kk,1,myblock),kind=db)
			 
#ifdef TWOCOMPONENT
		  phi_loc=real(phifields_s(ii,jj,kk,1,myblock),kind=db)
		  lap_phi_loc=real(locauxfields_s(ii,jj,kk,1,myblock),kind=db)
          ! interface mask 2 sqrt(phi(1-phi)), zero outside 0 < phi < 1: it is 1 at phi = 1/2 and vanishes
          ! in the bulk phases; it weights the pressure and viscous forces below
!          mymask=ONE-HALF*(ONE-tanh((phi_loc-0.015_db)/(HALF*0.005_db)))
!          mymask=    HALF*(ONE-tanh((phi_loc-0.995_db)/(HALF*0.005_db)))*mymask
          mymask=sqrt(max(ZERO,FOUR*phi_loc*(ONE-phi_loc)))
#endif
                  ! body force rho*(fx,fy,fz); with BUOYANCY_FORCING (rho - (rho_r+rho_b)/2)*(fx,fy,fz), i.e.
                  ! relative to the mean density of the two components; rho = 1 without DENSRATIO
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
		   
		  ! Selectable capillary force: geometric CSF or Jacqmin chemical potential.
                  ! grad(phi) = n |grad(phi)| (auxfields 1-3 and 4), also used for grad(rho) below
		  mytemp=real(auxfields_s(ii,jj,kk,4,myblock),kind=db) !modgrad
		  gradfix=real(auxfields_s(ii,jj,kk,1,myblock),kind=db)*mytemp !normx*modgrad
		  gradfiy=real(auxfields_s(ii,jj,kk,2,myblock),kind=db)*mytemp !normy*modgrad
		  gradfiz=real(auxfields_s(ii,jj,kk,3,myblock),kind=db)*mytemp !normz*modgrad

                  ! CSF: forces_s holds sigma*curvature*grad(phi), curvature = -div(n), written by
                  ! compute_csf_force_kernel. Otherwise mu*grad(phi), with the chemical potential
                  ! mu = 4 beta phi (phi-1)(phi-1/2) - kapp lap(phi), beta = 12 sigma/width,
                  ! kapp = 3 sigma width/2 (Lauricella et al. 2025, Eq. (3))
#ifdef CSF
		  forcex=forcex+real(forces_s(ii,jj,kk,1,myblock),kind=db)
		  forcey=forcey+real(forces_s(ii,jj,kk,2,myblock),kind=db)
		  forcez=forcez+real(forces_s(ii,jj,kk,3,myblock),kind=db)
#else
		  forcex = forcex + &
                   (4.0_db*beta*phi_loc*(phi_loc-1.0_db)*(phi_loc-0.5_db) - kapp*lap_phi_loc)*gradfix
		  forcey = forcey + &
                   (4.0_db*beta*phi_loc*(phi_loc-1.0_db)*(phi_loc-0.5_db) - kapp*lap_phi_loc)*gradfiy
		  forcez = forcez + &
                   (4.0_db*beta*phi_loc*(phi_loc-1.0_db)*(phi_loc-0.5_db) - kapp*lap_phi_loc)*gradfiz
#endif
				
				   				   

                  ! REPULSIVE_FLUX: force rho*(rho J), each component of rho J limited to magnitude 1e-3,
                  ! where J (locauxfields 6-8) is the repulsive flux between facing interfaces of a thin
                  ! film, computed by repulsive_flux_normal_kernel
#ifdef REPULSIVE_FLUX
                  mytemp=real(locauxfields_s(ii,jj,kk,6,myblock),kind=db)*rhophi_loc 
                  if(abs(mytemp)>1.0d-3) mytemp=1.0d-3*sign(1.0,mytemp)!mytemp*0.1_db
                  forcex=forcex + mytemp*rhophi_loc
		  
                  mytemp=real(locauxfields_s(ii,jj,kk,7,myblock),kind=db)*rhophi_loc 
                  if(abs(mytemp)>1.0d-3) mytemp=1.0d-3*sign(1.0,mytemp)!mytemp*0.1_db
                  forcey=forcey + mytemp*rhophi_loc
		  
                  mytemp=real(locauxfields_s(ii,jj,kk,8,myblock),kind=db)*rhophi_loc 
                  if(abs(mytemp)>1.0d-3) mytemp=1.0d-3*sign(1.0,mytemp)!mytemp*0.1_db
                  forcez=forcez + mytemp*rhophi_loc
#endif
#endif


                  ! MULTIHIT: forcing field (ABCx, ABCy, ABCz) per unit mass, multiplied by rho
#if defined(MULTIHIT)   
                  forcex=forcex + rhophi_loc*ABCx(i,j,k) !+ AAA*sin(k_zero*gk) + AAA*sin(k_zero*gj)  
                  forcey=forcey + rhophi_loc*ABCy(i,j,k) !+ AAA*sin(k_zero*gi) + AAA*sin(k_zero*gk)
                  forcez=forcez + rhophi_loc*ABCz(i,j,k) !+ AAA*sin(k_zero*gj) + AAA*sin(k_zero*gi) 	
#endif
                  
#ifdef DENSRATIO				  
		  ! pressure and viscous forces
				  
                  ! grad(rho) = (rho_r-rho_b) grad(phi). Pressure force -p* cs^2 grad(rho), Eq. (18) of the
                  ! article: grad(p) = rho cs^2 grad(p*) + p* cs^2 grad(rho), and the equilibrium written in
                  ! terms of p* accounts only for the first term
                  gradrhox=(rho_r-rho_b)*gradfix
                  gradrhoy=(rho_r-rho_b)*gradfiy
                  gradrhoz=(rho_r-rho_b)*gradfiz
				  
                  forcex=forcex - mymask* &
                   press_loc*cssq*gradrhox   !+ (rhophi_loc-(rho_r+rho_b)*HALF)*fx
                  forcey=forcey - mymask* &
                   press_loc*cssq*gradrhoy   !+ (rhophi_loc-(rho_r+rho_b)*HALF)*fy
                  forcez=forcez - mymask* &
                   press_loc*cssq*gradrhoz   !+ (rhophi_loc-(rho_r+rho_b)*HALF)*fz
		  !! from this point I compute the force terms that depend on the velocity
		  !! these terms should be not included in force arrays since they must be computed with the updated velocity
		  !! at the end of this subroutine
#endif
		  
                  ! tforce: velocity-independent part of the force; the velocity-dependent terms are added
                  ! to it below, evaluated with the updated velocity
                  tforcex=forcex
                  tforcey=forcey
                  tforcez=forcez

                  ! predictor: velocity-dependent terms evaluated with the velocity of the previous time step
                  u_loc=real(hfields_old(ii,jj,kk,2,myblock),kind=db) !velocity
                  v_loc=real(hfields_old(ii,jj,kk,3,myblock),kind=db)
                  w_loc=real(hfields_old(ii,jj,kk,4,myblock),kind=db)       
                 
				  
                  ! ELASTIC_FORCE: force rho*phi*k_elastic*lambda_rel*(u - u_ref), where u_ref is a running
                  ! average of the velocity (updated below)
#if defined(ELASTIC_FORCE)
                  forcex=forcex + &
                   rhophi_loc*(u_loc - u_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+   
                  forcey=forcey + &
                   rhophi_loc*(v_loc - v_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+ 
                  forcez=forcez + &
                   rhophi_loc*(w_loc - w_ref(i,j,k))*k_elastic*lambda_rel*phi_loc + rhophi_loc*fz !+  	
#endif

                  ! viscous force -mymask*(nu/(tau cs^2)) Pi_neq.grad(rho), with Pi_neq = a2 - u u the
                  ! non-equilibrium second-order moment (Eq. (19) of the article, Eq. (16) of Kumar et al.,
                  ! Phys. Fluids 31, 072103 (2019)). To leading order Pi_neq = -tau cs^2 (grad u + grad u^T),
                  ! so the force is nu (grad u + grad u^T).grad(rho). Local kinematic viscosity
                  ! nu = (rho_r visc1 phi + rho_b visc2 (1-phi))/rho and tau = nu/cs^2 + 1/2. Here a2 is the
                  ! moment after streaming and u the velocity of the previous step.
#ifdef DENSRATIO 
                  pxx=real(hfields_s(ii,jj,kk,5,myblock),kind=db)
                  pyy=real(hfields_s(ii,jj,kk,6,myblock),kind=db)
                  pzz=real(hfields_s(ii,jj,kk,7,myblock),kind=db)
                  pxy=real(hfields_s(ii,jj,kk,8,myblock),kind=db)
                  pxz=real(hfields_s(ii,jj,kk,9,myblock),kind=db)
                  pyz=real(hfields_s(ii,jj,kk,10,myblock),kind=db)   
                               
                  !1-2
                  !*1
                  ! 2nd order
                  pxx=pxx - u_loc*u_loc
                  pyy=pyy - v_loc*v_loc
                  pzz=pzz - w_loc*w_loc
                  pxy=pxy - u_loc*v_loc
                  pxz=pxz - u_loc*w_loc
                  pyz=pyz - v_loc*w_loc

		  visc_loc=(rho_r*visc1*phi_loc+(1.0_db-phi_loc)*visc2*rho_b)/rhophi_loc
				  
                  tau_loc=(visc_loc*invcssq + HALF) !this is a tau
				  
		  forcex=forcex - mymask* &
		   (visc_loc/(tau_loc*cssq))*(pxx*gradrhox + pxy*gradrhoy + pxz*gradrhoz)
		  forcey=forcey - mymask* &
		   (visc_loc/(tau_loc*cssq))*(pyy*gradrhoy + pxy*gradrhox + pyz*gradrhoz)
		  forcez=forcez - mymask* &
		   (visc_loc/(tau_loc*cssq))*(pzz*gradrhoz + pxz*gradrhox + pyz*gradrhoy)
#endif	
                  !I compute the new velocities
                  ! first moment after streaming plus half the predictor force per unit mass, u + F/(2 rho)
                  ! (Eq. (14) of the article)
                  u_loc=real(hfields_s(ii,jj,kk,2,myblock),kind=db) !velocity
                  v_loc=real(hfields_s(ii,jj,kk,3,myblock),kind=db)
                  w_loc=real(hfields_s(ii,jj,kk,4,myblock),kind=db)  
					 
                  ! INTERFACE_INCOMP (defined with DENSRATIO): interfacial force F_a = -q u of Kumar et al.,
                  ! Phys. Fluids 31, 072103 (2019), Eqs. (17)-(18), with q = (rho_r-rho_b)*(tau_diff lap(phi)
                  ! - sharp_c div(theta)) = drho/dt + div(rho u), the density change due to the diffusive and
                  ! sharpening terms of the phase-field equation. As in Eq. (18) there, the velocity is
                  ! (u + F/(2 rho))/(1 - q/(2 rho)); F_a is then evaluated with this velocity and added to the
                  ! predictor force and to tforce.
#if defined(INTERFACE_INCOMP) && defined(DENSRATIO)
                  mytemp= -sharp_c*real(locauxfields_s(ii,jj,kk,2,myblock),kind=db)
                  u_loc = u_loc + 0.5_db*forcex/(rhophi_loc)
                  v_loc = v_loc + 0.5_db*forcey/(rhophi_loc)
                  w_loc = w_loc + 0.5_db*forcez/(rhophi_loc)
				  
                  u_loc = u_loc/ &
                   (1.0_db - 0.5_db*(rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)/rhophi_loc )
				  
                  v_loc = v_loc/ &
                   (1.0_db - 0.5_db*(rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)/rhophi_loc ) 
				  
                  w_loc = w_loc/ &
                   (1.0_db - 0.5_db*(rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)/rhophi_loc ) 

                  tforcex= tforcex - &
                   (rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)*u_loc
                  tforcey= tforcey - &
                   (rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)*v_loc
                  tforcez= tforcez - &
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
                  
                  hfields_s(ii,jj,kk,2,myblock)=real(u_loc,kind=strdb) !put the new velocity in hfields_s
                  hfields_s(ii,jj,kk,3,myblock)=real(v_loc,kind=strdb)
                  hfields_s(ii,jj,kk,4,myblock)=real(w_loc,kind=strdb)           
                                    

!regularized 
                  pxx=real(hfields_s(ii,jj,kk,5,myblock),kind=db)
                  pyy=real(hfields_s(ii,jj,kk,6,myblock),kind=db)
                  pzz=real(hfields_s(ii,jj,kk,7,myblock),kind=db)
                  pxy=real(hfields_s(ii,jj,kk,8,myblock),kind=db)
                  pxz=real(hfields_s(ii,jj,kk,9,myblock),kind=db)
                  pyz=real(hfields_s(ii,jj,kk,10,myblock),kind=db)
                  
                  ! forcing term of the second-order moments, a2 = a2 + (u F + F u)/(2 rho), with the updated
                  ! velocity and the predictor force
                  pxx=pxx + forcex*u_loc/rhophi_loc
                  pyy=pyy + forcey*v_loc/rhophi_loc
                  pzz=pzz + forcez*w_loc/rhophi_loc
                  pxy=pxy + HALF*(forcey*u_loc+forcex*v_loc)/rhophi_loc
                  pxz=pxz + HALF*(forcez*u_loc+forcex*w_loc)/rhophi_loc
                  pyz=pyz + HALF*(forcez*v_loc+forcey*w_loc)/rhophi_loc

                  hfields_s(ii,jj,kk,5,myblock)=real(pxx,kind=strdb)
                  hfields_s(ii,jj,kk,6,myblock)=real(pyy,kind=strdb)
                  hfields_s(ii,jj,kk,7,myblock)=real(pzz,kind=strdb)
                  hfields_s(ii,jj,kk,8,myblock)=real(pxy,kind=strdb)
                  hfields_s(ii,jj,kk,9,myblock)=real(pxz,kind=strdb)
                  hfields_s(ii,jj,kk,10,myblock)=real(pyz,kind=strdb)              
                  
                  ! ELASTIC_FORCE: u_ref relaxes towards the updated velocity with rate lambda_rel, and the
                  ! elastic force evaluated with the updated u and u_ref is added to tforce
#if defined(ELASTIC_FORCE)
                  u_ref(i,j,k) = u_ref(i,j,k) + &
                   lambda_rel*(u_loc - u_ref(i,j,k))
                  v_ref(i,j,k) = v_ref(i,j,k) + &
                   lambda_rel*(v_loc - v_ref(i,j,k))
                  w_ref(i,j,k) = w_ref(i,j,k) + &
                   lambda_rel*(w_loc - w_ref(i,j,k))
				  
                  tforcex = tforcex + &
                   rhophi_loc*(u_loc - u_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+ 
                  tforcey = tforcey +&
                   rhophi_loc*(v_loc - v_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+ 
                  tforcez = tforcez + &
                   rhophi_loc*(w_loc - w_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+rhophi_loc*fz  
#endif               
                  
                  ! corrector: viscous force evaluated with the updated velocity and second-order moments,
                  ! added to tforce
#if defined(DENSRATIO)	
		  
                  pxx=pxx - u_loc*u_loc 
                  pyy=pyy - v_loc*v_loc
                  pzz=pzz - w_loc*w_loc 
                  pxy=pxy - u_loc*v_loc
                  pxz=pxz - u_loc*w_loc
                  pyz=pyz - v_loc*w_loc
                  
                  tforcex = tforcex - mymask* &
                   (visc_loc/(tau_loc*cssq))*(pxx*gradrhox + pxy*gradrhoy + pxz*gradrhoz)
                  tforcey = tforcey - mymask* &
                   (visc_loc/(tau_loc*cssq))*(pyy*gradrhoy + pxy*gradrhox + pyz*gradrhoz)
                  tforcez = tforcez - mymask* &
                   (visc_loc/(tau_loc*cssq))*(pzz*gradrhoz + pxz*gradrhox + pyz*gradrhoy)
#endif  
        
                  ! corrected total force for the fused collide-and-stream kernel, which divides it by rho;
                  ! with CSF it replaces the capillary force read above
                  forces_s(ii,jj,kk,1,myblock)=real(tforcex,kind=strdb)
		          forces_s(ii,jj,kk,2,myblock)=real(tforcey,kind=strdb)
		          forces_s(ii,jj,kk,3,myblock)=real(tforcez,kind=strdb)

   endsubroutine moments_LB_kernel  
   
   ! moments_LB_kernel_int: as moments_LB_kernel, for the interior tiles, i.e. the tiles that do not
   ! touch the subdomain boundary (tiles 2..n/TILE_DIM-1 in each direction). Launched on dimGridInt =
   ! (nx/TILE_DIMx-2, ny/TILE_DIMy-2, nz/TILE_DIMz-2) blocks when all three are positive.
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
       ,hfields_old,hfields_s &
#ifdef TWOCOMPONENT        
       ,auxfields_s,locauxfields_s &
#endif
       ,forces_s)
 

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
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields,nblocks_d) :: phifields_s 
#endif           
#if defined(ELASTIC_FORCE)
      real(kind=db) :: lambda_rel,k_elastic
      real(kind=db), dimension(1:nx,1:ny,1:nz) :: u_ref,v_ref,w_ref
#endif
      real(kind=db) :: fx,fy,fz

      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields,nblocks_d) :: hfields_old,hfields_s
#ifdef TWOCOMPONENT 
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nauxfields,nblocks_d) :: auxfields_s
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nlocauxfields,nblocks_d) :: locauxfields_s
#endif
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nforces,nblocks_d) :: forces_s
      

      real(kind=db) :: pxx,pyy,pzz,pxy,pxz,pyz
      real(kind=db) :: mytemp,forcex,forcey,forcez,rhophi_loc,press_loc,u_loc,v_loc,w_loc
      real(kind=db) :: tforcex,tforcey,tforcez
      
#ifdef TWOCOMPONENT
	  real(kind=db) ::gradfix,gradfiy,gradfiz,wet_loc,phi_loc,lap_phi_loc,mymask
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
      
      ! block b = 1..nx/TILE_DIMx-2 handles tile t = b+1, i.e. the nodes i = (t-1)*TILE_DIMx + ii
      i = (blockIdx%x) * TILE_DIMx + ii
      j = (blockIdx%y) * TILE_DIMy + jj
      k = (blockIdx%z) * TILE_DIMz + kk
      
      if(abs(isfluid(i,j,k)) .ne. 1)return
      
      !gi=nx*coords(1)+i
      !gj=ny*coords(2)+j
      !gk=nz*coords(3)+k
      
      ! tile t has block coordinate t (coordinate 0 is the halo tile), hence blockIdx+1
      myblock=(blockIdx%x+1)+(blockIdx%y+1)*nxblock_d+(blockIdx%z+1)*nxyblock_d+1
				 
		  press_loc=real(hfields_s(ii,jj,kk,1,myblock),kind=db)
			 
#ifdef TWOCOMPONENT
		  phi_loc=real(phifields_s(ii,jj,kk,1,myblock),kind=db)
		  lap_phi_loc=real(locauxfields_s(ii,jj,kk,1,myblock),kind=db)
!          mymask=ONE-HALF*(ONE-tanh((phi_loc-0.015_db)/(HALF*0.005_db)))
!          mymask=    HALF*(ONE-tanh((phi_loc-0.995_db)/(HALF*0.005_db)))*mymask
          mymask=sqrt(max(ZERO,FOUR*phi_loc*(ONE-phi_loc)))
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
		   
		  ! Selectable capillary force: geometric CSF or Jacqmin chemical potential.
		  mytemp=real(auxfields_s(ii,jj,kk,4,myblock),kind=db) !modgrad
		  gradfix=real(auxfields_s(ii,jj,kk,1,myblock),kind=db)*mytemp !normx*modgrad
		  gradfiy=real(auxfields_s(ii,jj,kk,2,myblock),kind=db)*mytemp !normy*modgrad
		  gradfiz=real(auxfields_s(ii,jj,kk,3,myblock),kind=db)*mytemp !normz*modgrad

#ifdef CSF
		  forcex=forcex+real(forces_s(ii,jj,kk,1,myblock),kind=db)
		  forcey=forcey+real(forces_s(ii,jj,kk,2,myblock),kind=db)
		  forcez=forcez+real(forces_s(ii,jj,kk,3,myblock),kind=db)
#else
		  forcex = forcex + &
                   (4.0_db*beta*phi_loc*(phi_loc-1.0_db)*(phi_loc-0.5_db) - kapp*lap_phi_loc)*gradfix
		  forcey = forcey + &
                   (4.0_db*beta*phi_loc*(phi_loc-1.0_db)*(phi_loc-0.5_db) - kapp*lap_phi_loc)*gradfiy
		  forcez = forcez + &
                   (4.0_db*beta*phi_loc*(phi_loc-1.0_db)*(phi_loc-0.5_db) - kapp*lap_phi_loc)*gradfiz
#endif
				   				   

#ifdef REPULSIVE_FLUX
                  mytemp=real(locauxfields_s(ii,jj,kk,6,myblock),kind=db)*rhophi_loc 
                  if(abs(mytemp)>1.0d-3) mytemp=1.0d-3*sign(1.0,mytemp)!mytemp*0.1_db
                  forcex=forcex + mytemp*rhophi_loc
		  
                  mytemp=real(locauxfields_s(ii,jj,kk,7,myblock),kind=db)*rhophi_loc 
                  if(abs(mytemp)>1.0d-3) mytemp=1.0d-3*sign(1.0,mytemp)!mytemp*0.1_db
                  forcey=forcey + mytemp*rhophi_loc
		  
                  mytemp=real(locauxfields_s(ii,jj,kk,8,myblock),kind=db)*rhophi_loc 
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
				  
          forcex=forcex - mymask* &
           press_loc*cssq*gradrhox   !+ (rhophi_loc-(rho_r+rho_b)*HALF)*fx
          forcey=forcey - mymask* &
           press_loc*cssq*gradrhoy   !+ (rhophi_loc-(rho_r+rho_b)*HALF)*fy
		  forcez=forcez - mymask* &
		   press_loc*cssq*gradrhoz   !+ (rhophi_loc-(rho_r+rho_b)*HALF)*fz
		  !! from this point I compute the force terms that depend on the velocity
		  !! these terms should be not included in force arrays since they must be computed with the updated velocity
		  !! at the end of this subroutine
#endif
                  
              
                  
                  tforcex=forcex
                  tforcey=forcey
                  tforcez=forcez
				  
                  u_loc=real(hfields_old(ii,jj,kk,2,myblock),kind=db) !velocity
                  v_loc=real(hfields_old(ii,jj,kk,3,myblock),kind=db)
                  w_loc=real(hfields_old(ii,jj,kk,4,myblock),kind=db)
                  
                 
				  
#if defined(ELASTIC_FORCE)
                  forcex=forcex + &
                   rhophi_loc*(u_loc - u_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+   
		  forcey=forcey + &
		   rhophi_loc*(v_loc - v_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+ 
		  forcez=forcez + &
		   rhophi_loc*(w_loc - w_ref(i,j,k))*k_elastic*lambda_rel*phi_loc + rhophi_loc*fz !+  	
#endif

#ifdef DENSRATIO 
			  
                  pxx=real(hfields_s(ii,jj,kk,5,myblock),kind=db)
                  pyy=real(hfields_s(ii,jj,kk,6,myblock),kind=db)
                  pzz=real(hfields_s(ii,jj,kk,7,myblock),kind=db)
                  pxy=real(hfields_s(ii,jj,kk,8,myblock),kind=db)
                  pxz=real(hfields_s(ii,jj,kk,9,myblock),kind=db)
                  pyz=real(hfields_s(ii,jj,kk,10,myblock),kind=db)
                  
                  !1-2
                  !*1
                  ! 2nd order
                  pxx=pxx - u_loc*u_loc
                  pyy=pyy - v_loc*v_loc
                  pzz=pzz - w_loc*w_loc
                  pxy=pxy - u_loc*v_loc
                  pxz=pxz - u_loc*w_loc
                  pyz=pyz - v_loc*w_loc

		  visc_loc=(rho_r*visc1*phi_loc+(1.0_db-phi_loc)*visc2*rho_b)/rhophi_loc
				  
                  tau_loc=(visc_loc*invcssq + HALF) !this is a tau
				  
		  forcex=forcex - mymask* &
		   (visc_loc/(tau_loc*cssq))*(pxx*gradrhox + pxy*gradrhoy + pxz*gradrhoz)
		  forcey=forcey - mymask* &
		   (visc_loc/(tau_loc*cssq))*(pyy*gradrhoy + pxy*gradrhox + pyz*gradrhoz)
		  forcez=forcez - mymask* &
		   (visc_loc/(tau_loc*cssq))*(pzz*gradrhoz + pxz*gradrhox + pyz*gradrhoy)
#endif	
                  !I compute the new velocities
                  u_loc=real(hfields_s(ii,jj,kk,2,myblock),kind=db) !velocity
                  v_loc=real(hfields_s(ii,jj,kk,3,myblock),kind=db)
                  w_loc=real(hfields_s(ii,jj,kk,4,myblock),kind=db)  
					 
#if defined(INTERFACE_INCOMP) && defined(DENSRATIO)
                  mytemp= -sharp_c*real(locauxfields_s(ii,jj,kk,2,myblock),kind=db)

		  u_loc = u_loc + 0.5_db*forcex/(rhophi_loc)
                  v_loc = v_loc + 0.5_db*forcey/(rhophi_loc)
                  w_loc = w_loc + 0.5_db*forcez/(rhophi_loc)
				  
		  u_loc = u_loc/ &
		   (1.0_db - 0.5_db*(rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)/rhophi_loc )
				  
		  v_loc = v_loc/ &
		   (1.0_db - 0.5_db*(rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)/rhophi_loc ) 
				  
		  w_loc = w_loc/ &
	           (1.0_db - 0.5_db*(rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)/rhophi_loc ) 

		  tforcex= tforcex - &
		   (rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)*u_loc
		  tforcey= tforcey - &
		   (rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)*v_loc
		  tforcez= tforcez - &
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
                  

                  hfields_s(ii,jj,kk,2,myblock)=real(u_loc,kind=strdb) !put the new velocity in hfields_s
                  hfields_s(ii,jj,kk,3,myblock)=real(v_loc,kind=strdb)
                  hfields_s(ii,jj,kk,4,myblock)=real(w_loc,kind=strdb)
                  
                                    

!regularized 
                  pxx=real(hfields_s(ii,jj,kk,5,myblock),kind=db)
                  pyy=real(hfields_s(ii,jj,kk,6,myblock),kind=db)
                  pzz=real(hfields_s(ii,jj,kk,7,myblock),kind=db)
                  pxy=real(hfields_s(ii,jj,kk,8,myblock),kind=db)
                  pxz=real(hfields_s(ii,jj,kk,9,myblock),kind=db)
                  pyz=real(hfields_s(ii,jj,kk,10,myblock),kind=db)
                  
                  
                  pxx=pxx + forcex*u_loc/rhophi_loc
                  pyy=pyy + forcey*v_loc/rhophi_loc
                  pzz=pzz + forcez*w_loc/rhophi_loc
                  pxy=pxy + HALF*(forcey*u_loc+forcex*v_loc)/rhophi_loc
                  pxz=pxz + HALF*(forcez*u_loc+forcex*w_loc)/rhophi_loc
                  pyz=pyz + HALF*(forcez*v_loc+forcey*w_loc)/rhophi_loc


                  hfields_s(ii,jj,kk,5,myblock)=real(pxx,kind=strdb)
                  hfields_s(ii,jj,kk,6,myblock)=real(pyy,kind=strdb)
                  hfields_s(ii,jj,kk,7,myblock)=real(pzz,kind=strdb)
                  hfields_s(ii,jj,kk,8,myblock)=real(pxy,kind=strdb)
                  hfields_s(ii,jj,kk,9,myblock)=real(pxz,kind=strdb)
                  hfields_s(ii,jj,kk,10,myblock)=real(pyz,kind=strdb)
                  
                  
#if defined(ELASTIC_FORCE)
		  u_ref(i,j,k) = u_ref(i,j,k) + &
		   lambda_rel*(u_loc - u_ref(i,j,k))
		  v_ref(i,j,k) = v_ref(i,j,k) + &
		   lambda_rel*(v_loc - v_ref(i,j,k))
		  w_ref(i,j,k) = w_ref(i,j,k) + &
		   lambda_rel*(w_loc - w_ref(i,j,k))
				  
                  tforcex= tforcex + &
                   rhophi_loc*(u_loc - u_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+ 
		  tforcey= tforcey +&
		   rhophi_loc*(v_loc - v_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+ 
		  tforcez= tforcez + &
		   rhophi_loc*(w_loc - w_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+rhophi_loc*fz  
#endif               
                  
#if defined(DENSRATIO)	
                  pxx=pxx - u_loc*u_loc 
                  pyy=pyy - v_loc*v_loc
                  pzz=pzz - w_loc*w_loc 
                  pxy=pxy - u_loc*v_loc
                  pxz=pxz - u_loc*w_loc
                  pyz=pyz - v_loc*w_loc
		  
		  tforcex= tforcex - mymask* &
		   (visc_loc/(tau_loc*cssq))*(pxx*gradrhox + pxy*gradrhoy + pxz*gradrhoz)
		  tforcey= tforcey - mymask* &
		   (visc_loc/(tau_loc*cssq))*(pyy*gradrhoy + pxy*gradrhox + pyz*gradrhoz)
		  tforcez= tforcez - mymask* &
		   (visc_loc/(tau_loc*cssq))*(pzz*gradrhoz + pxz*gradrhox + pyz*gradrhoy)
#endif               

                  forces_s(ii,jj,kk,1,myblock)=real(tforcex,kind=strdb)
		          forces_s(ii,jj,kk,2,myblock)=real(tforcey,kind=strdb)
		          forces_s(ii,jj,kk,3,myblock)=real(tforcez,kind=strdb)

   endsubroutine moments_LB_kernel_int
   
      ! moments_LB_kernel_xplus: as moments_LB_kernel, for the last tile layer along x (bx = nxblock_d-2),
      ! y and z tiles 2..n/TILE_DIM-1 (dimGridx); not launched when there is a single tile along x.
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
       ,hfields_old,hfields_s &
#ifdef TWOCOMPONENT        
       ,auxfields_s,locauxfields_s &
#endif
       ,forces_s)
 

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
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields,nblocks_d) :: phifields_s 
#endif           
#if defined(ELASTIC_FORCE)
      real(kind=db) :: lambda_rel,k_elastic
      real(kind=db), dimension(1:nx,1:ny,1:nz) :: u_ref,v_ref,w_ref
#endif
      real(kind=db) :: fx,fy,fz

      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields,nblocks_d) :: hfields_old,hfields_s
#ifdef TWOCOMPONENT 
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nauxfields,nblocks_d) :: auxfields_s
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nlocauxfields,nblocks_d) :: locauxfields_s
#endif
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nforces,nblocks_d) :: forces_s
      

      real(kind=db) :: pxx,pyy,pzz,pxy,pxz,pyz
      real(kind=db) :: mytemp,forcex,forcey,forcez,rhophi_loc,press_loc,u_loc,v_loc,w_loc
      real(kind=db) :: tforcex,tforcey,tforcez
      
#ifdef TWOCOMPONENT
	  real(kind=db) ::gradfix,gradfiy,gradfiz,wet_loc,phi_loc,lap_phi_loc,mymask
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
      
				 
		  press_loc=real(hfields_s(ii,jj,kk,1,myblock),kind=db)
			 
#ifdef TWOCOMPONENT
		  phi_loc=real(phifields_s(ii,jj,kk,1,myblock),kind=db)
		  lap_phi_loc=real(locauxfields_s(ii,jj,kk,1,myblock),kind=db)
!          mymask=ONE-HALF*(ONE-tanh((phi_loc-0.015_db)/(HALF*0.005_db)))
!          mymask=    HALF*(ONE-tanh((phi_loc-0.995_db)/(HALF*0.005_db)))*mymask
          mymask=sqrt(max(ZERO,FOUR*phi_loc*(ONE-phi_loc)))
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
		   
		  ! Selectable capillary force: geometric CSF or Jacqmin chemical potential.
		  mytemp=real(auxfields_s(ii,jj,kk,4,myblock),kind=db) !modgrad
		  gradfix=real(auxfields_s(ii,jj,kk,1,myblock),kind=db)*mytemp !normx*modgrad
		  gradfiy=real(auxfields_s(ii,jj,kk,2,myblock),kind=db)*mytemp !normy*modgrad
		  gradfiz=real(auxfields_s(ii,jj,kk,3,myblock),kind=db)*mytemp !normz*modgrad

#ifdef CSF
		  forcex=forcex+real(forces_s(ii,jj,kk,1,myblock),kind=db)
		  forcey=forcey+real(forces_s(ii,jj,kk,2,myblock),kind=db)
		  forcez=forcez+real(forces_s(ii,jj,kk,3,myblock),kind=db)
#else
		  forcex = forcex + &
                   (4.0_db*beta*phi_loc*(phi_loc-1.0_db)*(phi_loc-0.5_db) - kapp*lap_phi_loc)*gradfix
		  forcey = forcey + &
                   (4.0_db*beta*phi_loc*(phi_loc-1.0_db)*(phi_loc-0.5_db) - kapp*lap_phi_loc)*gradfiy
		  forcez = forcez + &
                   (4.0_db*beta*phi_loc*(phi_loc-1.0_db)*(phi_loc-0.5_db) - kapp*lap_phi_loc)*gradfiz
#endif
				   				   

#ifdef REPULSIVE_FLUX
                  mytemp=real(locauxfields_s(ii,jj,kk,6,myblock),kind=db)*rhophi_loc 
                  if(abs(mytemp)>1.0d-3) mytemp=1.0d-3*sign(1.0,mytemp)!mytemp*0.1_db
                  forcex=forcex + mytemp*rhophi_loc
		  
                  mytemp=real(locauxfields_s(ii,jj,kk,7,myblock),kind=db)*rhophi_loc 
                  if(abs(mytemp)>1.0d-3) mytemp=1.0d-3*sign(1.0,mytemp)!mytemp*0.1_db
                  forcey=forcey + mytemp*rhophi_loc
		  
                  mytemp=real(locauxfields_s(ii,jj,kk,8,myblock),kind=db)*rhophi_loc 
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
				  
          forcex=forcex - mymask* &
           press_loc*cssq*gradrhox   !+ (rhophi_loc-(rho_r+rho_b)*HALF)*fx
          forcey=forcey - mymask* &
           press_loc*cssq*gradrhoy   !+ (rhophi_loc-(rho_r+rho_b)*HALF)*fy
		  forcez=forcez - mymask* &
		   press_loc*cssq*gradrhoz   !+ (rhophi_loc-(rho_r+rho_b)*HALF)*fz
		  !! from this point I compute the force terms that depend on the velocity
		  !! these terms should be not included in force arrays since they must be computed with the updated velocity
		  !! at the end of this subroutine
#endif
                  
              
                  
                  tforcex=forcex
                  tforcey=forcey
                  tforcez=forcez
				  
                  u_loc=real(hfields_old(ii,jj,kk,2,myblock),kind=db) !velocity
                  v_loc=real(hfields_old(ii,jj,kk,3,myblock),kind=db)
                  w_loc=real(hfields_old(ii,jj,kk,4,myblock),kind=db)
                  
                 
				  
#if defined(ELASTIC_FORCE)
                  forcex=forcex + &
                   rhophi_loc*(u_loc - u_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+   
		  forcey=forcey + &
		   rhophi_loc*(v_loc - v_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+ 
		  forcez=forcez + &
		   rhophi_loc*(w_loc - w_ref(i,j,k))*k_elastic*lambda_rel*phi_loc + rhophi_loc*fz !+  	
#endif

#ifdef DENSRATIO 
			  
                  pxx=real(hfields_s(ii,jj,kk,5,myblock),kind=db)
                  pyy=real(hfields_s(ii,jj,kk,6,myblock),kind=db)
                  pzz=real(hfields_s(ii,jj,kk,7,myblock),kind=db)
                  pxy=real(hfields_s(ii,jj,kk,8,myblock),kind=db)
                  pxz=real(hfields_s(ii,jj,kk,9,myblock),kind=db)
                  pyz=real(hfields_s(ii,jj,kk,10,myblock),kind=db)
                  
                  !1-2
                  !*1
                  ! 2nd order
                  pxx=pxx - u_loc*u_loc
                  pyy=pyy - v_loc*v_loc
                  pzz=pzz - w_loc*w_loc
                  pxy=pxy - u_loc*v_loc
                  pxz=pxz - u_loc*w_loc
                  pyz=pyz - v_loc*w_loc

		  visc_loc=(rho_r*visc1*phi_loc+(1.0_db-phi_loc)*visc2*rho_b)/rhophi_loc
				  
                  tau_loc=(visc_loc*invcssq + HALF) !this is a tau
				  
		  forcex=forcex - mymask* &
		   (visc_loc/(tau_loc*cssq))*(pxx*gradrhox + pxy*gradrhoy + pxz*gradrhoz)
		  forcey=forcey - mymask* &
		   (visc_loc/(tau_loc*cssq))*(pyy*gradrhoy + pxy*gradrhox + pyz*gradrhoz)
		  forcez=forcez - mymask* &
		   (visc_loc/(tau_loc*cssq))*(pzz*gradrhoz + pxz*gradrhox + pyz*gradrhoy)
#endif	
                  !I compute the new velocities
                  u_loc=real(hfields_s(ii,jj,kk,2,myblock),kind=db) !velocity
                  v_loc=real(hfields_s(ii,jj,kk,3,myblock),kind=db)
                  w_loc=real(hfields_s(ii,jj,kk,4,myblock),kind=db)  

					 
#if defined(INTERFACE_INCOMP) && defined(DENSRATIO)
                  mytemp= -sharp_c*real(locauxfields_s(ii,jj,kk,2,myblock),kind=db)

		  u_loc = u_loc + 0.5_db*forcex/(rhophi_loc)
                  v_loc = v_loc + 0.5_db*forcey/(rhophi_loc)
                  w_loc = w_loc + 0.5_db*forcez/(rhophi_loc)
				  
		  u_loc = u_loc/ &
		   (1.0_db - 0.5_db*(rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)/rhophi_loc )
				  
		  v_loc = v_loc/ &
		   (1.0_db - 0.5_db*(rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)/rhophi_loc ) 
				  
		  w_loc = w_loc/ &
	           (1.0_db - 0.5_db*(rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)/rhophi_loc ) 

		  tforcex= tforcex - &
		   (rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)*u_loc
		  tforcey= tforcey - &
		   (rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)*v_loc
		  tforcez= tforcez - &
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
                  

                  hfields_s(ii,jj,kk,2,myblock)=real(u_loc,kind=strdb) !put the new velocity in hfields_s
                  hfields_s(ii,jj,kk,3,myblock)=real(v_loc,kind=strdb)
                  hfields_s(ii,jj,kk,4,myblock)=real(w_loc,kind=strdb)
                  
                                    

!regularized 
                  pxx=real(hfields_s(ii,jj,kk,5,myblock),kind=db)
                  pyy=real(hfields_s(ii,jj,kk,6,myblock),kind=db)
                  pzz=real(hfields_s(ii,jj,kk,7,myblock),kind=db)
                  pxy=real(hfields_s(ii,jj,kk,8,myblock),kind=db)
                  pxz=real(hfields_s(ii,jj,kk,9,myblock),kind=db)
                  pyz=real(hfields_s(ii,jj,kk,10,myblock),kind=db)
                  
                  
                  pxx=pxx + forcex*u_loc/rhophi_loc
                  pyy=pyy + forcey*v_loc/rhophi_loc
                  pzz=pzz + forcez*w_loc/rhophi_loc
                  pxy=pxy + HALF*(forcey*u_loc+forcex*v_loc)/rhophi_loc
                  pxz=pxz + HALF*(forcez*u_loc+forcex*w_loc)/rhophi_loc
                  pyz=pyz + HALF*(forcez*v_loc+forcey*w_loc)/rhophi_loc


                  hfields_s(ii,jj,kk,5,myblock)=real(pxx,kind=strdb)
                  hfields_s(ii,jj,kk,6,myblock)=real(pyy,kind=strdb)
                  hfields_s(ii,jj,kk,7,myblock)=real(pzz,kind=strdb)
                  hfields_s(ii,jj,kk,8,myblock)=real(pxy,kind=strdb)
                  hfields_s(ii,jj,kk,9,myblock)=real(pxz,kind=strdb)
                  hfields_s(ii,jj,kk,10,myblock)=real(pyz,kind=strdb)
                  
                  
#if defined(ELASTIC_FORCE)
		  u_ref(i,j,k) = u_ref(i,j,k) + &
		   lambda_rel*(u_loc - u_ref(i,j,k))
		  v_ref(i,j,k) = v_ref(i,j,k) + &
		   lambda_rel*(v_loc - v_ref(i,j,k))
		  w_ref(i,j,k) = w_ref(i,j,k) + &
		   lambda_rel*(w_loc - w_ref(i,j,k))
				  
          tforcex= tforcex + &
           rhophi_loc*(u_loc - u_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+ 
		  tforcey= tforcey +&
		   rhophi_loc*(v_loc - v_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+ 
		  tforcez= tforcez + &
		   rhophi_loc*(w_loc - w_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+rhophi_loc*fz  
#endif               
                  
#if defined(DENSRATIO)	
                  pxx=pxx - u_loc*u_loc 
                  pyy=pyy - v_loc*v_loc
                  pzz=pzz - w_loc*w_loc 
                  pxy=pxy - u_loc*v_loc
                  pxz=pxz - u_loc*w_loc
                  pyz=pyz - v_loc*w_loc
		  
		  tforcex= tforcex - mymask* &
		   (visc_loc/(tau_loc*cssq))*(pxx*gradrhox + pxy*gradrhoy + pxz*gradrhoz)
		  tforcey= tforcey - mymask* &
		   (visc_loc/(tau_loc*cssq))*(pyy*gradrhoy + pxy*gradrhox + pyz*gradrhoz)
		  tforcez= tforcez - mymask* &
		   (visc_loc/(tau_loc*cssq))*(pzz*gradrhoz + pxz*gradrhox + pyz*gradrhoy)
#endif               

                  forces_s(ii,jj,kk,1,myblock)=real(tforcex,kind=strdb)
		          forces_s(ii,jj,kk,2,myblock)=real(tforcey,kind=strdb)
		          forces_s(ii,jj,kk,3,myblock)=real(tforcez,kind=strdb)

   endsubroutine moments_LB_kernel_xplus
   
   ! moments_LB_kernel_xminus: as moments_LB_kernel, for the first tile layer along x (bx = 1),
   ! y and z tiles 2..n/TILE_DIM-1 (dimGridx).
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
       ,hfields_old,hfields_s &
#ifdef TWOCOMPONENT        
       ,auxfields_s,locauxfields_s &
#endif
       ,forces_s)
 

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
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields,nblocks_d) :: phifields_s 
#endif           
#if defined(ELASTIC_FORCE)
      real(kind=db) :: lambda_rel,k_elastic
      real(kind=db), dimension(1:nx,1:ny,1:nz) :: u_ref,v_ref,w_ref
#endif
      real(kind=db) :: fx,fy,fz

      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields,nblocks_d) :: hfields_old,hfields_s
#ifdef TWOCOMPONENT 
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nauxfields,nblocks_d) :: auxfields_s
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nlocauxfields,nblocks_d) :: locauxfields_s
#endif
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nforces,nblocks_d) :: forces_s
      

      real(kind=db) :: pxx,pyy,pzz,pxy,pxz,pyz
      real(kind=db) :: mytemp,forcex,forcey,forcez,rhophi_loc,press_loc,u_loc,v_loc,w_loc
      real(kind=db) :: tforcex,tforcey,tforcez
      
#ifdef TWOCOMPONENT
	  real(kind=db) ::gradfix,gradfiy,gradfiz,wet_loc,phi_loc,lap_phi_loc,mymask
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
      
				 
		  press_loc=real(hfields_s(ii,jj,kk,1,myblock),kind=db)
			 
#ifdef TWOCOMPONENT
		  phi_loc=real(phifields_s(ii,jj,kk,1,myblock),kind=db)
		  lap_phi_loc=real(locauxfields_s(ii,jj,kk,1,myblock),kind=db)
!          mymask=ONE-HALF*(ONE-tanh((phi_loc-0.015_db)/(HALF*0.005_db)))
!          mymask=    HALF*(ONE-tanh((phi_loc-0.995_db)/(HALF*0.005_db)))*mymask
          mymask=sqrt(max(ZERO,FOUR*phi_loc*(ONE-phi_loc)))
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
		   
		  ! Selectable capillary force: geometric CSF or Jacqmin chemical potential.
		  mytemp=real(auxfields_s(ii,jj,kk,4,myblock),kind=db) !modgrad
		  gradfix=real(auxfields_s(ii,jj,kk,1,myblock),kind=db)*mytemp !normx*modgrad
		  gradfiy=real(auxfields_s(ii,jj,kk,2,myblock),kind=db)*mytemp !normy*modgrad
		  gradfiz=real(auxfields_s(ii,jj,kk,3,myblock),kind=db)*mytemp !normz*modgrad

#ifdef CSF
		  forcex=forcex+real(forces_s(ii,jj,kk,1,myblock),kind=db)
		  forcey=forcey+real(forces_s(ii,jj,kk,2,myblock),kind=db)
		  forcez=forcez+real(forces_s(ii,jj,kk,3,myblock),kind=db)
#else
		  forcex = forcex + &
                   (4.0_db*beta*phi_loc*(phi_loc-1.0_db)*(phi_loc-0.5_db) - kapp*lap_phi_loc)*gradfix
		  forcey = forcey + &
                   (4.0_db*beta*phi_loc*(phi_loc-1.0_db)*(phi_loc-0.5_db) - kapp*lap_phi_loc)*gradfiy
		  forcez = forcez + &
                   (4.0_db*beta*phi_loc*(phi_loc-1.0_db)*(phi_loc-0.5_db) - kapp*lap_phi_loc)*gradfiz
#endif
				   				   

#ifdef REPULSIVE_FLUX
                  mytemp=real(locauxfields_s(ii,jj,kk,6,myblock),kind=db)*rhophi_loc 
                  if(abs(mytemp)>1.0d-3) mytemp=1.0d-3*sign(1.0,mytemp)!mytemp*0.1_db
                  forcex=forcex + mytemp*rhophi_loc
		  
                  mytemp=real(locauxfields_s(ii,jj,kk,7,myblock),kind=db)*rhophi_loc 
                  if(abs(mytemp)>1.0d-3) mytemp=1.0d-3*sign(1.0,mytemp)!mytemp*0.1_db
                  forcey=forcey + mytemp*rhophi_loc
		  
                  mytemp=real(locauxfields_s(ii,jj,kk,8,myblock),kind=db)*rhophi_loc 
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
				  
          forcex=forcex - mymask* &
           press_loc*cssq*gradrhox   !+ (rhophi_loc-(rho_r+rho_b)*HALF)*fx
          forcey=forcey - mymask* &
           press_loc*cssq*gradrhoy   !+ (rhophi_loc-(rho_r+rho_b)*HALF)*fy
		  forcez=forcez - mymask* &
		   press_loc*cssq*gradrhoz   !+ (rhophi_loc-(rho_r+rho_b)*HALF)*fz
		  !! from this point I compute the force terms that depend on the velocity
		  !! these terms should be not included in force arrays since they must be computed with the updated velocity
		  !! at the end of this subroutine
#endif
                  
              
                  
                  tforcex=forcex
                  tforcey=forcey
                  tforcez=forcez
				  
                  u_loc=real(hfields_old(ii,jj,kk,2,myblock),kind=db) !velocity
                  v_loc=real(hfields_old(ii,jj,kk,3,myblock),kind=db)
                  w_loc=real(hfields_old(ii,jj,kk,4,myblock),kind=db)
                  
                 
				  
#if defined(ELASTIC_FORCE)
                  forcex=forcex + &
                   rhophi_loc*(u_loc - u_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+   
		  forcey=forcey + &
		   rhophi_loc*(v_loc - v_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+ 
		  forcez=forcez + &
		   rhophi_loc*(w_loc - w_ref(i,j,k))*k_elastic*lambda_rel*phi_loc + rhophi_loc*fz !+  	
#endif

#ifdef DENSRATIO 
			  
                  pxx=real(hfields_s(ii,jj,kk,5,myblock),kind=db)
                  pyy=real(hfields_s(ii,jj,kk,6,myblock),kind=db)
                  pzz=real(hfields_s(ii,jj,kk,7,myblock),kind=db)
                  pxy=real(hfields_s(ii,jj,kk,8,myblock),kind=db)
                  pxz=real(hfields_s(ii,jj,kk,9,myblock),kind=db)
                  pyz=real(hfields_s(ii,jj,kk,10,myblock),kind=db)
                  
                  !1-2
                  !*1
                  ! 2nd order
                  pxx=pxx - u_loc*u_loc
                  pyy=pyy - v_loc*v_loc
                  pzz=pzz - w_loc*w_loc
                  pxy=pxy - u_loc*v_loc
                  pxz=pxz - u_loc*w_loc
                  pyz=pyz - v_loc*w_loc

		  visc_loc=(rho_r*visc1*phi_loc+(1.0_db-phi_loc)*visc2*rho_b)/rhophi_loc
				  
                  tau_loc=(visc_loc*invcssq + HALF) !this is a tau
				  
		  forcex=forcex - mymask* &
		   (visc_loc/(tau_loc*cssq))*(pxx*gradrhox + pxy*gradrhoy + pxz*gradrhoz)
		  forcey=forcey - mymask* &
		   (visc_loc/(tau_loc*cssq))*(pyy*gradrhoy + pxy*gradrhox + pyz*gradrhoz)
		  forcez=forcez - mymask* &
		   (visc_loc/(tau_loc*cssq))*(pzz*gradrhoz + pxz*gradrhox + pyz*gradrhoy)
#endif	
                  !I compute the new velocities
                  u_loc=real(hfields_s(ii,jj,kk,2,myblock),kind=db) !velocity
                  v_loc=real(hfields_s(ii,jj,kk,3,myblock),kind=db)
                  w_loc=real(hfields_s(ii,jj,kk,4,myblock),kind=db)  

					 
#if defined(INTERFACE_INCOMP) && defined(DENSRATIO)
                  mytemp= -sharp_c*real(locauxfields_s(ii,jj,kk,2,myblock),kind=db)
                  
		  u_loc = u_loc + 0.5_db*forcex/(rhophi_loc)
                  v_loc = v_loc + 0.5_db*forcey/(rhophi_loc)
                  w_loc = w_loc + 0.5_db*forcez/(rhophi_loc)
				  
		  u_loc = u_loc/ &
		   (1.0_db - 0.5_db*(rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)/rhophi_loc )
				  
		  v_loc = v_loc/ &
		   (1.0_db - 0.5_db*(rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)/rhophi_loc ) 
				  
		  w_loc = w_loc/ &
	           (1.0_db - 0.5_db*(rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)/rhophi_loc ) 

		  tforcex= tforcex - &
		   (rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)*u_loc
		  tforcey= tforcey - &
		   (rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)*v_loc
		  tforcez= tforcez - &
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
                  

                  hfields_s(ii,jj,kk,2,myblock)=real(u_loc,kind=strdb) !put the new velocity in hfields_s
                  hfields_s(ii,jj,kk,3,myblock)=real(v_loc,kind=strdb)
                  hfields_s(ii,jj,kk,4,myblock)=real(w_loc,kind=strdb)
                  
                                    

!regularized 
                  pxx=real(hfields_s(ii,jj,kk,5,myblock),kind=db)
                  pyy=real(hfields_s(ii,jj,kk,6,myblock),kind=db)
                  pzz=real(hfields_s(ii,jj,kk,7,myblock),kind=db)
                  pxy=real(hfields_s(ii,jj,kk,8,myblock),kind=db)
                  pxz=real(hfields_s(ii,jj,kk,9,myblock),kind=db)
                  pyz=real(hfields_s(ii,jj,kk,10,myblock),kind=db)
                  
                  
                  pxx=pxx + forcex*u_loc/rhophi_loc
                  pyy=pyy + forcey*v_loc/rhophi_loc
                  pzz=pzz + forcez*w_loc/rhophi_loc
                  pxy=pxy + HALF*(forcey*u_loc+forcex*v_loc)/rhophi_loc
                  pxz=pxz + HALF*(forcez*u_loc+forcex*w_loc)/rhophi_loc
                  pyz=pyz + HALF*(forcez*v_loc+forcey*w_loc)/rhophi_loc


                  hfields_s(ii,jj,kk,5,myblock)=real(pxx,kind=strdb)
                  hfields_s(ii,jj,kk,6,myblock)=real(pyy,kind=strdb)
                  hfields_s(ii,jj,kk,7,myblock)=real(pzz,kind=strdb)
                  hfields_s(ii,jj,kk,8,myblock)=real(pxy,kind=strdb)
                  hfields_s(ii,jj,kk,9,myblock)=real(pxz,kind=strdb)
                  hfields_s(ii,jj,kk,10,myblock)=real(pyz,kind=strdb)
                  
                  
#if defined(ELASTIC_FORCE)
		  u_ref(i,j,k) = u_ref(i,j,k) + &
		   lambda_rel*(u_loc - u_ref(i,j,k))
		  v_ref(i,j,k) = v_ref(i,j,k) + &
		   lambda_rel*(v_loc - v_ref(i,j,k))
		  w_ref(i,j,k) = w_ref(i,j,k) + &
		   lambda_rel*(w_loc - w_ref(i,j,k))
				  
                  tforcex= tforcex + &
                   rhophi_loc*(u_loc - u_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+ 
		  tforcey= tforcey +&
		   rhophi_loc*(v_loc - v_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+ 
		  tforcez= tforcez + &
		   rhophi_loc*(w_loc - w_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+rhophi_loc*fz  
#endif               
                  
#if defined(DENSRATIO)	
                  pxx=pxx - u_loc*u_loc 
                  pyy=pyy - v_loc*v_loc
                  pzz=pzz - w_loc*w_loc 
                  pxy=pxy - u_loc*v_loc
                  pxz=pxz - u_loc*w_loc
                  pyz=pyz - v_loc*w_loc
		  
		  tforcex= tforcex - mymask* &
		   (visc_loc/(tau_loc*cssq))*(pxx*gradrhox + pxy*gradrhoy + pxz*gradrhoz)
		  tforcey= tforcey - mymask* &
		   (visc_loc/(tau_loc*cssq))*(pyy*gradrhoy + pxy*gradrhox + pyz*gradrhoz)
		  tforcez= tforcez - mymask* &
		   (visc_loc/(tau_loc*cssq))*(pzz*gradrhoz + pxz*gradrhox + pyz*gradrhoy)
#endif               

                  forces_s(ii,jj,kk,1,myblock)=real(tforcex,kind=strdb)
		          forces_s(ii,jj,kk,2,myblock)=real(tforcey,kind=strdb)
		          forces_s(ii,jj,kk,3,myblock)=real(tforcez,kind=strdb)

   endsubroutine moments_LB_kernel_xminus
   
         ! moments_LB_kernel_yplus: as moments_LB_kernel, for the last tile layer along y
         ! (by = nyblock_d-2), all x tiles, z tiles 2..nz/TILE_DIMz-1 (dimGridy); not launched when
         ! there is a single tile along y.
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
       ,hfields_old,hfields_s &
#ifdef TWOCOMPONENT        
       ,auxfields_s,locauxfields_s &
#endif
       ,forces_s)
 

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
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields,nblocks_d) :: phifields_s 
#endif           
#if defined(ELASTIC_FORCE)
      real(kind=db) :: lambda_rel,k_elastic
      real(kind=db), dimension(1:nx,1:ny,1:nz) :: u_ref,v_ref,w_ref
#endif
      real(kind=db) :: fx,fy,fz

      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields,nblocks_d) :: hfields_old,hfields_s
#ifdef TWOCOMPONENT 
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nauxfields,nblocks_d) :: auxfields_s
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nlocauxfields,nblocks_d) :: locauxfields_s
#endif
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nforces,nblocks_d) :: forces_s
      

      real(kind=db) :: pxx,pyy,pzz,pxy,pxz,pyz
      real(kind=db) :: mytemp,forcex,forcey,forcez,rhophi_loc,press_loc,u_loc,v_loc,w_loc
      real(kind=db) :: tforcex,tforcey,tforcez
      
#ifdef TWOCOMPONENT
	  real(kind=db) ::gradfix,gradfiy,gradfiz,wet_loc,phi_loc,lap_phi_loc,mymask
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
      
				 
		  press_loc=real(hfields_s(ii,jj,kk,1,myblock),kind=db)
			 
#ifdef TWOCOMPONENT
		  phi_loc=real(phifields_s(ii,jj,kk,1,myblock),kind=db)
		  lap_phi_loc=real(locauxfields_s(ii,jj,kk,1,myblock),kind=db)
!          mymask=ONE-HALF*(ONE-tanh((phi_loc-0.015_db)/(HALF*0.005_db)))
!          mymask=    HALF*(ONE-tanh((phi_loc-0.995_db)/(HALF*0.005_db)))*mymask
          mymask=sqrt(max(ZERO,FOUR*phi_loc*(ONE-phi_loc)))
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
		   
		  ! Selectable capillary force: geometric CSF or Jacqmin chemical potential.
		  mytemp=real(auxfields_s(ii,jj,kk,4,myblock),kind=db) !modgrad
		  gradfix=real(auxfields_s(ii,jj,kk,1,myblock),kind=db)*mytemp !normx*modgrad
		  gradfiy=real(auxfields_s(ii,jj,kk,2,myblock),kind=db)*mytemp !normy*modgrad
		  gradfiz=real(auxfields_s(ii,jj,kk,3,myblock),kind=db)*mytemp !normz*modgrad

#ifdef CSF
		  forcex=forcex+real(forces_s(ii,jj,kk,1,myblock),kind=db)
		  forcey=forcey+real(forces_s(ii,jj,kk,2,myblock),kind=db)
		  forcez=forcez+real(forces_s(ii,jj,kk,3,myblock),kind=db)
#else
		  forcex = forcex + &
                   (4.0_db*beta*phi_loc*(phi_loc-1.0_db)*(phi_loc-0.5_db) - kapp*lap_phi_loc)*gradfix
		  forcey = forcey + &
                   (4.0_db*beta*phi_loc*(phi_loc-1.0_db)*(phi_loc-0.5_db) - kapp*lap_phi_loc)*gradfiy
		  forcez = forcez + &
                   (4.0_db*beta*phi_loc*(phi_loc-1.0_db)*(phi_loc-0.5_db) - kapp*lap_phi_loc)*gradfiz
#endif
				   				   

#ifdef REPULSIVE_FLUX
                  mytemp=real(locauxfields_s(ii,jj,kk,6,myblock),kind=db)*rhophi_loc 
                  if(abs(mytemp)>1.0d-3) mytemp=1.0d-3*sign(1.0,mytemp)!mytemp*0.1_db
                  forcex=forcex + mytemp*rhophi_loc
		  
                  mytemp=real(locauxfields_s(ii,jj,kk,7,myblock),kind=db)*rhophi_loc 
                  if(abs(mytemp)>1.0d-3) mytemp=1.0d-3*sign(1.0,mytemp)!mytemp*0.1_db
                  forcey=forcey + mytemp*rhophi_loc
		  
                  mytemp=real(locauxfields_s(ii,jj,kk,8,myblock),kind=db)*rhophi_loc 
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
				  
          forcex=forcex - mymask* &
           press_loc*cssq*gradrhox   !+ (rhophi_loc-(rho_r+rho_b)*HALF)*fx
          forcey=forcey - mymask* &
           press_loc*cssq*gradrhoy   !+ (rhophi_loc-(rho_r+rho_b)*HALF)*fy
		  forcez=forcez - mymask* &
		   press_loc*cssq*gradrhoz   !+ (rhophi_loc-(rho_r+rho_b)*HALF)*fz
		  !! from this point I compute the force terms that depend on the velocity
		  !! these terms should be not included in force arrays since they must be computed with the updated velocity
		  !! at the end of this subroutine
#endif
                  
              
                  
                  tforcex=forcex
                  tforcey=forcey
                  tforcez=forcez
				  
                  u_loc=real(hfields_old(ii,jj,kk,2,myblock),kind=db) !velocity
                  v_loc=real(hfields_old(ii,jj,kk,3,myblock),kind=db)
                  w_loc=real(hfields_old(ii,jj,kk,4,myblock),kind=db)
                  
                 
				  
#if defined(ELASTIC_FORCE)
                  forcex=forcex + &
                   rhophi_loc*(u_loc - u_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+   
		  forcey=forcey + &
		   rhophi_loc*(v_loc - v_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+ 
		  forcez=forcez + &
		   rhophi_loc*(w_loc - w_ref(i,j,k))*k_elastic*lambda_rel*phi_loc + rhophi_loc*fz !+  	
#endif

#ifdef DENSRATIO 
                  pxx=real(hfields_s(ii,jj,kk,5,myblock),kind=db)
                  pyy=real(hfields_s(ii,jj,kk,6,myblock),kind=db)
                  pzz=real(hfields_s(ii,jj,kk,7,myblock),kind=db)
                  pxy=real(hfields_s(ii,jj,kk,8,myblock),kind=db)
                  pxz=real(hfields_s(ii,jj,kk,9,myblock),kind=db)
                  pyz=real(hfields_s(ii,jj,kk,10,myblock),kind=db)
                  
                  !1-2
                  !*1
                  ! 2nd order
                  pxx=pxx - u_loc*u_loc
                  pyy=pyy - v_loc*v_loc
                  pzz=pzz - w_loc*w_loc
                  pxy=pxy - u_loc*v_loc
                  pxz=pxz - u_loc*w_loc
                  pyz=pyz - v_loc*w_loc

		  visc_loc=(rho_r*visc1*phi_loc+(1.0_db-phi_loc)*visc2*rho_b)/rhophi_loc
				  
                  tau_loc=(visc_loc*invcssq + HALF) !this is a tau
				  
		  forcex=forcex - mymask* &
		   (visc_loc/(tau_loc*cssq))*(pxx*gradrhox + pxy*gradrhoy + pxz*gradrhoz)
		  forcey=forcey - mymask* &
		   (visc_loc/(tau_loc*cssq))*(pyy*gradrhoy + pxy*gradrhox + pyz*gradrhoz)
		  forcez=forcez - mymask* &
		   (visc_loc/(tau_loc*cssq))*(pzz*gradrhoz + pxz*gradrhox + pyz*gradrhoy)
#endif	
                  !I compute the new velocities
                  u_loc=real(hfields_s(ii,jj,kk,2,myblock),kind=db) !velocity
                  v_loc=real(hfields_s(ii,jj,kk,3,myblock),kind=db)
                  w_loc=real(hfields_s(ii,jj,kk,4,myblock),kind=db)  

					 
#if defined(INTERFACE_INCOMP) && defined(DENSRATIO)
                  mytemp= -sharp_c*real(locauxfields_s(ii,jj,kk,2,myblock),kind=db)
                  
		  u_loc = u_loc + 0.5_db*forcex/(rhophi_loc)
                  v_loc = v_loc + 0.5_db*forcey/(rhophi_loc)
                  w_loc = w_loc + 0.5_db*forcez/(rhophi_loc)
				  
		  u_loc = u_loc/ &
		   (1.0_db - 0.5_db*(rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)/rhophi_loc )
				  
		  v_loc = v_loc/ &
		   (1.0_db - 0.5_db*(rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)/rhophi_loc ) 
				  
		  w_loc = w_loc/ &
	           (1.0_db - 0.5_db*(rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)/rhophi_loc ) 

		  tforcex= tforcex - &
		   (rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)*u_loc
		  tforcey= tforcey - &
		   (rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)*v_loc
		  tforcez= tforcez - &
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
                  

                  hfields_s(ii,jj,kk,2,myblock)=real(u_loc,kind=strdb) !put the new velocity in hfields_s
                  hfields_s(ii,jj,kk,3,myblock)=real(v_loc,kind=strdb)
                  hfields_s(ii,jj,kk,4,myblock)=real(w_loc,kind=strdb)
                  
                                    

!regularized 
                  pxx=real(hfields_s(ii,jj,kk,5,myblock),kind=db)
                  pyy=real(hfields_s(ii,jj,kk,6,myblock),kind=db)
                  pzz=real(hfields_s(ii,jj,kk,7,myblock),kind=db)
                  pxy=real(hfields_s(ii,jj,kk,8,myblock),kind=db)
                  pxz=real(hfields_s(ii,jj,kk,9,myblock),kind=db)
                  pyz=real(hfields_s(ii,jj,kk,10,myblock),kind=db)
                  
                  
                  pxx=pxx + forcex*u_loc/rhophi_loc
                  pyy=pyy + forcey*v_loc/rhophi_loc
                  pzz=pzz + forcez*w_loc/rhophi_loc
                  pxy=pxy + HALF*(forcey*u_loc+forcex*v_loc)/rhophi_loc
                  pxz=pxz + HALF*(forcez*u_loc+forcex*w_loc)/rhophi_loc
                  pyz=pyz + HALF*(forcez*v_loc+forcey*w_loc)/rhophi_loc


                  hfields_s(ii,jj,kk,5,myblock)=real(pxx,kind=strdb)
                  hfields_s(ii,jj,kk,6,myblock)=real(pyy,kind=strdb)
                  hfields_s(ii,jj,kk,7,myblock)=real(pzz,kind=strdb)
                  hfields_s(ii,jj,kk,8,myblock)=real(pxy,kind=strdb)
                  hfields_s(ii,jj,kk,9,myblock)=real(pxz,kind=strdb)
                  hfields_s(ii,jj,kk,10,myblock)=real(pyz,kind=strdb)
                  
                  
#if defined(ELASTIC_FORCE)
		  u_ref(i,j,k) = u_ref(i,j,k) + &
		   lambda_rel*(u_loc - u_ref(i,j,k))
		  v_ref(i,j,k) = v_ref(i,j,k) + &
		   lambda_rel*(v_loc - v_ref(i,j,k))
		  w_ref(i,j,k) = w_ref(i,j,k) + &
		   lambda_rel*(w_loc - w_ref(i,j,k))
				  
                  tforcex= tforcex + &
                   rhophi_loc*(u_loc - u_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+ 
		  tforcey= tforcey +&
		   rhophi_loc*(v_loc - v_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+ 
		  tforcez= tforcez + &
		   rhophi_loc*(w_loc - w_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+rhophi_loc*fz  
#endif               
                  
#if defined(DENSRATIO)	
                  pxx=pxx - u_loc*u_loc 
                  pyy=pyy - v_loc*v_loc
                  pzz=pzz - w_loc*w_loc 
                  pxy=pxy - u_loc*v_loc
                  pxz=pxz - u_loc*w_loc
                  pyz=pyz - v_loc*w_loc
		  
		  tforcex= tforcex - mymask* &
		   (visc_loc/(tau_loc*cssq))*(pxx*gradrhox + pxy*gradrhoy + pxz*gradrhoz)
		  tforcey= tforcey - mymask* &
		   (visc_loc/(tau_loc*cssq))*(pyy*gradrhoy + pxy*gradrhox + pyz*gradrhoz)
		  tforcez= tforcez - mymask* &
		   (visc_loc/(tau_loc*cssq))*(pzz*gradrhoz + pxz*gradrhox + pyz*gradrhoy)
#endif               

                  forces_s(ii,jj,kk,1,myblock)=real(tforcex,kind=strdb)
		          forces_s(ii,jj,kk,2,myblock)=real(tforcey,kind=strdb)
		          forces_s(ii,jj,kk,3,myblock)=real(tforcez,kind=strdb)

   endsubroutine moments_LB_kernel_yplus
   
   ! moments_LB_kernel_yminus: as moments_LB_kernel, for the first tile layer along y (by = 1),
   ! all x tiles, z tiles 2..nz/TILE_DIMz-1 (dimGridy).
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
       ,hfields_old,hfields_s &
#ifdef TWOCOMPONENT        
       ,auxfields_s,locauxfields_s &
#endif
       ,forces_s)
 

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
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields,nblocks_d) :: phifields_s 
#endif           
#if defined(ELASTIC_FORCE)
      real(kind=db) :: lambda_rel,k_elastic
      real(kind=db), dimension(1:nx,1:ny,1:nz) :: u_ref,v_ref,w_ref
#endif
      real(kind=db) :: fx,fy,fz

      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields,nblocks_d) :: hfields_old,hfields_s
#ifdef TWOCOMPONENT 
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nauxfields,nblocks_d) :: auxfields_s
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nlocauxfields,nblocks_d) :: locauxfields_s
#endif
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nforces,nblocks_d) :: forces_s
      

      real(kind=db) :: pxx,pyy,pzz,pxy,pxz,pyz
      real(kind=db) :: mytemp,forcex,forcey,forcez,rhophi_loc,press_loc,u_loc,v_loc,w_loc
      real(kind=db) :: tforcex,tforcey,tforcez
      
#ifdef TWOCOMPONENT
	  real(kind=db) ::gradfix,gradfiy,gradfiz,wet_loc,phi_loc,lap_phi_loc,mymask
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
      
				 
		  press_loc=real(hfields_s(ii,jj,kk,1,myblock),kind=db)
			 
#ifdef TWOCOMPONENT
		  phi_loc=real(phifields_s(ii,jj,kk,1,myblock),kind=db)
		  lap_phi_loc=real(locauxfields_s(ii,jj,kk,1,myblock),kind=db)
!          mymask=ONE-HALF*(ONE-tanh((phi_loc-0.015_db)/(HALF*0.005_db)))
!          mymask=    HALF*(ONE-tanh((phi_loc-0.995_db)/(HALF*0.005_db)))*mymask
          mymask=sqrt(max(ZERO,FOUR*phi_loc*(ONE-phi_loc)))
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
		   
		  ! Selectable capillary force: geometric CSF or Jacqmin chemical potential.
		  mytemp=real(auxfields_s(ii,jj,kk,4,myblock),kind=db) !modgrad
		  gradfix=real(auxfields_s(ii,jj,kk,1,myblock),kind=db)*mytemp !normx*modgrad
		  gradfiy=real(auxfields_s(ii,jj,kk,2,myblock),kind=db)*mytemp !normy*modgrad
		  gradfiz=real(auxfields_s(ii,jj,kk,3,myblock),kind=db)*mytemp !normz*modgrad

#ifdef CSF
		  forcex=forcex+real(forces_s(ii,jj,kk,1,myblock),kind=db)
		  forcey=forcey+real(forces_s(ii,jj,kk,2,myblock),kind=db)
		  forcez=forcez+real(forces_s(ii,jj,kk,3,myblock),kind=db)
#else
		  forcex = forcex + &
                   (4.0_db*beta*phi_loc*(phi_loc-1.0_db)*(phi_loc-0.5_db) - kapp*lap_phi_loc)*gradfix
		  forcey = forcey + &
                   (4.0_db*beta*phi_loc*(phi_loc-1.0_db)*(phi_loc-0.5_db) - kapp*lap_phi_loc)*gradfiy
		  forcez = forcez + &
                   (4.0_db*beta*phi_loc*(phi_loc-1.0_db)*(phi_loc-0.5_db) - kapp*lap_phi_loc)*gradfiz
#endif
				   				   

#ifdef REPULSIVE_FLUX
                  mytemp=real(locauxfields_s(ii,jj,kk,6,myblock),kind=db)*rhophi_loc 
                  if(abs(mytemp)>1.0d-3) mytemp=1.0d-3*sign(1.0,mytemp)!mytemp*0.1_db
                  forcex=forcex + mytemp*rhophi_loc
		  
                  mytemp=real(locauxfields_s(ii,jj,kk,7,myblock),kind=db)*rhophi_loc 
                  if(abs(mytemp)>1.0d-3) mytemp=1.0d-3*sign(1.0,mytemp)!mytemp*0.1_db
                  forcey=forcey + mytemp*rhophi_loc
		  
                  mytemp=real(locauxfields_s(ii,jj,kk,8,myblock),kind=db)*rhophi_loc 
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
				  
          forcex=forcex - mymask* &
           press_loc*cssq*gradrhox   !+ (rhophi_loc-(rho_r+rho_b)*HALF)*fx
          forcey=forcey - mymask* &
           press_loc*cssq*gradrhoy   !+ (rhophi_loc-(rho_r+rho_b)*HALF)*fy
		  forcez=forcez - mymask* &
		   press_loc*cssq*gradrhoz   !+ (rhophi_loc-(rho_r+rho_b)*HALF)*fz
		  !! from this point I compute the force terms that depend on the velocity
		  !! these terms should be not included in force arrays since they must be computed with the updated velocity
		  !! at the end of this subroutine
#endif
                  
              
                  
                  tforcex=forcex
                  tforcey=forcey
                  tforcez=forcez
				  
                  u_loc=real(hfields_old(ii,jj,kk,2,myblock),kind=db) !velocity
                  v_loc=real(hfields_old(ii,jj,kk,3,myblock),kind=db)
                  w_loc=real(hfields_old(ii,jj,kk,4,myblock),kind=db)
                  
                 
				  
#if defined(ELASTIC_FORCE)
                  forcex=forcex + &
                   rhophi_loc*(u_loc - u_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+   
		  forcey=forcey + &
		   rhophi_loc*(v_loc - v_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+ 
		  forcez=forcez + &
		   rhophi_loc*(w_loc - w_ref(i,j,k))*k_elastic*lambda_rel*phi_loc + rhophi_loc*fz !+  	
#endif

#ifdef DENSRATIO 
			  
                  pxx=real(hfields_s(ii,jj,kk,5,myblock),kind=db)
                  pyy=real(hfields_s(ii,jj,kk,6,myblock),kind=db)
                  pzz=real(hfields_s(ii,jj,kk,7,myblock),kind=db)
                  pxy=real(hfields_s(ii,jj,kk,8,myblock),kind=db)
                  pxz=real(hfields_s(ii,jj,kk,9,myblock),kind=db)
                  pyz=real(hfields_s(ii,jj,kk,10,myblock),kind=db)
                  
                  !1-2
                  !*1
                  ! 2nd order
                  pxx=pxx - u_loc*u_loc
                  pyy=pyy - v_loc*v_loc
                  pzz=pzz - w_loc*w_loc
                  pxy=pxy - u_loc*v_loc
                  pxz=pxz - u_loc*w_loc
                  pyz=pyz - v_loc*w_loc

		  visc_loc=(rho_r*visc1*phi_loc+(1.0_db-phi_loc)*visc2*rho_b)/rhophi_loc
				  
                  tau_loc=(visc_loc*invcssq + HALF) !this is a tau
				  
		  forcex=forcex - mymask* &
		   (visc_loc/(tau_loc*cssq))*(pxx*gradrhox + pxy*gradrhoy + pxz*gradrhoz)
		  forcey=forcey - mymask* &
		   (visc_loc/(tau_loc*cssq))*(pyy*gradrhoy + pxy*gradrhox + pyz*gradrhoz)
		  forcez=forcez - mymask* &
		   (visc_loc/(tau_loc*cssq))*(pzz*gradrhoz + pxz*gradrhox + pyz*gradrhoy)
#endif	
                  !I compute the new velocities
                  u_loc=real(hfields_s(ii,jj,kk,2,myblock),kind=db) !velocity
                  v_loc=real(hfields_s(ii,jj,kk,3,myblock),kind=db)
                  w_loc=real(hfields_s(ii,jj,kk,4,myblock),kind=db)  
					 
#if defined(INTERFACE_INCOMP) && defined(DENSRATIO)

                  mytemp= -sharp_c*real(locauxfields_s(ii,jj,kk,2,myblock),kind=db)
                  
		  u_loc = u_loc + 0.5_db*forcex/(rhophi_loc)
                  v_loc = v_loc + 0.5_db*forcey/(rhophi_loc)
                  w_loc = w_loc + 0.5_db*forcez/(rhophi_loc)
				  
		  u_loc = u_loc/ &
		   (1.0_db - 0.5_db*(rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)/rhophi_loc )
				  
		  v_loc = v_loc/ &
		   (1.0_db - 0.5_db*(rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)/rhophi_loc ) 
				  
		  w_loc = w_loc/ &
	           (1.0_db - 0.5_db*(rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)/rhophi_loc ) 

		  tforcex= tforcex - &
		   (rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)*u_loc
		  tforcey= tforcey - &
		   (rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)*v_loc
		  tforcez= tforcez - &
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
                  

                  hfields_s(ii,jj,kk,2,myblock)=real(u_loc,kind=strdb) !put the new velocity in hfields_s
                  hfields_s(ii,jj,kk,3,myblock)=real(v_loc,kind=strdb)
                  hfields_s(ii,jj,kk,4,myblock)=real(w_loc,kind=strdb)
                  
                                    

!regularized 
                  pxx=real(hfields_s(ii,jj,kk,5,myblock),kind=db)
                  pyy=real(hfields_s(ii,jj,kk,6,myblock),kind=db)
                  pzz=real(hfields_s(ii,jj,kk,7,myblock),kind=db)
                  pxy=real(hfields_s(ii,jj,kk,8,myblock),kind=db)
                  pxz=real(hfields_s(ii,jj,kk,9,myblock),kind=db)
                  pyz=real(hfields_s(ii,jj,kk,10,myblock),kind=db)
                  
                  
                  pxx=pxx + forcex*u_loc/rhophi_loc
                  pyy=pyy + forcey*v_loc/rhophi_loc
                  pzz=pzz + forcez*w_loc/rhophi_loc
                  pxy=pxy + HALF*(forcey*u_loc+forcex*v_loc)/rhophi_loc
                  pxz=pxz + HALF*(forcez*u_loc+forcex*w_loc)/rhophi_loc
                  pyz=pyz + HALF*(forcez*v_loc+forcey*w_loc)/rhophi_loc


                  hfields_s(ii,jj,kk,5,myblock)=real(pxx,kind=strdb)
                  hfields_s(ii,jj,kk,6,myblock)=real(pyy,kind=strdb)
                  hfields_s(ii,jj,kk,7,myblock)=real(pzz,kind=strdb)
                  hfields_s(ii,jj,kk,8,myblock)=real(pxy,kind=strdb)
                  hfields_s(ii,jj,kk,9,myblock)=real(pxz,kind=strdb)
                  hfields_s(ii,jj,kk,10,myblock)=real(pyz,kind=strdb)
                  
                  
#if defined(ELASTIC_FORCE)
		  u_ref(i,j,k) = u_ref(i,j,k) + &
		   lambda_rel*(u_loc - u_ref(i,j,k))
		  v_ref(i,j,k) = v_ref(i,j,k) + &
		   lambda_rel*(v_loc - v_ref(i,j,k))
		  w_ref(i,j,k) = w_ref(i,j,k) + &
		   lambda_rel*(w_loc - w_ref(i,j,k))
				  
                  tforcex= tforcex + &
                   rhophi_loc*(u_loc - u_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+ 
		  tforcey= tforcey +&
		   rhophi_loc*(v_loc - v_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+ 
		  tforcez= tforcez + &
		   rhophi_loc*(w_loc - w_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+rhophi_loc*fz  
#endif               
                  
#if defined(DENSRATIO)	
                  pxx=pxx - u_loc*u_loc 
                  pyy=pyy - v_loc*v_loc
                  pzz=pzz - w_loc*w_loc 
                  pxy=pxy - u_loc*v_loc
                  pxz=pxz - u_loc*w_loc
                  pyz=pyz - v_loc*w_loc
		  
		  tforcex= tforcex - mymask* &
		   (visc_loc/(tau_loc*cssq))*(pxx*gradrhox + pxy*gradrhoy + pxz*gradrhoz)
		  tforcey= tforcey - mymask* &
		   (visc_loc/(tau_loc*cssq))*(pyy*gradrhoy + pxy*gradrhox + pyz*gradrhoz)
		  tforcez= tforcez - mymask* &
		   (visc_loc/(tau_loc*cssq))*(pzz*gradrhoz + pxz*gradrhox + pyz*gradrhoy)
#endif               

                  forces_s(ii,jj,kk,1,myblock)=real(tforcex,kind=strdb)
		          forces_s(ii,jj,kk,2,myblock)=real(tforcey,kind=strdb)
		          forces_s(ii,jj,kk,3,myblock)=real(tforcez,kind=strdb)

   endsubroutine moments_LB_kernel_yminus
   
         ! moments_LB_kernel_zplus: as moments_LB_kernel, for the last tile layer along z
         ! (bz = nzblock_d-2), all x and y tiles (dimGridz); not launched when there is a single tile
         ! along z.
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
       ,hfields_old,hfields_s &
#ifdef TWOCOMPONENT        
       ,auxfields_s,locauxfields_s &
#endif
       ,forces_s)
 

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
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields,nblocks_d) :: phifields_s 
#endif           
#if defined(ELASTIC_FORCE)
      real(kind=db) :: lambda_rel,k_elastic
      real(kind=db), dimension(1:nx,1:ny,1:nz) :: u_ref,v_ref,w_ref
#endif
      real(kind=db) :: fx,fy,fz

      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields,nblocks_d) :: hfields_old,hfields_s
#ifdef TWOCOMPONENT 
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nauxfields,nblocks_d) :: auxfields_s
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nlocauxfields,nblocks_d) :: locauxfields_s
#endif
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nforces,nblocks_d) :: forces_s
      

      real(kind=db) :: pxx,pyy,pzz,pxy,pxz,pyz
      real(kind=db) :: mytemp,forcex,forcey,forcez,rhophi_loc,press_loc,u_loc,v_loc,w_loc
      real(kind=db) :: tforcex,tforcey,tforcez
      
#ifdef TWOCOMPONENT
	  real(kind=db) ::gradfix,gradfiy,gradfiz,wet_loc,phi_loc,lap_phi_loc,mymask
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
      
				 
		  press_loc=real(hfields_s(ii,jj,kk,1,myblock),kind=db)
			 
#ifdef TWOCOMPONENT
		  phi_loc=real(phifields_s(ii,jj,kk,1,myblock),kind=db)
		  lap_phi_loc=real(locauxfields_s(ii,jj,kk,1,myblock),kind=db)
!          mymask=ONE-HALF*(ONE-tanh((phi_loc-0.015_db)/(HALF*0.005_db)))
!          mymask=    HALF*(ONE-tanh((phi_loc-0.995_db)/(HALF*0.005_db)))*mymask
          mymask=sqrt(max(ZERO,FOUR*phi_loc*(ONE-phi_loc)))
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
		   
		  ! Selectable capillary force: geometric CSF or Jacqmin chemical potential.
		  mytemp=real(auxfields_s(ii,jj,kk,4,myblock),kind=db) !modgrad
		  gradfix=real(auxfields_s(ii,jj,kk,1,myblock),kind=db)*mytemp !normx*modgrad
		  gradfiy=real(auxfields_s(ii,jj,kk,2,myblock),kind=db)*mytemp !normy*modgrad
		  gradfiz=real(auxfields_s(ii,jj,kk,3,myblock),kind=db)*mytemp !normz*modgrad

#ifdef CSF
		  forcex=forcex+real(forces_s(ii,jj,kk,1,myblock),kind=db)
		  forcey=forcey+real(forces_s(ii,jj,kk,2,myblock),kind=db)
		  forcez=forcez+real(forces_s(ii,jj,kk,3,myblock),kind=db)
#else
		  forcex = forcex + &
                   (4.0_db*beta*phi_loc*(phi_loc-1.0_db)*(phi_loc-0.5_db) - kapp*lap_phi_loc)*gradfix
		  forcey = forcey + &
                   (4.0_db*beta*phi_loc*(phi_loc-1.0_db)*(phi_loc-0.5_db) - kapp*lap_phi_loc)*gradfiy
		  forcez = forcez + &
                   (4.0_db*beta*phi_loc*(phi_loc-1.0_db)*(phi_loc-0.5_db) - kapp*lap_phi_loc)*gradfiz
#endif
				   				   

#ifdef REPULSIVE_FLUX
                  mytemp=real(locauxfields_s(ii,jj,kk,6,myblock),kind=db)*rhophi_loc 
                  if(abs(mytemp)>1.0d-3) mytemp=1.0d-3*sign(1.0,mytemp)!mytemp*0.1_db
                  forcex=forcex + mytemp*rhophi_loc
		  
                  mytemp=real(locauxfields_s(ii,jj,kk,7,myblock),kind=db)*rhophi_loc 
                  if(abs(mytemp)>1.0d-3) mytemp=1.0d-3*sign(1.0,mytemp)!mytemp*0.1_db
                  forcey=forcey + mytemp*rhophi_loc
		  
                  mytemp=real(locauxfields_s(ii,jj,kk,8,myblock),kind=db)*rhophi_loc 
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
				  
          forcex=forcex - mymask* &
           press_loc*cssq*gradrhox   !+ (rhophi_loc-(rho_r+rho_b)*HALF)*fx
          forcey=forcey - mymask* &
           press_loc*cssq*gradrhoy   !+ (rhophi_loc-(rho_r+rho_b)*HALF)*fy
		  forcez=forcez - mymask* &
		   press_loc*cssq*gradrhoz   !+ (rhophi_loc-(rho_r+rho_b)*HALF)*fz
		  !! from this point I compute the force terms that depend on the velocity
		  !! these terms should be not included in force arrays since they must be computed with the updated velocity
		  !! at the end of this subroutine
#endif
                  
              
                  
                  tforcex=forcex
                  tforcey=forcey
                  tforcez=forcez
				  
                  u_loc=real(hfields_old(ii,jj,kk,2,myblock),kind=db) !velocity
                  v_loc=real(hfields_old(ii,jj,kk,3,myblock),kind=db)
                  w_loc=real(hfields_old(ii,jj,kk,4,myblock),kind=db)
                  
                 
				  
#if defined(ELASTIC_FORCE)
                  forcex=forcex + &
                   rhophi_loc*(u_loc - u_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+   
		  forcey=forcey + &
		   rhophi_loc*(v_loc - v_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+ 
		  forcez=forcez + &
		   rhophi_loc*(w_loc - w_ref(i,j,k))*k_elastic*lambda_rel*phi_loc + rhophi_loc*fz !+  	
#endif

#ifdef DENSRATIO 
			  
                  pxx=real(hfields_s(ii,jj,kk,5,myblock),kind=db)
                  pyy=real(hfields_s(ii,jj,kk,6,myblock),kind=db)
                  pzz=real(hfields_s(ii,jj,kk,7,myblock),kind=db)
                  pxy=real(hfields_s(ii,jj,kk,8,myblock),kind=db)
                  pxz=real(hfields_s(ii,jj,kk,9,myblock),kind=db)
                  pyz=real(hfields_s(ii,jj,kk,10,myblock),kind=db)
                  
                  !1-2
                  !*1
                  ! 2nd order
                  pxx=pxx - u_loc*u_loc
                  pyy=pyy - v_loc*v_loc
                  pzz=pzz - w_loc*w_loc
                  pxy=pxy - u_loc*v_loc
                  pxz=pxz - u_loc*w_loc
                  pyz=pyz - v_loc*w_loc

		  visc_loc=(rho_r*visc1*phi_loc+(1.0_db-phi_loc)*visc2*rho_b)/rhophi_loc
				  
                  tau_loc=(visc_loc*invcssq + HALF) !this is a tau
				  
		  forcex=forcex - mymask* &
		   (visc_loc/(tau_loc*cssq))*(pxx*gradrhox + pxy*gradrhoy + pxz*gradrhoz)
		  forcey=forcey - mymask* &
		   (visc_loc/(tau_loc*cssq))*(pyy*gradrhoy + pxy*gradrhox + pyz*gradrhoz)
		  forcez=forcez - mymask* &
		   (visc_loc/(tau_loc*cssq))*(pzz*gradrhoz + pxz*gradrhox + pyz*gradrhoy)
#endif	
                  !I compute the new velocities
                  u_loc=real(hfields_s(ii,jj,kk,2,myblock),kind=db) !velocity
                  v_loc=real(hfields_s(ii,jj,kk,3,myblock),kind=db)
                  w_loc=real(hfields_s(ii,jj,kk,4,myblock),kind=db)  

					 
#if defined(INTERFACE_INCOMP) && defined(DENSRATIO)

                  mytemp= -sharp_c*real(locauxfields_s(ii,jj,kk,2,myblock),kind=db)

		  u_loc = u_loc + 0.5_db*forcex/(rhophi_loc)
                  v_loc = v_loc + 0.5_db*forcey/(rhophi_loc)
                  w_loc = w_loc + 0.5_db*forcez/(rhophi_loc)
				  
		  u_loc = u_loc/ &
		   (1.0_db - 0.5_db*(rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)/rhophi_loc )
				  
		  v_loc = v_loc/ &
		   (1.0_db - 0.5_db*(rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)/rhophi_loc ) 
				  
		  w_loc = w_loc/ &
	           (1.0_db - 0.5_db*(rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)/rhophi_loc ) 

		  tforcex= tforcex - &
		   (rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)*u_loc
		  tforcey= tforcey - &
		   (rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)*v_loc
		  tforcez= tforcez - &
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
                  

                  hfields_s(ii,jj,kk,2,myblock)=real(u_loc,kind=strdb) !put the new velocity in hfields_s
                  hfields_s(ii,jj,kk,3,myblock)=real(v_loc,kind=strdb)
                  hfields_s(ii,jj,kk,4,myblock)=real(w_loc,kind=strdb)
                  
                                    

!regularized 
                  pxx=real(hfields_s(ii,jj,kk,5,myblock),kind=db)
                  pyy=real(hfields_s(ii,jj,kk,6,myblock),kind=db)
                  pzz=real(hfields_s(ii,jj,kk,7,myblock),kind=db)
                  pxy=real(hfields_s(ii,jj,kk,8,myblock),kind=db)
                  pxz=real(hfields_s(ii,jj,kk,9,myblock),kind=db)
                  pyz=real(hfields_s(ii,jj,kk,10,myblock),kind=db)
                  
                  
                  pxx=pxx + forcex*u_loc/rhophi_loc
                  pyy=pyy + forcey*v_loc/rhophi_loc
                  pzz=pzz + forcez*w_loc/rhophi_loc
                  pxy=pxy + HALF*(forcey*u_loc+forcex*v_loc)/rhophi_loc
                  pxz=pxz + HALF*(forcez*u_loc+forcex*w_loc)/rhophi_loc
                  pyz=pyz + HALF*(forcez*v_loc+forcey*w_loc)/rhophi_loc


                  hfields_s(ii,jj,kk,5,myblock)=real(pxx,kind=strdb)
                  hfields_s(ii,jj,kk,6,myblock)=real(pyy,kind=strdb)
                  hfields_s(ii,jj,kk,7,myblock)=real(pzz,kind=strdb)
                  hfields_s(ii,jj,kk,8,myblock)=real(pxy,kind=strdb)
                  hfields_s(ii,jj,kk,9,myblock)=real(pxz,kind=strdb)
                  hfields_s(ii,jj,kk,10,myblock)=real(pyz,kind=strdb) 
                  
                  
#if defined(ELASTIC_FORCE)
		  u_ref(i,j,k) = u_ref(i,j,k) + &
		   lambda_rel*(u_loc - u_ref(i,j,k))
		  v_ref(i,j,k) = v_ref(i,j,k) + &
		   lambda_rel*(v_loc - v_ref(i,j,k))
		  w_ref(i,j,k) = w_ref(i,j,k) + &
		   lambda_rel*(w_loc - w_ref(i,j,k))
				  
                  tforcex= tforcex + &
                   rhophi_loc*(u_loc - u_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+ 
		  tforcey= tforcey +&
		   rhophi_loc*(v_loc - v_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+ 
		  tforcez= tforcez + &
		   rhophi_loc*(w_loc - w_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+rhophi_loc*fz  
#endif               
                  
#if defined(DENSRATIO)	
                  pxx=pxx - u_loc*u_loc 
                  pyy=pyy - v_loc*v_loc
                  pzz=pzz - w_loc*w_loc 
                  pxy=pxy - u_loc*v_loc
                  pxz=pxz - u_loc*w_loc
                  pyz=pyz - v_loc*w_loc
		  
		  tforcex= tforcex - mymask* &
		   (visc_loc/(tau_loc*cssq))*(pxx*gradrhox + pxy*gradrhoy + pxz*gradrhoz)
		  tforcey= tforcey - mymask* &
		   (visc_loc/(tau_loc*cssq))*(pyy*gradrhoy + pxy*gradrhox + pyz*gradrhoz)
		  tforcez= tforcez - mymask* &
		   (visc_loc/(tau_loc*cssq))*(pzz*gradrhoz + pxz*gradrhox + pyz*gradrhoy)
#endif               

                  forces_s(ii,jj,kk,1,myblock)=real(tforcex,kind=strdb)
		          forces_s(ii,jj,kk,2,myblock)=real(tforcey,kind=strdb)
		          forces_s(ii,jj,kk,3,myblock)=real(tforcez,kind=strdb)

   endsubroutine moments_LB_kernel_zplus
   
   ! moments_LB_kernel_zminus: as moments_LB_kernel, for the first tile layer along z (bz = 1),
   ! all x and y tiles (dimGridz).
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
       ,hfields_old,hfields_s &
#ifdef TWOCOMPONENT        
       ,auxfields_s,locauxfields_s &
#endif
       ,forces_s)
 

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
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields,nblocks_d) :: phifields_s 
#endif           
#if defined(ELASTIC_FORCE)
      real(kind=db) :: lambda_rel,k_elastic
      real(kind=db), dimension(1:nx,1:ny,1:nz) :: u_ref,v_ref,w_ref
#endif
      real(kind=db) :: fx,fy,fz

      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields,nblocks_d) :: hfields_old,hfields_s
#ifdef TWOCOMPONENT 
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nauxfields,nblocks_d) :: auxfields_s
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nlocauxfields,nblocks_d) :: locauxfields_s
#endif
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nforces,nblocks_d) :: forces_s
      

      real(kind=db) :: pxx,pyy,pzz,pxy,pxz,pyz
      real(kind=db) :: mytemp,forcex,forcey,forcez,rhophi_loc,press_loc,u_loc,v_loc,w_loc
      real(kind=db) :: tforcex,tforcey,tforcez

#ifdef TWOCOMPONENT
	  real(kind=db) ::gradfix,gradfiy,gradfiz,wet_loc,phi_loc,lap_phi_loc,mymask
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
      
				 
		  press_loc=real(hfields_s(ii,jj,kk,1,myblock),kind=db)
			 
#ifdef TWOCOMPONENT
		  phi_loc=real(phifields_s(ii,jj,kk,1,myblock),kind=db)
		  lap_phi_loc=real(locauxfields_s(ii,jj,kk,1,myblock),kind=db)
!          mymask=ONE-HALF*(ONE-tanh((phi_loc-0.015_db)/(HALF*0.005_db)))
!          mymask=    HALF*(ONE-tanh((phi_loc-0.995_db)/(HALF*0.005_db)))*mymask
          mymask=sqrt(max(ZERO,FOUR*phi_loc*(ONE-phi_loc)))
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
		   
		  ! Selectable capillary force: geometric CSF or Jacqmin chemical potential.
		  mytemp=real(auxfields_s(ii,jj,kk,4,myblock),kind=db) !modgrad
		  gradfix=real(auxfields_s(ii,jj,kk,1,myblock),kind=db)*mytemp !normx*modgrad
		  gradfiy=real(auxfields_s(ii,jj,kk,2,myblock),kind=db)*mytemp !normy*modgrad
		  gradfiz=real(auxfields_s(ii,jj,kk,3,myblock),kind=db)*mytemp !normz*modgrad

#ifdef CSF
		  forcex=forcex+real(forces_s(ii,jj,kk,1,myblock),kind=db)
		  forcey=forcey+real(forces_s(ii,jj,kk,2,myblock),kind=db)
		  forcez=forcez+real(forces_s(ii,jj,kk,3,myblock),kind=db)
#else
		  forcex = forcex + &
                   (4.0_db*beta*phi_loc*(phi_loc-1.0_db)*(phi_loc-0.5_db) - kapp*lap_phi_loc)*gradfix
		  forcey = forcey + &
                   (4.0_db*beta*phi_loc*(phi_loc-1.0_db)*(phi_loc-0.5_db) - kapp*lap_phi_loc)*gradfiy
		  forcez = forcez + &
                   (4.0_db*beta*phi_loc*(phi_loc-1.0_db)*(phi_loc-0.5_db) - kapp*lap_phi_loc)*gradfiz
#endif
				   				   

#ifdef REPULSIVE_FLUX
                  mytemp=real(locauxfields_s(ii,jj,kk,6,myblock),kind=db)*rhophi_loc 
                  if(abs(mytemp)>1.0d-3) mytemp=1.0d-3*sign(1.0,mytemp)!mytemp*0.1_db
                  forcex=forcex + mytemp*rhophi_loc
		  
                  mytemp=real(locauxfields_s(ii,jj,kk,7,myblock),kind=db)*rhophi_loc 
                  if(abs(mytemp)>1.0d-3) mytemp=1.0d-3*sign(1.0,mytemp)!mytemp*0.1_db
                  forcey=forcey + mytemp*rhophi_loc
		  
                  mytemp=real(locauxfields_s(ii,jj,kk,8,myblock),kind=db)*rhophi_loc 
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
				  
          forcex=forcex - mymask* &
           press_loc*cssq*gradrhox   !+ (rhophi_loc-(rho_r+rho_b)*HALF)*fx
          forcey=forcey - mymask* &
           press_loc*cssq*gradrhoy   !+ (rhophi_loc-(rho_r+rho_b)*HALF)*fy
		  forcez=forcez - mymask* &
		   press_loc*cssq*gradrhoz   !+ (rhophi_loc-(rho_r+rho_b)*HALF)*fz
		  !! from this point I compute the force terms that depend on the velocity
		  !! these terms should be not included in force arrays since they must be computed with the updated velocity
		  !! at the end of this subroutine
#endif
                  
              
                  
                  tforcex=forcex
                  tforcey=forcey
                  tforcez=forcez
				  
                  u_loc=real(hfields_old(ii,jj,kk,2,myblock),kind=db) !velocity
                  v_loc=real(hfields_old(ii,jj,kk,3,myblock),kind=db)
                  w_loc=real(hfields_old(ii,jj,kk,4,myblock),kind=db)
                  
                 
				  
#if defined(ELASTIC_FORCE)
                  forcex=forcex + &
                   rhophi_loc*(u_loc - u_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+   
		  forcey=forcey + &
		   rhophi_loc*(v_loc - v_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+ 
		  forcez=forcez + &
		   rhophi_loc*(w_loc - w_ref(i,j,k))*k_elastic*lambda_rel*phi_loc + rhophi_loc*fz !+  	
#endif

#ifdef DENSRATIO 
			  
                  pxx=real(hfields_s(ii,jj,kk,5,myblock),kind=db)
                  pyy=real(hfields_s(ii,jj,kk,6,myblock),kind=db)
                  pzz=real(hfields_s(ii,jj,kk,7,myblock),kind=db)
                  pxy=real(hfields_s(ii,jj,kk,8,myblock),kind=db)
                  pxz=real(hfields_s(ii,jj,kk,9,myblock),kind=db)
                  pyz=real(hfields_s(ii,jj,kk,10,myblock),kind=db)
                  
                  !1-2
                  !*1
                  ! 2nd order
                  pxx=pxx - u_loc*u_loc
                  pyy=pyy - v_loc*v_loc
                  pzz=pzz - w_loc*w_loc
                  pxy=pxy - u_loc*v_loc
                  pxz=pxz - u_loc*w_loc
                  pyz=pyz - v_loc*w_loc

		  visc_loc=(rho_r*visc1*phi_loc+(1.0_db-phi_loc)*visc2*rho_b)/rhophi_loc
				  
                  tau_loc=(visc_loc*invcssq + HALF) !this is a tau
				  
		  forcex=forcex - mymask* &
		   (visc_loc/(tau_loc*cssq))*(pxx*gradrhox + pxy*gradrhoy + pxz*gradrhoz)
		  forcey=forcey - mymask* &
		   (visc_loc/(tau_loc*cssq))*(pyy*gradrhoy + pxy*gradrhox + pyz*gradrhoz)
		  forcez=forcez - mymask* &
		   (visc_loc/(tau_loc*cssq))*(pzz*gradrhoz + pxz*gradrhox + pyz*gradrhoy)
#endif	
                  !I compute the new velocities
                  u_loc=real(hfields_s(ii,jj,kk,2,myblock),kind=db) !velocity
                  v_loc=real(hfields_s(ii,jj,kk,3,myblock),kind=db)
                  w_loc=real(hfields_s(ii,jj,kk,4,myblock),kind=db)  
					 
#if defined(INTERFACE_INCOMP) && defined(DENSRATIO)

                  mytemp= -sharp_c*real(locauxfields_s(ii,jj,kk,2,myblock),kind=db)
                  
		  u_loc = u_loc + 0.5_db*forcex/(rhophi_loc)
                  v_loc = v_loc + 0.5_db*forcey/(rhophi_loc)
                  w_loc = w_loc + 0.5_db*forcez/(rhophi_loc)
				  
		  u_loc = u_loc/ &
		   (1.0_db - 0.5_db*(rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)/rhophi_loc )
				  
		  v_loc = v_loc/ &
		   (1.0_db - 0.5_db*(rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)/rhophi_loc ) 
				  
		  w_loc = w_loc/ &
	           (1.0_db - 0.5_db*(rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)/rhophi_loc ) 

		  tforcex= tforcex - &
		   (rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)*u_loc
		  tforcey= tforcey - &
		   (rho_r-rho_b)*(tau_diff*lap_phi_loc + mytemp)*v_loc
		  tforcez= tforcez - &
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
                  

                  hfields_s(ii,jj,kk,2,myblock)=real(u_loc,kind=strdb) !put the new velocity in hfields_s
                  hfields_s(ii,jj,kk,3,myblock)=real(v_loc,kind=strdb)
                  hfields_s(ii,jj,kk,4,myblock)=real(w_loc,kind=strdb)
                  
                                    

!regularized 
                  pxx=real(hfields_s(ii,jj,kk,5,myblock),kind=db)
                  pyy=real(hfields_s(ii,jj,kk,6,myblock),kind=db)
                  pzz=real(hfields_s(ii,jj,kk,7,myblock),kind=db)
                  pxy=real(hfields_s(ii,jj,kk,8,myblock),kind=db)
                  pxz=real(hfields_s(ii,jj,kk,9,myblock),kind=db)
                  pyz=real(hfields_s(ii,jj,kk,10,myblock),kind=db)
                  
                  
                  pxx=pxx + forcex*u_loc/rhophi_loc
                  pyy=pyy + forcey*v_loc/rhophi_loc
                  pzz=pzz + forcez*w_loc/rhophi_loc
                  pxy=pxy + HALF*(forcey*u_loc+forcex*v_loc)/rhophi_loc
                  pxz=pxz + HALF*(forcez*u_loc+forcex*w_loc)/rhophi_loc
                  pyz=pyz + HALF*(forcez*v_loc+forcey*w_loc)/rhophi_loc


                  hfields_s(ii,jj,kk,5,myblock)=real(pxx,kind=strdb)
                  hfields_s(ii,jj,kk,6,myblock)=real(pyy,kind=strdb)
                  hfields_s(ii,jj,kk,7,myblock)=real(pzz,kind=strdb)
                  hfields_s(ii,jj,kk,8,myblock)=real(pxy,kind=strdb)
                  hfields_s(ii,jj,kk,9,myblock)=real(pxz,kind=strdb)
                  hfields_s(ii,jj,kk,10,myblock)=real(pyz,kind=strdb)
                  
                  
#if defined(ELASTIC_FORCE)
		  u_ref(i,j,k) = u_ref(i,j,k) + &
		   lambda_rel*(u_loc - u_ref(i,j,k))
		  v_ref(i,j,k) = v_ref(i,j,k) + &
		   lambda_rel*(v_loc - v_ref(i,j,k))
		  w_ref(i,j,k) = w_ref(i,j,k) + &
		   lambda_rel*(w_loc - w_ref(i,j,k))
				  
                  tforcex= tforcex + &
                   rhophi_loc*(u_loc - u_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+ 
		  tforcey= tforcey +&
		   rhophi_loc*(v_loc - v_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+ 
		  tforcez= tforcez + &
		   rhophi_loc*(w_loc - w_ref(i,j,k))*k_elastic*lambda_rel*phi_loc !+rhophi_loc*fz  
#endif               
                  
#if defined(DENSRATIO)	
                  pxx=pxx - u_loc*u_loc 
                  pyy=pyy - v_loc*v_loc
                  pzz=pzz - w_loc*w_loc 
                  pxy=pxy - u_loc*v_loc
                  pxz=pxz - u_loc*w_loc
                  pyz=pyz - v_loc*w_loc
		  
		  tforcex= tforcex - mymask* &
		   (visc_loc/(tau_loc*cssq))*(pxx*gradrhox + pxy*gradrhoy + pxz*gradrhoz)
		  tforcey= tforcey - mymask* &
		   (visc_loc/(tau_loc*cssq))*(pyy*gradrhoy + pxy*gradrhox + pyz*gradrhoz)
		  tforcez= tforcez - mymask* &
		   (visc_loc/(tau_loc*cssq))*(pzz*gradrhoz + pxz*gradrhox + pyz*gradrhoy)
#endif               

                  forces_s(ii,jj,kk,1,myblock)=real(tforcex,kind=strdb)
		          forces_s(ii,jj,kk,2,myblock)=real(tforcey,kind=strdb)
		          forces_s(ii,jj,kk,3,myblock)=real(tforcez,kind=strdb)

   endsubroutine moments_LB_kernel_zminus
   
endmodule lb_cuda_moments
