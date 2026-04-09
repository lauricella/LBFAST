#include "defines.h"
module initial_condts

   use vars
   use mpi_template, only: coords,myoffset,sum_world_float,dostop,myrank
   use stat_module, only: Ekin0,Ekin, &
#if defined(LAMBTEST) && defined(TWOCOMPONENT)   
    lamb_A,lamb_visc, &
#endif
    probe_loc
   implicit none

contains

   subroutine initial_conditions_all
#ifdef MPI  
	  use mpi
#endif
      implicit none
      integer:: subchords(3)
      real(kind=db) :: dist1,dist2,sel1,sel2,dist
      real(kind=db),dimension(3) :: dist3d,dist3dout,invdim
      real(kind=db) :: fneq1,feq, rhophi_loc,visc_loc,lamb_visc_temp  
#if defined(POISEUILLE) || defined(TWOPOISEUILLE)        
      real(kind=db) :: H_pois,xc_pois,distabs, &
       coefL_pois,coefB_pois,coefR_pois,wleft,wright
#endif    
#ifdef LAMBTEST
      real(kind=db) :: myp2,xx,yy,zz,rr,costh,rloc,eta0  
      real(kind=db) :: myfreq,nrat,mu1,mu2,chi,myfreq_corr,myperiod
#endif 
#if defined(MULTIHIT)
	  real(kind=db) :: k_zero
#endif
#ifdef INTERNAL_OBSTACLES
	integer :: line_number
	integer :: val, n, count
	integer :: fid
	character(len=256) :: line
	integer, allocatable :: buffer(:),ierr
#endif
    
    integer :: xblock,yblock,zblock,myblock,ii,jj,kk,link_dist2
    integer :: oii,ojj,okk,link_id,nei_x,nei_y,nei_z
    integer :: oxblock,oyblock,ozblock,omyblock
    real(kind=db) :: tempphi,tempphi2,loc_u,loc_v,loc_w,loc_press,stdev, &
     link_scale
    real(db), parameter :: inv_sqrt2  = 0.70710678118654752440_db
    real(db), parameter :: inv_sqrt3  = 0.57735026918962576450_db
#if defined(TAYLORGREEN) && !defined(TWOCOMPONENT)     
    real(kind=db) :: c2x,c2y,c2z,x,y,z
#endif      

    
    invdim(1) = ONE/real(lx,kind=db)
    invdim(2) = ONE/real(ly,kind=db)
    invdim(3) = ONE/real(lz,kind=db)
       
      !*************************************initial conditions ************************
      
      do k = 1, nz
		   gk = nz*coords(3) + k
		   zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
		   do j = 1, ny
			  gj = ny*coords(2) + j
			  yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
			  do i = 1, nx
				 gi = nx*coords(1) + i
                 xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
      
                 myblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
                 ii=i-xblock*TILE_DIMx+2*TILE_DIMx
                 jj=j-yblock*TILE_DIMy+2*TILE_DIMy
                 kk=k-zblock*TILE_DIMz+2*TILE_DIMz               

                 hfields_flip(ii,jj,kk,1,myblock)=ZEROSTR
                 hfields_flip(ii,jj,kk,2,myblock)=ZEROSTR
                 hfields_flip(ii,jj,kk,3,myblock)=ZEROSTR
                 hfields_flip(ii,jj,kk,4,myblock)=ZEROSTR
                 hfields_flip(ii,jj,kk,5,myblock)=ZEROSTR
                 hfields_flip(ii,jj,kk,6,myblock)=ZEROSTR
                 hfields_flip(ii,jj,kk,7,myblock)=ZEROSTR
                 hfields_flip(ii,jj,kk,8,myblock)=ZEROSTR
                 hfields_flip(ii,jj,kk,9,myblock)=ZEROSTR
                 hfields_flip(ii,jj,kk,10,myblock)=ZEROSTR
                 

#ifdef TWOCOMPONENT
				 
				 phifields_flip(ii,jj,kk,1,myblock)=ZEROSTR
				 
				 auxfields(ii,jj,kk,1,myblock)=ZEROSTR
				 auxfields(ii,jj,kk,2,myblock)=ZEROSTR
				 auxfields(ii,jj,kk,3,myblock)=ZEROSTR
				 auxfields(ii,jj,kk,4,myblock)=ZEROSTR
#endif
              enddo
            enddo
        enddo

		do k = 1, nz
		   gk = nz*coords(3) + k
		   zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
		   do j = 1, ny
			  gj = ny*coords(2) + j
			  yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
			  do i = 1, nx
				 gi = nx*coords(1) + i
				 xblock=(i+2*TILE_DIMx-1)/TILE_DIMx

                 myblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
                 ii=i-xblock*TILE_DIMx+2*TILE_DIMx
                 jj=j-yblock*TILE_DIMy+2*TILE_DIMy
                 kk=k-zblock*TILE_DIMz+2*TILE_DIMz   
                 
                 loc_press=ZERO
                 loc_u=ZERO
                 loc_v=ZERO
                 loc_w=ZERO  
                 
				 ! distances to the two centers
				 dist1 = sqrt( (real(gi,kind=db) - ((float(lx) / TWO) - radius-7))**2 + &
							   (real(gj,kind=db) - float(ly)/two )**2 + &
							   (real(gk,kind=db) - float(lz)/two - radius + 4 )**2 )

				 dist2 = sqrt( (real(gi,kind=db) - ((float(lx) / TWO) + radius+6))**2 + &
							   (real(gj,kind=db) - float(ly)/two )**2 + &
							   (real(gk,kind=db) - float(lz)/two + radius - 3)**2 )

				 ! individual smooth indicator fields (φ ~ 0 inside, ~1 outside)
				 sel1 = 0.5_db + 0.5_db * tanh( (radius - dist1) / 2.0_db )
				 sel2 = 0.5_db + 0.5_db * tanh( (radius - dist2) / 2.0_db )

				 ! union of droplets: inside either -> take the larger (closer to 1)
				 
				 tempphi=max(sel1, sel2)
				 
#ifdef TWOCOMPONENT		
				 phifields_flip(ii,jj,kk,1,myblock)=real(tempphi,kind=strdb)               
#ifdef DENSRATIO
                 rhophi_loc=rho_r*tempphi+(ONE-tempphi)*rho_b
#else
                 rhophi_loc = 1.0_db
#endif	
#endif                
                tempphi2 = tempphi*(sigma*2.0_db)/radius/(rhophi_loc/3.0_db)
				loc_press = loc_press + tempphi2
				 
			    ! crisp velocity: uniform inside each core (phi~1 region)
				!if (sel1 > 0.1_db) then 
				loc_u = 0.02*sel1-0.02*sel2

                hfields_flip(ii,jj,kk,1,myblock)= real(loc_press,kind=strdb)
                hfields_flip(ii,jj,kk,2,myblock)=real(loc_u,kind=strdb)
                hfields_flip(ii,jj,kk,3,myblock)=real(loc_v,kind=strdb) 
                hfields_flip(ii,jj,kk,4,myblock)=real(loc_w,kind=strdb) 
                hfields_flip(ii,jj,kk,5,myblock)=real(loc_u*loc_u+cssq*loc_press,kind=strdb)
                hfields_flip(ii,jj,kk,6,myblock)=real(loc_v*loc_v+cssq*loc_press,kind=strdb)
                hfields_flip(ii,jj,kk,7,myblock)=real(loc_w*loc_w+cssq*loc_press,kind=strdb)
                hfields_flip(ii,jj,kk,8,myblock)=real(loc_u*loc_v,kind=strdb)
                hfields_flip(ii,jj,kk,9,myblock)=real(loc_u*loc_w,kind=strdb)
                hfields_flip(ii,jj,kk,10,myblock)=real(loc_v*loc_w,kind=strdb)                
                        
				!else if (sel2 > 0.1_db) then
				 ! u(i,j,k) = -0.01*sel2
				!else
				!  u(i,j,k) = 0.0_db
			    !endif
			  end do
		   end do
		end do
		


    

#if defined(DOBENCHMARK)


      do k = 1, nz
		   gk = nz*coords(3) + k
		   zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
		   do j = 1, ny
			  gj = ny*coords(2) + j
			  yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
			  do i = 1, nx
				 gi = nx*coords(1) + i
                 xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
      
                 myblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
                 ii=i-xblock*TILE_DIMx+2*TILE_DIMx
                 jj=j-yblock*TILE_DIMy+2*TILE_DIMy
                 kk=k-zblock*TILE_DIMz+2*TILE_DIMz         
                    

                 hfields_flip(ii,jj,kk,1,myblock)=ZEROSTR
                 hfields_flip(ii,jj,kk,2,myblock)=ZEROSTR
                 hfields_flip(ii,jj,kk,3,myblock)=ZEROSTR
                 hfields_flip(ii,jj,kk,4,myblock)=ZEROSTR
                 hfields_flip(ii,jj,kk,5,myblock)=ZEROSTR
                 hfields_flip(ii,jj,kk,6,myblock)=ZEROSTR
                 hfields_flip(ii,jj,kk,7,myblock)=ZEROSTR
                 hfields_flip(ii,jj,kk,8,myblock)=ZEROSTR
                 hfields_flip(ii,jj,kk,9,myblock)=ZEROSTR
                 hfields_flip(ii,jj,kk,10,myblock)=ZEROSTR
                 
                 
#ifdef TWOCOMPONENT
				 
				 phifields_flip(ii,jj,kk,1,myblock)=ZEROSTR
				 
				 auxfields(ii,jj,kk,1,myblock)=ZEROSTR
				 auxfields(ii,jj,kk,2,myblock)=ZEROSTR
				 auxfields(ii,jj,kk,3,myblock)=ZEROSTR
				 auxfields(ii,jj,kk,4,myblock)=ZEROSTR
#endif
              enddo
            enddo
        enddo

       ! ampiezza del TG: scegli tu, ma resta a Mach basso
       ! es. uwall = u0tg = 0.02_db oppure 0.05_db
       stdev=0.0e-2
       
#if defined(LAMBTEST) && defined(TWOCOMPONENT)
       lamb_A=uwall / radius
       lamb_visc=ZERO
       lamb_visc_temp=ZERO
       myfreq = sqrt(TWENTYFOUR*sigma / (radius**THREE * (TWO*rho_b + THREE*rho_r)))
	   nrat = rho_r/rho_b
       mu1  = rho_r*visc1
       mu2  = rho_b*visc2

       chi = ((TWO*nrat + ONE)**TWO * sqrt(mu1*mu2*rho_r*rho_b)) / &
        (TWO*radius * (nrat*rho_b + (nrat + ONE)*rho_r) * &
        (sqrt(mu1*rho_r) + sqrt(mu2*rho_b)))
	   myfreq_corr = myfreq - HALF*chi*sqrt(myfreq) + (ONE/FOUR)*chi**TWO

       myperiod = TWO*pi_greek / myfreq_corr
       
       if(myrank==0)then
          write(6,'(a,f20.10)')'LAMB: myfreq_corr',myfreq_corr
          write(6,'(a,f20.10)')'LAMB: myperiod',myperiod
       endif
#endif       
       
#if defined(POISEUILLE) || defined(TWOPOISEUILLE)
#ifdef BOUNCE_BACK
       H_pois  = 0.5_db * real(lx-2,db)
#else
       H_pois  = 0.5_db * real(lx-1,db)
#endif
       xc_pois = 0.5_db * real(lx+1,db)
#endif

#if defined(TWOCOMPONENT) && defined(TWOPOISEUILLE)
#ifdef DENSRATIO
       coefB_pois = fz * H_pois * H_pois * (rho_r + rho_b) / &
                   (2.0_db * (rho_r*visc1 + rho_b*visc2))

       coefL_pois = fz * H_pois * rho_b * (visc1 - visc2) / &
                   (2.0_db * visc1 * (rho_r*visc1 + rho_b*visc2))

       coefR_pois = fz * H_pois * rho_r * (visc1 - visc2) / &
                   (2.0_db * visc2 * (rho_r*visc1 + rho_b*visc2))
#else
       coefB_pois = fz * H_pois * H_pois / (visc1 + visc2)

       coefL_pois = fz * H_pois * (visc1 - visc2) / &
                   (2.0_db * visc1 * (visc1 + visc2))

       coefR_pois = fz * H_pois * (visc1 - visc2) / &
                   (2.0_db * visc2 * (visc1 + visc2))
#endif
#endif
       do k=1,nz
          gk = nz*coords(3) + k
		  zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
          do j=1,ny
             gj=ny*coords(2)+j
             yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
             do i=1,nx
                gi=nx*coords(1)+i
                xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
                
                myblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
                ii=i-xblock*TILE_DIMx+2*TILE_DIMx
                jj=j-yblock*TILE_DIMy+2*TILE_DIMy
                kk=k-zblock*TILE_DIMz+2*TILE_DIMz  
                 
                 loc_press=ZERO
                 loc_u=ZERO
                 loc_v=ZERO
                 loc_w=ZERO  
                
                
                if(abs(isfluid(i,j,k)).eq.1)then
#if defined(TAYLORGREEN) && !defined(TWOCOMPONENT)
                  rhophi_loc = 1.0_db
                  ! coordinate periodiche globali
                  x = TWO*pi_greek * real(gi-1, db) / real(lx, db)
                  y = TWO*pi_greek * real(gj-1, db) / real(ly, db)
                  z = TWO*pi_greek * real(gk-1, db) / real(lz, db)

                  loc_u =  uwall * sin(x) * cos(y) * cos(z)
                  loc_v = -uwall * cos(x) * sin(y) * cos(z)
                  loc_w =  ZERO
                  
                  c2x = cos(TWO*x)
                  c2y = cos(TWO*y)
                  c2z = cos(TWO*z)

                  ! pressione idrodinamica fluttuante, media nulla
                  loc_press = ((uwall*uwall/SIXTEEN) * (c2x + c2y) * (TWO + c2z))/(rhophi_loc*cssq)

#elif defined(TWOCOMPONENT) && defined(TWOPOISEUILLE)
               dist   = real(gi,db) - xc_pois
               distabs = abs(dist)
               tempphi=fcut_tanh(dist,ZERO,width)
               !tempphi = ONE - (0.5 + 0.5*TANH(2.0_db*(dist)/width))	
!               if (dist <= ZERO) then
!                  tempphi = ONE
!               else
!                  tempphi = ZERO
!               endif

#ifdef DENSRATIO
               rhophi_loc=rho_r*tempphi+(ONE-tempphi)*rho_b       
#endif

               if (distabs <= H_pois) then
                  !if (dist <= ZERO) then
                     wleft = -(fz/(2.0_db*visc1)) * dist*dist + &
                                 coefL_pois * dist + coefB_pois
                  !else
                     wright = -(fz/(2.0_db*visc2)) * dist*dist + &
                                 coefR_pois * dist + coefB_pois
                  !endif
                  loc_w=wleft*tempphi + (ONE-tempphi)*wright   
               else
                  loc_w = ZERO
               endif


#elif defined(POISEUILLE)
                  rhophi_loc = ONE     
                  dist = real(gi,db) - xc_pois
                  distabs = abs(dist)
                  if (distabs <= H_pois) then
                    loc_w=(fz)/(TWO*visc1) * (H_pois*H_pois - dist*dist)+ stdev*randgauss_CPU()
                  else
                    loc_w=ZERO
                  endif
                  loc_press=ZERO
#elif defined(LAMBTEST) && defined(TWOCOMPONENT)


                  dist3d(1)=real(gi,kind=db)-center(1)
                  dist3d(2)=real(gj,kind=db)-center(2)
                  dist3d(3)=real(gk,kind=db)-center(3)
                  call pbc_images(invdim,dist3d,dist3dout)
                   
                  dist=sqrt((dist3dout(1)/(ONE+lamb_eps))**TWO + dist3dout(2)**TWO + dist3dout(3)**TWO)
                   
                  tempphi=ONE*fcut_tanh(dist,radius,width)
                   
				  loc_u = ZERO
				  loc_v = ZERO
				  loc_w = ZERO
			  
!                 if (dist <= radius) then
!                    loc_u =  TWO *lamb_A*dist3dout(1)
!                    loc_v =  TWO *lamb_A*dist3dout(2)
!                    loc_w = -FOUR*lamb_A*dist3dout(3)
!                 else
!                    loc_u = ZERO
!                    loc_v = ZERO
!                    loc_w = ZERO
!                 endif 
#ifdef DENSRATIO
                  rhophi_loc=rho_r*tempphi+(ONE-tempphi)*rho_b
#else
                  rhophi_loc = 1.0_db
#endif	
                  
                  tempphi2 = tempphi*(sigma*TWO)/radius/(rhophi_loc*cssq)
				  loc_press = loc_press + tempphi2
                  
                  visc_loc=(rho_r*visc1*tempphi+(1.0_db-tempphi)*visc2*rho_b)/rhophi_loc
                  
                  lamb_visc=lamb_visc+visc_loc*tempphi*(ONE-tempphi) 
                  lamb_visc_temp=lamb_visc_temp+tempphi*(ONE-tempphi)

#else
                   !dist=sqrt((float(gi)-lx/TWO)**TWO + (float(gj)-ly/TWO)**TWO+(float(gk)-(lz/TWO)+1.5*radius)**TWO)
                   dist3d(1)=real(gi,kind=db)-center(1)
                   dist3d(2)=real(gj,kind=db)-center(2)
                   dist3d(3)=real(gk,kind=db)-center(3)
                   call pbc_images(invdim,dist3d,dist3dout)
                   
                   dist=sqrt(dist3dout(1)**TWO + dist3dout(2)**TWO + dist3dout(3)**TWO)
                   
                   tempphi=ONE*fcut_tanh(dist,radius,width)
                   
#ifdef TWOCOMPONENT
                   phifields_flip(ii,jj,kk,1,myblock)=real(tempphi,kind=strdb)
#ifdef DENSRATIO
                   rhophi_loc=rho_r*tempphi+(ONE-tempphi)*rho_b
#else
                   rhophi_loc = 1.0_db
#endif				  

                  tempphi2 = tempphi*(sigma*TWO)/radius/(rhophi_loc*cssq)
				  loc_press = loc_press + tempphi2
				  loc_w=ZERO!fcut(dist,radius-width*0.5,radius+width*0.5)*uwall !   - fcut(dist2,radius-width*0.5,radius+width*0.5)*HALF*uwall
 !                if(gi==lx/2 .and. gj==ly/2 .and. gk==lz/2)phi(i,j,k)=ONE		 
#else				  
 
                  !if(gi==lx/2 .and. gj==ly/2 .and. gk==lz/2)then
                  !if(gi.le.lx/4+2 .and. gi.ge.lx/4-2 .and. gj.le.ly/4+2 .and. gj.ge.ly/4-2 .and. gk.le.lz/3+2 .and. gk.ge.lz/3-2)then
                   loc_press = loc_press + 0.01_db*tempphi
                  !endif
                  loc_w=ZERO!fcut(dist,radius-width*0.5,radius+width*0.5)*uwall !   - fcut(dist2,radius-width*0.5,radius+width*0.5)*HALF*uwall

#endif       
#endif

#ifdef VELUNIFORMV
                loc_press=ZERO
                loc_w=uwall        
#endif

                hfields_flip(ii,jj,kk,1,myblock)= real(loc_press,kind=strdb)
                hfields_flip(ii,jj,kk,2,myblock)=real(loc_u,kind=strdb)
                hfields_flip(ii,jj,kk,3,myblock)=real(loc_v,kind=strdb) 
                hfields_flip(ii,jj,kk,4,myblock)=real(loc_w,kind=strdb) 
                hfields_flip(ii,jj,kk,5,myblock)=real(loc_u*loc_u,kind=strdb)
                hfields_flip(ii,jj,kk,6,myblock)=real(loc_v*loc_v,kind=strdb)
                hfields_flip(ii,jj,kk,7,myblock)=real(loc_w*loc_w,kind=strdb)
                hfields_flip(ii,jj,kk,8,myblock)=real(loc_u*loc_v,kind=strdb)
                hfields_flip(ii,jj,kk,9,myblock)=real(loc_u*loc_w,kind=strdb)
                hfields_flip(ii,jj,kk,10,myblock)=real(loc_v*loc_w,kind=strdb)  

#ifdef TWOCOMPONENT
				 
				 phifields_flip(ii,jj,kk,1,myblock)=real(tempphi,kind=strdb)
				 
				 auxfields(ii,jj,kk,1,myblock)=ZEROSTR
				 auxfields(ii,jj,kk,2,myblock)=ZEROSTR
				 auxfields(ii,jj,kk,3,myblock)=ZEROSTR
				 auxfields(ii,jj,kk,4,myblock)=ZEROSTR
#endif
           
                endif
                !if(j==jprobe .and. k==kprobe)write(6,*)'velocità iniziale', i,loc_w
             enddo
          enddo
       enddo
#if defined(LAMBTEST) && defined(TWOCOMPONENT)
       call sum_world_float(lamb_visc)
       call sum_world_float(lamb_visc_temp)
       lamb_visc=lamb_visc/lamb_visc_temp
#endif


#endif

#if defined(INTERNAL_OBSTACLES) && defined(TWOCOMPONENT)     
		global_phi_sum_ini=ZERO
		do k=1,nz
		   gk = nz*coords(3) + k
		   zblock=(k+2*TILE_DIMz-1)/TILE_DIMz		
           do j=1,ny
			  gj = ny*coords(2) + j
			  yblock=(j+2*TILE_DIMy-1)/TILE_DIMy         
              do i=1,nx
				 gi = nx*coords(1) + i
				 xblock=(i+2*TILE_DIMx-1)/TILE_DIMx

                 myblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
                 ii=i-xblock*TILE_DIMx+2*TILE_DIMx
                 jj=j-yblock*TILE_DIMy+2*TILE_DIMy
                 kk=k-zblock*TILE_DIMz+2*TILE_DIMz
                 if(abs(isfluid(i,j,k)).eq.1)then
			       global_phi_sum_ini=global_phi_sum_ini + real(phifields_flip(ii,jj,kk,1,myblock),kind=db)
                 endif
              enddo
           enddo
      enddo
      call sum_world_float(global_phi_sum_ini)
      
      do k=1,nz
		   gk = nz*coords(3) + k
		   zblock=(k+2*TILE_DIMz-1)/TILE_DIMz		
           do j=1,ny
			  gj = ny*coords(2) + j
			  yblock=(j+2*TILE_DIMy-1)/TILE_DIMy         
              do i=1,nx
				 gi = nx*coords(1) + i
				 xblock=(i+2*TILE_DIMx-1)/TILE_DIMx

                 myblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
                 ii=i-xblock*TILE_DIMx+2*TILE_DIMx
                 jj=j-yblock*TILE_DIMy+2*TILE_DIMy
                 kk=k-zblock*TILE_DIMz+2*TILE_DIMz
                 if(isfluid(i,j,k).eq.0)then
                   tempphi=ZERO
                   tempphi2=ZERO
                   do link_id = 1, nlinks
		             link_dist2 = ex(link_id)*ex(link_id) + ey(link_id)*ey(link_id) + ez(link_id)*ez(link_id)
		             if (link_dist2 == 0) cycle
		             nei_x = i + ex(link_id)
		             nei_y = j + ey(link_id)
		             nei_z = k + ez(link_id)
                     if (nei_x < 1 .or. nei_x > nx) cycle
		             if (nei_y < 1 .or. nei_y > ny) cycle
		             if (nei_z < 1 .or. nei_z > nz) cycle
		             if (isfluid(nei_x,nei_y,nei_z) == 0) cycle
		             if (link_dist2 == 1) then
		               link_scale = 1.0_db
		             elseif (link_dist2 == 2) then
		               link_scale = inv_sqrt2
		             else
		               link_scale = inv_sqrt3
		             end if
		             oxblock=(nei_x+2*TILE_DIMx-1)/TILE_DIMx   
	                 oyblock=(nei_y+2*TILE_DIMy-1)/TILE_DIMy     
	                 ozblock=(nei_z+2*TILE_DIMz-1)/TILE_DIMz 
	                 omyblock=(oxblock-1)+(oyblock-1)*nxblock+(ozblock-1)*nxyblock+1
	                 oii=nei_x-oxblock*TILE_DIMx+2*TILE_DIMx
	                 ojj=nei_y-oyblock*TILE_DIMy+2*TILE_DIMy
	                 okk=nei_z-ozblock*TILE_DIMz+2*TILE_DIMz
	                 tempphi = tempphi + link_scale*real(phifields_flip(oii,ojj,okk,1,omyblock),kind=db)
	                 tempphi2 = tempphi2 + link_scale
			       enddo
			       tempphi=tempphi/tempphi2
			       phifields_flip(ii,jj,kk,1,myblock)=real(tempphi,kind=strdb)
                 endif
              enddo
           enddo
      enddo
#endif

! #if defined(MULTIHIT)
	  
	  ! AAA=2.0_db*1.0d-7
	  ! k_zero=2.0_db*2.0_db*pi_greek
      ! do k=1,nz
         ! gk=nz*coords(3)+k
         ! do j=1,ny
            ! gj=ny*coords(2)+j
            ! do i=1,nx
               ! gi=nx*coords(1)+i
               ! if(abs(isfluid(i,j,k)).eq.1)then
				  
				  ! ABCx(i,j,k)= AAA*sin(k_zero*real(gk-1)/real(lz)) + AAA*sin(k_zero*real(gj-1)/real(ly))  
				  ! ABCy(i,j,k)= AAA*sin(k_zero*real(gi-1)/real(lx)) + AAA*sin(k_zero*real(gk-1)/real(lz))
				  ! ABCz(i,j,k)= AAA*sin(k_zero*real(gj-1)/real(ly)) + AAA*sin(k_zero*real(gi-1)/real(lx)) 
				  ! u(i,j,k)= ABCx(i,j,k)
				  ! v(i,j,k)= ABCy(i,j,k)
				  ! w(i,j,k)= ABCz(i,j,k)
! #if defined(TWOCOMPONENT)
				  ! dist=sqrt((float(gi)-lx/TWO)**TWO + (float(gj)-ly/TWO)**TWO+(float(gk)-lz/TWO)**TWO)
				  ! !phi(i,j,k)=1.0_db-ONE*(fcut(dist,radius-width*0.5,radius+width*0.5))
				  ! selphi(i,j,k,flip) = 0.5 + 0.5*TANH(2.0_db*(radius-dist)/width)	!droplet
				  
				  ! !phi(i,j,k) = 0.5 + 0.5*TANH(2.0_db*(radius-dist)/width)	!droplet
! #ifdef DENSRATIO
                  ! rhophi(i,j,k)=rho_r*selphi(i,j,k,flip)+(1.0_db-selphi(i,j,k,flip))*rho_b
                  ! rhophi_loc = rhophi(i,j,k) 
				 

! #endif				  
				  ! !rho(i,j,k) = rho(i,j,k) + phi(i,j,k)*(sigma*2.0_db)/radius/(rhophi_loc/3.0_db)!droplet
				  ! rho(i,j,k) = rho(i,j,k) - selphi(i,j,k,flip)*(sigma*2.0_db)/radius/(rhophi_loc/3.0_db)!bubble
! #endif				  				  				  
               ! endif
            ! enddo
         ! enddo
      ! enddo
! #endif


     hfields_flop=hfields_flip
#ifdef TWOCOMPONENT
     phifields_flop=phifields_flip
#endif     
     
   endsubroutine
   
    subroutine pbc_images(aaa,xxs,xout)
      real(kind=db), dimension(3), intent(in) :: aaa,xxs
      real(kind=db), dimension(3), intent(out) :: xout
      
      xout(1) = xxs(1) - real(lx,kind=db)*nint(aaa(1)*xxs(1)) 
      xout(2) = xxs(2) - real(ly,kind=db)*nint(aaa(2)*xxs(2)) 
      xout(3) = xxs(3) - real(lz,kind=db)*nint(aaa(3)*xxs(3)) 
      
    end subroutine
   

endmodule
