#include "defines.h"
module initial_condts

   use vars
   use mpi_template, only: coords,myoffset,sum_world_float

   implicit none

contains

   subroutine initial_conditions_all
#ifdef MPI  
	  use mpi
#endif
      implicit none
      integer:: subchords(3)
      real(kind=db) :: dist1,dist2,sel1,sel2,dist
      real(kind=db) :: fneq1,feq, rhophi_loc
	  
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
      step=0

      flip=mod(step,2)+1     
      flop = 3 - flip
       
      !*************************************initial conditions ************************
      u=ZERO
      v=ZERO
      w=ZERO
      rho=0.0_db  !tot dens
#ifdef DENSRATIO
      rhophi=ZERO  !tot dens
#endif
#ifdef TWOCOMPONENT
      selphi=0.0
      normx=0
      normy=0
      normz=0
#endif

		do k = 1, nz
		   gk = nz*coords(3) + k
		   do j = 1, ny
			  gj = ny*coords(2) + j
			  do i = 1, nx
				 gi = nx*coords(1) + i

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
				 selphi(i,j,k,flip) = max(sel1, sel2)
				 selphi(i,j,k,flop) = selphi(i,j,k,flip)

				 ! density interpolation and pressure correction
#ifdef DENSRATIO
				 rhophi(i,j,k) = rho_r * selphi(i,j,k,flip) + &
								(1.0_db - selphi(i,j,k,flip)) * rho_b
#endif
			    ! crisp velocity: uniform inside each core (phi~1 region)
				!if (sel1 > 0.1_db) then 
				  u(i,j,k) = 0.02*sel1-0.02*sel2
				!else if (sel2 > 0.1_db) then
				 ! u(i,j,k) = -0.01*sel2
				!else
				!  u(i,j,k) = 0.0_db
			    !endif
			  end do
		   end do
		end do
		
#ifdef INTERNAL_OBSTACLES     
		global_phi_sum_ini=ZERO
		do k=1,nz
         do j=1,ny
            do i=1,nx
			      global_phi_sum_ini=global_phi_sum_ini + selphi(i,j,k,flip)
            enddo
         enddo
      enddo
      call sum_world_float(global_phi_sum_ini)
#endif

    

#if defined(DOBENCHMARK)

      u=ZERO
      v=ZERO
      w=ZERO
      rho=0.0_db  !tot dens
#ifdef DENSRATIO
      rhophi=ZERO  !tot dens
#endif
#ifdef TWOCOMPONENT
      selphi=0.0
      normx=0
      normy=0
      normz=0
#endif


       do k=1,nz
          gk=nz*coords(3)+k
          do j=1,ny
             gj=ny*coords(2)+j
             do i=1,nx
                gi=nx*coords(1)+i
                w(i,j,k)=0.0
                if(abs(isfluid(i,j,k)).eq.1)then
#ifdef TWOCOMPONENT
                   !dist=sqrt((float(gi)-lx/TWO)**TWO + (float(gj)-ly/TWO)**TWO+(float(gk)-(lz/TWO)+1.5*radius)**TWO)
                   dist=sqrt((float(gi)-lx/TWO)**TWO + (float(gj)-ly/TWO)**TWO+(float(gk)-(lz/TWO))**TWO)
!                   selphi(i,j,k,flip)=ONE*(fcut(dist,radius-width*0.75_db,radius+width*0.75_db))! +fcut(dist2,radius-width*0.5,radius+width*0.5))
                   selphi(i,j,k,flip) = ONE*fcut_tanh(dist,radius,width)
#ifdef DENSRATIO
                   rhophi(i,j,k)=rho_r*selphi(i,j,k,flip)+(1.0_db-selphi(i,j,k,flip))*rho_b
                   rhophi_loc = rhophi(i,j,k)
#else
                   rhophi_loc = 1.0_db
#endif				  
				   rho(i,j,k) = rho(i,j,k) + selphi(i,j,k,flip)*(sigma*2.0_db)/radius/(rhophi_loc/3.0_db)
				   w(i,j,k)=ZERO!fcut(dist,radius-width*0.5,radius+width*0.5)*uwall !   - fcut(dist2,radius-width*0.5,radius+width*0.5)*HALF*uwall
 !                 if(gi==lx/2 .and. gj==ly/2 .and. gk==lz/2)phi(i,j,k)=ONE		 
#else				  
                  !rho(i,j,k) = 1.0_db
 
                  if(gi==lx/8 .and. gj==ly/8 .and. gk==lz/4)rho(i,j,k) = rho(i,j,k) + 0.05_db
                  w(i,j,k)=ZERO!fcut(dist,radius-width*0.5,radius+width*0.5)*uwall !   - fcut(dist2,radius-width*0.5,radius+width*0.5)*HALF*uwall
#endif                                
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


!!!!initialize non-eq flux tensor
      do k=1,nz
         do j=1,ny
            do i=1,nx
               if(abs(isfluid(i,j,k)).eq.1)then
                 pxx(i,j,k)=ZERO
                 pyy(i,j,k)=ZERO
                 pzz(i,j,k)=ZERO
                 pxy(i,j,k)=ZERO
                 pxz(i,j,k)=ZERO
                 pyz(i,j,k)=ZERO
               endif
            enddo
         enddo
      enddo

!!!!common to every LB
      do k=1,nz
         gk=nz*coords(3)+k
         do j=1,ny
            gj=ny*coords(2)+j
            do i=1,nx
               gi=nx*coords(1)+i
               if(abs(isfluid(i,j,k)).eq.1)then
				  uu=HALF*(u(i,j,k)*u(i,j,k)+v(i,j,k)*v(i,j,k)+w(i,j,k)*w(i,j,k))/cssq
				  do l=0,nlinks
					udotc=(u(i,j,k)*dex(l) + v(i,j,k)*dey(l)+ w(i,j,k)*dez(l))/cssq
					f(i,j,k,l)=p(l)*(rho(i,j,k) + udotc+HALF*udotc*udotc - uu)
				  enddo
               endif
            enddo
         enddo
      enddo
      
      



         selphi(:,:,:,flop)=selphi(:,:,:,flip)



   endsubroutine

endmodule
