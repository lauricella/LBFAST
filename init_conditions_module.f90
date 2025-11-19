#include "defines.h"
module initial_condts

   use vars
   use mpi_template, only: coords,myoffset,sum_world_float,dostop
   
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
    
    integer :: xblock,yblock,zblock,myblock,ii,jj,kk
    real(kind=db) :: tempphi,tempphi2

       
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

                 hfields_flip(idx5(ii,jj,kk,1,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields))=ZERO
                 hfields_flip(idx5(ii,jj,kk,2,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields))=ZERO
                 hfields_flip(idx5(ii,jj,kk,3,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields))=ZERO
                 hfields_flip(idx5(ii,jj,kk,4,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields))=ZERO
                               
                 
				 u(i,j,k)=ZERO
				 v(i,j,k)=ZERO
				 w(i,j,k)=ZERO
				 rho(i,j,k)=ZERO  !tot dens
#ifdef DENSRATIO
				 rhophi(i,j,k)=ZERO  !tot dens
				 !locauxfields(idx5(ii,jj,kk,5,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nauxfields))=ZERO
#endif
#ifdef TWOCOMPONENT
				 selphi(i,j,k,flip)=ZERO
				 
				 phifields_flip(idx5(ii,jj,kk,1,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields))=ZERO
				 
				 normx(i,j,k)=ZERO
				 normy(i,j,k)=ZERO
				 normz(i,j,k)=ZERO
				 auxfields(idx5(ii,jj,kk,1,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nauxfields))=ZERO
				 auxfields(idx5(ii,jj,kk,2,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nauxfields))=ZERO
				 auxfields(idx5(ii,jj,kk,3,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nauxfields))=ZERO
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
				 
				 selphi(i,j,k,flip) = tempphi
				 
				 phifields_flip(idx5(ii,jj,kk,1,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields))=tempphi

				 ! density interpolation and pressure correction
#ifdef DENSRATIO
				 rhophi(i,j,k) = rho_r * tempphi + &
								(ONE - tempphi) * rho_b
				 !locauxfields(idx5(ii,jj,kk,5,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nauxfields))=rho_r * tempphi + &
				!				(ONE - tempphi) * rho_b
#endif
			    ! crisp velocity: uniform inside each core (phi~1 region)
				!if (sel1 > 0.1_db) then 
				  u(i,j,k) = 0.02*sel1-0.02*sel2
                hfields_flip(idx5(ii,jj,kk,2,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields))=0.02*sel1-0.02*sel2
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
			     global_phi_sum_ini=global_phi_sum_ini + phifields_flip(idx5(ii,jj,kk,1,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields))
              enddo
           enddo
      enddo
      call sum_world_float(global_phi_sum_ini)
#endif

    

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

                 hfields_flip(idx5(ii,jj,kk,1,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields))=ZERO
                 hfields_flip(idx5(ii,jj,kk,2,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields))=ZERO
                 hfields_flip(idx5(ii,jj,kk,3,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields))=ZERO
                 hfields_flip(idx5(ii,jj,kk,4,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields))=ZERO
                               
                 
				 u(i,j,k)=ZERO
				 v(i,j,k)=ZERO
				 w(i,j,k)=ZERO
				 rho(i,j,k)=ZERO  !tot dens
#ifdef DENSRATIO
				 rhophi(i,j,k)=ZERO  !tot dens
				 !locauxfields(idx5(ii,jj,kk,5,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nauxfields))=ZERO
#endif
#ifdef TWOCOMPONENT
				 selphi(i,j,k,flip)=ZERO
				 
				 phifields_flip(idx5(ii,jj,kk,1,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields))=ZERO
				 
				 normx(i,j,k)=ZERO
				 normy(i,j,k)=ZERO
				 normz(i,j,k)=ZERO
				 auxfields(idx5(ii,jj,kk,1,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nauxfields))=ZERO
				 auxfields(idx5(ii,jj,kk,2,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nauxfields))=ZERO
				 auxfields(idx5(ii,jj,kk,3,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nauxfields))=ZERO
#endif
              enddo
            enddo
        enddo



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
                
                
                
                w(i,j,k)=0.0
                hfields_flip(idx5(ii,jj,kk,4,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields))=ZERO !w(i,j,k)
                
                if(abs(isfluid(i,j,k)).eq.1)then
#ifdef TWOCOMPONENT
                   !dist=sqrt((float(gi)-lx/TWO)**TWO + (float(gj)-ly/TWO)**TWO+(float(gk)-(lz/TWO)+1.5*radius)**TWO)
                   dist=sqrt((float(gi)-lx/TWO)**TWO + (float(gj)-ly/TWO)**TWO+(float(gk)-(lz/TWO))**TWO)
!                   selphi(i,j,k,flip)=ONE*(fcut(dist,radius-width*0.75_db,radius+width*0.75_db))! +fcut(dist2,radius-width*0.5,radius+width*0.5))
                   
                   tempphi=ONE*fcut_tanh(dist,radius,width)
                   
                   selphi(i,j,k,flip) = tempphi
                   phifields_flip(idx5(ii,jj,kk,1,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields))=tempphi
#ifdef DENSRATIO
                   rhophi(i,j,k)=rho_r*tempphi+(ONE-tempphi)*rho_b
                   rhophi_loc = rhophi(i,j,k)
                   
                   !locauxfields(idx5(ii,jj,kk,5,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nauxfields))=rho_r*tempphi+(ONE-tempphi)*rho_b !rhophi(i,j,k)
                   rhophi_loc = rho_r*tempphi+(ONE-tempphi)*rho_b !rhophi(i,j,k)
                   
#else
                   rhophi_loc = 1.0_db
#endif				  
                   tempphi2 = tempphi*(sigma*2.0_db)/radius/(rhophi_loc/3.0_db)
				   rho(i,j,k) = rho(i,j,k) + tempphi2
				   hfields_flip(idx5(ii,jj,kk,1,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields))= &
				    hfields_flip(idx5(ii,jj,kk,1,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields)) + tempphi2 !rho(i,j,k)
				   w(i,j,k)=ZERO!fcut(dist,radius-width*0.5,radius+width*0.5)*uwall !   - fcut(dist2,radius-width*0.5,radius+width*0.5)*HALF*uwall
                   hfields_flip(idx5(ii,jj,kk,4,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields))=ZERO !w(i,j,k)
 !                 if(gi==lx/2 .and. gj==ly/2 .and. gk==lz/2)phi(i,j,k)=ONE		 
#else				  
                  !rho(i,j,k) = 1.0_db
 
                  if(gi==lx/8 .and. gj==ly/8 .and. gk==lz/4)then
                    rho(i,j,k) = rho(i,j,k) + 0.05_db
                    hfields_flip(idx5(ii,jj,kk,1,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields)) = &
                     hfields_flip(idx5(ii,jj,kk,1,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields)) + 0.05_db
                  endif
                  w(i,j,k)=ZERO!fcut(dist,radius-width*0.5,radius+width*0.5)*uwall !   - fcut(dist2,radius-width*0.5,radius+width*0.5)*HALF*uwall
                  hfields_flip(idx5(ii,jj,kk,4,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields))=ZERO !w(i,j,k)
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
               
               tempphi=phifields_flip(idx5(ii,jj,kk,1,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields))
               rhophi_loc = rho_r*tempphi+(ONE-tempphi)*rho_b
               if(rhophi(i,j,k).ne.rhophi_loc)then
                 write(6,*)'cazzi1 in',i,j,k,nxblock,nxyblock
               endif
               
               if(abs(isfluid(i,j,k)).eq.1)then
                 pxx(i,j,k)=ZERO
                 pyy(i,j,k)=ZERO
                 pzz(i,j,k)=ZERO
                 pxy(i,j,k)=ZERO
                 pxz(i,j,k)=ZERO
                 pyz(i,j,k)=ZERO
                 
                 hfields_flip(idx5(ii,jj,kk,5,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields))=ZERO
                 hfields_flip(idx5(ii,jj,kk,6,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields))=ZERO
                 hfields_flip(idx5(ii,jj,kk,7,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields))=ZERO
                 hfields_flip(idx5(ii,jj,kk,8,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields))=ZERO
                 hfields_flip(idx5(ii,jj,kk,9,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields))=ZERO
                 hfields_flip(idx5(ii,jj,kk,10,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields))=ZERO
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

     hfields_flop=hfields_flip
     phifields_flop=phifields_flip
     

   endsubroutine

endmodule
