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
      real(kind=db),dimension(3) :: dist3d,dist3dout,invdim
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
    real(kind=db) :: tempphi,tempphi2,loc_u,loc_v,loc_w,loc_press
    
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

			     global_phi_sum_ini=global_phi_sum_ini + real(phifields_flip(ii,jj,kk,1,myblock),kind=db)

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

                   tempphi2 = tempphi*(sigma*2.0_db)/radius/(rhophi_loc/3.0_db)
				   loc_press = loc_press + tempphi2
				   loc_w=ZERO!fcut(dist,radius-width*0.5,radius+width*0.5)*uwall !   - fcut(dist2,radius-width*0.5,radius+width*0.5)*HALF*uwall
 !                 if(gi==lx/2 .and. gj==ly/2 .and. gk==lz/2)phi(i,j,k)=ONE		 
#else				  
 
                  !if(gi==lx/2 .and. gj==ly/2 .and. gk==lz/2)then
                  !if(gi.le.lx/4+2 .and. gi.ge.lx/4-2 .and. gj.le.ly/4+2 .and. gj.ge.ly/4-2 .and. gk.le.lz/3+2 .and. gk.ge.lz/3-2)then
                   loc_press = loc_press + 0.01_db*tempphi
                  !endif
                  loc_w=ZERO!fcut(dist,radius-width*0.5,radius+width*0.5)*uwall !   - fcut(dist2,radius-width*0.5,radius+width*0.5)*HALF*uwall

#endif       
#endif
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
