#include "defines.h"
module stat_module

   use vars
   use mpi_template, only : sum_world_farr,coords,myoffset,sum_world_float,dostop,myrank
   
   implicit none

contains

#if defined(FALLING_DROP_noslip) || defined(MULTI_DROPS)
	subroutine rising_vel
		
		implicit none
		
		real(kind=db),dimension(3) :: vel_com
		real(kind=db) :: phi_loc,w_loc
		integer :: ierr,i,j,k,gi,gj,gk,xblock,yblock,zblock,myblock
		integer :: ii,kk,jj
		

		  vel_com=0.0_db
		  !$acc update host(phifields_flip,hfields_flip)
		  !$acc wait
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

                phi_loc=real(phifields_flip(ii,jj,kk,1,myblock),kind=db)
                if(abs(isfluid(i,j,k))==1)vel_com(3)=vel_com(3)+phi_loc !storo la massa
			    if(phi_loc<0.5_db)then 
                  w_loc=real(hfields_flip(ii,jj,kk,4,myblock),kind=db)
				  vel_com(1)=vel_com(1) + w_loc*(1.0_db-phi_loc)
				  vel_com(2)=vel_com(2) + 1.0_db*(1.0_db-phi_loc) 
				endif
			  enddo
			enddo
		  enddo  
		 
		  call sum_world_farr(3,vel_com) 
           
		  write(193,*) (vel_com(1)/vel_com(2)),vel_com(3)
		  flush(193)
		 
		 
	endsubroutine
#endif

    subroutine compare_benchmark
    
      implicit none
      
      character(1), parameter :: bs = char(92)   ! backslash '\'
	  
#ifdef POISEUILLE
 
      real(kind=db) :: H_pois,xc_pois,dist1,rhophi_loc,mytemp
      real(kind=db), allocatable, dimension(:) :: w_pois,w_num,i_num
      integer :: xblock,yblock,zblock,myblock,ii,jj,kk


      allocate(w_pois(lx),w_num(lx),i_num(lx))
      w_num=ZERO;w_pois=ZERO;i_num=ZERO
      do gi = 1, lx
         rhophi_loc = 1.0_db
         H_pois  = 0.5_db * real(lx-2,db)
         xc_pois = 0.5_db * real(lx+1,db)
         dist1 = real(gi,db) - xc_pois
         dist1 = abs(dist1)
         if (dist1 <= H_pois) then
           w_pois(gi)=(fz)/(2.0_db*visc1) * (H_pois*H_pois - dist1*dist1)
         else
           w_pois(gi)=ZERO
         endif
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
                
               mytemp=real(hfields_flip(ii,jj,kk,4,myblock),kind=db)
               if(isfluid(i,j,k)==0)mytemp=ZERO
               gi = nx*coords(1) + i  
               w_num(gi)=w_num(gi) + mytemp
               i_num(gi)=i_num(gi)+ONE
            enddo
          enddo
       enddo
       call sum_world_farr(lx,w_num)
       call sum_world_farr(lx,i_num)
       
       do gi = 1, lx
         w_num(gi)=w_num(gi)/i_num(gi)
       enddo
       if(myrank==0)then
         open(unit=42,file='plot_poiseuille.dat',status='replace',action='write')
         do gi = 1, lx
           write(42,'(i8,2g16.8)')gi,w_num(gi),w_pois(gi)
         enddo
         close(42)
#ifdef USEGNUPLOT         
         open(unit=42,file='plot_poiseuille.gp',status='replace',action='write')
         write(42,'(a)')'# plot_poiseuille.gp'
         write(42,'(a)')'set terminal pngcairo size 1400,900 enhanced font "Helvetica,20"'
         write(42,'(a)')'set output "plot_poiseuille.png"'
         write(42,'(4a,es12.4,a,es12.4,a,es12.4,a)')'set title "Poiseuille: numerical vs analytical',bs,'n',&
         'Fz=',fz,'  nu=',visc1,'"'
         write(42,'(a)')'set xlabel "i"'
         write(42,'(a)')'set ylabel "u"'
         write(42,'(a)')'set grid'
         write(42,'(a)')'set key left top'
         write(42,'(a)')'set autoscale'
         write(42,'(a)') 'plot ' // char(92)
         write(42,'(a)') '"plot_poiseuille.dat" using 1:2 with lines lw 3 title "numerical", ' // char(92)
         write(42,'(a)')'"plot_poiseuille.dat" using 1:3 with lines lw 3 dt 2 title "analytical"'
         write(42,'(a)')'unset output'
         close(42)
         call execute_command_line('gnuplot plot_poiseuille.gp', wait=.true.)
#endif         
       endif
       deallocate(w_num,w_pois,i_num)
       
       
#endif    
    
    endsubroutine

endmodule
