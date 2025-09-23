#include "defines.h"
module stat_module

   use vars
#ifdef MPI
   use mpi_template, only : sum_world_farr
#endif
   implicit none

contains

#if defined(FALLING_DROP_noslip) || defined(MULTI_DROPS)
	subroutine rising_vel
		
		implicit none
		
		real(kind=db),dimension(3) :: vel_com
		integer :: ierr
		

		  vel_com=0.0_db
		  !$acc update host(phi,w,rhophi)
		  !$acc wait
		  do k=1,nz
			 do j=1,ny
				do i=1,nx
                     if(abs(isfluid(i,j,k))==1)vel_com(3)=selphi(i,j,k,flip)   !storo la massa
					 if(selphi(i,j,k,flip)<0.5_db)then 
						vel_com(1)=vel_com(1) + w(i,j,k)*(1.0_db-selphi(i,j,k,flip))
						vel_com(2)=vel_com(2) + 1.0_db*(1.0_db-selphi(i,j,k,flip)) 
					 endif
				enddo
			 enddo
		 enddo  
		 
		 call sum_world_farr(3,vel_com) 
           
		 write(193,*) (vel_com(1)/vel_com(2)),vel_com(3)
		 flush(193)
		 
		 
	endsubroutine
#endif

endmodule
