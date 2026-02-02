#include "defines.h"
module stat_module

   use vars
   use mpi_template, only : sum_world_farr,coords,myoffset,sum_world_float,dostop
   
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
#ifdef MIXEDPRC
                phi_loc=real(phifields_flip(ii,jj,kk,1,myblock),kind=db)
#else
                phi_loc=phifields_flip(ii,jj,kk,1,myblock)
#endif
                if(abs(isfluid(i,j,k))==1)vel_com(3)=vel_com(3)+phi_loc !storo la massa
			    if(phi_loc<0.5_db)then 
#ifdef MIXEDPRC
                  w_loc=real(hfields_flip(ii,jj,kk,4,myblock),kind=db)
#else
			      w_loc=hfields_flip(ii,jj,kk,4,myblock)
#endif
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

endmodule
