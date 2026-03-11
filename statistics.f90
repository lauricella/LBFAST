#include "defines.h"
module stat_module

   use vars
   use mpi_template, only : sum_world_farr,coords,myoffset,sum_world_float,dostop,myrank
   
   implicit none
   
   real(kind=db) :: Ekin,Ekin0

contains

   subroutine probe_loc(hfields_s &
#ifdef TWOCOMPONENT	       
     ,phifields_s &
#endif     
     ,mystring,mystring2,inumber)
      

      implicit none
      
      real(kind=strdb), allocatable, dimension(:,:,:,:,:) :: hfields_s
#ifdef TWOCOMPONENT	       
      real(kind=strdb), allocatable, dimension(:,:,:,:,:) :: phifields_s
#endif    
      character(len=*), optional :: mystring
      character(len=*), optional :: mystring2
      integer, optional :: inumber
      integer :: subchords(3),xblock,yblock,zblock,myblock,ii,jj,kk,ierr
      
      !$acc update host(hfields_s)
#ifdef TWOCOMPONENT	        
      !$acc update host(phifields_s)
#endif    
      !$acc wait

      gi=iprobe;gj=jprobe;gk=kprobe
      subchords(1)=(gi-1)/nx
      subchords(2)=(gj-1)/ny
      subchords(3)=(gk-1)/nz 
      if(all(subchords==coords))then      
        i=gi-myoffset(1)
        j=gj-myoffset(2)
        k=gk-myoffset(3)
        zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
	    yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
	    xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
        myblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
        ii=i-xblock*TILE_DIMx+2*TILE_DIMx
        jj=j-yblock*TILE_DIMy+2*TILE_DIMy
        kk=k-zblock*TILE_DIMz+2*TILE_DIMz 
        write(6,*)'at step :',step
        if(present(mystring))then
          if(present(inumber))then
            write(6,'(4a,i0)')trim(mystring),' - ',trim(mystring2),':',inumber
            call flush(6)
          else
            if(present(mystring2))then
              write(6,'(3a)')trim(mystring),' - ',trim(mystring2)
              call flush(6)
            else
              write(6,'(a)')mystring
              call flush(6)
            endif
          endif
        endif
        write(6,*)hfields_s(ii,jj,kk,1,myblock)     
		write(6,*)hfields_s(ii,jj,kk,2,myblock)    
		write(6,*)hfields_s(ii,jj,kk,3,myblock)  
		write(6,*)hfields_s(ii,jj,kk,4,myblock)  
		write(6,*)hfields_s(ii,jj,kk,5,myblock)  
		write(6,*)hfields_s(ii,jj,kk,6,myblock)  
		write(6,*)hfields_s(ii,jj,kk,7,myblock)  
		write(6,*)hfields_s(ii,jj,kk,8,myblock)
		write(6,*)hfields_s(ii,jj,kk,9,myblock)  
		write(6,*)hfields_s(ii,jj,kk,10,myblock)
#ifdef TWOCOMPONENT	        
        write(6,*)phifields_s(ii,jj,kk,10,myblock)
#endif              
      endif
      
      return
      
    end subroutine probe_loc

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

    subroutine linear_regression(x, y, n, a, b)
	  implicit none
	  integer, intent(in) :: n
	  real(kind=db), intent(in)  :: x(n), y(n)
	  real(kind=db), intent(out) :: a, b
	
	  integer :: i
	  real(kind=db) :: sx, sy, sxx, sxy, denom
	
	  sx  = 0.0_db
	  sy  = 0.0_db
	  sxx = 0.0_db
	  sxy = 0.0_db
	
	  do i = 1, n
	     sx  = sx  + x(i)
	     sy  = sy  + y(i)
	     sxx = sxx + x(i)*x(i)
	     sxy = sxy + x(i)*y(i)
	  end do
	
	  denom = real(n,kind=db)*sxx - sx*sx
	
	  if (abs(denom) < 1.0e-30_db) then
	     a = 0.0_db
	     b = 0.0_db
	  else
	     a = (real(n,kind=db)*sxy - sx*sy) / denom
	     b = (sy - a*sx) / real(n,kind=db)
	  end if
	end subroutine linear_regression
	
	
	subroutine fit_taylorgreen_nu(filename, Lbox, step_max_fit, nu_num, slope, intercept)
	  implicit none
	  character(len=*), intent(in) :: filename
	  real(kind=db),    intent(in) :: Lbox
	  integer,          intent(in) :: step_max_fit
	  real(kind=db),    intent(out):: nu_num, slope, intercept
	
	  integer, parameter :: nmax = 100000
	  integer :: ios, n, step_i
	  real(kind=db) :: dum, loge
	  real(kind=db) :: x(nmax), y(nmax)
	  real(kind=db) :: k2
	
	  n = 0
	
	  open(unit=77,file=filename,status='old',action='read',iostat=ios)
	  if (ios /= 0) then
	     nu_num    = -1.0_db
	     slope     = 0.0_db
	     intercept = 0.0_db
	     return
	  end if
	
	  do
	     read(77,*,iostat=ios) step_i, dum, loge
	     if (ios /= 0) exit
	
	     if (step_i <= step_max_fit) then
	        n = n + 1
	        if (n > nmax) exit
	        x(n) = real(step_i,kind=db)
	        y(n) = loge
	     end if
	  end do
	  close(77)
	
	  if (n < 2) then
	     nu_num    = -1.0_db
	     slope     = 0.0_db
	     intercept = 0.0_db
	     return
	  end if
	
	  call linear_regression(x(1:n), y(1:n), n, slope, intercept)
	
	  k2 = (TWO*pi_greek/Lbox)**2
	  nu_num = -slope / (6.0_db*k2)
	
	end subroutine fit_taylorgreen_nu

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
#ifdef BOUNCE_BACK
         gi=1
         write(42,'(3g16.8)')1.5_db,w_num(gi),w_pois(gi)
         do gi = 2, lx-1
           write(42,'(3g16.8)')real(gi,db),w_num(gi),w_pois(gi)
         enddo
         gi=lx
         write(42,'(3g16.8)')real(lx,db)-0.5_db,w_num(gi),w_pois(gi)
#else         
         do gi = 1, lx
           write(42,'(3g16.8)')real(gi,db),w_num(gi),w_pois(gi)
         enddo
#endif
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

#elif defined(TAYLORGREEN) && !defined(TWOCOMPONENT) 

       real(kind=db) :: mytime, eth, log_eth ,k2,nu_num,fit_slope,fit_intercept
       integer :: stepfit
 
       if(myrank==0)then   
         open(unit=42,file='taylorgreen_theory.dat',status='replace',action='write')
         do i = 0, nsteps/stamp_term
           mytime   = real(i*stamp_term,kind=db)
           k2 = (TWO*pi_greek/real(lx,db))**2
           eth     = exp(-6.0_db*visc1*k2*mytime)
           log_eth = -6.0_db*visc1*k2*mytime
           write(42,'(i12,2es20.10)') i*stamp_term, eth, log_eth
           if(log_eth>=-0.25)stepfit=i*stamp_term
         enddo
         close(42)
         call fit_taylorgreen_nu('taylorgreen.dat', real(lx,kind=db), stepfit, nu_num, fit_slope, fit_intercept)
#ifdef USEGNUPLOT
         open(unit=42,file='plot_taylorgreen.gp',status='replace',action='write')
         write(42,'(a)') '# plot_taylorgreen.gp'
         write(42,'(a)') 'set terminal pngcairo size 1400,900 enhanced font "Helvetica,20"'
         write(42,'(a)') 'set output "plot_taylorgreen.png"'
         write(42,'(a,es12.4,a,es12.4,a,es12.4,a)') &
              'set title "Taylor-Green: numerical vs theoretical log(E/E_{0})'//bs//'n'// &
              'u0=', uwall, '   nu_{th}=', visc1, '   nu_{fit}=', nu_num, '"'
         write(42,'(a)') 'set title offset 0,-0.5'
         write(42,'(a)') 'set xlabel "time step"'
         write(42,'(a)') 'set ylabel "log(E/E_{0})"'
         write(42,'(a)') 'set grid'
         write(42,'(a)') 'set key right top'
         write(42,'(a)') 'set autoscale'
         write(42,'(a)') 'plot ' // char(92)
         write(42,'(a)') '"taylorgreen.dat" using 1:3 with lines lw 3 title "numerical", ' // char(92)
         write(42,'(a)') '"taylorgreen_theory.dat" using 1:3 with lines lw 3 dt 2 title "theoretical"'
         write(42,'(a)') 'unset output'
         close(42)
         call execute_command_line('gnuplot plot_taylorgreen.gp', wait=.true.)
#endif
       endif
#endif    
    
    endsubroutine

endmodule
