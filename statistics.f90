#include "defines.h"
module stat_module

   use vars
   use mpi_template, only : sum_world_farr,coords,myoffset,sum_world_float, &
    dostop,myrank,doerror
   
   implicit none
   
   real(kind=db) :: Ekin,Ekin0
#if defined(LAMBTEST) && defined(TWOCOMPONENT)   
   real(kind=db) :: pos_x_int_left=ZERO
   real(kind=db) :: pos_x_int_node_left=ZERO
   real(kind=db) :: pos_x_int_right=ZERO
   real(kind=db) :: pos_x_int_node_right=ZERO
   real(kind=db) :: pos_z_int_left=ZERO
   real(kind=db) :: pos_z_int_node_left=ZERO
   real(kind=db) :: pos_z_int_right=ZERO
   real(kind=db) :: pos_z_int_node_right=ZERO
   real(kind=db) :: lamb_z,lamb_cm_z,lamb_x,lamb_cm_x
   real(kind=db) :: lamb_dosc,lamb_A,lamb_visc
#endif
#ifdef LAPLACE
   real(kind=db) :: pos_x_int_left=ZERO
   real(kind=db) :: pos_x_int_node_left=ZERO
   real(kind=db) :: pos_x_int_right=ZERO
   real(kind=db) :: pos_x_int_node_right=ZERO   
   real(kind=db) :: lamb_x,lamb_cm_x,laplace_rad
#endif
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
	  
#if defined(TWOCOMPONENT) && defined(TWOPOISEUILLE)

      real(kind=db) :: H_pois,xc_pois,dist1,rhophi_loc,mytemp,distabs
      real(kind=db) :: coefB_pois,coefL_pois,coefR_pois,phi_loc,wleft,wright
      real(kind=db), allocatable, dimension(:) :: w_pois,w_num,i_num
      integer :: xblock,yblock,zblock,myblock,ii,jj,kk

      allocate(w_pois(lx),w_num(lx),i_num(lx))
      w_num=ZERO;w_pois=ZERO;i_num=ZERO
      
#ifdef BOUNCE_BACK
      H_pois  = 0.5_db * real(lx-2,db)
#else
      H_pois  = 0.5_db * real(lx-1,db)
#endif
      xc_pois = 0.5_db * real(lx+1,db)      
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
      do gi = 1, lx
	    dist1   = real(gi,db) - xc_pois
	    distabs = abs(dist1)
	    phi_loc=fcut_tanh(dist1,ZERO,width)
!               if (dist1 <= ZERO) then
!                  selphi(i,j,k,flip) = ONE
!               else
!                  selphi(i,j,k,flip) = ZERO
!               endif

#ifdef DENSRATIO
	    rhophi_loc = rho_r*phi_loc + (ONE-phi_loc)*rho_b        
#endif

	    if (distabs <= H_pois) then
		  !if (dist1 <= ZERO) then
			 wleft = -(fz/(2.0_db*visc1)) * dist1*dist1 + &
						 coefL_pois * dist1 + coefB_pois
		  !else
			 wright = -(fz/(2.0_db*visc2)) * dist1*dist1 + &
						 coefR_pois * dist1 + coefB_pois
		  !endif
		  w_pois(gi)=wleft*phi_loc + (ONE-phi_loc)*wright   
	    else
		  w_pois(gi) = ZERO
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
         open(unit=42,file='plot_twopoiseuille.dat',status='replace',action='write')
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
         open(unit=42,file='plot_twopoiseuille.gp',status='replace',action='write')
         write(42,'(a)')'# plot_twopoiseuille.gp'
         write(42,'(a)')'set terminal pngcairo size 1400,900 enhanced font "Helvetica,20"'
         write(42,'(a)')'set output "plot_twopoiseuille.png"'
         write(42,'(4a,es12.4,a,es12.4,a,es12.4,3a,es12.4,a,es12.4,a)') &
          'set title "Poiseuille: numerical vs analytical',bs,'n', &
          'Fz=',fz,'  nu1=',visc1,'  nu2=',visc2,bs,'n', &
          'rho1=',rho_r,'  rho2=',rho_b,'"'
         write(42,'(a)')'set xlabel "i"'
         write(42,'(a)')'set ylabel "u"'
         write(42,'(a)')'set grid'
         write(42,'(a)')'set key left top'
         write(42,'(a)')'set autoscale'
         write(42,'(a)') 'plot ' // char(92)
         write(42,'(a)') '"plot_twopoiseuille.dat" using 1:2 with lines lw 3 title "numerical", ' // char(92)
         write(42,'(a)')'"plot_twopoiseuille.dat" using 1:3 with lines lw 3 dt 2 title "analytical"'
         write(42,'(a)')'unset output'
         close(42)
         call execute_command_line('gnuplot plot_twopoiseuille.gp', wait=.true.)
#endif         
       endif
       deallocate(w_num,w_pois,i_num) 
	  
#elif defined(POISEUILLE)
 
      real(kind=db) :: H_pois,xc_pois,dist1,rhophi_loc,mytemp
      real(kind=db), allocatable, dimension(:) :: w_pois,w_num,i_num
      integer :: xblock,yblock,zblock,myblock,ii,jj,kk


      allocate(w_pois(lx),w_num(lx),i_num(lx))
      w_num=ZERO;w_pois=ZERO;i_num=ZERO
#ifdef BOUNCE_BACK
      H_pois  = 0.5_db * real(lx-2,db)
#else
      H_pois  = 0.5_db * real(lx-1,db)
#endif
      xc_pois = 0.5_db * real(lx+1,db)         
      do gi = 1, lx
         rhophi_loc = 1.0_db
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


#elif defined(LAMBTEST) && defined(TWOCOMPONENT)

       real(kind=db) :: mytime, dosc_th,myfreq,fit_slope,fit_intercept,amp0, &
        period_sum,tmean,omega_num,nrat,mu1,mu2,chi,myfreq_corr,zint0_minus_zcm
       integer :: stepfit,maxn,npeaks,myn,ios,nfit
       real(kind=db) :: dt, xfit, yfit, beta_num
       real(kind=db) :: sx, sy, sxx, sxy
       real(kind=db) :: lamb_x_tmp, lamb_z_tmp, lamb_cm_x_tmp, lamb_cm_z_tmp,my_tmp
       real(kind=db) :: Aenv, t1
       integer, allocatable :: mystep(:), pstep(:)
       real(db), allocatable :: dosc(:), adosc(:), pval(:)
       integer, parameter :: min_peak_dist = 2000
 
		if(myrank==0)then
		
		  dt = ONE
		
		  myfreq = sqrt(TWENTYFOUR*sigma / (radius**THREE * (TWO*rho_b + THREE*rho_r)))
		  nrat = rho_r/rho_b
          mu1  = rho_r*visc1
          mu2  = rho_b*visc2

          chi = ((TWO*nrat + ONE)**TWO * sqrt(mu1*mu2*rho_r*rho_b)) / &
           (TWO*radius * (nrat*rho_b + (nrat + ONE)*rho_r) * &
           (sqrt(mu1*rho_r) + sqrt(mu2*rho_b)))
		  myfreq_corr = myfreq - HALF*chi*sqrt(myfreq) + (ONE/FOUR)*chi**TWO
		  
		  maxn = nsteps/stamp_term + 1
		
		  allocate(mystep(0:maxn), dosc(0:maxn), adosc(0:maxn))
		  allocate(pstep(maxn), pval(maxn))
		
		  mystep(:) = 0
		  pstep(:)  = 0
		  dosc(:)   = ZERO
		  adosc(:)  = ZERO
		  pval(:)   = ZERO
		
		  open(unit=142,file='lamb.dat',status='old',action='read',iostat=ios)
		  if(ios /= 0) then
		    write(6,*) 'Errore apertura lamb.dat'
		    return
		  endif
		
		  myn = 0
		  do
		    if(myn >= maxn) exit
		    myn = myn + 1
#ifdef MILLER 
		    read(142,'(i8,7g16.8)',iostat=ios) mystep(myn),my_tmp, &
		         lamb_x_tmp, lamb_z_tmp, lamb_cm_x_tmp, lamb_cm_z_tmp,dosc(myn)
#else
		    read(142,'(i8,6g16.8)',iostat=ios) mystep(myn), dosc(myn), &
		         lamb_x_tmp, lamb_z_tmp, lamb_cm_x_tmp, lamb_cm_z_tmp,my_tmp
#endif
		    if(ios /= 0) then
		      myn = myn - 1
		      exit
		    endif
		  enddo
		  close(142)
		
		  if(myn < 3) then
		    write(6,*) 'Troppi pochi punti in lamb.dat'
		    return
		  endif
		  
		  mystep(0)=0
		  dosc(0)=dosc(1)-ONE
		  zint0_minus_zcm=dosc(1)
          amp0 = zint0_minus_zcm - radius


		  
		  
		
		  do i = 0, myn
		    adosc(i) = abs(dosc(i))
		  enddo
		
#ifdef MILLER 
                  npeaks = 0
                  do i = 1, myn-1
                    if ( dosc(i) < dosc(i-1) .and. dosc(i) <= dosc(i+1) .and. dosc(i) > 1.1*radius ) then
                      if (npeaks == 0 .or. mystep(i)-pstep(npeaks) >= min_peak_dist) then
                        npeaks = npeaks + 1
                        pstep(npeaks) = mystep(i)
                        pval(npeaks)  = dosc(i)
                      endif
                    endif
                  enddo
                  write(6,*)'internal sanity check',pval(1),amp0+radius
#else                  
                  npeaks = 0
                  do i = 2, myn-1
                    if ( dosc(i) < dosc(i-1) .and. dosc(i) <= dosc(i+1) .and. dosc(i) < -1.0e-3_db ) then
                      if (npeaks == 0 .or. mystep(i)-pstep(npeaks) >= min_peak_dist) then
                        npeaks = npeaks + 1
                        pstep(npeaks) = mystep(i)
                        pval(npeaks)  = -dosc(i)
                      endif
                    endif
                  enddo
#endif                		
		  if (npeaks < 2) then
		    write(6,*) 'Non ho trovato abbastanza picchi per stimare omega.'
		    return
		  endif
		 
     
		  do i = 1, npeaks
            if (pval(i) <= 1.0e-8_db) then
              write(6,'(a,i8,a,g16.8,a,i12)') 'Picco troppo piccolo: indice picco = ', i, &
              '  ampiezza = ', pval(i), '  step = ', pstep(i)
               return
            endif
         enddo
		
		  ! Picchi di abs(dosc): distanza tra picchi consecutivi = mezzo periodo
		  period_sum = ZERO
                  do i = 1, npeaks-1
                    period_sum = period_sum + real(pstep(i+1) - pstep(i), db) * dt
                  enddo

		
		  tmean = period_sum / real(npeaks-1, db)
		  omega_num = TWO*pi_greek / tmean
#ifndef MILLER			
		  ! Fit lineare di log(pval) = intercept - beta * t
		  sx  = ZERO
		  sy  = ZERO
		  sxx = ZERO
		  sxy = ZERO
		  nfit = 0
		
		  do i = 1, npeaks
		    if (pval(i) > ZERO) then
		      xfit = real(pstep(i), db) * dt
		      yfit = log(pval(i))
		      sx   = sx  + xfit
		      sy   = sy  + yfit
		      sxx  = sxx + xfit*xfit
		      sxy  = sxy + xfit*yfit
		      nfit = nfit + 1
		    endif
		  enddo
		
		  if (nfit >= 2) then
		    fit_slope = ( real(nfit,db)*sxy - sx*sy ) / ( real(nfit,db)*sxx - sx*sx )
		    fit_intercept = ( sy - fit_slope*sx ) / real(nfit,db)
		    beta_num = -fit_slope
		  else
		    fit_slope = ZERO
		    fit_intercept = ZERO
		    beta_num = ZERO
		    write(6,*) 'Non abbastanza picchi validi per stimare beta.'
		  endif
#endif		
		  write(6,'(a,i10)')    'Numero punti attesi             = ', maxn
		  write(6,'(a,i10)')    'Numero punti letti              = ', myn
		  write(6,'(a,i10)')    'Numero picchi trovati           = ', npeaks
		  write(6,'(a,f20.10)') 'Omega teorica                   = ', myfreq_corr
		  write(6,'(a,f20.10)') 'Periodo medio T_num             = ', tmean
		  write(6,'(a,f20.10)') 'Frequenza angolare omega_num    = ', omega_num
#ifndef MILLER		  
		  write(6,'(a,f20.10)') 'Dissipazione beta_num           = ', beta_num
		  write(6,'(a,f20.10)') 'Fit slope                       = ', fit_slope
		  write(6,'(a,f20.10)') 'Fit intercept                   = ', fit_intercept
#endif		  
                  do i = 1, npeaks
                    write(6,'(a,i6,a,i10,a,g16.8)') 'peak ',i,' step=',pstep(i),' amp=',pval(i)
                  enddo
         t1   = real(pstep(1),db)*dt
#ifndef MILLER         
         Aenv = pval(1)*exp(beta_num*t1)
#endif
         open(unit=42,file='lamb_theory.dat',status='replace',action='write')
         
         do i = 1, myn
           mytime  = real(mystep(i),kind=db)*dt
#ifdef MILLER 
           dosc_th = radius + amp0 * cos(myfreq_corr * mytime)
#else
           dosc_th = -Aenv*sin(myfreq_corr*mytime)
#endif
           write(42,'(i8,g16.8)') mystep(i), dosc_th
         enddo
         close(42)
  
#ifdef USEGNUPLOT
         open(unit=42,file='plot_lamb.gp',status='replace',action='write')
         write(42,'(a)') '# plot_lamb.gp'
         write(42,'(a)') 'set terminal pngcairo size 1400,900 enhanced font "Helvetica,20"'
         write(42,'(a)') 'set output "plot_lamb.png"'
         write(42,'(a,es12.4,a,es12.4,a,es12.4,a)') &
              'set title "Lamb: numerical vs theoretical D'//bs//'n'// &
              'u0=', uwall, '   \omega_{th}=', myfreq_corr, '   \omega_{fit}=',omega_num, '"'
         write(42,'(a)') 'set title offset 0,-0.5'
         write(42,'(a)') 'set xlabel "time step"'
         write(42,'(a)') 'set ylabel "D"'
         write(42,'(a)') 'set grid'
         write(42,'(a)') 'set key right top'
         write(42,'(a)') 'set autoscale'
         write(42,'(a)') 'plot ' // char(92)
#ifdef MILLER
         write(42,'(a)') '"lamb.dat" using 1:7 with lines lw 3 title "numerical", ' // char(92)
#else         
         write(42,'(a)') '"lamb.dat" using 1:2 with lines lw 3 title "numerical", ' // char(92)
#endif
         write(42,'(a)') '"lamb_theory.dat" using 1:2 with lines lw 3 dt 2 title "theoretical"'
         write(42,'(a)') 'unset output'
         close(42)
         call execute_command_line('gnuplot plot_lamb.gp', wait=.true.)
#endif
         deallocate(mystep, dosc, adosc, pstep, pval)
       endif
#elif defined(LAPLACE) && defined(TWOCOMPONENT)

       real(kind=db) :: mytime
       integer :: stepfit,maxn,npeaks,myn,ios,nfit,mystep
       real(kind=db) :: dt,mysum,mysqsum,avg,std,avgP,stdP,avgR,stdR
       real(kind=db) :: mysumP,mysqsumP,mysumR,mysqsumR
       real(kind=db) :: sigma_eff,delta_p,laplace_rad,p_in,p_out
       
 
		if(myrank==0)then
		
		  dt = ONE

		
		  open(unit=142,file='laplace.dat',status='old',action='read',iostat=ios)
		  if(ios /= 0) then
		    write(6,*) 'Errore apertura laplace.dat'
		    return
		  endif
		
		  myn = 0
		  mysum=ZERO
		  mysqsum=ZERO
		  mysumP=ZERO;mysqsumP=ZERO;mysumR=ZERO;mysqsumR=ZERO
		  do
		    read(142,'(i8,5g16.8)',iostat=ios) mystep, sigma_eff,delta_p,laplace_rad,p_in,p_out
		    if(ios /= 0) then
		      exit
		    endif
		    myn = myn + 1
		    mysum=mysum+sigma_eff
		    mysqsum=mysqsum+sigma_eff*sigma_eff
		    mysumP=mysumP+delta_p
		    mysqsumP=mysqsumP+delta_p*delta_p
		    mysumR=mysumR+laplace_rad
		    mysqsumR=mysqsumR+laplace_rad*laplace_rad
		  enddo
		  close(142)
		
		  if(myn < 3) then
		    write(6,*) 'Troppi pochi punti in laplace.dat'
		    return
		  endif
		  
		  avg = mysum / real(myn, db)
          std = sqrt( real(myn,db)*mysqsum - mysum*mysum ) / real(myn,db)
          avgP = mysumP / real(myn, db)
          stdP = sqrt( real(myn,db)*mysqsumP - mysumP*mysumP ) / real(myn,db)
          avgR = mysumR / real(myn, db)
          stdR = sqrt( real(myn,db)*mysqsumR - mysumR*mysumR ) / real(myn,db)
		  
		  write(6,'(a,3i10)')   'posizione interna               = ', nint(center(1:3))
		  write(6,'(a,3i10)')   'posizione esterna               = ', iprobe,jprobe,kprobe
		  write(6,'(a,i10)')    'Numero punti letti              = ', myn
		  write(6,'(a,f20.10)') 'Sigma teorica                   = ', sigma
		  write(6,'(a,f20.10,a,f20.10)') 'Sigma numerica                  = ', avg,' +- ',std
		  write(6,'(a,f20.10,a,f20.10)') 'Press numerica                  = ', avgP,' +- ',stdP
		  write(6,'(a,f20.10,a,f20.10)') 'Radius numerica                 = ', avgR,' +- ',stdR

		  open(unit=142, file='laplace.dat', status='old', action='write', position='append', iostat=ios)
		  if(ios /= 0) then
		    write(6,*) 'Errore apertura laplace.dat'
		    return
		  endif
		  write(142,'(a,3i10)')   '#posizione interna               = ', nint(center(1:3))
		  write(142,'(a,3i10)')   '#posizione esterna               = ', iprobe,jprobe,kprobe
		  write(142,'(a,i10)')    '#Numero punti letti              = ', myn
		  write(142,'(a,f20.10)') '#Sigma teorica                   = ', sigma
		  write(142,'(a,f20.10,a,f20.10)') '#Sigma numerica                  = ', avg,' +- ',std
		  write(142,'(a,f20.10,a,f20.10)') '#Press numerica                  = ', avgP,' +- ',stdP
		  write(142,'(a,f20.10,a,f20.10)') '#Radius numerica                 = ', avgR,' +- ',stdR
		  close(142)
       endif
#endif  


    
    endsubroutine

endmodule
