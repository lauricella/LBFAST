#include "defines.h"
module integrator_module
!   use mpi_template
#ifdef _OPENACC
   use openacc
#endif
   use mpi_template, only : nbuff,coords,myoffset,myrank,nprocs,intpbc_dir, &
           num_links_pops,links_pops,f_datampi,fvec_datampi,b_datampi,c_datampi,i_datampi, &
           f_send_extr,f_recv_extr,fvec_send_extr,fvec_recv_extr, &
           b_send_extr,b_recv_extr,c_send_extr,c_recv_extr, &
           f_send_buffmpi,f_recv_buffmpi,f_nbuffmpi_send,f_nbuffmpi_recv,lbuff, &
           fvec_send_buffmpi,fvec_recv_buffmpi,fvec_nbuffmpi_send,fvec_nbuffmpi_recv, &
           b_send_buffmpi,b_recv_buffmpi,b_nbuffmpi_send,b_nbuffmpi_recv, &
           c_send_buffmpi,c_recv_buffmpi,c_nbuffmpi_send,c_nbuffmpi_recv, &
           exchange_phifields_sendrecv,exchange_phifields_intpbc,exchange_phifields_wait, &
           exchange_auxfields_sendrecv,exchange_auxfields_intpbc,exchange_auxfields_wait, &
           exchange_hfields_sendrecv,exchange_hfields_intpbc,exchange_hfields_wait, &
           exchange_forces_sendrecv,exchange_forces_intpbc,exchange_forces_wait, &
#ifdef CRAY 
           int_buffpbc,nbuffpbc_int, &
#endif
           skip_myoffset,or_world_l, &
           dostop
           
           
   use prints, only : get_memory_gpu,print_memory_registration_gpu, &
    driver_print_raw_sync,driver_print_raw_sync, &
    init_output,driver_print_vtk_sync,read_restart_2c,read_restart_1c, &
#ifdef _OPENACC
    printdeviceproperties, &
#endif
    write_restart_1c,write_restart_2c,driver_print_raw_sync2d, &
    driver_print_raw_isfluid,copy_print
   use vars
   use bcs3D, only : bcs_mesoscopic_hfields,bcs_mesoscopic_phifields
#if defined(_OPENACC)        
   use lb_cuda_driver, only : moments_LB_cuda,fused_lb_cuda,test_LB_cuda, &
    compute_norm_interface_cuda,thinfilm_scan_mark_cuda,repulsive_flux_normal_cuda, &
    compute_div_theta_n,update_phifields,fused_LB_cuda_int,fused_LB_cuda_ext, &
    moments_LB_cuda,moments_LB_cuda_int,moments_LB_cuda_ext, &
    update_phifields_int,update_phifields_ext, &
    compute_norm_interface_cuda_int,compute_norm_interface_cuda_ext
#endif
   use profiling_m,   only : timer_init,itime_start, &
      startPreprocessingTime,print_timing_partial, &
      reset_timing_partial,printSimulationTime, &
      print_timing_final,itime_counter,idiagnostic, &
      ldiagnostic,start_timing2,end_timing2, &
      startSimulationTime,print_memory_registration, &
      get_memory,current_time, &
#ifdef CUDA
      get_totram,get_memory_cuda,print_memory_registration_cuda
#else
   get_totram
#endif

   implicit none

contains

   subroutine integrator
      implicit none
      integer :: subchords(3)
      logical :: ltime_actual=.false.
      integer :: ii,jj,kk
      integer :: xblock,yblock,zblock,myblock,iii,jjj,kkk

      step=0

      flip=mod(step,2)+1     
      flop = 3 - flip
 
      !$acc data copy(step,lx,ly,lz,nx,ny,nz,coords,myoffset,isfluid,myrank, &
      !$acc& rhoprint,velprint,radius,iprobe,jprobe,kprobe, &
	  !$acc& tau1,visc1,rho_r,rho_b,invrho_r,invrho_b,omega,arr_3d, &
      !$acc& intpbc_dir,num_links_pops,links_pops,f_datampi,uwall,udotc,uu, &
      !$acc& f_send_extr,f_recv_extr,openbc,openbc_type_x,openbc_type_y,openbc_type_z, &
      !$acc& ntothfields,ntotphifields,ntotauxfields,ntotlocauxfields,ntotforces, &
      !$acc& hfields_flip,hfields_flop,auxfields,locauxfields,forces, &
      !$acc& nblocks,nxblock,nxyblock,openbc_press_x,openbc_press_y,openbc_press_z, &
      !$acc& openbc_u_x,openbc_u_y,openbc_u_z,openbc_v_x,openbc_v_y,openbc_v_z, &
      !$acc& openbc_w_x,openbc_w_y,openbc_w_z, &
#ifdef TWOCOMPONENT
      !$acc& phifields_flip,phifields_flop, &
      !$acc& sharp_c,beta,tau2,visc2, kapp, sigma,width,tau_diff,corr,global_phi_sum_ini,global_count,global_phi_sum,&
	  !$acc& global_phi_change, &
#ifdef MONOD
	  !$acc&  mu_max,Ks, &
#endif
#ifdef IMPOSED_PRESSURE_GRADIENT
      !$acc& rhoIN,rhoOUT, &
#endif
#ifdef REPULSIVE_FLUX
	 !$acc& rep_mask,q_th,cosOppT,pwr,A_rep,win, &
#endif
#endif
      !$acc& flip,flop,stepskip,nxskip,nyskip,nzskip, &
#ifdef WRITEPRESS
      !$acc& pressprint, &
#endif
#ifdef MULTIHIT
	  !$acc& ABCx,ABCy,ABCz, &
#endif
#ifdef ELASTIC_FORCE
	  !$acc& u_ref,v_ref,w_ref, &
#endif
      !$acc& b_send_extr,b_recv_extr, &
      !$acc& c_send_extr,c_recv_extr, &
      !$acc& f_send_extr,f_recv_extr, &
#ifdef CRAY 
      !$acc& int_buffpbc,nbuffpbc_int, &
#endif
      !$acc& fvec_send_extr,fvec_recv_extr,fx,fy,fz) &
      !$acc& create(f_send_buffmpi, &
      !$acc& f_recv_buffmpi,fvec_send_buffmpi &
      !$acc& ,b_send_buffmpi,b_recv_buffmpi &
      !$acc& ,c_send_buffmpi,c_recv_buffmpi &
      !$acc& ,fvec_recv_buffmpi)
      ! quali sono i buff effettivamente da tenere? servono tutti!
      !$acc wait
      
#ifdef _OPENACC
      call printDeviceProperties(ngpus,devNum,devType,6)
#endif
      iframe=0
      iframe2D=0
      
      if(myrank==0)then
         write(6,'(a,i8,a,i8,3f16.4)')'start simulation'
         call flush(6)
      endif

      if(lprint)then
         call init_output(1,lvtk,lraw,nplanes,ndir,npoint)
         call string_char(head1,nheadervtk(1),headervtk(1))
         call string_char(head2,nheadervtk(2),headervtk(2))
#if defined(WRITEPRESS)
         call string_char(head3,nheadervtk(3),headervtk(3))
#endif
         if(lraw)then
           call driver_print_raw_isfluid(iframe)
         endif
      endif
      
      !***********************************read restart************************
#ifdef TWOCOMPONENT
      if(lrestart)then
         call read_restart_2c(iframe,iframe2D,hfields_flip,phifields_flip)
      endif
#else
      if(lrestart)then
         call read_restart_1c(iframe,iframe2D,hfields_flip)
      endif
#endif      

      if(lprint .and. (.not. lrestart))then
	     call copy_print(iframe,hfields_flip &
#ifdef TWOCOMPONENT	      
	      ,phifields_flip &
#endif
	      )
	      
         if(lvtk)then
            call driver_print_vtk_sync(iframe)
         endif
         if(lraw)then
            call driver_print_raw_sync(iframe)
            !2d planes print
            call driver_print_raw_sync2D(iframe2D,nplanes,ndir,npoint)
         endif
      endif
       
#ifdef TWOCOMPONENT	 
!****************scambio phi: boundary condition periodiche su phi************************
      if(ldiagnostic)call start_timing2("LB","ex_phifields_sendrecv")
	  call exchange_phifields_sendrecv(phifields_flip)
	  if(ldiagnostic)call end_timing2("LB","ex_phifields_sendrecv")
      if(ldiagnostic)call start_timing2("LB","ex_phifields_intpbc")
	  call exchange_phifields_intpbc(phifields_flip)
	  if(ldiagnostic)call end_timing2("LB","ex_phifields_intpbc")
      if(ldiagnostic)call start_timing2("LB","ex_phifields_wait")
	  call exchange_phifields_wait(phifields_flip)
	  if(ldiagnostic)call end_timing2("LB","ex_phifields_wait")
#endif
 
	  !***********************************pbcs boundary conditions ********************************!
      if(ldiagnostic)call start_timing2("LB","ex_hfields_sendrecv")
	  call exchange_hfields_sendrecv(hfields_flip)
	  if(ldiagnostic)call end_timing2("LB","ex_hfields_sendrecv")
	  if(ldiagnostic)call start_timing2("LB","ex_hfields_intpbc")
	  call exchange_hfields_intpbc(hfields_flip)
	  if(ldiagnostic)call end_timing2("LB","ex_hfields_intpbc")
	  if(ldiagnostic)call start_timing2("LB","ex_hfields_wait")
	  call exchange_hfields_wait(hfields_flip)
	  if(ldiagnostic)call end_timing2("LB","ex_hfields_wait")	  
      !************ thread-safe boundary condition setup
#ifdef TWOCOMPONENT  
      if(ldiagnostic)call start_timing2("LB","bcs_phi")
      call bcs_mesoscopic_phifields(hfields_flop,phifields_flop)
      if(ldiagnostic)call end_timing2("LB","bcs_phi") 
#endif
      ! start diagnostic if requested
      if(ldiagnostic)then
         !call print_timing_partial(1,1,itime_start,6)
         !call reset_timing_partial()
         call startSimulationTime()
         call get_memory(smemory)
         call get_totram(sram)
         call print_memory_registration(6,&
            'Occupied memory after setup MPI','Total memory',smemory,sram)
      endif
      !***********************************moments initialization************************

      call get_memory_gpu(mymemory,totmemory)
      call print_memory_registration_gpu(6,'DEVICE memory occupied at the start', &
         'total DEVICE memory',mymemory,totmemory)

      !write(6,'(a,i8,1x,i8,1x,a)')' sono arrivato ',myrank,__LINE__,__FILE__
      !call dostop
      time_init=current_time()
      time_actual_old=time_init
      call cpu_time(ts1)
      
#if 0
      call dostop('ciao',__FILE__,__LINE__)
      !call test_LB_cuda
      goto 110
#endif      
      do step_flip=1,nsteps,2
!***********************************************************************
!***********************************FLIP********************************
!***********************************************************************
         step=step+1
         flip=mod(step,2)+1     
         flop = 3 - flip
         !$acc wait
         !$acc update device(step,flip,flop)
         !$acc wait

#ifdef TWOCOMPONENT	 
#ifdef ASYNCMPI
         if(ldiagnostic)call start_timing2("LB","compute_norm_ext")
         call compute_norm_interface_cuda_ext(phifields_flip)
         if(ldiagnostic)call end_timing2("LB","compute_norm_ext")
         if(ldiagnostic)call start_timing2("LB","ex_auxfields_sendrecv")
		 call exchange_auxfields_sendrecv
		 if(ldiagnostic)call end_timing2("LB","ex_auxfields_sendrecv")
		 if(ldiagnostic)call start_timing2("LB","ex_auxfields_intpbc")
		 call exchange_auxfields_intpbc
		 if(ldiagnostic)call end_timing2("LB","ex_auxfields_intpbc")
		 if(ldiagnostic)call start_timing2("LB","compute_norm_int")
         call compute_norm_interface_cuda_int(phifields_flip)
         if(ldiagnostic)call end_timing2("LB","compute_norm_int")
         if(ldiagnostic)call start_timing2("LB","ex_auxfields_wait")
		 call exchange_auxfields_wait
		 if(ldiagnostic)call end_timing2("LB","ex_auxfields_wait")
#else
         !***********************************ora che ho phi in cornice calcolo normx normy normz************************
         if(ldiagnostic)call start_timing2("LB","compute_norm")
         call compute_norm_interface_cuda(phifields_flip)
         if(ldiagnostic)call end_timing2("LB","compute_norm")
         !***********************************scambio normx normy normz************************
         if(ldiagnostic)call start_timing2("LB","ex_auxfields_sendrecv")
		 call exchange_auxfields_sendrecv
		 if(ldiagnostic)call end_timing2("LB","ex_auxfields_sendrecv")
		 if(ldiagnostic)call start_timing2("LB","ex_auxfields_intpbc")
		 call exchange_auxfields_intpbc
		 if(ldiagnostic)call end_timing2("LB","ex_auxfields_intpbc")
		 if(ldiagnostic)call start_timing2("LB","ex_auxfields_wait")
		 call exchange_auxfields_wait
		 if(ldiagnostic)call end_timing2("LB","ex_auxfields_wait")
#endif         
#endif

         !************ thread-safe boundary condition setup
#ifdef TWOCOMPONENT           
         if(ldiagnostic)call start_timing2("LB","bcs_phi")
         call bcs_mesoscopic_phifields(hfields_flip,phifields_flip)
         if(ldiagnostic)call end_timing2("LB","bcs_phi") 
#endif
         !***********************************moments************************

#ifdef TWOCOMPONENT
         if(ldiagnostic)call start_timing2("LB","compute_div_theta_n")
         call compute_div_theta_n(phifields_flip)
         if(ldiagnostic)call end_timing2("LB","compute_div_theta_n")
         
#ifdef REPULSIVE_FLUX
         if(ldiagnostic)call start_timing2("LB","repulsive_flux")
		 call thinfilm_scan_mark_cuda(phifields_flip)
		 call repulsive_flux_normal_cuda(phifields_flip)
		 if(ldiagnostic)call end_timing2("LB","repulsive_flux")
#endif
#endif

#ifdef ASYNCMPI
         if(ldiagnostic)call start_timing2("LB","moments_LB_cuda_ext")
         call moments_LB_cuda_ext(hfields_flop,hfields_flip &
#ifdef TWOCOMPONENT         
          ,phifields_flip &
#endif
          )
         if(ldiagnostic)call end_timing2("LB","moments_LB_cuda_ext")
         
         if(ldiagnostic)call start_timing2("LB","ex_hfields_sendrecv")
		 call exchange_hfields_sendrecv(hfields_flip)
		 if(ldiagnostic)call end_timing2("LB","ex_hfields_sendrecv")
		 if(ldiagnostic)call start_timing2("LB","ex_forces_sendrecv")
		 call exchange_forces_sendrecv
		 if(ldiagnostic)call end_timing2("LB","ex_forces_sendrecv")
		 if(ldiagnostic)call start_timing2("LB","ex_hfields_intpbc")
		 call exchange_hfields_intpbc(hfields_flip)
		 if(ldiagnostic)call end_timing2("LB","ex_hfields_intpbc")
		 if(ldiagnostic)call start_timing2("LB","ex_forces_intpbc")
		 call exchange_forces_intpbc
		 if(ldiagnostic)call end_timing2("LB","ex_forces_intpbc")
		 
		 if(ldiagnostic)call start_timing2("LB","moments_LB_cuda_int")
         call moments_LB_cuda_int(hfields_flop,hfields_flip &
#ifdef TWOCOMPONENT         
          ,phifields_flip &
#endif
          )
         if(ldiagnostic)call end_timing2("LB","moments_LB_cuda_int")
         
         if(ldiagnostic)call start_timing2("LB","ex_hfields_wait")
		 call exchange_hfields_wait(hfields_flip)
		 if(ldiagnostic)call end_timing2("LB","ex_hfields_wait")
		 if(ldiagnostic)call start_timing2("LB","ex_forces_wait")
		 call exchange_forces_wait
		 if(ldiagnostic)call end_timing2("LB","ex_forces_wait")
#else
         if(ldiagnostic)call start_timing2("LB","moments_LB_cuda")
         call moments_LB_cuda(hfields_flop,hfields_flip &
#ifdef TWOCOMPONENT         
          ,phifields_flip &
#endif
          )
         if(ldiagnostic)call end_timing2("LB","moments_LB_cuda")
		 !***********************************pbcs boundary conditions ********************************!
         if(ldiagnostic)call start_timing2("LB","ex_hfields_sendrecv")
		 call exchange_hfields_sendrecv(hfields_flip)
		 if(ldiagnostic)call end_timing2("LB","ex_hfields_sendrecv")
		 if(ldiagnostic)call start_timing2("LB","ex_hfields_intpbc")
		 call exchange_hfields_intpbc(hfields_flip)
		 if(ldiagnostic)call end_timing2("LB","ex_hfields_intpbc")
		 if(ldiagnostic)call start_timing2("LB","ex_hfields_wait")
		 call exchange_hfields_wait(hfields_flip)
		 if(ldiagnostic)call end_timing2("LB","ex_hfields_wait")
		 !***********************************force pbcs boundary conditions ********************************!
         if(ldiagnostic)call start_timing2("LB","ex_forces_sendrecv")
		 call exchange_forces_sendrecv
		 if(ldiagnostic)call end_timing2("LB","ex_forces_sendrecv")
		 if(ldiagnostic)call start_timing2("LB","ex_forces_intpbc")
		 call exchange_forces_intpbc
		 if(ldiagnostic)call end_timing2("LB","ex_forces_intpbc")
		 if(ldiagnostic)call start_timing2("LB","ex_forces_wait")
		 call exchange_forces_wait
		 if(ldiagnostic)call end_timing2("LB","ex_forces_wait")      
#endif        
         if(lprint)then 
           if(mod(step,stamp).eq.0 .or. mod(step,stamp2D).eq.0 .or. mod(step,stamp_term).eq.0)then
             if(ldiagnostic)call start_timing2("IO","print")
	         call copy_print(iframe,hfields_flip &
#ifdef TWOCOMPONENT	      
	          ,phifields_flip &
#endif
	         )
             !write(6,*)'vorrei stampare'
             if(mod(step,stamp).eq.0)then
               iframe=iframe+1
               if(lvtk)then
                 call driver_print_vtk_sync(iframe)
               endif
               if(lraw)then
                 call driver_print_raw_sync(iframe)
               endif
             endif
             !***********************************Print on files 2D************************
             if(mod(step,stamp2D).eq.0 .and. lraw)then
               iframe2D=iframe2D+1
               call driver_print_raw_sync2D(iframe2D,nplanes,ndir,npoint)
             endif

             if(mod(step,stamp_term).eq.0)then
               time_actual=current_time()
               gi=iprobe;gj=jprobe;gk=kprobe
               subchords(1)=(gi-1)/nx
               subchords(2)=(gj-1)/ny
               subchords(3)=(gk-1)/nz
               if(all(subchords==coords))then
                 i=gi/stepskip-skip_myoffset(1)
                 j=gj/stepskip-skip_myoffset(2)
                 k=gk/stepskip-skip_myoffset(3)
                 write(6,'(a,4i6,f10.2,a,i2,9g16.8)')'stamp step : ',step, &
                  (gi/stepskip)*stepskip,(gj/stepskip)*stepskip,(gk/stepskip)*stepskip, &
                  (time_actual-time_actual_old)/real(stamp_term,kind=db)*real(nsteps-step,kind=db), &
                  '; probe values : ',isfluid(i,j,k),rhoprint(i,j,k),velprint(1:3,i,j,k)!, &
!#if defined(INTERNAL_OBSTACLES)                       
!                   corr,dphi,global_phi_sum_ini,global_phi_sum_new-global_phi_change_new,global_phi_sum_new !real(global_count_new)
!#else
!                   0.0e0_db                       
!#endif                       
                 call flush(6)
               endif
               time_actual_old=time_actual
             endif
             if(ldiagnostic)call end_timing2("IO","print")
           endif
         endif
       
         !***********************************collision + no slip + forcing: fused implementation*********
		 if(ldiagnostic)call start_timing2("LB","fused")
         call fused_LB_cuda(hfields_flip,hfields_flop &
#ifdef TWOCOMPONENT	   
          ,phifields_flip &
#endif
         )
         if(ldiagnostic)call end_timing2("LB","fused")
         
#ifdef TWOCOMPONENT	       
#ifdef ASYNCMPI 
         if(ldiagnostic)call start_timing2("LB","update_phifields_ext")       
         call update_phifields_ext(hfields_flip,phifields_flip,phifields_flop)
         if(ldiagnostic)call end_timing2("LB","update_phifields_ext")
         
         if(ldiagnostic)call start_timing2("LB","ex_phifields_sendrecv")
		 call exchange_phifields_sendrecv(phifields_flop)
		 if(ldiagnostic)call end_timing2("LB","ex_phifields_sendrecv")
		 
         if(ldiagnostic)call start_timing2("LB","ex_phifields_intpbc")
		 call exchange_phifields_intpbc(phifields_flop)
		 if(ldiagnostic)call end_timing2("LB","ex_phifields_intpbc")
		 
		 if(ldiagnostic)call start_timing2("LB","update_phifields_int")       
         call update_phifields_int(hfields_flip,phifields_flip,phifields_flop)
         if(ldiagnostic)call end_timing2("LB","update_phifields_int")
		 
		 if(ldiagnostic)call start_timing2("LB","ex_phifields_wait")
		 call exchange_phifields_wait(phifields_flop)
		 if(ldiagnostic)call end_timing2("LB","ex_phifields_wait")
#else  
         if(ldiagnostic)call start_timing2("LB","update_phifields")       
         call update_phifields(hfields_flip,phifields_flip,phifields_flop)
         if(ldiagnostic)call end_timing2("LB","update_phifields")
         
!****************scambio phi: boundary condition periodiche su phi************************
		 if(ldiagnostic)call start_timing2("LB","ex_phifields_sendrecv")
		 call exchange_phifields_sendrecv(phifields_flop)
		 if(ldiagnostic)call end_timing2("LB","ex_phifields_sendrecv")
		 
         if(ldiagnostic)call start_timing2("LB","ex_phifields_intpbc")
		 call exchange_phifields_intpbc(phifields_flop)
		 if(ldiagnostic)call end_timing2("LB","ex_phifields_intpbc")
		 
         if(ldiagnostic)call start_timing2("LB","ex_phifields_wait")
		 call exchange_phifields_wait(phifields_flop)
		 if(ldiagnostic)call end_timing2("LB","ex_phifields_wait")
#endif
#endif
         
         !************ thread-safe boundary condition setup
#ifdef INTERNAL_OBSTACLES    
         if(ldiagnostic)call start_timing2("LB","bcs_hfields")
         call bcs_mesoscopic_hfields(hfields_flip,hfields_flop &
#ifdef TWOCOMPONENT	          
          ,phifields_flip &
#endif
         )
         if(ldiagnostic)call end_timing2("LB","bcs_hfields") 
#endif         
         if(time_limit>ZERO)then
           if(mod(step,every_time_check).eq.0)then
             !$acc wait
             if(ldiagnostic)call start_timing2("IO","time_limit")
             time_actual = (current_time()) - time_init
             ltime_actual=(time_actual>time_limit)
             call or_world_l(ltime_actual)
             if(ltime_actual)then
               if(myrank==0)then
                 write(6,'(a,f16.2,a,i16,a,f16.2,a)')'Time limit ',time_limit, &
                  ' sec reached at step ',step,' after ',time_actual,' sec'
                 call flush(6)
               endif
               lwriterestart=.true.
               if(ldiagnostic)call end_timing2("IO","time_limit")
               goto 110
             endif
             if(ldiagnostic)call end_timing2("IO","time_limit")
           endif
         endif
	 
!***********************************************************************
!***********************************FLOP********************************
!***********************************************************************
         step=step+1
         flip=mod(step,2)+1     
         flop = 3 - flip
         !$acc wait
         !$acc update device(step,flip,flop)
         !$acc wait
         
#ifdef TWOCOMPONENT	 
#ifdef ASYNCMPI
         if(ldiagnostic)call start_timing2("LB","compute_norm_ext")
         call compute_norm_interface_cuda_ext(phifields_flop)
         if(ldiagnostic)call end_timing2("LB","compute_norm_ext")
         if(ldiagnostic)call start_timing2("LB","ex_auxfields_sendrecv")
		 call exchange_auxfields_sendrecv
		 if(ldiagnostic)call end_timing2("LB","ex_auxfields_sendrecv")
		 if(ldiagnostic)call start_timing2("LB","ex_auxfields_intpbc")
		 call exchange_auxfields_intpbc
		 if(ldiagnostic)call end_timing2("LB","ex_auxfields_intpbc")
		 if(ldiagnostic)call start_timing2("LB","compute_norm_int")
         call compute_norm_interface_cuda_int(phifields_flop)
         if(ldiagnostic)call end_timing2("LB","compute_norm_int")
         if(ldiagnostic)call start_timing2("LB","ex_auxfields_wait")
		 call exchange_auxfields_wait
		 if(ldiagnostic)call end_timing2("LB","ex_auxfields_wait")
#else
         !***********************************ora che ho phi in cornice calcolo normx normy normz************************
         if(ldiagnostic)call start_timing2("LB","compute_norm")
         call compute_norm_interface_cuda(phifields_flop)
         if(ldiagnostic)call end_timing2("LB","compute_norm")
         !***********************************scambio normx normy normz************************
         if(ldiagnostic)call start_timing2("LB","ex_auxfields_sendrecv")
		 call exchange_auxfields_sendrecv
		 if(ldiagnostic)call end_timing2("LB","ex_auxfields_sendrecv")
		 if(ldiagnostic)call start_timing2("LB","ex_auxfields_intpbc")
		 call exchange_auxfields_intpbc
		 if(ldiagnostic)call end_timing2("LB","ex_auxfields_intpbc")
		 if(ldiagnostic)call start_timing2("LB","ex_auxfields_wait")
		 call exchange_auxfields_wait
		 if(ldiagnostic)call end_timing2("LB","ex_auxfields_wait")
#endif         
#endif
 
         !************ thread-safe boundary condition setup
#ifdef TWOCOMPONENT         
         if(ldiagnostic)call start_timing2("LB","bcs_phi")
         call bcs_mesoscopic_phifields(hfields_flop,phifields_flop)
         if(ldiagnostic)call end_timing2("LB","bcs_phi") 
#endif
         !***********************************moments************************
         
#ifdef TWOCOMPONENT
         if(ldiagnostic)call start_timing2("LB","compute_div_theta_n")
         call compute_div_theta_n(phifields_flop)
         if(ldiagnostic)call end_timing2("LB","compute_div_theta_n")
         
#ifdef REPULSIVE_FLUX
         if(ldiagnostic)call start_timing2("LB","repulsive_flux")
		 call thinfilm_scan_mark_cuda(phifields_flop)
		 call repulsive_flux_normal_cuda(phifields_flop)
		 if(ldiagnostic)call end_timing2("LB","repulsive_flux")
#endif
#endif

#ifdef ASYNCMPI
         if(ldiagnostic)call start_timing2("LB","moments_LB_cuda_ext")
         call moments_LB_cuda_ext(hfields_flip,hfields_flop &
#ifdef TWOCOMPONENT         
          ,phifields_flop &
#endif
          )
         if(ldiagnostic)call end_timing2("LB","moments_LB_cuda_ext")	
         
         if(ldiagnostic)call start_timing2("LB","ex_hfields_sendrecv")
		 call exchange_hfields_sendrecv(hfields_flop)
		 if(ldiagnostic)call end_timing2("LB","ex_hfields_sendrecv")
		 if(ldiagnostic)call start_timing2("LB","ex_forces_sendrecv")
		 call exchange_forces_sendrecv
		 if(ldiagnostic)call end_timing2("LB","ex_forces_sendrecv")
		 if(ldiagnostic)call start_timing2("LB","ex_hfields_intpbc")
		 call exchange_hfields_intpbc(hfields_flop)
		 if(ldiagnostic)call end_timing2("LB","ex_hfields_intpbc")
		 if(ldiagnostic)call start_timing2("LB","ex_forces_intpbc")
		 call exchange_forces_intpbc
		 if(ldiagnostic)call end_timing2("LB","ex_forces_intpbc")
		 
         if(ldiagnostic)call start_timing2("LB","moments_LB_cuda_int")
         call moments_LB_cuda_int(hfields_flip,hfields_flop &
#ifdef TWOCOMPONENT         
          ,phifields_flop &
#endif
          )
         if(ldiagnostic)call end_timing2("LB","moments_LB_cuda_int")
         
         if(ldiagnostic)call start_timing2("LB","ex_hfields_wait")
		 call exchange_hfields_wait(hfields_flop)
		 if(ldiagnostic)call end_timing2("LB","ex_hfields_wait")
		 if(ldiagnostic)call start_timing2("LB","ex_forces_wait")
		 call exchange_forces_wait
		 if(ldiagnostic)call end_timing2("LB","ex_forces_wait")
#else

         if(ldiagnostic)call start_timing2("LB","moments_LB_cuda")
         call moments_LB_cuda(hfields_flip,hfields_flop &
#ifdef TWOCOMPONENT         
          ,phifields_flop &
#endif
          )
         if(ldiagnostic)call end_timing2("LB","moments_LB_cuda")	 
		 !***********************************pbcs boundary conditions ********************************!
         if(ldiagnostic)call start_timing2("LB","ex_hfields_sendrecv")
		 call exchange_hfields_sendrecv(hfields_flop)
		 if(ldiagnostic)call end_timing2("LB","ex_hfields_sendrecv")
		 if(ldiagnostic)call start_timing2("LB","ex_hfields_intpbc")
		 call exchange_hfields_intpbc(hfields_flop)
		 if(ldiagnostic)call end_timing2("LB","ex_hfields_intpbc")
		 if(ldiagnostic)call start_timing2("LB","ex_hfields_wait")
		 call exchange_hfields_wait(hfields_flop)
		 if(ldiagnostic)call end_timing2("LB","ex_hfields_wait")
		 !***********************************force pbcs boundary conditions ********************************!
         if(ldiagnostic)call start_timing2("LB","ex_forces_sendrecv")
		 call exchange_forces_sendrecv
		 if(ldiagnostic)call end_timing2("LB","ex_forces_sendrecv")
		 if(ldiagnostic)call start_timing2("LB","ex_forces_intpbc")
		 call exchange_forces_intpbc
		 if(ldiagnostic)call end_timing2("LB","ex_forces_intpbc")
		 if(ldiagnostic)call start_timing2("LB","ex_forces_wait")
		 call exchange_forces_wait
		 if(ldiagnostic)call end_timing2("LB","ex_forces_wait")  
#endif        
         if(lprint)then 
           if(mod(step,stamp).eq.0 .or. mod(step,stamp2D).eq.0 .or. mod(step,stamp_term).eq.0)then
             if(ldiagnostic)call start_timing2("IO","print")
	         call copy_print(iframe,hfields_flop &
#ifdef TWOCOMPONENT	      
	          ,phifields_flop &
#endif
	         )
             !write(6,*)'vorrei stampare'
             if(mod(step,stamp).eq.0)then
               iframe=iframe+1
               if(lvtk)then
                 call driver_print_vtk_sync(iframe)
               endif
               if(lraw)then
                 call driver_print_raw_sync(iframe)
               endif
             endif
             !***********************************Print on files 2D************************
             if(mod(step,stamp2D).eq.0 .and. lraw)then
               iframe2D=iframe2D+1
               call driver_print_raw_sync2D(iframe2D,nplanes,ndir,npoint)
             endif

             if(mod(step,stamp_term).eq.0)then
               time_actual=current_time()
               gi=iprobe;gj=jprobe;gk=kprobe
               subchords(1)=(gi-1)/nx
               subchords(2)=(gj-1)/ny
               subchords(3)=(gk-1)/nz
               if(all(subchords==coords))then
                 i=gi/stepskip-skip_myoffset(1)
                 j=gj/stepskip-skip_myoffset(2)
                 k=gk/stepskip-skip_myoffset(3)
                 write(6,'(a,4i6,f10.2,a,i2,9g16.8)')'stamp step : ',step, &
                  (gi/stepskip)*stepskip,(gj/stepskip)*stepskip,(gk/stepskip)*stepskip, &
                  (time_actual-time_actual_old)/real(stamp_term,kind=db)*real(nsteps-step,kind=db), &
                  '; probe values : ',isfluid(i,j,k),rhoprint(i,j,k),velprint(1:3,i,j,k)!, &
!#if defined(INTERNAL_OBSTACLES)                       
!                   corr,dphi,global_phi_sum_ini,global_phi_sum_new-global_phi_change_new,global_phi_sum_new !real(global_count_new)
!#else
!                   0.0e0                       
!#endif                       
                 call flush(6)
               endif
               time_actual_old=time_actual
             endif
             if(ldiagnostic)call end_timing2("IO","print")
           endif
         endif
         
         !***********************************collision + no slip + forcing: fused implementation*********
		 if(ldiagnostic)call start_timing2("LB","fused")  
		 call fused_LB_cuda(hfields_flop,hfields_flip &
#ifdef TWOCOMPONENT	            
         ,phifields_flop &
#endif
         )
         if(ldiagnostic)call end_timing2("LB","fused")
         

#ifdef TWOCOMPONENT	       
#ifdef ASYNCMPI

         if(ldiagnostic)call start_timing2("LB","update_phifields_ext")  
         call update_phifields_ext(hfields_flop,phifields_flop,phifields_flip)
         if(ldiagnostic)call end_timing2("LB","update_phifields_ext")
         
         if(ldiagnostic)call start_timing2("LB","ex_phifields_sendrecv")
		 call exchange_phifields_sendrecv(phifields_flip)
		 if(ldiagnostic)call end_timing2("LB","ex_phifields_sendrecv")
		 
         if(ldiagnostic)call start_timing2("LB","ex_phifields_intpbc")
		 call exchange_phifields_intpbc(phifields_flip)
		 if(ldiagnostic)call end_timing2("LB","ex_phifields_intpbc")
         
         if(ldiagnostic)call start_timing2("LB","update_phifields_int")  
         call update_phifields_int(hfields_flop,phifields_flop,phifields_flip)
         if(ldiagnostic)call end_timing2("LB","update_phifields_int")
         
         if(ldiagnostic)call start_timing2("LB","ex_phifields_wait")
		 call exchange_phifields_wait(phifields_flip)
		 if(ldiagnostic)call end_timing2("LB","ex_phifields_wait")
#else 
         if(ldiagnostic)call start_timing2("LB","update_phifields")  
         call update_phifields(hfields_flop,phifields_flop,phifields_flip)
         if(ldiagnostic)call end_timing2("LB","update_phifields")
         
!****************scambio phi: boundary condition periodiche su phi************************
		 if(ldiagnostic)call start_timing2("LB","ex_phifields_sendrecv")
		 call exchange_phifields_sendrecv(phifields_flip)
		 if(ldiagnostic)call end_timing2("LB","ex_phifields_sendrecv")
		 
         if(ldiagnostic)call start_timing2("LB","ex_phifields_intpbc")
		 call exchange_phifields_intpbc(phifields_flip)
		 if(ldiagnostic)call end_timing2("LB","ex_phifields_intpbc")
		 
         if(ldiagnostic)call start_timing2("LB","ex_phifields_wait")
		 call exchange_phifields_wait(phifields_flip)
		 if(ldiagnostic)call end_timing2("LB","ex_phifields_wait")
#endif
#endif
         
         
         !************ thread-safe boundary condition setup
#ifdef INTERNAL_OBSTACLES          
         if(ldiagnostic)call start_timing2("LB","bcs_hfields")
         call bcs_mesoscopic_hfields(hfields_flop,hfields_flip &
#ifdef TWOCOMPONENT	          
          ,phifields_flop &
#endif
         )
         if(ldiagnostic)call end_timing2("LB","bcs_hfields") 
#endif         
         if(time_limit>ZERO)then
           if(mod(step,every_time_check).eq.0)then
             !$acc wait
             if(ldiagnostic)call start_timing2("IO","time_limit")
             time_actual = (current_time()) - time_init
             ltime_actual=(time_actual>time_limit)
             call or_world_l(ltime_actual)
             if(ltime_actual)then
               if(myrank==0)then
                 write(6,'(a,f16.2,a,i16,a,f16.2,a)')'Time limit ',time_limit, &
                  ' sec reached at step ',step,' after ',time_actual,' sec'
                 call flush(6)
               endif
               lwriterestart=.true.
               if(ldiagnostic)call end_timing2("IO","time_limit")
               goto 110
             endif
             if(ldiagnostic)call end_timing2("IO","time_limit")
           endif
         endif
	 
      enddo
110   continue      
      call cpu_time(ts2)
      !***********************************write restart************************
      if(lwriterestart)then
#ifdef TWOCOMPONENT
      call write_restart_2c(iframe,iframe2D,hfields_flop,phifields_flop)
#else      
      call write_restart_1c(iframe,iframe2D,hfields_flop)
#endif
      endif
120   continue
	  call get_memory_gpu(mymemory,totmemory)
	  call print_memory_registration_gpu(6,'DEVICE memory occupied at the end', &
      'total DEVICE memory',mymemory,totmemory)
      !$wait
      !$acc end data


   endsubroutine

endmodule
