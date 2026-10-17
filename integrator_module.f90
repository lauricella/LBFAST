#include "defines.h"
module integrator_module
!   use mpi_template
#ifdef _OPENACC
   use openacc
#endif
   use iso_c_binding
#ifdef _NVML
   use nvml_interface
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
           dostop,sum_world_float
           
   use prints, only : get_memory_gpu,print_memory_registration_gpu, &
    driver_print_raw_sync,driver_print_raw_sync, &
    init_output,driver_print_vtk_sync,read_restart_2c,read_restart_1c, &
#ifdef _OPENACC
    printdeviceproperties, &
#endif
    write_restart_1c,write_restart_2c,driver_print_raw_sync2d, &
    driver_print_raw_isfluid,copy_print
#ifdef FUSED_MOMENTS
   use prints, only : lprint_moments_update
#endif
   use vars
   use bcs3D, only : bcs_mesoscopic_hfields,bcs_mesoscopic_phifields
#if defined(_OPENACC)        
   use lb_cuda_driver, only : moments_LB_cuda,fused_lb_cuda,test_LB_cuda, &
    compute_norm_interface_cuda,thinfilm_scan_mark_cuda,repulsive_flux_normal_cuda, &
    compute_div_theta_n,update_phifields, &
#ifdef CSF
    compute_csf_force, &
#endif
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
   use stat_module, only: Ekin0,Ekin, &
#if defined(LAMBTEST) && defined(TWOCOMPONENT)
    lamb_z,lamb_cm_z,lamb_x,lamb_cm_x, &
    lamb_dosc,lamb_A,pos_x_int_left,pos_x_int_node_left,pos_x_int_right, &
    pos_x_int_node_right,pos_z_int_left,pos_z_int_node_left, &
    pos_z_int_right,pos_z_int_node_right, &
    open_lamb,print_lamb, &
#elif defined(LAPLACE)
    lamb_x,lamb_cm_x, &
    pos_x_int_left,pos_x_int_node_left,pos_x_int_right, &
    pos_x_int_node_right, &
    open_laplace,print_laplace, &
    rhophi_in,pstar_in,p_in, &
    rhophi_out,pstar_out,p_out, &
    delta_p,laplace_rad, &
    sigma_eff, &  
#elif defined(CAPILLARYWAVE)
    print_capillary_wave,print_capillary_force_projection, &
#endif
    probe_loc,open_taylorgreen,print_taylorgreen,close_benchmark
    
   implicit none

contains

   !---------------------------------------------------------------------------------------------
   ! integrator: runs the time loop. It opens an OpenACC data region that keeps the fields on the
   ! GPU, writes the initial output (or reads the restart files), fills the halos of the initial
   ! fields, performs the time steps (see the description before the loop), and writes the
   ! restart files when lwriterestart is set (in the input, or when time_limit is reached).
   !---------------------------------------------------------------------------------------------
   subroutine integrator
      implicit none
      integer :: subchords(3)
      logical :: ltime_actual=.false.
      integer :: ii,jj,kk
      integer :: xblock,yblock,zblock,myblock,iii,jjj,kkk
      integer(c_int) :: ierrc

      step=0

      flip=mod(step,2)+1     
      flop = 3 - flip
 
      ! OpenACC data region for the whole run: the arrays in copy() are copied to the GPU here and
      ! back to the host at the end of the region; the MPI buffers in create() are allocated on the
      ! GPU without a copy. Kernels and halo exchanges work on the device copies.
      !$acc data copy(step,lx,ly,lz,nx,ny,nz,coords,myoffset,isfluid,myrank, &
      !$acc& rhoprint,velprint,radius,width,center,iprobe,jprobe,kprobe, &
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
      !$acc& sharp_c,beta,tau2,visc2,kapp,sigma,tau_diff,corr,global_phi_sum_ini,global_count,global_phi_sum,&
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
      ! which buffers really have to be kept? All of them are needed!
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
      ! The restart files provide the moments (and phi) in the flip buffers and the frame counters
      ! iframe and iframe2D; the step counter starts again from zero.
#ifdef TWOCOMPONENT
      if(lrestart)then
         call read_restart_2c(iframe,iframe2D,hfields_flip,phifields_flip)
      endif
#else
      if(lrestart)then
         call read_restart_1c(iframe,iframe2D,hfields_flip)
      endif
#endif      

      !************ thread-safe boundary condition setup
#ifdef TWOCOMPONENT  
      if(ldiagnostic)call start_timing2("LB","bcs_phi")
      call bcs_mesoscopic_phifields(hfields_flip,phifields_flip)
      if(ldiagnostic)call end_timing2("LB","bcs_phi") 
#endif

      ! Initial output (frame 0, not after a restart) and opening of the benchmark files on unit 142.
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
#if defined(TAYLORGREEN) && !defined(TWOCOMPONENT)
         call open_taylorgreen(142)
#endif
#if defined(LAMBTEST) && defined(TWOCOMPONENT)
         call open_lamb(142)  
#endif
#ifdef LAPLACE
         call open_laplace(142)  
#endif
#ifdef CAPILLARYWAVE
         call print_capillary_wave(142,.true.)
#endif
      endif
       
      ! Fill the halos of the initial phase field and moments before the first time step.
#ifdef TWOCOMPONENT	 
!****************phi exchange: periodic boundary conditions on phi************************
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

      !write(6,'(a,i8,1x,i8,1x,a)')' reached ',myrank,__LINE__,__FILE__
      !call dostop
      ! time_init is the reference for time_limit, time_actual_old for the estimate of the
      ! remaining time printed with the probe; ts1 and ts2 bracket the time loop.
      time_init=current_time()
      time_actual_old=time_init
#if defined(_NVML) && !defined(GETPOWER)      
      ierrc=get_gpu_energy_mJ_u64(mydev_c,energy_1)
#endif  
      ts1=current_time()
      !call cpu_time(ts1)
    
#if 0
      call dostop('debug',__FILE__,__LINE__)
      !call test_LB_cuda
      goto 110
#endif      
#ifdef FUSED_MOMENTS
      ! from here on copy_print adds the update of moments_LB_kernel to the printed velocity
      lprint_moments_update=.true.
#endif
      !--------------------------------------------------------------------------------------------
      ! Time loop. Each iteration performs two time steps with double buffering:
      !  - FLIP half (odd step): the kernels read hfields_flip and phifields_flip, and the new
      !    moments and phase field are written into hfields_flop and phifields_flop;
      !  - FLOP half (even step): the same sequence with the roles of the two buffers exchanged.
      ! After a complete iteration the newest fields are again in the flip buffers. If nsteps is
      ! odd, the last iteration also performs step nsteps+1.
      !
      ! Sequence of one time step (FLIP half: 'in' = flip buffers, 'out' = flop buffers):
      !  1. TWOCOMPONENT: unit normal, |grad phi|, phi (1-phi) n and Laplacian of phi
      !     (compute_norm_interface) and halo exchange of auxfields; boundary conditions on phi
      !     (bcs_mesoscopic_phifields); div(phi (1-phi) n) (compute_div_theta_n); with CSF the
      !     surface-tension force (compute_csf_force); with REPULSIVE_FLUX the repulsive flux.
      !  2. Unless FUSED_MOMENTS is defined: moments_LB_cuda computes the force on every fluid node,
      !     adds the forcing terms to the velocity and the second-order moments of hfields 'in' in
      !     place, and stores the force in forces. The 'out' buffer, which still holds the moments of
      !     the previous time step, is passed only to read the velocity of that step.
      !  3. Halo exchange of hfields 'in' (and, unless FUSED_MOMENTS, of forces), since the fused
      !     kernel reads the halo nodes.
      !  4. Output (lprint) every stamp steps (3D fields), stamp2D steps (2D planes) and stamp_term
      !     steps (benchmark data and probe line), from the 'in' buffers.
      !  5. fused_LB_cuda: collision, forcing and streaming; the new moments go to hfields 'out'.
      !  6. TWOCOMPONENT: update_phifields computes phi 'out' from phi 'in' and the velocity of
      !     hfields 'in', followed by the halo exchange of phi 'out'.
      !  7. INTERNAL_OBSTACLES: boundary conditions on hfields 'out' (bcs_mesoscopic_hfields).
      !  8. If time_limit > 0, every every_time_check steps: when the time measured by current_time
      !     exceeds time_limit, lwriterestart is set and the loop is left (goto 110).
      !
      ! With ASYNCMPI, the halo exchanges of steps 1, 3 and 6 overlap with computation: the _ext
      ! kernels first update the tiles at the faces of the subdomain, the exchange is started
      ! (sendrecv packs the halo data and posts non-blocking MPI transfers, intpbc copies the
      ! periodic images that stay on the same rank), the _int kernels update the interior tiles,
      ! and wait completes the exchange and unpacks the received halo. With FUSED_MOMENTS there is
      ! no moment kernel to overlap with the exchange of step 3.
      !
      ! FUSED_MOMENTS (single-component flows without INTERNAL_OBSTACLES, MULTIHIT and
      ! ELASTIC_FORCE, see defines.h): moments_LB_cuda and the exchange of forces are skipped,
      ! because the fused kernels apply the moment update themselves when they read the moments of
      ! a node; the stored moments are therefore those before the update (docs/FUSED_KERNELS.md).
      !--------------------------------------------------------------------------------------------
      do step_flip=1,nsteps,2
!***********************************************************************
!***********************************FLIP********************************
!***********************************************************************
         step=step+1
         flip=mod(step,2)+1     
         flop = 3 - flip
         ! step, flip and flop are passed to the kernels with host_data use_device, so their device
         ! copies are updated at every step. The field buffers are not indexed by flip and flop: each
         ! call passes hfields_flip/hfields_flop (and the phi buffers) explicitly.
         !$acc wait
         !$acc update device(step,flip,flop)
         !$acc wait

         ! Step 1 (TWOCOMPONENT): interface geometry, halo exchange of auxfields, capillary terms.
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
         !***********************************now that phi is available in the halo, compute normx normy normz************************
         if(ldiagnostic)call start_timing2("LB","compute_norm")
         call compute_norm_interface_cuda(phifields_flip)
         if(ldiagnostic)call end_timing2("LB","compute_norm")
         !***********************************exchange normx normy normz************************
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
#ifdef CSF
         if(ldiagnostic)call start_timing2("LB","compute_csf_force")
         call compute_csf_force
         if(ldiagnostic)call end_timing2("LB","compute_csf_force")
#ifdef CAPILLARYWAVE
         if(step==1)call print_capillary_force_projection
#endif
#endif
         
#ifdef REPULSIVE_FLUX
         if(ldiagnostic)call start_timing2("LB","repulsive_flux")
		 call thinfilm_scan_mark_cuda(phifields_flip)
		 call repulsive_flux_normal_cuda(phifields_flip)
		 if(ldiagnostic)call end_timing2("LB","repulsive_flux")
#endif
#endif

         ! Steps 2-3: moment update (unless FUSED_MOMENTS) and halo exchange of hfields_flip and forces.
#ifdef ASYNCMPI
#ifndef FUSED_MOMENTS
         if(ldiagnostic)call start_timing2("LB","moments_LB_cuda_ext")
         call moments_LB_cuda_ext(hfields_flop,hfields_flip &
#ifdef TWOCOMPONENT         
          ,phifields_flip &
#endif
          )
         if(ldiagnostic)call end_timing2("LB","moments_LB_cuda_ext")
#endif
         
         if(ldiagnostic)call start_timing2("LB","ex_hfields_sendrecv")
		 call exchange_hfields_sendrecv(hfields_flip)
		 if(ldiagnostic)call end_timing2("LB","ex_hfields_sendrecv")
#ifndef FUSED_MOMENTS
		 if(ldiagnostic)call start_timing2("LB","ex_forces_sendrecv")
		 call exchange_forces_sendrecv
		 if(ldiagnostic)call end_timing2("LB","ex_forces_sendrecv")
#endif
		 if(ldiagnostic)call start_timing2("LB","ex_hfields_intpbc")
		 call exchange_hfields_intpbc(hfields_flip)
		 if(ldiagnostic)call end_timing2("LB","ex_hfields_intpbc")
#ifndef FUSED_MOMENTS
		 if(ldiagnostic)call start_timing2("LB","ex_forces_intpbc")
		 call exchange_forces_intpbc
		 if(ldiagnostic)call end_timing2("LB","ex_forces_intpbc")
#endif
		 
#ifndef FUSED_MOMENTS
		 if(ldiagnostic)call start_timing2("LB","moments_LB_cuda_int")
         call moments_LB_cuda_int(hfields_flop,hfields_flip &
#ifdef TWOCOMPONENT         
          ,phifields_flip &
#endif
          )
         if(ldiagnostic)call end_timing2("LB","moments_LB_cuda_int")
#endif
         
         if(ldiagnostic)call start_timing2("LB","ex_hfields_wait")
		 call exchange_hfields_wait(hfields_flip)
		 if(ldiagnostic)call end_timing2("LB","ex_hfields_wait")
#ifndef FUSED_MOMENTS
		 if(ldiagnostic)call start_timing2("LB","ex_forces_wait")
		 call exchange_forces_wait
		 if(ldiagnostic)call end_timing2("LB","ex_forces_wait")
#endif
#else
#ifndef FUSED_MOMENTS
         if(ldiagnostic)call start_timing2("LB","moments_LB_cuda")
         call moments_LB_cuda(hfields_flop,hfields_flip &
#ifdef TWOCOMPONENT         
          ,phifields_flip &
#endif
          )
         if(ldiagnostic)call end_timing2("LB","moments_LB_cuda")
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
		 !***********************************force pbcs boundary conditions ********************************!
#ifndef FUSED_MOMENTS
         if(ldiagnostic)call start_timing2("LB","ex_forces_sendrecv")
		 call exchange_forces_sendrecv
		 if(ldiagnostic)call end_timing2("LB","ex_forces_sendrecv")
#endif
#ifndef FUSED_MOMENTS
		 if(ldiagnostic)call start_timing2("LB","ex_forces_intpbc")
		 call exchange_forces_intpbc
		 if(ldiagnostic)call end_timing2("LB","ex_forces_intpbc")
#endif
#ifndef FUSED_MOMENTS
		 if(ldiagnostic)call start_timing2("LB","ex_forces_wait")
		 call exchange_forces_wait
		 if(ldiagnostic)call end_timing2("LB","ex_forces_wait")      
#endif
#endif        
         ! Step 4: output. copy_print copies the fields of the flip buffers to the print arrays.
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
#if defined(TAYLORGREEN) && !defined(TWOCOMPONENT)
                 call print_taylorgreen(142)
#endif
#if defined(LAMBTEST) && defined(TWOCOMPONENT)
                 call print_lamb(142)
#endif
#ifdef LAPLACE
                 call print_laplace(142)
#endif
#ifdef CAPILLARYWAVE
                 call print_capillary_wave(142,.false.)
#endif

               ! Probe line, printed by the rank that owns the probe node (iprobe,jprobe,kprobe): step,
               ! probe position on the output grid, estimated remaining time in seconds (time
               ! per step over the last stamp_term steps times the remaining steps), and isfluid, rhoprint
               ! and velprint at the probe. With _NVML and GETPOWER, the power sampled in fused_LB_cuda is
               ! integrated over the elapsed interval.
               time_actual=current_time()
               gi=iprobe;gj=jprobe;gk=kprobe
               subchords(1)=(gi-1)/nx
               subchords(2)=(gj-1)/ny
               subchords(3)=(gk-1)/nz
#if defined(_NVML) && defined(GETPOWER)              
               num_p_w=num_p_w+1
               p_w=real(p_mw,kind=db)*1.0e-3_db
               tot_energy=tot_energy+p_w*(time_actual-time_actual_old)
#endif 
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
         ! Step 5: reads hfields_flip (including its halo) and writes the new moments into hfields_flop.
		 if(ldiagnostic)call start_timing2("LB","fused")
		 !call probe_loc(hfields_flip,'hfields_flip',__FILE__,__LINE__)
         call fused_LB_cuda(hfields_flip,hfields_flop &
#ifdef TWOCOMPONENT	   
          ,phifields_flip &
#endif
         )
         if(ldiagnostic)call end_timing2("LB","fused")
         
         ! Step 6 (TWOCOMPONENT): new phase field in phifields_flop, then its halo exchange.
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
         
!****************phi exchange: periodic boundary conditions on phi************************
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
         ! Step 7 (INTERNAL_OBSTACLES): boundary conditions on the new moments in hfields_flop.
#ifdef INTERNAL_OBSTACLES    
         if(ldiagnostic)call start_timing2("LB","bcs_hfields")
         call bcs_mesoscopic_hfields(hfields_flip,hfields_flop &
#ifdef TWOCOMPONENT	          
          ,phifields_flip &
#endif
         )
         if(ldiagnostic)call end_timing2("LB","bcs_hfields") 
#endif         
         ! Step 8: time limit. or_world_l combines the test of all ranks, so that all of them
         ! leave the loop at the same step; the restart files are then written (lwriterestart).
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
         ! Second time step of the iteration: the same sequence as the FLIP half, with hfields_flop and
         ! phifields_flop as input and hfields_flip and phifields_flip as output.
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
         !***********************************now that phi is available in the halo, compute normx normy normz************************
         if(ldiagnostic)call start_timing2("LB","compute_norm")
         call compute_norm_interface_cuda(phifields_flop)
         if(ldiagnostic)call end_timing2("LB","compute_norm")
         !***********************************exchange normx normy normz************************
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
#ifdef CSF
         if(ldiagnostic)call start_timing2("LB","compute_csf_force")
         call compute_csf_force
         if(ldiagnostic)call end_timing2("LB","compute_csf_force")
#ifdef CAPILLARYWAVE
         if(step==1)call print_capillary_force_projection
#endif
#endif
         
#ifdef REPULSIVE_FLUX
         if(ldiagnostic)call start_timing2("LB","repulsive_flux")
		 call thinfilm_scan_mark_cuda(phifields_flop)
		 call repulsive_flux_normal_cuda(phifields_flop)
		 if(ldiagnostic)call end_timing2("LB","repulsive_flux")
#endif
#endif

#ifdef ASYNCMPI
#ifndef FUSED_MOMENTS
         if(ldiagnostic)call start_timing2("LB","moments_LB_cuda_ext")
         call moments_LB_cuda_ext(hfields_flip,hfields_flop &
#ifdef TWOCOMPONENT         
          ,phifields_flop &
#endif
          )
         if(ldiagnostic)call end_timing2("LB","moments_LB_cuda_ext")	
#endif
         
         if(ldiagnostic)call start_timing2("LB","ex_hfields_sendrecv")
		 call exchange_hfields_sendrecv(hfields_flop)
		 if(ldiagnostic)call end_timing2("LB","ex_hfields_sendrecv")
#ifndef FUSED_MOMENTS
		 if(ldiagnostic)call start_timing2("LB","ex_forces_sendrecv")
		 call exchange_forces_sendrecv
		 if(ldiagnostic)call end_timing2("LB","ex_forces_sendrecv")
#endif
		 if(ldiagnostic)call start_timing2("LB","ex_hfields_intpbc")
		 call exchange_hfields_intpbc(hfields_flop)
		 if(ldiagnostic)call end_timing2("LB","ex_hfields_intpbc")
#ifndef FUSED_MOMENTS
		 if(ldiagnostic)call start_timing2("LB","ex_forces_intpbc")
		 call exchange_forces_intpbc
		 if(ldiagnostic)call end_timing2("LB","ex_forces_intpbc")
#endif
		 
#ifndef FUSED_MOMENTS
         if(ldiagnostic)call start_timing2("LB","moments_LB_cuda_int")
         call moments_LB_cuda_int(hfields_flip,hfields_flop &
#ifdef TWOCOMPONENT         
          ,phifields_flop &
#endif
          )
         if(ldiagnostic)call end_timing2("LB","moments_LB_cuda_int")
#endif
         
         if(ldiagnostic)call start_timing2("LB","ex_hfields_wait")
		 call exchange_hfields_wait(hfields_flop)
		 if(ldiagnostic)call end_timing2("LB","ex_hfields_wait")
#ifndef FUSED_MOMENTS
		 if(ldiagnostic)call start_timing2("LB","ex_forces_wait")
		 call exchange_forces_wait
		 if(ldiagnostic)call end_timing2("LB","ex_forces_wait")
#endif
#else

#ifndef FUSED_MOMENTS
         if(ldiagnostic)call start_timing2("LB","moments_LB_cuda")
         call moments_LB_cuda(hfields_flip,hfields_flop &
#ifdef TWOCOMPONENT         
          ,phifields_flop &
#endif
          )
         if(ldiagnostic)call end_timing2("LB","moments_LB_cuda")	 
#endif
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
#ifndef FUSED_MOMENTS
         if(ldiagnostic)call start_timing2("LB","ex_forces_sendrecv")
		 call exchange_forces_sendrecv
		 if(ldiagnostic)call end_timing2("LB","ex_forces_sendrecv")
#endif
#ifndef FUSED_MOMENTS
		 if(ldiagnostic)call start_timing2("LB","ex_forces_intpbc")
		 call exchange_forces_intpbc
		 if(ldiagnostic)call end_timing2("LB","ex_forces_intpbc")
#endif
#ifndef FUSED_MOMENTS
		 if(ldiagnostic)call start_timing2("LB","ex_forces_wait")
		 call exchange_forces_wait
		 if(ldiagnostic)call end_timing2("LB","ex_forces_wait")  
#endif
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
#if defined(TAYLORGREEN) && !defined(TWOCOMPONENT)
               call print_taylorgreen(142)
#endif
#if defined(LAMBTEST) && defined(TWOCOMPONENT)
               call print_lamb(142)
#endif
#ifdef LAPLACE
               call print_laplace(142)
#endif
#ifdef CAPILLARYWAVE
               call print_capillary_wave(142,.false.)
#endif
               time_actual=current_time()
               gi=iprobe;gj=jprobe;gk=kprobe
               subchords(1)=(gi-1)/nx
               subchords(2)=(gj-1)/ny
               subchords(3)=(gk-1)/nz
#if defined(_NVML) && defined(GETPOWER)  
               num_p_w=num_p_w+1         
               p_w=real(p_mw,kind=db)*1.0e-3_db
               tot_energy=tot_energy+p_w*(time_actual-time_actual_old)
#endif 
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
         
!****************phi exchange: periodic boundary conditions on phi************************
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
      !call cpu_time(ts2)
      ts2=current_time()
#if defined(_NVML) && !defined(GETPOWER)      
      ierrc=get_gpu_energy_mJ_u64(mydev_c,energy_2)
#endif  
      !***********************************write restart************************
      ! The restart files are written from hfields_flop and phifields_flop. These hold the newest
      ! fields after a time-limit exit from the FLIP half; after a complete iteration the newest
      ! fields are in the flip buffers.
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
      ! step_energy: energy per time step summed over all ranks, from the power sampled once every
      ! stamp_term steps (GETPOWER) or from the energy counters read before and after the loop.
#ifdef _NVML 
#ifdef GETPOWER
      call sum_world_float(tot_energy)
      step_energy=tot_energy/real(num_p_w*stamp_term,kind=db)
#else
      tot_energy=u64_delta_J(energy_1,energy_2)
      call sum_world_float(tot_energy)
      step_energy=tot_energy/real(nsteps,kind=db)
#endif
#endif 
#if defined(TAYLORGREEN) || defined(LAMBTEST) || defined(LAPLACE) || defined(CAPILLARYWAVE)
      call close_benchmark(142) 
#endif
      !$wait
      !$acc end data


   endsubroutine

endmodule
