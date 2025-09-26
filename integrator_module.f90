#include "defines.h"
module integrator_module
!   use mpi_template
#ifdef _OPENACC
   use openacc
#endif
   use mpi_template, only : nbuff,coords,myoffset,myrank,nprocs,intpbc_dir, &
           num_links_pops,links_pops,datampi,f_datampi,fvec_datampi,b_datampi,i_datampi, &
           send_extr,recv_extr,f_send_extr,f_recv_extr,fvec_send_extr,fvec_recv_extr, &
           b_send_extr,b_recv_extr,send_buffmpi,recv_buffmpi,nbuffmpi_send,nbuffmpi_recv, &
           f_send_buffmpi,f_recv_buffmpi,f_nbuffmpi_send,f_nbuffmpi_recv,lbuff, &
           fvec_send_buffmpi,fvec_recv_buffmpi,fvec_nbuffmpi_send,fvec_nbuffmpi_recv, &
           b_send_buffmpi,b_recv_buffmpi,b_nbuffmpi_send,b_nbuffmpi_recv, &
           exchange_float_sendrecv,exchange_float_intpbc,exchange_float_wait, &
           exchange_floatvec_sendrecv,exchange_floatvec_intpbc,exchange_floatvec_wait, &
#ifdef REPULSIVE_FLUX 
           exchange_bvec_intpbc,exchange_bvec_wait,exchange_bvec_sendrecv, &
#endif
#ifdef CRAY 
           int_buffpbc,nbuffpbc_int, &
#endif
           intpbcsend_extr,intpbcrecv_extr,skip_myoffset,or_world_l, &
           exchange_pops_sendrecv,exchange_pops_intpbc,exchange_pops_wait,dostop
           
           
   use prints, only : get_memory_gpu,print_memory_registration_gpu, &
    driver_print_raw_sync,driver_print_raw_sync, &
    init_output,driver_print_vtk_sync,read_restart_2c,read_restart_1c, &
#ifdef _OPENACC
    printdeviceproperties, &
#endif
    write_restart_1c,write_restart_2c,driver_print_raw_sync2d, &
    driver_print_raw_isfluid
   use vars
   use bcs3D, only : bcs_mesoscopic_all
   use lb_kernels, only :  &
#ifdef TWOCOMPONENT  
    compute_norm_interface,compute_laplacian_phi, &
#endif
#ifdef REPULSIVE_FLUX
    thinfilm_scan_mark,repulsive_flux_tangential,repulsive_flux_normal, &
#endif
    moments_lb,compute_densityratio,fused_lb
#if defined(_KERNELCUDA) && defined(_OPENACC)        
   use lb_cuda_kernels, only : moments_LB_cuda,fused_lb_cuda
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

      step=0

      flip=mod(step,2)+1     
      flop = 3 - flip
      
      !$acc data copy(step,lx,ly,lz,nx,ny,nz,coords,myoffset,f,isfluid,myrank, &
      !$acc& pxx,pyy,pzz,pxy,pxz,pyz,rho,u,v,w,rhoprint,velprint,radius, &
	  !$acc& tau1,visc1,rho_r,rho_b,invrho_r,invrho_b,omega,lap_phi, &
      !$acc& intpbc_dir,num_links_pops,links_pops,datampi,f_datampi,uwall,udotc,uu, &
      !$acc& send_extr,recv_extr,f_send_extr,f_recv_extr,fux,fvy,fwz, &
#ifdef TWOCOMPONENT
      !$acc& sharp_c,beta,tau2,visc2, kapp, sigma,width,tau_diff,corr,global_phi_sum_ini,global_count,global_phi_sum,&
	  !$acc& arr_x, arr_y, arr_z,normx,normy,normz,modgrad,selphi,global_phi_change, &
#ifdef DENSRATIO
	  !$acc& rhophi, &
#endif
#ifdef MONOD
	  !$acc&  mu_max,Ks, &
#endif
#ifdef IMPOSED_PRESSURE_GRADIENT
      !$acc& rhoIN,rhoOUT, &
#endif
#ifdef REPULSIVE_FLUX
	 !$acc& Jx,Jy,Jz,pair_i,pair_j,pair_k,rep_mask,q_th,cosOppT,pwr,A_rep,win, &
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
#if defined REPULSIVE_FLUX
      !$acc& b_send_extr,b_recv_extr, &
#endif
#ifdef CRAY 
      !$acc& int_buffpbc,nbuffpbc_int, &
#endif
      !$acc& intpbcsend_extr,intpbcrecv_extr, &
      !$acc& fvec_send_extr,fvec_recv_extr,fx,fy,fz) &
      !$acc& create(send_buffmpi,recv_buffmpi,f_send_buffmpi, &
      !$acc& f_recv_buffmpi,fvec_send_buffmpi &
#if defined REPULSIVE_FLUX
      !$acc& ,b_send_buffmpi,b_recv_buffmpi &
#endif
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
         call read_restart_2c(iframe,iframe2D)
      endif
      if(ldiagnostic)call start_timing2("LB","moments_phi")
      call compute_densityratio
      if(ldiagnostic)call end_timing2("LB","moments_phi")
#else
      if(lrestart)then
         call read_restart_1c(iframe,iframe2D)
      endif
#endif      

      if(lprint)then
	 
#ifdef ACCNOKERNELS
#if defined(DENSRATIO) && defined(TWOCOMPONENT)
#ifdef WRITEPRESS
         !$acc parallel loop independent collapse(3) present(rhoprint,velprint,pressprint,rhophi,u,v,w)
#else
         !$acc parallel loop independent collapse(3) present(rhoprint,velprint,rhophi,u,v,w)
#endif
#endif

#if defined(TWOCOMPONENT) && !defined(DENSRATIO)
#ifdef WRITEPRESS
		 !$acc parallel loop independent collapse(3) present(rhoprint,velprint,,pressprint,selphi,rho,u,v,w)
#else
         !$acc parallel loop independent collapse(3) present(rhoprint,velprint,selphi,u,v,w)
#endif 
#endif                 
#ifndef TWOCOMPONENT 
         !$acc parallel loop independent collapse(3) present(rhoprint,velprint,rho,u,v,w)
#endif         
#else
#if defined(DENSRATIO) && defined(TWOCOMPONENT)
#ifdef WRITEPRESS
         !$acc kernels present(rhoprint,velprint,pressprint,rhophi,u,v,w)
#else
         !$acc kernels present(rhoprint,velprint,rhophi,u,v,w)
#endif
#endif
#if defined(TWOCOMPONENT) && !defined(DENSRATIO)
#ifdef WRITEPRESS
		 !$acc kernels present(rhoprint,velprint,pressprint,rho,selphi,u,v,w)
#else
         !$acc kernels present(rhoprint,velprint,selphi,u,v,w)
#endif
#endif
#ifndef TWOCOMPONENT
         !$acc kernels present(rhoprint,velprint,rho,u,v,w)
#endif
         !$acc loop independent collapse(3)  private(i,j,k)
#endif
         do k=1,nzskip
            do j=1,nyskip
               do i=1,nxskip
#if defined(DENSRATIO) && defined(TWOCOMPONENT)
                  rhoprint(i,j,k)=real(rhophi(i*stepskip,j*stepskip,k*stepskip),kind=printdb)
#ifdef WRITEPRESS                   
                  pressprint(i,j,k)=real(rho(i*stepskip,j*stepskip,k*stepskip),kind=printdb)
#endif
#endif

#if defined(TWOCOMPONENT) && !defined(DENSRATIO)
                  rhoprint(i,j,k)=real(selphi(i*stepskip,j*stepskip,k*stepskip,flip),kind=printdb)
#ifdef WRITEPRESS                  
                  pressprint(i,j,k)=real(rho(i*stepskip,j*stepskip,k*stepskip),kind=printdb)				  
#endif
#endif
#ifndef TWOCOMPONENT
				  rhoprint(i,j,k)=real(rho(i*stepskip,j*stepskip,k*stepskip),kind=printdb)
#endif
					velprint(1,i,j,k)=real(u(i*stepskip,j*stepskip,k*stepskip),kind=printdb)
					velprint(2,i,j,k)=real(v(i*stepskip,j*stepskip,k*stepskip),kind=printdb)
					velprint(3,i,j,k)=real(w(i*stepskip,j*stepskip,k*stepskip),kind=printdb)
               enddo
            enddo
         enddo
#ifdef ACCNOKERNELS
         !$acc end parallel loop
#else
         !$acc end kernels
#endif
#if defined(WRITEPRESS)
         !$acc update host(rhoprint,velprint,pressprint)
#else
         !$acc update host(rhoprint,velprint)
#endif
         !$acc wait
         if(lvtk)then
            call driver_print_vtk_sync(iframe)
         endif
         if(lraw)then
            call driver_print_raw_sync(iframe)
            !2d planes print
            call driver_print_raw_sync2D(iframe2D,nplanes,ndir,npoint)
         endif
      endif



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
  

      !***********************************scambio moments************************
#ifdef TWOCOMPONENT
	  if(ldiagnostic)call start_timing2("LB","exchange_float_sendrecv")         
          if(lbuff)call exchange_bvec_sendrecv
	  call exchange_float_sendrecv
	  if(ldiagnostic)call end_timing2("LB","exchange_float_sendrecv")

	  if(ldiagnostic)call start_timing2("LB","exchange_float_intpbc")
 	  if(lbuff)call exchange_bvec_intpbc
	  call exchange_float_intpbc
	  if(ldiagnostic)call end_timing2("LB","exchange_float_intpbc")
	  
	  if(ldiagnostic)call start_timing2("LB","exchange_float_wait")
          if(lbuff)call exchange_bvec_wait
	  call exchange_float_wait
	  if(ldiagnostic)call end_timing2("LB","exchange_float_wait")
      !***********************************ora che ho phi in cornice calcolo normx normy normz************************
      if(ldiagnostic)call start_timing2("LB","compute_norm")
      call compute_norm_interface
      if(ldiagnostic)call end_timing2("LB","compute_norm")
      !***********************************scambio normx normy normz************************
      if(ldiagnostic)call start_timing2("LB","exchange_floatvec_sendrecv")
      call exchange_floatvec_sendrecv
	  if(ldiagnostic)call end_timing2("LB","exchange_floatvec_sendrecv")
	  if(ldiagnostic)call start_timing2("LB","exchange_floatvec_intpbc")
      call exchange_floatvec_intpbc
	  if(ldiagnostic)call end_timing2("LB","exchange_floatvec_intpbc")
	  
      if(ldiagnostic)call start_timing2("LB","exchange_floatvec_wait")
	  call exchange_floatvec_wait
      if(ldiagnostic)call end_timing2("LB","exchange_floatvec_wait")
      !***********************************compute laplacian************************
      if(ldiagnostic)call start_timing2("LB","force_2c")
      call compute_laplacian_phi
      if(ldiagnostic)call end_timing2("LB","force_2c")
#endif	  
      !***********************************compute moments***********************
	  if(ldiagnostic)call start_timing2("LB","moments")

#if defined(_KERNELCUDA) && defined(_OPENACC)
      call moments_LB_cuda
#else
      call moments_LB
#endif
      if(ldiagnostic)call end_timing2("LB","moments")

      call get_memory_gpu(mymemory,totmemory)
      call print_memory_registration_gpu(6,'DEVICE memory occupied at the start', &
         'total DEVICE memory',mymemory,totmemory)

      !write(6,'(a,i8,1x,i8,1x,a)')' sono arrivato ',myrank,__LINE__,__FILE__
      !call dostop
      time_init=current_time()
      time_actual_old=time_init
      call cpu_time(ts1)
      do step=1,nsteps
         !***********************************Print on files 3D************************

         flip=mod(step,2)+1     
         flop = 3 - flip
         !$acc wait
         !$acc update device(step,flip,flop)
         !$acc wait
        
         if(lprint)then
            
            if(mod(step,stamp).eq.0 .or. mod(step,stamp2D).eq.0 .or. mod(step,stamp_term).eq.0)then
               if(ldiagnostic)call start_timing2("IO","print")
               !write(6,*)'vorrei stampare'
#ifdef ACCNOKERNELS
#if defined(DENSRATIO) && defined(TWOCOMPONENT)
#ifdef WRITEPRESS
               !$acc parallel loop independent collapse(3) present(rhoprint,velprint,pressprint,rhophi,u,v,w)
#else
               !$acc parallel loop independent collapse(3) present(rhoprint,velprint,rhophi,u,v,w)
#endif
#endif
#if defined(TWOCOMPONENT) && !defined(DENSRATIO)
#ifdef WRITEPRESS
			   !$acc parallel loop independent collapse(3) present(rhoprint,velprint,pressprint,rho,selphi,u,v,w)
#else
               !$acc parallel loop independent collapse(3) present(rhoprint,velprint,selphi,u,v,w)
#endif
#endif
#ifndef TWOCOMPONENT
                !$acc parallel loop independent collapse(3) present(rhoprint,velprint,rho,u,v,w)
#endif
#else
#if defined(DENSRATIO) && defined(TWOCOMPONENT)
#ifdef WRITEPRESS
               !$acc kernels present(rhoprint,velprint,pressprint,rhophi,u,v,w)
#else
               !$acc kernels present(rhoprint,velprint,rhophi,u,v,w)
#endif
#endif
#if defined(TWOCOMPONENT) && !defined(DENSRATIO)
#ifdef WRITEPRESS
		       !$acc kernels present(rhoprint,velprint,pressprint,rho,selphi,u,v,w)
#else
			   !$acc kernels present(rhoprint,velprint,selphi,u,v,w)
#endif
#endif
#ifndef TWOCOMPONENT
			   !$acc kernels present(rhoprint,velprint,rho,u,v,w)
#endif
               !$acc loop independent collapse(3)  private(i,j,k)
#endif
               do k=1,nzskip
                  do j=1,nyskip
                     do i=1,nxskip
#if defined(DENSRATIO) && defined(TWOCOMPONENT)
						rhoprint(i,j,k)=real(rhophi(i*stepskip,j*stepskip,k*stepskip),kind=printdb)
#ifdef WRITEPRESS                  
                        pressprint(i,j,k)=real(rho(i*stepskip,j*stepskip,k*stepskip),kind=printdb)
#endif
#endif

#if defined(TWOCOMPONENT) && !defined(DENSRATIO)
                        rhoprint(i,j,k)=real(selphi(i*stepskip,j*stepskip,k*stepskip,flip),kind=printdb)
#ifdef WRITEPRESS                  
                        pressprint(i,j,k)=real(rho(i*stepskip,j*stepskip,k*stepskip),kind=printdb)
#endif
#endif
#ifndef TWOCOMPONENT
						rhoprint(i,j,k)=real(rho(i*stepskip,j*stepskip,k*stepskip),kind=printdb)
#endif
                        velprint(1,i,j,k)=real(u(i*stepskip,j*stepskip,k*stepskip),kind=printdb)
						velprint(2,i,j,k)=real(v(i*stepskip,j*stepskip,k*stepskip),kind=printdb)
						velprint(3,i,j,k)=real(w(i*stepskip,j*stepskip,k*stepskip),kind=printdb)
                     enddo
                  enddo
               enddo
#ifdef ACCNOKERNELS
               !$acc end parallel loop
#else
               !$acc end kernels
#endif
#if defined(WRITEPRESS)
               !$acc update host(rhoprint,velprint,pressprint)
#else
               !$acc update host(rhoprint,velprint)
#endif
               !$acc wait
            !endif

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
               gi=18;gj=18;gk=2
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
                       '; probe values : ',isfluid(i,j,k),rhoprint(i,j,k),velprint(1:3,i,j,k), &
#if defined(INTERNAL_OBSTACLES)                       
                       corr,dphi,global_phi_sum_ini,global_phi_sum_new-global_phi_change_new,global_phi_sum_new !real(global_count_new)
#else
                       0.0e0                       
#endif                       
                 call flush(6)
               endif
               time_actual_old=time_actual
             endif
             if(ldiagnostic)call end_timing2("IO","print")
           endif
         endif
         
         !***********************************collision + no slip + forcing: fused implementation*********
		 if(ldiagnostic)call start_timing2("LB","fused")

#if defined(_KERNELCUDA) && defined(_OPENACC)      
         call fused_LB_cuda
#else
         call fused_LB
#endif
         if(ldiagnostic)call end_timing2("LB","fused")
#ifdef TWOCOMPONENT	 
!****************scambio phi: boundary condition periodiche su phi************************
		 if(ldiagnostic)call start_timing2("LB","exchange_float_sendrecv")
                 if(lbuff)call exchange_bvec_sendrecv
		 call exchange_float_sendrecv
		 if(ldiagnostic)call end_timing2("LB","exchange_float_sendrecv")
		 
         if(ldiagnostic)call start_timing2("LB","exchange_float_intpbc")
                 if(lbuff)call exchange_bvec_intpbc
		 call exchange_float_intpbc
		 if(ldiagnostic)call end_timing2("LB","exchange_float_intpbc")
		 
         if(ldiagnostic)call start_timing2("LB","exchange_float_wait")
                 if(lbuff)call exchange_bvec_wait
		 call exchange_float_wait
		 if(ldiagnostic)call end_timing2("LB","exchange_float_wait")
         !***********************************ora che ho phi in cornice calcolo normx normy normz************************
         if(ldiagnostic)call start_timing2("LB","compute_norm")
         call compute_norm_interface
         if(ldiagnostic)call end_timing2("LB","compute_norm")
         !***********************************scambio normx normy normz************************
         if(ldiagnostic)call start_timing2("LB","exchange_floatvec_sendrecv")
		 call exchange_floatvec_sendrecv
		 if(ldiagnostic)call end_timing2("LB","exchange_floatvec_sendrecv")
		 if(ldiagnostic)call start_timing2("LB","exchange_floatvec_intpbc")
		 call exchange_floatvec_intpbc
		 if(ldiagnostic)call end_timing2("LB","exchange_floatvec_intpbc")
		 if(ldiagnostic)call start_timing2("LB","exchange_floatvec_wait")
		 call exchange_floatvec_wait
		 if(ldiagnostic)call end_timing2("LB","exchange_floatvec_wait")
         
         if(ldiagnostic)call start_timing2("LB","force_2c")
         call compute_laplacian_phi
         if(ldiagnostic)call end_timing2("LB","force_2c")
#endif
 
		 !***********************************pbcs boundary conditions ********************************!
		 !call pbcs
         if(ldiagnostic)call start_timing2("LB","pbcs")
         call exchange_pops_sendrecv
         call exchange_pops_intpbc
         call exchange_pops_wait
         if(ldiagnostic)call end_timing2("LB","pbcs")
         !************ thread-safe boundary condition setup
         if(ldiagnostic)call start_timing2("LB","bcs_TSLB")
         call bcs_mesoscopic_all  
         if(ldiagnostic)call end_timing2("LB","bcs_TSLB") 
         !***********************************moments************************
#ifdef DENSRATIO
         if(ldiagnostic)call start_timing2("LB","moments_phi")
         call compute_densityratio
         if(ldiagnostic)call end_timing2("LB","moments_phi")
#endif
         if(ldiagnostic)call start_timing2("LB","moments")


 #ifdef REPULSIVE_FLUX
		  call thinfilm_scan_mark
		  call repulsive_flux_normal
 #endif


#if defined(_KERNELCUDA) && defined(_OPENACC) 
         call moments_LB_cuda
#else
         call moments_LB
#endif
         if(ldiagnostic)call end_timing2("LB","moments")
         
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
      !$acc update host(rho,u,v,w,pxx,pxy,pxz,pyy,pyz,pzz,selphi)
      !$wait
      call write_restart_2c(iframe,iframe2D)
#else      
      !$acc update host(rho,u,v,w,pxx,pxy,pxz,pyy,pyz,pzz)
      !$wait
      call write_restart_1c(iframe,iframe2D)
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
