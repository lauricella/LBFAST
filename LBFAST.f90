#include "defines.h"
program threadsafeLB
   use mpi_template, only : myrank,proc_x,proc_y,proc_z,pbc_x,pbc_y,pbc_z, &
    coords,start_mpi,setup_mpi,exchange_isf_sendrecv,exchange_isf_intpbc, &
    exchange_isf_wait,dostop,print_version_code,nprocs
#ifdef _OPENACC
   use openacc
#endif
#if defined(_OPENACC)            
   use lb_cuda_driver, only : setup_cuda
#endif
   use prints, only : intstr,copystring,driver_read_isfluid_raw, &
    driver_read_init_raw
   use vars
   use allocate_arrays, only : allocate_struct
   use integrator_module, only : integrator
   use profiling_m, only : ldiagnostic,idiagnostic,itime_start,itime_counter, &
    set_value_ldiagnostic,set_value_idiagnostic,get_memory,print_memory_registration, &
    timer_init,startpreprocessingtime,printsimulationtime,print_timing_final, &
    get_totram
   
   use initial_condts, only: initial_conditions_all
   use stat_module, only: compare_benchmark
   implicit none
   
   logical :: lexist
   integer, parameter :: inputio=24
   real(kind=db) :: dist,tot_time,Glups
#ifdef _NVML   
   real(kind=db) :: Pavg=ZERO,eff_glups_per_w,E_tot,J_per_LUP
#endif
   namelist /simulation/ nsteps,stamp,stamp2D,stepskip,lreadisfluid, &
    lreadinit,lx,ly,lz,lprint,lvtk,lraw,lrestart,pbc_x,pbc_y,pbc_z,openbc,lasync, &
    iprobe,jprobe,kprobe,openbc_type_x,openbc_type_y,openbc_type_z, &
    openbc_press_x,openbc_press_y,openbc_press_z, &
    openbc_u_x,openbc_u_y,openbc_u_z,openbc_v_x,openbc_v_y,openbc_v_z, &
    openbc_w_x,openbc_w_y,openbc_w_z,&
    nplanes,stamp_term,time_limit,every_time_check,lwriterestart,lweakscaling
    
   namelist /fluid/ fx,fy,fz,visc1, &
    radius,width,center, &
#ifdef TWOCOMPONENT
    tau_diff,sigma,wettab_r,wettab_b,visc2, &
#ifdef MONOD	
	mu_max,Ks, &
#endif
#ifdef DENSRATIO
    rho_r,rho_b, &
#endif
#endif    
#ifdef ELASTIC_FORCE
    k_elastic,lambda_rel, &
#endif 
#ifdef IMPOSED_PRESSURE_GRADIENT
    rhoIN,rhoOUT,&
#endif 
#ifdef REPULSIVE_FLUX
	A_rep,cosOppT,q_th &
#endif
    uwall,npoint,ndir


#ifdef _OPENACC
   devType = acc_get_device_type()
   devNum=acc_get_device_num(devType)
#endif

#ifdef _OPENACC
   ngpus=acc_get_num_devices(devType)
#else
   ngpus=0
#endif

   !!!!!!! START MPI!!!!!!!!!
   call start_mpi

   proc_x=1
   proc_y=1
   proc_z=1

#ifdef MPI
   !leggi decomposizione da riga di comando
   narg = command_argument_count()
   if (narg /= 3 .and. (.not. lreadinput)) then
      if(myrank==0)then
        write(6,*) 'error!'
        write(6,*) 'the command line should be'
        write(6,*) '[executable] [proc_x] [proc_y] [proc_z]'
        write(6,*) 'proc_x = decomposition along x'
        write(6,*) 'proc_y = decomposition along y'
        write(6,*) 'proc_z = decomposition along z'
        write(6,*) 'STOP!'
      endif
      call dostop
   endif
   
   if (narg /= 4 .and. lreadinput) then
      if(myrank==0)then 
        write(6,*) 'error!'
        write(6,*) 'the command line should be'
        write(6,*) '[executable] [proc_x] [proc_y] [proc_z] [input]'
        write(6,*) 'proc_x = decomposition along x'
        write(6,*) 'proc_y = decomposition along y'
        write(6,*) 'proc_z = decomposition along z'
        write(6,*) 'input  = name of input file'
        write(6,*) 'STOP!'
      endif
      call dostop
   endif

   do i = 1, narg
      arg=repeat(' ',mxln)
      call getarg(i, arg)
      if(i==1)then
         call copystring(arg,directive,mxln)
         proc_x=intstr(directive,mxln,inumchar)
         if(myrank==0)write(6,'(a,i8)') 'proc_x  = ',proc_x
      elseif(i==2)then
         call copystring(arg,directive,mxln)
         proc_y=intstr(directive,mxln,inumchar)
         if(myrank==0)write(6,'(a,i8)') 'proc_y  = ',proc_y
      elseif(i==3)then
         call copystring(arg,directive,mxln)
         proc_z=intstr(directive,mxln,inumchar)
         if(myrank==0)write(6,'(a,i8)') 'proc_z  = ',proc_z
      elseif(i==4 .and. lreadinput)then
         inipFile=repeat(' ',mxln)
         inipFile=trim(arg)
         if(myrank==0)write(6,'(2a)') 'file  = ',trim(inipFile)
      endif
   enddo
#else
   !leggi decomposizione da riga di comando
   narg = command_argument_count()
   if(lreadinput)then
     if(narg /= 1)then
        if(myrank==0)then 
          write(6,*) 'error!'
          write(6,*) 'the command line should be'
          write(6,*) '[executable] [input]'
          write(6,*) 'input  = name of input file'
          write(6,*) 'STOP!'
        endif
        call dostop
     endif
     
     do i = 1, narg
        arg=repeat(' ',mxln)
        call getarg(i, arg)
        if(i==1)then
          inipFile=repeat(' ',mxln)
          inipFile=trim(arg)
          if(myrank==0)write(6,*) 'file  = ',trim(inipFile)
        endif
     enddo
   endif
#endif

#ifdef DOBENCHMARK   
   nsteps=100
   stamp=100
   stamp2D=10
   stepskip=1
   lreadisfluid=.false.
   lreadinit=.false.
#ifndef MYSIDE
#define MYSIDE 512
#endif
   lx=MYSIDE
   ly=MYSIDE
   lz=MYSIDE
   if(lweakscaling)then
     !lz=MYSIDE*proc_x*proc_y*proc_z !weakscaling 1D
     lx = MYSIDE * proc_x
     ly = MYSIDE * proc_y
     lz = MYSIDE * proc_z
   endif
   radius=lx/4
   uwall=-0.0e0
   lprint=.true.
   lvtk=.false.
   lraw=.true.
   lasync=.false.
   lrestart=.false. !metti true se riparti da restart
   pbc_x=1  !(0=false 1=true)
   pbc_y=1
   pbc_z=1
   fx=0.0e0
   fy=0.0e0
   fz=0.0e0
#endif

   if(lreadinput)then
    
     inquire(file=trim(inipFile),exist=lexist)
      if(.not. lexist)then
        if(myrank==0)then
          write(6,*)'ERROR: file ',trim(inipFile),' does not exist!'
          write(6,*) 'STOP!'
        endif
        call dostop
      endif
      open(unit=inputio,file=trim(inipFile),status='old')
      read(inputio,nml=simulation)
      
      if(allocated(ndir))deallocate(ndir)
      if(allocated(npoint))deallocate(npoint)
      allocate(ndir(nplanes),npoint(nplanes))
      
      read(inputio,nml=fluid)
      close(inputio)
      tau1=visc1*invcssq + 0.5_db
      omega=1.0_db/tau1      
#ifdef TWOCOMPONENT      
      beta=12.0_db*sigma/width 
      kapp=1.5_db*sigma*width 
      sharp_c=4.0_db*tau_diff/width
      !visc2= 1.0_db*visc1 
      tau2=visc2*invcssq + 0.5_db
#ifdef REPULSIVE_FLUX
      win=3 
	  !q_th=0.125
      !cosOppT=-0.8	   
#endif 
#endif    
      if(lweakscaling)then
        !lz=lz*proc_x*proc_y*proc_z
        lx = MYSIDE * proc_x
        ly = MYSIDE * proc_y
        lz = MYSIDE * proc_z
      endif  
   endif
     

   ! start diagnostic if requested
   mydiagnostic=.true.
   tdiagnostic=1
   call set_value_ldiagnostic(mydiagnostic)
   call set_value_idiagnostic(tdiagnostic)
   if(ldiagnostic)then
      call timer_init()
      call startPreprocessingTime()
   endif
      
   !!!!!!!SETUP MPI!!!!!!!!!!!!!!!!!!!
   call setup_mpi()
   
#if defined(_OPENACC)
   call setup_cuda
#endif

   call allocate_struct
    
   !ex=(/0, 1, -1, 0,  0,  0,  0,  1,  -1,  1,  -1,  0,   0,  0,   0,  1,  -1,  -1,   1/)
   !ey=(/0, 0,  0, 1, -1,  0,  0,  1,  -1, -1,   1,  1,  -1,  1,  -1,  0,   0,   0,   0/)
   !ez=(/0, 0,  0, 0,  0,  1, -1,  0,   0,  0,   0,  1,  -1, -1,   1,  1,  -1,   1,  -1/)
   !*****************************geometry************************************************
   isfluid=1
   
   if(lreadisfluid)then
     if(myrank==0)write(6,'(a)') 'Reading initialization isfluid file.....'
     call driver_read_isfluid_raw(0)
     if(myrank==0)write(6,'(a)') 'Completed!'
   endif
   
  
#if defined(INTERNAL_OBSTACLES)
	   do k=1,nz
	      gk=nz*coords(3)+k
	      do j=1,ny
	         gj=ny*coords(2)+j
	         do i=1,nx
	            gi=nx*coords(1)+i
	             !*****some isfluid def
	             if(pbc_x==0)then
	               if(gi==1)isfluid(i,j,k)=0
	               if(gi==lx)isfluid(i,j,k)=0
	             endif
	             if(pbc_y==0)then
	               if(gj==1)isfluid(i,j,k)=0
	               if(gj==ly)isfluid(i,j,k)=0
	             endif
	             if(pbc_z==0)then
	               if(gk==1)isfluid(i,j,k)=0
	               if(gk==lz)isfluid(i,j,k)=0
	             endif
#define noKARMANN
#ifdef KARMANN
	             !dist=sqrt((float(gi)-lx/TWO)**TWO + (float(gj)-ly/TWO)**TWO+(float(gk)-(lz/TWO))**TWO)
	             dist=sqrt((float(gi)-lx/FOUR)**TWO + (float(gj)-ly/TWO)**TWO)
	             if(int(dist)<=8)isfluid(i,j,k)=0
#endif
	         enddo
	      enddo
	   enddo
#endif   
   
#ifdef DOBENCHMARK 
   isfluid=1
#endif

   !setup domain decomposition among MPI process
   !******************************(only once at the beginning)************************
   call exchange_isf_sendrecv
   call exchange_isf_intpbc
   call exchange_isf_wait

!  identify fluid nodes close to the boundary
   do k=1,nz
      gk=nz*coords(3)+k
      do j=1,ny
         gj=ny*coords(2)+j
         do i=1,nx
            gi=nx*coords(1)+i
            if(isfluid(i,j,k).eq.1)then
              do l=1,nlinks
                !if the fluid node is close to a boundary put -1
			    if(abs(isfluid(i+ex(l),j+ey(l),k+ez(l))).ne.1)then
			      isfluid(i,j,k)=-1
			    endif
			  enddo 
            endif
         enddo
      enddo
    enddo
   call exchange_isf_sendrecv
   call exchange_isf_intpbc
   call exchange_isf_wait
   
   
   !*************************************initial conditions ************************
   step=0

   flip=mod(step,2)+1     
   flop = 3 - flip
   call initial_conditions_all
   if(lreadinit)then
     if(myrank==0)write(6,'(a)') 'Reading initialization files.....'
     call driver_read_init_raw(0)
     if(myrank==0)write(6,'(a)') 'Completed!'
   endif
   
   !*************************************check data ************************
   if(myrank==0)then
      call print_version_code
#if defined(_OPENACC)
      write(6,'(a)') 'GPU VERSION COMPILED'
#endif
      write(6,'(a)') '*******************LB type*****************'
#if LATTICE == 27
#ifdef HIGHORDER
      write(6,'(a)') 'Using D3Q27 lattice with HIGH-ORDER'
#else
      write(6,'(a)') 'Using D3Q27 lattice'
#endif
#elif LATTICE == 19
      write(6,'(a)') 'Using D3Q19 lattice'
#elif LATTICE == 15
      write(6,'(a)') 'Using D3Q15 lattice'
#endif
#ifdef TWOCOMPONENT
      write(6,'(a)') 'Compiled for two components'
#else
      write(6,'(a)') 'Compiled for single component'
#endif
#ifdef MIXEDPRC
      write(6,'(a)') 'Compiled with mixed precision'
#endif
#if PRC==4
      write(6,'(a)') 'Compiled in single precision'
#elif PRC==8
      write(6,'(a)') 'Compiled in double precision'
#endif
#ifdef MIXEDPRC
#if STRPRC==2
      write(6,'(a)') 'Compiled in half precision for storaing'
#elif STRPRC==4
      write(6,'(a)') 'Compiled in single precision for storaing'
#elif STRPRC==8
      write(6,'(a)') 'Compiled in double precision for storaing'
#endif
#endif
#ifdef PRINTHALF
      write(6,'(a)') 'File written in half precision'
#else
      write(6,'(a)') 'File written in single precision'
#endif
      write(6,'(a)') '*******************LB data*****************'
      write(6,'(a,g16.8)') 'tau1',tau1
	  write(6,'(a,g16.8)') 'visc',visc1
	  write(6,'(a,g16.8)') 'omega',omega
#ifdef TWOCOMPONENT
      write(6,'(a,g16.8)') 'tau2',tau2
	  write(6,'(a,g16.8)') 'visc',visc2
	  write(6,'(a,g16.8)') 'diff',tau_diff
	  write(6,'(a,g16.8)') 'sharp_c',sharp_c
	  write(6,'(a,g16.8)') 'beta',beta
	  write(6,'(a,g16.8)') 'kapp',kapp
      write(6,'(a,g16.8)') 'sigma',sigma
#endif
      write(6,'(a,g16.8)') 'radius',radius
      write(6,'(a,g16.8)') 'width',width
      write(6,'(a,3g16.8)')'center',center
#ifdef REPULSIVE_FLUX
	  write(6,'(a,g16.8)') 'A_Rep', A_rep
	  write(6,'(a,g16.8)') 'qth', q_th
	  write(6,'(a,g16.8)') 'cosoppT', cosOppT
#endif
      write(6,'(a,g16.8)') 'fx',fx
      write(6,'(a,g16.8)') 'fy',fy
      write(6,'(a,g16.8)') 'fz',fz
      write(6,'(a,g16.8)') 'cssq',cssq
      
      write(6,'(a,g16.8)') 'rho_r',rho_r
      write(6,'(a,g16.8)') 'rho_b',rho_b

      write(6,'(a)') '*******************INPUT data*****************'
      write(6,'(a,i8)') 'lx',lx
      write(6,'(a,i8)') 'ly',ly
      write(6,'(a,i8)') 'lz',lz
      write(6,'(a,3i4)') 'pbc',pbc_x,pbc_y,pbc_z
      write(6,'(a,i4)') 'openbc',openbc
      write(6,'(a,2i4)') 'openbc_type_x',openbc_type_x
      write(6,'(a,2i4)') 'openbc_type_y',openbc_type_y
      write(6,'(a,2i4)') 'openbc_type_z',openbc_type_z
      write(6,'(a,2g16.8)') 'openbc_press_x',openbc_press_x
      write(6,'(a,2g16.8)') 'openbc_press_y',openbc_press_y
      write(6,'(a,2g16.8)') 'openbc_press_z',openbc_press_z
      write(6,'(a,2g16.8)') 'openbc_u_x',openbc_u_x
      write(6,'(a,2g16.8)') 'openbc_u_y',openbc_u_y
      write(6,'(a,2g16.8)') 'openbc_u_z',openbc_u_z
      write(6,'(a,2g16.8)') 'openbc_v_x',openbc_v_x
      write(6,'(a,2g16.8)') 'openbc_v_y',openbc_v_y
      write(6,'(a,2g16.8)') 'openbc_v_z',openbc_v_z
      write(6,'(a,2g16.8)') 'openbc_w_x',openbc_w_x
      write(6,'(a,2g16.8)') 'openbc_w_y',openbc_w_y
      write(6,'(a,2g16.8)') 'openbc_w_z',openbc_w_z
      write(6,'(a,3i4)') 'probe',iprobe,jprobe,kprobe
      write(6,'(a,l8)') 'lprint',lprint
      write(6,'(a,l8)') 'lvtk',lvtk
      write(6,'(a,l8)') 'lraw',lraw
      write(6,'(a,l8)') 'lrestart',lrestart
      write(6,'(a,l8)') 'lreadisfluid',lreadisfluid
      write(6,'(a,l8)') 'lreadinit',lreadinit
      write(6,'(a,l8)') 'lreadinput',lreadinput
      write(6,'(a,l8)') 'lwriterestart',lwriterestart
#ifdef IMPOSED_PRESSURE_GRADIENT
	  write(6,'(a,g16.8)') 'rhoIN', rhoIN
	  write(6,'(a,g16.8)') 'rhoIN', rhoOUT
#endif
      if(lweakscaling)write(6,'(a,l8)') 'lweakscaling',lweakscaling
      write(6,'(a,i16)') 'nsteps',nsteps
      if(time_limit>ZERO)then
      write(6,'(a,g16.8)') 'time_limit',time_limit
      write(6,'(a,i16)') 'every_time_check', every_time_check
      endif
      write(6,'(a,i16)') 'stamp',stamp
      write(6,'(a,i16)') 'stamp_term',stamp_term
      write(6,'(a,i16)') 'stamp2D',stamp2D
      write(6,'(a,i8)') 'nplanes',nplanes
      do i=1,nplanes
        write(6,'(a,3i8)') '    nplanes : ',i,ndir(i),npoint(i)
      enddo
#if defined(_OPENACC)
      write(6,'(a,i8)') 'TILE_DIMx',TILE_DIMx
      write(6,'(a,i8)') 'TILE_DIMy',TILE_DIMy
      write(6,'(a,i8)') 'TILE_DIMz',TILE_DIMz
      write(6,'(a,i8)') 'TILE_DIM',TILE_DIM
#endif
      write(6,'(a)') '*******************************************'
      ! info gpu
!      call get_memory_gpu(mymemory,totmemory)
!      call print_memory_registration_gpu(6,'DEVICE memory occupied at the start', &
!         'total DEVICE memory',mymemory,totmemory)
      call flush(6)
   endif
   !*************************************time loop************************
   call integrator


   if(ldiagnostic)then
      call printSimulationTime()
      call print_timing_final(idiagnostic,itime_counter, &
         itime_start,1,1,6)
      call get_memory(smemory)
      call get_totram(sram)
      call print_memory_registration(6,&
         'Occupied memory after setup MPI','Total memory',smemory,sram)
   endif


   
   tot_time=ts2-ts1
   Glups = real(lx,db)*real(ly,db)*real(lz,db)*real(nsteps,db)/tot_time*1.e-9_db
#ifdef _NVML               
   Pavg=step_energy*real(nsteps,db)/tot_time
   if (Pavg > ZERO) eff_glups_per_w = Glups / Pavg
   E_tot = step_energy * real(nsteps,db)
   J_per_LUP = step_energy / ( real(lx,db)*real(ly,db)*real(lz,db) )
#endif  

   if(myrank==0)then
      write(6,*) 'time elapsed: ', tot_time, ' s of your life time'
      write(6,*) 'glups: ',  Glups
#ifdef _NVML    
      write(6,*) 'Pavg: ',  Pavg
      write(6,*) 'eff_glups_per_w: ',  eff_glups_per_w
      write(6,*) 'J_per_LUP: ',  J_per_LUP
      write(6,*) 'E_tot: ',  E_tot
      write(6,*) 'E_tot_per_rank: ',  E_tot / real(nprocs,db)
      write(6,*) 'step_energy: ',  step_energy
      write(6,*) 'step_energy_per_rank: ', step_energy / real(nprocs,db)
#endif
   endif  

   call get_memory(smemory)
   call get_totram(sram)
   call print_memory_registration(6,&
      'Occupied memory on HOST at the end','Total HOST memory',smemory,sram)

!   call get_memory_gpu(mymemory,totmemory)
!   call print_memory_registration_gpu(6,'DEVICE memory occupied at the end', &
!      'total DEVICE memory',mymemory,totmemory)
   call compare_benchmark
   call dostop('program correctly closed')

end program
