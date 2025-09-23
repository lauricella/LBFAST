
program accLLB

#ifdef _OPENACC
   use openacc
#endif    
   use vars
   use mpi_template, only : myrank,proc_x,proc_y,proc_z,pbc_x,pbc_y,pbc_z, &
    coords,start_mpi,setup_mpi,exchange_isf_sendrecv,exchange_isf_intpbc, &
    exchange_isf_wait,dostop,print_version_code
   use prints, only : intstr,copystring,driver_read_isfluid_raw, &
    driver_read_init_raw
   use lb_cuda_kernels, only : setup_cuda,tile_dimx,tile_dimy,tile_dimz, &
    tile_dim
    
   implicit none
    
   
   logical :: lexist
   integer, parameter :: inputio=24
   
   namelist /simulation/ nsteps,stamp,stamp2D,stepskip,lreadisfluid, &
    lreadinit,lx,ly,lz,lprint,lvtk,lraw,lrestart,pbc_x,pbc_y,pbc_z,lasync, &
#if defined(_KERNELCUDA) && defined(_OPENACC)
    TILE_DIMx,TILE_DIMy,TILE_DIMz,TILE_DIM, &
#endif    
    nplanes,stamp_term,time_limit,every_time_check,lwriterestart,lweakscaling
    
   namelist /fluid/ fx,fy,fz,visc1, &
#ifdef TWOCOMPONENT
    tau_diff,radius,width,sigma,wettab_r,wettab_b, &
    visc2, &
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
    
end program accLLB
