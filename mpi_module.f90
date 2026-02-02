#include "defines.h"
module mpi_template
   
   use vars, only: db,strdb,isf,i,j,k,nx,ny,nz,lx,ly,lz,rho,f,isfluid,l,ll,opp,&
      ex,ey,ez,nlinks,filenamevtk,namevarvtk,sevt1,sevt2,dir_out,write_fmtnumb2, &
      write_fmtnumb,headervtk,nheadervtk,vtkoffset,ndatavtk,footervtk,printdb, &
      rhoprint,velprint,pressprint,space_fmtnumb,mxln,sevt3,arr_3d,ndir, &
      nxskip,nyskip,nzskip,lxskip,lyskip,lzskip,stepskip,nplanes,skip_npoint,npoint, &
#ifdef _OPENACC
      devType, &
#endif
      flip,flop,rho_r,rho_b, &
      physic_type,acc_device_radeon, &
      nhfields,nphifields,auxfields,nauxfields,forces,nforces, &
      TILE_DIMx,TILE_DIMy,TILE_DIMz,TILE_DIM,nxblock,nyblock,nzblock,nxyblock,nblocks
#ifdef _OPENACC
   use openacc
#endif
#ifdef MPI  
   use mpi
#endif
   implicit none
#ifdef MPI
   !include 'mpif.h'
#endif
   
   integer :: STRMPIREAL
   integer :: MYMPIREAL
   integer :: MYMPIINTS

   integer, save :: nprocs,myrank,lbecomm,localcomm

   integer :: mydev, ndev
   integer :: file_offset
   integer :: proc_x,proc_y,proc_z
   integer :: pbc_x,pbc_y,pbc_z
   integer :: mem_stop
   logical :: rreorder
   integer, parameter::  mpid=3     ! mpi dimension
   logical :: periodic(mpid)
   integer :: prgrid(mpid)
   integer:: coords(mpid)
   integer:: up(mpid),down(mpid),left(mpid)
   integer:: front(mpid),rear(mpid),right(mpid)
   integer, allocatable, dimension(:) :: xinidom,xfindom
   integer, allocatable, dimension(:) :: yinidom,yfindom
   integer, allocatable, dimension(:) :: zinidom,zfindom
   
   integer, allocatable, dimension(:) :: io_comm2d
   logical, allocatable, dimension(:) :: lnoparallel2d
   
   logical, allocatable, dimension(:,:) :: log_plane2d

   integer :: right_send_x,left_recv_x
   integer :: left_send_x,right_recv_x
   integer :: up_send_y,down_recv_y
   integer :: down_send_y, up_recv_y
   integer :: front_send_z,rear_recv_z
   integer :: rear_send_z,front_recv_z

   integer, dimension(mpid) :: myoffset,lsizes,start_idx,gsizes,end_idx
   integer, dimension(mpid) :: skip_myoffset,skip_lsizes,skip_start_idx,skip_gsizes,skip_end_idx
   integer, allocatable, dimension(:,:,:) :: myoffset_plane2d
   integer, allocatable, dimension(:,:,:) :: lsizes_plane2d
   integer, allocatable, dimension(:,:) :: lsizes_3d,myoffset_3d

   integer, parameter :: nlinksmpi=26

   integer, dimension(0:nlinksmpi), parameter :: &
   ! 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26
      exmpi=(/0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 0, 0, 0, 0, 1,-1,-1, 1, 1,-1, 1,-1,-1, 1, 1,-1/)
   integer, dimension(0:nlinksmpi), parameter :: &
      eympi=(/0, 0, 0, 1,-1, 0, 0, 1,-1,-1, 1, 1,-1, 1,-1, 0, 0, 0, 0, 1,-1,-1, 1,-1, 1,-1, 1/)
   integer, dimension(0:nlinksmpi), parameter :: &
      ezmpi=(/0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 1,-1,-1, 1, 1,-1, 1,-1, 1,-1, 1,-1, 1,-1,-1, 1/)
   integer, dimension(0:nlinksmpi), parameter ::&
      oppmpi=(/0, 2, 1, 4, 3, 6, 5, 8, 7,10, 9,12,11,14,13,16,15,18,17,20,19,22,21,24,23,26,25/)

   
   integer, save :: nlinks_faces,nlinks_edges,nlinks_corners,nlinks_max
   integer, save :: nfaces,nedges,ncorners
   logical, save :: lintbb=.false.
   integer, dimension(nlinksmpi) :: send_dir,recv_dir
   logical, dimension(nlinksmpi) :: lsendpop_dir,lrecvpop_dir
   logical, dimension(nlinksmpi) :: lsend_dir,lrecv_dir
   logical, dimension(nlinksmpi) :: lintpbc_dir
   
   logical, dimension(3,nlinksmpi) :: intpbc_dir
   integer, dimension(3,nlinksmpi) :: send_dir_coord,recv_dir_coord


   integer, allocatable, save, dimension(:,:) :: links_faces,links_edges, &
      links_corners,links_pops
   
   integer, allocatable, save, dimension(:,:) :: f_send_extr,f_recv_extr
   integer, allocatable, save, dimension(:,:) :: fvec_send_extr,fvec_recv_extr
   integer, allocatable, save, dimension(:,:) :: b_send_extr,b_recv_extr
   integer, allocatable, save, dimension(:,:) :: c_send_extr,c_recv_extr

   integer, save, dimension(nlinksmpi) :: f_num_extr,fvec_num_extr, &
      b_num_extr,c_num_extr,i_num_extr
   integer, save :: f_numtot_extr,fvec_numtot_extr, &
      b_numtot_extr,c_numtot_extr,i_numtot_extr
   integer, allocatable, save, dimension(:) :: num_links_pops

   integer, dimension(13) :: f_datampi,fvec_datampi,b_datampi,c_datampi,i_datampi
   
   integer, parameter :: num_hfields_datampi=10 
#ifdef TWOCOMPONENT
   integer, parameter :: num_auxfields_datampi=7    ! 3 norm unit vec ! 1 modgrad ! 3 arr_ 
#else
   integer, parameter :: num_auxfields_datampi=0
#endif
   integer, parameter :: num_phifields_datampi=1
   
   integer, parameter :: num_forces_datampi=3 

   real(kind=strdb), allocatable, save, dimension(:) :: f_send_buffmpi,f_recv_buffmpi
   integer, dimension(nlinksmpi), save :: f_nbuffmpi_send,f_nbuffmpi_recv

   real(kind=strdb), allocatable, save, dimension(:) :: fvec_send_buffmpi,fvec_recv_buffmpi
   integer, dimension(nlinksmpi), save :: fvec_nbuffmpi_send,fvec_nbuffmpi_recv

   real(kind=strdb), allocatable, save, dimension(:) :: b_send_buffmpi,b_recv_buffmpi
   integer, dimension(nlinksmpi), save :: b_nbuffmpi_send,b_nbuffmpi_recv
   
   real(kind=strdb), allocatable, save, dimension(:) :: c_send_buffmpi,c_recv_buffmpi
   integer, dimension(nlinksmpi), save :: c_nbuffmpi_send,c_nbuffmpi_recv

   integer(kind=isf), allocatable, save, dimension(:) :: i_send_buffmpi,i_recv_buffmpi
   integer, dimension(nlinksmpi), save :: i_nbuffmpi_send,i_nbuffmpi_recv

   integer, dimension(nlinksmpi), save :: mpitag

   integer, dimension(nlinksmpi), save :: f_mpitag
   integer, dimension(nlinksmpi), save :: fvec_mpitag
   integer, dimension(nlinksmpi), save :: b_mpitag
   integer, dimension(nlinksmpi), save :: c_mpitag
   integer, dimension(nlinksmpi), save :: i_mpitag

   integer, dimension(nlinksmpi*2), save :: reqs

   integer, dimension(nlinksmpi*2), save :: f_reqs
   integer, dimension(nlinksmpi*2), save :: fvec_reqs
   integer, dimension(nlinksmpi*2), save :: b_reqs
   integer, dimension(nlinksmpi*2), save :: c_reqs
   integer, dimension(nlinksmpi*2), save :: i_reqs 
#ifdef TWOCOMPONENT
   logical, parameter :: ltwocomp=.true.
#else   
   logical, parameter :: ltwocomp=.false.
#endif
   integer, save :: nreqs

   integer, save :: nf_reqs
   integer, save :: nfvec_reqs
   integer, save :: nb_reqs
   integer, save :: nc_reqs
   integer, save :: ni_reqs

#ifdef REPULSIVE_FLUX
   integer, parameter :: nbuff=3    !phifields
   integer, parameter :: nbuffvec=3 !auxfields
#else
   integer, parameter :: nbuff=1    !phifields
   integer, parameter :: nbuffvec=1 !auxfields
#endif
   
   integer, parameter :: nbuffbvec=1 !hfields


contains

   subroutine start_mpi

      implicit none

      integer:: ierr

#ifdef MPI
!
      call mpi_init(ierr)
      call MPI_comm_size(MPI_COMM_WORLD, nprocs, ierr)
      call MPI_comm_rank(MPI_COMM_WORLD, myrank, ierr)

#else
      nprocs=1
      myrank=0
      proc_x=1
      proc_y=1
      proc_z=1
#endif


   end subroutine start_mpi
!
   subroutine setup_mpi

      use iso_c_binding
      implicit none
!
      integer:: uni,lopp,idrank,oi,oj,ok
      integer:: ierr, ijlen,istat
! mpi variables
      integer,dimension(mpid) :: temp_coord
      logical :: lcheck=.false.
!
      real(db):: knorm
      character(len=256) :: gpu_env
      integer :: env_length,actual_dev 
      integer :: name_len,myierr
#ifdef MPI      
      character(len=MPI_MAX_PROCESSOR_NAME) :: hname
#else
      character(len=mxln) :: hname
      interface
      function gethostname(name, len) bind(C, name="gethostname")
        import :: c_char, c_int
        character(kind=c_char), dimension(*) :: name
        integer(c_int), value :: len
        integer(c_int) :: gethostname
      end function gethostname
      end interface      
#endif      
!
#ifdef MPI
      if(db==4)then
         MYMPIREAL = MPI_REAL
      elseif(db==8)then
         MYMPIREAL = MPI_DOUBLE_PRECISION
      else
         write(6,*)'ERROR db not defined'
         call dostop
      endif
      if(strdb==2)then
         STRMPIREAL = MPI_REAL2
      elseif(strdb==4)then
         STRMPIREAL = MPI_REAL
      elseif(strdb==8)then
         STRMPIREAL = MPI_DOUBLE_PRECISION
      else
         write(6,*)'ERROR strdb not defined'
         call dostop
      endif      
      if(isf==4)then
         MYMPIINTS=MPI_INTEGER
      elseif(isf==1)then
         MYMPIINTS=MPI_INTEGER1
      else
         write(6,*)'ERROR isf not defined'
         call dostop
      endif
#endif
      knorm = 1.0/1024.0


      nx = lx/proc_x
      ny = ly/proc_y
      nz = lz/proc_z

!
! some check
      lcheck=.false.
      if((nx*proc_x).NE.lx) then
         write(6,*) "ERROR: global and local size along x not valid!!" &
         &                      , lx, nx, proc_x
         lcheck=.true.
      endif
!
      if((ny*proc_y).NE.ly) then
         write(6,*) "ERROR: global and local size along y not valid!!" &
         &                      , ly, ny, proc_y
         lcheck=.true.
      endif
!
      if((nz*proc_z).NE.lz) then
         write(6,*) "ERROR: global and local size along z not valid!!" &
         &                      , lz, nz, proc_z
         lcheck=.true.
      endif

      if(lcheck)then
#ifdef MPI
         call MPI_finalize(ierr)
#endif
         stop
      endif
!
      if(stepskip==0)then
        call doerror(6,'ERROR: stepskip equal to zero')
      endif
      
      nxskip = nx/stepskip
      nyskip = ny/stepskip
      nzskip = nz/stepskip
      
      lcheck=.false.
      if((nxskip*stepskip).NE.nx) then
         write(6,*) "ERROR: stepskip along x not valid!!" &
         &                      , nx, nxskip, stepskip
         lcheck=.true.
      endif
!
      if((nyskip*stepskip).NE.ny) then
         write(6,*) "ERROR: stepskip along y not valid!!" &
         &                      , ny, nyskip, stepskip
         lcheck=.true.
      endif
!
      if((nzskip*stepskip).NE.nz) then
         write(6,*) "ERROR: stepskip along z not valid!!" &
         &                      , nz, nzskip, stepskip
         lcheck=.true.
      endif

      if(lcheck)then
#ifdef MPI
         call MPI_finalize(ierr)
#endif
         stop
      endif
      
      lxskip = lx/stepskip
      lyskip = ly/stepskip
      lzskip = lz/stepskip
      
      lcheck=.false.
      if((lxskip*stepskip).NE.lx) then
         write(6,*) "ERROR: stepskip along x not valid!!" &
         &                      , lx, lxskip, stepskip
         lcheck=.true.
      endif
!
      if((lyskip*stepskip).NE.ly) then
         write(6,*) "ERROR: stepskip along y not valid!!" &
         &                      , ly, lyskip, stepskip
         lcheck=.true.
      endif
!
      if((lzskip*stepskip).NE.lz) then
         write(6,*) "ERROR: stepskip along z not valid!!" &
         &                      , lz, lzskip, stepskip
         lcheck=.true.
      endif

      if(lcheck)then
#ifdef MPI
         call MPI_finalize(ierr)
#endif
         stop
      endif
      
      allocate(skip_npoint(nplanes))
      if(nplanes>0)then
        do i=1,nplanes
          skip_npoint(i)=nint(real(npoint(i),kind=db)/real(stepskip,kind=db))
          if(skip_npoint(i)<1)skip_npoint(i)=1
          if(ndir(i)==1 .and. skip_npoint(i)>lxskip)skip_npoint(i)=lxskip
          if(ndir(i)==2 .and. skip_npoint(i)>lyskip)skip_npoint(i)=lyskip
          if(ndir(i)==3 .and. skip_npoint(i)>lzskip)skip_npoint(i)=lzskip
          if(skip_npoint(i)*stepskip.ne.npoint(i))then
            npoint(i)=skip_npoint(i)*stepskip
            if(myrank==0)then
              write(6,'(a,i0,a,i0)')'PRINT 2D: npoint for the ',i, &
               ' plane was reset to ',npoint(i)
            endif
          endif
        enddo
      endif

#ifdef _OPENACC
    ndev = acc_get_num_devices(devType)
    write(6,'(a,2i8)') "number of device on each node: ", ndev, myrank
    if (ndev == 0) then
        if (myrank == 0) write(6,*) 'WARNING: No GPUs found:', ndev
        call dostop
    endif

#ifdef CRAY
    ! CRAY features: reading ROCR_VISIBLE_DEVICES
    if(myrank==0)then
      call get_environment_variable("ROCR_VISIBLE_DEVICES", gpu_env, env_length)
      write(6,'(a,i0,2a)')"Rank ", myrank, " sees ROCR_VISIBLE_DEVICES: ", trim(gpu_env)
    endif
 !   read(gpu_env, '(I0)') mydev
 !   write(6,'(a,2i8)') "Device to set for the rank on each node: ", mydev, myrank
#else
    call get_environment_variable("CUDA_VISIBLE_DEVICES", gpu_env, env_length)
     write(6,'(a,i0,2a)') "Rank ", myrank, " sees CUDA_VISIBLE_DEVICES=", trim(gpu_env)
#endif
    mydev = mod(myrank, ndev)
    call acc_set_device_num(mydev, devType)
    actual_dev = acc_get_device_num(devType)
#endif

#ifdef MPI
      hname=repeat(' ',MPI_MAX_PROCESSOR_NAME)
      call MPI_Get_processor_name(hname, name_len, myierr)
#else
      hname=repeat(' ',mxln)
      myierr = gethostname(hname, mxln)      
#endif

      call flush(6)
#ifdef MPI
      call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif
      do idrank=0,nprocs-1
         if(idrank==myrank)then  !
            write(6,'(a,i3,2a)') 'Rank ', myrank, ' on host: ', trim(hname)
#ifdef _OPENACC
            write(6,'(a,i4,a,i4)') 'I want to set the devide num ', mydev, ' for the mpi process', myrank
            ! DEBUG EXTRA
            write(6,'(a,i8,a,i8)') 'DEBUG: acc_get_num_devices = ', acc_get_num_devices(devType),' in rank ', myrank
            write(6,'(a,i4,a,i4)') 'The MPI process ', myrank, ' is using the GPU ', mydev
            write(6,'(a,i4,a,i4)') 'CHECK: Process ', myrank, ' is REALLY using GPU ', actual_dev
#endif               
            call flush(6)
         endif
#ifdef MPI
         call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif
      enddo
      call flush(6)


      rreorder=.false.

      periodic(1) = (pbc_x==1)
      periodic(2) = (pbc_y==1)
      periodic(3) = (pbc_z==1)

      prgrid(1) = proc_x!proc_x
      prgrid(2) = proc_y!proc_y
      prgrid(3) = proc_z

      if(myrank==0)write(6,'(a,3i4)')'MPI processes= ',proc_x,proc_y,proc_z
      if(myrank==0)write(6,'(a,3i4)')'pbc applied= ',pbc_x,pbc_y,pbc_z
      call flush(6)


#ifdef MPI
      !$acc wait
      call mpi_barrier(MPI_COMM_WORLD,ierr)
!
! set the gpu to the task id
      call MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, &
         MPI_INFO_NULL, localcomm, ierr)
      call MPI_Comm_rank(localcomm, mydev, ierr)
      !call MPI_get_processor_name(hname,ijlen,ierr)

!
! check
      if ((proc_x*proc_y*proc_z).ne.nprocs) then
         if (myrank.eq.0) then
            write(6,'(a,3i4)') 'ERROR: decomposed for x y z procs ', &
               proc_x,proc_y,proc_z
            write(6,'(a,i4,a)') 'ERROR: decomposed for', &
               proc_x*proc_y*proc_z, 'procs'
            write(6,'(a,i4,a)') 'ERROR: launched on', nprocs, 'processes'
         end if
         call dostop('ERROR mpi job lunched with wrong number of processes')
      end if
!
!
! building virtual topology
      call MPI_cart_create(mpi_comm_world, mpid, prgrid, &
         periodic,rreorder,lbecomm,ierr)

      call MPI_comm_rank(lbecomm, myrank, ierr)
      call MPI_cart_coords(lbecomm, myrank, mpid, &
         coords, ierr)

      call mpi_barrier(MPI_COMM_WORLD,ierr)
!

      !pause
      !********************************************************!
      !*******************neighbour processes******************!
      !********************************************************!
      !call MPI_ERRHANDLER_SET(lbecomm, MPI_ERRORS_RETURN, ierr)
      do l=1,nlinksmpi
         mpitag(l) = 400 + l
         f_mpitag(l) = 500 + l
         b_mpitag(l) = 600 + l
         c_mpitag(l) = 300 + l
         i_mpitag(l) = 700 + l
         fvec_mpitag(l) = 900 + l
         temp_coord(1) = coords(1) + exmpi(l)
         temp_coord(2) = coords(2) + eympi(l)
         temp_coord(3) = coords(3) + ezmpi(l)
         !call MPI_Cart_rank(lbecomm, temp_coord, send_dir(l),ierr)

         oi=temp_coord(1)
         oj=temp_coord(2)
         ok=temp_coord(3)
         !oi=mod(oi+nx-1,nx)+1
         if(periodic(1))then
            oi=mod(oi+proc_x,proc_x)
         else
            oi=min(max(oi,0),proc_x-1)
         endif
         if(periodic(2))then
            oj=mod(oj+proc_y,proc_y)
         else
            oj=min(max(oj,0),proc_y-1)
         endif
         if(periodic(3))then
            ok=mod(ok+proc_z,proc_z)
         else
            ok=min(max(ok,0),proc_z-1)
         endif

         send_dir_coord(1:3,l)=[oi,oj,ok]
         call MPI_Cart_rank(lbecomm, [oi,oj,ok], send_dir(l),ierr)


         lopp=oppmpi(l)
         temp_coord(1) = coords(1) + exmpi(lopp)
         temp_coord(2) = coords(2) + eympi(lopp)
         temp_coord(3) = coords(3) + ezmpi(lopp)
         !call MPI_Cart_rank(lbecomm, temp_coord, recv_dir(l),ierr)

         oi=temp_coord(1)
         oj=temp_coord(2)
         ok=temp_coord(3)
!        if(periodic(1))oi=mod(oi+proc_x,proc_x)
!        if(periodic(2))oj=mod(oj+proc_y,proc_y)
!        if(periodic(3))ok=mod(ok+proc_z,proc_z)
         if(periodic(1))then
            oi=mod(oi+proc_x,proc_x)
         else
            oi=min(max(oi,0),proc_x-1)
         endif
         if(periodic(2))then
            oj=mod(oj+proc_y,proc_y)
         else
            oj=min(max(oj,0),proc_y-1)
         endif
         if(periodic(3))then
            ok=mod(ok+proc_z,proc_z)
         else
            ok=min(max(ok,0),proc_z-1)
         endif

         recv_dir_coord(1:3,l)=[oi,oj,ok]
         call MPI_Cart_rank(lbecomm, [oi,oj,ok],recv_dir(l),ierr)

      enddo

      file_offset = 0    !to check

#else

      coords=0

      do l=1,nlinksmpi
         !vado sempre a me stesso perche sono il solo processo
         send_dir_coord(1:3,l)=coords(1:3)
         send_dir(l)=myrank

         recv_dir_coord(1:3,l)=coords(1:3)
         recv_dir(l)=myrank
      enddo

#endif
      !gestisci se fare o no send e receive (deve andare su nodi diversi e non sfondare il range coords se non periodico)
      do l=1,nlinksmpi

         lsend_dir(l)=(myrank .ne. send_dir(l))
         if(lsend_dir(l))then
            temp_coord(1) = coords(1) + exmpi(l)
            temp_coord(2) = coords(2) + eympi(l)
            temp_coord(3) = coords(3) + ezmpi(l)

            oi=temp_coord(1)
            oj=temp_coord(2)
            ok=temp_coord(3)

            if(periodic(1))then
               oi=mod(oi+proc_x,proc_x)
            endif
            if(periodic(2))then
               oj=mod(oj+proc_y,proc_y)
            endif
            if(periodic(3))then
               ok=mod(ok+proc_z,proc_z)
            endif
            !se sfondo allora non devo inviare
            if(oi<0 .or. oj<0 .or. ok<0 .or. &
               oi>=proc_x .or. oj>=proc_y .or. ok>=proc_z)then
               lsend_dir(l)=.false.
            endif
         endif

         lrecv_dir(l)=(myrank .ne. recv_dir(l))
         if(lrecv_dir(l))then
            lopp=oppmpi(l)
            temp_coord(1) = coords(1) + exmpi(lopp)
            temp_coord(2) = coords(2) + eympi(lopp)
            temp_coord(3) = coords(3) + ezmpi(lopp)

            oi=temp_coord(1)
            oj=temp_coord(2)
            ok=temp_coord(3)
            if(periodic(1))then
               oi=mod(oi+proc_x,proc_x)
            endif
            if(periodic(2))then
               oj=mod(oj+proc_y,proc_y)
            endif
            if(periodic(3))then
               ok=mod(ok+proc_z,proc_z)
            endif
            !se sfondo allora non devo ricevere
            if(oi<0 .or. oj<0 .or. ok<0 .or. &
               oi>=proc_x .or. oj>=proc_y .or. ok>=proc_z)then
               lrecv_dir(l)=.false.
            endif
         endif

      enddo


#ifdef MPI
      call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif

      !gestisci se fare o no le pbc interne (deve andare sullo stesso nodo in send e recv e non sfondare il range coords se non periodico)
      do l=1,nlinksmpi
         lintpbc_dir(l)=((myrank==send_dir(l)) .and. (myrank==recv_dir(l)))
         intpbc_dir(1,l)=.false.
         intpbc_dir(2,l)=.false.
         intpbc_dir(3,l)=.false.
         if(lintpbc_dir(l))then
            temp_coord(1) = coords(1) + exmpi(l)
            temp_coord(2) = coords(2) + eympi(l)
            temp_coord(3) = coords(3) + ezmpi(l)

            oi=temp_coord(1)
            oj=temp_coord(2)
            ok=temp_coord(3)

            if(periodic(1))then
               oi=mod(oi+proc_x,proc_x)
            endif
            if(periodic(2))then
               oj=mod(oj+proc_y,proc_y)
            endif
            if(periodic(3))then
               ok=mod(ok+proc_z,proc_z)
            endif
            !non devo sfondare se non periodico nella direzione specifica l
            if(oi<0 .or. oj<0 .or. ok<0 .or. &
               oi>=proc_x .or. oj>=proc_y .or. ok>=proc_z)then
               lintpbc_dir(l)=.false.
            endif
         endif
         !se sono periodico nel mio processo allora non devo sfondare neanche quando ricevo
         if(lintpbc_dir(l))then
            lopp=oppmpi(l)
            temp_coord(1) = coords(1) + exmpi(lopp)
            temp_coord(2) = coords(2) + eympi(lopp)
            temp_coord(3) = coords(3) + ezmpi(lopp)

            oi=temp_coord(1)
            oj=temp_coord(2)
            ok=temp_coord(3)
            if(periodic(1))then
               oi=mod(oi+proc_x,proc_x)
            endif
            if(periodic(2))then
               oj=mod(oj+proc_y,proc_y)
            endif
            if(periodic(3))then
               ok=mod(ok+proc_z,proc_z)
            endif
            !non devo sfondare se non periodico nella direzione specifica l
            if(oi<0 .or. oj<0 .or. ok<0 .or. &
               oi>=proc_x .or. oj>=proc_y .or. ok>=proc_z)then
               lintpbc_dir(l)=.false.
            endif
         endif
         !sto in ballo devo fare periodico interno se true
         !allora storo se sono periodico per direzione l e condizioni lungo i tre assi
         if(lintpbc_dir(l))then
            !se lungo il vettore l mi muovo lungo un asse e sono periodico allora intpbc_dir è true
            intpbc_dir(1,l)=(abs(exmpi(l))==1 .and.periodic(1))
            intpbc_dir(2,l)=(abs(eympi(l))==1 .and.periodic(2))
            intpbc_dir(3,l)=(abs(ezmpi(l))==1 .and.periodic(3))
         endif
      enddo


      !setto variabili utili in particolare per MPI-IO
      !offset in coordinate globali di ogni processo MPI
      myoffset(1) = coords(1)*nx
      myoffset(2) = coords(2)*ny
      myoffset(3) = coords(3)*nz
      !dimensione locale grid di ogni processo MPI
      lsizes(1)=nx
      lsizes(2)=ny
      lsizes(3)=nz
      !dimensione globale grid di ogni processo MPI
      gsizes(1)=lx
      gsizes(2)=ly
      gsizes(3)=lz
      !inizio dimensione locale grid in coordinate globali di ogni processo MPI
      start_idx(1)=myoffset(1)+1
      start_idx(2)=myoffset(2)+1
      start_idx(3)=myoffset(3)+1
      !fine dimensione locale grid in coordinate globali di ogni processo MPI
      end_idx(1) = myoffset(1)+lsizes(1)
      end_idx(2) = myoffset(2)+lsizes(2)
      end_idx(3) = myoffset(3)+lsizes(3)
      !!!!!!stepskip
      !offset in coordinate globali di ogni processo MPI
      skip_myoffset(1) = coords(1)*nxskip
      skip_myoffset(2) = coords(2)*nyskip
      skip_myoffset(3) = coords(3)*nzskip
      !dimensione locale grid di ogni processo MPI
      skip_lsizes(1)=nxskip
      skip_lsizes(2)=nyskip
      skip_lsizes(3)=nzskip
      !dimensione globale grid di ogni processo MPI
      skip_gsizes(1)=lxskip
      skip_gsizes(2)=lyskip
      skip_gsizes(3)=lzskip
      !inizio dimensione locale grid in coordinate globali di ogni processo MPI
      skip_start_idx(1)=skip_myoffset(1)+1
      skip_start_idx(2)=skip_myoffset(2)+1
      skip_start_idx(3)=skip_myoffset(3)+1
      !fine dimensione locale grid in coordinate globali di ogni processo MPI
      skip_end_idx(1) = skip_myoffset(1)+skip_lsizes(1)
      skip_end_idx(2) = skip_myoffset(2)+skip_lsizes(2)
      skip_end_idx(3) = skip_myoffset(3)+skip_lsizes(3)

      !questo mi serve per la funzione GET_COORD_POINT
      !funzione di debug che dalle coordinate generali x y z mi da le coordinate della decomposizione MPI
      !così so su quale nodo vengono lavorate le coordinate x y z
      allocate(xinidom(0:proc_x-1))
      allocate(xfindom(0:proc_x-1))
      xinidom(:)=0
      xfindom(:)=-1
      xinidom(0)=1
      xfindom(0)=nx+xinidom(0)-1
      do i=1,proc_x-1
         xinidom(i)=xfindom(i-1)+1
         xfindom(i)=xinidom(i)+nx-1
      enddo

      allocate(yinidom(0:proc_y-1))
      allocate(yfindom(0:proc_y-1))
      yinidom(:)=0
      yfindom(:)=-1
      yinidom(0)=1
      yfindom(0)=ny+yinidom(0)-1
      do i=1,proc_y-1
         yinidom(i)=yfindom(i-1)+1
         yfindom(i)=yinidom(i)+ny-1
      enddo

      allocate(zinidom(0:proc_z-1))
      allocate(zfindom(0:proc_z-1))
      zinidom(:)=0
      zfindom(:)=-1
      zinidom(0)=1
      zfindom(0)=nz+zinidom(0)-1
      do i=1,proc_z-1
         zinidom(i)=zfindom(i-1)+1
         zfindom(i)=zinidom(i)+nz-1
      enddo


      call flush(6)
#ifdef MPI
      call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif
      do idrank=0,nprocs
         if(idrank==myrank)then  !
            write(6,'(a,i4,a,9i8)')'DEC myrank ',myrank,' coords ',coords,&
               start_idx(1),end_idx(1),&
               start_idx(2),end_idx(2),&
               start_idx(3),end_idx(3)
            call flush(6)
         endif
#ifdef MPI
         call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif
      enddo
      call flush(6)
#ifdef VERBOSE

#ifdef MPI
      call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif

      do idrank=0,nprocs-1
         do l=1,nlinksmpi
            if(idrank==myrank)then
               write(6,'(a,i4,a,4i3,a,i4,a,i4,2l2,a,3l2)')'Myrank ',myrank,' l ',l,exmpi(l),eympi(l),ezmpi(l),&
                  ' recv ',recv_dir(l),' send ',send_dir(l),lrecv_dir(l),lsend_dir(l),' intpbc ',intpbc_dir(1:3,l)
               call flush(6)
            endif
#ifdef MPI
            call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif
         enddo
      enddo
      call flush(6)
#ifdef MPI
      call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif

#endif
      !trovo per ogni direzione l quali popolazioni devono essere inviate e le storo in links_faces
      !occhio che devo vedere quali popolazioni mandare ma il lattice può essere
      !tipo d3q15 o d3q19 e quindi non devo mandare nulla
      !exmpi eympi ezmpi sono le 26 direzioni di MPI
      !ex ey ez sono le direzioni del lattice d3q15 o d3q19 o d3q27
      allocate(num_links_pops(1:nlinksmpi))
      !faces
      nfaces=6
      nlinks_faces=9
      allocate(links_faces(1:nlinks_faces,1:6))
      do l=1,6
         nlinks_faces=0
         do ll=1,nlinks
            if((abs(exmpi(l))==1 .and. exmpi(l)==ex(ll)) .or. &
               (abs(eympi(l))==1 .and. eympi(l)==ey(ll)) .or. &
               (abs(ezmpi(l))==1 .and. ezmpi(l)==ez(ll)))then
               nlinks_faces=nlinks_faces+1
               links_faces(nlinks_faces,l)=ll
#ifdef VERBOSE
               if(myrank==0)write(6,'(a,i3,a,3i3,a,a2,a,3i3)')'dir l ',l,' disp ',&
                  exmpi(l),eympi(l),ezmpi(l),&
                  ' f',write_fmtnumb2(ll),' dir ',ex(ll),ey(ll),ez(ll)
               call flush(6)
#endif
            endif
         enddo
         num_links_pops(l)=nlinks_faces
      enddo
      !edges
      nedges=12
      nlinks_edges=3
      allocate(links_edges(1:nlinks_edges,7:18))
      do l=7,18
         nlinks_edges=0
         do ll=1,nlinks
            if((abs(exmpi(l))==1 .and. abs(eympi(l))==1 .and. exmpi(l)==ex(ll) .and. eympi(l)==ey(ll)) .or. &
               (abs(exmpi(l))==1 .and. abs(ezmpi(l))==1 .and. exmpi(l)==ex(ll) .and. ezmpi(l)==ez(ll)) .or. &
               (abs(eympi(l))==1 .and. abs(ezmpi(l))==1 .and. eympi(l)==ey(ll) .and. ezmpi(l)==ez(ll)))then
               nlinks_edges=nlinks_edges+1
               links_edges(nlinks_edges,l)=ll
#ifdef VERBOSE
               if(myrank==0)write(6,'(a,i3,a,3i3,a,a2,a,3i3)')'dir l ',l,' disp ',&
                  exmpi(l),eympi(l),ezmpi(l),&
                  ' f',write_fmtnumb2(ll),' dir ',ex(ll),ey(ll),ez(ll)
               call flush(6)
#endif
            endif
         enddo
         num_links_pops(l)=nlinks_edges
      enddo
      !corner
      ncorners=8
      nlinks_corners=1
      allocate(links_corners(1:nlinks_corners,19:nlinksmpi))
      do l=19,nlinksmpi
         nlinks_corners=0
         do ll=1,nlinks
            if(exmpi(l)==ex(ll) .and. eympi(l)==ey(ll) .and. ezmpi(l)==ez(ll))then
               nlinks_corners=nlinks_corners+1
               links_corners(nlinks_corners,l)=ll
#ifdef VERBOSE
               if(myrank==0)write(6,'(a,i3,a,3i3,a,a2,a,3i3)')'dir l ',l,' disp ',&
                  exmpi(l),eympi(l),ezmpi(l),&
                  ' f',write_fmtnumb2(ll),' dir ',ex(ll),ey(ll),ez(ll)
               call flush(6)
#endif
            endif
         enddo
         num_links_pops(l)=nlinks_corners
      enddo

#ifdef MPI
      call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif
      if(myrank==0)then
         write(6,'(a)')'dir nlinks: faces edges corners'
         write(6,'(13x,i4,2x,i4,4x,i4)')nlinks_faces,nlinks_edges,nlinks_corners
      endif

      !riempio la lista links_pops con le popolazioni da inviare per ogni direzione l
      !il primo indice di lista links_pops è preso dal massimo delle pops da mandare tra faces edges e corners (ovviamente è sempre faces)
      !nota che num_links_pops è il numero di poplazioni da inviare per direzione l
      nlinks_max=max(nlinks_faces,nlinks_edges,nlinks_corners)
      allocate(links_pops(1:nlinks_max,1:nlinksmpi))
      do l=1,6
         do ll=1,num_links_pops(l)
            links_pops(ll,l)=links_faces(ll,l)
         enddo
      enddo
      do l=7,18
         do ll=1,num_links_pops(l)
            links_pops(ll,l)=links_edges(ll,l)
         enddo
      enddo
      do l=19,nlinksmpi
         do ll=1,num_links_pops(l)
            links_pops(ll,l)=links_corners(ll,l)
         enddo
      enddo
#ifdef MPI
      call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif
      deallocate(links_faces,links_edges,links_corners)
      !solo se ci sono popolazioni da mandare
      !allora metti lsendpop_dir e lrecvpop_dir true
      !con d3q19 o d3q15 spesso non devo mandare nulla
      do l=1,nlinksmpi
         lsendpop_dir(l)=(lsend_dir(l) .and. num_links_pops(l)>0)
         lrecvpop_dir(l)=(lrecv_dir(l) .and. num_links_pops(l)>0)
      enddo
      !solo se ci sono popolazioni da fare pbc interno lo fai


      !mi storo gli estremi i j k che devono essere inviati e ricevuti lungo ogni direzione l
      allocate(f_send_extr(6,nlinksmpi))
      allocate(f_recv_extr(6,nlinksmpi))
      allocate(fvec_send_extr(6,nlinksmpi))
      allocate(fvec_recv_extr(6,nlinksmpi))
      allocate(b_send_extr(6,nlinksmpi))
      allocate(b_recv_extr(6,nlinksmpi))
      allocate(c_send_extr(6,nlinksmpi))
      allocate(c_recv_extr(6,nlinksmpi))

      !faces
      do l=1,6
         if(exmpi(l)==1)then

            f_send_extr(1,l)=nx-nbuff+1
            f_send_extr(2,l)=nx
            f_send_extr(3,l)=1
            f_send_extr(4,l)=ny
            f_send_extr(5,l)=1
            f_send_extr(6,l)=nz

            f_recv_extr(1,l)=1-nbuff
            f_recv_extr(2,l)=0
            f_recv_extr(3,l)=1
            f_recv_extr(4,l)=ny
            f_recv_extr(5,l)=1
            f_recv_extr(6,l)=nz
            
            fvec_send_extr(1,l)=nx-nbuffvec+1
            fvec_send_extr(2,l)=nx
            fvec_send_extr(3,l)=1
            fvec_send_extr(4,l)=ny
            fvec_send_extr(5,l)=1
            fvec_send_extr(6,l)=nz

            fvec_recv_extr(1,l)=1-nbuffvec
            fvec_recv_extr(2,l)=0
            fvec_recv_extr(3,l)=1
            fvec_recv_extr(4,l)=ny
            fvec_recv_extr(5,l)=1
            fvec_recv_extr(6,l)=nz

            b_send_extr(1,l)=nx-nbuffbvec+1
            b_send_extr(2,l)=nx
            b_send_extr(3,l)=1
            b_send_extr(4,l)=ny
            b_send_extr(5,l)=1
            b_send_extr(6,l)=nz

            b_recv_extr(1,l)=1-nbuffbvec
            b_recv_extr(2,l)=0
            b_recv_extr(3,l)=1
            b_recv_extr(4,l)=ny
            b_recv_extr(5,l)=1
            b_recv_extr(6,l)=nz
         endif
         if(exmpi(l)==-1)then

            f_send_extr(1,l)=1
            f_send_extr(2,l)=nbuff
            f_send_extr(3,l)=1
            f_send_extr(4,l)=ny
            f_send_extr(5,l)=1
            f_send_extr(6,l)=nz

            f_recv_extr(1,l)=nx+1
            f_recv_extr(2,l)=nx+nbuff
            f_recv_extr(3,l)=1
            f_recv_extr(4,l)=ny
            f_recv_extr(5,l)=1
            f_recv_extr(6,l)=nz
            
            fvec_send_extr(1,l)=1
            fvec_send_extr(2,l)=nbuffvec
            fvec_send_extr(3,l)=1
            fvec_send_extr(4,l)=ny
            fvec_send_extr(5,l)=1
            fvec_send_extr(6,l)=nz

            fvec_recv_extr(1,l)=nx+1
            fvec_recv_extr(2,l)=nx+nbuffvec
            fvec_recv_extr(3,l)=1
            fvec_recv_extr(4,l)=ny
            fvec_recv_extr(5,l)=1
            fvec_recv_extr(6,l)=nz

            b_send_extr(1,l)=1
            b_send_extr(2,l)=nbuffbvec
            b_send_extr(3,l)=1
            b_send_extr(4,l)=ny
            b_send_extr(5,l)=1
            b_send_extr(6,l)=nz

            b_recv_extr(1,l)=nx+1
            b_recv_extr(2,l)=nx+nbuffbvec
            b_recv_extr(3,l)=1
            b_recv_extr(4,l)=ny
            b_recv_extr(5,l)=1
            b_recv_extr(6,l)=nz
         endif
         if(eympi(l)==1)then

            f_send_extr(1,l)=1
            f_send_extr(2,l)=nx
            f_send_extr(3,l)=ny-nbuff+1
            f_send_extr(4,l)=ny
            f_send_extr(5,l)=1
            f_send_extr(6,l)=nz

            f_recv_extr(1,l)=1
            f_recv_extr(2,l)=nx
            f_recv_extr(3,l)=1-nbuff
            f_recv_extr(4,l)=0
            f_recv_extr(5,l)=1
            f_recv_extr(6,l)=nz
            
            fvec_send_extr(1,l)=1
            fvec_send_extr(2,l)=nx
            fvec_send_extr(3,l)=ny-nbuffvec+1
            fvec_send_extr(4,l)=ny
            fvec_send_extr(5,l)=1
            fvec_send_extr(6,l)=nz

            fvec_recv_extr(1,l)=1
            fvec_recv_extr(2,l)=nx
            fvec_recv_extr(3,l)=1-nbuffvec
            fvec_recv_extr(4,l)=0
            fvec_recv_extr(5,l)=1
            fvec_recv_extr(6,l)=nz

            b_send_extr(1,l)=1
            b_send_extr(2,l)=nx
            b_send_extr(3,l)=ny-nbuffbvec+1
            b_send_extr(4,l)=ny
            b_send_extr(5,l)=1
            b_send_extr(6,l)=nz

            b_recv_extr(1,l)=1
            b_recv_extr(2,l)=nx
            b_recv_extr(3,l)=1-nbuffbvec
            b_recv_extr(4,l)=0
            b_recv_extr(5,l)=1
            b_recv_extr(6,l)=nz
         endif
         if(eympi(l)==-1)then

            f_send_extr(1,l)=1
            f_send_extr(2,l)=nx
            f_send_extr(3,l)=1
            f_send_extr(4,l)=nbuff
            f_send_extr(5,l)=1
            f_send_extr(6,l)=nz

            f_recv_extr(1,l)=1
            f_recv_extr(2,l)=nx
            f_recv_extr(3,l)=ny+1
            f_recv_extr(4,l)=ny+nbuff
            f_recv_extr(5,l)=1
            f_recv_extr(6,l)=nz
            
            fvec_send_extr(1,l)=1
            fvec_send_extr(2,l)=nx
            fvec_send_extr(3,l)=1
            fvec_send_extr(4,l)=nbuffvec
            fvec_send_extr(5,l)=1
            fvec_send_extr(6,l)=nz

            fvec_recv_extr(1,l)=1
            fvec_recv_extr(2,l)=nx
            fvec_recv_extr(3,l)=ny+1
            fvec_recv_extr(4,l)=ny+nbuffvec
            fvec_recv_extr(5,l)=1
            fvec_recv_extr(6,l)=nz

            b_send_extr(1,l)=1
            b_send_extr(2,l)=nx
            b_send_extr(3,l)=1
            b_send_extr(4,l)=nbuffbvec
            b_send_extr(5,l)=1
            b_send_extr(6,l)=nz

            b_recv_extr(1,l)=1
            b_recv_extr(2,l)=nx
            b_recv_extr(3,l)=ny+1
            b_recv_extr(4,l)=ny+nbuffbvec
            b_recv_extr(5,l)=1
            b_recv_extr(6,l)=nz
         endif
         if(ezmpi(l)==1)then

            f_send_extr(1,l)=1
            f_send_extr(2,l)=nx
            f_send_extr(3,l)=1
            f_send_extr(4,l)=ny
            f_send_extr(5,l)=nz-nbuff+1
            f_send_extr(6,l)=nz

            f_recv_extr(1,l)=1
            f_recv_extr(2,l)=nx
            f_recv_extr(3,l)=1
            f_recv_extr(4,l)=ny
            f_recv_extr(5,l)=1-nbuff
            f_recv_extr(6,l)=0
            
            fvec_send_extr(1,l)=1
            fvec_send_extr(2,l)=nx
            fvec_send_extr(3,l)=1
            fvec_send_extr(4,l)=ny
            fvec_send_extr(5,l)=nz-nbuffvec+1
            fvec_send_extr(6,l)=nz

            fvec_recv_extr(1,l)=1
            fvec_recv_extr(2,l)=nx
            fvec_recv_extr(3,l)=1
            fvec_recv_extr(4,l)=ny
            fvec_recv_extr(5,l)=1-nbuffvec
            fvec_recv_extr(6,l)=0

            b_send_extr(1,l)=1
            b_send_extr(2,l)=nx
            b_send_extr(3,l)=1
            b_send_extr(4,l)=ny
            b_send_extr(5,l)=nz-nbuffbvec+1
            b_send_extr(6,l)=nz

            b_recv_extr(1,l)=1
            b_recv_extr(2,l)=nx
            b_recv_extr(3,l)=1
            b_recv_extr(4,l)=ny
            b_recv_extr(5,l)=1-nbuffbvec
            b_recv_extr(6,l)=0
         endif
         if(ezmpi(l)==-1)then

            f_send_extr(1,l)=1
            f_send_extr(2,l)=nx
            f_send_extr(3,l)=1
            f_send_extr(4,l)=ny
            f_send_extr(5,l)=1
            f_send_extr(6,l)=nbuff

            f_recv_extr(1,l)=1
            f_recv_extr(2,l)=nx
            f_recv_extr(3,l)=1
            f_recv_extr(4,l)=ny
            f_recv_extr(5,l)=nz+1
            f_recv_extr(6,l)=nz+nbuff
            
            fvec_send_extr(1,l)=1
            fvec_send_extr(2,l)=nx
            fvec_send_extr(3,l)=1
            fvec_send_extr(4,l)=ny
            fvec_send_extr(5,l)=1
            fvec_send_extr(6,l)=nbuffvec

            fvec_recv_extr(1,l)=1
            fvec_recv_extr(2,l)=nx
            fvec_recv_extr(3,l)=1
            fvec_recv_extr(4,l)=ny
            fvec_recv_extr(5,l)=nz+1
            fvec_recv_extr(6,l)=nz+nbuffvec

            b_send_extr(1,l)=1
            b_send_extr(2,l)=nx
            b_send_extr(3,l)=1
            b_send_extr(4,l)=ny
            b_send_extr(5,l)=1
            b_send_extr(6,l)=nbuffbvec

            b_recv_extr(1,l)=1
            b_recv_extr(2,l)=nx
            b_recv_extr(3,l)=1
            b_recv_extr(4,l)=ny
            b_recv_extr(5,l)=nz+1
            b_recv_extr(6,l)=nz+nbuffbvec
         endif
      enddo
      !edges
      do l=7,18
         !!!!   x   y
         if(exmpi(l)==1 .and. eympi(l)==1)then

            f_send_extr(1,l)=nx-nbuff+1
            f_send_extr(2,l)=nx
            f_send_extr(3,l)=ny-nbuff+1
            f_send_extr(4,l)=ny
            f_send_extr(5,l)=1
            f_send_extr(6,l)=nz

            f_recv_extr(1,l)=1-nbuff
            f_recv_extr(2,l)=0
            f_recv_extr(3,l)=1-nbuff
            f_recv_extr(4,l)=0
            f_recv_extr(5,l)=1
            f_recv_extr(6,l)=nz
            
            fvec_send_extr(1,l)=nx-nbuffvec+1
            fvec_send_extr(2,l)=nx
            fvec_send_extr(3,l)=ny-nbuffvec+1
            fvec_send_extr(4,l)=ny
            fvec_send_extr(5,l)=1
            fvec_send_extr(6,l)=nz

            fvec_recv_extr(1,l)=1-nbuffvec
            fvec_recv_extr(2,l)=0
            fvec_recv_extr(3,l)=1-nbuffvec
            fvec_recv_extr(4,l)=0
            fvec_recv_extr(5,l)=1
            fvec_recv_extr(6,l)=nz

            b_send_extr(1,l)=nx-nbuffbvec+1
            b_send_extr(2,l)=nx
            b_send_extr(3,l)=ny-nbuffbvec+1
            b_send_extr(4,l)=ny
            b_send_extr(5,l)=1
            b_send_extr(6,l)=nz

            b_recv_extr(1,l)=1-nbuffbvec
            b_recv_extr(2,l)=0
            b_recv_extr(3,l)=1-nbuffbvec
            b_recv_extr(4,l)=0
            b_recv_extr(5,l)=1
            b_recv_extr(6,l)=nz
         endif
         if(exmpi(l)==-1 .and. eympi(l)==-1)then

            f_send_extr(1,l)=1
            f_send_extr(2,l)=nbuff
            f_send_extr(3,l)=1
            f_send_extr(4,l)=nbuff
            f_send_extr(5,l)=1
            f_send_extr(6,l)=nz

            f_recv_extr(1,l)=nx+1
            f_recv_extr(2,l)=nx+nbuff
            f_recv_extr(3,l)=ny+1
            f_recv_extr(4,l)=ny+nbuff
            f_recv_extr(5,l)=1
            f_recv_extr(6,l)=nz
            
            fvec_send_extr(1,l)=1
            fvec_send_extr(2,l)=nbuffvec
            fvec_send_extr(3,l)=1
            fvec_send_extr(4,l)=nbuffvec
            fvec_send_extr(5,l)=1
            fvec_send_extr(6,l)=nz

            fvec_recv_extr(1,l)=nx+1
            fvec_recv_extr(2,l)=nx+nbuffvec
            fvec_recv_extr(3,l)=ny+1
            fvec_recv_extr(4,l)=ny+nbuffvec
            fvec_recv_extr(5,l)=1
            fvec_recv_extr(6,l)=nz

            b_send_extr(1,l)=1
            b_send_extr(2,l)=nbuffbvec
            b_send_extr(3,l)=1
            b_send_extr(4,l)=nbuffbvec
            b_send_extr(5,l)=1
            b_send_extr(6,l)=nz

            b_recv_extr(1,l)=nx+1
            b_recv_extr(2,l)=nx+nbuffbvec
            b_recv_extr(3,l)=ny+1
            b_recv_extr(4,l)=ny+nbuffbvec
            b_recv_extr(5,l)=1
            b_recv_extr(6,l)=nz
         endif
         if(exmpi(l)==1 .and. eympi(l)==-1)then

            f_send_extr(1,l)=nx-nbuff+1
            f_send_extr(2,l)=nx
            f_send_extr(3,l)=1
            f_send_extr(4,l)=nbuff
            f_send_extr(5,l)=1
            f_send_extr(6,l)=nz

            f_recv_extr(1,l)=1-nbuff
            f_recv_extr(2,l)=0
            f_recv_extr(3,l)=ny+1
            f_recv_extr(4,l)=ny+nbuff
            f_recv_extr(5,l)=1
            f_recv_extr(6,l)=nz
            
            fvec_send_extr(1,l)=nx-nbuffvec+1
            fvec_send_extr(2,l)=nx
            fvec_send_extr(3,l)=1
            fvec_send_extr(4,l)=nbuffvec
            fvec_send_extr(5,l)=1
            fvec_send_extr(6,l)=nz

            fvec_recv_extr(1,l)=1-nbuffvec
            fvec_recv_extr(2,l)=0
            fvec_recv_extr(3,l)=ny+1
            fvec_recv_extr(4,l)=ny+nbuffvec
            fvec_recv_extr(5,l)=1
            fvec_recv_extr(6,l)=nz

            b_send_extr(1,l)=nx-nbuffbvec+1
            b_send_extr(2,l)=nx
            b_send_extr(3,l)=1
            b_send_extr(4,l)=nbuffbvec
            b_send_extr(5,l)=1
            b_send_extr(6,l)=nz

            b_recv_extr(1,l)=1-nbuffbvec
            b_recv_extr(2,l)=0
            b_recv_extr(3,l)=ny+1
            b_recv_extr(4,l)=ny+nbuffbvec
            b_recv_extr(5,l)=1
            b_recv_extr(6,l)=nz
         endif
         if(exmpi(l)==-1 .and. eympi(l)==1)then

            f_send_extr(1,l)=1
            f_send_extr(2,l)=nbuff
            f_send_extr(3,l)=ny-nbuff+1
            f_send_extr(4,l)=ny
            f_send_extr(5,l)=1
            f_send_extr(6,l)=nz

            f_recv_extr(1,l)=nx+1
            f_recv_extr(2,l)=nx+nbuff
            f_recv_extr(3,l)=1-nbuff
            f_recv_extr(4,l)=0
            f_recv_extr(5,l)=1
            f_recv_extr(6,l)=nz
            
            fvec_send_extr(1,l)=1
            fvec_send_extr(2,l)=nbuffvec
            fvec_send_extr(3,l)=ny-nbuffvec+1
            fvec_send_extr(4,l)=ny
            fvec_send_extr(5,l)=1
            fvec_send_extr(6,l)=nz

            fvec_recv_extr(1,l)=nx+1
            fvec_recv_extr(2,l)=nx+nbuffvec
            fvec_recv_extr(3,l)=1-nbuffvec
            fvec_recv_extr(4,l)=0
            fvec_recv_extr(5,l)=1
            fvec_recv_extr(6,l)=nz

            b_send_extr(1,l)=1
            b_send_extr(2,l)=nbuffbvec
            b_send_extr(3,l)=ny-nbuffbvec+1
            b_send_extr(4,l)=ny
            b_send_extr(5,l)=1
            b_send_extr(6,l)=nz

            b_recv_extr(1,l)=nx+1
            b_recv_extr(2,l)=nx+nbuffbvec
            b_recv_extr(3,l)=1-nbuffbvec
            b_recv_extr(4,l)=0
            b_recv_extr(5,l)=1
            b_recv_extr(6,l)=nz
         endif
         !!!!   x   z
         if(exmpi(l)==1 .and. ezmpi(l)==1)then

            f_send_extr(1,l)=nx-nbuff+1
            f_send_extr(2,l)=nx
            f_send_extr(3,l)=1
            f_send_extr(4,l)=ny
            f_send_extr(5,l)=nz-nbuff+1
            f_send_extr(6,l)=nz

            f_recv_extr(1,l)=1-nbuff
            f_recv_extr(2,l)=0
            f_recv_extr(3,l)=1
            f_recv_extr(4,l)=ny
            f_recv_extr(5,l)=1-nbuff
            f_recv_extr(6,l)=0
            
            fvec_send_extr(1,l)=nx-nbuffvec+1
            fvec_send_extr(2,l)=nx
            fvec_send_extr(3,l)=1
            fvec_send_extr(4,l)=ny
            fvec_send_extr(5,l)=nz-nbuffvec+1
            fvec_send_extr(6,l)=nz

            fvec_recv_extr(1,l)=1-nbuffvec
            fvec_recv_extr(2,l)=0
            fvec_recv_extr(3,l)=1
            fvec_recv_extr(4,l)=ny
            fvec_recv_extr(5,l)=1-nbuffvec
            fvec_recv_extr(6,l)=0

            b_send_extr(1,l)=nx-nbuffbvec+1
            b_send_extr(2,l)=nx
            b_send_extr(3,l)=1
            b_send_extr(4,l)=ny
            b_send_extr(5,l)=nz-nbuffbvec+1
            b_send_extr(6,l)=nz

            b_recv_extr(1,l)=1-nbuffbvec
            b_recv_extr(2,l)=0
            b_recv_extr(3,l)=1
            b_recv_extr(4,l)=ny
            b_recv_extr(5,l)=1-nbuffbvec
            b_recv_extr(6,l)=0
         endif
         if(exmpi(l)==-1 .and. ezmpi(l)==-1)then

            f_send_extr(1,l)=1
            f_send_extr(2,l)=nbuff
            f_send_extr(3,l)=1
            f_send_extr(4,l)=ny
            f_send_extr(5,l)=1
            f_send_extr(6,l)=nbuff

            f_recv_extr(1,l)=nx+1
            f_recv_extr(2,l)=nx+nbuff
            f_recv_extr(3,l)=1
            f_recv_extr(4,l)=ny
            f_recv_extr(5,l)=nz+1
            f_recv_extr(6,l)=nz+nbuff
            
            fvec_send_extr(1,l)=1
            fvec_send_extr(2,l)=nbuffvec
            fvec_send_extr(3,l)=1
            fvec_send_extr(4,l)=ny
            fvec_send_extr(5,l)=1
            fvec_send_extr(6,l)=nbuffvec

            fvec_recv_extr(1,l)=nx+1
            fvec_recv_extr(2,l)=nx+nbuffvec
            fvec_recv_extr(3,l)=1
            fvec_recv_extr(4,l)=ny
            fvec_recv_extr(5,l)=nz+1
            fvec_recv_extr(6,l)=nz+nbuffvec

            b_send_extr(1,l)=1
            b_send_extr(2,l)=nbuffbvec
            b_send_extr(3,l)=1
            b_send_extr(4,l)=ny
            b_send_extr(5,l)=1
            b_send_extr(6,l)=nbuffbvec

            b_recv_extr(1,l)=nx+1
            b_recv_extr(2,l)=nx+nbuffbvec
            b_recv_extr(3,l)=1
            b_recv_extr(4,l)=ny
            b_recv_extr(5,l)=nz+1
            b_recv_extr(6,l)=nz+nbuffbvec
         endif
         if(exmpi(l)==1 .and. ezmpi(l)==-1)then

            f_send_extr(1,l)=nx-nbuff+1
            f_send_extr(2,l)=nx
            f_send_extr(3,l)=1
            f_send_extr(4,l)=ny
            f_send_extr(5,l)=1
            f_send_extr(6,l)=nbuff

            f_recv_extr(1,l)=1-nbuff
            f_recv_extr(2,l)=0
            f_recv_extr(3,l)=1
            f_recv_extr(4,l)=ny
            f_recv_extr(5,l)=nz+1
            f_recv_extr(6,l)=nz+nbuff
            
            fvec_send_extr(1,l)=nx-nbuffvec+1
            fvec_send_extr(2,l)=nx
            fvec_send_extr(3,l)=1
            fvec_send_extr(4,l)=ny
            fvec_send_extr(5,l)=1
            fvec_send_extr(6,l)=nbuffvec

            fvec_recv_extr(1,l)=1-nbuffvec
            fvec_recv_extr(2,l)=0
            fvec_recv_extr(3,l)=1
            fvec_recv_extr(4,l)=ny
            fvec_recv_extr(5,l)=nz+1
            fvec_recv_extr(6,l)=nz+nbuffvec

            b_send_extr(1,l)=nx-nbuffbvec+1
            b_send_extr(2,l)=nx
            b_send_extr(3,l)=1
            b_send_extr(4,l)=ny
            b_send_extr(5,l)=1
            b_send_extr(6,l)=nbuffbvec

            b_recv_extr(1,l)=1-nbuffbvec
            b_recv_extr(2,l)=0
            b_recv_extr(3,l)=1
            b_recv_extr(4,l)=ny
            b_recv_extr(5,l)=nz+1
            b_recv_extr(6,l)=nz+nbuffbvec
         endif
         if(exmpi(l)==-1 .and. ezmpi(l)==1)then

            f_send_extr(1,l)=1
            f_send_extr(2,l)=nbuff
            f_send_extr(3,l)=1
            f_send_extr(4,l)=ny
            f_send_extr(5,l)=nz-nbuff+1
            f_send_extr(6,l)=nz

            f_recv_extr(1,l)=nx+1
            f_recv_extr(2,l)=nx+nbuff
            f_recv_extr(3,l)=1
            f_recv_extr(4,l)=ny
            f_recv_extr(5,l)=1-nbuff
            f_recv_extr(6,l)=0
            
            fvec_send_extr(1,l)=1
            fvec_send_extr(2,l)=nbuffvec
            fvec_send_extr(3,l)=1
            fvec_send_extr(4,l)=ny
            fvec_send_extr(5,l)=nz-nbuffvec+1
            fvec_send_extr(6,l)=nz

            fvec_recv_extr(1,l)=nx+1
            fvec_recv_extr(2,l)=nx+nbuffvec
            fvec_recv_extr(3,l)=1
            fvec_recv_extr(4,l)=ny
            fvec_recv_extr(5,l)=1-nbuffvec
            fvec_recv_extr(6,l)=0

            b_send_extr(1,l)=1
            b_send_extr(2,l)=nbuffbvec
            b_send_extr(3,l)=1
            b_send_extr(4,l)=ny
            b_send_extr(5,l)=nz-nbuffbvec+1
            b_send_extr(6,l)=nz

            b_recv_extr(1,l)=nx+1
            b_recv_extr(2,l)=nx+nbuffbvec
            b_recv_extr(3,l)=1
            b_recv_extr(4,l)=ny
            b_recv_extr(5,l)=1-nbuffbvec
            b_recv_extr(6,l)=0
         endif
         !!!!   y   z
         if(eympi(l)==1 .and. ezmpi(l)==1)then

            f_send_extr(1,l)=1
            f_send_extr(2,l)=nx
            f_send_extr(3,l)=ny-nbuff+1
            f_send_extr(4,l)=ny
            f_send_extr(5,l)=nz-nbuff+1
            f_send_extr(6,l)=nz

            f_recv_extr(1,l)=1
            f_recv_extr(2,l)=nx
            f_recv_extr(3,l)=1-nbuff
            f_recv_extr(4,l)=0
            f_recv_extr(5,l)=1-nbuff
            f_recv_extr(6,l)=0
            
            fvec_send_extr(1,l)=1
            fvec_send_extr(2,l)=nx
            fvec_send_extr(3,l)=ny-nbuffvec+1
            fvec_send_extr(4,l)=ny
            fvec_send_extr(5,l)=nz-nbuffvec+1
            fvec_send_extr(6,l)=nz

            fvec_recv_extr(1,l)=1
            fvec_recv_extr(2,l)=nx
            fvec_recv_extr(3,l)=1-nbuffvec
            fvec_recv_extr(4,l)=0
            fvec_recv_extr(5,l)=1-nbuffvec
            fvec_recv_extr(6,l)=0

            b_send_extr(1,l)=1
            b_send_extr(2,l)=nx
            b_send_extr(3,l)=ny-nbuffbvec+1
            b_send_extr(4,l)=ny
            b_send_extr(5,l)=nz-nbuffbvec+1
            b_send_extr(6,l)=nz

            b_recv_extr(1,l)=1
            b_recv_extr(2,l)=nx
            b_recv_extr(3,l)=1-nbuffbvec
            b_recv_extr(4,l)=0
            b_recv_extr(5,l)=1-nbuffbvec
            b_recv_extr(6,l)=0
         endif
         if(eympi(l)==-1 .and. ezmpi(l)==-1)then

            f_send_extr(1,l)=1
            f_send_extr(2,l)=nx
            f_send_extr(3,l)=1
            f_send_extr(4,l)=nbuff
            f_send_extr(5,l)=1
            f_send_extr(6,l)=nbuff

            f_recv_extr(1,l)=1
            f_recv_extr(2,l)=nx
            f_recv_extr(3,l)=ny+1
            f_recv_extr(4,l)=ny+nbuff
            f_recv_extr(5,l)=nz+1
            f_recv_extr(6,l)=nz+nbuff
            
            fvec_send_extr(1,l)=1
            fvec_send_extr(2,l)=nx
            fvec_send_extr(3,l)=1
            fvec_send_extr(4,l)=nbuffvec
            fvec_send_extr(5,l)=1
            fvec_send_extr(6,l)=nbuffvec

            fvec_recv_extr(1,l)=1
            fvec_recv_extr(2,l)=nx
            fvec_recv_extr(3,l)=ny+1
            fvec_recv_extr(4,l)=ny+nbuffvec
            fvec_recv_extr(5,l)=nz+1
            fvec_recv_extr(6,l)=nz+nbuffvec

            b_send_extr(1,l)=1
            b_send_extr(2,l)=nx
            b_send_extr(3,l)=1
            b_send_extr(4,l)=nbuffbvec
            b_send_extr(5,l)=1
            b_send_extr(6,l)=nbuffbvec

            b_recv_extr(1,l)=1
            b_recv_extr(2,l)=nx
            b_recv_extr(3,l)=ny+1
            b_recv_extr(4,l)=ny+nbuffbvec
            b_recv_extr(5,l)=nz+1
            b_recv_extr(6,l)=nz+nbuffbvec
         endif
         if(eympi(l)==1 .and. ezmpi(l)==-1)then

            f_send_extr(1,l)=1
            f_send_extr(2,l)=nx
            f_send_extr(3,l)=ny-nbuff+1
            f_send_extr(4,l)=ny
            f_send_extr(5,l)=1
            f_send_extr(6,l)=nbuff

            f_recv_extr(1,l)=1
            f_recv_extr(2,l)=nx
            f_recv_extr(3,l)=1-nbuff
            f_recv_extr(4,l)=0
            f_recv_extr(5,l)=nz+1
            f_recv_extr(6,l)=nz+nbuff
            
            fvec_send_extr(1,l)=1
            fvec_send_extr(2,l)=nx
            fvec_send_extr(3,l)=ny-nbuffvec+1
            fvec_send_extr(4,l)=ny
            fvec_send_extr(5,l)=1
            fvec_send_extr(6,l)=nbuffvec

            fvec_recv_extr(1,l)=1
            fvec_recv_extr(2,l)=nx
            fvec_recv_extr(3,l)=1-nbuffvec
            fvec_recv_extr(4,l)=0
            fvec_recv_extr(5,l)=nz+1
            fvec_recv_extr(6,l)=nz+nbuffvec

            b_send_extr(1,l)=1
            b_send_extr(2,l)=nx
            b_send_extr(3,l)=ny-nbuffbvec+1
            b_send_extr(4,l)=ny
            b_send_extr(5,l)=1
            b_send_extr(6,l)=nbuffbvec

            b_recv_extr(1,l)=1
            b_recv_extr(2,l)=nx
            b_recv_extr(3,l)=1-nbuffbvec
            b_recv_extr(4,l)=0
            b_recv_extr(5,l)=nz+1
            b_recv_extr(6,l)=nz+nbuffbvec
         endif
         if(eympi(l)==-1 .and. ezmpi(l)==1)then

            f_send_extr(1,l)=1
            f_send_extr(2,l)=nx
            f_send_extr(3,l)=1
            f_send_extr(4,l)=nbuff
            f_send_extr(5,l)=nz-nbuff+1
            f_send_extr(6,l)=nz

            f_recv_extr(1,l)=1
            f_recv_extr(2,l)=nx
            f_recv_extr(3,l)=ny+1
            f_recv_extr(4,l)=ny+nbuff
            f_recv_extr(5,l)=1-nbuff
            f_recv_extr(6,l)=0
            
            fvec_send_extr(1,l)=1
            fvec_send_extr(2,l)=nx
            fvec_send_extr(3,l)=1
            fvec_send_extr(4,l)=nbuffvec
            fvec_send_extr(5,l)=nz-nbuffvec+1
            fvec_send_extr(6,l)=nz

            fvec_recv_extr(1,l)=1
            fvec_recv_extr(2,l)=nx
            fvec_recv_extr(3,l)=ny+1
            fvec_recv_extr(4,l)=ny+nbuffvec
            fvec_recv_extr(5,l)=1-nbuffvec
            fvec_recv_extr(6,l)=0

            b_send_extr(1,l)=1
            b_send_extr(2,l)=nx
            b_send_extr(3,l)=1
            b_send_extr(4,l)=nbuffbvec
            b_send_extr(5,l)=nz-nbuffbvec+1
            b_send_extr(6,l)=nz

            b_recv_extr(1,l)=1
            b_recv_extr(2,l)=nx
            b_recv_extr(3,l)=ny+1
            b_recv_extr(4,l)=ny+nbuffbvec
            b_recv_extr(5,l)=1-nbuffbvec
            b_recv_extr(6,l)=0
         endif

      enddo
      !corner
      do l=19,nlinksmpi
         if(exmpi(l)==1)then

            f_send_extr(1,l)=nx-nbuff+1
            f_send_extr(2,l)=nx

            f_recv_extr(1,l)=1-nbuff
            f_recv_extr(2,l)=0
            
            fvec_send_extr(1,l)=nx-nbuffvec+1
            fvec_send_extr(2,l)=nx

            fvec_recv_extr(1,l)=1-nbuffvec
            fvec_recv_extr(2,l)=0

            b_send_extr(1,l)=nx-nbuffbvec+1
            b_send_extr(2,l)=nx

            b_recv_extr(1,l)=1-nbuffbvec
            b_recv_extr(2,l)=0
         else

            f_send_extr(1,l)=1
            f_send_extr(2,l)=nbuff

            f_recv_extr(1,l)=nx+1
            f_recv_extr(2,l)=nx+nbuff
            
            fvec_send_extr(1,l)=1
            fvec_send_extr(2,l)=nbuffvec

            fvec_recv_extr(1,l)=nx+1
            fvec_recv_extr(2,l)=nx+nbuffvec

            b_send_extr(1,l)=1
            b_send_extr(2,l)=nbuffbvec

            b_recv_extr(1,l)=nx+1
            b_recv_extr(2,l)=nx+nbuffbvec
         endif
         if(eympi(l)==1)then

            f_send_extr(3,l)=ny-nbuff+1
            f_send_extr(4,l)=ny

            f_recv_extr(3,l)=1-nbuff
            f_recv_extr(4,l)=0
            
            fvec_send_extr(3,l)=ny-nbuffvec+1
            fvec_send_extr(4,l)=ny

            fvec_recv_extr(3,l)=1-nbuffvec
            fvec_recv_extr(4,l)=0

            b_send_extr(3,l)=ny-nbuffbvec+1
            b_send_extr(4,l)=ny

            b_recv_extr(3,l)=1-nbuffbvec
            b_recv_extr(4,l)=0
         else

            f_send_extr(3,l)=1
            f_send_extr(4,l)=nbuff

            f_recv_extr(3,l)=ny+1
            f_recv_extr(4,l)=ny+nbuff
            
            fvec_send_extr(3,l)=1
            fvec_send_extr(4,l)=nbuffvec

            fvec_recv_extr(3,l)=ny+1
            fvec_recv_extr(4,l)=ny+nbuffvec

            b_send_extr(3,l)=1
            b_send_extr(4,l)=nbuffbvec

            b_recv_extr(3,l)=ny+1
            b_recv_extr(4,l)=ny+nbuffbvec
         endif
         if(ezmpi(l)==1)then

            f_send_extr(5,l)=nz-nbuff+1
            f_send_extr(6,l)=nz

            f_recv_extr(5,l)=1-nbuff
            f_recv_extr(6,l)=0
            
            fvec_send_extr(5,l)=nz-nbuffvec+1
            fvec_send_extr(6,l)=nz

            fvec_recv_extr(5,l)=1-nbuffvec
            fvec_recv_extr(6,l)=0

            b_send_extr(5,l)=nz-nbuffbvec+1
            b_send_extr(6,l)=nz

            b_recv_extr(5,l)=1-nbuffbvec
            b_recv_extr(6,l)=0
         else

            f_send_extr(5,l)=1
            f_send_extr(6,l)=nbuff

            f_recv_extr(5,l)=nz+1
            f_recv_extr(6,l)=nz+nbuff
            
            fvec_send_extr(5,l)=1
            fvec_send_extr(6,l)=nbuffvec

            fvec_recv_extr(5,l)=nz+1
            fvec_recv_extr(6,l)=nz+nbuffvec

            b_send_extr(5,l)=1
            b_send_extr(6,l)=nbuffbvec

            b_recv_extr(5,l)=nz+1
            b_recv_extr(6,l)=nz+nbuffbvec
         endif
      enddo

      !gli estremi di c sono gli stessi di b
      c_send_extr(1:6,1:nlinksmpi)=b_send_extr(1:6,1:nlinksmpi)
      c_recv_extr(1:6,1:nlinksmpi)=b_recv_extr(1:6,1:nlinksmpi)
      
	  !gli estremi di bvec sono diversi da fvec perche seguono nbuffbvec

      !calcolo le quantita complessive da movimentare per ogni direzione l
      do l=1,nlinksmpi
         i_num_extr(l)=(f_recv_extr(2,l)-f_recv_extr(1,l)+1)* &
            (f_recv_extr(4,l)-f_recv_extr(3,l)+1)* &
            (f_recv_extr(6,l)-f_recv_extr(5,l)+1)
         f_num_extr(l)=(f_recv_extr(2,l)-f_recv_extr(1,l)+1)* &
            (f_recv_extr(4,l)-f_recv_extr(3,l)+1)* &
            (f_recv_extr(6,l)-f_recv_extr(5,l)+1)*num_phifields_datampi
         fvec_num_extr(l)=(fvec_recv_extr(2,l)-fvec_recv_extr(1,l)+1)* &
            (fvec_recv_extr(4,l)-fvec_recv_extr(3,l)+1)* &
            (fvec_recv_extr(6,l)-fvec_recv_extr(5,l)+1)*num_auxfields_datampi
         b_num_extr(l)=(b_recv_extr(2,l)-b_recv_extr(1,l)+1)* &
            (b_recv_extr(4,l)-b_recv_extr(3,l)+1)* &
            (b_recv_extr(6,l)-b_recv_extr(5,l)+1)*num_hfields_datampi
         c_num_extr(l)=(c_recv_extr(2,l)-c_recv_extr(1,l)+1)* &
            (c_recv_extr(4,l)-c_recv_extr(3,l)+1)* &
            (c_recv_extr(6,l)-c_recv_extr(5,l)+1)*num_forces_datampi
      enddo



#ifdef VERBOSE
      !stampo per debug
      if(ltwocomp)then
        if(myrank==0)write(6,'(a)')'####################### f_send_extr  f_recv_extr #######################'
        do l=1,nlinksmpi
           if(myrank==0)write(6,'(a,i3,a,3i3,a,6i4,a,6i4,a,i4)')'dir l ',l,' disp ',&
              exmpi(l),eympi(l),ezmpi(l),' extremes d',f_send_extr(1:6,l),&
              ' s ',f_recv_extr(1:6,l),' num ',f_num_extr(l)
           call flush(6)
#ifdef MPI
           call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif
        enddo
#ifdef MPI
        call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif
        if(myrank==0)write(6,'(a)')'####################### fvec_send_extr  fvec_recv_extr #######################'
        do l=1,nlinksmpi
           if(myrank==0)write(6,'(a,i3,a,3i3,a,6i4,a,6i4,a,i4)')'dir l ',l,' disp ',&
              exmpi(l),eympi(l),ezmpi(l),' extremes d',f_send_extr(1:6,l),&
              ' s ',f_recv_extr(1:6,l),' num ',fvec_num_extr(l)
           call flush(6)
#ifdef MPI
           call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif
        enddo
#ifdef MPI
        call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif
      endif
      if(myrank==0)write(6,'(a)')'####################### b_send_extr  b_recv_extr #######################'
         do l=1,nlinksmpi
            if(myrank==0)write(6,'(a,i3,a,3i3,a,6i4,a,6i4,a,i4)')'dir l ',l,' disp ',&
               exmpi(l),eympi(l),ezmpi(l),' extremes d',b_send_extr(1:6,l),&
               ' s ',b_recv_extr(1:6,l),' num ',b_num_extr(l)
            call flush(6)
#ifdef MPI
            call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif
         enddo
      
#ifdef MPI
      call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif

      if(myrank==0)write(6,'(a)')'####################### c_send_extr  c_recv_extr #######################'
         do l=1,nlinksmpi
            if(myrank==0)write(6,'(a,i3,a,3i3,a,6i4,a,6i4,a,i4)')'dir l ',l,' disp ',&
               exmpi(l),eympi(l),ezmpi(l),' extremes d',c_send_extr(1:6,l),&
               ' s ',c_recv_extr(1:6,l),' num ',c_num_extr(l)
            call flush(6)
#ifdef MPI
            call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif
         enddo
      
#ifdef MPI
      call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif

#endif

      !creo i tipi MPI contigui che mi servono per i send e receive
      !lo faccio su 13 direzioni perchè per l ed lopp le quantità da muovere sono uguali
#ifdef MPI

      !dir
      if(ltwocomp)then
        do l=1,nlinksmpi,2
          ll=(l+1)/2
          call MPI_type_contiguous(f_num_extr(l), STRMPIREAL, f_datampi(ll), ierr)
          call MPI_type_commit(f_datampi(ll),ierr)
#ifdef VERBOSE
          if(myrank.eq.0) then
            write(6,'(a,2i4,a,f16.8)') 'CREATE BUFFER: f_datampi',ll*2-1,ll*2,' (KB)-->', &
               real(f_num_extr(l),kind=db) *4 / 1024
            call flush(6)
          endif
#endif
        enddo

      !dir
        do l=1,nlinksmpi,2
          ll=(l+1)/2
          call MPI_type_contiguous(fvec_num_extr(l), STRMPIREAL, fvec_datampi(ll), ierr)
          call MPI_type_commit(fvec_datampi(ll),ierr)
#ifdef VERBOSE
          if(myrank.eq.0) then
            write(6,'(a,2i4,a,f16.8)') 'CREATE BUFFER: fvec_datampi',ll*2-1,ll*2,' (KB)-->', &
               real(fvec_num_extr(l),kind=db) *4 / 1024
            call flush(6)
          endif
#endif
        enddo
      endif
      
      !dir
      do l=1,nlinksmpi,2
         ll=(l+1)/2
         call MPI_type_contiguous(b_num_extr(l), STRMPIREAL, b_datampi(ll), ierr)
         call MPI_type_commit(b_datampi(ll),ierr)
#ifdef VERBOSE
         if(myrank.eq.0) then
            write(6,'(a,2i4,a,f16.8)') 'CREATE BUFFER: b_datampi',ll*2-1,ll*2,' (KB)-->', &
               real(b_num_extr(l),kind=db) *4 / 1024
            call flush(6)
         endif
#endif
      enddo

      !dir
      do l=1,nlinksmpi,2
         ll=(l+1)/2
         call MPI_type_contiguous(c_num_extr(l), STRMPIREAL, c_datampi(ll), ierr)
         call MPI_type_commit(c_datampi(ll),ierr)
#ifdef VERBOSE
         if(myrank.eq.0) then
            write(6,'(a,2i4,a,f16.8)') 'CREATE BUFFER: c_datampi',ll*2-1,ll*2,' (KB)-->', &
               real(c_num_extr(l),kind=db) *4 / 1024
            call flush(6)
         endif
#endif
      enddo   

      !dir
      do l=1,nlinksmpi,2
         ll=(l+1)/2
         call MPI_type_contiguous(i_num_extr(l), MYMPIINTS, i_datampi(ll), ierr)
         call MPI_type_commit(i_datampi(ll),ierr)
#ifdef VERBOSE
         if(myrank.eq.0) then
            write(6,'(a,2i4,a,f16.8)') 'CREATE BUFFER: i_datampi',ll*2-1,ll*2,' (KB)-->', &
               real(i_num_extr(l),kind=db) *1 / 1024
            call flush(6)
         endif
#endif
      enddo

      call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif

  
      



      !alloca i buffer per mandare e ricevere
      

      !alloco per f
    if(ltwocomp)then
      f_numtot_extr=sum(f_num_extr)
      ll=0
      do l=1,nlinksmpi
         f_nbuffmpi_send(l)=ll+1
         if(lsend_dir(l))ll=ll+f_num_extr(l)
      enddo
      allocate(f_send_buffmpi(ll))
      f_send_buffmpi=real(0.d0,kind=db)

      ll=0
      do l=1,nlinksmpi
         f_nbuffmpi_recv(l)=ll+1
         if(lrecv_dir(l))ll=ll+f_num_extr(l)
      enddo
      allocate(f_recv_buffmpi(ll))
      f_recv_buffmpi=real(0.d0,kind=db)

      !alloco per fvec
      fvec_numtot_extr=sum(fvec_num_extr)
      ll=0
      do l=1,nlinksmpi
         fvec_nbuffmpi_send(l)=ll+1
         if(lsend_dir(l))ll=ll+fvec_num_extr(l)
      enddo
      allocate(fvec_send_buffmpi(ll))
      fvec_send_buffmpi=real(0.d0,kind=db)

      ll=0
      do l=1,nlinksmpi
         fvec_nbuffmpi_recv(l)=ll+1
         if(lrecv_dir(l))ll=ll+fvec_num_extr(l)
      enddo
      allocate(fvec_recv_buffmpi(ll))
      fvec_recv_buffmpi=real(0.d0,kind=db)
    endif
      
      !alloco per b
      b_numtot_extr=sum(b_num_extr)
      ll=0
      do l=1,nlinksmpi
         b_nbuffmpi_send(l)=ll+1
         if(lsend_dir(l))ll=ll+b_num_extr(l)
      enddo
      allocate(b_send_buffmpi(ll))
      b_send_buffmpi=real(0.d0,kind=db)

      ll=0
      do l=1,nlinksmpi
         b_nbuffmpi_recv(l)=ll+1
         if(lrecv_dir(l))ll=ll+b_num_extr(l)
      enddo
      allocate(b_recv_buffmpi(ll))
      b_recv_buffmpi=real(0.d0,kind=db)
     
      !alloco per c
      c_numtot_extr=sum(c_num_extr)
      ll=0
      do l=1,nlinksmpi
         c_nbuffmpi_send(l)=ll+1
         if(lsend_dir(l))ll=ll+c_num_extr(l)
      enddo
      allocate(c_send_buffmpi(ll))
      c_send_buffmpi=real(0.d0,kind=db)

      ll=0
      do l=1,nlinksmpi
         c_nbuffmpi_recv(l)=ll+1
         if(lrecv_dir(l))ll=ll+c_num_extr(l)
      enddo
      allocate(c_recv_buffmpi(ll))
      c_recv_buffmpi=real(0.d0,kind=db)

      !alloco per i
      i_numtot_extr=sum(i_num_extr)
      ll=0
      do l=1,nlinksmpi
         i_nbuffmpi_send(l)=ll+1
         if(lsend_dir(l))ll=ll+i_num_extr(l)
      enddo
      allocate(i_send_buffmpi(ll))
      i_send_buffmpi=int(0,kind=isf)

      ll=0
      do l=1,nlinksmpi
         i_nbuffmpi_recv(l)=ll+1
         if(lrecv_dir(l))ll=ll+i_num_extr(l)
      enddo
      allocate(i_recv_buffmpi(ll))
      i_recv_buffmpi=int(0,kind=isf)

#ifdef MPI
      call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif

#ifdef MEM_CHECK
      if(myrank == 0) then
         mem_stop = get_mem();
         write(6,*) "MEM_CHECK: after sub. input mem =", mem_stop
      endif
#endif


   end subroutine setup_mpi
   
   subroutine setup_io_comm2d(ldowrite,l)
   
      implicit none
      
      logical, intent(in) :: ldowrite
      integer, intent(in) :: l
      
      integer :: color,ierr
#if defined(MPI)        
      io_comm2d(l) = MPI_COMM_NULL
      !qui seleziono solo i processi buoni, quelli che devono scrivere
      if (ldowrite)then
        color = 1  ! Gruppo che partecipa all'IO
      else
        color = MPI_UNDEFINED  ! Escludiamo gli altri processi
      endif
      
      ! Creiamo un nuovo comunicatore per i processi che scrivono
      call MPI_COMM_SPLIT(MPI_COMM_WORLD, color, myrank, io_comm2d(l), ierr)
#endif
      
   end subroutine setup_io_comm2d   
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!*******************************PHI********************************************************************!
   subroutine exchange_phifields_intpbc(phifields_s)

      implicit none
      
      real(kind=strdb), allocatable, dimension(:,:,:,:,:) :: phifields_s
      integer :: lmio
      integer :: oi,oj,ok
      integer :: ii,jj,kk
      integer :: oii,ojj,okk
      integer :: imin,imax,jmin,jmax,kmin,kmax
      integer :: oxblock,oyblock,ozblock,omyblock
      integer :: xblock,yblock,zblock,myblock

      do lmio=1,nlinksmpi
         if(.not. lintpbc_dir(lmio))cycle
         imin=f_recv_extr(1,lmio)
         imax=f_recv_extr(2,lmio)
         jmin=f_recv_extr(3,lmio)
         jmax=f_recv_extr(4,lmio)
         kmin=f_recv_extr(5,lmio)
         kmax=f_recv_extr(6,lmio)
#ifdef ACCNOKERNELS
         !$acc parallel loop independent collapse(3) present(intpbc_dir,phifields_s, &
         !$acc& ) private(i,j,k,ii,jj,kk,myblock, &
         !$acc& oi,oj,ok,oii,ojj,okk,omyblock)
#else
         !$acc kernels present(intpbc_dir,phifields_s)
         !$acc loop independent collapse(3) private(i,j,k,ii,jj,kk,myblock, &
         !$acc& oi,oj,ok,oii,ojj,okk,omyblock)
#endif
         do k=kmin,kmax
            do j=jmin,jmax
               do i=imin,imax
                  oi=i
                  oj=j
                  ok=k
                  if(intpbc_dir(1,lmio))oi=mod(oi+nx-1,nx)+1
                  if(intpbc_dir(2,lmio))oj=mod(oj+ny-1,ny)+1
                  if(intpbc_dir(3,lmio))ok=mod(ok+nz-1,nz)+1
                  
                  oxblock=(oi+2*TILE_DIMx-1)/TILE_DIMx   
                  oyblock=(oj+2*TILE_DIMy-1)/TILE_DIMy     
                  ozblock=(ok+2*TILE_DIMz-1)/TILE_DIMz 
                  omyblock=(oxblock-1)+(oyblock-1)*nxblock+(ozblock-1)*nxyblock+1
                  oii=oi-oxblock*TILE_DIMx+2*TILE_DIMx
                  ojj=oj-oyblock*TILE_DIMy+2*TILE_DIMy
                  okk=ok-ozblock*TILE_DIMz+2*TILE_DIMz
                  
                  xblock=(i+2*TILE_DIMx-1)/TILE_DIMx   
                  yblock=(j+2*TILE_DIMy-1)/TILE_DIMy     
                  zblock=(k+2*TILE_DIMz-1)/TILE_DIMz  
                  myblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
                  ii=i-xblock*TILE_DIMx+2*TILE_DIMx
                  jj=j-yblock*TILE_DIMy+2*TILE_DIMy
                  kk=k-zblock*TILE_DIMz+2*TILE_DIMz   
                  
                 
                  phifields_s(ii,jj,kk,1,myblock)= &
                   phifields_s(oii,ojj,okk,1,omyblock)
                                 
               enddo
            enddo
         enddo
#ifdef ACCNOKERNELS
         !$acc end parallel loop
#else
         !$acc end kernels
#endif
      enddo

   end subroutine exchange_phifields_intpbc

   subroutine exchange_phifields_sendrecv(phifields_s)

      implicit none
      
      real(kind=strdb), allocatable, dimension(:,:,:,:,:) :: phifields_s
      integer :: l,ll,myoffset,tag,ierr,mm
#ifndef MPI
      return
#endif
      do l=1,nlinksmpi
         if(lsend_dir(l))call packaging_phifields_buffmpi(l,phifields_s)
      enddo
      !$acc wait
      !f_recv_buffmpi=f_send_buffmpi
#ifdef MPI
      mm=0
      do l=1,nlinksmpi
         ll=(l+1)/2
         if(lsend_dir(l))then
            mm=mm+1
            myoffset=f_nbuffmpi_send(l)
            !$acc host_data use_device(f_send_buffmpi)
            call mpi_isend(f_send_buffmpi(myoffset),1,f_datampi(ll),send_dir(l), &
               f_mpitag(l),lbecomm,f_reqs(mm),ierr)
            !$acc end host_data
         endif
         if(lrecv_dir(l))then
            mm=mm+1
            myoffset=f_nbuffmpi_recv(l)
            !$acc host_data use_device(f_recv_buffmpi)
            call mpi_irecv(f_recv_buffmpi(myoffset),1,f_datampi(ll),recv_dir(l), &
               f_mpitag(l),lbecomm,f_reqs(mm),ierr)
            !$acc end host_data
         endif
      enddo
      nf_reqs=mm
#endif


   end subroutine exchange_phifields_sendrecv

   subroutine exchange_phifields_wait(phifields_s)

      implicit none
      
      real(kind=strdb), allocatable, dimension(:,:,:,:,:) :: phifields_s
      integer :: l,ll,myoffset,tag,ierr
#ifdef MPI
      integer, dimension(:), allocatable :: myierr
      integer, dimension(:,:), allocatable :: mystatus
#else
      return
#endif
#ifdef MPI

      allocate(myierr(nf_reqs))
      myierr=0

      allocate(mystatus(MPI_STATUS_SIZE,nf_reqs))
      !$acc wait
      call mpi_waitall(nf_reqs,f_reqs,mystatus,ierr)
      !$acc wait

      if(any(myierr.ne.0))call doerror(6,'ERROR in exchange_phifields_sendrecv')
#endif
      do l=1,nlinksmpi
         if(lrecv_dir(l))call depackaging_phifields_buffmpi(l,phifields_s)
      enddo

   end subroutine exchange_phifields_wait

   subroutine packaging_phifields_buffmpi(lmio,phifields_s)

      implicit none

      integer, intent(in) :: lmio
      real(kind=strdb), allocatable, dimension(:,:,:,:,:) :: phifields_s
      integer :: myoffset

      integer :: i,j,k,l,ll,m1,m2,m3
      integer :: imin,imax,jmin,jmax,kmin,kmax
      integer :: idx=0
      integer :: ii,jj,kk,xblock,yblock,zblock,myblock

      myoffset=f_nbuffmpi_send(lmio)
      m1=f_send_extr(2,lmio)-f_send_extr(1,lmio)+1
      m2=f_send_extr(4,lmio)-f_send_extr(3,lmio)+1
      m3=f_send_extr(6,lmio)-f_send_extr(5,lmio)+1
      !scorro sul numero di campi da prendere (per scalare = 1)
      ll=1
      imin=f_send_extr(1,lmio)
      imax=f_send_extr(2,lmio)
      jmin=f_send_extr(3,lmio)
      jmax=f_send_extr(4,lmio)
      kmin=f_send_extr(5,lmio)
      kmax=f_send_extr(6,lmio)
#ifdef ACCNOKERNELS
      !$acc parallel loop independent collapse(3) present(f_send_buffmpi,phifields_s,f_send_extr) &
      !$acc& private(idx,ii,jj,kk,xblock,yblock,zblock,myblock)
#else
      !$acc kernels present(f_send_buffmpi,phifields_s,f_send_extr)
      !$acc loop independent collapse(3)  private(i,j,k,idx,ii,jj,kk,xblock,yblock,zblock,myblock)
#endif
      do k=kmin,kmax
         do j=jmin,jmax
            do i=imin,imax
               xblock=(i+2*TILE_DIMx-1)/TILE_DIMx   
               yblock=(j+2*TILE_DIMy-1)/TILE_DIMy     
               zblock=(k+2*TILE_DIMz-1)/TILE_DIMz  
               myblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
               ii=i-xblock*TILE_DIMx+2*TILE_DIMx
               jj=j-yblock*TILE_DIMy+2*TILE_DIMy
               kk=k-zblock*TILE_DIMz+2*TILE_DIMz    
               !linearizzo con l'ordine naturale e metto nel buffer unico per tutte le direzioni
               !poi mandero solo i pezzi contigui che mi servono per la data direzione
               idx=myoffset+(i-f_send_extr(1,lmio))+(j-f_send_extr(3,lmio))*m1+(&
                  k-f_send_extr(5,lmio))*(m1*m2)+(ll-1)*(m1*m2*m3)
               f_send_buffmpi(idx)=phifields_s(ii,jj,kk,1,myblock)		   
            enddo
         enddo
      enddo
#ifdef ACCNOKERNELS
      !$acc end parallel loop
#else
      !$acc end kernels
#endif


   end subroutine packaging_phifields_buffmpi

   subroutine depackaging_phifields_buffmpi(lmio,phifields_s)

      implicit none

      integer, intent(in) :: lmio
      real(kind=strdb), allocatable, dimension(:,:,:,:,:) :: phifields_s
      integer :: myoffset

      integer :: i,j,k,l,ll,m1,m2,m3
      integer :: imin,imax,jmin,jmax,kmin,kmax
      integer :: idx=0
      integer :: ii,jj,kk,xblock,yblock,zblock,myblock

      myoffset=f_nbuffmpi_recv(lmio)
      m1=f_recv_extr(2,lmio)-f_recv_extr(1,lmio)+1
      m2=f_recv_extr(4,lmio)-f_recv_extr(3,lmio)+1
      m3=f_recv_extr(6,lmio)-f_recv_extr(5,lmio)+1
      !scorro sul numero di campi da prendere (per scalare = 1)
      ll=1
      imin=f_recv_extr(1,lmio)
      imax=f_recv_extr(2,lmio)
      jmin=f_recv_extr(3,lmio)
      jmax=f_recv_extr(4,lmio)
      kmin=f_recv_extr(5,lmio)
      kmax=f_recv_extr(6,lmio)
#ifdef ACCNOKERNELS
      !$acc parallel loop independent collapse(3) present(f_recv_buffmpi,phifields_s,f_recv_extr) &
      !$acc& private(idx,ii,jj,kk,xblock,yblock,zblock,myblock)
#else
      !$acc kernels present(f_recv_buffmpi,phifields_s,f_recv_extr)
      !$acc loop independent collapse(3)  private(i,j,k,idx,ii,jj,kk,xblock,yblock,zblock,myblock)
#endif
      do k=kmin,kmax
         do j=jmin,jmax
            do i=imin,imax
               xblock=(i+2*TILE_DIMx-1)/TILE_DIMx   
               yblock=(j+2*TILE_DIMy-1)/TILE_DIMy     
               zblock=(k+2*TILE_DIMz-1)/TILE_DIMz  
               myblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
               ii=i-xblock*TILE_DIMx+2*TILE_DIMx
               jj=j-yblock*TILE_DIMy+2*TILE_DIMy
               kk=k-zblock*TILE_DIMz+2*TILE_DIMz 
               !linearizzo con l'ordine naturale e metto nel buffer unico per tutte le direzioni
               !poi mandero solo i pezzi contigui che mi servono per la data direzione
               idx=myoffset+(i-f_recv_extr(1,lmio))+(j-f_recv_extr(3,lmio))*m1+(&
                  k-f_recv_extr(5,lmio))*(m1*m2)+(ll-1)*(m1*m2*m3)
               phifields_s(ii,jj,kk,1,myblock)=f_recv_buffmpi(idx)
            enddo
         enddo
      enddo
#ifdef ACCNOKERNELS
      !$acc end parallel loop
#else
      !$acc end kernels
#endif


   end subroutine depackaging_phifields_buffmpi

   !****************************** normx normy normz********************************************************************!
   subroutine exchange_auxfields_intpbc

      implicit none

      integer :: lmio,oi,oj,ok
      integer :: imin,imax,jmin,jmax,kmin,kmax
      integer :: ii,jj,kk
      integer :: oii,ojj,okk
      integer :: oxblock,oyblock,ozblock,omyblock
      integer :: xblock,yblock,zblock,myblock
      
      if(nauxfields==0)return

      do lmio=1,nlinksmpi
         if(.not. lintpbc_dir(lmio))cycle
         imin=fvec_recv_extr(1,lmio)
         imax=fvec_recv_extr(2,lmio)
         jmin=fvec_recv_extr(3,lmio)
         jmax=fvec_recv_extr(4,lmio)
         kmin=fvec_recv_extr(5,lmio)
         kmax=fvec_recv_extr(6,lmio)
#ifdef ACCNOKERNELS
         !$acc parallel loop independent collapse(3) present(intpbc_dir,auxfields) &
         !$acc& private(ii,jj,kk,oi,oj,ok,oii,ojj,okk,oxblock,oyblock,ozblock, &
         !$acc& omyblock,xblock,yblock,zblock,myblock)
#else
         !$acc kernels present(intpbc_dir,auxfields)
         !$acc loop independent collapse(3)  private(i,j,k,ii,jj,kk, &
         !$acc& oi,oj,ok,oii,ojj,okk,oxblock,oyblock,ozblock, &
         !$acc& omyblock,xblock,yblock,zblock,myblock)
#endif
         do k=kmin,kmax
            do j=jmin,jmax
               do i=imin,imax
                  oi=i
                  oj=j
                  ok=k
                  if(intpbc_dir(1,lmio))oi=mod(oi+nx-1,nx)+1
                  if(intpbc_dir(2,lmio))oj=mod(oj+ny-1,ny)+1
                  if(intpbc_dir(3,lmio))ok=mod(ok+nz-1,nz)+1
                  
                  oxblock=(oi+2*TILE_DIMx-1)/TILE_DIMx   
                  oyblock=(oj+2*TILE_DIMy-1)/TILE_DIMy     
                  ozblock=(ok+2*TILE_DIMz-1)/TILE_DIMz 
                  omyblock=(oxblock-1)+(oyblock-1)*nxblock+(ozblock-1)*nxyblock+1
                  oii=oi-oxblock*TILE_DIMx+2*TILE_DIMx
                  ojj=oj-oyblock*TILE_DIMy+2*TILE_DIMy
                  okk=ok-ozblock*TILE_DIMz+2*TILE_DIMz
                  
                  xblock=(i+2*TILE_DIMx-1)/TILE_DIMx   
                  yblock=(j+2*TILE_DIMy-1)/TILE_DIMy     
                  zblock=(k+2*TILE_DIMz-1)/TILE_DIMz  
                  myblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
                  ii=i-xblock*TILE_DIMx+2*TILE_DIMx
                  jj=j-yblock*TILE_DIMy+2*TILE_DIMy
                  kk=k-zblock*TILE_DIMz+2*TILE_DIMz 
                   
#ifdef TWOCOMPONENT
                  auxfields(ii,jj,kk,1,myblock)= &
                   auxfields(oii,ojj,okk,1,omyblock)
                   
                  auxfields(ii,jj,kk,2,myblock)= &
                   auxfields(oii,ojj,okk,2,omyblock)
                   
                  auxfields(ii,jj,kk,3,myblock)= &
                   auxfields(oii,ojj,okk,3,omyblock)
                   
                  auxfields(ii,jj,kk,4,myblock)= &
                   auxfields(oii,ojj,okk,4,omyblock)
                   
                  auxfields(ii,jj,kk,5,myblock)= &
                   auxfields(oii,ojj,okk,5,omyblock)
                  
                  auxfields(ii,jj,kk,6,myblock)= &
                   auxfields(oii,ojj,okk,6,omyblock)
                  
                  auxfields(ii,jj,kk,7,myblock)= &
                   auxfields(oii,ojj,okk,7,omyblock)
                    
#endif
               enddo
            enddo
         enddo
#ifdef ACCNOKERNELS
         !$acc end parallel loop
#else
         !$acc end kernels
#endif
      enddo

   end subroutine exchange_auxfields_intpbc

   subroutine exchange_auxfields_sendrecv

      implicit none

      integer :: l,ll,myoffset,tag,ierr,mm
      
#ifndef MPI
      return
#endif
      if(nauxfields==0)return
      
      do l=1,nlinksmpi
         if(lsend_dir(l))call packaging_auxfields_buffmpi(l)
      enddo
      !$acc wait
#ifdef MPI
      mm=0
      do l=1,nlinksmpi
         ll=(l+1)/2
         if(lsend_dir(l))then
            mm=mm+1
            myoffset=fvec_nbuffmpi_send(l)
            !$acc host_data use_device(fvec_send_buffmpi)
            call mpi_isend(fvec_send_buffmpi(myoffset),1,fvec_datampi(ll),send_dir(l), &
               fvec_mpitag(l),lbecomm,fvec_reqs(mm),ierr)
            !$acc end host_data
         endif
         if(lrecv_dir(l))then
            mm=mm+1
            myoffset=fvec_nbuffmpi_recv(l)
            !$acc host_data use_device(fvec_recv_buffmpi)
            call mpi_irecv(fvec_recv_buffmpi(myoffset),1,fvec_datampi(ll),recv_dir(l), &
               fvec_mpitag(l),lbecomm,fvec_reqs(mm),ierr)
            !$acc end host_data
         endif
      enddo
      nfvec_reqs=mm
#endif


   end subroutine exchange_auxfields_sendrecv

   subroutine exchange_auxfields_wait

      implicit none

      integer :: l,ll,myoffset,tag,ierr
#ifdef MPI
      integer, dimension(:), allocatable :: myierr
      integer, dimension(:,:), allocatable :: mystatus
#else
      return
#endif
      if(nauxfields==0)return
#ifdef MPI

      allocate(myierr(nfvec_reqs))
      myierr=0

      allocate(mystatus(MPI_STATUS_SIZE,nfvec_reqs))
      !$acc wait
      call mpi_waitall(nfvec_reqs,fvec_reqs,mystatus,ierr)
      !$acc wait

      if(any(myierr.ne.0))call doerror(6,'ERROR in exchange_auxfields_wait')
#endif
      do l=1,nlinksmpi
         if(lrecv_dir(l))call depackaging_auxfields_buffmpi(l)
      enddo

   end subroutine exchange_auxfields_wait

   subroutine packaging_auxfields_buffmpi(lmio)

      implicit none

      integer, intent(in) :: lmio
      integer :: myoffset

      integer :: i,j,k,l,ll=0,m1,m2,m3
      integer :: imin,imax,jmin,jmax,kmin,kmax
      integer :: idx=0
      integer :: ii,jj,kk,xblock,yblock,zblock,myblock

      myoffset=fvec_nbuffmpi_send(lmio)
      m1=fvec_send_extr(2,lmio)-fvec_send_extr(1,lmio)+1
      m2=fvec_send_extr(4,lmio)-fvec_send_extr(3,lmio)+1
      m3=fvec_send_extr(6,lmio)-fvec_send_extr(5,lmio)+1
      imin=fvec_send_extr(1,lmio)
      imax=fvec_send_extr(2,lmio)
      jmin=fvec_send_extr(3,lmio)
      jmax=fvec_send_extr(4,lmio)
      kmin=fvec_send_extr(5,lmio)
      kmax=fvec_send_extr(6,lmio)
#ifdef ACCNOKERNELS
      !$acc parallel loop independent collapse(3) present(fvec_send_buffmpi,fvec_send_extr,auxfields) &
      !$acc& private(ll,idx,ii,jj,kk,xblock,yblock,zblock,myblock)
#else
      !$acc kernels present(fvec_send_buffmpi,fvec_send_extr,auxfields)
      !$acc loop independent collapse(3) private(i,j,k,ll,idx,ii,jj,kk,xblock,yblock,zblock,myblock)
#endif
      do k=kmin,kmax
         do j=jmin,jmax
            do i=imin,imax
            
               xblock=(i+2*TILE_DIMx-1)/TILE_DIMx   
               yblock=(j+2*TILE_DIMy-1)/TILE_DIMy     
               zblock=(k+2*TILE_DIMz-1)/TILE_DIMz  
               myblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
               ii=i-xblock*TILE_DIMx+2*TILE_DIMx
               jj=j-yblock*TILE_DIMy+2*TILE_DIMy
               kk=k-zblock*TILE_DIMz+2*TILE_DIMz   
            
               !linearizzo con l'ordine naturale e metto nel buffer unico per tutte le direzioni
               !poi mandero solo i pezzi contigui che mi servono per la data direzione

               !scorro sul numero di campi da prendere il massimo previsto è num_auxfields_datampi

#ifdef TWOCOMPONENT

               ll=1
               idx=myoffset+(i-fvec_send_extr(1,lmio))+(j-fvec_send_extr(3,lmio))*m1+(&
                  k-fvec_send_extr(5,lmio))*(m1*m2)+(ll-1)*(m1*m2*m3)
               fvec_send_buffmpi(idx)=auxfields(ii,jj,kk,1,myblock)

               ll=2
               idx=myoffset+(i-fvec_send_extr(1,lmio))+(j-fvec_send_extr(3,lmio))*m1+(&
                  k-fvec_send_extr(5,lmio))*(m1*m2)+(ll-1)*(m1*m2*m3)
               fvec_send_buffmpi(idx)=auxfields(ii,jj,kk,2,myblock)

               ll=3
               idx=myoffset+(i-fvec_send_extr(1,lmio))+(j-fvec_send_extr(3,lmio))*m1+(&
                  k-fvec_send_extr(5,lmio))*(m1*m2)+(ll-1)*(m1*m2*m3)
               fvec_send_buffmpi(idx)=auxfields(ii,jj,kk,3,myblock)
			   
               ll=4
               idx=myoffset+(i-fvec_send_extr(1,lmio))+(j-fvec_send_extr(3,lmio))*m1+(&
                  k-fvec_send_extr(5,lmio))*(m1*m2)+(ll-1)*(m1*m2*m3)
               fvec_send_buffmpi(idx)=auxfields(ii,jj,kk,4,myblock)
			   
			   ll=5
               idx=myoffset+(i-fvec_send_extr(1,lmio))+(j-fvec_send_extr(3,lmio))*m1+(&
                  k-fvec_send_extr(5,lmio))*(m1*m2)+(ll-1)*(m1*m2*m3)
               fvec_send_buffmpi(idx)=auxfields(ii,jj,kk,5,myblock)
			   
			   ll=6
               idx=myoffset+(i-fvec_send_extr(1,lmio))+(j-fvec_send_extr(3,lmio))*m1+(&
                  k-fvec_send_extr(5,lmio))*(m1*m2)+(ll-1)*(m1*m2*m3)
               fvec_send_buffmpi(idx)=auxfields(ii,jj,kk,6,myblock)
			   
			   ll=7
               idx=myoffset+(i-fvec_send_extr(1,lmio))+(j-fvec_send_extr(3,lmio))*m1+(&
                  k-fvec_send_extr(5,lmio))*(m1*m2)+(ll-1)*(m1*m2*m3)
               fvec_send_buffmpi(idx)=auxfields(ii,jj,kk,7,myblock)

#endif
            enddo
         enddo
      enddo
#ifdef ACCNOKERNELS
      !$acc end parallel loop
#else
      !$acc end kernels
#endif

   end subroutine packaging_auxfields_buffmpi

   subroutine depackaging_auxfields_buffmpi(lmio)

      implicit none

      integer, intent(in) :: lmio
      integer :: myoffset

      integer :: i,j,k,l,ll=0,m1,m2,m3
      integer :: imin,imax,jmin,jmax,kmin,kmax
      integer :: idx=0
      integer :: ii,jj,kk,xblock,yblock,zblock,myblock

      myoffset=fvec_nbuffmpi_recv(lmio)
      m1=fvec_recv_extr(2,lmio)-fvec_recv_extr(1,lmio)+1
      m2=fvec_recv_extr(4,lmio)-fvec_recv_extr(3,lmio)+1
      m3=fvec_recv_extr(6,lmio)-fvec_recv_extr(5,lmio)+1
      imin=fvec_recv_extr(1,lmio)
      imax=fvec_recv_extr(2,lmio)
      jmin=fvec_recv_extr(3,lmio)
      jmax=fvec_recv_extr(4,lmio)
      kmin=fvec_recv_extr(5,lmio)
      kmax=fvec_recv_extr(6,lmio)
#ifdef ACCNOKERNELS
      !$acc parallel loop independent collapse(3) present(fvec_send_buffmpi,fvec_send_extr,auxfields) &
      !$acc& private(ll,idx,ii,jj,kk,xblock,yblock,zblock,myblock)
#else
      !$acc kernels present(fvec_send_buffmpi,fvec_send_extr,auxfields)
      !$acc loop independent collapse(3) private(i,j,k,ll,idx,ii,jj,kk,xblock,yblock,zblock,myblock)
#endif
      do k=kmin,kmax
         do j=jmin,jmax
            do i=imin,imax
            
               xblock=(i+2*TILE_DIMx-1)/TILE_DIMx   
               yblock=(j+2*TILE_DIMy-1)/TILE_DIMy     
               zblock=(k+2*TILE_DIMz-1)/TILE_DIMz  
               myblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
               ii=i-xblock*TILE_DIMx+2*TILE_DIMx
               jj=j-yblock*TILE_DIMy+2*TILE_DIMy
               kk=k-zblock*TILE_DIMz+2*TILE_DIMz   
               
               !linearizzo con l'ordine naturale e metto nel buffer unico per tutte le direzioni
               !poi mandero solo i pezzi contigui che mi servono per la data direzione

               !scorro sul numero di campi da prendere il massimo previsto è num_auxfields_datampi
#ifdef TWOCOMPONENT

               ll=1
               idx=myoffset+(i-fvec_recv_extr(1,lmio))+(j-fvec_recv_extr(3,lmio))*m1+(&
                  k-fvec_recv_extr(5,lmio))*(m1*m2)+(ll-1)*(m1*m2*m3)
               auxfields(ii,jj,kk,1,myblock)=fvec_recv_buffmpi(idx)

               ll=2
               idx=myoffset+(i-fvec_recv_extr(1,lmio))+(j-fvec_recv_extr(3,lmio))*m1+(&
                  k-fvec_recv_extr(5,lmio))*(m1*m2)+(ll-1)*(m1*m2*m3)
               auxfields(ii,jj,kk,2,myblock)=fvec_recv_buffmpi(idx)

               ll=3
               idx=myoffset+(i-fvec_recv_extr(1,lmio))+(j-fvec_recv_extr(3,lmio))*m1+(&
                  k-fvec_recv_extr(5,lmio))*(m1*m2)+(ll-1)*(m1*m2*m3)
               auxfields(ii,jj,kk,3,myblock)=fvec_recv_buffmpi(idx)

               ll=4
               idx=myoffset+(i-fvec_recv_extr(1,lmio))+(j-fvec_recv_extr(3,lmio))*m1+(&
                  k-fvec_recv_extr(5,lmio))*(m1*m2)+(ll-1)*(m1*m2*m3)
               auxfields(ii,jj,kk,4,myblock)=fvec_recv_buffmpi(idx)
			   
			   ll=5
               idx=myoffset+(i-fvec_recv_extr(1,lmio))+(j-fvec_recv_extr(3,lmio))*m1+(&
                  k-fvec_recv_extr(5,lmio))*(m1*m2)+(ll-1)*(m1*m2*m3)
               auxfields(ii,jj,kk,5,myblock)=fvec_recv_buffmpi(idx)
			   
			   ll=6
               idx=myoffset+(i-fvec_recv_extr(1,lmio))+(j-fvec_recv_extr(3,lmio))*m1+(&
                  k-fvec_recv_extr(5,lmio))*(m1*m2)+(ll-1)*(m1*m2*m3)
               auxfields(ii,jj,kk,6,myblock)=fvec_recv_buffmpi(idx)
			   
			   ll=7
               idx=myoffset+(i-fvec_recv_extr(1,lmio))+(j-fvec_recv_extr(3,lmio))*m1+(&
                  k-fvec_recv_extr(5,lmio))*(m1*m2)+(ll-1)*(m1*m2*m3)
               auxfields(ii,jj,kk,7,myblock)=fvec_recv_buffmpi(idx)

#endif
            enddo
         enddo
      enddo
#ifdef ACCNOKERNELS
      !$acc end parallel loop
#else
      !$acc end kernels
#endif

   end subroutine depackaging_auxfields_buffmpi
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!hfields!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine exchange_hfields_intpbc(hfields_s)

      implicit none
      
      real(kind=strdb), allocatable, dimension(:,:,:,:,:) :: hfields_s
      integer :: imin,imax,jmin,jmax,kmin,kmax
      integer :: lmio,oi,oj,ok
      integer :: ii,jj,kk,xblock,yblock,zblock,myblock
      integer :: oii,ojj,okk
      integer :: oxblock,oyblock,ozblock,omyblock
      
!!!!!b_send_extr   !!!!!!!!!!!!b_recv_extr
      do lmio=1,nlinksmpi
         if(.not. lintpbc_dir(lmio))cycle
         imin=b_recv_extr(1,lmio)
         imax=b_recv_extr(2,lmio)
         jmin=b_recv_extr(3,lmio)
         jmax=b_recv_extr(4,lmio)
         kmin=b_recv_extr(5,lmio)
         kmax=b_recv_extr(6,lmio)
#ifdef ACCNOKERNELS
         !$acc parallel loop independent collapse(3) present(intpbc_dir,hfields_s) &
         !$acc& private(ii,jj,kk,oi,oj,ok,oii,ojj,okk,oxblock,oyblock,ozblock, &
         !$acc& omyblock,xblock,yblock,zblock,myblock)
#else
         !$acc kernels present(intpbc_dir,hfields_s)
         !$acc loop independent collapse(3) &
         !$acc& private(i,j,k,ii,jj,kk,oi,oj,ok,oii,ojj,okk,oxblock,oyblock,ozblock, &
         !$acc& omyblock,xblock,yblock,zblock,myblock)
#endif
         do k=kmin,kmax
           do j=jmin,jmax
              do i=imin,imax
                  oi=i
                  oj=j
                  ok=k
                  if(intpbc_dir(1,lmio))oi=mod(oi+nx-1,nx)+1
                  if(intpbc_dir(2,lmio))oj=mod(oj+ny-1,ny)+1
                  if(intpbc_dir(3,lmio))ok=mod(ok+nz-1,nz)+1
                  
                  oxblock=(oi+2*TILE_DIMx-1)/TILE_DIMx   
                  oyblock=(oj+2*TILE_DIMy-1)/TILE_DIMy     
                  ozblock=(ok+2*TILE_DIMz-1)/TILE_DIMz 
                  omyblock=(oxblock-1)+(oyblock-1)*nxblock+(ozblock-1)*nxyblock+1
                  oii=oi-oxblock*TILE_DIMx+2*TILE_DIMx
                  ojj=oj-oyblock*TILE_DIMy+2*TILE_DIMy
                  okk=ok-ozblock*TILE_DIMz+2*TILE_DIMz
                  
                  xblock=(i+2*TILE_DIMx-1)/TILE_DIMx   
                  yblock=(j+2*TILE_DIMy-1)/TILE_DIMy     
                  zblock=(k+2*TILE_DIMz-1)/TILE_DIMz  
                  myblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
                  ii=i-xblock*TILE_DIMx+2*TILE_DIMx
                  jj=j-yblock*TILE_DIMy+2*TILE_DIMy
                  kk=k-zblock*TILE_DIMz+2*TILE_DIMz 
                  
                  
                  
                  hfields_s(ii,jj,kk,1,myblock)= &
                   hfields_s(oii,ojj,okk,1,omyblock)
                  
                  hfields_s(ii,jj,kk,2,myblock)= &
                   hfields_s(oii,ojj,okk,2,omyblock)
                  
                  hfields_s(ii,jj,kk,3,myblock)= &
                   hfields_s(oii,ojj,okk,3,omyblock)
                  
                  hfields_s(ii,jj,kk,4,myblock)= &
                   hfields_s(oii,ojj,okk,4,omyblock)
                  
                  hfields_s(ii,jj,kk,5,myblock)= &
                   hfields_s(oii,ojj,okk,5,omyblock)
                  
                  hfields_s(ii,jj,kk,6,myblock)= &
                   hfields_s(oii,ojj,okk,6,omyblock)
                  
                  hfields_s(ii,jj,kk,7,myblock)= &
                   hfields_s(oii,ojj,okk,7,omyblock)
                  
                  hfields_s(ii,jj,kk,8,myblock)= &
                   hfields_s(oii,ojj,okk,8,omyblock)
                  
                  hfields_s(ii,jj,kk,9,myblock)= &
                   hfields_s(oii,ojj,okk,9,omyblock)
                  
                  hfields_s(ii,jj,kk,10,myblock)= &
                   hfields_s(oii,ojj,okk,10,omyblock)
                  
               enddo
            enddo
         enddo
#ifdef ACCNOKERNELS
        !$acc end parallel loop
#else
        !$acc end kernels
#endif
      enddo

   end subroutine exchange_hfields_intpbc

   subroutine exchange_hfields_sendrecv(hfields_s)

      implicit none
      
      real(kind=strdb), allocatable, dimension(:,:,:,:,:) :: hfields_s

      integer :: l,ll,myoffset,tag,ierr,mm
#ifndef MPI
      return
#endif
      do l=1,nlinksmpi
         if(lsend_dir(l))call packaging_hfields_buffmpi(l,hfields_s)
      enddo
      !$acc wait
#ifdef MPI
      mm=0
      do l=1,nlinksmpi
         ll=(l+1)/2
         if(lsend_dir(l))then
            mm=mm+1
            myoffset=b_nbuffmpi_send(l)
            !$acc host_data use_device(b_send_buffmpi)
            call mpi_isend(b_send_buffmpi(myoffset),1,b_datampi(ll),send_dir(l), &
               b_mpitag(l),lbecomm,b_reqs(mm),ierr)
            !$acc end host_data
         endif
         if(lrecv_dir(l))then
            mm=mm+1
            myoffset=b_nbuffmpi_recv(l)
            !$acc host_data use_device(b_recv_buffmpi)
            call mpi_irecv(b_recv_buffmpi(myoffset),1,b_datampi(ll),recv_dir(l), &
               b_mpitag(l),lbecomm,b_reqs(mm),ierr)
            !$acc end host_data
         endif
      enddo
      nb_reqs=mm
#endif


   end subroutine exchange_hfields_sendrecv

   subroutine exchange_hfields_wait(hfields_s)

      implicit none
      
      real(kind=strdb), allocatable, dimension(:,:,:,:,:) :: hfields_s

      integer :: l,ll,myoffset,tag,ierr
      integer, dimension(:), allocatable :: myierr
      integer, dimension(:,:), allocatable :: mystatus
#ifndef MPI
      return
#endif
#ifdef MPI

      allocate(myierr(nb_reqs))
      myierr=0

      allocate(mystatus(MPI_STATUS_SIZE,nb_reqs))
      !$acc wait
      call mpi_waitall(nb_reqs,b_reqs,mystatus,ierr)
      !$acc wait

      if(any(myierr.ne.0))call doerror(6,'ERROR in exchange_hfields_wait')
#endif
      do l=1,nlinksmpi
         if(lrecv_dir(l))call depackaging_hfields_buffmpi(l,hfields_s)
      enddo

   end subroutine exchange_hfields_wait

   subroutine packaging_hfields_buffmpi(lmio,hfields_s)

      implicit none

      integer, intent(in) :: lmio
      real(kind=strdb), allocatable, dimension(:,:,:,:,:) :: hfields_s

      integer :: myoffset

      integer :: i,j,k,l,ll=0,m1,m2,m3
      integer :: imin,imax,jmin,jmax,kmin,kmax
      integer :: idx=0
      integer :: ii,jj,kk,xblock,yblock,zblock,myblock

      myoffset=b_nbuffmpi_send(lmio)
      m1=b_send_extr(2,lmio)-b_send_extr(1,lmio)+1
      m2=b_send_extr(4,lmio)-b_send_extr(3,lmio)+1
      m3=b_send_extr(6,lmio)-b_send_extr(5,lmio)+1
      imin=b_send_extr(1,lmio)
      imax=b_send_extr(2,lmio)
      jmin=b_send_extr(3,lmio)
      jmax=b_send_extr(4,lmio)
      kmin=b_send_extr(5,lmio)
      kmax=b_send_extr(6,lmio)
#ifdef ACCNOKERNELS
      !$acc parallel loop independent collapse(3) present(b_send_buffmpi,b_send_extr,hfields_s) &
      !$acc& private(ll,idx,ii,jj,kk,xblock,yblock,zblock,myblock)
#else
      !$acc kernels present(b_send_buffmpi,b_send_extr,hfields_s)
      !$acc loop independent collapse(3)  private(i,j,k,ll,idx,ii,jj,kk,xblock,yblock,zblock,myblock)
#endif
      do k=kmin,kmax
         do j=jmin,jmax
            do i=imin,imax
               
               xblock=(i+2*TILE_DIMx-1)/TILE_DIMx   
               yblock=(j+2*TILE_DIMy-1)/TILE_DIMy     
               zblock=(k+2*TILE_DIMz-1)/TILE_DIMz  
               myblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
               ii=i-xblock*TILE_DIMx+2*TILE_DIMx
               jj=j-yblock*TILE_DIMy+2*TILE_DIMy
               kk=k-zblock*TILE_DIMz+2*TILE_DIMz   
               !linearizzo con l'ordine naturale e metto nel buffer unico per tutte le direzioni
               !poi mandero solo i pezzi contigui che mi servono per la data direzione

               !scorro sul numero di campi da prendere il massimo previsto è num_hfields_datampi

               ll=1
               idx=myoffset+(i-b_send_extr(1,lmio))+(j-b_send_extr(3,lmio))*m1+(&
                  k-b_send_extr(5,lmio))*(m1*m2)+(ll-1)*(m1*m2*m3)
               b_send_buffmpi(idx)=hfields_s(ii,jj,kk,1,myblock)
               
               ll=2
               idx=myoffset+(i-b_send_extr(1,lmio))+(j-b_send_extr(3,lmio))*m1+(&
                  k-b_send_extr(5,lmio))*(m1*m2)+(ll-1)*(m1*m2*m3)
               b_send_buffmpi(idx)=hfields_s(ii,jj,kk,2,myblock)
               
               ll=3
               idx=myoffset+(i-b_send_extr(1,lmio))+(j-b_send_extr(3,lmio))*m1+(&
                  k-b_send_extr(5,lmio))*(m1*m2)+(ll-1)*(m1*m2*m3)
               b_send_buffmpi(idx)=hfields_s(ii,jj,kk,3,myblock)
               
               ll=4
               idx=myoffset+(i-b_send_extr(1,lmio))+(j-b_send_extr(3,lmio))*m1+(&
                  k-b_send_extr(5,lmio))*(m1*m2)+(ll-1)*(m1*m2*m3)
               b_send_buffmpi(idx)=hfields_s(ii,jj,kk,4,myblock)
               
               ll=5
               idx=myoffset+(i-b_send_extr(1,lmio))+(j-b_send_extr(3,lmio))*m1+(&
                  k-b_send_extr(5,lmio))*(m1*m2)+(ll-1)*(m1*m2*m3)
               b_send_buffmpi(idx)=hfields_s(ii,jj,kk,5,myblock)
               
               ll=6
               idx=myoffset+(i-b_send_extr(1,lmio))+(j-b_send_extr(3,lmio))*m1+(&
                  k-b_send_extr(5,lmio))*(m1*m2)+(ll-1)*(m1*m2*m3)
               b_send_buffmpi(idx)=hfields_s(ii,jj,kk,6,myblock)
               
               ll=7
               idx=myoffset+(i-b_send_extr(1,lmio))+(j-b_send_extr(3,lmio))*m1+(&
                  k-b_send_extr(5,lmio))*(m1*m2)+(ll-1)*(m1*m2*m3)
               b_send_buffmpi(idx)=hfields_s(ii,jj,kk,7,myblock)
               
               ll=8
               idx=myoffset+(i-b_send_extr(1,lmio))+(j-b_send_extr(3,lmio))*m1+(&
                  k-b_send_extr(5,lmio))*(m1*m2)+(ll-1)*(m1*m2*m3)
               b_send_buffmpi(idx)=hfields_s(ii,jj,kk,8,myblock)
               
               ll=9
               idx=myoffset+(i-b_send_extr(1,lmio))+(j-b_send_extr(3,lmio))*m1+(&
                  k-b_send_extr(5,lmio))*(m1*m2)+(ll-1)*(m1*m2*m3)
               b_send_buffmpi(idx)=hfields_s(ii,jj,kk,9,myblock)
               
               ll=10
               idx=myoffset+(i-b_send_extr(1,lmio))+(j-b_send_extr(3,lmio))*m1+(&
                  k-b_send_extr(5,lmio))*(m1*m2)+(ll-1)*(m1*m2*m3)
               b_send_buffmpi(idx)=hfields_s(ii,jj,kk,10,myblock)

            enddo
         enddo
      enddo
#ifdef ACCNOKERNELS
      !$acc end parallel loop
#else
      !$acc end kernels
#endif


   end subroutine packaging_hfields_buffmpi

   subroutine depackaging_hfields_buffmpi(lmio,hfields_s)

      implicit none

      integer, intent(in) :: lmio
      real(kind=strdb), allocatable, dimension(:,:,:,:,:) :: hfields_s
      
      integer :: myoffset

      integer :: i,j,k,l,ll=0,m1,m2,m3
      integer :: imin,imax,jmin,jmax,kmin,kmax
      integer :: idx=0
      integer :: ii,jj,kk,xblock,yblock,zblock,myblock

      myoffset=b_nbuffmpi_recv(lmio)
      m1=b_recv_extr(2,lmio)-b_recv_extr(1,lmio)+1
      m2=b_recv_extr(4,lmio)-b_recv_extr(3,lmio)+1
      m3=b_recv_extr(6,lmio)-b_recv_extr(5,lmio)+1
      imin=b_recv_extr(1,lmio)
      imax=b_recv_extr(2,lmio)
      jmin=b_recv_extr(3,lmio)
      jmax=b_recv_extr(4,lmio)
      kmin=b_recv_extr(5,lmio)
      kmax=b_recv_extr(6,lmio)
#ifdef ACCNOKERNELS
      !$acc parallel loop independent collapse(3) present(b_recv_buffmpi,b_recv_extr,hfields_s) &
      !$acc& private(ll,idx,ii,jj,kk,xblock,yblock,zblock,myblock)
#else
      !$acc kernels present(b_recv_buffmpi,b_recv_extr,hfields_s)
      !$acc loop independent collapse(3) private(i,j,k,ll,idx,ii,jj,kk,xblock,yblock,zblock,myblock)
#endif
      do k=kmin,kmax
         do j=jmin,jmax
            do i=imin,imax
               
               xblock=(i+2*TILE_DIMx-1)/TILE_DIMx   
               yblock=(j+2*TILE_DIMy-1)/TILE_DIMy     
               zblock=(k+2*TILE_DIMz-1)/TILE_DIMz  
               myblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
               ii=i-xblock*TILE_DIMx+2*TILE_DIMx
               jj=j-yblock*TILE_DIMy+2*TILE_DIMy
               kk=k-zblock*TILE_DIMz+2*TILE_DIMz   
               !linearizzo con l'ordine naturale e metto nel buffer unico per tutte le direzioni
               !poi mandero solo i pezzi contigui che mi servono per la data direzione

               !scorro sul numero di campi da prendere il massimo previsto è num_hfields_datampi

               ll=1
               idx=myoffset+(i-b_recv_extr(1,lmio))+(j-b_recv_extr(3,lmio))*m1+(&
                  k-b_recv_extr(5,lmio))*(m1*m2)+(ll-1)*(m1*m2*m3)
               hfields_s(ii,jj,kk,1,myblock)=b_recv_buffmpi(idx)
               
               ll=2
               idx=myoffset+(i-b_recv_extr(1,lmio))+(j-b_recv_extr(3,lmio))*m1+(&
                  k-b_recv_extr(5,lmio))*(m1*m2)+(ll-1)*(m1*m2*m3)
               hfields_s(ii,jj,kk,2,myblock)=b_recv_buffmpi(idx)
               
               ll=3
               idx=myoffset+(i-b_recv_extr(1,lmio))+(j-b_recv_extr(3,lmio))*m1+(&
                  k-b_recv_extr(5,lmio))*(m1*m2)+(ll-1)*(m1*m2*m3)
               hfields_s(ii,jj,kk,3,myblock)=b_recv_buffmpi(idx)
               
               ll=4
               idx=myoffset+(i-b_recv_extr(1,lmio))+(j-b_recv_extr(3,lmio))*m1+(&
                  k-b_recv_extr(5,lmio))*(m1*m2)+(ll-1)*(m1*m2*m3)
               hfields_s(ii,jj,kk,4,myblock)=b_recv_buffmpi(idx)
               
               ll=5
               idx=myoffset+(i-b_recv_extr(1,lmio))+(j-b_recv_extr(3,lmio))*m1+(&
                  k-b_recv_extr(5,lmio))*(m1*m2)+(ll-1)*(m1*m2*m3)
               hfields_s(ii,jj,kk,5,myblock)=b_recv_buffmpi(idx)
               
               ll=6
               idx=myoffset+(i-b_recv_extr(1,lmio))+(j-b_recv_extr(3,lmio))*m1+(&
                  k-b_recv_extr(5,lmio))*(m1*m2)+(ll-1)*(m1*m2*m3)
               hfields_s(ii,jj,kk,6,myblock)=b_recv_buffmpi(idx)
               
               ll=7
               idx=myoffset+(i-b_recv_extr(1,lmio))+(j-b_recv_extr(3,lmio))*m1+(&
                  k-b_recv_extr(5,lmio))*(m1*m2)+(ll-1)*(m1*m2*m3)
               hfields_s(ii,jj,kk,7,myblock)=b_recv_buffmpi(idx)
               
               ll=8
               idx=myoffset+(i-b_recv_extr(1,lmio))+(j-b_recv_extr(3,lmio))*m1+(&
                  k-b_recv_extr(5,lmio))*(m1*m2)+(ll-1)*(m1*m2*m3)
               hfields_s(ii,jj,kk,8,myblock)=b_recv_buffmpi(idx)
               
               ll=9
               idx=myoffset+(i-b_recv_extr(1,lmio))+(j-b_recv_extr(3,lmio))*m1+(&
                  k-b_recv_extr(5,lmio))*(m1*m2)+(ll-1)*(m1*m2*m3)
               hfields_s(ii,jj,kk,9,myblock)=b_recv_buffmpi(idx)
               
               ll=10
               idx=myoffset+(i-b_recv_extr(1,lmio))+(j-b_recv_extr(3,lmio))*m1+(&
                  k-b_recv_extr(5,lmio))*(m1*m2)+(ll-1)*(m1*m2*m3)
               hfields_s(ii,jj,kk,10,myblock)=b_recv_buffmpi(idx)

            enddo
         enddo
      enddo
#ifdef ACCNOKERNELS
      !$acc end parallel loop
#else
      !$acc end kernels
#endif


   end subroutine depackaging_hfields_buffmpi   
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!forces!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine exchange_forces_intpbc()

      implicit none
      
      integer :: imin,imax,jmin,jmax,kmin,kmax
      integer :: lmio,oi,oj,ok
      integer :: ii,jj,kk,xblock,yblock,zblock,myblock
      integer :: oii,ojj,okk
      integer :: oxblock,oyblock,ozblock,omyblock
      
!!!!!c_send_extr   !!!!!!!!!!!!c_recv_extr
      do lmio=1,nlinksmpi
         if(.not. lintpbc_dir(lmio))cycle
         imin=c_recv_extr(1,lmio)
         imax=c_recv_extr(2,lmio)
         jmin=c_recv_extr(3,lmio)
         jmax=c_recv_extr(4,lmio)
         kmin=c_recv_extr(5,lmio)
         kmax=c_recv_extr(6,lmio)
#ifdef ACCNOKERNELS
         !$acc parallel loop independent collapse(3) present(intpbc_dir,forces) &
         !$acc& private(ii,jj,kk,oi,oj,ok,oii,ojj,okk,oxblock,oyblock,ozblock, &
         !$acc& omyblock,xblock,yblock,zblock,myblock)
#else
         !$acc kernels present(intpbc_dir,forces)
         !$acc loop independent collapse(3) &
         !$acc& private(i,j,k,ii,jj,kk,oi,oj,ok,oii,ojj,okk,oxblock,oyblock,ozblock, &
         !$acc& omyblock,xblock,yblock,zblock,myblock)
#endif
         do k=kmin,kmax
           do j=jmin,jmax
              do i=imin,imax
                  oi=i
                  oj=j
                  ok=k
                  if(intpbc_dir(1,lmio))oi=mod(oi+nx-1,nx)+1
                  if(intpbc_dir(2,lmio))oj=mod(oj+ny-1,ny)+1
                  if(intpbc_dir(3,lmio))ok=mod(ok+nz-1,nz)+1
                  
                  oxblock=(oi+2*TILE_DIMx-1)/TILE_DIMx   
                  oyblock=(oj+2*TILE_DIMy-1)/TILE_DIMy     
                  ozblock=(ok+2*TILE_DIMz-1)/TILE_DIMz 
                  omyblock=(oxblock-1)+(oyblock-1)*nxblock+(ozblock-1)*nxyblock+1
                  oii=oi-oxblock*TILE_DIMx+2*TILE_DIMx
                  ojj=oj-oyblock*TILE_DIMy+2*TILE_DIMy
                  okk=ok-ozblock*TILE_DIMz+2*TILE_DIMz
                  
                  xblock=(i+2*TILE_DIMx-1)/TILE_DIMx   
                  yblock=(j+2*TILE_DIMy-1)/TILE_DIMy     
                  zblock=(k+2*TILE_DIMz-1)/TILE_DIMz  
                  myblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
                  ii=i-xblock*TILE_DIMx+2*TILE_DIMx
                  jj=j-yblock*TILE_DIMy+2*TILE_DIMy
                  kk=k-zblock*TILE_DIMz+2*TILE_DIMz 
                  
                  forces(ii,jj,kk,1,myblock)= &
                   forces(oii,ojj,okk,1,omyblock)
                  
                  forces(ii,jj,kk,2,myblock)= &
                   forces(oii,ojj,okk,2,omyblock)
                  
                  forces(ii,jj,kk,3,myblock)= &
                   forces(oii,ojj,okk,3,omyblock)
                  
               enddo
            enddo
         enddo
#ifdef ACCNOKERNELS
        !$acc end parallel loop
#else
        !$acc end kernels
#endif
      enddo

   end subroutine exchange_forces_intpbc

   subroutine exchange_forces_sendrecv()

      implicit none

      integer :: l,ll,myoffset,tag,ierr,mm
#ifndef MPI
      return
#endif
      do l=1,nlinksmpi
         if(lsend_dir(l))call packaging_forces_buffmpi(l)
      enddo
      !$acc wait
#ifdef MPI
      mm=0
      do l=1,nlinksmpi
         ll=(l+1)/2
         if(lsend_dir(l))then
            mm=mm+1
            myoffset=c_nbuffmpi_send(l)
            !$acc host_data use_device(c_send_buffmpi)
            call mpi_isend(c_send_buffmpi(myoffset),1,c_datampi(ll),send_dir(l), &
               c_mpitag(l),lbecomm,c_reqs(mm),ierr)
            !$acc end host_data
         endif
         if(lrecv_dir(l))then
            mm=mm+1
            myoffset=c_nbuffmpi_recv(l)
            !$acc host_data use_device(c_recv_buffmpi)
            call mpi_irecv(c_recv_buffmpi(myoffset),1,c_datampi(ll),recv_dir(l), &
               c_mpitag(l),lbecomm,c_reqs(mm),ierr)
            !$acc end host_data
         endif
      enddo
      nc_reqs=mm
#endif


   end subroutine exchange_forces_sendrecv

   subroutine exchange_forces_wait()

      implicit none

      integer :: l,ll,myoffset,tag,ierr
      integer, dimension(:), allocatable :: myierr
      integer, dimension(:,:), allocatable :: mystatus
#ifndef MPI
      return
#endif
#ifdef MPI

      allocate(myierr(nc_reqs))
      myierr=0

      allocate(mystatus(MPI_STATUS_SIZE,nc_reqs))
      !$acc wait
      call mpi_waitall(nc_reqs,c_reqs,mystatus,ierr)
      !$acc wait

      if(any(myierr.ne.0))call doerror(6,'ERROR in exchange_forces_wait')
#endif
      do l=1,nlinksmpi
         if(lrecv_dir(l))call depackaging_forces_buffmpi(l)
      enddo

   end subroutine exchange_forces_wait

   subroutine packaging_forces_buffmpi(lmio)

      implicit none

      integer, intent(in) :: lmio

      integer :: myoffset

      integer :: i,j,k,l,ll=0,m1,m2,m3
      integer :: imin,imax,jmin,jmax,kmin,kmax
      integer :: idx=0
      integer :: ii,jj,kk,xblock,yblock,zblock,myblock

      myoffset=c_nbuffmpi_send(lmio)
      m1=c_send_extr(2,lmio)-c_send_extr(1,lmio)+1
      m2=c_send_extr(4,lmio)-c_send_extr(3,lmio)+1
      m3=c_send_extr(6,lmio)-c_send_extr(5,lmio)+1
      imin=c_send_extr(1,lmio)
      imax=c_send_extr(2,lmio)
      jmin=c_send_extr(3,lmio)
      jmax=c_send_extr(4,lmio)
      kmin=c_send_extr(5,lmio)
      kmax=c_send_extr(6,lmio)
#ifdef ACCNOKERNELS
      !$acc parallel loop independent collapse(3) present(c_send_buffmpi,c_send_extr,forces) &
      !$acc& private(ll,idx,ii,jj,kk,xblock,yblock,zblock,myblock)
#else
      !$acc kernels present(c_send_buffmpi,c_send_extr,forces)
      !$acc loop independent collapse(3)  private(i,j,k,ll,idx,ii,jj,kk,xblock,yblock,zblock,myblock)
#endif
      do k=kmin,kmax
         do j=jmin,jmax
            do i=imin,imax
               
               xblock=(i+2*TILE_DIMx-1)/TILE_DIMx   
               yblock=(j+2*TILE_DIMy-1)/TILE_DIMy     
               zblock=(k+2*TILE_DIMz-1)/TILE_DIMz  
               myblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
               ii=i-xblock*TILE_DIMx+2*TILE_DIMx
               jj=j-yblock*TILE_DIMy+2*TILE_DIMy
               kk=k-zblock*TILE_DIMz+2*TILE_DIMz   
               !linearizzo con l'ordine naturale e metto nel buffer unico per tutte le direzioni
               !poi mandero solo i pezzi contigui che mi servono per la data direzione

               !scorro sul numero di campi da prendere il massimo previsto è num_forces_datampi

               ll=1
               idx=myoffset+(i-c_send_extr(1,lmio))+(j-c_send_extr(3,lmio))*m1+(&
                  k-c_send_extr(5,lmio))*(m1*m2)+(ll-1)*(m1*m2*m3)
               c_send_buffmpi(idx)=forces(ii,jj,kk,1,myblock)
               
               ll=2
               idx=myoffset+(i-c_send_extr(1,lmio))+(j-c_send_extr(3,lmio))*m1+(&
                  k-c_send_extr(5,lmio))*(m1*m2)+(ll-1)*(m1*m2*m3)
               c_send_buffmpi(idx)=forces(ii,jj,kk,2,myblock)
               
               ll=3
               idx=myoffset+(i-c_send_extr(1,lmio))+(j-c_send_extr(3,lmio))*m1+(&
                  k-c_send_extr(5,lmio))*(m1*m2)+(ll-1)*(m1*m2*m3)
               c_send_buffmpi(idx)=forces(ii,jj,kk,3,myblock)

            enddo
         enddo
      enddo
#ifdef ACCNOKERNELS
      !$acc end parallel loop
#else
      !$acc end kernels
#endif


   end subroutine packaging_forces_buffmpi

   subroutine depackaging_forces_buffmpi(lmio)

      implicit none

      integer, intent(in) :: lmio
      
      integer :: myoffset

      integer :: i,j,k,l,ll=0,m1,m2,m3
      integer :: imin,imax,jmin,jmax,kmin,kmax
      integer :: idx=0
      integer :: ii,jj,kk,xblock,yblock,zblock,myblock

      myoffset=c_nbuffmpi_recv(lmio)
      m1=c_recv_extr(2,lmio)-c_recv_extr(1,lmio)+1
      m2=c_recv_extr(4,lmio)-c_recv_extr(3,lmio)+1
      m3=c_recv_extr(6,lmio)-c_recv_extr(5,lmio)+1
      imin=c_recv_extr(1,lmio)
      imax=c_recv_extr(2,lmio)
      jmin=c_recv_extr(3,lmio)
      jmax=c_recv_extr(4,lmio)
      kmin=c_recv_extr(5,lmio)
      kmax=c_recv_extr(6,lmio)
#ifdef ACCNOKERNELS
      !$acc parallel loop independent collapse(3) present(c_recv_buffmpi,c_recv_extr,forces) &
      !$acc& private(ll,idx,ii,jj,kk,xblock,yblock,zblock,myblock)
#else
      !$acc kernels present(c_recv_buffmpi,c_recv_extr,forces)
      !$acc loop independent collapse(3) private(i,j,k,ll,idx,ii,jj,kk,xblock,yblock,zblock,myblock)
#endif
      do k=kmin,kmax
         do j=jmin,jmax
            do i=imin,imax
               
               xblock=(i+2*TILE_DIMx-1)/TILE_DIMx   
               yblock=(j+2*TILE_DIMy-1)/TILE_DIMy     
               zblock=(k+2*TILE_DIMz-1)/TILE_DIMz  
               myblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
               ii=i-xblock*TILE_DIMx+2*TILE_DIMx
               jj=j-yblock*TILE_DIMy+2*TILE_DIMy
               kk=k-zblock*TILE_DIMz+2*TILE_DIMz   
               !linearizzo con l'ordine naturale e metto nel buffer unico per tutte le direzioni
               !poi mandero solo i pezzi contigui che mi servono per la data direzione

               !scorro sul numero di campi da prendere il massimo previsto è num_forces_datampi

               ll=1
               idx=myoffset+(i-c_recv_extr(1,lmio))+(j-c_recv_extr(3,lmio))*m1+(&
                  k-c_recv_extr(5,lmio))*(m1*m2)+(ll-1)*(m1*m2*m3)
               forces(ii,jj,kk,1,myblock)=c_recv_buffmpi(idx)
               
               ll=2
               idx=myoffset+(i-c_recv_extr(1,lmio))+(j-c_recv_extr(3,lmio))*m1+(&
                  k-c_recv_extr(5,lmio))*(m1*m2)+(ll-1)*(m1*m2*m3)
               forces(ii,jj,kk,2,myblock)=c_recv_buffmpi(idx)
               
               ll=3
               idx=myoffset+(i-c_recv_extr(1,lmio))+(j-c_recv_extr(3,lmio))*m1+(&
                  k-c_recv_extr(5,lmio))*(m1*m2)+(ll-1)*(m1*m2*m3)
               forces(ii,jj,kk,3,myblock)=c_recv_buffmpi(idx)

            enddo
         enddo
      enddo
#ifdef ACCNOKERNELS
      !$acc end parallel loop
#else
      !$acc end kernels
#endif


   end subroutine depackaging_forces_buffmpi   
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!ISFLUID!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine exchange_isf_intpbc

      implicit none

      integer :: lmio,oi,oj,ok

      do lmio=1,nlinksmpi
         if(.not. lintpbc_dir(lmio))cycle
         do k=f_recv_extr(5,lmio),f_recv_extr(6,lmio)
            ok=k
            if(intpbc_dir(3,lmio))ok=mod(ok+nz-1,nz)+1
            do j=f_recv_extr(3,lmio),f_recv_extr(4,lmio)
               oj=j
               if(intpbc_dir(2,lmio))oj=mod(oj+ny-1,ny)+1
               do i=f_recv_extr(1,lmio),f_recv_extr(2,lmio)
                  oi=i
                  if(intpbc_dir(1,lmio))oi=mod(oi+nx-1,nx)+1
                  isfluid(i,j,k)=isfluid(oi,oj,ok)
               enddo
            enddo
         enddo
      enddo

   end subroutine exchange_isf_intpbc

   subroutine exchange_isf_sendrecv

      implicit none

      integer :: l,ll,myoffset,tag,ierr,mm

      do l=1,nlinksmpi
         if(lsend_dir(l))call packaging_isf_buffmpi(l)
      enddo
      !i_recv_buffmpi=i_send_buffmpi
#ifdef MPI
      mm=0
      do l=1,nlinksmpi
         ll=(l+1)/2
         if(lsend_dir(l))then
            mm=mm+1
            myoffset=i_nbuffmpi_send(l)
            call mpi_isend(i_send_buffmpi(myoffset),1,i_datampi(ll),send_dir(l), &
               i_mpitag(l),lbecomm,i_reqs(mm),ierr)
         endif
         if(lrecv_dir(l))then
            mm=mm+1
            myoffset=i_nbuffmpi_recv(l)
            call mpi_irecv(i_recv_buffmpi(myoffset),1,i_datampi(ll),recv_dir(l), &
               i_mpitag(l),lbecomm,i_reqs(mm),ierr)
         endif
      enddo
      ni_reqs=mm
#endif


   end subroutine exchange_isf_sendrecv

   subroutine exchange_isf_wait

      implicit none

      integer :: l,ll,myoffset,tag,ierr
#ifdef MPI
      integer, dimension(:), allocatable :: myierr
      integer, dimension(:,:), allocatable :: mystatus

      allocate(myierr(ni_reqs))
      myierr=0

      allocate(mystatus(MPI_STATUS_SIZE,ni_reqs))
      !$acc wait
      call mpi_waitall(ni_reqs,i_reqs,mystatus,ierr)

      if(any(myierr.ne.0))call doerror(6,'ERROR in mpi_wait send')
#endif
      do l=1,nlinksmpi
         if(lrecv_dir(l))call depackaging_isf_buffmpi(l)
      enddo

   end subroutine exchange_isf_wait

   subroutine packaging_isf_buffmpi(lmio)

      implicit none

      integer, intent(in) :: lmio
      integer :: myoffset

      integer :: i,j,k,l,ll,m1,m2,m3

      integer :: idx

      myoffset=i_nbuffmpi_send(lmio)
      m1=f_send_extr(2,lmio)-f_send_extr(1,lmio)+1
      m2=f_send_extr(4,lmio)-f_send_extr(3,lmio)+1
      m3=f_send_extr(6,lmio)-f_send_extr(5,lmio)+1
      !scorro sul numero di campi da prendere (per scalare = 1)
      ll=1
      do k=f_send_extr(5,lmio),f_send_extr(6,lmio)
         do j=f_send_extr(3,lmio),f_send_extr(4,lmio)
            do i=f_send_extr(1,lmio),f_send_extr(2,lmio)
               !linearizzo con l'ordine naturale e metto nel buffer unico per tutte le direzioni
               !poi mandero solo i pezzi contigui che mi servono per la data direzione
               idx=myoffset+(i-f_send_extr(1,lmio))+(j-f_send_extr(3,lmio))*m1+(&
                  k-f_send_extr(5,lmio))*(m1*m2)+(ll-1)*(m1*m2*m3)

               i_send_buffmpi(idx)=isfluid(i,j,k)
            enddo
         enddo
      enddo



   end subroutine packaging_isf_buffmpi

   subroutine depackaging_isf_buffmpi(lmio)

      implicit none

      integer, intent(in) :: lmio
      integer :: myoffset

      integer :: i,j,k,l,ll,m1,m2,m3

      integer :: idx

      myoffset=i_nbuffmpi_recv(lmio)
      m1=f_recv_extr(2,lmio)-f_recv_extr(1,lmio)+1
      m2=f_recv_extr(4,lmio)-f_recv_extr(3,lmio)+1
      m3=f_recv_extr(6,lmio)-f_recv_extr(5,lmio)+1
      !scorro sul numero di campi da prendere (per scalare = 1)
      ll=1
      do k=f_recv_extr(5,lmio),f_recv_extr(6,lmio)
         do j=f_recv_extr(3,lmio),f_recv_extr(4,lmio)
            do i=f_recv_extr(1,lmio),f_recv_extr(2,lmio)
               !linearizzo con l'ordine naturale e metto nel buffer unico per tutte le direzioni
               !poi mandero solo i pezzi contigui che mi servono per la data direzione
               idx=myoffset+(i-f_recv_extr(1,lmio))+(j-f_recv_extr(3,lmio))*m1+(&
                  k-f_recv_extr(5,lmio))*(m1*m2)+(ll-1)*(m1*m2*m3)

               isfluid(i,j,k)=i_recv_buffmpi(idx)
            enddo
         enddo
      enddo



   end subroutine depackaging_isf_buffmpi

   subroutine write_file_vtk_par(iframe,e_io)

!***********************************************************************
!
!     LBsoft subroutine for opening the vtk legacy file
!     in parallel IO
!
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2019
!
!***********************************************************************

      implicit none

      integer, intent(in) ::iframe
      integer, intent(out) :: e_io
#ifdef MPI
      integer(kind=MPI_OFFSET_KIND) :: tempoffset
#endif
      integer :: nns
      character(len=500) :: sheadervtk

      integer :: ioffset
      character(1), parameter :: end_rec = char(10)
      integer, parameter :: bytechar=kind(end_rec)
      integer, parameter :: byteint = 4
      integer, parameter :: byter4  = 4
      integer, parameter :: byter8  = 8
      integer, parameter :: nbuffsub = 0
      integer :: filetypesub,imemtype,filetypesubv,ierr

      integer, dimension(3) :: memDims,memOffs
      integer, dimension(4) :: velglobalDims,velldims,velmystarts, &
         velmemDims,velmemOffs
      integer :: elen,amode
      logical :: lexist
#ifdef MPI
      integer :: fdens=MPI_FILE_NULL
      integer :: fvel=MPI_FILE_NULL
      character(len=MPI_MAX_ERROR_STRING) :: emsg
      
     


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!density!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call MPI_Barrier(MPI_COMM_WORLD, e_io)
      sevt1 = trim(dir_out) // trim(filenamevtk)//'_'//trim(namevarvtk(1))// &
         '_'//trim(write_fmtnumb(iframe)) // '.vti'

      if (myrank == 0) then
        inquire(file=trim(sevt1), exist=lexist)
        if (lexist) then
          call MPI_File_delete(trim(sevt1), MPI_INFO_NULL, e_io)
        endif
      endif
      call MPI_Barrier(MPI_COMM_WORLD, e_io)
      
      amode = IOR(MPI_MODE_CREATE, MPI_MODE_WRONLY) 
      
      call MPI_FILE_OPEN(MPI_COMM_WORLD, trim(sevt1),amode, &
        MPI_INFO_NULL,fdens,e_io)
      
      if (e_io /= MPI_SUCCESS) then
        call MPI_Error_string(e_io, emsg, elen, ierr)
        write(6,'(A,I6,2A)') 'Rank', myrank, ' MPI_File_open error: ', trim(emsg)
        write(6,'(A)') 'Path tried: '//trim(sevt1)
        call MPI_Abort(MPI_COMM_WORLD, e_io, ierr)
      endif     
         
      call MPI_File_set_size(fdens, 0_MPI_OFFSET_KIND, e_io) 
      call MPI_Barrier(MPI_COMM_WORLD, e_io)
      tempoffset=int(0,kind=MPI_OFFSET_KIND)

      sheadervtk=repeat(' ',500)
      sheadervtk=headervtk(1)
      nns=nheadervtk(1)

      if(myrank==0)call MPI_File_write_at(fdens,tempoffset,sheadervtk(1:nns),nns, &
         MPI_CHARACTER,MPI_STATUS_IGNORE,e_io)

      ioffset=vtkoffset(1)
      tempoffset=int(ioffset,kind=MPI_OFFSET_KIND)

      if(myrank==0)call MPI_File_write_at(fdens,tempoffset,int(ndatavtk(1),kind=4),1, &
         MPI_INTEGER,MPI_STATUS_IGNORE,e_io)
#ifdef PRINTHALF
      call MPI_Type_create_subarray(3,skip_gsizes,skip_lsizes,skip_myoffset, &
         MPI_ORDER_FORTRAN,MPI_UNSIGNED_SHORT,filetypesub,e_io)
#else
      call MPI_Type_create_subarray(3,skip_gsizes,skip_lsizes,skip_myoffset, &
         MPI_ORDER_FORTRAN,MPI_REAL,filetypesub,e_io)
#endif
      call MPI_Type_commit(filetypesub, e_io)

      ioffset=vtkoffset(1)+byteint
      tempoffset=int(ioffset,kind=MPI_OFFSET_KIND)
#ifdef PRINTHALF
      call MPI_File_Set_View(fdens,tempoffset,MPI_UNSIGNED_SHORT,filetypesub, &
         "native",MPI_INFO_NULL,e_io)
#else
      call MPI_File_Set_View(fdens,tempoffset,MPI_REAL,filetypesub, &
         "native",MPI_INFO_NULL,e_io)
#endif
      ! We need full local sizes: memDims
      memDims = skip_lsizes + 2*nbuffsub
      memOffs = [ nbuffsub, nbuffsub, nbuffsub ]

#ifdef PRINTHALF
      call MPI_TYPE_CREATE_SUBARRAY(3,memDims,skip_lsizes,memOffs, &
         MPI_ORDER_FORTRAN,MPI_UNSIGNED_SHORT,imemtype,e_io)
#else
      call MPI_TYPE_CREATE_SUBARRAY(3,memDims,skip_lsizes,memOffs, &
         MPI_ORDER_FORTRAN,MPI_REAL,imemtype,e_io)
#endif
      call MPI_TYPE_COMMIT(imemtype,e_io)

      call MPI_FILE_WRITE_ALL(fdens,rhoprint,1,imemtype,MPI_STATUS_IGNORE,e_io)

      ioffset=vtkoffset(1)+byteint+ndatavtk(1)
      tempoffset=int(ioffset,kind=MPI_OFFSET_KIND)

      if(myrank==0)call MPI_File_write_at(fdens,tempoffset,footervtk(1),30, &
         MPI_CHARACTER,MPI_STATUS_IGNORE,e_io)

      call MPI_Type_free(imemtype,   e_io)
      call MPI_Type_free(filetypesub,e_io)
      call MPI_File_sync(fdens, e_io)
      
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
      call MPI_FILE_CLOSE(fdens, e_io)

      if (e_io /= MPI_SUCCESS) then
        call MPI_Error_string(e_io, emsg, elen, ierr)
        write(6,'(A,I6,2A)') 'Rank', myrank, ' MPI_File_close error: ', trim(emsg)
        write(6,'(A)') 'Path tried: '//trim(sevt1)
        call MPI_Abort(MPI_COMM_WORLD, e_io, ierr)
      endif  
      
      call mpi_barrier(MPI_COMM_WORLD,e_io)
      fdens = MPI_FILE_NULL
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!velocity!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
      sevt2 = trim(dir_out) // trim(filenamevtk)//'_'//trim(namevarvtk(2))// &
         '_'//trim(write_fmtnumb(iframe)) // '.vti'
      
      if (myrank == 0) then
        inquire(file=trim(sevt2), exist=lexist)
        if (lexist) then
          call MPI_File_delete(trim(sevt2), MPI_INFO_NULL, e_io)
        endif
      endif
      call MPI_Barrier(MPI_COMM_WORLD, e_io)
      
      amode = IOR(MPI_MODE_CREATE, MPI_MODE_WRONLY) 
      
      call MPI_FILE_OPEN(MPI_COMM_WORLD,trim(sevt2),amode, &
         MPI_INFO_NULL,fvel,e_io)

      if (e_io /= MPI_SUCCESS) then
        call MPI_Error_string(e_io, emsg, elen, ierr)
        write(6,'(A,I6,2A)') 'Rank', myrank, ' MPI_File_open error: ', trim(emsg)
        write(6,'(A)') 'Path tried: '//trim(sevt2)
        call MPI_Abort(MPI_COMM_WORLD, e_io, ierr)
      endif     
         
      call MPI_File_set_size(fvel, 0_MPI_OFFSET_KIND, e_io) 
      call MPI_Barrier(MPI_COMM_WORLD, e_io)
      tempoffset=int(0,kind=MPI_OFFSET_KIND)

      sheadervtk=repeat(' ',500)
      sheadervtk=headervtk(2)
      nns=nheadervtk(2)

      if(myrank==0)call MPI_File_write_at(fvel,tempoffset,sheadervtk(1:nns),nns, &
         MPI_CHARACTER,MPI_STATUS_IGNORE,e_io)

      velglobalDims(1)=3
      velglobalDims(2:4)=skip_gsizes(1:3)
      velldims(1)=3
      velldims(2:4)=skip_lsizes(1:3)
      velmystarts(1) = 0
      velmystarts(2:4) = skip_myoffset(1:3)

      ioffset=vtkoffset(2)
      tempoffset=int(ioffset,kind=MPI_OFFSET_KIND)

      if(myrank==0)call MPI_File_write_at(fvel,tempoffset,int(ndatavtk(2),kind=4),1, &
         MPI_INTEGER,MPI_STATUS_IGNORE,e_io)

#ifdef PRINTHALF
      call MPI_Type_create_subarray(4,velglobalDims,velldims,velmystarts, &
         MPI_ORDER_FORTRAN,MPI_UNSIGNED_SHORT,filetypesubv,e_io)
#else
      call MPI_Type_create_subarray(4,velglobalDims,velldims,velmystarts, &
         MPI_ORDER_FORTRAN,MPI_REAL,filetypesubv,e_io)
#endif
      call MPI_Type_commit(filetypesubv, e_io)

      ioffset=vtkoffset(2)+byteint
      tempoffset=int(ioffset,kind=MPI_OFFSET_KIND)
#ifdef PRINTHALF
      call MPI_File_Set_View(fvel,tempoffset,MPI_UNSIGNED_SHORT,filetypesubv, &
         "native",MPI_INFO_NULL,e_io)
#else
      call MPI_File_Set_View(fvel,tempoffset,MPI_REAL,filetypesubv, &
         "native",MPI_INFO_NULL,e_io)
#endif
      ! We need full local sizes: memDims
      velmemDims(1) = vellDims(1)
      velmemDims(2:4) = vellDims(2:4) + 2*nbuffsub
      velmemOffs = [ 0, nbuffsub, nbuffsub, nbuffsub ]
#ifdef PRINTHALF
      call MPI_TYPE_CREATE_SUBARRAY(4,velmemDims,velldims,velmemOffs, &
         MPI_ORDER_FORTRAN,MPI_UNSIGNED_SHORT,imemtype,e_io)
#else
      call MPI_TYPE_CREATE_SUBARRAY(4,velmemDims,velldims,velmemOffs, &
         MPI_ORDER_FORTRAN,MPI_REAL,imemtype,e_io)
#endif
      call MPI_TYPE_COMMIT(imemtype,e_io)

      call MPI_FILE_WRITE_ALL(fvel,velprint,1,imemtype,MPI_STATUS_IGNORE,e_io)

      ioffset=vtkoffset(2)+byteint+ndatavtk(2)
      tempoffset=int(ioffset,kind=MPI_OFFSET_KIND)

      if(myrank==0)call MPI_File_write_at(fvel,tempoffset,footervtk(2),30, &
         MPI_CHARACTER,MPI_STATUS_IGNORE,e_io)
      
      call MPI_Type_free(imemtype,   e_io)
      call MPI_Type_free(filetypesub,e_io)
      call MPI_File_sync(fvel, e_io)
      
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
      call MPI_FILE_CLOSE(fvel, e_io)
      
      if (e_io /= MPI_SUCCESS) then
        call MPI_Error_string(e_io, emsg, elen, ierr)
        write(6,'(A,I6,2A)') 'Rank', myrank, ' MPI_File_close error: ', trim(emsg)
        write(6,'(A)') 'Path tried: '//trim(sevt2)
        call MPI_Abort(MPI_COMM_WORLD, e_io, ierr)
      endif  
      
      call mpi_barrier(MPI_COMM_WORLD,e_io)
      fvel = MPI_FILE_NULL
      
#if defined(TWOCOMPONENT) && defined(WRITEPRESS)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!pressure!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      sevt3 = trim(dir_out) // trim(filenamevtk)//'_'//trim(namevarvtk(3))// &
         '_'//trim(write_fmtnumb(iframe)) // '.vti'


      if (myrank == 0) then
        inquire(file=trim(sevt3), exist=lexist)
        if (lexist) then
          call MPI_File_delete(trim(sevt3), MPI_INFO_NULL, e_io)
        endif
      endif
      call MPI_Barrier(MPI_COMM_WORLD, e_io)
      
      amode = IOR(MPI_MODE_CREATE, MPI_MODE_WRONLY) 

      call MPI_FILE_OPEN(MPI_COMM_WORLD, trim(sevt3),amode, &
         MPI_INFO_NULL,fdens,e_io)
      
      if (e_io /= MPI_SUCCESS) then
        call MPI_Error_string(e_io, emsg, elen, ierr)
        write(6,'(A,I6,2A)') 'Rank', myrank, ' MPI_File_open error: ', trim(emsg)
        write(6,'(A)') 'Path tried: '//trim(sevt3)
        call MPI_Abort(MPI_COMM_WORLD, e_io, ierr)
      endif     
         
      call MPI_File_set_size(fdens, 0_MPI_OFFSET_KIND, e_io) 
      call MPI_Barrier(MPI_COMM_WORLD, e_io)
      tempoffset=int(0,kind=MPI_OFFSET_KIND)

      sheadervtk=repeat(' ',500)
      sheadervtk=headervtk(3)
      nns=nheadervtk(3)

      if(myrank==0)call MPI_File_write_at(fdens,tempoffset,sheadervtk(1:nns),nns, &
         MPI_CHARACTER,MPI_STATUS_IGNORE,e_io)

      ioffset=vtkoffset(3)
      tempoffset=int(ioffset,kind=MPI_OFFSET_KIND)

      if(myrank==0)call MPI_File_write_at(fdens,tempoffset,int(ndatavtk(3),kind=4),1, &
         MPI_INTEGER,MPI_STATUS_IGNORE,e_io)
#ifdef PRINTHALF
      call MPI_Type_create_subarray(3,skip_gsizes,skip_lsizes,skip_myoffset, &
         MPI_ORDER_FORTRAN,MPI_UNSIGNED_SHORT,filetypesub,e_io)
#else
      call MPI_Type_create_subarray(3,skip_gsizes,skip_lsizes,skip_myoffset, &
         MPI_ORDER_FORTRAN,MPI_REAL,filetypesub,e_io)
#endif
      call MPI_Type_commit(filetypesub, e_io)

      ioffset=vtkoffset(3)+byteint
      tempoffset=int(ioffset,kind=MPI_OFFSET_KIND)
      
#ifdef PRINTHALF
      call MPI_File_Set_View(fdens,tempoffset,MPI_UNSIGNED_SHORT,filetypesub, &
         "native",MPI_INFO_NULL,e_io)
#else
      call MPI_File_Set_View(fdens,tempoffset,MPI_REAL,filetypesub, &
         "native",MPI_INFO_NULL,e_io)
#endif
      ! We need full local sizes: memDims
      memDims = skip_lsizes + 2*nbuffsub
      memOffs = [ nbuffsub, nbuffsub, nbuffsub ]

#ifdef PRINTHALF
      call MPI_TYPE_CREATE_SUBARRAY(3,memDims,skip_lsizes,memOffs, &
         MPI_ORDER_FORTRAN,MPI_UNSIGNED_SHORT,imemtype,e_io)
#else
      call MPI_TYPE_CREATE_SUBARRAY(3,memDims,lsizes,memOffs, &
         MPI_ORDER_FORTRAN,MPI_REAL,imemtype,e_io)
#endif
      call MPI_TYPE_COMMIT(imemtype,e_io)

      call MPI_FILE_WRITE_ALL(fdens,pressprint,1,imemtype,MPI_STATUS_IGNORE,e_io)

      ioffset=vtkoffset(3)+byteint+ndatavtk(3)
      tempoffset=int(ioffset,kind=MPI_OFFSET_KIND)

      if(myrank==0)call MPI_File_write_at(fdens,tempoffset,footervtk(3),30, &
         MPI_CHARACTER,MPI_STATUS_IGNORE,e_io)

      
      call MPI_Type_free(imemtype,   e_io)
      call MPI_Type_free(filetypesub,e_io)
      call MPI_File_sync(fdens, e_io)
      
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
      call MPI_FILE_CLOSE(fdens, e_io)
      
      if (e_io /= MPI_SUCCESS) then
        call MPI_Error_string(e_io, emsg, elen, ierr)
        write(6,'(A,I6,2A)') 'Rank', myrank, ' MPI_File_close error: ', trim(emsg)
        write(6,'(A)') 'Path tried: '//trim(sevt3)
        call MPI_Abort(MPI_COMM_WORLD, e_io, ierr)
      endif  
      
      call mpi_barrier(MPI_COMM_WORLD,e_io)
      fdens = MPI_FILE_NULL

#endif

#endif
      return

   end subroutine write_file_vtk_par

   subroutine write_restart_parallel_1c(iframe,iframe2D,e_io,hfields_s)

!***********************************************************************
!
!     LBsoft subroutine for opening the vtk legacy file
!     in parallel IO
!
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2019
!
!***********************************************************************

      implicit none

      integer(kind=4), intent(in) :: iframe,iframe2D
      real(kind=strdb), allocatable, dimension(:,:,:,:,:) :: hfields_s
      integer, intent(out) :: e_io
#ifdef MPI
      integer(kind=MPI_OFFSET_KIND) :: tempoffset,offset_final
#endif

      integer :: ioffset
      character(1), parameter :: end_rec = char(10)
      integer, parameter :: bytechar=kind(end_rec)
      integer, parameter :: byteint = 4
      integer, parameter :: byter4  = 4
      integer, parameter :: byter8  = 8
      integer, parameter :: nbuffsub = 0
      integer :: filetypesub,imemtype,filetypesubv,ierr
      logical :: lexist
      integer, dimension(3) :: memDims,memOffs
      integer :: elen,amode
      integer :: xblock,yblock,zblock,myblock
      integer :: ii,jj,kk

#ifdef MPI
      integer :: fdens=MPI_FILE_NULL
      character(len=MPI_MAX_ERROR_STRING) :: emsg
      
      call MPI_Barrier(MPI_COMM_WORLD, e_io)
      
      sevt1 = trim(dir_out) // 'restart.raw'
      
      if (myrank == 0) then
        inquire(file=trim(sevt1), exist=lexist)
        if (lexist) then
          call MPI_File_delete(trim(sevt1), MPI_INFO_NULL, e_io)
          !write(6,*)'replace ',trim(sevt1)
        endif
      endif
      call MPI_Barrier(MPI_COMM_WORLD, e_io)
      
      amode = IOR(MPI_MODE_CREATE, MPI_MODE_WRONLY) 
      
      call MPI_File_open(MPI_COMM_WORLD, trim(sevt1),amode, MPI_INFO_NULL, fdens, e_io)   
         
      if (e_io /= MPI_SUCCESS) then
        call MPI_Error_string(e_io, emsg, elen, ierr)
        write(6,'(A,I6,2A)') 'Rank', myrank, ' MPI_File_open error: ', trim(emsg)
        write(6,'(A)') 'Path tried: '//trim(sevt1)
        call MPI_Abort(MPI_COMM_WORLD, e_io, ierr)
      endif    
      
      call MPI_File_set_size(fdens, 0_MPI_OFFSET_KIND, e_io) 
      call MPI_Barrier(MPI_COMM_WORLD, e_io)
      tempoffset=int(0,kind=MPI_OFFSET_KIND)


      call MPI_Type_create_subarray(3,gsizes,lsizes,myoffset, &
         MPI_ORDER_FORTRAN,MPI_REAL,filetypesub,e_io)

      call MPI_Type_commit(filetypesub, e_io)

      call MPI_File_Set_View(fdens,tempoffset,MPI_REAL,filetypesub, &
         "native",MPI_INFO_NULL,e_io)
      
      ! We need full local sizes: memDims
      memDims = lsizes + 2*nbuffsub
      memOffs = [ nbuffsub, nbuffsub, nbuffsub ]

      call MPI_TYPE_CREATE_SUBARRAY(3,memDims,lsizes,memOffs, &
         MPI_ORDER_FORTRAN,MPI_REAL,imemtype,e_io)

      call MPI_TYPE_COMMIT(imemtype,e_io)

      
      do k=1,nz
  	     zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
         do j=1,ny
            yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
            do i=1,nx
               xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
               myblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
               ii=i-xblock*TILE_DIMx+2*TILE_DIMx
               jj=j-yblock*TILE_DIMy+2*TILE_DIMy
               kk=k-zblock*TILE_DIMz+2*TILE_DIMz  
               arr_3d(i,j,k)=real(hfields_s(ii,jj,kk,1,myblock),kind=4)
             enddo
         enddo
      enddo
      call MPI_FILE_WRITE_ALL(fdens,arr_3d,1,imemtype,MPI_STATUS_IGNORE,e_io)

      
      do k=1,nz
  	     zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
         do j=1,ny
            yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
            do i=1,nx
               xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
               myblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
               ii=i-xblock*TILE_DIMx+2*TILE_DIMx
               jj=j-yblock*TILE_DIMy+2*TILE_DIMy
               kk=k-zblock*TILE_DIMz+2*TILE_DIMz  
               arr_3d(i,j,k)=real(hfields_s(ii,jj,kk,2,myblock),kind=4)
             enddo
         enddo
      enddo 
      call MPI_FILE_WRITE_ALL(fdens,arr_3d,1,imemtype,MPI_STATUS_IGNORE,e_io)
           
      
      do k=1,nz
  	     zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
         do j=1,ny
            yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
            do i=1,nx
               xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
               myblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
               ii=i-xblock*TILE_DIMx+2*TILE_DIMx
               jj=j-yblock*TILE_DIMy+2*TILE_DIMy
               kk=k-zblock*TILE_DIMz+2*TILE_DIMz  
               arr_3d(i,j,k)=real(hfields_s(ii,jj,kk,3,myblock),kind=4)
             enddo
         enddo
      enddo 
      call MPI_FILE_WRITE_ALL(fdens,arr_3d,1,imemtype,MPI_STATUS_IGNORE,e_io)
      
      
      do k=1,nz
  	     zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
         do j=1,ny
            yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
            do i=1,nx
               xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
               myblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
               ii=i-xblock*TILE_DIMx+2*TILE_DIMx
               jj=j-yblock*TILE_DIMy+2*TILE_DIMy
               kk=k-zblock*TILE_DIMz+2*TILE_DIMz  
               arr_3d(i,j,k)=real(hfields_s(ii,jj,kk,4,myblock),kind=4)
             enddo
         enddo
      enddo 
      call MPI_FILE_WRITE_ALL(fdens,arr_3d,1,imemtype,MPI_STATUS_IGNORE,e_io)

      
      do k=1,nz
  	     zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
         do j=1,ny
            yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
            do i=1,nx
               xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
               myblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
               ii=i-xblock*TILE_DIMx+2*TILE_DIMx
               jj=j-yblock*TILE_DIMy+2*TILE_DIMy
               kk=k-zblock*TILE_DIMz+2*TILE_DIMz  
               arr_3d(i,j,k)=real(hfields_s(ii,jj,kk,5,myblock),kind=4)
             enddo
         enddo
      enddo 
      call MPI_FILE_WRITE_ALL(fdens,arr_3d,1,imemtype,MPI_STATUS_IGNORE,e_io)
      
      
      do k=1,nz
  	     zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
         do j=1,ny
            yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
            do i=1,nx
               xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
               myblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
               ii=i-xblock*TILE_DIMx+2*TILE_DIMx
               jj=j-yblock*TILE_DIMy+2*TILE_DIMy
               kk=k-zblock*TILE_DIMz+2*TILE_DIMz  
               arr_3d(i,j,k)=real(hfields_s(ii,jj,kk,8,myblock),kind=4)
             enddo
         enddo
      enddo 
      call MPI_FILE_WRITE_ALL(fdens,arr_3d,1,imemtype,MPI_STATUS_IGNORE,e_io)
      
      
      do k=1,nz
  	     zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
         do j=1,ny
            yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
            do i=1,nx
               xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
               myblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
               ii=i-xblock*TILE_DIMx+2*TILE_DIMx
               jj=j-yblock*TILE_DIMy+2*TILE_DIMy
               kk=k-zblock*TILE_DIMz+2*TILE_DIMz  
               arr_3d(i,j,k)=real(hfields_s(ii,jj,kk,9,myblock),kind=4)
             enddo
         enddo
      enddo 
      call MPI_FILE_WRITE_ALL(fdens,arr_3d,1,imemtype,MPI_STATUS_IGNORE,e_io)
      
       
      do k=1,nz
  	     zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
         do j=1,ny
            yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
            do i=1,nx
               xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
               myblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
               ii=i-xblock*TILE_DIMx+2*TILE_DIMx
               jj=j-yblock*TILE_DIMy+2*TILE_DIMy
               kk=k-zblock*TILE_DIMz+2*TILE_DIMz  
               arr_3d(i,j,k)=real(hfields_s(ii,jj,kk,6,myblock),kind=4)
             enddo
         enddo
      enddo 
      call MPI_FILE_WRITE_ALL(fdens,arr_3d,1,imemtype,MPI_STATUS_IGNORE,e_io)

      
      do k=1,nz
  	     zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
         do j=1,ny
            yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
            do i=1,nx
               xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
               myblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
               ii=i-xblock*TILE_DIMx+2*TILE_DIMx
               jj=j-yblock*TILE_DIMy+2*TILE_DIMy
               kk=k-zblock*TILE_DIMz+2*TILE_DIMz  
               arr_3d(i,j,k)=real(hfields_s(ii,jj,kk,10,myblock),kind=4)
             enddo
         enddo
      enddo 
      call MPI_FILE_WRITE_ALL(fdens,arr_3d,1,imemtype,MPI_STATUS_IGNORE,e_io)
      
       
      do k=1,nz
  	     zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
         do j=1,ny
            yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
            do i=1,nx
               xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
               myblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
               ii=i-xblock*TILE_DIMx+2*TILE_DIMx
               jj=j-yblock*TILE_DIMy+2*TILE_DIMy
               kk=k-zblock*TILE_DIMz+2*TILE_DIMz  
               arr_3d(i,j,k)=real(hfields_s(ii,jj,kk,7,myblock),kind=4)
             enddo
         enddo
      enddo 
      call MPI_FILE_WRITE_ALL(fdens,arr_3d,1,imemtype,MPI_STATUS_IGNORE,e_io)
      
      call MPI_Type_free(imemtype,   e_io)
      call MPI_Type_free(filetypesub,e_io)
      call MPI_File_sync(fdens, e_io)
      
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
      call MPI_FILE_CLOSE(fdens, e_io)
      
      if (e_io /= MPI_SUCCESS) then
        call MPI_Error_string(e_io, emsg, elen, ierr)
        write(6,'(A,I6,2A)') 'Rank', myrank, ' MPI_File_close error: ', trim(emsg)
        write(6,'(A)') 'Path tried: '//trim(sevt1)
        call MPI_Abort(MPI_COMM_WORLD, e_io, ierr)
      endif  
      
      call mpi_barrier(MPI_COMM_WORLD,e_io)
      fdens = MPI_FILE_NULL
      if (myrank == 0) then
        call MPI_FILE_OPEN(MPI_COMM_SELF, trim(sevt1), MPI_MODE_WRONLY, MPI_INFO_NULL, fdens, e_io)
        if (e_io /= MPI_SUCCESS) then
          call MPI_Error_string(e_io, emsg, elen, ierr)
          write(6,'(A,I6,2A)') 'Rank', myrank, ' MPI_File_open error: ', trim(emsg)
          write(6,'(A)') 'Path tried: '//trim(sevt1)
          call MPI_Abort(MPI_COMM_WORLD, e_io, ierr)
        endif    
        call MPI_File_set_view(fdens, 0_MPI_OFFSET_KIND, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL, e_io) 
        offset_final = int(lx * ly * lz * 10 * db, kind=MPI_OFFSET_KIND)
        call MPI_FILE_WRITE_AT(fdens, offset_final, iframe, 1, MPI_INTEGER, MPI_STATUS_IGNORE, e_io)
        call MPI_FILE_WRITE_AT(fdens, offset_final + byteint, iframe2D, 1, MPI_INTEGER, MPI_STATUS_IGNORE, e_io)
        call MPI_File_sync(fdens, e_io)
        call MPI_FILE_CLOSE(fdens, e_io)
      endif
      fdens = MPI_FILE_NULL
      call mpi_barrier(MPI_COMM_WORLD,e_io)
#endif
      return

   end subroutine write_restart_parallel_1c

   subroutine read_restart_parallel_1c(iframe,iframe2D,e_io,hfields_s)

!***********************************************************************
!
!     LBsoft subroutine for opening the vtk legacy file
!     in parallel IO
!
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2019
!
!***********************************************************************

      implicit none

      integer(kind=4), intent(out) ::  iframe,iframe2D
      real(kind=strdb), allocatable, dimension(:,:,:,:,:) :: hfields_s
      integer, intent(out) :: e_io
#ifdef MPI
      integer(kind=MPI_OFFSET_KIND) :: tempoffset,offset_final
#endif

      integer :: ioffset
      character(1), parameter :: end_rec = char(10)
      integer, parameter :: bytechar=kind(end_rec)
      integer, parameter :: byteint = 4
      integer, parameter :: byter4  = 4
      integer, parameter :: byter8  = 8
      integer, parameter :: nbuffsub = 0
      integer :: filetypesub,imemtype,filetypesubv,ierr
      logical :: lexist
      integer, dimension(3) :: memDims,memOffs
      integer :: elen,amode

#ifdef MPI
      integer :: fdens=MPI_FILE_NULL
      character(len=MPI_MAX_ERROR_STRING) :: emsg
      logical :: file_exists
      integer :: ii,jj,kk,myblock,xblock,yblock,zblock

      sevt1 = trim(dir_out) // 'restart.raw'
      
      inquire(file=trim(sevt1), exist=file_exists)
      if (.not. file_exists) then
         call doerror(6,'ERROR: unable to find the file: '//trim(sevt1))
      end if    

      call MPI_FILE_OPEN(MPI_COMM_WORLD, trim(sevt1), &
         MPI_MODE_RDONLY, &
         MPI_INFO_NULL,fdens,e_io)
         
      if (e_io /= MPI_SUCCESS) then
        call MPI_Error_string(e_io, emsg, elen, ierr)
        write(6,'(A,I6,2A)') 'Rank', myrank, ' MPI_File_open error: ', trim(emsg)
        write(6,'(A)') 'Path tried: '//trim(sevt1)
        call MPI_Abort(MPI_COMM_WORLD, e_io, ierr)
      endif   
      call MPI_Barrier(MPI_COMM_WORLD, e_io)

      call MPI_Type_create_subarray(3,gsizes,lsizes,myoffset, &
         MPI_ORDER_FORTRAN,MPI_REAL,filetypesub,e_io)

      call MPI_Type_commit(filetypesub, e_io)
      
      tempoffset=int(0,kind=MPI_OFFSET_KIND)
      
      call MPI_File_Set_View(fdens,tempoffset,MPI_REAL,filetypesub, &
         "native",MPI_INFO_NULL,e_io)
         
      call MPI_Barrier(MPI_COMM_WORLD, e_io)
      
      ! We need full local sizes: memDims
      memDims = lsizes + 2*nbuffsub
      memOffs = [ nbuffsub, nbuffsub, nbuffsub ]

      call MPI_TYPE_CREATE_SUBARRAY(3,memDims,lsizes,memOffs, &
         MPI_ORDER_FORTRAN,MPI_REAL,imemtype,e_io)

      call MPI_TYPE_COMMIT(imemtype,e_io)

      call MPI_FILE_READ_ALL(fdens,arr_3d,1,imemtype,MPI_STATUS_IGNORE,e_io)
      do k=1,nz
  	     zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
         do j=1,ny
            yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
            do i=1,nx
               xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
               myblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
               ii=i-xblock*TILE_DIMx+2*TILE_DIMx
               jj=j-yblock*TILE_DIMy+2*TILE_DIMy
               kk=k-zblock*TILE_DIMz+2*TILE_DIMz  
               hfields_s(ii,jj,kk,1,myblock)=real(arr_3d(i,j,k),kind=strdb)
             enddo
         enddo
      enddo   
      

      call MPI_FILE_READ_ALL(fdens,arr_3d,1,imemtype,MPI_STATUS_IGNORE,e_io)
      do k=1,nz
  	     zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
         do j=1,ny
            yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
            do i=1,nx
               xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
               myblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
               ii=i-xblock*TILE_DIMx+2*TILE_DIMx
               jj=j-yblock*TILE_DIMy+2*TILE_DIMy
               kk=k-zblock*TILE_DIMz+2*TILE_DIMz  
               hfields_s(ii,jj,kk,2,myblock)=real(arr_3d(i,j,k) ,kind=strdb)  
             enddo
         enddo
      enddo   
       
      
      call MPI_FILE_READ_ALL(fdens,arr_3d,1,imemtype,MPI_STATUS_IGNORE,e_io)
      do k=1,nz
  	     zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
         do j=1,ny
            yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
            do i=1,nx
               xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
               myblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
               ii=i-xblock*TILE_DIMx+2*TILE_DIMx
               jj=j-yblock*TILE_DIMy+2*TILE_DIMy
               kk=k-zblock*TILE_DIMz+2*TILE_DIMz  
               hfields_s(ii,jj,kk,3,myblock)=real(arr_3d(i,j,k),kind=strdb) 
             enddo
         enddo
      enddo   
      
      
      call MPI_FILE_READ_ALL(fdens,arr_3d,1,imemtype,MPI_STATUS_IGNORE,e_io)
      do k=1,nz
  	     zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
         do j=1,ny
            yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
            do i=1,nx
               xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
               myblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
               ii=i-xblock*TILE_DIMx+2*TILE_DIMx
               jj=j-yblock*TILE_DIMy+2*TILE_DIMy
               kk=k-zblock*TILE_DIMz+2*TILE_DIMz  
               hfields_s(ii,jj,kk,4,myblock)=real(arr_3d(i,j,k),kind=strdb)
             enddo
         enddo
      enddo   
      

      call MPI_FILE_READ_ALL(fdens,arr_3d,1,imemtype,MPI_STATUS_IGNORE,e_io)
      do k=1,nz
  	     zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
         do j=1,ny
            yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
            do i=1,nx
               xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
               myblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
               ii=i-xblock*TILE_DIMx+2*TILE_DIMx
               jj=j-yblock*TILE_DIMy+2*TILE_DIMy
               kk=k-zblock*TILE_DIMz+2*TILE_DIMz  
               hfields_s(ii,jj,kk,5,myblock)=real(arr_3d(i,j,k),kind=strdb)
             enddo
         enddo
      enddo   
      
      
      call MPI_FILE_READ_ALL(fdens,arr_3d,1,imemtype,MPI_STATUS_IGNORE,e_io)
      do k=1,nz
  	     zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
         do j=1,ny
            yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
            do i=1,nx
               xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
               myblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
               ii=i-xblock*TILE_DIMx+2*TILE_DIMx
               jj=j-yblock*TILE_DIMy+2*TILE_DIMy
               kk=k-zblock*TILE_DIMz+2*TILE_DIMz  
               hfields_s(ii,jj,kk,8,myblock)=real(arr_3d(i,j,k),kind=strdb)
             enddo
         enddo
      enddo   
      
      
      call MPI_FILE_READ_ALL(fdens,arr_3d,1,imemtype,MPI_STATUS_IGNORE,e_io)
      do k=1,nz
  	     zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
         do j=1,ny
            yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
            do i=1,nx
               xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
               myblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
               ii=i-xblock*TILE_DIMx+2*TILE_DIMx
               jj=j-yblock*TILE_DIMy+2*TILE_DIMy
               kk=k-zblock*TILE_DIMz+2*TILE_DIMz  
               hfields_s(ii,jj,kk,9,myblock)=real(arr_3d(i,j,k),kind=strdb) 
             enddo
         enddo
      enddo   
      
      
      call MPI_FILE_READ_ALL(fdens,arr_3d,1,imemtype,MPI_STATUS_IGNORE,e_io)
      do k=1,nz
  	     zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
         do j=1,ny
            yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
            do i=1,nx
               xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
               myblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
               ii=i-xblock*TILE_DIMx+2*TILE_DIMx
               jj=j-yblock*TILE_DIMy+2*TILE_DIMy
               kk=k-zblock*TILE_DIMz+2*TILE_DIMz  
               hfields_s(ii,jj,kk,6,myblock)=real(arr_3d(i,j,k),kind=strdb)
             enddo
         enddo
      enddo   
      
      
      call MPI_FILE_READ_ALL(fdens,arr_3d,1,imemtype,MPI_STATUS_IGNORE,e_io)
      do k=1,nz
  	     zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
         do j=1,ny
            yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
            do i=1,nx
               xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
               myblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
               ii=i-xblock*TILE_DIMx+2*TILE_DIMx
               jj=j-yblock*TILE_DIMy+2*TILE_DIMy
               kk=k-zblock*TILE_DIMz+2*TILE_DIMz  
               hfields_s(ii,jj,kk,10,myblock)=real(arr_3d(i,j,k),kind=strdb)
             enddo
         enddo
      enddo   
      
      
      call MPI_FILE_READ_ALL(fdens,arr_3d,1,imemtype,MPI_STATUS_IGNORE,e_io)
      do k=1,nz
  	     zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
         do j=1,ny
            yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
            do i=1,nx
               xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
               myblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
               ii=i-xblock*TILE_DIMx+2*TILE_DIMx
               jj=j-yblock*TILE_DIMy+2*TILE_DIMy
               kk=k-zblock*TILE_DIMz+2*TILE_DIMz  
               hfields_s(ii,jj,kk,7,myblock)=real(arr_3d(i,j,k),kind=strdb)  
             enddo
         enddo
      enddo   
      
      call MPI_Type_free(imemtype, e_io)
      call MPI_Type_free(filetypesub, e_io)
      
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
      call MPI_FILE_CLOSE(fdens, e_io)
      if (e_io /= MPI_SUCCESS) then
        call MPI_Error_string(e_io, emsg, elen, ierr)
        write(6,'(A,I6,2A)') 'Rank', myrank, ' MPI_File_close error: ', trim(emsg)
        write(6,'(A)') 'Path tried: '//trim(sevt1)
        call MPI_Abort(MPI_COMM_WORLD, e_io, ierr)
      endif     
      fdens = MPI_FILE_NULL
      call mpi_barrier(MPI_COMM_WORLD,e_io)
      
      if(myrank == 0) then
        call MPI_FILE_OPEN(MPI_COMM_SELF, trim(sevt1), MPI_MODE_RDONLY, MPI_INFO_NULL, fdens, e_io)
        if (e_io /= MPI_SUCCESS) then
          call MPI_Error_string(e_io, emsg, elen, ierr)
          write(6,'(A,I6,2A)') 'Rank', myrank, ' MPI_File_open error: ', trim(emsg)
          write(6,'(A)') 'Path tried: '//trim(sevt1)
          call MPI_Abort(MPI_COMM_WORLD, e_io, ierr)
        endif  
        call MPI_File_set_view(fdens, 0_MPI_OFFSET_KIND, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL, e_io)

        offset_final = int(lx * ly * lz * 10 * db, kind=MPI_OFFSET_KIND)
        call MPI_FILE_READ_AT(fdens, offset_final, iframe, 1, MPI_INTEGER, MPI_STATUS_IGNORE, e_io)
        call MPI_FILE_READ_AT(fdens, offset_final + byteint, iframe2D, 1, MPI_INTEGER, MPI_STATUS_IGNORE, e_io)

        call MPI_FILE_CLOSE(fdens, e_io)
      endif
      fdens = MPI_FILE_NULL
      call mpi_barrier(MPI_COMM_WORLD,e_io)

#endif
      return

   end subroutine read_restart_parallel_1c

   subroutine write_restart_parallel_2c(iframe,iframe2D,e_io,hfields_s,phifields_s)

!***********************************************************************
!
!     LBsoft subroutine for opening the vtk legacy file
!     in parallel IO
!
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2019
!
!***********************************************************************

      implicit none

      integer(kind=4), intent(in) :: iframe,iframe2D
      real(kind=strdb), allocatable, dimension(:,:,:,:,:) :: hfields_s,phifields_s
      integer, intent(out) :: e_io
#ifdef MPI
      integer(kind=MPI_OFFSET_KIND) :: tempoffset,offset_final
#endif

      integer :: ioffset
      character(1), parameter :: end_rec = char(10)
      integer, parameter :: bytechar=kind(end_rec)
      integer, parameter :: byteint = 4
      integer, parameter :: byter4  = 4
      integer, parameter :: byter8  = 8
      integer, parameter :: nbuffsub = 0
      integer :: filetypesub,imemtype,filetypesubv,ierr
      logical :: lexist
      integer, dimension(3) :: memDims,memOffs
      integer :: elen,amode
      integer :: xblock,yblock,zblock,myblock
      integer :: ii,jj,kk

#ifdef MPI
      integer :: fdens=MPI_FILE_NULL
      character(len=MPI_MAX_ERROR_STRING) :: emsg
      
      call MPI_Barrier(MPI_COMM_WORLD, e_io)
      
      sevt1 = trim(dir_out) // 'restart.raw'
      
      if (myrank == 0) then
        inquire(file=trim(sevt1), exist=lexist)
        if (lexist) then
          call MPI_File_delete(trim(sevt1), MPI_INFO_NULL, e_io)
          !write(6,*)'replace ',trim(sevt1)
        endif
      endif
      call MPI_Barrier(MPI_COMM_WORLD, e_io)
      
      amode = IOR(MPI_MODE_CREATE, MPI_MODE_WRONLY) 
      
      call MPI_File_open(MPI_COMM_WORLD, trim(sevt1),amode, MPI_INFO_NULL, fdens, e_io)   
         
      if (e_io /= MPI_SUCCESS) then
        call MPI_Error_string(e_io, emsg, elen, ierr)
        write(6,'(A,I6,2A)') 'Rank', myrank, ' MPI_File_open error: ', trim(emsg)
        write(6,'(A)') 'Path tried: '//trim(sevt1)
        call MPI_Abort(MPI_COMM_WORLD, e_io, ierr)
      endif     
         
      call MPI_File_set_size(fdens, 0_MPI_OFFSET_KIND, e_io) 
      call MPI_Barrier(MPI_COMM_WORLD, e_io)
      tempoffset=int(0,kind=MPI_OFFSET_KIND)


      call MPI_Type_create_subarray(3,gsizes,lsizes,myoffset, &
         MPI_ORDER_FORTRAN,MPI_REAL,filetypesub,e_io)

      call MPI_Type_commit(filetypesub, e_io)

      call MPI_File_Set_View(fdens,tempoffset,MPI_REAL,filetypesub, &
         "native",MPI_INFO_NULL,e_io)
      ! We need full local sizes: memDims
      memDims = lsizes + 2*nbuffsub
      memOffs = [ nbuffsub, nbuffsub, nbuffsub ]

      call MPI_TYPE_CREATE_SUBARRAY(3,memDims,lsizes,memOffs, &
         MPI_ORDER_FORTRAN,MPI_REAL,imemtype,e_io)

      call MPI_TYPE_COMMIT(imemtype,e_io)

      
      do k=1,nz
  	     zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
         do j=1,ny
            yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
            do i=1,nx
               xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
               myblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
               ii=i-xblock*TILE_DIMx+2*TILE_DIMx
               jj=j-yblock*TILE_DIMy+2*TILE_DIMy
               kk=k-zblock*TILE_DIMz+2*TILE_DIMz  
               arr_3d(i,j,k)=real(hfields_s(ii,jj,kk,1,myblock),kind=4)
             enddo
         enddo
      enddo 
      call MPI_FILE_WRITE_ALL(fdens,arr_3d,1,imemtype,MPI_STATUS_IGNORE,e_io)

      
      do k=1,nz
  	     zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
         do j=1,ny
            yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
            do i=1,nx
               xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
               myblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
               ii=i-xblock*TILE_DIMx+2*TILE_DIMx
               jj=j-yblock*TILE_DIMy+2*TILE_DIMy
               kk=k-zblock*TILE_DIMz+2*TILE_DIMz  
               arr_3d(i,j,k)=real(hfields_s(ii,jj,kk,2,myblock),kind=4)
             enddo
         enddo
      enddo 
      call MPI_FILE_WRITE_ALL(fdens,arr_3d,1,imemtype,MPI_STATUS_IGNORE,e_io)
      
      
      do k=1,nz
  	     zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
         do j=1,ny
            yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
            do i=1,nx
               xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
               myblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
               ii=i-xblock*TILE_DIMx+2*TILE_DIMx
               jj=j-yblock*TILE_DIMy+2*TILE_DIMy
               kk=k-zblock*TILE_DIMz+2*TILE_DIMz  
               arr_3d(i,j,k)=real(hfields_s(ii,jj,kk,3,myblock),kind=4)
             enddo
         enddo
      enddo 
      call MPI_FILE_WRITE_ALL(fdens,arr_3d,1,imemtype,MPI_STATUS_IGNORE,e_io)
      
      
      do k=1,nz
  	     zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
         do j=1,ny
            yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
            do i=1,nx
               xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
               myblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
               ii=i-xblock*TILE_DIMx+2*TILE_DIMx
               jj=j-yblock*TILE_DIMy+2*TILE_DIMy
               kk=k-zblock*TILE_DIMz+2*TILE_DIMz  
               arr_3d(i,j,k)=real(hfields_s(ii,jj,kk,4,myblock),kind=4)
             enddo
         enddo
      enddo 
      call MPI_FILE_WRITE_ALL(fdens,arr_3d,1,imemtype,MPI_STATUS_IGNORE,e_io)

      
      do k=1,nz
  	     zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
         do j=1,ny
            yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
            do i=1,nx
               xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
               myblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
               ii=i-xblock*TILE_DIMx+2*TILE_DIMx
               jj=j-yblock*TILE_DIMy+2*TILE_DIMy
               kk=k-zblock*TILE_DIMz+2*TILE_DIMz  
               arr_3d(i,j,k)=real(hfields_s(ii,jj,kk,5,myblock),kind=4)
             enddo
         enddo
      enddo 
      call MPI_FILE_WRITE_ALL(fdens,arr_3d,1,imemtype,MPI_STATUS_IGNORE,e_io)
      
      
      do k=1,nz
  	     zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
         do j=1,ny
            yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
            do i=1,nx
               xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
               myblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
               ii=i-xblock*TILE_DIMx+2*TILE_DIMx
               jj=j-yblock*TILE_DIMy+2*TILE_DIMy
               kk=k-zblock*TILE_DIMz+2*TILE_DIMz  
               arr_3d(i,j,k)=real(hfields_s(ii,jj,kk,8,myblock),kind=4)
             enddo
         enddo
      enddo 
      call MPI_FILE_WRITE_ALL(fdens,arr_3d,1,imemtype,MPI_STATUS_IGNORE,e_io)
      
      
      do k=1,nz
  	     zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
         do j=1,ny
            yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
            do i=1,nx
               xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
               myblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
               ii=i-xblock*TILE_DIMx+2*TILE_DIMx
               jj=j-yblock*TILE_DIMy+2*TILE_DIMy
               kk=k-zblock*TILE_DIMz+2*TILE_DIMz  
               arr_3d(i,j,k)=real(hfields_s(ii,jj,kk,9,myblock),kind=4)
             enddo
         enddo
      enddo 
      call MPI_FILE_WRITE_ALL(fdens,arr_3d,1,imemtype,MPI_STATUS_IGNORE,e_io)
      
      
      do k=1,nz
  	     zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
         do j=1,ny
            yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
            do i=1,nx
               xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
               myblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
               ii=i-xblock*TILE_DIMx+2*TILE_DIMx
               jj=j-yblock*TILE_DIMy+2*TILE_DIMy
               kk=k-zblock*TILE_DIMz+2*TILE_DIMz  
               arr_3d(i,j,k)=real(hfields_s(ii,jj,kk,6,myblock),kind=4)
             enddo
         enddo
      enddo 
      call MPI_FILE_WRITE_ALL(fdens,arr_3d,1,imemtype,MPI_STATUS_IGNORE,e_io)
      
      
      do k=1,nz
  	     zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
         do j=1,ny
            yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
            do i=1,nx
               xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
               myblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
               ii=i-xblock*TILE_DIMx+2*TILE_DIMx
               jj=j-yblock*TILE_DIMy+2*TILE_DIMy
               kk=k-zblock*TILE_DIMz+2*TILE_DIMz  
               arr_3d(i,j,k)=real(hfields_s(ii,jj,kk,10,myblock),kind=4)
             enddo
         enddo
      enddo 
      call MPI_FILE_WRITE_ALL(fdens,arr_3d,1,imemtype,MPI_STATUS_IGNORE,e_io)
      
      
      do k=1,nz
  	     zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
         do j=1,ny
            yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
            do i=1,nx
               xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
               myblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
               ii=i-xblock*TILE_DIMx+2*TILE_DIMx
               jj=j-yblock*TILE_DIMy+2*TILE_DIMy
               kk=k-zblock*TILE_DIMz+2*TILE_DIMz
               arr_3d(i,j,k)=real(hfields_s(ii,jj,kk,7,myblock),kind=4)
             enddo
         enddo
      enddo
      call MPI_FILE_WRITE_ALL(fdens,arr_3d,1,imemtype,MPI_STATUS_IGNORE,e_io)
      
      
      do k=1,nz
  	     zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
         do j=1,ny
            yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
            do i=1,nx
               xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
               myblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
               ii=i-xblock*TILE_DIMx+2*TILE_DIMx
               jj=j-yblock*TILE_DIMy+2*TILE_DIMy
               kk=k-zblock*TILE_DIMz+2*TILE_DIMz  
               arr_3d(i,j,k)=real(phifields_s(ii,jj,kk,1,myblock),kind=4)
             enddo
         enddo
      enddo
      call MPI_FILE_WRITE_ALL(fdens,arr_3d,1,imemtype,MPI_STATUS_IGNORE,e_io)
      
      call MPI_Type_free(imemtype,   e_io)
      call MPI_Type_free(filetypesub,e_io)
      call MPI_File_sync(fdens, e_io)
      
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
      call MPI_FILE_CLOSE(fdens, e_io)
      
      if (e_io /= MPI_SUCCESS) then
        call MPI_Error_string(e_io, emsg, elen, ierr)
        write(6,'(A,I6,2A)') 'Rank', myrank, ' MPI_File_close error: ', trim(emsg)
        write(6,'(A)') 'Path tried: '//trim(sevt1)
        call MPI_Abort(MPI_COMM_WORLD, e_io, ierr)
      endif  
      
      call mpi_barrier(MPI_COMM_WORLD,e_io)
      fdens = MPI_FILE_NULL
      
      if (myrank == 0) then

        call MPI_FILE_OPEN(MPI_COMM_SELF, trim(sevt1), MPI_MODE_WRONLY, MPI_INFO_NULL, fdens, e_io)
        if (e_io /= MPI_SUCCESS) then
          call MPI_Error_string(e_io, emsg, elen, ierr)
          write(6,'(A,I6,2A)') 'Rank', myrank, ' MPI_File_open error: ', trim(emsg)
          write(6,'(A)') 'Path tried: '//trim(sevt1)
          call MPI_Abort(MPI_COMM_WORLD, e_io, ierr)
        endif    
        call MPI_File_set_view(fdens, 0_MPI_OFFSET_KIND, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL, e_io) 
        offset_final = int(lx * ly * lz * 11 * db, kind=MPI_OFFSET_KIND)
        call MPI_FILE_WRITE_AT(fdens, offset_final, iframe, 1, MPI_INTEGER, MPI_STATUS_IGNORE, e_io)
        call MPI_FILE_WRITE_AT(fdens, offset_final + byteint, iframe2D, 1, MPI_INTEGER, MPI_STATUS_IGNORE, e_io)
        call MPI_File_sync(fdens, e_io)
        call MPI_FILE_CLOSE(fdens, e_io)
      endif
      fdens = MPI_FILE_NULL
      call mpi_barrier(MPI_COMM_WORLD,e_io)

#endif
      return

   end subroutine write_restart_parallel_2c

   subroutine read_restart_parallel_2c(iframe,iframe2D,e_io,hfields_s,phifields_s)

!***********************************************************************
!
!     LBsoft subroutine for opening the vtk legacy file
!     in parallel IO
!
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2019
!
!***********************************************************************

      implicit none

      integer(kind=4), intent(out) :: iframe,iframe2D
      real(kind=strdb), allocatable, dimension(:,:,:,:,:) :: hfields_s,phifields_s
      integer, intent(out) :: e_io
#ifdef MPI
      integer(kind=MPI_OFFSET_KIND) :: tempoffset,offset_final
#endif

      integer :: ioffset
      character(1), parameter :: end_rec = char(10)
      integer, parameter :: bytechar=kind(end_rec)
      integer, parameter :: byteint = 4
      integer, parameter :: byter4  = 4
      integer, parameter :: byter8  = 8
      integer, parameter :: nbuffsub = 0
      integer :: filetypesub,imemtype,filetypesubv,ierr

      logical :: lexist
      integer, dimension(3) :: memDims,memOffs
      integer :: elen,amode

#ifdef MPI
      integer :: fdens = MPI_FILE_NULL
      logical :: file_exists
      integer :: ii,jj,kk,myblock,xblock,yblock,zblock
      character(len=MPI_MAX_ERROR_STRING) :: emsg
      
      call MPI_Barrier(MPI_COMM_WORLD, e_io)
      sevt1 = trim(dir_out) // 'restart.raw'
      
      inquire(file=trim(sevt1), exist=file_exists)
      if (.not. file_exists) then
         call doerror(6,'ERROR: unable to find the file: '//trim(sevt1))
      end if   

      call MPI_FILE_OPEN(MPI_COMM_WORLD, trim(sevt1), &
         MPI_MODE_RDONLY, &
         MPI_INFO_NULL,fdens,e_io)
      
      if (e_io /= MPI_SUCCESS) then
        call MPI_Error_string(e_io, emsg, elen, ierr)
        write(6,'(A,I6,2A)') 'Rank', myrank, ' MPI_File_open error: ', trim(emsg)
        write(6,'(A)') 'Path tried: '//trim(sevt1)
        call MPI_Abort(MPI_COMM_WORLD, e_io, ierr)
      endif     
      
      call MPI_Barrier(MPI_COMM_WORLD, e_io)
      
      call MPI_Type_create_subarray(3,gsizes,lsizes,myoffset, &
         MPI_ORDER_FORTRAN,MPI_REAL,filetypesub,e_io)

      call MPI_Type_commit(filetypesub, e_io)
      
      tempoffset=int(0,kind=MPI_OFFSET_KIND)

      call MPI_File_Set_View(fdens,tempoffset,MPI_REAL,filetypesub, &
         "native",MPI_INFO_NULL,e_io)
      
      call MPI_Barrier(MPI_COMM_WORLD, e_io)
      
      ! We need full local sizes: memDims
      memDims = lsizes + 2*nbuffsub
      memOffs = [ nbuffsub, nbuffsub, nbuffsub ]

      call MPI_TYPE_CREATE_SUBARRAY(3,memDims,lsizes,memOffs, &
         MPI_ORDER_FORTRAN,MPI_REAL,imemtype,e_io)

      call MPI_TYPE_COMMIT(imemtype,e_io)

      call MPI_FILE_READ_ALL(fdens,arr_3d,1,imemtype,MPI_STATUS_IGNORE,e_io)
      do k=1,nz
  	     zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
         do j=1,ny
            yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
            do i=1,nx
               xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
               myblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
               ii=i-xblock*TILE_DIMx+2*TILE_DIMx
               jj=j-yblock*TILE_DIMy+2*TILE_DIMy
               kk=k-zblock*TILE_DIMz+2*TILE_DIMz  
               hfields_s(ii,jj,kk,1,myblock)=real(arr_3d(i,j,k),kind=strdb)  
             enddo
         enddo
      enddo   
      

      call MPI_FILE_READ_ALL(fdens,arr_3d,1,imemtype,MPI_STATUS_IGNORE,e_io)
      do k=1,nz
  	     zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
         do j=1,ny
            yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
            do i=1,nx
               xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
               myblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
               ii=i-xblock*TILE_DIMx+2*TILE_DIMx
               jj=j-yblock*TILE_DIMy+2*TILE_DIMy
               kk=k-zblock*TILE_DIMz+2*TILE_DIMz  
               hfields_s(ii,jj,kk,2,myblock)=real(arr_3d(i,j,k),kind=strdb) 
             enddo
         enddo
      enddo   
      
      
      call MPI_FILE_READ_ALL(fdens,arr_3d,1,imemtype,MPI_STATUS_IGNORE,e_io)
      do k=1,nz
  	     zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
         do j=1,ny
            yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
            do i=1,nx
               xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
               myblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
               ii=i-xblock*TILE_DIMx+2*TILE_DIMx
               jj=j-yblock*TILE_DIMy+2*TILE_DIMy
               kk=k-zblock*TILE_DIMz+2*TILE_DIMz  
               hfields_s(ii,jj,kk,3,myblock)=real(arr_3d(i,j,k),kind=strdb) 
             enddo
         enddo
      enddo   
      
       
      call MPI_FILE_READ_ALL(fdens,arr_3d,1,imemtype,MPI_STATUS_IGNORE,e_io)
      do k=1,nz
  	     zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
         do j=1,ny
            yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
            do i=1,nx
               xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
               myblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
               ii=i-xblock*TILE_DIMx+2*TILE_DIMx
               jj=j-yblock*TILE_DIMy+2*TILE_DIMy
               kk=k-zblock*TILE_DIMz+2*TILE_DIMz  
               hfields_s(ii,jj,kk,4,myblock)=real(arr_3d(i,j,k),kind=strdb) 
             enddo
         enddo
      enddo   
      

      call MPI_FILE_READ_ALL(fdens,arr_3d,1,imemtype,MPI_STATUS_IGNORE,e_io)
      do k=1,nz
  	     zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
         do j=1,ny
            yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
            do i=1,nx
               xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
               myblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
               ii=i-xblock*TILE_DIMx+2*TILE_DIMx
               jj=j-yblock*TILE_DIMy+2*TILE_DIMy
               kk=k-zblock*TILE_DIMz+2*TILE_DIMz  
               hfields_s(ii,jj,kk,5,myblock)=real(arr_3d(i,j,k),kind=strdb) 
             enddo
         enddo
      enddo   
       
      
      call MPI_FILE_READ_ALL(fdens,arr_3d,1,imemtype,MPI_STATUS_IGNORE,e_io)
      do k=1,nz
  	     zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
         do j=1,ny
            yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
            do i=1,nx
               xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
               myblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
               ii=i-xblock*TILE_DIMx+2*TILE_DIMx
               jj=j-yblock*TILE_DIMy+2*TILE_DIMy
               kk=k-zblock*TILE_DIMz+2*TILE_DIMz  
               hfields_s(ii,jj,kk,8,myblock)=real(arr_3d(i,j,k),kind=strdb) 
             enddo
         enddo
      enddo   
      
      
      call MPI_FILE_READ_ALL(fdens,arr_3d,1,imemtype,MPI_STATUS_IGNORE,e_io)
      do k=1,nz
  	     zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
         do j=1,ny
            yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
            do i=1,nx
               xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
               myblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
               ii=i-xblock*TILE_DIMx+2*TILE_DIMx
               jj=j-yblock*TILE_DIMy+2*TILE_DIMy
               kk=k-zblock*TILE_DIMz+2*TILE_DIMz  
               hfields_s(ii,jj,kk,9,myblock)=real(arr_3d(i,j,k),kind=strdb) 
             enddo
         enddo
      enddo   
      
      
      call MPI_FILE_READ_ALL(fdens,arr_3d,1,imemtype,MPI_STATUS_IGNORE,e_io)
      do k=1,nz
  	     zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
         do j=1,ny
            yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
            do i=1,nx
               xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
               myblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
               ii=i-xblock*TILE_DIMx+2*TILE_DIMx
               jj=j-yblock*TILE_DIMy+2*TILE_DIMy
               kk=k-zblock*TILE_DIMz+2*TILE_DIMz  
               hfields_s(ii,jj,kk,6,myblock)=real(arr_3d(i,j,k),kind=strdb) 
             enddo
         enddo
      enddo   
      
      
      call MPI_FILE_READ_ALL(fdens,arr_3d,1,imemtype,MPI_STATUS_IGNORE,e_io)
      do k=1,nz
  	     zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
         do j=1,ny
            yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
            do i=1,nx
               xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
               myblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
               ii=i-xblock*TILE_DIMx+2*TILE_DIMx
               jj=j-yblock*TILE_DIMy+2*TILE_DIMy
               kk=k-zblock*TILE_DIMz+2*TILE_DIMz  
               hfields_s(ii,jj,kk,10,myblock)=real(arr_3d(i,j,k),kind=strdb) 
             enddo
         enddo
      enddo   
      
      
      call MPI_FILE_READ_ALL(fdens,arr_3d,1,imemtype,MPI_STATUS_IGNORE,e_io)
      do k=1,nz
  	     zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
         do j=1,ny
            yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
            do i=1,nx
               xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
               myblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
               ii=i-xblock*TILE_DIMx+2*TILE_DIMx
               jj=j-yblock*TILE_DIMy+2*TILE_DIMy
               kk=k-zblock*TILE_DIMz+2*TILE_DIMz  
               hfields_s(ii,jj,kk,7,myblock)=real(arr_3d(i,j,k),kind=strdb) 
             enddo
         enddo
      enddo   
      

      call MPI_FILE_READ_ALL(fdens,arr_3d,1,imemtype,MPI_STATUS_IGNORE,e_io)
      do k=1,nz
  	     zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
         do j=1,ny
            yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
            do i=1,nx
               xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
               myblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
               ii=i-xblock*TILE_DIMx+2*TILE_DIMx
               jj=j-yblock*TILE_DIMy+2*TILE_DIMy
               kk=k-zblock*TILE_DIMz+2*TILE_DIMz   
               phifields_s(ii,jj,kk,1,myblock)=real(arr_3d(i,j,k),kind=strdb) 
             enddo
         enddo
      enddo    
      
      call MPI_Type_free(imemtype, e_io)
      call MPI_Type_free(filetypesub, e_io)
      
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
      call MPI_FILE_CLOSE(fdens, e_io)
      if (e_io /= MPI_SUCCESS) then
        call MPI_Error_string(e_io, emsg, elen, ierr)
        write(6,'(A,I6,2A)') 'Rank', myrank, ' MPI_File_close error: ', trim(emsg)
        write(6,'(A)') 'Path tried: '//trim(sevt1)
        call MPI_Abort(MPI_COMM_WORLD, e_io, ierr)
      endif     
      fdens = MPI_FILE_NULL
      call mpi_barrier(MPI_COMM_WORLD,e_io)
      if(myrank == 0) then
        call MPI_FILE_OPEN(MPI_COMM_SELF, trim(sevt1), MPI_MODE_RDONLY, MPI_INFO_NULL, fdens, e_io)
        if (e_io /= MPI_SUCCESS) then
          call MPI_Error_string(e_io, emsg, elen, ierr)
          write(6,'(A,I6,2A)') 'Rank', myrank, ' MPI_File_open error: ', trim(emsg)
          write(6,'(A)') 'Path tried: '//trim(sevt1)
          call MPI_Abort(MPI_COMM_WORLD, e_io, ierr)
        endif  
        call MPI_File_set_view(fdens, 0_MPI_OFFSET_KIND, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL, e_io)

        offset_final = int(lx * ly * lz * 11 * db, kind=MPI_OFFSET_KIND)
        call MPI_FILE_READ_AT(fdens, offset_final, iframe, 1, MPI_INTEGER, MPI_STATUS_IGNORE, e_io)
        call MPI_FILE_READ_AT(fdens, offset_final + byteint, iframe2D, 1, MPI_INTEGER, MPI_STATUS_IGNORE, e_io)

        call MPI_FILE_CLOSE(fdens, e_io)
        if (e_io /= MPI_SUCCESS) then
          call MPI_Error_string(e_io, emsg, elen, ierr)
          write(6,'(A,I6,2A)') 'Rank', myrank, ' MPI_File_close error: ', trim(emsg)
          write(6,'(A)') 'Path tried: '//trim(sevt1)
          call MPI_Abort(MPI_COMM_WORLD, e_io, ierr)
        endif 
      endif
      fdens = MPI_FILE_NULL
      call mpi_barrier(MPI_COMM_WORLD,e_io)
#endif
      return

   end subroutine read_restart_parallel_2c
   
   subroutine read_isfluid_parallel(iframe,e_io)

!***********************************************************************
!
!     LBsoft subroutine for opening the vtk legacy file
!     in parallel IO
!
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2019
!
!***********************************************************************

      implicit none

      integer(kind=4), intent(in) ::  iframe
      integer, intent(out) :: e_io
#ifdef MPI
      integer(kind=MPI_OFFSET_KIND) :: tempoffset,offset_final
#endif

      integer :: ioffset
      character(1), parameter :: end_rec = char(10)
      integer, parameter :: bytechar=kind(end_rec)
      integer, parameter :: byteint = 4
      integer, parameter :: byter4  = 4
      integer, parameter :: byter8  = 8
      integer, parameter :: nbuffsub = 0
      integer :: filetypesub,imemtype,filetypesubv,ierr
      logical :: lexist
      integer, dimension(3) :: memDims,memOffs
      integer :: elen,amode

#ifdef MPI
      integer :: fdens=MPI_FILE_NULL
      character(len=MPI_MAX_ERROR_STRING) :: emsg
      logical :: file_exists

      sevt1 = 'isfluid.raw'
      
      inquire(file=trim(sevt1), exist=file_exists)
      if (.not. file_exists) then
         call doerror(6,'ERROR: unable to find the file: '//trim(sevt1))
      end if    

      call MPI_FILE_OPEN(MPI_COMM_WORLD, trim(sevt1), &
         MPI_MODE_RDONLY, &
         MPI_INFO_NULL,fdens,e_io)
      if (e_io /= MPI_SUCCESS) then
        call MPI_Error_string(e_io, emsg, elen, ierr)
        write(6,'(A,I6,2A)') 'Rank', myrank, ' MPI_File_open error: ', trim(emsg)
        write(6,'(A)') 'Path tried: '//trim(sevt1)
        call MPI_Abort(MPI_COMM_WORLD, e_io, ierr)
      endif   
      call MPI_Barrier(MPI_COMM_WORLD, e_io)

      call MPI_Type_create_subarray(3,gsizes,lsizes,myoffset, &
         MPI_ORDER_FORTRAN,MPI_REAL,filetypesub,e_io)
      
      call MPI_Type_commit(filetypesub, e_io)
      
      tempoffset=int(0,kind=MPI_OFFSET_KIND)

      call MPI_File_Set_View(fdens,tempoffset,MPI_REAL,filetypesub, &
         "native",MPI_INFO_NULL,e_io)
      
      call MPI_Barrier(MPI_COMM_WORLD, e_io)
      
      ! We need full local sizes: memDims
      memDims = lsizes + 2*nbuffsub
      memOffs = [ nbuffsub, nbuffsub, nbuffsub ]

      call MPI_TYPE_CREATE_SUBARRAY(3,memDims,lsizes,memOffs, &
         MPI_ORDER_FORTRAN,MPI_REAL,imemtype,e_io)

      call MPI_TYPE_COMMIT(imemtype,e_io)

      call MPI_FILE_READ_ALL(fdens,arr_3d,1,imemtype,MPI_STATUS_IGNORE,e_io)
      isfluid(1:nx,1:ny,1:nz)= int(arr_3d(1:nx,1:ny,1:nz),kind=isf) 
      
      call MPI_Type_free(imemtype, e_io)
      call MPI_Type_free(filetypesub, e_io)
      
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
      call MPI_FILE_CLOSE(fdens, e_io)
      if (e_io /= MPI_SUCCESS) then
        call MPI_Error_string(e_io, emsg, elen, ierr)
        write(6,'(A,I6,2A)') 'Rank', myrank, ' MPI_File_close error: ', trim(emsg)
        write(6,'(A)') 'Path tried: '//trim(sevt1)
        call MPI_Abort(MPI_COMM_WORLD, e_io, ierr)
      endif     
      fdens = MPI_FILE_NULL
      call mpi_barrier(MPI_COMM_WORLD,e_io)
      
#endif
      return

   end subroutine read_isfluid_parallel
   
   subroutine read_init_parallel(iframe,e_io,hfields_s,phifields_s)

!***********************************************************************
!
!     LBsoft subroutine for opening the vtk legacy file
!     in parallel IO
!
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2019
!
!***********************************************************************

      implicit none

      integer(kind=4), intent(in) ::  iframe
      integer, intent(out) :: e_io
      real(kind=strdb), allocatable, dimension(:,:,:,:,:) :: hfields_s,phifields_s
#ifdef MPI
      integer(kind=MPI_OFFSET_KIND) :: tempoffset,offset_final
#endif

      integer :: ioffset
      character(1), parameter :: end_rec = char(10)
      integer, parameter :: bytechar=kind(end_rec)
      integer, parameter :: byteint = 4
      integer, parameter :: byter4  = 4
      integer, parameter :: byter8  = 8
      integer, parameter :: nbuffsub = 0
      integer :: filetypesub,imemtype,filetypesubv,ierr
      logical :: lexist
      integer, dimension(3) :: memDims,memOffs
      integer :: elen,amode

#ifdef MPI
      integer :: fdens=MPI_FILE_NULL
      character(len=MPI_MAX_ERROR_STRING) :: emsg
      logical :: file_exists
      integer :: ii,jj,kk,myblock,xblock,yblock,zblock
!!!!!!!!!!!!!!!!!!!!!!!!!!rho!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      sevt1 = 'rho.raw'
      
      inquire(file=trim(sevt1), exist=file_exists)
      if (.not. file_exists) then
         call doerror(6,'ERROR: unable to find the file: '//trim(sevt1))
      end if    

      call MPI_FILE_OPEN(MPI_COMM_WORLD, trim(sevt1), &
         MPI_MODE_RDONLY, &
         MPI_INFO_NULL,fdens,e_io)
      if (e_io /= MPI_SUCCESS) then
        call MPI_Error_string(e_io, emsg, elen, ierr)
        write(6,'(A,I6,2A)') 'Rank', myrank, ' MPI_File_open error: ', trim(emsg)
        write(6,'(A)') 'Path tried: '//trim(sevt1)
        call MPI_Abort(MPI_COMM_WORLD, e_io, ierr)
      endif   
      call MPI_Barrier(MPI_COMM_WORLD, e_io)

      call MPI_Type_create_subarray(3,gsizes,lsizes,myoffset, &
         MPI_ORDER_FORTRAN,MPI_REAL,filetypesub,e_io)

      call MPI_Type_commit(filetypesub, e_io)
      
      tempoffset=int(0,kind=MPI_OFFSET_KIND)

      call MPI_File_Set_View(fdens,tempoffset,MPI_REAL,filetypesub, &
         "native",MPI_INFO_NULL,e_io)
      
      call MPI_Barrier(MPI_COMM_WORLD, e_io)
      
      ! We need full local sizes: memDims
      memDims = lsizes + 2*nbuffsub
      memOffs = [ nbuffsub, nbuffsub, nbuffsub ]

      call MPI_TYPE_CREATE_SUBARRAY(3,memDims,lsizes,memOffs, &
         MPI_ORDER_FORTRAN,MPI_REAL,imemtype,e_io)

      call MPI_TYPE_COMMIT(imemtype,e_io)

      call MPI_FILE_READ_ALL(fdens,arr_3d,1,imemtype,MPI_STATUS_IGNORE,e_io)
      
      call MPI_Type_free(imemtype, e_io)
      call MPI_Type_free(filetypesub, e_io)
      
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
      call MPI_FILE_CLOSE(fdens, e_io)
      if (e_io /= MPI_SUCCESS) then
        call MPI_Error_string(e_io, emsg, elen, ierr)
        write(6,'(A,I6,2A)') 'Rank', myrank, ' MPI_File_close error: ', trim(emsg)
        write(6,'(A)') 'Path tried: '//trim(sevt1)
        call MPI_Abort(MPI_COMM_WORLD, e_io, ierr)
      endif     
      fdens = MPI_FILE_NULL
      call mpi_barrier(MPI_COMM_WORLD,e_io)
      
      do k=1,nz
  	     zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
         do j=1,ny
            yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
            do i=1,nx
               xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
               myblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
               ii=i-xblock*TILE_DIMx+2*TILE_DIMx
               jj=j-yblock*TILE_DIMy+2*TILE_DIMy
               kk=k-zblock*TILE_DIMz+2*TILE_DIMz  
               hfields_s(ii,jj,kk,1,myblock)=real(arr_3d(i,j,k),kind=strdb)   
             enddo
         enddo
      enddo 
      
!!!!!!!!!!!!!!!!!!!!!!!!!!u!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      sevt1 = 'u.raw'
      
      inquire(file=trim(sevt1), exist=file_exists)
      if (.not. file_exists) then
         call doerror(6,'ERROR: unable to find the file: '//trim(sevt1))
      end if    

      call MPI_FILE_OPEN(MPI_COMM_WORLD, trim(sevt1), &
         MPI_MODE_RDONLY, &
         MPI_INFO_NULL,fdens,e_io)
      if (e_io /= MPI_SUCCESS) then
        call MPI_Error_string(e_io, emsg, elen, ierr)
        write(6,'(A,I6,2A)') 'Rank', myrank, ' MPI_File_open error: ', trim(emsg)
        write(6,'(A)') 'Path tried: '//trim(sevt1)
        call MPI_Abort(MPI_COMM_WORLD, e_io, ierr)
      endif   
      call MPI_Barrier(MPI_COMM_WORLD, e_io)

      call MPI_Type_create_subarray(3,gsizes,lsizes,myoffset, &
         MPI_ORDER_FORTRAN,MPI_REAL,filetypesub,e_io)

      call MPI_Type_commit(filetypesub, e_io)
      
      tempoffset=int(0,kind=MPI_OFFSET_KIND)

      call MPI_File_Set_View(fdens,tempoffset,MPI_REAL,filetypesub, &
         "native",MPI_INFO_NULL,e_io)
      
      call MPI_Barrier(MPI_COMM_WORLD, e_io)
      
      ! We need full local sizes: memDims
      memDims = lsizes + 2*nbuffsub
      memOffs = [ nbuffsub, nbuffsub, nbuffsub ]

      call MPI_TYPE_CREATE_SUBARRAY(3,memDims,lsizes,memOffs, &
         MPI_ORDER_FORTRAN,MPI_REAL,imemtype,e_io)

      call MPI_TYPE_COMMIT(imemtype,e_io)

      call MPI_FILE_READ_ALL(fdens,arr_3d,1,imemtype,MPI_STATUS_IGNORE,e_io)
      
      call MPI_Type_free(imemtype, e_io)
      call MPI_Type_free(filetypesub, e_io)
      
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
      call MPI_FILE_CLOSE(fdens, e_io)
      if (e_io /= MPI_SUCCESS) then
        call MPI_Error_string(e_io, emsg, elen, ierr)
        write(6,'(A,I6,2A)') 'Rank', myrank, ' MPI_File_close error: ', trim(emsg)
        write(6,'(A)') 'Path tried: '//trim(sevt1)
        call MPI_Abort(MPI_COMM_WORLD, e_io, ierr)
      endif     
      fdens = MPI_FILE_NULL
      call mpi_barrier(MPI_COMM_WORLD,e_io)
          
      do k=1,nz
  	     zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
         do j=1,ny
            yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
            do i=1,nx
               xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
               myblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
               ii=i-xblock*TILE_DIMx+2*TILE_DIMx
               jj=j-yblock*TILE_DIMy+2*TILE_DIMy
               kk=k-zblock*TILE_DIMz+2*TILE_DIMz  
               hfields_s(ii,jj,kk,2,myblock)=real(arr_3d(i,j,k),kind=strdb) 
             enddo
         enddo
      enddo 
      
!!!!!!!!!!!!!!!!!!!!!!!!!!v!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      sevt1 = 'v.raw'
      
      inquire(file=trim(sevt1), exist=file_exists)
      if (.not. file_exists) then
         call doerror(6,'ERROR: unable to find the file: '//trim(sevt1))
      end if    

      call MPI_FILE_OPEN(MPI_COMM_WORLD, trim(sevt1), &
         MPI_MODE_RDONLY, &
         MPI_INFO_NULL,fdens,e_io)
      if (e_io /= MPI_SUCCESS) then
        call MPI_Error_string(e_io, emsg, elen, ierr)
        write(6,'(A,I6,2A)') 'Rank', myrank, ' MPI_File_open error: ', trim(emsg)
        write(6,'(A)') 'Path tried: '//trim(sevt1)
        call MPI_Abort(MPI_COMM_WORLD, e_io, ierr)
      endif   
      call MPI_Barrier(MPI_COMM_WORLD, e_io)

      call MPI_Type_create_subarray(3,gsizes,lsizes,myoffset, &
         MPI_ORDER_FORTRAN,MPI_REAL,filetypesub,e_io)

      call MPI_Type_commit(filetypesub, e_io)
      
      tempoffset=int(0,kind=MPI_OFFSET_KIND)

      call MPI_File_Set_View(fdens,tempoffset,MPI_REAL,filetypesub, &
         "native",MPI_INFO_NULL,e_io)
      
      call MPI_Barrier(MPI_COMM_WORLD, e_io)
      
      ! We need full local sizes: memDims
      memDims = lsizes + 2*nbuffsub
      memOffs = [ nbuffsub, nbuffsub, nbuffsub ]

      call MPI_TYPE_CREATE_SUBARRAY(3,memDims,lsizes,memOffs, &
         MPI_ORDER_FORTRAN,MPI_REAL,imemtype,e_io)

      call MPI_TYPE_COMMIT(imemtype,e_io)

      call MPI_FILE_READ_ALL(fdens,arr_3d,1,imemtype,MPI_STATUS_IGNORE,e_io)
      
      call MPI_Type_free(imemtype, e_io)
      call MPI_Type_free(filetypesub, e_io)
      
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
      call MPI_FILE_CLOSE(fdens, e_io)
      if (e_io /= MPI_SUCCESS) then
        call MPI_Error_string(e_io, emsg, elen, ierr)
        write(6,'(A,I6,2A)') 'Rank', myrank, ' MPI_File_close error: ', trim(emsg)
        write(6,'(A)') 'Path tried: '//trim(sevt1)
        call MPI_Abort(MPI_COMM_WORLD, e_io, ierr)
      endif     
      fdens = MPI_FILE_NULL
      call mpi_barrier(MPI_COMM_WORLD,e_io)
      
      do k=1,nz
  	     zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
         do j=1,ny
            yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
            do i=1,nx
               xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
               myblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
               ii=i-xblock*TILE_DIMx+2*TILE_DIMx
               jj=j-yblock*TILE_DIMy+2*TILE_DIMy
               kk=k-zblock*TILE_DIMz+2*TILE_DIMz  
               hfields_s(ii,jj,kk,3,myblock)=real(arr_3d(i,j,k),kind=strdb) 
             enddo
         enddo
      enddo 
      
!!!!!!!!!!!!!!!!!!!!!!!!!!w!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      sevt1 = 'w.raw'
      
      inquire(file=trim(sevt1), exist=file_exists)
      if (.not. file_exists) then
         call doerror(6,'ERROR: unable to find the file: '//trim(sevt1))
      end if    

      call MPI_FILE_OPEN(MPI_COMM_WORLD, trim(sevt1), &
         MPI_MODE_RDONLY, &
         MPI_INFO_NULL,fdens,e_io)
      if (e_io /= MPI_SUCCESS) then
        call MPI_Error_string(e_io, emsg, elen, ierr)
        write(6,'(A,I6,2A)') 'Rank', myrank, ' MPI_File_open error: ', trim(emsg)
        write(6,'(A)') 'Path tried: '//trim(sevt1)
        call MPI_Abort(MPI_COMM_WORLD, e_io, ierr)
      endif   
      call MPI_Barrier(MPI_COMM_WORLD, e_io)

      call MPI_Type_create_subarray(3,gsizes,lsizes,myoffset, &
         MPI_ORDER_FORTRAN,MPI_REAL,filetypesub,e_io)

      call MPI_Type_commit(filetypesub, e_io)
      
      tempoffset=int(0,kind=MPI_OFFSET_KIND)

      call MPI_File_Set_View(fdens,tempoffset,MPI_REAL,filetypesub, &
         "native",MPI_INFO_NULL,e_io)
      
      call MPI_Barrier(MPI_COMM_WORLD, e_io)
      
      ! We need full local sizes: memDims
      memDims = lsizes + 2*nbuffsub
      memOffs = [ nbuffsub, nbuffsub, nbuffsub ]

      call MPI_TYPE_CREATE_SUBARRAY(3,memDims,lsizes,memOffs, &
         MPI_ORDER_FORTRAN,MPI_REAL,imemtype,e_io)

      call MPI_TYPE_COMMIT(imemtype,e_io)

      call MPI_FILE_READ_ALL(fdens,arr_3d,1,imemtype,MPI_STATUS_IGNORE,e_io)
      
      call MPI_Type_free(imemtype, e_io)
      call MPI_Type_free(filetypesub, e_io)
      
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
      call MPI_FILE_CLOSE(fdens, e_io)
      if (e_io /= MPI_SUCCESS) then
        call MPI_Error_string(e_io, emsg, elen, ierr)
        write(6,'(A,I6,2A)') 'Rank', myrank, ' MPI_File_close error: ', trim(emsg)
        write(6,'(A)') 'Path tried: '//trim(sevt1)
        call MPI_Abort(MPI_COMM_WORLD, e_io, ierr)
      endif     
      fdens = MPI_FILE_NULL
      call mpi_barrier(MPI_COMM_WORLD,e_io)
      
      do k=1,nz
  	     zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
         do j=1,ny
            yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
            do i=1,nx
               xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
               myblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
               ii=i-xblock*TILE_DIMx+2*TILE_DIMx
               jj=j-yblock*TILE_DIMy+2*TILE_DIMy
               kk=k-zblock*TILE_DIMz+2*TILE_DIMz  
               hfields_s(ii,jj,kk,4,myblock)=real(arr_3d(i,j,k),kind=strdb) 
             enddo
         enddo
      enddo 

#ifdef TWOCOMPONENT      
!!!!!!!!!!!!!!!!!!!!!!!!!!phi!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      sevt1 = 'phi.raw'
      
      inquire(file=trim(sevt1), exist=file_exists)
      if (.not. file_exists) then
         call doerror(6,'ERROR: unable to find the file: '//trim(sevt1))
      end if    

      call MPI_FILE_OPEN(MPI_COMM_WORLD, trim(sevt1), &
         MPI_MODE_RDONLY, &
         MPI_INFO_NULL,fdens,e_io)
      if (e_io /= MPI_SUCCESS) then
        call MPI_Error_string(e_io, emsg, elen, ierr)
        write(6,'(A,I6,2A)') 'Rank', myrank, ' MPI_File_open error: ', trim(emsg)
        write(6,'(A)') 'Path tried: '//trim(sevt1)
        call MPI_Abort(MPI_COMM_WORLD, e_io, ierr)
      endif   
      call MPI_Barrier(MPI_COMM_WORLD, e_io)
      
      call MPI_Type_create_subarray(3,gsizes,lsizes,myoffset, &
         MPI_ORDER_FORTRAN,MPI_REAL,filetypesub,e_io)

      call MPI_Type_commit(filetypesub, e_io)
      
      tempoffset=int(0,kind=MPI_OFFSET_KIND)

      call MPI_File_Set_View(fdens,tempoffset,MPI_REAL,filetypesub, &
         "native",MPI_INFO_NULL,e_io)
      
      call MPI_Barrier(MPI_COMM_WORLD, e_io)
      
      ! We need full local sizes: memDims
      memDims = lsizes + 2*nbuffsub
      memOffs = [ nbuffsub, nbuffsub, nbuffsub ]

      call MPI_TYPE_CREATE_SUBARRAY(3,memDims,lsizes,memOffs, &
         MPI_ORDER_FORTRAN,MPI_REAL,imemtype,e_io)

      call MPI_TYPE_COMMIT(imemtype,e_io)

      call MPI_FILE_READ_ALL(fdens,arr_3d,1,imemtype,MPI_STATUS_IGNORE,e_io)
      
      call MPI_Type_free(imemtype, e_io)
      call MPI_Type_free(filetypesub, e_io)
      
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
      call MPI_FILE_CLOSE(fdens, e_io)
      if (e_io /= MPI_SUCCESS) then
        call MPI_Error_string(e_io, emsg, elen, ierr)
        write(6,'(A,I6,2A)') 'Rank', myrank, ' MPI_File_close error: ', trim(emsg)
        write(6,'(A)') 'Path tried: '//trim(sevt1)
        call MPI_Abort(MPI_COMM_WORLD, e_io, ierr)
      endif     
      fdens = MPI_FILE_NULL
      call mpi_barrier(MPI_COMM_WORLD,e_io)

      do k=1,nz
  	     zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
         do j=1,ny
            yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
            do i=1,nx
               xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
               myblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
               ii=i-xblock*TILE_DIMx+2*TILE_DIMx
               jj=j-yblock*TILE_DIMy+2*TILE_DIMy
               kk=k-zblock*TILE_DIMz+2*TILE_DIMz   
               phifields_s(ii,jj,kk,1,myblock)=real(arr_3d(i,j,k),kind=strdb) 
             enddo
         enddo
      enddo    
            
#endif
#endif
      return

   end subroutine read_init_parallel
   
   subroutine write_file_raw_par_isfluid(iframe,e_io)

!***********************************************************************
!
!     LBsoft subroutine for opening the vtk legacy file
!     in parallel IO
!
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2019
!
!***********************************************************************

      implicit none

      integer, intent(in) ::iframe
      integer, intent(out) :: e_io
      character(len=8) :: namevarvtk_sub='isfluid '
#ifdef MPI
      integer(kind=MPI_OFFSET_KIND) :: tempoffset
#endif

      integer :: ioffset
      character(1), parameter :: end_rec = char(10)
      integer, parameter :: bytechar=kind(end_rec)
      integer, parameter :: byteint = 4
      integer, parameter :: byter4  = 4
      integer, parameter :: byter8  = 8
      integer, parameter :: nbuffsub = 0
      integer :: filetypesub,imemtype,filetypesubv,ierr
      logical :: lexist,lexit

      integer :: elen,amode
      
      integer, dimension(3) :: memDims,memOffs
      integer, dimension(4) :: velglobalDims,velldims,velmystarts, &
         velmemDims,velmemOffs
#ifdef MPI
      integer :: fdens= MPI_FILE_NULL
      integer :: fvel= MPI_FILE_NULL
      character(len=MPI_MAX_ERROR_STRING) :: emsg
#ifdef DOXDMF      
      integer, parameter :: xml_file=734
      character(len=32) :: strx,stry,strz
#endif

      call MPI_Barrier(MPI_COMM_WORLD, e_io)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!density!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 

      sevt1 = trim(dir_out) // trim(filenamevtk)//'_'//trim(namevarvtk_sub)// &
         '_'//trim(write_fmtnumb(iframe)) // '.raw'

#ifdef DOXDMF
      if(myrank==0)then
        write(strx, '(I0)') skip_gsizes(3)
        write(stry, '(I0)') skip_gsizes(2)
        write(strz, '(I0)') skip_gsizes(1)
        
        sevt2=repeat(' ',mxln)
        sevt2 =  trim(dir_out) // trim(filenamevtk)//'_'//trim(namevarvtk_sub)// &
           '_'//trim(write_fmtnumb(iframe)) // '.xdmf'
        open(unit=xml_file, file=trim(sevt2), status="replace", action="write")
        
        sevt2=repeat(' ',mxln)
        sevt2 =  trim(filenamevtk)//'_'//trim(namevarvtk_sub)// &
           '_'//trim(write_fmtnumb(iframe)) // '.raw'
        
        write(xml_file,'(a)') '<?xml version="1.0" ?>'
        write(xml_file,'(a)') '<Xdmf Version="3.0">'
        write(xml_file,'(a)') '  <Domain>'
        write(xml_file,'(a)') '    <!-- Uniform Grid -->'
        write(xml_file,'(a)') '    <Grid Name="ScalarField" GridType="Uniform">'
        write(xml_file,'(a)') ''
        write(xml_file,'(a)') '      <!-- Topology: 3D grid -->'
        write(xml_file,'(7a)') '      <Topology TopologyType="3DCoRectMesh" Dimensions="' // &
         trim(adjustl(strx)) // ' ' // trim(adjustl(stry)) // ' ' // trim(adjustl(strz)) // '"/>'
        write(xml_file,'(a)') ''
        write(xml_file,'(a)') '      <!-- Geometry: origin and spacing (dx, dy, dz) -->'
        write(xml_file,'(a)') '      <Geometry GeometryType="ORIGIN_DXDYDZ">'
        write(xml_file,'(a)') '        <DataItem Dimensions="3" Format="XML" Precision="4">0.0 0.0 0.0</DataItem>'
        write(xml_file,'(a)') '        <DataItem Dimensions="3" Format="XML" Precision="4">1.0 1.0 1.0</DataItem>'
        write(xml_file,'(a)') '      </Geometry>'
        write(xml_file,'(a)') ''
        write(xml_file,'(a)') '      <!-- Scalar field associated with the grid -->'
        write(xml_file,'(3a)') '      <Attribute Name="'//trim(namevarvtk_sub)//'" AttributeType="Scalar" Center="Node">'
        write(xml_file,'(7a)') '        <DataItem Format="Binary" NumberType="Float" Precision="4" Dimensions="' // &
         trim(adjustl(strx)) // ' ' // trim(adjustl(stry)) // ' ' // trim(adjustl(strz)) // '">'
        write(xml_file,'(2a)') '          ',trim(sevt2)
        write(xml_file,'(a)') '        </DataItem>'
        write(xml_file,'(a)') '      </Attribute>'
        write(xml_file,'(a)') ''
        write(xml_file,'(a)') '    </Grid>'
        write(xml_file,'(a)') '  </Domain>'
        write(xml_file,'(a)') '</Xdmf>'

        close(xml_file)
        
      endif
#endif
      
      do k=1,nzskip
        do j=1,nyskip
          do i=1,nxskip
            rhoprint(i,j,k)=real(isfluid(i*stepskip,j*stepskip,k*stepskip),kind=printdb)
          enddo
        enddo
      enddo
      
      lexit=.false.
      if (myrank == 0) then
        inquire(file=trim(sevt1), exist=lexist)
        if (lexist) then
          lexit=.true.
        endif
      endif
      call MPI_Barrier(MPI_COMM_WORLD, e_io)
      call or_world_l(lexit)
      if(lexit)return
      
      amode = IOR(MPI_MODE_CREATE, MPI_MODE_WRONLY) 

      call MPI_File_open(MPI_COMM_WORLD, trim(sevt1),amode, MPI_INFO_NULL, fdens, e_io)

      if (e_io /= MPI_SUCCESS) then
        call MPI_Error_string(e_io, emsg, elen, ierr)
        write(6,'(A,I6,2A)') 'Rank', myrank, ' MPI_File_open error: ', trim(emsg)
        write(6,'(A)') 'Path tried: '//trim(sevt1)
        call MPI_Abort(MPI_COMM_WORLD, e_io, ierr)
      endif     
      
        call MPI_File_set_size(fdens, 0_MPI_OFFSET_KIND, e_io)  
      
      call MPI_Barrier(MPI_COMM_WORLD, e_io)     ! o il communicator usato per aprire il file
      tempoffset=int(0,kind=MPI_OFFSET_KIND)

#ifdef PRINTHALF
      call MPI_Type_create_subarray(3,skip_gsizes,skip_lsizes,skip_myoffset, &
         MPI_ORDER_FORTRAN,MPI_UNSIGNED_SHORT,filetypesub,e_io)
#else
      call MPI_Type_create_subarray(3,skip_gsizes,skip_lsizes,skip_myoffset, &
         MPI_ORDER_FORTRAN,MPI_REAL,filetypesub,e_io)
#endif
      call MPI_Type_commit(filetypesub, e_io)

#ifdef PRINTHALF
      call MPI_File_Set_View(fdens,tempoffset,MPI_UNSIGNED_SHORT,filetypesub, &
         "native",MPI_INFO_NULL,e_io)
#else
      call MPI_File_Set_View(fdens,tempoffset,MPI_REAL,filetypesub, &
         "native",MPI_INFO_NULL,e_io)
#endif
      ! We need full local sizes: memDims
      memDims = skip_lsizes + 2*nbuffsub
      memOffs = [ nbuffsub, nbuffsub, nbuffsub ]

#ifdef PRINTHALF
      call MPI_TYPE_CREATE_SUBARRAY(3,memDims,skip_lsizes,memOffs, &
         MPI_ORDER_FORTRAN,MPI_UNSIGNED_SHORT,imemtype,e_io)
#else
      call MPI_TYPE_CREATE_SUBARRAY(3,memDims,skip_lsizes,memOffs, &
         MPI_ORDER_FORTRAN,MPI_REAL,imemtype,e_io)
#endif
      call MPI_TYPE_COMMIT(imemtype,e_io)

      call MPI_FILE_WRITE_ALL(fdens,rhoprint,1,imemtype,MPI_STATUS_IGNORE,e_io)
      call MPI_Type_free(imemtype,   e_io)
      call MPI_Type_free(filetypesub,e_io)
      call MPI_File_sync(fdens, e_io)

      call MPI_Barrier(MPI_COMM_WORLD, ierr)
      call MPI_FILE_CLOSE(fdens,e_io)
      if (e_io /= MPI_SUCCESS) then
        call MPI_Error_string(e_io, emsg, elen, ierr)
        write(6,'(A,I6,2A)') 'Rank', myrank, ' MPI_File_close error: ', trim(emsg)
        write(6,'(A)') 'Path tried: '//trim(sevt1)
        call MPI_Abort(MPI_COMM_WORLD, e_io, ierr)
      endif     
      fdens = MPI_FILE_NULL
#endif

      return

   end subroutine write_file_raw_par_isfluid

   subroutine write_file_raw_par(iframe,e_io)

!***********************************************************************
!
!     LBsoft subroutine for opening the vtk legacy file
!     in parallel IO
!
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2019
!
!***********************************************************************

      implicit none

      integer, intent(in) ::iframe
      integer, intent(out) :: e_io
#ifdef MPI
      integer(kind=MPI_OFFSET_KIND) :: tempoffset
#endif

      integer :: ioffset
      character(1), parameter :: end_rec = char(10)
      integer, parameter :: bytechar=kind(end_rec)
      integer, parameter :: byteint = 4
      integer, parameter :: byter4  = 4
      integer, parameter :: byter8  = 8
      integer, parameter :: nbuffsub = 0
      integer :: filetypesub,imemtype,filetypesubv,ierr

      integer, dimension(3) :: memDims,memOffs
      integer, dimension(4) :: velglobalDims,velldims,velmystarts, &
         velmemDims,velmemOffs
      integer :: elen,amode
      logical :: lexist
#ifdef MPI
      integer :: fdens=MPI_FILE_NULL
      integer :: fvel=MPI_FILE_NULL
      character(len=MPI_MAX_ERROR_STRING) :: emsg
#ifdef DOXDMF      
      integer, parameter :: xml_file=734
      character(len=32) :: strx,stry,strz
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!density!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      sevt1 = trim(dir_out) // trim(filenamevtk)//'_'//trim(namevarvtk(1))// &
         '_'//trim(write_fmtnumb(iframe)) // '.raw'

#ifdef DOXDMF
      if(myrank==0)then
        write(strx, '(I0)') skip_gsizes(3)
        write(stry, '(I0)') skip_gsizes(2)
        write(strz, '(I0)') skip_gsizes(1)
        
        sevt2=repeat(' ',mxln)
        sevt2 =  trim(dir_out) // trim(filenamevtk)//'_'//trim(namevarvtk(1))// &
           '_'//trim(write_fmtnumb(iframe)) // '.xdmf'
        open(unit=xml_file, file=trim(sevt2), status="replace", action="write")
        
        sevt2=repeat(' ',mxln)
        sevt2 =  trim(filenamevtk)//'_'//trim(namevarvtk(1))// &
           '_'//trim(write_fmtnumb(iframe)) // '.raw'
        
        write(xml_file,'(a)') '<?xml version="1.0" ?>'
        write(xml_file,'(a)') '<Xdmf Version="3.0">'
        write(xml_file,'(a)') '  <Domain>'
        write(xml_file,'(a)') '    <!-- Uniform Grid -->'
        write(xml_file,'(a)') '    <Grid Name="ScalarField" GridType="Uniform">'
        write(xml_file,'(a)') ''
        write(xml_file,'(a)') '      <!-- Topology: 3D grid -->'
        write(xml_file,'(7a)') '      <Topology TopologyType="3DCoRectMesh" Dimensions="' // &
         trim(adjustl(strx)) // ' ' // trim(adjustl(stry)) // ' ' // trim(adjustl(strz)) // '"/>'
        write(xml_file,'(a)') ''
        write(xml_file,'(a)') '      <!-- Geometry: origin and spacing (dx, dy, dz) -->'
        write(xml_file,'(a)') '      <Geometry GeometryType="ORIGIN_DXDYDZ">'
        write(xml_file,'(a)') '        <DataItem Dimensions="3" Format="XML" Precision="4">0.0 0.0 0.0</DataItem>'
        write(xml_file,'(a)') '        <DataItem Dimensions="3" Format="XML" Precision="4">1.0 1.0 1.0</DataItem>'
        write(xml_file,'(a)') '      </Geometry>'
        write(xml_file,'(a)') ''
        write(xml_file,'(a)') '      <!-- Scalar field associated with the grid -->'
        write(xml_file,'(3a)') '      <Attribute Name="'//trim(namevarvtk(1))//'" AttributeType="Scalar" Center="Node">'
        write(xml_file,'(7a)') '        <DataItem Format="Binary" NumberType="Float" Precision="4" Dimensions="' // &
         trim(adjustl(strx)) // ' ' // trim(adjustl(stry)) // ' ' // trim(adjustl(strz)) // '">'
        write(xml_file,'(2a)') '          ',trim(sevt2)
        write(xml_file,'(a)') '        </DataItem>'
        write(xml_file,'(a)') '      </Attribute>'
        write(xml_file,'(a)') ''
        write(xml_file,'(a)') '    </Grid>'
        write(xml_file,'(a)') '  </Domain>'
        write(xml_file,'(a)') '</Xdmf>'

        close(xml_file)
        
      endif
#endif
      
      call MPI_Barrier(MPI_COMM_WORLD, e_io)
      if (myrank == 0) then
        inquire(file=trim(sevt1), exist=lexist)
        if (lexist) then
          call MPI_File_delete(trim(sevt1), MPI_INFO_NULL, e_io)
          !write(6,*)'replace ',trim(sevt1)
        endif
      endif
      call MPI_Barrier(MPI_COMM_WORLD, e_io)

      call MPI_FILE_OPEN(MPI_COMM_WORLD, trim(sevt1), &
         MPI_MODE_CREATE + MPI_MODE_WRONLY, &
         MPI_INFO_NULL,fdens,e_io)
      call MPI_File_set_size(fdens, 0_MPI_OFFSET_KIND, e_io)
      tempoffset=int(0,kind=MPI_OFFSET_KIND)

#ifdef PRINTHALF
      call MPI_Type_create_subarray(3,skip_gsizes,skip_lsizes,skip_myoffset, &
         MPI_ORDER_FORTRAN,MPI_UNSIGNED_SHORT,filetypesub,e_io)
#else
      call MPI_Type_create_subarray(3,skip_gsizes,skip_lsizes,skip_myoffset, &
         MPI_ORDER_FORTRAN,MPI_REAL,filetypesub,e_io)
#endif
      call MPI_Type_commit(filetypesub, e_io)

#ifdef PRINTHALF
      call MPI_File_Set_View(fdens,tempoffset,MPI_UNSIGNED_SHORT,filetypesub, &
         "native",MPI_INFO_NULL,e_io)
#else
      call MPI_File_Set_View(fdens,tempoffset,MPI_REAL,filetypesub, &
         "native",MPI_INFO_NULL,e_io)
#endif
      ! We need full local sizes: memDims
      memDims = skip_lsizes + 2*nbuffsub
      memOffs = [ nbuffsub, nbuffsub, nbuffsub ]

#ifdef PRINTHALF
      call MPI_TYPE_CREATE_SUBARRAY(3,memDims,skip_lsizes,memOffs, &
         MPI_ORDER_FORTRAN,MPI_UNSIGNED_SHORT,imemtype,e_io)
#else
      call MPI_TYPE_CREATE_SUBARRAY(3,memDims,skip_lsizes,memOffs, &
         MPI_ORDER_FORTRAN,MPI_REAL,imemtype,e_io)
#endif
      call MPI_TYPE_COMMIT(imemtype,e_io)

      call MPI_FILE_WRITE_ALL(fdens,rhoprint,1,imemtype,MPI_STATUS_IGNORE,e_io)
      call MPI_Type_free(imemtype,   e_io)
      call MPI_Type_free(filetypesub,e_io)
      call MPI_File_sync(fdens, e_io)  
      call MPI_FILE_CLOSE(fdens,e_io)

      if (e_io /= MPI_SUCCESS) then
        call MPI_Error_string(e_io, emsg, elen, ierr)
        write(6,'(A,I6,2A)') 'Rank', myrank, ' MPI_File_close error: ', trim(emsg)
        write(6,'(A)') 'Path tried: '//trim(sevt1)
        call MPI_Abort(MPI_COMM_WORLD, e_io, ierr)
      endif  
      
      call mpi_barrier(MPI_COMM_WORLD,e_io)
      fdens = MPI_FILE_NULL


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!velocity!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      sevt2 = trim(dir_out) // trim(filenamevtk)//'_'//trim(namevarvtk(2))// &
         '_'//trim(write_fmtnumb(iframe)) // '.raw'


#ifdef DOXDMF
      if(myrank==0)then
        sevt1=repeat(' ',mxln)
        sevt1 =  trim(dir_out) // trim(filenamevtk)//'_'//trim(namevarvtk(2))// &
         '_'//trim(write_fmtnumb(iframe)) // '.xdmf'
        open(unit=xml_file, file=trim(sevt1), status="replace", action="write")
        
        sevt1=repeat(' ',mxln)
        sevt1 =  trim(filenamevtk)//'_'//trim(namevarvtk(2))// &
         '_'//trim(write_fmtnumb(iframe)) // '.raw'
        
        write(xml_file,'(a)') '<?xml version="1.0" ?>'
        write(xml_file,'(a)') '<Xdmf Version="3.0">'
        write(xml_file,'(a)') '  <Domain>'
        write(xml_file,'(a)') '    <!-- Uniform Grid -->'
        write(xml_file,'(a)') '    <Grid Name="VelocityField" GridType="Uniform">'
        write(xml_file,'(a)') ''
        write(xml_file,'(a)') '      <!-- Topology: 3D grid -->'
        write(xml_file,'(7a)') '      <Topology TopologyType="3DCoRectMesh" Dimensions="' // &
         trim(adjustl(strx)) // ' ' // trim(adjustl(stry)) // ' ' // trim(adjustl(strz)) // '"/>'
        write(xml_file,'(a)') ''
        write(xml_file,'(a)') '      <!-- Geometry: origin and spacing (dx, dy, dz) -->'
        write(xml_file,'(a)') '      <Geometry GeometryType="ORIGIN_DXDYDZ">'
        write(xml_file,'(a)') '        <DataItem Dimensions="3" Format="XML" Precision="4">0.0 0.0 0.0</DataItem>'
        write(xml_file,'(a)') '        <DataItem Dimensions="3" Format="XML" Precision="4">1.0 1.0 1.0</DataItem>'
        write(xml_file,'(a)') '      </Geometry>'
        write(xml_file,'(a)') ''
        write(xml_file,'(a)') '      <!-- Vector field associated with the grid -->'
        write(xml_file,'(3a)') '      <Attribute Name="'//trim(namevarvtk(2))//'" AttributeType="Vector" Center="Node">'
        write(xml_file,'(7a)') '        <DataItem Format="Binary" NumberType="Float" Precision="4" Dimensions="' // &
         trim(adjustl(strx)) // ' ' // trim(adjustl(stry)) // ' ' // trim(adjustl(strz)) // ' 3">'
        write(xml_file,'(2a)') '          ',trim(sevt1)
        write(xml_file,'(a)') '        </DataItem>'
        write(xml_file,'(a)') '      </Attribute>'
        write(xml_file,'(a)') ''
        write(xml_file,'(a)') '    </Grid>'
        write(xml_file,'(a)') '  </Domain>'
        write(xml_file,'(a)') '</Xdmf>'

        close(xml_file)
        
      endif
#endif

      call MPI_Barrier(MPI_COMM_WORLD, e_io)
      if (myrank == 0) then
        inquire(file=trim(sevt2), exist=lexist)
        if (lexist) then
          call MPI_File_delete(trim(sevt2), MPI_INFO_NULL, e_io)
          !write(6,*)'replace ',trim(sevt2)
        endif
      endif
      call MPI_Barrier(MPI_COMM_WORLD, e_io)
 
      call MPI_FILE_OPEN(MPI_COMM_WORLD,trim(sevt2), &
         MPI_MODE_WRONLY + MPI_MODE_CREATE, &
         MPI_INFO_NULL,fvel,e_io)
      call MPI_File_set_size(fvel, 0_MPI_OFFSET_KIND, e_io)

      tempoffset=int(0,kind=MPI_OFFSET_KIND)


      velglobalDims(1)=3
      velglobalDims(2:4)=skip_gsizes(1:3)
      velldims(1)=3
      velldims(2:4)=skip_lsizes(1:3)
      velmystarts(1) = 0
      velmystarts(2:4) = skip_myoffset(1:3)
#ifdef PRINTHALF
      call MPI_Type_create_subarray(4,velglobalDims,velldims,velmystarts, &
         MPI_ORDER_FORTRAN,MPI_UNSIGNED_SHORT,filetypesubv,e_io)
#else
      call MPI_Type_create_subarray(4,velglobalDims,velldims,velmystarts, &
         MPI_ORDER_FORTRAN,MPI_REAL,filetypesubv,e_io)
#endif
      call MPI_Type_commit(filetypesubv, e_io)
#ifdef PRINTHALF
      call MPI_File_Set_View(fvel,tempoffset,MPI_UNSIGNED_SHORT,filetypesubv, &
         "native",MPI_INFO_NULL,e_io)
#else
      call MPI_File_Set_View(fvel,tempoffset,MPI_REAL,filetypesubv, &
         "native",MPI_INFO_NULL,e_io)
#endif
      ! We need full local sizes: memDims
      velmemDims(1) = vellDims(1)
      velmemDims(2:4) = vellDims(2:4) + 2*nbuffsub
      velmemOffs = [ 0, nbuffsub, nbuffsub, nbuffsub ]
#ifdef PRINTHALF
      call MPI_TYPE_CREATE_SUBARRAY(4,velmemDims,velldims,velmemOffs, &
         MPI_ORDER_FORTRAN,MPI_UNSIGNED_SHORT,imemtype,e_io)
#else
      call MPI_TYPE_CREATE_SUBARRAY(4,velmemDims,velldims,velmemOffs, &
         MPI_ORDER_FORTRAN,MPI_REAL,imemtype,e_io)
#endif
      call MPI_TYPE_COMMIT(imemtype,e_io)

      call MPI_FILE_WRITE_ALL(fvel,velprint,1,imemtype,MPI_STATUS_IGNORE,e_io)
      call MPI_Type_free(imemtype,    e_io)
      call MPI_Type_free(filetypesubv, e_io)
      call MPI_File_sync(fvel, e_io)   
      call MPI_FILE_CLOSE(fvel, e_io)

      if (e_io /= MPI_SUCCESS) then
        call MPI_Error_string(e_io, emsg, elen, ierr)
        write(6,'(A,I6,2A)') 'Rank', myrank, ' MPI_File_close error: ', trim(emsg)
        write(6,'(A)') 'Path tried: '//trim(sevt2)
        call MPI_Abort(MPI_COMM_WORLD, e_io, ierr)
      endif  
      
      call mpi_barrier(MPI_COMM_WORLD,e_io)
      fvel = MPI_FILE_NULL
      
      
#if defined(TWOCOMPONENT) && defined(WRITEPRESS)      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!pressure!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      sevt1 = trim(dir_out) // trim(filenamevtk)//'_'//trim(namevarvtk(3))// &
         '_'//trim(write_fmtnumb(iframe)) // '.raw'

#ifdef DOXDMF
      if(myrank==0)then
        write(strx, '(I0)') skip_gsizes(3)
        write(stry, '(I0)') skip_gsizes(2)
        write(strz, '(I0)') skip_gsizes(1)
        
        sevt2=repeat(' ',mxln)
        sevt2 =  trim(dir_out) // trim(filenamevtk)//'_'//trim(namevarvtk(3))// &
           '_'//trim(write_fmtnumb(iframe)) // '.xdmf'
        open(unit=xml_file, file=trim(sevt2), status="replace", action="write")
        
        sevt2=repeat(' ',mxln)
        sevt2 =  trim(filenamevtk)//'_'//trim(namevarvtk(3))// &
           '_'//trim(write_fmtnumb(iframe)) // '.raw'
        
        write(xml_file,'(a)') '<?xml version="1.0" ?>'
        write(xml_file,'(a)') '<Xdmf Version="3.0">'
        write(xml_file,'(a)') '  <Domain>'
        write(xml_file,'(a)') '    <!-- Uniform Grid -->'
        write(xml_file,'(a)') '    <Grid Name="ScalarField" GridType="Uniform">'
        write(xml_file,'(a)') ''
        write(xml_file,'(a)') '      <!-- Topology: 3D grid -->'
        write(xml_file,'(7a)') '      <Topology TopologyType="3DCoRectMesh" Dimensions="' // &
         trim(adjustl(strx)) // ' ' // trim(adjustl(stry)) // ' ' // trim(adjustl(strz)) // '"/>'
        write(xml_file,'(a)') ''
        write(xml_file,'(a)') '      <!-- Geometry: origin and spacing (dx, dy, dz) -->'
        write(xml_file,'(a)') '      <Geometry GeometryType="ORIGIN_DXDYDZ">'
        write(xml_file,'(a)') '        <DataItem Dimensions="3" Format="XML" Precision="4">0.0 0.0 0.0</DataItem>'
        write(xml_file,'(a)') '        <DataItem Dimensions="3" Format="XML" Precision="4">1.0 1.0 1.0</DataItem>'
        write(xml_file,'(a)') '      </Geometry>'
        write(xml_file,'(a)') ''
        write(xml_file,'(a)') '      <!-- Scalar field associated with the grid -->'
        write(xml_file,'(3a)') '      <Attribute Name="'//trim(namevarvtk(3))//'" AttributeType="Scalar" Center="Node">'
        write(xml_file,'(7a)') '        <DataItem Format="Binary" NumberType="Float" Precision="4" Dimensions="' // &
         trim(adjustl(strx)) // ' ' // trim(adjustl(stry)) // ' ' // trim(adjustl(strz)) // '">'
        write(xml_file,'(2a)') '          ',trim(sevt2)
        write(xml_file,'(a)') '        </DataItem>'
        write(xml_file,'(a)') '      </Attribute>'
        write(xml_file,'(a)') ''
        write(xml_file,'(a)') '    </Grid>'
        write(xml_file,'(a)') '  </Domain>'
        write(xml_file,'(a)') '</Xdmf>'

        close(xml_file)
        
      endif
#endif
      
      call MPI_Barrier(MPI_COMM_WORLD, e_io)
      if (myrank == 0) then
        inquire(file=trim(sevt1), exist=lexist)
        if (lexist) then
          call MPI_File_delete(trim(sevt1), MPI_INFO_NULL, e_io)
          !write(6,*)'replace ',trim(sevt1)
        endif
      endif
      call MPI_Barrier(MPI_COMM_WORLD, e_io)

      call MPI_FILE_OPEN(MPI_COMM_WORLD, trim(sevt1), &
         MPI_MODE_CREATE + MPI_MODE_WRONLY, &
         MPI_INFO_NULL,fdens,e_io)
      call MPI_File_set_size(fdens, 0_MPI_OFFSET_KIND, e_io)
      tempoffset=int(0,kind=MPI_OFFSET_KIND)

#ifdef PRINTHALF
      call MPI_Type_create_subarray(3,skip_gsizes,skip_lsizes,skip_myoffset, &
         MPI_ORDER_FORTRAN,MPI_UNSIGNED_SHORT,filetypesub,e_io)
#else
      call MPI_Type_create_subarray(3,skip_gsizes,skip_lsizes,skip_myoffset, &
         MPI_ORDER_FORTRAN,MPI_REAL,filetypesub,e_io)
#endif
      call MPI_Type_commit(filetypesub, e_io)
#ifdef PRINTHALF
      call MPI_File_Set_View(fdens,tempoffset,MPI_UNSIGNED_SHORT,filetypesub, &
         "native",MPI_INFO_NULL,e_io)
#else
      call MPI_File_Set_View(fdens,tempoffset,MPI_REAL,filetypesub, &
         "native",MPI_INFO_NULL,e_io)
#endif
      ! We need full local sizes: memDims
      memDims = skip_lsizes + 2*nbuffsub
      memOffs = [ nbuffsub, nbuffsub, nbuffsub ]
#ifdef PRINTHALF
      call MPI_TYPE_CREATE_SUBARRAY(3,memDims,skip_lsizes,memOffs, &
         MPI_ORDER_FORTRAN,MPI_UNSIGNED_SHORT,imemtype,e_io)
#else
      call MPI_TYPE_CREATE_SUBARRAY(3,memDims,skip_lsizes,memOffs, &
         MPI_ORDER_FORTRAN,MPI_REAL,imemtype,e_io)
#endif
      call MPI_TYPE_COMMIT(imemtype,e_io)

      call MPI_FILE_WRITE_ALL(fdens,pressprint,1,imemtype,MPI_STATUS_IGNORE,e_io)
      call MPI_Type_free(imemtype,   e_io)
      call MPI_Type_free(filetypesub,e_io)
      call MPI_File_sync(fdens, e_io) 
      call MPI_FILE_CLOSE(fdens,e_io)    

      if (e_io /= MPI_SUCCESS) then
        call MPI_Error_string(e_io, emsg, elen, ierr)
        write(6,'(A,I6,2A)') 'Rank', myrank, ' MPI_File_close error: ', trim(emsg)
        write(6,'(A)') 'Path tried: '//trim(sevt1)
        call MPI_Abort(MPI_COMM_WORLD, e_io, ierr)
      endif  
      
      call mpi_barrier(MPI_COMM_WORLD,e_io)
      fdens = MPI_FILE_NULL
#endif


#endif
      return

   end subroutine write_file_raw_par
   
   subroutine write_file_raw_par_nompiio(iframe,e_io)

!***********************************************************************
!
!     LBsoft subroutine for opening the vtk legacy file
!     in parallel IO
!
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2019
!
!***********************************************************************

      implicit none

      integer, intent(in) ::iframe
      integer, intent(out) :: e_io
#ifdef MPI
      integer(kind=MPI_OFFSET_KIND) :: tempoffset
#endif

      integer :: ioffset
      character(1), parameter :: end_rec = char(10)
      integer, parameter :: bytechar=kind(end_rec)
      integer, parameter :: byteint = 4
      integer, parameter :: byter4  = 4
      integer, parameter :: byter8  = 8
      integer, parameter :: nbuffsub = 0
      integer :: filetypesub,imemtype,filetypesubv,ierr

      integer, dimension(3) :: memDims,memOffs
      integer, dimension(4) :: velglobalDims,velldims,velmystarts, &
         velmemDims,velmemOffs
#ifdef MPI
      integer :: fdens
      integer :: fvel
#ifdef DOXDMF      
      integer, parameter :: xml_file=734
      character(len=32) :: strx,stry,strz
      character(len=32) :: strx1,stry1,strz1
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!density!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      

#ifdef DOXDMF
      if(myrank==0)then
        
        sevt2=repeat(' ',mxln)
        sevt2 =  trim(dir_out) // trim(filenamevtk)//'_'//trim(namevarvtk(1))// &
           '_'//trim(write_fmtnumb(iframe)) // '.xdmf'
        open(unit=xml_file, file=trim(sevt2), status="replace", action="write")
        
        write(xml_file, '(a)') '<?xml version="1.0" ?>'
        write(xml_file, '(a)') '<Xdmf Version="3.0">'
        write(xml_file, '(a)') '  <Domain>'
        write(xml_file, '(a)') '    <Grid Name="'//trim(namevarvtk(1))//'" GridType="Collection" CollectionType="Spatial">'

        do i = 0,nprocs-1
          
          sevt2=repeat(' ',mxln)
          sevt2 = trim(filenamevtk)//'_'//trim(namevarvtk(1))// &
         '_'//trim(write_fmtnumb(iframe)) // '_p'//trim(write_fmtnumb(i)) //'.raw'
          
          write(xml_file, '(a)') '      <Grid Name="Piece' // trim(write_fmtnumb(i)) // '">' 
          !invert index since fortran is column-ordered
          write(strx, '(I0)') lsizes_3d(3,i)
          write(stry, '(I0)') lsizes_3d(2,i)
          write(strz, '(I0)') lsizes_3d(1,i)
          write(xml_file, '(a)') '        <Topology TopologyType="3DCoRectMesh" Dimensions="' // &
           trim(adjustl(strx)) // ' ' // trim(adjustl(stry)) // ' ' // trim(adjustl(strz)) // '"/>'
          write(xml_file, '(a)') '        <Geometry GeometryType="ORIGIN_DXDYDZ">'
          write(strx1, '(I0)') myoffset_3d(3,i)
          write(stry1, '(I0)') myoffset_3d(2,i)
          write(strz1, '(I0)') myoffset_3d(1,i)
          write(xml_file, '(a)') '          <DataItem Dimensions="3" Format="XML">', &
           trim(adjustl(strx1)) // ' ' // trim(adjustl(stry1)) // ' ' // trim(adjustl(strz1)) // '</DataItem>'
          write(xml_file, '(a)') '          <DataItem Dimensions="3" Format="XML">1 1 1</DataItem>'
          write(xml_file, '(a)') '        </Geometry>'
          write(xml_file, '(a)') '        <Attribute Name="'//trim(namevarvtk(1))//'" AttributeType="Scalar" Center="Node">'
          write(xml_file, '(a)') '          <DataItem Format="Binary" NumberType="Float" Precision="4" Dimensions="' // &
           trim(adjustl(strx)) // ' ' // trim(adjustl(stry)) // ' ' // trim(adjustl(strz)) // '">'
          write(xml_file, '(2a)') '            ',trim(sevt2)
          write(xml_file, '(a)') '          </DataItem>'
          write(xml_file, '(a)') '        </Attribute>'
          write(xml_file, '(a)') '      </Grid>'
    end do

    write(xml_file, '(a)') '    </Grid>'
    write(xml_file, '(a)') '  </Domain>'
    write(xml_file, '(a)') '</Xdmf>'
    close(xml_file)
        
        
    endif
    
#endif

      sevt1 = trim(dir_out) // trim(filenamevtk)//'_'//trim(namevarvtk(1))// &
         '_'//trim(write_fmtnumb(iframe)) // '_p'//trim(write_fmtnumb(myrank)) //'.raw'
      
      ! Each process writes its own RAW file if ldowrite is true.
      open(unit=91, file=trim(sevt1), form='unformatted', access='stream', status='replace', iostat=e_io)
      write(91) rhoprint
      close(91)
      
      
      


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!velocity!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      


#ifdef DOXDMF
      if(myrank==0)then
        sevt1=repeat(' ',mxln)
        sevt1 =  trim(dir_out) // trim(filenamevtk)//'_'//trim(namevarvtk(2))// &
         '_'//trim(write_fmtnumb(iframe)) // '.xdmf'
        open(unit=xml_file, file=trim(sevt1), status="replace", action="write")
        
        
        
        write(xml_file, '(a)') '<?xml version="1.0" ?>'
        write(xml_file, '(a)') '<Xdmf Version="3.0">'
        write(xml_file, '(a)') '  <Domain>'
        write(xml_file, '(a)') '    <Grid Name="'//trim(namevarvtk(2))//'" GridType="Collection" CollectionType="Spatial">'

        do i = 0,nprocs-1

         
         sevt1=repeat(' ',mxln)
         sevt1 =  trim(filenamevtk)//'_'//trim(namevarvtk(2))// &
         '_'//trim(write_fmtnumb(iframe)) // '_p'//trim(write_fmtnumb(i)) //'.raw'
          
          write(xml_file, '(a)') '      <Grid Name="Piece' // trim(write_fmtnumb(i)) // '">' 
          !invert index since fortran is column-ordered
          write(strx, '(I0)') lsizes_3d(3,i)
          write(stry, '(I0)') lsizes_3d(2,i)
          write(strz, '(I0)') lsizes_3d(1,i)
          write(xml_file, '(a)') '        <Topology TopologyType="3DCoRectMesh" Dimensions="' // &
           trim(adjustl(strx)) // ' ' // trim(adjustl(stry)) // ' ' // trim(adjustl(strz)) // '"/>'
          write(xml_file, '(a)') '        <Geometry GeometryType="ORIGIN_DXDYDZ">'
          write(strx1, '(I0)') myoffset_3d(3,i)
          write(stry1, '(I0)') myoffset_3d(2,i)
          write(strz1, '(I0)') myoffset_3d(1,i)
          write(xml_file, '(a)') '          <DataItem Dimensions="3" Format="XML">', &
           trim(adjustl(strx1)) // ' ' // trim(adjustl(stry1)) // ' ' // trim(adjustl(strz1)) // '</DataItem>'
          write(xml_file, '(a)') '          <DataItem Dimensions="3" Format="XML">1 1 1</DataItem>'
          write(xml_file, '(a)') '        </Geometry>'
          write(xml_file, '(a)') '        <Attribute Name="'//trim(namevarvtk(2))//'" AttributeType="Vector" Center="Node">'
          write(xml_file, '(a)') '          <DataItem Format="Binary" NumberType="Float" Precision="4" Dimensions="' // &
           trim(adjustl(strx)) // ' ' // trim(adjustl(stry)) // ' ' // trim(adjustl(strz)) // ' 3">'
          write(xml_file, '(2a)') '            ',trim(sevt1)
          write(xml_file, '(a)') '          </DataItem>'
          write(xml_file, '(a)') '        </Attribute>'
          write(xml_file, '(a)') '      </Grid>'
    end do

    write(xml_file, '(a)') '    </Grid>'
    write(xml_file, '(a)') '  </Domain>'
    write(xml_file, '(a)') '</Xdmf>'
    close(xml_file)        
        
    endif
#endif

      sevt2 = trim(dir_out) // trim(filenamevtk)//'_'//trim(namevarvtk(2))// &
         '_'//trim(write_fmtnumb(iframe)) // '_p'//trim(write_fmtnumb(myrank)) //'.raw'
      ! Each process writes its own RAW file if ldowrite is true.
      open(unit=92, file=trim(sevt2), form='unformatted', access='stream', status='replace', iostat=e_io)
      write(92) velprint
      close(92)
      
      
      
#if defined(TWOCOMPONENT) && defined(WRITEPRESS)      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!pressure!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      

#ifdef DOXDMF
      if(myrank==0)then
        write(strx, '(I0)') gsizes(3)
        write(stry, '(I0)') gsizes(2)
        write(strz, '(I0)') gsizes(1)
        
        sevt2=repeat(' ',mxln)
        sevt2 =  trim(dir_out) // trim(filenamevtk)//'_'//trim(namevarvtk(3))// &
           '_'//trim(write_fmtnumb(iframe)) // '.xdmf'
        open(unit=xml_file, file=trim(sevt2), status="replace", action="write")
                write(xml_file, '(a)') '<?xml version="1.0" ?>'
        write(xml_file, '(a)') '<Xdmf Version="3.0">'
        write(xml_file, '(a)') '  <Domain>'
        write(xml_file, '(a)') '    <Grid Name="'//trim(namevarvtk(3))//'" GridType="Collection" CollectionType="Spatial">'

        do i = 0,nprocs-1
          
          sevt2=repeat(' ',mxln)
          sevt2 = trim(filenamevtk)//'_'//trim(namevarvtk(3))// &
           '_'//trim(write_fmtnumb(iframe)) // '_p'//trim(write_fmtnumb(i)) //'.raw'
          
          write(xml_file, '(a)') '      <Grid Name="Piece' // trim(write_fmtnumb(i)) // '">' 
          !invert index since fortran is column-ordered
          write(strx, '(I0)') lsizes_3d(3,i)
          write(stry, '(I0)') lsizes_3d(2,i)
          write(strz, '(I0)') lsizes_3d(1,i)
          write(xml_file, '(a)') '        <Topology TopologyType="3DCoRectMesh" Dimensions="' // &
           trim(adjustl(strx)) // ' ' // trim(adjustl(stry)) // ' ' // trim(adjustl(strz)) // '"/>'
          write(xml_file, '(a)') '        <Geometry GeometryType="ORIGIN_DXDYDZ">'
          write(strx1, '(I0)') myoffset_3d(3,i)
          write(stry1, '(I0)') myoffset_3d(2,i)
          write(strz1, '(I0)') myoffset_3d(1,i)
          write(xml_file, '(a)') '          <DataItem Dimensions="3" Format="XML">', &
           trim(adjustl(strx1)) // ' ' // trim(adjustl(stry1)) // ' ' // trim(adjustl(strz1)) // '</DataItem>'
          write(xml_file, '(a)') '          <DataItem Dimensions="3" Format="XML">1 1 1</DataItem>'
          write(xml_file, '(a)') '        </Geometry>'
          write(xml_file, '(a)') '        <Attribute Name="'//trim(namevarvtk(3))//'" AttributeType="Scalar" Center="Node">'
          write(xml_file, '(a)') '          <DataItem Format="Binary" NumberType="Float" Precision="4" Dimensions="' // &
           trim(adjustl(strx)) // ' ' // trim(adjustl(stry)) // ' ' // trim(adjustl(strz)) // '">'
          write(xml_file, '(2a)') '            ',trim(sevt2)
          write(xml_file, '(a)') '          </DataItem>'
          write(xml_file, '(a)') '        </Attribute>'
          write(xml_file, '(a)') '      </Grid>'
    end do

    write(xml_file, '(a)') '    </Grid>'
    write(xml_file, '(a)') '  </Domain>'
    write(xml_file, '(a)') '</Xdmf>'
    close(xml_file)
        
    endif
#endif

      sevt1 = trim(dir_out) // trim(filenamevtk)//'_'//trim(namevarvtk(3))// &
         '_'//trim(write_fmtnumb(iframe)) // '_p'//trim(write_fmtnumb(myrank)) //'.raw'
      ! Each process writes its own RAW file if ldowrite is true.
      open(unit=93, file=trim(sevt1), form='unformatted', access='stream', status='replace', iostat=e_io)
      write(93) pressprint
      close(93)
  
#endif


#endif
      return

   end subroutine write_file_raw_par_nompiio
   
   subroutine write_file_raw_par2D(ldowrite,iframe,myid,mydir,mypoint,service1,service3, &
      myoffset_plane,lsizes_plane,gsizes_plane,e_io,service2)

!***********************************************************************
!
!     LBsoft subroutine for opening the vtk legacy file
!     in parallel IO
!
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2019
!
!***********************************************************************

      implicit none
      
      logical, intent(in) :: ldowrite
      integer, intent(in) :: myid,iframe,mydir,mypoint
      real(4), dimension(:,:,:), allocatable :: service1
      real(4), dimension(:,:,:,:), allocatable :: service3
      integer, dimension(mpid), intent(in) :: myoffset_plane,lsizes_plane,gsizes_plane
      integer, intent(out) :: e_io
      real(4), dimension(:,:,:), allocatable, optional :: service2
#ifdef MPI
      integer(kind=MPI_OFFSET_KIND) :: tempoffset
#endif

      integer :: ioffset
      character(1), parameter :: end_rec = char(10)
      integer, parameter :: bytechar=kind(end_rec)
      integer, parameter :: byteint = 4
      integer, parameter :: byter4  = 4
      integer, parameter :: byter8  = 8
      integer, parameter :: nbuffsub = 0
      integer :: filetypesub,imemtype,filetypesubv,ierr

      integer, dimension(3) :: memDims,memOffs
      integer, dimension(4) :: velglobalDims,velldims,velmystarts, &
         velmemDims,velmemOffs
      integer :: elen,amode
      logical :: lexist
#ifdef MPI
      integer :: fdens=MPI_FILE_NULL
      integer :: fvel=MPI_FILE_NULL
      character(len=MPI_MAX_ERROR_STRING) :: emsg

#ifdef DOXDMF      
      integer, parameter :: xml_file=734
      character(len=32) :: strx,stry,strz
#endif
    
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!density!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      sevt1=repeat(' ',mxln)
      sevt1 = trim(dir_out) // trim(filenamevtk)//'_'//trim(adjustl(space_fmtnumb(mydir)))// &
         '_'//trim(adjustl(space_fmtnumb(mypoint)))//'_'//trim(namevarvtk(1))// &
         '_'//trim(write_fmtnumb(iframe)) // '.raw'
         
#ifdef DOXDMF
      if(myrank==0)then
        write(strx, '(I0)') gsizes_plane(3)
        write(stry, '(I0)') gsizes_plane(2)
        write(strz, '(I0)') gsizes_plane(1)
        
        sevt2=repeat(' ',mxln)
        sevt2 =  trim(dir_out) // trim(filenamevtk)//'_'//trim(adjustl(space_fmtnumb(mydir)))// &
         '_'//trim(adjustl(space_fmtnumb(mypoint)))//'_'//trim(namevarvtk(1))// &
         '_'//trim(write_fmtnumb(iframe)) // '.xdmf'
        open(unit=xml_file, file=trim(sevt2), status="replace", action="write")
        
        sevt2=repeat(' ',mxln)
        sevt2 =  trim(filenamevtk)//'_'//trim(adjustl(space_fmtnumb(mydir)))// &
         '_'//trim(adjustl(space_fmtnumb(mypoint)))//'_'//trim(namevarvtk(1))// &
         '_'//trim(write_fmtnumb(iframe)) // '.raw'
        
        write(xml_file,'(a)') '<?xml version="1.0" ?>'
        write(xml_file,'(a)') '<Xdmf Version="3.0">'
        write(xml_file,'(a)') '  <Domain>'
        write(xml_file,'(a)') '    <!-- Uniform Grid -->'
        write(xml_file,'(a)') '    <Grid Name="ScalarField" GridType="Uniform">'
        write(xml_file,'(a)') ''
        write(xml_file,'(a)') '      <!-- Topology: 3D grid -->'
        write(xml_file,'(7a)') '      <Topology TopologyType="3DCoRectMesh" Dimensions="' // &
         trim(adjustl(strx)) // ' ' // trim(adjustl(stry)) // ' ' // trim(adjustl(strz)) // '"/>'
        write(xml_file,'(a)') ''
        write(xml_file,'(a)') '      <!-- Geometry: origin and spacing (dx, dy, dz) -->'
        write(xml_file,'(a)') '      <Geometry GeometryType="ORIGIN_DXDYDZ">'
        write(xml_file,'(a)') '        <DataItem Dimensions="3" Format="XML" Precision="4">0.0 0.0 0.0</DataItem>'
        write(xml_file,'(a)') '        <DataItem Dimensions="3" Format="XML" Precision="4">1.0 1.0 1.0</DataItem>'
        write(xml_file,'(a)') '      </Geometry>'
        write(xml_file,'(a)') ''
        write(xml_file,'(a)') '      <!-- Scalar field associated with the grid -->'
        write(xml_file,'(3a)') '      <Attribute Name="'//trim(namevarvtk(1))//'" AttributeType="Scalar" Center="Node">'
        write(xml_file,'(7a)') '        <DataItem Format="Binary" NumberType="Float" Precision="4" Dimensions="' // &
         trim(adjustl(strx)) // ' ' // trim(adjustl(stry)) // ' ' // trim(adjustl(strz)) // '">'
        write(xml_file,'(2a)') '          ',trim(sevt2)
        write(xml_file,'(a)') '        </DataItem>'
        write(xml_file,'(a)') '      </Attribute>'
        write(xml_file,'(a)') ''
        write(xml_file,'(a)') '    </Grid>'
        write(xml_file,'(a)') '  </Domain>'
        write(xml_file,'(a)') '</Xdmf>'

        close(xml_file)
        
      endif
#endif

      call MPI_Barrier(MPI_COMM_WORLD, e_io)
      if (myrank == 0) then
        inquire(file=trim(sevt1), exist=lexist)
        if (lexist) then
          call MPI_File_delete(trim(sevt1), MPI_INFO_NULL, e_io)
          !write(6,*)'replace ',trim(sevt1)
        endif
      endif
      call MPI_Barrier(MPI_COMM_WORLD, e_io)
     
    if(ldowrite)then
      
      call MPI_Barrier(io_comm2d(myid), e_io)
      
      amode = IOR(MPI_MODE_CREATE, MPI_MODE_WRONLY) 

      call MPI_FILE_OPEN(io_comm2d(myid), trim(sevt1),amode, &
         MPI_INFO_NULL,fdens,e_io)
         
      if (e_io /= MPI_SUCCESS) then
        call MPI_Error_string(e_io, emsg, elen, ierr)
        write(6,'(A,I6,2A)') 'Rank', myrank, ' MPI_File_open error: ', trim(emsg)
        write(6,'(A)') 'Path tried: '//trim(sevt1)
        call MPI_Abort(io_comm2d(myid), e_io, ierr)
      endif
         
      call MPI_File_set_size(fdens, 0_MPI_OFFSET_KIND, e_io)
      call MPI_Barrier(io_comm2d(myid), e_io)  
      tempoffset=int(0,kind=MPI_OFFSET_KIND)

#ifdef PRINTHALF
      call MPI_Type_create_subarray(3,gsizes_plane,lsizes_plane,myoffset_plane, &
         MPI_ORDER_FORTRAN,MPI_UNSIGNED_SHORT,filetypesub,e_io)
#else
      call MPI_Type_create_subarray(3,gsizes_plane,lsizes_plane,myoffset_plane, &
         MPI_ORDER_FORTRAN,MPI_REAL,filetypesub,e_io)
#endif
      call MPI_Type_commit(filetypesub, e_io)
#ifdef PRINTHALF
      call MPI_File_Set_View(fdens,tempoffset,MPI_UNSIGNED_SHORT,filetypesub, &
         "native",MPI_INFO_NULL,e_io)
#else
      call MPI_File_Set_View(fdens,tempoffset,MPI_REAL,filetypesub, &
         "native",MPI_INFO_NULL,e_io)
#endif
      ! We need full local sizes: memDims
      memDims = lsizes_plane + 2*nbuffsub
      memOffs = [ nbuffsub, nbuffsub, nbuffsub ]
#ifdef PRINTHALF
      call MPI_TYPE_CREATE_SUBARRAY(3,memDims,lsizes_plane,memOffs, &
         MPI_ORDER_FORTRAN,MPI_UNSIGNED_SHORT,imemtype,e_io)
#else
      call MPI_TYPE_CREATE_SUBARRAY(3,memDims,lsizes_plane,memOffs, &
         MPI_ORDER_FORTRAN,MPI_REAL,imemtype,e_io)
#endif
      call MPI_TYPE_COMMIT(imemtype,e_io)

      call MPI_FILE_WRITE_ALL(fdens,service1,1,imemtype,MPI_STATUS_IGNORE,e_io)
      call MPI_Type_free(imemtype,   e_io)
      call MPI_Type_free(filetypesub,e_io)
      call MPI_File_sync(fdens, e_io)
      
      call MPI_Barrier(io_comm2d(myid), ierr)
      call MPI_FILE_CLOSE(fdens, e_io)
      
      if (e_io /= MPI_SUCCESS) then
        call MPI_Error_string(e_io, emsg, elen, ierr)
        write(6,'(A,I6,2A)') 'Rank', myrank, ' MPI_File_close error: ', trim(emsg)
        write(6,'(A)') 'Path tried: '//trim(sevt1)
        call MPI_Abort(io_comm2d(myid), e_io, ierr)
      endif  
      
      call mpi_barrier(io_comm2d(myid),e_io)
      fdens = MPI_FILE_NULL

    endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!velocity!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call MPI_Barrier(MPI_COMM_WORLD, e_io)
      if (myrank == 0) then
        inquire(file=trim(sevt2), exist=lexist)
        if (lexist) then
          call MPI_File_delete(trim(sevt2), MPI_INFO_NULL, e_io)
          !write(6,*)'replace ',trim(sevt2)
        endif
      endif
      call MPI_Barrier(MPI_COMM_WORLD, e_io)
      
      sevt2=repeat(' ',mxln)
      sevt2 = trim(dir_out) // trim(filenamevtk)//'_'//trim(adjustl(space_fmtnumb(mydir)))// &
         '_'//trim(adjustl(space_fmtnumb(mypoint)))//'_'//trim(namevarvtk(2))// &
         '_'//trim(write_fmtnumb(iframe)) // '.raw'

#ifdef DOXDMF
      if(myrank==0)then
        sevt1=repeat(' ',mxln)
        sevt1 =  trim(dir_out) // trim(filenamevtk)//'_'//trim(adjustl(space_fmtnumb(mydir)))// &
         '_'//trim(adjustl(space_fmtnumb(mypoint)))//'_'//trim(namevarvtk(2))// &
         '_'//trim(write_fmtnumb(iframe)) // '.xdmf'
        open(unit=xml_file, file=trim(sevt1), status="replace", action="write")
        
        sevt1=repeat(' ',mxln)
        sevt1 =  trim(filenamevtk)//'_'//trim(adjustl(space_fmtnumb(mydir)))// &
         '_'//trim(adjustl(space_fmtnumb(mypoint)))//'_'//trim(namevarvtk(2))// &
         '_'//trim(write_fmtnumb(iframe)) // '.raw'
        
        write(xml_file,'(a)') '<?xml version="1.0" ?>'
        write(xml_file,'(a)') '<Xdmf Version="3.0">'
        write(xml_file,'(a)') '  <Domain>'
        write(xml_file,'(a)') '    <!-- Uniform Grid -->'
        write(xml_file,'(a)') '    <Grid Name="VelocityField" GridType="Uniform">'
        write(xml_file,'(a)') ''
        write(xml_file,'(a)') '      <!-- Topology: 3D grid -->'
        write(xml_file,'(7a)') '      <Topology TopologyType="3DCoRectMesh" Dimensions="' // &
         trim(adjustl(strx)) // ' ' // trim(adjustl(stry)) // ' ' // trim(adjustl(strz)) // '"/>'
        write(xml_file,'(a)') ''
        write(xml_file,'(a)') '      <!-- Geometry: origin and spacing (dx, dy, dz) -->'
        write(xml_file,'(a)') '      <Geometry GeometryType="ORIGIN_DXDYDZ">'
        write(xml_file,'(a)') '        <DataItem Dimensions="3" Format="XML" Precision="4">0.0 0.0 0.0</DataItem>'
        write(xml_file,'(a)') '        <DataItem Dimensions="3" Format="XML" Precision="4">1.0 1.0 1.0</DataItem>'
        write(xml_file,'(a)') '      </Geometry>'
        write(xml_file,'(a)') ''
        write(xml_file,'(a)') '      <!-- Vector field associated with the grid -->'
        write(xml_file,'(3a)') '      <Attribute Name="'//trim(namevarvtk(2))//'" AttributeType="Vector" Center="Node">'
        write(xml_file,'(7a)') '        <DataItem Format="Binary" NumberType="Float" Precision="4" Dimensions="' // &
         trim(adjustl(strx)) // ' ' // trim(adjustl(stry)) // ' ' // trim(adjustl(strz)) // ' 3">'
        write(xml_file,'(2a)') '          ',trim(sevt1)
        write(xml_file,'(a)') '        </DataItem>'
        write(xml_file,'(a)') '      </Attribute>'
        write(xml_file,'(a)') ''
        write(xml_file,'(a)') '    </Grid>'
        write(xml_file,'(a)') '  </Domain>'
        write(xml_file,'(a)') '</Xdmf>'

        close(xml_file)
        
      endif
#endif
    if(ldowrite)then
      call MPI_FILE_OPEN(io_comm2d(myid),trim(sevt2), &
         MPI_MODE_WRONLY + MPI_MODE_CREATE, &
         MPI_INFO_NULL,fvel,e_io)

      call MPI_File_set_size(fvel, 0_MPI_OFFSET_KIND, e_io)
      tempoffset=int(0,kind=MPI_OFFSET_KIND)


      velglobalDims(1)=3
      velglobalDims(2:4)=gsizes_plane(1:3)
      velldims(1)=3
      velldims(2:4)=lsizes_plane(1:3)
      velmystarts(1) = 0
      velmystarts(2:4) = myoffset_plane(1:3)
#ifdef PRINTHALF
      call MPI_Type_create_subarray(4,velglobalDims,velldims,velmystarts, &
         MPI_ORDER_FORTRAN,MPI_UNSIGNED_SHORT,filetypesubv,e_io)
#else
      call MPI_Type_create_subarray(4,velglobalDims,velldims,velmystarts, &
         MPI_ORDER_FORTRAN,MPI_REAL,filetypesubv,e_io)
#endif
      call MPI_Type_commit(filetypesubv, e_io)
#ifdef PRINTHALF
      call MPI_File_Set_View(fvel,tempoffset,MPI_UNSIGNED_SHORT,filetypesubv, &
         "native",MPI_INFO_NULL,e_io)
#else
      call MPI_File_Set_View(fvel,tempoffset,MPI_REAL,filetypesubv, &
         "native",MPI_INFO_NULL,e_io)
#endif
      ! We need full local sizes: memDims
      velmemDims(1) = vellDims(1)
      velmemDims(2:4) = vellDims(2:4) + 2*nbuffsub
      velmemOffs = [ 0, nbuffsub, nbuffsub, nbuffsub ]
#ifdef PRINTHALF
      call MPI_TYPE_CREATE_SUBARRAY(4,velmemDims,velldims,velmemOffs, &
         MPI_ORDER_FORTRAN,MPI_UNSIGNED_SHORT,imemtype,e_io)
#else
      call MPI_TYPE_CREATE_SUBARRAY(4,velmemDims,velldims,velmemOffs, &
         MPI_ORDER_FORTRAN,MPI_REAL,imemtype,e_io)
#endif
      call MPI_TYPE_COMMIT(imemtype,e_io)

      call MPI_FILE_WRITE_ALL(fvel,service3,1,imemtype,MPI_STATUS_IGNORE,e_io)
      call MPI_Type_free(imemtype, e_io) 
      call MPI_Type_free(filetypesubv,e_io) 
      call MPI_File_sync(fvel, e_io)

      call MPI_FILE_CLOSE(fvel, e_io)

      if (e_io /= MPI_SUCCESS) then
        call MPI_Error_string(e_io, emsg, elen, ierr)
        write(6,'(A,I6,2A)') 'Rank', myrank, ' MPI_File_close error: ', trim(emsg)
        write(6,'(A)') 'Path tried: '//trim(sevt2)
        call MPI_Abort(io_comm2d(myid), e_io, ierr)
      endif  
      call mpi_barrier(io_comm2d(myid),e_io)
      fvel = MPI_FILE_NULL  
    endif
#if defined(TWOCOMPONENT) && defined(WRITEPRESS)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!pressure!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      sevt1=repeat(' ',mxln)
      sevt1 = trim(dir_out) // trim(filenamevtk)//'_'//trim(adjustl(space_fmtnumb(mydir)))// &
         '_'//trim(adjustl(space_fmtnumb(mypoint)))//'_'//trim(namevarvtk(3))// &
         '_'//trim(write_fmtnumb(iframe)) // '.raw'
         
#ifdef DOXDMF
      if(myrank==0)then
        write(strx, '(I0)') gsizes_plane(3)
        write(stry, '(I0)') gsizes_plane(2)
        write(strz, '(I0)') gsizes_plane(1)
        
        sevt2=repeat(' ',mxln)
        sevt2 =  trim(dir_out) // trim(filenamevtk)//'_'//trim(adjustl(space_fmtnumb(mydir)))// &
         '_'//trim(adjustl(space_fmtnumb(mypoint)))//'_'//trim(namevarvtk(3))// &
         '_'//trim(write_fmtnumb(iframe)) // '.xdmf'
        open(unit=xml_file, file=trim(sevt2), status="replace", action="write")
        
        sevt2=repeat(' ',mxln)
        sevt2 =  trim(filenamevtk)//'_'//trim(adjustl(space_fmtnumb(mydir)))// &
         '_'//trim(adjustl(space_fmtnumb(mypoint)))//'_'//trim(namevarvtk(3))// &
         '_'//trim(write_fmtnumb(iframe)) // '.raw'
        
        write(xml_file,'(a)') '<?xml version="1.0" ?>'
        write(xml_file,'(a)') '<Xdmf Version="3.0">'
        write(xml_file,'(a)') '  <Domain>'
        write(xml_file,'(a)') '    <!-- Uniform Grid -->'
        write(xml_file,'(a)') '    <Grid Name="ScalarField" GridType="Uniform">'
        write(xml_file,'(a)') ''
        write(xml_file,'(a)') '      <!-- Topology: 3D grid -->'
        write(xml_file,'(7a)') '      <Topology TopologyType="3DCoRectMesh" Dimensions="' // &
         trim(adjustl(strx)) // ' ' // trim(adjustl(stry)) // ' ' // trim(adjustl(strz)) // '"/>'
        write(xml_file,'(a)') ''
        write(xml_file,'(a)') '      <!-- Geometry: origin and spacing (dx, dy, dz) -->'
        write(xml_file,'(a)') '      <Geometry GeometryType="ORIGIN_DXDYDZ">'
        write(xml_file,'(a)') '        <DataItem Dimensions="3" Format="XML" Precision="4">0.0 0.0 0.0</DataItem>'
        write(xml_file,'(a)') '        <DataItem Dimensions="3" Format="XML" Precision="4">1.0 1.0 1.0</DataItem>'
        write(xml_file,'(a)') '      </Geometry>'
        write(xml_file,'(a)') ''
        write(xml_file,'(a)') '      <!-- Scalar field associated with the grid -->'
        write(xml_file,'(3a)') '      <Attribute Name="'//trim(namevarvtk(3))//'" AttributeType="Scalar" Center="Node">'
        write(xml_file,'(7a)') '        <DataItem Format="Binary" NumberType="Float" Precision="4" Dimensions="' // &
         trim(adjustl(strx)) // ' ' // trim(adjustl(stry)) // ' ' // trim(adjustl(strz)) // '">'
        write(xml_file,'(2a)') '          ',trim(sevt2)
        write(xml_file,'(a)') '        </DataItem>'
        write(xml_file,'(a)') '      </Attribute>'
        write(xml_file,'(a)') ''
        write(xml_file,'(a)') '    </Grid>'
        write(xml_file,'(a)') '  </Domain>'
        write(xml_file,'(a)') '</Xdmf>'

        close(xml_file)
        
      endif
#endif
    
      call MPI_Barrier(MPI_COMM_WORLD, e_io)
      if (myrank == 0) then
        inquire(file=trim(sevt1), exist=lexist)
        if (lexist) then
          call MPI_File_delete(trim(sevt1), MPI_INFO_NULL, e_io)
          !write(6,*)'replace ',trim(sevt1)
        endif
      endif
      call MPI_Barrier(MPI_COMM_WORLD, e_io)
    
    if(ldowrite)then
      call MPI_FILE_OPEN(io_comm2d(myid), trim(sevt1), &
         MPI_MODE_CREATE + MPI_MODE_WRONLY, &
         MPI_INFO_NULL,fdens,e_io)
      call MPI_File_set_size(fdens, 0_MPI_OFFSET_KIND, e_io)
      tempoffset=int(0,kind=MPI_OFFSET_KIND)

#ifdef PRINTHALF
      call MPI_Type_create_subarray(3,gsizes_plane,lsizes_plane,myoffset_plane, &
         MPI_ORDER_FORTRAN,MPI_UNSIGNED_SHORT,filetypesub,e_io)
#else
      call MPI_Type_create_subarray(3,gsizes_plane,lsizes_plane,myoffset_plane, &
         MPI_ORDER_FORTRAN,MPI_REAL,filetypesub,e_io)
#endif
      call MPI_Type_commit(filetypesub, e_io)
#ifdef PRINTHALF
      call MPI_File_Set_View(fdens,tempoffset,MPI_UNSIGNED_SHORT,filetypesub, &
         "native",MPI_INFO_NULL,e_io)
#else
      call MPI_File_Set_View(fdens,tempoffset,MPI_REAL,filetypesub, &
         "native",MPI_INFO_NULL,e_io)
#endif
      ! We need full local sizes: memDims
      memDims = lsizes_plane + 2*nbuffsub
      memOffs = [ nbuffsub, nbuffsub, nbuffsub ]
#ifdef PRINTHALF
      call MPI_TYPE_CREATE_SUBARRAY(3,memDims,lsizes_plane,memOffs, &
         MPI_ORDER_FORTRAN,MPI_UNSIGNED_SHORT,imemtype,e_io)
#else
      call MPI_TYPE_CREATE_SUBARRAY(3,memDims,lsizes_plane,memOffs, &
         MPI_ORDER_FORTRAN,MPI_REAL,imemtype,e_io)
#endif
      call MPI_TYPE_COMMIT(imemtype,e_io)

      call MPI_FILE_WRITE_ALL(fdens,service2,1,imemtype,MPI_STATUS_IGNORE,e_io)
      call MPI_Type_free(imemtype, e_io) 
      call MPI_Type_free(filetypesub,e_io) 
      call MPI_File_sync(fdens, e_io)
      call MPI_FILE_CLOSE(fdens,e_io)
      if (e_io /= MPI_SUCCESS) then
        call MPI_Error_string(e_io, emsg, elen, ierr)
        write(6,'(A,I6,2A)') 'Rank', myrank, ' MPI_File_close error: ', trim(emsg)
        write(6,'(A)') 'Path tried: '//trim(sevt1)
        call MPI_Abort(io_comm2d(myid), e_io, ierr)
      endif
      call mpi_barrier(io_comm2d(myid),e_io)
      fdens = MPI_FILE_NULL  
    endif
#endif

#endif
      return

   end subroutine write_file_raw_par2D
   
   subroutine write_file_raw_par2D_nompiio(ldowrite,iframe,myid,mydir,mypoint,service1,service3, &
      myoffset_plane,lsizes_plane,gsizes_plane,e_io,service2)

!***********************************************************************
!
!     LBsoft subroutine for opening the vtk legacy file
!     in parallel IO
!
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2019
!
!***********************************************************************

      implicit none
      
      logical, intent(in) :: ldowrite
      integer, intent(in) :: myid,iframe,mydir,mypoint
      real(4), dimension(:,:,:), allocatable :: service1
      real(4), dimension(:,:,:,:), allocatable :: service3
      integer, dimension(mpid), intent(in) :: myoffset_plane,lsizes_plane,gsizes_plane
      integer, intent(out) :: e_io
      real(4), dimension(:,:,:), allocatable, optional :: service2
#ifdef MPI
      integer(kind=MPI_OFFSET_KIND) :: tempoffset
#endif

      integer :: ioffset
      character(1), parameter :: end_rec = char(10)
      integer, parameter :: bytechar=kind(end_rec)
      integer, parameter :: byteint = 4
      integer, parameter :: byter4  = 4
      integer, parameter :: byter8  = 8
      integer, parameter :: nbuffsub = 0
      integer :: filetypesub,imemtype,filetypesubv,ierr

      integer, dimension(3) :: memDims,memOffs
      integer, dimension(4) :: velglobalDims,velldims,velmystarts, &
         velmemDims,velmemOffs
#ifdef MPI
      integer :: fdens
      integer :: fvel

#ifdef DOXDMF      
      integer, parameter :: xml_file=734
      character(len=32) :: strx,stry,strz
      character(len=32) :: strx1,stry1,strz1
#endif

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!density!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         
#ifdef DOXDMF
      if(myrank==0)then
        
        sevt2=repeat(' ',mxln)
        sevt2 =  trim(dir_out) // trim(filenamevtk)//'_'//trim(adjustl(space_fmtnumb(mydir)))// &
         '_'//trim(adjustl(space_fmtnumb(mypoint)))//'_'//trim(namevarvtk(1))// &
         '_'//trim(write_fmtnumb(iframe)) // '.xdmf'
        open(unit=xml_file, file=trim(sevt2), status="replace", action="write")
        
        
        write(xml_file, '(a)') '<?xml version="1.0" ?>'
        write(xml_file, '(a)') '<Xdmf Version="3.0">'
        write(xml_file, '(a)') '  <Domain>'
        write(xml_file, '(a)') '    <Grid Name="'//trim(namevarvtk(1))//'" GridType="Collection" CollectionType="Spatial">'

        do i = 0,nprocs-1

          if (.not. log_plane2d(i,myid))cycle
          
          sevt2=repeat(' ',mxln)
          sevt2 = trim(filenamevtk)//'_'//trim(adjustl(space_fmtnumb(mydir)))// &
         '_'//trim(adjustl(space_fmtnumb(mypoint)))//'_'//trim(namevarvtk(1))// &
         '_'//trim(write_fmtnumb(iframe)) // '_p'//trim(write_fmtnumb(i)) //'.raw'
          
          write(xml_file, '(a)') '      <Grid Name="Piece' // trim(write_fmtnumb(i)) // '">' 
          !invert index since fortran is column-ordered
          write(strx, '(I0)') lsizes_plane2d(3,i,myid)
          write(stry, '(I0)') lsizes_plane2d(2,i,myid)
          write(strz, '(I0)') lsizes_plane2d(1,i,myid)
          write(xml_file, '(a)') '        <Topology TopologyType="3DCoRectMesh" Dimensions="' // &
           trim(adjustl(strx)) // ' ' // trim(adjustl(stry)) // ' ' // trim(adjustl(strz)) // '"/>'
          write(xml_file, '(a)') '        <Geometry GeometryType="ORIGIN_DXDYDZ">'
          write(strx1, '(I0)') myoffset_plane2d(3,i,myid)
          write(stry1, '(I0)') myoffset_plane2d(2,i,myid)
          write(strz1, '(I0)') myoffset_plane2d(1,i,myid)
          write(xml_file, '(a)') '          <DataItem Dimensions="3" Format="XML">', &
           trim(adjustl(strx1)) // ' ' // trim(adjustl(stry1)) // ' ' // trim(adjustl(strz1)) // '</DataItem>'
          write(xml_file, '(a)') '          <DataItem Dimensions="3" Format="XML">1 1 1</DataItem>'
          write(xml_file, '(a)') '        </Geometry>'
          write(xml_file, '(a)') '        <Attribute Name="'//trim(namevarvtk(1))//'" AttributeType="Scalar" Center="Node">'
          write(xml_file, '(a)') '          <DataItem Format="Binary" NumberType="Float" Precision="4" Dimensions="' // &
           trim(adjustl(strx)) // ' ' // trim(adjustl(stry)) // ' ' // trim(adjustl(strz)) // '">'
          write(xml_file, '(2a)') '            ',trim(sevt2)
          write(xml_file, '(a)') '          </DataItem>'
          write(xml_file, '(a)') '        </Attribute>'
          write(xml_file, '(a)') '      </Grid>'
    end do

    write(xml_file, '(a)') '    </Grid>'
    write(xml_file, '(a)') '  </Domain>'
    write(xml_file, '(a)') '</Xdmf>'
    close(xml_file)
        
    endif
#endif
     
    if(ldowrite)then
      sevt1=repeat(' ',mxln)
      sevt1 = trim(dir_out) // trim(filenamevtk)//'_'//trim(adjustl(space_fmtnumb(mydir)))// &
         '_'//trim(adjustl(space_fmtnumb(mypoint)))//'_'//trim(namevarvtk(1))// &
         '_'//trim(write_fmtnumb(iframe)) // '_p'//trim(write_fmtnumb(myrank)) //'.raw'
      ! Each process writes its own RAW file if ldowrite is true.
      open(unit=81, file=trim(sevt1), form='unformatted', access='stream', status='replace', iostat=e_io)
      write(81) service1
      close(81)
    endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!velocity!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef DOXDMF
      if(myrank==0)then
        sevt1=repeat(' ',mxln)
        sevt1 =  trim(dir_out) // trim(filenamevtk)//'_'//trim(adjustl(space_fmtnumb(mydir)))// &
         '_'//trim(adjustl(space_fmtnumb(mypoint)))//'_'//trim(namevarvtk(2))// &
         '_'//trim(write_fmtnumb(iframe)) // '.xdmf'
        open(unit=xml_file, file=trim(sevt1), status="replace", action="write")
        
        
        
        write(xml_file, '(a)') '<?xml version="1.0" ?>'
        write(xml_file, '(a)') '<Xdmf Version="3.0">'
        write(xml_file, '(a)') '  <Domain>'
        write(xml_file, '(a)') '    <Grid Name="'//trim(namevarvtk(2))//'" GridType="Collection" CollectionType="Spatial">'

        do i = 0,nprocs-1

          if (.not. log_plane2d(i,myid))cycle
         
         
         sevt1=repeat(' ',mxln)
         sevt1 =  trim(filenamevtk)//'_'//trim(adjustl(space_fmtnumb(mydir)))// &
         '_'//trim(adjustl(space_fmtnumb(mypoint)))//'_'//trim(namevarvtk(2))// &
         '_'//trim(write_fmtnumb(iframe)) // '_p'//trim(write_fmtnumb(i)) //'.raw'
          
          write(xml_file, '(a)') '      <Grid Name="Piece' // trim(write_fmtnumb(i)) // '">' 
          !invert index since fortran is column-ordered
          write(strx, '(I0)') lsizes_plane2d(3,i,myid)
          write(stry, '(I0)') lsizes_plane2d(2,i,myid)
          write(strz, '(I0)') lsizes_plane2d(1,i,myid)
          write(xml_file, '(a)') '        <Topology TopologyType="3DCoRectMesh" Dimensions="' // &
           trim(adjustl(strx)) // ' ' // trim(adjustl(stry)) // ' ' // trim(adjustl(strz)) // '"/>'
          write(xml_file, '(a)') '        <Geometry GeometryType="ORIGIN_DXDYDZ">'
          write(strx1, '(I0)') myoffset_plane2d(3,i,myid)
          write(stry1, '(I0)') myoffset_plane2d(2,i,myid)
          write(strz1, '(I0)') myoffset_plane2d(1,i,myid)
          write(xml_file, '(a)') '          <DataItem Dimensions="3" Format="XML">', &
           trim(adjustl(strx1)) // ' ' // trim(adjustl(stry1)) // ' ' // trim(adjustl(strz1)) // '</DataItem>'
          write(xml_file, '(a)') '          <DataItem Dimensions="3" Format="XML">1 1 1</DataItem>'
          write(xml_file, '(a)') '        </Geometry>'
          write(xml_file, '(a)') '        <Attribute Name="'//trim(namevarvtk(2))//'" AttributeType="Vector" Center="Node">'
          write(xml_file, '(a)') '          <DataItem Format="Binary" NumberType="Float" Precision="4" Dimensions="' // &
           trim(adjustl(strx)) // ' ' // trim(adjustl(stry)) // ' ' // trim(adjustl(strz)) // ' 3">'
          write(xml_file, '(2a)') '            ',trim(sevt1)
          write(xml_file, '(a)') '          </DataItem>'
          write(xml_file, '(a)') '        </Attribute>'
          write(xml_file, '(a)') '      </Grid>'
    end do

    write(xml_file, '(a)') '    </Grid>'
    write(xml_file, '(a)') '  </Domain>'
    write(xml_file, '(a)') '</Xdmf>'
    close(xml_file)        
        
        
    endif
#endif
    if(ldowrite)then
      sevt2=repeat(' ',mxln)
      sevt2 = trim(dir_out) // trim(filenamevtk)//'_'//trim(adjustl(space_fmtnumb(mydir)))// &
         '_'//trim(adjustl(space_fmtnumb(mypoint)))//'_'//trim(namevarvtk(2))// &
         '_'//trim(write_fmtnumb(iframe)) // '_p'//trim(write_fmtnumb(myrank)) //'.raw'
      ! Each process writes its own RAW file if ldowrite is true.
      open(unit=82, file=trim(sevt2), form='unformatted', access='stream', status='replace', iostat=e_io)
      write(82) service3
      close(82)
    endif

#if defined(TWOCOMPONENT) && defined(WRITEPRESS)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!pressure!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
         
#ifdef DOXDMF
      if(myrank==0)then
        
        
                
        sevt2=repeat(' ',mxln)
        sevt2 =  trim(dir_out) // trim(filenamevtk)//'_'//trim(adjustl(space_fmtnumb(mydir)))// &
         '_'//trim(adjustl(space_fmtnumb(mypoint)))//'_'//trim(namevarvtk(3))// &
         '_'//trim(write_fmtnumb(iframe)) // '.xdmf'
        open(unit=xml_file, file=trim(sevt2), status="replace", action="write")
        
        
        write(xml_file, '(a)') '<?xml version="1.0" ?>'
        write(xml_file, '(a)') '<Xdmf Version="3.0">'
        write(xml_file, '(a)') '  <Domain>'
        write(xml_file, '(a)') '    <Grid Name="'//trim(namevarvtk(3))//'" GridType="Collection" CollectionType="Spatial">'

        do i = 0,nprocs-1

          if (.not. log_plane2d(i,myid))cycle
          
         
         sevt2=repeat(' ',mxln)
         sevt2 =  trim(filenamevtk)//'_'//trim(adjustl(space_fmtnumb(mydir)))// &
          '_'//trim(adjustl(space_fmtnumb(mypoint)))//'_'//trim(namevarvtk(3))// &
          '_'//trim(write_fmtnumb(iframe)) // '_p'//trim(write_fmtnumb(i)) //'.raw'
          
          write(xml_file, '(a)') '      <Grid Name="Piece' // trim(write_fmtnumb(i)) // '">' 
          !invert index since fortran is column-ordered
          write(strx, '(I0)') lsizes_plane2d(3,i,myid)
          write(stry, '(I0)') lsizes_plane2d(2,i,myid)
          write(strz, '(I0)') lsizes_plane2d(1,i,myid)
          write(xml_file, '(a)') '        <Topology TopologyType="3DCoRectMesh" Dimensions="' // &
           trim(adjustl(strx)) // ' ' // trim(adjustl(stry)) // ' ' // trim(adjustl(strz)) // '"/>'
          write(xml_file, '(a)') '        <Geometry GeometryType="ORIGIN_DXDYDZ">'
          write(strx1, '(I0)') myoffset_plane2d(3,i,myid)
          write(stry1, '(I0)') myoffset_plane2d(2,i,myid)
          write(strz1, '(I0)') myoffset_plane2d(1,i,myid)
          write(xml_file, '(a)') '          <DataItem Dimensions="3" Format="XML">', &
           trim(adjustl(strx1)) // ' ' // trim(adjustl(stry1)) // ' ' // trim(adjustl(strz1)) // '</DataItem>'
          write(xml_file, '(a)') '          <DataItem Dimensions="3" Format="XML">1 1 1</DataItem>'
          write(xml_file, '(a)') '        </Geometry>'
          write(xml_file, '(a)') '        <Attribute Name="'//trim(namevarvtk(3))//'" AttributeType="Scalar" Center="Node">'
          write(xml_file, '(a)') '          <DataItem Format="Binary" NumberType="Float" Precision="4" Dimensions="' // &
           trim(adjustl(strx)) // ' ' // trim(adjustl(stry)) // ' ' // trim(adjustl(strz)) // '">'
          write(xml_file, '(2a)') '            ',trim(sevt2)
          write(xml_file, '(a)') '          </DataItem>'
          write(xml_file, '(a)') '        </Attribute>'
          write(xml_file, '(a)') '      </Grid>'
    end do

    write(xml_file, '(a)') '    </Grid>'
    write(xml_file, '(a)') '  </Domain>'
    write(xml_file, '(a)') '</Xdmf>'
    close(xml_file)
        
    endif
#endif
    if(ldowrite)then
      sevt1=repeat(' ',mxln)
      sevt1 = trim(dir_out) // trim(filenamevtk)//'_'//trim(adjustl(space_fmtnumb(mydir)))// &
         '_'//trim(adjustl(space_fmtnumb(mypoint)))//'_'//trim(namevarvtk(3))// &
         '_'//trim(write_fmtnumb(iframe)) // '_p'//trim(write_fmtnumb(myrank)) //'.raw'
      ! Each process writes its own RAW file if ldowrite is true.
      open(unit=83, file=trim(sevt1), form='unformatted', access='stream', status='replace', iostat=e_io)
      write(83) service2
      close(83)
      
    endif
    
#endif
#endif
      return

   end subroutine write_file_raw_par2D_nompiio

   function GET_COORD_POINT(ii,jj,kk)

      implicit none

      integer, intent(in) :: ii,jj,kk

      integer :: i,j,k
      integer, dimension(mpid) :: GET_COORD_POINT

      do i=0,proc_x-1
         if(ii<=xfindom(i))then
            GET_COORD_POINT(1)=i
            exit
         endif
      enddo

      do j=0,proc_y-1
         if(jj<=yfindom(j))then
            GET_COORD_POINT(2)=j
            exit
         endif
      enddo

      do k=0,proc_z-1
         if(kk<=zfindom(k))then
            GET_COORD_POINT(3)=k
            exit
         endif
      enddo


   end function GET_COORD_POINT

   function GET_RANK_POINT(ii,jj,kk,ierr)

      implicit none

      integer, intent(in) :: ii,jj,kk
      integer, dimension(mpid) :: temp
      integer :: ierr

      integer :: GET_RANK_POINT

      temp=GET_COORD_POINT(ii,jj,kk)
#ifdef MPI
      call MPI_Cart_rank(lbecomm, temp, GET_RANK_POINT,ierr)
#else
      GET_RANK_POINT=0
#endif

   end function GET_RANK_POINT

   subroutine dostop(mystring,mystring2,inumber)

      implicit none

      integer :: ierr

      character(len=*), optional :: mystring
      character(len=*), optional :: mystring2
      integer, optional :: inumber
      
      if(myrank==0)then
        if(present(mystring))then
          if(present(inumber))then
            write(6,'(4a,i0)')trim(mystring),' - ',trim(mystring2),':',inumber
            call flush(6)
          else
            write(6,'(a)')mystring
            call flush(6)
          endif
        endif
      endif

#ifdef MPI
      !$acc wait
      call MPI_Barrier(MPI_COMM_WORLD,ierr)
      call MPI_finalize(ierr)
#endif
      stop

   end subroutine dostop

   subroutine doerror(errcode,mystring)

      implicit none

      integer, value :: errcode
      integer :: ierr

      character(len=*), optional :: mystring

      if(present(mystring))then
         write(6,'(a)')mystring
         call flush(6)
      endif

#ifdef MPI
      call MPI_Abort(MPI_COMM_WORLD, errcode, ierr)
#endif
      stop

   end subroutine doerror
   
   subroutine or_world_l(argument)

!***********************************************************************
!
!     LBsoft global 'logical or' subroutine for a logical array
!     originally written in JETSPIN by M. Lauricella et al.
!
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification March 2015
!
!***********************************************************************

      implicit none



      logical, intent(inout) :: argument
      logical, dimension(1) :: buffersub,lbuffer

      integer ierr

#ifdef MPI
      buffersub(1)=argument

      call MPI_ALLREDUCE(buffersub,lbuffer,1,MPI_LOGICAL, &
         MPI_LOR,MPI_COMM_WORLD,ierr)

      argument=lbuffer(1)

      !$acc wait
      call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif

      return

   end subroutine or_world_l

   subroutine or_world_larr(narr,argument)

!***********************************************************************
!
!     LBsoft global 'logical or' subroutine for a logical array
!     originally written in JETSPIN by M. Lauricella et al.
!
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification March 2015
!
!***********************************************************************

      implicit none

      integer, intent(in) :: narr
      logical, intent(inout), dimension(narr) :: argument
      logical, dimension(narr) :: lbuffer

      integer ierr

#ifdef MPI

      call MPI_ALLREDUCE(argument,lbuffer,narr,MPI_LOGICAL, &
         MPI_LOR,MPI_COMM_WORLD,ierr)

      argument(1:narr)=lbuffer(1:narr)

      !$acc wait
      call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif

      return

   end subroutine or_world_larr
   
   subroutine sum_world_iarr(narr,argument)
 
!***********************************************************************
!     
!     LBsoft global summation subroutine for a integer array
!     originally written in JETSPIN by M. Lauricella et al.
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(in) :: narr
  integer, intent(inout), dimension(narr) :: argument
  integer, dimension(narr) :: ibuffer
  
  integer ier
  
#ifdef MPI
      
      call MPI_ALLREDUCE(argument,ibuffer,narr,MPI_INTEGER, &
       MPI_SUM,MPI_COMM_WORLD,ier)
    
      argument(1:narr)=ibuffer(1:narr)
  
      !$acc wait
      call MPI_Barrier(MPI_COMM_WORLD,ier)
#endif
  
  return
  
 end subroutine sum_world_iarr
 
 subroutine sum_world_farr(narr,argument)
 
!***********************************************************************
!     
!     LBsoft global summation subroutine for a floating point array
!     originally written in JETSPIN by M. Lauricella et al.
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(in) :: narr
  real(kind=db), intent(inout), dimension(narr) :: argument
  real(kind=db), dimension(narr) :: ibuffer
  
  integer ier
  
#ifdef MPI
      
      call MPI_ALLREDUCE(argument,ibuffer,narr,MYMPIREAL, &
       MPI_SUM,MPI_COMM_WORLD,ier)
    
      argument(1:narr)=ibuffer(1:narr)
  
      !$acc wait
      call MPI_Barrier(MPI_COMM_WORLD,ier)
#endif
  
  return
  
 end subroutine sum_world_farr
 
  subroutine sum_world_int(argument)
 
!***********************************************************************
!     
!     LBsoft global summation subroutine for a floating point scalar
!     originally written in JETSPIN by M. Lauricella et al.
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(inout) :: argument
  integer, dimension(1) :: ibuffer1,ibuffer2
  
  integer ier
  
#ifdef MPI
      ibuffer1(1)=argument
      
      call MPI_ALLREDUCE(ibuffer1,ibuffer2,1,MPI_INTEGER, &
       MPI_SUM,MPI_COMM_WORLD,ier)
    
      argument=ibuffer2(1)
  
      !$acc wait
      call MPI_Barrier(MPI_COMM_WORLD,ier)
#endif
  
  return
  
 end subroutine sum_world_int
 
 subroutine sum_world_float(argument)
 
!***********************************************************************
!     
!     LBsoft global summation subroutine for a floating point scalar
!     originally written in JETSPIN by M. Lauricella et al.
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  real(kind=db), intent(inout) :: argument
  real(kind=db), dimension(1) :: ibuffer1,ibuffer2
  
  integer ier
  
#ifdef MPI
      ibuffer1(1)=argument
      
      call MPI_ALLREDUCE(ibuffer1,ibuffer2,1,MYMPIREAL, &
       MPI_SUM,MPI_COMM_WORLD,ier)
    
      argument=ibuffer2(1)
  
      !$acc wait
      call MPI_Barrier(MPI_COMM_WORLD,ier)
#endif
  
  return
  
 end subroutine sum_world_float
 
 subroutine print_version_code
 
  implicit none
  
  if(myrank.ne.0)return
   
#ifdef MPI
      write(6,'(a)') 'MPI VERSION COMPILED'
#else
      write(6,'(a)') 'SERIAL VERSION COMPILED'
#endif
  
  return
  
 end subroutine print_version_code
 
end module mpi_template
