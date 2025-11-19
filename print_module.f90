#include "defines.h"
module prints

   use vars
   use mpi_template, only : myrank,nprocs,mpid,coords,myoffset,gsizes,dostop, &
    doerror,write_file_vtk_par,write_file_raw_par, &
    write_restart_parallel_2c,read_restart_parallel_2c,io_comm2d, &
    write_restart_parallel_1c,read_restart_parallel_1c,myoffset_plane2d, &
#ifdef MPI
    mpi_comm_world, &
#endif
    setup_io_comm2d,write_file_raw_par2D,or_world_larr,lnoparallel2d, &
    sum_world_iarr,log_plane2d,or_world_l,lsizes_plane2d,lsizes_3d, &
    myoffset_3d,write_file_raw_par2D_nompiio,write_file_raw_par_nompiio, &
    lsizes,skip_lsizes,skip_myoffset,skip_gsizes,read_isfluid_parallel, &
    write_file_raw_par_isfluid,read_init_parallel

   implicit none

contains

   subroutine header_vtk(nxs,nys,nzs,mystring500,namevar,extent,ncomps,iinisub,iend,myoffset, &
      new_myoffset,indent)

      implicit none

      integer, intent(in) :: nxs,nys,nzs
      character(len=8),intent(in) :: namevar
      character(len=120),intent(in) :: extent
      integer, intent(in) :: ncomps,iinisub,myoffset
      integer, intent(out) :: iend,new_myoffset
      integer, intent(inout) :: indent

      !namevar='density1'

      character(len=500), intent(out) :: mystring500
      ! End-character for binary-record finalize.
      character(1), parameter:: end_rec = char(10)
      character(1) :: string1
      character(len=*),parameter :: topology='ImageData'
      integer :: ioffset,nele,bytechar,byteint,byter2,byter4,byter8,iini

      iini=iinisub
      bytechar=kind(end_rec)
      byteint=kind(iini)
      byter2  = 2
      byter4  = 4
      byter8  = 8

      mystring500=repeat(' ',500)

      iend=iini

      iini=iend+1
      nele=22
      iend=iend+nele
      mystring500(iini:iend)='<?xml version="1.0"?>'//end_rec

      new_myoffset=myoffset
      new_myoffset = new_myoffset + nele * bytechar


      iini=iend+1
      nele=67
      iend=iend+nele
      if(lelittle)then
         mystring500(iini:iend) = '<VTKFile type="'//trim(topology)// &
            '" version="0.1" byte_order="LittleEndian">'//end_rec
      else
         mystring500(iini:iend) = '<VTKFile type="'//trim(topology)// &
            '" version="0.1" byte_order="BigEndian">   '//end_rec
      endif

      new_myoffset = new_myoffset + 67 * bytechar


      indent = indent + 2
      iini=iend+1
      nele=70
      iend=iend+nele
      mystring500(iini:iend) = repeat(' ',indent)//'<'//trim(topology)//' WholeExtent="'//&
         trim(extent)//'">'//end_rec


      new_myoffset = new_myoffset + 70 * bytechar


      indent = indent + 2
      iini=iend+1
      nele=63
      iend=iend+nele
      mystring500(iini:iend) = repeat(' ',indent)//'<Piece Extent="'//trim(extent)//'">'//end_rec

      new_myoffset = new_myoffset + 63 * bytechar


      ! initializing offset pointer
      ioffset = 0

      indent = indent + 2
      iini=iend+1
      nele=18
      iend=iend+nele
      mystring500(iini:iend)=repeat(' ',indent)//'<PointData>'//end_rec

      new_myoffset = new_myoffset + 18 * bytechar

      indent = indent + 2
      iini=iend+1
#ifdef PRINTHALF
      nele=121
      iend=iend+nele

      if(ncomps/=1 .and. ncomps/=3)then
         write(6,'(a)')'ERROR in header_vtk'
         stop
      endif
      write(string1,'(i1)')ncomps
      mystring500(iini:iend)=repeat(' ',indent)//'<DataArray type="UnsignedShort" Name="'// &
         namevar//'" NumberOfComponents="'//string1// '" '//&
         'format="appended" offset="'//space_fmtnumb12(ioffset)//'"/>'//end_rec

      new_myoffset = new_myoffset + 121 * bytechar
#else      
      nele=115
      iend=iend+nele

      if(ncomps/=1 .and. ncomps/=3)then
         write(6,'(a)')'ERROR in header_vtk'
         stop
      endif
      write(string1,'(i1)')ncomps
      mystring500(iini:iend)=repeat(' ',indent)//'<DataArray type="Float32" Name="'// &
         namevar//'" NumberOfComponents="'//string1// '" '//&
         'format="appended" offset="'//space_fmtnumb12(ioffset)//'"/>'//end_rec

      new_myoffset = new_myoffset + 115 * bytechar
#endif

      indent = indent - 2
      iini=iend+1
      nele=19
      iend=iend+nele
      mystring500(iini:iend)=repeat(' ',indent)//'</PointData>'//end_rec

      new_myoffset = new_myoffset + 19 * bytechar


      indent = indent - 2
      iini=iend+1
      nele=13
      iend=iend+nele
      mystring500(iini:iend)=repeat(' ',indent)//'</Piece>'//end_rec


      new_myoffset = new_myoffset + 13 * bytechar


      indent = indent - 2
      iini=iend+1
      nele=15
      iend=iend+nele
      mystring500(iini:iend)=repeat(' ',indent)//'</'//trim(topology)//'>'//end_rec

      new_myoffset = new_myoffset + 15 * bytechar


      iini=iend+1
      nele=32
      iend=iend+nele
      mystring500(iini:iend)=repeat(' ',indent)//'<AppendedData encoding="raw">'//end_rec

      new_myoffset = new_myoffset + 32 * bytechar

      iini=iend+1
      nele=1
      iend=iend+nele
      mystring500(iini:iend)='_'

      new_myoffset = new_myoffset + 1 * bytechar

      return

   end subroutine header_vtk

   subroutine footer_vtk(nxs,nys,nzs,mystring30,iinisub,iend,myoffset, &
      new_myoffset,indent)

      implicit none

      integer, intent(in) :: nxs,nys,nzs
      integer, intent(in) :: iinisub,myoffset
      integer, intent(out) :: iend,new_myoffset
      integer, intent(inout) :: indent


      character(len=30), intent(out) :: mystring30
      ! End-character for binary-record finalize.
      character(1), parameter:: end_rec = char(10)
      character(1) :: string1
      character(len=*),parameter :: topology='ImageData'
      integer :: ioffset,nele,bytechar,byteint,byter2,byter4,byter8,iini

      iini=iinisub
      bytechar=kind(end_rec)
      byteint=kind(iini)
      byter2  = 2
      byter4  = 4
      byter8  = 8

      mystring30=repeat(' ',30)

      iend=iini

      iini=iend+1
      nele=1
      iend=iend+nele
      mystring30(iini:iend)=end_rec

      new_myoffset = myoffset
      new_myoffset = new_myoffset + 1 * bytechar



      iini=iend+1
      nele=18
      iend=iend+nele
      mystring30(iini:iend)=repeat(' ',indent)//'</AppendedData>'//end_rec

      new_myoffset = new_myoffset + 18 * bytechar

      iini=iend+1
      nele=11
      iend=iend+nele
      mystring30(iini:iend)='</VTKFile>'//end_rec

      if(iend/=30)then
         write(6,'(a)')'ERROR in footer_vtk'
         stop
      endif

      return

   end subroutine footer_vtk

   subroutine test_little_endian(ltest)

      !***********************************************************************
      !
      !     LBsoft subroutine for checking if the computing architecture
      !     is working in little-endian or big-endian
      !
      !     licensed under Open Software License v. 3.0 (OSL-3.0)
      !     author: M. Lauricella
      !     last modification October 2019
      !
      !***********************************************************************

      implicit none
      integer, parameter :: ik1 = selected_int_kind(2)
      integer, parameter :: ik4 = selected_int_kind(9)

      logical, intent(out) :: ltest

      if(btest(transfer(int((/1,0,0,0/),ik1),1_ik4),0)) then
         !it is little endian
         ltest=.true.
      else
         !it is big endian
         ltest=.false.
      end if

      return

   end subroutine test_little_endian

   subroutine init_output(ncomp,lvtk,lraw,nplanes,ndir,npoint)

      !***********************************************************************
      !
      !     LBsoft subroutine for creating the folders containing the files
      !     in image VTK legacy binary format in parallel IO
      !
      !     licensed under Open Software License v. 3.0 (OSL-3.0)
      !     author: M. Lauricella
      !     last modification October 2018
      !
      !***********************************************************************


      implicit none

      integer, intent(in) :: ncomp
      logical, intent(in) :: lvtk,lraw
      integer, intent(in) :: nplanes
      integer, dimension(nplanes), intent(in) :: ndir,npoint
      character(len=255) :: path,makedirectory
      logical :: lexist

      integer :: i,j,k,nn,indent,temp_offset,new_myoffset,iend
      integer, parameter :: byter2=2
      integer, parameter :: byter4=4
      integer, parameter :: byteint=4
      integer, dimension(3) :: printlistvtk
      integer, parameter :: ioxyz=54
      character(len=*), parameter :: filexyz='isfluid.xyz'
      character(len=120) :: mystring120

      if((.not. lvtk) .and. (.not. lraw) .and. (.not. lprint))return

      call test_little_endian(lelittle)

      sevt1=repeat(' ',mxln)
      sevt2=repeat(' ',mxln)

      path = repeat(' ',255)
      call getcwd(path)

      !call get_environment_variable('DELIMITER',delimiter)
      path = trim(path)
      delimiter = path(1:1)
      if (delimiter==' ') delimiter='/'



      makedirectory=repeat(' ',255)
      makedirectory = 'output'//delimiter
      dir_out=trim(makedirectory)
#ifdef _INTEL
      inquire(directory=trim(makedirectory),exist=lexist)
#else
      inquire(file=trim(makedirectory),exist=lexist)
#endif

      if(.not. lexist)then
         makedirectory=repeat(' ',255)
         makedirectory = 'mkdir output'
         if(lvtk .or. lraw)then
           if(myrank==0)call system(makedirectory)
         endif
      endif
      mystring120=repeat(' ',120)


      makedirectory=repeat(' ',255)
      makedirectory=trim(path)//delimiter//'output'//delimiter

      extentvtk =  space_fmtnumb(1) // ' ' // space_fmtnumb(lxskip) // ' ' &
         // space_fmtnumb(1) // ' ' // space_fmtnumb(lyskip) // ' ' &
         // space_fmtnumb(1) // ' ' // space_fmtnumb(lzskip)


#if defined(TWOCOMPONENT) && defined(WRITEPRESS)
       nfilevtk=3
#else
       nfilevtk=2
#endif

      do i=1,nfilevtk
         printlistvtk(i)=i
      enddo

      allocate(varlistvtk(nfilevtk))
      allocate(namevarvtk(nfilevtk))
      allocate(ndimvtk(nfilevtk))
      allocate(headervtk(nfilevtk))
      allocate(footervtk(nfilevtk))
      allocate(nheadervtk(nfilevtk))
      allocate(vtkoffset(nfilevtk))
      allocate(ndatavtk(nfilevtk))
      varlistvtk(1:nfilevtk)=printlistvtk(1:nfilevtk)

#if defined(TWOCOMPONENT) && defined(WRITEPRESS)
       do i=1,nfilevtk
          select case(printlistvtk(i))
           case(1)
             namevarvtk(i)='rho     '
             ndimvtk(i)=1
           case(2)
             namevarvtk(i)='vel     '
             ndimvtk(i)=3
           case(3)
             namevarvtk(i)='press   '
             ndimvtk(i)=1
           case default
             write(6,'(a)')'ERROR in init_output'
             stop
          end select
       enddo
#else
       do i=1,nfilevtk
          select case(printlistvtk(i))
           case(1)
             namevarvtk(i)='rho     '
             ndimvtk(i)=1
           case(2)
             namevarvtk(i)='vel     '
             ndimvtk(i)=3
           case default
             write(6,'(a)')'ERROR in init_output'
             stop
          end select
       enddo
#endif

      nn=lxskip*lyskip*lzskip

      do i=1,nfilevtk
         temp_offset=0
         indent=0
         call header_vtk(lxskip,lyskip,lzskip,headervtk(i),namevarvtk(i),extentvtk,ndimvtk(i),0,iend,temp_offset, &
            new_myoffset,indent)
         vtkoffset(i)=new_myoffset
#ifdef PRINTHALF
         temp_offset=new_myoffset+byteint+ndimvtk(i)*nn*byter2
         ndatavtk(i)=ndimvtk(i)*nn*byter2
#else
         temp_offset=new_myoffset+byteint+ndimvtk(i)*nn*byter4
         ndatavtk(i)=ndimvtk(i)*nn*byter4
#endif
         nheadervtk(i)=iend
         call footer_vtk(lxskip,lyskip,lzskip,footervtk(i),0,iend,temp_offset, &
            new_myoffset,indent)
      enddo
      
      allocate(lsizes_3d(mpid,0:nprocs-1))
      lsizes_3d=0
      
      lsizes_3d(1:mpid,myrank)=skip_lsizes
      iend=mpid*nprocs
      call sum_world_iarr(iend,lsizes_3d)
      
      allocate(myoffset_3d(mpid,0:nprocs-1))
      myoffset_3d=0
      myoffset_3d(1:mpid,myrank)=skip_myoffset
      
      call sum_world_iarr(iend,myoffset_3d)
      
      call  init_output_2D(nplanes,ndir,npoint,skip_npoint)

      return

   end subroutine init_output


   subroutine get_memory_gpu(fout,fout2)

      !***********************************************************************
      !
      !     LBsoft subroutine for register the memory usage
      !
      !     licensed under the 3-Clause BSD License (BSD-3-Clause)
      !     modified by: M. Lauricella
      !     last modification July 2018
      !
      !***********************************************************************
#ifdef _OPENACC
!      use openacc
#if defined(_OPENACC) && !defined(CRAY)
      use accel_lib
#endif
      use iso_c_binding
#elif defined(CUDA)
      use cudafor
#endif

      implicit none

      real(kind=db), intent(out) :: fout,fout2
      real(kind=db) :: myd(2),myd2(2)
      integer :: istat
#ifdef _OPENACC
      integer(c_size_t) :: myfree, total
#elif defined(CUDA)
      integer(kind=cuda_count_kind) :: myfree, total
#else
      integer :: myfree, total
#endif

#if defined(_OPENACC) && !defined(CRAY)
      myfree=acc_get_free_memory()
      total=acc_get_memory()
#elif defined(CUDA)
      istat = cudaMemGetInfo( myfree, total )
#else
      myfree=0
      total=0
#endif
      fout = real(total-myfree,kind=4)/(1024.0**3.0)
      fout2 = real(total,kind=4)/(1024.0**3.0)

      return

   end subroutine get_memory_gpu

   subroutine print_memory_registration_gpu(iu,mybanner,mybanner2,&
      mymemory,totmem)

      !***********************************************************************
      !
      !     LBcuda subroutine for printing the memory registration
      !
      !     licensed under the 3-Clause BSD License (BSD-3-Clause)
      !     author: M. Lauricella
      !     last modification April 2022
      !
      !***********************************************************************

      implicit none

      integer, intent(in) :: iu
      character(len=*), intent(in) :: mybanner,mybanner2
      real(kind=db), intent(in) :: mymemory,totmem

      character(len=12) :: r_char,r_char2

      character(len=*),parameter :: of='(a)'



      if(myrank/=0)return
#if defined(CRAY) && defined(_OPENACC)
      call system("rocm-smi --showmeminfo all -d 0")
      return
#endif      
      write (r_char,'(f12.4)')mymemory
      write (r_char2,'(f12.4)')totmem
      write(iu,of)"                                                                               "
      write(iu,of)"******************************GPU MEMORY MONITOR*******************************"
      write(iu,of)"                                                                               "
      write(iu,'(4a)')trim(mybanner)," = ",trim(adjustl(r_char))," (GB)"
      write(iu,'(4a)')trim(mybanner2)," = ",trim(adjustl(r_char2))," (GB)"
      write(iu,of)"                                                                               "
      write(iu,of)"*******************************************************************************"
      write(iu,of)"                                                                               "

      return

   end subroutine print_memory_registration_gpu
#if 0
   subroutine print_memory_registration(iu,mybanner,mybanner2,mymemory,totmem)

!***********************************************************************
!
!     LBsoft subroutine for printing the memory registration
!
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!
!***********************************************************************

      implicit none

      integer, intent(in) :: iu
      character(len=*), intent(in) :: mybanner,mybanner2
      real(kind=PRC), intent(in) :: mymemory,totmem

      character(len=12) :: r_char,r_char2

      character(len=*),parameter :: of='(a)'


      if(myrank/=0)return
      write (r_char,'(f12.4)')mymemory
      write (r_char2,'(f12.4)')totmem/real(1024.0,kind=db)
      write(iu,of)"                                                                               "
      write(iu,of)"********************************MEMORY MONITOR*********************************"
      write(iu,of)"                                                                               "
      write(iu,'(4a)')trim(mybanner)," = ",trim(adjustl(r_char))," (MB)"
      write(iu,'(4a)')trim(mybanner2)," = ",trim(adjustl(r_char2))," (GB)"
      write(iu,of)"                                                                               "
      write(iu,of)"*******************************************************************************"
      write(iu,of)"                                                                               "

      return

   end subroutine print_memory_registration
#endif
   !******************************************************************************************************!
#ifdef _OPENACC
   subroutine printDeviceProperties(ngpus,dev_Num,dev_Type,iu)


      use openacc
      use iso_c_binding

      integer :: ngpus,dev_Num
      integer(acc_device_kind) :: dev_Type

      integer,intent(in) :: iu
      integer(c_size_t) :: tot_mem
      character(len=255) :: myname,myvendor,mydriver
#if !defined(CRAY) && defined(_OPENACC) 
      call acc_get_property_string(dev_num,dev_Type,acc_property_name,myname)
      tot_mem = acc_get_property(dev_num,dev_Type,acc_property_memory)
      call acc_get_property_string(dev_num,dev_Type,acc_property_vendor,myvendor)
      call acc_get_property_string(dev_num,dev_Type,acc_property_driver,mydriver)
#else
      myname=repeat(' ',255)
      myvendor=repeat(' ',255)
      mydriver=repeat(' ',255)
      tot_mem=int(0,kind=c_size_t)
#endif
      if(myrank/=0)return
#if defined(CRAY) && defined(_OPENACC)      
      call system("rocm-smi --showproductname -d 0")
      return
#endif
      write(iu,907)"                                                                               "
      write(iu,907)"*****************************GPU FEATURE MONITOR*******************************"
      write(iu,907)"                                                                               "

      write (iu,900) "Device Number: "      ,ngpus,' per node'
      write (iu,901) "Device Name: "        ,trim(strip_null(myname))
      write (iu,903) "Total Global Memory: ",real(tot_mem)/(1024.0**3.0)," Gbytes"
      write (iu,901) "Vendor: "        ,trim(strip_null(myvendor))
      write (iu,901) "Driver: "        ,trim(strip_null(mydriver))

      write(iu,907)"                                                                               "
      write(iu,907)"*******************************************************************************"
      write(iu,907)"                                                                               "

900   format (a,i0,a)
901   format (a,a)
902   format (a,i0,a)
903   format (a,f16.8,a)
904   format (a,2(i0,1x,'x',1x),i0)
905   format (a,i0,'.',i0)
906   format (a,l0)
907   format (a)

      return

   end subroutine printDeviceProperties
#endif

   subroutine copy_print(iframe,hfields_s,phifields_s)

      implicit none

      integer, intent(in) :: iframe
      real(kind=db), allocatable, dimension(:) :: hfields_s,phifields_s
      
      integer :: ii,jj,kk
      integer :: iii,jjj,kkk
      integer :: xblock,yblock,zblock,myblock
      
      real(kind=db) :: rhophi_loc

#ifdef ACCNOKERNELS
#if defined(DENSRATIO) && defined(TWOCOMPONENT)
#ifdef WRITEPRESS
         !$acc parallel loop independent collapse(3) present(rhoprint,velprint,pressprint, &
         !$acc& rhophi,u,v,w,hfields_s,phifields_s,TILE_DIMx,TILE_DIMy,TILE_DIMz) &
         !$acc& private(i,j,k,ii,jj,kk,iii,jjj,kkk,xblock,yblock,zblock,myblock,rhophi_loc)
#else
         !$acc parallel loop independent collapse(3) present(rhoprint,velprint,rhophi, &
         !$acc& u,v,w,hfields_s,phifields_s,TILE_DIMx,TILE_DIMy,TILE_DIMz) &
         !$acc& private(i,j,k,ii,jj,kk,iii,jjj,kkk,xblock,yblock,zblock,myblock,rhophi_loc)
#endif
#endif

#if defined(TWOCOMPONENT) && !defined(DENSRATIO)
#ifdef WRITEPRESS
		 !$acc parallel loop independent collapse(3) present(rhoprint,velprint,pressprint, &
		 !$acc& selphi,rho,u,v,w,hfields_s,phifields_s,TILE_DIMx,TILE_DIMy,TILE_DIMz) &
		 !$acc& private(i,j,k,ii,jj,kk,iii,jjj,kkk,xblock,yblock,zblock,myblock,rhophi_loc)
#else
         !$acc parallel loop independent collapse(3) present(rhoprint,velprint,selphi, &
         !$acc& u,v,w,hfields_s,phifields_s,TILE_DIMx,TILE_DIMy,TILE_DIMz) &
         !$acc& private(i,j,k,ii,jj,kk,iii,jjj,kkk,xblock,yblock,zblock,myblock,rhophi_loc)
#endif 
#endif                 
#ifndef TWOCOMPONENT 
         !$acc parallel loop independent collapse(3) present(rhoprint,velprint, &
         !$acc& rho,u,v,w,hfields_s,phifields_s,TILE_DIMx,TILE_DIMy,TILE_DIMz) &
         !$acc& private(i,j,k,ii,jj,kk,iii,jjj,kkk,xblock,yblock,zblock,myblock,rhophi_loc)
#endif         
#else
#if defined(DENSRATIO) && defined(TWOCOMPONENT)
#ifdef WRITEPRESS
         !$acc kernels present(rhoprint,velprint,pressprint,rhophi,u,v,w, &
         !$acc& hfields_s,phifields_s,TILE_DIMx,TILE_DIMy,TILE_DIMz,rhophi_loc)
#else
         !$acc kernels present(rhoprint,velprint,rhophi,u,v,w, &
         !$acc& hfields_s,phifields_s,TILE_DIMx,TILE_DIMy,TILE_DIMz,rhophi_loc)
#endif
#endif
#if defined(TWOCOMPONENT) && !defined(DENSRATIO)
#ifdef WRITEPRESS
		 !$acc kernels present(rhoprint,velprint,pressprint,rho,selphi,u,v,w, &
		 !$acc& hfields_s,phifields_s,TILE_DIMx,TILE_DIMy,TILE_DIMz,rhophi_loc)
#else
         !$acc kernels present(rhoprint,velprint,selphi,u,v,w, &
         !$acc& hfields_s,phifields_s,TILE_DIMx,TILE_DIMy,TILE_DIMz,rhophi_loc)
#endif
#endif
#ifndef TWOCOMPONENT
         !$acc kernels present(rhoprint,velprint,rho,u,v,w,hfields_s, &
         !$acc& phifields_s,TILE_DIMx,TILE_DIMy,TILE_DIMz)
#endif
         !$acc loop independent collapse(3)  private(i,j,k,ii,jj,kk,iii,jjj,kkk, &
         !$acc& xblock,yblock,zblock,myblock,rhophi_loc)
#endif
         do k=1,nzskip
            do j=1,nyskip
               do i=1,nxskip
                  ii=i*stepskip
                  jj=j*stepskip
                  kk=k*stepskip

                  xblock=(ii+2*TILE_DIMx-1)/TILE_DIMx   
                  yblock=(jj+2*TILE_DIMy-1)/TILE_DIMy     
                  zblock=(kk+2*TILE_DIMz-1)/TILE_DIMz   
                  
                  myblock=(xblock-1)+(yblock-1)*nxblock+(zblock-1)*nxyblock+1
                  iii=ii-xblock*TILE_DIMx+2*TILE_DIMx
                  jjj=jj-yblock*TILE_DIMy+2*TILE_DIMy
                  kkk=kk-zblock*TILE_DIMz+2*TILE_DIMz                            
                  
#if defined(DENSRATIO) && defined(TWOCOMPONENT)
                  !rhoprint(i,j,k)=real(rhophi(i*stepskip,j*stepskip,k*stepskip),kind=printdb)
                  rhophi_loc=phifields_s(idx5(iii,jjj,kkk,1,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields))
                  rhophi_loc=rho_r*rhophi_loc+(ONE-rhophi_loc)*rho_b
                  rhoprint(i,j,k)=real(rhophi_loc,kind=printdb)
#ifdef WRITEPRESS 
                  !pressprint(i,j,k)=real(rho(i*stepskip,j*stepskip,k*stepskip),kind=printdb)                  
                  pressprint(i,j,k)=real(hfields_s(idx5(iii,jjj,kkk,1,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields)),kind=printdb)
#endif
#endif

#if defined(TWOCOMPONENT) && !defined(DENSRATIO)
                  !rhoprint(i,j,k)=real(selphi(i*stepskip,j*stepskip,k*stepskip,flip),kind=printdb)
                  rhoprint(i,j,k)=real(phifields_s(idx5(iii,jjj,kkk,1,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields)),kind=printdb)
#ifdef WRITEPRESS   
                  !pressprint(i,j,k)=real(rho(i*stepskip,j*stepskip,k*stepskip),kind=printdb)               
                  pressprint(i,j,k)=real(hfields_s(idx5(iii,jjj,kkk,1,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields)),kind=printdb)				  
#endif
#endif
#ifndef TWOCOMPONENT
				  !rhoprint(i,j,k)=real(rho(i*stepskip,j*stepskip,k*stepskip),kind=printdb)
				  rhoprint(i,j,k)=real(hfields_s(idx5(iii,jjj,kkk,1,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields)),kind=printdb)
#endif
				  velprint(1,i,j,k)=real(hfields_s(idx5(iii,jjj,kkk,2,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields)),kind=printdb)
				  velprint(2,i,j,k)=real(hfields_s(idx5(iii,jjj,kkk,3,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields)),kind=printdb)
				  velprint(3,i,j,k)=real(hfields_s(idx5(iii,jjj,kkk,4,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields)),kind=printdb)
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

   end subroutine copy_print

   subroutine driver_read_isfluid_raw(iframe)

      implicit none

      integer, intent(in) :: iframe
      integer :: e_io

      if(nprocs==1)then
         call read_isfluid_serial(iframe)
      else
         call read_isfluid_parallel(iframe,e_io)
      endif

   end subroutine driver_read_isfluid_raw
   
   subroutine driver_read_init_raw(iframe)

      implicit none

      integer, intent(in) :: iframe
      integer :: e_io
      real(kind=db) :: fneq1,feq, rhophi_loc

      if(nprocs==1)then
        if(flop==1)then
           call read_init_serial(iframe,hfields_flip,phifields_flip)
        else
           call read_init_serial(iframe,hfields_flop,phifields_flop)
        endif
      else
        if(flop==1)then
          call read_init_parallel(iframe,e_io,hfields_flip,phifields_flip)
        else
          call read_init_parallel(iframe,e_io,hfields_flop,phifields_flop)
        endif
      endif
      
      if(flop==1)then
        hfields_flop=hfields_flip
        phifields_flop=phifields_flip
      else
        hfields_flip=hfields_flop
        phifields_flip=phifields_flop
      endif
      
!!!!common to every LB
      do k=1,nz
         gk=nz*coords(3)+k
         do j=1,ny
            gj=ny*coords(2)+j
            do i=1,nx
               gi=nx*coords(1)+i
               if(abs(isfluid(i,j,k)).eq.1)then
				  uu=HALF*(u(i,j,k)*u(i,j,k)+v(i,j,k)*v(i,j,k)+w(i,j,k)*w(i,j,k))/cssq
				  do l=0,nlinks
					udotc=(u(i,j,k)*dex(l) + v(i,j,k)*dey(l)+ w(i,j,k)*dez(l))/cssq
					feq=p(l)*(rho(i,j,k) + udotc+0.5_db*udotc*udotc - uu)
                    f(i,j,k,l)=feq
				  enddo
               endif
            enddo
         enddo
      enddo      

   end subroutine driver_read_init_raw
   
   subroutine read_isfluid_serial(iframe)

      implicit none

      integer, intent(in) :: iframe
      logical :: file_exists

      sevt1 = 'isfluid.raw'
      
      inquire(file=trim(sevt1), exist=file_exists)
      if (.not. file_exists) then
         call doerror(6,'ERROR: unable to find the file: '//trim(sevt1))
      end if    
      
      open(unit=345,file=trim(sevt1), &
         status='old',action='read',access='stream',form='unformatted')

      read(345)lap_phi
      isfluid(1:nx,1:ny,1:nz)= int(lap_phi(1:nx,1:ny,1:nz),kind=isf) 

      close(345)

   end subroutine read_isfluid_serial
   
   subroutine read_init_serial(iframe,hfields_s,phifields_s)

      implicit none

      integer, intent(in) :: iframe
      real(kind=db), allocatable, dimension(:) :: hfields_s,phifields_s
      logical :: file_exists
      integer :: ii,jj,kk,myblock,xblock,yblock,zblock
      
      sevt1 = 'rho.raw'
      
      inquire(file=trim(sevt1), exist=file_exists)
      if (.not. file_exists) then
         call doerror(6,'ERROR: unable to find the file: '//trim(sevt1))
      end if    
      
      open(unit=345,file=trim(sevt1), &
         status='old',action='read',access='stream',form='unformatted')

      read(345)lap_phi
      rho(1:nx,1:ny,1:nz)= lap_phi(1:nx,1:ny,1:nz)
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
               hfields_s(idx5(ii,jj,kk,1,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields))=lap_phi(i,j,k)   
             enddo
         enddo
      enddo 

      close(345)
      
      sevt1 = 'u.raw'
      
      inquire(file=trim(sevt1), exist=file_exists)
      if (.not. file_exists) then
         call doerror(6,'ERROR: unable to find the file: '//trim(sevt1))
      end if    
      
      open(unit=345,file=trim(sevt1), &
         status='old',action='read',access='stream',form='unformatted')

      read(345)lap_phi
      u(1:nx,1:ny,1:nz)= lap_phi(1:nx,1:ny,1:nz)
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
               hfields_s(idx5(ii,jj,kk,2,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields))=lap_phi(i,j,k)   
             enddo
         enddo
      enddo 
      
      close(345)      
      
      sevt1 = 'v.raw'
      
      inquire(file=trim(sevt1), exist=file_exists)
      if (.not. file_exists) then
         call doerror(6,'ERROR: unable to find the file: '//trim(sevt1))
      end if    
      
      open(unit=345,file=trim(sevt1), &
         status='old',action='read',access='stream',form='unformatted')

      read(345)lap_phi
      v(1:nx,1:ny,1:nz)= lap_phi(1:nx,1:ny,1:nz)
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
               hfields_s(idx5(ii,jj,kk,3,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields))=lap_phi(i,j,k)   
             enddo
         enddo
      enddo 

      close(345)   
      
      sevt1 = 'w.raw'
      
      inquire(file=trim(sevt1), exist=file_exists)
      if (.not. file_exists) then
         call doerror(6,'ERROR: unable to find the file: '//trim(sevt1))
      end if    
      
      open(unit=345,file=trim(sevt1), &
         status='old',action='read',access='stream',form='unformatted')

      read(345)lap_phi
      w(1:nx,1:ny,1:nz)= lap_phi(1:nx,1:ny,1:nz)
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
               hfields_s(idx5(ii,jj,kk,4,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields))=lap_phi(i,j,k)   
             enddo
         enddo
      enddo 
      
      close(345)      

#ifdef TWOCOMPONENT
      sevt1 = 'phi.raw'
      
      inquire(file=trim(sevt1), exist=file_exists)
      if (.not. file_exists) then
         call doerror(6,'ERROR: unable to find the file: '//trim(sevt1))
      end if    
      
      open(unit=345,file=trim(sevt1), &
         status='old',action='read',access='stream',form='unformatted')

      read(345)lap_phi
      selphi(1:nx,1:ny,1:nz,flip)= lap_phi(1:nx,1:ny,1:nz)
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
               phifields_s(idx5(ii,jj,kk,1,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields))=lap_phi(i,j,k) 
             enddo
         enddo
      enddo 

      close(345)
      
      selphi(:,:,:,flop)=selphi(:,:,:,flip)
#ifdef DENSRATIO
      rhophi(:,:,:)=rho_r*selphi(:,:,:,flip)+(1.0_db-selphi(:,:,:,flip))*rho_b   
#endif         
#endif         

   end subroutine read_init_serial

   subroutine read_restart_1c(iframe,iframe2D,hfields_s)

      implicit none

      integer, intent(out) :: iframe,iframe2D
      real(kind=db), allocatable, dimension(:) :: hfields_s

      integer :: e_io,l,i,j,k
      
      real(kind=db) :: feq,fneq1

      if(nprocs==1)then
        call read_restart_serial_1c(iframe,iframe2D,hfields_s)
      else
        call read_restart_parallel_1c(iframe,iframe2D,e_io,hfields_s)
      endif
      !$acc update device(hfields_s)
      !$acc wait
      
   end subroutine read_restart_1c

   subroutine write_restart_1c(iframe,iframe2D,hfields_s)

      implicit none

      integer, intent(in) :: iframe,iframe2D
      real(kind=db), allocatable, dimension(:) :: hfields_s

      integer :: e_io
      
      !$acc update host(hfields_s)
      !$wait
      if(nprocs==1)then
        call write_restart_serial_1c(iframe,iframe2D,hfields_s)
      else
        call write_restart_parallel_1c(iframe,iframe2D,e_io,hfields_s)
      endif

   end subroutine write_restart_1c

   subroutine write_restart_serial_1c(iframe,iframe2D,hfields_s)

      implicit none

      integer, intent(in) :: iframe,iframe2D
      real(kind=db), allocatable, dimension(:) :: hfields_s
      integer :: ii,jj,kk,xblock,yblock,zblock,myblock

      sevt1 = trim(dir_out) // 'restart.raw'
      open(unit=345,file=trim(sevt1), &
         status='replace',action='write',access='stream',form='unformatted')

      lap_phi(1:nx,1:ny,1:nz)= rho(1:nx,1:ny,1:nz) 
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
               lap_phi(i,j,k)=hfields_s(idx5(ii,jj,kk,1,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields))
             enddo
         enddo
      enddo 
      write(345)lap_phi

      lap_phi(1:nx,1:ny,1:nz)= u(1:nx,1:ny,1:nz) 
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
               lap_phi(i,j,k)=hfields_s(idx5(ii,jj,kk,2,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields))
             enddo
         enddo
      enddo 
      write(345)lap_phi
      
      lap_phi(1:nx,1:ny,1:nz)= v(1:nx,1:ny,1:nz)
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
               lap_phi(i,j,k)=hfields_s(idx5(ii,jj,kk,3,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields))
             enddo
         enddo
      enddo  
      write(345)lap_phi
      
      
      lap_phi(1:nx,1:ny,1:nz)= w(1:nx,1:ny,1:nz) 
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
               lap_phi(i,j,k)=hfields_s(idx5(ii,jj,kk,4,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields))
             enddo
         enddo
      enddo 
      write(345)lap_phi

      lap_phi(1:nx,1:ny,1:nz)= pxx(1:nx,1:ny,1:nz) 
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
               lap_phi(i,j,k)=hfields_s(idx5(ii,jj,kk,5,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields))
             enddo
         enddo
      enddo 
      write(345)lap_phi
      
      lap_phi(1:nx,1:ny,1:nz)= pxy(1:nx,1:ny,1:nz) 
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
               lap_phi(i,j,k)=hfields_s(idx5(ii,jj,kk,8,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields))
             enddo
         enddo
      enddo 
      write(345)lap_phi
      
      lap_phi(1:nx,1:ny,1:nz)= pxz(1:nx,1:ny,1:nz) 
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
               lap_phi(i,j,k)=hfields_s(idx5(ii,jj,kk,9,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields))
             enddo
         enddo
      enddo 
      write(345)lap_phi
      
      lap_phi(1:nx,1:ny,1:nz)= pyy(1:nx,1:ny,1:nz) 
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
               lap_phi(i,j,k)=hfields_s(idx5(ii,jj,kk,6,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields))
             enddo
         enddo
      enddo 
      write(345)lap_phi
      
      lap_phi(1:nx,1:ny,1:nz)= pyz(1:nx,1:ny,1:nz) 
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
               lap_phi(i,j,k)=hfields_s(idx5(ii,jj,kk,10,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields))
             enddo
         enddo
      enddo 
      write(345)lap_phi
      
      lap_phi(1:nx,1:ny,1:nz)= pzz(1:nx,1:ny,1:nz)
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
               lap_phi(i,j,k)=hfields_s(idx5(ii,jj,kk,7,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields))
             enddo
         enddo
      enddo  
      write(345)lap_phi
      
      write(345)iframe,iframe2D

      close(345)

   end subroutine write_restart_serial_1c

   subroutine read_restart_serial_1c(iframe,iframe2D,hfields_s)

      implicit none

      integer, intent(out) :: iframe,iframe2D
      real(kind=db), allocatable, dimension(:) :: hfields_s
      logical :: file_exists
      integer :: ii,jj,kk,xblock,yblock,zblock,myblock

      sevt1 = trim(dir_out) // 'restart.raw'
      
      inquire(file=trim(sevt1), exist=file_exists)
      if (.not. file_exists) then
         call doerror(6,'ERROR: unable to find the file: '//trim(sevt1))
      end if    
      
      open(unit=345,file=trim(sevt1), &
         status='old',action='read',access='stream',form='unformatted')

      read(345)lap_phi
      rho(1:nx,1:ny,1:nz)= lap_phi(1:nx,1:ny,1:nz) 
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
               hfields_s(idx5(ii,jj,kk,1,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields))=lap_phi(i,j,k)   
             enddo
         enddo
      enddo       

      read(345)lap_phi
      u(1:nx,1:ny,1:nz)= lap_phi(1:nx,1:ny,1:nz) 
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
               hfields_s(idx5(ii,jj,kk,2,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields))=lap_phi(i,j,k)   
             enddo
         enddo
      enddo 
            
      read(345)lap_phi
      v(1:nx,1:ny,1:nz)= lap_phi(1:nx,1:ny,1:nz) 
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
               hfields_s(idx5(ii,jj,kk,3,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields))=lap_phi(i,j,k)   
             enddo
         enddo
      enddo 
      
      read(345)lap_phi
      w(1:nx,1:ny,1:nz)= lap_phi(1:nx,1:ny,1:nz) 
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
               hfields_s(idx5(ii,jj,kk,4,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields))=lap_phi(i,j,k)   
             enddo
         enddo
      enddo       

      read(345)lap_phi
      pxx(1:nx,1:ny,1:nz)= lap_phi(1:nx,1:ny,1:nz) 
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
               hfields_s(idx5(ii,jj,kk,5,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields))=lap_phi(i,j,k)   
             enddo
         enddo
      enddo 
            
      read(345)lap_phi
      pxy(1:nx,1:ny,1:nz)= lap_phi(1:nx,1:ny,1:nz) 
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
               hfields_s(idx5(ii,jj,kk,8,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields))=lap_phi(i,j,k)   
             enddo
         enddo
      enddo 
      
      read(345)lap_phi
      pxz(1:nx,1:ny,1:nz)= lap_phi(1:nx,1:ny,1:nz) 
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
               hfields_s(idx5(ii,jj,kk,9,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields))=lap_phi(i,j,k)   
             enddo
         enddo
      enddo 
      
      read(345)lap_phi
      pyy(1:nx,1:ny,1:nz)= lap_phi(1:nx,1:ny,1:nz) 
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
               hfields_s(idx5(ii,jj,kk,6,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields))=lap_phi(i,j,k)   
             enddo
         enddo
      enddo 
      
      read(345)lap_phi
      pyz(1:nx,1:ny,1:nz)= lap_phi(1:nx,1:ny,1:nz) 
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
               hfields_s(idx5(ii,jj,kk,10,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields))=lap_phi(i,j,k)   
             enddo
         enddo
      enddo 
      
      read(345)lap_phi
      pzz(1:nx,1:ny,1:nz)= lap_phi(1:nx,1:ny,1:nz) 
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
               hfields_s(idx5(ii,jj,kk,7,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields))=lap_phi(i,j,k)   
             enddo
         enddo
      enddo 
      
      read(345)iframe,iframe2D

      close(345)

   end subroutine read_restart_serial_1c

   subroutine read_restart_2c(iframe,iframe2D,hfields_s,phifields_s)

      implicit none

      integer, intent(out) :: iframe,iframe2D
      real(kind=db), allocatable, dimension(:) :: hfields_s,phifields_s

      integer :: e_io,l,i,j,k
      
      real(kind=db) :: feq,fneq1

      if(nprocs==1)then
        call read_restart_serial_2c(iframe,iframe2D,hfields_s,phifields_s)
      else
        call read_restart_parallel_2c(iframe,iframe2D,e_io,hfields_s,phifields_s)
      endif
      
      !$acc update device(hfields_s,phifields_s)
      !$acc wait
      
   end subroutine read_restart_2c

   subroutine write_restart_2c(iframe,iframe2D,hfields_s,phifields_s)

      implicit none

      integer, intent(in) :: iframe,iframe2D
      real(kind=db), allocatable, dimension(:) :: hfields_s,phifields_s

      integer :: e_io
      
      !$acc update host(hfields_s,phifields_s)
      !$wait
      if(nprocs==1)then
        call write_restart_serial_2c(iframe,iframe2D,hfields_s,phifields_s)
      else
        call write_restart_parallel_2c(iframe,iframe2D,e_io,hfields_s,phifields_s)
      endif

   end subroutine write_restart_2c

   subroutine write_restart_serial_2c(iframe,iframe2D,hfields_s,phifields_s)

      implicit none

      integer, intent(in) :: iframe,iframe2D
      real(kind=db), allocatable, dimension(:) :: hfields_s,phifields_s
      integer :: ii,jj,kk,xblock,yblock,zblock,myblock

      sevt1 = trim(dir_out) // 'restart.raw'
      open(unit=345,file=trim(sevt1), &
         status='replace',action='write',access='stream',form='unformatted')

      lap_phi(1:nx,1:ny,1:nz)= rho(1:nx,1:ny,1:nz) 
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
               lap_phi(i,j,k)=hfields_s(idx5(ii,jj,kk,1,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields))
             enddo
         enddo
      enddo 
      write(345)lap_phi

      lap_phi(1:nx,1:ny,1:nz)= u(1:nx,1:ny,1:nz) 
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
               lap_phi(i,j,k)=hfields_s(idx5(ii,jj,kk,2,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields))
             enddo
         enddo
      enddo 
      write(345)lap_phi
      
      lap_phi(1:nx,1:ny,1:nz)= v(1:nx,1:ny,1:nz) 
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
               lap_phi(i,j,k)=hfields_s(idx5(ii,jj,kk,3,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields))
             enddo
         enddo
      enddo 
      write(345)lap_phi
      
      lap_phi(1:nx,1:ny,1:nz)= w(1:nx,1:ny,1:nz) 
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
               lap_phi(i,j,k)=hfields_s(idx5(ii,jj,kk,4,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields))
             enddo
         enddo
      enddo 
      write(345)lap_phi

      lap_phi(1:nx,1:ny,1:nz)= pxx(1:nx,1:ny,1:nz) 
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
               lap_phi(i,j,k)=hfields_s(idx5(ii,jj,kk,5,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields))
             enddo
         enddo
      enddo 
      write(345)lap_phi
      
      lap_phi(1:nx,1:ny,1:nz)= pxy(1:nx,1:ny,1:nz) 
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
               lap_phi(i,j,k)=hfields_s(idx5(ii,jj,kk,8,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields))
             enddo
         enddo
      enddo 
      write(345)lap_phi
      
      lap_phi(1:nx,1:ny,1:nz)= pxz(1:nx,1:ny,1:nz) 
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
               lap_phi(i,j,k)=hfields_s(idx5(ii,jj,kk,9,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields))
             enddo
         enddo
      enddo 
      write(345)lap_phi
      
      lap_phi(1:nx,1:ny,1:nz)= pyy(1:nx,1:ny,1:nz) 
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
               lap_phi(i,j,k)=hfields_s(idx5(ii,jj,kk,6,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields))
             enddo
         enddo
      enddo 
      write(345)lap_phi
      
      lap_phi(1:nx,1:ny,1:nz)= pyz(1:nx,1:ny,1:nz) 
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
               lap_phi(i,j,k)=hfields_s(idx5(ii,jj,kk,10,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields))
             enddo
         enddo
      enddo 
      write(345)lap_phi
      
      lap_phi(1:nx,1:ny,1:nz)= pzz(1:nx,1:ny,1:nz) 
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
               lap_phi(i,j,k)=hfields_s(idx5(ii,jj,kk,7,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields))
             enddo
         enddo
      enddo 
      write(345)lap_phi

      lap_phi(1:nx,1:ny,1:nz)= selphi(1:nx,1:ny,1:nz,flop)
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
               lap_phi(i,j,k)=phifields_s(idx5(ii,jj,kk,1,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields))
             enddo
         enddo
      enddo         
      write(345)lap_phi
      
      write(345)iframe,iframe2D

      close(345)

   end subroutine write_restart_serial_2c

   subroutine read_restart_serial_2c(iframe,iframe2D,hfields_s,phifields_s)

      implicit none

      integer(kind=4), intent(out) :: iframe,iframe2D
      real(kind=db), allocatable, dimension(:) :: hfields_s,phifields_s
      logical :: file_exists
      integer :: ii,jj,kk,myblock,xblock,yblock,zblock

      sevt1 = trim(dir_out) // 'restart.raw'
      
      inquire(file=trim(sevt1), exist=file_exists)
      if (.not. file_exists) then
         call doerror(6,'ERROR: unable to find the file: '//trim(sevt1))
      end if       
      
      write(6,*)'sto leggendo il restart'
      
      open(unit=345,file=trim(sevt1), &
         status='old',action='read',access='stream',form='unformatted')

      read(345)lap_phi
      rho(1:nx,1:ny,1:nz)= lap_phi(1:nx,1:ny,1:nz)
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
               hfields_s(idx5(ii,jj,kk,1,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields))=lap_phi(i,j,k)   
             enddo
         enddo
      enddo 
                
      read(345)lap_phi
      u(1:nx,1:ny,1:nz)= lap_phi(1:nx,1:ny,1:nz) 
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
               hfields_s(idx5(ii,jj,kk,2,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields))=lap_phi(i,j,k)   
             enddo
         enddo
      enddo       
      
      read(345)lap_phi
      v(1:nx,1:ny,1:nz)= lap_phi(1:nx,1:ny,1:nz) 
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
               hfields_s(idx5(ii,jj,kk,3,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields))=lap_phi(i,j,k)   
             enddo
         enddo
      enddo 
      
      read(345)lap_phi
      w(1:nx,1:ny,1:nz)= lap_phi(1:nx,1:ny,1:nz) 
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
               hfields_s(idx5(ii,jj,kk,4,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields))=lap_phi(i,j,k)   
             enddo
         enddo
      enddo 

      read(345)lap_phi
      pxx(1:nx,1:ny,1:nz)= lap_phi(1:nx,1:ny,1:nz) 
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
               hfields_s(idx5(ii,jj,kk,5,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields))=lap_phi(i,j,k)   
             enddo
         enddo
      enddo 
            
      read(345)lap_phi
      pxy(1:nx,1:ny,1:nz)= lap_phi(1:nx,1:ny,1:nz) 
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
               hfields_s(idx5(ii,jj,kk,8,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields))=lap_phi(i,j,k)   
             enddo
         enddo
      enddo 
      
      read(345)lap_phi
      pxz(1:nx,1:ny,1:nz)= lap_phi(1:nx,1:ny,1:nz) 
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
               hfields_s(idx5(ii,jj,kk,9,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields))=lap_phi(i,j,k)   
             enddo
         enddo
      enddo 
      
      read(345)lap_phi
      pyy(1:nx,1:ny,1:nz)= lap_phi(1:nx,1:ny,1:nz) 
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
               hfields_s(idx5(ii,jj,kk,6,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields))=lap_phi(i,j,k)   
             enddo
         enddo
      enddo 
      
      read(345)lap_phi
      pyz(1:nx,1:ny,1:nz)= lap_phi(1:nx,1:ny,1:nz) 
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
               hfields_s(idx5(ii,jj,kk,10,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields))=lap_phi(i,j,k)   
             enddo
         enddo
      enddo 
      
      read(345)lap_phi
      pzz(1:nx,1:ny,1:nz)= lap_phi(1:nx,1:ny,1:nz) 
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
               hfields_s(idx5(ii,jj,kk,7,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nhfields))=lap_phi(i,j,k)   
             enddo
         enddo
      enddo 

      read(345)lap_phi
      selphi(1:nx,1:ny,1:nz,flop)= lap_phi(1:nx,1:ny,1:nz) 
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
               phifields_s(idx5(ii,jj,kk,1,myblock,TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields))=lap_phi(i,j,k) 
             enddo
         enddo
      enddo       
      
      read(345)iframe,iframe2D

      close(345)

   end subroutine read_restart_serial_2c
   
   subroutine driver_print_raw_isfluid(iframe)

      implicit none

      integer, intent(in) :: iframe

      if(nprocs==1)then
         call print_raw_isfluid(iframe)
      else
         call print_parraw_isfluid(iframe)
      endif

   end subroutine driver_print_raw_isfluid
   
   subroutine print_raw_isfluid(iframe)

      implicit none

      integer, intent(in) :: iframe
      character(len=8) :: namevarvtk_sub='isfluid '
      logical :: lexit,lexist
#ifdef DOXDMF      
      integer, parameter :: xml_file=734
      character(len=32) :: strx,stry,strz,bitorder
      bitorder=repeat(' ',32)
      bitorder='Big'
      if(lelittle)then
        bitorder='Little'
      endif
#endif
      sevt1 = trim(dir_out) // trim(filenamevtk)//'_'//trim(namevarvtk_sub)// &
         '_'//trim(write_fmtnumb(iframe)) // '.raw'
#ifdef DOXDMF
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
#endif      

      lexit=.false.
      if (myrank == 0) then
        inquire(file=trim(sevt1), exist=lexist)
        if (lexist) then
          lexit=.true.
        endif
      endif
      call or_world_l(lexit)
      if(lexit)return
      
      do k=1,nzskip
        do j=1,nyskip
          do i=1,nxskip
            rhoprint(i,j,k)=real(isfluid(i*stepskip,j*stepskip,k*stepskip),kind=printdb)
          enddo
        enddo
      enddo
      
      open(unit=345,file=trim(sevt1), &
         status='replace',action='write',access='stream',form='unformatted')
      write(345)rhoprint
      close(345)

   end subroutine print_raw_isfluid

   subroutine driver_print_raw_sync(iframe)

      implicit none

      integer, intent(in) :: iframe

      if(nprocs==1)then
         call print_raw_sync(iframe)
      else
         call print_parraw_sync(iframe)
      endif

   end subroutine driver_print_raw_sync


   subroutine print_raw_sync(iframe)

      implicit none

      integer, intent(in) :: iframe
#ifdef DOXDMF      
      integer, parameter :: xml_file=734
      character(len=32) :: strx,stry,strz,bitorder
      bitorder=repeat(' ',32)
      bitorder='Big'
      if(lelittle)then
        bitorder='Little'
      endif
#endif
      sevt1 = trim(dir_out) // trim(filenamevtk)//'_'//trim(namevarvtk(1))// &
         '_'//trim(write_fmtnumb(iframe)) // '.raw'
#ifdef DOXDMF
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
      write(xml_file,'(7a)') '      <Topology TopologyType="3DCoRectMesh" Dimensions="' , &
       trim(adjustl(strx)) , ' ' , trim(adjustl(stry)) , ' ' , trim(adjustl(strz)) , '"/>'
      write(xml_file,'(a)') ''
      write(xml_file,'(a)') '      <!-- Geometry: origin and spacing (dx, dy, dz) -->'
      write(xml_file,'(a)') '      <Geometry GeometryType="ORIGIN_DXDYDZ">'
      write(xml_file,'(a)') '        <DataItem Dimensions="3" Format="XML" Precision="4">0.0 0.0 0.0</DataItem>'
      write(xml_file,'(a)') '        <DataItem Dimensions="3" Format="XML" Precision="4">1.0 1.0 1.0</DataItem>'
      write(xml_file,'(a)') '      </Geometry>'
      write(xml_file,'(a)') ''
      write(xml_file,'(a)') '      <!-- Scalar field associated with the grid -->'
      write(xml_file,'(3a)') '      <Attribute Name="',trim(namevarvtk(1)),'" AttributeType="Scalar" Center="Node">'
#ifdef PRINTHALF
      write(xml_file,'(9a)') '        <DataItem Format="Binary" NumberType="UInt" Precision="2" Dimensions="' , &
       trim(adjustl(strx)) , ' ' , trim(adjustl(stry)) , ' ' , trim(adjustl(strz)) , '" Endian="',trim(adjustl(bitorder)),'">'
#else
      write(xml_file,'(9a)') '        <DataItem Format="Binary" NumberType="Float" Precision="4" Dimensions="' , &
       trim(adjustl(strx)) , ' ' , trim(adjustl(stry)) , ' ' , trim(adjustl(strz)) , '" Endian="',trim(adjustl(bitorder)),'">'
#endif
      write(xml_file,'(2a)') '          ',trim(sevt2)
      write(xml_file,'(a)') '        </DataItem>'
      write(xml_file,'(a)') '      </Attribute>'
      write(xml_file,'(a)') ''
      write(xml_file,'(a)') '    </Grid>'
      write(xml_file,'(a)') '  </Domain>'
      write(xml_file,'(a)') '</Xdmf>'
      close(xml_file)     
#endif      
      open(unit=345,file=trim(sevt1), &
         status='replace',action='write',access='stream',form='unformatted')
      write(345)rhoprint
      close(345)
      
      sevt2 = trim(dir_out) // trim(filenamevtk)//'_'//trim(namevarvtk(2))// &
         '_'//trim(write_fmtnumb(iframe)) // '.raw'
#ifdef DOXDMF
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
      write(xml_file,'(7a)') '      <Topology TopologyType="3DCoRectMesh" Dimensions="' , &
       trim(adjustl(strx)) , ' ' , trim(adjustl(stry)) , ' ' , trim(adjustl(strz)) , '"/>'
      write(xml_file,'(a)') ''
      write(xml_file,'(a)') '      <!-- Geometry: origin and spacing (dx, dy, dz) -->'
      write(xml_file,'(a)') '      <Geometry GeometryType="ORIGIN_DXDYDZ">'
      write(xml_file,'(a)') '        <DataItem Dimensions="3" Format="XML" Precision="4">0.0 0.0 0.0</DataItem>'
      write(xml_file,'(a)') '        <DataItem Dimensions="3" Format="XML" Precision="4">1.0 1.0 1.0</DataItem>'
      write(xml_file,'(a)') '      </Geometry>'
      write(xml_file,'(a)') ''
      write(xml_file,'(a)') '      <!-- Vector field associated with the grid -->'
      write(xml_file,'(3a)') '      <Attribute Name="',trim(namevarvtk(2)),'" AttributeType="Vector" Center="Node">'
#ifdef PRINTHALF
      write(xml_file,'(9a)') '        <DataItem Format="Binary" NumberType="UInt" Precision="2" Dimensions="' , &
       trim(adjustl(strx)),' ',trim(adjustl(stry)),' ',trim(adjustl(strz)),' 3" Endian="',trim(adjustl(bitorder)),'">'
#else
      write(xml_file,'(9a)') '        <DataItem Format="Binary" NumberType="Float" Precision="4" Dimensions="' , &
       trim(adjustl(strx)),' ',trim(adjustl(stry)),' ',trim(adjustl(strz)),' 3" Endian="',trim(adjustl(bitorder)),'">'
#endif
      write(xml_file,'(2a)') '          ',trim(sevt1)
      write(xml_file,'(a)') '        </DataItem>'
      write(xml_file,'(a)') '      </Attribute>'
      write(xml_file,'(a)') ''
      write(xml_file,'(a)') '    </Grid>'
      write(xml_file,'(a)') '  </Domain>'
      write(xml_file,'(a)') '</Xdmf>'
      close(xml_file)
#endif

      open(unit=346,file=trim(sevt2), &
         status='replace',action='write',access='stream',form='unformatted')
      write(346) velprint
      close(346)
      
#if defined(TWOCOMPONENT) && defined(WRITEPRESS)

      sevt3 = trim(dir_out) // trim(filenamevtk)//'_'//trim(namevarvtk(3))// &
         '_'//trim(write_fmtnumb(iframe)) // '.raw'
#ifdef DOXDMF
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
      write(xml_file,'(7a)') '      <Topology TopologyType="3DCoRectMesh" Dimensions="' , &
       trim(adjustl(strx)) , ' ' , trim(adjustl(stry)) , ' ' , trim(adjustl(strz)) , '"/>'
      write(xml_file,'(a)') ''
      write(xml_file,'(a)') '      <!-- Geometry: origin and spacing (dx, dy, dz) -->'
      write(xml_file,'(a)') '      <Geometry GeometryType="ORIGIN_DXDYDZ">'
      write(xml_file,'(a)') '        <DataItem Dimensions="3" Format="XML" Precision="4">0.0 0.0 0.0</DataItem>'
      write(xml_file,'(a)') '        <DataItem Dimensions="3" Format="XML" Precision="4">1.0 1.0 1.0</DataItem>'
      write(xml_file,'(a)') '      </Geometry>'
      write(xml_file,'(a)') ''
      write(xml_file,'(a)') '      <!-- Scalar field associated with the grid -->'
      write(xml_file,'(3a)') '      <Attribute Name="',trim(namevarvtk(3)),'" AttributeType="Scalar" Center="Node">'
#ifdef PRINTHALF
      write(xml_file,'(9a)') '        <DataItem Format="Binary" NumberType="UInt" Precision="2" Dimensions="' , &
       trim(adjustl(strx)) , ' ' , trim(adjustl(stry)) , ' ' , trim(adjustl(strz)) , '" Endian="',trim(adjustl(bitorder)),'">'
#else
      write(xml_file,'(9a)') '        <DataItem Format="Binary" NumberType="Float" Precision="4" Dimensions="' , &
       trim(adjustl(strx)) , ' ' , trim(adjustl(stry)) , ' ' , trim(adjustl(strz)) , '" Endian="',trim(adjustl(bitorder)),'">'
#endif
      write(xml_file,'(2a)') '          ',trim(sevt2)
      write(xml_file,'(a)') '        </DataItem>'
      write(xml_file,'(a)') '      </Attribute>'
      write(xml_file,'(a)') ''
      write(xml_file,'(a)') '    </Grid>'
      write(xml_file,'(a)') '  </Domain>'
      write(xml_file,'(a)') '</Xdmf>'
      close(xml_file)     
#endif      
      open(unit=347,file=trim(sevt3), &
         status='replace',action='write',access='stream',form='unformatted')
      write(347)pressprint
      close(347)      
#endif      

   end subroutine print_raw_sync

   subroutine print_parraw_sync(iframe)

      implicit none

      integer, intent(in) :: iframe
      integer :: e_io

#ifdef AVOIDMPIIO
      call write_file_raw_par_nompiio(iframe,e_io)
#else
      call write_file_raw_par(iframe,e_io)
#endif
   end subroutine print_parraw_sync
   
   subroutine print_parraw_isfluid(iframe)

      implicit none

      integer, intent(in) :: iframe
      integer :: e_io

      call write_file_raw_par_isfluid(iframe,e_io)

   end subroutine print_parraw_isfluid
   
   subroutine init_output_2D(nplanes,ndir,npoint,skip_npoint)

      implicit none

      integer, intent(in) :: nplanes
      integer, dimension(nplanes), intent(in) :: ndir,npoint,skip_npoint
      integer :: subchords(3),gi,gj,gk,ierr,nele
      integer, dimension(mpid) :: myoffset_plane,lsizes_plane,gsizes_plane
      logical :: ldoserial,ldowrite,lnoparallel
      

      allocate(io_comm2d(nplanes))
      allocate(lnoparallel2d(nplanes))
      allocate(myoffset_plane2d(mpid,0:nprocs-1,nplanes))
      myoffset_plane2d=0 
      
      allocate(lsizes_plane2d(mpid,0:nprocs-1,nplanes))
      lsizes_plane2d=0
      
      allocate(log_plane2d(0:nprocs-1,nplanes))
      log_plane2d=.false.

      if(nprocs.ne.1)then
         !questo è il caso MPI rognoso
         !devo travare quali processi hanno in carico il piano ed escludere gli altri
         !devo far caricare solo dai processi giusti i dati sugli array service1 e service3
         do l=1,nplanes
            ldoserial=.false.
            ldowrite=.false.
            select case(ndir(l))
             case(1)
               !caso perpendicolare a x
               !devo trovare quali processi hanno in carico i nodi lungo il piano gi=npoint(l)
               gi=npoint(l)
               subchords(1)=(gi-1)/nx
               !sono buoni tutti i processi che hanno subchords(1)==coords(1)
               !ciclo su tutti gli altri così iniziano a lavorare al prossimo file
               if(subchords(1).ne.coords(1))goto 240
               !myoffset(1) è il mio offset lungo x del mio sottodominio MPI e mi ridà il valore di i nel sottodominio
               i=gi-myoffset(1)
               !setto variabili utili in particolare per MPI-IO
               !offset in coordinate globali di ogni processo MPI
               myoffset_plane = [ 0, coords(2)*nyskip, coords(3)*nzskip ]
               lsizes_plane=[ 1, nyskip, nzskip ]
               gsizes_plane=[ 1, lyskip, lzskip ]
               !se le locali e globali dimensioni lungo y e z sono le stesse
               !allora un solo processo ha in capo il piano quindi tanto vale fare la stampa seriale
               ldoserial=(ny==ly .and. nz==lz)
               ldowrite=.true.
             case(2)
               !caso perpendicolare a y
               gj=npoint(l)
               subchords(2)=(gj-1)/ny
               if(subchords(2).ne.coords(2))goto 240
               j=gj-myoffset(2)
               myoffset_plane = [ coords(1)*nxskip, 0, coords(3)*nzskip ]
               lsizes_plane=[ nxskip, 1, nzskip ]
               gsizes_plane=[ lxskip, 1, lzskip ]
               !se le locali e globali dimensioni lungo x e z sono le stesse
               !allora un solo processo ha in capo il piano quindi tanto vale fare la stampa seriale
               ldoserial=(nx==lx .and. nz==lz)
               ldowrite=.true.
             case(3)
               !caso perpendicolare a z
               gk=npoint(l)
               subchords(3)=(gk-1)/nz
               if(subchords(3).ne.coords(3))goto 240 
               k=gk-myoffset(3)
               myoffset_plane = [ coords(1)*nxskip, coords(2)*nyskip, 0 ]
               lsizes_plane=[ nxskip, nyskip, 1 ]
               gsizes_plane=[ lxskip, lyskip, 1 ]
               !se le locali e globali dimensioni lungo x e y sono le stesse
               !allora un solo processo ha in capo il piano quindi tanto vale fare la stampa seriale
               ldoserial=(nx==lx .and. ny==ly)
               ldowrite=.true.
             case default
               call dostop('wrong argument in ndir for 2d plane print: only 1 to 3 is permitted')
            end select
240         continue
            myoffset_plane2d(1:mpid,myrank,l) = myoffset_plane
            lsizes_plane2d(1:mpid,myrank,l) = lsizes_plane
            log_plane2d(myrank,l)=ldowrite
            lnoparallel=ldoserial
            call or_world_l(lnoparallel)
            lnoparallel2d(l)=lnoparallel
            if(.not. lnoparallel2d(l))then
                call setup_io_comm2d(ldowrite,l)
            endif
         enddo
      endif
      
      nele=mpid*nprocs*nplanes
      call sum_world_iarr(nele,myoffset_plane2d)
      call sum_world_iarr(nele,lsizes_plane2d)
      
      nele=nprocs*nplanes
      call or_world_larr(nele,log_plane2d)

   end subroutine init_output_2D   
   
   subroutine driver_print_raw_sync2D(iframe,nplanes,ndir,npoint)

      implicit none

      integer, intent(in) :: iframe,nplanes
      integer, dimension(nplanes), intent(in) :: ndir,npoint
      real(4), dimension(:,:,:), allocatable :: service1
#if defined(TWOCOMPONENT) && defined(WRITEPRESS)
      real(4), dimension(:,:,:), allocatable :: service2
#endif
      real(4), dimension(:,:,:,:), allocatable :: service3
      integer :: subchords(3),gi,gj,gk,ierr
      integer, dimension(mpid) :: myoffset_plane,lsizes_plane,gsizes_plane
      logical :: ldoserial,ldowrite,lnoparallel
      

      if(nprocs==1)then
         do l=1,nplanes
            select case(ndir(l))
             case(1)
               gsizes_plane=[ 1, lyskip, lzskip ]
               allocate(service1(1,1:nyskip,1:nzskip))
               do k=1,nzskip
                  do j=1,nyskip
                     service1(1,j,k)=rhoprint(skip_npoint(l),j,k)
                  enddo
               enddo
#if defined(TWOCOMPONENT) && defined(WRITEPRESS)
               allocate(service2(1,1:nyskip,1:nzskip))
               do k=1,nzskip
                  do j=1,nyskip
                     service2(1,j,k)=pressprint(skip_npoint(l),j,k)
                  enddo
               enddo
#endif
               allocate(service3(3,1,1:nyskip,1:nzskip))
               do k=1,nzskip
                  do j=1,nyskip
                     service3(1:3,1,j,k)=velprint(1:3,skip_npoint(l),j,k)
                  enddo
               enddo
             case(2)
               gsizes_plane=[ lxskip, 1, lzskip ]
               allocate(service1(1:nxskip,1,1:nzskip))
               do k=1,nzskip
                  do i=1,nxskip
                     service1(i,1,k)=rhoprint(i,skip_npoint(l),k)
                  enddo
               enddo
#if defined(TWOCOMPONENT) && defined(WRITEPRESS)
               allocate(service2(1:nxskip,1,1:nzskip))
               do k=1,nzskip
                  do i=1,nxskip
                     service2(i,1,k)=pressprint(i,skip_npoint(l),k)
                  enddo
               enddo
#endif
               allocate(service3(3,1:nxskip,1,1:nzskip))
               do k=1,nzskip
                  do i=1,nxskip
                     service3(1:3,i,1,k)=velprint(1:3,i,skip_npoint(l),k)
                  enddo
               enddo
             case(3)
               gsizes_plane=[ lxskip, lyskip, 1 ]
               allocate(service1(1:nxskip,1:nyskip,1))
               do j=1,nyskip
                  do i=1,nxskip
                     service1(i,j,1)=rhoprint(i,j,skip_npoint(l))
                  enddo
               enddo
#if defined(TWOCOMPONENT) && defined(WRITEPRESS)
               allocate(service2(1:nxskip,1:nyskip,1))
               do j=1,nyskip
                  do i=1,nxskip
                     service2(i,j,1)=pressprint(i,j,skip_npoint(l))
                  enddo
               enddo
#endif
               allocate(service3(3,1:nxskip,1:nyskip,1))
               do j=1,nyskip
                  do i=1,nxskip
                     service3(1:3,i,j,1)=velprint(1:3,i,j,skip_npoint(l))
                  enddo
               enddo
             case default
               call dostop('wrong argument in ndir for 2d plane print: only 1 to 3 is permitted')
            end select
#if defined(TWOCOMPONENT) && defined(WRITEPRESS)
            call print_raw_sync2D(iframe,ndir(l),npoint(l),gsizes_plane,service1,service3,service2)
            deallocate(service2)
#else   
            call print_raw_sync2D(iframe,ndir(l),npoint(l),gsizes_plane,service1,service3)
#endif
            deallocate(service1,service3)
         enddo
      else
         !questo è il caso MPI rognoso
         !devo travare quali processi hanno in carico il piano ed escludere gli altri
         !devo far caricare solo dai processi giusti i dati sugli array service1 e service3
         do l=1,nplanes
            ldoserial=.false.
            ldowrite=.false.
            select case(ndir(l))
             case(1)
               !caso perpendicolare a x
               !devo trovare quali processi hanno in carico i nodi lungo il piano gi=npoint(l)
               gi=npoint(l)
               subchords(1)=(gi-1)/nx
               !sono buoni tutti i processi che hanno subchords(1)==coords(1)
               !ciclo su tutti gli altri così iniziano a lavorare al prossimo file
               if(subchords(1).ne.coords(1))goto 240
               !myoffset(1) è il mio offset lungo x del mio sottodominio MPI e mi ridà il valore di i nel sottodominio
               i=skip_npoint(l)-skip_myoffset(1)
               !setto variabili utili in particolare per MPI-IO
               !offset in coordinate globali di ogni processo MPI
               myoffset_plane = [ 0, coords(2)*nyskip, coords(3)*nzskip ]
               lsizes_plane=[ 1, nyskip, nzskip ]
               gsizes_plane=[ 1, lyskip, lzskip ]
               allocate(service1(1,1:nyskip,1:nzskip))
               do k=1,nzskip
                  do j=1,nyskip
                     service1(1,j,k)=rhoprint(i,j,k)
                  enddo
               enddo
#if defined(TWOCOMPONENT) && defined(WRITEPRESS)
               allocate(service2(1,1:nyskip,1:nzskip))
               do k=1,nzskip
                  do j=1,nyskip
                     service2(1,j,k)=pressprint(i,j,k)
                  enddo
               enddo
#endif
               allocate(service3(3,1,1:nyskip,1:nzskip))
               do k=1,nzskip
                  do j=1,nyskip
                     service3(1:3,1,j,k)=velprint(1:3,i,j,k)
                  enddo
               enddo
               !se le locali e globali dimensioni lungo y e z sono le stesse
               !allora un solo processo ha in capo il piano quindi tanto vale fare la stampa seriale
               ldoserial=(ny==ly .and. nz==lz)
               ldowrite=.true.
             case(2)

               !caso perpendicolare a y
               gj=npoint(l)
               subchords(2)=(gj-1)/ny
               if(subchords(2).ne.coords(2))goto 240
               j=skip_npoint(l)-skip_myoffset(2)
               myoffset_plane = [ coords(1)*nxskip, 0, coords(3)*nzskip ]
               lsizes_plane=[ nxskip, 1, nzskip ]
               gsizes_plane=[ lxskip, 1, lzskip ]
               allocate(service1(1:nxskip,1,1:nzskip))
               do k=1,nzskip
                  do i=1,nxskip
                     service1(i,1,k)=rhoprint(i,j,k)
                  enddo
               enddo
#if defined(TWOCOMPONENT) && defined(WRITEPRESS)
               allocate(service2(1:nxskip,1,1:nzskip))
               do k=1,nzskip
                  do i=1,nxskip
                     service2(i,1,k)=pressprint(i,j,k)
                  enddo
               enddo
#endif
               allocate(service3(3,1:nxskip,1,1:nzskip))
               do k=1,nzskip
                  do i=1,nxskip
                     service3(1:3,i,1,k)=velprint(1:3,i,j,k)
                  enddo
               enddo
               !se le locali e globali dimensioni lungo x e z sono le stesse
               !allora un solo processo ha in capo il piano quindi tanto vale fare la stampa seriale
               ldoserial=(nx==lx .and. nz==lz)
               ldowrite=.true.
             case(3)
               !caso perpendicolare a z
               gk=npoint(l)
               subchords(3)=(gk-1)/nz
               if(subchords(3).ne.coords(3))goto 240 
               k=skip_npoint(l)-skip_myoffset(3)
               myoffset_plane = [ coords(1)*nxskip, coords(2)*nyskip, 0 ]
               lsizes_plane=[ nxskip, nyskip, 1 ]
               gsizes_plane=[ lxskip, lyskip, 1 ]
               allocate(service1(1:nxskip,1:nyskip,1))
               do j=1,nyskip
                  do i=1,nxskip
                     service1(i,j,1)=rhoprint(i,j,k)
                  enddo
               enddo
#if defined(TWOCOMPONENT) && defined(WRITEPRESS)
               allocate(service2(1:nxskip,1:nyskip,1))
               do j=1,nyskip
                  do i=1,nxskip
                     service2(i,j,1)=pressprint(i,j,k)
                  enddo
               enddo
#endif
               allocate(service3(3,1:nxskip,1:nyskip,1))
               do j=1,nyskip
                  do i=1,nxskip
                     service3(1:3,i,j,1)=velprint(1:3,i,j,k)
                  enddo
               enddo
               !se le locali e globali dimensioni lungo x e y sono le stesse
               !allora un solo processo ha in capo il piano quindi tanto vale fare la stampa seriale
               ldoserial=(nx==lx .and. ny==ly)
               ldowrite=.true.
             case default
               call dostop('wrong argument in ndir for 2d plane print: only 1 to 3 is permitted')
            end select
240         continue
            if(lnoparallel2d(l))then
              if(ldoserial)then
#if defined(TWOCOMPONENT) && defined(WRITEPRESS)
                call print_raw_sync2D(iframe,ndir(l),npoint(l),gsizes_plane,service1,service3,service2)
#else   
                call print_raw_sync2D(iframe,ndir(l),npoint(l),gsizes_plane,service1,service3)
#endif
              endif
            else
               !ldowrite è true solo nei processi buoni 
#if defined(TWOCOMPONENT) && defined(WRITEPRESS)
               call print_parraw_sync2D(ldowrite,iframe,l,ndir(l),npoint(l),service1,service3, &
                  myoffset_plane,lsizes_plane,gsizes_plane,service2)
#else
               call print_parraw_sync2D(ldowrite,iframe,l,ndir(l),npoint(l),service1,service3, &
                  myoffset_plane,lsizes_plane,gsizes_plane)
#endif
            endif
            if(ldowrite)then
              deallocate(service1,service3)
#if defined(TWOCOMPONENT) && defined(WRITEPRESS)
              deallocate(service2)
#endif
            endif
         enddo
      endif

   end subroutine driver_print_raw_sync2D

   subroutine print_raw_sync2D(iframe,mydir,mypoint,gsizes_plane,service1,service3,service2)

      implicit none

      integer, intent(in) :: iframe,mydir,mypoint
      integer, intent(in), dimension(mpid) :: gsizes_plane
      real(4), dimension(:,:,:), allocatable :: service1
      real(4), dimension(:,:,:,:), allocatable :: service3

      real(4), dimension(:,:,:), allocatable, optional :: service2

#ifdef DOXDMF      
      integer, parameter :: xml_file=734
      character(len=32) :: strx,stry,strz,bitorder
      bitorder=repeat(' ',32)
      bitorder='Big'
      if(lelittle)then
        bitorder='Little'
      endif
#endif

      sevt1 = trim(dir_out) // trim(filenamevtk)//'_'//trim(adjustl(space_fmtnumb(mydir)))// &
         '_'//trim(adjustl(space_fmtnumb(mypoint)))//'_'//trim(namevarvtk(1))// &
         '_'//trim(write_fmtnumb(iframe)) // '.raw'
      
#ifdef DOXDMF
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
      write(xml_file,'(7a)') '      <Topology TopologyType="3DCoRectMesh" Dimensions="' , &
       trim(adjustl(strx)) , ' ' , trim(adjustl(stry)) , ' ' , trim(adjustl(strz)) , '"/>'
      write(xml_file,'(a)') ''
      write(xml_file,'(a)') '      <!-- Geometry: origin and spacing (dx, dy, dz) -->'
      write(xml_file,'(a)') '      <Geometry GeometryType="ORIGIN_DXDYDZ">'
      write(xml_file,'(a)') '        <DataItem Dimensions="3" Format="XML" Precision="4">0.0 0.0 0.0</DataItem>'
      write(xml_file,'(a)') '        <DataItem Dimensions="3" Format="XML" Precision="4">1.0 1.0 1.0</DataItem>'
      write(xml_file,'(a)') '      </Geometry>'
      write(xml_file,'(a)') ''
      write(xml_file,'(a)') '      <!-- Scalar field associated with the grid -->'
      write(xml_file,'(3a)') '      <Attribute Name="',trim(namevarvtk(1)),'" AttributeType="Scalar" Center="Node">'
#ifdef PRINTHALF
      write(xml_file,'(7a)') '        <DataItem Format="Binary" NumberType="UInt" Precision="2" Dimensions="' , &
       trim(adjustl(strx)),' ',trim(adjustl(stry)),' ',trim(adjustl(strz)),'" Endian="',trim(adjustl(bitorder)),'">'
#else
      write(xml_file,'(7a)') '        <DataItem Format="Binary" NumberType="Float" Precision="4" Dimensions="' , &
       trim(adjustl(strx)),' ',trim(adjustl(stry)),' ',trim(adjustl(strz)),'" Endian="',trim(adjustl(bitorder)),'">'
#endif
      write(xml_file,'(2a)') '          ',trim(sevt2)
      write(xml_file,'(a)') '        </DataItem>'
      write(xml_file,'(a)') '      </Attribute>'
      write(xml_file,'(a)') ''
      write(xml_file,'(a)') '    </Grid>'
      write(xml_file,'(a)') '  </Domain>'
      write(xml_file,'(a)') '</Xdmf>'
      close(xml_file)     
#endif      
      open(unit=345,file=trim(sevt1), &
         status='replace',action='write',access='stream',form='unformatted')
      write(345)service1
      close(345)
      
      sevt2 = trim(dir_out) // trim(filenamevtk)//'_'//trim(adjustl(space_fmtnumb(mydir)))// &
         '_'//trim(adjustl(space_fmtnumb(mypoint)))//'_'//trim(namevarvtk(2))// &
         '_'//trim(write_fmtnumb(iframe)) // '.raw'
#ifdef DOXDMF
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
      write(xml_file,'(7a)') '      <Topology TopologyType="3DCoRectMesh" Dimensions="' , &
       trim(adjustl(strx)) , ' ' , trim(adjustl(stry)) , ' ' , trim(adjustl(strz)) , '"/>'
      write(xml_file,'(a)') ''
      write(xml_file,'(a)') '      <!-- Geometry: origin and spacing (dx, dy, dz) -->'
      write(xml_file,'(a)') '      <Geometry GeometryType="ORIGIN_DXDYDZ">'
      write(xml_file,'(a)') '        <DataItem Dimensions="3" Format="XML" Precision="4">0.0 0.0 0.0</DataItem>'
      write(xml_file,'(a)') '        <DataItem Dimensions="3" Format="XML" Precision="4">1.0 1.0 1.0</DataItem>'
      write(xml_file,'(a)') '      </Geometry>'
      write(xml_file,'(a)') ''
      write(xml_file,'(a)') '      <!-- Vector field associated with the grid -->'
      write(xml_file,'(3a)') '      <Attribute Name="',trim(namevarvtk(2)),'" AttributeType="Vector" Center="Node">'
#ifdef PRINTHALF
      write(xml_file,'(7a)') '        <DataItem Format="Binary" NumberType="UInt" Precision="2" Dimensions="' , &
       trim(adjustl(strx)),' ',trim(adjustl(stry)),' ',trim(adjustl(strz)),' 3" Endian="',trim(adjustl(bitorder)),'">'
#else
      write(xml_file,'(7a)') '        <DataItem Format="Binary" NumberType="Float" Precision="4" Dimensions="' , &
       trim(adjustl(strx)),' ',trim(adjustl(stry)),' ',trim(adjustl(strz)),' 3" Endian="',trim(adjustl(bitorder)),'">'
#endif
      write(xml_file,'(2a)') '          ',trim(sevt1)
      write(xml_file,'(a)') '        </DataItem>'
      write(xml_file,'(a)') '      </Attribute>'
      write(xml_file,'(a)') ''
      write(xml_file,'(a)') '    </Grid>'
      write(xml_file,'(a)') '  </Domain>'
      write(xml_file,'(a)') '</Xdmf>'
      close(xml_file)
#endif
      open(unit=346,file=trim(sevt2), &
         status='replace',action='write',access='stream',form='unformatted')
      write(346)service3
      close(346)
      
#if defined(TWOCOMPONENT) && defined(WRITEPRESS)
            sevt3 = trim(dir_out) // trim(filenamevtk)//'_'//trim(adjustl(space_fmtnumb(mydir)))// &
         '_'//trim(adjustl(space_fmtnumb(mypoint)))//'_'//trim(namevarvtk(3))// &
         '_'//trim(write_fmtnumb(iframe)) // '.raw'
      
#ifdef DOXDMF
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
      write(xml_file,'(7a)') '      <Topology TopologyType="3DCoRectMesh" Dimensions="' , &
       trim(adjustl(strx)) , ' ' , trim(adjustl(stry)) , ' ' , trim(adjustl(strz)) , '"/>'
      write(xml_file,'(a)') ''
      write(xml_file,'(a)') '      <!-- Geometry: origin and spacing (dx, dy, dz) -->'
      write(xml_file,'(a)') '      <Geometry GeometryType="ORIGIN_DXDYDZ">'
      write(xml_file,'(a)') '        <DataItem Dimensions="3" Format="XML" Precision="4">0.0 0.0 0.0</DataItem>'
      write(xml_file,'(a)') '        <DataItem Dimensions="3" Format="XML" Precision="4">1.0 1.0 1.0</DataItem>'
      write(xml_file,'(a)') '      </Geometry>'
      write(xml_file,'(a)') ''
      write(xml_file,'(a)') '      <!-- Scalar field associated with the grid -->'
      write(xml_file,'(3a)') '      <Attribute Name="',trim(namevarvtk(3)),'" AttributeType="Scalar" Center="Node">'
#ifdef PRINTHALF
      write(xml_file,'(7a)') '        <DataItem Format="Binary" NumberType="UInt" Precision="2" Dimensions="' , &
       trim(adjustl(strx)),' ',trim(adjustl(stry)),' ',trim(adjustl(strz)),'" Endian="',trim(adjustl(bitorder)),'">'
#else
      write(xml_file,'(7a)') '        <DataItem Format="Binary" NumberType="Float" Precision="4" Dimensions="' , &
       trim(adjustl(strx)),' ',trim(adjustl(stry)),' ',trim(adjustl(strz)),'" Endian="',trim(adjustl(bitorder)),'">'
#endif
      write(xml_file,'(2a)') '          ',trim(sevt2)
      write(xml_file,'(a)') '        </DataItem>'
      write(xml_file,'(a)') '      </Attribute>'
      write(xml_file,'(a)') ''
      write(xml_file,'(a)') '    </Grid>'
      write(xml_file,'(a)') '  </Domain>'
      write(xml_file,'(a)') '</Xdmf>'
      close(xml_file)     
#endif      
      open(unit=347,file=trim(sevt3), &
         status='replace',action='write',access='stream',form='unformatted')
      write(347)service2
      close(347)
#endif

   end subroutine print_raw_sync2D
   

   subroutine print_parraw_sync2D(ldowrite,iframe,myid,mydir,mypoint,service1,service3,&
      myoffset_plane,lsizes_plane,gsizes_plane,service2)

      implicit none
      
      logical, intent(in) :: ldowrite
      integer, intent(in) :: iframe,myid,mydir,mypoint
      real(4), dimension(:,:,:), allocatable :: service1
      real(4), dimension(:,:,:,:), allocatable :: service3
      integer, dimension(mpid), intent(in) :: myoffset_plane,lsizes_plane,gsizes_plane
      real(4), dimension(:,:,:), allocatable, optional :: service2
      integer :: e_io
      
#if defined(TWOCOMPONENT) && defined(WRITEPRESS)
#ifdef AVOIDMPIIO
      call write_file_raw_par2D_nompiio(ldowrite,iframe,myid,mydir,mypoint,service1,service3, &
         myoffset_plane,lsizes_plane,gsizes_plane,e_io,service2)
#else
      call write_file_raw_par2D(ldowrite,iframe,myid,mydir,mypoint,service1,service3, &
         myoffset_plane,lsizes_plane,gsizes_plane,e_io,service2)
#endif
#else
#ifdef AVOIDMPIIO
      call write_file_raw_par2D_nompiio(ldowrite,iframe,myid,mydir,mypoint,service1,service3, &
         myoffset_plane,lsizes_plane,gsizes_plane,e_io)
#else
      call write_file_raw_par2D(ldowrite,iframe,myid,mydir,mypoint,service1,service3, &
         myoffset_plane,lsizes_plane,gsizes_plane,e_io)
#endif
#endif

   end subroutine print_parraw_sync2D

   subroutine print_raw_slice_sync(iframe)

      implicit none

      integer, intent(in) :: iframe
      !rho
      sevt1 = trim(dir_out) // 'out'//'_'//'rhoxy'// &
         '_'//trim(write_fmtnumb(iframe)) // '.raw'
      open(unit=745,file=trim(sevt1), &
         status='replace',action='write',access='stream',form='unformatted')
      write(745) rho(:,:,nz/2)
      close(745)
      sevt1 = trim(dir_out) // 'out'//'_'//'rhoxz'// &
         '_'//trim(write_fmtnumb(iframe)) // '.raw'
      open(unit=746,file=trim(sevt1), &
         status='replace',action='write',access='stream',form='unformatted')
      write(746) rho(:,ny/2,:)
      close(746)
      sevt1 = trim(dir_out) // 'out'//'_'//'rhoyz'// &
         '_'//trim(write_fmtnumb(iframe)) // '.raw'
      open(unit=747,file=trim(sevt1), &
         status='replace',action='write',access='stream',form='unformatted')
      write(747) rho(nx/2,:,:)
      close(747)

      !u
      sevt1 = trim(dir_out) // 'out'//'_'//'uxy'// &
         '_'//trim(write_fmtnumb(iframe)) // '.raw'
      open(unit=845,file=trim(sevt1), &
         status='replace',action='write',access='stream',form='unformatted')
      write(845) u(:,:,nz/2)
      close(845)
      sevt1 = trim(dir_out) // 'out'//'_'//'uxz'// &
         '_'//trim(write_fmtnumb(iframe)) // '.raw'
      open(unit=846,file=trim(sevt1), &
         status='replace',action='write',access='stream',form='unformatted')
      write(846) u(:,ny/2,:)
      close(846)
      sevt1 = trim(dir_out) // 'out'//'_'//'uyz'// &
         '_'//trim(write_fmtnumb(iframe)) // '.raw'
      open(unit=847,file=trim(sevt1), &
         status='replace',action='write',access='stream',form='unformatted')
      write(847) u(nx/2,:,:)
      close(847)

      !v
      sevt1 = trim(dir_out) // 'out'//'_'//'vxy'// &
         '_'//trim(write_fmtnumb(iframe)) // '.raw'
      open(unit=848,file=trim(sevt1), &
         status='replace',action='write',access='stream',form='unformatted')
      write(848) v(:,:,nz/2)
      close(848)
      sevt1 = trim(dir_out) // 'out'//'_'//'vyz'// &
         '_'//trim(write_fmtnumb(iframe)) // '.raw'
      open(unit=849,file=trim(sevt1), &
         status='replace',action='write',access='stream',form='unformatted')
      write(849) v(nx/2,:,:)
      close(849)
      sevt1 = trim(dir_out) // 'out'//'_'//'vxz'// &
         '_'//trim(write_fmtnumb(iframe)) // '.raw'
      open(unit=850,file=trim(sevt1), &
         status='replace',action='write',access='stream',form='unformatted')
      write(850) v(:,ny/2,:)
      close(850)

      !w
      sevt1 = trim(dir_out) // 'out'//'_'//'wxz'// &
         '_'//trim(write_fmtnumb(iframe)) // '.raw'
      open(unit=851,file=trim(sevt1), &
         status='replace',action='write',access='stream',form='unformatted')
      write(851) w(:,ny/2,:)
      close(851)
      sevt1 = trim(dir_out) // 'out'//'_'//'wyz'// &
         '_'//trim(write_fmtnumb(iframe)) // '.raw'
      open(unit=852,file=trim(sevt1), &
         status='replace',action='write',access='stream',form='unformatted')
      write(852) w(nx/2,:,:)
      close(852)
      sevt1 = trim(dir_out) // 'out'//'_'//'wxy'// &
         '_'//trim(write_fmtnumb(iframe)) // '.raw'
      open(unit=853,file=trim(sevt1), &
         status='replace',action='write',access='stream',form='unformatted')
      write(853) w(:,:,nz/2)
      close(853)

   end subroutine print_raw_slice_sync

   subroutine print_raw_slice_2c_sync(iframe)

      implicit none

      integer, intent(in) :: iframe
      !phi
      sevt1 = trim(dir_out) // 'out'//'_'//'phixy'// &
         '_'//trim(write_fmtnumb(iframe)) // '.raw'
      open(unit=745,file=trim(sevt1), &
         status='replace',action='write',access='stream',form='unformatted')
      write(745) selphi(:,:,nz/2,flip)
      close(745)
      sevt1 = trim(dir_out) // 'out'//'_'//'phixz'// &
         '_'//trim(write_fmtnumb(iframe)) // '.raw'
      open(unit=746,file=trim(sevt1), &
         status='replace',action='write',access='stream',form='unformatted')
      write(746) selphi(:,ny/2,:,flip)
      close(746)
      sevt1 = trim(dir_out) // 'out'//'_'//'phiyz'// &
         '_'//trim(write_fmtnumb(iframe)) // '.raw'
      open(unit=747,file=trim(sevt1), &
         status='replace',action='write',access='stream',form='unformatted')
      write(747) selphi(nx/2,:,:,flip)
      close(747)

      !u
      sevt1 = trim(dir_out) // 'out'//'_'//'uxy'// &
         '_'//trim(write_fmtnumb(iframe)) // '.raw'
      open(unit=845,file=trim(sevt1), &
         status='replace',action='write',access='stream',form='unformatted')
      write(845) u(:,:,nz/2)
      close(845)
      sevt1 = trim(dir_out) // 'out'//'_'//'uxz'// &
         '_'//trim(write_fmtnumb(iframe)) // '.raw'
      open(unit=846,file=trim(sevt1), &
         status='replace',action='write',access='stream',form='unformatted')
      write(846) u(:,ny/2,:)
      close(846)
      sevt1 = trim(dir_out) // 'out'//'_'//'uyz'// &
         '_'//trim(write_fmtnumb(iframe)) // '.raw'
      open(unit=847,file=trim(sevt1), &
         status='replace',action='write',access='stream',form='unformatted')
      write(847) u(nx/2,:,:)
      close(847)

      !v
      sevt1 = trim(dir_out) // 'out'//'_'//'vxy'// &
         '_'//trim(write_fmtnumb(iframe)) // '.raw'
      open(unit=848,file=trim(sevt1), &
         status='replace',action='write',access='stream',form='unformatted')
      write(848) v(:,:,nz/2)
      close(848)
      sevt1 = trim(dir_out) // 'out'//'_'//'vyz'// &
         '_'//trim(write_fmtnumb(iframe)) // '.raw'
      open(unit=849,file=trim(sevt1), &
         status='replace',action='write',access='stream',form='unformatted')
      write(849) v(nx/2,:,:)
      close(849)
      sevt1 = trim(dir_out) // 'out'//'_'//'vxz'// &
         '_'//trim(write_fmtnumb(iframe)) // '.raw'
      open(unit=850,file=trim(sevt1), &
         status='replace',action='write',access='stream',form='unformatted')
      write(850) v(:,ny/2,:)
      close(850)

      !w
      sevt1 = trim(dir_out) // 'out'//'_'//'wxz'// &
         '_'//trim(write_fmtnumb(iframe)) // '.raw'
      open(unit=851,file=trim(sevt1), &
         status='replace',action='write',access='stream',form='unformatted')
      write(851) w(:,ny/2,:)
      close(851)
      sevt1 = trim(dir_out) // 'out'//'_'//'wyz'// &
         '_'//trim(write_fmtnumb(iframe)) // '.raw'
      open(unit=852,file=trim(sevt1), &
         status='replace',action='write',access='stream',form='unformatted')
      write(852) w(nx/2,:,:)
      close(852)
      sevt1 = trim(dir_out) // 'out'//'_'//'wxy'// &
         '_'//trim(write_fmtnumb(iframe)) // '.raw'
      open(unit=853,file=trim(sevt1), &
         status='replace',action='write',access='stream',form='unformatted')
      write(853) w(:,:,nz/2)
      close(853)

   end subroutine print_raw_slice_2c_sync

   subroutine driver_print_vtk_sync(iframe)

      implicit none

      integer, intent(in) :: iframe

      if(nprocs==1)then
         call print_vtk_sync(iframe)
      else
         call print_parvtk_sync(iframe)
      endif

   end subroutine driver_print_vtk_sync

   subroutine print_vtk_sync(iframe)
      implicit none

      integer, intent(in) :: iframe

      sevt1 = trim(dir_out) // trim(filenamevtk)//'_'//trim(namevarvtk(1))// &
         '_'//trim(write_fmtnumb(iframe)) // '.vti'
      sevt2 = trim(dir_out) // trim(filenamevtk)//'_'//trim(namevarvtk(2))// &
         '_'//trim(write_fmtnumb(iframe)) // '.vti'
      open(unit=345,file=trim(sevt1), &
         status='replace',action='write',access='stream',form='unformatted')
      write(345)head1,ndatavtk(1),rhoprint,footervtk(1)
      close(345)
      open(unit=346,file=trim(sevt2), &
         status='replace',action='write',access='stream',form='unformatted')
      write(346)head2,ndatavtk(2),velprint,footervtk(2)
      close(346)
#if defined(TWOCOMPONENT) && defined(WRITEPRESS)
      sevt3 = trim(dir_out) // trim(filenamevtk)//'_'//trim(namevarvtk(3))// &
         '_'//trim(write_fmtnumb(iframe)) // '.vti'
      open(unit=347,file=trim(sevt3), &
         status='replace',action='write',access='stream',form='unformatted')
      write(347)head3,ndatavtk(3),pressprint,footervtk(3)
      close(347)
#endif
   end subroutine print_vtk_sync

   subroutine print_parvtk_sync(iframe)
      implicit none

      integer, intent(in) :: iframe

      integer :: e_io

      call write_file_vtk_par(iframe,e_io)

   end subroutine print_parvtk_sync

   subroutine print_raw_async(iframe)

      implicit none

      integer, intent(in) :: iframe

      sevt1 = trim(dir_out) // trim(filenamevtk)//'_'//trim(namevarvtk(1))// &
         '_'//trim(write_fmtnumb(iframe)) // '.raw'
      sevt2 = trim(dir_out) // trim(filenamevtk)//'_'//trim(namevarvtk(2))// &
         '_'//trim(write_fmtnumb(iframe)) // '.raw'
      open(unit=345,file=trim(sevt1), &
         status='replace',action='write',access='stream',form='unformatted',&
         asynchronous='yes')
      write(345,asynchronous='yes')rhoprint

      open(unit=346,file=trim(sevt2), &
         status='replace',action='write',access='stream',form='unformatted',&
         asynchronous='yes')
      write(346,asynchronous='yes')velprint
      
#if defined(TWOCOMPONENT) && defined(WRITEPRESS)
      sevt3 = trim(dir_out) // trim(filenamevtk)//'_'//trim(namevarvtk(3))// &
         '_'//trim(write_fmtnumb(iframe)) // '.raw'
      open(unit=37,file=trim(sevt3), &
         status='replace',action='write',access='stream',form='unformatted',&
         asynchronous='yes')
      write(347,asynchronous='yes')pressprint
#endif
   end subroutine print_raw_async

   subroutine print_vtk_async(iframe)
      implicit none

      integer, intent(in) :: iframe

      sevt1 = trim(dir_out) // trim(filenamevtk)//'_'//trim(namevarvtk(1))// &
         '_'//trim(write_fmtnumb(iframe)) // '.vti'
      sevt2 = trim(dir_out) // trim(filenamevtk)//'_'//trim(namevarvtk(2))// &
         '_'//trim(write_fmtnumb(iframe)) // '.vti'

      open(unit=345,file=trim(sevt1), &
         status='replace',action='write',access='stream',form='unformatted',&
         asynchronous='yes')
      write(345,asynchronous='yes')head1,ndatavtk(1),rhoprint


      open(unit=780,file=trim(sevt2), &
         status='replace',action='write',access='stream',form='unformatted',&
         asynchronous='yes')
      write(780,asynchronous='yes')head2,ndatavtk(2),velprint

#if defined(TWOCOMPONENT) && defined(WRITEPRESS)
      sevt3 = trim(dir_out) // trim(filenamevtk)//'_'//trim(namevarvtk(3))// &
         '_'//trim(write_fmtnumb(iframe)) // '.vti'

      open(unit=347,file=trim(sevt3), &
         status='replace',action='write',access='stream',form='unformatted',&
         asynchronous='yes')
      write(347,asynchronous='yes')head3,ndatavtk(3),pressprint
#endif
   end subroutine print_vtk_async

   subroutine close_print_async(lvtk)

      implicit none
      logical, intent(in) :: lvtk

      wait(345)
      if(lvtk)write(345)footervtk(1)
      close(345)


      wait(780)
      if(lvtk)write(780)footervtk(2)
      close(780)
      
#if defined(TWOCOMPONENT) && defined(WRITEPRESS)
      wait(347)
      if(lvtk)write(347)footervtk(3)
      close(347)
#endif      

   end subroutine close_print_async

   subroutine copystring(oldstring,newstring,lenstring)

!***********************************************************************
!
!     LBsoft subroutine to copy one character string into another
!     originally written in JETSPIN by M. Lauricella et al.
!
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!
!***********************************************************************

      implicit none

      character(len=*), intent(in) :: oldstring
      character(len=*), intent(out) :: newstring
      integer, intent(in) :: lenstring

      integer :: i

      do i=1,lenstring
         newstring(i:i)=oldstring(i:i)
      enddo

      return

   end subroutine copystring

   function intstr(string,lenstring,laststring)

!***********************************************************************
!
!     LBsoft subroutine for extracting integers from a character
!     string
!     originally written in JETSPIN by M. Lauricella et al.
!
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!
!***********************************************************************

      implicit none

      character(len=*), intent(inout) :: string
      integer, intent(in) :: lenstring
      integer, intent(out) :: laststring

      integer :: intstr

      integer :: j,isn
      character*1, parameter, dimension(0:9) :: &
         n=(/'0','1','2','3','4','5','6','7','8','9'/)
      logical :: flag,lcount,final
      character*1 :: ksn
      character*1, dimension(lenstring) :: word

      do j=1,lenstring
         word(j)=string(j:j)
      enddo

      isn=1
      laststring=0
      ksn='+'
      intstr=0
      flag=.false.
      final=.false.
      lcount=.false.


      do while(laststring<lenstring.and.(.not.final))

         laststring=laststring+1
         flag=.false.

         do j=0,9

            if(n(j)==word(laststring))then

               intstr=10*intstr+j
               lcount=.true.
               flag=.true.

            endif

         enddo

         if(lcount.and.(.not.flag))final=.true.
         if(flag .and. ksn=='-')isn=-1
         ksn=word(laststring)

      enddo

      intstr=isn*intstr

      do j=laststring,lenstring
         word(j-laststring+1)=word(j)
      enddo
      do j=lenstring-laststring+2,lenstring
         word(j)=' '
      enddo

      do j=1,lenstring
         string(j:j)=word(j)
      enddo

      return

   end function intstr

   function dblstr(string,lenstring,laststring)

!***********************************************************************
!
!     LBsoft subroutine for extracting double precisions from a
!     character string
!     originally written in JETSPIN by M. Lauricella et al.
!
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!
!***********************************************************************

      implicit none

      character(len=*), intent(inout) :: string
      integer, intent(in) :: lenstring
      integer, intent(out) :: laststring

      real(kind=db) :: dblstr

      logical :: flag,ldot,start,final
      integer :: iexp,idum,i,j,fail
      real(kind=db) :: sn,sten,sone

      character*1, parameter, dimension(0:9) :: &
         n=(/'0','1','2','3','4','5','6','7','8','9'/)
      character*1, parameter :: dot='.'
      character*1, parameter :: d='d'
      character*1, parameter :: e='e'

      character*1 :: ksn
      character*1, dimension(lenstring) :: word
      character(len=lenstring) :: work

      do j=1,lenstring
         word(j)=string(j:j)
      enddo

      laststring=0
      sn= ONE
      ksn='+'
      sten= TEN
      sone= ONE

      dblstr = ZERO
      iexp=0
      idum=0
      start=.false.
      ldot=.false.
      final=.false.

      do while(laststring<lenstring .and. (.not.final))

         laststring=laststring+1
         flag=.false.

         do j=0,9

            if(n(j)==word(laststring))then

               dblstr=sten*dblstr+sone*real(j,kind=db)
               flag=.true.
               start=.true.
            endif

         enddo


         if(dot==word(laststring))then

            flag=.true.
            sten= ONE
            ldot=.true.
            start=.true.

         endif

         if(flag .and. ksn=='-') sn=- ONE
         if(ldot)sone= real(sone,kind=db)/ TEN
         ksn=word(laststring)
         if(ksn=="D")ksn="d"
         if(ksn=="E")ksn="e"

         if(start)then
            if(d==ksn .or. e==ksn)then
               do i=1,lenstring-laststring
                  work(i:i)=word(i+laststring)
               enddo
               iexp=intstr(work,lenstring-laststring,idum)
               final=.true.
            endif
            if(.not.flag)final=.true.
         endif
      enddo

      dblstr=sn*dblstr*( TEN ** iexp)
      laststring=laststring+idum

      do j=laststring,lenstring
         word(j-laststring+1)=word(j)
      enddo
      do j=lenstring-laststring+2,lenstring
         word(j)=' '
      enddo

      do j=1,lenstring
         string(j:j)=word(j)
      enddo

      return

   end function dblstr

   pure function strip_null(s) result(r)

      implicit none
      character(len=*), intent(in) :: s
      character(len=len(s)) :: r
      integer :: i
      i = index(s, char(0))
      if (i > 0) then
        r = s(:i-1)
      else
        r = trim(s)
      end if
end function strip_null

endmodule
