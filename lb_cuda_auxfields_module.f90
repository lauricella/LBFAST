#include "defines.h"
#if !defined(_OPENACC)  
#error "To use this module the macros _OPENACC should be defined."
#endif

! CUDA kernels for the phase-field auxiliary quantities of the two-component model (TWOCOMPONENT).
! All fields are stored tile by tile, field(ii,jj,kk,ifield,myblock), with tiles of
! TILE_DIMx x TILE_DIMy x TILE_DIMz nodes (8x8x8 with the default GPUTILEX/Y/Z). The sub-domain
! node i (1-TILE_DIMx <= i <= nx+TILE_DIMx) lies in the x tile
!   xblock = (i+2*TILE_DIMx-1)/TILE_DIMx = floor((i-1)/TILE_DIMx)+2
! at the local index ii = i-xblock*TILE_DIMx+2*TILE_DIMx = i-(xblock-2)*TILE_DIMx (1..TILE_DIMx);
! the 2*TILE_DIMx offset keeps the numerator positive for the halo nodes i <= 0. Tiles 1 and
! nxblock_d = nx/TILE_DIMx+2 are halo tiles, whose layers next to the sub-domain are filled by the
! MPI/periodic halo exchange. The same holds in y and z, and the linear tile index is
!   myblock = (xblock-1) + (yblock-1)*nxblock_d + (zblock-1)*nxyblock_d + 1.
! Field indices: phifields 1 = phi; auxfields 1-3 = unit normal n, 4 = |grad phi|,
! 5-7 = theta = phi(1-phi) n; locauxfields 1 = lap phi, 2 = div theta; forces 1-3 = CSF force.
! Order in a half step (integrator_module): compute_norm_interface kernels -> auxfields halo
! exchange -> bcs_mesoscopic_phifields -> compute_div_theta_n_kernel -> compute_csf_force_kernel.
module lb_cuda_auxfields

   use vars
   use iso_c_binding
   use cudafor
   use mpi_template, only: coords,dostop,doerror,mydev,myrank,nprocs,nbuff,nbuffbvec
   use lb_cuda_vars

   implicit none
   

contains
#ifdef TWOCOMPONENT
  ! compute_norm_interface_kernel: interface geometry on all the tiles of the sub-domain (grid
  ! nx/TILE_DIMx x ny/TILE_DIMy x nz/TILE_DIMz, one thread per node), used without ASYNCMPI.
  ! For each fluid node (|isfluid| = 1) it reads phi (phifields 1) on the node and its 26
  ! neighbours and writes auxfields 1-3 (n), 4 (|grad phi|), 5-7 (theta = phi(1-phi) n) and
  ! locauxfields 1 (lap phi); other nodes are not written. The kernels _int, _xplus, ... below
  ! perform the same computation split into interior and face tiles for ASYNCMPI.
  attributes(global) subroutine compute_norm_interface_kernel(step, &
   iprobe,jprobe,kprobe,flop,nx,ny,nz,coords,isfluid,rho_r,rho_b, &
   ntotphifields,ntotauxfields,ntotlocauxfields,phifields_s,auxfields_s,locauxfields_s)

      implicit none
      integer :: step,iprobe,jprobe,kprobe,flop,nx,ny,nz
      integer, dimension(3) :: coords
      integer(kind=isf), dimension(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff) :: isfluid
      real(kind=db) :: rho_r,rho_b
      integer :: ntotphifields,ntotauxfields,ntotlocauxfields
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields,nblocks_d) :: phifields_s
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nauxfields,nblocks_d) :: auxfields_s
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nlocauxfields,nblocks_d) :: locauxfields_s
      
      real(kind=db), shared :: myphi(0:TILE_DIMx+1,0:TILE_DIMy+1,0:TILE_DIMz+1)
      real(kind=db):: grad_fix,grad_fiy,grad_fiz,mod_grad,normal_denom
      
      integer :: i,j,k,gi,gj,gk,myblock,intblock
      integer :: ii,jj,kk
      integer :: li,lj,lk
      integer :: xblock,yblock,zblock
      integer, parameter :: nhx=TILE_DIMx+2, nhy=TILE_DIMy+2, nhz=TILE_DIMz+2
      integer, parameter :: nthreads=TILE_DIMx*TILE_DIMy*TILE_DIMz
      integer :: ithread,ih
      
      ! one thread per interior node; the threads first load the tile and its halo into
      ! shared memory, nodes ih=ithread, ithread+nthreads, ... in x-fastest order
      ! With 8x8x8 tiles the shared array myphi(0:9,0:9,0:9) holds 1000 nodes for 512 threads, so
      ! each thread loads one or two of them. (li,lj,lk) are shared-memory coordinates, 0 and
      ! TILE_DIM+1 being the one-node halo; the sub-domain node (i,j,k) of a halo point lies in a
      ! neighbouring tile, possibly a halo tile, hence the full tile and local-index computation.
      ithread=(threadIdx%x-1)+TILE_DIMx*(threadIdx%y-1)+TILE_DIMx*TILE_DIMy*(threadIdx%z-1)
      do ih=ithread,nhx*nhy*nhz-1,nthreads
        lk=ih/(nhx*nhy)
        lj=(ih-lk*nhx*nhy)/nhx
        li=ih-lk*nhx*nhy-lj*nhx
        i = (blockIdx%x-1) * TILE_DIMx + li
        j = (blockIdx%y-1) * TILE_DIMy + lj
        k = (blockIdx%z-1) * TILE_DIMz + lk
        xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
        yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
        zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
        myblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
        ii=i-xblock*TILE_DIMx+2*TILE_DIMx
        jj=j-yblock*TILE_DIMy+2*TILE_DIMy
        kk=k-zblock*TILE_DIMz+2*TILE_DIMz
        myphi(li,lj,lk)=real(phifields_s(ii,jj,kk,1,myblock),kind=db)
      enddo
      call syncthreads

      ! each thread now handles the tile node (threadIdx%x,threadIdx%y,threadIdx%z), i.e. the
      ! shared-memory point (li,lj,lk) with the same indices
      li=threadIdx%x
      lj=threadIdx%y
      lk=threadIdx%z
      i = (blockIdx%x-1) * TILE_DIMx + li
      j = (blockIdx%y-1) * TILE_DIMy + lj
      k = (blockIdx%z-1) * TILE_DIMz + lk

      if(abs(isfluid(i,j,k)) .ne. 1)return

      ! the node is not in a halo tile: xblock = blockIdx%x+1 (y, z alike) and ii=li, jj=lj, kk=lk
      !internal-node block of the tile
      myblock=blockIdx%x+blockIdx%y*nxblock_d+blockIdx%z*nxyblock_d+1
      ii=li
      jj=lj
      kk=lk

      ! isotropic D3Q27 gradient, grad phi = (1/cs^2) sum_i w_i phi(x+c_i) c_i (factor 3 = 1/cs^2);
      ! for each component only the links with c_ia = +-1 contribute: one face pair (p1d3q27 = 2/27),
      ! four edge pairs (p2d3q27 = 1/54) and four corner pairs (p3d3q27 = 1/216)
	  grad_fix=3.0_db*(p1d3q27*(myphi(li+1,lj,lk)-myphi(li-1,lj,lk)) + &
		 p2d3q27*( (myphi(li+1,lj+1,lk)-myphi(li-1,lj-1,lk))+ &
		 (myphi(li+1,lj-1,lk)-myphi(li-1,lj+1,lk))+ &
		 (myphi(li+1,lj,lk+1)-myphi(li-1,lj,lk-1))+ &
		 (myphi(li+1,lj,lk-1)-myphi(li-1,lj,lk+1)) )  + &
		 p3d3q27*((myphi(li+1,lj+1,lk+1)-myphi(li-1,lj-1,lk-1))+ &
		 (myphi(li+1,lj-1,lk-1)-myphi(li-1,lj+1,lk+1))+ &
		 (myphi(li+1,lj-1,lk+1)-myphi(li-1,lj+1,lk-1))+ &
		 (myphi(li+1,lj+1,lk-1)-myphi(li-1,lj-1,lk+1))))

	  grad_fiy=3.0_db*(p1d3q27*(myphi(li,lj+1,lk)-myphi(li,lj-1,lk)) + &
		 p2d3q27*((myphi(li+1,lj+1,lk)-myphi(li-1,lj-1,lk))+ &
		 (myphi(li-1,lj+1,lk)-myphi(li+1,lj-1,lk))+ &
		 (myphi(li,lj+1,lk+1)-myphi(li,lj-1,lk-1))+ &
		 (myphi(li,lj+1,lk-1)-myphi(li,lj-1,lk+1)) ) + &
		 p3d3q27*((myphi(li+1,lj+1,lk+1)-myphi(li-1,lj-1,lk-1))+ &
		 (myphi(li-1,lj+1,lk-1)-myphi(li+1,lj-1,lk+1))+ &
		 (myphi(li+1,lj+1,lk-1)-myphi(li-1,lj-1,lk+1))+ &
		 (myphi(li-1,lj+1,lk+1)-myphi(li+1,lj-1,lk-1))))

	  grad_fiz=3.0_db*(p1d3q27*(myphi(li,lj,lk+1)-myphi(li,lj,lk-1)) + &
		 p2d3q27*((myphi(li+1,lj,lk+1)-myphi(li-1,lj,lk-1))+ &
		 (myphi(li-1,lj,lk+1)-myphi(li+1,lj,lk-1))+ &
		 (myphi(li,lj+1,lk+1)-myphi(li,lj-1,lk-1))+ &
		 (myphi(li,lj-1,lk+1)-myphi(li,lj+1,lk-1)) ) + &
		 p3d3q27*((myphi(li+1,lj+1,lk+1)-myphi(li-1,lj-1,lk-1)) &
		 +(myphi(li-1,lj-1,lk+1)-myphi(li+1,lj+1,lk-1))+ &
		 (myphi(li+1,lj-1,lk+1)-myphi(li-1,lj+1,lk-1))+ &
		 (myphi(li-1,lj+1,lk+1)-myphi(li+1,lj-1,lk-1))))
      
	  mod_grad=sqrt(grad_fix**TWO+grad_fiy**TWO+grad_fiz**TWO)
      ! regularized normalization, finite also where grad phi = 0 (bulk phases): with CSF
      ! n = grad phi/sqrt(|grad phi|^2+eps^2), eps = 1e-12 (docs/CSF_LAPLACE_VALIDATION.md),
      ! otherwise n = grad phi/(|grad phi|+1e-9)
#ifdef CSF
	  normal_denom=sqrt(mod_grad*mod_grad+1.0e-24_db)
#else
	  normal_denom=mod_grad+1.0e-9_db
#endif

	  auxfields_s(ii,jj,kk,1,myblock)= &
	   real(grad_fix/normal_denom,kind=strdb)
	  auxfields_s(ii,jj,kk,2,myblock)= &
	   real(grad_fiy/normal_denom,kind=strdb)
	  auxfields_s(ii,jj,kk,3,myblock)= &
	   real(grad_fiz/normal_denom,kind=strdb)
	  
	  auxfields_s(ii,jj,kk,4,myblock)=mod_grad 

      ! theta = phi(1-phi) n; its divergence (locauxfields 2, compute_div_theta_n_kernel) enters the
      ! phi update and the moment kernels as -sharp_c*div theta
	  auxfields_s(ii,jj,kk,5,myblock)=real( &
	   myphi(li,lj,lk)*(1.0_db-myphi(li,lj,lk))*(grad_fix/normal_denom),kind=strdb)
	  auxfields_s(ii,jj,kk,6,myblock)=real( &
	   myphi(li,lj,lk)*(1.0_db-myphi(li,lj,lk))*(grad_fiy/normal_denom),kind=strdb)
	  auxfields_s(ii,jj,kk,7,myblock)=real( &
	   myphi(li,lj,lk)*(1.0_db-myphi(li,lj,lk))*(grad_fiz/normal_denom),kind=strdb)
	   
      ! isotropic D3Q27 Laplacian, lap phi = (2/cs^2) sum_{i/=0} w_i (phi(x+c_i)-phi(x)), written
      ! with sum_{i/=0} w_i = 1-w_0, hence the central term phi*(p0d3q27-1)
      !lap_phi here
      locauxfields_s(ii,jj,kk,1,myblock)=real( &
                   (2.0_db*invcssq)*(myphi(li,lj,lk)*(p0d3q27-1.0_db) + &
                   ( p1d3q27*(myphi(li+1,lj,lk)+myphi(li-1,lj,lk) + &
                   myphi(li,lj+1,lk)+myphi(li,lj-1,lk) + &
                   myphi(li,lj,lk+1)+myphi(li,lj,lk-1)) + &
                   p2d3q27*( (myphi(li+1,lj+1,lk)+myphi(li-1,lj-1,lk))+ &
                   (myphi(li+1,lj-1,lk)+myphi(li-1,lj+1,lk))+ &
                   (myphi(li+1,lj,lk+1)+myphi(li-1,lj,lk-1))+ &
                   (myphi(li+1,lj,lk-1)+myphi(li-1,lj,lk+1)) + &
                   (myphi(li,lj+1,lk+1)+myphi(li,lj-1,lk-1))+ &
                   (myphi(li,lj+1,lk-1)+myphi(li,lj-1,lk+1)) )  + &
                   p3d3q27*((myphi(li+1,lj+1,lk+1)+myphi(li-1,lj-1,lk-1))+ &
                   (myphi(li+1,lj-1,lk-1)+myphi(li-1,lj+1,lk+1))+ &
                   (myphi(li+1,lj-1,lk+1)+myphi(li-1,lj+1,lk-1))+ &
                   (myphi(li+1,lj+1,lk-1)+myphi(li-1,lj-1,lk+1))))),kind=strdb)
  
      return
      
   endsubroutine compute_norm_interface_kernel
   
   ! ASYNCMPI variants of compute_norm_interface_kernel: the body is identical, only the tiles
   ! covered differ (nxb = nx/TILE_DIMx = nxblock_d-2, nyb and nzb alike). The six face kernels,
   ! launched by compute_norm_interface_cuda_ext, cover the tiles touching the sub-domain faces
   ! without overlap: _zminus/_zplus all x and y tiles of the first/last z tile layer (edges and
   ! corners included), _yminus/_yplus all x tiles and z tiles 2..nzb-1 of the first/last y layer,
   ! _xminus/_xplus y tiles 2..nyb-1 and z tiles 2..nzb-1 of the first/last x layer. The
   ! integrator then starts the auxfields halo exchange and runs the _int kernel on the remaining
   ! tiles while the messages are in flight.
   ! compute_norm_interface_kernel_int: interior tiles 2..nxb-1, 2..nyb-1, 2..nzb-1
   ! (grid dimGridInt = (nxb-2,nyb-2,nzb-2), xblock = blockIdx%x+2, y and z alike).
     attributes(global) subroutine compute_norm_interface_kernel_int(step, &
      iprobe,jprobe,kprobe,flop,nx,ny,nz,coords,isfluid,rho_r,rho_b, &
      ntotphifields,ntotauxfields,ntotlocauxfields,phifields_s,auxfields_s,locauxfields_s)

      implicit none
      integer :: step,iprobe,jprobe,kprobe,flop,nx,ny,nz
      integer, dimension(3) :: coords
      integer(kind=isf), dimension(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff) :: isfluid
      real(kind=db) :: rho_r,rho_b
      integer :: ntotphifields,ntotauxfields,ntotlocauxfields
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields,nblocks_d) :: phifields_s
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nauxfields,nblocks_d) :: auxfields_s
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nlocauxfields,nblocks_d) :: locauxfields_s
      
      real(kind=db), shared :: myphi(0:TILE_DIMx+1,0:TILE_DIMy+1,0:TILE_DIMz+1)
      real(kind=db):: grad_fix,grad_fiy,grad_fiz,mod_grad,normal_denom
      
      integer :: i,j,k,myblock,intblock
      integer :: ii,jj,kk
      integer :: li,lj,lk
      integer :: xblock,yblock,zblock
      integer, parameter :: nhx=TILE_DIMx+2, nhy=TILE_DIMy+2, nhz=TILE_DIMz+2
      integer, parameter :: nthreads=TILE_DIMx*TILE_DIMy*TILE_DIMz
      integer :: ithread,ih
!      integer :: gi,gj,gk
      
      ! one thread per interior node; the threads first load the tile and its halo into
      ! shared memory, nodes ih=ithread, ithread+nthreads, ... in x-fastest order
      ithread=(threadIdx%x-1)+TILE_DIMx*(threadIdx%y-1)+TILE_DIMx*TILE_DIMy*(threadIdx%z-1)
      do ih=ithread,nhx*nhy*nhz-1,nthreads
        lk=ih/(nhx*nhy)
        lj=(ih-lk*nhx*nhy)/nhx
        li=ih-lk*nhx*nhy-lj*nhx
        i = (blockIdx%x) * TILE_DIMx + li
        j = (blockIdx%y) * TILE_DIMy + lj
        k = (blockIdx%z) * TILE_DIMz + lk
        xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
        yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
        zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
        myblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
        ii=i-xblock*TILE_DIMx+2*TILE_DIMx
        jj=j-yblock*TILE_DIMy+2*TILE_DIMy
        kk=k-zblock*TILE_DIMz+2*TILE_DIMz
        myphi(li,lj,lk)=real(phifields_s(ii,jj,kk,1,myblock),kind=db)
      enddo
      call syncthreads

      li=threadIdx%x
      lj=threadIdx%y
      lk=threadIdx%z
      i = (blockIdx%x) * TILE_DIMx + li
      j = (blockIdx%y) * TILE_DIMy + lj
      k = (blockIdx%z) * TILE_DIMz + lk

      if(abs(isfluid(i,j,k)) .ne. 1)return

      !internal-node block of the tile
      myblock=(blockIdx%x+1)+(blockIdx%y+1)*nxblock_d+(blockIdx%z+1)*nxyblock_d+1
      ii=li
      jj=lj
      kk=lk

	  grad_fix=3.0_db*(p1d3q27*(myphi(li+1,lj,lk)-myphi(li-1,lj,lk)) + &
		 p2d3q27*( (myphi(li+1,lj+1,lk)-myphi(li-1,lj-1,lk))+ &
		 (myphi(li+1,lj-1,lk)-myphi(li-1,lj+1,lk))+ &
		 (myphi(li+1,lj,lk+1)-myphi(li-1,lj,lk-1))+ &
		 (myphi(li+1,lj,lk-1)-myphi(li-1,lj,lk+1)) )  + &
		 p3d3q27*((myphi(li+1,lj+1,lk+1)-myphi(li-1,lj-1,lk-1))+ &
		 (myphi(li+1,lj-1,lk-1)-myphi(li-1,lj+1,lk+1))+ &
		 (myphi(li+1,lj-1,lk+1)-myphi(li-1,lj+1,lk-1))+ &
		 (myphi(li+1,lj+1,lk-1)-myphi(li-1,lj-1,lk+1))))

	  grad_fiy=3.0_db*(p1d3q27*(myphi(li,lj+1,lk)-myphi(li,lj-1,lk)) + &
		 p2d3q27*((myphi(li+1,lj+1,lk)-myphi(li-1,lj-1,lk))+ &
		 (myphi(li-1,lj+1,lk)-myphi(li+1,lj-1,lk))+ &
		 (myphi(li,lj+1,lk+1)-myphi(li,lj-1,lk-1))+ &
		 (myphi(li,lj+1,lk-1)-myphi(li,lj-1,lk+1)) ) + &
		 p3d3q27*((myphi(li+1,lj+1,lk+1)-myphi(li-1,lj-1,lk-1))+ &
		 (myphi(li-1,lj+1,lk-1)-myphi(li+1,lj-1,lk+1))+ &
		 (myphi(li+1,lj+1,lk-1)-myphi(li-1,lj-1,lk+1))+ &
		 (myphi(li-1,lj+1,lk+1)-myphi(li+1,lj-1,lk-1))))

	  grad_fiz=3.0_db*(p1d3q27*(myphi(li,lj,lk+1)-myphi(li,lj,lk-1)) + &
		 p2d3q27*((myphi(li+1,lj,lk+1)-myphi(li-1,lj,lk-1))+ &
		 (myphi(li-1,lj,lk+1)-myphi(li+1,lj,lk-1))+ &
		 (myphi(li,lj+1,lk+1)-myphi(li,lj-1,lk-1))+ &
		 (myphi(li,lj-1,lk+1)-myphi(li,lj+1,lk-1)) ) + &
		 p3d3q27*((myphi(li+1,lj+1,lk+1)-myphi(li-1,lj-1,lk-1)) &
		 +(myphi(li-1,lj-1,lk+1)-myphi(li+1,lj+1,lk-1))+ &
		 (myphi(li+1,lj-1,lk+1)-myphi(li-1,lj+1,lk-1))+ &
		 (myphi(li-1,lj+1,lk+1)-myphi(li+1,lj-1,lk-1))))
      
	  mod_grad=sqrt(grad_fix**TWO+grad_fiy**TWO+grad_fiz**TWO)
#ifdef CSF
	  normal_denom=sqrt(mod_grad*mod_grad+1.0e-24_db)
#else
	  normal_denom=mod_grad+1.0e-9_db
#endif

	  auxfields_s(ii,jj,kk,1,myblock)= &
	   real(grad_fix/normal_denom,kind=strdb)
	  auxfields_s(ii,jj,kk,2,myblock)= &
	   real(grad_fiy/normal_denom,kind=strdb)
	  auxfields_s(ii,jj,kk,3,myblock)= &
	   real(grad_fiz/normal_denom,kind=strdb)
	  
	  auxfields_s(ii,jj,kk,4,myblock)=mod_grad 

	  auxfields_s(ii,jj,kk,5,myblock)=real( &
	   myphi(li,lj,lk)*(1.0_db-myphi(li,lj,lk))*(grad_fix/normal_denom),kind=strdb)
	  auxfields_s(ii,jj,kk,6,myblock)=real( &
	   myphi(li,lj,lk)*(1.0_db-myphi(li,lj,lk))*(grad_fiy/normal_denom),kind=strdb)
	  auxfields_s(ii,jj,kk,7,myblock)=real( &
	   myphi(li,lj,lk)*(1.0_db-myphi(li,lj,lk))*(grad_fiz/normal_denom),kind=strdb)
	   
      !lap_phi here
      locauxfields_s(ii,jj,kk,1,myblock)=real( &
                   (2.0_db*invcssq)*(myphi(li,lj,lk)*(p0d3q27-1.0_db) + &
                   ( p1d3q27*(myphi(li+1,lj,lk)+myphi(li-1,lj,lk) + &
                   myphi(li,lj+1,lk)+myphi(li,lj-1,lk) + &
                   myphi(li,lj,lk+1)+myphi(li,lj,lk-1)) + &
                   p2d3q27*( (myphi(li+1,lj+1,lk)+myphi(li-1,lj-1,lk))+ &
                   (myphi(li+1,lj-1,lk)+myphi(li-1,lj+1,lk))+ &
                   (myphi(li+1,lj,lk+1)+myphi(li-1,lj,lk-1))+ &
                   (myphi(li+1,lj,lk-1)+myphi(li-1,lj,lk+1)) + &
                   (myphi(li,lj+1,lk+1)+myphi(li,lj-1,lk-1))+ &
                   (myphi(li,lj+1,lk-1)+myphi(li,lj-1,lk+1)) )  + &
                   p3d3q27*((myphi(li+1,lj+1,lk+1)+myphi(li-1,lj-1,lk-1))+ &
                   (myphi(li+1,lj-1,lk-1)+myphi(li-1,lj+1,lk+1))+ &
                   (myphi(li+1,lj-1,lk+1)+myphi(li-1,lj+1,lk-1))+ &
                   (myphi(li+1,lj+1,lk-1)+myphi(li-1,lj-1,lk+1))))),kind=strdb)
      
      return
      
   endsubroutine compute_norm_interface_kernel_int
   
   ! _xplus: last x tile layer (xblock=nxblock_d-1), y,z tiles 2..nyb-1,2..nzb-1; grid (1,nyb-2,nzb-2)
   attributes(global) subroutine compute_norm_interface_kernel_xplus(step, &
      iprobe,jprobe,kprobe,flop,nx,ny,nz,coords,isfluid,rho_r,rho_b, &
      ntotphifields,ntotauxfields,ntotlocauxfields,phifields_s,auxfields_s,locauxfields_s)

      implicit none
      integer :: step,iprobe,jprobe,kprobe,flop,nx,ny,nz
      integer, dimension(3) :: coords
      integer(kind=isf), dimension(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff) :: isfluid
      real(kind=db) :: rho_r,rho_b
      integer :: ntotphifields,ntotauxfields,ntotlocauxfields
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields,nblocks_d) :: phifields_s
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nauxfields,nblocks_d) :: auxfields_s
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nlocauxfields,nblocks_d) :: locauxfields_s
      
      real(kind=db), shared :: myphi(0:TILE_DIMx+1,0:TILE_DIMy+1,0:TILE_DIMz+1)
      real(kind=db):: grad_fix,grad_fiy,grad_fiz,mod_grad,normal_denom
      
      integer :: i,j,k,myblock,intblock
      integer :: ii,jj,kk
      integer :: li,lj,lk
      integer :: xblock,yblock,zblock
      integer, parameter :: nhx=TILE_DIMx+2, nhy=TILE_DIMy+2, nhz=TILE_DIMz+2
      integer, parameter :: nthreads=TILE_DIMx*TILE_DIMy*TILE_DIMz
      integer :: ithread,ih
!      integer :: gi,gj,gk
      
      
      ! one thread per interior node; the threads first load the tile and its halo into
      ! shared memory, nodes ih=ithread, ithread+nthreads, ... in x-fastest order
      ithread=(threadIdx%x-1)+TILE_DIMx*(threadIdx%y-1)+TILE_DIMx*TILE_DIMy*(threadIdx%z-1)
      do ih=ithread,nhx*nhy*nhz-1,nthreads
        lk=ih/(nhx*nhy)
        lj=(ih-lk*nhx*nhy)/nhx
        li=ih-lk*nhx*nhy-lj*nhx
        i = ((nxblock_d-2)-1) * TILE_DIMx + li
        j = (blockIdx%y) * TILE_DIMy + lj
        k = (blockIdx%z) * TILE_DIMz + lk
        xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
        yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
        zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
        myblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
        ii=i-xblock*TILE_DIMx+2*TILE_DIMx
        jj=j-yblock*TILE_DIMy+2*TILE_DIMy
        kk=k-zblock*TILE_DIMz+2*TILE_DIMz
        myphi(li,lj,lk)=real(phifields_s(ii,jj,kk,1,myblock),kind=db)
      enddo
      call syncthreads

      li=threadIdx%x
      lj=threadIdx%y
      lk=threadIdx%z
      i = ((nxblock_d-2)-1) * TILE_DIMx + li
      j = (blockIdx%y) * TILE_DIMy + lj
      k = (blockIdx%z) * TILE_DIMz + lk

      if(abs(isfluid(i,j,k)) .ne. 1)return

      !internal-node block of the tile
      myblock=(nxblock_d-2)+(blockIdx%y+1)*nxblock_d+(blockIdx%z+1)*nxyblock_d+1
      ii=li
      jj=lj
      kk=lk

	  grad_fix=3.0_db*(p1d3q27*(myphi(li+1,lj,lk)-myphi(li-1,lj,lk)) + &
		 p2d3q27*( (myphi(li+1,lj+1,lk)-myphi(li-1,lj-1,lk))+ &
		 (myphi(li+1,lj-1,lk)-myphi(li-1,lj+1,lk))+ &
		 (myphi(li+1,lj,lk+1)-myphi(li-1,lj,lk-1))+ &
		 (myphi(li+1,lj,lk-1)-myphi(li-1,lj,lk+1)) )  + &
		 p3d3q27*((myphi(li+1,lj+1,lk+1)-myphi(li-1,lj-1,lk-1))+ &
		 (myphi(li+1,lj-1,lk-1)-myphi(li-1,lj+1,lk+1))+ &
		 (myphi(li+1,lj-1,lk+1)-myphi(li-1,lj+1,lk-1))+ &
		 (myphi(li+1,lj+1,lk-1)-myphi(li-1,lj-1,lk+1))))

	  grad_fiy=3.0_db*(p1d3q27*(myphi(li,lj+1,lk)-myphi(li,lj-1,lk)) + &
		 p2d3q27*((myphi(li+1,lj+1,lk)-myphi(li-1,lj-1,lk))+ &
		 (myphi(li-1,lj+1,lk)-myphi(li+1,lj-1,lk))+ &
		 (myphi(li,lj+1,lk+1)-myphi(li,lj-1,lk-1))+ &
		 (myphi(li,lj+1,lk-1)-myphi(li,lj-1,lk+1)) ) + &
		 p3d3q27*((myphi(li+1,lj+1,lk+1)-myphi(li-1,lj-1,lk-1))+ &
		 (myphi(li-1,lj+1,lk-1)-myphi(li+1,lj-1,lk+1))+ &
		 (myphi(li+1,lj+1,lk-1)-myphi(li-1,lj-1,lk+1))+ &
		 (myphi(li-1,lj+1,lk+1)-myphi(li+1,lj-1,lk-1))))

	  grad_fiz=3.0_db*(p1d3q27*(myphi(li,lj,lk+1)-myphi(li,lj,lk-1)) + &
		 p2d3q27*((myphi(li+1,lj,lk+1)-myphi(li-1,lj,lk-1))+ &
		 (myphi(li-1,lj,lk+1)-myphi(li+1,lj,lk-1))+ &
		 (myphi(li,lj+1,lk+1)-myphi(li,lj-1,lk-1))+ &
		 (myphi(li,lj-1,lk+1)-myphi(li,lj+1,lk-1)) ) + &
		 p3d3q27*((myphi(li+1,lj+1,lk+1)-myphi(li-1,lj-1,lk-1)) &
		 +(myphi(li-1,lj-1,lk+1)-myphi(li+1,lj+1,lk-1))+ &
		 (myphi(li+1,lj-1,lk+1)-myphi(li-1,lj+1,lk-1))+ &
		 (myphi(li-1,lj+1,lk+1)-myphi(li+1,lj-1,lk-1))))
      
	  mod_grad=sqrt(grad_fix**TWO+grad_fiy**TWO+grad_fiz**TWO)
#ifdef CSF
	  normal_denom=sqrt(mod_grad*mod_grad+1.0e-24_db)
#else
	  normal_denom=mod_grad+1.0e-9_db
#endif

	  auxfields_s(ii,jj,kk,1,myblock)= &
	   real(grad_fix/normal_denom,kind=strdb)
	  auxfields_s(ii,jj,kk,2,myblock)= &
	   real(grad_fiy/normal_denom,kind=strdb)
	  auxfields_s(ii,jj,kk,3,myblock)= &
	   real(grad_fiz/normal_denom,kind=strdb)
	  
	  auxfields_s(ii,jj,kk,4,myblock)=mod_grad 

	  auxfields_s(ii,jj,kk,5,myblock)=real( &
	   myphi(li,lj,lk)*(1.0_db-myphi(li,lj,lk))*(grad_fix/normal_denom),kind=strdb)
	  auxfields_s(ii,jj,kk,6,myblock)=real( &
	   myphi(li,lj,lk)*(1.0_db-myphi(li,lj,lk))*(grad_fiy/normal_denom),kind=strdb)
	  auxfields_s(ii,jj,kk,7,myblock)=real( &
	   myphi(li,lj,lk)*(1.0_db-myphi(li,lj,lk))*(grad_fiz/normal_denom),kind=strdb)
	   
      !lap_phi here
      locauxfields_s(ii,jj,kk,1,myblock)=real( &
                   (2.0_db*invcssq)*(myphi(li,lj,lk)*(p0d3q27-1.0_db) + &
                   ( p1d3q27*(myphi(li+1,lj,lk)+myphi(li-1,lj,lk) + &
                   myphi(li,lj+1,lk)+myphi(li,lj-1,lk) + &
                   myphi(li,lj,lk+1)+myphi(li,lj,lk-1)) + &
                   p2d3q27*( (myphi(li+1,lj+1,lk)+myphi(li-1,lj-1,lk))+ &
                   (myphi(li+1,lj-1,lk)+myphi(li-1,lj+1,lk))+ &
                   (myphi(li+1,lj,lk+1)+myphi(li-1,lj,lk-1))+ &
                   (myphi(li+1,lj,lk-1)+myphi(li-1,lj,lk+1)) + &
                   (myphi(li,lj+1,lk+1)+myphi(li,lj-1,lk-1))+ &
                   (myphi(li,lj+1,lk-1)+myphi(li,lj-1,lk+1)) )  + &
                   p3d3q27*((myphi(li+1,lj+1,lk+1)+myphi(li-1,lj-1,lk-1))+ &
                   (myphi(li+1,lj-1,lk-1)+myphi(li-1,lj+1,lk+1))+ &
                   (myphi(li+1,lj-1,lk+1)+myphi(li-1,lj+1,lk-1))+ &
                   (myphi(li+1,lj+1,lk-1)+myphi(li-1,lj-1,lk+1))))),kind=strdb) 
      
      return
      
   endsubroutine compute_norm_interface_kernel_xplus
   
   ! _xminus: first x tile layer (xblock=2), y,z tiles 2..nyb-1,2..nzb-1; grid (1,nyb-2,nzb-2)
   attributes(global) subroutine compute_norm_interface_kernel_xminus(step, &
      iprobe,jprobe,kprobe,flop,nx,ny,nz,coords,isfluid,rho_r,rho_b, &
      ntotphifields,ntotauxfields,ntotlocauxfields,phifields_s,auxfields_s,locauxfields_s)

      implicit none
      integer :: step,iprobe,jprobe,kprobe,flop,nx,ny,nz
      integer, dimension(3) :: coords
      integer(kind=isf), dimension(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff) :: isfluid
      real(kind=db) :: rho_r,rho_b
      integer :: ntotphifields,ntotauxfields,ntotlocauxfields
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields,nblocks_d) :: phifields_s
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nauxfields,nblocks_d) :: auxfields_s
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nlocauxfields,nblocks_d) :: locauxfields_s
      
      real(kind=db), shared :: myphi(0:TILE_DIMx+1,0:TILE_DIMy+1,0:TILE_DIMz+1)
      real(kind=db):: grad_fix,grad_fiy,grad_fiz,mod_grad,normal_denom
      
      integer :: i,j,k,myblock,intblock
      integer :: ii,jj,kk
      integer :: li,lj,lk
      integer :: xblock,yblock,zblock
      integer, parameter :: nhx=TILE_DIMx+2, nhy=TILE_DIMy+2, nhz=TILE_DIMz+2
      integer, parameter :: nthreads=TILE_DIMx*TILE_DIMy*TILE_DIMz
      integer :: ithread,ih
!      integer :: gi,gj,gk

      
      ! one thread per interior node; the threads first load the tile and its halo into
      ! shared memory, nodes ih=ithread, ithread+nthreads, ... in x-fastest order
      ithread=(threadIdx%x-1)+TILE_DIMx*(threadIdx%y-1)+TILE_DIMx*TILE_DIMy*(threadIdx%z-1)
      do ih=ithread,nhx*nhy*nhz-1,nthreads
        lk=ih/(nhx*nhy)
        lj=(ih-lk*nhx*nhy)/nhx
        li=ih-lk*nhx*nhy-lj*nhx
        i = li
        j = (blockIdx%y) * TILE_DIMy + lj
        k = (blockIdx%z) * TILE_DIMz + lk
        xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
        yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
        zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
        myblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
        ii=i-xblock*TILE_DIMx+2*TILE_DIMx
        jj=j-yblock*TILE_DIMy+2*TILE_DIMy
        kk=k-zblock*TILE_DIMz+2*TILE_DIMz
        myphi(li,lj,lk)=real(phifields_s(ii,jj,kk,1,myblock),kind=db)
      enddo
      call syncthreads

      li=threadIdx%x
      lj=threadIdx%y
      lk=threadIdx%z
      i = li
      j = (blockIdx%y) * TILE_DIMy + lj
      k = (blockIdx%z) * TILE_DIMz + lk

      if(abs(isfluid(i,j,k)) .ne. 1)return

      !internal-node block of the tile
      myblock=1+(blockIdx%y+1)*nxblock_d+(blockIdx%z+1)*nxyblock_d+1
      ii=li
      jj=lj
      kk=lk

	  grad_fix=3.0_db*(p1d3q27*(myphi(li+1,lj,lk)-myphi(li-1,lj,lk)) + &
		 p2d3q27*( (myphi(li+1,lj+1,lk)-myphi(li-1,lj-1,lk))+ &
		 (myphi(li+1,lj-1,lk)-myphi(li-1,lj+1,lk))+ &
		 (myphi(li+1,lj,lk+1)-myphi(li-1,lj,lk-1))+ &
		 (myphi(li+1,lj,lk-1)-myphi(li-1,lj,lk+1)) )  + &
		 p3d3q27*((myphi(li+1,lj+1,lk+1)-myphi(li-1,lj-1,lk-1))+ &
		 (myphi(li+1,lj-1,lk-1)-myphi(li-1,lj+1,lk+1))+ &
		 (myphi(li+1,lj-1,lk+1)-myphi(li-1,lj+1,lk-1))+ &
		 (myphi(li+1,lj+1,lk-1)-myphi(li-1,lj-1,lk+1))))

	  grad_fiy=3.0_db*(p1d3q27*(myphi(li,lj+1,lk)-myphi(li,lj-1,lk)) + &
		 p2d3q27*((myphi(li+1,lj+1,lk)-myphi(li-1,lj-1,lk))+ &
		 (myphi(li-1,lj+1,lk)-myphi(li+1,lj-1,lk))+ &
		 (myphi(li,lj+1,lk+1)-myphi(li,lj-1,lk-1))+ &
		 (myphi(li,lj+1,lk-1)-myphi(li,lj-1,lk+1)) ) + &
		 p3d3q27*((myphi(li+1,lj+1,lk+1)-myphi(li-1,lj-1,lk-1))+ &
		 (myphi(li-1,lj+1,lk-1)-myphi(li+1,lj-1,lk+1))+ &
		 (myphi(li+1,lj+1,lk-1)-myphi(li-1,lj-1,lk+1))+ &
		 (myphi(li-1,lj+1,lk+1)-myphi(li+1,lj-1,lk-1))))

	  grad_fiz=3.0_db*(p1d3q27*(myphi(li,lj,lk+1)-myphi(li,lj,lk-1)) + &
		 p2d3q27*((myphi(li+1,lj,lk+1)-myphi(li-1,lj,lk-1))+ &
		 (myphi(li-1,lj,lk+1)-myphi(li+1,lj,lk-1))+ &
		 (myphi(li,lj+1,lk+1)-myphi(li,lj-1,lk-1))+ &
		 (myphi(li,lj-1,lk+1)-myphi(li,lj+1,lk-1)) ) + &
		 p3d3q27*((myphi(li+1,lj+1,lk+1)-myphi(li-1,lj-1,lk-1)) &
		 +(myphi(li-1,lj-1,lk+1)-myphi(li+1,lj+1,lk-1))+ &
		 (myphi(li+1,lj-1,lk+1)-myphi(li-1,lj+1,lk-1))+ &
		 (myphi(li-1,lj+1,lk+1)-myphi(li+1,lj-1,lk-1))))
      
	  mod_grad=sqrt(grad_fix**TWO+grad_fiy**TWO+grad_fiz**TWO)
#ifdef CSF
	  normal_denom=sqrt(mod_grad*mod_grad+1.0e-24_db)
#else
	  normal_denom=mod_grad+1.0e-9_db
#endif

	  auxfields_s(ii,jj,kk,1,myblock)= &
	   real(grad_fix/normal_denom,kind=strdb)
	  auxfields_s(ii,jj,kk,2,myblock)= &
	   real(grad_fiy/normal_denom,kind=strdb)
	  auxfields_s(ii,jj,kk,3,myblock)= &
	   real(grad_fiz/normal_denom,kind=strdb)
	  
	  auxfields_s(ii,jj,kk,4,myblock)=mod_grad 

	  auxfields_s(ii,jj,kk,5,myblock)=real( &
	   myphi(li,lj,lk)*(1.0_db-myphi(li,lj,lk))*(grad_fix/normal_denom),kind=strdb)
	  auxfields_s(ii,jj,kk,6,myblock)=real( &
	   myphi(li,lj,lk)*(1.0_db-myphi(li,lj,lk))*(grad_fiy/normal_denom),kind=strdb)
	  auxfields_s(ii,jj,kk,7,myblock)=real( &
	   myphi(li,lj,lk)*(1.0_db-myphi(li,lj,lk))*(grad_fiz/normal_denom),kind=strdb)
	   
      !lap_phi here
      locauxfields_s(ii,jj,kk,1,myblock)=real( &
                   (2.0_db*invcssq)*(myphi(li,lj,lk)*(p0d3q27-1.0_db) + &
                   ( p1d3q27*(myphi(li+1,lj,lk)+myphi(li-1,lj,lk) + &
                   myphi(li,lj+1,lk)+myphi(li,lj-1,lk) + &
                   myphi(li,lj,lk+1)+myphi(li,lj,lk-1)) + &
                   p2d3q27*( (myphi(li+1,lj+1,lk)+myphi(li-1,lj-1,lk))+ &
                   (myphi(li+1,lj-1,lk)+myphi(li-1,lj+1,lk))+ &
                   (myphi(li+1,lj,lk+1)+myphi(li-1,lj,lk-1))+ &
                   (myphi(li+1,lj,lk-1)+myphi(li-1,lj,lk+1)) + &
                   (myphi(li,lj+1,lk+1)+myphi(li,lj-1,lk-1))+ &
                   (myphi(li,lj+1,lk-1)+myphi(li,lj-1,lk+1)) )  + &
                   p3d3q27*((myphi(li+1,lj+1,lk+1)+myphi(li-1,lj-1,lk-1))+ &
                   (myphi(li+1,lj-1,lk-1)+myphi(li-1,lj+1,lk+1))+ &
                   (myphi(li+1,lj-1,lk+1)+myphi(li-1,lj+1,lk-1))+ &
                   (myphi(li+1,lj+1,lk-1)+myphi(li-1,lj-1,lk+1))))),kind=strdb)
      
      return
      
   endsubroutine compute_norm_interface_kernel_xminus
   
   ! _yplus: last y tile layer (yblock=nyblock_d-1), all x tiles, z tiles 2..nzb-1; grid (nxb,1,nzb-2)
   attributes(global) subroutine compute_norm_interface_kernel_yplus(step, &
      iprobe,jprobe,kprobe,flop,nx,ny,nz,coords,isfluid,rho_r,rho_b, &
      ntotphifields,ntotauxfields,ntotlocauxfields,phifields_s,auxfields_s,locauxfields_s)

      implicit none
      integer :: step,iprobe,jprobe,kprobe,flop,nx,ny,nz
      integer, dimension(3) :: coords
      integer(kind=isf), dimension(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff) :: isfluid
      real(kind=db) :: rho_r,rho_b
      integer :: ntotphifields,ntotauxfields,ntotlocauxfields
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields,nblocks_d) :: phifields_s
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nauxfields,nblocks_d) :: auxfields_s
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nlocauxfields,nblocks_d) :: locauxfields_s
      
      real(kind=db), shared :: myphi(0:TILE_DIMx+1,0:TILE_DIMy+1,0:TILE_DIMz+1)
      real(kind=db):: grad_fix,grad_fiy,grad_fiz,mod_grad,normal_denom
      
      integer :: i,j,k,myblock,intblock
      integer :: ii,jj,kk
      integer :: li,lj,lk
      integer :: xblock,yblock,zblock
      integer, parameter :: nhx=TILE_DIMx+2, nhy=TILE_DIMy+2, nhz=TILE_DIMz+2
      integer, parameter :: nthreads=TILE_DIMx*TILE_DIMy*TILE_DIMz
      integer :: ithread,ih
!      integer :: gi,gj,gk

      
      ! one thread per interior node; the threads first load the tile and its halo into
      ! shared memory, nodes ih=ithread, ithread+nthreads, ... in x-fastest order
      ithread=(threadIdx%x-1)+TILE_DIMx*(threadIdx%y-1)+TILE_DIMx*TILE_DIMy*(threadIdx%z-1)
      do ih=ithread,nhx*nhy*nhz-1,nthreads
        lk=ih/(nhx*nhy)
        lj=(ih-lk*nhx*nhy)/nhx
        li=ih-lk*nhx*nhy-lj*nhx
        i = (blockIdx%x-1) * TILE_DIMx + li
        j = ((nyblock_d-2)-1) * TILE_DIMy + lj
        k = (blockIdx%z) * TILE_DIMz + lk
        xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
        yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
        zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
        myblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
        ii=i-xblock*TILE_DIMx+2*TILE_DIMx
        jj=j-yblock*TILE_DIMy+2*TILE_DIMy
        kk=k-zblock*TILE_DIMz+2*TILE_DIMz
        myphi(li,lj,lk)=real(phifields_s(ii,jj,kk,1,myblock),kind=db)
      enddo
      call syncthreads

      li=threadIdx%x
      lj=threadIdx%y
      lk=threadIdx%z
      i = (blockIdx%x-1) * TILE_DIMx + li
      j = ((nyblock_d-2)-1) * TILE_DIMy + lj
      k = (blockIdx%z) * TILE_DIMz + lk

      if(abs(isfluid(i,j,k)) .ne. 1)return

      !internal-node block of the tile
      myblock=blockIdx%x+(nyblock_d-2)*nxblock_d+(blockIdx%z+1)*nxyblock_d+1
      ii=li
      jj=lj
      kk=lk

	  grad_fix=3.0_db*(p1d3q27*(myphi(li+1,lj,lk)-myphi(li-1,lj,lk)) + &
		 p2d3q27*( (myphi(li+1,lj+1,lk)-myphi(li-1,lj-1,lk))+ &
		 (myphi(li+1,lj-1,lk)-myphi(li-1,lj+1,lk))+ &
		 (myphi(li+1,lj,lk+1)-myphi(li-1,lj,lk-1))+ &
		 (myphi(li+1,lj,lk-1)-myphi(li-1,lj,lk+1)) )  + &
		 p3d3q27*((myphi(li+1,lj+1,lk+1)-myphi(li-1,lj-1,lk-1))+ &
		 (myphi(li+1,lj-1,lk-1)-myphi(li-1,lj+1,lk+1))+ &
		 (myphi(li+1,lj-1,lk+1)-myphi(li-1,lj+1,lk-1))+ &
		 (myphi(li+1,lj+1,lk-1)-myphi(li-1,lj-1,lk+1))))

	  grad_fiy=3.0_db*(p1d3q27*(myphi(li,lj+1,lk)-myphi(li,lj-1,lk)) + &
		 p2d3q27*((myphi(li+1,lj+1,lk)-myphi(li-1,lj-1,lk))+ &
		 (myphi(li-1,lj+1,lk)-myphi(li+1,lj-1,lk))+ &
		 (myphi(li,lj+1,lk+1)-myphi(li,lj-1,lk-1))+ &
		 (myphi(li,lj+1,lk-1)-myphi(li,lj-1,lk+1)) ) + &
		 p3d3q27*((myphi(li+1,lj+1,lk+1)-myphi(li-1,lj-1,lk-1))+ &
		 (myphi(li-1,lj+1,lk-1)-myphi(li+1,lj-1,lk+1))+ &
		 (myphi(li+1,lj+1,lk-1)-myphi(li-1,lj-1,lk+1))+ &
		 (myphi(li-1,lj+1,lk+1)-myphi(li+1,lj-1,lk-1))))

	  grad_fiz=3.0_db*(p1d3q27*(myphi(li,lj,lk+1)-myphi(li,lj,lk-1)) + &
		 p2d3q27*((myphi(li+1,lj,lk+1)-myphi(li-1,lj,lk-1))+ &
		 (myphi(li-1,lj,lk+1)-myphi(li+1,lj,lk-1))+ &
		 (myphi(li,lj+1,lk+1)-myphi(li,lj-1,lk-1))+ &
		 (myphi(li,lj-1,lk+1)-myphi(li,lj+1,lk-1)) ) + &
		 p3d3q27*((myphi(li+1,lj+1,lk+1)-myphi(li-1,lj-1,lk-1)) &
		 +(myphi(li-1,lj-1,lk+1)-myphi(li+1,lj+1,lk-1))+ &
		 (myphi(li+1,lj-1,lk+1)-myphi(li-1,lj+1,lk-1))+ &
		 (myphi(li-1,lj+1,lk+1)-myphi(li+1,lj-1,lk-1))))
      
	  mod_grad=sqrt(grad_fix**TWO+grad_fiy**TWO+grad_fiz**TWO)
#ifdef CSF
	  normal_denom=sqrt(mod_grad*mod_grad+1.0e-24_db)
#else
	  normal_denom=mod_grad+1.0e-9_db
#endif

	  auxfields_s(ii,jj,kk,1,myblock)= &
	   real(grad_fix/normal_denom,kind=strdb)
	  auxfields_s(ii,jj,kk,2,myblock)= &
	   real(grad_fiy/normal_denom,kind=strdb)
	  auxfields_s(ii,jj,kk,3,myblock)= &
	   real(grad_fiz/normal_denom,kind=strdb)
	  
	  auxfields_s(ii,jj,kk,4,myblock)=mod_grad 

	  auxfields_s(ii,jj,kk,5,myblock)=real( &
	   myphi(li,lj,lk)*(1.0_db-myphi(li,lj,lk))*(grad_fix/normal_denom),kind=strdb)
	  auxfields_s(ii,jj,kk,6,myblock)=real( &
	   myphi(li,lj,lk)*(1.0_db-myphi(li,lj,lk))*(grad_fiy/normal_denom),kind=strdb)
	  auxfields_s(ii,jj,kk,7,myblock)=real( &
	   myphi(li,lj,lk)*(1.0_db-myphi(li,lj,lk))*(grad_fiz/normal_denom),kind=strdb)
	   
      !lap_phi here
      locauxfields_s(ii,jj,kk,1,myblock)=real( &
                   (2.0_db*invcssq)*(myphi(li,lj,lk)*(p0d3q27-1.0_db) + &
                   ( p1d3q27*(myphi(li+1,lj,lk)+myphi(li-1,lj,lk) + &
                   myphi(li,lj+1,lk)+myphi(li,lj-1,lk) + &
                   myphi(li,lj,lk+1)+myphi(li,lj,lk-1)) + &
                   p2d3q27*( (myphi(li+1,lj+1,lk)+myphi(li-1,lj-1,lk))+ &
                   (myphi(li+1,lj-1,lk)+myphi(li-1,lj+1,lk))+ &
                   (myphi(li+1,lj,lk+1)+myphi(li-1,lj,lk-1))+ &
                   (myphi(li+1,lj,lk-1)+myphi(li-1,lj,lk+1)) + &
                   (myphi(li,lj+1,lk+1)+myphi(li,lj-1,lk-1))+ &
                   (myphi(li,lj+1,lk-1)+myphi(li,lj-1,lk+1)) )  + &
                   p3d3q27*((myphi(li+1,lj+1,lk+1)+myphi(li-1,lj-1,lk-1))+ &
                   (myphi(li+1,lj-1,lk-1)+myphi(li-1,lj+1,lk+1))+ &
                   (myphi(li+1,lj-1,lk+1)+myphi(li-1,lj+1,lk-1))+ &
                   (myphi(li+1,lj+1,lk-1)+myphi(li-1,lj-1,lk+1))))),kind=strdb) 
      
      return
      
   endsubroutine compute_norm_interface_kernel_yplus
   
   ! _yminus: first y tile layer (yblock=2), all x tiles, z tiles 2..nzb-1; grid (nxb,1,nzb-2)
   attributes(global) subroutine compute_norm_interface_kernel_yminus(step, &
      iprobe,jprobe,kprobe,flop,nx,ny,nz,coords,isfluid,rho_r,rho_b, &
      ntotphifields,ntotauxfields,ntotlocauxfields,phifields_s,auxfields_s,locauxfields_s)

      implicit none
      integer :: step,iprobe,jprobe,kprobe,flop,nx,ny,nz
      integer, dimension(3) :: coords
      integer(kind=isf), dimension(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff) :: isfluid
      real(kind=db) :: rho_r,rho_b
      integer :: ntotphifields,ntotauxfields,ntotlocauxfields
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields,nblocks_d) :: phifields_s
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nauxfields,nblocks_d) :: auxfields_s
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nlocauxfields,nblocks_d) :: locauxfields_s
      
      real(kind=db), shared :: myphi(0:TILE_DIMx+1,0:TILE_DIMy+1,0:TILE_DIMz+1)
      real(kind=db):: grad_fix,grad_fiy,grad_fiz,mod_grad,normal_denom
      
      integer :: i,j,k,myblock,intblock
      integer :: ii,jj,kk
      integer :: li,lj,lk
      integer :: xblock,yblock,zblock
      integer, parameter :: nhx=TILE_DIMx+2, nhy=TILE_DIMy+2, nhz=TILE_DIMz+2
      integer, parameter :: nthreads=TILE_DIMx*TILE_DIMy*TILE_DIMz
      integer :: ithread,ih
!      integer :: gi,gj,gk

      
      ! one thread per interior node; the threads first load the tile and its halo into
      ! shared memory, nodes ih=ithread, ithread+nthreads, ... in x-fastest order
      ithread=(threadIdx%x-1)+TILE_DIMx*(threadIdx%y-1)+TILE_DIMx*TILE_DIMy*(threadIdx%z-1)
      do ih=ithread,nhx*nhy*nhz-1,nthreads
        lk=ih/(nhx*nhy)
        lj=(ih-lk*nhx*nhy)/nhx
        li=ih-lk*nhx*nhy-lj*nhx
        i = (blockIdx%x-1) * TILE_DIMx + li
        j = lj
        k = (blockIdx%z) * TILE_DIMz + lk
        xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
        yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
        zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
        myblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
        ii=i-xblock*TILE_DIMx+2*TILE_DIMx
        jj=j-yblock*TILE_DIMy+2*TILE_DIMy
        kk=k-zblock*TILE_DIMz+2*TILE_DIMz
        myphi(li,lj,lk)=real(phifields_s(ii,jj,kk,1,myblock),kind=db)
      enddo
      call syncthreads

      li=threadIdx%x
      lj=threadIdx%y
      lk=threadIdx%z
      i = (blockIdx%x-1) * TILE_DIMx + li
      j = lj
      k = (blockIdx%z) * TILE_DIMz + lk

      if(abs(isfluid(i,j,k)) .ne. 1)return

      !internal-node block of the tile
      myblock=blockIdx%x+1*nxblock_d+(blockIdx%z+1)*nxyblock_d+1
      ii=li
      jj=lj
      kk=lk

	  grad_fix=3.0_db*(p1d3q27*(myphi(li+1,lj,lk)-myphi(li-1,lj,lk)) + &
		 p2d3q27*( (myphi(li+1,lj+1,lk)-myphi(li-1,lj-1,lk))+ &
		 (myphi(li+1,lj-1,lk)-myphi(li-1,lj+1,lk))+ &
		 (myphi(li+1,lj,lk+1)-myphi(li-1,lj,lk-1))+ &
		 (myphi(li+1,lj,lk-1)-myphi(li-1,lj,lk+1)) )  + &
		 p3d3q27*((myphi(li+1,lj+1,lk+1)-myphi(li-1,lj-1,lk-1))+ &
		 (myphi(li+1,lj-1,lk-1)-myphi(li-1,lj+1,lk+1))+ &
		 (myphi(li+1,lj-1,lk+1)-myphi(li-1,lj+1,lk-1))+ &
		 (myphi(li+1,lj+1,lk-1)-myphi(li-1,lj-1,lk+1))))

	  grad_fiy=3.0_db*(p1d3q27*(myphi(li,lj+1,lk)-myphi(li,lj-1,lk)) + &
		 p2d3q27*((myphi(li+1,lj+1,lk)-myphi(li-1,lj-1,lk))+ &
		 (myphi(li-1,lj+1,lk)-myphi(li+1,lj-1,lk))+ &
		 (myphi(li,lj+1,lk+1)-myphi(li,lj-1,lk-1))+ &
		 (myphi(li,lj+1,lk-1)-myphi(li,lj-1,lk+1)) ) + &
		 p3d3q27*((myphi(li+1,lj+1,lk+1)-myphi(li-1,lj-1,lk-1))+ &
		 (myphi(li-1,lj+1,lk-1)-myphi(li+1,lj-1,lk+1))+ &
		 (myphi(li+1,lj+1,lk-1)-myphi(li-1,lj-1,lk+1))+ &
		 (myphi(li-1,lj+1,lk+1)-myphi(li+1,lj-1,lk-1))))

	  grad_fiz=3.0_db*(p1d3q27*(myphi(li,lj,lk+1)-myphi(li,lj,lk-1)) + &
		 p2d3q27*((myphi(li+1,lj,lk+1)-myphi(li-1,lj,lk-1))+ &
		 (myphi(li-1,lj,lk+1)-myphi(li+1,lj,lk-1))+ &
		 (myphi(li,lj+1,lk+1)-myphi(li,lj-1,lk-1))+ &
		 (myphi(li,lj-1,lk+1)-myphi(li,lj+1,lk-1)) ) + &
		 p3d3q27*((myphi(li+1,lj+1,lk+1)-myphi(li-1,lj-1,lk-1)) &
		 +(myphi(li-1,lj-1,lk+1)-myphi(li+1,lj+1,lk-1))+ &
		 (myphi(li+1,lj-1,lk+1)-myphi(li-1,lj+1,lk-1))+ &
		 (myphi(li-1,lj+1,lk+1)-myphi(li+1,lj-1,lk-1))))
      
	  mod_grad=sqrt(grad_fix**TWO+grad_fiy**TWO+grad_fiz**TWO)
#ifdef CSF
	  normal_denom=sqrt(mod_grad*mod_grad+1.0e-24_db)
#else
	  normal_denom=mod_grad+1.0e-9_db
#endif

	  auxfields_s(ii,jj,kk,1,myblock)= &
	   real(grad_fix/normal_denom,kind=strdb)
	  auxfields_s(ii,jj,kk,2,myblock)= &
	   real(grad_fiy/normal_denom,kind=strdb)
	  auxfields_s(ii,jj,kk,3,myblock)= &
	   real(grad_fiz/normal_denom,kind=strdb)
	  
	  auxfields_s(ii,jj,kk,4,myblock)=mod_grad 

	  auxfields_s(ii,jj,kk,5,myblock)=real( &
	   myphi(li,lj,lk)*(1.0_db-myphi(li,lj,lk))*(grad_fix/normal_denom),kind=strdb)
	  auxfields_s(ii,jj,kk,6,myblock)=real( &
	   myphi(li,lj,lk)*(1.0_db-myphi(li,lj,lk))*(grad_fiy/normal_denom),kind=strdb)
	  auxfields_s(ii,jj,kk,7,myblock)=real( &
	   myphi(li,lj,lk)*(1.0_db-myphi(li,lj,lk))*(grad_fiz/normal_denom),kind=strdb)
	   
      !lap_phi here
      locauxfields_s(ii,jj,kk,1,myblock)=real( &
                   (2.0_db*invcssq)*(myphi(li,lj,lk)*(p0d3q27-1.0_db) + &
                   ( p1d3q27*(myphi(li+1,lj,lk)+myphi(li-1,lj,lk) + &
                   myphi(li,lj+1,lk)+myphi(li,lj-1,lk) + &
                   myphi(li,lj,lk+1)+myphi(li,lj,lk-1)) + &
                   p2d3q27*( (myphi(li+1,lj+1,lk)+myphi(li-1,lj-1,lk))+ &
                   (myphi(li+1,lj-1,lk)+myphi(li-1,lj+1,lk))+ &
                   (myphi(li+1,lj,lk+1)+myphi(li-1,lj,lk-1))+ &
                   (myphi(li+1,lj,lk-1)+myphi(li-1,lj,lk+1)) + &
                   (myphi(li,lj+1,lk+1)+myphi(li,lj-1,lk-1))+ &
                   (myphi(li,lj+1,lk-1)+myphi(li,lj-1,lk+1)) )  + &
                   p3d3q27*((myphi(li+1,lj+1,lk+1)+myphi(li-1,lj-1,lk-1))+ &
                   (myphi(li+1,lj-1,lk-1)+myphi(li-1,lj+1,lk+1))+ &
                   (myphi(li+1,lj-1,lk+1)+myphi(li-1,lj+1,lk-1))+ &
                   (myphi(li+1,lj+1,lk-1)+myphi(li-1,lj-1,lk+1))))),kind=strdb)
      
      return
      
   endsubroutine compute_norm_interface_kernel_yminus
   
   ! _zplus: last z tile layer (zblock=nzblock_d-1), all x and y tiles; grid (nxb,nyb,1)
   attributes(global) subroutine compute_norm_interface_kernel_zplus(step, &
      iprobe,jprobe,kprobe,flop,nx,ny,nz,coords,isfluid,rho_r,rho_b, &
      ntotphifields,ntotauxfields,ntotlocauxfields,phifields_s,auxfields_s,locauxfields_s)

      implicit none
      integer :: step,iprobe,jprobe,kprobe,flop,nx,ny,nz
      integer, dimension(3) :: coords
      integer(kind=isf), dimension(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff) :: isfluid
      real(kind=db) :: rho_r,rho_b
      integer :: ntotphifields,ntotauxfields,ntotlocauxfields
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields,nblocks_d) :: phifields_s
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nauxfields,nblocks_d) :: auxfields_s
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nlocauxfields,nblocks_d) :: locauxfields_s
      
      real(kind=db), shared :: myphi(0:TILE_DIMx+1,0:TILE_DIMy+1,0:TILE_DIMz+1)
      real(kind=db):: grad_fix,grad_fiy,grad_fiz,mod_grad,normal_denom
      
      integer :: i,j,k,myblock,intblock
      integer :: ii,jj,kk
      integer :: li,lj,lk
      integer :: xblock,yblock,zblock
      integer, parameter :: nhx=TILE_DIMx+2, nhy=TILE_DIMy+2, nhz=TILE_DIMz+2
      integer, parameter :: nthreads=TILE_DIMx*TILE_DIMy*TILE_DIMz
      integer :: ithread,ih
!      integer :: gi,gj,gk
      
            
      ! one thread per interior node; the threads first load the tile and its halo into
      ! shared memory, nodes ih=ithread, ithread+nthreads, ... in x-fastest order
      ithread=(threadIdx%x-1)+TILE_DIMx*(threadIdx%y-1)+TILE_DIMx*TILE_DIMy*(threadIdx%z-1)
      do ih=ithread,nhx*nhy*nhz-1,nthreads
        lk=ih/(nhx*nhy)
        lj=(ih-lk*nhx*nhy)/nhx
        li=ih-lk*nhx*nhy-lj*nhx
        i = (blockIdx%x-1) * TILE_DIMx + li
        j = (blockIdx%y-1) * TILE_DIMy + lj
        k = ((nzblock_d-2)-1) * TILE_DIMz + lk
        xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
        yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
        zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
        myblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
        ii=i-xblock*TILE_DIMx+2*TILE_DIMx
        jj=j-yblock*TILE_DIMy+2*TILE_DIMy
        kk=k-zblock*TILE_DIMz+2*TILE_DIMz
        myphi(li,lj,lk)=real(phifields_s(ii,jj,kk,1,myblock),kind=db)
      enddo
      call syncthreads

      li=threadIdx%x
      lj=threadIdx%y
      lk=threadIdx%z
      i = (blockIdx%x-1) * TILE_DIMx + li
      j = (blockIdx%y-1) * TILE_DIMy + lj
      k = ((nzblock_d-2)-1) * TILE_DIMz + lk

      if(abs(isfluid(i,j,k)) .ne. 1)return

      !internal-node block of the tile
      myblock=blockIdx%x+blockIdx%y*nxblock_d+(nzblock_d-2)*nxyblock_d+1
      ii=li
      jj=lj
      kk=lk

	  grad_fix=3.0_db*(p1d3q27*(myphi(li+1,lj,lk)-myphi(li-1,lj,lk)) + &
		 p2d3q27*( (myphi(li+1,lj+1,lk)-myphi(li-1,lj-1,lk))+ &
		 (myphi(li+1,lj-1,lk)-myphi(li-1,lj+1,lk))+ &
		 (myphi(li+1,lj,lk+1)-myphi(li-1,lj,lk-1))+ &
		 (myphi(li+1,lj,lk-1)-myphi(li-1,lj,lk+1)) )  + &
		 p3d3q27*((myphi(li+1,lj+1,lk+1)-myphi(li-1,lj-1,lk-1))+ &
		 (myphi(li+1,lj-1,lk-1)-myphi(li-1,lj+1,lk+1))+ &
		 (myphi(li+1,lj-1,lk+1)-myphi(li-1,lj+1,lk-1))+ &
		 (myphi(li+1,lj+1,lk-1)-myphi(li-1,lj-1,lk+1))))

	  grad_fiy=3.0_db*(p1d3q27*(myphi(li,lj+1,lk)-myphi(li,lj-1,lk)) + &
		 p2d3q27*((myphi(li+1,lj+1,lk)-myphi(li-1,lj-1,lk))+ &
		 (myphi(li-1,lj+1,lk)-myphi(li+1,lj-1,lk))+ &
		 (myphi(li,lj+1,lk+1)-myphi(li,lj-1,lk-1))+ &
		 (myphi(li,lj+1,lk-1)-myphi(li,lj-1,lk+1)) ) + &
		 p3d3q27*((myphi(li+1,lj+1,lk+1)-myphi(li-1,lj-1,lk-1))+ &
		 (myphi(li-1,lj+1,lk-1)-myphi(li+1,lj-1,lk+1))+ &
		 (myphi(li+1,lj+1,lk-1)-myphi(li-1,lj-1,lk+1))+ &
		 (myphi(li-1,lj+1,lk+1)-myphi(li+1,lj-1,lk-1))))

	  grad_fiz=3.0_db*(p1d3q27*(myphi(li,lj,lk+1)-myphi(li,lj,lk-1)) + &
		 p2d3q27*((myphi(li+1,lj,lk+1)-myphi(li-1,lj,lk-1))+ &
		 (myphi(li-1,lj,lk+1)-myphi(li+1,lj,lk-1))+ &
		 (myphi(li,lj+1,lk+1)-myphi(li,lj-1,lk-1))+ &
		 (myphi(li,lj-1,lk+1)-myphi(li,lj+1,lk-1)) ) + &
		 p3d3q27*((myphi(li+1,lj+1,lk+1)-myphi(li-1,lj-1,lk-1)) &
		 +(myphi(li-1,lj-1,lk+1)-myphi(li+1,lj+1,lk-1))+ &
		 (myphi(li+1,lj-1,lk+1)-myphi(li-1,lj+1,lk-1))+ &
		 (myphi(li-1,lj+1,lk+1)-myphi(li+1,lj-1,lk-1))))
      
	  mod_grad=sqrt(grad_fix**TWO+grad_fiy**TWO+grad_fiz**TWO)
#ifdef CSF
	  normal_denom=sqrt(mod_grad*mod_grad+1.0e-24_db)
#else
	  normal_denom=mod_grad+1.0e-9_db
#endif

	  auxfields_s(ii,jj,kk,1,myblock)= &
	   real(grad_fix/normal_denom,kind=strdb)
	  auxfields_s(ii,jj,kk,2,myblock)= &
	   real(grad_fiy/normal_denom,kind=strdb)
	  auxfields_s(ii,jj,kk,3,myblock)= &
	   real(grad_fiz/normal_denom,kind=strdb)
	  
	  auxfields_s(ii,jj,kk,4,myblock)=mod_grad 

	  auxfields_s(ii,jj,kk,5,myblock)=real( &
	   myphi(li,lj,lk)*(1.0_db-myphi(li,lj,lk))*(grad_fix/normal_denom),kind=strdb)
	  auxfields_s(ii,jj,kk,6,myblock)=real( &
	   myphi(li,lj,lk)*(1.0_db-myphi(li,lj,lk))*(grad_fiy/normal_denom),kind=strdb)
	  auxfields_s(ii,jj,kk,7,myblock)=real( &
	   myphi(li,lj,lk)*(1.0_db-myphi(li,lj,lk))*(grad_fiz/normal_denom),kind=strdb)
	   
      !lap_phi here
      locauxfields_s(ii,jj,kk,1,myblock)=real( &
                   (2.0_db*invcssq)*(myphi(li,lj,lk)*(p0d3q27-1.0_db) + &
                   ( p1d3q27*(myphi(li+1,lj,lk)+myphi(li-1,lj,lk) + &
                   myphi(li,lj+1,lk)+myphi(li,lj-1,lk) + &
                   myphi(li,lj,lk+1)+myphi(li,lj,lk-1)) + &
                   p2d3q27*( (myphi(li+1,lj+1,lk)+myphi(li-1,lj-1,lk))+ &
                   (myphi(li+1,lj-1,lk)+myphi(li-1,lj+1,lk))+ &
                   (myphi(li+1,lj,lk+1)+myphi(li-1,lj,lk-1))+ &
                   (myphi(li+1,lj,lk-1)+myphi(li-1,lj,lk+1)) + &
                   (myphi(li,lj+1,lk+1)+myphi(li,lj-1,lk-1))+ &
                   (myphi(li,lj+1,lk-1)+myphi(li,lj-1,lk+1)) )  + &
                   p3d3q27*((myphi(li+1,lj+1,lk+1)+myphi(li-1,lj-1,lk-1))+ &
                   (myphi(li+1,lj-1,lk-1)+myphi(li-1,lj+1,lk+1))+ &
                   (myphi(li+1,lj-1,lk+1)+myphi(li-1,lj+1,lk-1))+ &
                   (myphi(li+1,lj+1,lk-1)+myphi(li-1,lj-1,lk+1))))),kind=strdb)
      
      return
      
   endsubroutine compute_norm_interface_kernel_zplus
   
   ! _zminus: first z tile layer (zblock=2), all x and y tiles; grid (nxb,nyb,1)
   attributes(global) subroutine compute_norm_interface_kernel_zminus(step, &
      iprobe,jprobe,kprobe,flop,nx,ny,nz,coords,isfluid,rho_r,rho_b, &
      ntotphifields,ntotauxfields,ntotlocauxfields,phifields_s,auxfields_s,locauxfields_s)

      implicit none
      integer :: step,iprobe,jprobe,kprobe,flop,nx,ny,nz
      integer, dimension(3) :: coords
      integer(kind=isf), dimension(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff) :: isfluid
      real(kind=db) :: rho_r,rho_b
      integer :: ntotphifields,ntotauxfields,ntotlocauxfields
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields,nblocks_d) :: phifields_s
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nauxfields,nblocks_d) :: auxfields_s
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nlocauxfields,nblocks_d) :: locauxfields_s
      
      real(kind=db), shared :: myphi(0:TILE_DIMx+1,0:TILE_DIMy+1,0:TILE_DIMz+1)
      real(kind=db):: grad_fix,grad_fiy,grad_fiz,mod_grad,normal_denom
      
      integer :: i,j,k,myblock,intblock
      integer :: ii,jj,kk
      integer :: li,lj,lk
      integer :: xblock,yblock,zblock
      integer, parameter :: nhx=TILE_DIMx+2, nhy=TILE_DIMy+2, nhz=TILE_DIMz+2
      integer, parameter :: nthreads=TILE_DIMx*TILE_DIMy*TILE_DIMz
      integer :: ithread,ih
!      integer :: gi,gj,gk

      
      ! one thread per interior node; the threads first load the tile and its halo into
      ! shared memory, nodes ih=ithread, ithread+nthreads, ... in x-fastest order
      ithread=(threadIdx%x-1)+TILE_DIMx*(threadIdx%y-1)+TILE_DIMx*TILE_DIMy*(threadIdx%z-1)
      do ih=ithread,nhx*nhy*nhz-1,nthreads
        lk=ih/(nhx*nhy)
        lj=(ih-lk*nhx*nhy)/nhx
        li=ih-lk*nhx*nhy-lj*nhx
        i = (blockIdx%x-1) * TILE_DIMx + li
        j = (blockIdx%y-1) * TILE_DIMy + lj
        k = lk
        xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
        yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
        zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
        myblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
        ii=i-xblock*TILE_DIMx+2*TILE_DIMx
        jj=j-yblock*TILE_DIMy+2*TILE_DIMy
        kk=k-zblock*TILE_DIMz+2*TILE_DIMz
        myphi(li,lj,lk)=real(phifields_s(ii,jj,kk,1,myblock),kind=db)
      enddo
      call syncthreads

      li=threadIdx%x
      lj=threadIdx%y
      lk=threadIdx%z
      i = (blockIdx%x-1) * TILE_DIMx + li
      j = (blockIdx%y-1) * TILE_DIMy + lj
      k = lk

      if(abs(isfluid(i,j,k)) .ne. 1)return

      !internal-node block of the tile
      myblock=blockIdx%x+blockIdx%y*nxblock_d+1*nxyblock_d+1
      ii=li
      jj=lj
      kk=lk

	  grad_fix=3.0_db*(p1d3q27*(myphi(li+1,lj,lk)-myphi(li-1,lj,lk)) + &
		 p2d3q27*( (myphi(li+1,lj+1,lk)-myphi(li-1,lj-1,lk))+ &
		 (myphi(li+1,lj-1,lk)-myphi(li-1,lj+1,lk))+ &
		 (myphi(li+1,lj,lk+1)-myphi(li-1,lj,lk-1))+ &
		 (myphi(li+1,lj,lk-1)-myphi(li-1,lj,lk+1)) )  + &
		 p3d3q27*((myphi(li+1,lj+1,lk+1)-myphi(li-1,lj-1,lk-1))+ &
		 (myphi(li+1,lj-1,lk-1)-myphi(li-1,lj+1,lk+1))+ &
		 (myphi(li+1,lj-1,lk+1)-myphi(li-1,lj+1,lk-1))+ &
		 (myphi(li+1,lj+1,lk-1)-myphi(li-1,lj-1,lk+1))))

	  grad_fiy=3.0_db*(p1d3q27*(myphi(li,lj+1,lk)-myphi(li,lj-1,lk)) + &
		 p2d3q27*((myphi(li+1,lj+1,lk)-myphi(li-1,lj-1,lk))+ &
		 (myphi(li-1,lj+1,lk)-myphi(li+1,lj-1,lk))+ &
		 (myphi(li,lj+1,lk+1)-myphi(li,lj-1,lk-1))+ &
		 (myphi(li,lj+1,lk-1)-myphi(li,lj-1,lk+1)) ) + &
		 p3d3q27*((myphi(li+1,lj+1,lk+1)-myphi(li-1,lj-1,lk-1))+ &
		 (myphi(li-1,lj+1,lk-1)-myphi(li+1,lj-1,lk+1))+ &
		 (myphi(li+1,lj+1,lk-1)-myphi(li-1,lj-1,lk+1))+ &
		 (myphi(li-1,lj+1,lk+1)-myphi(li+1,lj-1,lk-1))))

	  grad_fiz=3.0_db*(p1d3q27*(myphi(li,lj,lk+1)-myphi(li,lj,lk-1)) + &
		 p2d3q27*((myphi(li+1,lj,lk+1)-myphi(li-1,lj,lk-1))+ &
		 (myphi(li-1,lj,lk+1)-myphi(li+1,lj,lk-1))+ &
		 (myphi(li,lj+1,lk+1)-myphi(li,lj-1,lk-1))+ &
		 (myphi(li,lj-1,lk+1)-myphi(li,lj+1,lk-1)) ) + &
		 p3d3q27*((myphi(li+1,lj+1,lk+1)-myphi(li-1,lj-1,lk-1)) &
		 +(myphi(li-1,lj-1,lk+1)-myphi(li+1,lj+1,lk-1))+ &
		 (myphi(li+1,lj-1,lk+1)-myphi(li-1,lj+1,lk-1))+ &
		 (myphi(li-1,lj+1,lk+1)-myphi(li+1,lj-1,lk-1))))
      
	  mod_grad=sqrt(grad_fix**TWO+grad_fiy**TWO+grad_fiz**TWO)
#ifdef CSF
	  normal_denom=sqrt(mod_grad*mod_grad+1.0e-24_db)
#else
	  normal_denom=mod_grad+1.0e-9_db
#endif

	  auxfields_s(ii,jj,kk,1,myblock)= &
	   real(grad_fix/normal_denom,kind=strdb)
	  auxfields_s(ii,jj,kk,2,myblock)= &
	   real(grad_fiy/normal_denom,kind=strdb)
	  auxfields_s(ii,jj,kk,3,myblock)= &
	   real(grad_fiz/normal_denom,kind=strdb)
	  
	  auxfields_s(ii,jj,kk,4,myblock)=mod_grad 

	  auxfields_s(ii,jj,kk,5,myblock)=real( &
	   myphi(li,lj,lk)*(1.0_db-myphi(li,lj,lk))*(grad_fix/normal_denom),kind=strdb)
	  auxfields_s(ii,jj,kk,6,myblock)=real( &
	   myphi(li,lj,lk)*(1.0_db-myphi(li,lj,lk))*(grad_fiy/normal_denom),kind=strdb)
	  auxfields_s(ii,jj,kk,7,myblock)=real( &
	   myphi(li,lj,lk)*(1.0_db-myphi(li,lj,lk))*(grad_fiz/normal_denom),kind=strdb)
	   
      !lap_phi here
      locauxfields_s(ii,jj,kk,1,myblock)=real( &
                   (2.0_db*invcssq)*(myphi(li,lj,lk)*(p0d3q27-1.0_db) + &
                   ( p1d3q27*(myphi(li+1,lj,lk)+myphi(li-1,lj,lk) + &
                   myphi(li,lj+1,lk)+myphi(li,lj-1,lk) + &
                   myphi(li,lj,lk+1)+myphi(li,lj,lk-1)) + &
                   p2d3q27*( (myphi(li+1,lj+1,lk)+myphi(li-1,lj-1,lk))+ &
                   (myphi(li+1,lj-1,lk)+myphi(li-1,lj+1,lk))+ &
                   (myphi(li+1,lj,lk+1)+myphi(li-1,lj,lk-1))+ &
                   (myphi(li+1,lj,lk-1)+myphi(li-1,lj,lk+1)) + &
                   (myphi(li,lj+1,lk+1)+myphi(li,lj-1,lk-1))+ &
                   (myphi(li,lj+1,lk-1)+myphi(li,lj-1,lk+1)) )  + &
                   p3d3q27*((myphi(li+1,lj+1,lk+1)+myphi(li-1,lj-1,lk-1))+ &
                   (myphi(li+1,lj-1,lk-1)+myphi(li-1,lj+1,lk+1))+ &
                   (myphi(li+1,lj-1,lk+1)+myphi(li-1,lj+1,lk-1))+ &
                   (myphi(li+1,lj+1,lk-1)+myphi(li-1,lj-1,lk+1))))),kind=strdb)
      
      return
      
   endsubroutine compute_norm_interface_kernel_zminus
  
   ! compute_div_theta_n_kernel: div theta on all the tiles, with the shared-memory tiling of
   ! compute_norm_interface_kernel (three shared arrays, one per component of theta). It reads
   ! auxfields 5-7 on the node and its 26 neighbours (the auxfields halo exchange must be
   ! complete) and writes, for fluid nodes, the isotropic D3Q27 divergence
   !   locauxfields 2 = div theta = (1/cs^2) sum_i w_i theta(x+c_i).c_i,
   ! used as -sharp_c*div theta (sharp_c = 4*tau_diff/width) in the phi update and moment kernels.
#if defined(GPUTILEX) && defined(GPUTILEY) && defined(GPUTILEZ) && (GPUTILEX*GPUTILEY*GPUTILEZ) == 512
   ! four blocks per multiprocessor: 32 registers, no spills in double precision
   attributes(global) launch_bounds(512,4) subroutine compute_div_theta_n_kernel(flop,nx,ny,nz,coords,isfluid, &
#else
   attributes(global) subroutine compute_div_theta_n_kernel(flop,nx,ny,nz,coords,isfluid, &
#endif
     rho_r,rho_b,ntotphifields,ntotauxfields,ntotlocauxfields,phifields_s,auxfields_s,locauxfields_s)

      implicit none
      integer :: flop,nx,ny,nz
      integer, dimension(3) :: coords
      integer(kind=isf), dimension(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff) :: isfluid
      real(kind=db) :: rho_r,rho_b
      integer :: ntotphifields,ntotauxfields,ntotlocauxfields
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields,nblocks_d) :: phifields_s
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nauxfields,nblocks_d) :: auxfields_s
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nlocauxfields,nblocks_d) :: locauxfields_s
      
      real(kind=db), shared :: myarrx(0:TILE_DIMx+1,0:TILE_DIMy+1,0:TILE_DIMz+1)
      real(kind=db), shared :: myarry(0:TILE_DIMx+1,0:TILE_DIMy+1,0:TILE_DIMz+1)
      real(kind=db), shared :: myarrz(0:TILE_DIMx+1,0:TILE_DIMy+1,0:TILE_DIMz+1)
      
      integer :: i,j,k,gi,gj,gk,myblock,idblock
      integer :: ii,jj,kk
      integer :: li,lj,lk
      integer :: xblock,yblock,zblock
      integer, parameter :: nhx=TILE_DIMx+2, nhy=TILE_DIMy+2, nhz=TILE_DIMz+2
      integer, parameter :: nthreads=TILE_DIMx*TILE_DIMy*TILE_DIMz
      integer :: ithread,ih
      
      ! one thread per interior node; the threads first load the tile and its halo into
      ! shared memory, nodes ih=ithread, ithread+nthreads, ... in x-fastest order
      ithread=(threadIdx%x-1)+TILE_DIMx*(threadIdx%y-1)+TILE_DIMx*TILE_DIMy*(threadIdx%z-1)
      do ih=ithread,nhx*nhy*nhz-1,nthreads
        lk=ih/(nhx*nhy)
        lj=(ih-lk*nhx*nhy)/nhx
        li=ih-lk*nhx*nhy-lj*nhx
        i = (blockIdx%x-1) * TILE_DIMx + li
        j = (blockIdx%y-1) * TILE_DIMy + lj
        k = (blockIdx%z-1) * TILE_DIMz + lk
        xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
        yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
        zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
        myblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
        ii=i-xblock*TILE_DIMx+2*TILE_DIMx
        jj=j-yblock*TILE_DIMy+2*TILE_DIMy
        kk=k-zblock*TILE_DIMz+2*TILE_DIMz
        myarrx(li,lj,lk)=real(auxfields_s(ii,jj,kk,5,myblock),kind=db)
        myarry(li,lj,lk)=real(auxfields_s(ii,jj,kk,6,myblock),kind=db)
        myarrz(li,lj,lk)=real(auxfields_s(ii,jj,kk,7,myblock),kind=db)
      enddo
      call syncthreads

      li=threadIdx%x
      lj=threadIdx%y
      lk=threadIdx%z
      i = (blockIdx%x-1) * TILE_DIMx + li
      j = (blockIdx%y-1) * TILE_DIMy + lj
      k = (blockIdx%z-1) * TILE_DIMz + lk

      if(abs(isfluid(i,j,k)) .ne. 1)return

      !internal-node block of the tile
      myblock=blockIdx%x+blockIdx%y*nxblock_d+blockIdx%z*nxyblock_d+1
      ii=li
      jj=lj
      kk=lk


	   
      !div_thetan here
      locauxfields_s(ii,jj,kk,2,myblock)=real( &
       (( p1d3q27*(myarrx(li+1,lj,lk)-myarrx(li-1,lj,lk)) + &
       p2d3q27*((myarrx(li+1,lj+1,lk)-myarrx(li-1,lj-1,lk))+(myarrx(li+1,lj-1,lk)-myarrx(li-1,lj+1,lk))+ &
       (myarrx(li+1,lj,lk+1)-myarrx(li-1,lj,lk-1))+(myarrx(li+1,lj,lk-1)-myarrx(li-1,lj,lk+1)))+ &
       p3d3q27*((myarrx(li+1,lj+1,lk+1)-myarrx(li-1,lj-1,lk-1))+(myarrx(li+1,lj-1,lk-1)-myarrx(li-1,lj+1,lk+1))+ &
       (myarrx(li+1,lj-1,lk+1)-myarrx(li-1,lj+1,lk-1))+(myarrx(li+1,lj+1,lk-1)-myarrx(li-1,lj-1,lk+1))))+ &
       (p1d3q27*(myarry(li,lj+1,lk)-myarry(li,lj-1,lk)) + &
       p2d3q27*((myarry(li+1,lj+1,lk)-myarry(li-1,lj-1,lk))+(myarry(li-1,lj+1,lk)-myarry(li+1,lj-1,lk))+ &
       (myarry(li,lj+1,lk+1)-myarry(li,lj-1,lk-1))+(myarry(li,lj+1,lk-1)-myarry(li,lj-1,lk+1)))+ &
       p3d3q27*((myarry(li+1,lj+1,lk+1)-myarry(li-1,lj-1,lk-1))+(myarry(li-1,lj+1,lk-1)-myarry(li+1,lj-1,lk+1))+ &
       (myarry(li+1,lj+1,lk-1)-myarry(li-1,lj-1,lk+1))+(myarry(li-1,lj+1,lk+1)-myarry(li+1,lj-1,lk-1))))+ &
       (p1d3q27*(myarrz(li,lj,lk+1)-myarrz(li,lj,lk-1)) + &
       p2d3q27*((myarrz(li+1,lj,lk+1)-myarrz(li-1,lj,lk-1))+(myarrz(li-1,lj,lk+1)-myarrz(li+1,lj,lk-1))+ &
       (myarrz(li,lj+1,lk+1)-myarrz(li,lj-1,lk-1))+(myarrz(li,lj-1,lk+1)-myarrz(li,lj+1,lk-1)))+ &
       p3d3q27*((myarrz(li+1,lj+1,lk+1)-myarrz(li-1,lj-1,lk-1))+(myarrz(li-1,lj-1,lk+1)-myarrz(li+1,lj+1,lk-1))+ &
       (myarrz(li+1,lj-1,lk+1)-myarrz(li-1,lj+1,lk-1))+(myarrz(li-1,lj+1,lk+1)-myarrz(li+1,lj-1,lk-1)))))*invcssq &
       ,kind=strdb)
  
      return
      
   endsubroutine compute_div_theta_n_kernel

#ifdef CSF
   ! compute_csf_force_kernel (CSF): capillary force of the continuum-surface-force model on all
   ! the tiles, with the shared-memory tiling of compute_norm_interface_kernel. It reads n
   ! (auxfields 1-3) on the node and its 26 neighbours and |grad phi| (auxfields 4) on the node,
   ! and writes, for fluid nodes, forces 1-3 = sigma*kappa*n*|grad phi| (= sigma*kappa*grad phi up
   ! to the regularization of n), with the curvature kappa = -div n from the isotropic D3Q27
   ! divergence. kappa is not stored; the moment kernels add forces 1-3 to the hydrodynamic force.
   ! With phi = 1 inside a droplet n points inward and kappa > 0 (docs/CSF_LAPLACE_VALIDATION.md).
#if defined(GPUTILEX) && defined(GPUTILEY) && defined(GPUTILEZ) && (GPUTILEX*GPUTILEY*GPUTILEZ) == 512
   ! four blocks per multiprocessor: 32 registers, no spills in double precision
   attributes(global) launch_bounds(512,4) subroutine compute_csf_force_kernel(nx,ny,nz,isfluid,sigma, &
#else
   attributes(global) subroutine compute_csf_force_kernel(nx,ny,nz,isfluid,sigma, &
#endif
     ntotauxfields,ntotforces,auxfields_s,forces_s)

      implicit none
      integer :: nx,ny,nz,ntotauxfields,ntotforces
      real(kind=db) :: sigma
      integer(kind=isf), dimension(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff) :: isfluid
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nauxfields,nblocks_d) :: auxfields_s
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nforces,nblocks_d) :: forces_s

      real(kind=db), shared :: normal_x(0:TILE_DIMx+1,0:TILE_DIMy+1,0:TILE_DIMz+1)
      real(kind=db), shared :: normal_y(0:TILE_DIMx+1,0:TILE_DIMy+1,0:TILE_DIMz+1)
      real(kind=db), shared :: normal_z(0:TILE_DIMx+1,0:TILE_DIMy+1,0:TILE_DIMz+1)
      real(kind=db) :: div_normal,curvature,mod_grad
      integer :: i,j,k,myblock,idblock,ii,jj,kk,li,lj,lk,xblock,yblock,zblock
      integer, parameter :: nhx=TILE_DIMx+2, nhy=TILE_DIMy+2, nhz=TILE_DIMz+2
      integer, parameter :: nthreads=TILE_DIMx*TILE_DIMy*TILE_DIMz
      integer :: ithread,ih

      ! one thread per interior node; the threads first load the tile and its halo into
      ! shared memory, nodes ih=ithread, ithread+nthreads, ... in x-fastest order
      ithread=(threadIdx%x-1)+TILE_DIMx*(threadIdx%y-1)+TILE_DIMx*TILE_DIMy*(threadIdx%z-1)
      do ih=ithread,nhx*nhy*nhz-1,nthreads
        lk=ih/(nhx*nhy)
        lj=(ih-lk*nhx*nhy)/nhx
        li=ih-lk*nhx*nhy-lj*nhx
        i = (blockIdx%x-1)*TILE_DIMx+li
        j = (blockIdx%y-1)*TILE_DIMy+lj
        k = (blockIdx%z-1)*TILE_DIMz+lk
        xblock=(i+2*TILE_DIMx-1)/TILE_DIMx
        yblock=(j+2*TILE_DIMy-1)/TILE_DIMy
        zblock=(k+2*TILE_DIMz-1)/TILE_DIMz
        myblock=(xblock-1)+(yblock-1)*nxblock_d+(zblock-1)*nxyblock_d+1
        ii=i-xblock*TILE_DIMx+2*TILE_DIMx
        jj=j-yblock*TILE_DIMy+2*TILE_DIMy
        kk=k-zblock*TILE_DIMz+2*TILE_DIMz
        normal_x(li,lj,lk)=real(auxfields_s(ii,jj,kk,1,myblock),kind=db)
        normal_y(li,lj,lk)=real(auxfields_s(ii,jj,kk,2,myblock),kind=db)
        normal_z(li,lj,lk)=real(auxfields_s(ii,jj,kk,3,myblock),kind=db)
      enddo
      call syncthreads

      li=threadIdx%x
      lj=threadIdx%y
      lk=threadIdx%z
      i = (blockIdx%x-1)*TILE_DIMx+li
      j = (blockIdx%y-1)*TILE_DIMy+lj
      k = (blockIdx%z-1)*TILE_DIMz+lk

      if(abs(isfluid(i,j,k)) .ne. 1)return

      !internal-node block of the tile
      myblock=blockIdx%x+blockIdx%y*nxblock_d+blockIdx%z*nxyblock_d+1
      ii=li
      jj=lj
      kk=lk

      ! div_normal = sum_i w_i n(x+c_i).c_i, so that div n = div_normal/cs^2 (invcssq below)
      div_normal= &
       p1d3q27*(normal_x(li+1,lj,lk)-normal_x(li-1,lj,lk) + &
                 normal_y(li,lj+1,lk)-normal_y(li,lj-1,lk) + &
                 normal_z(li,lj,lk+1)-normal_z(li,lj,lk-1)) + &
       p2d3q27*( &
        normal_x(li+1,lj+1,lk)-normal_x(li-1,lj-1,lk) + &
        normal_x(li+1,lj-1,lk)-normal_x(li-1,lj+1,lk) + &
        normal_x(li+1,lj,lk+1)-normal_x(li-1,lj,lk-1) + &
        normal_x(li+1,lj,lk-1)-normal_x(li-1,lj,lk+1) + &
        normal_y(li+1,lj+1,lk)-normal_y(li-1,lj-1,lk) + &
        normal_y(li-1,lj+1,lk)-normal_y(li+1,lj-1,lk) + &
        normal_y(li,lj+1,lk+1)-normal_y(li,lj-1,lk-1) + &
        normal_y(li,lj+1,lk-1)-normal_y(li,lj-1,lk+1) + &
        normal_z(li+1,lj,lk+1)-normal_z(li-1,lj,lk-1) + &
        normal_z(li-1,lj,lk+1)-normal_z(li+1,lj,lk-1) + &
        normal_z(li,lj+1,lk+1)-normal_z(li,lj-1,lk-1) + &
        normal_z(li,lj-1,lk+1)-normal_z(li,lj+1,lk-1)) + &
       p3d3q27*( &
        normal_x(li+1,lj+1,lk+1)-normal_x(li-1,lj-1,lk-1) + &
        normal_x(li+1,lj-1,lk-1)-normal_x(li-1,lj+1,lk+1) + &
        normal_x(li+1,lj-1,lk+1)-normal_x(li-1,lj+1,lk-1) + &
        normal_x(li+1,lj+1,lk-1)-normal_x(li-1,lj-1,lk+1) + &
        normal_y(li+1,lj+1,lk+1)-normal_y(li-1,lj-1,lk-1) + &
        normal_y(li-1,lj+1,lk-1)-normal_y(li+1,lj-1,lk+1) + &
        normal_y(li+1,lj+1,lk-1)-normal_y(li-1,lj-1,lk+1) + &
        normal_y(li-1,lj+1,lk+1)-normal_y(li+1,lj-1,lk-1) + &
        normal_z(li+1,lj+1,lk+1)-normal_z(li-1,lj-1,lk-1) + &
        normal_z(li-1,lj-1,lk+1)-normal_z(li+1,lj+1,lk-1) + &
        normal_z(li+1,lj-1,lk+1)-normal_z(li-1,lj+1,lk-1) + &
        normal_z(li-1,lj+1,lk+1)-normal_z(li+1,lj-1,lk-1))

      curvature=-invcssq*div_normal
      mod_grad=real(auxfields_s(ii,jj,kk,4,myblock),kind=db)
      forces_s(ii,jj,kk,1,myblock)=real(sigma*curvature*normal_x(li,lj,lk)*mod_grad,kind=strdb)
      forces_s(ii,jj,kk,2,myblock)=real(sigma*curvature*normal_y(li,lj,lk)*mod_grad,kind=strdb)
      forces_s(ii,jj,kk,3,myblock)=real(sigma*curvature*normal_z(li,lj,lk)*mod_grad,kind=strdb)
   endsubroutine compute_csf_force_kernel
#endif
#endif 
endmodule lb_cuda_auxfields
