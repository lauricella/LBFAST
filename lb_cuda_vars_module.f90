#include "defines.h"
#if !defined(_OPENACC)  
#error "To use this module the macros _OPENACC should be defined."
#endif


! Launch configuration and device constants shared by the CUDA modules; they are set in
! setup_cuda (lb_cuda_driver_module.f90).
module lb_cuda_vars

   !use iso_c_binding
   use cudafor

   implicit none
   
   integer :: istat
   
   ! Sizes in constant memory: local and global lattice sizes and numbers of array blocks
   ! (nxblock_d = nx/TILE_DIMx + 2, nxyblock_d = nxblock_d*nyblock_d, nblocks_d), used by the
   ! kernels for the tiled indexing.
   integer, constant :: nx_d,ny_d,nz_d
   integer, constant :: lx_d,ly_d,lz_d
   integer, constant :: nxblock_d,nyblock_d,nzblock_d
   integer, constant :: nxyblock_d,nblocks_d
   ! Launch configurations (see the header of lb_cuda_driver_module.f90): dimBlock, one thread per
   ! tile node; dimGrid, all tiles; dimGridInt, interior tiles; dimGridx/y/z, face tiles (ASYNCMPI).
   ! dimGridhalo (all tiles plus the halo tiles), dimBlockhalo and dimBlockshared (a tile plus its
   ! one-node halo) are used only by the test kernels; dimBlock2 is set but not used in any launch.
   type (dim3) :: dimGrid,dimBlock
   type (dim3) :: dimGridhalo,dimBlockhalo
   type (dim3) :: dimBlockshared
   type (dim3) :: dimGridInt
   type (dim3) :: dimGridx,dimGridy,dimGridz
   type (dim3) :: dimBlock2
   ! Device properties and shared memory per block, read in setup_cuda.
   type(cudaDeviceProp) :: prop
   integer(8) :: mshared
   
   ! ldodimGridInt, ldodimGridx, ldodimGridy: false when the corresponding grid has no tiles.
   ! ldodimGrid*_twice: false with a single tile along the direction, so that only the minus face
   ! kernel is launched.
   logical, save :: ldodimGridInt,ldodimGridx,ldodimGridy
   logical, save :: ldodimGridx_twice,ldodimGridy_twice,ldodimGridz_twice
   
   contains
   
   

endmodule lb_cuda_vars
