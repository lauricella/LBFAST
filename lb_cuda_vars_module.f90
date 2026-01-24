#include "defines.h"
#if !defined(_OPENACC)  
#error "To use this module the macros _OPENACC should be defined."
#endif


module lb_cuda_vars

   !use iso_c_binding
   use cudafor

   implicit none
   
   integer :: istat
   
   integer, constant :: nx_d,ny_d,nz_d
   integer, constant :: lx_d,ly_d,lz_d
   integer, constant :: nxblock_d,nyblock_d,nzblock_d
   integer, constant :: nxyblock_d,nblocks_d
   type (dim3) :: dimGrid,dimBlock
   type (dim3) :: dimGridhalo,dimBlockhalo
   type (dim3) :: dimBlockshared
   type (dim3) :: dimGridInt
   type (dim3) :: dimGridx,dimGridy,dimGridz
   type (dim3) :: dimBlock2
   type(cudaDeviceProp) :: prop
   integer(8) :: mshared
   
   logical, save :: ldodimGridInt,ldodimGridx,ldodimGridy
   
   contains
   
   

endmodule lb_cuda_vars
