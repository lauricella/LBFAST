#include "defines.h"
#if !defined(_OPENACC)  
#error "To use this module the macros _OPENACC should be defined."
#endif


module lb_cuda_vars

   !use iso_c_binding
   use cudafor

   implicit none
   
   integer :: istat

   integer, constant :: TILE_DIMx_d,TILE_DIMy_d,TILE_DIMz_d,TILE_DIM_d
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
   
   contains
   
  attributes(device) elemental function idx5d(ind1,ind2,ind3,ind4,ind5,m1,m2,m3,m4)
 
  implicit none
  
  integer, intent(in) :: ind1,ind2,ind3,ind4,ind5,m1,m2,m3,m4
  
  integer :: idx5d
  
  idx5d=1+(ind1-1)+(ind2-1)*m1+(ind3-1)*(m1*m2)+(ind4-1)*(m1*m2*m3)+ &
   (ind5-1)*(m1*m2*m3*m4)
  
  return
  
 end function idx5d
   

endmodule lb_cuda_vars
