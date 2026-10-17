#include "defines.h"
module nvml_interface
  ! Fortran interfaces to the NVML helpers of nvml_wrapper.c, available only when _NVML is
  ! defined (Makefile targets nvfortran-nvml and nvfortran-nvml-mpi, linked with -lnvidia-ml),
  ! and the conversion of two readings of the GPU energy counter into joules.
  use iso_c_binding
  implicit none
#ifdef _NVML
  ! Power draw of device idx in milliwatts, or a negative error code. With GETPOWER,
  ! fused_LB_cuda samples it every stamp_term steps for the device mydev_c of the rank.
  interface
    function get_gpu_power_index(idx) bind(C, name="get_gpu_power_index")
      import :: c_int
      integer(c_int), value :: idx
      integer(c_int) :: get_gpu_power_index
    end function
  end interface
  
  ! Energy consumed by device idx since the driver was loaded, in millijoules, stored in e_mJ;
  ! returns 0, or a negative error code. NVML returns an unsigned 64-bit value, held here in a
  ! signed c_long_long, so values of 2**63 mJ or more appear negative (see u64_delta_J).
  ! Without GETPOWER, integrator_module reads it before and after the time loop.
  interface
    integer(c_int) function get_gpu_energy_mJ_u64(idx, e_mJ) bind(C, name="get_gpu_energy_mJ_u64")
      import :: c_int, c_long_long
      integer(c_int), value              :: idx
      integer(c_long_long)      :: e_mJ   ! output (passed by reference)
    end function get_gpu_energy_mJ_u64
  end interface  
  
  contains
  
  ! Energy in joules between two readings e0 and e1 of the NVML energy counter. Both readings
  ! are unsigned 64-bit integers stored in signed variables (Fortran has no unsigned type), so
  ! the difference e1 - e0 is computed modulo 2**64 from the two 32-bit halves. The split assumes
  ! little-endian byte order, in which element 1 of the transferred array is the low word.
  pure function u64_delta_J(e0, e1) result(dJ)
    implicit none
    integer(c_long_long), intent(in) :: e0, e1
    real(kind=PRC) :: dJ

    integer(c_int32_t) :: a0(2), a1(2)      ! [lo, hi] as int32 (bitwise)
    integer(c_long_long) :: lo0, hi0, lo1, hi1
    integer(c_long_long) :: dlo, dhi, borrow

    ! Split into two 32-bit words (transfer preserves the bits)
    a0 = transfer(e0, a0)
    a1 = transfer(e1, a1)

    ! Interpret each word as unsigned 32-bit: add 2^32 if negative
    lo0 = int(a0(1), c_long_long); if (lo0 < 0_c_long_long) lo0 = lo0 + 4294967296_c_long_long
    hi0 = int(a0(2), c_long_long); if (hi0 < 0_c_long_long) hi0 = hi0 + 4294967296_c_long_long
    lo1 = int(a1(1), c_long_long); if (lo1 < 0_c_long_long) lo1 = lo1 + 4294967296_c_long_long
    hi1 = int(a1(2), c_long_long); if (hi1 < 0_c_long_long) hi1 = hi1 + 4294967296_c_long_long

    ! Unsigned 64-bit subtraction: (hi1:lo1) - (hi0:lo0) modulo 2^64
    if (lo1 >= lo0) then
      dlo = lo1 - lo0
      borrow = 0_c_long_long
    else
      dlo = (lo1 + 4294967296_c_long_long) - lo0
      borrow = 1_c_long_long
    end if

    dhi = hi1 - hi0 - borrow
    if (dhi < 0_c_long_long) dhi = dhi + 4294967296_c_long_long

    ! de_mJ = dhi*2^32 + dlo, then mJ->J
    dJ = real(1.0d-3 * ( real(dhi, c_double) * 4294967296.0d0 + real(dlo, c_double) ),kind=PRC)
  end function u64_delta_J
  
#endif
end module nvml_interface

