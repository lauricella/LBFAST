#include "defines.h"
module nvml_interface
  use iso_c_binding
  implicit none
#ifdef MONITORENERGY
  interface
    function get_gpu_power_index(idx) bind(C, name="get_gpu_power_index")
      import :: c_int
      integer(c_int), value :: idx
      integer(c_int) :: get_gpu_power_index
    end function
  end interface
  
  interface
    integer(c_int) function get_gpu_energy_mJ_u64(idx, e_mJ) bind(C, name="get_gpu_energy_mJ_u64")
      import :: c_int, c_long_long
      integer(c_int), value              :: idx
      integer(c_long_long)      :: e_mJ   ! output (passato per riferimento)
    end function get_gpu_energy_mJ_u64
  end interface  
  
  contains
  
  pure function u64_delta_J(e0, e1) result(dJ)
    integer(c_long_long), intent(in) :: e0, e1
    real(kind=PRC) :: dJ

    integer(c_int32_t) :: a0(2), a1(2)      ! [lo, hi] come int32 (bitwise)
    integer(c_long_long) :: lo0, hi0, lo1, hi1
    integer(c_long_long) :: dlo, dhi, borrow

    ! Spezza in due parole da 32 bit (transfer conserva i bit)
    a0 = transfer(e0, a0)
    a1 = transfer(e1, a1)

    ! Interpreta ciascuna parola come unsigned 32: se negativa aggiungi 2^32
    lo0 = int(a0(1), c_long_long); if (lo0 < 0_c_long_long) lo0 = lo0 + 4294967296_c_long_long
    hi0 = int(a0(2), c_long_long); if (hi0 < 0_c_long_long) hi0 = hi0 + 4294967296_c_long_long
    lo1 = int(a1(1), c_long_long); if (lo1 < 0_c_long_long) lo1 = lo1 + 4294967296_c_long_long
    hi1 = int(a1(2), c_long_long); if (hi1 < 0_c_long_long) hi1 = hi1 + 4294967296_c_long_long

    ! Sottrazione unsigned 64: (hi1:lo1) - (hi0:lo0) modulo 2^64
    if (lo1 >= lo0) then
      dlo = lo1 - lo0
      borrow = 0_c_long_long
    else
      dlo = (lo1 + 4294967296_c_long_long) - lo0
      borrow = 1_c_long_long
    end if

    dhi = hi1 - hi0 - borrow
    if (dhi < 0_c_long_long) dhi = dhi + 4294967296_c_long_long

    ! de_mJ = dhi*2^32 + dlo, poi mJ->J
    dJ = real(1.0d-3 * ( real(dhi, c_double) * 4294967296.0d0 + real(dlo, c_double) ),kind=PRC)
  end function u64_delta_J

  
#endif
end module nvml_interface

