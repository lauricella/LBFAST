module nvml_interface
  use iso_c_binding
  implicit none

  interface
    function get_gpu_power_index(idx) bind(C, name="get_gpu_power_index")
      import :: c_int
      integer(c_int), value :: idx
      integer(c_int) :: get_gpu_power_index
    end function
  end interface

end module nvml_interface

