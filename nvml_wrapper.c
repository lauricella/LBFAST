#include "defines.h"

/*
 * NVML helpers, compiled only with _NVML (Makefile targets nvfortran-nvml and
 * nvfortran-nvml-mpi, which link -lnvidia-ml) and called from Fortran through the bind(C)
 * interfaces of nvml_interface_module.f90. Each function initializes NVML at its first call
 * and never shuts it down. idx is an NVML device index; the Fortran callers pass mydev_c, the
 * OpenACC device number of the rank. NVML numbers the GPUs independently of
 * CUDA_VISIBLE_DEVICES, so the two numbers can refer to different GPUs when that variable
 * restricts or reorders the devices.
 * Return codes: -1 NVML initialization failed, -2 device lookup failed, -3 query failed.
 */
#ifdef _NVML
#include <nvml.h>
  
/* Current power draw of device idx in milliwatts (nvmlDeviceGetPowerUsage). */
int get_gpu_power_index(int idx) {
    static int inited = 0;
    if(!inited){
        if (nvmlInit_v2() != NVML_SUCCESS) return -1;
        inited = 1;
    }

    nvmlDevice_t dev;
    if (nvmlDeviceGetHandleByIndex_v2(idx, &dev) != NVML_SUCCESS) return -2;

    unsigned int p_mw = 0;
    if (nvmlDeviceGetPowerUsage(dev, &p_mw) != NVML_SUCCESS) return -3;

    return (int)p_mw;
}

/*
 * Total energy consumed by device idx since the driver was loaded, in millijoules
 * (nvmlDeviceGetTotalEnergyConsumption), stored in *e_mJ; returns 0 on success.
 */
int get_gpu_energy_mJ_u64(int idx, unsigned long long *e_mJ){
  static int inited = 0;
  if(!inited){
    if(nvmlInit_v2() != NVML_SUCCESS) return -1;
    inited = 1;
  }
  nvmlDevice_t dev;
  if(nvmlDeviceGetHandleByIndex_v2(idx, &dev) != NVML_SUCCESS) return -2;

  nvmlReturn_t r = nvmlDeviceGetTotalEnergyConsumption(dev, e_mJ);
  if(r != NVML_SUCCESS) return -3;
  return 0;
}

#endif
