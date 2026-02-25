#include "defines.h"

#ifdef MONITORENERGY
#include <nvml.h>
  
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
