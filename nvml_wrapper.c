#include <nvml.h>

unsigned int get_gpu_power_index(int idx) {
    static int inited = 0;
    if(!inited){ nvmlInit(); inited = 1; }

    nvmlDevice_t dev;
    if (nvmlDeviceGetHandleByIndex(idx, &dev) != NVML_SUCCESS) return 0;

    unsigned int p_mw = 0;
    if (nvmlDeviceGetPowerUsage(dev, &p_mw) != NVML_SUCCESS) return 0;

    return p_mw;
}
