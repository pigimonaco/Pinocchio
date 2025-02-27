#pragma once

#include "energy_pmt_methods.h"

#if defined(_ENERGY_RAPL_) || defined(_ENERGY_NVIDIA_) || defined(_ENERGY_AMD_)
   #define PMT_CREATE(devID, numGPUs, Task, NTasks) Create_PMT((devID), (numGPUs), (Task), (NTasks))
#else
   #define PMT_CREATE(devID, numGPUs, Task, NTasks)
#endif // defined(_ENERGY_RAPL_) || defined(_ENERGY_NVIDIA_) || defined(_ENERGY_AMD_)

#if defined(_ENERGY_RAPL_)
   #define PMT_CPU_START(string, task) Start_PMT_CPU((string), (task))
   #define PMT_CPU_STOP(string, task)  Stop_PMT_CPU((string), (task))
   #define PMT_CPU_SHOW(string, task)  Show_PMT_CPU((string), (task))
#else
   #define PMT_CPU_START(string, task)
   #define PMT_CPU_STOP(string, task)
   #define PMT_CPU_SHOW(string, task)
#endif // _ENERGY_RAPL_

#if defined(_ENERGY_NVIDIA_) || defined(_ENERGY_AMD_)
   #define PMT_GPU_START(string, devID, task) Start_PMT_GPU((string), (devID), (task))
   #define PMT_GPU_STOP(string, devID, task)  Stop_PMT_GPU((string), (devID), (task))
   #define PMT_GPU_SHOW(string, task)         Show_PMT_GPU((string), (task))
#else
   #define PMT_GPU_START(string, devID, task)
   #define PMT_GPU_STOP(string, devID, task)
   #define PMT_GPU_SHOW(string, task)
#endif // defined(_ENERGY_NVIDIA_) || defined(_ENERGY_AMD_)
