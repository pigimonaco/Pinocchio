#pragma once

#if defined(_ENERGY_RAPL_) || defined(_ENERGY_NVIDIA_) || defined(_ENERGY_AMD_)
   void Create_PMT(const int * const devID, const int numGPUs, const int Task, const int NTasks);
#endif // defined(_ENERGY_RAPL_) || defined(_ENERGY_NVIDIA_) || defined(_ENERGY_AMD_)

#if defined(_ENERGY_RAPL_)
   void Start_PMT_CPU(const char *string, const int Task);
   void Stop_PMT_CPU(const char *string, const int Task);
   void Show_PMT_CPU(const char *string, const int Task);
#endif // _ENERGY_RAPL_

#if defined(_ENERGY_NVIDIA_) || defined(_ENERGY_AMD_)
   void Start_PMT_GPU(const char *string, const int devID, const int Task);
   void Stop_PMT_GPU(const char *string, const int devID, const int Task);
   void Show_PMT_GPU(const char *string, const int task);
#endif // defined(_ENERGY_NVIDIA_) || defined(_ENERGY_AMD_)


