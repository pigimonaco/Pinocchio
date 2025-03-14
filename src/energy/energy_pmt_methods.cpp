#if defined(_ENERGY_PMT_)

#include <pmt.h>
#include <assert.h>
#include <iostream>
#include <vector>
#include <string>
#include <map>

struct EnergyState
{
  pmt::State   start;
  pmt::State   stop;
  double       joules;
  double       watts;
  double       seconds;
  unsigned int count;
};

static std::vector<bool> PMT_ERROR;

void PMT_err(const int task)
{
  std::cout << "\n\t Task: " << task << ": PMT Error \n" << std::endl;
  
  return;
}

#if defined(_ENERGY_RAPL_)
   static std::vector<std::unique_ptr<pmt::PMT>> sensor_cpu;
   static std::vector<std::map<std::string, EnergyState>> state_cpu;
#endif // _ENERGY_RAPL_

#if defined(_ENERGY_NVIDIA_) || defined(_ENERGY_AMD_)
   static std::vector<std::map<int, std::unique_ptr<pmt::PMT>>> sensor_gpu;
   static std::vector<std::map<int, std::map<std::string, EnergyState>>> state_gpu;
#endif // defined(_ENERGY_NVIDIA_) || defined(_ENERGY_AMD_)

//#if !defined(__NVCC__) && !defined(__NVCOMPILER)
   extern "C"
      {
         #include "energy_pmt_methods.h"
      }
//#endif // !defined(__NVCC__) && !defined(__NVCOMPILER)

#if defined(_ENERGY_RAPL_) || defined(_ENERGY_NVIDIA_) || defined(_ENERGY_AMD_)
void Create_PMT(const int * const  devID,
		const int          numGPUs,
		const int          task,
		const int          ntasks)
{
  PMT_ERROR.resize(ntasks);
  PMT_ERROR.at(task) = ((ntasks > 0) ? false : true);

  if (PMT_ERROR.at(task))
    {
      PMT_err(task);
      return;
    }

#if defined(_ENERGY_RAPL_)

  sensor_cpu.resize(ntasks);
  
  sensor_cpu.at(task) = pmt::rapl::Rapl::Create();

  state_cpu.resize(ntasks);
  
#endif // _ENERGY_RAPL_

#if defined(_ENERGY_NVIDIA_) || defined(_ENERGY_AMD_)
  if ((numGPUs > 0) && (devID != nullptr))
    {
      sensor_gpu.resize(ntasks);
      state_gpu.resize(ntasks);
      
      for (int dev=0 ; dev<numGPUs ; dev++)
	{
#if defined(_ENERGY_NVIDIA_)
	  sensor_gpu.at(task).insert(std::pair<int, std::unique_ptr<pmt::PMT>>(devID[dev],
									       pmt::nvml::NVML::Create(devID[dev])));
#elif defined(_ENERGY_AMD_)
	  sensor_gpu.at(task).insert(std::pair<int, std::unique_ptr<pmt::PMT>>(devID[dev],
									       pmt::rocm::ROCM::Create(devID[dev])));
#endif
	} // numGPUs
    }
#endif // defined(_ENERGY_NVIDIA_) || defined(_ENERGY_AMD_)
  
  return;
}
#endif // defined(_ENERGY_RAPL_) || defined(_ENERGY_NVIDIA_) || defined(_ENERGY_AMD_)

#if defined(_ENERGY_RAPL_)
void Start_PMT_CPU(const char *label,
		   const int   task)
{
  if (PMT_ERROR.at(task))
    {
      PMT_err(task);
      return;
    }

  const std::string tag{std::string{label}};

  // check if the label already exists
  if (state_cpu.at(task).count(tag))
    {
      state_cpu.at(task)[tag].start = sensor_cpu.at(task)->Read();
    }
  else
    {
      // create new EnergyState
      const EnergyState newState{sensor_cpu.at(task)->Read(),
                                 static_cast<pmt::State>(0),
                                 static_cast<double>(0),
                                 static_cast<double>(0),
                                 static_cast<double>(0),
                                 static_cast<unsigned int>(0)};

      // insert the key and initialize the counters
      state_cpu.at(task).insert(std::pair<std::string, EnergyState>(tag, newState));
    }

  return;
}

void Stop_PMT_CPU(const char *label,
		  const int   task)
{
  if (PMT_ERROR.at(task))
    {
      PMT_err(task);
      return;
    }

  const std::string tag{std::string{label}};
  
  // check if the label already exists
  // if not error
  if (!state_cpu.at(task).count(tag))
    {
      PMT_ERROR.at(task) = true;
      PMT_err(task);
      return;
    }
  else
    {
      // get the energy state
      EnergyState &State = state_cpu.at(task)[tag];
      
      // read the counter
      State.stop = sensor_cpu.at(task)->Read();

      // update quantities
      State.seconds += sensor_cpu.at(task)->seconds(State.start, State.stop);
      
      State.joules  += sensor_cpu.at(task)->joules(State.start, State.stop);

      State.watts   += sensor_cpu.at(task)->watts(State.start, State.stop);

      State.count++;
    }
  
  return;
}

void Show_PMT_CPU(const char *label,
		  const int   task)
{
  const std::string tag{std::string{label}};
  
  if (PMT_ERROR.at(task))
    {
      PMT_err(task);
      return;
    }
  else if (state_cpu.at(task).count(tag))
    {
      std::cout << "\n\t Task: " << task << std::endl;
      std::cout << "\n\t\t CPU Kernel:" << tag << ":" << std::endl;
      std::cout << "\t\t\t" << state_cpu.at(task)[tag].seconds << " [S]" << std::endl;
      std::cout << "\t\t\t" << state_cpu.at(task)[tag].joules  << " [J]" << std::endl;
      std::cout << "\t\t\t" << state_cpu.at(task)[tag].watts / state_cpu.at(task)[tag].count  << " [W]" << "\n" << std::endl;
    }
  
  return;
}
#endif // _ENERGY_RAPL_

#if defined(_ENERGY_NVIDIA_) || defined(_ENERGY_AMD_)
void Start_PMT_GPU(const char *label,
		   const int   devID,
		   const int   task)
{
  if (PMT_ERROR.at(task) || !sensor_gpu.at(task).count(devID))
    {
      PMT_err(task);
      return;
    }

  const std::string tag{std::string{label}};
  
  // check if the devID already exists
  if (!state_gpu.at(task).count(devID))
    {
      // initialize the counters
      const EnergyState newState{sensor_gpu.at(task)[devID]->Read(),
                                 static_cast<pmt::State>(0),
                                 static_cast<double>(0),
                                 static_cast<double>(0),
                                 static_cast<double>(0),
                                 static_cast<unsigned int>(0)};

      // insert devID and the associated pair< label, EnergyState >
      std::map<std::string, EnergyState> tag_EnergyState;
      tag_EnergyState.insert(std::pair<std::string, EnergyState>(tag, newState));
      state_gpu.at(task).insert(std::pair<int, std::map<std::string, EnergyState>>(devID, tag_EnergyState));
    }
  else
    {
      // the tag exists
      if (state_gpu.at(task)[devID].count(tag))
	{
	  // read the sensor
	  state_gpu.at(task)[devID][tag].start = sensor_gpu.at(task)[devID]->Read();
	}
      else
	{
	  // insert the label and initialize the counters
	  const EnergyState newState{sensor_gpu.at(task)[devID]->Read(),
				     static_cast<pmt::State>(0),
				     static_cast<double>(0),
				     static_cast<double>(0),
				     static_cast<double>(0),
				     static_cast<unsigned int>(0)};

	  state_gpu.at(task)[devID].insert(std::pair<std::string, EnergyState>(tag, newState));
	}
    }

  return;
}

void Stop_PMT_GPU(const char *label,
		  const int   devID,
		  const int   task)
{
  // check if the devID already exists
  // if not error
  if (!state_gpu.at(task).count(devID) || PMT_ERROR.at(task) || !sensor_gpu.at(task).count(devID))
    {
      PMT_ERROR.at(task) = true;
      PMT_err(task);
      return;
    }
  else
    {
      const std::string tag{std::string{label}};
      
      // check if the label already exists
      // if not error
      if (!state_gpu.at(task)[devID].count(tag))
	{
	  PMT_ERROR.at(task) = true;
	  PMT_err(task);
	  return;
	}
      else
	{
	  EnergyState &State = state_gpu.at(task)[devID][tag];
	  
	  // read the counter
	  State.stop = sensor_gpu.at(task)[devID]->Read();

	  // update quantities
	  State.seconds +=
	    sensor_gpu.at(task)[devID]->seconds(State.start,
						State.stop);
      
	  State.joules +=
	    sensor_gpu.at(task)[devID]->joules(State.start,
					       State.stop);

	  State.watts +=
	    sensor_gpu.at(task)[devID]->watts(State.start,
					      State.stop);
      
	  State.count++;
	}
    }
  
  return;
}

void Show_PMT_GPU(const char *label,
		  const int   task)
{
  if (PMT_ERROR.at(task))
    {
      PMT_err(task);
      return;
    }
  else
    {
      const std::string tag{std::string{label}};

      std::cout << "\n\t Task: " << task << std::endl;
      
      for (const auto &[key, value]: state_gpu.at(task))
	{
	  if (value.count(tag))
	    {
	      std::cout << "\n\t GPU [" << key << "] kernel:" << tag << ":" << std::endl;
	      std::cout << "\t\t" << value.at(tag).seconds << " [s]" << std::endl;
	      std::cout << "\t\t" << value.at(tag).joules  << " [J]" << std::endl;
	      std::cout << "\t\t" << value.at(tag).watts / value.at(tag).count  << " [W]" << "\n" << std::endl;
	    }
	}
    }
  
  return;
}

#endif // defined(_ENERGY_NVIDIA_) || defined(_ENERGY_AMD_)

#endif // _ENERGY_PMT_
