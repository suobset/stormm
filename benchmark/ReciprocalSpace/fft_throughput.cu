// -*-c++-*-
#include <cuda.h>
#include <cuda_runtime.h>
#include <cufft.h>
#include <string>
#include <vector>
#include "copyright.h"
#include "../../src/Accelerator/gpu_details.h"
#include "../../src/Accelerator/hpc_config.h"
#include "../../src/Accelerator/hybrid.h"
#include "../../src/Constants/behavior.h"
#include "../../src/DataTypes/common_types.h"
#include "../../src/Random/random.h"
#include "../../src/Reporting/summary_file.h"
#include "../../src/UnitTesting/test_environment.h"
#include "../../src/UnitTesting/test_system_manager.h"
#include "../../src/UnitTesting/stopwatch.h"

using namespace stormm::card;
using namespace stormm::constants;
using namespace stormm::review;
using namespace stormm::testing;

//-------------------------------------------------------------------------------------------------
// main
//-------------------------------------------------------------------------------------------------
int main(const int argc, const char* argv[]) {

  // Baseline variables
  TestEnvironment oe(argc, argv, ExceptionResponse::SILENT);
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmSplash();
  }
  StopWatch timer;
#ifdef STORMM_USE_HPC
  const HpcConfig gpu_config(ExceptionResponse::WARN);
  const std::vector<int> my_gpus = gpu_config.getGpuDevice(1);
  const GpuDetails gpu = gpu_config.getGpuInfo(my_gpus[0]);
  Hybrid<int> force_gpu_to_engage(1);
#endif

  // Perform a simple FFT, 8 x 8 x 8, and check the result.
  Xoshiro256ppGenerator xrs;
  Hybrid<double> trial(512, "trial_fft");
  const std::vector<double> trial_fill = gaussianRand(&xrs, trial.size(), 1.0);
  trial.putHost(trial_fill);
  trial.upload();

  // Lay out cuFFT essentials
  cufftHandle trial_plan;
#if 0
  cufftPlan3d(&trial_plan, 8, 8, 8, CUFFT_R2C);
#endif
  
  // Summary evaluation
  if (oe.getDisplayTimingsOrder()) {
    timer.assignTime(0);
    timer.printResults();
  }
  return 0;
}
