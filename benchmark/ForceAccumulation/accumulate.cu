// -*-c++-*-
#include <cuda_runtime.h>
#include <cusolverDn.h>
#include <nvml.h>
#include "../../src/Accelerator/hpc_config.cuh"
#include "../../src/Accelerator/hybrid.h"
#include "../../src/Constants/behavior.h"
#include "../../src/Constants/scaling.h"
#include "../../src/Constants/fixed_precision.h"
#include "../../src/DataTypes/common_types.h"
#include "../../src/UnitTesting/stopwatch.h"

using omni::card::GpuDetails;
using omni::card::HpcConfig;
using omni::card::Hybrid;
using omni::card::HybridTargetLevel;
using omni::constants::warp_size_int;
using omni::constants::ExceptionResponse;
using omni::data_types::llint;
using omni::data_types::ullint;
using omni::numerics::max_int_accumulation_f;
using omni::numerics::max_int_accumulation_ll;
using omni::testing::StopWatch;

// Copy the inline __device__ functions
#include "../../src/Potential/accumulation.cui"

//-------------------------------------------------------------------------------------------------
// This kernel will accumulate values using split int32 accumulators.
//-------------------------------------------------------------------------------------------------
__global__ void __launch_bounds__(512, 2) kAddSplit(int* overflow) {
  __shared__ int primary[512], overflow_active[16];
  float contrib = 0.0f;
  float incr = 1.3f;
  for (int i = 0; i < 1000000; i++) {
    if (contrib > 57.5f) {
      incr = -1.2f;
    }
    else if (contrib < 31.5f) {
      incr = 1.3f;
    }
    contrib += incr;
    splitForceContribution(contrib * 64.0f, threadIdx.x, primary, overflow_active, overflow);
  }
}

//-------------------------------------------------------------------------------------------------
// This kernel will accumulate values using unified (standard) int64 accumulators.
//-------------------------------------------------------------------------------------------------
__global__ void __launch_bounds__(512, 2) kAddUnified() {
  __shared__ llint primary[512];
  float contrib = 0.0f;
  float incr = 1.3f;
  for (int i = 0; i < 1000000; i++) {
    if (contrib > 57.5f) {
      incr = -1.2f;
    }
    else if (contrib < 31.5f) {
      incr = 1.3f;
    }
    contrib += incr;
    atomicAdd((ullint*)&primary[threadIdx.x], llitoulli(__float2ll_rn(contrib * 64.0f)));
  }
}

//-------------------------------------------------------------------------------------------------
// main
//-------------------------------------------------------------------------------------------------
int main(const int argc, const char* argv[]) {

  // Prime the timings device
  StopWatch timer;
  timer.addCategory("Split accumulation");
  timer.addCategory("Unified accumulation");

  // Select a GPU
  HpcConfig gpu_config(ExceptionResponse::WARN);
  std::vector<int> my_gpus = gpu_config.getGpuDevice(1);
  GpuDetails gpu = gpu_config.getGpuInfo(my_gpus[0]);

  // Prepare a buffer to hold overflow data
  Hybrid<int> overflow(131072, "overflow");
  int* devc_overflow = overflow.data(HybridTargetLevel::DEVICE);
  cudaFuncSetSharedMemConfig(kAddUnified, cudaSharedMemBankSizeEightByte);
  timer.assignTime(0);

  // Launch each kernel 100x
  for (int i = 0; i < 100; i++) {
    kAddSplit<<<56, 512>>>(devc_overflow);
  }
  cudaDeviceSynchronize();
  timer.assignTime(1);
  for (int i = 0; i < 100; i++) {
    kAddUnified<<<56, 512>>>();
  }
  cudaDeviceSynchronize();
  timer.assignTime(2);

  // Report the results
  timer.printResults();
}
