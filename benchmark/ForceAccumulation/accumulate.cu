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
#include "../../src/UnitTesting/unit_test.h"

using omni::constants::warp_size_int;
using omni::constants::ExceptionResponse;
using omni::data_types::llint;
using omni::data_types::ullint;
using omni::numerics::max_int_accumulation_f;
using omni::numerics::max_int_accumulation_ll;
using namespace omni::card;
using namespace omni::testing;

// Copy the inline __device__ functions
#include "../../src/Potential/accumulation.cui"

//-------------------------------------------------------------------------------------------------
// This kernel will accumulate values using split int32 accumulators.
//-------------------------------------------------------------------------------------------------
__global__ void __launch_bounds__(512, 2) kAddSplit(int* overflow, llint* result) {
  __shared__ int primary[512], overflow_active[512 / warp_size_int];
  float contrib = (float)(threadIdx.x);
  float incr = 1.3f;
  primary[threadIdx.x] = 0;
  const size_t pos = (blockIdx.x * blockDim.x) + threadIdx.x;
  overflow[pos] = 0;
  int* ovrf_ptr = &overflow[blockIdx.x * blockDim.x];
  for (int i = 0; i < 1000000; i++) {
    if (contrib > 57.5f) {
      incr = -1.2f;
    }
    else if (contrib < 31.5f) {
      incr = 1.3f;
    }
    contrib += incr;
    splitForceContribution(contrib * 64.0f, threadIdx.x, primary, overflow_active,
                           ovrf_ptr);
  }
  if (overflow_active[threadIdx.x / warp_size_int]) {
    const llint ovrf_val = overflow[pos];
    result[pos] = (ovrf_val * max_int_accumulation_ll) + (llint)(primary[threadIdx.x]);
  }
  else {
    result[pos] = (llint)(primary[threadIdx.x]);
  }
}

//-------------------------------------------------------------------------------------------------
// This kernel will accumulate values using unified (standard) int64 accumulators.
//-------------------------------------------------------------------------------------------------
__global__ void __launch_bounds__(512, 2) kAddUnified(llint* result) {
  __shared__ llint primary[512];
  float contrib = (float)(threadIdx.x);
  float incr = 1.3f;
  primary[threadIdx.x] = 0LL;
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
  result[(blockIdx.x * blockDim.x) + threadIdx.x] = primary[threadIdx.x];
}

//-------------------------------------------------------------------------------------------------
// main
//-------------------------------------------------------------------------------------------------
int main(const int argc, const char* argv[]) {

  // Prime the timings device
  StopWatch timer;
  timer.addCategory("Split accumulation");
  timer.addCategory("Unified accumulation");

  // Some baseline initialization
  TestEnvironment oe(argc, argv);

  // Select a GPU
  HpcConfig gpu_config(ExceptionResponse::WARN);
  std::vector<int> my_gpus = gpu_config.getGpuDevice(1);
  GpuDetails gpu = gpu_config.getGpuInfo(my_gpus[0]);
  const int nblocks  = gpu.getSMPCount() * 2;
  const int nthreads = 512;
  const size_t buffer_size = nblocks * nthreads;
  
  // Prepare a buffer to hold overflow data
  Hybrid<int> overflow(buffer_size, "overflow");
  Hybrid<llint> result_split(HybridKind::ARRAY, "res_split", HybridFormat::HOST_ONLY, buffer_size);
  Hybrid<llint> result_unified(HybridKind::ARRAY, "res_unified", HybridFormat::HOST_ONLY,
                               buffer_size);
  int* devc_overflow = overflow.data(HybridTargetLevel::DEVICE);
  llint* split_ptr = result_split.data();
  llint* unified_ptr = result_unified.data();
  cudaFuncSetSharedMemConfig(kAddUnified, cudaSharedMemBankSizeEightByte);
  timer.assignTime(0);
  
  // Launch each kernel 100x
  for (int i = 0; i < 100; i++) {
    kAddSplit<<<nblocks, nthreads>>>(devc_overflow, split_ptr);
  }
  cudaDeviceSynchronize();
  timer.assignTime(1);
  for (int i = 0; i < 100; i++) {
    kAddUnified<<<nblocks, nthreads>>>(unified_ptr);
  }
  cudaDeviceSynchronize();
  timer.assignTime(2);
  check(result_split.readHost(), RelationalOperator::EQUAL, result_unified.readHost(), "Split "
        "accumulation does not produce the same outcomes as unified int64 accumulation.");
  
  // Report the timings
  if (oe.getDisplayTimingsOrder()) {
    timer.assignTime(0);
    timer.printResults();
  }

  // Print results
  printTestSummary(oe.getVerbosity());
}
