// -*-c++-*-
#include <cuda_runtime.h>
#include <cusolverDn.h>
#include <nvml.h>
#include "../../src/Accelerator/hpc_config.h"
#include "../../src/Accelerator/hybrid.h"
#include "../../src/Accelerator/ptx_macros.h"
#include "../../src/Constants/behavior.h"
#include "../../src/Constants/scaling.h"
#include "../../src/Constants/fixed_precision.h"
#include "../../src/DataTypes/common_types.h"
#include "../../src/DataTypes/stormm_vector_types.h"
#include "../../src/UnitTesting/stopwatch.h"
#include "../../src/UnitTesting/unit_test.h"

using stormm::constants::warp_size_int;
using stormm::constants::ExceptionResponse;
using stormm::data_types::int95_t;
using stormm::data_types::llint;
using stormm::data_types::ullint;
using stormm::numerics::max_int_accumulation;
using stormm::numerics::max_int_accumulation_f;
using stormm::numerics::max_int_accumulation_ll;
using stormm::numerics::max_llint_accumulation;
using stormm::numerics::max_llint_accumulation_f;
using namespace stormm::card;
using namespace stormm::testing;

// Copy the inline __device__ functions
#include "../../src/Numerics/accumulation.cui"

//-------------------------------------------------------------------------------------------------
// Perform arithmetic with +, -, and * using int32 numbers.
//-------------------------------------------------------------------------------------------------
__global__ void __launch_bounds__(512, 2) kWorkInt32(llint* result) {
  int quant_i  = threadIdx.x;
  int quant_ii = blockIdx.x;
  int incr_i  = 1;
  int incr_ii = 3;
  int mult = -1;
  for (int i = 0; i < 1000000; i++) {
    quant_i += quant_ii + incr_i;
    incr_i += 1 - 2 * (quant_i > 10000000);
    incr_ii += 1 + (quant_i > 10000000) + (quant_i < -10000000);
    quant_i -= incr_ii;
    quant_i -= 100000 * (quant_i > 10000000);
    quant_ii += i * mult;
    mult *= -1;
  }
  result[blockIdx.x * blockDim.x + threadIdx.x] = quant_i;
}

//-------------------------------------------------------------------------------------------------
// Perform arithmetic with +, -, and * using int64 numbers.
//-------------------------------------------------------------------------------------------------
__global__ void __launch_bounds__(512, 2) kWorkInt64(llint* result) {
  llint quant_i  = threadIdx.x;
  llint quant_ii = blockIdx.x;
  llint incr_i  = 1;
  llint incr_ii = 3;
  llint mult = -1;
  for (int i = 0; i < 1000000; i++) {
    quant_i += quant_ii + incr_i;
    incr_i += 1LL - 2LL * (quant_i > 10000000LL);
    incr_ii += 1LL + (quant_i > 10000000LL) + (quant_i < -10000000LL);
    quant_i -= incr_ii;
    quant_i -= 100000LL * (quant_i > 10000000LL);
    quant_ii += i * mult;
    mult *= -1LL;
  }
  result[blockIdx.x * blockDim.x + threadIdx.x] = quant_i;
}

//-------------------------------------------------------------------------------------------------
// This kernel will accumulate values using split int32 accumulators.
//-------------------------------------------------------------------------------------------------
__global__ void __launch_bounds__(512, 2) kAddSplit(llint* result) {
  __shared__ int primary[512], overflow[512];
  float contrib = (float)(threadIdx.x);
  float incr = 1.3f;
  primary[threadIdx.x] = 0;
  const size_t pos = (blockIdx.x * blockDim.x) + threadIdx.x;
  overflow[threadIdx.x] = 0;
  for (int i = 0; i < 250000; i++) {
    if (contrib > 57.5f) {
      incr = -1.2f;
    }
    else if (contrib < 31.5f) {
      incr = 1.3f;
    }
    contrib += incr;
    atomicSplit(contrib * 256.0f, threadIdx.x, primary, overflow);
    contrib += incr;
    atomicSplit(contrib * 256.0f, threadIdx.x, primary, overflow);
    contrib += incr;
    atomicSplit(contrib * 256.0f, threadIdx.x, primary, overflow);
    contrib += incr;
    atomicSplit(contrib * 256.0f, threadIdx.x, primary, overflow);
  }
  const llint ovrf_val = overflow[threadIdx.x];
  result[pos] = (ovrf_val * max_int_accumulation_ll) + (llint)(primary[threadIdx.x]);
}

//-------------------------------------------------------------------------------------------------
// This kernel will accumulate values using unified (standard) int64 accumulators.
//-------------------------------------------------------------------------------------------------
__global__ void __launch_bounds__(512, 2) kAddUnified(llint* result) {
  __shared__ llint primary[512];
  float contrib = (float)(threadIdx.x);
  float incr = 1.3f;
  primary[threadIdx.x] = 0LL;
  for (int i = 0; i < 250000; i++) {
    if (contrib > 57.5f) {
      incr = -1.2f;
    }
    else if (contrib < 31.5f) {
      incr = 1.3f;
    }
    contrib += incr;
    atomicAdd((ullint*)&primary[threadIdx.x], lliToUlli(__float2ll_rn(contrib * 256.0f)));
    contrib += incr;
    atomicAdd((ullint*)&primary[threadIdx.x], lliToUlli(__float2ll_rn(contrib * 256.0f)));
    contrib += incr;
    atomicAdd((ullint*)&primary[threadIdx.x], lliToUlli(__float2ll_rn(contrib * 256.0f)));
    contrib += incr;
    atomicAdd((ullint*)&primary[threadIdx.x], lliToUlli(__float2ll_rn(contrib * 256.0f)));
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
    kAddSplit<<<nblocks, nthreads>>>(split_ptr);
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

  // Make the GPU do int32 and int64-based calculations with +, -, and *
  timer.addCategory("int32 Work");
  timer.addCategory("int64 Work");
  Hybrid<llint> testll(buffer_size);
  llint* testll_data = testll.data(HybridTargetLevel::DEVICE);
  for (int i = 0; i < 100; i++) {
    kWorkInt32<<<nblocks, nthreads>>>(testll_data);
  }
  cudaDeviceSynchronize();
  timer.assignTime(3);
  for (int i = 0; i < 100; i++) {
    kWorkInt64<<<nblocks, nthreads>>>(testll_data);
  }
  cudaDeviceSynchronize();
  timer.assignTime(4);  
    
  // Report the timings
  if (oe.getDisplayTimingsOrder()) {
    timer.assignTime(0);
    timer.printResults();
  }

  // Print results
  printTestSummary(oe.getVerbosity());
}
