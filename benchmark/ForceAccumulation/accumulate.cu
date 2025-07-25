// -*-c++-*-
#include <cuda_runtime.h>
#include <cusolverDn.h>
#include <nvml.h>
#include "../../src/Accelerator/hpc_config.h"
#include "../../src/Accelerator/hybrid.h"
#include "../../src/Accelerator/ptx_macros.h"
#include "../../src/Constants/behavior.h"
#include "../../src/Constants/hpc_bounds.h"
#include "../../src/Constants/scaling.h"
#include "../../src/DataTypes/common_types.h"
#include "../../src/DataTypes/stormm_vector_types.h"
#include "../../src/Namelists/command_line_parser.h"
#include "../../src/Namelists/namelist_emulator.h"
#include "../../src/Namelists/namelist_enumerators.h"
#include "../../src/Numerics/split_fixed_precision.h"
#include "../../src/Parsing/parse.h"
#include "../../src/Parsing/parsing_enumerators.h"
#include "../../src/UnitTesting/stopwatch.h"
#include "../../src/UnitTesting/unit_test.h"

using stormm::constants::CaseSensitivity;
using stormm::constants::ExceptionResponse;
using stormm::constants::warp_size_int;
using stormm::constants::warp_size_zu;
using stormm::constants::warp_bits;
using stormm::data_types::int95_t;
using stormm::data_types::llint;
using stormm::data_types::ullint;
using stormm::errors::rtErr;
using stormm::namelist::NamelistEmulator;
using stormm::namelist::NamelistType;
using stormm::numerics::max_int_accumulation;
using stormm::numerics::max_int_accumulation_f;
using stormm::numerics::max_int_accumulation_ll;
using stormm::numerics::max_llint_accumulation;
using stormm::numerics::max_llint_accumulation_f;
using stormm::numerics::max_short_accumulation;
using stormm::parse::NumberFormat;
using stormm::parse::strcmpCased;
using stormm::parse::verifyNumberFormat;
using namespace stormm::card;
using namespace stormm::testing;

// Copy the inline __device__ functions
#include "../../src/Numerics/accumulation.cui"

// Define the cycle counts for various kernels
#define SMALL_CYCLE_COUNT  250000
#define LARGE_CYCLE_COUNT  1000000

//-------------------------------------------------------------------------------------------------
// Perform arithmetic with +, -, and * using int32 numbers.
//-------------------------------------------------------------------------------------------------
__global__ void __launch_bounds__(512, 2) kWorkInt32(llint* result) {
  int quant_i  = threadIdx.x;
  int quant_ii = blockIdx.x;
  int incr_i  = 1;
  int incr_ii = 3;
  int mult = -1;
  for (int i = 0; i < LARGE_CYCLE_COUNT; i++) {
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
  for (int i = 0; i < LARGE_CYCLE_COUNT; i++) {
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
__global__ void __launch_bounds__(512, 2) kAddSplit(llint* result, const float scale_factor) {
  __shared__ int primary[512], overflow[512];
  float contrib = (float)(threadIdx.x);
  float incr = 1.3f;
  primary[threadIdx.x] = 0;
  const size_t pos = (blockIdx.x * blockDim.x) + threadIdx.x;
  overflow[threadIdx.x] = 0;
  for (int i = 0; i < SMALL_CYCLE_COUNT; i++) {
    if (contrib > 57.5f) {
      incr = -1.2f;
    }
    else if (contrib < 31.5f) {
      incr = 1.3f;
    }
    contrib += incr;
    atomicSplit(contrib * scale_factor, threadIdx.x, primary, overflow);
    contrib += incr;
    atomicSplit(contrib * scale_factor, threadIdx.x, primary, overflow);
    contrib += incr;
    atomicSplit(contrib * scale_factor, threadIdx.x, primary, overflow);
    contrib += incr;
    atomicSplit(contrib * scale_factor, threadIdx.x, primary, overflow);
  }
  const llint ovrf_val = overflow[threadIdx.x];
  result[pos] = (ovrf_val * max_int_accumulation_ll) + (llint)(primary[threadIdx.x]);
}

//-------------------------------------------------------------------------------------------------
// This kernel will accumulate values using unified (standard) int64 accumulators.
//-------------------------------------------------------------------------------------------------
__global__ void __launch_bounds__(512, 2) kAddUnified(llint* result, const float scale_factor) {
  __shared__ llint primary[512];
  float contrib = (float)(threadIdx.x);
  float incr = 1.3f;
  primary[threadIdx.x] = 0LL;
  for (int i = 0; i < SMALL_CYCLE_COUNT; i++) {
    if (contrib > 57.5f) {
      incr = -1.2f;
    }
    else if (contrib < 31.5f) {
      incr = 1.3f;
    }
    contrib += incr;
    atomicAdd((ullint*)&primary[threadIdx.x], (ullint)(__float2ll_rn(contrib * scale_factor)));
    contrib += incr;
    atomicAdd((ullint*)&primary[threadIdx.x], (ullint)(__float2ll_rn(contrib * scale_factor)));
    contrib += incr;
    atomicAdd((ullint*)&primary[threadIdx.x], (ullint)(__float2ll_rn(contrib * scale_factor)));
    contrib += incr;
    atomicAdd((ullint*)&primary[threadIdx.x], (ullint)(__float2ll_rn(contrib * scale_factor)));
  }
  result[(blockIdx.x * blockDim.x) + threadIdx.x] = primary[threadIdx.x];
}

//-------------------------------------------------------------------------------------------------
// This kernel will accumulate values using split int32 accumulators.
//-------------------------------------------------------------------------------------------------
__global__ void __launch_bounds__(512, 2) kAddSplitGbl(llint* result, int* primary,
                                                       int* overflow, const float scale_factor) {
  float contrib = (float)(threadIdx.x);
  float incr = 1.3f;
  const size_t pos = (blockIdx.x * blockDim.x) + threadIdx.x;
  primary[pos] = 0;
  overflow[pos] = 0;
  for (int i = 0; i < SMALL_CYCLE_COUNT; i++) {
    if (contrib > 57.5f) {
      incr = -1.2f;
    }
    else if (contrib < 31.5f) {
      incr = 1.3f;
    }
    contrib += incr;
    atomicSplit(contrib * scale_factor, pos, primary, overflow);
    contrib += incr;
    atomicSplit(contrib * scale_factor, pos, primary, overflow);
    contrib += incr;
    atomicSplit(contrib * scale_factor, pos, primary, overflow);
    contrib += incr;
    atomicSplit(contrib * scale_factor, pos, primary, overflow);
  }
  const llint ovrf_val = overflow[pos];
  result[pos] = (ovrf_val * max_int_accumulation_ll) + (llint)(primary[pos]);
}

//-------------------------------------------------------------------------------------------------
// This kernel will accumulate values using split int32 accumulators, with the overflow stored in
// __global__ memory and the primary accumulator in __shared__.
//-------------------------------------------------------------------------------------------------
__global__ void __launch_bounds__(512, 2) kAddSplitMix(llint* result, int* overflow,
                                                       const float scale_factor) {
  __shared__ int primary[512];
  float contrib = (float)(threadIdx.x);
  float incr = 1.3f;
  const size_t pos = (blockIdx.x * blockDim.x) + threadIdx.x;
  const size_t base = blockIdx.x * blockDim.x;
  primary[threadIdx.x] = 0;
  overflow[pos] = 0;
  for (int i = 0; i < SMALL_CYCLE_COUNT; i++) {
    if (contrib > 57.5f) {
      incr = -1.2f;
    }
    else if (contrib < 31.5f) {
      incr = 1.3f;
    }
    contrib += incr;
    atomicSplit(contrib * scale_factor, threadIdx.x, primary, &overflow[base]);
    contrib += incr;
    atomicSplit(contrib * scale_factor, threadIdx.x, primary, &overflow[base]);
    contrib += incr;
    atomicSplit(contrib * scale_factor, threadIdx.x, primary, &overflow[base]);
    contrib += incr;
    atomicSplit(contrib * scale_factor, threadIdx.x, primary, &overflow[base]);
  }
  const llint ovrf_val = overflow[pos];
  result[pos] = (ovrf_val * max_int_accumulation_ll) + (llint)(primary[threadIdx.x]);
}

//-------------------------------------------------------------------------------------------------
// This kernel will accumulate values using unified (standard) int64 accumulators.
//-------------------------------------------------------------------------------------------------
__global__ void __launch_bounds__(512, 2) kAddUnifiedGbl(llint* result, const float scale_factor) {
  float contrib = (float)(threadIdx.x);
  float incr = 1.3f;
  const size_t pos = (blockIdx.x * blockDim.x) + threadIdx.x;
  result[pos] = 0LL;
  for (int i = 0; i < SMALL_CYCLE_COUNT; i++) {
    if (contrib > 57.5f) {
      incr = -1.2f;
    }
    else if (contrib < 31.5f) {
      incr = 1.3f;
    }
    contrib += incr;
    atomicAdd((ullint*)&result[pos], (ullint)(__float2ll_rn(contrib * scale_factor)));
    contrib += incr;
    atomicAdd((ullint*)&result[pos], (ullint)(__float2ll_rn(contrib * scale_factor)));
    contrib += incr;
    atomicAdd((ullint*)&result[pos], (ullint)(__float2ll_rn(contrib * scale_factor)));
    contrib += incr;
    atomicAdd((ullint*)&result[pos], (ullint)(__float2ll_rn(contrib * scale_factor)));
  }
}

//-------------------------------------------------------------------------------------------------
// This kernel will accumulate values using unified (standard) int632 accumulators.
//-------------------------------------------------------------------------------------------------
__global__ void __launch_bounds__(512, 2) kAddUnifiedGbl(int* result, const float scale_factor) {
  float contrib = (float)(threadIdx.x);
  float incr = 1.3f;
  const size_t pos = (blockIdx.x * blockDim.x) + threadIdx.x;
  result[pos] = 0;
  for (int i = 0; i < SMALL_CYCLE_COUNT; i++) {
    if (contrib > 57.5f) {
      incr = -1.2f;
    }
    else if (contrib < 31.5f) {
      incr = 1.3f;
    }
    contrib += incr;
    atomicAdd(&result[pos], __float2int_rn(contrib * scale_factor));
    contrib += incr;
    atomicAdd(&result[pos], __float2int_rn(contrib * scale_factor));
    contrib += incr;
    atomicAdd(&result[pos], __float2int_rn(contrib * scale_factor));
    contrib += incr;
    atomicAdd(&result[pos], __float2int_rn(contrib * scale_factor));
  }
}

//-------------------------------------------------------------------------------------------------
// main
//-------------------------------------------------------------------------------------------------
int main(const int argc, const char* argv[]) {

  // Prime the timings device
  StopWatch timer;
  timer.addCategory("Split 2x int32_t accumulation");
  timer.addCategory("Unified int64_taccumulation");

  // Parse command line instructions
  CommandLineParser clip("accumulate", "A benchmarking program to test the performance of split "
                         "fixed-precision accumulation.", { "-timings" });
  clip.addStandardBenchmarkingInputs("-trials");
  NamelistEmulator *t_nml = clip.getNamelistPointer();
  t_nml->addKeyword("-min_bits", NamelistType::INTEGER, std::to_string(12));
  t_nml->addHelp("-min_bits", "The minimum number of bits after the point (in the fraction part "
                 "of each number representation) to test");
  t_nml->addKeyword("-max_bits", NamelistType::INTEGER, std::to_string(36));
  t_nml->addHelp("-max_bits", "The maximum number of bits after the point (in the fraction part "
                 "of each number representation) to test");
  TestEnvironment oe(argc, argv, &clip, TmpdirStatus::NOT_REQUIRED, ExceptionResponse::SILENT);
  clip.parseUserInput(argc, argv);
  const int min_bits = t_nml->getIntValue("-min_bits");
  const int max_bits = t_nml->getIntValue("-max_bits");
  const int trials = t_nml->getIntValue("-trials");

  // Select a GPU
  HpcConfig gpu_config(ExceptionResponse::WARN);
  std::vector<int> my_gpus = gpu_config.getGpuDevice(1);
  GpuDetails gpu = gpu_config.getGpuInfo(my_gpus[0]);
  const int nblocks  = gpu.getSMPCount() * 2;
  const int nthreads = 512;
  const size_t buffer_size = nblocks * nthreads;
  
  // Prepare a buffer to hold overflow data
  Hybrid<int> primary(buffer_size, "primary");
  Hybrid<int> overflow(buffer_size, "overflow");
  Hybrid<llint> result_split(HybridKind::ARRAY, "res_split", HybridFormat::EXPEDITED, buffer_size);
  Hybrid<llint> result_unified(HybridKind::ARRAY, "res_unified", HybridFormat::EXPEDITED,
                               buffer_size);
  int* primary_ptr = primary.data(HybridTargetLevel::DEVICE);
  int* overflow_ptr = overflow.data(HybridTargetLevel::DEVICE);
  llint* split_ptr = result_split.data(HybridTargetLevel::DEVICE);
  llint* unified_ptr = result_unified.data(HybridTargetLevel::DEVICE);
  timer.assignTime(0);
  
  // Launch each kernel 10x and perform multiple trials.  Benchmark __shared__ memory accumulation
  // first.
  for (int bits = min_bits; bits <= max_bits; bits++) {
    const int ti_split = timer.addCategory("Split accumulation, " + std::to_string(bits));
    const float scale = pow(2.0, bits);
    timer.assignTime(0);
    for (int n = 0; n < trials; n++) {
      for (int i = 0; i < 10; i++) {
        kAddSplit<<<nblocks, nthreads>>>(split_ptr, scale);
      }
      cudaDeviceSynchronize();
      timer.assignTime(ti_split);
    }
  }
  for (int bits = min_bits; bits <= max_bits; bits++) {
    const int ti_unite = timer.addCategory("Unified accumulation, " + std::to_string(bits));
    const float scale = pow(2.0, bits);
    timer.assignTime(0);
    for (int n = 0; n < trials; n++) {
      for (int i = 0; i < 10; i++) {
        kAddUnified<<<nblocks, nthreads>>>(unified_ptr, scale);
      }
      cudaDeviceSynchronize();
      timer.assignTime(ti_unite);
    }
  }
  result_split.download();
  result_unified.download();
  check(result_split.readHost(), RelationalOperator::EQUAL, result_unified.readHost(), "Split "
        "accumulation does not produce the same outcomes as unified int64 accumulation when "
        "performing operations in __shared__ memory.");

  // Benchmark mixed __shared__ and __global__ memory accumulation.
  for (int bits = min_bits; bits <= max_bits; bits++) {
    const int ti_split = timer.addCategory("Split accumulation, Mix " + std::to_string(bits));
    const float scale = pow(2.0, bits);
    for (int n = 0; n < trials; n++) {
      timer.assignTime(0);
      for (int i = 0; i < 10; i++) {
        kAddSplitMix<<<nblocks, nthreads>>>(split_ptr, overflow_ptr, scale);
      }
      cudaDeviceSynchronize();
      timer.assignTime(ti_split);
    }
  }
  result_split.download();
  check(result_split.readHost(), RelationalOperator::EQUAL, result_unified.readHost(), "Split "
        "accumulation does not produce the same outcomes as unified int64 accumulation when "
        "performing operations with primary accumulators in __shared__ memory and overflow in "
        "__global__.");

  // Benchmark __global__ memory accumulation
  for (int bits = min_bits; bits <= max_bits; bits++) {
    const int ti_split = timer.addCategory("Split accumulation, L2 " + std::to_string(bits));
    const float scale = pow(2.0, bits);
    for (int n = 0; n < trials; n++) {
      timer.assignTime(0);
      for (int i = 0; i < 10; i++) {
        kAddSplitGbl<<<nblocks, nthreads>>>(split_ptr, primary_ptr, overflow_ptr, scale);
      }
      cudaDeviceSynchronize();
      timer.assignTime(ti_split);
    }
  }
  for (int bits = min_bits; bits <= max_bits; bits++) {
    const int ti_unite = timer.addCategory("Unified accumulation, L2 " + std::to_string(bits));
    const float scale = pow(2.0, bits);
    for (int n = 0; n < trials; n++) {
      timer.assignTime(0);
      for (int i = 0; i < 10; i++) {
        kAddUnifiedGbl<<<nblocks, nthreads>>>(unified_ptr, scale);
      }
      cudaDeviceSynchronize();
      timer.assignTime(ti_unite);
    }
    const int ti_short = timer.addCategory("Simple int32_t accumulation, L2 " +
                                           std::to_string(bits));
    for (int n = 0; n < trials; n++) {
      timer.assignTime(0);
      for (int i = 0; i < 10; i++) {
        kAddUnifiedGbl<<<nblocks, nthreads>>>(primary_ptr, scale);
      }
      cudaDeviceSynchronize();
      timer.assignTime(ti_short);
    }
  }
  result_split.download();
  result_unified.download();
  check(result_split.readHost(), RelationalOperator::EQUAL, result_unified.readHost(), "Split "
        "accumulation does not produce the same outcomes as unified int64 accumulation when "
        "performing operations in global memory.");
  
  // Make the GPU do int32 and int64-based calculations with +, -, and *
  const int ti_short = timer.addCategory("int32 Work");
  const int ti_long  = timer.addCategory("int64 Work");
  Hybrid<llint> testll(buffer_size);
  llint* testll_data = testll.data(HybridTargetLevel::DEVICE);
  timer.assignTime(0);
  for (int i = 0; i < 10; i++) {
    kWorkInt32<<<nblocks, nthreads>>>(testll_data);
  }
  cudaDeviceSynchronize();
  timer.assignTime(ti_short);
  for (int i = 0; i < 10; i++) {
    kWorkInt64<<<nblocks, nthreads>>>(testll_data);
  }
  cudaDeviceSynchronize();
  timer.assignTime(ti_long);  
    
  // Report the timings
  if (oe.getDisplayTimingsOrder()) {
    timer.assignTime(0);
    timer.printResults();
  }

  // Print results
  printTestSummary(oe.getVerbosity());
}
