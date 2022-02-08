// -*-c++-*-
#include "../../src/Constants/behavior.h"
#include "../../src/Constants/hpc_bounds.h"
#include "../../src/Accelerator/hybrid.h"
#include "../../src/Accelerator/hpc_config.cuh"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/Math/summation.h"
#include "../../src/Math/hpc_summation.cuh"
#include "../../src/Reporting/error_format.h"
#include "../../src/Random/random.h"
#include "../../src/Random/hpc_random.cuh"
#include "../../src/UnitTesting/unit_test.h"

// CHECK
#include "../../src/UnitTesting/benchmark.h"
// END CHECK

using omni::constants::ExceptionResponse;
using omni::constants::tiny;
using omni::constants::large_block_size;
using omni::constants::warp_bits_mask_int;
using omni::card::HpcConfig;
using omni::card::Hybrid;
using omni::card::HybridFormat;
using omni::card::HybridTargetLevel;
using omni::card::GpuDetails;
using omni::data_types::llint;
using omni::data_types::ullint;
using omni::data_types::ullint2;
using omni::data_types::ullint4;
using omni::diskutil::DrivePathType;
using omni::diskutil::getDrivePathType;
using omni::diskutil::osSeparator;
using omni::errors::rtWarn;
using omni::hpc_random::max_xo_long_jumps;
using omni::hpc_random::seedXoroshiro128p;
using omni::hpc_random::seedXoshiro256pp;
using omni::random::Ran2Generator;
using omni::random::Xoroshiro128pGenerator;
using omni::random::Xoshiro256ppGenerator;
using namespace omni::math;
using namespace omni::hpc_math;
using namespace omni::testing;

//-------------------------------------------------------------------------------------------------
// Load a vector of random numbers using a CPU-based Ran2 generator.
//
// Arguments:
//   va:    The Hybrid object to load (modified and returned)
//   prng:  Ran2 generator object
//-------------------------------------------------------------------------------------------------
void loadRan2ByCPU(Hybrid<double> *va, Ran2Generator *prng) {
  double* vdata = va->data();
  const int nval = va->size();
  for (int i = 0; i < nval; i++) {
    vdata[i] = (prng->uniformRandomNumber() - 0.5) * sqrt(static_cast<double>(i + 1));
  }
}

//-------------------------------------------------------------------------------------------------
// Load a vector of random numbers using a CPU-based xoroshiro128+ generator.
//
// Arguments:
//   va:    The Hybrid object to load (modified and returned)
//   prng:  Xoroshiro128+ generator object
//-------------------------------------------------------------------------------------------------
void loadXoroshiro128pByCPU(Hybrid<double> *va, Xoroshiro128pGenerator *prng) {
  double* vdata = va->data();
  const int nval = va->size();
  for (int i = 0; i < nval; i++) {
    vdata[i] = (prng->uniformRandomNumber() - 0.5) * sqrt(static_cast<double>(i + 1));
  }
}

//-------------------------------------------------------------------------------------------------
// Load a vector of random numbers using a CPU-based xoshiro256++ generator.
//
// Arguments:
//   va:    The Hybrid object to load (modified and returned)
//   prng:  Xoroshiro128+ generator object
//-------------------------------------------------------------------------------------------------
void loadXoshiro256ppByCPU(Hybrid<double> *va, Xoshiro256ppGenerator *prng) {
  double* vdata = va->data();
  const int nval = va->size();
  for (int i = 0; i < nval; i++) {
    vdata[i] = (prng->uniformRandomNumber() - 0.5) * sqrt(static_cast<double>(i + 1));
  }
}

//-------------------------------------------------------------------------------------------------
// Evaluate random numbers from a plethora of xoroshiro128+ generator states.
//
// Arguments:
//   state_vector:  The vector of all states
//   n_generators:  Total number of states
//   samples:       Vector of states and iteration counts to watch out for
//   random_ouput:  Random numbers produced by the requested states and iterations
//   n_samples:     The number of samples (length of samples and random_output)
//-------------------------------------------------------------------------------------------------
__global__ void __launch_bounds__(large_block_size, 1)
kEvalXoroshiro128p(ullint2* state_vector, const int n_generators, const int2* samples,
                   double* random_output, const int n_samples) {
  __shared__ volatile int smax_iter;
  
  // Find the maximum iteration and cache the samples
  if (threadIdx.x == 0) {
    smax_iter = 0;
  }
  int max_iter = 0;
  for (int pos = threadIdx.x; pos < n_samples; pos += blockDim.x) {
    const int2 tmp_sample = samples[pos];
    if (tmp_sample.y > max_iter) {
      max_iter = tmp_sample.y;
    }
  }
#ifdef OMNI_USE_HIP
  const int init_shfl_stride = 32;
#else
  const int init_shfl_stride = 16;
#endif
  for (int i = init_shfl_stride; i > 0; i >>= 1) {
    const int max_comp = SHFL_DOWN(max_iter, i);
    if (max_comp > max_iter) {
      max_iter = max_comp;
    }
  }
  __syncthreads();
  const int lane_idx = (threadIdx.x & warp_bits_mask_int);
  if (lane_idx == 0 && max_iter > 0) {
    atomicMax((int*)&smax_iter, max_iter);
  }
  __syncthreads();
  const int ginc = blockDim.x * gridDim.x;
  for (int pos = threadIdx.x + (blockIdx.x * blockDim.x); pos < n_generators; pos += ginc) {
    ullint2 tmp_state = state_vector[pos];
    for (int i = 0; i <= smax_iter; i++) {

      // Get a uniform random number in the range [0, 1).  Casting the unsigned long long int to
      // a (signed) long long int works as anything great than 2^3 just becomes negative.  The
      // bit string is unchanged.
      const ullint s0 = tmp_state.x;
      ullint       s1 = tmp_state.y;
      const ullint rndbits = s0 + s1;
      const llint work = (((rndbits >> 12) & 0xfffffffffffff) | 0x3ff0000000000000);
      const double rn_out = __longlong_as_double(work) - 1.0;

      // Horribly inefficient loop and access pattern, but this is just
      // a test program.  All of this is L1-cached, too.
      for (int j = 0; j < n_samples; j++) {
        const int2 tmp_sample = samples[j];
        if (pos == tmp_sample.x && i == tmp_sample.y) {
          random_output[j] = rn_out;
        }
      }

      // Push the state forward
      s1 ^= s0;
      tmp_state.x = (((s0 << 24) | (s0 >> (64 - 24))) ^ s1 ^ (s1 << 16));
      tmp_state.y =  ((s1 << 37) | (s1 >> (64 - 37)));      
    }

    // In a simulation, the state vector would be updated after computations
    // for this thread or atom are complete.
    state_vector[pos] = tmp_state;    
  }
}

//-------------------------------------------------------------------------------------------------
// Evaluate random numbers from a plethora of xoshiro256++ generator states.
//
// Arguments:
//   state_vector:  The vector of all states
//   n_generators:  Total number of states
//   samples:       Vector of states and iteration counts to watch out for
//   random_ouput:  Random numbers produced by the requested states and iterations
//   n_samples:     The number of samples (length of samples and random_output)
//-------------------------------------------------------------------------------------------------
__global__ void __launch_bounds__(large_block_size, 1)
kEvalXoshiro256pp(ullint4* state_vector, const int n_generators, const int2* samples,
                  double* random_output, const int n_samples) {
  __shared__ volatile int smax_iter;
  
  // Find the maximum iteration and cache the samples
  if (threadIdx.x == 0) {
    smax_iter = 0;
  }
  int max_iter = 0;
  for (int pos = threadIdx.x; pos < n_samples; pos += blockDim.x) {
    const int2 tmp_sample = samples[pos];
    if (tmp_sample.y > max_iter) {
      max_iter = tmp_sample.y;
    }
  }
#ifdef OMNI_USE_HIP
  const int init_shfl_stride = 32;
#else
  const int init_shfl_stride = 16;
#endif
  for (int i = init_shfl_stride; i > 0; i >>= 1) {
    const int max_comp = SHFL_DOWN(max_iter, i);
    if (max_comp > max_iter) {
      max_iter = max_comp;
    }
  }
  __syncthreads();
  const int lane_idx = (threadIdx.x & warp_bits_mask_int);
  if (lane_idx == 0 && max_iter > 0) {
    atomicMax((int*)&smax_iter, max_iter);
  }
  __syncthreads();
  const int ginc = blockDim.x * gridDim.x;
  for (int pos = threadIdx.x + (blockIdx.x * blockDim.x); pos < n_generators; pos += ginc) {
    ullint4 tmp_state = state_vector[pos];
    for (int i = 0; i <= smax_iter; i++) {

      // Get a uniform random number in the range [0, 1).  Casting the unsigned long long int to
      // a (signed) long long int works as anything great than 2^3 just becomes negative.  The
      // bit string is unchanged.
      const ullint sxsw = tmp_state.x + tmp_state.w;
      const ullint rndbits = tmp_state.x + ((sxsw << 23) | (sxsw >> (64 - 23)));
      const llint work = (((rndbits >> 12) & 0xfffffffffffff) | 0x3ff0000000000000);
      const double rn_out = __longlong_as_double(work) - 1.0;

      // Horribly inefficient loop and access pattern, but this is just
      // a test program.  All of this is L1-cached, too.
      for (int j = 0; j < n_samples; j++) {
        const int2 tmp_sample = samples[j];
        if (pos == tmp_sample.x && i == tmp_sample.y) {
          random_output[j] = rn_out;
        }
      }

      // Push the state forward
      const ullint t = (tmp_state.y << 17);
      tmp_state.z ^= tmp_state.x;
      tmp_state.w ^= tmp_state.y;
      tmp_state.y ^= tmp_state.z;
      tmp_state.x ^= tmp_state.w;
      tmp_state.z ^= t;
      tmp_state.w = ((tmp_state.w << 45) | (tmp_state.w >> (64 - 45)));
    }

    // In a simulation, the state vector would be updated after computations
    // for this thread or atom are complete.
    state_vector[pos] = tmp_state;    
  }
}

//-------------------------------------------------------------------------------------------------
// Reproduce the random number generated by one of a series of xoroshiro128+ state vectors.
//
// Arguments:
//   rng_states:     List of xoroshiro128+ long-jump generator states (the dimension of this array
//                   implies the number of long-jumps taken when seeding the various states)
//   generator_idx:  Index of the generator state to query.  In a simulation, this might correspond
//                   to the atom index, or perhaps to the thread index within a particular launch
//                   grid.
//   iteration:      Produce the pseudo-random number for this point in the sequence of the
//                   particular atom or thread.
//-------------------------------------------------------------------------------------------------
double pinpointXoroshiro128p(std::vector<ullint2> &rng_states, const int generator_idx,
                             const uint iteration) {
  
  // Determine the initial generator state, possibly after having taken some short jumps.
  // The rng_states vector was created by taking up to max_xo_long_jumps long jumps of a single
  // state.  A given generator state is then created by tiling this series of initial, long jump
  // states with 1, 2, ..., n additional jumps, up to the total number of generators needed.  The
  // total number of generators is irrelevant.  This calculation requires the index for just one.
  const int n_seeds       = rng_states.size();
  const int n_short_jumps = generator_idx / n_seeds;
  const int seed_idx      = generator_idx - (n_short_jumps * n_seeds);

  // Recover the generator initial state of interest.  In a simulation, an initial generator X can
  // be jumped forward by P long jumps and N short jumps to arrive at a subsidiary generator X'.
  // Advancing X' by K iterations will produce the same output random number as advancing X by K
  // iterations, then jumping forward by P long jumps and N short jumps.
  Xoroshiro128pGenerator tgen(rng_states[seed_idx]);
  for (int i = 0; i < n_short_jumps; i++) {
    tgen.jump();
  }
  
  // Get the double-precision random number resulting from the requested iteration
  double result;
  for (uint i = 0; i <= iteration; i++) {
    result = tgen.uniformRandomNumber();
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
// Reproduce the random number generated by one of a series of xoshiro256++ state vectors.
//
// Arguments:
//   rng_states:     List of xoshiro256++ long-jump generator states (the dimension of this array
//                   implies the number of long-jumps taken when seeding the various states)
//   generator_idx:  Index of the generator state to query.  In a simulation, this might correspond
//                   to the atom index, or perhaps to the thread index within a particular launch
//                   grid.
//   iteration:      Produce the pseudo-random number for this point in the sequence of the
//                   particular atom or thread.
//-------------------------------------------------------------------------------------------------
double pinpointXoshiro256pp(std::vector<ullint4> &rng_states, const int generator_idx,
                            const uint iteration) {
  
  // Determine the initial generator state, possibly after having taken some short jumps.
  // The rng_states vector was created by taking up to max_xo_long_jumps long jumps of a single
  // state.  A given generator state is then created by tiling this series of initial, long jump
  // states with 1, 2, ..., n additional jumps, up to the total number of generators needed.  The
  // total number of generators is irrelevant.  This calculation requires the index for just one.
  const int n_seeds       = rng_states.size();
  const int n_short_jumps = generator_idx / n_seeds;
  const int seed_idx      = generator_idx - (n_short_jumps * n_seeds);

  // Recover the generator initial state of interest.  In a simulation, an initial generator X can
  // be jumped forward by P long jumps and N short jumps to arrive at a subsidiary generator X'.
  // Advancing X' by K iterations will produce the same output random number as advancing X by K
  // iterations, then jumping forward by P long jumps and N short jumps.
  Xoshiro256ppGenerator tgen(rng_states[seed_idx]);
  for (int i = 0; i < n_short_jumps; i++) {
    tgen.jump();
  }
  
  // Get the double-precision random number resulting from the requested iteration
  double result;
  for (uint i = 0; i <= iteration; i++) {
    result = tgen.uniformRandomNumber();
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
// main
//-------------------------------------------------------------------------------------------------
int main(int argc, char* argv[]) {

  // Some baseline initialization
  TestEnvironment oe(argc, argv);
  HpcConfig gpu_config(ExceptionResponse::WARN);
  std::vector<int> my_gpus = gpu_config.getGpuDevice(1);
  GpuDetails gpu = gpu_config.getGpuInfo(my_gpus[0]);
  
  // Section 1
  section("Vector processing capabilities");
  
  // Section 2
  section("GPU-based Xoroshiro128+ PRNG");
  
  // Perform a summation over a double-precision real vector using the GPU
  section(1);
  const int n_tiny  = 16;
  const int n_small = 128;
  const int n_chunk = 517;
  const int n_block = 1024;
  const int n_large = 23552;
  const int n_giant = 2500000;
  Hybrid<double> tiny_set(n_tiny, "tiny_vector");
  Hybrid<double> small_set(n_small, "small_vector");
  Hybrid<double> chunk_set(n_chunk, "chunk_vector");
  Hybrid<double> block_set(n_block, "block_vector");
  Hybrid<double> large_set(n_large, "large_vector");
  Hybrid<double> giant_set(n_giant, "giant_vector");
  Hybrid<double> tb_buffer(gpu.getSMPCount(), "sum_accumulators", HybridFormat::HOST_ONLY);
  Ran2Generator prng(oe.getRandomSeed());
  loadRan2ByCPU(&tiny_set, &prng);
  loadRan2ByCPU(&small_set, &prng);
  loadRan2ByCPU(&chunk_set, &prng);
  loadRan2ByCPU(&block_set, &prng);
  tiny_set.upload();
  loadRan2ByCPU(&tiny_set, &prng);
  const double gpu_tiny_sum  = sum(tiny_set, &tb_buffer, gpu);
  const double cpu_tiny_sum  = sum(tiny_set, &tb_buffer, gpu, HybridTargetLevel::HOST);
  check(gpu_tiny_sum, RelationalOperator::EQUAL, Approx(-0.4835516929).margin(1.0e-8),
        "The tiniest vector did not sum correctly on the GPU.\n");
  check(cpu_tiny_sum, RelationalOperator::EQUAL, Approx(-1.1269689474).margin(1.0e-8),
        "The tiniest vector did not sum correctly on the CPU.\n");
  small_set.upload();
  chunk_set.upload();
  block_set.upload();
  check(sum(small_set, &tb_buffer, gpu), RelationalOperator::EQUAL, sum<double>(small_set),
        "The small vector did not sum correctly on the GPU.\n");
  check(sum(chunk_set, &tb_buffer, gpu), RelationalOperator::EQUAL, sum<double>(chunk_set),
        "The medium, odd-sized vector did not sum correctly on the GPU.\n");
  check(sum(block_set, &tb_buffer, gpu), RelationalOperator::EQUAL, sum<double>(block_set),
        "The full block-sized vector did not sum correctly on the GPU.\n");
  Xoroshiro128pGenerator fast_prng(78172);
  loadXoroshiro128pByCPU(&large_set, &fast_prng);
  loadXoroshiro128pByCPU(&giant_set, &fast_prng);
  large_set.upload();
  giant_set.upload();
  check(sum(large_set, &tb_buffer, gpu), RelationalOperator::EQUAL, sum<double>(large_set),
        "The large vector did not sum correctly on the GPU.\n");
  check(sum(giant_set, &tb_buffer, gpu), RelationalOperator::EQUAL, sum<double>(giant_set),
        "The giant vector did not sum correctly on the GPU.\n");
  
  // Test the GPU-base random number seeding and synchronized CPU/GPU generation
  const int n_cellulose_atoms = 408609;
  Hybrid<ullint2> rng128p_states(n_cellulose_atoms, "xoroshiro128p_state");
  Hybrid<ullint4> rng256pp_states(n_cellulose_atoms, "xoroshiro256pp_state");
  seedXoroshiro128p(&rng128p_states, 8773925, gpu);
  seedXoshiro256pp(&rng256pp_states, 4091832, gpu);
  rng128p_states.download();
  rng256pp_states.download();
  Xoroshiro128pGenerator xrs128p_check(8773925);
  Xoshiro256ppGenerator xrs256pp_check(4091832);
  const int n_seeds_made = std::min(max_xo_long_jumps, n_cellulose_atoms);
  std::vector<ullint2> cpu_128p_seeds(n_seeds_made);
  std::vector<ullint4> cpu_256pp_seeds(n_seeds_made);
  for (int i = 0; i < n_seeds_made; i++) {
    cpu_128p_seeds[i] = xrs128p_check.revealState();
    cpu_256pp_seeds[i] = xrs256pp_check.revealState();
    xrs128p_check.longJump();
    xrs256pp_check.longJump();
  }

  // Create a smattering of generator indices (within the bounds of the rng128p_states above) and
  // some (low) iteration counts at which to test each of them.  Predict the results on the CPU,
  // then compute them on the GPU.
  const int n_samples = 16;
  Hybrid<int2> samples(n_samples, "assorted_points");
  int2* samp_ptr = samples.data();
  for (int i = 0; i < n_samples; i++) {
    samp_ptr[i] = { 918 * i, (2 * i) + 5 };
  }
  Hybrid<double> random_output(n_samples, "random_pluckings");
  double* ro_ptr = random_output.data();
  for (int i = 0; i < n_samples; i++) {
    ro_ptr[i] = pinpointXoroshiro128p(cpu_128p_seeds, samp_ptr[i].x, samp_ptr[i].y);
  }
  samples.upload();
  const int nsmp = gpu.getSMPCount();
  const int nthr = gpu.getMaxThreadsPerBlock();
  kEvalXoroshiro128p<<<nsmp, nthr>>>(rng128p_states.data(HybridTargetLevel::DEVICE),
                                     n_cellulose_atoms, samples.data(HybridTargetLevel::DEVICE),
                                     random_output.data(HybridTargetLevel::DEVICE), n_samples);
  check(random_output.readHost(), RelationalOperator::EQUAL, random_output.readDevice(), "Random "
        "numbers from an array of Xoroshiro128+ generators computed on the CPU and GPU do not "
        "agree.");
  kEvalXoshiro256pp<<<nsmp, nthr>>>(rng256pp_states.data(HybridTargetLevel::DEVICE),
                                    n_cellulose_atoms, samples.data(HybridTargetLevel::DEVICE),
                                    random_output.data(HybridTargetLevel::DEVICE), n_samples);
  for (int i = 0; i < n_samples; i++) {
    ro_ptr[i] = pinpointXoshiro256pp(cpu_256pp_seeds, samp_ptr[i].x, samp_ptr[i].y);
  }
  check(random_output.readHost(), RelationalOperator::EQUAL, random_output.readDevice(), "Random "
        "numbers from an array of Xoshiro256++ generators computed on the CPU and GPU do not "
        "agree.");
  
  // Print results
  printTestSummary(oe.getVerbosity());
  
  return 0;
}
