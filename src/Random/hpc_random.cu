// -*-c++-*-
#include "copyright.h"
#include "Constants/hpc_bounds.h"
#include "hpc_random.h"
#include "random.h"

namespace stormm {
namespace random {

using card::HybridFormat;
using card::HybridTargetLevel;
  
#include "xor_shift_rng.cui"

//-------------------------------------------------------------------------------------------------
__global__ void __launch_bounds__(large_block_size, 1)
kInitXoroshiro128pArray(ullint2* state_vector, const int n_seeds, const int n_generators) {
  
  // The only threads that seed additional generators are those that can take a seed
  int seed_pos = threadIdx.x + (blockIdx.x * blockDim.x);
  if (seed_pos < n_seeds) {

    // Read the seed from the array
    ullint2 my_state = state_vector[seed_pos];
    state_vector[seed_pos] = my_state;
    for (int gen_pos = seed_pos + n_seeds; gen_pos < n_generators; gen_pos += n_seeds) {
      my_state = xoroshiro128p_jump(my_state);
      state_vector[gen_pos] = my_state;
    }
  }
}

//-------------------------------------------------------------------------------------------------
extern void initXoroshiro128pArray(Hybrid<ullint2> *state_vector, const int igseed,
                                   const int scrub_cycles, const GpuDetails &gpu) {

  // Sanity check on the type of hybrid object
  ullint2* svdata;
  std::vector<ullint2> staging_space;
  const int n_generators = state_vector->size();
  const int n_seeds = std::min(max_xo_long_jumps, n_generators);
  const HybridFormat svfmt = state_vector->getFormat();
  switch (svfmt) {
  case HybridFormat::EXPEDITED:
  case HybridFormat::DECOUPLED:
  case HybridFormat::UNIFIED:
  case HybridFormat::HOST_ONLY:
    svdata = state_vector->data();
    break;
  case HybridFormat::DEVICE_ONLY:

    // If there is no allocated host data, create a vector for staging the work
    staging_space.resize(n_seeds);
    svdata = staging_space.data();
    break;
  }
  Xoroshiro128pGenerator prng(igseed, scrub_cycles);
  svdata[0] = prng.revealState();
  for (int i = 1; i < n_seeds; i++) {
    prng.longJump();
    svdata[i] = prng.revealState();
  }
  switch (svfmt) {
  case HybridFormat::EXPEDITED:
  case HybridFormat::DECOUPLED:
    state_vector->upload();
    break;
  case HybridFormat::UNIFIED:
  case HybridFormat::HOST_ONLY:
    break;
  case HybridFormat::DEVICE_ONLY:
    state_vector->putDevice(staging_space);
    break;
  }
  const int nsmp = gpu.getSMPCount();
  const int nthr = gpu.getMaxThreadsPerBlock();
  kInitXoroshiro128pArray<<<nsmp, nthr>>>(state_vector->data(HybridTargetLevel::DEVICE), n_seeds,
                                          n_generators);
  switch (svfmt) {
  case HybridFormat::EXPEDITED:
  case HybridFormat::DECOUPLED:
    state_vector->download();
    break;
  case HybridFormat::UNIFIED:
  case HybridFormat::HOST_ONLY:
  case HybridFormat::DEVICE_ONLY:
    break;
  }
}

//-------------------------------------------------------------------------------------------------
__global__ void __launch_bounds__(large_block_size, 1)
kInitXoshiro256ppArray(ullint4* state_vector, const int n_seeds, const int n_generators) {
  
  // The only threads that seed additional generators are those that can take a seed
  int seed_pos = threadIdx.x + (blockIdx.x * blockDim.x);
  if (seed_pos < n_seeds) {
    
    // Read the seed from the array
    ullint4 my_state = state_vector[seed_pos];
    state_vector[seed_pos] = my_state;
    for (int gen_pos = seed_pos + n_seeds; gen_pos < n_generators; gen_pos += n_seeds) {
      my_state = xoshiro256pp_jump(my_state);
      state_vector[gen_pos] = my_state;
    }
  }
}

//-------------------------------------------------------------------------------------------------
extern void initXoshiro256ppArray(Hybrid<ullint4> *state_vector, const int igseed,
                                  const int scrub_cycles, const GpuDetails &gpu) {

  // Sanity check on the type of hybrid object
  ullint4* svdata;
  std::vector<ullint4> staging_space;
  const int n_generators = state_vector->size();
  const int n_seeds = std::min(max_xo_long_jumps, n_generators);
  const HybridFormat svfmt = state_vector->getFormat();
  switch (svfmt) {
  case HybridFormat::EXPEDITED:
  case HybridFormat::DECOUPLED:
  case HybridFormat::UNIFIED:
  case HybridFormat::HOST_ONLY:
    svdata = state_vector->data();
    break;
  case HybridFormat::DEVICE_ONLY:

    // If there is no allocated host data, create a vector for staging the work
    staging_space.resize(n_seeds);
    svdata = staging_space.data();
    break;
  }
  Xoshiro256ppGenerator prng(igseed, scrub_cycles);
  svdata[0] = prng.revealState();
  for (int i = 1; i < n_seeds; i++) {
    prng.longJump();
    svdata[i] = prng.revealState();
  }
  switch (svfmt) {
  case HybridFormat::EXPEDITED:
  case HybridFormat::DECOUPLED:
    state_vector->upload();
    break;
  case HybridFormat::UNIFIED:
  case HybridFormat::HOST_ONLY:
    break;
  case HybridFormat::DEVICE_ONLY:
    state_vector->putDevice(staging_space);
    break;
  }
  const int nsmp = gpu.getSMPCount();
  const int nthr = gpu.getMaxThreadsPerBlock();
  kInitXoshiro256ppArray<<<nsmp, nthr>>>(state_vector->data(HybridTargetLevel::DEVICE), n_seeds,
                                         n_generators);
  switch (svfmt) {
  case HybridFormat::EXPEDITED:
  case HybridFormat::DECOUPLED:
    state_vector->download();
    break;
  case HybridFormat::UNIFIED:
  case HybridFormat::HOST_ONLY:
  case HybridFormat::DEVICE_ONLY:
    break;
  }
}

} // namespace random
} // namespace stormm
