// -*-c++-*-
#include "copyright.h"
#include "Constants/hpc_bounds.h"
#include "hpc_random.cuh"
#include "random.h"

namespace stormm {
namespace hpc_random {

using card::HybridFormat;
using card::HybridTargetLevel;
using random::Xoroshiro128pGenerator;
using random::Xoshiro256ppGenerator;
using random::xrs128p_jump_i;
using random::xrs128p_jump_ii;
using random::xrs256pp_jump_i;
using random::xrs256pp_jump_ii;
using random::xrs256pp_jump_iii;
using random::xrs256pp_jump_iv;
  
//-------------------------------------------------------------------------------------------------
__global__ void __launch_bounds__(large_block_size, 1)
kSeedXoroshiro128p(ullint2* state_vector, const int n_seeds, const int n_generators) {

  // Useful constants
  const ullint2 jump = { xrs128p_jump_i, xrs128p_jump_ii };
  
  // The only threads that seed additional generators are those that can take a seed
  int seed_pos = threadIdx.x + (blockIdx.x * blockDim.x);
  if (seed_pos < n_seeds) {

    // Read the seed from the array
    ullint2 my_state = state_vector[seed_pos];
    for (int gen_pos = seed_pos + n_seeds; gen_pos < n_generators; gen_pos += n_seeds) {

      // Execute the Xoroshiro128+ jump function on this seed.
      ullint s0 = 0LLU;
      ullint s1 = 0LLU;
      for (int b = 0; b < 64; b++) {
        if (jump.x & (0x1LLU << b)) {
          s0 ^= my_state.x;
          s1 ^= my_state.y;
        }
        const ullint ns0 = my_state.x;
        const ullint ns1 = my_state.y ^ ns0;
        my_state.x = (((ns0 << 24) | (ns0 >> (64 - 24))) ^ ns1 ^ (ns1 << 16));
        my_state.y =  ((ns1 << 37) | (ns1 >> (64 - 37)));
      }
      for (int b = 0; b < 64; b++) {
        if (jump.y & (0x1LLU << b)) {
          s0 ^= my_state.x;
          s1 ^= my_state.y;
        }
        const ullint ns0 = my_state.x;
        const ullint ns1 = my_state.y ^ ns0;
        my_state.x = (((ns0 << 24) | (ns0 >> (64 - 24))) ^ ns1 ^ (ns1 << 16));
        my_state.y =  ((ns1 << 37) | (ns1 >> (64 - 37)));
      }
      my_state.x = s0;
      my_state.y = s1;
      state_vector[gen_pos] = my_state;
    }
  }
}

//-------------------------------------------------------------------------------------------------
extern void seedXoroshiro128p(Hybrid<ullint2> *state_vector, const int igseed,
                              const GpuDetails &gpu) {

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
  Xoroshiro128pGenerator prng(igseed);
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
  kSeedXoroshiro128p<<<nsmp, nthr>>>(state_vector->data(HybridTargetLevel::DEVICE), n_seeds,
                                     n_generators);
}

//-------------------------------------------------------------------------------------------------
__global__ void __launch_bounds__(large_block_size, 1)
kSeedXoshiro256pp(ullint4* state_vector, const int n_seeds, const int n_generators) {

  // Useful constants
  const ullint4 jump = { xrs256pp_jump_i, xrs256pp_jump_ii, xrs256pp_jump_iii, xrs256pp_jump_iv };
  
  // The only threads that seed additional generators are those that can take a seed
  int seed_pos = threadIdx.x + (blockIdx.x * blockDim.x);
  if (seed_pos < n_seeds) {

    // Read the seed from the array
    ullint4 my_state = state_vector[seed_pos];
    for (int gen_pos = seed_pos + n_seeds; gen_pos < n_generators; gen_pos += n_seeds) {

      // Execute the Xoroshiro256++ jump function on this seed.
      ullint s0 = 0LLU;
      ullint s1 = 0LLU;
      ullint s2 = 0LLU;
      ullint s3 = 0LLU;
      for (int b = 0; b < 64; b++) {
        if (jump.x & (0x1LLU << b)) {
          s0 ^= my_state.x;
          s1 ^= my_state.y;
          s2 ^= my_state.z;
          s3 ^= my_state.w;
        }
        const ullint t = (my_state.y << 17);
        my_state.z ^= my_state.x;
        my_state.w ^= my_state.y;
        my_state.y ^= my_state.z;
        my_state.x ^= my_state.w;
        my_state.z ^= t;
        my_state.w = ((my_state.w << 45) | (my_state.w >> (64 - 45)));
      }
      for (int b = 0; b < 64; b++) {
        if (jump.y & (0x1LLU << b)) {
          s0 ^= my_state.x;
          s1 ^= my_state.y;
          s2 ^= my_state.z;
          s3 ^= my_state.w;
        }
        const ullint t = (my_state.y << 17);
        my_state.z ^= my_state.x;
        my_state.w ^= my_state.y;
        my_state.y ^= my_state.z;
        my_state.x ^= my_state.w;
        my_state.z ^= t;
        my_state.w = ((my_state.w << 45) | (my_state.w >> (64 - 45)));
      }
      for (int b = 0; b < 64; b++) {
        if (jump.z & (0x1LLU << b)) {
          s0 ^= my_state.x;
          s1 ^= my_state.y;
          s2 ^= my_state.z;
          s3 ^= my_state.w;
        }
        const ullint t = (my_state.y << 17);
        my_state.z ^= my_state.x;
        my_state.w ^= my_state.y;
        my_state.y ^= my_state.z;
        my_state.x ^= my_state.w;
        my_state.z ^= t;
        my_state.w = ((my_state.w << 45) | (my_state.w >> (64 - 45)));
      }
      for (int b = 0; b < 64; b++) {
        if (jump.w & (0x1LLU << b)) {
          s0 ^= my_state.x;
          s1 ^= my_state.y;
          s2 ^= my_state.z;
          s3 ^= my_state.w;
        }
        const ullint t = (my_state.y << 17);
        my_state.z ^= my_state.x;
        my_state.w ^= my_state.y;
        my_state.y ^= my_state.z;
        my_state.x ^= my_state.w;
        my_state.z ^= t;
        my_state.w = ((my_state.w << 45) | (my_state.w >> (64 - 45)));
      }
      my_state.x = s0;
      my_state.y = s1;
      my_state.z = s2;
      my_state.w = s3;
      state_vector[gen_pos] = my_state;
    }
  }
}

//-------------------------------------------------------------------------------------------------
extern void seedXoshiro256pp(Hybrid<ullint4> *state_vector, const int igseed,
                             const GpuDetails &gpu) {

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
  Xoshiro256ppGenerator prng(igseed);
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
  kSeedXoshiro256pp<<<nsmp, nthr>>>(state_vector->data(HybridTargetLevel::DEVICE), n_seeds,
                                    n_generators);
}

} // namespace hpc_random
} // namespace stormm
