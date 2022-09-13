// -*-c++-*-
#ifndef STORMM_HPC_RANDOM_H
#define STORMM_HPC_RANDOM_H

#include "copyright.h"
#include "Accelerator/gpu_details.h"
#include "Accelerator/hybrid.h"
#include "DataTypes/stormm_vector_types.h"

namespace stormm {
namespace hpc_random {

using card::GpuDetails;
using card::Hybrid;

/// \brief The maximum number of long jumps to take with any of the scrambled linear pseudo-random
///        number generators.  The CPU will take the long jumps, while individual GPU threads make
///        standard jumps ahead in the pseudo-random sequence.  The CPU is 10-30 times faster
///        than a single GPU thread, so this number could in fact scale with the number of state
///        vectors being requested, but the optimization would be on the order of milliseconds
///        and making this a constant makes it easier to reproduce the sequence at a given point.
constexpr int max_xo_long_jumps = 1024;
  
/// \brief Launch the above kernel to seed Xoroshiro128+ random number generator states.  The
///        algorithm's long jump is executed using CPU code as this work must be serialized.
///
/// \param state_vector  Blank vector allocated to hold the requested number of random number
///                      generator states.  Filled and returned.
/// \param igseed        Random number seed for the first of the generator states (all others
///                      follow from this)
/// \param gpu           Details of the GPU on which to launch this operation
void seedXoroshiro128p(Hybrid<ullint2> *state_vector, int igseed, const GpuDetails &gpu);

/// \brief Launch the above kernel to seed Xoshiro256++ random number generator states.  The
///        algorithm's long jump is executed using CPU code as this work must be serialized.
///
/// \param state_vector  Blank vector allocated to hold the requested number of random number
///                      generator states.  Filled and returned.
/// \param igseed        Random number seed for the first of the generator states (all others
///                      follow from this)
/// \param gpu           Details of the GPU on which to launch this operation
void seedXoshiro256pp(Hybrid<ullint4> *state_vector, int igseed, const GpuDetails &gpu);

} // namespace hpc_random
} // namespace stormm

#endif
