// -*-c++-*-
#ifndef STORMM_HPC_RANDOM_CUH
#define STORMM_HPC_RANDOM_CUH

#include "Accelerator/gpu_details.h"
#include "Accelerator/hybrid.h"
#include "DataTypes/common_types.h"
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
  
/// \brief Seeds a specified number of Xoroshiro128+ random number generator states in an array
///        using device resources to accelerate the production.  A preset number of non-zero
///        states created by the algorithm's long jump method are propagated with the algorithm's
///        short jump method.
///
/// \param state_vector  As input, contains the states of long jump seeds.  Modified and output
///                      containing n_generators states set at regular, well-spaced intervals along
///                      the algorithm's sequence of (2^128) - 1 random numbers.
/// \param n_seeds       The initial number of seeds in state_vector
/// \param n_generators  The requested number of distinct pseudo-random number generators
__global__ void __launch_bounds__(large_block_size, 1)
kSeedXoroshiro128p(ullint2* state_vector, const int n_seeds, const int n_generators);

/// \brief Launch the above kernel to seed Xoroshiro128+ random number generator states.  The
///        algorithm's long jump is executed using CPU code as this work must be serialized.
///
/// \param state_vector  Blank vector allocated to hold the requested number of random number
///                      generator states.  Filled and returned.
/// \param igseed        Random number seed for the first of the generator states (all others
///                      follow from this)
/// \param gpu           Details of the GPU on which to launch this operation
void seedXoroshiro128p(Hybrid<ullint2> *state_vector, int igseed, const GpuDetails &gpu);

/// \brief Seeds a specified number of Xoshiro256+ random number generator states in an array
///        using device resources to accelerate the production.  A preset number of non-zero
///        states created by the algorithm's long jump method are propagated with the algorithm's
///        short jump method.
///
/// \param state_vector  As input, contains the states of long jump seeds.  Modified and output
///                      containing n_generators states set at regular, well-spaced intervals along
///                      the algorithm's sequence of (2^128) - 1 random numbers.
/// \param n_seeds       The initial number of seeds in state_vector
/// \param n_generators  The requested number of distinct pseudo-random number generators
__global__ void __launch_bounds__(large_block_size, 1)
kSeedXoshiro256pp(ullint4* state_vector, const int n_seeds, const int n_generators);

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
