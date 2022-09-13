// -*-c++-*-
#ifndef STORMM_HPC_RANDOM_CUH
#define STORMM_HPC_RANDOM_CUH

#include "copyright.h"
#include "Constants/hpc_bounds.h"
#include "DataTypes/stormm_vector_types.h"

namespace stormm {
namespace hpc_random {

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

} // namespace hpc_random
} // namespace stormm

#endif
