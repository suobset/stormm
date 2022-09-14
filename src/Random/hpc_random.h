// -*-c++-*-
#ifndef STORMM_HPC_RANDOM_H
#define STORMM_HPC_RANDOM_H

#include "copyright.h"
#include "Accelerator/gpu_details.h"
#include "Accelerator/hybrid.h"
#include "DataTypes/stormm_vector_types.h"

namespace stormm {
namespace random {

using card::GpuDetails;
using card::Hybrid;

/// \brief Launch the eponymous kernel to seed Xoroshiro128+ random number generator states.  The
///        algorithm's long jump is executed using CPU code as this work must be serialized.
///
/// \param state_vector  Blank vector allocated to hold the requested number of random number
///                      generator states.  Filled and returned.
/// \param igseed        Random number seed for the first of the generator states (all others
///                      follow from this)
/// \param scrub_cycles  Number of times to run the random number generator and discard results in
///                      order to raise the quality of random numbers coming out of it
/// \param gpu           Details of the GPU on which to launch this operation
void initXoroshiro128pArray(Hybrid<ullint2> *state_vector, int igseed, int scrub_cycles,
                            const GpuDetails &gpu);

/// \brief Launch the eponymous kernel to seed Xoshiro256++ random number generator states.  The
///        algorithm's long jump is executed using CPU code as this work must be serialized.
///
/// \param state_xy      Blank vector allocated to hold the first half of all random number
///                      generator states, in the requested number.  Filled and returned.
/// \param state_zw      Blank vector allocated to hold the second half of all random number
///                      generator states, in the requested number.  Filled and returned.
/// \param igseed        Random number seed for the first of the generator states (all others
///                      follow from this)
/// \param scrub_cycles  Number of times to run the random number generator and discard results in
///                      order to raise the quality of random numbers coming out of it
/// \param gpu           Details of the GPU on which to launch this operation
void initXoshiro256ppArray(Hybrid<ullint2> *state_xy, Hybrid<ullint2> *state_zw, int igseed,
                           int scrub_cycles, const GpuDetails &gpu);

} // namespace random
} // namespace stormm

#endif
