// -*-c++-*-
#ifndef OMNI_HPC_SCORECARD_H
#define OMNI_HPC_SCORECARD_H

#include <vector>
#include "Accelerator/gpu_details.h"
#include "energy_enumerators.h"

namespace omni {
namespace energy {

using card::GpuDetails;

/// \brief Launch the appropriate kernel to initialize some or all of a scorecard on the device.
///
/// Overloaded:
///   - Accept a single state variable
///   - Accept a list of state variables
///
/// \param var           One or more state variables to initialize
/// \param system_index  Index of the system to initialize (-1 will initialize all systems)
/// \param accumulators  Pointer to the energy tracking object's instantaneous accumulators on
///                      the HPC device
/// \param system_count  Total number of systems traked by the ScoreCard (energy tracking object)
/// \{
void launchScoreCardInitialization(StateVariable var, int system_index, llint* accumulators,
                                   int system_count, const GpuDetails &gpu);

void launchScoreCardInitialization(const std::vector<StateVariable> &var, int system_index,
                                   llint* accumulators, int system_count, const GpuDetails &gpu);
/// \}

} // namespace energy
} // namespace omni

#endif
