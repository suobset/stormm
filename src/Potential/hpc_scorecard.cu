// -*-c++-*-
#include "Constants/hpc_bounds.h"
#include "Constants/scaling.h"
#include "DataTypes/common_types.h"
#include "hpc_scorecard.h"

namespace omni {
namespace energy {

//-------------------------------------------------------------------------------------------------
__global__ void __launch_bounds__(large_block_size, 1)
kInitializeScoreCard(const ullint state_mask, llint* accumulators, const int system_index,
                     const int system_count) {

  // Get the extent of the StateVariable enumerator
  const int max_states = (int)(StateVariable::ALL_STATES);
  const int padded_max_states = (((max_states + warp_size_int - 1) >> warp_bits) << warp_bits);
  const int warp_idx = (threadIdx.x >> warp_bits);
  const int lane_idx = (threadIdx.x & warp_bits_mask_int);
  if (system_index >= 0) {
    if (warp_idx == 0) {
      for (int i = lane_idx; i < max_states; i += warp_size_int) {
        if ((state_mask >> i) & 0x1LLU) {
          accumulators[(padded_max_states * system_index) + i] = 0LL;
        }
      }
    }
  }
  else {
    int syspos = warp_idx + (blockIdx.x * (blockDim.x >> warp_bits));
    while (syspos < system_count) {
      for (int i = lane_idx; i < max_states; i += warp_size_int) {
        if ((state_mask >> i) & 0x1LLU) {
          accumulators[(padded_max_states * syspos) + i] = 0LL;
        }
      }
      syspos += gridDim.x * (blockDim.x >> warp_bits);
    }
  }
}

//-------------------------------------------------------------------------------------------------
void launchScoreCardInitialization(const StateVariable var, const int system_index,
                                   llint* accumulators, const int system_count,
                                   const GpuDetails &gpu) {
  const ullint state_mask = (var == StateVariable::ALL_STATES) ? 0xffffffffffffffffLLU :
                                                                 (0x1 << static_cast<int>(var));
  kInitializeScoreCard<<<gpu.getSMPCount(), large_block_size>>>(state_mask, accumulators,
                                                                system_index, system_count);
}

//-------------------------------------------------------------------------------------------------
void launchScoreCardInitialization(const std::vector<StateVariable> &var, const int system_index,
                                   llint* accumulators, const int system_count,
                                   const GpuDetails &gpu) {
  const size_t nvar = var.size();
  ullint state_mask = 0LLU;
  for (size_t i = 0; i < nvar; i++) {
    state_mask |= (0x1 << static_cast<int>(var[i]));
  }
  kInitializeScoreCard<<<gpu.getSMPCount(), large_block_size>>>(state_mask, accumulators,
                                                                system_index, system_count);
}

} // namespace energy
} // namespace omni
