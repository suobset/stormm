// -*-c++-*-
#include "Constants/hpc_bounds.h"
#include "Constants/scaling.h"
#include "DataTypes/common_types.h"
#include "hpc_scorecard.h"

namespace stormm {
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
    if (warp_idx == 0 && blockIdx.x == 0) {
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
extern void launchScoreCardInitialization(const StateVariable var, const int system_index,
                                          llint* accumulators, const int system_count,
                                          const GpuDetails &gpu) {
  const ullint state_mask = (var == StateVariable::ALL_STATES) ? 0xffffffffffffffffLLU :
                                                                 (0x1 << static_cast<int>(var));
  kInitializeScoreCard<<<gpu.getSMPCount(), large_block_size>>>(state_mask, accumulators,
                                                                system_index, system_count);
}

//-------------------------------------------------------------------------------------------------
extern void launchScoreCardInitialization(const std::vector<StateVariable> &var,
                                          const int system_index, llint* accumulators,
                                          const int system_count, const GpuDetails &gpu) {
  const size_t nvar = var.size();
  ullint state_mask = 0LLU;
  for (size_t i = 0; i < nvar; i++) {
    state_mask |= (0x1 << static_cast<int>(var[i]));
  }
  kInitializeScoreCard<<<gpu.getSMPCount(), large_block_size>>>(state_mask, accumulators,
                                                                system_index, system_count);
}

//-------------------------------------------------------------------------------------------------
__global__ void __launch_bounds__(large_block_size, 1)
kScoreCardCommit(const ullint state_mask, const int system_index, const int system_count,
                 const size_t sample_count, const llint* inst_acc, double* run_acc,
                 double* sqd_acc, llint* time_ser_acc) {

  // Get the extent of the StateVariable enumerator
  const int max_states = (int)(StateVariable::ALL_STATES);
  const int padded_max_states = (((max_states + warp_size_int - 1) >> warp_bits) << warp_bits);
  const int warp_idx = (threadIdx.x >> warp_bits);
  const int lane_idx = (threadIdx.x & warp_bits_mask_int);
  if (system_index >= 0) {
    if (warp_idx == 0 && blockIdx.x == 0) {
      for (int i = lane_idx; i < max_states; i += warp_size_int) {
        if ((state_mask >> i) & 0x1LLU) {
          const int size_t slot = (padded_max_states * system_index) + i;
          const llint ll_amt = accumulators[slot];
          const double d_amt = (double)(ll_amt);
          run_acc[slot] += d_amt;
          sqd_acc[slot] += d_amt * d_amt;
          const size_t ts_slot = slot + (static_cast<size_t>(padded_max_states * system_count) *
                                         sample_count);
          time_ser_acc[ts_slot] = ll_amt;
        }
      }
    }    
  }
  else {
    int syspos = warp_idx + (blockIdx.x * (blockDim.x >> warp_bits));
    while (syspos < system_count) {
      for (int i = lane_idx; i < max_states; i += warp_size_int) {
        if ((state_mask >> i) & 0x1LLU) {
          const int size_t slot = (padded_max_states * syspos) + i;
          const llint ll_amt = accumulators[slot];
          const double d_amt = (double)(ll_amt);
          run_acc[slot] += d_amt;
          sqd_acc[slot] += d_amt * d_amt;
          const size_t ts_slot = slot + (static_cast<size_t>(padded_max_states * system_count) *
                                         sample_count);
          time_ser_acc[ts_slot] = ll_amt;
        }
      }
      syspos += gridDim.x * (blockDim.x >> warp_bits);
    }
  }
}

//-------------------------------------------------------------------------------------------------
extern void launchScoreCardCommit(std::vector<StateVariable> &var, const int system_index,
                                  const int system_count, const size_t sample_count,
                                  llint* inst_acc, double* run_acc, double* sqd_acc,
                                  llint* time_ser_acc, const GpuDetails &gpu) {
  const size_t nvar = var.size();
  ullint state_mask = 0LLU;
  for (size_t i = 0; i < nvar; i++) {
    state_mask |= (0x1 << static_cast<int>(var[i]));
  }
  const HybridTargetLevel tier = HybridTargetLevel::DEVICE;
  kScoreCardCommit<<<gpu.getSMPCount(), large_block_size>>>(state_mask, system_index, system_count,
                                                            sample_count, inst_acc, run_acc,
                                                            sqd_acc, time_ser_acc);
}

//-------------------------------------------------------------------------------------------------
extern void launchScoreCardCommit(StateVariable var, const int system_index,
                                  const int system_count, const size_t sample_count,
                                  llint* inst_acc, double* run_acc, double* sqd_acc,
                                  llint* time_ser_acc, const GpuDetails &gpu) {
  const ullint state_mask = (var == StateVariable::ALL_STATES) ? 0xffffffffffffffffLLU :
                                                                 (0x1 << static_cast<int>(var));
  const HybridTargetLevel tier = HybridTargetLevel::DEVICE;
  kScoreCardCommit<<<gpu.getSMPCount(), large_block_size>>>(state_mask, system_index, system_count,
                                                            sample_count, inst_acc, run_acc,
                                                            sqd_acc, time_ser_acc);
}

} // namespace energy
} // namespace stormm
