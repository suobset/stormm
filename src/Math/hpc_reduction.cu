// -*-c++-*-
#include "hpc_reduction.cuh"

namespace omni {
namespace synthesis {

//-------------------------------------------------------------------------------------------------
__global__ void __launch_bounds__(small_block_size, 4)
kConjGradGather(ReductionKit redk, ConjGradKit cgk) {

  // Two small arrays will store the double-precision accumulants for squared force magnitudes
  // and force differential measurements.
  __shared__ double gg_collector[small_block_size >> warp_bits];
  __shared__ double dgg_collector[small_block_size >> warp_bits];
  int rdwu_idx = blockIdx.x;
  while (rdwu_idx < redk.nrdwu) {
    const int wabs_pos   = (rdwu_idx * rdwu_abstract_length);
    const int start_pos  = redk.rdwu_abstracts[wabs_pos + (int)(RdwuAbstractMap::ATOM_START)];
    const int end_pos    = redk.rdwu_abstracts[wabs_pos + (int)(RdwuAbstractMap::ATOM_END)];
    const int result_pos = redk.rdwu_abstracts[wabs_pos + (int)(RdwuAbstractMap::RESULT_INDEX)];

    // A for-loop manages the reading for this work unit, as the per-thread workload is consistent,
    // whereas a while loop still manages the work unit progression as the sizes of different work
    // units may not be consistent.  The advancement through work units is therefore asynchronous.
    for (int i = start_pos; i < end_pos; i += blockDim.x) {
const double dfx = cgk.
    }
  }
}

} // namespace synthesis
} // namespace omni
