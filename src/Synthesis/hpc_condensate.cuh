// -*-c++-*-
#ifndef STORMM_HPC_CONDENSATE_CUH
#define STORMM_HPC_CONDENSATE_CUH

#include "copyright.h"
#include "Constants/hpc_bounds.h"
#include "Trajectory/coordinate_series.h"
#include "condensate.h"

namespace stormm {
namespace synthesis {

using trajectory::CoordinateSeries;
using trajectory::CoordinateSeriesReader;

//-------------------------------------------------------------------------------------------------
// Transfer data from a CoordinateSeries object to the COndensate's coordinate arrays.
//
// Arguments:
//   cdw:  Writeable abstract for the Condensate object
//   csr:  Read-only abstract for the source coordinate series
//-------------------------------------------------------------------------------------------------
template <typename T>
__global__ void __launch_bounds__(large_block_size, 1)
kCondensateUpdate(CondensateWriter cdw, const CoordinateSeriesReader<T> csr) {    
  const size_t nframe_zu = csr.nframe - 1;

  // Explicitly type out the rounding function to avoid including it in a header-level
  // environment, where it could cause a mess for the compiler.
  const size_t natom_zu  = ((csr.natom + warp_bits_mask_int) / warp_size_int) * warp_size_int;
  const size_t padded_atoms = nframe_zu * natom_zu;
  const size_t istride = gridDim.x * blockDim.x;
  const float inv_gpos_scale_f = csr.inv_gpos_scale;
  for (size_t i = threadIdx.x + (blockIdx.x * blockDim.x); i < padded_atoms; i += istride) {
    switch (cdw.mode) {
    case CondensationLevel::DOUBLE:
      cdw.xcrd[i] = (double)(csr.xcrd[i]) * csr.inv_gpos_scale;
      cdw.ycrd[i] = (double)(csr.ycrd[i]) * csr.inv_gpos_scale;
      cdw.zcrd[i] = (double)(csr.zcrd[i]) * csr.inv_gpos_scale;
      break;
    case CondensationLevel::FLOAT:
      cdw.xcrd_sp[i] = (float)(csr.xcrd[i]) * inv_gpos_scale_f;
      cdw.ycrd_sp[i] = (float)(csr.ycrd[i]) * inv_gpos_scale_f;
      cdw.zcrd_sp[i] = (float)(csr.zcrd[i]) * inv_gpos_scale_f;
      break;
    }
  }
}

} // namespace synthesis
} // namespace stormm

#endif
