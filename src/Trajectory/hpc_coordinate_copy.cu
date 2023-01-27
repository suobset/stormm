// -*-c++-*-
#include "copyright.h"
#include "Constants/hpc_bounds.h"
#include "DataTypes/common_types.h"
#include "coordinate_series.h"
#include "hpc_coordinate_copy.h"

namespace stormm {
namespace trajectory {

#include "Math/rounding.cui"

//-------------------------------------------------------------------------------------------------
__global__ void __launch_bounds__(medium_block_size, 2)
kCopyCoordinateXYZ(void* xdest, void* ydest, void* zdest, const double dest_scale,
                   const size_t dest_offset, const size_t ct_dest, const void* xorig,
                   const void* yorig, const void* zorig, const double orig_scale,
                   const size_t orig_offset, const size_t ct_orig, const int natom,
                   const ValidCoordinateTypes vct) {
  __shared__ double *shd_xdest, *shd_ydest, *shd_zdest, *shd_xorig, *shd_yorig, *shd_zorig;
  __shared__ float  *shf_xdest, *shf_ydest, *shf_zdest, *shf_xorig, *shf_yorig, *shf_zorig;
  __shared__ short  *shs_xdest, *shs_ydest, *shs_zdest, *shs_xorig, *shs_yorig, *shs_zorig;
  __shared__ int    *shi_xdest, *shi_ydest, *shi_zdest, *shi_xorig, *shi_yorig, *shi_zorig;
  __shared__ llint  *shl_xdest, *shl_ydest, *shl_zdest, *shl_xorig, *shl_yorig, *shl_zorig;

  // Use each warp's first thread to set a different collection of the __shared__ pointers
  if (threadIdx.x == 0) {
    shd_xdest = (double*)xdest;
    shd_ydest = (double*)ydest;
    shd_zdest = (double*)zdest;
    shd_xorig = (double*)xorig;
    shd_yorig = (double*)yorig;
    shd_zorig = (double*)zorig;
  }
  else if (threadIdx.x == warp_size_int) {
    shf_xdest = (float*)xdest;
    shf_ydest = (float*)ydest;
    shf_zdest = (float*)zdest;
    shf_xorig = (float*)xorig;
    shf_yorig = (float*)yorig;
    shf_zorig = (float*)zorig;
  }
  else if (threadIdx.x == twice_warp_size_int) {
    shs_xdest = (short int*)xdest;
    shs_ydest = (short int*)ydest;
    shs_zdest = (short int*)zdest;
    shs_xorig = (short int*)xorig;
    shs_yorig = (short int*)yorig;
    shs_zorig = (short int*)zorig;
  }
  else if (threadIdx.x == 3 * warp_size_int) {
    shi_xdest = (int*)xdest;
    shi_ydest = (int*)ydest;
    shi_zdest = (int*)zdest;
    shi_xorig = (int*)xorig;
    shi_yorig = (int*)yorig;
    shi_zorig = (int*)zorig;
  }
  else if (threadIdx.x == 4 * warp_size_int) {
    shl_xdest = (llint*)xdest;
    shl_ydest = (llint*)ydest;
    shl_zdest = (llint*)zdest;
    shl_xorig = (llint*)xorig;
    shl_yorig = (llint*)yorig;
    shl_zorig = (llint*)zorig;
  }
  __syncthreads();

  // For transfers between long long ints, the double type may not be sufficient to preserve all
  // bits.  Otherwise, use double-precision numbers as the medium of exchange.
  const size_t padded_natom = devcRoundUp(natom, warp_size_int);
  const int gdim = gridDim.x * blockDim.x;
  const double inv_orig_scale = 1.0 / orig_scale;
  if (ct_dest == vct.llint_id && ct_orig == vct.llint_id) {
    const bool gets_bigger = (dest_scale > orig_scale);
    const llint conv_factor = (gets_bigger) ? __double2ll_rn(dest_scale / orig_scale) :
                                              __double2ll_rn(orig_scale / dest_scale);
    for (int i = threadIdx.x + (blockIdx.x * blockDim.x); i < natom; i += gdim) {
      const size_t orig_idx = orig_offset + (size_t)(i);
      const size_t dest_idx = dest_offset + (size_t)(i);
      if (gets_bigger) {
        shl_xdest[dest_idx] = shl_xorig[orig_idx] * conv_factor;
        shl_ydest[dest_idx] = shl_yorig[orig_idx] * conv_factor;
        shl_zdest[dest_idx] = shl_zorig[orig_idx] * conv_factor;
      }
      else {
        shl_xdest[dest_idx] = shl_xorig[orig_idx] / conv_factor;
        shl_ydest[dest_idx] = shl_yorig[orig_idx] / conv_factor;
        shl_zdest[dest_idx] = shl_zorig[orig_idx] / conv_factor;
      }
    }
  }
  else {
    for (int i = threadIdx.x + (blockIdx.x * blockDim.x); i < natom; i += gdim) {
      const size_t orig_idx = orig_offset + (size_t)(i);
      const size_t dest_idx = dest_offset + (size_t)(i);
      double tmp_xorig, tmp_yorig, tmp_zorig;
      if (ct_orig == vct.double_id) {
        tmp_xorig = (double)(shd_xorig[orig_idx]) * inv_orig_scale;
        tmp_yorig = (double)(shd_yorig[orig_idx]) * inv_orig_scale;
        tmp_zorig = (double)(shd_zorig[orig_idx]) * inv_orig_scale;
      }
      else if (ct_orig == vct.float_id) {
        tmp_xorig = (double)(shf_xorig[orig_idx]) * inv_orig_scale;
        tmp_yorig = (double)(shf_yorig[orig_idx]) * inv_orig_scale;
        tmp_zorig = (double)(shf_zorig[orig_idx]) * inv_orig_scale;
      }
      else if (ct_orig == vct.short_id) {
        tmp_xorig = (double)(shs_xorig[orig_idx]) * inv_orig_scale;
        tmp_yorig = (double)(shs_yorig[orig_idx]) * inv_orig_scale;
        tmp_zorig = (double)(shs_zorig[orig_idx]) * inv_orig_scale;
      }
      else if (ct_orig == vct.int_id) {
        tmp_xorig = (double)(shi_xorig[orig_idx]) * inv_orig_scale;
        tmp_yorig = (double)(shi_yorig[orig_idx]) * inv_orig_scale;
        tmp_zorig = (double)(shi_zorig[orig_idx]) * inv_orig_scale;
      }
      else if (ct_orig == vct.llint_id) {
        tmp_xorig = (double)(shl_xorig[orig_idx]) * inv_orig_scale;
        tmp_yorig = (double)(shl_yorig[orig_idx]) * inv_orig_scale;
        tmp_zorig = (double)(shl_zorig[orig_idx]) * inv_orig_scale;
      }
      if (ct_dest == vct.double_id) {
        shd_xdest[dest_idx] = tmp_xorig;
        shd_ydest[dest_idx] = tmp_yorig;
        shd_zdest[dest_idx] = tmp_zorig;
      }
      else if (ct_dest == vct.float_id) {
        shf_xdest[dest_idx] = tmp_xorig;
        shf_ydest[dest_idx] = tmp_yorig;
        shf_zdest[dest_idx] = tmp_zorig;
      }
      else if (ct_dest == vct.short_id) {
        shs_xdest[dest_idx] = __double2int_rn(tmp_xorig * dest_scale);
        shs_ydest[dest_idx] = __double2int_rn(tmp_yorig * dest_scale);
        shs_zdest[dest_idx] = __double2int_rn(tmp_zorig * dest_scale);
      }
      else if (ct_dest == vct.int_id) {
        shi_xdest[dest_idx] = __double2int_rn(tmp_xorig * dest_scale);
        shi_ydest[dest_idx] = __double2int_rn(tmp_yorig * dest_scale);
        shi_zdest[dest_idx] = __double2int_rn(tmp_zorig * dest_scale);
      }
      else if (ct_dest == vct.llint_id) {
        shl_xdest[dest_idx] = __double2ll_rn(tmp_xorig * dest_scale);
        shl_ydest[dest_idx] = __double2ll_rn(tmp_yorig * dest_scale);
        shl_zdest[dest_idx] = __double2ll_rn(tmp_zorig * dest_scale);
      }
    }
  }
}
  
//-------------------------------------------------------------------------------------------------
void launchCopyCoordinateXYZ(void* xdest, void* ydest, void* zdest, const double dest_scale,
		             const size_t dest_offset, const size_t ct_dest, const void* xorig,
                             const void* yorig, const void* zorig, const double orig_scale,
                             const size_t orig_offset, const size_t ct_orig, const int natom,
                             const GpuDetails &gpu) {
  const ValidCoordinateTypes vct(double_type_index, float_type_index, short_type_index,
                                 int_type_index, llint_type_index);
  const int nblock = 2 * gpu.getSMPCount();
  kCopyCoordinateXYZ<<<nblock, medium_block_size>>>(xdest, ydest, zdest, dest_scale, dest_offset,
                                                    ct_dest, xorig, yorig, zorig, orig_scale,
                                                    orig_offset, ct_orig, natom, vct);
}
  
} // namespace trajectory
} // namespace stormm
