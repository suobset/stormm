// -*-c++-*-
#include "copyright.h"
#include "Constants/hpc_bounds.h"
#include "DataTypes/common_types.h"
#include "DataTypes/stormm_vector_types.h"
#include "coordinate_series.h"
#include "hpc_coordinate_copy.cuh"
#include "hpc_coordinate_copy.h"

namespace stormm {
namespace trajectory {

#include "Math/rounding.cui"
  
//-------------------------------------------------------------------------------------------------
__global__ void __launch_bounds__(medium_block_size, 2)
kCopyCoordinateXYZ(void* xdest, void* ydest, void* zdest, const double dest_scale,
                   const size_t ct_dest, const void* xorig, const void* yorig, const void* zorig,
                   const double orig_scale, const size_t ct_orig, const int natom,
                   double* umat_dest, double* invu_dest, double* bdim_dest,
                   const double* umat_orig, const double* invu_orig, const double* bdim_orig, 
                   const ValidCoordinateTypes vct, const int dest_atom_offset,
                   const int orig_atom_offset, const int dest_xfrm_offset,
                   const int orig_xfrm_offset, const int dest_bdim_offset,
                   const int orig_bdim_offset) {
  __shared__ double *shd_xdest, *shd_ydest, *shd_zdest, *shd_xorig, *shd_yorig, *shd_zorig;
  __shared__ float  *shf_xdest, *shf_ydest, *shf_zdest, *shf_xorig, *shf_yorig, *shf_zorig;
  __shared__ short  *shs_xdest, *shs_ydest, *shs_zdest, *shs_xorig, *shs_yorig, *shs_zorig;
  __shared__ int    *shi_xdest, *shi_ydest, *shi_zdest, *shi_xorig, *shi_yorig, *shi_zorig;
  __shared__ llint  *shl_xdest, *shl_ydest, *shl_zdest, *shl_xorig, *shl_yorig, *shl_zorig;

  // Use each warp's first thread to set a different collection of the __shared__ pointers.  This
  // will cast away the const-ness of the origin coordinate arrays, but no changes will be imparted
  // to that data by this kernel.
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

  // Copy the box information, if provided.  This assumes that there are at least eight warps in
  // the block, which is true for all known architectures.
  if (umat_dest != nullptr && blockIdx.x == 0) {
    const int warp_idx = (threadIdx.x >> warp_bits);
    if (warp_idx == 5) {
      int pos = (threadIdx.x & warp_bits_mask_int);
      while (pos < 9) {
        umat_dest[dest_xfrm_offset + pos] = umat_orig[orig_xfrm_offset + pos];
        pos += warp_size_int;
      }
    }
    else if (warp_idx == 6) {
      int pos = (threadIdx.x & warp_bits_mask_int);
      while (pos < 9) {
        invu_dest[dest_xfrm_offset + pos] = invu_orig[orig_xfrm_offset + pos];
        pos += warp_size_int;
      }
    }
    else if (warp_idx == 7) {
      int pos = (threadIdx.x & warp_bits_mask_int);
      while (pos < 6) {
        bdim_dest[dest_bdim_offset + pos] = bdim_orig[orig_bdim_offset + pos];
        pos += warp_size_int;
      }
    }
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
      const size_t orig_idx = orig_atom_offset + i;
      const size_t dest_idx = dest_atom_offset + i;
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
      const size_t orig_idx = orig_atom_offset + i;
      const size_t dest_idx = dest_atom_offset + i;
      double tmp_xorig, tmp_yorig, tmp_zorig;
      if (ct_orig == vct.double_id) {
        tmp_xorig = shd_xorig[orig_idx];
        tmp_yorig = shd_yorig[orig_idx];
        tmp_zorig = shd_zorig[orig_idx];
      }
      else if (ct_orig == vct.float_id) {
        tmp_xorig = (double)(shf_xorig[orig_idx]);
        tmp_yorig = (double)(shf_yorig[orig_idx]);
        tmp_zorig = (double)(shf_zorig[orig_idx]);
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
                             const size_t ct_dest, const void* xorig, const void* yorig,
                             const void* zorig, const double orig_scale, const size_t ct_orig,
                             const int natom, const int dest_atom_offset,
                             const int orig_atom_offset, const GpuDetails &gpu) {
  const ValidCoordinateTypes vct(double_type_index, float_type_index, short_type_index,
                                 int_type_index, llint_type_index);
  const int nblock = 2 * gpu.getSMPCount();
  kCopyCoordinateXYZ<<<nblock, medium_block_size>>>(xdest, ydest, zdest, dest_scale, ct_dest,
                                                    xorig, yorig, zorig, orig_scale, ct_orig,
                                                    natom, nullptr, nullptr, nullptr, nullptr,
                                                    nullptr, nullptr, vct, dest_atom_offset,
                                                    orig_atom_offset, 0, 0, 0, 0);
}

//-------------------------------------------------------------------------------------------------
void launchCopyCoordinateXYZBox(void* xdest, void* ydest, void* zdest, const double dest_scale,
                                const size_t ct_dest, const void* xorig, const void* yorig,
                                const void* zorig, const double orig_scale, const size_t ct_orig,
                                const int natom, double* umat_dest, double* invu_dest,
                                double* bdim_dest, const double* umat_orig,
                                const double* invu_orig, const double* bdim_orig,
                                const int dest_atom_offset, const int orig_atom_offset,
                                const int dest_xfrm_offset, const int orig_xfrm_offset,
                                const int dest_bdim_offset, const int orig_bdim_offset,
                                const GpuDetails &gpu) {
  const ValidCoordinateTypes vct(double_type_index, float_type_index, short_type_index,
                                 int_type_index, llint_type_index);
  const int nblock = 2 * gpu.getSMPCount();
  kCopyCoordinateXYZ<<<nblock, medium_block_size>>>(xdest, ydest, zdest, dest_scale, ct_dest,
                                                    xorig, yorig, zorig, orig_scale, ct_orig,
                                                    natom, umat_dest, invu_dest, bdim_dest,
                                                    umat_orig, invu_orig, bdim_orig, vct,
                                                    dest_atom_offset, orig_atom_offset,
                                                    dest_xfrm_offset, orig_xfrm_offset,
                                                    dest_bdim_offset, orig_bdim_offset);
}

//-------------------------------------------------------------------------------------------------
__global__ void __launch_bounds__(medium_block_size, 2)
kCopyCoordinateXYZ(void* xdest, void* ydest, void* zdest, const double dest_scale,
                   const int dest_scale_bits, const size_t ct_dest, const llint* xorig,
                   const llint* yorig, const llint* zorig, const int* xorig_ovrf,
                   const int* yorig_ovrf, const int* zorig_ovrf, const double orig_scale,
                   const int orig_scale_bits, const int natom, double* umat_dest,
                   double* invu_dest, double* bdim_dest, const double* umat_orig,
                   const double* invu_orig, const double* bdim_orig,
                   const ValidCoordinateTypes vct, const int dest_atom_offset,
                   const int orig_atom_offset, const int dest_xfrm_offset,
                   const int orig_xfrm_offset, const int dest_bdim_offset,
                   const int orig_bdim_offset) {
  __shared__ double *shd_xdest, *shd_ydest, *shd_zdest;
  __shared__ float  *shf_xdest, *shf_ydest, *shf_zdest;
  __shared__ short  *shs_xdest, *shs_ydest, *shs_zdest;
  __shared__ int    *shi_xdest, *shi_ydest, *shi_zdest;
  __shared__ llint  *shl_xdest, *shl_ydest, *shl_zdest;

  // Use each warp's first thread to set a different collection of the __shared__ pointers.  This
  // will cast away the const-ness of the origin coordinate arrays, but no changes will be imparted
  // to that data by this kernel.
  if (threadIdx.x == 0) {
    shd_xdest = (double*)xdest;
    shd_ydest = (double*)ydest;
    shd_zdest = (double*)zdest;
  }
  else if (threadIdx.x == warp_size_int) {
    shf_xdest = (float*)xdest;
    shf_ydest = (float*)ydest;
    shf_zdest = (float*)zdest;
  }
  else if (threadIdx.x == twice_warp_size_int) {
    shs_xdest = (short int*)xdest;
    shs_ydest = (short int*)ydest;
    shs_zdest = (short int*)zdest;
  }
  else if (threadIdx.x == 3 * warp_size_int) {
    shi_xdest = (int*)xdest;
    shi_ydest = (int*)ydest;
    shi_zdest = (int*)zdest;
  }
  else if (threadIdx.x == 4 * warp_size_int) {
    shl_xdest = (llint*)xdest;
    shl_ydest = (llint*)ydest;
    shl_zdest = (llint*)zdest;
  }
  
  // Copy the box information, if provided.  This assumes that there are at least eight warps in
  // the block, which is true for all known architectures.
  if (umat_dest != nullptr && blockIdx.x == 0) {
    const int warp_idx = (threadIdx.x >> warp_bits);
    if (warp_idx == 5) {
      int pos = (threadIdx.x & warp_bits_mask_int);
      while (pos < 9) {
        umat_dest[dest_xfrm_offset + pos] = umat_orig[orig_xfrm_offset + pos];
        pos += warp_size_int;
      }
    }
    else if (warp_idx == 6) {
      int pos = (threadIdx.x & warp_bits_mask_int);
      while (pos < 9) {
        invu_dest[dest_xfrm_offset + pos] = invu_orig[orig_xfrm_offset + pos];
        pos += warp_size_int;
      }
    }
    else if (warp_idx == 7) {
      int pos = (threadIdx.x & warp_bits_mask_int);
      while (pos < 6) {
        bdim_dest[dest_bdim_offset + pos] = bdim_orig[orig_bdim_offset + pos];
        pos += warp_size_int;
      }
    }
  }
  __syncthreads();

  // For transfers between long long ints, the double type may not be sufficient to preserve all
  // bits.  Otherwise, use double-precision numbers as the medium of exchange.
  const size_t padded_natom = devcRoundUp(natom, warp_size_int);
  const int gdim = gridDim.x * blockDim.x;
  if (ct_dest == vct.llint_id) {
    if (dest_scale >= orig_scale) {
      const llint conv_factor = dest_scale / orig_scale;
      for (int i = threadIdx.x + (blockIdx.x * blockDim.x); i < natom; i += gdim) {
        const size_t orig_idx = orig_atom_offset + i;
        const size_t dest_idx = dest_atom_offset + i;

        // A 64-bit signed integer destination cannot hold more than the 64-bit part of an int95_t
        // coordinate if the scaling of the destination is equal to or greater than the original.
        shl_xdest[dest_idx] = xorig[orig_idx] * conv_factor;
        shl_ydest[dest_idx] = yorig[orig_idx] * conv_factor;
        shl_zdest[dest_idx] = zorig[orig_idx] * conv_factor;
      }
    }
    else {
      const double dconv_factor = orig_scale / dest_scale;
      const int conv_bits = orig_scale_bits - dest_scale_bits;
      const double dovrf_factor = max_llint_accumulation / dconv_factor;
      for (int i = threadIdx.x + (blockIdx.x * blockDim.x); i < natom; i += gdim) {
        const size_t orig_idx = orig_atom_offset + i;
        const size_t dest_idx = dest_atom_offset + i;
        shl_xdest[dest_idx] = (xorig[orig_idx] >> conv_bits);
        shl_ydest[dest_idx] = (yorig[orig_idx] >> conv_bits);
        shl_zdest[dest_idx] = (zorig[orig_idx] >> conv_bits);
        shl_xdest[dest_idx] += __double2ll_rn((double)(xorig_ovrf[orig_idx]) * dovrf_factor);
        shl_ydest[dest_idx] += __double2ll_rn((double)(yorig_ovrf[orig_idx]) * dovrf_factor);
        shl_zdest[dest_idx] += __double2ll_rn((double)(zorig_ovrf[orig_idx]) * dovrf_factor);
      }
    }
  }
  else {
    const double inv_oscale = 1.0 / orig_scale;
    for (int i = threadIdx.x + (blockIdx.x * blockDim.x); i < natom; i += gdim) {
      const size_t orig_idx = orig_atom_offset + i;
      const size_t dest_idx = dest_atom_offset + i;
      const double tmp_xorig = int95ToDouble(xorig[orig_idx], xorig_ovrf[orig_idx]) * inv_oscale;
      const double tmp_yorig = int95ToDouble(yorig[orig_idx], yorig_ovrf[orig_idx]) * inv_oscale;
      const double tmp_zorig = int95ToDouble(zorig[orig_idx], zorig_ovrf[orig_idx]) * inv_oscale;
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
    }
  }
}

//-------------------------------------------------------------------------------------------------
void launchCopyCoordinateXYZ(void* xdest, void* ydest, void* zdest, const double dest_scale,
                             const int dest_scale_bits, const size_t ct_dest, const llint* xorig,
                             const llint* yorig, const llint* zorig, const int* xorig_ovrf,
                             const int* yorig_ovrf, const int* zorig_ovrf, const double orig_scale,
                             const int orig_scale_bits, const int natom,
                             const int dest_atom_offset, const int orig_atom_offset,
                             const GpuDetails &gpu) {
  const ValidCoordinateTypes vct(double_type_index, float_type_index, short_type_index,
                                 int_type_index, llint_type_index);
  const int nblock = 2 * gpu.getSMPCount();
  kCopyCoordinateXYZ<<<nblock, medium_block_size>>>(xdest, ydest, zdest, dest_scale,
                                                    dest_scale_bits, ct_dest, xorig, yorig, zorig,
                                                    xorig_ovrf, yorig_ovrf, zorig_ovrf, orig_scale,
                                                    orig_scale_bits, natom, nullptr, nullptr,
                                                    nullptr, nullptr, nullptr, nullptr, vct,
                                                    dest_atom_offset, orig_atom_offset, 0, 0,
                                                    0, 0);
}

//-------------------------------------------------------------------------------------------------
void launchCopyCoordinateXYZBox(void* xdest, void* ydest, void* zdest, const double dest_scale,
                                const int dest_scale_bits, const size_t ct_dest,
                                const llint* xorig, const llint* yorig, const llint* zorig,
                                const int* xorig_ovrf, const int* yorig_ovrf,
                                const int* zorig_ovrf, const double orig_scale,
                                const int orig_scale_bits, const int natom, double* umat_dest,
                                double* invu_dest, double* bdim_dest, const double* umat_orig,
                                const double* invu_orig, const double* bdim_orig,
                                const int dest_atom_offset, const int orig_atom_offset,
                                const int dest_xfrm_offset, const int orig_xfrm_offset,
                                const int dest_bdim_offset, const int orig_bdim_offset,
                                const GpuDetails &gpu) {
  const ValidCoordinateTypes vct(double_type_index, float_type_index, short_type_index,
                                 int_type_index, llint_type_index);
  const int nblock = 2 * gpu.getSMPCount();
  kCopyCoordinateXYZ<<<nblock, medium_block_size>>>(xdest, ydest, zdest, dest_scale,
                                                    dest_scale_bits, ct_dest, xorig, yorig, zorig,
                                                    xorig_ovrf, yorig_ovrf, zorig_ovrf, orig_scale,
                                                    orig_scale_bits, natom, umat_dest, invu_dest,
                                                    bdim_dest, umat_orig, invu_orig, bdim_orig,
                                                    vct, dest_atom_offset, orig_atom_offset,
                                                    dest_xfrm_offset, orig_xfrm_offset,
                                                    dest_bdim_offset, orig_bdim_offset);
}

//-------------------------------------------------------------------------------------------------
__global__ void __launch_bounds__(medium_block_size, 2)
kCopyCoordinateXYZ(llint* xdest, llint* ydest, llint* zdest, int* xdest_ovrf, int* ydest_ovrf,
                   int* zdest_ovrf, const double dest_scale, const int dest_scale_bits,
                   const void* xorig, const void* yorig, const void* zorig,
                   const double orig_scale, const int orig_scale_bits, const size_t ct_orig,
                   const int natom, double* umat_dest, double* invu_dest, double* bdim_dest,
                   const double* umat_orig, const double* invu_orig, const double* bdim_orig,
                   llint* boxvecs_dest, int* boxvec_ovrf_dest, const ValidCoordinateTypes vct,
                   const int dest_atom_offset, const int orig_atom_offset,
                   const int dest_xfrm_offset, const int orig_xfrm_offset,
                   const int dest_bdim_offset, const int orig_bdim_offset) {
  __shared__ double *shd_xorig, *shd_yorig, *shd_zorig;
  __shared__ float  *shf_xorig, *shf_yorig, *shf_zorig;
  __shared__ short  *shs_xorig, *shs_yorig, *shs_zorig;
  __shared__ int    *shi_xorig, *shi_yorig, *shi_zorig;
  __shared__ llint  *shl_xorig, *shl_yorig, *shl_zorig;

  // Use each warp's first thread to set a different collection of the __shared__ pointers.  This
  // will cast away the const-ness of the origin coordinate arrays, but no changes will be imparted
  // to that data by this kernel.
  if (threadIdx.x == 0) {
    shd_xorig = (double*)xorig;
    shd_yorig = (double*)yorig;
    shd_zorig = (double*)zorig;
  }
  else if (threadIdx.x == warp_size_int) {
    shf_xorig = (float*)xorig;
    shf_yorig = (float*)yorig;
    shf_zorig = (float*)zorig;
  }
  else if (threadIdx.x == twice_warp_size_int) {
    shs_xorig = (short int*)xorig;
    shs_yorig = (short int*)yorig;
    shs_zorig = (short int*)zorig;
  }
  else if (threadIdx.x == 3 * warp_size_int) {
    shi_xorig = (int*)xorig;
    shi_yorig = (int*)yorig;
    shi_zorig = (int*)zorig;
  }
  else if (threadIdx.x == 4 * warp_size_int) {
    shl_xorig = (llint*)xorig;
    shl_yorig = (llint*)yorig;
    shl_zorig = (llint*)zorig;
  }

  // Copy the box information, if provided.  This assumes that there are at least eight warps in
  // the block, which is true for all known architectures.
  if (umat_dest != nullptr && blockIdx.x == 0) {
    const int warp_idx = (threadIdx.x >> warp_bits);
    if (warp_idx == 5) {
      int pos = (threadIdx.x & warp_bits_mask_int);
      while (pos < 9) {
        umat_dest[dest_xfrm_offset + pos] = umat_orig[orig_xfrm_offset + pos];
        pos += warp_size_int;
      }
    }
    else if (warp_idx == 6) {
      int pos = (threadIdx.x & warp_bits_mask_int);
      while (pos < 9) {
        invu_dest[dest_xfrm_offset + pos] = invu_orig[orig_xfrm_offset + pos];

        // With the inverse transformation matrix at hand, compute the fixed-precision box vectors.
        const int95_t bv = doubleToInt95(invu_dest[orig_xfrm_offset + pos] * dest_scale);
        boxvecs_dest[dest_xfrm_offset + pos] = bv.x;
        boxvec_ovrf_dest[dest_xfrm_offset + pos] = bv.x;
        pos += warp_size_int;
      }
    }
    else if (warp_idx == 7) {
      int pos = (threadIdx.x & warp_bits_mask_int);
      while (pos < 6) {
        bdim_dest[dest_bdim_offset + pos] = bdim_orig[orig_xfrm_offset + pos];
        pos += warp_size_int;
      }
    }
  }
  __syncthreads();

  // For transfers between long long ints, the double type may not be sufficient to preserve all
  // bits.  Otherwise, use double-precision numbers as the medium of exchange.
  const size_t padded_natom = devcRoundUp(natom, warp_size_int);
  const int gdim = gridDim.x * blockDim.x;
  if (ct_orig == vct.llint_id) {
    if (dest_scale <= orig_scale) {
      const llint conv_bits = orig_scale_bits - dest_scale_bits;
      for (int i = threadIdx.x + (blockIdx.x * blockDim.x); i < natom; i += gdim) {
        const size_t orig_idx = orig_atom_offset + i;
        const size_t dest_idx = dest_atom_offset + i;

        // The 64-bit part of the int95_t format is sufficient to hold all of the information
        // coming from the original 64-bit integer representation.
        if (dest_scale_bits == orig_scale_bits) {
          xdest[dest_idx] = shl_xorig[orig_idx];
          ydest[dest_idx] = shl_yorig[orig_idx];
          zdest[dest_idx] = shl_zorig[orig_idx];
        }
        else {
          xdest[dest_idx] = shl_xorig[orig_idx] >> conv_bits;
          ydest[dest_idx] = shl_yorig[orig_idx] >> conv_bits;
          zdest[dest_idx] = shl_zorig[orig_idx] >> conv_bits;
        }
      }
    }
    else {
      for (int i = threadIdx.x + (blockIdx.x * blockDim.x); i < natom; i += gdim) {
        const size_t orig_idx = orig_atom_offset + i;
        const size_t dest_idx = dest_atom_offset + i;
        const int95_t xcurr = { shl_xorig[orig_idx], 0 };
        const int95_t ycurr = { shl_yorig[orig_idx], 0 };
        const int95_t zcurr = { shl_zorig[orig_idx], 0 };
        const int95_t xoutp = changeFPBits(xcurr, orig_scale_bits, dest_scale_bits);
        const int95_t youtp = changeFPBits(ycurr, orig_scale_bits, dest_scale_bits);
        const int95_t zoutp = changeFPBits(zcurr, orig_scale_bits, dest_scale_bits);
        xdest[dest_idx] = xoutp.x;
        ydest[dest_idx] = youtp.x;
        zdest[dest_idx] = zoutp.x;
        xdest_ovrf[dest_idx] = xoutp.y;
        ydest_ovrf[dest_idx] = youtp.y;
        zdest_ovrf[dest_idx] = zoutp.y;
      }
    }
  }
  else {
    const double dconv_factor = dest_scale / orig_scale;
    for (int i = threadIdx.x + (blockIdx.x * blockDim.x); i < natom; i += gdim) {
      const size_t orig_idx = orig_atom_offset + i;
      const size_t dest_idx = dest_atom_offset + i;
      double dx, dy, dz;
      if (ct_orig == vct.double_id) {
        dx = shd_xorig[orig_idx];
        dy = shd_yorig[orig_idx];
        dz = shd_zorig[orig_idx];
      }
      else if (ct_orig == vct.float_id) {
        dx = (double)(shf_xorig[orig_idx]);
        dy = (double)(shf_yorig[orig_idx]);
        dz = (double)(shf_zorig[orig_idx]);
      }
      else if (ct_orig == vct.short_id) {
        dx = (double)(shs_xorig[orig_idx]);
        dy = (double)(shs_yorig[orig_idx]);
        dz = (double)(shs_zorig[orig_idx]);
      }
      else if (ct_orig == vct.int_id) {
        dx = (double)(shi_xorig[orig_idx]);
        dy = (double)(shi_yorig[orig_idx]);
        dz = (double)(shi_zorig[orig_idx]);
      }
      const int95_t idx = doubleToInt95(dx * dconv_factor);
      const int95_t idy = doubleToInt95(dy * dconv_factor);
      const int95_t idz = doubleToInt95(dz * dconv_factor);
      xdest[dest_idx] = idx.x;
      ydest[dest_idx] = idy.x;
      zdest[dest_idx] = idz.x;
      xdest_ovrf[dest_idx] = idx.y;
      ydest_ovrf[dest_idx] = idy.y;
      zdest_ovrf[dest_idx] = idz.y;
    }
  }
}

//-------------------------------------------------------------------------------------------------
void launchCopyCoordinateXYZ(llint* xdest, llint* ydest, llint* zdest, int* xdest_ovrf,
                             int* ydest_ovrf, int* zdest_ovrf, const double dest_scale,
                             const int dest_scale_bits, const void* xorig, const void* yorig,
                             const void* zorig, const double orig_scale, const int orig_scale_bits,
                             const size_t ct_orig, const int natom, const int dest_atom_offset,
                             const int orig_atom_offset, const GpuDetails &gpu) {
  const ValidCoordinateTypes vct(double_type_index, float_type_index, short_type_index,
                                 int_type_index, llint_type_index);
  const int nblock = 2 * gpu.getSMPCount();
  kCopyCoordinateXYZ<<<nblock, medium_block_size>>>(xdest, ydest, zdest, xdest_ovrf, ydest_ovrf,
                                                    zdest_ovrf, dest_scale, dest_scale_bits,
                                                    xorig, yorig, zorig, orig_scale,
                                                    orig_scale_bits, ct_orig, natom, nullptr,
                                                    nullptr, nullptr, nullptr, nullptr, nullptr,
                                                    nullptr, nullptr, vct, dest_atom_offset,
                                                    orig_atom_offset, 0, 0, 0, 0);
}

//-------------------------------------------------------------------------------------------------
void launchCopyCoordinateXYZBox(llint* xdest, llint* ydest, llint* zdest, int* xdest_ovrf,
                                int* ydest_ovrf, int* zdest_ovrf, const double dest_scale,
                                const int dest_scale_bits, const void* xorig, const void* yorig,
                                const void* zorig, const double orig_scale,
                                const int orig_scale_bits, const size_t ct_orig, const int natom,
                                double* umat_dest, double* invu_dest, double* bdim_dest,
                                const double* umat_orig, const double* invu_orig,
                                const double* bdim_orig, llint* boxvecs_dest,
                                int* boxvec_ovrf_dest, const int dest_atom_offset,
                                const int orig_atom_offset, const int dest_xfrm_offset,
                                const int orig_xfrm_offset, const int dest_bdim_offset,
                                const int orig_bdim_offset, const GpuDetails &gpu) {
  const ValidCoordinateTypes vct(double_type_index, float_type_index, short_type_index,
                                 int_type_index, llint_type_index);
  const int nblock = 2 * gpu.getSMPCount();
  kCopyCoordinateXYZ<<<nblock, medium_block_size>>>(xdest, ydest, zdest, xdest_ovrf, ydest_ovrf,
                                                    zdest_ovrf, dest_scale, dest_scale_bits,
                                                    xorig, yorig, zorig, orig_scale,
                                                    orig_scale_bits, ct_orig, natom, umat_dest,
                                                    invu_dest, bdim_dest, umat_orig, invu_orig,
                                                    bdim_orig, boxvecs_dest, boxvec_ovrf_dest, vct,
                                                    dest_atom_offset, orig_atom_offset,
                                                    dest_xfrm_offset, orig_xfrm_offset,
                                                    dest_bdim_offset, orig_bdim_offset);
}

//-------------------------------------------------------------------------------------------------
__global__ void __launch_bounds__(medium_block_size, 2)
kCopyCoordinates(PhaseSpaceWriter dest, const PhaseSpaceReader orig) {

  // Copy the box information
  if (blockIdx.x == 0) {
    const int warp_idx = (threadIdx.x >> warp_bits);
    const int lane_idx = (threadIdx.x & warp_bits_mask_int);
    int i = lane_idx;
    if (warp_idx == 0) {
      while (i < 9) {
        dest.umat[i] = orig.umat[i];
        i += warp_size_int;
      }
    }
    else if (warp_idx == 1) {
      while (i < 9) {
        dest.invu[i] = orig.invu[i];
        i += warp_size_int;
      }
    }
    else if (warp_idx == 2) {
      while (i < 9) {
        dest.umat_alt[i] = orig.umat_alt[i];
        i += warp_size_int;
      }
    }
    else if (warp_idx == 3) {
      while (i < 9) {
        dest.invu_alt[i] = orig.invu_alt[i];
        i += warp_size_int;
      }
    }
    else if (warp_idx == 4) {
      while (i < 6) {
        dest.boxdim[i] = orig.boxdim[i];
        i += warp_size_int;
      }
    }
    else if (warp_idx == 5) {
      while (i < 6) {
        dest.boxdim_alt[i] = orig.boxdim_alt[i];
        i += warp_size_int;
      }
    }
  }
  size_t pos = threadIdx.x + (blockIdx.x * blockDim.x);
  const size_t stride = devcRoundUp(dest.natom, warp_size_int);
  const size_t advance = (blockDim.x * gridDim.x);
  pos = copyCoordinateSet(dest.xcrd, 0, 1.0, 0, orig.xcrd, 0, 1.0, 0, orig.natom, pos,  0, stride,
                          advance);
  pos = copyCoordinateSet(dest.ycrd, 0, 1.0, 0, orig.ycrd, 0, 1.0, 0, orig.natom, pos,  1, stride,
                          advance);
  pos = copyCoordinateSet(dest.zcrd, 0, 1.0, 0, orig.zcrd, 0, 1.0, 0, orig.natom, pos,  2, stride,
                          advance);
  pos = copyCoordinateSet(dest.xvel, 0, 1.0, 0, orig.xvel, 0, 1.0, 0, orig.natom, pos,  3, stride,
                          advance);
  pos = copyCoordinateSet(dest.yvel, 0, 1.0, 0, orig.yvel, 0, 1.0, 0, orig.natom, pos,  4, stride,
                          advance);
  pos = copyCoordinateSet(dest.zvel, 0, 1.0, 0, orig.zvel, 0, 1.0, 0, orig.natom, pos,  5, stride,
                          advance);
  pos = copyCoordinateSet(dest.xfrc, 0, 1.0, 0, orig.xfrc, 0, 1.0, 0, orig.natom, pos,  6, stride,
                          advance);
  pos = copyCoordinateSet(dest.yfrc, 0, 1.0, 0, orig.yfrc, 0, 1.0, 0, orig.natom, pos,  7, stride,
                          advance);
  pos = copyCoordinateSet(dest.zfrc, 0, 1.0, 0, orig.zfrc, 0, 1.0, 0, orig.natom, pos,  8, stride,
                          advance);
  pos = copyCoordinateSet(dest.xalt, 0, 1.0, 0, orig.xalt, 0, 1.0, 0, orig.natom, pos,  9, stride,
                          advance);
  pos = copyCoordinateSet(dest.yalt, 0, 1.0, 0, orig.yalt, 0, 1.0, 0, orig.natom, pos, 10, stride,
                          advance);
  pos = copyCoordinateSet(dest.zalt, 0, 1.0, 0, orig.zalt, 0, 1.0, 0, orig.natom, pos, 11, stride,
                          advance);
  pos = copyCoordinateSet(dest.vxalt, 0, 1.0, 0, orig.vxalt, 0, 1.0, 0, orig.natom, pos, 12,
                          stride, advance);
  pos = copyCoordinateSet(dest.vyalt, 0, 1.0, 0, orig.vyalt, 0, 1.0, 0, orig.natom, pos, 13,
                          stride, advance);
  pos = copyCoordinateSet(dest.vzalt, 0, 1.0, 0, orig.vzalt, 0, 1.0, 0, orig.natom, pos, 14,
                          stride, advance);
  pos = copyCoordinateSet(dest.fxalt, 0, 1.0, 0, orig.fxalt, 0, 1.0, 0, orig.natom, pos, 15,
                          stride, advance);
  pos = copyCoordinateSet(dest.fyalt, 0, 1.0, 0, orig.fyalt, 0, 1.0, 0, orig.natom, pos, 16,
                          stride, advance);
  pos = copyCoordinateSet(dest.fzalt, 0, 1.0, 0, orig.fzalt, 0, 1.0, 0, orig.natom, pos, 17,
                          stride, advance);
}

//-------------------------------------------------------------------------------------------------
void launchCopyCoordinates(PhaseSpaceWriter *destination, const PhaseSpaceReader &origin,
                           const GpuDetails &gpu) {
  const int nblock = 2 * gpu.getSMPCount();
  kCopyCoordinates<<<nblock, medium_block_size>>>(*destination, origin);
}

//-------------------------------------------------------------------------------------------------
__global__ void __launch_bounds__(medium_block_size, 2)
kCopyCoordinates(PhaseSpaceWriter dest, const PsSynthesisReader orig,
                 const size_t orig_atom_offset, const size_t orig_xfrm_offset,
                 const size_t orig_bdim_offset) {

  // Copy the box information
  if (blockIdx.x == 0) {
    const int warp_idx = (threadIdx.x >> warp_bits);
    const int lane_idx = (threadIdx.x & warp_bits_mask_int);
    size_t i = lane_idx;
    if (warp_idx == 0) {
      while (i < 9) {
        dest.umat[i] = orig.umat[orig_xfrm_offset + i];
        i += warp_size_zu;
      }
    }
    else if (warp_idx == 1) {
      while (i < 9) {
        dest.invu[i] = orig.invu[orig_xfrm_offset + i];
        i += warp_size_zu;
      }
    }
    else if (warp_idx == 2) {
      while (i < 9) {
        dest.umat_alt[i] = orig.umat_alt[orig_xfrm_offset + i];
        i += warp_size_zu;
      }
    }
    else if (warp_idx == 3) {
      while (i < 9) {
        dest.invu_alt[i] = orig.invu_alt[orig_xfrm_offset + i];
        i += warp_size_zu;
      }
    }
    else if (warp_idx == 4) {
      while (i < 6) {
        dest.boxdim[i] = orig.boxdims[orig_bdim_offset + i];
        i += warp_size_zu;
      }
    }
    else if (warp_idx == 5) {
      while (i < 6) {
        dest.boxdim_alt[i] = orig.alt_boxdims[orig_bdim_offset + i];
        i += warp_size_zu;
      }
    }
  }
  size_t pos = threadIdx.x + (blockIdx.x * blockDim.x);
  const size_t stride = devcRoundUp(dest.natom, warp_size_int);
  const size_t advance = (blockDim.x * gridDim.x);
  pos = copyCoordinateSet(dest.xcrd, 0, 0, orig.xcrd, orig.xcrd_ovrf, orig_atom_offset,
                          orig.gpos_scale, orig.gpos_bits, dest.natom, pos,  0, stride, advance);
  pos = copyCoordinateSet(dest.ycrd, 0, 0, orig.ycrd, orig.ycrd_ovrf, orig_atom_offset,
                          orig.gpos_scale, orig.gpos_bits, dest.natom, pos,  1, stride, advance);
  pos = copyCoordinateSet(dest.zcrd, 0, 0, orig.zcrd, orig.zcrd_ovrf, orig_atom_offset,
                          orig.gpos_scale, orig.gpos_bits, dest.natom, pos,  2, stride, advance);
  pos = copyCoordinateSet(dest.xvel, 0, 0, orig.xvel, orig.xvel_ovrf, orig_atom_offset,
                          orig.vel_scale, orig.gpos_bits, dest.natom, pos,   3, stride, advance);
  pos = copyCoordinateSet(dest.yvel, 0, 0, orig.yvel, orig.yvel_ovrf, orig_atom_offset,
                          orig.vel_scale, orig.gpos_bits, dest.natom, pos,   4, stride, advance);
  pos = copyCoordinateSet(dest.zvel, 0, 0, orig.zvel, orig.zvel_ovrf, orig_atom_offset,
                          orig.vel_scale, orig.gpos_bits, dest.natom, pos,   5, stride, advance);
  pos = copyCoordinateSet(dest.xfrc, 0, 0, orig.xfrc, orig.xfrc_ovrf, orig_atom_offset,
                          orig.frc_scale, orig.gpos_bits, dest.natom, pos,   6, stride, advance);
  pos = copyCoordinateSet(dest.yfrc, 0, 0, orig.yfrc, orig.yfrc_ovrf, orig_atom_offset,
                          orig.frc_scale, orig.gpos_bits, dest.natom, pos,   7, stride, advance);
  pos = copyCoordinateSet(dest.zfrc, 0, 0, orig.zfrc, orig.zfrc_ovrf, orig_atom_offset,
                          orig.frc_scale, orig.gpos_bits, dest.natom, pos,   8, stride, advance);
  pos = copyCoordinateSet(dest.xalt, 0, 0, orig.xalt, orig.xalt_ovrf, orig_atom_offset,
                          orig.gpos_scale, orig.gpos_bits, dest.natom, pos,  9, stride, advance);
  pos = copyCoordinateSet(dest.yalt, 0, 0, orig.yalt, orig.yalt_ovrf, orig_atom_offset,
                          orig.gpos_scale, orig.gpos_bits, dest.natom, pos, 10, stride, advance);
  pos = copyCoordinateSet(dest.zalt, 0, 0, orig.zalt, orig.zalt_ovrf, orig_atom_offset,
                          orig.gpos_scale, orig.gpos_bits, dest.natom, pos, 11, stride, advance);
  pos = copyCoordinateSet(dest.vxalt, 0, 0, orig.vxalt, orig.vxalt_ovrf, orig_atom_offset,
                          orig.vel_scale, orig.gpos_bits, dest.natom, pos,  12, stride, advance);
  pos = copyCoordinateSet(dest.vyalt, 0, 0, orig.vyalt, orig.vyalt_ovrf, orig_atom_offset,
                          orig.vel_scale, orig.gpos_bits, dest.natom, pos,  13, stride, advance);
  pos = copyCoordinateSet(dest.vzalt, 0, 0, orig.vzalt, orig.vzalt_ovrf, orig_atom_offset,
                          orig.vel_scale, orig.gpos_bits, dest.natom, pos,  14, stride, advance);
  pos = copyCoordinateSet(dest.fxalt, 0, 0, orig.fxalt, orig.fxalt_ovrf, orig_atom_offset,
                          orig.frc_scale, orig.gpos_bits, dest.natom, pos,  15, stride, advance);
  pos = copyCoordinateSet(dest.fyalt, 0, 0, orig.fyalt, orig.fyalt_ovrf, orig_atom_offset,
                          orig.frc_scale, orig.gpos_bits, dest.natom, pos,  16, stride, advance);
  pos = copyCoordinateSet(dest.fzalt, 0, 0, orig.fzalt, orig.fzalt_ovrf, orig_atom_offset,
                          orig.frc_scale, orig.gpos_bits, dest.natom, pos,  17, stride, advance);
}

//-------------------------------------------------------------------------------------------------
void launchCopyCoordinates(PhaseSpaceWriter *destination, const PsSynthesisReader &origin,
                           const size_t orig_atom_offset, const int index_orig,
                           const GpuDetails &gpu) {
  const int nblock = 2 * gpu.getSMPCount();

  // Compute the offsets here and pass them to the kernel.  Host memory will have lower latency
  // device memory, and computations of the offsets would simply have to be computed on every
  // thread, so shift the work to the CPU.  The API call for launchCopyCoordinates() remains as
  // simple as possible and compute time is optimized.
  const size_t orig_xfrm_offset = roundUp(9, warp_size_int) * index_orig;
  const size_t orig_bdim_offset = roundUp(6, warp_size_int) * index_orig;
  kCopyCoordinates<<<nblock, medium_block_size>>>(*destination, origin, orig_atom_offset,
                                                  orig_xfrm_offset, orig_bdim_offset);
}

//-------------------------------------------------------------------------------------------------
__global__ void __launch_bounds__(medium_block_size, 2)
kCopyCoordinates(PsSynthesisWriter dest, const size_t dest_atom_offset,
                 const size_t dest_xfrm_offset, const size_t dest_bdim_offset,
                 const PhaseSpaceReader orig) {

  // Copy the box information
  if (blockIdx.x == 0) {
    const int warp_idx = (threadIdx.x >> warp_bits);
    const int lane_idx = (threadIdx.x & warp_bits_mask_int);
    size_t i = lane_idx;
    if (warp_idx == 0) {
      while (i < 9) {
        dest.umat[dest_xfrm_offset + i] = orig.umat[i];
        i += warp_size_zu;
      }
    }
    else if (warp_idx == 1) {
      while (i < 9) {
        dest.invu[dest_xfrm_offset + i] = orig.invu[i];
        i += warp_size_zu;
      }
    }
    else if (warp_idx == 2) {
      while (i < 9) {
        dest.umat_alt[dest_xfrm_offset + i] = orig.umat_alt[i];
        i += warp_size_zu;
      }
    }
    else if (warp_idx == 3) {
      while (i < 9) {
        dest.invu_alt[dest_xfrm_offset + i] = orig.invu_alt[i];
        i += warp_size_zu;
      }
    }
    else if (warp_idx == 4) {
      while (i < 6) {
        dest.boxdims[dest_bdim_offset + i] = orig.boxdim[i];
        i += warp_size_zu;
      }
    }
    else if (warp_idx == 5) {
      while (i < 6) {
        dest.alt_boxdims[dest_bdim_offset + i] = orig.boxdim[i];
        i += warp_size_zu;
      }
    }
  }
  size_t pos = threadIdx.x + (blockIdx.x * blockDim.x);
  const size_t stride = devcRoundUp(orig.natom, warp_size_int);
  const size_t advance = (blockDim.x * gridDim.x);
  pos = copyCoordinateSet(dest.xcrd, dest.xcrd_ovrf, dest_atom_offset, dest.gpos_scale,
                          dest.gpos_bits, orig.xcrd, 0, 0, orig.natom, pos,  0, stride, advance);
  pos = copyCoordinateSet(dest.ycrd, dest.ycrd_ovrf, dest_atom_offset, dest.gpos_scale,
                          dest.gpos_bits, orig.ycrd, 0, 0, orig.natom, pos,  1, stride, advance);
  pos = copyCoordinateSet(dest.zcrd, dest.zcrd_ovrf, dest_atom_offset, dest.gpos_scale,
                          dest.gpos_bits, orig.zcrd, 0, 0, orig.natom, pos,  2, stride, advance);
  pos = copyCoordinateSet(dest.xvel, dest.xvel_ovrf, dest_atom_offset, dest.vel_scale,
                          dest.vel_bits, orig.xvel, 0, 0, orig.natom, pos,   3, stride, advance);
  pos = copyCoordinateSet(dest.yvel, dest.yvel_ovrf, dest_atom_offset, dest.vel_scale,
                          dest.vel_bits, orig.yvel, 0, 0, orig.natom, pos,   4, stride, advance);
  pos = copyCoordinateSet(dest.zvel, dest.zvel_ovrf, dest_atom_offset, dest.vel_scale,
                          dest.vel_bits, orig.zvel, 0, 0, orig.natom, pos,   5, stride, advance);
  pos = copyCoordinateSet(dest.xfrc, dest.xfrc_ovrf, dest_atom_offset, dest.frc_scale,
                          dest.frc_bits, orig.xfrc, 0, 0, orig.natom, pos,   6, stride, advance);
  pos = copyCoordinateSet(dest.yfrc, dest.yfrc_ovrf, dest_atom_offset, dest.frc_scale,
                          dest.frc_bits, orig.yfrc, 0, 0, orig.natom, pos,   7, stride, advance);
  pos = copyCoordinateSet(dest.zfrc, dest.zfrc_ovrf, dest_atom_offset, dest.frc_scale,
                          dest.frc_bits, orig.zfrc, 0, 0, orig.natom, pos,   8, stride, advance);
  pos = copyCoordinateSet(dest.xalt, dest.xalt_ovrf, dest_atom_offset, dest.gpos_scale,
                          dest.gpos_bits, orig.xalt, 0, 0, orig.natom, pos,  9, stride, advance);
  pos = copyCoordinateSet(dest.yalt, dest.yalt_ovrf, dest_atom_offset, dest.gpos_scale,
                          dest.gpos_bits, orig.yalt, 0, 0, orig.natom, pos, 10, stride, advance);
  pos = copyCoordinateSet(dest.zalt, dest.zalt_ovrf, dest_atom_offset, dest.gpos_scale,
                          dest.gpos_bits, orig.zalt, 0, 0, orig.natom, pos, 11, stride, advance);
  pos = copyCoordinateSet(dest.vxalt, dest.vxalt_ovrf, dest_atom_offset, dest.vel_scale,
                          dest.vel_bits, orig.vxalt, 0, 0, orig.natom, pos,  12, stride, advance);
  pos = copyCoordinateSet(dest.vyalt, dest.vyalt_ovrf, dest_atom_offset, dest.vel_scale,
                          dest.vel_bits, orig.vyalt, 0, 0, orig.natom, pos,  13, stride, advance);
  pos = copyCoordinateSet(dest.vzalt, dest.vzalt_ovrf, dest_atom_offset, dest.vel_scale,
                          dest.vel_bits, orig.vzalt, 0, 0, orig.natom, pos,  14, stride, advance);
  pos = copyCoordinateSet(dest.fxalt, dest.fxalt_ovrf, dest_atom_offset, dest.frc_scale,
                          dest.frc_bits, orig.fxalt, 0, 0, orig.natom, pos,  15, stride, advance);
  pos = copyCoordinateSet(dest.fyalt, dest.fyalt_ovrf, dest_atom_offset, dest.frc_scale,
                          dest.frc_bits, orig.fyalt, 0, 0, orig.natom, pos,  16, stride, advance);
  pos = copyCoordinateSet(dest.fzalt, dest.fzalt_ovrf, dest_atom_offset, dest.frc_scale,
                          dest.frc_bits, orig.fzalt, 0, 0, orig.natom, pos,  17, stride, advance);
}

//-------------------------------------------------------------------------------------------------
void launchCopyCoordinates(PsSynthesisWriter *destination, const size_t dest_atom_offset,
                           const int index_dest, const PhaseSpaceReader &origin,
                           const GpuDetails &gpu) {
  const int nblock = 2 * gpu.getSMPCount();

  // Compute the offsets here and pass them to the kernel.  Host memory will have lower latency
  // device memory, and computations of the offsets would simply have to be computed on every
  // thread, so shift the work to the CPU.  The API call for launchCopyCoordinates() remains as
  // simple as possible and compute time is optimized.
  const size_t dest_xfrm_offset = roundUp(9, warp_size_int) * index_dest;
  const size_t dest_bdim_offset = roundUp(6, warp_size_int) * index_dest;
  kCopyCoordinates<<<nblock, medium_block_size>>>(*destination, dest_atom_offset, dest_xfrm_offset,
                                                  dest_bdim_offset, origin);
}

//-------------------------------------------------------------------------------------------------
// Copy one coordinate set into another when both coordinate sets are based on 95-bit
// fixed-precision representations.  The arugments to this function follow from additional,
// templated overloads found in hpc_coordinate_copy.cuh. 
//-------------------------------------------------------------------------------------------------
__device__ __forceinline__
size_t copyCoordinateSet(llint* dest_crd, int* dest_crd_ovrf, const int dest_start_idx,
                         const int dest_bits, const llint* orig_crd, const int* orig_crd_ovrf,
                         const int orig_start_idx, const int orig_bits, const int count,
                         const size_t pos, const size_t iter, const size_t stride,
                         const size_t advance) {
  size_t ipos = pos;
  const size_t pos_limit = (iter + 1) * stride;
  while (ipos < pos_limit) {
    const int rel_pos = ipos - (iter * stride);
    if (rel_pos < count) {
      const size_t orp_idx = orig_start_idx + rel_pos;
      const int95_t orig_val = { orig_crd[orp_idx], orig_crd_ovrf[orp_idx] };
      const int95_t dest_val = changeFPBits(orig_val, orig_bits, dest_bits);
      const size_t drp_idx = dest_start_idx + rel_pos;
      dest_crd[drp_idx]      = dest_val.x;
      dest_crd_ovrf[drp_idx] = dest_val.y;
    }
    ipos += advance;
  }
  return ipos;
}
  
//-------------------------------------------------------------------------------------------------
__global__ void __launch_bounds__(medium_block_size, 2)
kCopyCoordinates(PsSynthesisWriter dest, const size_t dest_atom_offset,
                 const size_t dest_xfrm_offset, const size_t dest_bdim_offset,
                 const PsSynthesisReader orig, const size_t orig_atom_offset,
                 const size_t orig_xfrm_offset, const size_t orig_bdim_offset, const int natom) {

  // Copy the box information
  if (blockIdx.x == 0) {
    const int warp_idx = (threadIdx.x >> warp_bits);
    const int lane_idx = (threadIdx.x & warp_bits_mask_int);
    size_t i = lane_idx;
    if (warp_idx == 0) {
      while (i < 9) {
        dest.umat[dest_xfrm_offset + i] = orig.umat[orig_xfrm_offset + i];
        i += warp_size_zu;
      }
    }
    else if (warp_idx == 1) {
      while (i < 9) {
        dest.invu[dest_xfrm_offset + i] = orig.invu[orig_xfrm_offset + i];
        i += warp_size_zu;
      }
    }
    else if (warp_idx == 2) {
      while (i < 9) {
        dest.umat_alt[dest_xfrm_offset + i] = orig.umat_alt[orig_xfrm_offset + i];
        i += warp_size_zu;
      }
    }
    else if (warp_idx == 3) {
      while (i < 9) {
        dest.invu_alt[dest_xfrm_offset + i] = orig.invu_alt[orig_xfrm_offset + i];
        i += warp_size_zu;
      }
    }
    else if (warp_idx == 4) {
      while (i < 6) {
        dest.boxdims[dest_bdim_offset + i] = orig.boxdims[orig_bdim_offset + i];
        i += warp_size_zu;
      }
    }
    else if (warp_idx == 5) {
      while (i < 6) {
        dest.alt_boxdims[dest_bdim_offset + i] = orig.alt_boxdims[orig_bdim_offset + i];
        i += warp_size_zu;
      }
    }
  }
  size_t pos = threadIdx.x + (blockIdx.x * blockDim.x);
  const size_t stride = devcRoundUp(natom, warp_size_int);
  const size_t advance = (blockDim.x * gridDim.x);
  pos = copyCoordinateSet(dest.xcrd, dest.xcrd_ovrf, dest_atom_offset, dest.gpos_bits, orig.xcrd,
                          orig.xcrd_ovrf, orig_atom_offset, orig.gpos_bits, natom, pos,  0,
                          stride, advance);
  pos = copyCoordinateSet(dest.ycrd, dest.ycrd_ovrf, dest_atom_offset, dest.gpos_bits, orig.ycrd,
                          orig.ycrd_ovrf, orig_atom_offset, orig.gpos_bits, natom, pos,  1,
                          stride, advance);
  pos = copyCoordinateSet(dest.zcrd, dest.zcrd_ovrf, dest_atom_offset, dest.gpos_bits, orig.zcrd,
                          orig.zcrd_ovrf, orig_atom_offset, orig.gpos_bits, natom, pos,  2,
                          stride, advance);
  pos = copyCoordinateSet(dest.xvel, dest.xvel_ovrf, dest_atom_offset, dest.vel_bits, orig.xvel,
                          orig.xvel_ovrf, orig_atom_offset, orig.vel_bits, natom, pos,   3,
                          stride, advance);
  pos = copyCoordinateSet(dest.yvel, dest.yvel_ovrf, dest_atom_offset, dest.vel_bits, orig.yvel,
                          orig.yvel_ovrf, orig_atom_offset, orig.vel_bits, natom, pos,   4,
                          stride, advance);
  pos = copyCoordinateSet(dest.zvel, dest.zvel_ovrf, dest_atom_offset, dest.vel_bits, orig.zvel,
                          orig.zvel_ovrf, orig_atom_offset, orig.vel_bits, natom, pos,   5,
                          stride, advance);
  pos = copyCoordinateSet(dest.xfrc, dest.xfrc_ovrf, dest_atom_offset, dest.frc_bits, orig.xfrc,
                          orig.xfrc_ovrf, orig_atom_offset, orig.frc_bits, natom, pos,   6,
                          stride, advance);
  pos = copyCoordinateSet(dest.yfrc, dest.yfrc_ovrf, dest_atom_offset, dest.frc_bits, orig.yfrc,
                          orig.yfrc_ovrf, orig_atom_offset, orig.frc_bits, natom, pos,   7,
                          stride, advance);
  pos = copyCoordinateSet(dest.zfrc, dest.zfrc_ovrf, dest_atom_offset, dest.frc_bits, orig.zfrc,
                          orig.zfrc_ovrf, orig_atom_offset, orig.frc_bits, natom, pos,   8,
                          stride, advance);
  pos = copyCoordinateSet(dest.xalt, dest.xalt_ovrf, dest_atom_offset, dest.gpos_bits, orig.xalt,
                          orig.xalt_ovrf, orig_atom_offset, orig.gpos_bits, natom, pos,  9,
                          stride, advance);
  pos = copyCoordinateSet(dest.yalt, dest.yalt_ovrf, dest_atom_offset, dest.gpos_bits, orig.yalt,
                          orig.yalt_ovrf, orig_atom_offset, orig.gpos_bits, natom, pos, 10,
                          stride, advance);
  pos = copyCoordinateSet(dest.zalt, dest.zalt_ovrf, dest_atom_offset, dest.gpos_bits, orig.zalt,
                          orig.zalt_ovrf, orig_atom_offset, orig.gpos_bits, natom, pos, 11,
                          stride, advance);
  pos = copyCoordinateSet(dest.vxalt, dest.vxalt_ovrf, dest_atom_offset, dest.vel_bits, orig.vxalt,
                          orig.vxalt_ovrf, orig_atom_offset, orig.vel_bits, natom, pos,  12,
                          stride, advance);
  pos = copyCoordinateSet(dest.vyalt, dest.vyalt_ovrf, dest_atom_offset, dest.vel_bits, orig.vyalt,
                          orig.vyalt_ovrf, orig_atom_offset, orig.vel_bits, natom, pos,  13,
                          stride, advance);
  pos = copyCoordinateSet(dest.vzalt, dest.vzalt_ovrf, dest_atom_offset, dest.vel_bits, orig.vzalt,
                          orig.vzalt_ovrf, orig_atom_offset, orig.vel_bits, natom, pos,  14,
                          stride, advance);
  pos = copyCoordinateSet(dest.fxalt, dest.fxalt_ovrf, dest_atom_offset, dest.frc_bits, orig.fxalt,
                          orig.fxalt_ovrf, orig_atom_offset, orig.frc_bits, natom, pos,  15,
                          stride, advance);
  pos = copyCoordinateSet(dest.fyalt, dest.fyalt_ovrf, dest_atom_offset, dest.frc_bits, orig.fyalt,
                          orig.fyalt_ovrf, orig_atom_offset, orig.frc_bits, natom, pos,  16,
                          stride, advance);
  pos = copyCoordinateSet(dest.fzalt, dest.fzalt_ovrf, dest_atom_offset, dest.frc_bits, orig.fzalt,
                          orig.fzalt_ovrf, orig_atom_offset, orig.frc_bits, natom, pos,  17,
                          stride, advance);
}

//-------------------------------------------------------------------------------------------------
void launchCopyCoordinates(PsSynthesisWriter *destination, const size_t dest_atom_offset,
                           const int index_dest, const PsSynthesisReader &origin,
                           const size_t orig_atom_offset, const int index_orig,
                           const GpuDetails &gpu) {
  const int nblock = 2 * gpu.getSMPCount();

  // Compute the offsets here and pass them to the kernel.  Host memory will have lower latency
  // device memory, and computations of the offsets would simply have to be computed on every
  // thread, so shift the work to the CPU.  The API call for launchCopyCoordinates() remains as
  // simple as possible and compute time is optimized.
  const int xfrm_w = roundUp(9, warp_size_int);
  const int bdim_w = roundUp(6, warp_size_int);
  const size_t dest_xfrm_offset = xfrm_w * index_dest;
  const size_t dest_bdim_offset = bdim_w * index_dest;
  const size_t orig_xfrm_offset = xfrm_w * index_orig;
  const size_t orig_bdim_offset = bdim_w * index_orig;
  const int natom = origin.atom_counts[index_orig];
  kCopyCoordinates<<<nblock, medium_block_size>>>(*destination, dest_atom_offset, dest_xfrm_offset,
                                                  dest_bdim_offset, origin, orig_atom_offset,
                                                  orig_xfrm_offset, orig_bdim_offset, natom);
}


} // namespace trajectory
} // namespace stormm
