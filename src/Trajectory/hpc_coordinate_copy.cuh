// -*-c++-*-
#ifndef STORMM_HPC_COORDINATE_COPY_CUH
#define STORMM_HPC_COORDINATE_COPY_CUH

#include "DataTypes/common_types.h"
#include "DataTypes/stormm_vector_types.h"
#include "Numerics/split_fixed_precision.h"

namespace stormm {
namespace trajectory {

using numerics::max_int_accumulation;
using numerics::max_int_accumulation_f;
using numerics::max_int_accumulation_ll;
using numerics::max_llint_accumulation;
using numerics::max_llint_accumulation_f;

#include "Numerics/accumulation.cui"

/// \brief Copy a set of coordinates from one array into another, inferring the conversion
///        operation based on the origin and destination data types.
///
/// Overloaded:
///   - Convert one general type into another (here, a "general type" is a coordinate format
///     which expresses the X, Y, or Z component of a position, velocity, or force as a scalar
///     type, e.g. float, llint, double)
///   - Convert a general type into 95-bit fixed-precision (int95_t)
///   - Convert 95-bit fixed-precision into a general type
///
/// \param dest_crd         Destination array to hold coordinates
/// \param dest_crd_ovrf    Overflow bits for dest_crd (its presence implies dest_crd is llint)
/// \param dest_start_idx   Starting index of the destination array at which to begin writing
///                         (likely to be zero for single-system coordinate objects)
/// \param dest_scale       Scaling factor for coordinates in the destination array
/// \param dest_scale_bits  Number of bits after the decimal for fixed-precision coordinates in
///                         the destination array (the base-2 logarithm of dest_scale)
/// \param orig_crd         Array holding the original coordinates
/// \param orig_crd_ovrf    Overflow bits for orig_crd (its presence implies dest_crd is llint)
/// \param orig_start_idx   Starting index of the origin array at which to begin reading (likely
///                         to be zero for single-system coordinate objects)
/// \param orig_scale       Scaling factor for coordinates in the origin array
/// \param orig_scale_bits  Number of bits after the decimal for fixed-precision coordinates in
///                         the origin array (the base-2 logarithm of dest_scale)
/// \param count            The number of particle coordinates to copy
/// \param pos              Internal counter passed in from the calling function
/// \param iter             The number of passes through this array, each incrementing the counter
///                         behind pos
/// \param stride           The padded form of count (pre-computed and passed in for convenience
///                         and probably reduced register pressure)
/// \para, advance          The amount to increment an internal variant of the pos counter (this
///                         can be the width of a warp, the width of a block, or the width of the
///                         entire kernel launch grid depending on the context in which this copy
///                         routine is called)
/// \{
template <typename Tdest, typename Torig> __device__ __forceinline__
size_t copyCoordinateSet(Tdest* dest_crd, const int dest_start_idx, const double dest_scale,
                         const int dest_scale_bits, const Torig* orig_crd,
                         const int orig_start_idx, const double orig_scale,
                         const int orig_scale_bits, const int count, const size_t pos,
                         const size_t iter, const size_t stride, const size_t advance) {
  const size_t pos_limit = (iter + 1) * stride;
  size_t ipos = pos;
  while (ipos < pos_limit) {
    const int rel_pos = ipos - (iter * stride);
    if (rel_pos < count) {

      const size_t drp_idx = dest_start_idx + rel_pos;
      const size_t orp_idx = orig_start_idx + rel_pos;
      if (orig_scale_bits == 0) {
        if (dest_scale_bits == 0) {
          dest_crd[drp_idx] = orig_crd[orp_idx];
        }
        else {
          const llint ival = __double2ll_rn(orig_crd[orp_idx] * dest_scale);
          dest_crd[drp_idx] = ival;
        }
      }
      else {
        if (dest_scale_bits == 0) {
          dest_crd[drp_idx] = (double)(orig_crd[orp_idx]) / orig_scale;
        }
        else {
          if (dest_scale_bits > orig_scale_bits) {
            const llint dconv = orig_crd[orp_idx];
            dest_crd[drp_idx] = (dconv << (dest_scale_bits - orig_scale_bits));
          }
          else if (dest_scale_bits == orig_scale_bits) {
            dest_crd[drp_idx] = orig_crd[orp_idx];
          }
          else {
            const llint dconv = orig_crd[orp_idx];
            dest_crd[drp_idx] = (dconv >> (orig_scale_bits - dest_scale_bits));
          }
        }
      }
    }
    ipos += advance;
  }
  return ipos;
}

template <typename Torig> __device__ __forceinline__
size_t copyCoordinateSet(llint* dest_crd, int* dest_crd_ovrf, const int dest_start_idx,
                         const double dest_scale, const int dest_scale_bits, const Torig* orig_crd,
                         const int orig_start_idx, const int orig_scale_bits, const int count,
                         const size_t pos, const size_t iter, const size_t stride,
                         const size_t advance) {
  const size_t pos_limit = (iter + 1) * stride;
  size_t ipos = pos;
  while (ipos < pos_limit) {
    const int rel_pos = ipos - (iter * stride);
    if (rel_pos < count) {
      const size_t drp_idx = dest_start_idx + rel_pos;

      // Detect the presence of a signed integer, fixed-precision representation in the origin
      // by a nonzero number of scaling bits.  Copy that directly into the long long int portion
      // of a 95-bit representation and upgrade (or downgrade) the bit count as appropriate.
      if (orig_scale_bits > 0) {
        const int95_t oval = { (long long int)(orig_crd[orig_start_idx + rel_pos]), 0 };
        const int95_t dval = changeFPBits(oval, orig_scale_bits, dest_scale_bits);
        dest_crd[drp_idx]      = dval.x;
        dest_crd_ovrf[drp_idx] = dval.y;
      }
      else {

        // Ignore limits on various types of coordinates which could support an assumption that
        // the high 32 bits are unnecessary.  The memory bandwidth will be the limiting factor.
        // Treat every conversion as a double-to-int95_t operation.
        const int95_t dval = doubleToInt95(orig_crd[orig_start_idx + rel_pos] * dest_scale);
        dest_crd[drp_idx]      = dval.x;
        dest_crd_ovrf[drp_idx] = dval.y;
      }
    }
    ipos += advance;
  }
  return ipos;
}

template <typename Tdest> __device__ __forceinline__
size_t copyCoordinateSet(Tdest* dest_crd, const int dest_start_idx, const int dest_scale_bits,
                         const llint* orig_crd, const int* orig_crd_ovrf, const int orig_start_idx,
                         const double orig_scale, const int orig_scale_bits, const int count,
                         const size_t pos, const size_t iter, const size_t stride,
                         const size_t advance) {
  const size_t pos_limit = (iter + 1) * stride;
  size_t ipos = pos;
  while (ipos < pos_limit) {
    const int rel_pos = ipos - (iter * stride);
    if (rel_pos < count) {
      const size_t drp_idx = dest_start_idx + rel_pos;
      const size_t orp_idx = orig_start_idx + rel_pos;

      // Detect the presence of a signed integer, fixed-precision representation in the destination
      // by a nonzero number of scaling bits.  The result will emerge as, at most, the low 64 bits
      // of the rescaled 95-bit quantity.
      if (dest_scale_bits > 0) {
        const int95_t oval = { orig_crd[orp_idx], orig_crd_ovrf[orp_idx] };
        const int95_t dval = changeFPBits(oval, orig_scale_bits, dest_scale_bits);
        dest_crd[drp_idx] = dval.x;
      }
      else {
        const double dval = int95ToDouble(orig_crd[orp_idx], orig_crd_ovrf[orp_idx]) / orig_scale;
        dest_crd[drp_idx] = dval;
      }
    }
    ipos += advance;
  }
  return ipos;  
}
/// \}

} // namespace trajectory
} // namespace stormm

#endif
