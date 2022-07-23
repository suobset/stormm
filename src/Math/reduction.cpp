#include "Constants/scaling.h"
#include "Math/rounding.h"
#include "reduction.h"

namespace omni {
namespace math {

using card::HybridKind;
using math::roundUp;

//-------------------------------------------------------------------------------------------------
ReductionBridge::ReductionBridge(const size_t n_values) :
  x_buffer{HybridKind::POINTER, "bridge_xbuff"},
  y_buffer{HybridKind::POINTER, "bridge_ybuff"},
  z_buffer{HybridKind::POINTER, "bridge_zbuff"},
  storage{3LLU * roundUp(n_values, warp_size_zu)}
{
  const size_t padded_nval = roundUp(n_values, warp_size_zu);
  x_buffer.setPointer(&storage,                  0, n_values);
  y_buffer.setPointer(&storage,        padded_nval, n_values);
  z_buffer.setPointer(&storage, 2LLU * padded_nval, n_values);
}

//-------------------------------------------------------------------------------------------------
const double* ReductionBridge::getPointer(const CartesianDimension cdim,
                                          const HybridTargetLevel tier) const {
  switch (cdim) {
  case CartesianDimension::X:
    return x_buffer.data(tier);
  case CartesianDimension::Y:
    return y_buffer.data(tier);
  case CartesianDimension::Z:
    return z_buffer.data(tier);    
  }
  __builtin_unreachable();
}
  
//-------------------------------------------------------------------------------------------------
double* ReductionBridge::getPointer(const CartesianDimension cdim, const HybridTargetLevel tier) {
  switch (cdim) {
  case CartesianDimension::X:
    return x_buffer.data(tier);
  case CartesianDimension::Y:
    return y_buffer.data(tier);
  case CartesianDimension::Z:
    return z_buffer.data(tier);    
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
ReductionKit::ReductionKit(const int nrdwu_in, const RdwuPerSystem rps_in,
                           const int* rdwu_abstracts_in, const int* atom_counts_in) :
  nrdwu{nrdwu_in}, rps{rps_in}, rdwu_abstracts{rdwu_abstracts_in}, atom_counts{atom_counts_in}
{}

//-------------------------------------------------------------------------------------------------
ReductionKit::ReductionKit(const AtomGraphSynthesis &poly_ag, const HybridTargetLevel tier) :
  nrdwu{poly_ag.getReductionWorkUnitCount()},
  rps{poly_ag.getRdwuPerSystem()},
  rdwu_abstracts{poly_ag.getReductionWorkUnitAbstracts().data(tier)},
  atom_counts{poly_ag.getSystemAtomCounts().data(tier)}
{}

//-------------------------------------------------------------------------------------------------
double gatherNormalization(const ReductionSubstrate<llint> rsbs, const int start_pos,
                           const int end_pos) {

  // Scale down the values of vector elements if fixed precision is in use (if it is not, the
  // value of the scaling factor will just one 1.0).  While the double-precision format could
  // likely handle the squared values of large forces scaled up by, say, 2^72, that is 2^144 extra
  // stress on the number format that is not needed.  If this were ever to go to single-precision,
  // it would be impossible.
  const bool extended_precision = (rsbs.x_read_ovrf != nullptr);
  double tsum = 0.0;
  if (extended_precision) {
    if (rsbs.y_read == nullptr && rsbs.z_read == nullptr) {
      for (int j = start_pos; j < end_pos; j++) {
        const double dx = ((static_cast<double>(rsbs.x_read_ovrf[j]) *
                            max_llint_accumulation) + static_cast<double>(rsbs.x_read[j])) *
                          rsbs.inv_fp_scaling;
        tsum += (dx * dx);
      }
    }
    else if (rsbs.y_read != nullptr && rsbs.z_read == nullptr) {
      for (int j = start_pos; j < end_pos; j++) {
        const double dx = ((static_cast<double>(rsbs.x_read_ovrf[j]) *
                            max_llint_accumulation) + static_cast<double>(rsbs.x_read[j])) *
                          rsbs.inv_fp_scaling;
        const double dy = ((static_cast<double>(rsbs.y_read_ovrf[j]) *
                            max_llint_accumulation) + static_cast<double>(rsbs.y_read[j])) *
                          rsbs.inv_fp_scaling;
        tsum += (dx * dx) + (dy * dy);
      }
    }
    else {
      for (int j = start_pos; j < end_pos; j++) {
        const double dx = ((static_cast<double>(rsbs.x_read_ovrf[j]) *
                            max_llint_accumulation) + static_cast<double>(rsbs.x_read[j])) *
                          rsbs.inv_fp_scaling;
        const double dy = ((static_cast<double>(rsbs.y_read_ovrf[j]) *
                            max_llint_accumulation) + static_cast<double>(rsbs.y_read[j])) *
                          rsbs.inv_fp_scaling;
        const double dz = ((static_cast<double>(rsbs.z_read_ovrf[j]) *
                            max_llint_accumulation) + static_cast<double>(rsbs.z_read[j])) *
                          rsbs.inv_fp_scaling;
        tsum += (dx * dx) + (dy * dy) + (dz * dz);
      }
    }
  }
  else {
    double tsum = 0.0;
    if (rsbs.y_read == nullptr && rsbs.z_read == nullptr) {
      for (int j = start_pos; j < end_pos; j++) {
        const double dx = static_cast<double>(rsbs.x_read[j]) * rsbs.inv_fp_scaling;
        tsum += (dx * dx);
      }
    }
    else if (rsbs.y_read != nullptr && rsbs.z_read == nullptr) {
      for (int j = start_pos; j < end_pos; j++) {
        const double dx = static_cast<double>(rsbs.x_read[j]) * rsbs.inv_fp_scaling;
        const double dy = static_cast<double>(rsbs.x_read[j]) * rsbs.inv_fp_scaling;
        tsum += (dx * dx) + (dy * dy);
      }
    }
    else {
      for (int j = start_pos; j < end_pos; j++) {
        const double dx = static_cast<double>(rsbs.x_read[j]) * rsbs.inv_fp_scaling;
        const double dy = static_cast<double>(rsbs.x_read[j]) * rsbs.inv_fp_scaling;
        const double dz = static_cast<double>(rsbs.x_read[j]) * rsbs.inv_fp_scaling;
        tsum += (dx * dx) + (dy * dy) + (dz * dz);
      }
    }
  }
  return tsum;
}

//-------------------------------------------------------------------------------------------------
double3 gatherCenterOnZero(const ReductionSubstrate<llint> rsbs, const int start_pos,
                           const int end_pos) {
  const bool extended_precision = (rsbs.x_read_ovrf != nullptr);
  double tsum_x = 0.0;
  double tsum_y = 0.0;
  double tsum_z = 0.0;
  const int nval = end_pos - start_pos;
  if (extended_precision) {
    tsum_x = sum<double>(&rsbs.x_read_ovrf[start_pos], nval);
    tsum_x *= max_llint_accumulation;
    tsum_x += sum<double>(&rsbs.x_read[start_pos], nval);
    if (rsbs.y_read != nullptr) {
      tsum_y = sum<double>(&rsbs.y_read_ovrf[start_pos], nval);
      tsum_y *= max_llint_accumulation;
      tsum_y += sum<double>(&rsbs.y_read[start_pos], nval);
    }
    if (rsbs.z_read != nullptr) {
      tsum_z = sum<double>(&rsbs.z_read_ovrf[start_pos], nval);
      tsum_z *= max_llint_accumulation;
      tsum_z += sum<double>(&rsbs.z_read[start_pos], nval);
    }
  }
  else {
    tsum_x = sum<double>(&rsbs.x_read[start_pos], nval);
    if (rsbs.y_read != nullptr) {
      tsum_y = sum<double>(&rsbs.y_read[start_pos], nval);
    }
    if (rsbs.z_read != nullptr) {
      tsum_z = sum<double>(&rsbs.z_read[start_pos], nval);
    }
  }
  return { -tsum_x, -tsum_y, -tsum_z };
}

//-------------------------------------------------------------------------------------------------
void scatterNormalization(ReductionSubstrate<llint> rsbs, const double tsum, const int start_pos,
                          const int end_pos) {

  // When normalizing extended fixed-precision values, it is not necessary to scale them down with
  // whatever fixed-precision scaling factor only to scale them right back up.  Just divide by the
  // magnitude of the vector found earlier.
  const bool extended_precision = (rsbs.x_read_ovrf != nullptr);
  const double nfactor = 1.0 / sqrt(tsum);
  if (extended_precision) {
    if (rsbs.y_read == nullptr && rsbs.z_read == nullptr) {
      for (int j = start_pos; j < end_pos; j++) {
        const double dx = ((static_cast<double>(rsbs.x_write_ovrf[j]) *
                            max_llint_accumulation) + static_cast<double>(rsbs.x_write[j]));
        splitRealAccumulation(dx * nfactor, &rsbs.x_write[j], &rsbs.x_write_ovrf[j]);
      }
    }
    else if (rsbs.y_read != nullptr && rsbs.z_read == nullptr) {
      for (int j = start_pos; j < end_pos; j++) {
        const double dx = ((static_cast<double>(rsbs.x_write_ovrf[j]) *
                            max_llint_accumulation) + static_cast<double>(rsbs.x_write[j]));
        const double dy = ((static_cast<double>(rsbs.y_write_ovrf[j]) *
                            max_llint_accumulation) + static_cast<double>(rsbs.y_write[j]));
        splitRealAccumulation(dx * nfactor, &rsbs.x_write[j], &rsbs.x_write_ovrf[j]);
        splitRealAccumulation(dy * nfactor, &rsbs.y_write[j], &rsbs.y_write_ovrf[j]);
      }
    }
    else {
      for (int j = start_pos; j < end_pos; j++) {
        const double dx = ((static_cast<double>(rsbs.x_write_ovrf[j]) *
                            max_llint_accumulation) + static_cast<double>(rsbs.x_write[j]));
        const double dy = ((static_cast<double>(rsbs.y_write_ovrf[j]) *
                            max_llint_accumulation) + static_cast<double>(rsbs.y_write[j]));
        const double dz = ((static_cast<double>(rsbs.z_write_ovrf[j]) *
                            max_llint_accumulation) + static_cast<double>(rsbs.z_write[j]));
        splitRealAccumulation(dx * nfactor, &rsbs.x_write[j], &rsbs.x_write_ovrf[j]);
        splitRealAccumulation(dy * nfactor, &rsbs.y_write[j], &rsbs.y_write_ovrf[j]);
        splitRealAccumulation(dz * nfactor, &rsbs.z_write[j], &rsbs.z_write_ovrf[j]);
      }
    }
  }
  else {
    if (rsbs.y_read == nullptr && rsbs.z_read == nullptr) {
      for (int j = start_pos; j < end_pos; j++) {
        const double dx = static_cast<double>(rsbs.x_write[j]);
        rsbs.x_write[j] = llround(dx * nfactor);
      }
    }
    else if (rsbs.y_read != nullptr && rsbs.z_read == nullptr) {
      for (int j = start_pos; j < end_pos; j++) {
        const double dx = static_cast<double>(rsbs.x_write[j]);
        const double dy = static_cast<double>(rsbs.y_write[j]);
        rsbs.x_write[j] = llround(dx * nfactor);
        rsbs.y_write[j] = llround(dy * nfactor);
      }
    }
    else {
      for (int j = start_pos; j < end_pos; j++) {
        const double dx = static_cast<double>(rsbs.x_write[j]);
        const double dy = static_cast<double>(rsbs.y_write[j]);
        const double dz = static_cast<double>(rsbs.z_write[j]);
        rsbs.x_write[j] = llround(dx * nfactor);
        rsbs.y_write[j] = llround(dy * nfactor);
        rsbs.z_write[j] = llround(dz * nfactor);
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------
void scatterCenterOnZero(ReductionSubstrate<llint> rsbs, const double tsum_x, const double tsum_y,
                         const double tsum_z, const int natom, const int start_pos,
                         const int end_pos) {

  // Values entering the calculation of the sum were never scaled down from their fixed-precision
  // values, so there is no need to rescale any fixed-precision representations here.
  const bool extended_precision = (rsbs.x_read_ovrf != nullptr);
  const double inv_norm = 1.0 / static_cast<double>(natom);
  if (extended_precision) {
    const double center_x = tsum_x * inv_norm;
    for (int j = start_pos; j < end_pos; j++) {
      splitRealAccumulation(center_x, &rsbs.x_write[j], &rsbs.x_write_ovrf[j]);
    }
    if (rsbs.y_read != nullptr) {
      const double center_y = tsum_y * inv_norm;
      for (int j = start_pos; j < end_pos; j++) {
        splitRealAccumulation(center_y, &rsbs.y_write[j], &rsbs.y_write_ovrf[j]);
      }
    }
    if (rsbs.z_read != nullptr) {
      const double center_z = tsum_z * inv_norm;
      for (int j = start_pos; j < end_pos; j++) {
        splitRealAccumulation(center_z, &rsbs.z_write[j], &rsbs.z_write_ovrf[j]);
      }
    }
  }
  else {
    const llint center_x = llround(tsum_x * inv_norm);
    addScalarToVector(&rsbs.x_write[start_pos], end_pos - start_pos, center_x);
    if (rsbs.y_read != nullptr) {
      const llint center_y = llround(tsum_y * inv_norm);
      addScalarToVector(&rsbs.y_write[start_pos], end_pos - start_pos, center_y);
    }
    if (rsbs.z_read != nullptr) {
      const llint center_z = llround(tsum_z * inv_norm);
      addScalarToVector(&rsbs.z_write[start_pos], end_pos - start_pos, center_z);
    }
  }
}

} // namespace math
} // namespace omni
