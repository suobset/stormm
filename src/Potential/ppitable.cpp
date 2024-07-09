#include "copyright.h"
#include "ppitable.h"

namespace stormm {
namespace energy {

//-------------------------------------------------------------------------------------------------
PPITable::PPITable(const NonbondedTheme theme_in, const BasisFunctions basis_set_in,
                   const TableIndexing indexing_method_in, const double cutoff_in,
                   const double argument_offset_in, const double dsum_tol_in,
                   const int mantissa_bits_in, const double coulomb_in,
                   const double min_range_in) :
    theme{theme_in}, basis_set{basis_set_in}, indexing_method{indexing_method_in},
    cutoff{cutoff_in}, max_range{cutoff_in}, min_range{min_range_in},
    argument_offset{argument_offset_in}, dsum_tol{dsum_tol_in},
    ew_coeff{ewaldCoefficient(cutoff_in, dsum_tol_in)}, mantissa_bits{mantissa_bits_in},
    coulomb{coulomb_in}, exclusion_offset{0},
    dp_detail_bitmask{doublePrecisionSplineDetailMask(mantissa_bits)},
    sp_detail_bitmask{singlePrecisionSplineDetailMask(mantissa_bits)},
    energy{HybridKind::POINTER, "etab_energy"},
    force{HybridKind::POINTER, "etab_force"},
    energy_with_exclusions{HybridKind::POINTER, "etab_energy_excl"},
    force_with_exclusions{HybridKind::POINTER, "etab_force_excl"},
    sp_energy{HybridKind::POINTER, "etab_energy"},
    sp_force{HybridKind::POINTER, "etab_force"},
    sp_energy_with_exclusions{HybridKind::POINTER, "etab_energy_excl"},
    sp_force_with_exclusions{HybridKind::POINTER, "etab_force_excl"},
    coeffs{HybridKind::ARRAY, "etab_coeffs"},
    sp_coeffs{HybridKind::ARRAY, "etab_sp_coeffs"}
{
  // Create the four spline tables, then feed into the same templated functions for re-arranging
  // them as the constructors taking spline table inputs.
  max_range = exp2(ceil(log2(cutoff)));
  switch (indexing_method) {
  case TableIndexing::ARG:
  case TableIndexing::ARG_OFFSET:
    break;
  case TableIndexing::SQUARED_ARG:
  case TableIndexing::SQ_ARG_OFFSET:
    argument_offset *= argument_offset;
    min_range *= min_range;
    break;
  }
  const std::vector<LogScaleSpline<double4>> splv = buildAllSplineTables<double4>();
  const std::vector<LogScaleSpline<float4>> sp_splv = buildAllSplineTables<float4>();
  populateCoefficients<double4>(splv[0], splv[1], splv[2], splv[3]);
  populateCoefficients<float4>(sp_splv[0], sp_splv[1], sp_splv[2], sp_splv[3]);
}

//-------------------------------------------------------------------------------------------------
NonbondedTheme PPITable::getTheme() const {
  return theme;
}

//-------------------------------------------------------------------------------------------------
double PPITable::getCutoff() const {
  return cutoff;
}

//-------------------------------------------------------------------------------------------------
double PPITable::getMaximumRange() const {
  return max_range;
}

//-------------------------------------------------------------------------------------------------
double PPITable::getIndexingOffset() const {
  return argument_offset;
}

//-------------------------------------------------------------------------------------------------
double PPITable::getDirectSumTolerance() const {
  return dsum_tol;
}

//-------------------------------------------------------------------------------------------------
int PPITable::getBitStride() const {
  return mantissa_bits;
}

//-------------------------------------------------------------------------------------------------
double PPITable::getEwaldCoefficient() const {
  return ew_coeff;
}

//-------------------------------------------------------------------------------------------------
double PPITable::getGaussianWidth() const {
  return dsum_tol;
}

//-------------------------------------------------------------------------------------------------
const PPIKit<double, double4> PPITable::dpData(const HybridTargetLevel tier) const {
  const int index_bit_count = mantissa_bits + 12;
  const int index_shift = 64 - index_bit_count;
  return PPIKit<double, double4>(theme, basis_set, indexing_method, exclusion_offset / 2,
                                 exclusion_offset, index_shift, dp_detail_bitmask,
                                 sp_detail_bitmask, argument_offset, energy.data(tier),
                                 force.data(tier), energy_with_exclusions.data(tier),
                                 force_with_exclusions.data(tier));
}

//-------------------------------------------------------------------------------------------------
const PPIKit<float, float4> PPITable::spData(const HybridTargetLevel tier) const {
  const int index_bit_count = mantissa_bits + 9;
  const int index_shift = 32 - index_bit_count;
  return PPIKit<float, float4>(theme, basis_set, indexing_method, exclusion_offset / 2,
                               exclusion_offset, index_shift, dp_detail_bitmask, sp_detail_bitmask,
                               argument_offset, sp_energy.data(tier), sp_force.data(tier),
                               sp_energy_with_exclusions.data(tier),
                               sp_force_with_exclusions.data(tier));
}

//-------------------------------------------------------------------------------------------------
LogSplineForm PPITable::findMissingForm(const uint holdings) const {
  if ((holdings & 0x1) == 0) {
    return LogSplineForm::ELEC_PME_DIRECT;
  }
  else if ((holdings & 0x2) == 0) {
    return LogSplineForm::ELEC_PME_DIRECT_EXCL;
  }
  else if ((holdings & 0x4) == 0) {
    return LogSplineForm::DELEC_PME_DIRECT;
  }
  else if ((holdings & 0x8) == 0) {
    return LogSplineForm::DELEC_PME_DIRECT_EXCL;
  }
  __builtin_unreachable();
}

} // namespace energy
} // namespace stormm
