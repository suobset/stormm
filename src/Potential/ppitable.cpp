#include "copyright.h"
#include "ppitable.h"

namespace stormm {
namespace energy {

//-------------------------------------------------------------------------------------------------
PPITable::PPITable(const NonbondedTheme theme_in, const PrecisionModel precision_in,
                   const BasisFunctions basis_set_in, const TableIndexing indexing_method_in,
                   const double cutoff_in, const double dsum_tol_in, const int mantissa_bits_in,
                   const double coulomb_in, const double min_spl_compute_range_in,
                   const double min_offset_in) :
    theme{theme_in}, precision{precision_in}, basis_set{basis_set_in},
    indexing_method{indexing_method_in}, cutoff{cutoff_in}, max_range{cutoff_in},
    dsum_tol{dsum_tol_in}, ew_coeff{ewaldCoefficient(cutoff_in, dsum_tol_in)},
    mantissa_bits{mantissa_bits_in}, coulomb{coulomb_in}, exclusion_offset{0},
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
  double min_spl_compute_range = min_spl_compute_range_in;
  double min_offset = min_offset_in;
  switch (indexing_method) {
  case TableIndexing::ARG:
  case TableIndexing::ARG_OFFSET:
    break;
  case TableIndexing::SQUARED_ARG:
  case TableIndexing::SQ_ARG_OFFSET:
    min_offset *= min_offset;
    min_spl_compute_range *= min_spl_compute_range;
    break;
  }
  switch (precision) {
  case PrecisionModel::DOUBLE:
    {
      std::vector<LogScaleSpline<double4>> splv;
      splv.reserve(4);
      switch (theme) {
      case NonbondedTheme::ELECTROSTATIC:
        {
          const std::vector<LogSplineForm> all_forms = { LogSplineForm::ELEC_PME_DIRECT,
                                                         LogSplineForm::ELEC_PME_DIRECT_EXCL,
                                                         LogSplineForm::DELEC_PME_DIRECT,
                                                         LogSplineForm::DELEC_PME_DIRECT_EXCL };
          for (size_t i = 0; i < all_forms.size(); i++) {
            splv.emplace_back(all_forms[i], ew_coeff, coulomb, mantissa_bits, max_range,
                              min_spl_compute_range, indexing_method, basis_set, 2, min_offset,
                              ExceptionResponse::DIE);
          }
        }
      case NonbondedTheme::VAN_DER_WAALS:
      case NonbondedTheme::ALL:
        break;
      }
    }
    break;
  case PrecisionModel::SINGLE:
    {
      std::vector<LogScaleSpline<float4>> splv;
      splv.reserve(4);
      switch (theme) {
      case NonbondedTheme::ELECTROSTATIC:
        {
          const std::vector<LogSplineForm> all_forms = { LogSplineForm::ELEC_PME_DIRECT,
                                                         LogSplineForm::ELEC_PME_DIRECT_EXCL,
                                                         LogSplineForm::DELEC_PME_DIRECT,
                                                         LogSplineForm::DELEC_PME_DIRECT_EXCL };
          for (size_t i = 0; i < all_forms.size(); i++) {
            splv.emplace_back(all_forms[i], ew_coeff, coulomb, mantissa_bits, max_range, 0.015625,
                              indexing_method, basis_set, 2, min_spl_compute_range,
                              ExceptionResponse::DIE);
          }
        }
      case NonbondedTheme::VAN_DER_WAALS:
      case NonbondedTheme::ALL:
        break;
      }
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
NonbondedTheme PPITable::getTheme() const {
  return theme;
}

//-------------------------------------------------------------------------------------------------
PrecisionModel PPITable::getPrecision() const {
  return precision;
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

} // namespace energy
} // namespace stormm
