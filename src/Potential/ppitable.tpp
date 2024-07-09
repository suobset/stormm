// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace energy {

//-------------------------------------------------------------------------------------------------
template <typename T, typename T4>
PPIKit<T, T4>::PPIKit(const NonbondedTheme theme_in, const BasisFunctions basis_in,
                      const TableIndexing lookup_in, const int index_bound_in,
                      const int excl_offset_in, const int index_shift_bits_in,
                      const ullint dp_detail_mask_in, const uint sp_detail_mask_in,
                      const T arg_offset_in, const T4* energy_in, const T4* force_in,
                      const T4* energy_excl_in, const T4* force_excl_in) :
    theme{theme_in}, basis{basis_in}, lookup{lookup_in}, index_bound{index_bound_in},
    excl_offset{excl_offset_in}, index_shift_bits{index_shift_bits_in},
    dp_detail_mask{dp_detail_mask_in}, sp_detail_mask{sp_detail_mask_in},
    arg_offset{arg_offset_in}, energy{energy_in}, force{force_in}, energy_excl{energy_excl_in},
    force_excl{force_excl_in}
{}
  
//-------------------------------------------------------------------------------------------------
template <typename T4>
PPITable::PPITable(const LogScaleSpline<T4> &spl_a, const LogScaleSpline<T4> &spl_b,
                   const LogScaleSpline<T4> &spl_c, const LogScaleSpline<T4> &spl_d,
                   const double cutoff_in) :
    theme{findTheme(spl_a)},
    basis_set{spl_a.getBasisSet()},
    indexing_method{spl_a.getIndexingMethod()},
    cutoff{cutoff_in},
    max_range{exp2(ceil(log2(cutoff)))},
    argument_offset{spl_a.getIndexingOffset()},
    dsum_tol{recoverDirectSumTolerance(spl_a.getEwaldCoefficient(), cutoff_in)},
    ew_coeff{ewaldCoefficient(cutoff_in, dsum_tol)},
    mantissa_bits{spl_a.getBitStride()},
    coulomb{spl_a.getCoulombConstant()},
    exclusion_offset{0},
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
  checkSplineCompatibility(spl_a, spl_b);
  checkSplineCompatibility(spl_a, spl_c);
  checkSplineCompatibility(spl_a, spl_d);

  // Determine which potential is which and populate the appropriate coefficients array.
  const LogScaleSpline<T4> *u_ptr = findNonExclPotential(spl_a, spl_b, spl_c, spl_d);
  const LogScaleSpline<T4> *du_ptr = findNonExclForce(spl_a, spl_b, spl_c, spl_d);
  const LogScaleSpline<T4> *u_excl_ptr = findExclPotential(spl_a, spl_b, spl_c, spl_d);
  const LogScaleSpline<T4> *du_excl_ptr = findExclForce(spl_a, spl_b, spl_c, spl_d);
  populateCoefficients<T4>(*u_ptr, *du_ptr, *u_excl_ptr, *du_excl_ptr);

  // Populate the remaining coefficients array.
  const size_t ct = std::type_index(typeid(T4)).hash_code();
  if (ct == double4_type_index) {
    const std::vector<LogScaleSpline<float4>> splv = buildAllSplineTables<float4>();
    populateCoefficients<float4>(splv[0], splv[1], splv[2], splv[3]);
  }
  else {
    const std::vector<LogScaleSpline<double4>> splv = buildAllSplineTables<double4>();
    populateCoefficients<double4>(splv[0], splv[1], splv[2], splv[3]);
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T4>
PPITable::PPITable(const LogScaleSpline<T4> &spl_a, const LogScaleSpline<T4> &spl_b,
                   const LogScaleSpline<T4> &spl_c) :
    PPITable(spl_a, spl_b, spl_c, getTablePriority(spl_a, spl_b, spl_c))
{}

//-------------------------------------------------------------------------------------------------
template <typename T4>
PPITable::PPITable(const LogScaleSpline<T4> &spl_a, const LogScaleSpline<T4> &spl_b) :
    PPITable(spl_a, spl_b, getTablePriority(spl_a, spl_b))
{}

//-------------------------------------------------------------------------------------------------
template <typename T4>
PPITable::PPITable(const LogScaleSpline<T4> &spl_a) :
    PPITable(spl_a, getTablePriority(spl_a))
{}

//-------------------------------------------------------------------------------------------------
template <typename T4> NonbondedTheme PPITable::findTheme(const LogScaleSpline<T4> &spl) const {
  switch (spl.getForm()) {
  case LogSplineForm::ELEC_PME_DIRECT:
  case LogSplineForm::ELEC_PME_DIRECT_EXCL:
  case LogSplineForm::DELEC_PME_DIRECT:
  case LogSplineForm::DELEC_PME_DIRECT_EXCL:
    return NonbondedTheme::ELECTROSTATIC;
  case LogSplineForm::CUSTOM:
    return NonbondedTheme::ALL;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename T4> std::vector<LogScaleSpline<T4>> PPITable::buildAllSplineTables() const {
  std::vector<LogScaleSpline<T4>> result;
  result.reserve(4);
  switch (theme) {
  case NonbondedTheme::ELECTROSTATIC:
    {
      const std::vector<LogSplineForm> all_forms = { LogSplineForm::ELEC_PME_DIRECT,
                                                     LogSplineForm::DELEC_PME_DIRECT,
                                                     LogSplineForm::ELEC_PME_DIRECT_EXCL,
                                                     LogSplineForm::DELEC_PME_DIRECT_EXCL };
      for (size_t i = 0; i < all_forms.size(); i++) {
        result.emplace_back(all_forms[i], ew_coeff, coulomb, mantissa_bits, max_range, min_range,
                            indexing_method, basis_set, 2, argument_offset,
                            ExceptionResponse::DIE);
      } 
    }
    break;
  case NonbondedTheme::VAN_DER_WAALS:
  case NonbondedTheme::ALL:
    break;
  }
  return result;
}
  
//-------------------------------------------------------------------------------------------------
template <typename T4>
const LogScaleSpline<T4>& PPITable::findNonExclPotential(const LogScaleSpline<T4> &spl_a,
                                                         const LogScaleSpline<T4> &spl_b,
                                                         const LogScaleSpline<T4> &spl_c,
                                                         const LogScaleSpline<T4> &spl_d) const {
  const std::vector<const LogScaleSpline<T4>*> lss = { spl_a.getSelfPointer(),
                                                       spl_b.getSelfPointer(),
                                                       spl_c.getSelfPointer(),
                                                       spl_d.getSelfPointer() };
  for (int i = 0; i < 4; i++) {
    switch (lss[i]->getForm()) {
    case LogSplineForm::ELEC_PME_DIRECT:
      return *lss[i];
    case LogSplineForm::ELEC_PME_DIRECT_EXCL:
    case LogSplineForm::DELEC_PME_DIRECT:
    case LogSplineForm::DELEC_PME_DIRECT_EXCL:
    case LogSplineForm::CUSTOM:
      break;
    }
  }
  rtErr("No potential function was provided.", "PPITable", "findNonExclPotential");
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename T4>
const LogScaleSpline<T4>& PPITable::findExclPotential(const LogScaleSpline<T4> &spl_a,
                                                      const LogScaleSpline<T4> &spl_b,
                                                      const LogScaleSpline<T4> &spl_c,
                                                      const LogScaleSpline<T4> &spl_d) const {
  const std::vector<const LogScaleSpline<T4>*> lss = { spl_a.getSelfPointer(),
                                                       spl_b.getSelfPointer(),
                                                       spl_c.getSelfPointer(),
                                                       spl_d.getSelfPointer() };
  for (int i = 0; i < 4; i++) {
    switch (lss[i]->getForm()) {
    case LogSplineForm::ELEC_PME_DIRECT_EXCL:
      return *lss[i];
    case LogSplineForm::ELEC_PME_DIRECT:
    case LogSplineForm::DELEC_PME_DIRECT:
    case LogSplineForm::DELEC_PME_DIRECT_EXCL:
    case LogSplineForm::CUSTOM:
      break;
    }
  }
  rtErr("No potential function (with exclusions) was provided.", "PPITable",
        "findNonExclPotential");
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename T4>
const LogScaleSpline<T4>& PPITable::findNonExclForce(const LogScaleSpline<T4> &spl_a,
                                                     const LogScaleSpline<T4> &spl_b,
                                                     const LogScaleSpline<T4> &spl_c,
                                                     const LogScaleSpline<T4> &spl_d) const {
  const std::vector<const LogScaleSpline<T4>*> lss = { spl_a.getSelfPointer(),
                                                       spl_b.getSelfPointer(),
                                                       spl_c.getSelfPointer(),
                                                       spl_d.getSelfPointer() };
  for (int i = 0; i < 4; i++) {
    switch (lss[i]->getForm()) {
    case LogSplineForm::DELEC_PME_DIRECT:
      return *lss[i];
    case LogSplineForm::ELEC_PME_DIRECT:
    case LogSplineForm::ELEC_PME_DIRECT_EXCL:
    case LogSplineForm::DELEC_PME_DIRECT_EXCL:
    case LogSplineForm::CUSTOM:
      break;
    }
  }
  rtErr("No force function was provided.", "PPITable", "findNonExclPotential");
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename T4>
const LogScaleSpline<T4>& PPITable::findExclForce(const LogScaleSpline<T4> &spl_a,
                                                  const LogScaleSpline<T4> &spl_b,
                                                  const LogScaleSpline<T4> &spl_c,
                                                  const LogScaleSpline<T4> &spl_d) const {
  const std::vector<const LogScaleSpline<T4>*> lss = { spl_a.getSelfPointer(),
                                                       spl_b.getSelfPointer(),
                                                       spl_c.getSelfPointer(),
                                                       spl_d.getSelfPointer() };
  for (int i = 0; i < 4; i++) {
    switch (lss[i]->getForm()) {
    case LogSplineForm::DELEC_PME_DIRECT_EXCL:
      return *lss[i];
    case LogSplineForm::ELEC_PME_DIRECT:
    case LogSplineForm::ELEC_PME_DIRECT_EXCL:
    case LogSplineForm::DELEC_PME_DIRECT:
    case LogSplineForm::CUSTOM:
      break;
    }
  }
  rtErr("No force function (with exclusions) was provided.", "PPITable", "findNonExclPotential");
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename T4> uint PPITable::checkPriority(const LogScaleSpline<T4> &spl_x) const {
  switch (spl_x.getForm()) {
  case LogSplineForm::ELEC_PME_DIRECT:
    return 0x1;
  case LogSplineForm::ELEC_PME_DIRECT_EXCL:
    return 0x2;
  case LogSplineForm::DELEC_PME_DIRECT:
    return 0x4;
  case LogSplineForm::DELEC_PME_DIRECT_EXCL:
    return 0x8;
  case LogSplineForm::CUSTOM:
    rtErr("Custom potentials cannot be compiled into this object.", "PPITable", "checkPriority");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename T4> void PPITable::checkSplineCompatibility(const LogScaleSpline<T4> &spl_a,
                                                               const LogScaleSpline<T4> &spl_b) {
  if (spl_a.getIndexingMethod() != spl_b.getIndexingMethod()) {
    rtErr("Spline tables with different indexing methods (" +
          getEnumerationName(spl_a.getIndexingMethod()) + ", " +
          getEnumerationName(spl_b.getIndexingMethod()) + ") cannot be compiled together.",
          "PPITable", "checkSplineCompatibility");
  }
  if (spl_a.getBasisSet() != spl_b.getBasisSet()) {
    rtErr("Spline tables with different basis sets (" +
          getEnumerationName(spl_a.getBasisSet()) + ", " +
          getEnumerationName(spl_b.getBasisSet()) + ") cannot be compiled together.",
          "PPITable", "checkSplineCompatibility");
  }
  if (spl_a.getBitStride() != spl_b.getBitStride()) {
    rtErr("Spline tables with different bit strides (" + std::to_string(spl_a.getBitStride()) +
          ", " + std::to_string(spl_b.getBitStride()) + ") cannot be compiled together.",
          "PPITable", "checkSplineCompatibility");
  }
  if (fabs(spl_a.getEwaldCoefficient() - spl_b.getEwaldCoefficient()) > constants::tiny) {
    rtErr("Spline tables with different Gaussian particle widths (" +
          realToString(0.5 / spl_a.getEwaldCoefficient(), 9, 6, NumberFormat::STANDARD_REAL) +
          ", " + 
          realToString(0.5 / spl_b.getEwaldCoefficient(), 9, 6, NumberFormat::STANDARD_REAL) +
          ") cannot be compiled together.", "PPITable", "checkSplineCompatibility");
  }
  if (fabs(spl_a.getCoulombConstant() - spl_b.getCoulombConstant()) > constants::tiny) {
    rtErr("Spline tables with different variations of Coulomb's constant (" +
          realToString(spl_a.getCoulombConstant(), 10, 6, NumberFormat::STANDARD_REAL) + ", " + 
          realToString(spl_b.getCoulombConstant(), 10, 6, NumberFormat::STANDARD_REAL) +
          ") cannot be compiled together.", "PPITable", "checkSplineCompatibility");
  }
  if (fabs(spl_a.getMaximumRange() - spl_b.getMaximumRange()) > constants::tiny) {
    rtErr("Spline tables with different ranges (" + 
          realToString(spl_a.getMaximumRange(), 10, 6, NumberFormat::STANDARD_REAL) + ", " + 
          realToString(spl_b.getMaximumRange(), 10, 6, NumberFormat::STANDARD_REAL) +
          ") cannot be compiled together.", "PPITable", "checkSplineCompatibility");
  }
  if (spl_a.getIndexingOffset() != spl_b.getIndexingOffset()) {
    rtErr("Spline tables with different indexing offsets (" +
          std::to_string(spl_a.getIndexingOffset()) + ", " +
          std::to_string(spl_a.getIndexingOffset()) + ") cannot be compiled together.",
          "PPITable", "checkSplineCompatibility");
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T4>
LogScaleSpline<T4> PPITable::getTablePriority(const LogScaleSpline<T4> &spl_a,
                                              const LogScaleSpline<T4> &spl_b,
                                              const LogScaleSpline<T4> &spl_c) const {
  checkSplineCompatibility(spl_a, spl_b);
  checkSplineCompatibility(spl_b, spl_c);
  const uint holdings = checkPriority(spl_a) | checkPriority(spl_b) | checkPriority(spl_c);
  const LogSplineForm missing_form = findMissingForm(holdings);
  return LogScaleSpline<T4>(missing_form, spl_a.getEwaldCoefficient(), spl_a.getCoulombConstant(),
                            spl_a.getBitStride(), spl_a.getMaximumRange(), spl_a.getMinimumRange(),
                            spl_a.getIndexingMethod(), spl_a.getBasisSet(),
                            spl_a.getOptimization_depth, spl_a.getIndexingOffset());
}

//-------------------------------------------------------------------------------------------------
template <typename T4>
LogScaleSpline<T4> PPITable::getTablePriority(const LogScaleSpline<T4> &spl_a,
                                              const LogScaleSpline<T4> &spl_b) const {
  checkSplineCompatibility(spl_a, spl_b);
  const uint holdings = checkPriority(spl_a) | checkPriority(spl_b);
  const LogSplineForm missing_form = findMissingForm(holdings);
  return LogScaleSpline<T4>(missing_form, spl_a.getEwaldCoefficient(), spl_a.getCoulombConstant(),
                            spl_a.getBitStride(), spl_a.getMaximumRange(), spl_a.getMinimumRange(),
                            spl_a.getIndexingMethod(), spl_a.getBasisSet(),
                            spl_a.getOptimization_depth, spl_a.getIndexingOffset());
}

//-------------------------------------------------------------------------------------------------
template <typename T4>
LogScaleSpline<T4> PPITable::getTablePriority(const LogScaleSpline<T4> &spl_a) const {
  const uint holdings = checkPriority(spl_a);
  const LogSplineForm missing_form = findMissingForm(holdings);
  return LogScaleSpline<T4>(missing_form, spl_a.getEwaldCoefficient(), spl_a.getCoulombConstant(),
                            spl_a.getBitStride(), spl_a.getMaximumRange(), spl_a.getMinimumRange(),
                            spl_a.getIndexingMethod(), spl_a.getBasisSet(),
                            spl_a.getOptimization_depth, spl_a.getIndexingOffset());
}

//-------------------------------------------------------------------------------------------------
template <typename T4>
void PPITable::populateCoefficients(const LogScaleSpline<T4> &u, const LogScaleSpline<T4> &du,
                                    const LogScaleSpline<T4> &ux, const LogScaleSpline<T4> &dux) {
  const LogSplineTable tbl_u   = u.data();
  const LogSplineTable tbl_ux  = ux.data();
  const LogSplineTable tbl_du  = du.data();
  const LogSplineTable tbl_dux = dux.data();

  // Allocate a common space for all four potentials.
  const size_t ct = std::type_index(typeid(T4)).hash_code();
  const size_t tbl_len = u.getSplineIndex(u.getMaximumRange());
  coeffs.resize(4 * exclusion_offset);
  double4* coef_ptr = coeffs.data();
  float4* sp_coef_ptr = sp_coeffs.data();
  if (ct == double4_type_index) {
    for (int i = 0; i < tbl_len; i++) {
      coef_ptr[                i].x = tbl_u.table[i].x;
      coef_ptr[                i].y = tbl_u.table[i].y;
      coef_ptr[                i].z = tbl_u.table[i].z;
      coef_ptr[                i].w = tbl_u.table[i].w;
      coef_ptr[      tbl_len + i].x = tbl_du.table[i].x;
      coef_ptr[      tbl_len + i].y = tbl_du.table[i].y;
      coef_ptr[      tbl_len + i].z = tbl_du.table[i].z;
      coef_ptr[      tbl_len + i].w = tbl_du.table[i].w;
      coef_ptr[(2 * tbl_len) + i].x = tbl_ux.table[i].x;
      coef_ptr[(2 * tbl_len) + i].y = tbl_ux.table[i].y;
      coef_ptr[(2 * tbl_len) + i].z = tbl_ux.table[i].z;
      coef_ptr[(2 * tbl_len) + i].w = tbl_ux.table[i].w;
      coef_ptr[(3 * tbl_len) + i].x = tbl_dux.table[i].x;
      coef_ptr[(3 * tbl_len) + i].y = tbl_dux.table[i].y;
      coef_ptr[(3 * tbl_len) + i].z = tbl_dux.table[i].z;
      coef_ptr[(3 * tbl_len) + i].w = tbl_dux.table[i].w;
    }
  }
  else if (ct == float4_type_index) {
    for (int i = 0; i < tbl_len; i++) {
      sp_coef_ptr[                i].x = tbl_u.table[i].x;
      sp_coef_ptr[                i].y = tbl_u.table[i].y;
      sp_coef_ptr[                i].z = tbl_u.table[i].z;
      sp_coef_ptr[                i].w = tbl_u.table[i].w;
      sp_coef_ptr[      tbl_len + i].x = tbl_du.table[i].x;
      sp_coef_ptr[      tbl_len + i].y = tbl_du.table[i].y;
      sp_coef_ptr[      tbl_len + i].z = tbl_du.table[i].z;
      sp_coef_ptr[      tbl_len + i].w = tbl_du.table[i].w;
      sp_coef_ptr[(2 * tbl_len) + i].x = tbl_ux.table[i].x;
      sp_coef_ptr[(2 * tbl_len) + i].y = tbl_ux.table[i].y;
      sp_coef_ptr[(2 * tbl_len) + i].z = tbl_ux.table[i].z;
      sp_coef_ptr[(2 * tbl_len) + i].w = tbl_ux.table[i].w;
      sp_coef_ptr[(3 * tbl_len) + i].x = tbl_dux.table[i].x;
      sp_coef_ptr[(3 * tbl_len) + i].y = tbl_dux.table[i].y;
      sp_coef_ptr[(3 * tbl_len) + i].z = tbl_dux.table[i].z;
      sp_coef_ptr[(3 * tbl_len) + i].w = tbl_dux.table[i].w;
    }
  }
  exclusion_offset = 2 * tbl_len;
  energy.setPointer(&coeffs,                           0, tbl_len);
  force.setPointer(&coeffs,                      tbl_len, tbl_len);
  energy_with_exclusions.setPointer(&coeffs, 2 * tbl_len, tbl_len);
  force_with_exclusions.setPointer(&coeffs,  3 * tbl_len, tbl_len);
}

} // namespace energy
} // namespace stormm
