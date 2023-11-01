// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace energy {

//-------------------------------------------------------------------------------------------------
template <typename T4>
PPITable::PPITable(const LogScaleSpline<T4> &spl_a, const LogScaleSpline<T4> &spl_b,
                   const LogScaleSpline<T4> &spl_c, const LogScaleSpline<T4> &spl_d,
                   const double cutoff_in) :
    theme{NonbondedTheme::ELECTROSTATIC},
    precision{(std::type_index(typeid).hash_code() == double_type_index) ?
              PrecisionModel::DOUBLE : PrecisionModel::SINGLE},
    basis_set{spl_a.getBasisSet()},
    indexing_method{spl_a.getIndexingMethod()},
    cutoff{cutoff_in},
    dsum_tol{recoverDirectSumTolerance(spl_a.getEwaldCoefficient(), cutoff_in)},
    mantissa_bits{spl_a.getBitStride()},
    coulomb{spl_a.getCoulombConstant()}
{
  checkSplineCompatibility(spl_a, spl_b);
  checkSplineCompatibility(spl_a, spl_c);
  checkSplineCompatibility(spl_a, spl_d);

  // Determine which potential is which
  const LogScaleSpline &u   = findNonExclPotential(spl_a, spl_b, spl_c, spl_d);
  const LogScaleSpline &ux  = findExclPotential(spl_a, spl_b, spl_c, spl_d);
  const LogScaleSpline &du  = findNonExclForce(spl_a, spl_b, spl_c, spl_d);
  const LogScaleSpline &dux = findExclForce(spl_a, spl_b, spl_c, spl_d);
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
template <typename T4>
void PPITable::checkSplineCompatibility(const LogScaleSpline<T4> &spl_a,
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
template <typename T4> LogSplineForm PPITable::findMissingForm(const uint holdings) const {
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
                                              const LogScaleSpline<T4> &spl_b) {
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
LogScaleSpline<T4> PPITable::getTablePriority(const LogScaleSpline<T4> &spl_a) {
  const uint holdings = checkPriority(spl_a);
  const LogSplineForm missing_form = findMissingForm(holdings);
  return LogScaleSpline<T4>(missing_form, spl_a.getEwaldCoefficient(), spl_a.getCoulombConstant(),
                            spl_a.getBitStride(), spl_a.getMaximumRange(), spl_a.getMinimumRange(),
                            spl_a.getIndexingMethod(), spl_a.getBasisSet(),
                            spl_a.getOptimization_depth, spl_a.getIndexingOffset());
}

} // namespace energy
} // namespace stormm
