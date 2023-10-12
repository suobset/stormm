// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace stmath {

//-------------------------------------------------------------------------------------------------
template <typename T4>
LogSplineTable<T4>::LogSplineTable(const int detail_bits_in, const int index_bound_in,
                                   const uint sp_detail_mask_in, const ullint dp_detail_mask_in,
                                   const T4* table_in) :
    detail_bits{detail_bits_in}, index_bound{index_bound_in}, sp_detail_mask{sp_detail_mask_in},
    dp_detail_mask{dp_detail_mask_in}, table{table_in}
{}

//-------------------------------------------------------------------------------------------------
template <typename T4>
LogScaleSpline<T4>::LogScaleSpline(const TableIndexing indexing_method_in,
                                   const BasisFunctions basis_set_in,
                                   const float indexing_offset_in) :
    target_form{LogSplineForm::CUSTOM},
    indexing_method{indexing_method_in},
    basis_set{basis_set_in},
    mantissa_bits{default_logtab_mantissa_bits},
    segments_per_stride{static_cast<int>(round(pow(2.0, mantissa_bits)))},
    ewald_coefficient{0.0},
    coulomb_constant{amber_ancient_bioq},
    maximum_range{default_logtab_max_range * default_logtab_max_range},
    minimum_absolute_range{default_logtab_min_range * default_logtab_min_range},
    minimum_significant_range{0.5},
    indexing_offset{checkIndexingOffset(indexing_offset_in)},
    ulp_optimization_depth{2},
    precision{PrecisionModel::SINGLE},
    table{},
    mean_overall_error{0.0}, stdev_overall_error{0.0}, max_overall_error{0.0},
    mean_segment_error{},
    stdev_segment_error{},
    max_segment_error{}
{
  // Set the precision model and check the data type
  const size_t ct = std::type_index(typeid(T4)).hash_code();
  if (ct == float4_type_index) {
    precision = PrecisionModel::SINGLE;
  }
  else if (ct == double4_type_index) {
    precision = PrecisionModel::DOUBLE;
  }
  else {
    rtErr("The only allowed table types are float and double four-tuples.", "LogScaleSpline");
  }
  
}

//-------------------------------------------------------------------------------------------------
template <typename T4>
LogScaleSpline<T4>::LogScaleSpline(const LogSplineForm target_form_in,
                                   const double ewald_coefficient_in,
                                   const double coulomb_constant_in, const int mantissa_bits_in,
                                   const double max_range_in, const double min_range_in,
                                   const TableIndexing indexing_method_in,
                                   const BasisFunctions basis_set_in,
                                   const float indexing_offset_in,
                                   const ExceptionResponse policy) :
    LogScaleSpline(indexing_method_in, basis_set_in, indexing_offset_in)
{
  // Set additional constants.  Coulomb's constant is accepted, whereas the Ewald coefficient will
  // be checked for validity. 
  setTargetForm(target_form_in);
  setEwaldCoefficient(ewald_coefficient_in, policy);
  setMantissaBits(mantissa_bits_in, policy);
  setMinimumAbsoluteRange(min_range_in);
  allocate(max_range_in);
  switch (target_form) {
  case LogSplineForm::ELEC_PME_DIRECT:
  case LogSplineForm::ELEC_PME_DIRECT_EXCL:
  case LogSplineForm::DELEC_PME_DIRECT:
  case LogSplineForm::DELEC_PME_DIRECT_EXCL:
    coulomb_constant = coulomb_constant_in;
    evaluatePrescribedTable();
    break;
  case LogSplineForm::CUSTOM:

    // The case of custom potentials, which cannot behandled by this constructor, will have been
    // trapped in setTargetForm().
    break;
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T4>
LogScaleSpline<T4>::LogScaleSpline(const LogSplineForm target_form_in,
                                   const std::vector<std::vector<double2>> &custom_form_in,
                                   const int mantissa_bits_in, const double max_range_in,
                                   const double min_range_in,
                                   const TableIndexing indexing_method_in,
                                   const BasisFunctions basis_set_in,
                                   const float indexing_offset_in,
                                   const ExceptionResponse policy) :
    LogScaleSpline(indexing_method_in, basis_set_in, indexing_offset_in)
{
  setTargetForm(target_form_in, custom_form_in);
  setMantissaBits(mantissa_bits_in, policy);
  setMinimumAbsoluteRange(min_range_in);
  allocate(max_range_in, custom_form_in);
  evaluateCustomTable();
}

//-------------------------------------------------------------------------------------------------
template <typename T4>
int LogScaleSpline<T4>::getBitStride() const {
  return mantissa_bits;
}

//-------------------------------------------------------------------------------------------------
template <typename T4>
int LogScaleSpline<T4>::getSplineDensity() const {
  return segments_per_stride;
}

//-------------------------------------------------------------------------------------------------
template <typename T4>
double LogScaleSpline<T4>::getEwaldCoefficient() const {
  return ewald_coefficient;
}

//-------------------------------------------------------------------------------------------------
template <typename T4>
double LogScaleSpline<T4>::getCoulombConstant() const {
  return coulomb_constant;
}

//-------------------------------------------------------------------------------------------------
template <typename T4>
double LogScaleSpline<T4>::getMaximumRange() const {
  return maximum_range;
}

//-------------------------------------------------------------------------------------------------
template <typename T4>
float LogScaleSpline<T4>::getIndexingOffset() const {
  return indexing_offset;
}

//-------------------------------------------------------------------------------------------------
template <typename T4>
int LogScaleSpline<T4>::getSplineIndex(const double r) const {
  int result;
  switch (precision) {
  case PrecisionModel::DOUBLE:
    {
      const Ecumenical8 workspc = { .d = r };
      result = (workspc.ulli >> (52 - mantissa_bits));
    }
    break;
  case PrecisionModel::SINGLE:
    {
      const Ecumenical4 workspc = { .f = static_cast<float>(r) };
      result = (workspc.ui >> (23 - mantissa_bits));
    }
    break;
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename T4>
int LogScaleSpline<T4>::getSplineIndexByRealArg(const double r) const {
  switch (indexing_method) {
  case TableIndexing::ARG:
    return getSplineIndex(r);
  case TableIndexing::SQUARED_ARG:
    switch (precision) {
    case PrecisionModel::DOUBLE:
      return getSplineIndex(r * r);
    case PrecisionModel::SINGLE:
      {
        const float rf = r;
        return getSplineIndex(rf * rf);
      }
      break;
    }
    break;
  case TableIndexing::ARG_OFFSET:
    switch (precision) {
    case PrecisionModel::DOUBLE:
      return getSplineIndex((r * r) + static_cast<double>(indexing_offset));
    case PrecisionModel::SINGLE:
      {
        const float rf = r;
        return getSplineIndex((rf * rf) + indexing_offset);
      }
      break;
    }    
    break;
  case TableIndexing::SQ_ARG_OFFSET:
    switch (precision) {
    case PrecisionModel::DOUBLE:
      return getSplineIndex(r + static_cast<double>(indexing_offset));
    case PrecisionModel::SINGLE:
      {
        const float rf = r;
        return getSplineIndex(rf + indexing_offset);
      }
      break;
    }    
    break;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename T4>
double3 LogScaleSpline<T4>::getErrorEstimate() const {
  return { mean_overall_error, stdev_overall_error, max_overall_error };
}

//-------------------------------------------------------------------------------------------------
template <typename T4>
double3 LogScaleSpline<T4>::getErrorEstimate(const double r) const {
  const size_t tbl_idx = getSplineIndex(r);
  return { mean_segment_error[tbl_idx], stdev_segment_error[tbl_idx], max_segment_error[tbl_idx] };
}

//-------------------------------------------------------------------------------------------------
template <typename T4>
double LogScaleSpline<T4>::evaluate(const double r) const {
  const size_t coef_idx = getSplineIndex(r);
  const T4 coefs = table.readHost(coef_idx);
  double result;
  switch (basis_set) {
  case BasisFunctions::MIXED_FRACTIONS:
    switch (precision) {
    case PrecisionModel::DOUBLE:
      {
        const double invr = 1.0 / r;
        result = (coefs.x * r) + coefs.y + ((coefs.z + (coefs.w * invr)) * invr);
      }
      break;
    case PrecisionModel::SINGLE:
      {
        const float rf = r;
        const float invrf = 1.0f / rf;
        result = (coefs.x * rf) + coefs.y + ((coefs.z + (coefs.w * invrf)) * invrf);
      }
      break;
    }
    break;
  case BasisFunctions::POLYNOMIAL:
    switch (precision) {
    case PrecisionModel::DOUBLE:
      {
        Ecumenical8 xfrm = { .d = r };
        xfrm.ulli = (((xfrm.ulli & dp_detail_bitmask) << mantissa_bits) | 0x3ff0000000000000ULL);
        double dr = xfrm.d;
        result = coefs.x + ((coefs.y + ((coefs.z + (coefs.w * dr)) * dr)) * dr);
      }
      break;
    case PrecisionModel::SINGLE:
      {
        // This function is used to evaluate the accuracy of the spline and perform unit tests.
        // As such, it is important to consider whether this conversion of the double-precision
        // input to single-precision will impose an artifact that would differ from conditions in
        // real applications.  However, if the original range argument was computed in
        // single-precision and then converted to double-precision, as it would be in situations
        // where the single-precision spline table is relevant, the double type number preserves
        // all information present in the original float, and reproduces the original float here.
        Ecumenical4 xfrm = { .f = static_cast<float>(r) };
        xfrm.ui = (((xfrm.ui & sp_detail_bitmask) << mantissa_bits) | 0x3f800000U);
        float dr = xfrm.f;
        result = coefs.x + ((coefs.y + ((coefs.z + (coefs.w * dr)) * dr)) * dr);
      }
      break;
    }
    break;
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename T4>
double LogScaleSpline<T4>::evaluateByRealArg(const double r) const {
  switch (indexing_method) {
  case TableIndexing::ARG:
    return evaluate(r);
  case TableIndexing::SQUARED_ARG:
    switch (precision) {
    case PrecisionModel::DOUBLE:
      return evaluate(r * r);
    case PrecisionModel::SINGLE:
      {
        const float rf = r;
        return evaluate(rf * rf);
      }
      break;
    }
    break;
  case TableIndexing::ARG_OFFSET:
    switch (precision) {
    case PrecisionModel::DOUBLE:
      return evaluate(r + static_cast<double>(indexing_offset));
    case PrecisionModel::SINGLE:
      {
        const float rf = r;
        return evaluate(rf + indexing_offset);
      }
      break;
    }
    break;
  case TableIndexing::SQ_ARG_OFFSET:
    switch (precision) {
    case PrecisionModel::DOUBLE:
      return evaluate((r * r) + static_cast<double>(indexing_offset));
    case PrecisionModel::SINGLE:
      {
        const float rf = r;
        return evaluate((rf * rf) + indexing_offset);
      }
      break;
    }
    break;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename T4>
const LogSplineTable<T4> LogScaleSpline<T4>::data(const HybridTargetLevel tier) const {

  // The number of bits in the exponent is getting one added to it, to account for the sign bit.
  int expsgn_bits;
  switch (precision) {
  case PrecisionModel::DOUBLE:
    expsgn_bits = 12;
  case PrecisionModel::SINGLE:
    expsgn_bits = 9;
  }
  return LogSplineTable<T4>(mantissa_bits + expsgn_bits, table.size(), sp_detail_bitmask,
                            dp_detail_bitmask, table.data(tier));
}

//-------------------------------------------------------------------------------------------------
template <typename T4> const LogSplineTable<void>
LogScaleSpline<T4>::templateFreeData(const HybridTargetLevel tier) const {
  int expsgn_bits;
  switch (precision) {
  case PrecisionModel::DOUBLE:
    expsgn_bits = 12;
  case PrecisionModel::SINGLE:
    expsgn_bits = 9;
  }
  return LogSplineTable<T4>(mantissa_bits + expsgn_bits, table.size(), sp_detail_bitmask,
                            dp_detail_bitmask, reinterpret_cast<const void*>(table.data(tier)));
}

#ifdef STORMM_USE_HPC
//-------------------------------------------------------------------------------------------------
template <typename T4>
void LogScaleSpline<T4>::upload() {
  table.upload();
}

//-------------------------------------------------------------------------------------------------
template <typename T4>
void LogScaleSpline<T4>::download() {
  table.download();
}
#endif

//-------------------------------------------------------------------------------------------------
template <typename T4>
void LogScaleSpline<T4>::setMinimumSignificantRange(const double min_range_in) {
  minimum_significant_range = min_range_in;

}

//-------------------------------------------------------------------------------------------------
template <typename T4>
void LogScaleSpline<T4>::setTargetForm(const LogSplineForm target_form_in,
                                       const std::vector<std::vector<double2>> &custom_form_in) {
  switch (target_form_in) {
  case LogSplineForm::ELEC_PME_DIRECT:
  case LogSplineForm::ELEC_PME_DIRECT_EXCL:
  case LogSplineForm::DELEC_PME_DIRECT:
  case LogSplineForm::DELEC_PME_DIRECT_EXCL:
    if (custom_form_in.size() > 0) {
      rtErr("A " + getEnumerationName(target_form) + " functional form is hard-coded and should "
            "not be created in the context of custom array data.", "LogScaleSpline");
    }
    break;
  case LogSplineForm::CUSTOM:
    if (custom_form_in.size() == 0) {
      rtErr("A " + getEnumerationName(target_form) + " functional form requires samples of the "
            "function and its values in order to generate a spline.", "LogScaleSpline");
    }
    break;
  }
  target_form = target_form_in;
}

//-------------------------------------------------------------------------------------------------
template <typename T4>
void LogScaleSpline<T4>::setMantissaBits(const int mantissa_bits_in,
                                         const ExceptionResponse policy) {
  mantissa_bits = mantissa_bits_in;
  if (mantissa_bits < 1 || mantissa_bits > 12) {
    switch (policy) {
    case ExceptionResponse::DIE:
      if (mantissa_bits < 1) {
        rtErr("A table of " + std::to_string(mantissa_bits) + " would be too coarse to give "
              "accurate representations of the function.  Only one spline segment would cover the "
              "range [ 4, 8 ), for example, and there is no reason to go lower.", "LogScaleSpline",
              "setMantissaBits");
      }
      else {
        rtErr("A table of " + std::to_string(mantissa_bits) + " would be so fine as to make every "
              "array access a cache miss on many architectures.  It could be extremely accurate, "
              "but 2048 or more spline segments would cover the range [ 1, 2 ) for example.",
              "LogScaleSpline", "setMantissaBits");
      }
      break;
    case ExceptionResponse::WARN:
      if (mantissa_bits < 1) {
        rtErr("A table of " + std::to_string(mantissa_bits) + " would be too coarse to give "
              "accurate representations of the function.  The default of " +
              std::to_string(default_logtab_mantissa_bits) + " will be restored.",
              "LogScaleSpline", "setMantissaBits");
      }
      else {
        rtErr("A table of " + std::to_string(mantissa_bits) + " would be so large as to defeat "
              "the purpose of generating splines.  The default of " +
              std::to_string(default_logtab_mantissa_bits) + " will be restored.",
              "LogScaleSpline", "setMantissaBits");
      }
      break;
    case ExceptionResponse::SILENT:
      break;
    }      
    mantissa_bits = default_logtab_mantissa_bits;
  }

  // The number of segments per stride is implied by the number of mantissa bits
  segments_per_stride = round(pow(2, mantissa_bits));

  // Compute the detail masks for single- and double-precision
  dp_detail_bitmask = 0LLU;
  for (int i = 0; i < 52 - mantissa_bits; i++) {
    dp_detail_bitmask |= (0x1LLU << i);
  }
  sp_detail_bitmask = 0U;
  for (int i = 0; i < 23 - mantissa_bits; i++) {
    sp_detail_bitmask |= (0x1U << i);
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T4>
void LogScaleSpline<T4>::setMinimumAbsoluteRange(const double min_range_in) {

  // The minimum applicable range of the spline is a speed optimization.  If it is too small then
  // simply set it to zero.
  const Ecumenical8 lowest_bar = { .ulli = (static_cast<ullint>(mantissa_bits + 4) << 52) };
  minimum_absolute_range = (min_range_in < lowest_bar.d) ? lowest_bar.d : min_range_in;
}

//-------------------------------------------------------------------------------------------------
template <typename T4>
float LogScaleSpline<T4>::checkIndexingOffset(float indexing_offset_in) {

  // Check that the number is positive, and greater than zero.
  Ecumenical4 nl = { .f = indexing_offset_in };
  if (((nl.ui >> 31) & 0x1) == 1 || (nl.ui & 0x7fffff) > 0) {
    rtErr("The indexing offset must be a positive number and some small power (or small "
          "negative power) of two.", "LogScaleSpline", "checkIndexingOffset");
  }
  return indexing_offset_in;
}

//-------------------------------------------------------------------------------------------------
template <typename T4>
void LogScaleSpline<T4>::setEwaldCoefficient(const double ewald_coefficient_in,
                                             const ExceptionResponse policy) {
  if (ewald_coefficient_in < minimum_ewald_coefficient ||
      ewald_coefficient >= maximum_ewald_coefficient) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("An Ewald coefficient of " +
            realToString(ewald_coefficient_in, 9, 5, NumberFormat::STANDARD_REAL) +
            " represents a Gaussian spread of " +
            realToString(0.5 / ewald_coefficient_in, 9, 5, NumberFormat::STANDARD_REAL) + ".",
            "LogScaleSpline", "setEwaldCoefficient");
    case ExceptionResponse::WARN:
      rtWarn("An Ewald coefficient of " +
             realToString(ewald_coefficient_in, 9, 5, NumberFormat::STANDARD_REAL) +
             " represents a Gaussian spread of " +
             realToString(0.5 / ewald_coefficient_in, 9, 5, NumberFormat::STANDARD_REAL) +
             ".  This seems unreasonable but will be tolerated.", "LogScaleSpline",
             "setEwaldCoefficient");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
  }
  ewald_coefficient = ewald_coefficient_in;
}

//-------------------------------------------------------------------------------------------------
template <typename T4>
void LogScaleSpline<T4>::allocate(const double max_range_in,
                                  const std::vector<std::vector<double2>> &custom_form_in,
                                  const ExceptionResponse policy) {

  // Set the maximum range according to the indexing method.
  switch (indexing_method) {
  case TableIndexing::ARG:
    maximum_range = max_range_in;
    break;
  case TableIndexing::SQUARED_ARG:
    maximum_range = max_range_in * max_range_in;
    break;
  case TableIndexing::ARG_OFFSET:
    maximum_range = max_range_in + indexing_offset;
    break;
  case TableIndexing::SQ_ARG_OFFSET:
    maximum_range = (max_range_in * max_range_in) + indexing_offset;
    break;
  }
  
  // Determine the size of the table.  This may modify the maximum_range member variable if using
  // a custom potential.  Allocate memory for the table.
  switch (target_form) {
  case LogSplineForm::ELEC_PME_DIRECT:
  case LogSplineForm::ELEC_PME_DIRECT_EXCL:
  case LogSplineForm::DELEC_PME_DIRECT:
  case LogSplineForm::DELEC_PME_DIRECT_EXCL:
    maximum_range = pow(2.0, ceil(log2(maximum_range)));
    break;
  case LogSplineForm::CUSTOM:
    {
      double observed_range = 0.0;
      for (size_t i = 0; i < custom_form_in.size(); i++) {
        const size_t tbl_length = custom_form_in[i].size();
        for (size_t j = 0; j < tbl_length; j++) {
          observed_range = std::max(custom_form_in[i][j].x, observed_range);
        }
      }
      if (observed_range < 0.0) {
        rtErr("The interpolated range in a logarithmically indexed table must be greater than "
              "zero.", "LogScaleSpline");
      }
      const double old_maximum_range = maximum_range;
      maximum_range = std::min(observed_range, maximum_range);
      if (maximum_range / old_maximum_range < 0.95) {
        switch (policy) {
        case ExceptionResponse::DIE:
          rtErr("The data provided are insufficient to cover a range of " +
                realToString(old_maximum_range, 9, 4, NumberFormat::STANDARD_REAL) + ".",
                "LogScaleSpline");
        case ExceptionResponse::WARN:
          rtWarn("The data provided are insufficient to cover a range of " +
                 realToString(old_maximum_range, 9, 4, NumberFormat::STANDARD_REAL) +
                 ".  A maximum range of " +
                 realToString(maximum_range, 9, 4, NumberFormat::STANDARD_REAL) + " will be "
                 "applied instead.", "LogScaleSpline");
          break;
        case ExceptionResponse::SILENT:
          break;
        }
      }
    }
    break;
  }
  const size_t max_rng_idx = getSplineIndex(maximum_range);
  table.resize(max_rng_idx + 1);
  mean_segment_error.resize(max_rng_idx + 1);
  stdev_segment_error.resize(max_rng_idx + 1);
  max_segment_error.resize(max_rng_idx + 1);
}

//-------------------------------------------------------------------------------------------------
template <typename T4>
void LogScaleSpline<T4>::evaluatePrescribedTable() {

  // Pre-determine whether a Coulombic exclusion term will be added to the potential and its
  // derivative.
  bool add_clmb_exclusion;
  switch (target_form) {
  case LogSplineForm::ELEC_PME_DIRECT:
  case LogSplineForm::DELEC_PME_DIRECT:
    add_clmb_exclusion = false;
    break;
  case LogSplineForm::ELEC_PME_DIRECT_EXCL:
  case LogSplineForm::DELEC_PME_DIRECT_EXCL:
    add_clmb_exclusion = true;
    break;
  case LogSplineForm::CUSTOM:
    rtErr("The spline table coefficients must be evaluated using a separate method for a " +
          getEnumerationName(target_form) + " functional form.", "LogScaleSpline",
          "evaluatePrescribedTable");
    break;
  }

  // Pre-allocate workspace arrays to hold the function value and its derivatives.
  std::vector<double> amat(32), bwrk(8), bref(8), xvec(4);

  // Establish the matrix for solving the spline coefficients in the case of polynomial basis
  // functions.  Lay out the matrix such that the nth column controls the coefficient for x^n.
  // Thus, the expression Ax^3 + Bx^2 + Cx + D evaluates such that D is solved in the first column
  // and A in the last.  Note that the matrix is written out in column format below, despite
  // appearances.  The first row of the matrix is [ 1, 0, 0, 0 ].  In addition, because the spline
  // in indexed by the square root and targets only values of the target function, not its
  // derivatives (the target function can be the derivative of some function of interest, but
  // derivatives of the spline itself are not expected to be used), the table is over-determined
  // with values of the target function.
  std::vector<double> amat_ref(32);
  switch (basis_set) {
  case BasisFunctions::MIXED_FRACTIONS:
    break;
  case BasisFunctions::POLYNOMIAL:
    for (int i = 0; i < 8; i++) {
      const double di = 1.0625 + (0.125 * static_cast<double>(i));
      amat_ref[     i] = 1.0;
      amat_ref[ 8 + i] = di;
      amat_ref[16 + i] = di * di;
      amat_ref[24 + i] = di * di * di;
    }
    break;
  }

  // Compute the absolute minimum range index.
  const int min_abs_rng_idx = getSplineIndex(minimum_absolute_range);

  // Loop over all segments.
  const double dseg_per_stride = segments_per_stride;
  int pos = 0;
  const int tbl_length = table.size();
  const double bfac = 2.0 * ewald_coefficient / sqrt(symbols::pi);
  T4* table_ptr = table.data();
  while (pos < tbl_length) {

    // In double-precision mode, there can be nothing smaller than the smallest exponent, and
    // so the first mantissa_bits strides, plus an additional four bits to account for the need
    // to discretize distances at 1 part in 16 of the spline segment width, will be left as zeros.
    if (pos < min_abs_rng_idx) {
      mean_segment_error[pos] = 0.0;
      stdev_segment_error[pos] = 0.0;
      max_segment_error[pos] = 0.0;
      pos++;
      continue;
    }

    // Determine the length scale.  Evaluate all splines for the segment and log them to the
    // extent possible.
    double stride_start, segment_width;
    switch (precision) {
    case PrecisionModel::DOUBLE:
      {
        const ullint expn_setting = pos / segments_per_stride;
        Ecumenical8 width_work = { .ulli = (expn_setting << 52) };
        stride_start = width_work.d;
        segment_width = stride_start / dseg_per_stride;
      }
      break;
    case PrecisionModel::SINGLE:
      {
        const uint expn_setting = pos / segments_per_stride;
        Ecumenical4 width_work = { .ui = (expn_setting << 23) };
        stride_start = width_work.f;
        segment_width = stride_start / dseg_per_stride;
      }
      break;
    }

    // In some table indexing modes, the vast majority of the table will never be accessed due to
    // the indexing value being raised by an arbitrary, small amount (indexing_offset).  Skip this 
    // segment if the real value of the function under such an indexing offset would be less than
    // or equal to zero.
    switch (indexing_method) {
    case TableIndexing::ARG:
    case TableIndexing::SQUARED_ARG:
      break;
    case TableIndexing::ARG_OFFSET:
    case TableIndexing::SQ_ARG_OFFSET:
      if (stride_start < indexing_offset) {
        mean_segment_error[pos] = 0.0;
        stdev_segment_error[pos] = 0.0;
        max_segment_error[pos] = 0.0;
        pos++;
        continue;
      }
      break;
    }

    // The table is indexed by the square of the function argument.  Find this and take its
    // square root to obtain the function argument for samples at 0.0625, 0.1875, ..., 0.9375
    // over the spline segment.
    const float segment_width_f = segment_width;
    double r2_i = stride_start;
    for (int i = 0; i < segments_per_stride; i++) {
      if (pos + i >= tbl_length) {
        continue;
      }
      
      // Evaluate the function at eight points to make an over-determined system and obtain the
      // best answer in a least-squares sense.
      for (int j = 0; j < 8; j++) {
        double r_ij = r2_i + ((0.0625 + (0.125 * static_cast<double>(j))) * segment_width);
        
        // Revise the fitting matrix for the medley of fractional basis functions.
        switch (basis_set) {
        case BasisFunctions::MIXED_FRACTIONS:
          amat[j     ] = r_ij;
          amat[j +  8] = 1.0;
          amat[j + 16] = 1.0 / r_ij;
          amat[j + 24] = 1.0 / (r_ij * r_ij);
          break;
        case BasisFunctions::POLYNOMIAL:
          break;
        }

        // Adjust the argument as needed to evaluate the actual target function.
        switch (indexing_method) {
        case TableIndexing::ARG:
          break;
        case TableIndexing::SQUARED_ARG:
          r_ij = sqrt(r_ij);
          break;
        case TableIndexing::ARG_OFFSET:
          r_ij = r_ij - static_cast<double>(indexing_offset);
          break;
        case TableIndexing::SQ_ARG_OFFSET:
          if (r_ij - static_cast<double>(indexing_offset) < constants::tiny) {
            r_ij = constants::tiny;
          }
          else {
            r_ij = sqrt(r_ij - static_cast<double>(indexing_offset));
          }
          break;
        }

        // Evaluate the target function for this spline index.
        switch (target_form) {
        case LogSplineForm::ELEC_PME_DIRECT:
          bwrk[j] = coulomb_constant * erfc(ewald_coefficient * r_ij) / r_ij;
          break;
        case LogSplineForm::ELEC_PME_DIRECT_EXCL:
          bwrk[j] = coulomb_constant * (erfc(ewald_coefficient * r_ij) - 1.0) / r_ij;
          break;
        case LogSplineForm::DELEC_PME_DIRECT:
          {
            const double ewr = ewald_coefficient * r_ij;
            bwrk[j] = -coulomb_constant *
                      ((bfac * exp(-ewr * ewr)) + (erfc(ewr) / r_ij)) / (r_ij * r_ij);
          }
          break;
        case LogSplineForm::DELEC_PME_DIRECT_EXCL:
          {
            const double ewr = ewald_coefficient * r_ij;
            bwrk[j] = -coulomb_constant *
                      ((bfac * exp(-ewr * ewr)) + ((erfc(ewr) - 1.0) / r_ij)) / (r_ij * r_ij);
          }
          break;
        case LogSplineForm::CUSTOM:
          break;
        }
        bref[j] = bwrk[j];
      }

      // Refresh the fitting matrix for polynomial basis functions (which work in a relative
      // coordinate frame), or store the reference form of the fitting matrix for the mixed
      // fractional basis functions (which work in an absolute frame of reference).
      switch (basis_set) {
      case BasisFunctions::MIXED_FRACTIONS:
        for (int j = 0; j < 32; j++) {
          amat_ref[j] = amat[j];
        }
        break;
      case BasisFunctions::POLYNOMIAL:
        for (int j = 0; j < 32; j++) {
          amat[j] = amat_ref[j];
        }
        break;
      }
      qrSolver<double>(&amat, &xvec, &bwrk);

      // While it would be safest to initialize this tuple as const and use brace-enclosed
      // initialization, it would require that the function itself take on an extra template
      // parameter, or generate unnecessary compiler warnings.
      T4 spl_tmp;
      spl_tmp.x = xvec[0];
      spl_tmp.y = xvec[1];
      spl_tmp.z = xvec[2];
      spl_tmp.w = xvec[3];

      // Evaluate the quality of the spline over its segment.
      switch (precision) {
      case PrecisionModel::DOUBLE:
        matrixVectorMultiply(amat_ref, xvec, &bwrk, 8, 4, 1.0, 1.0, 0.0);
        break;
      case PrecisionModel::SINGLE:

        // For single-precision spline tables, optimize a number of ULPs of the developer's
        // choosing in order to try and find better fits of the splines, when evaluated in single
        // precision, to the data at hand.
        if (r2_i >= minimum_significant_range) {
          
        }
        
        // Rather than the matrix-vector multiplication, evaluate each point as would occur in a
        // single-precision calculation setting.
        switch (basis_set) {
        case BasisFunctions::MIXED_FRACTIONS:
          for (int j = 0; j < 8; j++) {
            const float r_ij2 = r2_i + ((0.0625f + (0.125f * static_cast<float>(j))) *
                                        segment_width_f);
            bwrk[j] = (spl_tmp.x * r_ij2) + spl_tmp.y +
                      (((spl_tmp.w / r_ij2) + spl_tmp.z) / r_ij2);
          }
          break;
        case BasisFunctions::POLYNOMIAL:

          // This number will still be exact due to its composition from various powers of two
          for (int j = 0; j < 8; j++) {
            const float dj = 1.0625f + 0.125f * static_cast<float>(j);
            bwrk[j] = (((((spl_tmp.w * dj) + spl_tmp.z) * dj) + spl_tmp.y) * dj) + spl_tmp.x;
          }
          break;
        }
        break;
      }
      table_ptr[pos + i] = spl_tmp;
      
      // Analyze the error
      mean_segment_error[pos + i] = meanUnsignedError(bwrk, bref);
      max_segment_error[pos + i] = maxAbsoluteDifference(bwrk, bref);
      for (int j = 0; j < 8; j++) {
        bwrk[j] -= bref[j];
      }
      stdev_segment_error[pos + i] = variance(bwrk, VarianceMethod::STANDARD_DEVIATION);

      // Increment the squared function argument
      r2_i += segment_width;
    }

    // Leap forward to the next stride
    pos += segments_per_stride;
  }

  // Apply the value for the minimum explicitly computed interval to all intervals beneath it.
  for (int i = 0; i < min_abs_rng_idx; i++) {
    table_ptr[i] = table_ptr[min_abs_rng_idx];
  }
  
  // Compute the overall errors
  evaluateOverallError();
}

//-------------------------------------------------------------------------------------------------
template <typename T4>
void LogScaleSpline<T4>::evaluateCustomTable() {
}

//-------------------------------------------------------------------------------------------------
template <typename T4>
void LogScaleSpline<T4>::evaluateOverallError() {

  // Compute the stride in which the minimum significant values begin
  const size_t min_sig_idx = getSplineIndex(minimum_significant_range);
  const size_t npts = table.size() - min_sig_idx;
  mean_overall_error = mean(&mean_segment_error[min_sig_idx], npts);
  stdev_overall_error = variance(&mean_segment_error[min_sig_idx], npts,
                                 VarianceMethod::STANDARD_DEVIATION);
  max_overall_error = maxValue(&max_segment_error[min_sig_idx], npts);
}

//-------------------------------------------------------------------------------------------------
template <typename T4>
LogSplineTable<T4> restoreType(const LogSplineTable<void> *rasa) {
  return LogSplineTable<T4>(rasa->detail_bits, rasa->index_bound, rasa->sp_detail_mask,
                            rasa->dp_detail_mask, reinterpret_cast<const T4*>(rasa->table));
}

//-------------------------------------------------------------------------------------------------
template <typename T4>
LogSplineTable<T4> restoreType(const LogSplineTable<void> &rasa) {
  return LogSplineTable<T4>(rasa.detail_bits, rasa.index_bound, rasa.sp_detail_mask,
                            rasa.dp_detail_mask, reinterpret_cast<const T4*>(rasa.table));
}

} // namespace stmath
} // namespace stormm
