// -*-c++-*-
namespace omni {
namespace math {

using parse::NumberFormat;
using parse::realToString;

//-------------------------------------------------------------------------------------------------
template <typename T> void vectorComparisonCheck(const std::vector<T> &va,
                                                 const std::vector<T> &vb, const char* caller) {

  // Trap comparisons of non-scalar types
  if (isScalarType<T>() == false) {
    rtErr("Comparison between vectors of type " + std::string(typeid(T).name()) +
          " is not permitted.", caller);
  }

  // Trap vector size mismatch
  if (va.size() != vb.size()) {
    rtErr("Comparison requires vectors of identical sizes (" + std::to_string(va.size()) +
          " and " + std::to_string(vb.size()) + " given).", caller);
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T> void vectorComparisonCheck(const Hybrid<T> &va, const Hybrid<T> &vb,
                                                 const char* caller) {

  // Trap comparisons of non-scalar types
  if (isScalarType<T>() == false) {
    rtErr("Comparison between vectors of type " + std::string(typeid(T).name()) +
          " is not permitted.", caller);
  }

  // Trap vector size mismatch
  if (va.size() != vb.size()) {
    rtErr("Comparison requires vectors of identical sizes (" + std::to_string(va.size()) +
          " and " + std::to_string(vb.size()) + " given).", caller);
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T> T maxAbsoluteDifference(const T* va, const T* vb, const size_t length) {
  if (isSignedIntegralScalarType<T>()) {
    T mdev = 0;
    for (size_t i = 0; i < length; i++) {
      T tdev = abs(va[i] - vb[i]);
      if (tdev > mdev) {
        mdev = tdev;
      }
    }
    return mdev;
  }
  else if (isUnsignedIntegralScalarType<T>()) {
    T mdev = 0;
    for (size_t i = 0; i < length; i++) {
      T tdev = (va[i] > vb[i]) ? va[i] - vb[i] : vb[i] - va[i];
      if (tdev > mdev) {
        mdev = tdev;
      }
    }
    return mdev;
  }
  else if (isFloatingPointScalarType<T>()) {
    T mdev = 0;
    for (size_t i = 0; i < length; i++) {
      T tdev = fabs(va[i] - vb[i]);
      if (tdev > mdev) {
        mdev = tdev;
      }
    }
    return mdev;
  }
  else {
    rtErr("Data of type " + std::string(typeid(T).name()) + " is not suitable for comparisons.  "
          "Update isScalarType() and its related subclassification functions.",
          "maxAbsoluteDifference");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename T> T maxAbsoluteDifference(const std::vector<T> &va, const std::vector<T> &vb) {
  vectorComparisonCheck(va, vb, "maxAbsoluteDifference");
  return maxAbsoluteDifference(va.data(), vb.data(), va.size());
}

//-------------------------------------------------------------------------------------------------
template <typename T> T maxAbsoluteDifference(const Hybrid<T> &va, const Hybrid<T> &vb) {
  vectorComparisonCheck(va, vb, "maxAbsoluteDifference");
  return maxAbsoluteDifference(va.data(), vb.data(), va.size());
}

//-------------------------------------------------------------------------------------------------
template <typename T> double maxRelativeDifference(const T* va, const T* vb, const size_t length) {
  double mdev = 0.0;
  double abs_mdev = 0.0;
  for (size_t i = 0; i < length; i++) {
    double dvai = static_cast<double>(va[i]);
    double dvbi = static_cast<double>(vb[i]);
    double tdev = (dvai - dvbi) / dvbi;
    double abs_tdev = fabs(tdev);
    if (abs_tdev > abs_mdev) {
      mdev = tdev;
      abs_mdev = abs_tdev;
    }
  }
  return mdev;
}

//-------------------------------------------------------------------------------------------------
template <typename T> double maxRelativeDifference(const std::vector<T> &va,
                                                   const std::vector<T> &vb) {

  // Compute the maximum relative deviation in double precision
  vectorComparisonCheck(va, vb, "maxRelativeDifference");
  return maxRelativeDifference(va.data(), vb.data(), va.size());
}

//-------------------------------------------------------------------------------------------------
template <typename T> double maxRelativeDifference(const Hybrid<T> &va, const Hybrid<T> &vb) {

  // Compute the maximum relative deviation in double precision
  vectorComparisonCheck(va, vb, "maxRelativeDifference");
  return maxRelativeDifference(va.data(), vb.data(), va.size());
}

//-------------------------------------------------------------------------------------------------
template <typename T> double meanUnsignedError(const T* va, const T* vb, const size_t length) {

  // Accumulate the mean unsigned error in double precision
  // after computing differences in the vectors' native type
  double mue = 0.0;
  if (isSignedIntegralScalarType<T>() || isFloatingPointScalarType<T>()) {
    for (size_t i = 0; i < length; i++) {
      T tdev = va[i] - vb[i];
      mue += fabs(static_cast<double>(tdev));
    }
  }
  else if (isUnsignedIntegralScalarType<T>()) {
    for (size_t i = 0; i < length; i++) {
      T tdev = (va[i] > vb[i]) ? va[i] - vb[i] : vb[i] - va[i];
      mue += static_cast<double>(tdev);
    }
  }
  else {
    rtErr("Data of type " + std::string(typeid(T).name()) + " is not suitable for comparisons.  "
          "Update isScalarType() and its related subclassification functions.",
          "meanUnsignedError");
  }
  return mue / static_cast<double>(length);
}

//-------------------------------------------------------------------------------------------------
template <typename T> double meanUnsignedError(const std::vector<T> &va,
                                               const std::vector<T> &vb) {
  vectorComparisonCheck(va, vb, "meanUnsignedError");
  return meanUnsignedError(va.data(), vb.data(), va.size());
}

//-------------------------------------------------------------------------------------------------
template <typename T> double meanUnsignedError(const Hybrid<T> &va, const Hybrid<T> &vb) {
  vectorComparisonCheck(va, vb, "meanUnsignedError");
  return meanUnsignedError(va.data(), vb.data(), va.size());
}

//-------------------------------------------------------------------------------------------------
template <typename T> double relativeRmsError(const T* va, const T* vb, const size_t length) {

  // Accumulate the root mean squared error and mean unsigned
  // value of the reference vector, all in double precision
  double rmse = 0.0;
  double mvalue = 0.0;
  for (size_t i = 0; i < length; i++) {
    double dvai = static_cast<double>(va[i]);
    double dvbi = static_cast<double>(vb[i]);
    double tdev = (dvai - dvbi);
    rmse += tdev * tdev;
    mvalue += fabs(dvbi);
  }
  rmse = sqrt(rmse / static_cast<double>(length));
  return rmse / mvalue;
}

//-------------------------------------------------------------------------------------------------
template <typename T> double relativeRmsError(const std::vector<T> &va,
                                              const std::vector<T> &vb) {
  vectorComparisonCheck(va, vb, "relativeRmsError");
  return relativeRmsError(va.data(), vb.data(), va.size());
}

//-------------------------------------------------------------------------------------------------
template <typename T> double relativeRmsError(const Hybrid<T> &va, const Hybrid<T> &vb) {
  vectorComparisonCheck(va, vb, "relativeRmsError");
  return relativeRmsError(va.data(), vb.data(), va.size());
}

//-------------------------------------------------------------------------------------------------
template <typename T> double mean(const T* va, const size_t length) {
  const double inv_length = 1.0 / static_cast<double>(length);
  double mvalue = 0.0;
  if (isScalarType<T>()) {
    for (size_t i = 0; i < length; i++) {
      mvalue += static_cast<double>(va[i]);
    }
  }
  else {
    rtErr("Data type " + std::string(typeid(T).name()) + " is not suitable for computing a mean "
          "of scalar values.", "mean");
  }
  return mvalue * inv_length;
}

//-------------------------------------------------------------------------------------------------
template <typename T> double mean(const std::vector<T> &va) {
  return mean(va.data(), va.size());
}

//-------------------------------------------------------------------------------------------------
template <typename T> double mean(const Hybrid<T> &va) {
  return mean(va.data(), va.size());
}

//-------------------------------------------------------------------------------------------------
template <typename T> double variance(const T* va, const size_t length,
                                      const VarianceMethod method) {
  if (isScalarType<T>() == false) {
    rtErr("Data type " + std::string(typeid(T).name()) + " is not suitable for computing the "
          "variance of scalar values.", "mean");
  }
  if (length == 0) {
    rtErr("No data was present for computation.", "variance");
  }
  double s1 = 0.0;
  double s2 = 0.0;
  double mvalue = 0.0;
  const double dlength = static_cast<double>(length);
  switch (method) {
  case VarianceMethod::VARIANCE:
  case VarianceMethod::STANDARD_DEVIATION:
  case VarianceMethod::ROOT_MEAN_SQUARED_DEVIATION:
    for (size_t i = 0; i < length; i++) {
      double tval = static_cast<double>(va[i]);
      s1 += tval;
      s2 += tval * tval;
    }
    break;
  case VarianceMethod::COEFFICIENT_OF_VARIATION:
  case VarianceMethod::NORMALIZED_RMSD:
    for (size_t i = 0; i < length; i++) {
      double tval = static_cast<double>(va[i]);
      s1 += tval;
      s2 += tval * tval;
      mvalue += tval;
    }
    mvalue /= dlength;
    break;
  }
  switch (method) {
  case VarianceMethod::VARIANCE:
    return ((dlength * s2) - (s1 * s1)) / (dlength * dlength);
  case VarianceMethod::STANDARD_DEVIATION:
    if (length == 1) {
      rtErr("Standard deviation is undefined for a single number.");
    }
    return sqrt(((dlength * s2) - (s1 * s1)) / (dlength * (dlength - 1.0)));
  case VarianceMethod::ROOT_MEAN_SQUARED_DEVIATION:
    return sqrt(((dlength * s2) - (s1 * s1))) / dlength;
  case VarianceMethod::COEFFICIENT_OF_VARIATION:
    return sqrt(((dlength * s2) - (s1 * s1)) / (dlength * (dlength - 1.0))) / fabs(mvalue);
  case VarianceMethod::NORMALIZED_RMSD:
    return sqrt(((dlength * s2) - (s1 * s1))) / (dlength * fabs(mvalue));
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename T> double variance(const std::vector<T> &va, const VarianceMethod method) {
  return variance(va.data(), va.size(), method);
}

//-------------------------------------------------------------------------------------------------
template <typename T> double variance(const Hybrid<T> &va, const VarianceMethod method) {
  return variance(va.data(), va.size(), method);
}

//-------------------------------------------------------------------------------------------------
template <typename T> double2 bivariateMean(const T* va, const size_t length) {
  if (isHpcVectorType<T>() && getHpcVectorTypeSize<T>() == 2) {
    const double inv_length = 1.0 / static_cast<double>(length);
    double2 mvalue = { 0.0, 0.0 };
    for (size_t i = 0; i < length; i++) {
      mvalue.x += static_cast<double>(va[i].x);
      mvalue.y += static_cast<double>(va[i].y);
    }
    mvalue.x *= inv_length;
    mvalue.y *= inv_length;
    return mvalue;
  }
  else {
    rtErr("Data type " + std::string(typeid(T).name()) + " is not a suitable 2-tuple.",
          "bivariateMean");
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T> double2 bivariateMean(const std::vector<T> &va) {
  return bivariateMean(va.data(), va.size());
}

//-------------------------------------------------------------------------------------------------
template <typename T> double2 bivariateMean(const Hybrid<T> &va) {
  return bivariateMean(va.data(), va.size());
}

//-------------------------------------------------------------------------------------------------
template <typename T> double3 trivariateMean(const T* va, const size_t length) {
  if (isHpcVectorType<T>() && getHpcVectorTypeSize<T>() == 3) {
    const double inv_length = 1.0 / static_cast<double>(length);
    double3 mvalue = { 0.0, 0.0, 0.0 };
    for (size_t i = 0; i < length; i++) {
      mvalue.x += static_cast<double>(va[i].x);
      mvalue.y += static_cast<double>(va[i].y);
      mvalue.z += static_cast<double>(va[i].z);
    }
    mvalue.x *= inv_length;
    mvalue.y *= inv_length;
    mvalue.z *= inv_length;
    return mvalue;
  }
  else {
    rtErr("Data type " + std::string(typeid(T).name()) + " is not a suitable 3-tuple.",
          "trivariateMean");
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T> double3 trivariateMean(const std::vector<T> &va) {
  return trivariateMean(va.data(), va.size());
}

//-------------------------------------------------------------------------------------------------
template <typename T> double3 trivariateMean(const Hybrid<T> &va) {
  return trivariateMean(va.data(), va.size());
}

//-------------------------------------------------------------------------------------------------
template <typename T> double4 quadrivariateMean(const T* va, const size_t length) {
  if (isHpcVectorType<T>() && getHpcVectorTypeSize<T>() == 4) {
    const double inv_length = 1.0 / static_cast<double>(length);
    double4 mvalue = { 0.0, 0.0, 0.0, 0.0 };
    for (size_t i = 0; i < length; i++) {
      mvalue.x += static_cast<double>(va[i].x);
      mvalue.y += static_cast<double>(va[i].y);
      mvalue.z += static_cast<double>(va[i].z);
      mvalue.w += static_cast<double>(va[i].w);
    }
    mvalue.x *= inv_length;
    mvalue.y *= inv_length;
    mvalue.z *= inv_length;
    mvalue.w *= inv_length;
    return mvalue;
  }
  else {
    rtErr("Data type " + std::string(typeid(T).name()) + " is not a suitable 4-tuple.",
          "quadrivariateMean");
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T> double4 quadrivariateMean(const std::vector<T> &va) {
  return quadrivariateMean(va.data(), va.size());
}

//-------------------------------------------------------------------------------------------------
template <typename T> double4 quadrivariateMean(const Hybrid<T> &va) {
  return quadrivariateMean(va.data(), va.size());
}

//-------------------------------------------------------------------------------------------------
template <typename T> T maxValue(const T* va, const size_t length) {
  if (isScalarType<T>() == false) {
    rtErr("Data type " + std::string(typeid(T).name()) + " is not suitable for computing a "
          "maximum value.", "maxValue");
  }
  if (length == 0) {
    rtErr("There is no data to scan for a maximum value.", "maxValue");
  }
  T max_value = va[0];
  for (size_t i = 1; i < length; i++) {
    max_value = std::max(max_value, va[i]);
  }
  return max_value;
}

//-------------------------------------------------------------------------------------------------
template <typename T> T maxValue(const std::vector<T> &va) {
  return maxValue(va.data(), va.size());
}

//-------------------------------------------------------------------------------------------------
template <typename T> T maxValue(const Hybrid<T> &va) {
  return maxValue(va.data(), va.size());
}

//-------------------------------------------------------------------------------------------------
template <typename T> T minValue(const T* va, const size_t length) {
  if (isScalarType<T>() == false) {
    rtErr("Data type " + std::string(typeid(T).name()) + " is not suitable for computing a "
          "minimum value.", "minValue");
  }
  if (length == 0) {
    rtErr("There is no data to scan for a minimum value.", "minValue");
  }
  T min_value = va[0];
  for (size_t i = 1; i < length; i++) {
    min_value = std::min(min_value, va[i]);
  }
  return min_value;
}

//-------------------------------------------------------------------------------------------------
template <typename T> T minValue(const std::vector<T> &va) {
  return minValue(va.data(), va.size());
}

//-------------------------------------------------------------------------------------------------
template <typename T> T minValue(const Hybrid<T> &va) {
  return minValue(va.data(), va.size());
}

//-------------------------------------------------------------------------------------------------
template <typename T> T maxAbsValue(const T* va, const size_t length) {
  if (isScalarType<T>() == false) {
    rtErr("Data type " + std::string(typeid(T).name()) + " is not suitable for computing a "
          "maximum absolute value.", "maxAbsValue");
  }
  if (length == 0) {
    rtErr("There is no data to scan for a minimum value.", "minValue");
  }
  T maxabs_value = (va[0] < static_cast<T>(0)) ? -va[0] : va[0];
  for (size_t i = 1; i < length; i++) {
    maxabs_value = std::max(maxabs_value, (va[i] < static_cast<T>(0)) ? -va[i] : va[i]);
  }
  return maxabs_value;
}

//-------------------------------------------------------------------------------------------------
template <typename T> T maxAbsValue(const std::vector<T> &va) {
  return maxAbsValue(va.data(), va.size());
}

//-------------------------------------------------------------------------------------------------
template <typename T> T maxAbsValue(const Hybrid<T> &va) {
  return maxAbsValue(va.data(), va.size());
}

//-------------------------------------------------------------------------------------------------
template <typename T> double pearson(const T* va, const T* vb, const size_t length) {
  if (length < 2) {
    rtErr("There is insufficient data to correlate.", "minValue");
  }
  double sum_a  = 0.0;
  double sum_b  = 0.0;
  double sum_ab = 0.0;
  double sum_aa = 0.0;
  double sum_bb = 0.0;
  for (size_t i = 0; i < length; i++) {
    const double value_a = va[i];
    const double value_b = vb[i];
    sum_a += value_a;
    sum_b += value_b;
    sum_ab += value_a * value_b;
    sum_aa += value_a * value_a;
    sum_bb += value_b * value_b;
  }
  const double dlength = static_cast<double>(length);
  const double sq = (dlength*sum_aa - sum_a*sum_a) * (dlength*sum_bb - sum_b*sum_b);
  if (sq < constants::tiny) {
    return 0.0;
  }
  return (dlength*sum_ab - sum_a*sum_b) / sqrt(sq);
}

//-------------------------------------------------------------------------------------------------
template <typename T> double pearson(const std::vector<T> &va, const std::vector<T> &vb) {
  vectorComparisonCheck(va, vb, "pearson");
  return pearson(va.data(), vb.data(), va.size());
}

//-------------------------------------------------------------------------------------------------
template <typename T> double pearson(const Hybrid<T> &va, const Hybrid<T> &vb) {
  vectorComparisonCheck(va, vb, "pearson");
  return pearson(va.data(), vb.data(), va.size());
}

//-------------------------------------------------------------------------------------------------
template <typename T> void addScalarToVector(T* va, const size_t length, const T inc) {
  for (size_t i = 0; i < length; i++) {
    va[i] += inc;
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T> void addScalarToVector(std::vector<T> *va, const T inc) {
  addScalarToVector(va->data(), va->size(), inc);
}

//-------------------------------------------------------------------------------------------------
template <typename T> void addScalarToVector(Hybrid<T> *va, const T inc) {
  addScalarToVector(va->data(), va->size(), inc);
}

//-------------------------------------------------------------------------------------------------
template <typename T> void elementwiseMultiply(T* va, const size_t length, const T factor) {
  for (size_t i = 0; i < length; i++) {
    va[i] *= factor;
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T> void elementwiseMultiply(std::vector<T> *va, const T factor) {
  elementwiseMultiply(va->data(), va->size(), factor);
}

//-------------------------------------------------------------------------------------------------
template <typename T> void elementwiseMultiply(Hybrid<T> *va, const T factor) {
  elementwiseMultiply(va->data(), va->size(), factor);
}

//-------------------------------------------------------------------------------------------------
template <typename T> void elementwiseDivide(T* va, const size_t length, const T factor) {
  for (size_t i = 0; i < length; i++) {
    va[i] /= factor;
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T> void elementwiseDivide(std::vector<T> *va, const T factor) {
  elementwiseDivide(va->data(), va->size(), factor);
}

//-------------------------------------------------------------------------------------------------
template <typename T> void elementwiseDivide(Hybrid<T> *va, const T factor) {
  elementwiseDivide(va->data(), va->size(), factor);
}

//-------------------------------------------------------------------------------------------------
template <typename T> void crossProduct(const T* va, const T* vb, T* vc) {
  vc[0] = (va[1] * vb[2]) - (va[2] * vb[1]);
  vc[1] = (va[2] * vb[0]) - (va[0] * vb[2]);
  vc[2] = (va[0] * vb[1]) - (va[1] * vb[0]);
}

//-------------------------------------------------------------------------------------------------
template <typename T> void crossProduct(const std::vector<T> &va, const std::vector<T> &vb,
                                        std::vector<T> *vc) {

  // Bounds check
  if (va.size() < 3 || vb.size() < 3) {
    rtErr("Vectors of size " + std::to_string(va.size()) + " and " + std::to_string(vb.size()) +
          "cannot produce a cross product.", "crossProduct");
  }
  else if (vc->size() < 3) {
    rtErr("A vector of size " + std::to_string(vc->size()) + " cannot hold a cross product.",
          "crossProduct");
  }

  // Set a pointer to the data of vector C and compute
  T* c_ptr = vc->data();
  c_ptr[0] = (va[1] * vb[2]) - (va[2] * vb[1]);
  c_ptr[1] = (va[2] * vb[0]) - (va[0] * vb[2]);
  c_ptr[2] = (va[0] * vb[1]) - (va[1] * vb[0]);
}

//-------------------------------------------------------------------------------------------------
template <typename T> T crossProduct(const T va, const T vb) {
  T vc;
  vc.x = (va.y * vb.z) - (va.z * vb.y);
  vc.y = (va.z * vb.x) - (va.x * vb.z);
  vc.z = (va.x * vb.y) - (va.y * vb.x);
  return vc;
}

//-------------------------------------------------------------------------------------------------
template <typename T> double magnitude(const T* va, const size_t length) {
  double result = 0.0;
  for (size_t i = 0; i < length; i++) {
    result += va[i] * va[i];
  }
  return sqrt(result);
}

//-------------------------------------------------------------------------------------------------
template <typename T> double magnitude(const std::vector<T> &va) {
  return magnitude(va.data(), va.size());
}

//-------------------------------------------------------------------------------------------------
template <typename T> double magnitude(const Hybrid<T> &va) {
  return magnitude(va.data(), va.size());
}

//-------------------------------------------------------------------------------------------------
template <typename T> double normalize(T* va, const size_t length) {
  const double magvec = magnitude(va, length);
  const double invmag = 1.0 / magvec;
  for (size_t i = 0; i < length; i++) {
    va[i] *= invmag;
  }
  return magvec;
}

//-------------------------------------------------------------------------------------------------
template <typename T> double normalize(std::vector<T> *va) {
  return normalize(va->data(), va->size());
}

//-------------------------------------------------------------------------------------------------
template <typename T> double normalize(Hybrid<T> *va) {
  return normalize(va->data(), va->size());
}

//-------------------------------------------------------------------------------------------------
template <typename T> double dot(const T* va, const T* vb, const size_t length) {
  double result = 0.0;
  for (size_t i = 0; i < length; i++) {
    result += va[i] * vb[i];
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename T> double dot(const std::vector<T> &va, const std::vector<T> &vb) {
  vectorComparisonCheck(va, vb, "dot");
  return dot(va.data(), vb.data(), va.size());
}

//-------------------------------------------------------------------------------------------------
template <typename T> double dot(const Hybrid<T> &va, const Hybrid<T> &vb) {
  vectorComparisonCheck(va, vb, "dot");
  return dot(va.data(), vb.data(), va.size());
}

//-------------------------------------------------------------------------------------------------
template <typename T> double angleBetweenVectors(const T* va, const T* vb, const size_t length) {

  // Do this in slightly more optimized fashion
  double dotp_acc = 0.0;
  double maga_acc = 0.0;
  double magb_acc = 0.0;
  for (size_t i = 0; i < length; i++) {
    const double vai = va[i];
    const double vbi = vb[i];
    dotp_acc += vai * vbi;
    maga_acc += vai * vai;
    magb_acc += vbi * vbi;
  }
  const double mag2_ab = maga_acc * magb_acc;
  if (fabs(mag2_ab) < constants::tiny) {
    rtErr("One or both vectors are of close to zero length.", "angleBetweenVectors");
  }
  const double acos_arg = dotp_acc / (sqrt(maga_acc) * sqrt(magb_acc));
  if (acos_arg >= 1.0) {
    return 0.0;
  }
  else if (acos_arg <= -1.0) {
    return symbols::pi;
  }
  else {
    return acos(acos_arg);
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T> double angleBetweenVectors(const std::vector<T> &va,
                                                 const std::vector<T> &vb) {
  if (va.size() != vb.size()) {
    rtErr("Vectors of differing sizes, " + std::to_string(va.size()) + " and " +
          std::to_string(vb.size()) + ", cannot produce an angle value.", "angleBetweenVectors");
  }
  return angleBetweenVectors(va.data(), vb.data(), va.size());
}

//-------------------------------------------------------------------------------------------------
template <typename T> double angleBetweenVectors(const Hybrid<T> &va, const Hybrid<T> &vb) {
  if (va.size() != vb.size()) {
    rtErr("Hybrid objects " + std::string(va.getLabel().name) + " and " +
          std::string(vb.getLabel().name) + " have differing sizes, " + std::to_string(va.size()) +
          " and " + std::to_string(vb.size()) + ", and therefore cannot produce an angle value.",
          "angleBetweenVectors");
  }
  return angleBetweenVectors(va.data(), vb.data(), va.size());
}

//-------------------------------------------------------------------------------------------------
template <typename T> void project(const T* va, const T* vb, T* vc, const size_t length) {
  const double mag_vb = magnitude(vb, length);
  const double dp_val = dot(va, vb, length) / (mag_vb * mag_vb);
  for (size_t i = 0; i < length; i++) {
    vc[i] = vb[i] * dp_val;
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T> void project(const std::vector<T> &va, const std::vector<T> &vb,
                                   std::vector<T> *vc) {
  const size_t va_len = va.size();
  if (vb.size() != va_len || vc->size() != va_len) {
    rtErr("Vectors of length " + std::to_string(va_len) + " (vA), " + std::to_string(vb.size()) +
          " (vB), and " + std::to_string(vc->size()) + " (vC) are invalid for computing the "
          "projection of vA onto vB and storing it as vC.", "project");
  }
  project(va.data(), vb.data(), vc->data(), va_len);
}
  
//-------------------------------------------------------------------------------------------------
template <typename T> void project(const Hybrid<T> &va, const Hybrid<T> &vb, Hybrid<T> *vc) {
  const size_t va_len = va.size();
  if (vb.size() != va_len || vc->size() != va_len) {
    rtErr("Vectors of length " + std::to_string(va_len) + " (vA), " + std::to_string(vb.size()) +
          " (vB), and " + std::to_string(vc->size()) + " (vC) are invalid for computing the "
          "projection of vA (" + std::string(va.getLabel().name) + ") onto vB (" +
          std::string(vb.getLabel().name) + ") and storing it as vC (" +
          std::string(vc->getLabel().name) + ").", "project");
  }
  project(va.data(), vb.data(), vc->data(), va_len);
}

//-------------------------------------------------------------------------------------------------
template <typename T> std::vector<T> project(const std::vector<T> &va, const std::vector<T> &vb) {
  vectorComparisonCheck(va, vb, "project");
  const size_t va_len = va.size();
  const double mag_vb = magnitude(vb.data(), va_len);
  const double dp_val = dot(va.data(), vb.data(), va_len) / (mag_vb * mag_vb);
  std::vector<T> result(va_len);
  for (size_t i = 0; i < va_len; i++) {
    result[i] = vb[i] * dp_val;
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename T> void perpendicularComponent(const T* va, const T* vb, T* result,
                                                  const size_t length) {
  project(va, vb, result, length);
  for (size_t i = 0; i < length; i++) {
    result[i] = va[i] - result[i];
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T> void perpendicularComponent(const std::vector<T> &va,
                                                  const std::vector<T> &vb,
                                                  std::vector<T> *result) {
  vectorComparisonCheck(va, vb, "perpendicularComponent");
  vectorComparisonCheck(va, result, "perpendicularComponent");
  perpendicularComponent(va.data(), vb.data(), result->data(), va.size());
}

//-------------------------------------------------------------------------------------------------
template <typename T> std::vector<T> perpendicularComponent(const std::vector<T> &va,
                                                            const std::vector<T> &vb) {
  vectorComparisonCheck(va, vb, "perpendicularComponent");
  std::vector<double> result = project(va, vb);
  const size_t va_dim = va.size();
  for (size_t i = 0; i < va_dim; i++) {
    result[i] = va[i] - result[i];
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename T> double pointPlaneDistance(const T* va, const T* vb, const T* pt_pos) {
  double unit_normal[3], pt_plane_displacement[3];
  crossProduct(va, vb, unit_normal);
  normalize(unit_normal, 3);
  project(unit_normal, pt_pos, pt_plane_displacement, 3);
  return magnitude(pt_plane_displacement, 3);
}

//-------------------------------------------------------------------------------------------------
template <typename T> double pointPlaneDistance(const std::vector<T> &va, const std::vector<T> &vb,
                                                const std::vector<T> &pt_pos) {
  if (va.size() != 3LLU || vb.size() != 3LLU || pt_pos.size() != 3LLU) {
    rtErr("The distance from a point to a plane is computed in three dimensions.  Dimensions of "
          "vectors provided = [ " + std::to_string(va.size()) + ", " + std::to_string(vb.size()) +
          ", " + std::to_string(pt_pos.size()) + " ].", "pointPlaneDistance");
  }
  return pointPlaneDistance(va.data(), vb.data(), pt_pos.data());
}

//-------------------------------------------------------------------------------------------------
template <typename T> double pointPlaneDistance(const Hybrid<T> &va, const Hybrid<T> &vb,
                                                const Hybrid<T> &pt_pos) {
  if (va.size() != 3LLU || vb.size() != 3LLU || pt_pos.size() != 3LLU) {
    rtErr("The distance from a point to a plane is computed in three dimensions.  Dimensions of "
          "Hybrid objects " + std::string(va.getLabel().name) + " and " +
          std::string(vb.getLabel().name) + " are [ " + std::to_string(va.size()) + ", " +
          std::to_string(vb.size()) + ", " + std::to_string(pt_pos.size()) + " ].",
          "pointPlaneDistance");
  }
  return pointPlaneDistance(va.data(), vb.data(), pt_pos.data());
}

//-------------------------------------------------------------------------------------------------
template <typename T> void reportBinLimitError(const std::string &desc, const T value,
                                               const ExceptionResponse policy) {
  switch (policy) {
  case ExceptionResponse::DIE:
    if (isSignedIntegralScalarType<T>()) {
      rtErr("A value of " + std::to_string(static_cast<llint>(value)) + 
            " is off the " + desc + " end of the range.", "findBin");
    }
    else if (isUnsignedIntegralScalarType<T>()) {
      rtErr("A value of " + std::to_string(static_cast<ullint>(value)) + 
            " is off the " + desc + " end of the range.", "findBin");
    }
    else if (isFloatingPointScalarType<T>()) {
      rtErr("A value of " +
            realToString(static_cast<double>(value), 11, 4, NumberFormat::SCIENTIFIC) +
            " is off the " + desc + " end of the range.", "findBin");
    }
    break;
  case ExceptionResponse::WARN:
    if (isSignedIntegralScalarType<T>()) {
      rtWarn("A value of " + std::to_string(static_cast<llint>(value)) + 
             " is off the " + desc + " end of the range.", "findBin");
    }
    else if (isUnsignedIntegralScalarType<T>()) {
      rtWarn("A value of " + std::to_string(static_cast<ullint>(value)) + 
             " is off the " + desc + " end of the range.", "findBin");
    }
    else if (isFloatingPointScalarType<T>()) {
      rtWarn("A value of " +
             realToString(static_cast<double>(value), 11, 4, NumberFormat::SCIENTIFIC) +
             " is off the " + desc + " end of the range.", "findBin");
    }
    break;
  case ExceptionResponse::SILENT:
    break;
  }
}
  
//-------------------------------------------------------------------------------------------------
template <typename T> int findBin(const T* limits, const T value, const int length,
                                  const ExceptionResponse policy) {
  if (length == 0) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("A zero-length set of bin limits was supplied.", "findBin");
    case ExceptionResponse::WARN:
      rtWarn("A zero-length set of bin limits was supplied.", "findBin");
      return -1;
    case ExceptionResponse::SILENT:
      return -1;
    }
  }
  if (value < limits[0]) {
    reportBinLimitError("left", value, policy);
    return -1;
  }
  if (value >= limits[length]) {
    reportBinLimitError("right", value, policy);
    return length + 1;
  }
  int lguess = 0;
  int hguess = length;
  while (lguess < hguess - 1LLU) {

    // Choose a residue intermediate between the lower and upper bounds
    const int mguess = lguess + ((hguess - lguess) / 2);
    if (value >= limits[mguess + 1LLU]) {
      lguess = mguess;
    }
    else if (value < limits[mguess]) {
      hguess = mguess;
    }
    else {
      return mguess;
    }
  }
  return lguess;
}

//-------------------------------------------------------------------------------------------------
template <typename T> int findBin(const std::vector<T> &limits, const T value,
                                  const ExceptionResponse policy) {
  return findBin(limits.data(), value, limits.size(), policy);
}

//-------------------------------------------------------------------------------------------------
template <typename T> int findBin(const Hybrid<T> &limits, const T value,
                                  const ExceptionResponse policy) {
  return findBin(limits.data(), value, limits.size(), policy);
}

//-------------------------------------------------------------------------------------------------
template <typename T> size_t locateValue(const T* vdata, const T value, const size_t length,
                                         const DataOrder format) {
  if (length == 0LLU) {
    rtErr("Unable to search an array of zero length.", "locateValue");
  }
  switch (format) {
  case DataOrder::ASCENDING:
    if (value >= vdata[0] && value <= vdata[length - 1]) {
      size_t min_pos = 0LLU;
      size_t max_pos = length;
      while (max_pos - min_pos > 1) {
        size_t mid_pos = (max_pos + min_pos) / 2;
        if (vdata[mid_pos] > value) {
          max_pos = mid_pos;
        }
        else {
          min_pos = mid_pos;
        }
      }
      return (vdata[min_pos] == value) ? min_pos : length;
    }
    break;
  case DataOrder::DESCENDING:
    if (value <= vdata[0] && value >= vdata[length - 1]) {
      size_t min_pos = 0LLU;
      size_t max_pos = length;
      while (max_pos - min_pos > 1) {
        size_t mid_pos = (max_pos + min_pos) / 2;
        if (vdata[mid_pos] < value) {
          max_pos = mid_pos;
        }
        else {
          min_pos = mid_pos;
        }
      }
      return (vdata[min_pos] == value) ? min_pos : length;
    }
    break;
  case DataOrder::NONE:
    for (size_t i = 0; i < length; i++) {
      if (vdata[i] == value) {
        return i;
      }
    }
    break;
  }
  return length;
}

//-------------------------------------------------------------------------------------------------
template <typename T> size_t locateValue(const std::vector<T> &vdata, const T value,
                                         const DataOrder format) {
  return locateValue(vdata.data(), value, vdata.size(), format);
}
  
//-------------------------------------------------------------------------------------------------
template <typename T> size_t locateValue(const Hybrid<T> &vdata, const T value,
                                         const DataOrder format) {
  return locateValue(vdata.data(), value, vdata.size(), format);
}

//-------------------------------------------------------------------------------------------------
template <typename T> size_t locateValue(const Approx &value, const T* vdata, const size_t length,
                                         const DataOrder format) {
  if (length == 0LLU) {
    rtErr("Unable to search an array of zero length.", "locateValue");
  }
  switch (format) {
  case DataOrder::ASCENDING:
    if (value >= vdata[0] && value <= vdata[length - 1]) {
      size_t min_pos = 0LLU;
      size_t max_pos = length;
      while (max_pos - min_pos > 1) {
        size_t mid_pos = (max_pos + min_pos) / 2;
        if (vdata[mid_pos] > value) {
          max_pos = mid_pos;
        }
        else {
          min_pos = mid_pos;
        }
      }
      return (vdata[min_pos] == value) ? min_pos : length;
    }
    break;
  case DataOrder::DESCENDING:
    if (value <= vdata[0] && value >= vdata[length - 1]) {
      size_t min_pos = 0LLU;
      size_t max_pos = length;
      while (max_pos - min_pos > 1) {
        size_t mid_pos = (max_pos + min_pos) / 2;
        if (vdata[mid_pos] < value) {
          max_pos = mid_pos;
        }
        else {
          min_pos = mid_pos;
        }
      }
      return (vdata[min_pos] == value) ? min_pos : length;
    }
    break;
  case DataOrder::NONE:
    for (size_t i = 0; i < length; i++) {
      if (vdata[i] == value) {
        return i;
      }
    }
    break;
  }
  return length;
}
  
//-------------------------------------------------------------------------------------------------
template <typename T> size_t locateValue(const Approx &value, const std::vector<T> &vdata,
                                         const DataOrder format) {
  return locateValue(value, vdata.data(), vdata.size(), format);
}
  
//-------------------------------------------------------------------------------------------------
template <typename T> size_t locateValue(const Approx &value, const Hybrid<T> &vdata,
                                         const DataOrder format) {
  return locateValue(value, vdata.data(), vdata.size(), format);
}

//-------------------------------------------------------------------------------------------------
template <typename T> std::vector<T> reduceUniqueValues(const std::vector<T> &va) {
  std::vector<T> result(va);
  reduceUniqueValues(&result);
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename T> void reduceUniqueValues(std::vector<T> *va) {
  std::sort(va->begin(), va->end(), [](T a, T b) { return a < b; });
  size_t nunique = 0LLU;
  const size_t nval = va->size();
  T* va_ptr = va->data();
  const size_t nvalm1 = nval - 1LLU;
  T last_unique;
  for (size_t i = 0LLU; i < nval; i++) {
    if (i == 0LLU || va_ptr[i] != last_unique) {
      last_unique = va_ptr[i];
      va_ptr[nunique] = last_unique;
      nunique++;
    }
  }
  va->resize(nunique);
  va->shrink_to_fit();
}

//-------------------------------------------------------------------------------------------------
template <typename T> std::vector<ValueWithCounter<T>>
findUnmatchedValues(const std::vector<T> &va, const std::vector<T> &vb,
                    const UniqueValueHandling check_repeats) {
  const size_t length_a = va.size();
  const size_t length_b = vb.size();
  std::vector<T> va_duplicate(va);
  std::vector<T> vb_duplicate(vb);
  std::sort(va_duplicate.begin(), va_duplicate.end(), [](T a, T b) { return a < b; });
  std::sort(vb_duplicate.begin(), vb_duplicate.end(), [](T a, T b) { return a < b; });
  std::vector<ValueWithCounter<T>> unique_a;
  std::vector<ValueWithCounter<T>> unique_b;
  for (size_t i = 0; i < length_a; i++) {
    if (i == 0 || va_duplicate[i - 1] != va_duplicate[i]) {
      unique_a.push_back({va_duplicate[i], 1});
    }
    else {
      unique_a[unique_a.size() - 1LLU].count += 1;
    }
  }
  for (size_t i = 0; i < length_b; i++) {
    if (i == 0 || vb_duplicate[i - 1] != vb_duplicate[i]) {
      unique_b.push_back({vb_duplicate[i], 1});
    }
    else {
      unique_b[unique_b.size() - 1LLU].count += 1;
    }
  }
  const size_t nua = unique_a.size();
  const size_t nub = unique_b.size();
  size_t min_j = 0;
  std::vector<ValueWithCounter<T>> result;
  const bool counts_matter = (check_repeats == UniqueValueHandling::CONFIRM_ALL_COPIES);
  for (size_t i = 0; i < nua; i++) {
    const T aval = unique_a[i].value;
    const int acount = unique_a[i].count;
    bool found = false;
    for (size_t j = min_j; j < nub; j++) {
      if (unique_b[j].value == aval) {
        min_j += (j == min_j);
        if (counts_matter && unique_b[j].count != acount) {
          result.push_back({aval, acount - unique_b[j].count});
        }
        found = true;
      }
    }
    if (! found) {
      result.push_back({aval, acount});
    }
  }
  min_j = 0;
  for (size_t i = 0; i < nub; i++) {
    const T bval = unique_b[i].value;
    const int bcount = unique_b[i].count;
    bool found = false;
    for (size_t j = min_j; j < nua; j++) {
      if (unique_a[j].value == bval) {
        min_j += (j == min_j);
        found = true;
      }
    }
    if (! found) {
      result.push_back({bval, -bcount});
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename T> std::vector<ValueWithCounter<T>>
findUnmatchedValues(const T* va, const T* vb, const size_t length_a, const size_t length_b,
                    const UniqueValueHandling check_repeats) {
  std::vector<T> tva(length_a);
  std::vector<T> tvb(length_b);
  return findUnmatchedValues(tva, tvb, check_repeats);
}

//-------------------------------------------------------------------------------------------------
template <typename T> std::vector<ValueWithCounter<T>>
findUnmatchedValues(const Hybrid<T> &va, const Hybrid<T> &vb,
                    const UniqueValueHandling check_repeats) {
  std::vector<T> tva = va.readHost();
  std::vector<T> tvb = vb.readHost();
  return findUnmatchedValues(tva, tvb, check_repeats);
}

} // namespace math
} // namespace omni
