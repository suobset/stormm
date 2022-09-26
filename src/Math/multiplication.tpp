// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace math {

//-------------------------------------------------------------------------------------------------
template <typename T> double logProduct(const T* values, const size_t length) {
  double esum = 0.0;
  for (size_t i = 0LLU; i < length; i++) {
    if (values[i] <= (T)(0)) {
      rtErr("The logarithm of a negative number, or zero, is undefined.  " +
            realToString(values[i], 11, 4, NumberFormat::STANDARD_REAL) + " was encountered in "
            "position " + std::to_string(i) + ".", "logProduct");
    }
    esum += log(values[i]);
  }
  return esum;
}

//-------------------------------------------------------------------------------------------------
template <typename T> double logProduct(const std::vector<T> &values) {
  return logProduct(values.data(), values.size());
}

//-------------------------------------------------------------------------------------------------
template <typename T> double logProduct(const Hybrid<T> &values) {
  return logProduct(values.data(), values.size());
}

//-------------------------------------------------------------------------------------------------
template <typename T> T seriesProduct(const T* va, const size_t length) {
  if (isScalarType<T>() == false) {
    rtErr("Data type " + std::string(typeid(T).name()) + " is not suitable for computing the "
          "product series of scalar values.", "productSeries");
  }
  if (isSignedIntegralScalarType<T>() || isUnsignedIntegralScalarType<T>()) {
    T result = 1;
    for (size_t i = 0LLU; i < length; i++) {
      result *= va[i];
    }
    return result;
  }
  else if (isFloatingPointScalarType<T>()) {

    // Compute the logarithm of the product to ensure that it stays within the limits of the
    // hardware.
    double log_p = 0.0;
    const double base_e_limit = (sizeof(T) == 4) ? log(2.0) * 128.0 : log(2.0) * 1024.0;
    bool violation = false;
    double sign_mult = 1.0;
    for (size_t i = 0LLU; i < length; i++) {
      if (va[i] == 0.0) {
        return 0.0;
      }
      else if (va[i] < 0.0) {
        log_p += log(-va[i]);
        sign_mult *= -1.0;
      }
      else {
        log_p += log(va[i]);
      }
      violation = (violation || (fabs(log_p) >= base_e_limit));
    }
    if (violation) {
      return sign_mult * exp(log_p);
    }
    else {
      T result = static_cast<T>(1.0);
      for (size_t i = 0LLU; i < length; i++) {
        result *= va[i];
      }
      return result;
    }
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename T> T seriesProduct(const std::vector<T> &va) {
  return seriesProduct(va.data(), va.size());
}

//-------------------------------------------------------------------------------------------------
template <typename T> T seriesProduct(const Hybrid<T> &va) {
  return seriesProduct(va.data(), va.size());
}

} // namespace math
} // namespace stormm
