#include <cmath>
#include "copyright.h"
#include "formulas.h"

namespace stormm {
namespace stmath {

//-------------------------------------------------------------------------------------------------
int factorial(const int x) {
  if (x > 12) {
    rtErr("Factorials of integers greater than 12 would overflow a 32-bit integer format (signed "
          "or unsigned) and must be computed using the incomplete gamma function of the argument, "
          "plus one.", "factorial");
  }
  if (x < 0) {
    rtErr("Factorials of negative numbers are undefined.", "factorial");
  }
  int result = 1;
  for (int i = 1; i <= x; i++) {
    result *= i;
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
llint factorial(const llint x) {
  if (x > 20) {
    rtErr("Factorials of integers greater than 20 would overflow a 64-bit integer format (signed "
          "or unsigned) and must be computed using the incomplete gamma function of the argument, "
          "plus one.", "factoriall");
  }
  if (x < 0) {
    rtErr("Factorials of negative numbers are undefined.", "factoriall");
  }
  llint result = 1LL;
  for (llint i = 1; i <= x; i++) {
    result *= i;
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
double factoriald(const int x) {
  if (x > 170) {
    rtErr("Factorials of integers greater than 170 would overflow a 64-bit floating point format "
          "and cannot be computed in practical machine arithmetic.", "factoriald");
  }
  return tgamma(x + 1);
}

//-------------------------------------------------------------------------------------------------
int ipow(const int x, const int p) {
  if (p < 0) {
    rtErr("Negative powers are not allow in integer exponentiation.", "ipow");
  }
  int result = 1;
  int pwork = p;
  int base = x;
  while (pwork) {
    if (pwork & 1) {
      result *= base;
    }
    pwork >>= 1;
    if (pwork) {
      base *= base;
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
llint ipowl(const llint x, const int p) {
  if (p < 0) {
    rtErr("Negative powers are not allow in integer exponentiation.", "ipow");
  }
  llint result = 1;
  int pwork = p;
  llint base = x;
  while (pwork) {
    if (pwork & 1LL) {
      result *= base;
    }
    pwork >>= 1;
    if (pwork) {
      base *= base;
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
double4 sigmoid(const double r, const double crossover, const double intensity) {
  const double g = exp(intensity * (r - crossover));
  double4 result;
  result.x = 1.0 / (g + 1.0);
  result.y = -intensity * g / ((g + 1.0) * (g + 1.0));
  result.z = -intensity * intensity * g * (1.0 - g) / pow(g + 1.0, 3.0);
  result.w = -pow(intensity, 3.0) * g * ((g * g) - (4.0 * g) + 1.0) / pow(g + 1.0, 4.0);
  return result;
}

//-------------------------------------------------------------------------------------------------
double sigmoid(const double r, const double crossover, const double intensity, const int order) {
  const double g = exp(intensity * (r - crossover));
  switch (order) {
  case 0:
    return 1.0 / (g + 1.0);
  case 1:
    return -intensity * g / ((g + 1.0) * (g + 1.0));
  case 2:
    return -intensity * intensity * g * (1.0 - g) / pow(g + 1.0, 3.0);
  case 3:
    return -pow(intensity, 3.0) * g * ((g * g) - (4.0 * g) + 1.0) / pow(g + 1.0, 4.0);
  default:
    rtErr("Provide 0 for the function value or an integer 1-3 for the first, second, or third "
          "derivative.");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
float4 sigmoidf(const float r, const float crossover, const float intensity) {
  const float g = expf(intensity * (r - crossover));
  float4 result;
  result.x = 1.0f / (g + 1.0f);
  result.y = -intensity * g / ((g + 1.0f) * (g + 1.0f));
  result.z = -intensity * intensity * g * (1.0f - g) / powf(g + 1.0f, 3.0f);
  result.w = -powf(intensity, 3.0f) * g * ((g * g) - (4.0f * g) + 1.0f) / powf(g + 1.0f, 4.0f);
  return result;
}

//-------------------------------------------------------------------------------------------------
float sigmoidf(const float r, const float crossover, const float intensity, const int order) {
  const float g = expf(intensity * (r - crossover));
  switch (order) {
  case 0:
    return 1.0f / (g + 1.0f);
  case 1:
    return -intensity * g / ((g + 1.0f) * (g + 1.0f));
  case 2:
    return -intensity * intensity * g * (1.0f - g) / powf(g + 1.0f, 3.0f);
  case 3:
    return -powf(intensity, 3.0f) * g * ((g * g) - (4.0f * g) + 1.0f) / powf(g + 1.0f, 4.0f);
  default:
    rtErr("Provide 0 for the function value or an integer 1-3 for the first, second, or third "
          "derivative.");
  }
  __builtin_unreachable();
}

} // namespace stmath
} // namespace stormm
