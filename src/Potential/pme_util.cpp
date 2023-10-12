#include "copyright.h"
#include "Constants/scaling.h"
#include "Parsing/parse.h"
#include "Parsing/parsing_enumerators.h"
#include "Reporting/error_format.h"
#include "pme_util.h"

namespace stormm {
namespace energy {

using parse::NumberFormat;
using parse::realToString;
  
//-------------------------------------------------------------------------------------------------
double ewaldCoefficient(const double cutoff, const double direct_sum_tol) {
  if (cutoff < constants::tiny) {
    rtErr("A cutoff of " + realToString(cutoff, 9, 4, NumberFormat::STANDARD_REAL) +
          " is invalid.", "ewaldCoefficient");
  }
  const double effective_dsum_tol = (direct_sum_tol < 0.0) ? 0.0 : direct_sum_tol;
  double result = 1.0;
  double best_error = fabs((erfc(result * cutoff) / cutoff) - effective_dsum_tol);
  double amin = 0.0;
  double amax = 10.0;
  for (int i = 0; i < 12; i++) {
    const double delta_a = 0.005 * (amax - amin);
    for (int j = 0; j < 200; j++) {
      const double a_update = amin + (delta_a * static_cast<double>(j));
      const double test_error = fabs((erfc(a_update * cutoff) / cutoff) - effective_dsum_tol);
      if (test_error < best_error) {
        best_error = test_error;
        result = a_update;
      }
    }
    amin = result - (5.0 * delta_a);
    amax = result + (5.0 * delta_a);
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
double pmeGaussianSpread(const double cutoff, const double direct_sum_tol) {
  return 0.5 / ewaldCoefficient(cutoff, direct_sum_tol);
}

} // namespace energy
} // namespace stormm
