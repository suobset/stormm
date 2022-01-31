#include <cmath>
#include "Math/vector_ops.h"
#include "Parsing/parse.h"
#include "vector_report.h"
#include "approx.h"

namespace omni {
namespace testing {

using math::maxValue;
using math::mean;
using math::meanUnsignedError;
using math::pearson;
using math::variance;
using math::VarianceMethod;
using parse::realDecimalPlaces;
using parse::realToString;

//-------------------------------------------------------------------------------------------------
std::string vectorAlignmentReport(const std::vector<PolyNumeric> &va,
                                  const std::vector<PolyNumeric> &vb, NumberFormat data_format,
                                  const double tol) {

  // Initialize the output
  std::string result;

  // Different numerical formats (plus char4)
  switch (data_format) {
  case NumberFormat::SCIENTIFIC:
  case NumberFormat::STANDARD_REAL:
  case NumberFormat::INTEGER:
  case NumberFormat::LONG_LONG_INTEGER:
  case NumberFormat::UNSIGNED_INTEGER:
  case NumberFormat::UNSIGNED_LONG_LONG_INTEGER:
    {
      // Convert both vectors to double precision real for analysis
      const std::vector<double> dva = doubleFromPolyNumeric(va);
      const std::vector<double> dvb = doubleFromPolyNumeric(vb);
      std::vector<double> work(std::max(va.size(), vb.size()), 0.0);
      const int n_va = va.size();
      const int n_vb = vb.size();

      // Look for a common multiple if the vectors are of the same size
      if (n_va == n_vb) {
        for (int i = 0; i < n_va; i++) {
          work[i] = dva[i] / dvb[i];
        }
        if (variance(work, VarianceMethod::NORMALIZED_RMSD) < tol) {
          result += "The vectors appear to be multiples of one another: first vector = " +
                    realToString(mean(work)) + " x second vector.";
        }
        else {
          int n_match = 0;
          int n_mismatch = 0;
          std::vector<double> vec_differences;
          for (int i = 0; i < n_va; i++) {
            if (dva[i] == Approx(dvb[i], ComparisonType::ABSOLUTE, tol)) {
              n_match++;
            }
            else {
              n_mismatch++;
              vec_differences.push_back(fabs(dva[i] - dvb[i]));
            }
          }
          const NumberFormat scifm = NumberFormat::SCIENTIFIC;
          if (n_match > 0) {
            result += "The vectors fail to match in " + std::to_string(n_mismatch) +
                      " indices out of " + std::to_string(n_va) + ", with a mean unsigned error "
                      "(among the deviating entries) of " +
                      realToString(mean(vec_differences), 11, 4, scifm) +
                      " and maximum unsigned deviation " +
                      realToString(maxValue(vec_differences), 11, 4, scifm) + ".  ";
          }
          else {
            result += "Deviations occur throughout the data.  ";
          }
          result += "Pearson correlation between the vectors is " +
                    realToString(pearson(dva, dvb), 4) + ".";
          if (n_mismatch < 5) {
            result += "Mismatched entries:\n";
            const int ndec = realDecimalPlaces(tol);
            for (int i = 0; i < n_va; i++) {
              if (dva[i] != Approx(dvb, ComparisonType::ABSOLUTE, tol)) {
                result += "    " + realToString(dva[i], ndec + 7, ndec, scifm) +
                          " != " + realToString(dvb[i], ndec + 7, ndec, scifm) +
                          " (error " + realToString(fabs(dvb[i] - dva[i]), ndec + 7, ndec, scifm) +
                          ", entry " + std::to_string(i) + ")\n";
              }
            }
          }
	}
      }
    }
    break;
  case NumberFormat::CHAR4:
    rtErr("Vector analysis method not yet implemented.", "vectorAlignmentReport");
    break;
  }
  return result;
}

} // namespace testing
} // namespace omni
