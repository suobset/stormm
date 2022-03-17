#include <cmath>
#include "approx.h"

namespace omni {
namespace testing {

//-------------------------------------------------------------------------------------------------
Approx::Approx(const double value_in, ComparisonType style_in, const double tol_in) :
    Approx(std::vector<double>(1, value_in), style_in, tol_in)
{}

//-------------------------------------------------------------------------------------------------
Approx::Approx(const double value_in, const double tol_in, ComparisonType style_in) :
    Approx(std::vector<double>(1, value_in), style_in, tol_in)
{}

//-------------------------------------------------------------------------------------------------
int Approx::size() const {
  return values.size();
}

//-------------------------------------------------------------------------------------------------
double Approx::getValue() const {
  if (values.size() == 1) {
    return values[0];
  }
  else {
    rtErr("A single value was requested for an approximate comparison containing a vector of " +
          std::to_string(values.size()) + " values.  Use getValues() to retrieve them all.");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::vector<double> Approx::getValues() const {
  return values;
}

//-------------------------------------------------------------------------------------------------
ComparisonType Approx::getStyle() const {
  return style;
}

//-------------------------------------------------------------------------------------------------
double Approx::getMargin() const {
  return dtol;
}

//-------------------------------------------------------------------------------------------------
double Approx::getTolerance() const {
  return dtol;
}

//-------------------------------------------------------------------------------------------------
double Approx::getTol() const {
  return dtol;
}

//-------------------------------------------------------------------------------------------------
Approx Approx::margin(const double dtol_in) {
  return Approx(values, style, dtol_in);
}

//-------------------------------------------------------------------------------------------------
Approx Approx::tolerance(const double dtol_in) {
  return Approx(values, style, dtol_in);
}

//-------------------------------------------------------------------------------------------------
Approx Approx::tol(const double dtol_in) {
  return Approx(values, style, dtol_in);
}

//-------------------------------------------------------------------------------------------------
bool Approx::test(const double test_value) const {
  if (values.size() != 1) {
    return false;
  }
  switch (style) {
  case ComparisonType::ABSOLUTE:
  case ComparisonType::MEAN_UNSIGNED_ERROR:
    return (std::abs(test_value - values[0]) <= dtol);
  case ComparisonType::RELATIVE:
  case ComparisonType::RELATIVE_RMS_ERROR:
    if (std::abs(values[0]) > constants::tiny) {
      return (std::abs((test_value - values[0]) / values[0]) <= dtol);
    }
    else {
      return (std::abs((test_value - values[0]) / constants::tiny) <= dtol);
    }
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
bool operator==(const double d, const Approx &cr) {
  return cr.test(d);
}

//-------------------------------------------------------------------------------------------------
bool operator==(const Approx &cr, const double d) {
  return cr.test(d);
}

//-------------------------------------------------------------------------------------------------
bool operator!=(const double d, const Approx &cr) {
  return (cr.test(d) == false);
}

//-------------------------------------------------------------------------------------------------
bool operator!=(const Approx &cr, const double d) {
  return (cr.test(d) == false);
}

//-------------------------------------------------------------------------------------------------
bool operator>(const double d, const Approx &cr) {
  return (d > cr.getValue() + cr.getMargin());
}

//-------------------------------------------------------------------------------------------------
bool operator>(const Approx &cr, const double d) {
  return (cr.getValue() - cr.getMargin() > d);
}

//-------------------------------------------------------------------------------------------------
bool operator<(const double d, const Approx &cr) {
  return (d < cr.getValue() - cr.getMargin());
}

//-------------------------------------------------------------------------------------------------
bool operator<(const Approx &cr, const double d) {
  return (cr.getValue() + cr.getMargin() < d);
}

//-------------------------------------------------------------------------------------------------
bool operator>=(const double d, const Approx &cr) {
  return (d >= cr.getValue() - cr.getMargin());
}

//-------------------------------------------------------------------------------------------------
bool operator>=(const Approx &cr, const double d) {
  return (cr.getValue() + cr.getMargin() >= d);
}

//-------------------------------------------------------------------------------------------------
bool operator<=(const double d, const Approx &cr) {
  return (d <= cr.getValue() + cr.getMargin());
}

//-------------------------------------------------------------------------------------------------
bool operator<=(const Approx &cr, const double d) {
  return (cr.getValue() - cr.getMargin() <= d);
}

} // namespace testing
} // namespace omni
