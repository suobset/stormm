// -*-c++-*-
#ifndef STORMM_VECTOR_REPORT_H
#define STORMM_VECTOR_REPORT_H

#include <vector>
#include "copyright.h"
#include "Parsing/polynumeric.h"
#include "Reporting/error_format.h"

namespace stormm {
namespace testing {

using errors::listSeparator;
using parse::NumberFormat;
using parse::PolyNumeric;

/// \brief Report on the alignment of two vectors.  This should be called once the vectors have
///        been found to be different by some criterion, to check whether shifting them some
///        number of indices relative to one another would lead to a better match.
///
/// \param va           The first vector
/// \param vb           The second vector
/// \param data_format  The numericla format of each vector
/// \param tol          Comparison tolerance for the two vectors
std::string vectorAlignmentReport(const std::vector<PolyNumeric> &va,
                                  const std::vector<PolyNumeric> &vb, NumberFormat data_format,
                                  double tol);

/// \brief Create a string from a limited number of items from a list.  The items will be separated
///        by commas, using the Oxford comma convention.
///
///
/// Overloaded:
///   - Provide a list of numbers
///   - Provide a list of strings (other data types that cannot be interpreted by std::to_string()
///     should also receive their own special overloads)
///
/// \param list           The list of items (e.g. failing case indices)
/// \param maximum_items  The maximum number of items to present in the string.  If less than or
///                       equal to zero, all items of the list will be presented.
/// \param spread         Indicate whether to select items from throughout the list.  If FALSE,
///                       items will be pulled from the beginning of the list, up to the maximum
///                       number.
/// \{
template <typename T>
std::string listItemsAsString(const std::vector<T> &list, size_t maximum_items = 0,
                              bool spread = false);

std::string listItemsAsString(const std::vector<std::string> &list, size_t maximum_items = 0,
                              bool spread = false);
/// \}

} // namespace testing
} // namespace stormm

#include "vector_report.tpp"

#endif
