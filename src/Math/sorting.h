// -*-c++-*-
#ifndef STORMM_SORTING_H
#define STORMM_SORTING_H

#include "DataTypes/stormm_vector_types.h"

namespace stormm {
namespace math {

/// \brief Enumerate the directions that a sort can go.
enum class SortDirection {
  ASCENDING,  ///< Values will be ordered in ascending order
  DESCENDING, ///< Values will be ordered in descending order
  AUTOMATIC   ///< Triggers locally defined "default" behavior on a case-by-case basis
};

/// \brief Enumerate the ways to handle multiple copies of a single value occurring in an array
enum class UniqueValueHandling {
  UNIQUE_VALUES_ONLY,  ///< Unique values will register once, regardless of the number of copies
  CONFIRM_ALL_COPIES   ///< Each copy of any given value will be confirmed, i.e. if two vectors
                       ///<   contain four entries of 5.5 and one contains two entries of 2.7
                       ///<   while the other contains three entries of 2.7, they differ by a
                       ///<   single entry of 2.7.
};

} // namespace math
} // namespace stormm

#endif
