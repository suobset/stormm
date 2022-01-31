// -*-c++-*-
#ifndef OMNI_SORTING_H
#define OMNI_SORTING_H

#include "DataTypes/omni_vector_types.h"

namespace omni {
namespace math {

/// \brief Enumerate the directions that a sort can go.
enum class SortDirection {
  ASCENDING,  ///< Values will be ordered in ascending order
  DESCENDING, ///< Values will be ordered in descending order
  AUTOMATIC   ///< Triggers locally defined "default" behavior on a case-by-case basis
};

} // namespace math
} // namespace omni

#endif
