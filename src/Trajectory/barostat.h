// -*-c++-*-
#ifndef OMNI_BAROSTAT_H
#define OMNI_BAROSTAT_H

#include <string>

namespace omni {
namespace trajectory {

/// \brief Enumerate the available barostat types
enum class BarostatKind {
  NONE, MONTE_CARLO, BERENDSEN
};

/// \brief Store the parameters for a simulation barostat.  Includes Monte-Carlo and Berendesen
///        barostats with atomic virials.
struct Barostat {
  BarostatKind kind;
};

/// \brief Return the name of the barostat choice (another enumerator string conversion function)
///
/// \param kind  The type of barostat
std::string getBarostatName(BarostatKind kind);

} // namespace trajectory
} // namespace omni

#endif
