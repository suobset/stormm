// -*-c++-*-
#ifndef EMULATE_ENUMERATORS_H
#define EMULATE_ENUMERATORS_H

#include <string>
#include "../../../src/copyright.h"

namespace emulation {

// Different rows of the fitting matrix may include energy or force contributions.  This enumerator
// will list the various options.
enum class FittingContribution {
  ENERGY,  // The energy or relative energy of a group of atoms is the objective
  FORCE    // The force (e.g. zero force) exerted on atoms selected by a mask is of interest
};

// Return a human-readable string based on an enumerator's value.  Various overloads cover each
// enum class in the namespace, as is the policy in the libraries.
//
// Arguments:
//   input:    The enumerator of interest
std::string getEnumerationName(FittingContribution input);
  
} // namespace emulation

#endif
