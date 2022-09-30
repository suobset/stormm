// -*-c++-*-
#ifndef STORMM_BEHAVIOR_H
#define STORMM_BEHAVIOR_H

#include <string>
#include "copyright.h"

namespace stormm {
namespace constants {

/// \brief Enumerate the ways in which the code can respond to an exception under a variety of
///        circumstances: die, warn (and continue), or stay silent and continue
enum class ExceptionResponse {
  DIE, WARN, SILENT
};

/// \brief Express choices for the case-sensitivity of various inputs
enum class CaseSensitivity {
  YES,      ///< Case matters
  NO,       ///< Evaluate without case sensitivity--convert everything to uppercase before parsing
  AUTOMATIC ///< Do not impose one style or another--defer to local, default behavior for
            ///<   individual inputs
};

/// \brief Enumerate precision models
enum class PrecisionModel {
  SINGLE,  ///< Evaluates most interactions in 32-bit integer and floating point arithmetic.
           ///< References single-precision data arrays if available.
  DOUBLE,  ///< Evaluates most interactions in 64-bit integer and floating point arithmetic.
           ///< References double-precision data arrays if available.
};

/// \brief Enumerate the Cartesian dimensions
enum class CartesianDimension {
  X, Y, Z
};

/// \brief Convert a string object, possibly supplied by the user, into an ExceptionResponse
///        enumeration to indicate the program's intended behavior.
///
/// \param policy  Human-readable indication of the exception response behavior
ExceptionResponse translateExceptionResponse(const std::string &policy);
  
/// \brief Translate the exception response behavior into a human-readable string.
///
/// \param choice  The named exception response behavior
std::string getEnumerationName(ExceptionResponse policy);
  
/// \brief Translate a string into a known precision level enumeration.
///
/// \param choice  The named precision model (will be checked for validity)
PrecisionModel translatePrecisionModel(const std::string &choice, ExceptionResponse policy);

/// \brief Get a descriptive string corresponding to each enumerated compute precision model.
///
/// \param pmodel  The precision model to name
std::string getPrecisionModelName(PrecisionModel pmodel);

} // namespace constants
} // namespace stormm

#endif
