// -*-c++-*-
#ifndef OMNI_ENERGY_ENUMERATORS_H
#define OMNI_ENERGY_ENUMERATORS_H

namespace omni {
namespace energy {

/// \brief Enumerate the choices on whether to evaluate the force... yes or no.
enum class EvaluateForce {
  NO, YES
};

/// \brief Enumerate the choices on whether to evaluate the virial... yes or no.
enum class EvaluateVirial {
  NO, YES
};

} // namespace energy
} // namespace omni

#endif
