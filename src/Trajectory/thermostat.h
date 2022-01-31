// -*-c++-*-
#ifndef OMNI_THERMOSTAT_H
#define OMNI_THERMOSTAT_H

#include <string>

namespace omni {
namespace trajectory {

/// \brief Enumerate the various thermostats available for simulations
enum class ThermostatKind {
  NONE, ANDERSEN, LANGEVIN, BERENDSEN
};

/// \brief Store the parameters for a simulation thermostat.  Includes Berendsen, Andersen, and
///        Langevin methods.
struct Thermostat {
  ThermostatKind kind;
};

/// \brief Return the name of the thermostat choice (an enumerator string conversion function)
///
/// \param kind  The type of thermostat
std::string getThermostatName(ThermostatKind kind);

} // namespace trajectory
} // namespace omni

#endif
