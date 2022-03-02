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
class Thermostat {
public:

  /// \brief Constructors can be blank (implying a thermostat of kind NONE), take a specific kind
  ///        (implying default values for that kind), or take a kind and specific settings.
  ///
  /// \param kind_in  The type of thermostat to implement
  /// \{
  Thermostat();
  Thermostat(ThermostatKind kind_in);
  /// \}

  /// \brief Get the kind of thermostat
  ThermostatKind getKind() const;
  
private:
  ThermostatKind kind;      ///< The type of thermostat
};

/// \brief Return the name of the thermostat choice (an enumerator string conversion function)
///
/// \param kind  The type of thermostat
std::string getThermostatName(ThermostatKind kind);

} // namespace trajectory
} // namespace omni

#endif
