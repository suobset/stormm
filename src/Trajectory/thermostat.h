// -*-c++-*-
#ifndef STORMM_THERMOSTAT_H
#define STORMM_THERMOSTAT_H

#include <string>
#include "copyright.h"
#include "Accelerator/hybrid.h"

namespace stormm {
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
  ThermostatKind kind;             ///< The type of thermostat
  int total_atoms;                 ///< The total number of atoms, and thus the number of random
                                   ///<   state vectors, for which this thermostat is responsible
                                   ///<   (this gives the length of random_sv_bank)
  Hybrid<ullint4> random_sv_bank;  ///< Bank of Xoshiro256++ state vectors for creating random
                                   ///<   numbers.  Each atom has its own generator.
};

/// \brief Return the name of the thermostat choice (an enumerator string conversion function)
///
/// \param kind  The type of thermostat
std::string getThermostatName(ThermostatKind kind);

} // namespace trajectory
} // namespace stormm

#endif
