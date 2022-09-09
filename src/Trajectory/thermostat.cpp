#include "copyright.h"
#include "thermostat.h"

namespace stormm {
namespace trajectory {

//-------------------------------------------------------------------------------------------------
Thermostat::Thermostat() :
    kind{ThermostatKind::NONE}
{}

//-------------------------------------------------------------------------------------------------
Thermostat::Thermostat(const ThermostatKind kind_in) :
    kind{kind_in}
{}

//-------------------------------------------------------------------------------------------------
ThermostatKind Thermostat::getKind() const {
  return kind;
}

//-------------------------------------------------------------------------------------------------
std::string getThermostatName(const ThermostatKind kind) {
  switch (kind) {
  case ThermostatKind::NONE:
    return std::string("NONE");
  case ThermostatKind::ANDERSEN:
    return std::string("ANDERSEN");
  case ThermostatKind::LANGEVIN:
    return std::string("LANGEVIN");
  case ThermostatKind::BERENDSEN:
    return std::string("BERENDSEN");
  }
  __builtin_unreachable();
}

} // namespace trajectory
} // namespace stormm
