#include "thermostat.h"

namespace omni {
namespace trajectory {

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
} // namespace omni
