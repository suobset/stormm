#include "barostat.h"

namespace omni {
namespace trajectory {

std::string getBarostatName(const BarostatKind kind) {
  switch (kind) {
  case BarostatKind::NONE:
    return std::string("NONE");
  case BarostatKind::MONTE_CARLO:
    return std::string("MONTE_CARLO");
  case BarostatKind::BERENDSEN:
    return std::string("BERENDSEN");
  }
  __builtin_unreachable();
}

} // namespace trajectory
} // namespace omni
