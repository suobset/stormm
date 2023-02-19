#include "copyright.h"
#include "Constants/behavior.h"
#include "Parsing/parse.h"
#include "synthesis_enumerators.h"

namespace stormm {
namespace synthesis {

using constants::CaseSensitivity;
using parse::strcmpCased;

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const NbwuKind input) {
  switch (input) {
  case NbwuKind::TILE_GROUPS:
    return std::string("TILE_GROUPS");
  case NbwuKind::SUPERTILES:
    return std::string("SUPERTILES");
  case NbwuKind::HONEYCOMB:
    return std::string("HONEYCOMB");
  case NbwuKind::UNKNOWN:
    return std::string("UNKNOWN");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const VwuTask input) {
  switch (input) {
  case VwuTask::BOND:
    return std::string("BOND");
  case VwuTask::ANGL:
    return std::string("ANGL");
  case VwuTask::DIHE:
    return std::string("DIHE");
  case VwuTask::UBRD:
    return std::string("UBRD");
  case VwuTask::CBND:
    return std::string("CBND");
  case VwuTask::CIMP:
    return std::string("CIMP");
  case VwuTask::CDHE:
    return std::string("CDHE");
  case VwuTask::CMAP:
    return std::string("CMAP");
  case VwuTask::INFR14:
    return std::string("INFR14");
  case VwuTask::RPOSN:
    return std::string("RPOSN");
  case VwuTask::RBOND:
    return std::string("RBOND");
  case VwuTask::RANGL:
    return std::string("RANGL");
  case VwuTask::RDIHE:
    return std::string("RDIHE");
  case VwuTask::VSITE:
    return std::string("VSITE");
  case VwuTask::SETTLE:
    return std::string("SETTLE");
  case VwuTask::CGROUP:
    return std::string("CGROUP");
  case VwuTask::ALL_TASKS:
    return std::string("ALL_TASKS");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const VwuGoal input) {
  switch (input) {
  case VwuGoal::ACCUMULATE:
    return std::string("ACCUMULATE");
  case VwuGoal::MOVE_PARTICLES:
    return std::string("MOVE_PARTICLES");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const InitializationTask input) {
  switch (input) {
  case InitializationTask::GENERAL_DYNAMICS:
    return std::string("GENERAL_DYNAMICS");
  case InitializationTask::GB_DYNAMICS:
    return std::string("GB_DYNAMICS");
  case InitializationTask::GENERAL_MINIMIZATION:
    return std::string("GENERAL_MINIMIZATION");
  case InitializationTask::GB_MINIMIZATION:
    return std::string("GB_MINIMIZATION");
  case InitializationTask::LANGEVIN_DYNAMICS:
    return std::string("LANGEVIN_DYNAMICS");
  case InitializationTask::GB_LANGEVIN_DYNAMICS:
    return std::string("GB_LANGEVIN_DYNAMICS");
  case InitializationTask::NONE:
    return std::string("NONE");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const CondensateSource input) {
  switch (input) {
  case CondensateSource::SYNTHESIS:
    return std::string("SYNTHESIS");
  case CondensateSource::SERIES:
    return std::string("SERIES");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const SystemGrouping input) {
  switch (input) {
  case SystemGrouping::SOURCE:
    return std::string("SOURCE");
  case SystemGrouping::TOPOLOGY:
    return std::string("TOPOLOGY");
  case SystemGrouping::LABEL:
    return std::string("LABEL");
  }
  __builtin_unreachable();
}
  
//-------------------------------------------------------------------------------------------------
SystemGrouping translateSystemGrouping(const std::string &input_string) {
  if (strcmpCased(input_string, "system", CaseSensitivity::NO) ||
      strcmpCased(input_string, "source", CaseSensitivity::NO) ||
      strcmpCased(input_string, "sys", CaseSensitivity::NO)) {
    return SystemGrouping::SOURCE;
  }
  else if (strcmpCased(input_string, "topology", CaseSensitivity::NO)) {
    return SystemGrouping::TOPOLOGY;
  }
  else if (strcmpCased(input_string, "label", CaseSensitivity::NO)) {
    return SystemGrouping::LABEL;
  }
  else {
    rtErr("The input \"" + input_string + "\" is not a valid system grouping.",
          "translateSystemGrouping");
  }
  __builtin_unreachable();
}

} // namespace synthesis
} // namespace stormm
