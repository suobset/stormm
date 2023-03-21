#include "copyright.h"
#include "Constants/behavior.h"
#include "Parsing/parse.h"
#include "Reporting/error_format.h"
#include "structure_enumerators.h"

namespace stormm {
namespace structure {

using constants::CaseSensitivity;
using parse::strcmpCased;


//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const ImagingMethod input) {
  switch (input) {
  case ImagingMethod::PRIMARY_UNIT_CELL:
    return std::string("PRIMARY_UNIT_CELL");
  case ImagingMethod::MINIMUM_IMAGE:
    return std::string("MINIMUM_IMAGE");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const RMSDMethod input) {
  switch (input) {
  case RMSDMethod::ALIGN_MASS:
    return std::string("ALIGN_MASS");
  case RMSDMethod::ALIGN_GEOM:
    return std::string("ALIGN_GEOM");
  case RMSDMethod::NO_ALIGN_MASS:
    return std::string("NO_ALIGN_MASS");
  case RMSDMethod::NO_ALIGN_GEOM:
    return std::string("NO_ALIGN_GEOM");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const RMSDAlignmentProtocol input) {
  switch (input) {
  case RMSDAlignmentProtocol::BUILD_CORE:
    return std::string("BUILD_CORE");
  case RMSDAlignmentProtocol::ALIGN_CORE:
    return std::string("ALIGN_CORE");
  case RMSDAlignmentProtocol::ALIGN_ALL:
    return std::string("ALIGN_ALL");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const RMSDTask input) {
  switch (input) {
  case RMSDTask::REFERENCE:
    return std::string("REFERENCE");
  case RMSDTask::MATRIX:
    return std::string("MATRIX");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const VirtualSiteActivity input) {
  switch (input) {
  case VirtualSiteActivity::PLACEMENT:
    return std::string("PLACEMENT");
  case VirtualSiteActivity::TRANSMIT_FORCES:
    return std::string("TRANSMIT_FORCES");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const GridDetail input) {
  switch (input) {
  case GridDetail::OCCLUSION:
    return std::string("OCCLUSION");
  case GridDetail::NONBONDED_FIELD:
    return std::string("NONBONDED_FIELD");
  case GridDetail::NONBONDED_ATOMIC:
    return std::string("NONBONDED_ATOMIC");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const ClashKind input) {
  switch (input) {
  case ClashKind::VAN_DER_WAALS:
    return std::string("VAN_DER_WAALS");
  case ClashKind::PURE_DISTANCE:
    return std::string("PURE_DISTANCE");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const SamplingIntensity input) {
  switch (input) {
  case SamplingIntensity::MINIMAL:
    return std::string("MINIMAL");
  case SamplingIntensity::LIGHT:
    return std::string("LIGHT");
  case SamplingIntensity::HEAVY:
    return std::string("HEAVY");
  case SamplingIntensity::EXHAUSTIVE:
    return std::string("EXHAUSTIVE");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const BoundaryCondition input) {
  switch (input) {
  case BoundaryCondition::ISOLATED:
    return std::string("ISOLATED");
  case BoundaryCondition::PERIODIC:
    return std::string("PERIODIC");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
SamplingIntensity translateSamplingIntensity(const std::string &input) {
  if (strcmpCased(input, std::string("minimal"), CaseSensitivity::NO)) {
    return SamplingIntensity::MINIMAL;
  }
  else if (strcmpCased(input, std::string("light"), CaseSensitivity::NO)) {
    return SamplingIntensity::LIGHT;
  }
  else if (strcmpCased(input, std::string("heavy"), CaseSensitivity::NO)) {
    return SamplingIntensity::HEAVY;
  }
  else if (strcmpCased(input, std::string("exhaustive"), CaseSensitivity::NO)) {
    return SamplingIntensity::EXHAUSTIVE;
  }
  else {
    rtErr("\"" + input + "\" is not a valid sampling intensity.", "translateSamplingIntensity");
  }
  __builtin_unreachable();
}

} // namespace structure
} // namespace stormm
