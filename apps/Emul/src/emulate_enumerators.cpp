#include "../../../src/copyright.h"
#include "emulate_enumerators.h"

namespace emulation {

//-------------------------------------------------------------------------------------------------
std::string getEnumeratorName(const FittingContribution input) {
  switch (input) {
  case FittingContribution::ENERGY:
    return std::string("ENERGY");
  case FittingContribution::FORCE:
    return std::string("FORCE");
  }
  __builtin_unreachable();
}

} // namespace emulation
