#include "copyright.h"
#include <string>
#include "unit_test_enumerators.h"

namespace stormm {
namespace testing {

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const TestPriority x) {
  switch (x) {
  case TestPriority::CRITICAL:
    return std::string("CRITICAL");
  case TestPriority::NON_CRITICAL:
    return std::string("NON_CRITICAL");
  case TestPriority::ABORT:
    return std::string("ABORT");
  }
  __builtin_unreachable();
}

} // namespace testing
} // namespace stormm
