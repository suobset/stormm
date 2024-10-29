#include "copyright.h"
#include "Parsing/parse.h"
#include "Parsing/parsing_enumerators.h"
#include "numeric_enumerators.h"

namespace stormm {
namespace numerics {

using parse::CaseSensitivity;
using parse::strcmpCased;
  
//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const AccumulationMethod input) {
  switch (input) {
  case AccumulationMethod::SPLIT:
    return std::string("SPLIT");
  case AccumulationMethod::WHOLE:
    return std::string("WHOLE");
  case AccumulationMethod::AUTOMATIC:
    return std::string("AUTOMATIC");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
AccumulationMethod translateAccumulationMethod(const std::string &input) {
  if (strcmpCased(input, "SPLIT", CaseSensitivity::NO)) {
    return AccumulationMethod::SPLIT;
  }
  else if (strcmpCased(input, "WHOLE", CaseSensitivity::NO) ||
           strcmpCased(input, "FUSED", CaseSensitivity::NO)) {
    return AccumulationMethod::WHOLE;
  }
  else if (strcmpCased(input, "AUTO", CaseSensitivity::NO) ||
           strcmpCased(input, "AUTOMATIC", CaseSensitivity::NO)) {
    return AccumulationMethod::AUTOMATIC;
  }
  else {
    rtErr("Invalid specifier \"" + input + "\".  Valid inputs include \"SPLIT\", \"WHOLE\", and "
          "\"AUTO\".", "translateAccumulationMethod");
  }
  __builtin_unreachable();
}

} // namespace numerics
} // namespace stormm
