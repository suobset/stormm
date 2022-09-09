#include "copyright.h"
#include "Parsing/parse.h"
#include "Reporting/error_format.h"
#include "behavior.h"

namespace stormm {
namespace constants {

using parse::strcmpCased;

//-------------------------------------------------------------------------------------------------
PrecisionModel translatePrecisionModel(const std::string &choice, const ExceptionResponse policy) {
  if (strcmpCased(choice, std::string("single"))) {
    return PrecisionModel::SINGLE;
  }
  else if (strcmpCased(choice, std::string("double"))) {
    return PrecisionModel::DOUBLE;
  }
  else {
    rtErr("Invalid request for precision level " + choice + ".", "translatePrecisionModel");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getPrecisionModelName(const PrecisionModel plevel) {
  switch (plevel) {
  case PrecisionModel::SINGLE:
    return std::string("SINGLE");
  case PrecisionModel::DOUBLE:
    return std::string("DOUBLE");
  }
  __builtin_unreachable();
}

} // namespace constants
} // namespace stormm
