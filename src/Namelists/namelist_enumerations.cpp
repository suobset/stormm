#include "copyright.h"
#include "namelist_enumerations.h"

namespace stormm {
namespace namelist {

//-------------------------------------------------------------------------------------------------
std::string getNamelistTypeName(const NamelistType param_type) {
  switch (param_type) {
  case NamelistType::INTEGER:
    return std::string("INTEGER");
    break;
  case NamelistType::REAL:
    return std::string("REAL");
    break;
  case NamelistType::STRING:
    return std::string("STRING");
    break;
  case NamelistType::STRUCT:
    return std::string("STRUCT");
    break;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getInputStatusString(const InputStatus stt) {
  switch (stt) {
  case InputStatus::USER_SPECIFIED:
    return std::string("USER SPECIFIED");
  case InputStatus::DEFAULT:
    return std::string("DEFAULT");
  case InputStatus::MISSING:
    return std::string("MISSING");
  }
  __builtin_unreachable();
}

} // namespace namelist
} // namespace stormm
