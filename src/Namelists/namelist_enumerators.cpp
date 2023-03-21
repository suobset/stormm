#include "copyright.h"
#include "namelist_enumerators.h"

namespace stormm {
namespace namelist {

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const NamelistType input) {
  switch (input) {
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
std::string getEnumerationName(const InputRepeats input) {
  switch (input) {
  case InputRepeats::NO:
    return std::string("NO");
  case InputRepeats::YES:
    return std::string("YES");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const DefaultIsObligatory input) {
  switch (input) {
  case DefaultIsObligatory::NO:
    return std::string("NO");
  case DefaultIsObligatory::YES:
    return std::string("YES");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const InputStatus input) {
  switch (input) {
  case InputStatus::USER_SPECIFIED:
    return std::string("USER_SPECIFIED");
  case InputStatus::DEFAULT:
    return std::string("DEFAULT");
  case InputStatus::MISSING:
    return std::string("MISSING");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const RestraintAnchoring input) {
  switch (input) {
  case RestraintAnchoring::ATOMMASK:
    return std::string("ATOMMASK");
  case RestraintAnchoring::INDICES:
    return std::string("INDICES");
  case RestraintAnchoring::MIXED:
    return std::string("MIXED");
  case RestraintAnchoring::UNKNOWN:
    return std::string("UNKNOWN");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const SubkeyRequirement input) {
  switch (input) {
  case SubkeyRequirement::OPTIONAL:
    return std::string("OPTIONAL");
  case SubkeyRequirement::REQUIRED:
    return std::string("REQUIRED");
  case SubkeyRequirement::BOGUS:
    return std::string("BOGUS");
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
