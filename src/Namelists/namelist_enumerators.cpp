#include "copyright.h"
#include "Constants/behavior.h"
#include "Parsing/parse.h"
#include "namelist_enumerators.h"

namespace stormm {
namespace namelist {

using constants::CaseSensitivity;
using parse::strcmpCased;

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const NamelistType input) {
  switch (input) {
  case NamelistType::BOOLEAN:
    return std::string("BOOLEAN");
    break;
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
std::string getEnumerationName(const KeyRequirement input) {
  switch (input) {
  case KeyRequirement::OPTIONAL:
    return std::string("OPTIONAL");
  case KeyRequirement::REQUIRED:
    return std::string("REQUIRED");
  case KeyRequirement::BOGUS:
    return std::string("BOGUS");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const NamelistIntroduction input) {
  switch (input) {
  case NamelistIntroduction::HEADER:
    return std::string("HEADER");
  case NamelistIntroduction::COMPACT_HEADER:
    return std::string("COMPACT_HEADER");
  case NamelistIntroduction::BLANK_LINE:
    return std::string("BLANK_LINE");
  case NamelistIntroduction::NONE:
    return std::string("NONE");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const DayOfTheWeek input) {
  switch (input) {
  case DayOfTheWeek::SUN:
    return std::string("SUN");
  case DayOfTheWeek::MON:
    return std::string("MON");
  case DayOfTheWeek::TUE:
    return std::string("TUE");
  case DayOfTheWeek::WED:
    return std::string("WED");
  case DayOfTheWeek::THU:
    return std::string("THU");
  case DayOfTheWeek::FRI:
    return std::string("FRI");
  case DayOfTheWeek::SAT:
    return std::string("SAT");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
KeyRequirement translateKeyRequirement(const std::string &input) {
  if (strcmpCased(input, "optional", CaseSensitivity::NO)) {
    return KeyRequirement::OPTIONAL;
  }
  else if (strcmpCased(input, "required", CaseSensitivity::NO)) {
    return KeyRequirement::REQUIRED;
  }
  else if (strcmpCased(input, "bogus", CaseSensitivity::NO)) {
    return KeyRequirement::BOGUS;
  }
  else {
    rtErr("No translation of KeyRequirement \"" + input + "\" is available.",
          "translateKeyRequirement");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
DayOfTheWeek translateDayOfTheWeek(const std::string &input) {
  if (strcmpCased(input, "sun", CaseSensitivity::NO) ||	
      strcmpCased(input, "sunday", CaseSensitivity::NO) ||
      strcmpCased(input, "s", CaseSensitivity::NO)) {
    return DayOfTheWeek::SUN;
  }
  else if (strcmpCased(input, "mon", CaseSensitivity::NO) ||
           strcmpCased(input, "monday", CaseSensitivity::NO) ||
           strcmpCased(input, "m", CaseSensitivity::NO)) {      
    return DayOfTheWeek::MON;
  }
  else if (strcmpCased(input, "tue", CaseSensitivity::NO) ||	
           strcmpCased(input, "tuesday", CaseSensitivity::NO) ||
           strcmpCased(input, "t", CaseSensitivity::NO)) {
    return DayOfTheWeek::TUE;
  }
  else if (strcmpCased(input, "wed", CaseSensitivity::NO) ||	
           strcmpCased(input, "wednesday", CaseSensitivity::NO) ||
           strcmpCased(input, "w", CaseSensitivity::NO)) {
    return DayOfTheWeek::WED;
  }
  else if (strcmpCased(input, "thu", CaseSensitivity::NO) ||	
           strcmpCased(input, "thursday", CaseSensitivity::NO) ||
           strcmpCased(input, "th", CaseSensitivity::NO)) {
    return DayOfTheWeek::THU;
  }
  else if (strcmpCased(input, "fri", CaseSensitivity::NO) ||	
           strcmpCased(input, "friday", CaseSensitivity::NO) ||
           strcmpCased(input, "f", CaseSensitivity::NO)) {
    return DayOfTheWeek::FRI;
  }
  else if (strcmpCased(input, "sat", CaseSensitivity::NO) ||	
           strcmpCased(input, "saturday", CaseSensitivity::NO) ||
           strcmpCased(input, "sa", CaseSensitivity::NO)) {
    return DayOfTheWeek::SAT;
  }
  else {
    rtErr("Invalid day \"" + input + ".\"", "translateDayOfTheWeek");
  }
  __builtin_unreachable();
}

} // namespace namelist
} // namespace stormm
