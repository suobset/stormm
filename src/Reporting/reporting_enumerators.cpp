#include "copyright.h"
#include "Constants/behavior.h"
#include "Parsing/parse.h"
#include "reporting_enumerators.h"
#include "error_format.h"

namespace stormm {
namespace display {

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const HelpSignalKind input) {
  switch (input) {
  case HelpSignalKind::NO_ARGS:
    return std::string("NO_ARGS");
  case HelpSignalKind::NO_ARGS_ONLY:
    return std::string("NO_ARGS_ONLY");
  case HelpSignalKind::KEYWORD:
    return std::string("KEYWORD");
  case HelpSignalKind::KEYWORD_ONLY:
    return std::string("KEYWORD_ONLY");
  }
  __builtin_unreachable();
}

} // namespace display
  
namespace review {

using constants::CaseSensitivity;
using parse::strcmpCased;

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const OutputScope input) {
  switch (input) {
  case OutputScope::AVERAGES:
    return std::string("AVERAGES");
  case OutputScope::OUTLIERS:
    return std::string("OUTLIERS");
  case OutputScope::CLUSTER_AVERAGES:
    return std::string("CLUSTER_AVERAGES");
  case OutputScope::CLUSTER_OUTLIERS:
    return std::string("CLUSTER_OUTLIERS");
  case OutputScope::FULL:
    return std::string("FULL");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const OutputSyntax input) {
  switch (input) {
  case OutputSyntax::MATPLOTLIB:
    return std::string("MATPLOTLIB");
  case OutputSyntax::MATRIX_PKG:
    return std::string("MATRIX_PKG");
  case OutputSyntax::STANDALONE:
    return std::string("STANDALONE");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const GridFileSyntax input) {
  switch (input) {
  case GridFileSyntax::MATPLOTLIB:
    return std::string("MATPLOTLIB");
  case GridFileSyntax::MATRIX_PKG:
    return std::string("MATRIX_PKG");
  case GridFileSyntax::OPEN_DX:
    return std::string("OPEN_DX");
  case GridFileSyntax::CUBEGEN:
    return std::string("CUBEGEN");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const SectionComponent input) {
  switch (input) {
  case SectionComponent::NARRATIVE:
    return std::string("NARRATIVE");
  case SectionComponent::LIST:
    return std::string("LIST");
  case SectionComponent::TABLE:
    return std::string("TABLE");
  case SectionComponent::SCRIPT:
    return std::string("SCRIPT");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const ListEnumeration input) {
  switch (input) {
  case ListEnumeration::BULLET:
    return std::string("BULLET");
  case ListEnumeration::NUMBERED:
    return std::string("NUMBERED");
  case ListEnumeration::ALPHABETIC:
    return std::string("ALPHABETIC");
  case ListEnumeration::ROMAN:
    return std::string("ROMAN");
  case ListEnumeration::NONE:
    return std::string("NONE");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const TextEnds input) {
  switch (input) {
  case TextEnds::AS_IS:
    return std::string("AS_IS");
  case TextEnds::NEWLINE:
    return std::string("NEWLINE");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const TableContentKind input) {
  switch (input) {
  case TableContentKind::INTEGER:
    return std::string("INTEGER");
  case TableContentKind::REAL:
    return std::string("REAL");
  case TableContentKind::STRING:
    return std::string("STRING");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const SurfaceRender input) {
  switch (input) {
  case SurfaceRender::WIRE:
    return std::string("WIRE");
  case SurfaceRender::SOLID:
    return std::string("SOLID");
  case SurfaceRender::SCAFFOLD:
    return std::string("SCAFFOLD");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const LinePlotStyle input) {
  switch (input) {
  case LinePlotStyle::SOLID:
    return std::string("SOLID");
  case LinePlotStyle::DASHED:
    return std::string("DASHED");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
SurfaceRender translateSurfaceRender(const std::string &input) {
  if (strcmpCased(input, "wire", CaseSensitivity::NO) ||
      strcmpCased(input, "wire mesh", CaseSensitivity::NO) ||
      strcmpCased(input, "wire_mesh", CaseSensitivity::NO)) {
    return SurfaceRender::WIRE;
  }
  else if (strcmpCased(input, "solid", CaseSensitivity::NO)) {
    return SurfaceRender::SOLID;
  }
  else if (strcmpCased(input, "scaffold", CaseSensitivity::NO)) {
    return SurfaceRender::SCAFFOLD;    
  }
  else {
    rtErr("No SurfaceRender enumeration matches \"" + input + "\".", "translateSurfaceRender");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
LinePlotStyle translateLinePlotStyle(const std::string &input) {
  if (strcmpCased(input, "solid", CaseSensitivity::NO) ||
      strcmpCased(input, "continuous", CaseSensitivity::NO)) {
    return LinePlotStyle::SOLID;
  }
  else if (strcmpCased(input, "dashed", CaseSensitivity::NO) ||
           strcmpCased(input, "broken", CaseSensitivity::NO)) {
    return LinePlotStyle::DASHED;
  }
  else {
    rtErr("No LinePlotStyle enumeration matches \"" + input + "\".", "translateLinePlotStyle");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
const std::string& toRoman(const int x, const ExceptionResponse policy) {
  if (x <= 0 || x > maximum_roman_numeral) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("No Roman numeral is available for " + std::to_string(x) + ".", "toRoman");
      break;
    case ExceptionResponse::WARN:
      rtWarn("No Roman numeral is available for " + std::to_string(x) + ".  White space will be "
             "returned instead.", "toRoman");
      return roman_numerals[maximum_roman_numeral];
    case ExceptionResponse::SILENT:
      return roman_numerals[maximum_roman_numeral];
    }
  }
  else {
    return roman_numerals[x - 1];
  }
  __builtin_unreachable();
}
  
} // namespace review
} // namespace stormm
