#include "copyright.h"
#include "reporting_enumerators.h"

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

} // namespace review
} // namespace stormm
