// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace testing {

//-------------------------------------------------------------------------------------------------
template <typename T>
std::string listItemsAsString(const std::vector<T> &list, const size_t maximum_items,
                              const bool spread) {
  std::string result;
  const size_t n_display = std::min(maximum_items, list.size());
  if (n_display == maximum_items) {
    for (size_t i = 0; i < n_display; i++) {
      result += std::to_string(list[i]) + listSeparator(i, n_display);
    }
  }
  else {
    if (spread) {
      double prog = 0.0;
      const double frac = static_cast<double>(list.size()) / static_cast<double>(n_display);
      for (int i = 0; i < n_display; i++) {
        const size_t idx = static_cast<size_t>(prog);
        result += std::to_string(list[idx]) + listSeparator(i, n_display);
        prog += frac;
      }
    }
  }
  return result;
}

} // namespace testing
} // namespace stormm
