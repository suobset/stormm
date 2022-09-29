// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace math {

//-------------------------------------------------------------------------------------------------
template <typename T> void TickCounter::randomize(T *xrs) {
  randomize(xrs, std::vector<bool>(state_count, true));
}

//-------------------------------------------------------------------------------------------------
template <typename T> void TickCounter::randomize(T *xrs, const std::vector<bool> &apply_rng) {
  for (int i = 0; i < state_count; i++) {
    if (apply_rng[i]) {
      settings[i] = static_cast<int>(xrs->uniformRandomNumber() * state_limits[i]);
      settings[i] -= ((settings[i] < 0) * settings[i]) +
                     ((settings[i] >= state_limits[i]) * (settings[i] - state_limits[i] + 1));
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T> void TickCounter::randomize(T *xrs, const std::vector<int> &apply_rng) {
  std::vector<bool> bapply_rng(state_count, false);
  const	int nr	= apply_rng.size();
  for (int i = 0; i < nr; i++) {
    bapply_rng[apply_rng[i]] = true;
  }
  randomize(xrs, bapply_rng);
}

} // namespace math
} // namespace stormm
