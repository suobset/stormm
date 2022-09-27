#include "copyright.h"
#include "Reporting/error_format.h"
#include "multiplication.h"
#include "tickcounter.h"

namespace stormm {
namespace math {

//-------------------------------------------------------------------------------------------------
TickCounter::TickCounter(const std::vector<int> &state_limits_in) :
    state_count{static_cast<int>(state_limits_in.size())},
    log_permutations{logProduct(state_limits_in)},
    settings{std::vector<int>(state_limits_in.size(), 0)},
    state_limits{state_limits_in}
{
  if (state_count == 0) {
    rtErr("No variables were indicated.", "TickCounter");
  }
  int minloc = 0;
  int minval = state_limits[0];
  for (int i = 1; i < state_count; i++) {
    if (state_limits[i] < minval) {
      minval = state_limits[i];
      minloc = i;
    }
  }
  if (minval < 0) {
    rtErr("A variable with " + std::to_string(minval) + " possible states is invalid at "
          "position " + std::to_string(minloc) + ".", "TickCounter");
  }
}

//-------------------------------------------------------------------------------------------------
TickCounter::TickCounter(const std::vector<int> &state_limits_in,
                         const std::vector<int> &settings_in) :
    TickCounter(state_limits_in)
{
  for (int i = 0; i < state_count; i++) {
    settings[i] = settings_in[i];
  }
  validateSettings();
}

//-------------------------------------------------------------------------------------------------
int TickCounter::getStateCount() const {
  return state_count;
}

//-------------------------------------------------------------------------------------------------
double TickCounter::getLogPermutations() const {
  return log_permutations;
}

//-------------------------------------------------------------------------------------------------
const std::vector<int>& TickCounter::getSettings() const {
  return settings;
}

//-------------------------------------------------------------------------------------------------
const std::vector<int>& TickCounter::getStateLimits() const {
  return state_limits;
}

//-------------------------------------------------------------------------------------------------
llint TickCounter::getExactPermutationCount() const {
  if (getLogPermutationCount() >= static_cast<double>(sizeof(llint) - 1LLU) * log(2.0)) {
    rtErr("There are too many permutations to represent as a long long integer.", "TickCounter",
          "getExactPermutationCount");
  }
  return seriesProduct<llint>(state_limits);
}

//-------------------------------------------------------------------------------------------------
double TickCounter::getApproximatePermutationCount() const {
  return seriesProduct<double>(state_limits);
}

//-------------------------------------------------------------------------------------------------
double TickCounter::getLogPermutationCount() const {
  return logProduct(state_limits);
}

//-------------------------------------------------------------------------------------------------
void TickCounter::advance(const int steps) {
  if (steps == 1) {
    int wheel = 0;
    settings[wheel] += 1;
    while (wheel < state_count && settings[wheel] == state_limits[wheel]) {
      settings[wheel] = 0;
      wheel++;
      if (wheel < state_count) {
        settings[wheel] += 1;
      }
    }
  }
  else if (steps > 0) {
    int wheel = 0;
    int balance = steps;
    while (wheel < state_count && balance > 0) {

      // The current wheel ticks forward by the modulo of the remaining balance with the wheel's
      // maximum number of ticks.
      settings[wheel] += balance % state_limits[wheel];

      // Clean up overflows resulting from the current wheel ticking forward.
      int wheel_b = wheel;
      while (wheel_b < state_count && settings[wheel_b] >= state_limits[wheel_b]) {
        settings[wheel_b] -= state_limits[wheel_b];
        wheel_b++;
        if (wheel_b < state_count) {
          settings[wheel_b] += 1;
        }
      }

      // Update the balance that later wheels will need to incorporate.
      balance /= state_limits[wheel];
      wheel++;
    }
  }
  else if (steps == 0) {
    return;
  }
  else {
    rtErr("A forward move of " + std::to_string(steps) + " is invalid.", "TickCounter", "advance");
  }
}

//-------------------------------------------------------------------------------------------------
void TickCounter::reverse(const int steps) {
  if (steps == 1) {
    int wheel = 0;
    settings[wheel] -= 1;
    while (wheel < state_count && settings[wheel] < 0) {
      settings[wheel] = state_limits[wheel] - 1;
      wheel++;
      if (wheel < state_count) {
        settings[wheel] -= 1;
      }
    }
  }
  else if (steps > 0) {
    int wheel = 0;
    int balance = steps;
    while (wheel < state_count && balance > 0) {

      // The current wheel ticks forward by the modulo of the remaining balance with the wheel's
      // maximum number of ticks.
      settings[wheel] -= balance % state_limits[wheel];

      // Clean up overflows resulting from the current wheel ticking forward.
      int wheel_b = wheel;
      while (wheel_b < state_count && settings[wheel_b] < 0) {
        settings[wheel_b] += state_limits[wheel_b];
        wheel_b++;
        if (wheel_b < state_count) {
          settings[wheel_b] -= 1;
        }
      }

      // Update the balance that later wheels will need to incorporate.
      balance /= state_limits[wheel];
      wheel++;
    }
  }
  else if (steps == 0) {
    return;
  }
  else {
    rtErr("A reverse move of " + std::to_string(steps) + " is invalid.", "TickCounter", "reverse");
  }
}

//-------------------------------------------------------------------------------------------------
void TickCounter::set(const std::vector<int> &new_state) {
  if (static_cast<int>(new_state.size()) != state_count) {
    rtErr("The length of the input state (" + std::to_string(new_state.size()) + ") must equal "
          "the number of variables (" + std::to_string(state_count) + ").", "TickCounter", "set");
  }
  for (int i = 0; i < state_count; i++) {
    settings[i] = new_state[i];
  }
  validateSettings();
}

//-------------------------------------------------------------------------------------------------
void TickCounter::set(const int new_state, const int var_index) {
  if (var_index >= state_count || var_index < 0) {
    rtErr("Variable index " + std::to_string(var_index) + " is invalid for " +
          std::to_string(state_count) + " overall counters.", "TickCounter", "set");
  }
  settings[var_index] = new_state;
}

//-------------------------------------------------------------------------------------------------
void TickCounter::reset() {
  for (int i = 0; i < state_count; i++) {
    settings[i] = 0;
  }
}

//-------------------------------------------------------------------------------------------------
void TickCounter::validateSettings() const {
  for (int i = 0; i < state_count; i++) {
    if (settings[i] >= state_limits[i] || settings[i] < 0) {
      rtErr("Counter variable " + std::to_string(i) + " has a maximum of " +
            std::to_string(state_limits[i]) + " possible values.  A setting of " +
            std::to_string(settings[i]) + " is invalid.", "TickCounter", "validateSettings");
    }
  }
}

} // namespace math
} // namespace stormm
