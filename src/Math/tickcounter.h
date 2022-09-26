// -*-c++-*-
#ifndef STORMM_TICKCOUNTER_H
#define STORMM_TICKCOUNTER_H

#include <vector>

namespace stormm {
namespace math {

/// \brief Make an array of integral settings and tick it forward (or backward).  The analogy is a
///        series of wheels A - N, each with n_A, n_B, ... , n_N settings.  With each tick, the
///        lowest wheel (position 0 of the array) advances one tick.  When each wheel completes one
///        revolution, it triggers an advance of the next wheel.  This collapses to a simple base-K
///        counting system for a series of wheels with equal numbers of K settings apiece.
class TickCounter {
public:

  /// \brief The constructor takes a series of state counts for each "wheel" in the counter.  The
  ///        initial state can also be provided.
  /// \{
  TickCounter(const std::vector<int> &state_limits_in);
  TickCounter(const std::vector<int> &state_limits_in, const std::vector<int> &settings_in);
  /// \}

  /// \brief The copy and move constructors as well as assignment operators are all valid.
  /// \{
  TickCounter(const TickCounter &original) = default;
  TickCounter(TickCounter &&original) = default;
  TickCounter& operator=(const TickCounter &original) = default;
  TickCounter& operator=(TickCounter &&original) = default;
  /// \}

  /// \brief Get the number of independent settings, the number of variables.
  int getStateCount() const;

  /// \brief Get the logarithm of the number of overall permutations.
  double getLogPermutations() const;
  
  /// \brief Get a const reference to the vector of current settings.
  const std::vector<int>& getSettings() const;

  /// \brief Get a const reference to the vector of numbers of options for each state.
  const std::vector<int>& getStateLimits() const;

  /// \brief Advance the counter by a specified number of steps.
  void advance(int steps = 1);

  /// \brief Reverse the counter by a specified number of steps.
  void reverse(int steps = 1);

  /// \brief Set the object's state to an arbitrary series of values.
  ///
  /// Overloaded:
  ///   - Set all variables
  ///   - Set a specific variable
  ///
  /// \param new_state  The value or values to set
  /// \param var_index  Index of a particular coutner variable to set
  /// \{
  void set(const std::vector<int> &new_state);

  void set(int new_state, int var_index);
  /// \}
  
  /// \brief Reset all states of the counter to zero.
  void reset();
  
private:
  int state_count;
  double log_permutations;
  std::vector<int> settings;
  std::vector<int> state_limits;

  /// \brief Validate the settings to ensure that they are within the stated limits and above zero
  ///        in all counter variables.
  void validateSettings() const;
};

} // namespace math
} // namespace stormm

#endif
