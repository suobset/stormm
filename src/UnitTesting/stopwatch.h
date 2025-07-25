// -*-c++-*-
#ifndef STORMM_BENCHMARK_TEST_H
#define STORMM_BENCHMARK_TEST_H

#include <string>
#include <vector>
#include "copyright.h"
#include "Constants/behavior.h"

namespace stormm {
namespace testing {

using constants::ExceptionResponse;
  
/// \brief Object for managing calls to the C-standard function gettimeofday(), calculating deltas
///        and categorizing time spent according to a developer's wishes.  Calls to gettimeofday()
///        are precise to microseconds, and take considerably less time than that.  With the option
///        to lookup a section by numerical index, the timing can have no impact whatsoever on any
///        serious computation, and even lookup by section name is relatively fast thanks to the
///        findStringInVector() function.  The relative time can be computed in a number of ways,
///        so that the timer can be called less frequently but not accrue all of the time since
///        it last sampled the computation.
class StopWatch {
public:

  /// \brief The basic constructor records the time at which initialization was called and creates
  ///        a section called "Miscellaneous" to capture any times not assigned to a particular
  ///        section.
  ///
  /// \param title_in  Title to give this timing apparatus (default "Timer")
  StopWatch(const std::string &title_in = std::string("Timer"));
  
  /// \brief Default destructor
  ~StopWatch() = default;

  /// \brief With no const members or points to repair, the default copy and move constructors, as
  ///        well as copy and move assignment operators, will be applicable.
  ///
  /// \param original  The original object to copy or move
  /// \param other     Another object placed on the right hand side of the assignment statement
  /// \{
  StopWatch(const StopWatch &original) = default;
  StopWatch(StopWatch &&original) = default;
  StopWatch& operator=(const StopWatch &original) = default;
  StopWatch& operator=(StopWatch &&original) = default;
  /// \}
  
  /// \brief Get the time at which this StopWatch was first started
  double getTimeAtStart() const;

  /// \brief Get the time at which this StopWatch last recorded a test (this is not necessarily the
  ///        last time this StopWatch fired off translateCurrentTime() which contains the actual
  ///        call to the underlying gettimeofday(), as there are cases where one can simply query
  ///        how long ago the StopWatch did anything).
  double getTimeAtLastTest() const;

  /// \brief Get the time since this StopWatch was first started.
  double getTimeSinceStart() const;

  /// \brief Get the time since this StopWatch was last tested (see above for lengthier
  ///        explanation).
  double getTimeSinceLastTest() const;

  /// \brief Get the duration recorded under any one category
  ///
  /// Overloaded:
  ///   - Identify the category by name
  ///   - Identify the category by the order in which it was assigned
  ///
  /// \param query_index  Index of the category (the first added category has index 1)
  /// \param query_name   Name of the category
  /// \{
  double getCategoryDuration(int query_index) const;
  double getCategoryDuration(const std::string &query_name) const;
  /// \}

  /// \brief Get the average interval of time observed for any given category
  ///
  /// Overloaded:
  ///   - Identify the category by name
  ///   - Identify the category by the order in which it was assigned
  ///
  /// \param query_index  Index of the category (the first added category has index 1)
  /// \param query_name   Name of the category
  /// \{
  double getCategoryAverageInterval(int query_index) const;
  double getCategoryAverageInterval(const std::string &query_name) const;
  /// \}

  /// \brief Get the minimum stretch of time recorded under any one category
  ///
  /// Overloaded:
  ///   - Identify the category by name
  ///   - Identify the category by the order in which it was assigned
  ///
  /// \param query_index  Index of the category (the first added category has index 1)
  /// \param query_name   Name of the category
  /// \{
  double getCategoryMinimumTime(int query_index) const;
  double getCategoryMinimumTime(const std::string &query_name) const;
  /// \}

  /// \brief Get the maximum stretch of time recorded under any one category
  ///
  /// Overloaded:
  ///   - Identify the category by name
  ///   - Identify the category by the order in which it was assigned
  ///
  /// \param query_index  Index of the category (the first added category has index 1)
  /// \param query_name   Name of the category
  /// \{
  double getCategoryMaximumTime(int query_index) const;
  double getCategoryMaximumTime(const std::string &query_name) const;
  /// \}

  /// \brief Get the maximum stretch of time recorded under any one category
  ///
  /// Overloaded:
  ///   - Identify the category by name
  ///   - Identify the category by the order in which it was assigned
  ///
  /// \param query_index  Index of the category (the first added category has index 1)
  /// \param query_name   Name of the category
  /// \{
  int getCategorySamples(int query_index) const;
  int getCategorySamples(const std::string &query_name) const;
  /// \}

  /// \brief Get the name of a timing category based on its index in the program
  ///
  /// \param query_index  Index of the category (the first added category has index 1)
  std::string getCategoryName(int query_index) const;

  /// \brief Get the index of category based on its name.  Returns -1 or an error if the named
  ///        category is not present in the object.
  ///
  /// \param query   Name of the timings category to seek out
  /// \param policy  Course of action if the named category is not found
  int getCategoryIndex(const std::string &query,
                       ExceptionResponse policy = ExceptionResponse::DIE) const;  
  
  /// \brief Report the total duration recorded by this stopwatch under all sections.
  double getTotalDuration() const;
  
  /// \brief Assign an amount of time to this StopWatch.
  ///
  /// Overloaded:
  ///   - Assign time to the "Miscellaneous" section (no arguments)
  ///   - Assign time to a particular section, by index or by name, since the last time recorded by
  ///     that section
  ///   - Assign time to a particular section based on an arbitrary reference time (this could be
  ///     an independent call to gettimeofday() stored earlier in some variable, or the last time
  ///     recorded by some other section of this or even another StopWatch).
  ///
  /// \param target          The index or name of the section that will get time added
  /// \param reference_time  Arbitrary reference time from some other source
  /// \{
  void assignTime(int target, double reference_time);
  void assignTime(int target = 0);
  void assignTime(const std::string &target);
  void assignTime(const std::string &target, double reference_time);
  /// \}

  /// \brief Add a section to the current StopWatch, if no such section already exists.  Return
  ///        the index of the section, whether newly added or not.
  ///
  /// \param name  Proposed name of the new section
  int addCategory(const std::string &name);

  /// \brief Print the timings for this StopWatch.
  ///
  /// \param precision  Results will be reported such that the largest printed interval can be
  ///                   known to within one part in at least this amount. However, the maximum
  ///                   number of decimals that will be printed is 4, one 10,000th of a second.
  void printResults(double precision = 1.0e6);

  /// \brief Obtain the time taken as the difference of two categories, catg_a - catg_b.  The
  ///        result is returned as the mean differential, which is exact, and an estimate of the
  ///        standard deviation of the differential (based on the standard deviations of each
  ///        category individually).
  ///
  /// Overloaded:
  ///   - Provide the categories by their string identifiers
  ///   - Provide the categories by their integer identifiers
  ///
  /// \param catg_a              The first category, likely the longer of two categories or a
  ///                            category which comprises some or all of the activity in the second
  ///                            category
  /// \param catg_b              The second category, likely to be composed of a subset of the
  ///                            activities in the first category
  /// \param referring_function  The function calling for a differential (for error tracing)
  /// \{
  double2 timeDifferential(const std::string &catg_a, const std::string &catg_b,
                           const std::string &referring_function = std::string("")) const;

  double2 timeDifferential(int catg_a, int catg_b,
                           const std::string &referring_function = std::string("")) const;
  /// \}
  
private:
  int category_count;                         ///< Number of categories tracked by this StopWatch
  double initial_time;                        ///< The time at which this StopWatch was initialized
  double time_at_last_test;                   ///< The time at which any section of this StopWatch
                                              ///<   was last tested
  std::vector<double> category_start_times;   ///< The initial times at which each category began
  std::vector<double> category_total_times;   ///< Total accumulated times of each category 
  std::vector<double> category_squared_times; ///< Total squared times of samples in each category
                                              ///<   thus far (kept for running standard deviation)
  std::vector<double> category_last_times;    ///< Most recent time recorded under each category
  std::vector<double> category_min_interval;  ///< Shortest length of time clocked in each category
  std::vector<double> category_max_interval;  ///< Longest length of time clocked to each category
  std::vector<int> category_samples;          ///< The number of occasions on which time has been
                                              ///<   assigned to any given category
  std::vector<std::string> category_names;    ///< Names of each category
  std::string label;                          ///< Label for the timing data set as a whole

  /// \brief Convert the system time into a more convenient real number.
  double translateCurrentTime() const;

  /// \brief Compute the differential between two categories with trusted indices.  The public
  ///        member function timeDifferential() calls this private member function.
  ///
  /// \param catg_a  The first of the two categories
  /// \param catg_b  The second of the two categories
  double2 computeDifferential(int catg_a, int catg_b) const;
  
  /// \brief Determine if a section index is valid for this StopWatch.
  ///
  /// \param index               Index of the section in question
  /// \param referring_function  Function that called this check
  void validateCategoryIndex(int index, const std::string &referring_function) const;

  /// \brief Determine if a section name is valid for this StopWatch.  Return the index if so, or
  ///         print an error otherwise.
  ///
  /// \param name                Name of the section in question
  /// \param referring_function  Function that called this check
  int validateCategoryName(const std::string &query, const std::string &referring_function) const;
};
  
} // namespace testing
} // namespace stormm

#endif
