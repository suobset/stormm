// -*-c++-*-
#ifndef STORMM_NML_NICE_H
#define STORMM_NML_NICE_H

#include <sys/time.h>
#include "copyright.h"
#include "Parsing/textfile.h"
#include "namelist_emulator.h"
#include "namelist_enumerators.h"

namespace stormm {
namespace namelist {

using parse::TextFile;
using parse::WrapTextSearch;

/// \brief The default working hours (local time according to the computer is assumed)
/// \{
constexpr int default_workday_start_hour = 9;
constexpr int default_workday_end_hour = 17;
constexpr int default_workday_start_minute = 0;
constexpr int default_workday_end_minute = 0;
/// \}

/// \brief The default work week extends from Monday to Friday and does not include holidays.
///        The work week extends through the end of the last day indicated.
/// \{
constexpr char default_workweek_start[] = "MON";
constexpr char default_workweek_end[] = "FRI";
/// \}

/// \brief The default workday GPU level, barring any competing process, is 50%.
constexpr double default_workday_gpu_effort = 0.50;

/// \brief Bounds on the time interval between checks on GPU activity while the STORMM process is
///        running, in wall clock seconds
/// \{
constexpr double default_nice_running_check_interval = 10.0;
constexpr double minimum_running_check_interval = 0.1;
constexpr double maximum_running_check_interval = 60.0;
/// \}

/// \brief Bounds on the time interval between checks on GPU activity while the STORMM process is
///        waiting to resume after detecting a competing process, in wall clock seconds
/// \{
constexpr double default_nice_resume_check_interval = 60.0;
constexpr double maximum_resume_check_interval = 3600.0;
constexpr double minimum_resume_check_interval = 60.0;
/// \}

/// \brief Object to encapsulate job control and restraint information.  Like other namelist
///        encapsualtors, this object can take input file data as part of its construction, or
///        by a series of setters.  Validation of each piece of data is handled as it appears
///        either in the contructor or via setters.  Getter functions dispense the internal
///        information to any application using STORMM libraries.
class NiceControls {
public:

  /// \brief The constructor can prepare an object with default settings or read the corresponding
  ///        namelist to accept user input.
  ///
  /// \param tf          Input file translated into RAM
  /// \param start_line  Line of the input file to begin searching for the &solvent namelist
  /// \param found_nml   Indicator that the namelist was found in the input file
  /// \param policy_in   Requested error handling behavior
  /// \param wrap        Indicate that the search for a &dynamics namelist should carry on from
  ///                    the beginning of an input file if no such namelist is found starting
  ///                    from the original starting point
  /// \{
  NiceControls(ExceptionResponse policy_in = ExceptionResponse::DIE);

  NiceControls(const TextFile &tf, int *start_line, bool *found_nml,
               ExceptionResponse policy_in = ExceptionResponse::DIE,
               WrapTextSearch wrap = WrapTextSearch::NO);
  /// \}

  /// \brief As with other control objects, copy and move constructors, plus copy and move
  ///        assignment operators, can all take their default forms.
  /// \{
  NiceControls(const NiceControls &original) = default;
  NiceControls(NiceControls &&original) = default;
  NiceControls& operator=(const NiceControls &original) = default;
  NiceControls& operator=(NiceControls &&original) = default;
  /// \}

  /// \brief Get the workday start time as string recording hours:minutes:00 (workdays are rounded
  //         to the minute).
  std::string getWorkdayStart() const;

  /// \brief Get the workday end time as a real value in seconds after midnight.
  std::string getWorkdayEnd() const;

  /// \brief Get a list of the working days of the week, codified into a string
  std::string getWorkWeek() const;
  
  /// \brief Get a verdict on whether the workday is happening now.
  bool isWorkTimeNow() const;

  /// \brief Set the GPU effort level to undertake during working hours.
  ///
  /// \param workday_effort_level_in  The proportion of the GPU to utilize
  void setGpuEffortLevel(double workday_effort_level_in);

  /// \brief Set the interval at which to check for competing processes while the STORMM process
  ///        is in operation.
  ///
  /// \param running_check_interval_in  The chosen interval in wall clock seconds
  void setRunningCheckInterval(double running_check_interval_in);

  /// \brief Set the interval at which to check for competing processes after pausing the STORMM
  ///        process..
  ///
  /// \param resume_check_interval_in  The chosen interval in wall clock seconds
  void setResumeCheckInterval(double resume_check_interval_in);

  /// \brief Set the start of the workday.  The string will be interpreted first as an hour, or
  ///        as an hour and a minute if there is a ':' character separating two numbers.  Modifiers
  ///        such as 'am' and 'pm' will also be detected and handled appropriately.
  void setWorkTimeStart(const std::string &work_time_start_in);

  /// \brief Set the end of the workday.  The string will be interpreted as it is in the function
  ///        setWorkTimeStart(), above.
  void setWorkTimeEnd(const std::string &work_time_end_in);

private:

  ExceptionResponse policy;             ///< Set the behavior when bad inputs are encountered.
  double workday_effort_level;          ///< Set the level of GPU effort (out of the total, 1.0)
                                        ///<   to call for during working hours of the day
  double running_check_interval;        ///< The interval at which to check for other, competing
                                        ///<   GPU activity when running the STORMM process, in
                                        ///<   wall clock seconds
  double resume_check_interval;         ///< The frequency at which to check for other, competing
                                        ///<   GPU activity after pausing the STORMM process, in
                                        ///<   wall clock seconds
  bool stash_gpu_data;                  ///< Indicate that the STORMM process should not only
                                        ///<   pause but also create additional arrays on the CPU
                                        ///<   host to temporarily hold data from the GPU device
                                        ///<   allocations of certain large objects.
  int worktime_start;                   ///< The start of the workday, expressed in minutes after
                                        ///<   midnight
  int worktime_end;                     ///< The end of the workday, expressed in minutes after
                                        ///<   midnight
  std::vector<DayOfTheWeek> work_days;  ///< List of working days of the week
  std::vector<bool> working_day_mask;   ///< Mask for all work days of the week, for easy lookup
  
  /// Store a deep copy of the original namelist emulator as read from the input file.
  NamelistEmulator nml_transcript;

  /// \brief Validate a level of GPU effort during working hours.
  ///
  /// \param workday_effort_level_in  The chosen level of GPU effort
  bool validateGpuEffortLevel(const double workday_effort_level_in) const;

  /// \brief Validate the interval at which to check for competing GPU usage while running the
  ///        STORMM process (at any level of GPU utilization).
  ///
  /// \param running_check_interval_in  The chosen interval
  bool validateRunningCheckInterval(const double running_check_interval_in) const;

  /// \brief Validate the interval at which to check for competing GPU usage after pausing the
  ///        STORMM process to prioritize something competing with it.
  ///
  /// \param resume_check_interval_in  The chosen interval
  bool validateResumeCheckInterval(const double resume_check_interval_in) const;
  
  /// \brief Translate a string, e.g. 10:05PM, to a time in hours and minutes.  Return the number
  ///        of minutes after the previous midnight.
  ///
  /// \param time_statement  The string to translate
  int translateTimeString(const std::string &time_statement) const;

  /// \brief Translate one of the integers tored internally in this class's objects into a time
  ///        of day (appended with AM / PM, in "civilian" time).
  ///
  /// \param time_value  The integer to translate
  std::string translateTimeInteger(int time_value) const;
};

/// \brief Produce a namelist for specifying job control directives.  This feature provides STORMM
///        with an ability to run at reduced levels during working hours or even to pause itself
///        while another user or process is attempting to utilize the GPU.  This will allow users
///        to run STORMM on desktop resources without interfering regular work.
///
/// \param tf          Input text file to scan immediately after the namelist has been created
/// \param start_line  Line at which to begin scanning the input file for the namelist (this
///                    function will not wrap back to the beginning of the TextFile object, as the
///                    &rst namelist is intended to be repeatable)
/// \param found       Indicator that the namelist was present in the input file
/// \param policy      Reaction to exceptions encountered during namelist reading
/// \param wrap        Indicate that the search for an &dynamics namelist should carry on from the
///                    beginning of an input file if no such namelist is found starting from the
///                    original starting point 
NamelistEmulator niceInput(const TextFile &tf, int *start_line, bool *found,
                           ExceptionResponse policy = ExceptionResponse::DIE,
                           WrapTextSearch wrap = WrapTextSearch::NO);

} // namespace namelist
} // namespace stormm

#endif
