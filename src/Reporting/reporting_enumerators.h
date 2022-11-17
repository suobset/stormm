// -*-c++-*-
#ifndef STORMM_REPORTING_ENUMERATORS_H
#define STORMM_REPORTING_ENUMERATORS_H

#include <string>
#include "copyright.h"

namespace stormm {
namespace display {

/// \brief Possible modes in which command line program invocations can call for help messages
enum class HelpSignalKind {
  NO_ARGS,       ///< A command line call with no arguments will trigger a help message display
  NO_ARGS_ONLY,  ///< Only a call with no arguments, not a keyword, can trigger a help message
  KEYWORD,       ///< A keyword, by default any case-insensitive version of "help", "-help", or
                 ///<   or "--help" will trigger a help message display
  KEYWORD_ONLY   ///< Help messages are triggered only by an explicit keyword on the command line
};

} // namespace display

namespace review {

/// \brief Possible formats for organized numerical output.
enum class OutputSyntax {
  MATPLOTLIB,  ///< Results will be written to a script which the Python interpreter can run to
               ///<   display graphical results.  The script contains functions to simplify
               ///<   interactive plotting and rendering.
  MATRIX_PKG,  ///< A versatile format that works well with a number of matrix algebra packages,
               ///<   interspersing the data with commented blocks which serve the same purpose as
               ///<   narration in a typical output file.
  STANDALONE   ///< Results will be written to a unique format more typical of MD program output,
               ///<   with no syntax or other commands associated with plotting programs
};

/// \brief Options for reporting results from one or more systems.
enum class OutputScope {
  AVERAGES,          ///< Report the average energy of all systems, with standard deviations to
                     ///<   impart error bars.
  OUTLIERS,          ///< Report the average energy of all systems, with error bars, and individual
                     ///<   results for up to three outlier systems
  CLUSTER_AVERAGES,  ///< Report the average energy and error bars of systems governed by each
                     ///<   unique topology
  CLUSTER_OUTLIERS,  ///< Report the average energy and error bars of systems governed by each
                     ///<   unique topology, plus individual results for at most one outlier under
                     ///<   each topology
  FULL               ///< Report complete energies for all systems
};
  
/// \brief Produce strings detailing each of the enumerations above.
///
/// \param input  The enumeration to describe
/// \{
std::string getEnumerationName(OutputSyntax input);
std::string getEnumerationName(OutputScope input);
/// \}
  
} // namespace review
} // namespace stormm

#endif
