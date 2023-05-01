// -*-c++-*-
#ifndef STORMM_REPORTING_ENUMERATORS_H
#define STORMM_REPORTING_ENUMERATORS_H

#include <string>
#include <vector>
#include "copyright.h"
#include "Constants/behavior.h"

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

/// \brief Produce strings detailing each of the enumerations above.
///
/// \param input  The enumeration to describe
std::string getEnumerationName(HelpSignalKind input);

} // namespace display

namespace review {

using constants::ExceptionResponse;

/// \brief Options for reporting results from one or more systems
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
  
/// \brief Possible formats for organized numerical output
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

/// \brief Possible components of each output report file section
enum class SectionComponent {
  NARRATIVE,  ///< A block of narration suitable for processing by the
  LIST,       ///< A bulleted or ordered list within an OrderedList object (see ordered_list.h)
  TABLE,      ///< A formatted table, referencing the ReportTable object (see report_table.h)
  SCRIPT      ///< Additional script to print to the output file with no comment protection
};

/// \brief Possible methods of numbering or identifying separate items in an ordered list
enum class ListEnumeration {
  BULLET,      ///< Use a bulleted list, with every item taking the same single-character symbol
  NUMBERED,    ///< Number each item and place a ")" after the marker
  ALPHABETIC,  ///< Assign a (lowercase) letter to each item and place a ")" after the marker.  If
               ///<   the number of items rises above 26, the next items will have symbol aa, ab,
               ///<   (...).  This implies that there are 26 letters in the fastest incrementing
               ///<   position and 27 (add blank space) in the second, but the way to calculate the
               ///<   number of available combinations is as 26 + 26^2 + 26^3 + ...
  ROMAN,       ///< Assing a (lowercase) Roman numeral to each list item.
  NONE         ///< No enumeration symbol is used.
};

/// \brief Options for terminating a segment of text, such as a pargraph.  In some cases, a
///        new line is desirable, but in others the formatted text should just stop on its last
///        character.
enum class TextEnds {
  AS_IS,   ///< Stop printing with the last submitted character.
  NEWLINE  ///< Ensure that any line containing characters other than white space and a protective
           ///<   marker ends with a carriage return.
};

/// \brief Differtiate the types of ReportTable content with a specific enumeration.
enum class TableContentKind {
  INTEGER,  ///< The table contains data that can be read as integers (does not distinguish 32-bit
            ///<   from 64-bit integers)
  REAL,     ///< The table contains data that can be read as real numbers
  STRING    ///< The table data can only be interpreted as strings, possibly char4
};

/// \brief An array of Roman numerals, queried when using the ListEnumeration ROMAN setting.
const std::vector<std::string> roman_numerals = { "i", "ii", "iii", "iv", "v", "vi", "vii", "viii",
                                                  "ix", "x", "xi", "xii", "xiii", "xiv", "xv",
                                                  "xvi", "xvii", "xviii", "xix", "xx", "xxi",
                                                  "xxii", "xxiii", "xxiv", "xv", " " };

/// \brief The maximum value that STORMM will convert to Roman numerals.  Any integer greater than
///        or equal to this will produce an error.
constexpr int maximum_roman_numeral = 25;

/// \brief Produce strings detailing each of the enumerations above.
///
/// \param input  The enumeration to describe
/// \{
std::string getEnumerationName(OutputScope input);
std::string getEnumerationName(OutputSyntax input);
std::string getEnumerationName(SectionComponent input);
std::string getEnumerationName(ListEnumeration input);
std::string getEnumerationName(TextEnds input);
std::string getEnumerationName(TableContentKind input);
/// \}

/// \brief Perform a bounds check on an integer to see if it can be represented with one of the
///        hard-wired Roman numerals in the roman_numerals array above, then convert it if
///        possible.  Produces an error if the integer cannot be converted, or returns a single
///        whitespace character (the final entry of the array) with any degree of fault
///        tolerance.
///
/// \param x       The integer to validate
/// \param policy  The degree of fault tolerance
const std::string& toRoman(int x, ExceptionResponse policy = ExceptionResponse::DIE);
  
} // namespace review
} // namespace stormm

#endif
