// -*-c++-*-
#ifndef OMNI_PARSE_H
#define OMNI_PARSE_H

#include <vector>
#include <string>
#include "Constants/behavior.h"
#include "DataTypes/common_types.h"
#include "DataTypes/omni_vector_types.h"
#include "polynumeric.h"
#include "textfile.h"
#include "textguard.h"

namespace omni {
namespace parse {

using constants::ExceptionResponse;
using constants::CaseSensitivity;

/// \brief Specify the way in which to print a standard real number or integer (with or without
///        leading zeros)
enum class NumberPrintStyle {
  STANDARD,      ///< Print without leading zeros
  LEADING_ZEROS  ///< Print with leading zeros
};

/// \brief Enumerate the types of wildcard characters.  The exact nature of the characters,
///        i.e. "*" or ".", may be defined elsewhere.
enum class WildCardKind {
  NONE,            ///< Not a wildcard: direct match is the only way to success
  FREE_CHARACTER,  ///< Matches any one character, but must match a character
  FREE_STRETCH     ///< Matches any number of consecutive characters, including no characters
};

/// \brief Constant to define the default formatting for integer and real number representations.
///        If these values are seen then the program knows to do free formatting.
constexpr int free_number_format = -32788;

/// \brief Overload the + operator to extend a string with a char4 object.  The x, y, z, and w
///        components of the char4 will continue to extend the string until one of them is the
///        terminating null character.
///
/// \param lhs  The string to extend
/// \param rhs  The char tuple to add
std::string operator+(const std::string &lhs, const char4 rhs);

/// \brief Convert a char4 to a string.  This is slightly faster than creating an empty string
///        and then using the += operator.  As above, the x, y, z, and w components of the char4
///        will continue to extend the string unless a null character is reached.
///
/// \param value  The char4 to convert
std::string char4ToString(const char4 value);

/// \brief Convert the first four characters of a string to char4, ignoring subsequent characters
///        and padding with whitespace if the string is shorter than four characters.
///
/// \param value  The string to convert
char4 stringToChar4(const std::string &value);

/// \brief Overload == to handle two char4 tuples
///
/// \param lhs  The first value to compare
/// \param rhs  The second value to compare
bool operator==(const char4 lhs, const char4 rhs);

/// \brief Determine whether a character string can qualify as an integer, real, or char4
///
/// \param a           The character string
/// \param cform       The expected format of the number
/// \param read_begin  Point in the string at which to being reading the number (defaults to the
///                    start of the string)
/// \param len         Expected length of the number (default 0, continues until the end of the
///                    string)
bool verifyNumberFormat(const char* a, NumberFormat cform, int read_begin = 0, int len = 0);

/// \brief Convert one or more characters to uppercase, if they are letters between 'a' and 'z'.
///        Otherwise, return the original character(s).
///
/// Overloaded:
///   - Convert a single character
///   - Convert a C++ string
///   - Convert a C-style string (given a limit), in place
///   - Convert a constant C-style string, returning a C++ string
///
/// \param tc      The character in question
/// \param tcs     The character string in question
/// \param n_char  The number of characters to convert in the C-style string (trusted that the
///                string has this many characters)
/// \{
char uppercase(const char tc);
std::string uppercase(const std::string &ts);
void uppercase(char* tcs, size_t n_char = 0);
std::string uppercase(const char* tcs);
/// \}

/// \brief Convert one or more characters to lowercase, if they are letters between 'a' and 'z'.
///        Otherwise, return the original character(s).
///
/// Overloaded:
///   - Convert a single character
///   - Convert a C++ string
///   - Convert a C-style string (given a limit), in place
///   - Convert a constant C-style string, returning a C++ string
///
/// \param tc   The character in question
/// \param tcs     The character string in question
/// \param n_char  The number of characters to convert in the C-style string (trusted that the
///                string has this many characters)
/// \{
char lowercase(const char tc);
std::string lowercase(const std::string &ts);
void lowercase(char* tcs, size_t n_char = 0);
std::string lowercase(const char* tcs);
/// \}

/// \brief Match two strings based on a particular case sensitivity setting
///
/// Overloaded:
///   - Take two C-strings (the fundamental case)
///   - Take a C-string and a C++ string (in either order)
///   - Take two C++ strings
///
/// \param sa      The first string
/// \param sb      The second string
/// \param csen    Case sensitivity setting
/// \param length  The number of characters to compare (if 0 and a C++ string is provided, the
///                length of the string will be taken)
/// \{
bool strncmpCased(const char* sa, const char* sb, int length,
                  CaseSensitivity csen = CaseSensitivity::YES);

bool strncmpCased(const std::string &sa, const char* sb,
                  CaseSensitivity csen = CaseSensitivity::YES, int length = -1);

bool strncmpCased(const char* sa, const std::string &sb,
                  CaseSensitivity csen = CaseSensitivity::YES, int length = -1);

bool strncmpCased(const std::string &sa, const std::string &sb,
                  CaseSensitivity csen = CaseSensitivity::YES, int length = -1);
/// \}

/// \brief Match two strings, the second of which may contain wildcards.
///
/// \param target     The literal string, containing no wildcards, which must be matched
/// \param query      The match string, which may contain wildcards
/// \param wildcards  Array of indicators as to whether the indices in query can be treated as
///                   wildcards, and if so of what type
bool strcmpWildCard(const std::string &target, const std::string &query,
                    const std::vector<WildCardKind> &wildcards);

/// \brief Determine the appropriate number of decimal places in which to display a real number
///        without superfluous zeros.
///
/// \param value  The number to display
/// \param limit  The maximum number of decimal places to report (optional)
int realDecimalPlaces(double value, int limit = 10);

/// \brief Convert a real number to a formatted string.
///
/// \param value     The number to convert
/// \param format_a  First of up to two format specifiers.  If this is given and the second format
///                  specifier is blank, it refers to the number of decimal places to print for
///                  the number.  If the second format specifier is also present, this becomes the
///                  total width of the general format number.
/// \param format_b  If given, this will convert the first format specifier into the total width
///                  and then specify the number of digits after the decimal.
/// \param style     Whether to print leading zeros when there is more format space than digits
std::string realToString(double value, int format_a = free_number_format,
                         int format_b = free_number_format,
                         NumberFormat method = NumberFormat::STANDARD_REAL,
                         NumberPrintStyle = NumberPrintStyle::STANDARD);

/// Convert an integer to a formatted string.
///
/// \param value  The number to convert
/// \param width  The format specifier, showing how many characters to give the integer
/// \param style  Whether to print leading zeros when there is more format space than digits
std::string intToString(llint value, int width = free_number_format,
                        NumberPrintStyle = NumberPrintStyle::STANDARD);

/// \brief Overload the + operator to make it easier to append sequences of TextGuard objects in
///        std::vectors.
const std::vector<TextGuard> operator+(const std::vector<TextGuard> &lhs,
                                       const std::vector<TextGuard> &rhs);

/// \brief Decide whether the sequence of characters in a TextFile beginning at a given position
///        completes a guarded section of the text, given the guard sequence that started the
///        section.
///
/// Overloaded:
///   - Detect guards within a text file
///   - Detect guards within a character array (C-string, must have sufficient length)
///
/// \param tfr        Text file reading abstract
/// \param line_idx   Index at which to seek the putative termination sequence
/// \param pos_idx    Index at which to seek the putative termination sequence
/// \param guard_seq  Text guard containing the termination sequence
/// \{
bool detectGuard(const TextFile::Reader &tfr, const int line_idx, const int pos_idx,
                 const std::string &guard_seq);
bool detectGuard(const char* line, const int pos_idx, const std::string &guard_seq);
/// \}

/// \brief Test all left-hand guard sequences to determine whether one of them matches the present
///        text sequence.
///
/// \param tfr       Text file reading abstract
/// \param line_idx  Index at which to seek the putative termination sequence
/// \param pos_idx   Index at which to seek the putative termination sequence
/// \param markers   List of all text guards
int applyGuard(const TextFile::Reader &tfr, const int line_idx, const int pos_idx,
               const std::vector<TextGuard> &markers);

/// \brief Apply escape sequences to the text.  Escape sequences bind before any other guards such
///        as comments or quotations.
///
/// Overloaded:
///   - Mark all escaped characters (not escape sequences, but escaped characters) in a TextFile
///   - Mark all escaped characters in a string
///
/// \param tfr      Text file reading abstract
/// \param escapes  Vector of escape sequences
/// \{
std::vector<bool> markEscapedCharacters(const TextFile::Reader &tfr,
                                        const std::vector<TextGuard> &escapes);

std::vector<bool> markEscapedCharacters(const char* textstr, const int n_char,
                                        const std::vector<TextGuard> &escapes);
/// \}

/// \brief Scan a text file and mark off text that is in quotations.
///
/// \param tfr           Read-only pointers to the text file of interest once it is read into RAM
/// \param markers       List of markers denoting guards for the text of interest
/// \param alternatives  List of markers denoting guards for other text, binding at the same level
///                      as markers
std::vector<bool> markGuardedText(const TextFile::Reader &tfr,
                                  const std::vector<TextGuard> &markers,
                                  const std::vector<TextGuard> &alternatives,
                                  const std::vector<TextGuard> &escapes = { TextGuard("\\") });

/// \brief Bearing some similarities to markGuardedText, this function operates on strings only,
///        marking each guarded segment as one level "deeper" than the one enclosing it.  The scope
///        depth starts at zero, incrementing with each opening guard encountered and decrementing
///        when it reaches a relevant closing guard for the most recent opening guard.
///
/// Overloaded:
///   - Determine scopes in a C-style string, with the option to specify length
///   - Determine scopes in a C++ string
///
/// \param input_text    The text to parse for scope depth
/// \param scope_guards  Scope definitions, all must have closing guards.  Line limits are
///                      irrelevant as the sbstrate is merely a string.  The input text can include
///                      Quoted multi-line strings from input files, having been 
/// \{
std::vector<int> resolveScopes(const char* input_text, int length,
                               const std::vector<TextGuard> scope_guards,
                               ExceptionResponse policy = ExceptionResponse::DIE);

std::vector<int> resolveScopes(const std::string &input_text,
                               const std::vector<TextGuard> scope_guards,
                               ExceptionResponse policy = ExceptionResponse::DIE);
/// \}

/// \brief A quick diagnostic to assess whether a string contains characters that might separate or
///        otherwise annotate text in special ways.
///
/// \param text    The text to analyze
/// \param delims  A list of delimiters (one character each)
int countDelimiters(const std::string &text, const std::vector<char> &delms);

/// \brief Separate a sequece of ascii characters into words, using whitespace as the delimiter.
///
/// Overloaded:
///   - Operate on a C-string
///   - Operate on a C++ string
///   - Operate on a C++ string with custom delimiters and no quotations or comments
///   - Operate on a TextFile object
///
/// \param text            The text to parse
/// \param n_char          Number of characters in the text
/// \param comment_mask    A mask of booleans indicating whether text is protected by comment
///                        symbols (optional).  Commented text will be ignored and count as a
///                        delimiter when parsing text into strings.
/// \param quotation_mask  A mask of booleans indicating whether text is protected by quotation
///                        marks (optional).  Quoter text will count as a single string in the
///                        final output.
/// \param quote_marks     A list of tex tguards that represent quotes.  All are given as SINGLE-
///                        line guards as the string is interpreted as only a single line (even if
///                        it contains carriage returns--the multi-line detail pertains to TextFile
///                        objects, where lines are separated by an auxiliary array of limits and
///                        there are no carrieage returns).  The quote marks are needed in order to
///                        shave quote characters off of any substrings that come from stretches of
///                        quoted text.
/// \param delimiters      List of extra delimiters (common choices would include "=" and ",")
/// \{
std::vector<std::string>
separateText(const char* text, int n_char, const std::vector<bool> &comment_mask = {},
             const std::vector<bool> &quotation_mask = {},
             const std::vector<TextGuard> &quote_marks = { TextGuard("\"", "\""),
                                                           TextGuard("'", "'") },
             const std::vector<std::string> &delimiters = {},
             const std::vector<TextGuard> &escapes = { TextGuard("\\") });

std::vector<std::string>
separateText(const std::string &text, const std::vector<bool> &comment_mask = {},
             const std::vector<bool> &quotation_mask = {},
             const std::vector<TextGuard> &quote_marks = { TextGuard("\"", "\""),
                                                           TextGuard("'", "'") },
             const std::vector<std::string> &delimiters = {},
             const std::vector<TextGuard> &escapes = { TextGuard("\\") });

std::vector<std::string>
separateText(const std::string &text, const std::vector<std::string> &delimiters,
             const std::vector<TextGuard> &escapes = { TextGuard("\\") });

std::vector<std::string>
separateText(const TextFile &text, const std::vector<bool> &comment_mask = {},
             const std::vector<bool> &quotation_mask = {},
             const std::vector<TextGuard> &quote_marks = { TextGuard("\"", "\""),
                                                           TextGuard("'", "'") },
             const std::vector<std::string> &delimiters = {},
             const std::vector<TextGuard> &escapes = { TextGuard("\\") });
/// \}

/// \brief Find one string within a vector of other strings.  Return the index, or the length of
///        the input vector (akind to std::find) if the index is not found.
///
/// \param vec  The vector of strings to search
/// \param str  The string to search for
int findStringInVector(const std::vector<std::string> &vec, const std::string &query);

/// \brief Convert a vector of strings into a vector of integers, checking each for validity and
///        reporting the result or substituting a zero if it looks suspicious.
///
/// \param sv      Vector of strings to convert
/// \param policy  The course of action if a malformed number is found
std::vector<int> vectorStrtol(const std::vector<std::string> &sv,
                              ExceptionResponse policy = ExceptionResponse::DIE);

/// \brief Convert a vector of strings into a vector of double-precision reals, checking each for
///        validity and reporting the result or substituting a zero if it looks suspicious.
///
/// \param sv      Vector of strings to convert
/// \param policy  The course of action if a malformed number is found
std::vector<double> vectorStrtod(const std::vector<std::string> &sv,
                                 ExceptionResponse policy = ExceptionResponse::DIE);
  
} // namespace parse
} // namespace omni

#endif
