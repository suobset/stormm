#include <cmath>
#include <cstring>
#include <vector>
#include "copyright.h"
#include "Reporting/error_format.h"
#include "ascii_numbers.h"
#include "parse.h"

namespace stormm {
namespace parse {

//-------------------------------------------------------------------------------------------------
std::string operator+(const std::string &lhs, const char4 rhs) {
  std::string result = lhs;
  if (rhs.x == '\0') return result;
  result += rhs.x;
  if (rhs.y == '\0') return result;
  result += rhs.y;
  if (rhs.z == '\0') return result;
  result += rhs.z;
  if (rhs.w == '\0') return result;
  result += rhs.w;
  return result;
}

//-------------------------------------------------------------------------------------------------
std::string char4ToString(const char4 value) {
  std::string result;
  if (value.x == '\0') return result;
  result += value.x;
  if (value.y == '\0') return result;
  result += value.y;
  if (value.z == '\0') return result;
  result += value.z;
  if (value.w == '\0') return result;
  result += value.w;
  return result;
}

//-------------------------------------------------------------------------------------------------
char4 stringToChar4(const std::string &value) {
  char4 result;
  const int vlen = value.size();
  result.x = (vlen > 0) ? value[0] : ' ';
  result.y = (vlen > 1) ? value[1] : ' ';
  result.z = (vlen > 2) ? value[2] : ' ';
  result.w = (vlen > 3) ? value[3] : ' ';
  return result;
}

//-------------------------------------------------------------------------------------------------
bool operator==(const char4 lhs, const char4 rhs) {
  return (lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z && lhs.w == rhs.w);
}

//-------------------------------------------------------------------------------------------------
bool operator!=(const char4 lhs, const char4 rhs) {
  return (lhs.x != rhs.x || lhs.y != rhs.y || lhs.z != rhs.z || lhs.w != rhs.w);
}

//-------------------------------------------------------------------------------------------------
bool verifyNumberFormat(const char* a, const NumberFormat cform, const int read_begin,
                        const int len) {

  // Pre-allocate a buffer and fill it with string termination characters
  const int width = (len > 0) ? len : strlen(a) - read_begin;
  std::vector<char> buffer(width, ' ');

  // The CHAR4 type can accept any characters
  if (cform == NumberFormat::CHAR4) {
    return true;
  }
  
  // Fill the buffer with data from the input character array
  bool problem = false;
  bool number_begins = false;
  bool number_ends = false;
  for (int i = 0; i < width; i++) {
    buffer[i] = a[read_begin + i];
    number_begins = (number_begins || buffer[i] != ' ');
    number_ends = (number_ends || (number_begins && buffer[i] == ' '));
    problem = (problem || (buffer[i] != ' ' && number_ends));
  }

  // Check each character
  bool e_found = false;
  bool dot_found = false;
  int signs_found = 0;
  int minus_found = 0;
  for (int i = 0; i < width; i++) {
    if (buffer[i] == 'E' || buffer[i] == 'e') {
      problem = (problem || e_found);
      e_found = true;
    }
    else if (buffer[i] == '.') {
      problem = (problem || dot_found || e_found);
      dot_found = true;
    }
    else if (buffer[i] == '+' || buffer[i] == '-') {
      problem = (problem || (i > 0 && buffer[i - 1] != ' ' && e_found == false));
      minus_found += (buffer[i] == '-');
      signs_found++;
    }
    else {
      switch(cform) {
      case NumberFormat::SCIENTIFIC:
      case NumberFormat::STANDARD_REAL:
      case NumberFormat::INTEGER:
      case NumberFormat::LONG_LONG_INTEGER:
      case NumberFormat::UNSIGNED_INTEGER:
      case NumberFormat::UNSIGNED_LONG_LONG_INTEGER:
        problem = (problem || ((buffer[i] < '0' || buffer[i] > '9') && buffer[i] != ' '));
        break;
      case NumberFormat::CHAR4:
        break;
      }
    }
  }
  problem = (problem || number_begins == false);

  // Note any malformed numbers
  switch(cform) {
  case NumberFormat::SCIENTIFIC:
    problem = (problem || (e_found == false || dot_found == false || signs_found > 2));
    break;
  case NumberFormat::STANDARD_REAL:
    problem = (problem || (signs_found > 1 || e_found));
    break;
  case NumberFormat::INTEGER:
  case NumberFormat::LONG_LONG_INTEGER:
    problem = (problem || (signs_found > 1 || dot_found || e_found));
    break;
  case NumberFormat::UNSIGNED_INTEGER:
  case NumberFormat::UNSIGNED_LONG_LONG_INTEGER:
    problem = (problem || (signs_found > 1 || minus_found > 0 || dot_found || e_found));
    break;
  case NumberFormat::CHAR4:
    problem = (problem || width != 4);
    break;
  }

  return (problem == false);
}

//-------------------------------------------------------------------------------------------------
char uppercase(const char tc) {
  return tc - (tc >= 97 && tc <= 122) * 32;
}

//-------------------------------------------------------------------------------------------------
std::string uppercase(const std::string &ts) {
  std::string result = ts;
  const size_t n_char = ts.size();
  for (size_t i = 0; i < n_char; i++) {
    const char tc = result[i];
    result[i] = tc - (tc >= 97 && tc <= 122) * 32;
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
void uppercase(char* tcs, const size_t n_char) {
  const size_t scan_char = (n_char == 0) ? strlen(tcs) : n_char;
  for (size_t i = 0; i < scan_char; i++) {
    const char tc = tcs[i];
    tcs[i] = tc - (tc >= 97 && tc <= 122) * 32;
  }
}

//-------------------------------------------------------------------------------------------------
std::string uppercase(const char* tcs) {
  std::string result(tcs);
  const size_t n_char = result.size();
  for (size_t i = 0; i < n_char; i++) {
    const char tc = tcs[i];
    result[i] = tc - (tc >= 97 && tc <= 122) * 32;
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
char lowercase(const char tc) {
  return tc + (tc >= 65 && tc <= 90) * 32;
}

//-------------------------------------------------------------------------------------------------
std::string lowercase(const std::string &ts) {
  std::string result = ts;
  const size_t n_char = ts.size();
  for (size_t i = 0; i < n_char; i++) {
    const char tc = result[i];
    result[i] = tc + (tc >= 65 && tc <= 90) * 32;
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
void lowercase(char* tcs, const size_t n_char) {
  const size_t scan_char = (n_char == 0) ? strlen(tcs) : n_char;
  for (size_t i = 0; i < scan_char; i++) {
    const char tc = tcs[i];
    tcs[i] = tc + (tc >= 65 && tc <= 90) * 32;
  }
}

//-------------------------------------------------------------------------------------------------
std::string lowercase(const char* tcs) {
  std::string result(tcs);
  const size_t n_char = result.size();
  for (size_t i = 0; i < n_char; i++) {
    const char tc = tcs[i];
    result[i] = tc + (tc >= 65 && tc <= 90) * 32;
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
bool strcmpCased(const char* sa, const char* sb, const CaseSensitivity csen) {
  int i = 0;
  switch (csen) {
  case CaseSensitivity::YES:
    while (sa[i] == sb[i]) {
      if (sa[i] == '\0') {
        return true;
      }
      i++;
    }
    return false;
  case CaseSensitivity::NO:
    while (uppercase(sa[i]) == uppercase(sb[i])) {
      if (sa[i] == '\0') {
        return true;
      }
      i++;
    }
    return false;
  case CaseSensitivity::AUTOMATIC:
    rtErr("No AUTOMATIC behavior is defined for case-based string comparison.  AUTOMATIC "
          "settings for case sensitivity are defined at higher levels for specific situations, "
          "not the low-level implementation.", "strcmpCased");
    break;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
bool strcmpCased(const std::string &sa, const char* sb, const CaseSensitivity csen) {
  return strcmpCased(sa.c_str(), sb, csen);
}

//-------------------------------------------------------------------------------------------------
bool strcmpCased(const char* sa, const std::string &sb, const CaseSensitivity csen) {
  return strcmpCased(sa, sb.c_str(), csen);
}

//-------------------------------------------------------------------------------------------------
bool strcmpCased(const std::string &sa, const std::string &sb, const CaseSensitivity csen) {
  return strcmpCased(sa.c_str(), sb.c_str(), csen);
}

//-------------------------------------------------------------------------------------------------
bool strncmpCased(const char* sa, const char* sb, const int length, const CaseSensitivity csen) {
  bool match = true;
  switch (csen) {
  case CaseSensitivity::YES:
    for (int i = 0; i < length; i++) {
      match = (match && (sa[i] == sb[i]));
      if (sa[i] == '\0' && match) {
        return true;
      }
    }
    break;
  case CaseSensitivity::NO:
    for (int i = 0; i < length; i++) {
      match = (match && (uppercase(sa[i]) == uppercase(sb[i])));
      if (sa[i] == '\0' && match) {
        return true;
      }
    }
    break;
  case CaseSensitivity::AUTOMATIC:
    rtErr("No AUTOMATIC behavior is defined for case-based string comparison.  AUTOMATIC "
          "settings for case sensitivity are defined at higher levels for specific situations, "
          "not the low-level implementation.", "strncmpCased");
    break;
  }
  return match;
}

//-------------------------------------------------------------------------------------------------
bool strncmpCased(const std::string &sa, const char* sb, const CaseSensitivity csen,
                  const int length) {
  const int actual_length = (length < 0) ? sa.size() : length;
  if (length > static_cast<int>(sa.size())) {
    return false;
  }
  return strncmpCased(sa.c_str(), sb, actual_length, csen);
}

//-------------------------------------------------------------------------------------------------
bool strncmpCased(const char* sa, const std::string &sb, const CaseSensitivity csen,
                  const int length) {
  const int actual_length = (length < 0) ? sb.size() : length;
  if (length > static_cast<int>(sb.size())) {
    return false;
  }
  return strncmpCased(sa, sb.c_str(), actual_length, csen);
}

//-------------------------------------------------------------------------------------------------
bool strncmpCased(const std::string &sa, const std::string &sb, const CaseSensitivity csen,
                  const int length) {
  if (length < 0 && sa.size() != sb.size()) {
    rtErr("The extent of a comparison cannot be inferred from two string arguments of different "
          "lengths (" + std::to_string(sa.size()) + " and " + std::to_string(sb.size()) + ", " +
          sa + " and " + sb + ").", "strncmpCased");
  }
  if (length > static_cast<int>(sa.size()) || length > static_cast<int>(sb.size())) {
    rtErr("String of length " + std::to_string(sa.size()) + " and " + std::to_string(sb.size()) +
          " cannot be subjected to a comparison of " + std::to_string(length) + " characters.",
          "strncmpCased");
  }
  const int actual_length = (length < 0) ? sa.size() : length;
  return strncmpCased(sa.c_str(), sb.c_str(), actual_length, csen);
}

//-------------------------------------------------------------------------------------------------
bool strcmpWildCard(const std::string &target, const std::string &query,
                    const std::vector<WildCardKind> &wildcards) {

  // Check that the wildcards array is of the appropriate length
  if (wildcards.size() != query.size()) {
    rtErr("A wildcard string comparison requires that the second string, with " +
          std::to_string(query.size()) + " characters, be supported by the same number of "
          "wildcards (currently " + std::to_string(wildcards.size()) + ").", "strncmpWildCard");
  }
  
  // Ignore white space at the beginning of the target
  int target_con = 0;
  const int target_len = target.size();
  while (target_con < target_len && target[target_con] == ' ') {
    target_con++;
  }
  int query_con = 0;
  const int query_len = query.size();
  bool star_power_buff = false;
  const int n_wild = wildcards.size();
  for (int j = 0; j < query_len; j++) {

    // White space in the query does nothing to advance matching in the target
    if (query[j] == ' ') {
      query_con++;
      continue;
    }
    switch (wildcards[j]) {
    case WildCardKind::NONE:
      if (star_power_buff) {
        while (target_con < target_len && query[j] != target[target_con]) {
          target_con++;
        }
        star_power_buff = false;
      }
      if (target_con < target_len && query[j] == target[target_con]) {
        query_con++;
        target_con++;
      }
      break;
    case WildCardKind::FREE_CHARACTER:
      if (target_con < target_len) {
        query_con++;
        target_con++;
      }
      break;
    case WildCardKind::FREE_STRETCH:
      query_con++;
      star_power_buff = true;
      break;
    }
  }

  // If the query ends with a free stretch wildcard active, keep it going until the end
  target_con = (star_power_buff) ? target_len : target_con;

  // Ignore white space at the end of the target
  while (target_con < target_len && target[target_con] == ' ') {
    target_con++;
  }
  
  return (target_con == target_len && query_con == query_len);
}

//-------------------------------------------------------------------------------------------------
int realDecimalPlaces(const double value, const int limit) {
  
  // Determine the fractional component
  double abs_val = fabs(value);
  double frac = abs_val - floor(abs_val);
  double smallest_significant_amount = pow(0.1, limit);
  int n_place = 0;
  while (n_place < limit && 1.01 * frac > smallest_significant_amount &&
         1.01 * fabs(1.0 - frac) > smallest_significant_amount) {
    frac *= 10.0;
    frac -= floor(frac);
    n_place++;
  }
  return n_place;
}

//-------------------------------------------------------------------------------------------------
std::string realToString(const double value, const int format_a, const int format_b,
                         const NumberFormat method, const NumberPrintStyle style) {

  // Check the overall format
  switch (method) {
  case NumberFormat::SCIENTIFIC:
    if (format_a != free_number_format && format_b != free_number_format &&
        format_b > format_a - 7) {
      rtErr("A real number of format %" + std::to_string(format_a) + "." +
            std::to_string(format_b) + "e has too many decimal places for its overall length.",
            "realToString");
    }
    break;
  case NumberFormat::STANDARD_REAL:
    if (format_a != free_number_format && format_b != free_number_format &&
        format_b > format_a - 2) {
      rtErr("A real number of format %" + std::to_string(format_a) + "." +
            std::to_string(format_b) + "f has too many decimal places for its overall length.",
            "realToString");
    }
    break;
  case NumberFormat::INTEGER:
  case NumberFormat::LONG_LONG_INTEGER:
  case NumberFormat::UNSIGNED_INTEGER:
  case NumberFormat::UNSIGNED_LONG_LONG_INTEGER:
  case NumberFormat::CHAR4:
    rtErr("The printing method for a real number must be either SCIENTIFIC or STANDARD_REAL.",
          "realToString");
  }

  // Examine the given dimensions
  if (format_b < 0 && format_b != free_number_format) {
    rtErr("A nonsensical number of decimal places (" + std::to_string(format_b) +
          ") was specified.", "realToString");
  }
  if (abs(format_a) >= 62 && format_a != free_number_format) {
    rtErr("The requested number exceeds format limits (maximum 64 characters, %" +
          std::to_string(format_a) + "." + std::to_string(format_b) + "f requested).",
          "realToString");
  }

  // Print the formatted number and return it as a string
  char buffer[64];
  if (format_a != free_number_format && format_b != free_number_format) {
    switch (style) {
    case NumberPrintStyle::STANDARD:
      if (method == NumberFormat::SCIENTIFIC) {
        sprintf(buffer, "%*.*e", format_a, format_b, value);
      }
      else {
        sprintf(buffer, "%*.*f", format_a, format_b, value);
      }
      break;
    case NumberPrintStyle::LEADING_ZEROS:
      if (method == NumberFormat::SCIENTIFIC) {
        sprintf(buffer, "%0*.*e", format_a, format_b, value);
      }
      else {
        sprintf(buffer, "%0*.*f", format_a, format_b, value);
      }
      break;
    }
  }
  else if (format_b != free_number_format) {
    switch (style) {
    case NumberPrintStyle::STANDARD:
      if (method == NumberFormat::SCIENTIFIC) {
        sprintf(buffer, "%.*e", format_b, value);
      }
      else {
        sprintf(buffer, "%.*f", format_b, value);
      }
      break;
    case NumberPrintStyle::LEADING_ZEROS:
      if (method == NumberFormat::SCIENTIFIC) {
        sprintf(buffer, "%0.*e", format_b, value);
      }
      else {
        sprintf(buffer, "%0.*f", format_b, value);
      }
      break;
    }
  }
  else if (format_a != free_number_format) {
    switch (style) {
    case NumberPrintStyle::STANDARD:
      if (method == NumberFormat::SCIENTIFIC) {
        sprintf(buffer, "%.*e", format_a, value);
      }
      else {
        sprintf(buffer, "%.*f", format_a, value);

      }
      break;
    case NumberPrintStyle::LEADING_ZEROS:
      if (method == NumberFormat::SCIENTIFIC) {
        sprintf(buffer, "%0.*e", format_a, value);
      }
      else {
        sprintf(buffer, "%0.*f", format_a, value);
      }
      break;
    }
  }
  else {
    sprintf(buffer, "%.*f", realDecimalPlaces(value, 8), value);
  }
  return std::string(buffer);
}

//-------------------------------------------------------------------------------------------------
std::string realToString(const double value, const int format_a, const NumberFormat method,
                         const NumberPrintStyle style) {
  return realToString(value, format_a, free_number_format, method, style);
}

//-------------------------------------------------------------------------------------------------
std:: string intToString(const llint value, const int width, const NumberPrintStyle style) {

  // Examine the given dimensions
  if (abs(width) >= 62) {
    rtErr("The requested number exceeds format limits (maximum 64 characters, " +
          std::to_string(width) + " requested).", "intToString");
  }

  char buffer[64];
  if (width != free_number_format) {
    switch (style) {
    case NumberPrintStyle::STANDARD:
      sprintf(buffer, "%*lld", width, value);
      break;
    case NumberPrintStyle::LEADING_ZEROS:
      sprintf(buffer, "%0*lld", width, value);
      break;
    }
  }
  else {
    sprintf(buffer, "%lld", value);
  }
  return std::string(buffer);
}

//-------------------------------------------------------------------------------------------------
const std::vector<TextGuard> operator+(const std::vector<TextGuard> &lhs,
                                       const std::vector<TextGuard> &rhs) {
  std::vector<TextGuard> result(lhs.begin(), lhs.end());
  result.insert(result.end(), rhs.begin(), rhs.end());
  return result;
}

//-------------------------------------------------------------------------------------------------
bool detectGuard(const TextFileReader &tfr, const int line_idx, const int pos_idx,
                 const std::string &guard_seq) {
  const int gsq_length = guard_seq.size();
  if (gsq_length == 0) {
    return false;
  }
  if (tfr.line_limits[line_idx] + pos_idx + gsq_length <= tfr.line_limits[line_idx + 1]) {
    bool matched = true;
    const int readbase = tfr.line_limits[line_idx] + pos_idx;
    for (int i = 0; i < gsq_length; i++) {
      matched = (matched && guard_seq[i] == tfr.text[readbase + i]);
    }
    return matched;
  }
  return false;
}

//-------------------------------------------------------------------------------------------------
bool detectGuard(const char* line, const int pos_idx, const std::string &guard_seq) {
  const int gsq_length = guard_seq.size();
  if (gsq_length == 0) {
    return false;
  }
  bool matched = true;
  for (int i = 0; i < gsq_length; i++) {
    matched = (matched && guard_seq[i] == line[pos_idx + i]);
  }
  return matched;
}

//-------------------------------------------------------------------------------------------------
int applyGuard(const TextFileReader &tfr, const int line_idx, const int pos_idx,
               const std::vector<TextGuard> &markers) {
  const int n_markers = markers.size();
  for (int i = 0; i < n_markers; i++) {
    if (detectGuard(tfr, line_idx, pos_idx, markers[i].getLeft())) {
      return i;
    }
  }
  return -1;
}

//-------------------------------------------------------------------------------------------------
std::vector<bool> markEscapedCharacters(const TextFileReader &tfr,
                                        const std::vector<TextGuard> &escapes) {
  std::vector<bool> result(tfr.line_limits[tfr.line_count], false);
  for (int i = 0; i < tfr.line_count; i++) {
    for (int j = tfr.line_limits[i]; j < tfr.line_limits[i + 1]; j++) {
      int guarded = applyGuard(tfr, i, j - tfr.line_limits[i], escapes);
      if (guarded >= 0 && j + escapes[guarded].leftSize() < tfr.line_limits[i + 1]) {
        result[j + escapes[guarded].leftSize()] = true;
      }
    }
  }

  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<bool> markEscapedCharacters(const char* textstr, const int n_char,
                                        const std::vector<TextGuard> &escapes) {
  const int n_escape = escapes.size();
  std::vector<bool> result(n_char, false);
  for (int i = 0; i < n_char; i++) {
    for (int j = 0; j < n_escape; j++) {
      if (detectGuard(textstr, i, escapes[j].getLeft())) {
        result[i + escapes[j].leftSize()] = true;
      }
    }
  }

  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<bool> markGuardedText(const TextFileReader &tfr,
                                  const std::vector<TextGuard> &markers,
                                  const std::vector<TextGuard> &alternatives,
                                  const std::vector<TextGuard> &escapes) {
  const int n_markers = markers.size();
  std::vector<bool> result(tfr.line_limits[tfr.line_count], false);

  // Make all characters affected by an escape sequence.  Escapes are TextGuards with only a
  // left-hand sequence (by default, '\') and affect only the character directly behind them.
  std::vector<bool> literals = markEscapedCharacters(tfr, escapes);

  // Search each line of the text file
  std::vector<TextGuard> all_markers(markers);
  const int n_all_markers = all_markers.size();
  all_markers.insert(all_markers.end(), alternatives.begin(), alternatives.end());
  int guarded = -1;
  for (int i = 0; i < tfr.line_count; i++) {
    for (int j = tfr.line_limits[i]; j < tfr.line_limits[i + 1]; j++) {

      // Process literals to have the same guarded state the escape sequences that preceded them
      if (literals[j]) {
        result[j] = result[j-1];
        continue;
      }

      // Check for any single-character text guard markers
      if (guarded >= 0) {
        if (detectGuard(tfr, i, j - tfr.line_limits[i], all_markers[guarded].getRight())) {
          if (guarded < n_markers) {

            // Characters in the guard sequence count as part of the guarded section of text,
            // and this section was guarded in one of the manners we want to detect.
            for (int k = 0; k < all_markers[guarded].rightSize(); k++) {
              result[j + k] = true;
            }
          }
          j += all_markers[guarded].rightSize() - 1;
          guarded = -1;
        }
        else {
          if (guarded < n_markers) {
            result[j] = true;
          }
        }
      }
      else {

        // Try detecting a right guard in the absence of a guarded state.  If the left guard
        // of that guard pair is the same as the right guard, it doesn't count (false alarm),
        // as it means a sequence is just beginning.
        for (int k = 0; k < n_all_markers; k++) {
          if (detectGuard(tfr, i, j - tfr.line_limits[i], all_markers[k].getRight()) &&
              detectGuard(tfr, i, j - tfr.line_limits[i], all_markers[k].getLeft()) == false) {
            rtErr("A mismatched right-hand text guard was detected at line " +
                  std::to_string(i + 1) + ", position " + std::to_string(j + 1) + " in file " +
                  tfr.file_name + ".", "markGuardedText");
          }
        }

        // Try detecting a guard
        guarded = applyGuard(tfr, i, j - tfr.line_limits[i], all_markers);
        if (guarded >= 0) {

          // A guarded sequence of one of the types we want to detect just began.  The left
          // delimiter is a part of the sequence.
          if (guarded < n_markers) {
            for (int k = 0; k < all_markers[guarded].leftSize(); k++) {
              result[j + k] = true;
            }
          }
          j += all_markers[guarded].leftSize();

          // If the guard extends to the end of the line, take it there and no further
          if (all_markers[guarded].rightSize() == 0) {
            if (guarded < n_markers) {
              for (int k = j; k < tfr.line_limits[i + 1]; k++) {
                result[k] = true;
              }
            }
            j = tfr.line_limits[i + 1];
            guarded = -1;
          }

          // Decrement j before the loop iteration increments it
          j--;
        }

      }
    }

    // If we reach the end of a line and remain in a guarded state, make sure that the guard
    // can span multiple lines.
    if (guarded >= 0 && all_markers[guarded].getTerminationRequirement() &&
        all_markers[guarded].getSpan() != LineSpan::MULTIPLE) {
      rtErr("Guard sequence (\"" + all_markers[guarded].getLeft() + "\", \"" +
            all_markers[guarded].getRight() + "\") begins on line " + std::to_string(i + 1) +
            " of file " + tfr.file_name + ".  It requires a terminating sequence and cannot span "
            "multiple lines.", "markGuardedText");
    }
  }

  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<int> resolveScopes(const char* input_text, int length,
                               const std::vector<TextGuard> scope_guards,
                               const ExceptionResponse policy) {
  const int n_char = (length == 0) ? strlen(input_text) : length;
  int level = 0;
  std::vector<int> result(n_char, level);
  const int n_guard = scope_guards.size();
  std::vector<int> guard_layers;
  int i = 0;
  while (i < n_char) {

    // Check for any opening guard
    bool guard_found = false;
    for (int j = 0; j < n_guard; j++) {
      if (detectGuard(input_text, i, scope_guards[j].getLeft())) {
        guard_layers.push_back(j);
        const int n_left = scope_guards[j].leftSize();
        for (int k = 0; k < n_left; k++) {
          result[i + k] = level;
        }
        level++;
        i += n_left;
        guard_found = true;
        break;
      }
    }
    if (guard_found) {
      continue;
    }

    // Check for any closing guard 
    for (int j = 0; j < n_guard; j++) {
      if (detectGuard(input_text, i, scope_guards[j].getRight())) {
        if (j == guard_layers.back()) {
          level--;
          guard_layers.pop_back();
          const int n_right = scope_guards[j].rightSize();
          for (int k = 0; k < n_right; k++) {
            result[i + k] = level;
          }
          i += n_right;
          guard_found = true;
          break;
        }
        else {
          rtWarn("Closing guard \"" + scope_guards[j].getRight() + "\" found when expecting \"" +
                 scope_guards[guard_layers.back()].getRight() + "\" to close \"" +
                 scope_guards[guard_layers.back()].getLeft() + "\".", "resolveScopes");
        }
      }
    }
    if (guard_found) {
      continue;
    }

    // Final possibility: mark the one character as being at the current level and increment
    // the counter by one.
    result[i] = level;
    i++;
  }

  // Check that all scopes have been resolved
  if (guard_layers.size() > 0) {
    std::string all_remaining_scopes;
    const int n_rem = guard_layers.size();
    for (int i = 0; i < n_rem; i++) {
      all_remaining_scopes += "\"" + scope_guards[guard_layers[i]].getLeft() + "\"";
      if (i < n_rem - 1) {
        all_remaining_scopes += ", ";
      }
    }
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("String " + std::string(input_text) + " leaves " +
            std::to_string(guard_layers.size()) + " scopes unresolved: " + all_remaining_scopes +
            ".", "resolvedScopes");
    case ExceptionResponse::WARN:
      rtWarn("String " + std::string(input_text) + " leaves " +
             std::to_string(guard_layers.size()) + " scopes unresolved: " + all_remaining_scopes +
             ".", "resolvedScopes");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
  }

  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<int> resolveScopes(const std::string &input_text,
                               const std::vector<TextGuard> scope_guards,
                               const ExceptionResponse policy) {
  return resolveScopes(input_text.c_str(), input_text.size(), scope_guards, policy);
}

//-------------------------------------------------------------------------------------------------
int countDelimiters(const std::string &text, const std::vector<char> &delms) {
  int nd = 0;
  const int n_char = text.size();
  const int n_delm = delms.size();
  for (int i = 0; i < n_char; i++) {
    const char tmpc = text[i];
    for (int j = 0; j < n_delm; j++) {
      nd += (tmpc == delms[j]);
    }
  }
  return nd;
}

//-------------------------------------------------------------------------------------------------
std::vector<std::string> separateText(const char* text, int n_char,
                                      const std::vector<bool> &comment_mask,
                                      const std::vector<bool> &quotation_mask,
                                      const std::vector<TextGuard> &quote_marks,
                                      const std::vector<std::string> &delimiters,
                                      const std::vector<TextGuard> &escapes) {

  // Check inputs
  const bool comments_exist = (comment_mask.size() > 0);
  const bool quotes_exist = (quotation_mask.size() > 0);
  const bool escapes_exist = (escapes.size() > 0);
  const int n_quote_guards = quote_marks.size();
  if (quotes_exist && n_quote_guards == 0) {
    rtErr("Separating text with masked quotations requires a list of quotation symbols.",
          "separateText");
  }
  if ((comments_exist && comment_mask.size() != n_char) ||
      (quotes_exist && quotation_mask.size() != n_char)) {
    std::string sample;
    for (int i = 0; i < std::min(n_char, 24); i++) {
      sample += text[i];
    }
    const std::string problem_area = (comments_exist && comment_mask.size() != n_char) ?
                                     "comment" : "quotation";
    rtErr("The " + problem_area + " mask for text beginning \"" + sample +
          "\" does not conform to the overall size of the text.", "separateText");
  }

  // Mark escaped characters
  const std::vector<bool> escape_mask = (escapes_exist) ?
                                        markEscapedCharacters(text, n_char, escapes) :
                                        std::vector<bool>(n_char, false);

  // Convert all additional delimiters to whitespace
  const int n_delimiters = delimiters.size();
  std::vector<bool> delimiter_mask(n_char, false);
  for (int i = 0; i < n_char; i++) {
    delimiter_mask[i] = (text[i] == ' ');
  }
  for (int i = 0; i < n_delimiters; i++) {
    const int delimiter_size = delimiters[i].size();
    if (delimiter_size == 1) {
      const char idel = delimiters[i][0];
      if (comments_exist && quotes_exist) {
        for (int j = 0; j < n_char; j++) {
          delimiter_mask[j] = ((comment_mask[j] == false && quotation_mask[j] == false &&
                                escape_mask[j] == false && text[j] == idel) || delimiter_mask[j]);
        }
      }
      else if (comments_exist) {
        for (int j = 0; j < n_char; j++) {
          delimiter_mask[j] = ((comment_mask[j] == false && escape_mask[j] == false &&
                                text[j] == idel) || delimiter_mask[j]);
        }
      }
      else if (quotes_exist) {
        for (int j = 0; j < n_char; j++) {
          delimiter_mask[j] = ((quotation_mask[j] == false && escape_mask[j] == false &&
                                text[j] == idel) || delimiter_mask[j]);
        }
      }
      else {
        for (int j = 0; j < n_char; j++) {
          delimiter_mask[j] = ((escape_mask[j] == false && text[j] == idel) || delimiter_mask[j]);
        }
      }
    }
    else {
      const char* idel = delimiters[i].c_str();
      for (int j = 0; j < n_char; j++) {
        const int dlim = j + delimiters[i].size();
        if (dlim >= n_char) {
          continue;
        }
        const char* idel = delimiters[i].c_str();
        bool mismatch = false;
        if (comments_exist && quotes_exist) {
          for (int k = j; k < dlim; k++) {
            mismatch = ((comment_mask[k] == false && quotation_mask[k] == false &&
                         escape_mask[k] == false && text[k] != idel[k - j]) || mismatch);
          }
        }
        else if (comments_exist) {
          for (int k = j; k < dlim; k++) {
            mismatch = ((comment_mask[k] == false && escape_mask[k] == false &&
                         text[k] != idel[k - j]) || mismatch);
          }
        }
        else if (quotes_exist) {
          for (int k = j; k < dlim; k++) {
            mismatch = ((quotation_mask[k] == false && escape_mask[k] == false &&
                         text[k] != idel[k - j]) || mismatch);
          }
        }
        else {
          for (int k = j; k < dlim; k++) {
            mismatch = ((escape_mask[k] == false && text[k] != idel[k - j]) || mismatch);
          }
        }
        if (mismatch == false) {
          for (int k = j; k < dlim; k++) {
            delimiter_mask[k] = true;
          }
        }
      }
    }
  }

  // With the delimiter mask, comments mask, quotations mask, and
  // escaped characters mask in hand, find distinct strings.
  std::string growing_word;
  std::vector<std::string> result;
  for (int i = 0; i < n_char; i++) {

    // As escaped character just gets added to the growing word string, preventing anything
    // else from happening.
    if (escape_mask[i]) {
      growing_word += text[i];
      continue;
    }

    // Continue until a non-delimiter, non-comment character is encountered
    if (delimiter_mask[i] || (comments_exist && comment_mask[i])) {
      if (growing_word.size() > 0) {
        result.push_back(growing_word);
        growing_word.resize(0);
      }
      continue;
    }

    // A quotation is one word, unless it's multiple quotes back-to-back
    if (quotes_exist) {
      int n_quoted_char = 0;
      while (i + n_quoted_char < n_char && quotation_mask[i + n_quoted_char]) {
        n_quoted_char++;
      }
      if (n_quoted_char > 0) {

        // Search all delimiters and determine which quotes are responsible for this string
        int guarded = -1;
        for (int j = i; j < i + n_quoted_char; j++) {

          // Again, escaped characters just get added and pass over any further processing.
          // The guarded state (some sort of quotation) does not change.
          if (escape_mask[j] && guarded >= 0) {
            growing_word += text[j];
            continue;
          }

          if (guarded == -1) {
            for (int k = 0; k < n_quote_guards; k++) {
              if (j + quote_marks[k].leftSize() < i + n_quoted_char &&
                  detectGuard(text, j, quote_marks[k].getLeft())) {
                guarded = k;
                j += quote_marks[k].leftSize();
                break;
              }
            }
          }

          // Some form of guard must be in place by now
          if (guarded < 0) {
            std::string sample;
            for (int k = 0; k < std::min(n_quoted_char, 240); k++) {
              sample += text[i + k];
            }
            rtErr("Quotation beginning <<" + sample + ">> does not conform to any of the quote "
                  "marks provided.", "separateText");
          }

          // Check for a right-hand quotation terminator
          if (j + quote_marks[guarded].rightSize() <= i + n_quoted_char &&
              detectGuard(text, j, quote_marks[guarded].getRight())) {
            if (growing_word.size() > 0) {
              result.push_back(growing_word);
              growing_word.resize(0);
            }
            j += quote_marks[guarded].rightSize() - 1;
            guarded = -1;
          }

          // Grow the current word if this is still under a quotation guard
          if (guarded >= 0) {
            growing_word += text[j];
          }
        }

        // The quotation should end with a closing mark
        if (guarded >= 0) {
          std::string sample;
          for (int k = 0; k < std::min(n_quoted_char, 24); k++) {
            sample += text[i + k];
          }
          rtErr("Quotation beginning <<" + sample + ">> terminates without using any of the quote "
                "marks provided.", "separateText");
        }

        // Increment the counter to reflect progress through this quoted stretch of text
        i += n_quoted_char - 1;
      }
    }
    if (delimiter_mask[i] == false && (comments_exist == false || comment_mask[i] == false) &&
        (quotes_exist == false || quotation_mask[i] == false)) {

      // Grow the current word
      growing_word += text[i];
    }
  }

  // Append the final word 
  if (growing_word.size() > 0) {
    result.push_back(growing_word);
  }

  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<std::string> separateText(const std::string &text,
                                      const std::vector<bool> &comment_mask,
                                      const std::vector<bool> &quotation_mask,
                                      const std::vector<TextGuard> &quote_marks,
                                      const std::vector<std::string> &delimiters,
                                      const std::vector<TextGuard> &escapes) {
  return separateText(text.c_str(), text.size(), comment_mask, quotation_mask, quote_marks,
                      delimiters, escapes);
}

//-------------------------------------------------------------------------------------------------
std::vector<std::string> separateText(const std::string &text,
                                      const std::vector<std::string> &delimiters,
                                      const std::vector<TextGuard> &escapes) {
  return separateText(text.c_str(), text.size(), std::vector<bool>(), std::vector<bool>(),
                      std::vector<TextGuard>(), delimiters, escapes);
}

//-------------------------------------------------------------------------------------------------
std::vector<std::string> separateText(const TextFile &tf,
                                      const std::vector<bool> &comment_mask,
                                      const std::vector<bool> &quotation_mask,
                                      const std::vector<TextGuard> &quote_marks,
                                      const std::vector<std::string> &delimiters,
                                      const std::vector<TextGuard> &escapes) {
  TextFileReader tfr = tf.data();
  return separateText(tfr.text, tfr.line_limits[tfr.line_count], comment_mask, quotation_mask,
                      quote_marks, delimiters, escapes);
}

//-------------------------------------------------------------------------------------------------
int findStringInVector(const std::vector<std::string> &vec, const std::string &query) {

  // Check for simple problems
  const size_t qsize = query.size();
  const int vsize = vec.size();
  if (vsize == 0) {
    return 0;
  }

  // Peek at the first few characters to see if a match is possible.
  // There is at least one character in the query string.  Do a deeper search if warranted.
  const char first_char = query[0];  
  for (int i = 0; i < vsize; i++) {
    if (qsize == 0LLU && vec[i].size() == qsize) {
      return i;
    }
    if (vec[i].size() == qsize && vec[i][0] == first_char && vec[i] == query) {
      return i;
    }
  }

  // The query was not found.  Return the length of the vector.
  return vsize;
}

//-------------------------------------------------------------------------------------------------
std::vector<int> vectorStrtol(const std::vector<std::string> &sv, const ExceptionResponse policy) {
  const int nval = sv.size();
  std::vector<int> result(nval, 0.0);
  for (int i = 0; i < nval; i++) {
    if (verifyNumberFormat(sv[i].c_str(), NumberFormat::INTEGER)) {
      result[i] = strtol(sv[i].c_str(), nullptr, 10);
    }
    else {
      switch (policy) {
      case ExceptionResponse::DIE:
        rtErr("The string " + sv[i] + " is not a valid integer.", "vectorStrtol");
      case ExceptionResponse::WARN:
        rtWarn("The string " + sv[i] + " is not a valid integer.  Substituting 0.",
               "vectorStrotol");
	break;
      case ExceptionResponse::SILENT:
	break;
      }
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<double> vectorStrtod(const std::vector<std::string> &sv,
                                 const ExceptionResponse policy) {
  const int nval = sv.size();
  std::vector<double> result(nval, 0.0);
  for (int i = 0; i < nval; i++) {
    if (verifyNumberFormat(sv[i].c_str(), NumberFormat::STANDARD_REAL) ||
        verifyNumberFormat(sv[i].c_str(), NumberFormat::SCIENTIFIC)) {
      result[i] = strtod(sv[i].c_str(), nullptr);
    }
    else {
      switch (policy) {
      case ExceptionResponse::DIE:
        rtErr("The string " + sv[i] + " is not a valid real number.", "vectorStrtod");
      case ExceptionResponse::WARN:
        rtWarn("The string " + sv[i] + " is not a valid real number.  Substituting 0.0.",
               "vectorStrtod");
	break;
      case ExceptionResponse::SILENT:
	break;
      }
    }
  }
  return result;
}

} // namespace parse
} // namespace stormm
