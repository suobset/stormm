#include <string>
#include <vector>
#include "copyright.h"
#include "Parsing/parse.h"
#include "Parsing/tabulation.h"
#include "Reporting/display.h"
#include "Reporting/error_format.h"
#include "namelist_emulator.h"

namespace stormm {
namespace namelist {

using errors::RTMessageKind;
using display::terminalHorizontalRule;
using parse::findStringInVector;
using parse::lowercase;
using parse::NumberFormat;
  
//-------------------------------------------------------------------------------------------------
NamelistEmulator::NamelistEmulator(const std::string &title_in, const CaseSensitivity casing_in,
                                   const ExceptionResponse unknown_keyword_policy,
                                   const std::string &help_in) :
  title{title_in},
  keywords{},
  casing{casing_in},
  policy{unknown_keyword_policy},
  help_message{help_in},
  category_names{},
  categories{}
{}

//-------------------------------------------------------------------------------------------------
const std::string& NamelistEmulator::getTitle() const {
  return title;
}

//-------------------------------------------------------------------------------------------------
int NamelistEmulator::getKeywordCount() const {
  return keywords.size();
}

//-------------------------------------------------------------------------------------------------
CaseSensitivity NamelistEmulator::getCaseSensitivity() const {
  return casing;
}

//-------------------------------------------------------------------------------------------------
ExceptionResponse NamelistEmulator::getPolicy() const {
  return policy;
}

//-------------------------------------------------------------------------------------------------
const std::string& NamelistEmulator::getKeyword(const size_t index) const {
  if (index > keywords.size()) {
    rtErr("Namelist \"" + title + "\" has " + std::to_string(keywords.size()) + " keywords, " +
          "not enough to return index " + std::to_string(index) + ".", "getKeyword");
  }
  return keywords[index].getLabel();
}

//-------------------------------------------------------------------------------------------------
NamelistType NamelistEmulator::getKeywordKind(const std::string &keyword_query) const {
  const size_t p_index = findIndexByKeyword(keyword_query);
  if (p_index >= keywords.size()) {
    rtErr("Namelist \"" + title + "\" has no keyword \"" + keyword_query + "\".",
          "NamelistEmulator", "getKeywordKind");
  }
  return keywords[p_index].getKind();
}

//-------------------------------------------------------------------------------------------------
int NamelistEmulator::getKeywordEntries(const std::string &keyword_query) const {
  const size_t p_index = findIndexByKeyword(keyword_query);
  if (p_index >= keywords.size()) {
    rtErr("Namelist \"" + title + "\" has no keyword \"" + keyword_query + "\".",
          "NamelistEmulator", "getKeywordEntries");
  }
  switch (keywords[p_index].kind) {
  case NamelistType::INTEGER:
  case NamelistType::REAL:
  case NamelistType::STRING:
    return keywords[p_index].getEntryCount();
  case NamelistType::STRUCT:
    return keywords[p_index].getEntryCount() *
           (keywords[p_index].getEstablishment() != InputStatus::MISSING);
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
int NamelistEmulator::getSubKeyCount(const std::string &keyword_query) const {
  const size_t p_index = findIndexByKeyword(keyword_query);
  if (p_index >= keywords.size()) {
    rtErr("Namelist \"" + title + "\" has no keyword \"" + keyword_query + "\".",
          "NamelistEmulator", "getSubKeyCount");
  }
  return keywords[p_index].getTemplateSize();
}

//-------------------------------------------------------------------------------------------------
InputStatus NamelistEmulator::getKeywordStatus(const std::string &keyword_query) const {
  const size_t p_index = findIndexByKeyword(keyword_query);
  if (p_index >= keywords.size()) {
    rtErr("Namelist \"" + title + "\" has no keyword \"" + keyword_query + "\".",
          "NamelistEmulator", "getKeywordEntries");
  }
  return keywords[p_index].getEstablishment();
}

//-------------------------------------------------------------------------------------------------
InputStatus NamelistEmulator::getKeywordStatus(const std::string &keyword_query,
                                               const std::string &sub_key,
                                               const int repeat_no) const {
  const size_t p_index = findIndexByKeyword(keyword_query);
  if (p_index >= keywords.size()) {
    rtErr("Namelist \"" + title + "\" has no keyword \"" + keyword_query + "\".",
          "NamelistEmulator", "getKeywordEntries");
  }
  return keywords[p_index].getEstablishment(sub_key, repeat_no);
}

//-------------------------------------------------------------------------------------------------
int NamelistEmulator::getIntValue(const std::string &keyword_query, const int index) const {
  const size_t p_index = findIndexByKeyword(keyword_query);
  if (p_index >= keywords.size()) {
    rtErr("Namelist \"" + title + "\" has no keyword \"" + keyword_query + "\".",
          "NamelistEmulator", "getIntValue");
  }
  verifyEstablishment(keyword_query, p_index, "getIntValue");
  return keywords[p_index].getIntValue(index);
}

//-------------------------------------------------------------------------------------------------
int NamelistEmulator::getIntValue(const std::string &keyword_query, const std::string &sub_key,
                                  const int index) const {
  const size_t p_index = findIndexByKeyword(keyword_query);
  if (p_index >= keywords.size()) {
    rtErr("Namelist \"" + title + "\" has no keyword \"" + keyword_query + "\".",
          "NamelistEmulator", "getIntValue");
  }
  verifyEstablishment(keyword_query, p_index, "getIntValue");
  return keywords[p_index].getIntValue(sub_key, index);
}

//-------------------------------------------------------------------------------------------------
double NamelistEmulator::getRealValue(const std::string &keyword_query, const int index) const {
  const size_t p_index = findIndexByKeyword(keyword_query);
  if (p_index >= keywords.size()) {
    rtErr("Namelist \"" + title + "\" has no keyword \"" + keyword_query + "\".",
          "NamelistEmulator", "getRealValue");
  }
  verifyEstablishment(keyword_query, p_index, "getRealValue");
  return keywords[p_index].getRealValue(index);
}

//-------------------------------------------------------------------------------------------------
double NamelistEmulator::getRealValue(const std::string &keyword_query, const std::string &sub_key,
                                      const int index) const {
  const size_t p_index = findIndexByKeyword(keyword_query);
  if (p_index >= keywords.size()) {
    rtErr("Namelist \"" + title + "\" has no keyword \"" + keyword_query + "\".",
          "NamelistEmulator", "getRealValue");
  }
  verifyEstablishment(keyword_query, p_index, "getRealValue");
  return keywords[p_index].getRealValue(sub_key, index);
}

//-------------------------------------------------------------------------------------------------
const std::string& NamelistEmulator::getStringValue(const std::string &keyword_query,
                                                    const int index) const {
  const size_t p_index = findIndexByKeyword(keyword_query);
  if (p_index >= keywords.size()) {
    rtErr("Namelist \"" + title + "\" has no keyword \"" + keyword_query + "\".",
          "NamelistEmulator", "getStringValue");
  }
  verifyEstablishment(keyword_query, p_index, "getStringValue");
  return keywords[p_index].getStringValue(index);
}

//-------------------------------------------------------------------------------------------------
const std::string& NamelistEmulator::getStringValue(const std::string &keyword_query,
                                                    const std::string &sub_key,
                                                    const int index) const {
  const size_t p_index = findIndexByKeyword(keyword_query);
  if (p_index >= keywords.size()) {
    rtErr("Namelist \"" + title + "\" has no keyword \"" + keyword_query + "\".",
          "NamelistEmulator", "getStringValue");
  }
  verifyEstablishment(keyword_query, p_index, "getStringValue");
  return keywords[p_index].getStringValue(sub_key, index);
}

//-------------------------------------------------------------------------------------------------
std::vector<int> NamelistEmulator::getAllIntValues(const std::string &keyword_query,
                                                   const std::string &sub_key) const {
  const size_t p_index = findIndexByKeyword(keyword_query);
  if (p_index >= keywords.size()) {
    rtErr("Namelist \"" + title + "\" has no keyword \"" + keyword_query + "\".",
          "NamelistEmulator", "getAllIntValues");
  }
  verifyEstablishment(keyword_query, p_index, "getAllIntValues");
  return keywords[p_index].getIntValue(sub_key);
}

//-------------------------------------------------------------------------------------------------
std::vector<double> NamelistEmulator::getAllRealValues(const std::string &keyword_query,
                                                       const std::string &sub_key) const {
  const size_t p_index = findIndexByKeyword(keyword_query);
  if (p_index >= keywords.size()) {
    rtErr("Namelist \"" + title + "\" has no keyword \"" + keyword_query + "\".",
          "NamelistEmulator", "getAllRealValues");
  }
  verifyEstablishment(keyword_query, p_index, "getAllRealValues");
  return keywords[p_index].getRealValue(sub_key);
}

//-------------------------------------------------------------------------------------------------
std::vector<std::string> NamelistEmulator::getAllStringValues(const std::string &keyword_query,
                                                              const std::string &sub_key) const {
  const size_t p_index = findIndexByKeyword(keyword_query);
  if (p_index >= keywords.size()) {
    rtErr("Namelist \"" + title + "\" has no keyword \"" + keyword_query + "\".",
          "NamelistEmulator", "getAllStringValues");
  }
  verifyEstablishment(keyword_query, p_index, "getAllStringValues");
  return keywords[p_index].getStringValue(sub_key);
}

//-------------------------------------------------------------------------------------------------
const std::string& NamelistEmulator::getHelp() const {
  return help_message;
}

//-------------------------------------------------------------------------------------------------
const std::string& NamelistEmulator::getHelp(const std::string &keyword_query) const {
  const size_t p_index = findIndexByKeyword(keyword_query);
  if (p_index >= keywords.size()) {
    rtErr("Namelist \"" + title + "\" has no keyword \"" + keyword_query + "\".",
          "NamelistEmulator", "getHelp");
  }
  return keywords[p_index].help_message;
}

//-------------------------------------------------------------------------------------------------
const std::string& NamelistEmulator::getHelp(const std::string &keyword_query,
                                      const std::string &sub_key) const {
  const size_t p_index = findIndexByKeyword(keyword_query);
  if (p_index >= keywords.size()) {
    rtErr("Namelist \"" + title + "\" has no keyword \"" + keyword_query + "\".",
          "NamelistEmulator", "getHelp");
  }
  if (keywords[p_index].kind != NamelistType::STRUCT) {
    rtErr("The keyword \"" + keyword_query + "\" is a " +
          getNamelistTypeName(keywords[p_index].kind) + " keyword.  Only STRUCTs can produce user "
          "documentation on sub-keys such as \"" + sub_key + "\".", "NamelistEmulator", "getHelp");
  }
  const size_t member_index = findStringInVector(keywords[p_index].sub_keys, sub_key);
  if (member_index >= keywords[p_index].sub_help_messages.size()) {
    rtErr("No sub-key \"" + sub_key + "\" is present in namelist \"" + title + "\" keyword \"" +
          keyword_query + "\".", "NamelistEmulator", "getHelp");
  }
  return keywords[p_index].sub_help_messages[member_index];
}

//-------------------------------------------------------------------------------------------------
void NamelistEmulator::addKeyword(const std::vector<NamelistElement> &new_keys) {
  const int n_new_key = new_keys.size();
  for (int i = 0; i < n_new_key; i++) {

    // Check that the namelist does not already contain this keyword.  Such would be an outright
    // error, as no developer should be doing this when making an actual program.
    const int n_stored = keywords.size();
    for (int j = 0; j < n_stored; j++) {
      if (new_keys[i].label == keywords[j].label) {
        rtErr("Namelist \"" + title + "\" already has a " +
              getNamelistTypeName(keywords[j].kind) + " keyword " + keywords[j].label +
              ".  A " + getNamelistTypeName(new_keys[i].kind) + " keyword of the same name cannot "
              "be added.", "NamelistEmulator", "addKeyword");
      }
    }
    keywords.push_back(new_keys[i]);
    keywords.back().setPolicy(policy);
  }
}

//-------------------------------------------------------------------------------------------------
void NamelistEmulator::addKeywords(const std::vector<NamelistElement> &new_keys) {
  addKeyword(new_keys);
}

//-------------------------------------------------------------------------------------------------
void NamelistEmulator::addKeyword(const NamelistElement &new_key) {
  std::vector<NamelistElement> new_keys(1, new_key);
  addKeyword(new_keys);
}

//-------------------------------------------------------------------------------------------------
int NamelistEmulator::assignElement(const std::string &key, const std::string &value) {

  // The inputs presumably contain a keyword and value pair.  The question is what to do with the
  // value, and the answer is to find out what sort of input the namelist element expects.
  const size_t param_index = findIndexByKeyword(key);
  if (param_index >= keywords.size()) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("Namelist \"" + title + "\" has no keyword \"" + key + "\".", "NamelistEmulator",
            "assignElement");
    case ExceptionResponse::WARN:
      rtWarn("Namelist \"" + title + "\" has no keyword \"" + key + "\".", "NamelistEmulator",
             "assignElement");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
    return 0;
  }
  const NamelistType param_type = keywords[param_index].kind;
  bool problem = false;
  switch (param_type) {
  case NamelistType::STRUCT:
    rtErr("STRUCT keyword \"" + key + "\" should never be handled in this context.",
          "NamelistEmulator", "assignElement");
  case NamelistType::INTEGER:
    if (verifyNumberFormat(value.c_str(), NumberFormat::INTEGER)) {
      keywords[param_index].setIntValue(stol(value));
      keywords[param_index].establishment = InputStatus::USER_SPECIFIED;
    }
    else {
      problem = true;
    }
    break;
  case NamelistType::REAL:
    if (verifyNumberFormat(value.c_str(), NumberFormat::STANDARD_REAL) ||
        verifyNumberFormat(value.c_str(), NumberFormat::SCIENTIFIC)) {
      keywords[param_index].setRealValue(stod(value));
      keywords[param_index].establishment = InputStatus::USER_SPECIFIED;
    }
    else {
      problem = true;
    }
    break;
  case NamelistType::STRING:
    keywords[param_index].setStringValue(value);
    keywords[param_index].establishment = InputStatus::USER_SPECIFIED;
    break;
  }

  // Respond to input errors
  if (problem) {
    switch(policy) {
    case ExceptionResponse::DIE:
      rtErr("Keyword \"" + key + "\" in namelist \"" + title + "\" accepts " +
            getNamelistTypeName(param_type) + " values.  " + value + " is invalid.",
            "NamelistEmulator", "assignElement");
    case ExceptionResponse::WARN:
      rtWarn("Keyword \"" + key + "\" in namelist \"" + title + "\" accepts " +
             getNamelistTypeName(param_type) + " values.  " + value + " is invalid and no new "
             "value will be assigned.", "NamelistEmulator", "assignElement");
      return 0;
    case ExceptionResponse::SILENT:
      return 0;
    }
  }
  return 1;
}

//-------------------------------------------------------------------------------------------------
int NamelistEmulator::assignElement(const std::string &key, const std::string &sub_key,
                                    const std::string &value) {

  // The inputs contain a keyword, a sub-key, and a value.  This is only valid for STRUCTs.
  const size_t param_index = findIndexByKeyword(key);
  if (param_index >= keywords.size()) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("Namelist \"" + title + "\" has no keyword \"" + key + "\".", "NamelistEmulator",
            "assignElement");
    case ExceptionResponse::WARN:
      rtWarn("Namelist \"" + title + "\" has no keyword \"" + key + "\".", "NamelistEmulator",
             "assignElement");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
    return 0;
  }
  const NamelistType param_type = keywords[param_index].kind;
  bool problem = false;
  switch (param_type) {
  case NamelistType::STRUCT:
    {
      // Search within this keyword.  Use the getter function for encapsulated functionality,
      // even though friendship lets this NamelistEmulator go right in.
      NamelistType subtype = keywords[param_index].getKind(sub_key);
      bool problem = false;
      switch (subtype) {
      case NamelistType::STRUCT:
        break;
      case NamelistType::INTEGER:
        if (verifyNumberFormat(value.c_str(), NumberFormat::INTEGER)) {
          keywords[param_index].setIntValue(sub_key, stol(value));
        }
        else {
          problem = true;
        }
        break;
      case NamelistType::REAL:
        if (verifyNumberFormat(value.c_str(), NumberFormat::STANDARD_REAL) ||
            verifyNumberFormat(value.c_str(), NumberFormat::SCIENTIFIC)) {
          keywords[param_index].setRealValue(sub_key, stod(value));
        }
        else {
          problem = true;
        }
        break;
      case NamelistType::STRING:
        keywords[param_index].setStringValue(sub_key, value);
      }

      // Respond to input errors
      if (problem) {
        switch(policy) {
        case ExceptionResponse::DIE:
          rtErr("In namelist \"" + title + "\", keyword \"" + key + "\", sub-key \"" +
                sub_key + "\" accepts " + getNamelistTypeName(subtype) + " values.  " +
                value + " is invalid.", "NamelistEmulator", "assignElement");
        case ExceptionResponse::WARN:
          rtWarn("In namelist \"" + title + "\", keyword \"" + key + "\", sub-key \"" +
                 sub_key + "\" accepts " + getNamelistTypeName(subtype) + " values.  " +
                 value + " is invalid and no new value will be assigned.", "NamelistEmulator",
                 "assignElement");
          return 0;
        case ExceptionResponse::SILENT:
          return 0;
        }
      }
    }
    break;
  case NamelistType::INTEGER:
  case NamelistType::REAL:
  case NamelistType::STRING:
    rtErr(getNamelistTypeName(param_type) + " keyword \"" + key +
          "\" should never be handled in this context.", "NamelistEmulator", "assignElement");
  }
  return 1;
}

//-------------------------------------------------------------------------------------------------
void NamelistEmulator::triggerResizeBuffer(const std::string &key) {
  const size_t param_index = findIndexByKeyword(key);

  // No exception handling is needed here.  Problems will be caught earlier by assignElement().
  // Check that the keyword is part of the namelist, however.  If the keyword is associated with
  // a STRUCT, which is the purpose of providing this function for NamelistEmulators to control
  // advancement of their unerlying NamelistElement entry counts, ensure that all sub-keys in the
  // STRUCT either have default values or have been specified by the user.
  if (param_index < keywords.size()) {
    if (keywords[param_index].kind == NamelistType::STRUCT) {
      bool advance = false;
      bool problem = false;
      const int est_offset = keywords[param_index].next_entry_index *
                             keywords[param_index].template_size;
      for (int i = 0; i < keywords[param_index].template_size; i++) {
        const size_t eopi = est_offset + i;
        if (keywords[param_index].sub_key_found[i] == false) {

          // The subkey was not specified by the user.  Apply the template default, if available,
          // or issue an alert.  Continue issuing alerts until the entire STRUCT has been examined,
          // then trigger a runtime error.
          if (keywords[param_index].template_establishment[i] == InputStatus::MISSING) {
            keywords[param_index].instance_establishment[eopi] = InputStatus::MISSING;
            if (keywords[param_index].template_requirements[i] == SubkeyRequirement::REQUIRED) {
              rtAlert("Incomplete STRUCT entry for keyword \"" + key + "\" of namelist \"" +
                      title + "\": sub-key \"" + keywords[param_index].sub_keys[i] +
                      "\" does not have a user specification or a default value.",
                      "NamelistEmulator", "triggerResizeBuffer");              
              problem = true;
            }
          }
          else {
            keywords[param_index].instance_establishment[eopi] = InputStatus::DEFAULT;

            // The STRUCT keyword, as a whole, will be described as USER_SPECIFIED if any subkey
            // for at least one instance of the keyword has been specified by the user.
            if (keywords[param_index].establishment == InputStatus::MISSING) {
              keywords[param_index].establishment = InputStatus::DEFAULT;
            }
            
            // Without checking which type of variable this subkey pertains to, apply all of the
            // template values for the ith position to their respective arrays.
            const int tmp_int_default         = keywords[param_index].template_ints[i];
            const double tmp_real_default     = keywords[param_index].template_reals[i];
            const std::string tmp_str_default = keywords[param_index].template_strings[i];
            keywords[param_index].sub_int_values[eopi]    = tmp_int_default;
            keywords[param_index].sub_real_values[eopi]   = tmp_real_default;
            keywords[param_index].sub_string_values[eopi] = tmp_str_default;
            advance = true;
          }
        }
        else {

          // Check that the subkey is valid, then mark the user-specified input
          if (keywords[param_index].template_requirements[i] == SubkeyRequirement::BOGUS) {
            rtAlert("Situational violation in STRUCT keyword \"" + key + "\" of namelist \"" +
                    title + "\": sub-key \"" + keywords[param_index].sub_keys[i] +
                    "\" is bogus in this context.", "NamelistEmulator", "triggerResizeBuffer");
            problem = true;
          }
          else {
            keywords[param_index].establishment = InputStatus::USER_SPECIFIED;
            keywords[param_index].instance_establishment[eopi] = InputStatus::USER_SPECIFIED;
            advance = true;
          }
        }
      }
      if (problem) {
        rtErr("Incomplete specification of data associated with a STRUCT keyword.",
              "NamelistEmulator", "triggerResizeBuffer");
      }
      else if (advance) {
        keywords[param_index].resizeBuffer();
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------
void NamelistEmulator::addHelp(const std::string &blurb) {
  help_message = blurb;
}

//-------------------------------------------------------------------------------------------------
void NamelistEmulator::addHelp(const std::string &key, const std::string &blurb) {
  const size_t p_index = findIndexByKeyword(key);
  if (p_index >= keywords.size()) {
    rtErr("No keyword \"" + key + "\" exists in namelist \"" + title + "\".",
          "NamelistEmulator", "addHelp");
  }
  keywords[p_index].help_message = blurb;
}

//-------------------------------------------------------------------------------------------------
void NamelistEmulator::addHelp(const std::string &key, const std::string &sub_key,
                               const std::string &blurb) {
  const size_t p_index = findIndexByKeyword(key);
  if (p_index >= keywords.size()) {
    rtErr("No keyword \"" + key + "\" is present in namelist \"" + title + "\".",
          "NamelistEmulator", "addHelp");
  }
  if (keywords[p_index].kind != NamelistType::STRUCT) {
    rtErr("The keyword \"" + key + "\" accepts an " +
          getNamelistTypeName(keywords[p_index].kind) + " input.  It must be a STRUCT in "
          "order to incorporate documentation for a member variable's documentation.");
  }
  const size_t member_index = findStringInVector(keywords[p_index].sub_keys, sub_key);
  if (member_index >= keywords[p_index].sub_help_messages.size()) {
    rtErr("No sub-key \"" + sub_key + "\" is present in namelist \"" + title + "\" keyword \"" +
          key + "\".", "NamelistEmulator", "addHelp");
  }
  keywords[p_index].sub_help_messages[member_index] = blurb;
}

//-------------------------------------------------------------------------------------------------
void NamelistEmulator::addCategory(const std::string &new_category) {
  category_names.push_back(new_category);
}
  
//-------------------------------------------------------------------------------------------------
void NamelistEmulator::categorizeKeyword(const std::string &key,
                                         const std::string &category_label) {

  // Ensure that the keyword and category labels are valid
  const size_t category_index = findStringInVector(category_names, category_label);
  if (category_index == category_names.size()) {
    rtErr("Namelist \"" + title + "\" has no category \"" + category_label + "\".",
          "NamelistEmulator", "categorizeKeyword");
  }
  const size_t p_index = findIndexByKeyword(key);
  if (p_index < keywords.size()) {
    categories[category_index].push_back(key);
  }
  else {
    rtErr("Namelist \"" + title + "\" has no keyword \"" + key + "\".", "NamelistEmulator",
          "categorizeKeyword");
  }
}

//-------------------------------------------------------------------------------------------------
void NamelistEmulator::printKeywordDocumentation(const int p_idx, const int name_width) const {
  printf(" - %-*.*s : %s\n", name_width, name_width, keywords[p_idx].label.c_str(),
         terminalFormat("[" + getNamelistTypeName(keywords[p_idx].kind) + "] " +
                        keywords[p_idx].help_message, "", "", 6 + name_width,
                        0, 6 + name_width, 0, RTMessageKind::TABULAR).c_str());
  if (keywords[p_idx].kind == NamelistType::STRUCT) {
    int member_width = 0;
    for (int i = 0; i < keywords[p_idx].template_size; i++) {
      member_width = std::max(member_width,
                              static_cast<int>(keywords[p_idx].sub_keys[i].size()));
    }
    for (int i = 0; i < keywords[p_idx].template_size; i++) {
      printf("   * %-*.*s : %s\n", member_width, member_width,
             keywords[p_idx].sub_keys[i].c_str(),
             terminalFormat("[" + getNamelistTypeName(keywords[p_idx].sub_kinds[i]) + "] " +
                            keywords[p_idx].sub_help_messages[i], "", "", 8 + member_width, 0,
                            8 + member_width, 0, RTMessageKind::TABULAR).c_str());
    }
  }
}

//-------------------------------------------------------------------------------------------------
void NamelistEmulator::printHelp() const {
  terminalHorizontalRule();
  printf("%s\n", terminalFormat("&" + title + ": " + help_message, "", "", 0, 0, title.size() + 3,
                                0, RTMessageKind::TABULAR).c_str());
  terminalHorizontalRule();
  if (categories.size() == 0) {
    printf("\n Keywords:\n");
    const int n_params = keywords.size();

    // Determine the maximum width of any keyword
    int param_width = 0;
    for (int i = 0; i < n_params; i++) {
      param_width = std::max(param_width, static_cast<int>(keywords[i].label.size()));
    }
    for (int i = 0; i < n_params; i++) {
      printKeywordDocumentation(i, param_width);
    }
  }
  else {
    const int n_categories = categories.size();
    for (int i = 0; i < n_categories; i++) {
      printf("\n %s Keywords:\n", category_names[i].c_str());
      const int n_cati = categories[i].size();

      // Determine the maximum width of any keyword
      int param_width = 0;
      for (int j = 0; j < n_cati; j++) {
        param_width = std::max(param_width, static_cast<int>(categories[i][j].size()));
      }
      for (int j = 0; j < n_cati; j++) {
        // Safe to use the output of findIndexByKeyword() as an array index because the contents of
        // every category have already been checked as valid members of the namelist
        printKeywordDocumentation(findIndexByKeyword(categories[i][j]), param_width);
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------
void NamelistEmulator::printContents() const {
  int key_width = 0;
  int kind_width = 0;
  int source_width = 0;
  int value_width = 20;
  for (size_t i = 0; i < keywords.size(); i++) {
    NamelistType i_kind = keywords[i].kind;
    key_width = std::max(key_width, static_cast<int>(keywords[i].label.size()) + 2);
    kind_width = std::max(kind_width, static_cast<int>(getNamelistTypeName(i_kind).size()));
    const int i_entry_count = keywords[i].getEntryCount();
    if (i_kind == NamelistType::STRING) {
      for (int j = 0; j < i_entry_count; j++) {
        value_width = std::max(value_width,
                               static_cast<int>(keywords[i].getStringValue(j).size()) + 2);
      }
    }
    if (i_kind == NamelistType::STRUCT) {
      for (size_t j = 0; j < keywords[i].getTemplateSize(); j++) {
        NamelistType j_kind = keywords[i].sub_kinds[j];
        kind_width = std::max(kind_width, static_cast<int>(getNamelistTypeName(j_kind).size()));
        key_width = std::max(key_width, static_cast<int>(keywords[i].sub_keys[j].size()) + 3);
        if (j_kind == NamelistType::STRING) {
          for (int k = 0; k < i_entry_count; k++) {
            const int slen = keywords[i].getStringValue(keywords[i].sub_keys[j], k).size();
            value_width = std::max(value_width, slen + 2);
          }
        }
      }
    }
    switch (keywords[i].getEstablishment()) {
    case InputStatus::DEFAULT:
    case InputStatus::MISSING:
      source_width = std::max(source_width, 7);
      break;
    case InputStatus::USER_SPECIFIED:
      source_width = std::max(source_width, 14);
      break;
    }
  }

  // Don't let excessively long values make the table ridiculously wide
  value_width = (value_width > 32) ? 32 : value_width;

  // Print the table headings
  const int total_width = key_width + kind_width + value_width + source_width + 10;
  terminalHorizontalRule("-", "-", total_width);
  printf("Contents of &%s:\n", title.c_str());
  terminalHorizontalRule("-", "-", total_width);
  printf(" %-*.*s | %*.*s | %*.*s | %*.*s\n", key_width, key_width, "Keyword", kind_width,
         kind_width, "Type", value_width, value_width, "Value", source_width, source_width,
         "Source");
  terminalHorizontalRule("-", "-", total_width);
  for (size_t i = 0; i < keywords.size(); i++) {

    // Print the first value with a full description
    printf(" %-*.*s | %*.*s | ", key_width, key_width, keywords[i].getLabel().c_str(), kind_width,
           kind_width, getNamelistTypeName(keywords[i].getKind()).c_str());
    switch (keywords[i].getKind()) {
    case NamelistType::INTEGER:
      if (keywords[i].getEstablishment() == InputStatus::MISSING) {
        printf("%*.*s | ", value_width, value_width, "[ Needs user input ]");
      }
      else {
        printf("%*d | ", value_width, keywords[i].getIntValue(0));
      }
      break;
    case NamelistType::REAL:
      if (keywords[i].getEstablishment() == InputStatus::MISSING) {
        printf("%*.*s | ", value_width, value_width, "[ Needs user input ]");
      }
      else {
        printf("%*.4lf | ", value_width, keywords[i].getRealValue(0));
      }
      break;
    case NamelistType::STRING:
      {
        std::string outstring = "\"" + keywords[i].getStringValue(0) + "\"";
        if (outstring.size() > 32) {
          outstring = outstring.substr(0, 28) + "...\"";
        }
        printf("%*.*s | ", value_width, value_width, outstring.c_str());
      }
      break;
    case NamelistType::STRUCT:
      printf("%*.*s | ", value_width, value_width, "");
      break;
    }
    std::string source_tmp;
    switch (keywords[i].getEstablishment()) {
    case InputStatus::DEFAULT:
      source_tmp = "Default";
      break;
    case InputStatus::MISSING:
      source_tmp = "MISSING";
      break;
    case InputStatus::USER_SPECIFIED:
      source_tmp = "User Specified";
      break;
    }
    printf("%*.*s\n", source_width, source_width, source_tmp.c_str());

    // If the keyword is a STRUCT, print its sub-keys
    if (keywords[i].getKind() == NamelistType::STRUCT) {
      for (int j = 0; j < keywords[i].getTemplateSize(); j++) {
        const std::string ij_key = keywords[i].getSubLabel(j);
        const std::string ij_kind = lowercase(getNamelistTypeName(keywords[i].sub_kinds[j]));
        printf("  - %-*.*s | %*.*s | ", key_width - 3, key_width - 3, ij_key.c_str(), kind_width,
               kind_width, ij_kind.c_str());
        switch (keywords[i].sub_kinds[j]) {
        case NamelistType::INTEGER:
          if (keywords[i].getEstablishment() == InputStatus::MISSING) {
            printf("%*.*s | ", value_width, value_width, "[ Needs user input ]");
          }
          else {
            printf("%*d | ", value_width, keywords[i].getIntValue(ij_key, 0));
          }
          break;
        case NamelistType::REAL:
          if (keywords[i].getEstablishment() == InputStatus::MISSING) {
            printf("%*.*s | ", value_width, value_width, "[ Needs user input ]");
          }
          else {
            printf("%*.4lf | ", value_width, keywords[i].getRealValue(ij_key, 0));
          }
          break;
        case NamelistType::STRING:
          {
            std::string outstring = "\"" + keywords[i].getStringValue(ij_key, 0) + "\"";
            if (outstring.size() > 32) {
              outstring = outstring.substr(0, 28) + "...\"";
            }
            printf("%*.*s | ", value_width, value_width, outstring.c_str());
          }
          break;
        case NamelistType::STRUCT:
          break;
        }
        printf("\n");
      }
    }

    // Print information about additional values
    const int nval = keywords[i].getEntryCount();
    if (nval > 1) {
      std::string more_values = "[ +" + std::to_string(nval) + " more values ]";
      printf(" %-*.*s | %*.*s | %*.*s | %*.*s\n", key_width, key_width, "", kind_width, kind_width,
             "", value_width, value_width, more_values.c_str(), source_width, source_width, "");
    }
  }
}

//-------------------------------------------------------------------------------------------------
int NamelistEmulator::findIndexByKeyword(const std::string &query) const {
  const int n_element = keywords.size();
  if (n_element == 0) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("Namelist \"" + title + "\" contains no keywords with which to match keyword " +
            query + ".", "NamelistEmulator", "findIndexByKeyword");
    case ExceptionResponse::WARN:
      rtWarn("Namelist \"" + title + "\" contains no keywords with which to match keyword " +
             query + ".", "NamelistEmulator", "findIndexByKeyword");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
  }
  const int n_letters = query.size();
  if (n_letters == 0) {
    rtErr("A blank namelist (\"" + title + "\") is not searchable.", "NamelistEmulator",
          "findIndexByKeyword");
  }
  for (int i = 0; i < n_element; i++) {
    if (keywords[i].label.size() == n_letters && keywords[i].label[0] == query[0] &&
        keywords[i].label == query) {
      return i;
    }
  }
  
  // If this point is reached, return the end of the list (code that accepts this value will need
  // to check for this possibility and not act on it--there are cases where the code must always
  // throw an exception lest it cause a segmentation fault, or cases where it is OK to warn the
  // user (or even remain silent) and go on.
  return n_element;
}

//-------------------------------------------------------------------------------------------------
void NamelistEmulator::verifyEstablishment(const std::string &keyword_query, const size_t p_index,
                                           const char* caller) const {
  if (keywords[p_index].getEstablishment() == InputStatus::MISSING) {
    rtErr("Namelist \"" + title + "\" keyword \"" + keyword_query + "\" has not been set by "
          "default or by the user.", "NamelistEmulator", caller);
  }
}
  
} // namespace namelist
} // namespace stormm
