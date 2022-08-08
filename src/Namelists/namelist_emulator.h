// -*-c++-*-
#ifndef STORMM_NAMELIST_EMULATOR_H
#define STORMM_NAMELIST_EMULATOR_H

#include <string>
#include <vector>
#include "copyright.h"
#include "Constants/behavior.h"
#include "namelist_element.h"

namespace stormm {
namespace namelist {

using constants::CaseSensitivity;
using constants::ExceptionResponse;

/// \brief Collection of variables to transcribe information contained within a namelist
class NamelistEmulator {
public:

  /// \brief Construct an object to emulate Fortran namelist functionality, with improvements.
  ///
  /// \param title_in   The title of the namelist
  /// \param casing_in  Case sensitivity to abide (default "AUTOMATIC", which in this context means
  ///                   that namelist titles and keywords are case insensitive but values are
  ///                   case sensitive)
  NamelistEmulator(const std::string &title_in,
                   CaseSensitivity casing_in = CaseSensitivity::AUTOMATIC,
		   ExceptionResponse unknown_keyword_policy = ExceptionResponse::WARN,
                   const std::string &help_in = std::string("No description provided"));

  /// \brief Obtain the title of this namelist (i.e. &cntrl or &dock)
  std::string getTitle() const;

  /// \brief Obtain the number of parameters catalogged within this namelist emulator.
  int getKeywordCount() const;

  /// \brief Obtain the case sensitivity setting for this namelist.
  CaseSensitivity getCaseSensitivity() const;

  /// \brief Relay the exception handling policy for this namelist.
  ExceptionResponse getPolicy() const;

  /// \brief Get a keyword from this namelist based on an index.  This is for retrieving the
  ///        keyword itself, not a value associated with a keyword.
  ///
  /// \param index  Index of the keyword in the list held by this namelist.  In fact, this is
  ///               best used to step through the list of keywords in the order they were added,
  ///               not much more.
  std::string getKeyword(size_t index) const;

  /// \brief Get the type of a specific keyword within this namelist.
  ///
  /// \param keyword_query  The keyword of interest
  NamelistType getKeywordKind(const std::string &keyword_query) const;

  /// \brief Get the number of entries associated with a specific keyword.
  ///
  /// \param keyword_query  The keyword of interest
  int getKeywordEntries(const std::string &keyword_query) const;

  /// \brief Get the template size of a keyword, the number of sub-keys it contains.  For
  ///        non-struct keywords the template size is reported as 0.
  int getSubKeyCount(const std::string &keyword_query) const;
  
  /// \brief Test whether a keyword has been set, be that by default or user input.
  ///
  /// Overloaded:
  ///   - Get the status of the first instance of a keyword
  ///   - Get the status of a sub-key within a keyword
  ///   - Get the status of a specific, repeated specification of a keyword
  ///
  /// \param keyword_query  The keyword of interest
  /// \param sub_key        The keyword of interest
  /// \param repeat_no      The number of the repetition to check for its status.  This is needed
  ///                       only in the case of STRUCT-type keywords, as a plain INTEGER, REAL,
  ///                       or STRING keyword will be "established" if there is but one default
  ///                       or user-specified value given--subsequent applications of such a
  ///                       keyword are, by construction, established (status cannot be missing).
  ///                       In contrast, a STRUCT keyword can be specified many times, but the
  ///                       status of its individual members may still be inquestion if not all
  ///                       components of the STRUCT are given in each application.
  /// \{
  InputStatus getKeywordStatus(const std::string &keyword_query) const;
  InputStatus getKeywordStatus(const std::string &keyword_query, const std::string &sub_key,
                               int repeat_no = 0) const;
  /// \}
  
  /// \brief Get a labeled integer value from within the namelist.
  ///
  /// Overloaded:
  ///   - Get an INTEGER value for a non-STRUCT keyword
  ///   - Get an INTEGER value from within a STRUCT keyword
  ///
  /// \param keyword_query  Identifier of the keyword of interest
  /// \param sub_key        Identifier for the member variable within the STRUCT of interest
  /// \param index          For keywords that store multiple values, retrieve this value
  /// \{
  int getIntValue(const std::string &keyword_query, int index = 0) const;
  int getIntValue(const std::string &keyword_query, const std::string &sub_key,
                  int index = 0) const;
  /// \}

  /// \brief Get a labeled real number value from within the namelist.
  ///
  /// Overloaded:
  ///   - Get a REAL value associated with a non-STRUCT keyword
  ///   - Get a REAL value from within a STRUCT keyword
  ///
  /// \param keyword_query  Identifier of the keyword of interest
  /// \param sub_key        Identifier for the member variable within the STRUCT of interest
  /// \param index          For keywords that store multiple values, retrieve this value
  /// \{
  double getRealValue(const std::string &keyword_query, int index = 0) const;
  double getRealValue(const std::string &keyword_query, const std::string &sub_key,
                      int index = 0) const;
  /// \}

  /// \brief Get a labeled string value from within the namelist.
  ///
  /// Overloaded:
  ///   - Get a STRING value associated with a non-STRUCT keyword
  ///   - Get a STRING value from within a STRUCT keyword
  ///
  /// \param keyword_query  Identifier of the keyword of interest
  /// \param sub_key        Identifier for the member variable within the STRUCT of interest
  /// \param index          For keywords that store multiple values, retrieve this value
  /// \{
  std::string getStringValue(const std::string &keyword_query, int index = 0) const;
  std::string getStringValue(const std::string &keyword_query, const std::string &sub_key,
                             int index = 0) const;
  /// \}

  /// \brief Get all integer values assigned to a particular keyword
  ///
  /// \param keyword_query  Identifier of the keyword of interest
  std::vector<int> getAllIntValues(const std::string &keyword_query,
                                   const std::string &sub_key = std::string("")) const;

  /// \brief Get all real values assigned to a particular keyword
  ///
  /// \param keyword_query  Identifier of the keyword of interest
  std::vector<double> getAllRealValues(const std::string &keyword_query,
                                       const std::string &sub_key = std::string("")) const;

  /// \brief Get all string values assigned to a particular keyword
  ///
  /// \param keyword_query  Identifier of the keyword of interest
  std::vector<std::string>
  getAllStringValues(const std::string &keyword_query,
                     const std::string &sub_key = std::string("")) const;

  /// \brief Report the help message associated with a keyword or sub-key.  This can be useful for
  ///        developers who wish to alert users to erroneous input.
  ///
  /// Overloaded:
  ///   - Report user documentation for the namelist as a whole
  ///   - Report user documentation for any keyword (including STRUCTs) within the namelist
  ///   - Report user documentation for a sub-key within a STRUCT keyword in the namelist
  ///
  /// \param keyword_query  Identifier of the keyword of interest
  /// \param sub_key        Identifier for the member variable within the STRUCT of interest
  std::string getHelp() const;
  std::string getHelp(const std::string &keyword_query) const;
  std::string getHelp(const std::string &keyword_query, const std::string &sub_key) const;
  /// \}

  /// \brief Add a keyword to the namelist.
  ///
  /// Overloaded:
  ///   - Add a single keyword (be it a INTEGER, REAL, STRING, or STRUCT namelist element)
  ///   - Add multiple keywords
  ///
  /// \param new_key   The keyword to add
  /// \param new_keys  The keywords to add
  /// \{
  void addKeyword(const std::vector<NamelistElement> &new_keys);
  void addKeywords(const std::vector<NamelistElement> &new_keys);
  void addKeyword(const NamelistElement &new_key);
  /// \}

  /// \brief Assign values to elements of each particular NamelistType.  These overloaded functions
  ///        can be called from anywhere, but constructors making control objects for programs
  ///        using the STORMM libraries are the ideal place to use them.  They in turn call the
  ///        set(...)Value member functions of the target NamelistElement object in the
  ///        NamelistEmulator.  Returns 1 if the given value was successfully assigned to the
  ///        label or 0 if not.
  ///
  /// Overloaded:
  ///   - Assign a single integer, real, or string value to an INTEGER, REAL, or STRING namelist
  ///     element, respectively
  ///   - Assign a single integer, real, or string value to the INTEGER, REAL, or STRING member of
  ///     a STRUCT namelist element
  /// \{
  int assignElement(const std::string &key, const std::string &value);
  int assignElement(const std::string &key, const std::string &sub_key, const std::string &value);
  /// \}

  /// \brief When loading data for STRUCT-type keywords, the decision to increment the number of
  ///        entries on file cannot be made with the first sub-key assignment.  Instead, the entire
  ///        struct must be read from input before the number of entries can be incremented.
  ///        Because a function (readNamelist(), see input.h) that manages a NamelistEmulator
  ///        loops over the input that can be associated with each list of subk-eys for a given
  ///        STRUCT, that function must go through the NamelistEmulator in order to increment
  ///        the keyword's entry count.
  ///
  /// \param key  The STRUCT-type keyword of interest (its membership in the namelist will be
  ///             verified)
  void triggerResizeBuffer(const std::string &key);

  /// \brief Attach a help message to the namelist itself, to a keyword within the namelist, or
  ///        even to a member variable of a STRUCT-associated keyword in the namelist.  This is
  ///        provided so that developers do not have to include help messages at the initialization
  ///        of each namelist keyword.  This will overwrite existing help messages.
  ///
  /// Overloaded:
  ///   - Attach a help message to the namelist as a whole
  ///   - Attach a help message to any keyword (including STRUCTs) within the namelist
  ///   - Attach a help message to a sub-key within a STRUCT keyword in the namelist
  ///
  /// \param blurb      The help message to attach
  /// \param key        Label of the namelist keyword to which the help message gets attached
  /// \param sub_key    Label of the member variable of the STRUCT namelist keyword to which the
  ///                   help message gets attached
  /// \{
  void addHelp(const std::string &blurb);
  void addHelp(const std::string &key, const std::string &blurb);
  void addHelp(const std::string &key, const std::string &sub_key, const std::string &blurb);
  /// \}

  /// \brief Add a category to a namelist to group its keywords for user documentation
  ///
  /// \param new_category  The name of the new keyword category
  void addCategory(const std::string &new_category);

  /// \brief Place a namelist keyword into one of a list of arbitrary categories defined by the
  ///        developer.  This is for organizing the user documentation.
  ///
  /// \param key             The namelist keyword to find and catagorize
  /// \param category_label  The category to put it in
  void categorizeKeyword(const std::string &key, const std::string &category_label);

  /// \brief Print the documentation for a specific keyword.  The format is fixed in the sense
  ///        that it will have a set indentation, a dash for a bullet point, and the keyword
  ///        printed in a space large enough for a series of related keywords in a column.
  void printKeywordDocumentation(int p_idx, int name_width) const;

  /// \brief Print a detailed message concerning the user documentation for keywords in this
  ///        namelist.  This function can be called from the main program, for any namelists that
  ///        accept its inputs, and is itself called by the printProgramDocumentation() function
  ///        in the docs namespace (see Reporting/custom_help.h)
  void printHelp() const;

  /// \brief Print a complete table of the values for all parameters in this namelist, starting
  ///        including their sources (input statuses, i.e. DEFAULT, MISSING, or USER-SPECIFIED).
  void printContents() const;
  
private:
  std::string title;                       ///< Title of this namelist, i.e. &cntrl
  std::vector<NamelistElement> keywords;   ///< List of all keywords stored in this namelist
  CaseSensitivity casing;                  ///< Case sensitivity of the namelist
  ExceptionResponse policy;                ///< Reaction to unusual inputs: more than one instance
                                           ///<   of a non-repeating keyword, an unknown keyword,
                                           ///<   an invalid value
  std::string help_message;                ///< Help message for this namelist, as a whole
  std::vector<std::string> category_names; ///< Names of categories for keywords / labels in this
                                           ///<   namelist (this is for end-user documentation)

  /// Categories containing lists of labels (for sorting the end-user documentation output)
  std::vector<std::vector<std::string>> categories;

  /// \brief Find the index of a keyword based on its label.  This encapsulates throwing an
  ///        exception if the label cannot be found.  While this is similar to findStringInVector()
  ///        in the nearby parse library (see parse.h), the necessity of making some response if
  ///        the label cannot be found and the unusual shape of the vector elements makes it
  ///        reasonable to keep a separate function for this.
  ///
  /// \param query   The label to seek inside the namelist
  int findIndexByKeyword(const std::string &query) const;

  /// \brief Verify the establishment of a keyword, either by user specification or default.
  ///        This is called before returning any information collected from a namelist and will
  ///        trigger an error if the requested information does not exist.
  ///
  /// \param keyword_query  The label that was sought inside the namelist
  /// \param p_index        The keyword (parameter) index within the larger namelist emulator
  /// \param caller         Calling function that is seeking some keyword's information, i.e.
  ///                       getIntValue()
  void verifyEstablishment(const std::string &keyword_query, const size_t p_index,
                           const char* caller) const;
};

} // namespace namelist
} // namespace stormm

#endif
