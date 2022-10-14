// -*-c++-*-
#ifndef STORMM_MOLOBJ_PROPERTY_H
#define STORMM_MOLOBJ_PROPERTY_H

#include <string>
#include <vector>
#include "copyright.h"
#include "Parsing/textfile.h"
#include "molecule_format_enumerators.h"

namespace stormm {
namespace structure {

using parse::TextFile;
  
/// \brief A molecular or atomic property read from an MDL .mol or SDF file
class MolObjProperty {
public:

  /// \brief The constructor takes all member variable inputs, or the original text and a line
  ///        number within it.
  ///
  /// \param tf            Text of the original .sdf or .mol file, read into RAM
  /// \param line_number   Number of the line on which to read the data
  /// \param line_advance  The line advancement position, which will account for data lines in
  ///                      some properties
  /// \param title         The title of the structure, if known, for error tracing purposes
  /// \{
  MolObjProperty(const char4 code_in = { ' ', ' ', ' ', ' ' }, int substrate_in = -1,
                 int entry_count_in = 0, int entry_depth_in = 0,
                 const std::vector<MolObjPropField> &entry_detail_in = {},
                 const std::vector<MolObjIndexKind> &entry_adjustment_in = {},
                 const std::vector<int> &int_data_in = {},
                 const std::vector<double> &real_data_in = {},
                 const std::vector<std::string> &str_data_in = {},
                 const std::vector<std::string> &data_lines_in = {});

  MolObjProperty(const TextFile &tf, int line_number, int *line_advance,
                 const std::string &title = std::string(""));
  /// \}

  /// \brief Get the kind of property according to the internal enumerator.
  MolObjPropertyKind getKind() const;

  /// \brief Get the property code.  This will indicate whether to obliterate certain types of
  ///        information from the atom block of the V2000 MDL MOL format entry.
  char4 getCode() const;

  /// \brief Get the substrate atom or S-group.  The nature of the property will indicate which
  ///        type of substrate this is.
  int getSubstrate() const;

  /// \brief Get the number of entries for this property.
  int getEntryCount() const;

  /// \brief Get a value from the property for a specific entry.  Inputs to these functions
  ///        will be checked by checkAttributeValidity().
  ///
  /// \param entry_index      The index of the entry to access
  /// \param attribute_index  The index of the attribute to retrieve
  /// \{
  int getIntegerValue(int entry_index, int attribute_index) const;
  double getRealValue(int entry_index, int attribute_index) const;
  char4 getChar4Value(int entry_index, int attribute_index) const;
  std::string getStringValue(int entry_index, int attribute_index) const;
  /// \}

  /// \brief Get a const pointer to one of the data lines.
  const std::string& getDataLine(int index) const;
  
  /// \brief Define the property code.
  ///
  /// Overloaded:
  ///   - Provide the first letter and the final three
  ///   - Provide all four letters as a pre-made char tuple
  ///
  /// \param x        The first letter of the three-letter tuple in the property code
  /// \param y        The second letter of the three-letter tuple defining the property code
  /// \param z        The last letter of the three-letter tuple defining the property code
  /// \param major    The first letter that will appear on the line if the property is written to
  ///                 a file
  /// \param code_in  A pre-made code to apply, ordered 3-4-5-0 in terms of the indices that each
  ///                 character would occupy on a property line of an MDL MOL format entry
  /// \{
  void setCode(char x, char y, char z, char major = 'M');

  void setCode(char4 code_in);
  /// \}

  /// \brief Define the substrate to be used, whether an atom or an S-group.
  ///
  /// \param index  Atom or S-group index of the substrate.  This should be provided for the
  ///               C / C++ array element, not the file format (+1 will be added if and when the
  ///               information is committed to a file).
  void setSubstrate(int index);

  /// \brief Define the entry format for the property.  If applied to an existing MolObjProperty
  ///        with nonzero depth, this will create an error.
  ///
  /// \param entry_detail_in  List of data types of the entry elements
  /// \param entry_adjustment_in  List of index adjustments for each element of the entry.  Atom
  ///        and bond elements will have their indices adjusted by -1, but only if the information
  ///        is read from a file--programmatic input is still expected to occur in the C / C++
  ///        array numbering.
  void setEntryFormat(const std::vector<MolObjPropField> &entry_detail_in,
                      const std::vector <MolObjIndexKind> &entry_adjustment_in);

  /// \brief Add an entry to the MDL MOL V2000 format property.  The layout must match that
  ///        established in the detail array.  Adjustments will not be applied to the indexing for
  ///        input originating within the program.  Other methods will adjust the indexing of
  ///        atoms, bonds, and S-groups named in auxiliary user input when making properties out
  ///        of such information.
  ///
  /// \brief int_data_in     All integer data, given in the order that INTEGER MolObjPropField
  ///                        types appear in the detail array.  No char4, real or string data
  ///                        indices are included in this array.
  /// \param char4_data_in   List of char4 data elements, given back-to-back in the order that
  ///                        CHAR4 MolObjPropField types appear in the detail array.
  /// \param real_data_in    List of real data, ordered as the other arrays in this input.
  /// \param string_data_in  List of string data, ordered as the other arrays in this input.
  void addEntry(const std::vector <int> &int_data_in, const std::vector <char4> &char4_data_in,
                const std::vector<double> &real_data_in,
                const std::vector<std::string> &str_data_in);

private:
  char4 code;               ///< A three-letter code indicating what the property is.  The "w"
                            ///<   member stores the first letter on the line, which is usually
                            ///<   but not always 'M'.
  MolObjPropertyKind kind;  ///< The type of property
  int substrate;            ///< One atom or S-group that is central to all entries, relevant to
                            ///<   some properties
  int entry_count;          ///< Number of entries (some properties have maximum numbers of entries
                            ///<   hard-wired into the format, and thus into the constructor)
  int entry_depth;          ///< The number of fields in each entry

  /// Nature of each field in each entry.  The most important classifications are INTEGER and
  ///   CHAR4, although some properties contain longer strings or real numbers.
  std::vector<MolObjPropField> entry_detail;

  /// Indications of adjustments that must be made to convert from the MDL MOL format conventions
  /// into C / C++ array indexing.
  std::vector<MolObjIndexKind> entry_adjustment;

  std::vector<int> int_data;                 ///< Data for all entries, ordered for entries A, B,
                                             ///<   and C with depth 3: { A1, A2, A3, B1, B2, B3,
                                             ///<   C1, C2, C3 }
  std::vector<double> real_data;             ///< Real-valued information for the property.  If a
                                             ///<   component of an entry is REAL, the
                                             ///<   corresponding value in int_data will refer to
                                             ///<   the index of data_str at which to find the
                                             ///<   information.
  std::vector<std::string> str_data;         ///< String data for all entries, covering the
                                             ///<   STRING enumerations of the MolObjPropField
                                             ///<   entry details.  If a component of an entry is a
                                             ///<   STRING, the corresponding value in int_data
                                             ///<   will refer to the index of data_str at which to
                                             ///<   find the information.
  std::vector<std::string> data_lines;       ///< Lines of data, stored one line per string, that
                                             ///<   complement the property

  /// \brief Extract the number of entries for the property.  Returns FALSE if there is no error
  ///        encountered, TRUE if there is a problem.
  ///
  /// \param line_ptr   The line containing the property text
  /// \param start_pos  The starting position at which to read the value (default 6)
  /// \param length     The expected length of the entry count (default 3)
  bool readEntryCount(const char* line_ptr, int start_pos = 6, int length = 3);

  /// \brief Extract the index of the substrate atom or group for the property.  Returns FALSE if
  ///        there is no error encountered, TRUE if there is a problem.
  ///
  /// \param line_ptr   The line containing the property text
  /// \param start_pos  The starting position at which to read the value (default 7)
  /// \param length     The expected length of the entry count (default 3)
  bool readSubstrateIndex(const char* line_ptr, int start_pos = 7, int length = 3);

  /// \brief Parse the list of entries for the property, filling out the arrays of integer, real,
  ///        and string data in the process.
  ///
  /// \param tf           The original text file (passed in for error-tracing purposes)
  /// \param line_number  Number of the line containing the property text
  /// \param start_pos    Starting position for data pertaining to the first entry
  /// \param limits       Bounds array for reading each entry.  Formatting present in entry_detail
  ///                     will impose some limits on what is read (white space can be implied by
  ///                     using a char4 to read between limits 5 characters apart, for example),
  ///                     but the critical feature is that the character at which to begin reading
  ///                     the nth data element is given in the nth position of limits.
  void parseEntries(const TextFile &tf, int line_number, int start_pos,
                    const std::vector<int> &limits);

  /// \brief Check the validity of an attribute request.  Getting information from properties can
  ///        be a tedious and error-prone process.  These guardrails will help developers find
  ///        information safely.
  ///
  /// \param entry_index      The index of the entry to access (receives a bounds check)
  /// \param attribute_index  The index of the attribute to retrieve (receives a bounds check)
  /// \param expectation      The expected nature of the attribute (checked for sanity)
  void checkAttributeValidity(int entry_index, int attribute_index,
                              MolObjPropField expectation) const;
};

} // namespace structure
} // namespace stormm

#endif
