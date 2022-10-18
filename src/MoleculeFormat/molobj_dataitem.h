// -*-c++-*-
#ifndef STORMM_MOLOBJ_DATAITEM_H
#define STORMM_MOLOBJ_DATAITEM_H

#include <string>
#include <vector>
#include "copyright.h"
#include "DataTypes/stormm_vector_types.h"
#include "Parsing/textfile.h"
#include "Potential/energy_enumerators.h"
#include "molecule_format_enumerators.h"

namespace stormm {
namespace structure {

using energy::StateVariable;
using parse::TextFile;

/// \brief A request for one of a sanctioned list of information types to be included in a data
///        item of an SD file.  When printing an SD file, MolObjDataItem objects will be created
///        based on these requests.
class MolObjDataRequest {
public:

  /// \brief The constructor takes arguments corresponding to each member variable.  Various
  ///        overloads allow the object ot be constructed for a specific type of request, ignoring
  ///        member variables that are not relevant.
  /// \{
  MolObjDataRequest(const std::string &title_in = std::string(""),
                    const std::string &label_in = std::string(""));

  MolObjDataRequest(const std::string &title_in, StateVariable energy_component_in,
                    const std::string &label_in);

  MolObjDataRequest(DataRequestKind kind_in, const std::string &title_in,
                    const std::string &message_in, const std::string &label_in);

  MolObjDataRequest(const std::string &title_in, StateVariable valence_kind_in,
                    const std::vector<char4> &atom_types_in, const std::string &label_in);
  /// \}
  
  /// \brief The default copy and move constructors, as well as copy and move assignment operators,
  ///        are applicable to this object which has no const members or pointers to repair.
  /// \{
  MolObjDataRequest(const MolObjDataRequest &original) = default;
  MolObjDataRequest(MolObjDataRequest &&original) = default;
  MolObjDataRequest& operator=(const MolObjDataRequest &other) = default;
  MolObjDataRequest& operator=(MolObjDataRequest &&other) = default;
  /// \}

  /// \brief Get the kind of request.
  DataRequestKind getKind() const;

  /// \brief Get the title for the printed SD file's data item.
  const std::string& getTitle() const;
  
  /// \brief Get the type of energy requested, or produce a runtime error if this request is not
  ///        of the appropriate kind.
  StateVariable getEnergyComponent() const;

  /// \brief Get the atom mask string.  Raise a runtime error if this request is not of the
  ///        appropriate kind.
  const std::string& getAtomMask() const;

  /// \brief Get the valence parameter type.  Raise a runtime error if this request is not of the
  ///        appropriate kind.
  StateVariable getValenceParameter() const;

  /// \brief Get the custom string, again raising a runtime error if the request is not of the
  ///        appropriate kind.
  const std::string& getMessage() const;

  /// \brief Get the vector of atom types defining a valence parameter of interest, again raising
  ///        a runtime error if the request is not of the appropriate kind.
  const std::vector<char4>& getAtomTypes() const;

  /// \brief Get the system label to which the data item should be applied.
  const std::string& getSystemLabel() const;

private:
  DataRequestKind kind;           ///< Define the type of data request
  std::string title;              ///< Title for the data item, to be displayed as "> <title>" in
                                  ///<   the resulting SD file.
  StateVariable energy_component; ///< The type of energy quantity requested, i.e. PROPER_DIHEDRAL
  std::string atom_mask;          ///< Mask defining atoms for which valence interactions are to be
                                  ///<   listed in the SD file's data item, if the request is for
                                  ///<   atom influences
  StateVariable valence_kind;     ///< A type of valence parameter (including restraint terms)
                                  ///<   which, in conjunction with the appropriate number of named
                                  ///<   atom types, will cause the program to output all instances
                                  ///<   of the parameter in the model, listing the parameter's
                                  ///<   settings, atom indices to which it applies and the
                                  ///<   energies of each interaction.
  std::string message;            ///< A custom string provided by the user
  std::vector<char4> atom_types;  ///< Vector of atom types defining a valence
  std::string system_label;       ///< When applying the request to create data items in SD file
                                  ///<   output, this label will be checked against the molecular
                                  ///<   system's label to ensure that the information goes only
                                  ///<   into the intended places.

  /// \brief Check that the data request is of the right kind to return a particular member
  ///        variable, or possibly to perform other functions.
  ///
  /// \param accepted_kind  The kind of request that this object must be making in order for it to
  ///                       function as the developer wants
  void checkKind(DataRequestKind accepted_kind) const;
};
  
/// \brief Store a data item from within an SD file.  Data items begin with a line of the form
///        "> <ITEM_NAME>" (additional specifications are possible), end with a single blank line,
///        and can contain customized format information.  The information will be read as a series
///        of raw strings, one per line and each representing a whole line.  Once read, member
///        functions of the class can extract specific integer, real, char4, or string coded data
///        from within the item.
class MolObjDataItem {
public:

  /// \brief The constructor takes the original text file and the number of the first line of the
  ///        data item, or the information pertinent to the header line and a series of strings
  ///        that will form the body of data in the item.  Any information submitted as part of a
  ///        string such as an atom or bond index, which could then go to an output file, must be
  ///        in the Fortran array indexing (starting at 1).  Integer information read from the
  ///        data body can be adjusted to meet internal array indexing.
  ///
  /// \param item_name_in        Name of the data item (an item_name of "ABCD" will appear as
  ///                            "<ABCD>" on the data item's first line)
  /// \param external_regno_in   An identification number for the compound referencing an external
  ///                            database.  An external identification number of "PD-3647" will be
  ///                            represented as "(PD-3647)" on the data item's header line.
  /// \param maccs_ii_number_in  An integer representing the corresponding field in a MACCS-II
  ///                            database.  A field number such as "6519" will be represented as
  ///                            "DT6519" in the data item header line.
  /// \param header_info         A bit-packed unsigned integer containing, in its low to high bits,
  ///                            whether to display the internal identification number (based on
  ///                            the order in the SD file), the external identification number,
  ///                            the item name, the field number, and finally a "FROM ARCHIVES"
  ///                            declaration.
  /// \param tf                  The original text of the SD file, committed to RAM
  /// \param line_number         Line of the file at which to begin reading the data item
  /// \param line_advance        Line of the file that the construction of the data item leads to
  ///                            (after taking in lines of the item)
  /// \param compound_line_end   Index of the last line of the compound within the SD file (the
  ///                            line contains $$$$)
  /// \param title               The title of the structure, if known, for error tracing purposes
  /// \{
  MolObjDataItem(const TextFile &tf, int line_number, int *line_advance,
                 int compound_line_end = -1,
                 const std::string &title = std::string(""));

  MolObjDataItem(const std::string &item_name_in = std::string(""),
                 const std::string &external_regno_in = std::string(""),
                 int internal_regno_in = -1, int maccs_ii_number_in = -1, uint header_info = 0U,
                 const std::vector<std::string> &body_in = {});
  /// \}

  /// \brief The default copy and move constructors, as well as copy and move assignment operators,
  ///        are applicable to this object which has no const members or pointers to repair.
  /// \{
  MolObjDataItem(const MolObjDataItem &original) = default;
  MolObjDataItem(MolObjDataItem &&original) = default;
  MolObjDataItem& operator=(const MolObjDataItem &other) = default;
  MolObjDataItem& operator=(MolObjDataItem &&other) = default;
  /// \}
  
  /// \brief Get a const reference to the item name, if it exists.  If there is no item name, a
  ///        const reference to a blank string will be returned.
  const std::string& getItemName() const;

  /// \brief Get the external registry number.
  const std::string& getExternalRegistryNumber() const;

  /// \brief Get the internal registry number.  This will return -1 if no such registry number
  ///        exists.
  int getInternalRegistryNumber() const;

  /// \brief Get the MACCS-II database field number.
  int getMaccsFieldNumber() const;

  /// \brief Retrieve a string from the data lines of a data item, assuming that individual words
  ///        on each line are separated by one or more white space characters and that the
  ///        quotation marks "" and '' collect everything between them into a single word.
  ///
  /// Overloaded:
  ///   - Extract a free-format string based on words separated by white space, assuming that
  ///     individual words on each line are separated by one or more white space characters and
  ///     that the quotation marks "" and '' collect everything between them into a single word.
  ///   - Extract a column-formatted string based on strict character indices
  ///
  /// \param element_number  Number of the word on the line (words are separated by white space)
  /// \param start_pos       The starting position within the line at which to begin reading
  /// \param start_pos       Length of the fixed-column reading to perform
  /// \param line_number     Number of the data line on which to find the string
  /// \{
  std::string parseString(int element_number, int line_number) const;
  std::string parseString(int start_pos, int length, int line_number) const;
  /// \}

  /// \brief Retrieve a signed integer from the data lines of a data item.  The assumptions of
  ///        parseString() above apply, as do the overloads and descriptions of formal arguments.
  /// \{
  llint parseInteger(int element_number, int line_number) const;
  llint parseInteger(int start_pos, int length, int line_number) const;
  /// \}

  /// \brief Retrieve an unsigned integer from the data lines of a data item.  The assumptions of
  ///        parseString() above apply, as do the overloads and descriptions of formal arguments.
  /// \{
  ullint parseUnsigned(int element_number, int line_number) const;
  ullint parseUnsigned(int start_pos, int length, int line_number) const;
  /// \}

  /// \brief Retrieve a real number from the data lines of a data item.  The assumptions of
  ///        parseString() above apply, as do the overloads and descriptions of formal arguments.
  /// \{
  double parseReal(int element_number, int line_number) const;
  double parseReal(int start_pos, int length, int line_number) const;
  /// \}

  /// \brief Retrieve a tuple of four characters from the data lines of a data item.  Assumptions
  ///        from parseString() above apply, as do overloads and descriptions of formal arguments.
  /// \{
  char4 parseChar4(int element_number, int line_number) const;
  char4 parseChar4(int start_pos, int length, int line_number) const;
  /// \}

  /// \brief Match this data item with a series of identification tags.
  ///
  /// \param item_name_comp  Item name for comparison
  /// \param ext_regno_comp  External registry number for comparison
  /// \param maccs_ii_no_comp  MACCS-II database field number for comparison (omit leading "DT")
  /// \{
  bool matchItemName(const std::string &item_name_comp, const std::string &ext_regno_comp,
                     int maccs_ii_no_comp = -1) const;

  bool matchItemName(const std::string &item_name_comp, int maccs_ii_no_comp = -1) const;

  bool matchRegistryNumber(const std::string &ext_regno_comp, int maccs_ii_no_comp = -1) const;

  bool matchMaccsField(int maccs_ii_no_comp) const;
  /// \}

  /// \brief Set the item name and apply checks to the result.
  ///
  /// \param item_name_in  The item name to assign
  void setItemName(const std::string &item_name_in);

private:
  std::string item_name;       ///< Name of the data item (optional, but an effective means of
                               ///<   distiction)
  std::string external_regno;  ///< External registry number
  int internal_regno;          ///< Internal registry number, based on the structure index within
                               ///<   the SD file and a pure integer
  int maccs_ii_number;         ///< MACCS-II database field number
  bool use_internal_regno;     ///< Flag to have the header line use the SD file's internal
                               ///<   structure numbering in the data item header line
  bool use_external_regno;     ///< Flag to have the header line use the external registry number
  bool use_item_name;          ///< Flag to use the item name in the header between < > marks
  bool use_maccs_ii_number;    ///< Flag to use the MACCS-II field number in the header after "DT"
  bool note_archives;          ///< Flag to have the header line note "FROM ARCHIVES"

  /// Data lines from the item, one per array element
  std::vector<std::string> body;

  /// \brief Validate the choice of item name.
  void validateItemName() const;
};

} // namespace structure
} // namespace stormm

#endif
