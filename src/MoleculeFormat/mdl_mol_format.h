// -*-c++-*-
#ifndef STORMM_MDL_FILE_H
#define STORMM_MDL_FILE_H

#include <string>
#include <vector>
#include "copyright.h"
#include "Chemistry/znumber.h"
#include "Constants/behavior.h"
#include "DataTypes/stormm_vector_types.h"
#include "Parsing/ascii_numbers.h"
#include "Parsing/parse.h"
#include "Trajectory/coordinateframe.h"
#include "Trajectory/phasespace.h"
#include "molecule_file_io.h"
#include "molecule_format_enumerators.h"
#include "molobj_atomlist.h"
#include "molobj_bond.h"
#include "molobj_property.h"

namespace stormm {
namespace structure {

using chemistry::symbolToZNumber;
using constants::CartesianDimension;
using constants::CaseSensitivity;
using constants::ExceptionResponse;
using parse::TextFile;
using trajectory::CoordinateFrame;
using trajectory::PhaseSpace;

/// \brief Default settings for the MDL MOL object atom initializations
/// \{
constexpr char4 default_mdl_atomic_symbol = { ' ', ' ', ' ', ' ' };
constexpr int default_mdl_atomic_number = -1;
constexpr int default_mdl_formal_charge = 0;
constexpr RadicalState default_mdl_radical_state = RadicalState::NONE;
constexpr int default_mdl_isotopic_shift = 0;
constexpr MolObjAtomStereo default_mdl_stereo_parity = MolObjAtomStereo::NOT_STEREO;
constexpr HydrogenAssignment default_hydrogenation = HydrogenAssignment::DO_NOT_HYDROGENATE;
constexpr int default_mdl_implicit_hydrogen = 0;
constexpr int default_mdl_valence_connections = 0;
constexpr int default_mdl_map_count = 0;
constexpr bool default_mdl_stereo_considerations = false;
constexpr bool default_mdl_exact_change = false;
constexpr StereoRetention default_mdl_stereo_retention = StereoRetention::NOT_APPLIED;
/// \}

/// \brief An SText group from an MDL MOL (.mol) or SDF file.  This information is used only by
///        old ISIS / Desktop programs and is otherwise deprecated.
class MolObjSTextGroup {
public:
private:
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

/// \brief A molecular three-dimensional feature.  This special class of MOL object properies has
///        its own data lines.
  
/// \brief A molecule read from an MDL .mol file, or one of many read from a concatenated SDF file.
///        Many of the enumerators above are translated according to member functions of this
///        object based on the documentation in the ctfileformats.pdf file in this library's
///        directory, also available at:
///
///        http://help.accelrysonline.com/ (...)
///          ulm/onelab/1.0/content/ulm_pdfs/direct/reference/ctfileformats2016.pdf
class MdlMolObj {
public:

  /// \brief Constructors for the MDL molecule format object (known as MolObj in RDKit)
  ///
  /// Overloaded:
  ///   - Basic constructor for creating a blank MdlMolObj, referenced by all other constructors'
  ///     initializer lists
  ///   - Constructors based on a file name (string or const char* array)
  ///   - Constructor based on a TextFile object from an SDF container previously committed to RAM
  ///
  /// \param filename        Name of the file to read, containing one MDL file
  /// \param tf              Text file data previously read into memory (this is the way to handle
  ///                        an SDF file)
  /// \param line_start      The irst line at which to begin reading the TextFile object
  /// \param line_end        Last relevant line of the TextFile object at which reading will stop
  /// \param capitalization  Indicate whether atomic symbol capitalization can be ignored when
  ///                        assigning elements to each atom
  /// \param policy          Course of action to take if errors are encountered when inferring
  ///                        atomic elements
  /// \{
  MdlMolObj();
  MdlMolObj(const std::string &filename);
  MdlMolObj(const char* filename);
  MdlMolObj(const TextFile &tf, int line_start = 0, int line_end = -1,
            CaseSensitivity capitalization = CaseSensitivity::YES,
            ExceptionResponse policy = ExceptionResponse::WARN);
  /// \}

  /// \brief Get the system's atom count.
  int getAtomCount() const;
  
  /// \brief Get the number of bonds in the system.
  int getBondCount() const;

  /// \brief Get the { X, Y, Z } coordinate tuple for a particular atom, or for all atoms.
  ///
  /// Overloaded:
  ///   - Get a unique, modifiable tuple for a single atom.
  ///   - Get a const reference to the array of tuples for coordinates of all atoms.
  ///   - Get a modifiable vector of all Cartesian X, Y, or Z coordinates.
  ///
  /// \param index  Index of the atom of interest
  /// \param dim    The Cartesian dimension of interest
  /// \{
  double3 getCoordinates(int index) const;
  const std::vector<double3>& getCoordinates() const;
  std::vector<double> getCoordinates(CartesianDimension dim) const;
  /// \}

  /// \brief Export the coordinates as a PhaseSpace object suitable for basic molecular mechanics
  ///        force computations.
  PhaseSpace exportPhaseSpace() const ;

  /// \brief Export the coordinates as a stripped-down CoordinateFrame object suitable for
  ///        molecular mechanics energy computations.
  CoordinateFrame exportCoordinateFrame() const;

  /// \brief Get the atomic symbol for a particular atom.
  ///
  /// \param index  Index of the atom of interest
  char4 getAtomSymbol(int index) const;

  /// \brief Get a const reference to the atomic symbols for all atoms.
  const std::vector<char4>& getAtomSymbols() const;
  
  /// \brief Get the atomic number of a particular atom.
  ///
  /// \param index  Index number of the atom of interest
  int getAtomicNumber(int index) const;

  /// \brief Get a const reference to the vector of all atomic numbers in the system.
  const std::vector<int>& getAtomicNumbers() const;

  /// \brief Get the formal charge on a particular atom.
  ///
  /// \param index  Index number of the atom of interest
  int getFormalCharge(int index) const;

  /// \brief Get a const reference to the vector of all formal charges.
  const std::vector<int>& getFormalCharges() const;
  
  /// \brief Get the number of properties found in the MDL MOL entry
  int getPropertiesCount() const;
  
private:

  // Items describing quantities of information (most of them from the counts line)
  MdlMolVersion version_no;   ///< The format in which this entry was read (does not necessarily
                              ///<   dictate the version in which it will be written)
  int atom_count;             ///< The number of atoms in the molecule
  int bond_count;             ///< Number of bonds of all types between atoms in the system
  int list_count;             ///< The number of atom (element) lists
  int stext_entry_count;      ///< The number of S-text entries
  int properties_count;       ///< The number of additional properties (the default, and expected,
                              ///<   value is 999)
  int sgroup_count;           ///< The number of S-groups
  int constraint_count;       ///< The number of three-dimensional constraints
  MolObjChirality chirality;  ///< The molecule's chirality (assumes only one significant center)
  int registry_number;        ///< The molecule registry number
 
  // Atomic properties
  std::vector<double3> coordinates;         ///< Cartesian coordinates of all atoms
  std::vector<char4> atomic_symbols;        ///< Symbols for all atoms
  std::vector<int> atomic_numbers;          ///< Atomic numbers for all atoms
  std::vector<int> formal_charges;          ///< Formal charges for all atoms
  std::vector<RadicalState> radicals;       ///< Indications that any atom contains a radical in
                                            ///<   singlet, doublet, or triplet excitation state
  std::vector<int> isotopic_shifts;         ///< Isotope numbers shifting each atoms nuclear mass
                                            ///<   (~1 Dalton from the most common isotopic mass)
  std::vector<MolObjAtomStereo> parities;   ///< Stereochemical parity of each atom
  std::vector<int> implicit_hydrogens;      ///< The number of implicit hydrogens that may be
                                            ///<   considered to reside bonded to a particular atom
                                            ///<   (in addition to those that are explicitly drawn)
  std::vector<bool> stereo_considerations;  ///< Indications of whether each atom has steroisomeric
                                            ///<   considerations.  Having stereochemistry at both
                                            ///<   ends of a double bond indicates that the
                                            ///<   double-bond affects the structural properties
                                            ///<   of the molecule.
  std::vector<int> valence_connections;     ///< The number of bonds made by each atom to others
                                            ///<   in the same molecule
  std::vector<int> atom_atom_mapping_count; ///< The number of atoms mapped to this one in a
                                            ///<   chemical reaction
  std::vector<bool> exact_change_enforced;  ///< Flags to indicate that the changes on an atom
                                            ///<   must be exactly as described in the reaction

  /// The assigned method of adding implicit hydrogens: add at least the number defined in the
  /// field (up to the amount needed to satisfy the anticipated valence shell electron content),
  /// add the number needed to satisfy the valence shell, or do not add hydrogens. 
  std::vector<HydrogenAssignment> hydrogenation_protocol;

  /// Indication of whether stereochemistry is inverted or retained in a chemical reaction
  std::vector<StereoRetention> orientation_stability;

  /// Bonds between atoms
  std::vector<MolObjBond> bonds;

  /// Lists of atomic elements to be used in arbitrary operations
  std::vector<MolObjAtomList> element_lists;

  /// Stext entries (these are deprecated and used by ISIS / Desktop applications only)
  std::vector<MolObjSTextGroup> stext_entries;
  
  /// Properties
  std::vector<MolObjProperty> properties;

  /// Data items: this is what makes the SD file (.sdf) format so extensible.  While technically
  /// not part of the MDL MOL format, they are stored in the MdlMolObj for association.  Data items
  /// will be written back to new .sdf files but not .mol files based on the data in this object.
  std::vector<MolObjDataItem> data_items;
  
  /// Title (first line from the file header)
  std::string title;

  /// Details of the software that generated this MDL .mol file (second line of the header)
  std::string software_details;

  /// General comments for the file (third line of the file header)
  std::string general_comment;

  /// \brief Allocate space for information in amounts described on the MOL entry's counts line.
  void allocate();

  /// \brief Interpret the formal charge of an atom based on an integer code.  A doublet radical
  ///        can also emerge from this analysis.
  ///
  /// \param charge_in   The charge code to interpret
  /// \param atom_index  Index of the atom to which the charge applies (the radical state may also
  ///                    be set by this function)
  void interpretFormalCharge(int charge_in, int atom_index);

  /// \brief Interpret the stereochemical parity of an atom based on an integral numeric code.
  ///
  /// \param setting_in  The code to parse
  MolObjAtomStereo interpretStereoParity(int setting_in);

  /// \brief Interpret the implicit hydrogen count for an atom.  At least this many additional
  ///        hydrogens are implied to be present around an atom center, in addition to any hydrogen
  ///        atoms explicitly placed in the structure.
  ///
  /// \param nh_in       The number of hydrogens that are to be inferred around the atom's
  ///                    structure.  This number is decremented by one, meaning that 1 implies no
  ///                    implicit hydrogen content.
  /// \param atom_index  Atom to which the hydrogen content applies.  This is provided in order to
  ///                    set both the actual number and the action fields for the atom.
  void interpretImplicitHydrogenContent(int nh_in, int atom_index);

  /// \brief Interpret a boolean value from an integer.
  ///
  /// \param code_in  The code to translate.  This one is simple: 1 = TRUE, 0 = FALSE
  /// \param desc     Description of the activity, for error tracing purposes
  bool interpretBooleanValue(int value_in, const std::string &desc);

  /// \brief Translate the count of valence interactions for a particular atom.
  ///
  /// \param count_in  The count to translate.  The input is typically unchanged, but for the fact
  ///                  that 0 and 15 both translate to no valence interactions.
  int interpretValenceNumber(int count_in);

  /// \brief Interpret the stability of stereochemical arrangements listed for each atom.
  ///
  /// \param code_in  A numeric code to be translated into inversion or retention options
  StereoRetention interpretStereoStability(int code_in);

  /// \brief Scan properties (of a V2000-format molecule) and update information that may have been
  ///        read from the atoms block.
  void updateV2kAtomAttributes();
  
  /// \brief Add hydrogens to fill out valence shells based on the stated bond orders, formal
  ///        charges, and elements of the molecule.
  void hydrogenate();
};

/// \brief Read a structure data file (.sdf extension) containing one or more MDL MOL entries.
///        Return the results as a Standard Template Library vector of MDL MOL objects.
///
/// Overloaded:
///   - Take a file name
///   - Take a pre-converted TextFile object
///
/// \param tf       
/// \param file_name
/// \param capitalization  Indicate whether atomic symbol capitalization can be ignored when
///                        assigning elements to each atom
/// \param policy          Course of action to take if errors are encountered when inferring
///                        atomic elements
/// \{
std::vector<MdlMolObj> readStructureDataFile(const TextFile &tf, int low_frame_limit,
                                             int high_frame_limit,
                                             CaseSensitivity capitalization = CaseSensitivity::YES,
                                             ExceptionResponse policy = ExceptionResponse::WARN);
std::vector<MdlMolObj> readStructureDataFile(const std::string &file_name,
                                             CaseSensitivity capitalization = CaseSensitivity::YES,
                                             ExceptionResponse policy = ExceptionResponse::WARN);

std::vector<MdlMolObj> readStructureDataFile(const TextFile &tf,
                                             CaseSensitivity capitalization = CaseSensitivity::YES,
                                             ExceptionResponse policy = ExceptionResponse::WARN);
/// \}
  
} // namespace structure
} // namespace stormm

#endif
