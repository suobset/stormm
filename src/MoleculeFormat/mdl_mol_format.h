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
constexpr double default_mdl_formal_charge = 0.0;
constexpr bool default_mdl_doublet_radical_state = false;
constexpr int default_mdl_isotopic_shift = 0;
constexpr MolObjAtomStereo default_mdl_stereo_parity = MolObjAtomStereo::NOT_STEREO;
constexpr int default_mdl_implicit_hydrogen = 0;
constexpr int default_mdl_valence_connections = 0;
constexpr int default_mdl_map_count = 0;
constexpr bool default_mdl_stereo_considerations = false;
constexpr bool default_mdl_exact_change = false;
constexpr StereoRetention default_mdl_stereo_retention = StereoRetention::NOT_APPLIED;
constexpr MolObjBondOrder default_mdl_bond_order = MolObjBondOrder::SINGLE;
constexpr MolObjBondStereo default_mdl_bond_stereochemistry = MolObjBondStereo::NOT_STEREO;
constexpr MolObjRingState default_mdl_ring_status = MolObjRingState::EITHER;
constexpr MolObjReactionCenter default_mdl_bond_reactivity = MolObjReactionCenter::NON_CENTER;
/// \}

/// \brief A bond, as presented in the MDL molecule file format.  This unguarded struct will be
///        returned to the developer from a private array inside of the MdlMolObj object, so
///        further protection would be a hindrance.
class MolObjBond {
public:

  /// \brief The constructor can take all member variables, or just the atoms so that more
  ///        information can be filled in later.
  ///
  /// Overloaded:
  ///   - Construct a blank object with -1 atom indices to indicate its invalid nature
  ///   - Construct an object with just the atom indices
  ///   - Construct a complete object with all details pre-loaded
  ///   - Construct the object based on text input of trusted bounds
  ///
  /// \param tf           Text of the original .sdf or .mol file, read into RAM
  /// \param line_number  Number of the line on which to read the data
  /// \param title        The title of the structure, if known, for error tracing purposes
  /// \{
  MolObjBond();
  MolObjBond(int i_atom_in, int j_atom_in);
  MolObjBond(int i_atom_in, int j_atom_in, MolObjBondOrder order_in, MolObjBondStereo stereo_in,
             MolObjRingState ring_state_in, MolObjReactionCenter reactivity_in);
  MolObjBond(const TextFile &tf, int line_number, const std::string &title = std::string(""));
  /// \}

  /// \brief The default copy and move constructors as well as assignment operators are adequate.
  /// \{
  MolObjBond(const MolObjBond &original) = default;
  MolObjBond(MolObjBond &&original) = default;
  MolObjBond& operator=(const MolObjBond &original) = default;
  MolObjBond& operator=(MolObjBond &&original) = default;
  /// \}

  /// \brief Get the first atom in the bond.  Having separate functions for each atom is a more
  ///        intuitive way to offer the getter functions, whereas in a BoundedRestraint object the
  ///        atom index getter takes an argument for the first, second, third, or fourth atom.
  int getFirstAtom() const;

  /// \brief Get the second atom in the bond.
  int getSecondAtom() const;

  /// \brief Get the order of the bond.
  MolObjBondOrder getOrder() const;

  /// \brief Get the order of the bond.
  MolObjBondStereo getStereochemistry() const;

  /// \brief Get the ring status--is the bond known to be part of a ring?
  MolObjRingState getRingStatus() const;

  /// \brief Get the reactive potential of the bond.
  MolObjReactionCenter getReactivity() const;

  /// \brief Set the index of the first atom in the bond.
  void setFirstAtom(int index_in);

  /// \brief Set the index of the second atom in the bond.
  void setSecondAtom(int index_in);

  /// \brief Set the order of the bond, perhaps after computations with an associated
  ///        ChemicalFeatures object.
  void setOrder(MolObjBondOrder order_in);

  /// \brief Set the stereochemical details of the bond.
  void setStereochemistry(MolObjBondStereo stereo_in);

  /// \brief Mark the status of the bond with respect to any ring features.
  void setRingStatus(MolObjRingState status_in);

  /// \brief Mark the reactive potential of the bond.
  void setReactivity(MolObjReactionCenter potential_in);
  
private:
  int i_atom;                       ///< The first atom in the bond
  int j_atom;                       ///< The second atom in the bond
  MolObjBondOrder order;            ///< The bond order (single, double, aromatic, etc.)
  MolObjBondStereo stereo;          ///< Indicator of the bond stereochemistry
  MolObjRingState ring_state;       ///< Indicator of whether the atom is part of a ring
  MolObjReactionCenter reactivity;  ///< Indicator of a bond as a center of reactivity

  /// \brief Interpret a code for the order of a bond in the structure.
  ///
  /// \param code_in  A numeric code to be translated into the bond order
  /// \param title    Title of the parent structure, for error tracing purposes
  MolObjBondOrder interpretBondOrder(int code_in, const std::string &title);

  /// \brief Interpret a code for the stereochemistry of a bond in the structure.
  ///
  /// \param code_in  A numeric code to be translated into the bond stereochemistry
  /// \param title    Title of the parent structure, for error tracing purposes
  MolObjBondStereo interpretBondStereochemistry(int code_in, const std::string &title);
  
  /// \brief Interpret a code for the ring status of a bond in the structure.
  ///
  /// \param code_in  A numeric code to be translated into the ring status
  /// \param title    Title of the parent structure, for error tracing purposes
  MolObjRingState interpretRingState(int code_in, const std::string &title);

  /// \brief Interpret a code for the reactive potential of a bond in the structure.
  ///
  /// \param code_in  A numeric code to be translated into the reactive potential
  /// \param title    Title of the parent structure, for error tracing purposes
  MolObjReactionCenter interpretBondReactivePotential(int code_in, const std::string &title);
};

/// \brief An atom list entry (this object can be assembled either from one of the deprecated
///        V2000 format lines after the bonds block, or from one of the "M  ALS" properties)
class MolObjAtomList {
public:

  /// \brief The constructor can take all member variables (and all come with default values to
  ///        let this form of the constructor serve as the blank object constructor), or a pointer
  ///        to the line of a text file from which the information shall come.
  ///
  /// \param tf           Text of the original .sdf or .mol file, read into RAM
  /// \param line_number  Number of the line on which to read the data
  /// \param title        The title of the structure, if known, for error tracing purposes
  /// \{
  MolObjAtomList(const std::vector<int> &atomic_numbers_in = {}, bool exclusions_in = false,
                 int atom_attachment_in = 0);

  MolObjAtomList(const TextFile &tf, int line_number, const std::string &title = std::string(""));
  /// \}

private:
  int entry_count;                  ///< The number of atomic elements (identified by Z-numbers)
                                    ///<   in this list
  std::vector<int> atomic_numbers;  ///< Atomic (Z-) numbers of atoms that are to be excluded or
                                    ///<   included by processing this list
  bool exclusions;                  ///< Indicate whether this list covers atomic numbers which
                                    ///<   are to be excluded (TRUE) or included (FALSE)
  int atom_attachment;              ///< Attachment point of the list, an index of an atom in the
                                    ///<   molecule itself
};

/// \brief An SText group from an MDL MOL (.mol) or SDF file.  This information is used only by
///        old ISIS / Desktop programs and is otherwise deprecated.
class MolObjSTextGroup {
public:
private:
};
  
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
                 const std::vector<int> &int_data_in = {},
                 const std::vector<double> &real_data_in = {},
                 const std::vector<std::string> &str_data_in = {},
                 const std::vector<std::string> &data_lines_in = {});

  MolObjProperty(const TextFile &tf, int line_number, int *line_advance,
                 const std::string &title = std::string(""));
  /// \}
  
private:
  char4 code;                                ///< A three-letter code indicating what the property
                                             ///<   is.  The "w" member stores the first letter on
                                             ///<   the line, which is usually but not always 'M'.
  MolObjPropertyKind kind;                   ///< The type of property
  int substrate;                             ///< One atom or S-group that is central to all
                                             ///<   entries, relevant to some properties 
  int entry_count;                           ///< Number of entries (some properties have maximum
                                             ///<   numbers of entries hard-wired into the format,
                                             ///<   and thus into the constructor)
  int entry_depth;                           ///< The number of fields in each entry
  std::vector<MolObjPropField> entry_detail; ///< Nature of each field in each entry.  The most
                                             ///<   important classifications are INTEGER and
                                             ///<   CHAR4, although some properties contain longer
                                             ///<   strings.
                                             ///<   real numbers in the 3-character column format.
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
  std::vector<double> formal_charges;       ///< Formal charges for all atoms
  std::vector<bool> doublet_radicals;       ///< Indications that any atom contains a radical in
                                            ///<   doublet excitation state
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
  /// \param charge_in           The charge code to interpret
  /// \param is_doublet_radical  Flag to indicate code 4, no net charge but a double radical
  double interpretFormalCharge(int charge_in, bool *is_doublet_radical);

  /// \brief Interpret the stereochemical parity of an atom based on an integral numeric code.
  ///
  /// \param setting_in  The code to parse
  MolObjAtomStereo interpretStereoParity(int setting_in);

  /// \brief Interpret the implicit hydrogen count for an atom.  At least this many additional
  ///        hydrogens are implied to be present around an atom center, in addition to any hydrogen
  ///        atoms explicitly placed in the structure.
  ///
  /// \param nh_in  The number of hydrogens that are to be inferred around the atom's structure.
  ///               This number is decremented by one, meaning that 1 implies no implicit hydrogen
  ///               content.
  int interpretImplicitHydrogenContent(int nh_in);

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
};

/// \brief Overload the + operator to concatenate vectors of MDL and SDF bonds.
std::vector<MolObjBond> operator+(const std::vector<MolObjBond> &lhs,
				  const std::vector<MolObjBond> &rhs);

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
