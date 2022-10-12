// -*-c++-*-
#ifndef STORMM_MDL_FILE_H
#define STORMM_MDL_FILE_H

#include <fstream>
#include <string>
#include <vector>
#include "copyright.h"
#include "Chemistry/znumber.h"
#include "Constants/behavior.h"
#include "DataTypes/stormm_vector_types.h"
#include "FileManagement/file_util.h"
#include "Parsing/ascii_numbers.h"
#include "Parsing/parse.h"
#include "Trajectory/coordinateframe.h"
#include "Trajectory/phasespace.h"
#include "molecule_file_io.h"
#include "molecule_format_enumerators.h"
#include "molobj_atomlist.h"
#include "molobj_bond.h"
#include "molobj_dataitem.h"
#include "molobj_property.h"

namespace stormm {
namespace structure {

using chemistry::symbolToZNumber;
using constants::CartesianDimension;
using constants::CaseSensitivity;
using constants::ExceptionResponse;
using diskutil::PrintSituation;
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

  /// \brief Write a set of molecular coordinates in MDL MOL format.  Apply all properties already
  ///        stored in the object, in a 1:1 reversal of the operations for reading such a file.
  ///
  /// Overloaded:
  ///   - Write to an output file
  ///   - Write to a string
  ///
  /// \param foutp        Output file stream
  /// \param fname        Name of the output file
  /// \param vformat      Format version of the MDL MOL standard to use.  If provided as either
  ///                     V2000 or V3000, this will override information stored in the object.
  /// \param expectation  Anticipated (or required) condition of the output file that is to be
  ///                     opened, if only its name is provided
  /// \{
  void write(std::ofstream *foutp, MdlMolVersion vformat, PrintSituation expectation) const;
  void write(const std::string &fname, MdlMolVersion vformat, PrintSituation expectation) const;
  std::string write(MdlMolVersion vformat) const;
  /// \}

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
