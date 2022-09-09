// -*-c++-*-
#ifndef STORMM_MDL_FILE_H
#define STORMM_MDL_FILE_H

#include <string>
#include <vector>
#include "copyright.h"
#include "DataTypes/stormm_vector_types.h"
#include "Parsing/parse.h"

namespace stormm {
namespace structure {

using parse::TextFile;

/// \brief A molecule's chirality
enum class MolObjChirality {
  ACHIRAL = 0, CHIRAL = 1
};

/// \brief Enumerate possible bond orders, including aromatics or variable bonds
enum class MolObjBondOrder {
  SINGLE = 1, DOUBLE = 2, TRIPLE = 3, AROMATIC = 4, SINGLE_OR_DOUBLE = 5,
  SINGLE_OR_AROMATIC = 6, DOUBLE_OR_AROMATIC = 7, ANY = 8
};

/// \brief Enumerate possible stereochemistries arising from a bond
enum class MolObjStereo {
  NOT_STEREO = 0, UP = 1, CIS_OR_TRANS = 3, EITHER = 4, DOWN = 6
};

/// \brief Enumerate states in which a bond participates in a ring system or just a chain
enum class MolObjRingState {
  EITHER = 0, RING = 1, CHAIN = 2
};

/// \brief Enumerate ways in which a _bond_ can participate in reactions
enum class MolObjReactionCenter {
  NON_CENTER = -1, UNMARKED = 0, CENTER = 1, UNREACTIVE = 2, BOND_MADE_OR_BROKEN = 4,
  CENTER_WITH_FORMATION = 5, BOND_ORDER_CHANGE = 8, CENTER_WITH_ORDER_CHANGE = 9,
  BOND_FORMATION_AND_ORDER_CHANGE = 12, CENTER_WITH_FORMATION_AND_ORDER_CHANGE = 13
};

/// \brief A bond, as presented in the MDL molecule file format.  This unguarded struct will be
///        returned to the developer from a private array inside of the MdlMolObj object, so
///        further protection would be a hindrance.
struct MolObjBond {
  int i_atom;                       ///< The first atom in the bond
  int j_atom;                       ///< The second atom in the bond
  MolObjBondOrder order;            ///< The bond order (single, double, aromatic, etc.)
  MolObjStereo stereo;              ///< Indicator of the bond stereochemistry
  MolObjRingState ring_state;       ///< Indicator of whether the atom is part of a ring
  MolObjReactionCenter reactivity;  ///< Indicator of a bond as a center of reactivity
  int xxx_data;                     ///< Currently unused
};

/// \brief A molecular or atomic property read from an MDL .mol or SDF file
struct MolObjProperty {
  
};

/// \brief A molecule read from an MDL .mol file, or one of many read from a concatenated SDF file.
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
  /// \param filename    Name of the file to read, containing one MDL file
  /// \param tf          Text file data previously read into memory (this is the way to handle an
  ///                    SDF file)
  /// \param line_start  The first line at which to begin reading the TextFile object
  /// \param line_end    The last relevant line of the TextFile object at which reading will stop
  /// \{
  MdlMolObj();
  MdlMolObj(const std::string &filename);
  MdlMolObj(const char* filename);
  MdlMolObj(const TextFile &tf, int line_start = 0, int line_end = -1);
  /// \}
  
  // Getter functions, some of which trigger complex filtering behavior
  double3 getCoordinate(int index) const;
  std::vector<double3> getCoordinate() const;
  std::vector<double3> getCoordinates() const;
  int getAtomicNumber(int index) const;
  std::vector<int> getAtomicNumber() const;
  std::vector<int> getAtomicNumbers(int index) const;
  std::vector<MolObjBond> getBonds() const;
  std::vector<MolObjBond> getRotatableBonds() const;
  std::vector<MolObjBond> getBonds(int atom_index) const;
  std::vector<MolObjBond> getBonds(MolObjBondOrder query_order) const;
  std::vector<MolObjProperty> getProperties() const;

private:

  /// The number of atoms in the molecule
  int atom_count;

  /// The number of bonds of all types between atoms in the system
  int bond_count;

  /// The number of S-groups
  int sgroup_count;

  /// The number of three-dimensional constraints
  int constraint_count;
 
  /// The molecule's chirality (assumes only one significant center)
  MolObjChirality chirality;

  /// The molecule registry number
  int registry_number;
 
  /// Cartesian coordinates of all atoms
  std::vector<double3> coordinates;

  /// Symbols for all atoms
  std::vector<char4> atomic_symbols;

  /// Atomic numbers for all atoms
  std::vector<int> atomic_numbers;

  /// Formal charges for all atoms
  std::vector<double> formal_charges;

  /// Isotope numbers shifting each atoms nuclear mass ~1 Dalton from the most common isotopic mass
  std::vector<int> isotopic_shifts;

  /// Bonds between atoms
  std::vector<MolObjBond> bonds;

  /// Properties
  std::vector<MolObjProperty> properties;

  /// Title (first line from the file)
  std::string title;

  /// Details of the software that generated this MDL .mol file
  std::string software_details;

  /// General comments for the file
  std::string general_comment;
};

/// \brief Overload the + operator to concatenate vectors of MDL and SDF bonds.
std::vector<MolObjBond> operator+(const std::vector<MolObjBond> &lhs,
				  const std::vector<MolObjBond> &rhs);

} // namespace structure
} // namespace stormm

#endif
