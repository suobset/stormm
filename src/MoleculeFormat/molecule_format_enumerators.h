// -*-c++-*-
#ifndef STORMM_MOLECULE_FORMAT_ENUMERATORS_H
#define STORMM_MOLECULE_FORMAT_ENUMERATORS_H

#include "copyright.h"
#include "DataTypes/stormm_vector_types.h"

namespace stormm {
namespace structure {

/// \brief The possible version numbers of an MDL MOL entry
enum class MdlMolVersion {
  V2000, V3000, UNKNOWN
};

/// \brief Enumerate the types of indices that may be encountered.  The chief purpose of this is
///        to signal that certain indices such as atoms should have 1 subtracted from their values
///        in the original MDL MOL format, to shift into C-array indexing from the original format
///        which is likely based on Fortran.
enum class MolObjIndexKind {
  ATOM,  ///< Subtract 1 from the index found in the file.  Add 1 to the stored index when writing
         ///<   results to a new file.
  BOND,  ///< Subtract or add 1 to read / write this index, like ATOM
  OTHER  ///< No shift in the index is needed
};
  
/// \brief Enumerate the possible radical states for an MDL MOL entry
enum class RadicalState {
  NONE, SINGLET, DOUBLET, TRIPLET
};
  
/// \brief A molecule's chirality
enum class MolObjChirality {
  ACHIRAL = 0, CHIRAL = 1
};

/// \brief Different ways to go about adding hydrogens to a minimal representation of a molecule
enum class HydrogenAssignment {
  VALENCE_SHELL,      ///< Add at least the number of hydrogens given in a secondary array, until
                      ///<   the anticipated valence shell electron content is satisfied.  The
                      ///<   anticipated content is eight electrons for C, N, O, and F, 10 for P,
                      ///<   and 8 or 12 for S depending on the hybridization.
  DO_NOT_HYDROGENATE  ///< Do not attempt to add hydrogens to the atom.
};
  
/// \brief Enumerate possible bond orders, including aromatics or variable bonds
enum class MolObjBondOrder {
  SINGLE = 1, DOUBLE = 2, TRIPLE = 3, AROMATIC = 4, SINGLE_OR_DOUBLE = 5,
  SINGLE_OR_AROMATIC = 6, DOUBLE_OR_AROMATIC = 7, ANY = 8
};

/// \brief Enumerate possible stereochemistries arising from a bond
enum class MolObjBondStereo {
  NOT_STEREO = 0, UP = 1, CIS_OR_TRANS = 3, EITHER = 4, DOWN = 6
};

/// \brief Enumerate different steroe parity settings for an atom
enum class MolObjAtomStereo {
  NOT_STEREO = 0, ODD = 1, EVEN = 2, UNMARKED = 3
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

/// \brief Enumerate the outcomes for stereochemistry during a reaction
enum class StereoRetention {
  NOT_APPLIED = 0, INVERTED = 1, RETAINED = 2
};

/// \brief Enumerate types of data that could be found in the fields on an MDL MOL format property.
enum class MolObjPropField {
  INTEGER, CHAR4, REAL, STRING
};

/// \brief Enumerate the types of data that property fields could contain
enum class MolObjPropertyKind {
  ATOM_ALIAS,
  ATOM_VALUE,
  GROUP_ABBREVIATION,
  CHARGE,
  RADICAL,
  ISOTOPE,
  RING_BOND_COUNT,
  SUBSTITUTION_COUNT,
  UNSATURATED_COUNT,
  LINK_ATOM,
  ATOM_LIST,
  ATTACHMENT_POINT,
  ATTACHMENT_ORDER,
  RGROUP_LABEL_LOCATION,
  RGROUP_LOGIC,
  SGROUP_TYPE,
  SGROUP_SUBTYPE,
  SGROUP_LABELS,
  SGROUP_CONNECTIVITY,
  SGROUP_EXPANSION,
  SGROUP_ATOM_LIST,
  SGROUP_BOND_LIST,
  MG_PARENT_ATOM_LIST,
  SGROUP_SUBSCRIPT,
  SGROUP_CORRESPONENCE,
  SGROUP_DISPLAY_INFO,
  SGROUP_BOND_VECTOR,
  SGROUP_FIELD,
  SGROUP_DISPLAY,
  SGROUP_DATA,
  SGROUP_HIERARCHY,
  SGROUP_COMP_NUMBER,
  SPATIAL_FEATURE,
  PHANTOM_ATOM,
  SGROUP_ATTACH_POINT,
  SGROUP_CLASS,
  LARGE_REGNO,
  SGROUP_BRACKET_STYLE,
  SKIP
};

/// \brief Translate a four-character tuple into one of the known MDL MOL format properties.
///
/// \param input  The code to translate
MolObjPropertyKind translateMolObjPropertyKind(char4 input);

/// \brief Get a string corresponding to the name of the MolObjPropField enumeration.  Overloaded
///        various enumerations as input in this and other libraries.
///
/// \param input  The enumeration to translate
/// \{
std::string getEnumerationName(MolObjIndexKind input);

std::string getEnumerationName(MolObjPropField input);
/// \}

} // namespace structure
} // namespace stormm

#endif
