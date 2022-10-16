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
  ATOM,    ///< Subtract 1 from the index found in the file.  Add 1 to the stored index when
           ///<   writing results to a new file.  S-groups are also counted as ATOM for indexing
           ///<   purposes.
  BOND,    ///< Subtract or add 1 to read / write this index, like ATOM
  OTHER    ///< No shift in the index is needed
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
  VALENCE_SHELL = 0,      ///< Add at least the number of hydrogens given in a secondary array,
                          ///<   until the anticipated valence shell electron content is satisfied.
                          ///<   The anticipated content is eight electrons for C, N, O, and F, 10
                          ///<   for P, and 8 or 12 for S depending on the hybridization.
  DO_NOT_HYDROGENATE = 1  ///< Do not attempt to add hydrogens to the atom.
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

/// \brief Various kinds of information that can be encoded in an SD file data item upon request
enum class DataRequestKind {
  STATE_VARIABLE,      ///< Information from a state variable (see Potential/energy_enumerators.h)
                       ///<   will br presented as a single, real number in a data item.
  ATOM_INFLUENCES,     ///< The origins and effects of all valence interactions, including
                       ///<   restraint terms, affecting the specified atoms will be presented in
                       ///<   a formatted data item.
  TOPOLOGY_PARAMETER,  ///< Information to be extracted from the topology guiding the structure
  STRING               ///< Information in the form of a custom string issued by the user
};

/// \brief Enumerate the types of data that property fields could contain
enum class MolObjPropertyKind {
  ATOM_ALIAS,             ///< Deprecated ISIS / desktop concept
  ATOM_VALUE,             ///< Deprecated ISIS / desktop concept
  GROUP_ABBREVIATION,     ///< Deprecated ISIS / desktop concept
  CHARGE,                 ///< Atom formal charge (one or more properties of this type will
                          ///<   invalidate all charges and radical states imparted through the
                          ///<   block)
  RADICAL,                ///< Atom radical state (one or more properties of this type will
                          ///<   invalidate all charges and radical states imparted through the
                          ///<   block)
  ISOTOPE,                ///< Atom isotopic shift (one or more properties of this type will
                          ///<   invalidate all charges and radical states imparted through the
                          ///<   atom block)
  RING_BOND_COUNT,        ///< The number of ring bonds that an atom may make
  SUBSTITUTION_COUNT,     ///< The number of substitutions allowed for an atom
  UNSATURATED_COUNT,      ///< Indicate that an atom makes at least one double- or triple-bond
  LINK_ATOM,              ///< Designate a link atom
  ATOM_LIST,              ///< Give details about an atom (element) list
  ATTACHMENT_POINT,       ///< Indicate whether an atom is the first, second, or both first and
                          ///<   second attachment point for an R-group
  ATTACHMENT_ORDER,       ///< Orders of attachment bonds for up to two attachment points
                          ///<   connecting an R-group to a particular atom
  RGROUP_LABEL_LOCATION,  ///< Atoms at which one or more R-groups are labeled
  RGROUP_LOGIC,           ///< Details on an R-group
  SGROUP_TYPE,            ///< Type of S-group
  SGROUP_SUBTYPE,         ///< Subtype of one or more S-groups
  SGROUP_LABELS,          ///< Numeric labels for one or more S-groups
  SGROUP_CONNECTIVITY,    ///< Alphanumeric codes for the conenctions between S-groups
  SGROUP_EXPANSION,       ///< S-group indices of expanded abbreviation S-groups
  SGROUP_ATOM_LIST,       ///< Indices of atoms found in a particular S-group (the substrate
                          ///<   is the S-group index, the entry contents are the atom list
                          ///<   indices).  This is not related to ATOM_LIST, which pertains to
                          ///<   elements.
  SGROUP_BOND_LIST,       ///< Indices of bonds in an S-group
  MG_PARENT_ATOM_LIST,    ///< Multiple group parent atom list, with entries being atoms in a
                          ///<   repeating unit of some multiple group 
  SGROUP_SUBSCRIPT,       ///< S-group subscript with text
  SGROUP_CORRESPONENCE,   ///< Bonds that define correspondence between S-groups
  SGROUP_DISPLAY_INFO,    ///< Information to be passed to a graphical program for displaying an
                          ///<   S-group (followed by two types of S-group data, both of which are
                          ///<   tied to the display information property)
  SGROUP_BOND_VECTOR,     ///< A bond connecting to some contracted abbreviation S-group
  SGROUP_FIELD,           ///< Description of an S-group
  SGROUP_DISPLAY,         ///< Codes for graphical display information attached to an S-group 
  SGROUP_DATA,            ///< Information tied to the S-grop display
  SGROUP_HIERARCHY,       ///< Parent-child relationships among S-groups
  SGROUP_COMP_NUMBER,     ///< Indices an dorders of S-group components
  SPATIAL_FEATURE,        ///< A special type of property tracked by the MolObjFeature object
  PHANTOM_ATOM,           ///< A sort of virtual site encoded by the MDL MOL format
  SGROUP_ATTACH_POINT,    ///< Indexing for atoms by which a S-group attaches and atoms that can
                          ///<   leave if the S-group is cleaved off
  SGROUP_CLASS,           ///< The class of an S-group
  LARGE_REGNO,            ///< The registry number
  SGROUP_BRACKET_STYLE,   ///< Style of the S-group bracket
  SKIP                    ///< Indicator that a certain number of lines in the file should be
                          ///<   skipped (the lines will be retained as data lines in the property)
};

/// \brief Enumerate the types of data that property fields could contain
enum class MolObjFeatureKind {
  PT_FD_2,
  PT_FLEX_2,
  PT_PN_NORMAL,
  LN_BEST_FIT,
  PN_BEST_FIT,
  PN_PT_LN,
  CENTROID,
  LN_PN_NORMAL,
  DISTANCE_PT_PT,
  DISTANCE_PT_LN,
  DISTANCE_PT_PN,
  ANGLE_PT_PT_PT,
  ANGLE_LN_LN,
  ANGLE_PN_PN,
  SPHERE_PT,
  FIXED_ATOMS,
  ATOM_CONSTRAINT,
  PAIR_CONSTRAINT
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
std::string getEnumerationName(DataRequestKind input);

std::string getEnumerationName(MolObjIndexKind input);

std::string getEnumerationName(MolObjPropField input);
/// \}

} // namespace structure
} // namespace stormm

#endif
