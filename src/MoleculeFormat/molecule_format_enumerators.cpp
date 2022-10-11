#include <string>
#include "copyright.h"
#include "Parsing/parse.h"
#include "molecule_format_enumerators.h"

namespace stormm {
namespace structure {

using parse::operator==;
  
//-------------------------------------------------------------------------------------------------
MolObjPropertyKind translateMolObjPropertyKind(char4 input) {
  if      (input == char4({ ' ', ' ', ' ', 'A' })) return MolObjPropertyKind::ATOM_ALIAS;
  else if (input == char4({ ' ', ' ', ' ', 'V' })) return MolObjPropertyKind::ATOM_VALUE;
  else if (input == char4({ ' ', ' ', ' ', 'G' })) return MolObjPropertyKind::GROUP_ABBREVIATION;
  else if (input == char4({ 'C', 'H', 'G', 'M' })) return MolObjPropertyKind::CHARGE;
  else if (input == char4({ 'R', 'A', 'D', 'M' })) return MolObjPropertyKind::RADICAL;
  else if (input == char4({ 'I', 'S', 'O', 'M' })) return MolObjPropertyKind::ISOTOPE;
  else if (input == char4({ 'R', 'B', 'C', 'M' })) return MolObjPropertyKind::RING_BOND_COUNT;
  else if (input == char4({ 'S', 'U', 'B', 'M' })) return MolObjPropertyKind::SUBSTITUTION_COUNT;
  else if (input == char4({ 'U', 'N', 'S', 'M' })) return MolObjPropertyKind::UNSATURATED_COUNT;
  else if (input == char4({ 'L', 'I', 'N', 'M' })) return MolObjPropertyKind::LINK_ATOM;
  else if (input == char4({ 'A', 'L', 'S', 'M' })) return MolObjPropertyKind::ATOM_LIST;
  else if (input == char4({ 'A', 'P', 'O', 'M' })) return MolObjPropertyKind::ATTACHMENT_POINT;
  else if (input == char4({ 'A', 'A', 'L', 'M' })) return MolObjPropertyKind::ATTACHMENT_ORDER;
  else if (input == char4({ 'R', 'G', 'P', 'M' })) {
    return MolObjPropertyKind::RGROUP_LABEL_LOCATION;
  }
  else if (input == char4({ 'L', 'O', 'G', 'M' })) return MolObjPropertyKind::RGROUP_LOGIC;
  else if (input == char4({ 'S', 'T', 'Y', 'M' })) return MolObjPropertyKind::SGROUP_TYPE;
  else if (input == char4({ 'S', 'S', 'T', 'M' })) return MolObjPropertyKind::SGROUP_SUBTYPE;
  else if (input == char4({ 'S', 'L', 'B', 'M' })) return MolObjPropertyKind::SGROUP_LABELS;
  else if (input == char4({ 'S', 'C', 'N', 'M' })) return MolObjPropertyKind::SGROUP_CONNECTIVITY;
  else if (input == char4({ 'S', 'D', 'S', 'M' })) return MolObjPropertyKind::SGROUP_EXPANSION;
  else if (input == char4({ 'S', 'A', 'L', 'M' })) return MolObjPropertyKind::SGROUP_ATOM_LIST;
  else if (input == char4({ 'S', 'B', 'L', 'M' })) return MolObjPropertyKind::SGROUP_BOND_LIST;
  else if (input == char4({ 'S', 'P', 'A', 'M' })) return MolObjPropertyKind::MG_PARENT_ATOM_LIST;
  else if (input == char4({ 'S', 'M', 'T', 'M' })) return MolObjPropertyKind::SGROUP_SUBSCRIPT;
  else if (input == char4({ 'C', 'R', 'S', 'M' })) return MolObjPropertyKind::SGROUP_CORRESPONENCE;
  else if (input == char4({ 'S', 'D', 'I', 'M' })) return MolObjPropertyKind::SGROUP_DISPLAY_INFO;
  else if (input == char4({ 'S', 'B', 'V', 'M' })) return MolObjPropertyKind::SGROUP_BOND_VECTOR;
  else if (input == char4({ 'S', 'D', 'T', 'M' })) return MolObjPropertyKind::SGROUP_FIELD;
  else if (input == char4({ 'S', 'D', 'D', 'M' })) return MolObjPropertyKind::SGROUP_DISPLAY;
  else if (input == char4({ 'S', 'C', 'D', 'M' }) || input == char4({ 'S', 'E', 'D', 'M' })) {
    return MolObjPropertyKind::SGROUP_DATA;
  }
  else if (input == char4({ 'S', 'P', 'L', 'M' })) return MolObjPropertyKind::SGROUP_HIERARCHY;
  else if (input == char4({ 'S', 'N', 'C', 'M' })) return MolObjPropertyKind::SGROUP_COMP_NUMBER;
  else if (input == char4({ '$', '3', 'D', 'M' })) return MolObjPropertyKind::SPATIAL_FEATURE;
  else if (input == char4({ 'P', 'X', 'A', 'M' })) return MolObjPropertyKind::PHANTOM_ATOM;
  else if (input == char4({ 'S', 'A', 'P', 'M' })) return MolObjPropertyKind::SGROUP_ATTACH_POINT;
  else if (input == char4({ 'S', 'C', 'L', 'M' })) return MolObjPropertyKind::SGROUP_CLASS;
  else if (input == char4({ 'R', 'E', 'G', 'M' })) return MolObjPropertyKind::LARGE_REGNO;
  else if (input == char4({ 'S', 'B', 'T', 'M' })) return MolObjPropertyKind::SGROUP_BRACKET_STYLE;
  else if (input == char4({ 'S', 'K', 'P', 'S' })) return MolObjPropertyKind::SKIP;
  else {
    const std::string str_code = std::to_string(input.w) + "  " + input.x + input.y + input.z;
    rtErr("The MDL MOL property code " + str_code + " is unrecognized.",
          "translateMolObjPropertyKind");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const MolObjIndexKind input) {
  switch (input) {
  case MolObjIndexKind::ATOM:
    return std::string("ATOM");
  case MolObjIndexKind::BOND:
    return std::string("BOND");
  case MolObjIndexKind::OTHER:
    return std::string("OTHER");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const MolObjPropField input) {
  switch (input) {
  case MolObjPropField::INTEGER:
    return std::string("INTEGER");
  case MolObjPropField::REAL:
    return std::string("REAL");
  case MolObjPropField::CHAR4:
    return std::string("CHAR4");
  case MolObjPropField::STRING:
    return std::string("STRING");
  }
  __builtin_unreachable();
}

} // namespace structure
} // namespace stormm
