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
  else if (input == char4({ 'M', 'C', 'H', 'G' })) return MolObjPropertyKind::CHARGE;
  else if (input == char4({ 'M', 'R', 'A', 'D' })) return MolObjPropertyKind::RADICAL;
  else if (input == char4({ 'M', 'I', 'S', 'O' })) return MolObjPropertyKind::ISOTOPE;
  else if (input == char4({ 'M', 'R', 'B', 'C' })) return MolObjPropertyKind::RING_BOND_COUNT;
  else if (input == char4({ 'M', 'S', 'U', 'B' })) return MolObjPropertyKind::SUBSTITUTION_COUNT;
  else if (input == char4({ 'M', 'U', 'N', 'S' })) return MolObjPropertyKind::UNSATURATED_COUNT;
  else if (input == char4({ 'M', 'L', 'I', 'N' })) return MolObjPropertyKind::LINK_ATOM;
  else if (input == char4({ 'M', 'A', 'L', 'S' })) return MolObjPropertyKind::ATOM_LIST;
  else if (input == char4({ 'M', 'A', 'P', 'O' })) return MolObjPropertyKind::ATTACHMENT_POINT;
  else if (input == char4({ 'M', 'A', 'A', 'L' })) return MolObjPropertyKind::ATTACHMENT_ORDER;
  else if (input == char4({ 'M', 'R', 'G', 'P' })) return MolObjPropertyKind::RGROUP_LABEL_LOCATION;
  else if (input == char4({ 'M', 'L', 'O', 'G' })) return MolObjPropertyKind::RGROUP_LOGIC;
  else if (input == char4({ 'M', 'S', 'T', 'Y' })) return MolObjPropertyKind::SGROUP_TYPE;
  else if (input == char4({ 'M', 'S', 'S', 'T' })) return MolObjPropertyKind::SGROUP_SUBTYPE;
  else if (input == char4({ 'M', 'S', 'L', 'B' })) return MolObjPropertyKind::SGROUP_LABELS;
  else if (input == char4({ 'M', 'S', 'C', 'N' })) return MolObjPropertyKind::SGROUP_CONNECTIVITY;
  else if (input == char4({ 'M', 'S', 'D', 'S' })) return MolObjPropertyKind::SGROUP_EXPANSION;
  else if (input == char4({ 'M', 'S', 'A', 'L' })) return MolObjPropertyKind::SGROUP_ATOM_LIST;
  else if (input == char4({ 'M', 'S', 'B', 'L' })) return MolObjPropertyKind::SGROUP_BOND_LIST;
  else if (input == char4({ 'M', 'S', 'P', 'A' })) return MolObjPropertyKind::MG_PARENT_ATOM_LIST;
  else if (input == char4({ 'M', 'S', 'M', 'T' })) return MolObjPropertyKind::SGROUP_SUBSCRIPT;
  else if (input == char4({ 'M', 'C', 'R', 'S' })) return MolObjPropertyKind::SGROUP_CORRESPONENCE;
  else if (input == char4({ 'M', 'S', 'D', 'I' })) return MolObjPropertyKind::SGROUP_DISPLAY_INFO;
  else if (input == char4({ 'M', 'S', 'B', 'V' })) return MolObjPropertyKind::SGROUP_BOND_VECTOR;
  else if (input == char4({ 'M', 'S', 'D', 'T' })) return MolObjPropertyKind::SGROUP_FIELD;
  else if (input == char4({ 'M', 'S', 'D', 'D' })) return MolObjPropertyKind::SGROUP_DISPLAY;
  else if (input == char4({ 'M', 'S', 'C', 'D' })) return MolObjPropertyKind::SGROUP_DATA;
  else if (input == char4({ 'M', 'S', 'P', 'L' })) return MolObjPropertyKind::SGROUP_HIERARCHY;
  else if (input == char4({ 'M', 'S', 'N', 'C' })) return MolObjPropertyKind::SGROUP_COMP_NUMBER;
  else if (input == char4({ 'M', '$', '3', 'D' })) return MolObjPropertyKind::SPATIAL_FEATURE;
  else if (input == char4({ 'M', 'P', 'X', 'A' })) return MolObjPropertyKind::PHANTOM_ATOM;
  else if (input == char4({ 'M', 'S', 'A', 'P' })) return MolObjPropertyKind::SGROUP_ATTACH_POINT;
  else if (input == char4({ 'M', 'S', 'C', 'L' })) return MolObjPropertyKind::SGROUP_CLASS;
  else if (input == char4({ 'M', 'R', 'E', 'G' })) return MolObjPropertyKind::LARGE_REGNO;
  else if (input == char4({ 'M', 'S', 'B', 'T' })) return MolObjPropertyKind::SGROUP_BRACKET_STYLE;
  else if (input == char4({ 'S', 'S', 'K', 'P' })) return MolObjPropertyKind::SKIP;
  else {
    const std::string str_code = std::to_string(input.w) + "  " + input.x + input.y + input.z;
    rtErr("The MDL MOL property code " + str_code + " is unrecognized.", "translateMolObjPropertyKind");
  }
  __builtin_unreachable();
}

} // namespace structure
} // namespace stormm
