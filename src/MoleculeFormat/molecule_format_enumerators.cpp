#include <string>
#include "copyright.h"
#include "Parsing/parse.h"
#include "molecule_format_enumerators.h"

namespace stormm {
namespace structure {

using parse::operator==;
  
//-------------------------------------------------------------------------------------------------
MdlMolPropertyKind translateMdlMolPropertyKind(char4 input) {
  if      (input == char4({ ' ', ' ', ' ', 'A' })) return MdlMolPropertyKind::ATOM_ALIAS;
  else if (input == char4({ ' ', ' ', ' ', 'V' })) return MdlMolPropertyKind::ATOM_VALUE;
  else if (input == char4({ ' ', ' ', ' ', 'G' })) return MdlMolPropertyKind::GROUP_ABBREVIATION;
  else if (input == char4({ 'C', 'H', 'G', 'M' })) return MdlMolPropertyKind::CHARGE;
  else if (input == char4({ 'R', 'A', 'D', 'M' })) return MdlMolPropertyKind::RADICAL;
  else if (input == char4({ 'I', 'S', 'O', 'M' })) return MdlMolPropertyKind::ISOTOPE;
  else if (input == char4({ 'R', 'B', 'C', 'M' })) return MdlMolPropertyKind::RING_BOND_COUNT;
  else if (input == char4({ 'S', 'U', 'B', 'M' })) return MdlMolPropertyKind::SUBSTITUTION_COUNT;
  else if (input == char4({ 'U', 'N', 'S', 'M' })) return MdlMolPropertyKind::UNSATURATED_COUNT;
  else if (input == char4({ 'L', 'I', 'N', 'M' })) return MdlMolPropertyKind::LINK_ATOM;
  else if (input == char4({ 'A', 'L', 'S', 'M' })) return MdlMolPropertyKind::ATOM_LIST;
  else if (input == char4({ 'A', 'P', 'O', 'M' })) return MdlMolPropertyKind::ATTACHMENT_POINT;
  else if (input == char4({ 'A', 'A', 'L', 'M' })) return MdlMolPropertyKind::ATTACHMENT_ORDER;
  else if (input == char4({ 'R', 'G', 'P', 'M' })) {
    return MdlMolPropertyKind::RGROUP_LABEL_LOCATION;
  }
  else if (input == char4({ 'L', 'O', 'G', 'M' })) return MdlMolPropertyKind::RGROUP_LOGIC;
  else if (input == char4({ 'S', 'T', 'Y', 'M' })) return MdlMolPropertyKind::SGROUP_TYPE;
  else if (input == char4({ 'S', 'S', 'T', 'M' })) return MdlMolPropertyKind::SGROUP_SUBTYPE;
  else if (input == char4({ 'S', 'L', 'B', 'M' })) return MdlMolPropertyKind::SGROUP_LABELS;
  else if (input == char4({ 'S', 'C', 'N', 'M' })) return MdlMolPropertyKind::SGROUP_CONNECTIVITY;
  else if (input == char4({ 'S', 'D', 'S', 'M' })) return MdlMolPropertyKind::SGROUP_EXPANSION;
  else if (input == char4({ 'S', 'A', 'L', 'M' })) return MdlMolPropertyKind::SGROUP_ATOM_LIST;
  else if (input == char4({ 'S', 'B', 'L', 'M' })) return MdlMolPropertyKind::SGROUP_BOND_LIST;
  else if (input == char4({ 'S', 'P', 'A', 'M' })) return MdlMolPropertyKind::MG_PARENT_ATOM_LIST;
  else if (input == char4({ 'S', 'M', 'T', 'M' })) return MdlMolPropertyKind::SGROUP_SUBSCRIPT;
  else if (input == char4({ 'C', 'R', 'S', 'M' })) return MdlMolPropertyKind::SGROUP_CORRESPONENCE;
  else if (input == char4({ 'S', 'D', 'I', 'M' })) return MdlMolPropertyKind::SGROUP_DISPLAY_INFO;
  else if (input == char4({ 'S', 'B', 'V', 'M' })) return MdlMolPropertyKind::SGROUP_BOND_VECTOR;
  else if (input == char4({ 'S', 'D', 'T', 'M' })) return MdlMolPropertyKind::SGROUP_FIELD;
  else if (input == char4({ 'S', 'D', 'D', 'M' })) return MdlMolPropertyKind::SGROUP_DISPLAY;
  else if (input == char4({ 'S', 'C', 'D', 'M' }) || input == char4({ 'S', 'E', 'D', 'M' })) {
    return MdlMolPropertyKind::SGROUP_DATA;
  }
  else if (input == char4({ 'S', 'P', 'L', 'M' })) return MdlMolPropertyKind::SGROUP_HIERARCHY;
  else if (input == char4({ 'S', 'N', 'C', 'M' })) return MdlMolPropertyKind::SGROUP_COMP_NUMBER;
  else if (input == char4({ '$', '3', 'D', 'M' })) return MdlMolPropertyKind::SPATIAL_FEATURE;
  else if (input == char4({ 'P', 'X', 'A', 'M' })) return MdlMolPropertyKind::PHANTOM_ATOM;
  else if (input == char4({ 'S', 'A', 'P', 'M' })) return MdlMolPropertyKind::SGROUP_ATTACH_POINT;
  else if (input == char4({ 'S', 'C', 'L', 'M' })) return MdlMolPropertyKind::SGROUP_CLASS;
  else if (input == char4({ 'R', 'E', 'G', 'M' })) return MdlMolPropertyKind::LARGE_REGNO;
  else if (input == char4({ 'S', 'B', 'T', 'M' })) return MdlMolPropertyKind::SGROUP_BRACKET_STYLE;
  else if (input == char4({ 'S', 'K', 'P', 'S' })) return MdlMolPropertyKind::SKIP;
  else if (input == char4({ ' ', ' ', ' ', ' ' })) return MdlMolPropertyKind::NONE;
  else {
    const std::string str_code = std::to_string(input.w) + "  " + std::to_string(input.x) +
                                 std::to_string(input.y) + std::to_string(input.z);
    rtErr("The MDL MOL property code " + str_code + " is unrecognized.",
          "translateMdlMolPropertyKind");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const DataRequestKind input) {
  switch (input) {
  case DataRequestKind::STATE_VARIABLE:
    return std::string("STATE_VARIABLE");
  case DataRequestKind::ATOM_INFLUENCES:
    return std::string("ATOM_INLFUENCES");
  case DataRequestKind::TOPOLOGY_PARAMETER:
    return std::string("TOPOLOGY_PARAMETER");
  case DataRequestKind::STRING:
    return std::string("STRING");
  case DataRequestKind::ALL_KINDS:
    return std::string("ALL_KINDS");
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
