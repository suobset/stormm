#include "copyright.h"
#include "Constants/behavior.h"
#include "FileManagement/file_listing.h"
#include "Parsing/parse.h"
#include "molobj_dataitem.h"

namespace stormm {
namespace structure {

using constants::CaseSensitivity;
using diskutil::getBaseName;
using parse::NumberFormat;
using parse::strcmpCased;
using parse::stringToChar4;
using parse::strncmpCased;
using parse::separateText;
using parse::verifyContents;

//-------------------------------------------------------------------------------------------------
MolObjDataRequest::MolObjDataRequest(const std::string &title_in, const std::string &label_in) :
    kind{DataRequestKind::STRING}, title{title_in}, energy_component{StateVariable::BOND},
    atom_mask{std::string("")}, valence_kind{StateVariable::BOND}, message{std::string("")},
    atom_types{}, system_label{label_in}, use_maccs_ii_number{false}, maccs_ii_number{0},
    use_internal_registry{false}, external_regno{std::string("")}
{}

//-------------------------------------------------------------------------------------------------
MolObjDataRequest::MolObjDataRequest(const std::string &title_in,
                                     const StateVariable energy_component_in,
                                     const std::string &label_in) :
    MolObjDataRequest(title_in, label_in)
{
  kind = DataRequestKind::STATE_VARIABLE;
  energy_component = energy_component_in;
}

//-------------------------------------------------------------------------------------------------
MolObjDataRequest::MolObjDataRequest(const DataRequestKind kind_in, const std::string &title_in,
                                     const std::string &message_in, const std::string &label_in) :
    MolObjDataRequest(title_in, label_in)
{
  kind = kind_in;
  switch (kind) {
  case DataRequestKind::STATE_VARIABLE:
    rtErr("In order to construct a data request for a particular energy component, use "
          "MolObjDataRequest(<title>, <energy component>, <label>).", "MolObjDataRequest");
    break;
  case DataRequestKind::TOPOLOGY_PARAMETER:
    rtErr("In order to construct a data request for the actions of a topological energy "
          "parameter, use MolObjDataRequest(<title>, <interaction type>, <list of atom types>, "
          "<label>).", "MolObjDataRequest");
    break;
  case DataRequestKind::ATOM_INFLUENCES:
    atom_mask = message_in;
    break;
  case DataRequestKind::STRING:
    message = message_in;
    break;
  case DataRequestKind::ALL_KINDS:
    rtErr("The ALL_KINDS enumeration does not indicate an actual data request and no program "
          "should invoke it here.", "MolObjDataRequest");
    break;
  }
}

//-------------------------------------------------------------------------------------------------
MolObjDataRequest::MolObjDataRequest(const std::string &title_in,
                                     const StateVariable valence_kind_in,
                                     const std::vector<char4> &atom_types_in,
                                     const std::string &label_in) :
    MolObjDataRequest(title_in, label_in)
{
  kind = DataRequestKind::TOPOLOGY_PARAMETER;
  valence_kind = valence_kind_in;
  atom_types = atom_types_in;

  // Check the type of valence interaction and the corresponding number of atom types supplied.
  const int ntypes = atom_types.size();
  int nreq;
  switch (valence_kind) {
  case StateVariable::BOND:
  case StateVariable::UREY_BRADLEY:
    nreq = 2;
    break;
  case StateVariable::ANGLE:
    nreq = 3;
    break;
  case StateVariable::PROPER_DIHEDRAL:
  case StateVariable::IMPROPER_DIHEDRAL:
  case StateVariable::CHARMM_IMPROPER:
    nreq = 4;
    break;
  case StateVariable::CMAP:
    nreq = 5;
    break;
  case StateVariable::VDW:
  case StateVariable::VDW_ONE_FOUR:
  case StateVariable::ELECTROSTATIC:
  case StateVariable::ELECTROSTATIC_ONE_FOUR:
  case StateVariable::GENERALIZED_BORN:
  case StateVariable::RESTRAINT:
  case StateVariable::KINETIC:
  case StateVariable::PRESSURE:
  case StateVariable::VIRIAL_11:
  case StateVariable::VIRIAL_12:
  case StateVariable::VIRIAL_22:
  case StateVariable::VIRIAL_13:
  case StateVariable::VIRIAL_23:
  case StateVariable::VIRIAL_33:
  case StateVariable::VOLUME:
  case StateVariable::TEMPERATURE_ALL:
  case StateVariable::TEMPERATURE_PROTEIN:
  case StateVariable::TEMPERATURE_LIGAND:
  case StateVariable::TEMPERATURE_SOLVENT:
  case StateVariable::DU_DLAMBDA:
  case StateVariable::POTENTIAL_ENERGY:
  case StateVariable::TOTAL_ENERGY:
  case StateVariable::ALL_STATES:
    rtErr("The accepted topology parameter types for printing in data items of an SD file "
          "include BOND, ANGLE, PROPER_DIHEDRAL, IMPROPER_DIHEDRAL, UREY_BRADLEY, "
          "CHARMM_IMPROPER, and CMAP.", "MolObjDataRequest");
    break;
  }
  if (ntypes != nreq) {
    rtErr("A " + getEnumerationName(valence_kind) + " requires specification of " +
          std::to_string(nreq) + " atom types, but " + std::to_string(ntypes) + " were provided.",
          "MolObjDataRequest");
  }
}

//-------------------------------------------------------------------------------------------------
DataRequestKind MolObjDataRequest::getKind() const {
  return kind;
}

//-------------------------------------------------------------------------------------------------
const std::string& MolObjDataRequest::getTitle() const {
  return title;
}

//-------------------------------------------------------------------------------------------------
StateVariable MolObjDataRequest::getEnergyComponent() const {
  checkKind(DataRequestKind::STATE_VARIABLE);
  return energy_component;
}

//-------------------------------------------------------------------------------------------------
const std::string& MolObjDataRequest::getAtomMask() const {
  checkKind(DataRequestKind::ATOM_INFLUENCES);
  return atom_mask;
}

//-------------------------------------------------------------------------------------------------
StateVariable MolObjDataRequest::getValenceParameter() const {
  checkKind(DataRequestKind::TOPOLOGY_PARAMETER);
  return valence_kind;
}

//-------------------------------------------------------------------------------------------------
const std::string& MolObjDataRequest::getMessage() const {
  checkKind(DataRequestKind::STRING);
  return message;
}

//-------------------------------------------------------------------------------------------------
const std::vector<char4>& MolObjDataRequest::getAtomTypes() const {
  return atom_types;
}

//-------------------------------------------------------------------------------------------------
const std::string& MolObjDataRequest::getSystemLabel() const {
  return system_label;
}

//-------------------------------------------------------------------------------------------------
const std::string& MolObjDataRequest::getExternalRegistryNumber() const {
  return external_regno;
}

//-------------------------------------------------------------------------------------------------
bool MolObjDataRequest::placeMaccsFieldInHeader() const {
  return use_maccs_ii_number;
}

//-------------------------------------------------------------------------------------------------
bool MolObjDataRequest::placeInternalRegistryInHeader() const {
  return use_internal_registry;
}

//-------------------------------------------------------------------------------------------------
int MolObjDataRequest::getMaccsFieldNumber() const {
  return maccs_ii_number;
}

//-------------------------------------------------------------------------------------------------
void MolObjDataRequest::setExternalRegistryNumber(const std::string &regno_in) {
  external_regno = regno_in;
}

//-------------------------------------------------------------------------------------------------
void MolObjDataRequest::setMaccsFieldNumber(int maccs_in) {
  maccs_ii_number = maccs_in;
  use_maccs_ii_number = true;
}

//-------------------------------------------------------------------------------------------------
void MolObjDataRequest::setInternalRegistryUsage(const std::string &input) {
  if (strcmpCased(input, "on", CaseSensitivity::NO) ||
      strcmpCased(input, "active", CaseSensitivity::NO)) {
    use_internal_registry = true;
  }
  else if (strcmpCased(input, "off", CaseSensitivity::NO)) {
    use_internal_registry = false;
  }
  else {
    rtErr("Invalid directive " + input + " to an SD file data item request, pertaining to "
          "internal registry number usage.  Use ON or OFF.", "MolObjDataRequest",
          "setInternalRegistryUsage");
  }
}

//-------------------------------------------------------------------------------------------------
void MolObjDataRequest::checkKind(const DataRequestKind accepted_kind) const {
  if (kind != accepted_kind) {
    rtErr("A data request of type " + getEnumerationName(kind) + " cannot function as a request "
          "for " + getEnumerationName(accepted_kind) + " data.", "MolObjDataRequest",
          "checkKind");
  }
}

//-------------------------------------------------------------------------------------------------
MolObjDataItem::MolObjDataItem(const std::string &item_name_in,
                               const std::string &external_regno_in,
                               const int internal_regno_in, const int maccs_ii_number_in,
                               const uint header_info, const std::vector<std::string> &body_in) :
    item_name{item_name_in}, external_regno{external_regno_in}, internal_regno{internal_regno_in},
    maccs_ii_number{maccs_ii_number_in}, use_internal_regno{((header_info & 0x1) == 1U)},
    use_external_regno{((header_info & 0x2) == 2U)}, use_item_name{((header_info & 0x4) == 4U)},
    use_maccs_ii_number{((header_info & 0x8) == 8U)}, note_archives{((header_info & 0x10) == 16U)},
    body{body_in}
{
  // Check the sanity of the new object
  validateItemName();
}

//-------------------------------------------------------------------------------------------------
MolObjDataItem::MolObjDataItem(const TextFile &tf, const int line_number, int *line_advance,
                               const int compound_line_end, const std::string &title) :
    MolObjDataItem()
{
  // Find an item name in the header line
  const char* line_ptr = tf.getLinePointer(line_number);
  const int lnlength = tf.getLineLength(line_number);
  bool on_item_name = false;
  bool on_regno = false;
  for (int i = 1; i < lnlength; i++) {
    if (line_ptr[i] == '<') {
      if (on_item_name) {
        rtErr("The reserved character '<' appears twice on line " + std::to_string(line_number) +
              " of file " + getBaseName(tf.getFileName()) + ", compound " + title + ".",
              "MolObjDataItem");
      }
      on_item_name = true;
      use_item_name = true;
    }
    else if (line_ptr[i] == '>') {
      if (on_item_name == false) {
        rtErr("The reserved character '>' appears before its counterpart '<' on line " +
              std::to_string(line_number) + " of file " + getBaseName(tf.getFileName()) + ", "
              "failing to define a data item name in compound " + title + ".", "MolObjDataItem");
      }
      on_item_name = false;
    }
    else if (line_ptr[i] == '(') {
      if (on_regno) {
        rtErr("The reserved character '(' appears twice on line " + std::to_string(line_number) +
              " of file " + getBaseName(tf.getFileName()) + ", in compound " + title + ".",
              "MolObjDataItem");
      }
      on_regno = true;
      use_external_regno = true;
    }
    else if (line_ptr[i] == ')') {
      if (on_regno == false) {
        rtErr("The reserved character ')' appears before its counterpart '(' on line " +
              std::to_string(line_number) + " of file " + getBaseName(tf.getFileName()) + ", "
              "failing to define an external registry number in compound " + title + ".",
              "MolObjDataItem");
      }
      on_regno = false;
    }
    else if (on_item_name == false && on_regno == false) {
      if (line_ptr[i] == 'D' && i < lnlength - 2 && line_ptr[i + 1] == 'T') {
        maccs_ii_number = 0;
        int j = i + 2;
        while (j < lnlength && line_ptr[j] >= '0' && line_ptr[j] <= '9') {
          maccs_ii_number *= 10;
          maccs_ii_number += static_cast<int>(line_ptr[j]) - static_cast<int>('0');
          j++;
        }
        i = j - 1;
        use_maccs_ii_number = true;
      }
      else if (line_ptr[i] == 'F' && i <= lnlength - 13 &&
               strncmpCased(std::string("FROM ARCHIVES"), &line_ptr[i])) {
        note_archives = true;
        i += 12;
      }
      else if (line_ptr[i] >= '0' && line_ptr[i] <= '9') {
        internal_regno = 0;
        int j = i;
        while (j < lnlength && line_ptr[j] >= '0' && line_ptr[j] <= '9') {
          internal_regno *= 10;
          internal_regno += static_cast<int>(line_ptr[j]) - static_cast<int>('0');
          j++;
        }
        i = j - 1;
        use_internal_regno = true;
      }
    }
    else {
      if (on_item_name) {
        item_name += line_ptr[i];
      }
      else if (on_regno) {
        external_regno += line_ptr[i];
      }
    }
  }

  // Validate the header line information
  validateItemName();

  // Read the data lines
  int tmp_advance = line_number + 1;
  bool search_on = true;
  const int actual_compound_line_end = (compound_line_end == -1) ? tf.getLineCount() :
                                                                   compound_line_end;
  while (tmp_advance < actual_compound_line_end && search_on) {
    const int dl_length = tf.getLineLength(tmp_advance);
    if (dl_length == 0) {
      search_on = false;
    }
    else {
      search_on = false;
      const char* dl_ptr = tf.getLinePointer(tmp_advance);
      for (int i = 0; i < dl_length; i++) {
        search_on = (search_on || dl_ptr[i] != ' ');
      }
    }
    if (search_on) {
      tmp_advance++;
    }
  }
  tmp_advance--;
  body.reserve(tmp_advance - line_number);
  *line_advance = tmp_advance;
  for (int pos = line_number + 1; pos <= tmp_advance; pos++) {
    body.push_back(tf.extractString(pos));
  }
}

//-------------------------------------------------------------------------------------------------
MolObjDataItem::MolObjDataItem(const MolObjDataRequest &ask,
                               const std::vector<std::string> &body_in) :
    MolObjDataItem(ask.getTitle(), ask.getExternalRegistryNumber(), -1, ask.getMaccsFieldNumber(),
                   getDataItemHeaderCode(ask), body_in)
{}

//-------------------------------------------------------------------------------------------------
const std::string& MolObjDataItem::getItemName() const {
  return item_name;
}

//-------------------------------------------------------------------------------------------------
const std::string& MolObjDataItem::getExternalRegistryNumber() const {
  return external_regno;
}

//-------------------------------------------------------------------------------------------------
int MolObjDataItem::getInternalRegistryNumber() const {
  return internal_regno;
}

//-------------------------------------------------------------------------------------------------
int MolObjDataItem::getMaccsFieldNumber() const {
  return maccs_ii_number;
}

//-------------------------------------------------------------------------------------------------
std::string MolObjDataItem::parseString(const int element_number, const int line_number) const {
  if (line_number >= static_cast<int>(body.size())) {
    const std::string identifier = (item_name.size() > 0LLU) ? "Item name : <" + item_name + ">" :
                                                               std::string("");
    rtErr("Line index " + std::to_string(line_number) + " is inaccessible in a data item with " +
          std::to_string(body.size()) + " data lines.  " + identifier + ".", "MolObjDataItem",
          "parseString");
  }
  const std::vector<std::string> line_components = separateText(body[line_number]);
  if (element_number < static_cast<int>(line_components.size()) || element_number < 0) {
    const std::string identifier = (item_name.size() > 0LLU) ? "Item name : <" + item_name + ">" :
                                                               std::string("");
    rtErr("Line index " + std::to_string(line_number) + " has " +
          std::to_string(line_components.size()) + " components.  Index " +
          std::to_string(element_number) + " is invalid.", "MolObjDataItem", "parseString");
  }
  return line_components[element_number];
}

//-------------------------------------------------------------------------------------------------
std::string MolObjDataItem::parseString(const int start_pos, const int length,
                                        const int line_number) const {
  if (line_number >= static_cast<int>(body.size())) {
    const std::string identifier = (item_name.size() > 0LLU) ? "Item name : <" + item_name + ">" :
                                                               std::string("");
    rtErr("Line index " + std::to_string(line_number) + " is inaccessible in a data item with " +
          std::to_string(body.size()) + " data lines.  " + identifier + ".", "MolObjDataItem",
          "parseString");
  }
  if (start_pos < 0) {
    rtErr("Starting position " + std::to_string(start_pos) + " is invalid.", "MolObjDataItem",
          "parseString");
  }
  const int lnlength = body[line_number].size();
  const int actual_length = (length < 0) ? lnlength - start_pos : length;
  if (start_pos + length > lnlength) {
    rtErr("Starting position " + std::to_string(start_pos) + " and length " +
          std::to_string(actual_length) + " combine to make an invalid read of a string with "
          "length " + std::to_string(body[line_number].size()) + " characters.", "MolObjDataItem",
          "parseString");
  }
  return body[line_number].substr(start_pos, actual_length);
}

//-------------------------------------------------------------------------------------------------
llint MolObjDataItem::parseInteger(const int element_number, const int line_number) const {
  const std::string proto = parseString(element_number, line_number);
  if (verifyContents(proto, NumberFormat::LONG_LONG_INTEGER)) {
    return stoll(proto);
  }
  else {
    rtErr("Element " + std::to_string(element_number) + " of line " + std::to_string(line_number) +
          ", " + proto + ", could not be parsed as an integer.", "MolObjDataItem", "parseString");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
llint MolObjDataItem::parseInteger(const int start_pos, const int length,
                                   const int line_number) const {
  const std::string proto = parseString(start_pos, length, line_number);
  if (verifyContents(proto, NumberFormat::LONG_LONG_INTEGER)) {
    return stoll(proto);
  }
  else {
    rtErr("Column-formatted characters \"" + proto + "\" at position " +
          std::to_string(start_pos) + " of line " + std::to_string(line_number) + " could not be "
          "parsed as an integer.", "MolObjDataItem", "parseString");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
ullint MolObjDataItem::parseUnsigned(const int element_number, const int line_number) const {
  const std::string proto = parseString(element_number, line_number);
  if (verifyContents(proto, NumberFormat::UNSIGNED_LONG_LONG_INTEGER)) {
    return stoull(proto);
  }
  else {
    rtErr("Element " + std::to_string(element_number) + " of line " + std::to_string(line_number) +
          ", " + proto + ", could not be parsed as an unsigned integer.", "MolObjDataItem",
          "parseString");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
ullint MolObjDataItem::parseUnsigned(const int start_pos, const int length,
                                     const int line_number) const {
  const std::string proto = parseString(start_pos, length, line_number);
  if (verifyContents(proto, NumberFormat::UNSIGNED_LONG_LONG_INTEGER)) {
    return stoll(proto);
  }
  else {
    rtErr("Column-formatted characters \"" + proto + "\" at position " +
          std::to_string(start_pos) + " of line " + std::to_string(line_number) + " could not be "
          "parsed as an unsigned integer.", "MolObjDataItem", "parseString");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
double MolObjDataItem::parseReal(const int element_number, const int line_number) const {
  const std::string proto = parseString(element_number, line_number);
  if (verifyContents(proto, NumberFormat::STANDARD_REAL) ||
      verifyContents(proto, NumberFormat::SCIENTIFIC)) {
    return stod(proto);
  }
  else {
    rtErr("Element " + std::to_string(element_number) + " of line " + std::to_string(line_number) +
          ", " + proto + ", could not be parsed as a real number.", "MolObjDataItem",
          "parseString");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
double MolObjDataItem::parseReal(const int start_pos, const int length,
                                 const int line_number) const {
  const std::string proto = parseString(start_pos, length, line_number);
  if (verifyContents(proto, NumberFormat::STANDARD_REAL) ||
      verifyContents(proto, NumberFormat::SCIENTIFIC)) {
    return stod(proto);
  }
  else {
    rtErr("Column-formatted characters \"" + proto + "\" at position " +
          std::to_string(start_pos) + " of line " + std::to_string(line_number) + " could not be "
          "parsed as a real number.", "MolObjDataItem", "parseString");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
char4 MolObjDataItem::parseChar4(const int element_number, const int line_number) const {
  const std::string proto = parseString(element_number, line_number);
  if (verifyContents(proto, NumberFormat::CHAR4)) {
    return stringToChar4(proto);
  }
  else {
    rtErr("Element " + std::to_string(element_number) + " of line " + std::to_string(line_number) +
          ", " + proto + ", could not be parsed as a four-character tuple, having " +
          std::to_string(proto.size()) + " characters in all.", "MolObjDataItem",
          "parseString");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
char4 MolObjDataItem::parseChar4(const int start_pos, const int length,
                                 const int line_number) const {
  const std::string proto = parseString(start_pos, length, line_number);
  if (verifyContents(proto, NumberFormat::CHAR4)) {
    return stringToChar4(proto);
  }
  else {
    rtErr("Column-formatted characters \"" + proto + "\" at position " +
          std::to_string(start_pos) + " of line " + std::to_string(line_number) + " could not be "
          "parsed as a four-character tuple.", "MolObjDataItem", "parseString");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
bool MolObjDataItem::matchItemName(const std::string &item_name_comp,
                                   const std::string &ext_regno_comp,
                                   const int maccs_ii_no_comp) const {
  return ((use_item_name && item_name_comp == item_name) &&
          (use_external_regno && ext_regno_comp == external_regno) &&
          (use_maccs_ii_number && maccs_ii_no_comp == maccs_ii_number));
}

//-------------------------------------------------------------------------------------------------
bool MolObjDataItem::matchItemName(const std::string &item_name_comp,
                                   const int maccs_ii_no_comp) const {
  return ((use_item_name && item_name_comp == item_name) &&
          (use_maccs_ii_number && maccs_ii_no_comp == maccs_ii_number));
}

//-------------------------------------------------------------------------------------------------
bool MolObjDataItem::matchRegistryNumber(const std::string &ext_regno_comp,
                                         const int maccs_ii_no_comp) const {
  return ((use_external_regno && ext_regno_comp == external_regno) &&
          (use_maccs_ii_number && maccs_ii_no_comp == maccs_ii_number));
}

//-------------------------------------------------------------------------------------------------
bool MolObjDataItem::matchMaccsField(const int maccs_ii_no_comp) const {
  return (use_maccs_ii_number && maccs_ii_no_comp == maccs_ii_number);
}

//-------------------------------------------------------------------------------------------------
void MolObjDataItem::setItemName(const std::string &item_name_in) {
  item_name = item_name_in;
  validateItemName();
}

//-------------------------------------------------------------------------------------------------
void MolObjDataItem::validateItemName() const {
  const int nchar = item_name.size();
  bool problem = false;
  if (nchar > 0) {
    if ((item_name[0] >= 'a' && item_name[0] <= 'z') ||
        (item_name[0] >= 'A' && item_name[0] <= 'Z') || item_name[0] == '_') {
      for (int i = 1; i < nchar; i++) {
        problem = (problem ||
                   item_name[i] == '-' || item_name[i] == '.' || item_name[i] == '<' ||
                   item_name[i] == '>' || item_name[i] == '=' || item_name[i] == '%' ||
                   item_name[i] == ' ');
      }
    }
    else {
      problem = true;
    }
  }
  if (problem) {
    rtErr("Data item name " + item_name + " is invalid.  An item name must begin with an "
          "alphabetical character and thereafter contain alphanumeric characters and "
          "underscores, with no white space.", "MolObjDataItem", "setItemName");
  }
}

//-------------------------------------------------------------------------------------------------
uint getDataItemHeaderCode(const bool use_internal_regno, const bool use_external_regno,
                           const bool use_item_name, const bool use_maccs_ii_field,
                           const bool state_from_archives) {
  return static_cast<int>(use_internal_regno) + (static_cast<int>(use_external_regno) * 2) +
         (static_cast<int>(use_item_name) * 4) + (static_cast<int>(use_maccs_ii_field) * 8) +
         (static_cast<int>(state_from_archives) * 16);
}

//-------------------------------------------------------------------------------------------------
uint getDataItemHeaderCode(const MolObjDataRequest &ask) {

  // It is obligatory to display the item name for requested data items, and "FROM ARCHIVES" is
  // not displayed in such data items.
  uint result = 4U;
  if (ask.placeInternalRegistryInHeader()) {
    result += 1U;
  }
  if (ask.getExternalRegistryNumber().size() > 0LLU) { 
    result += 2U;
  }
  if (ask.placeMaccsFieldInHeader()) {
    result += 8U;
  }
  return result;
}

} // namespace structure
} // namespace stormm
