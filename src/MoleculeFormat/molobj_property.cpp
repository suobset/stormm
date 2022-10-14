#include "copyright.h"
#include "Constants/behavior.h"
#include "FileManagement/file_listing.h"
#include "Parsing/ascii_numbers.h"
#include "Parsing/parse.h"
#include "molobj_property.h"

namespace stormm {
namespace structure {

using constants::CaseSensitivity;
using constants::ExceptionResponse;
using diskutil::getBaseName;
using parse::NumberFormat;
using parse::readIntegerValue;
using parse::readRealValue;
using parse::strncmpCased;
using parse::verifyContents;
  
//-------------------------------------------------------------------------------------------------
MolObjProperty::MolObjProperty(const char4 code_in, const int substrate_in,
                               const int entry_count_in, const int entry_depth_in,
                               const std::vector<MolObjPropField> &entry_detail_in,
                               const std::vector<MolObjIndexKind> &entry_adjustment_in,
                               const std::vector<int> &int_data_in,
                               const std::vector<double> &real_data_in,
                               const std::vector<std::string> &str_data_in,
                               const std::vector<std::string> &data_lines_in) :
    code{code_in}, substrate{substrate_in}, entry_count{entry_count_in},
    entry_depth{entry_depth_in}, entry_detail{entry_detail_in}, int_data{int_data_in},
    real_data{real_data_in}, str_data{str_data_in}, data_lines{data_lines_in}
{}

//-------------------------------------------------------------------------------------------------
MolObjProperty::MolObjProperty(const TextFile &tf, const int line_number, int *line_advance,
                               const std::string &title) :
    MolObjProperty()
{
  const char* line_ptr = tf.getLinePointer(line_number);
  const int lnlength = tf.getLineLength(line_number);
  if (lnlength < 6) {
    rtErr("Line " + std::to_string(line_number) + " of " + getBaseName(tf.getFileName()) +
          " cannot contain MDL MOL property information due to its length being only " +
          std::to_string(lnlength) + ".");
  }

  // Read the initial letter and three-letter code,
  code.w = line_ptr[0];
  code.x = line_ptr[3];
  code.y = line_ptr[4];
  code.z = line_ptr[5];
  kind = translateMolObjPropertyKind(code);

  // Interpret the code according to known options, starting with the maximum number of entries.
  int max_entries = 10000;
  switch (kind) {
  case MolObjPropertyKind::ATOM_ALIAS:
  case MolObjPropertyKind::ATOM_VALUE:
  case MolObjPropertyKind::GROUP_ABBREVIATION:
  case MolObjPropertyKind::ATOM_LIST:
  case MolObjPropertyKind::SGROUP_SUBSCRIPT:
  case MolObjPropertyKind::SGROUP_BOND_VECTOR:
  case MolObjPropertyKind::SGROUP_FIELD:
  case MolObjPropertyKind::SGROUP_DISPLAY:
  case MolObjPropertyKind::SGROUP_DATA:
  case MolObjPropertyKind::SPATIAL_FEATURE:
  case MolObjPropertyKind::PHANTOM_ATOM:
  case MolObjPropertyKind::SGROUP_CLASS:
  case MolObjPropertyKind::LARGE_REGNO:
    break;
  case MolObjPropertyKind::SGROUP_EXPANSION:
  case MolObjPropertyKind::SGROUP_ATOM_LIST:
  case MolObjPropertyKind::SGROUP_BOND_LIST:
  case MolObjPropertyKind::MG_PARENT_ATOM_LIST:
    max_entries = 15;
    break;
  case MolObjPropertyKind::CHARGE:
  case MolObjPropertyKind::RADICAL:
  case MolObjPropertyKind::ISOTOPE:
  case MolObjPropertyKind::RING_BOND_COUNT:
  case MolObjPropertyKind::SUBSTITUTION_COUNT:
  case MolObjPropertyKind::UNSATURATED_COUNT:
  case MolObjPropertyKind::RGROUP_LABEL_LOCATION:
  case MolObjPropertyKind::SGROUP_TYPE:
  case MolObjPropertyKind::SGROUP_SUBTYPE:
  case MolObjPropertyKind::SGROUP_LABELS:
  case MolObjPropertyKind::SGROUP_CONNECTIVITY:
  case MolObjPropertyKind::SGROUP_HIERARCHY:
  case MolObjPropertyKind::SGROUP_COMP_NUMBER:
  case MolObjPropertyKind::SGROUP_BRACKET_STYLE:
    max_entries = 8;
    break;
  case MolObjPropertyKind::SGROUP_CORRESPONENCE:
  case MolObjPropertyKind::SGROUP_ATTACH_POINT:
    max_entries = 6;
    break;
  case MolObjPropertyKind::LINK_ATOM:
  case MolObjPropertyKind::SGROUP_DISPLAY_INFO:
    max_entries = 4;
    break;
  case MolObjPropertyKind::ATTACHMENT_POINT:
  case MolObjPropertyKind::ATTACHMENT_ORDER:
    max_entries = 2;
    break;
  case MolObjPropertyKind::RGROUP_LOGIC:
    max_entries = 1;
    break;
  case MolObjPropertyKind::SKIP:
    break;
  }

  // Determine the layout and allocate the necessary space for each property.  Set the line
  // advancement, if appropriate.
  int tmp_advance = 0;
  bool entry_count_unrecognized = false;
  bool substrate_unrecognized = false;
  int entry_read_start_pos;
  std::vector<int> entry_data_bounds;
  switch (kind) {
  case MolObjPropertyKind::ATOM_ALIAS:
  case MolObjPropertyKind::ATOM_VALUE:
    substrate_unrecognized = readSubstrateIndex(line_ptr, 3);
    entry_count = 1;
    entry_depth = 0;
    tmp_advance = 1;
    entry_read_start_pos = 3;
    entry_data_bounds = { 0, 3 };
    break;
  case MolObjPropertyKind::GROUP_ABBREVIATION:
    substrate_unrecognized = readSubstrateIndex(line_ptr, 3);
    entry_count = 1;
    entry_depth = 1;
    entry_detail = { MolObjPropField::INTEGER };
    entry_adjustment = { MolObjIndexKind::ATOM };
    tmp_advance = 1;
    entry_read_start_pos = 3;
    entry_data_bounds = { 0, 3, 6 };
    break;
  case MolObjPropertyKind::CHARGE:
  case MolObjPropertyKind::RADICAL:
  case MolObjPropertyKind::ISOTOPE:
  case MolObjPropertyKind::RING_BOND_COUNT:
  case MolObjPropertyKind::SUBSTITUTION_COUNT:
  case MolObjPropertyKind::UNSATURATED_COUNT:
  case MolObjPropertyKind::ATTACHMENT_POINT:
  case MolObjPropertyKind::RGROUP_LABEL_LOCATION:
  case MolObjPropertyKind::SGROUP_TYPE:
  case MolObjPropertyKind::SGROUP_SUBTYPE:
  case MolObjPropertyKind::SGROUP_LABELS:
  case MolObjPropertyKind::SGROUP_CONNECTIVITY:
  case MolObjPropertyKind::SGROUP_HIERARCHY:
  case MolObjPropertyKind::SGROUP_COMP_NUMBER:
  case MolObjPropertyKind::SGROUP_BRACKET_STYLE:
    entry_count_unrecognized = readEntryCount(line_ptr);
    entry_depth = 2;
    entry_detail = { MolObjPropField::INTEGER, MolObjPropField::INTEGER };
    entry_adjustment = { MolObjIndexKind::ATOM, MolObjIndexKind::OTHER };
    entry_read_start_pos = 9;
    entry_data_bounds = { 0, 4, 8 };
    break;
  case MolObjPropertyKind::LINK_ATOM:
  case MolObjPropertyKind::RGROUP_LOGIC:
    entry_count_unrecognized = readEntryCount(line_ptr);
    entry_depth = 4;
    entry_detail = std::vector<MolObjPropField>(4, MolObjPropField::INTEGER);
    entry_adjustment = std::vector<MolObjIndexKind>(4, MolObjIndexKind::OTHER);
    entry_read_start_pos = 9;
    entry_data_bounds = { 0, 4, 8, 12, 16 };
    break;
  case MolObjPropertyKind::ATOM_LIST:
    break;
  case MolObjPropertyKind::ATTACHMENT_ORDER:
    entry_count_unrecognized = readEntryCount(line_ptr, 10);
    substrate_unrecognized = readSubstrateIndex(line_ptr);
    entry_depth = 4;
    entry_detail = std::vector<MolObjPropField>(4, MolObjPropField::INTEGER);
    entry_adjustment = { MolObjIndexKind::ATOM, MolObjIndexKind::OTHER, MolObjIndexKind::ATOM,
                         MolObjIndexKind::OTHER };
    entry_read_start_pos = 13;
    entry_data_bounds = { 0, 4, 8, 12, 16 };
    break;
  case MolObjPropertyKind::SGROUP_EXPANSION:
    if (strncmpCased("EXP", tf.extractString(line_number, 7, 3), CaseSensitivity::YES) == false) {
      rtErr("A malformed S-Group expansion entry was encountered on line " +
            std::to_string(line_number) + " of file " + getBaseName(tf.getFileName()) + ".",
            "MolObjProperty");
    }
    entry_count_unrecognized = readEntryCount(line_ptr, 10);
    entry_depth = 1;
    entry_detail = { MolObjPropField::INTEGER };
    entry_adjustment = { MolObjIndexKind::BOND };
    entry_read_start_pos = 13;
    entry_data_bounds = { 0, 4 };
    break;
  case MolObjPropertyKind::SGROUP_ATOM_LIST:
  case MolObjPropertyKind::SGROUP_BOND_LIST:
  case MolObjPropertyKind::MG_PARENT_ATOM_LIST:
    substrate_unrecognized = readSubstrateIndex(line_ptr);
    entry_count_unrecognized = readEntryCount(line_ptr, 10);
    entry_depth = 1;
    entry_detail = { MolObjPropField::INTEGER };
    entry_adjustment = { MolObjIndexKind::ATOM };
    entry_read_start_pos = 13;
    entry_data_bounds = { 0, 4 };
    break;
  case MolObjPropertyKind::SGROUP_SUBSCRIPT:
    entry_count = 1;
    substrate_unrecognized = readSubstrateIndex(line_ptr);
    entry_depth = 1;
    entry_detail = { MolObjPropField::STRING };
    entry_adjustment = { MolObjIndexKind::OTHER };
    entry_read_start_pos = 10;
    entry_data_bounds = { 0, tf.getLineLength(line_number) - 11 };
    break;
  case MolObjPropertyKind::SGROUP_CORRESPONENCE:
    substrate_unrecognized = readSubstrateIndex(line_ptr);
    entry_count_unrecognized = readEntryCount(line_ptr, 10);
    entry_depth = 3;
    entry_detail = std::vector<MolObjPropField>(3, MolObjPropField::INTEGER);
    entry_adjustment = std::vector<MolObjIndexKind>(3, MolObjIndexKind::BOND);
    entry_read_start_pos = 13;
    entry_data_bounds = { 0, 4, 8, 12 };
    break;
  case MolObjPropertyKind::SGROUP_DISPLAY_INFO:
    substrate_unrecognized = readSubstrateIndex(line_ptr);
    entry_count_unrecognized = readEntryCount(line_ptr, 10);
    entry_depth = 4;
    entry_detail = std::vector<MolObjPropField>(4, MolObjPropField::INTEGER);
    entry_adjustment = std::vector<MolObjIndexKind>(4, MolObjIndexKind::BOND);
    entry_read_start_pos = 13;
    entry_data_bounds = { 0, 3, 6, 9, 12 };
    break;
  case MolObjPropertyKind::SGROUP_BOND_VECTOR:
    substrate_unrecognized = readSubstrateIndex(line_ptr);
    entry_count = 1;
    entry_depth = 3;
    entry_detail = std::vector<MolObjPropField>(3, MolObjPropField::INTEGER);
    entry_adjustment = { MolObjIndexKind::BOND, MolObjIndexKind::OTHER, MolObjIndexKind::OTHER };
    entry_read_start_pos = 10;
    entry_data_bounds = { 0, 4, 7, 10 };
    break;
  case MolObjPropertyKind::SGROUP_FIELD:
    substrate_unrecognized = readSubstrateIndex(line_ptr);
    entry_count = 1;
    entry_depth = 1;
    entry_detail = { MolObjPropField::STRING, MolObjPropField::CHAR4, MolObjPropField::STRING,
                     MolObjPropField::CHAR4, MolObjPropField::STRING };
    entry_adjustment = std::vector<MolObjIndexKind>(5, MolObjIndexKind::OTHER);
    entry_read_start_pos = 10;
    entry_data_bounds = { 0, 30, 32, 52, 54, tf.getLineLength(line_number) - 65 };
    break;
  case MolObjPropertyKind::SGROUP_DISPLAY:
    substrate_unrecognized = readSubstrateIndex(line_ptr);
    entry_count = 1;
    entry_depth = 13;
    entry_detail = { MolObjPropField::REAL, MolObjPropField::REAL, MolObjPropField::INTEGER,
                     MolObjPropField::CHAR4, MolObjPropField::CHAR4, MolObjPropField::CHAR4,
                     MolObjPropField::INTEGER, MolObjPropField::INTEGER, MolObjPropField::INTEGER,
                     MolObjPropField::INTEGER, MolObjPropField::CHAR4, MolObjPropField::INTEGER,
                     MolObjPropField::INTEGER };
    entry_adjustment = std::vector<MolObjIndexKind>(13, MolObjIndexKind::OTHER);
    entry_read_start_pos = 10;
    entry_data_bounds = { 0, 10, 20, 24, 25, 26, 27, 29, 33, 36, 39, 41, 43, 45 };

    // Each "M  SDD" entry will be followed by zero or more "M  SCD" entries and a final "M  SED"
    // entry.  Scan ahead to determine the number of subsequent lines to skip.
    tmp_advance = line_number + 1;
    while (tmp_advance < tf.getLineCount() && tf.getLineLength(tmp_advance) >= 6 &&
           strncmpCased(tf.extractString(tmp_advance, 0, 6), "M  SCD")) {
      tmp_advance++;
    }
    if (strncmpCased(tf.extractString(tmp_advance, 0, 6), "M  SED")) {
      tmp_advance++;
    }
    else {
      rtErr("An MDL MOL property S-group display (\"M  SDD\") entry beginning on line " +
            std::to_string(line_number) + " of file " + getBaseName(tf.getFileName()) + " must "
            "be followed by S-group data line (optional \"M  SCD\" lines, and a terminating "
            "\"M  SED\" line).", "MolObjProperty");
    }
    tmp_advance -= line_number;
    break;
  case MolObjPropertyKind::SGROUP_DATA:

    // The S-group data will be read as part of the parent S-group display entry.
    break;
  case MolObjPropertyKind::SPATIAL_FEATURE:
    break;
  case MolObjPropertyKind::PHANTOM_ATOM:
    substrate_unrecognized = readSubstrateIndex(line_ptr);
    entry_count = 1;
    entry_depth = 5;
    entry_detail = { MolObjPropField::REAL, MolObjPropField::REAL, MolObjPropField::REAL,
                     MolObjPropField::CHAR4, MolObjPropField::STRING };
    entry_adjustment = std::vector<MolObjIndexKind>(5, MolObjIndexKind::OTHER);
    entry_read_start_pos = 9;
    entry_data_bounds = { 0, 10, 20, 31, 35, tf.getLineLength(line_number) - 45 };
    break;
  case MolObjPropertyKind::SGROUP_ATTACH_POINT:
    substrate_unrecognized = readSubstrateIndex(line_ptr);
    entry_count_unrecognized = readEntryCount(line_ptr, 10);
    entry_depth = 5;
    entry_detail = { MolObjPropField::INTEGER, MolObjPropField::INTEGER, MolObjPropField::CHAR4 };
    entry_adjustment = { MolObjIndexKind::ATOM, MolObjIndexKind::ATOM, MolObjIndexKind::OTHER };
    entry_read_start_pos = 13;
    entry_data_bounds = { 0, 4, 8, 10 };
    break;
  case MolObjPropertyKind::SGROUP_CLASS:
    substrate_unrecognized = readSubstrateIndex(line_ptr);
    entry_count = 1;
    entry_depth = 1;
    entry_detail = { MolObjPropField::STRING };
    entry_adjustment = { MolObjIndexKind::OTHER };
    entry_read_start_pos = 10;
    entry_data_bounds = { 0, tf.getLineLength(line_number) - 11 };
    break;
  case MolObjPropertyKind::LARGE_REGNO:
    entry_count = 1;
    entry_depth = 1;
    entry_detail = { MolObjPropField::INTEGER };
    entry_adjustment = { MolObjIndexKind::OTHER };
    entry_read_start_pos = 6;
    entry_data_bounds = { 0, tf.getLineLength(line_number) - 7 };
    break;
  case MolObjPropertyKind::SKIP:
    entry_count = 0;
    entry_depth = 0;
    if (verifyContents(line_ptr, 6, 3, NumberFormat::INTEGER)) {
      tmp_advance = readIntegerValue(line_ptr, 6, 3);
    }
    else {
      rtErr("The S  SKP property cannot skip " + tf.extractString(line_number, 6, 3) + " lines.",
            "MolObjProperty");
    }
    break;
  }
  *line_advance = line_number + tmp_advance;
  if (*line_advance >= tf.getLineCount()) {
    rtErr("The advancement due to property " + std::to_string(code.w) + "  " + code.x + code.y +
          code.z + " overruns the length of file " + getBaseName(tf.getFileName()) + ".",
          "MolObjProperty");
  }
  if (entry_count_unrecognized) {
    rtErr("The entry count on line " + std::to_string(line_number) + " of file " +
          getBaseName(tf.getFileName()) + " was not recognizable as an integer for property \"" +
          tf.extractString(line_number, 0, 6) + "\".", "MolObjProperty");
  }
  if (substrate_unrecognized) {
    rtErr("The substrate index on line " + std::to_string(line_number) + " of file " +
          getBaseName(tf.getFileName()) + " was not recognizable as an integer for property \"" +
          tf.extractString(line_number, 0, 6) + "\".", "MolObjProperty");
  }
  if (entry_count < 0 || entry_count > max_entries) {
    rtErr("Property \"" + tf.extractString(line_number, 0, 6) + "\" on line " +
          std::to_string(line_number) + " of file " + getBaseName(tf.getFileName()) + " has an "
          "invalid number of entries (" + std::to_string(entry_count) + ") (max " +
          std::to_string(max_entries) + ").", "MolObjProperty");
  }

  // Read the line for each type of property.  Advance the auxiliary counter line_advance if there
  // are additional data lines associated with the property, to avoid later trying to read such
  // information as new properties.
  switch (kind) {
  case MolObjPropertyKind::ATOM_ALIAS:
  case MolObjPropertyKind::ATOM_VALUE:
  case MolObjPropertyKind::GROUP_ABBREVIATION:
    data_lines.push_back(tf.extractString(tmp_advance, 0, tf.getLineLength(line_number)));
    break;
  case MolObjPropertyKind::SGROUP_DISPLAY:
    for (int i = line_number + 1; i < tmp_advance; i++) {
      data_lines.push_back(tf.extractString(i, 12, std::min(tf.getLineLength(i) - 12, 69)));
    }
    break;
  case MolObjPropertyKind::CHARGE:
  case MolObjPropertyKind::RADICAL:
  case MolObjPropertyKind::ISOTOPE:
    parseEntries(tf, line_number, entry_read_start_pos, entry_data_bounds);
    break;
  case MolObjPropertyKind::RING_BOND_COUNT:
  case MolObjPropertyKind::SUBSTITUTION_COUNT:
  case MolObjPropertyKind::UNSATURATED_COUNT:
  case MolObjPropertyKind::LINK_ATOM:
  case MolObjPropertyKind::ATOM_LIST:
  case MolObjPropertyKind::ATTACHMENT_POINT:
  case MolObjPropertyKind::ATTACHMENT_ORDER:
  case MolObjPropertyKind::RGROUP_LABEL_LOCATION:
  case MolObjPropertyKind::RGROUP_LOGIC:
  case MolObjPropertyKind::SGROUP_TYPE:
  case MolObjPropertyKind::SGROUP_SUBTYPE:
  case MolObjPropertyKind::SGROUP_LABELS:
  case MolObjPropertyKind::SGROUP_CONNECTIVITY:
  case MolObjPropertyKind::SGROUP_EXPANSION:
  case MolObjPropertyKind::SGROUP_ATOM_LIST:
  case MolObjPropertyKind::SGROUP_BOND_LIST:
  case MolObjPropertyKind::MG_PARENT_ATOM_LIST:
  case MolObjPropertyKind::SGROUP_SUBSCRIPT:
  case MolObjPropertyKind::SGROUP_CORRESPONENCE:
  case MolObjPropertyKind::SGROUP_DISPLAY_INFO:
  case MolObjPropertyKind::SGROUP_BOND_VECTOR:
  case MolObjPropertyKind::SGROUP_FIELD:
  case MolObjPropertyKind::SGROUP_DATA:
  case MolObjPropertyKind::SGROUP_HIERARCHY:
  case MolObjPropertyKind::SGROUP_COMP_NUMBER:
  case MolObjPropertyKind::SPATIAL_FEATURE:
  case MolObjPropertyKind::PHANTOM_ATOM:
  case MolObjPropertyKind::SGROUP_ATTACH_POINT:
  case MolObjPropertyKind::SGROUP_CLASS:
  case MolObjPropertyKind::LARGE_REGNO:
  case MolObjPropertyKind::SGROUP_BRACKET_STYLE:
  case MolObjPropertyKind::SKIP:
    break;
  }
}

//-------------------------------------------------------------------------------------------------
MolObjPropertyKind MolObjProperty::getKind() const {
  return kind;
}

//-------------------------------------------------------------------------------------------------
char4 MolObjProperty::getCode() const {
  return code;
}

//-------------------------------------------------------------------------------------------------
int MolObjProperty::getSubstrate() const {
  return substrate;
}

//-------------------------------------------------------------------------------------------------
int MolObjProperty::getEntryCount() const {
  return entry_count;
}

//-------------------------------------------------------------------------------------------------
int MolObjProperty::getIntegerValue(const int entry_index, const int attribute_index) const {
  checkAttributeValidity(entry_index, attribute_index, MolObjPropField::INTEGER);
  return int_data[(entry_depth * entry_index) + attribute_index];
}

//-------------------------------------------------------------------------------------------------
double MolObjProperty::getRealValue(const int entry_index, const int attribute_index) const {
  checkAttributeValidity(entry_index, attribute_index, MolObjPropField::REAL);
  return real_data[int_data[(entry_depth * entry_index) + attribute_index]];
}

//-------------------------------------------------------------------------------------------------
char4 MolObjProperty::getChar4Value(const int entry_index, const int attribute_index) const {
  checkAttributeValidity(entry_index, attribute_index, MolObjPropField::CHAR4);
  const uint ui_result = int_data[(entry_depth * entry_index) + attribute_index];
  return { static_cast<char>(ui_result & 0xff),
           static_cast<char>((ui_result >>  8) & 0xff),
           static_cast<char>((ui_result >> 16) & 0xff),
           static_cast<char>((ui_result >> 24) & 0xff) };
}

//-------------------------------------------------------------------------------------------------
std::string MolObjProperty::getStringValue(const int entry_index,
                                           const int attribute_index) const {
  checkAttributeValidity(entry_index, attribute_index, MolObjPropField::STRING);
  return str_data[int_data[(entry_depth * entry_index) + attribute_index]];
}

//-------------------------------------------------------------------------------------------------
const std::string& MolObjProperty::getDataLine(const int index) const {
  if (index < 0 || index > static_cast<int>(data_lines.size())) {
    rtErr("A property (code \"" + std::to_string(code.w) + "  " + code.x + code.y + code.z +
          "\") with " + std::to_string(data_lines.size()) + " cannot produce a data line with "
          "index " + std::to_string(index) + ".", "MolObjProperty", "getDataLine");
  }
  return data_lines[index];
}

//-------------------------------------------------------------------------------------------------
void MolObjProperty::setCode(const char x, const char y, const char z, const char major) {
  code.x = x;
  code.y = y;
  code.z = z;
  code.w = major;
}

//-------------------------------------------------------------------------------------------------
void MolObjProperty::setCode(const char4 code_in) {
  code = code_in;
}

//-------------------------------------------------------------------------------------------------
void MolObjProperty::setSubstrate(const int index) {
  substrate = index;
}

//-------------------------------------------------------------------------------------------------
void MolObjProperty::setEntryFormat(const std::vector<MolObjPropField> &entry_detail_in,
                                    const std::vector <MolObjIndexKind> &entry_adjustment_in) {
  if (entry_depth != 0) {
    rtErr("A property with defined fields cannot be redesigned.", "MolObjProperty",
          "setEntryFormat");
  }
  if (entry_detail_in.size() != entry_adjustment_in.size()) {
    rtErr("Details and adjustment instructions for each property field must have a one-to-one "
          "correspondence.", "MolObjProperty", "setEntryFormat");
  }
  entry_detail = entry_detail_in;
  entry_adjustment = entry_adjustment_in;
  entry_depth = entry_detail.size();
  for (int i = 0; i < entry_depth; i++) {
    if (entry_detail[i] != MolObjPropField::INTEGER &&
        (entry_adjustment[i] == MolObjIndexKind::ATOM ||
         entry_adjustment[i] == MolObjIndexKind::BOND)) {
      rtErr("Atom and bond indices must have integer type.  Property field " +
            std::to_string(i) + " violates convention by combining a " +
            getEnumerationName(entry_detail[i]) + " data type with " +
            getEnumerationName(entry_adjustment[i]) + " index adjustment.", "MolObjProperty",
            "setEntryFormat");
    }
  }
}

//-------------------------------------------------------------------------------------------------
bool MolObjProperty::readEntryCount(const char* line_ptr, const int start_pos, const int length) {
  if (verifyContents(line_ptr, start_pos, length, NumberFormat::INTEGER)) {
    entry_count = readIntegerValue(line_ptr, start_pos, length);
    return false;
  }
  else {
    return true;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
bool MolObjProperty::readSubstrateIndex(const char* line_ptr, const int start_pos,
                                        const int length) {
  if (verifyContents(line_ptr, start_pos, length, NumberFormat::INTEGER)) {
    substrate = readIntegerValue(line_ptr, start_pos, length) - 1;
    return false;
  }
  else {
    return true;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
void MolObjProperty::parseEntries(const TextFile &tf, const int line_number, const int start_pos,
                                  const std::vector<int> &limits) {

  // There are integers for every piece of information: char4 data gets converted to a bit-packed
  // integer, pieces of string and double-precision data log integers for their index in the
  // corresponding string or double-precision lists.
  int_data.resize(entry_depth * entry_count);
  int n_real   = 0;
  int n_string = 0;
  for (int i = 0; i < entry_depth; i++) {
    n_real   += (entry_detail[i] == MolObjPropField::REAL);
    n_string += (entry_detail[i] == MolObjPropField::STRING);
  }
  real_data.reserve(n_real * entry_count);
  str_data.reserve(n_string * entry_count);
  const char* line_ptr = tf.getLinePointer(line_number);
  const int lnlength = tf.getLineLength(line_number);
  n_real = 0;
  n_string = 0;
  for (int i = 0; i < entry_count; i++) {
    for (int j = 0; j < entry_depth; j++) {
      const int llim = start_pos + (i * limits[entry_depth]) + limits[j];
      int slen = limits[j + 1] - limits[j];
      slen = (entry_detail[i] == MolObjPropField::CHAR4) ? std::min(slen, 4) : slen;
      if (llim + slen > lnlength) {
        rtErr("Reading entry " + std::to_string(i) + " field " + std::to_string(j) +
              " of property \"" + tf.extractString(line_number, 0, 6) + "\" at line " +
              std::to_string(line_number) + " of file " + getBaseName(tf.getFileName()) +
              "would overrun the line length (" + std::to_string(lnlength) + ").",
              "MolObjProperty", "parseEntries");
      }
      switch (entry_detail[i]) {
      case MolObjPropField::INTEGER:
        if (verifyContents(line_ptr, llim, slen, NumberFormat::INTEGER)) {
          int_data[(i * entry_depth) + j] = readIntegerValue(line_ptr, llim, slen);
          switch (entry_adjustment[j]) {
          case MolObjIndexKind::ATOM:
          case MolObjIndexKind::BOND:
            int_data[(i * entry_depth) + j] -= 1;
            break;
          case MolObjIndexKind::OTHER:
            break;
          }
        }
        else {
          rtErr("Failed to parse an integer from characters " + std::to_string(llim) + " - " +
                std::to_string(llim + slen - 1) + " of line " + std::to_string(line_number) +
                "(an \"" + tf.extractString(line_number, 0, 6) + "\" property) of file " +
                getBaseName(tf.getFileName()) + ".", "MolObjProperty", "parseEntries");
        }
        break;
      case MolObjPropField::CHAR4:
        if (verifyContents(line_ptr, llim, slen, NumberFormat::CHAR4)) {
          char4 result;
          result.x = (slen > 0) ? line_ptr[llim    ] : ' ';
          result.y = (slen > 1) ? line_ptr[llim + 1] : ' ';
          result.z = (slen > 2) ? line_ptr[llim + 2] : ' ';
          result.w = (slen > 3) ? line_ptr[llim + 3] : ' ';
          uint uiresult = ((static_cast<uint>(result.w) << 24) |
                           (static_cast<uint>(result.z) << 16) |
                           (static_cast<uint>(result.y) <<  8) |
                           (static_cast<uint>(result.x)));
          int_data[(i * entry_depth) + j] = uiresult;
        }
        break;
      case MolObjPropField::REAL:
        if (verifyContents(line_ptr, llim, slen, NumberFormat::STANDARD_REAL) ||
            verifyContents(line_ptr, llim, slen, NumberFormat::SCIENTIFIC)) {
          int_data[(i * entry_depth) + j] = n_real;
          real_data.push_back(readRealValue(line_ptr, llim, slen));
        }
        break;
      case MolObjPropField::STRING:
        int_data[(i * entry_depth) + j] = n_string;
        str_data.push_back(tf.extractString(line_number, llim, slen));
        break;
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------
void MolObjProperty::checkAttributeValidity(const int entry_index, const int attribute_index,
                                            const MolObjPropField expectation) const {
  if (entry_index < 0 || entry_index >= entry_count) {
    rtErr("Entry " + std::to_string(entry_index) + " does not exist for an \"" + code.w + "  " +
          code.x + code.y + code.z + "\" property with " + std::to_string(entry_count) +
          " entries.", "MolObjProperty", "checkAttributeValidity");
  }
  if (attribute_index < 0 || attribute_index >= entry_depth) {
    rtErr("Attribute index " + std::to_string(attribute_index) + " is invalid for an \"" +
          code.w + "  " + code.x + code.y + code.z + "\" property.", "MolObjProperty",
          "checkAttributeValidity");
  }
  if (entry_detail[attribute_index] != expectation) {
    rtErr("Attribute " + std::to_string(attribute_index) + " is of type " +
          getEnumerationName(entry_detail[attribute_index]).c_str() + ", not " +
          getEnumerationName(expectation).c_str() + ".", "MolObjProperty",
          "checkAttributeValidity");
  }
}

} // namespace structure
} // namespace stormm
