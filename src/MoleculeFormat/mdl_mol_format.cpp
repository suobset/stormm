#include "copyright.h"
#include "FileManagement/file_listing.h"
#include "mdl_mol_format.h"

namespace stormm {
namespace structure {

using diskutil::getBaseName;
using parse::NumberFormat;
using parse::readIntegerValue;
using parse::readRealValue;
using parse::strncmpCased;
using parse::TextFileReader;
using parse::verifyContents;
using trajectory::CoordinateFrameWriter;
using trajectory::PhaseSpaceWriter;
  
//-------------------------------------------------------------------------------------------------
MolObjBond::MolObjBond() :
    i_atom{-1}, j_atom{-1}, order{MolObjBondOrder::SINGLE}, stereo{MolObjBondStereo::NOT_STEREO},
    ring_state{MolObjRingState::EITHER}, reactivity{MolObjReactionCenter::NON_CENTER}
{}

//-------------------------------------------------------------------------------------------------
MolObjBond::MolObjBond(const int i_atom_in, const int j_atom_in) :
    i_atom{i_atom_in}, j_atom{j_atom_in}, order{MolObjBondOrder::SINGLE},
    stereo{MolObjBondStereo::NOT_STEREO}, ring_state{MolObjRingState::EITHER},
    reactivity{MolObjReactionCenter::NON_CENTER}
{}

//-------------------------------------------------------------------------------------------------
MolObjBond::MolObjBond(const int i_atom_in, const int j_atom_in, const MolObjBondOrder order_in,
                       const MolObjBondStereo stereo_in, const MolObjRingState ring_state_in,
                       const MolObjReactionCenter reactivity_in) :
    i_atom{i_atom_in}, j_atom{j_atom_in}, order{order_in}, stereo{stereo_in},
    ring_state{ring_state_in}, reactivity{reactivity_in}
{}

//-------------------------------------------------------------------------------------------------
MolObjBond::MolObjBond(const TextFile &tf, const int line_number, const std::string &title) :
    MolObjBond()
{
  const char* bond_line_ptr = tf.getLinePointer(line_number);
  if (tf.getLineLength(line_number) < 6) {
    rtErr("Line " + std::to_string(line_number) + " of " + getBaseName(tf.getFileName()) +
          " cannot contain MDL MOL bond information due to its length being only " +
          std::to_string(tf.getLineLength(line_number)));
  }
  i_atom = readIntegerValue(bond_line_ptr, 0, 3);
  j_atom = readIntegerValue(bond_line_ptr, 3, 3);
  order = (verifyContents(tf, line_number, 6, 3)) ?
          interpretBondOrder(readIntegerValue(bond_line_ptr, 6, 3), title) :
          default_mdl_bond_order;
  stereo = (verifyContents(tf, line_number, 9, 3)) ?
           interpretBondStereochemistry(readIntegerValue(bond_line_ptr, 9, 3), title) :
           default_mdl_bond_stereochemistry;
  ring_state = (verifyContents(tf, line_number, 12, 3)) ?
               interpretRingState(readIntegerValue(bond_line_ptr, 12, 3), title) :
               default_mdl_ring_status;
  reactivity = (verifyContents(tf, line_number, 15, 3)) ?
               interpretBondReactivePotential(readIntegerValue(bond_line_ptr, 15, 3), title) :
               default_mdl_bond_reactivity;
}

//-------------------------------------------------------------------------------------------------
MolObjBondOrder MolObjBond::interpretBondOrder(const int code_in, const std::string &title) {
  switch (code_in) {
  case 1:
    return MolObjBondOrder::SINGLE;
  case 2:
    return MolObjBondOrder::DOUBLE;
  case 3:
    return MolObjBondOrder::TRIPLE;
  case 4:
    return MolObjBondOrder::AROMATIC;
  case 5:
    return MolObjBondOrder::SINGLE_OR_DOUBLE;
  case 6:
    return MolObjBondOrder::SINGLE_OR_AROMATIC;
  case 7:
    return MolObjBondOrder::DOUBLE_OR_AROMATIC;
  case 8:
    return MolObjBondOrder::ANY;
  default:
    rtErr("An invalid bond order code " + std::to_string(code_in) + " was specified in an MDL MOL "
          "entry titled \"" + title + "\".", "MdlMolObj", "interpretBondOrder");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
MolObjBondStereo MolObjBond::interpretBondStereochemistry(const int code_in,
                                                          const std::string &title) {
  switch (code_in) {
  case 0:
    return MolObjBondStereo::NOT_STEREO;
  case 1:
    return MolObjBondStereo::UP;
  case 3:
    return MolObjBondStereo::CIS_OR_TRANS;
  case 4:
    return MolObjBondStereo::EITHER;
  case 6:
    return MolObjBondStereo::DOWN;
  default:
    rtErr("An invalid bond stereochemistry code " + std::to_string(code_in) + " was specified in "
          "an MDL MOL entry titled \"" + title + "\".", "MdlMolObj",
          "interpretBondStereoChemistry");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
MolObjRingState MolObjBond::interpretRingState(const int code_in, const std::string &title) {
  switch (code_in) {
  case 0:
    return MolObjRingState::EITHER;
  case 1:
    return MolObjRingState::RING;
  case 2:
    return MolObjRingState::CHAIN;
  default:
    rtErr("An invalid bond ring status code " + std::to_string(code_in) + " was specified in an "
          "MDL MOL entry titled \"" + title + "\".", "MdlMolObj", "interpretRingState");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
MolObjReactionCenter MolObjBond::interpretBondReactivePotential(const int code_in,
                                                                const std::string &title) {
  switch (code_in) {
  case -1:
    return MolObjReactionCenter::NON_CENTER;
  case 0:
    return MolObjReactionCenter::UNMARKED;
  case 1:
    return MolObjReactionCenter::CENTER;
  case 2:
    return MolObjReactionCenter::UNREACTIVE;
  case 4:
    return MolObjReactionCenter::BOND_MADE_OR_BROKEN;
  case 5:
    return MolObjReactionCenter::CENTER_WITH_FORMATION;
  case 8:
    return MolObjReactionCenter::BOND_ORDER_CHANGE;
  case 9:
    return MolObjReactionCenter::CENTER_WITH_ORDER_CHANGE;
  case 12:
    return MolObjReactionCenter::BOND_FORMATION_AND_ORDER_CHANGE;
  case 13:
    return MolObjReactionCenter::CENTER_WITH_FORMATION_AND_ORDER_CHANGE;
  default:
    rtErr("An invalid bond reactive potential code " + std::to_string(code_in) + " was specified "
          "in an MDL MOL entry titled \"" + title + "\".", "MdlMolObj",
          "interpretBondReactivePotential");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
int MolObjBond::getFirstAtom() const {
  return i_atom;
}

//-------------------------------------------------------------------------------------------------
int MolObjBond::getSecondAtom() const {
  return j_atom;
}

//-------------------------------------------------------------------------------------------------
MolObjBondOrder MolObjBond::getOrder() const {
  return order;
}

//-------------------------------------------------------------------------------------------------
MolObjBondStereo MolObjBond::getStereochemistry() const {
  return stereo;
}

//-------------------------------------------------------------------------------------------------
MolObjRingState MolObjBond::getRingStatus() const {
  return ring_state;
}

//-------------------------------------------------------------------------------------------------
MolObjReactionCenter MolObjBond::getReactivity() const {
  return reactivity;
}

//-------------------------------------------------------------------------------------------------
void MolObjBond::setFirstAtom(const int index_in) {
  i_atom = index_in;
}

//-------------------------------------------------------------------------------------------------
void MolObjBond::setSecondAtom(const int index_in) {
  j_atom = index_in;
}

//-------------------------------------------------------------------------------------------------
void MolObjBond::setOrder(const MolObjBondOrder order_in) {
  order = order_in;
}

//-------------------------------------------------------------------------------------------------
void MolObjBond::setStereochemistry(const MolObjBondStereo stereo_in) {
  stereo = stereo_in;
}

//-------------------------------------------------------------------------------------------------
void MolObjBond::setRingStatus(const MolObjRingState status_in) {
  ring_state = status_in;
}

//-------------------------------------------------------------------------------------------------
void MolObjBond::setReactivity(const MolObjReactionCenter potential_in) {
  reactivity = potential_in;
}

//-------------------------------------------------------------------------------------------------
MolObjAtomList::MolObjAtomList(const std::vector<int> &atomic_numbers_in, bool const exclusions_in,
                               const int atom_attachment_in) :
  entry_count{static_cast<int>(atomic_numbers_in.size())},
  atomic_numbers{atomic_numbers_in},
  exclusions{exclusions_in},
  atom_attachment{atom_attachment_in}
{}

//-------------------------------------------------------------------------------------------------
MolObjAtomList::MolObjAtomList(const TextFile &tf, const int line_number,
                               const std::string &title) :
    MolObjAtomList()
{
  // Determine whether this originates in an atom list entry or a property of the MDL section
  const char* line_ptr = tf.getLinePointer(line_number);
  const int lnlength = tf.getLineLength(line_number);
  if (lnlength < 6) {
    rtErr("Line " + std::to_string(line_number) + " of " + getBaseName(tf.getFileName()) +
          " cannot contain MDL MOL atom list information due to its length being only " +
          std::to_string(lnlength) + ".");
  }
  if (strncmpCased(line_ptr, std::string("M  ALS"))) {

    // The entry originates as a property.  This would invalidate any preceding atom list entries
    // of the deprecated format, but such a contingency has already been taken care of in the
    // parent MdlMolObj reader.
    if (lnlength < 15) {
      rtErr("There is insufficient information on line " + std::to_string(line_number) + " of " +
            tf.getFileName() + " to read a property-based MDL MOL atom list.", "MolObjAtomList");
    }
    if (verifyContents(line_ptr, 7, 3, NumberFormat::INTEGER)) {
      atom_attachment = readIntegerValue(line_ptr, 7, 3);
    }
    if (verifyContents(line_ptr, 10, 3, NumberFormat::INTEGER)) {
      entry_count = readIntegerValue(line_ptr, 10, 3);
    }
    if (entry_count < 0 || entry_count > 16) {
      rtErr("An invalid number of entries, " + std::to_string(entry_count) + ", was recorded at "
            "line " + std::to_string(line_number) + " of " + tf.getFileName() + ".",
            "MolObjAtomList");
    }
    switch (line_ptr[14]) {
    case 'T':
      exclusions = true;
      break;
    case 'F':
      exclusions = false;
      break;
    default:
      rtErr("The exclusions flag at line " + std::to_string(line_number) + " of " +
            tf.getFileName() + ", " + line_ptr[14] + ", is invalid.", "MolObjAtomList");
    }
    if (lnlength < 16 + (4 * entry_count)) {
      rtErr("There is insufficient information on line " + std::to_string(line_number) + " of " +
            tf.getFileName() + " to read an atom list with " + std::to_string(entry_count) +
            "items.", "MolObjAtomList");      
    }
    std::vector<char4> tmp_symbols;
    for (int i = 0; i < entry_count; i++) {
      tmp_symbols[i].x = line_ptr[16 + (4 * i)];
      tmp_symbols[i].y = line_ptr[17 + (4 * i)];
      tmp_symbols[i].z = line_ptr[18 + (4 * i)];
      tmp_symbols[i].w = line_ptr[19 + (4 * i)];
    }
    atomic_numbers = symbolToZNumber(tmp_symbols, CaseSensitivity::YES, ExceptionResponse::DIE);
  }
  else {
    if (lnlength < 10) {
      rtErr("There is insufficient information on line " + std::to_string(line_number) + " of " +
            tf.getFileName() + " to read one of the (deprecated) MDL MOL atom list entries.",
            "MolObjAtomList");
    }

    // Read the deprecated atom list entry format
    if (verifyContents(line_ptr, 0, 3, NumberFormat::INTEGER)) {
      atom_attachment = readIntegerValue(line_ptr, 0, 3);
    }
    switch (line_ptr[4]) {
    case 'T':
      exclusions = true;
      break;
    case 'F':
      exclusions = false;
      break;
    default:
      rtErr("The exclusions flag at line " + std::to_string(line_number) + " of " +
            tf.getFileName() + ", " + line_ptr[4] + ", is invalid.", "MolObjAtomList");
    }
    if (verifyContents(line_ptr, 5, 5, NumberFormat::INTEGER)) {
      entry_count = readIntegerValue(line_ptr, 5, 5);
    }
    if (entry_count < 0 || entry_count > 5) {
      rtErr("An invalid number of entries, " + std::to_string(entry_count) + ", was recorded at "
            "line " + std::to_string(line_number) + " of " + tf.getFileName() + ".",
            "MolObjAtomList");
    }
    if (lnlength < 10 + (4 * entry_count)) {
      rtErr("There is insufficient information on line " + std::to_string(line_number) + " of " +
            tf.getFileName() + " to read an atom list with " + std::to_string(entry_count) +
            "items.", "MolObjAtomList");
    }
    atomic_numbers.resize(entry_count);
    for (int i = 0; i < entry_count; i++) {
      if (verifyContents(line_ptr, 11 + (4 * i), 3, NumberFormat::INTEGER)) {
        atomic_numbers[i] = readIntegerValue(line_ptr, 11 + (4 * i), 3);
      }
      else {
        rtErr("An invalid entry was found on line " + std::to_string(line_number) + " of " +
              tf.getFileName() + ": " + tf.extractString(line_number, 11 + (4 * i), 3) + ".  The "
              "entries of one of the (deprecated) MDL MOL format atom lists must be integers "
              "corresponding to atomic numbers.  Use IUPAC element symbols in the new format, "
              "placed on MOL property lines beginning \"M  ALS\".", "MolObjAtomList");
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------
MolObjProperty::MolObjProperty(const char4 code_in, const int substrate_in,
                               const int entry_count_in, const int entry_depth_in,
                               const std::vector<MolObjPropField> &entry_detail_in,
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
  code.x = line_ptr[4];
  code.x = line_ptr[5];
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
    entry_count = 1;
    entry_depth = 1;
    entry_detail = { MolObjPropField::INTEGER };
    tmp_advance = 1;
    entry_read_start_pos = 3;
    entry_data_bounds = { 0, 3 };
    break;
  case MolObjPropertyKind::GROUP_ABBREVIATION:
    entry_count = 1;
    entry_depth = 2;
    entry_detail = { MolObjPropField::INTEGER, MolObjPropField::INTEGER };
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
    entry_read_start_pos = 10;
    entry_data_bounds = { 0, 4, 8 };
    break;
  case MolObjPropertyKind::LINK_ATOM:
  case MolObjPropertyKind::RGROUP_LOGIC:
    entry_count_unrecognized = readEntryCount(line_ptr);
    entry_depth = 4;
    entry_detail = { MolObjPropField::INTEGER, MolObjPropField::INTEGER, MolObjPropField::INTEGER,
                     MolObjPropField::INTEGER };
    entry_read_start_pos = 10;
    entry_data_bounds = { 0, 4, 8, 12, 16 };
    break;
  case MolObjPropertyKind::ATOM_LIST:
    break;
  case MolObjPropertyKind::ATTACHMENT_ORDER:
    entry_count_unrecognized = readEntryCount(line_ptr, 10);
    substrate_unrecognized = readSubstrateIndex(line_ptr);
    entry_depth = 4;
    entry_detail = { MolObjPropField::INTEGER, MolObjPropField::INTEGER, MolObjPropField::INTEGER,
                     MolObjPropField::INTEGER };
    entry_read_start_pos = 14;
    entry_data_bounds = { 0, 4, 8, 12, 16 };
    break;
  case MolObjPropertyKind::SGROUP_EXPANSION:
    if (strncmpCased("EXP", tf.extractString(line_number, 7, 3), CaseSensitivity::YES) == false) {
      rtErr("A malformed S-Group expansion entry was encountered on line " +
            std::to_string(line_number) + " of file " + tf.getFileName() + ".", "MolObjProperty");
    }
    entry_count_unrecognized = readEntryCount(line_ptr, 10);
    entry_depth = 1;
    entry_detail = { MolObjPropField::INTEGER };
    entry_read_start_pos = 14;
    entry_data_bounds = { 0, 4 };
    break;
  case MolObjPropertyKind::SGROUP_ATOM_LIST:
  case MolObjPropertyKind::SGROUP_BOND_LIST:
  case MolObjPropertyKind::MG_PARENT_ATOM_LIST:
    substrate_unrecognized = readSubstrateIndex(line_ptr);
    entry_count_unrecognized = readEntryCount(line_ptr, 10);
    entry_depth = 1;
    entry_detail = { MolObjPropField::INTEGER };
    entry_read_start_pos = 14;
    entry_data_bounds = { 0, 4 };
    break;
  case MolObjPropertyKind::SGROUP_SUBSCRIPT:
    entry_count = 1;
    substrate_unrecognized = readSubstrateIndex(line_ptr);
    entry_depth = 1;
    entry_detail = { MolObjPropField::STRING };
    entry_read_start_pos = 11;
    entry_data_bounds = { 0, tf.getLineLength(line_number) - 11 };
    break;
  case MolObjPropertyKind::SGROUP_CORRESPONENCE:
    substrate_unrecognized = readSubstrateIndex(line_ptr);
    entry_count_unrecognized = readEntryCount(line_ptr, 10);
    entry_depth = 3;
    entry_detail = { MolObjPropField::INTEGER, MolObjPropField::INTEGER,
                     MolObjPropField::INTEGER };
    entry_read_start_pos = 14;
    entry_data_bounds = { 0, 4, 8, 12 };
    break;
  case MolObjPropertyKind::SGROUP_DISPLAY_INFO:
    substrate_unrecognized = readSubstrateIndex(line_ptr);
    entry_count_unrecognized = readEntryCount(line_ptr, 10);
    entry_depth = 4;
    entry_detail = { MolObjPropField::INTEGER, MolObjPropField::INTEGER, MolObjPropField::INTEGER,
                     MolObjPropField::INTEGER };    
    entry_read_start_pos = 14;
    entry_data_bounds = { 0, 3, 6, 9, 12 };
    break;
  case MolObjPropertyKind::SGROUP_BOND_VECTOR:
    substrate_unrecognized = readSubstrateIndex(line_ptr);
    entry_count = 1;
    entry_depth = 3;
    entry_detail = { MolObjPropField::INTEGER, MolObjPropField::INTEGER,
                     MolObjPropField::INTEGER };    
    entry_read_start_pos = 11;
    entry_data_bounds = { 0, 4, 7, 10 };
    break;
  case MolObjPropertyKind::SGROUP_FIELD:
    substrate_unrecognized = readSubstrateIndex(line_ptr);
    entry_count = 1;
    entry_depth = 1;
    entry_detail = { MolObjPropField::STRING, MolObjPropField::CHAR4, MolObjPropField::STRING,
                     MolObjPropField::CHAR4, MolObjPropField::STRING };
    entry_read_start_pos = 11;
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
    entry_read_start_pos = 11;
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
    entry_read_start_pos = 10;
    entry_data_bounds = { 0, 10, 20, 31, 35, tf.getLineLength(line_number) - 45 }; 
    break;
  case MolObjPropertyKind::SGROUP_ATTACH_POINT:
    substrate_unrecognized = readSubstrateIndex(line_ptr);
    entry_count_unrecognized = readEntryCount(line_ptr, 10);
    entry_depth = 5;
    entry_detail = { MolObjPropField::INTEGER, MolObjPropField::INTEGER, MolObjPropField::CHAR4 };
    entry_read_start_pos = 14;
    entry_data_bounds = { 0, 4, 8, 10 };
    break;
  case MolObjPropertyKind::SGROUP_CLASS:
    substrate_unrecognized = readSubstrateIndex(line_ptr);
    entry_count = 1;
    entry_depth = 1;
    entry_detail = { MolObjPropField::STRING };
    entry_read_start_pos = 11;
    entry_data_bounds = { 0, tf.getLineLength(line_number) - 11 };
    break;
  case MolObjPropertyKind::LARGE_REGNO:
    entry_count = 1;
    entry_depth = 1;
    entry_detail = { MolObjPropField::INTEGER };
    entry_read_start_pos = 7;
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
          code.z + " overruns the length of file " + tf.getFileName() + ".", "MolObjProperty");
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
    parseEntry(line_ptr, entry_read_start_pos, entry_data_bounds);
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
    entry_count = readIntegerValue(line_ptr, start_pos, length);
    return false;
  }
  else {
    return true;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
void MolObjProperty::parseEntry(const char* line_ptr, const int start_pos,
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
  for (int i = 0; i < entry_count; i++) {
    for (int j = 0; j < entry_depth; j++) {
      switch (entry_detail[i]) {
      case MolObjPropField::INTEGER:
        
        break;
      case MolObjPropField::CHAR4:
        break;
      case MolObjPropField::REAL:
        break;
      case MolObjPropField::STRING:
        break;
      }
    }
  }
}
                           

//-------------------------------------------------------------------------------------------------
MdlMolObj::MdlMolObj():
    version_no{MdlMolVersion::V2000}, atom_count{0}, bond_count{0}, list_count{0}, sgroup_count{0},
    constraint_count{0}, chirality{MolObjChirality::ACHIRAL}, registry_number{-1}, coordinates{},
    atomic_symbols{}, atomic_numbers{}, formal_charges{}, isotopic_shifts{}, parities{},
    implicit_hydrogens{}, stereo_considerations{}, valence_connections{},
    atom_atom_mapping_count{}, orientation_stability{}, bonds{}, element_lists{}, stext_entries{},
    properties{}, title{""}, software_details{""}, general_comment{""}
{}

//-------------------------------------------------------------------------------------------------
MdlMolObj::MdlMolObj(const TextFile &tf, const int line_start, const int line_end_in,
                     const CaseSensitivity capitalization, const ExceptionResponse policy):
  MdlMolObj()
{
  const TextFileReader tfr = tf.data();

  // Default line end of -1 indicates reading to the end of the file.  Otherwise, identify the
  // end of the formatting ("M  END").
  const int mdl_section_end = getMdlMolSectionEnd(tfr, line_start, line_end_in);
                       
  // The range of data now extends from line_start to mdl_section_end.  Sift through that
  // information for a V2000 or V3000 specification.  This should be found on the fourth line.
  version_no = findMolObjVersion(tf, line_start + 3);
  
  // Begin by reading the molecule name (title), generating software details, and any general
  // comment (always three and only three distinct lines, even if left blank).
  if (line_start + 2 < tfr.line_count) {
    title            = tf.extractString(line_start);
    software_details = tf.extractString(line_start + 1);
    general_comment  = tf.extractString(line_start + 2);
  }

  // Read the counts line
  switch (version_no) {
  case MdlMolVersion::V2000:
    {
      const int counts_line_idx = line_start + 3;
      const char* counts_line_ptr = &tfr.text[tfr.line_limits[line_start + 3]];
      atom_count   = readIntegerValue(counts_line_ptr, 0, 3);
      bond_count   = readIntegerValue(counts_line_ptr, 3, 3);
      if (verifyContents(tf, counts_line_idx, 6, 3, NumberFormat::INTEGER)) {
        list_count   = readIntegerValue(counts_line_ptr, 6, 3);
      }
      if (verifyContents(tf, counts_line_idx, 12, 3, NumberFormat::INTEGER)) {
        const int chir_num = readIntegerValue(counts_line_ptr, 12, 3);
        if (chir_num == 0) {
          chirality = MolObjChirality::ACHIRAL;
        }
        else if (chir_num == 1) {
          chirality = MolObjChirality::CHIRAL;
        }
        else {
          rtErr("Invalid chirality setting detected at line " + std::to_string(counts_line_idx) +
                " in .sdf or MDL MOL file " + tf.getFileName() + ".", "MdlMolObj");
        }
      }
      if (verifyContents(tf, counts_line_idx, 15, 3, NumberFormat::INTEGER)) {
        stext_entry_count = readIntegerValue(counts_line_ptr, 15, 3);
      }
      if (verifyContents(tf, counts_line_idx, 30, 3, NumberFormat::INTEGER)) {
        properties_count  = readIntegerValue(counts_line_ptr, 30, 3);
      }
    
      // Validation
      if (atom_count > 255) {
        rtErr("A V2000 MOL format entry cannot contain more than 255 atoms.", "MdlMolObj");
      }
    }
    break;
  case MdlMolVersion::V3000:
    break;
  case MdlMolVersion::UNKNOWN:

    // This case would be unprocessable.
    rtErr("No valid MDL MOL version was detected.  Parsing in " + tf.getFileName() + " cannot "
          "proceed.");
  }

  // Allocate space for information to be read
  allocate();
    
  // Read the atoms block
  int iatm = 0;
  switch (version_no) {
  case MdlMolVersion::V2000:
    for (int i = line_start + 4; i < line_start + 4 + atom_count; i++) {
      const char* atom_line_ptr = &tfr.text[tfr.line_limits[i]];
      coordinates[iatm] = { readRealValue(atom_line_ptr,  0, 10),
                            readRealValue(atom_line_ptr, 10, 10),
                            readRealValue(atom_line_ptr, 20, 10) };
      if (verifyContents(tf, i, 31, 3, NumberFormat::CHAR4)) {
        atomic_symbols[iatm] = tf.extractChar4(i, 31, 3);
      }
      if (verifyContents(tf, i, 34, 2, NumberFormat::INTEGER)) {
        isotopic_shifts[iatm] = readIntegerValue(atom_line_ptr, 34, 2);
        if (isotopic_shifts[iatm] > 4 || isotopic_shifts[iatm] < -3) {
          rtErr("A V2000 MOL format entry should not describe an isotopic shift outside the range "
                "[-3, 4].  Shift found: " + std::to_string(isotopic_shifts[iatm]) +
                ".  Title of entry: \"" + title + "\".", "MdlMolObj");
        }
      }
      
      // Standard Template Library vector<bool> works differently from other vectors.  Set its
      // contents in a different manner.
      if (verifyContents(tf, i, 36, 3, NumberFormat::INTEGER)) {
        bool dblt_flag;
        formal_charges[iatm] = interpretFormalCharge(readIntegerValue(atom_line_ptr, 36, 3),
                                                     &dblt_flag);
        doublet_radicals[iatm] = dblt_flag;
      }
      if (verifyContents(tf, i, 39, 3, NumberFormat::INTEGER)) {
        parities[iatm] = interpretStereoParity(readIntegerValue(atom_line_ptr, 39, 3));
      }
      if (verifyContents(tf, i, 42, 3, NumberFormat::INTEGER)) {
        implicit_hydrogens[iatm] =
          interpretImplicitHydrogenContent(readIntegerValue(atom_line_ptr, 42, 3));
      }
      if (verifyContents(tf, i, 45, 3, NumberFormat::INTEGER)) {
        stereo_considerations[iatm] =
          interpretBooleanValue(readIntegerValue(atom_line_ptr, 45, 3),
                                "interpreting stereochemical considerations");
      }
      if (verifyContents(tf, i, 48, 3, NumberFormat::INTEGER)) {
        valence_connections[iatm] = interpretValenceNumber(readIntegerValue(atom_line_ptr, 48, 3));
      }
      if (verifyContents(tf, i, 51, 3, NumberFormat::INTEGER)) {
        if (readIntegerValue(atom_line_ptr, 51, 3) == 1 && implicit_hydrogens[iatm] > 0) {
          rtErr("The H0 designation, indicating that implicit hydrogens are not allowed on atom " +
                std::to_string(iatm + 1) + " of MDL MOL entry \"" +  title + "\", is present but "
                "the number of implicit hydrogens has also been indicated as " +
                std::to_string(implicit_hydrogens[iatm]) + ".", "MdlMolObj");
        }
      }
      if (verifyContents(tf, i, 60, 3, NumberFormat::INTEGER)) {
        atom_atom_mapping_count[iatm] = readIntegerValue(atom_line_ptr, 60, 3);
      }
      if (verifyContents(tf, i, 63, 3, NumberFormat::INTEGER)) {
        orientation_stability[iatm] =
          interpretStereoStability(readIntegerValue(atom_line_ptr, 63, 3));
      }
      if (verifyContents(tf, i, 66, 3, NumberFormat::INTEGER)) {
        exact_change_enforced[iatm] = interpretBooleanValue(readIntegerValue(atom_line_ptr, 66, 3),
                                                            "interpreting exact change flag");
      }
      iatm++;
    }
    break;
  case MdlMolVersion::V3000:
    break;
  case MdlMolVersion::UNKNOWN:
    break;
  }

  // Read the bonds block
  switch (version_no) {
  case MdlMolVersion::V2000:
    {
      const int bond_line_start = line_start + 4 + atom_count;
      for (int pos = 0; pos < bond_count; pos++) {
        bonds.emplace_back(tf, bond_line_start + pos, title);
      }
    }
    break;
  case MdlMolVersion::V3000:
    break;
  case MdlMolVersion::UNKNOWN:
    break;
  }

  // Read the atom lists (this information is superceded by the presence of "M  ALS" properties)
  switch (version_no) {
  case MdlMolVersion::V2000:
    {
      // Scan for atom lists imbedded in properties
      const std::string als_tag("M  ALS");
      for (int pos = line_start; pos < mdl_section_end; pos++) {
        if (tfr.line_lengths[pos] >= 6 && strncmpCased(tf.getLinePointer(pos), als_tag)) {
          element_lists.emplace_back(tf, pos, title);
        }
      }
      if (element_lists.size() == 0LLU) {
        const int alst_line_start = line_start + 4 + atom_count + bond_count;
        for (int pos = 0; pos < list_count; pos++) {
          element_lists.emplace_back(tf, alst_line_start + pos, title);          
        }
      }
    }
    break;
  case MdlMolVersion::V3000:
    break;
  case MdlMolVersion::UNKNOWN:
    break;
  }

  // Read the stext entries
  switch (version_no) {
  case MdlMolVersion::V2000:
    {
      const int stxt_line_start = line_start + 4 + atom_count + bond_count + list_count;
      for (int pos = 0; pos < stext_entry_count; pos += 2) {
      }
    }
    break;
  case MdlMolVersion::V3000:
    break;
  case MdlMolVersion::UNKNOWN:
    break;
  }

  // Read various properties
  switch (version_no) {
  case MdlMolVersion::V2000:
    {
      const int prop_line_start = line_start + 4 + atom_count + bond_count + list_count +
                                  (2 * stext_entry_count);
    }
    break;
  case MdlMolVersion::V3000:
    break;
  case MdlMolVersion::UNKNOWN:
    break;
  }

  // Make some basic inferences
  const std::vector<int> tmp_znum = symbolToZNumber(atomic_symbols, capitalization, policy);
  int nvs = 0;
  for (int i = 0; i < atom_count; i++) {
    atomic_numbers[i] = tmp_znum[i];
    nvs += (atomic_numbers[i] == 0);
  }
  if (nvs > 0) {
    
  }
  
}

//-------------------------------------------------------------------------------------------------
MdlMolObj::MdlMolObj(const std::string &filename):
  MdlMolObj(TextFile(filename))
{}

//-------------------------------------------------------------------------------------------------
MdlMolObj::MdlMolObj(const char* filename):
  MdlMolObj(std::string(filename))
{}

//-------------------------------------------------------------------------------------------------
int MdlMolObj::getAtomCount() const {
  return atom_count;
}

//-------------------------------------------------------------------------------------------------
int MdlMolObj::getBondCount() const {
  return bond_count;
}

//-------------------------------------------------------------------------------------------------
double3 MdlMolObj::getCoordinates(const int index) const {
  return coordinates[index];
}

//-------------------------------------------------------------------------------------------------
const std::vector<double3>& MdlMolObj::getCoordinates() const {
  return coordinates;
}

//-------------------------------------------------------------------------------------------------
std::vector<double> MdlMolObj::getCoordinates(const CartesianDimension dim) const {
  std::vector<double> result(atom_count);
  for (int i = 0; i < atom_count; i++) {
    switch (dim) {
    case CartesianDimension::X:
      result[i] = coordinates[i].x;
      break;
    case CartesianDimension::Y:
      result[i] = coordinates[i].y;
      break;
    case CartesianDimension::Z:
      result[i] = coordinates[i].z;
      break;
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
PhaseSpace MdlMolObj::exportPhaseSpace() const {
  PhaseSpace result(atom_count);
  PhaseSpaceWriter rw = result.data();
  for (int i = 0; i < atom_count; i++) {
    rw.xcrd[i] = coordinates[i].x;
    rw.ycrd[i] = coordinates[i].y;
    rw.zcrd[i] = coordinates[i].z;
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
CoordinateFrame MdlMolObj::exportCoordinateFrame() const {
  CoordinateFrame result(atom_count);
  CoordinateFrameWriter rw = result.data();
  for (int i = 0; i < atom_count; i++) {
    rw.xcrd[i] = coordinates[i].x;
    rw.ycrd[i] = coordinates[i].y;
    rw.zcrd[i] = coordinates[i].z;
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
char4 MdlMolObj::getAtomSymbol(const int index) const {
  return atomic_symbols[index];
}

//-------------------------------------------------------------------------------------------------
const std::vector<char4>& MdlMolObj::getAtomSymbols() const {
  return atomic_symbols;
}

//-------------------------------------------------------------------------------------------------
int MdlMolObj::getAtomicNumber(const int index) const {
  return atomic_numbers[index];
}

//-------------------------------------------------------------------------------------------------
const std::vector<int>& MdlMolObj::getAtomicNumbers() const {
  return atomic_numbers;
}

//-------------------------------------------------------------------------------------------------
void MdlMolObj::allocate() {

  // Atom property fields are resized and then set as part of a loop in the parent MdlMolObj
  // constructor.
  coordinates.resize(atom_count);
  atomic_symbols.resize(atom_count, default_mdl_atomic_symbol);
  atomic_numbers.resize(atom_count, default_mdl_atomic_number);
  formal_charges.resize(atom_count, default_mdl_formal_charge);
  doublet_radicals.resize(atom_count, default_mdl_doublet_radical_state);
  isotopic_shifts.resize(atom_count, default_mdl_isotopic_shift);
  parities.resize(atom_count, default_mdl_stereo_parity);
  implicit_hydrogens.resize(atom_count, default_mdl_implicit_hydrogen);
  stereo_considerations.resize(atom_count, default_mdl_stereo_considerations);
  valence_connections.resize(atom_count, default_mdl_valence_connections);
  atom_atom_mapping_count.resize(atom_count, default_mdl_map_count);
  exact_change_enforced.resize(atom_count, default_mdl_exact_change);
  orientation_stability.resize(atom_count, default_mdl_stereo_retention);

  // Other arrays are reserved and built with emplace_back().  The MDL MOL properties array is
  // not reserved to any specific length, however, as the number of properties is not known from
  // the counts line, where it is set to 999 by default.
  bonds.reserve(bond_count);
  element_lists.reserve(list_count);
  stext_entries.reserve(stext_entry_count);
}

//-------------------------------------------------------------------------------------------------
double MdlMolObj::interpretFormalCharge(const int charge_in, bool *is_doublet_radical) {
  *is_doublet_radical = false;
  switch (charge_in) {
  case 0:
    return 0.0;
  case 1:
    return 3.0;
  case 2:
    return 2.0;
  case 3:
    return 1.0;
  case 4:
    *is_doublet_radical = true;
    return 0.0;
  case 5:
    return -1.0;
  case 6:
    return -2.0;
  case 7:
    return -3.0;
  default:
    rtErr("A formal charge code of " + std::to_string(charge_in) + " is invalid for an MDL MOL "
          "entry.  Title of entry: \"" + title + "\".", "MdlMolObj", "interpretFormalCharge");    
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
MolObjAtomStereo MdlMolObj::interpretStereoParity(const int setting_in) {
  switch (setting_in) {
  case 0:
    return MolObjAtomStereo::NOT_STEREO;
  case 1:
    return MolObjAtomStereo::ODD;
  case 2:
    return MolObjAtomStereo::EVEN;
  case 3:
    return MolObjAtomStereo::UNMARKED;
  default:
    rtErr("A stereochemical parity setting of " + std::to_string(setting_in) + " is invalid.  "
          "Title of entry: \"" + title + "\".", "MdlMolObj", "interpretStereoParity");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
int MdlMolObj::interpretImplicitHydrogenContent(const int nh_in) {
  if (nh_in >= 5 || nh_in < 0) {
    rtErr("An implicit hydrogen content of " + std::to_string(nh_in) + " would imply " +
          std::to_string(nh_in - 1) + " hydrogens can be inferred around an atom in MDL MOL "
          "entry \"" + title + "\".", "MdlMolObj", "interpretImplicitHydrogenContent");
  }
  return std::max(nh_in - 1, 0);
}

//-------------------------------------------------------------------------------------------------
bool MdlMolObj::interpretBooleanValue(const int value_in, const std::string &desc) {
  if (value_in != 0 && value_in != 1) {
    rtErr("A directive of " + std::to_string(value_in) + " is invalid when " + desc + ".  Title "
          "of entry: \"" + title + "\".", "MdlMolObj", "interpretBooleanValue");
  }
  return (value_in == 1);
}

//-------------------------------------------------------------------------------------------------
int MdlMolObj::interpretValenceNumber(const int count_in) {
  if (count_in == 15) {
    return 0;
  }
  else if (count_in >= 0 && count_in < 15) {
    return count_in;
  }
  else {
    rtErr("An atom cannot have " + std::to_string(count_in) + " valence connections, as is the "
          "case for one atom in MDL MOL entry \"" + title + "\".", "MdlMolObj",
          "interpretValenceNumber");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
StereoRetention MdlMolObj::interpretStereoStability(const int code_in) {
  switch (code_in) {
  case 0:
    return StereoRetention::NOT_APPLIED;
  case 1:
    return StereoRetention::INVERTED;
  case 2:
    return StereoRetention::RETAINED;
  default:
    rtErr("A stereochemistry retention code of " + std::to_string(code_in) + " is invalid in MDL "
          "MOL entry \"" + title + "\".", "MdlMolObj", "interpretStereoStability");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::vector<MolObjBond> operator+(const std::vector<MolObjBond> &lhs,
				  const std::vector<MolObjBond> &rhs) {
  std::vector<MolObjBond> result(lhs.begin(), lhs.end());
  result.insert(result.end(), rhs.begin(), rhs.end());
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<MdlMolObj> readStructureDataFile(const TextFile &tf, const int low_frame_limit,
                                             const int high_frame_limit,
                                             const CaseSensitivity capitalization,
                                             const ExceptionResponse policy) {
  std::vector<MdlMolObj> result;

  // Find the limits for different MDL MOL entries
  const std::vector<int2> mol_entry_limits = findSdfMolEntryLimits(tf);
  const int nsection = mol_entry_limits.size();
  int actual_low_limit, actual_high_limit;
  if (low_frame_limit >= nsection) {
    rtErr("An SD file with " + std::to_string(nsection) + " frames cannot be read starting at "
          "frame index " + std::to_string(low_frame_limit) + ".", "readStructureDataFile");
  }
  else if (low_frame_limit < 0 || high_frame_limit < 0 || high_frame_limit >= nsection) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("The frame range " + std::to_string(low_frame_limit) + " to " +
            std::to_string(high_frame_limit) + " is invalid for a file with " +
            std::to_string(nsection) + " frames.", "readStructureDataFile");
    case ExceptionResponse::WARN:
      rtWarn("The frame range " + std::to_string(low_frame_limit) + " to " +
             std::to_string(high_frame_limit) + " is invalid for a file with " +
             std::to_string(nsection) + " frames.  Only the valid range will be taken",
             "readStructureDataFile");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
    actual_low_limit = std::max(0, low_frame_limit);
    actual_high_limit = (high_frame_limit < low_frame_limit || high_frame_limit >= nsection) ?
                        nsection - 1 : high_frame_limit;
  }
  else {
    actual_low_limit = low_frame_limit;
    actual_high_limit = (high_frame_limit < low_frame_limit) ? nsection - 1: high_frame_limit;
  }
  result.reserve(actual_high_limit - actual_low_limit + 1);
  for (int i = actual_low_limit; i <= actual_high_limit; i++) {
    result.emplace_back(tf, mol_entry_limits[i].x, mol_entry_limits[i].y, capitalization, policy);
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<MdlMolObj> readStructureDataFile(const TextFile &tf,
                                             const CaseSensitivity capitalization,
                                             const ExceptionResponse policy) {
  std::vector<MdlMolObj> result;

  // Find the limits for different MDL MOL entries
  const std::vector<int2> mol_entry_limits = findSdfMolEntryLimits(tf);

  // Parse each MDL MOL entry
  const int nsection = mol_entry_limits.size();
  result.reserve(nsection);
  for (int i = 0; i < nsection; i++) {
    result.emplace_back(tf, mol_entry_limits[i].x, mol_entry_limits[i].y, capitalization, policy);
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<MdlMolObj> readStructureDataFile(const std::string &file_name,
                                             const CaseSensitivity capitalization,
                                             const ExceptionResponse policy) {
  const TextFile tf(file_name);
  return readStructureDataFile(tf, capitalization, policy);
}

} // namespace structure
} // namespace stormm
