#include "copyright.h"
#include "mdl_mol_format.h"

namespace stormm {
namespace structure {

using parse::NumberFormat;
using parse::readIntegerValue;
using parse::readRealValue;
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
MdlMolObj::MdlMolObj():
    atom_count{0}, bond_count{0}, list_count{0}, sgroup_count{0}, constraint_count{0},
    chirality{MolObjChirality::ACHIRAL}, registry_number{-1}, coordinates{}, atomic_symbols{},
    atomic_numbers{}, formal_charges{}, isotopic_shifts{}, parities{}, implicit_hydrogens{},
    stereo_considerations{}, valence_connections{}, atom_atom_mapping_count{},
    orientation_stability{}, bonds{}, properties{}, title{""}, software_details{""},
    general_comment{""}
{}

//-------------------------------------------------------------------------------------------------
MdlMolObj::MdlMolObj(const TextFile &tf, const int line_start, const int line_end_in,
                     const CaseSensitivity capitalization, const ExceptionResponse policy):
  MdlMolObj()
{
  // Default line end of -1 indicates reading to the end of the file
  int line_end = (line_end_in < 0) ? tf.getLineCount() : line_end_in;

  // Identify the next delimiter (a $$$$ marking at the beginning of the line)
  const TextFileReader tfr = tf.data();
  int next_delim_loc = line_start;
  bool found = false;
  while (found == false && next_delim_loc < line_end) {
    const char* lptr = &tfr.text[tfr.line_limits[next_delim_loc]];
    if (tfr.line_limits[next_delim_loc + 1] - tfr.line_limits[next_delim_loc] >= 6 &&
        lptr[0] == 'M' && lptr[1] == ' ' && lptr[2] == ' ' && lptr[3] == 'E' && lptr[4] == 'N' &&
        lptr[5] == 'D') {
      found = true;
    }
    else {
      next_delim_loc++;
    }
  }
  line_end = next_delim_loc;

  // The range of data now extends from line_start to line_end.  Sift through that information
  // for a V2000 or V3000 specification.  This should be found on the fourth line.
  int version = findMolObjVersion(&tfr.text[tfr.line_limits[line_start + 3]],
                                  tfr.line_limits[line_start + 4] -
                                  tfr.line_limits[line_start + 3]);
  
  // Begin by reading the molecule name (title), generating software details, and any general
  // comment (always three and only three distinct lines, even if left blank).
  title            = tf.extractString(line_start);
  software_details = tf.extractString(line_start + 1);
  general_comment  = tf.extractString(line_start + 2);

  // Read the counts line
  if (version == 2000) {
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

    // Allocate space for information to be read
    allocate();
  }
    
  // Read the atoms block
  if (version == 2000) {
    int iatm = 0;
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
  }

  // Read the bonds block
  if (version == 2000) {
    const int bond_line_start = line_start + 4 + atom_count;
    for (int pos = 0; pos < bond_count; pos++) {
      const char* bond_line_ptr = &tfr.text[tfr.line_limits[bond_line_start + pos]];
      const MolObjBondOrder bnd_order = (verifyContents(tf, bond_line_start + pos, 6, 3)) ?
                                        interpretBondOrder(readIntegerValue(bond_line_ptr, 6, 3)) :
                                        default_mdl_bond_order;
      const MolObjBondStereo bnd_stereo = (verifyContents(tf, bond_line_start + pos, 9, 3)) ?
        interpretBondStereochemistry(readIntegerValue(bond_line_ptr, 9, 3)) :
        default_mdl_bond_stereochemistry;
      const MolObjRingState ring_status = (verifyContents(tf, bond_line_start + pos, 12, 3)) ?
        interpretRingState(readIntegerValue(bond_line_ptr, 12, 3)) :
        default_mdl_ring_status;
      const MolObjReactionCenter reax_potential = (verifyContents(tf, bond_line_start + pos,
                                                                  15, 3)) ?
        interpretBondReactivePotential(readIntegerValue(bond_line_ptr, 15, 3)) :
        default_mdl_bond_reactivity;
      bonds[pos] = MolObjBond(readIntegerValue(bond_line_ptr, 0, 3),
                              readIntegerValue(bond_line_ptr, 3, 3), bnd_order, bnd_stereo,
                              ring_status, reax_potential);
    }
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
  
  // Read various properties

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
PhaseSpace MdlMolObj::exportPhaseSpace() {
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
CoordinateFrame MdlMolObj::exportCoordinateFrame() {
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
  bonds.resize(bond_count);
  properties.resize(properties_count);
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
MolObjBondOrder MdlMolObj::interpretBondOrder(const int code_in) {
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
MolObjBondStereo MdlMolObj::interpretBondStereochemistry(const int code_in) {
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
MolObjRingState MdlMolObj::interpretRingState(const int code_in) {
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
MolObjReactionCenter MdlMolObj::interpretBondReactivePotential(const int code_in) {
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
std::vector<MolObjBond> operator+(const std::vector<MolObjBond> &lhs,
				  const std::vector<MolObjBond> &rhs) {
  std::vector<MolObjBond> result(lhs.begin(), lhs.end());
  result.insert(result.end(), rhs.begin(), rhs.end());
  return result;
}

//-------------------------------------------------------------------------------------------------
int findMolObjVersion(const char* text, const int nchar) {
  for (int i = 0; i < nchar; i++) {
    if (text[i] == 'V' || text[i] == 'v') {
      if (i < nchar - 4 && text[i + 2] == '0' && text[i + 3] == '0' && text[i + 4] == '0') {
        if (text[i + 1] == '2') {
          return 2000;
        }
        else if (text[i + 1] == '3') {
          return 3000;
        }
      }
    }
  }
  return 2000;
}

//-------------------------------------------------------------------------------------------------
std::vector<MdlMolObj> readStructureDataFile(const TextFile &tf,
                                             const CaseSensitivity capitalization,
                                             const ExceptionResponse policy) {

}
                                             
} // namespace structure
} // namespace stormm
