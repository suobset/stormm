#include "copyright.h"
#include "FileManagement/file_listing.h"
#include "mdl_mol_format.h"

namespace stormm {
namespace structure {

using diskutil::getBaseName;
using parse::NumberFormat;
using parse::readIntegerValue;
using parse::readRealValue;
using parse::separateText;
using parse::stringToChar4;
using parse::strncmpCased;
using parse::TextFileReader;
using parse::verifyContents;
using parse::operator==;
using trajectory::CoordinateFrameWriter;
using trajectory::PhaseSpaceWriter;
  
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
  const int sd_compound_end = getCompoundSectionEnd(tfr, line_start, line_end_in);

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
                " in .sdf or MDL MOL file " + getBaseName(tf.getFileName()) + ".", "MdlMolObj");
        }
      }
      if (verifyContents(tf, counts_line_idx, 15, 3, NumberFormat::INTEGER)) {
        stext_entry_count = readIntegerValue(counts_line_ptr, 15, 3);
      }
      if (verifyContents(tf, counts_line_idx, 30, 3, NumberFormat::INTEGER)) {
        properties_count  = readIntegerValue(counts_line_ptr, 30, 3);
      }

      // Validation of the atom counts line
      if (atom_count <= 0 || bond_count < 0) {
        switch (policy) {
        case ExceptionResponse::DIE:
          rtErr("A V2000 MOL format entry in file " + getBaseName(tf.getFileName()) + " at line " +
                std::to_string(line_start) + " contains invalid numbers of atoms (" +
                std::to_string(atom_count) + ") and/or bonds (" + std::to_string(bond_count) +
                ").", "MdlMolObj");
        case ExceptionResponse::WARN:
          rtWarn("A V2000 MOL format entry in file " + getBaseName(tf.getFileName()) + " at " +
                 "line " + std::to_string(line_start) + " contains invalid numbers of atoms (" +
                 std::to_string(atom_count) + ") and/or bonds (" + std::to_string(bond_count) +
                 ").  This will become a blank entry", "MdlMolObj");
          break;
        case ExceptionResponse::SILENT:
          break;
        }
        atom_count = 0;
        bond_count = 0;
        list_count = 0;
        chirality = MolObjChirality::ACHIRAL;
        stext_entry_count = 0;
        properties_count = 0;
        return;
      }
    }
    break;
  case MdlMolVersion::V3000:
    break;
  case MdlMolVersion::UNKNOWN:

    // This case would be unprocessable.
    rtErr("No valid MDL MOL version was detected.  Parsing in " + getBaseName(tf.getFileName()) +
          " cannot proceed.");
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
        interpretFormalCharge(readIntegerValue(atom_line_ptr, 36, 3), iatm);
      }
      if (verifyContents(tf, i, 39, 3, NumberFormat::INTEGER)) {
        parities[iatm] = interpretStereoParity(readIntegerValue(atom_line_ptr, 39, 3));
      }
      if (verifyContents(tf, i, 42, 3, NumberFormat::INTEGER)) {
        interpretImplicitHydrogenContent(readIntegerValue(atom_line_ptr, 42, 3), iatm);
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
      for (int pos = prop_line_start; pos < mdl_section_end; pos++) {
        int adv_pos;
        properties.emplace_back(tf, pos, &adv_pos, title);
      }
      
      // Update the properties count
      properties_count = properties.size();
    }
    break;
  case MdlMolVersion::V3000:
    break;
  case MdlMolVersion::UNKNOWN:
    break;
  }

  // Read the data items
  for (int pos = mdl_section_end; pos < sd_compound_end; pos++) {
    if (tf.getLineLength(pos) >= 2 && tf.getChar(tf.getLineLimits(pos)) == '>') {
      int adv_pos;
      data_items.emplace_back(tf, pos, &adv_pos, sd_compound_end, title);
      pos = adv_pos;
    }
  }
  if (data_items.size() == 0LLU && sd_compound_end - mdl_section_end > 1) {
    rtErr("If there are no data items, the compound section must terminate immediately after the "
          "MDL MOL format section.  File " + getBaseName(tf.getFileName()) + " violates SD file "
          "conventions at lines " + std::to_string(mdl_section_end) + " - " +
          std::to_string(sd_compound_end));
  }

  // CHECK
  printf("There are %zu data items.\n", data_items.size());
  for (size_t i = 0; i < data_items.size(); i++) {
    printf("Data item %2zu : > <%s>(%s)\n", i, data_items[i].getItemName().c_str(),
           data_items[i].getExternalRegistryNumber().c_str());
  }
  // END CHECK
  
  // Update the atom attributes based on properties data.  This provides V3000 functionality and
  // backwards compatibility for the V2000 format.
  updateV2kAtomAttributes();
  
  // Make some basic inferences.  The version number is irrelevant by now: all information has been
  // converted into version-agnostic internal representations and the version number is kept only
  // as a footnote for reference when re-printing the file later.
  const std::vector<int> tmp_znum = symbolToZNumber(atomic_symbols, capitalization, policy);
  int nvs = 0;
  for (int i = 0; i < atom_count; i++) {
    atomic_numbers[i] = tmp_znum[i];
    nvs += (atomic_numbers[i] == 0);
  }
  if (nvs > 0) {
    
  }

  // Add implicit hydrogens.  This may re-allocate the data arrays and extend the bonding patterns.
  hydrogenate();
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
int MdlMolObj::getFormalCharge(const int index) const {
  return formal_charges[index];
}

//-------------------------------------------------------------------------------------------------
const std::vector<int>& MdlMolObj::getFormalCharges() const {
  return formal_charges;
}

//-------------------------------------------------------------------------------------------------
int MdlMolObj::getPropertiesCount() const {
  return properties_count;
}

//-------------------------------------------------------------------------------------------------
void MdlMolObj::allocate() {

  // Atom property fields are resized and then set as part of a loop in the parent MdlMolObj
  // constructor.
  coordinates.resize(atom_count);
  atomic_symbols.resize(atom_count, default_mdl_atomic_symbol);
  atomic_numbers.resize(atom_count, default_mdl_atomic_number);
  formal_charges.resize(atom_count, default_mdl_formal_charge);
  radicals.resize(atom_count, default_mdl_radical_state);
  isotopic_shifts.resize(atom_count, default_mdl_isotopic_shift);
  parities.resize(atom_count, default_mdl_stereo_parity);
  implicit_hydrogens.resize(atom_count, default_mdl_implicit_hydrogen);
  stereo_considerations.resize(atom_count, default_mdl_stereo_considerations);
  valence_connections.resize(atom_count, default_mdl_valence_connections);
  atom_atom_mapping_count.resize(atom_count, default_mdl_map_count);
  exact_change_enforced.resize(atom_count, default_mdl_exact_change);
  hydrogenation_protocol.resize(atom_count, default_hydrogenation);
  orientation_stability.resize(atom_count, default_mdl_stereo_retention);

  // Other arrays are reserved and built with emplace_back().  The MDL MOL properties array is
  // not reserved to any specific length, however, as the number of properties is not known from
  // the counts line, where it is set to 999 by default.
  bonds.reserve(bond_count);
  element_lists.reserve(list_count);
  stext_entries.reserve(stext_entry_count);
}

//-------------------------------------------------------------------------------------------------
void MdlMolObj::interpretFormalCharge(const int charge_in, const int atom_index) {
  switch (charge_in) {
  case 0:
    formal_charges[atom_index] = 0;
    break;
  case 1:
    formal_charges[atom_index] = 3;
    break;
  case 2:
    formal_charges[atom_index] = 2;
    break;
  case 3:
    formal_charges[atom_index] = 1;
    break;
  case 4:
    radicals[atom_index] = RadicalState::DOUBLET;
    break;
  case 5:
    formal_charges[atom_index] = -1;
    break;
  case 6:
    formal_charges[atom_index] = -2;
    break;
  case 7:
    formal_charges[atom_index] = -3;
    break;
  default:
    rtErr("A formal charge code of " + std::to_string(charge_in) + " is invalid for an MDL MOL "
          "entry.  Title of entry: \"" + title + "\".", "MdlMolObj", "interpretFormalCharge");    
  }
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
void MdlMolObj::interpretImplicitHydrogenContent(const int nh_in, const int atom_index) {
  if (nh_in > 5 || nh_in < 0) {
    rtErr("An implicit hydrogen content of " + std::to_string(nh_in) + " would imply " +
          std::to_string(nh_in - 1) + " hydrogens can be inferred around an atom in MDL MOL "
          "entry \"" + title + "\".", "MdlMolObj", "interpretImplicitHydrogenContent");
  }
  if (nh_in > 0) {
    implicit_hydrogens[atom_index] = nh_in - 1;
    hydrogenation_protocol[atom_index] = (nh_in == 1) ? HydrogenAssignment::DO_NOT_HYDROGENATE :
                                                        HydrogenAssignment::VALENCE_SHELL;
  }
  else {

    // An implicit hydrogen indicator of 0 does not correspond to H0, H1, H2, H3, or H4, but it is
    // very common.  This final possibility implies "free hydrogen content." While the actual
    // number of hydrogens will read 0, the flag will be set to apply as many as are needed to
    // fill the valence shell given the bonding considerations.
    implicit_hydrogens[atom_index] = 0;
    hydrogenation_protocol[atom_index] = HydrogenAssignment::VALENCE_SHELL;
  }
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
void MdlMolObj::updateV2kAtomAttributes() {

  // Return immediately if the version is not V2000
  switch (version_no) {
  case MdlMolVersion::V2000:
    break;
  case MdlMolVersion::V3000:
  case MdlMolVersion::UNKNOWN:
    return;
  }
  
  // Scan all properties for items that would invalidate atom-block information.  Wipe the
  // relevant arrays.
  for (int i = 0; i < properties_count; i++) {
    if (properties[i].getCode() == char4({ 'C', 'H', 'G', 'M' })) {
      for (int j = 0; j < atom_count; j++) {
        formal_charges[j] = 0;
      }
    }
    if (properties[i].getCode() == char4({ 'R', 'A', 'D', 'M' })) {
      for (int j = 0; j < atom_count; j++) {
        formal_charges[j] = 0;
        radicals[j] = RadicalState::NONE;
      }
    }
    if (properties[i].getCode() == char4({ 'I', 'S', 'O', 'M' })) {
      for (int j = 0; j < atom_count; j++) {
        isotopic_shifts[j] = 0;
      }
    }
  }

  // Scan the properties again and add details to the atoms.
  for (int i = 0; i < properties_count; i++) {
    const int n_entry = properties[i].getEntryCount();
    if (properties[i].getCode() == char4({ 'C', 'H', 'G', 'M' })) {
      for (int j = 0; j < n_entry; j++) {
        formal_charges[properties[i].getIntegerValue(j, 0)] = properties[i].getIntegerValue(j, 1);
      }
    }
    if (properties[i].getCode() == char4({ 'R', 'A', 'D', 'M' })) {
      for (int j = 0; j < n_entry; j++) {
        const int atom_idx = properties[i].getIntegerValue(j, 0);
        switch (properties[i].getIntegerValue(j, 1)) {
        case 0:
          radicals[atom_idx] = RadicalState::NONE;
          break;
        case 1:
          radicals[atom_idx] = RadicalState::SINGLET;
          break;
        case 2:
          radicals[atom_idx] = RadicalState::DOUBLET;
          break;
        case 3:
          radicals[atom_idx] = RadicalState::TRIPLET;
          break;
        }
      }
    }
    if (properties[i].getCode() == char4({ 'I', 'S', 'O', 'M' })) {
      for (int j = 0; j < n_entry; j++) {
        isotopic_shifts[properties[i].getIntegerValue(j, 0)] = properties[i].getIntegerValue(j, 1);
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------
void MdlMolObj::hydrogenate() {

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
