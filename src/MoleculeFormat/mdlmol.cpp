#include <cmath>
#include "copyright.h"
#include "FileManagement/file_listing.h"
#include "Topology/atomgraph_abstracts.h"
#include "Topology/atomgraph_enumerators.h"
#include "Trajectory/write_annotated_frame.h"
#include "mdlmol.h"

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
using parse::TextOrigin;
using parse::verifyContents;
using parse::operator==;
using topology::ChemicalDetailsKit;
using topology::TorsionKind;
using topology::ValenceKit;
using trajectory::writeFrame;

//-------------------------------------------------------------------------------------------------
MdlMol::MdlMol(const ExceptionResponse policy_in):
    policy{policy_in}, version_no{MdlMolVersion::V2000}, atom_count{0}, bond_count{0},
    list_count{0}, sgroup_count{0}, constraint_count{0}, chirality{MolObjChirality::ACHIRAL},
    registry_number{-1}, data_item_count{0}, property_formal_charges{false},
    property_radicals{false}, property_isotopes{false}, property_element_lists{false},
    coordinates{}, atomic_symbols{}, atomic_numbers{}, formal_charges{}, isotopic_shifts{},
    parities{}, implicit_hydrogens{}, stereo_considerations{}, valence_connections{},
    atom_atom_mapping_count{}, orientation_stability{}, bonds{}, element_lists{}, stext_entries{},
    properties{}, title{""}, software_details{""}, general_comment{""}
{}

//-------------------------------------------------------------------------------------------------
MdlMol::MdlMol(const TextFile &tf, const int line_start, const int line_end_in,
               const CaseSensitivity capitalization, const ExceptionResponse policy_in):
  MdlMol(policy_in)
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
                " in .sdf or MDL MOL file " + getBaseName(tf.getFileName()) + ".", "MdlMol");
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
                ").", "MdlMol");
        case ExceptionResponse::WARN:
          rtWarn("A V2000 MOL format entry in file " + getBaseName(tf.getFileName()) + " at " +
                 "line " + std::to_string(line_start) + " contains invalid numbers of atoms (" +
                 std::to_string(atom_count) + ") and/or bonds (" + std::to_string(bond_count) +
                 ").  This will become a blank entry", "MdlMol");
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
    
  // Read the atom block
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
                ".  Title of entry: \"" + title + "\".", "MdlMol");
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
                std::to_string(implicit_hydrogens[iatm]) + ".", "MdlMol");
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
      compareExternalRegistryNumbers(data_items.back().getExternalRegistryNumber());
      pos = adv_pos;
    }
  }
  data_item_count = data_items.size();
  if (data_item_count == 0 && sd_compound_end - mdl_section_end > 1) {
    rtErr("If there are no data items, the compound section must terminate immediately after the "
          "MDL MOL format section.  File " + getBaseName(tf.getFileName()) + " violates SD file "
          "conventions at lines " + std::to_string(mdl_section_end) + " - " +
          std::to_string(sd_compound_end));
  }
  
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
MdlMol::MdlMol(const std::string &filename, const ExceptionResponse policy_in):
    MdlMol(TextFile(filename), 0, -1, CaseSensitivity::YES, policy_in)
{}

//-------------------------------------------------------------------------------------------------
MdlMol::MdlMol(const char* filename, const ExceptionResponse policy_in):
    MdlMol(std::string(filename), policy_in)
{}

//-------------------------------------------------------------------------------------------------
int MdlMol::getAtomCount() const {
  return atom_count;
}

//-------------------------------------------------------------------------------------------------
int MdlMol::getBondCount() const {
  return bond_count;
}

//-------------------------------------------------------------------------------------------------
double3 MdlMol::getCoordinates(const int index) const {
  return coordinates[index];
}

//-------------------------------------------------------------------------------------------------
const std::vector<double3>& MdlMol::getCoordinates() const {
  return coordinates;
}

//-------------------------------------------------------------------------------------------------
std::vector<double> MdlMol::getCoordinates(const CartesianDimension dim) const {
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
PhaseSpace MdlMol::exportPhaseSpace() const {
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
CoordinateFrame MdlMol::exportCoordinateFrame() const {
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
char4 MdlMol::getAtomSymbol(const int index) const {
  return atomic_symbols[index];
}

//-------------------------------------------------------------------------------------------------
const std::vector<char4>& MdlMol::getAtomSymbols() const {
  return atomic_symbols;
}

//-------------------------------------------------------------------------------------------------
int MdlMol::getAtomicNumber(const int index) const {
  return atomic_numbers[index];
}

//-------------------------------------------------------------------------------------------------
const std::vector<int>& MdlMol::getAtomicNumbers() const {
  return atomic_numbers;
}

//-------------------------------------------------------------------------------------------------
int MdlMol::getFormalCharge(const int index) const {
  return formal_charges[index];
}

//-------------------------------------------------------------------------------------------------
const std::vector<int>& MdlMol::getFormalCharges() const {
  return formal_charges;
}

//-------------------------------------------------------------------------------------------------
int MdlMol::getPropertiesCount() const {
  return properties_count;
}

//-------------------------------------------------------------------------------------------------
void MdlMol::impartCoordinates(const PhaseSpaceReader &psr) {
  checkAtomCount(psr.natom);
  impartCoordinates(psr.xcrd, psr.ycrd, psr.zcrd, 1.0);
}

//-------------------------------------------------------------------------------------------------
void MdlMol::impartCoordinates(const PhaseSpaceWriter &psw) {
  checkAtomCount(psw.natom);
  impartCoordinates(psw.xcrd, psw.ycrd, psw.zcrd, 1.0);
}

//-------------------------------------------------------------------------------------------------
void MdlMol::impartCoordinates(const PhaseSpace &ps, const CoordinateCycle orientation,
                               const HybridTargetLevel tier) {
  impartCoordinates(ps.data(orientation, tier));
}

//-------------------------------------------------------------------------------------------------
void MdlMol::impartCoordinates(const PhaseSpace &ps, const HybridTargetLevel tier) {
  impartCoordinates(ps, ps.getCyclePosition(), tier);
}

//-------------------------------------------------------------------------------------------------
void MdlMol::impartCoordinates(const CoordinateFrameReader &cfr) {
  checkAtomCount(cfr.natom);
  impartCoordinates(cfr.xcrd, cfr.ycrd, cfr.zcrd, 1.0);
}

//-------------------------------------------------------------------------------------------------
void MdlMol::impartCoordinates(const CoordinateFrameWriter &cfw) {
  checkAtomCount(cfw.natom);
  impartCoordinates(cfw.xcrd, cfw.ycrd, cfw.zcrd, 1.0);
}

//-------------------------------------------------------------------------------------------------
void MdlMol::impartCoordinates(const CoordinateFrame &cf, const HybridTargetLevel tier) {
  impartCoordinates(cf.data());
}

//-------------------------------------------------------------------------------------------------
void MdlMol::addDataItem(const MolObjDataRequest &ask, const AtomGraph &ag,
                         const RestraintApparatus &ra) {
  std::vector<std::string> di_body;
  switch (ask.getKind()) {
  case DataRequestKind::STATE_VARIABLE:

    // The text will a single number, supplied by an external function or re-evaluated.
    break;
  case DataRequestKind::ATOM_INFLUENCES:

    // The text will be constructed during the energy re-evaluation when the data item is printed.
    break;
  case DataRequestKind::TOPOLOGY_PARAMETER:
    {
      // Search for all instances of a particular parameter type affecting the stated atom types
      // in the given order (or the reverse order).  Each instance of the parameter in the topology
      // will become one data line of the data item.  If no instances are found, a single data line
      // is entered to indicate this fact.  The first string in the body will be reserved to hold
      // a message indicating however many instances were found.
      int nfound = 0;
      di_body.push_back("");
      char buffer_data[256];
      const ValenceKit<double> vk  = ag.getDoublePrecisionValenceKit();
      const ChemicalDetailsKit cdk = ag.getChemicalDetailsKit();
      const std::vector<char4> atom_types = ask.getAtomTypes();
      switch (ask.getValenceParameter()) {
      case StateVariable::BOND:
        {
          const char4 i_tref = atom_types[0];
          const char4 j_tref = atom_types[1];
          for (int pos = 0; pos < vk.nbond; pos++) {
            const int i_atom = vk.bond_i_atoms[pos];
            const int j_atom = vk.bond_j_atoms[pos];
            const char4 i_type = cdk.atom_types[i_atom];
            const char4 j_type = cdk.atom_types[j_atom];
            if ((i_type == i_tref && j_type == j_tref) || (i_type == j_tref && j_type == i_tref)) {
              const bool forward = (i_type == i_tref);
              const char4 i_name = (forward) ? cdk.atom_names[i_atom] : cdk.atom_names[j_atom];
              const char4 j_name = (forward) ? cdk.atom_names[j_atom] : cdk.atom_names[i_atom];
              const int param_idx = vk.bond_param_idx[pos];
              if (nfound == 0) {
                sprintf(buffer_data, " Atom Names   Indices     Eq. L0    Stiffness\n");
                di_body.push_back(buffer_data);
                sprintf(buffer_data, "  ---- ----  ---- ----  ---------- ----------\n");
                di_body.push_back(buffer_data);
              }
              const int i_index = (forward) ? i_atom + 1 : j_atom + 1;
              const int j_index = (forward) ? j_atom + 1 : i_atom + 1;
              sprintf(buffer_data, "  %c%c%c%c %c%c%c%c  %4d %4d  %10.4lf %10.4lf\n", i_name.x,
                      i_name.y, i_name.z, i_name.w, j_name.x, j_name.y, j_name.z, j_name.w,
                      i_index, j_index, vk.bond_leq[param_idx], vk.bond_keq[param_idx]);
              di_body.push_back(buffer_data);
              nfound++;
            }
          }
        }
        break;
      case StateVariable::ANGLE:
        {
          const char4 i_tref = atom_types[0];
          const char4 j_tref = atom_types[1];
          const char4 k_tref = atom_types[2];
          for (int pos = 0; pos < vk.nangl; pos++) {
            const int i_atom = vk.angl_i_atoms[pos];
            const int j_atom = vk.angl_j_atoms[pos];
            const int k_atom = vk.angl_k_atoms[pos];
            const char4 i_type = cdk.atom_types[i_atom];
            const char4 j_type = cdk.atom_types[j_atom];
            const char4 k_type = cdk.atom_types[k_atom];
            if ((i_type == i_tref && j_type == j_tref && k_type == k_tref) ||
                (i_type == k_tref && j_type == j_tref && k_type == i_tref)) {
              const bool forward = (i_type == i_tref);
              const char4 i_name = (forward) ? cdk.atom_names[i_atom] : cdk.atom_names[k_atom];
              const char4 j_name = cdk.atom_names[j_atom];
              const char4 k_name = (forward) ? cdk.atom_names[k_atom] : cdk.atom_names[i_atom];
              const int param_idx = vk.angl_param_idx[pos];
              if (nfound == 0) {
                sprintf(buffer_data, "    Atom Names        Indices       Eq. L0    Stiffness\n");
                di_body.push_back(buffer_data);
                sprintf(buffer_data, "  ---- ---- ----  ---- ---- ----  ---------- ----------\n");
                di_body.push_back(buffer_data);
              }
              const int i_index = (forward) ? i_atom + 1 : k_atom + 1;
              const int j_index = j_atom + 1;
              const int k_index = (forward) ? k_atom + 1 : i_atom + 1;
              sprintf(buffer_data, "  %c%c%c%c %c%c%c%c %c%c%c%c  %4d %4d %4d  %10.4lf %10.4lf\n",
                      i_name.x, i_name.y, i_name.z, i_name.w, j_name.x, j_name.y, j_name.z,
                      j_name.w, k_name.x, k_name.y, k_name.z, k_name.w, i_index, j_index, k_index,
                      vk.angl_theta[param_idx] * 180.0 / symbols::pi, vk.angl_keq[param_idx]);
              di_body.push_back(buffer_data);
              nfound++;
            }
          }
        }
        break;
      case StateVariable::PROPER_DIHEDRAL:
      case StateVariable::IMPROPER_DIHEDRAL:
        {
          const char4 i_tref = atom_types[0];
          const char4 j_tref = atom_types[1];
          const char4 k_tref = atom_types[2];
          const char4 l_tref = atom_types[3];
          const bool is_proper = (ask.getValenceParameter() == StateVariable::PROPER_DIHEDRAL);
          for (int pos = 0; pos < vk.ndihe; pos++) {
            switch (static_cast<TorsionKind>(vk.dihe_modifiers[pos].w)) {
            case TorsionKind::PROPER:
            case TorsionKind::PROPER_NO_14:
              if (is_proper == false) {
                continue;
              }
              break;
            case TorsionKind::IMPROPER:
            case TorsionKind::IMPROPER_NO_14:
              if (is_proper) {
                continue;
              }
            }
            const int i_atom = vk.dihe_i_atoms[pos];
            const int j_atom = vk.dihe_j_atoms[pos];
            const int k_atom = vk.dihe_k_atoms[pos];
            const int l_atom = vk.dihe_l_atoms[pos];
            const char4 i_type = cdk.atom_types[i_atom];
            const char4 j_type = cdk.atom_types[j_atom];
            const char4 k_type = cdk.atom_types[k_atom];
            const char4 l_type = cdk.atom_types[l_atom];
            if ((i_type == i_tref && j_type == j_tref && k_type == k_tref && l_type == l_tref) ||
                (i_type == l_tref && j_type == k_tref && k_type == j_tref && l_type == i_tref)) {
              const bool forward = (i_type == i_tref);
              const char4 i_name = (forward) ? cdk.atom_names[i_atom] : cdk.atom_names[l_atom];
              const char4 j_name = (forward) ? cdk.atom_names[j_atom] : cdk.atom_names[k_atom];
              const char4 k_name = (forward) ? cdk.atom_names[k_atom] : cdk.atom_names[j_atom];
              const char4 l_name = (forward) ? cdk.atom_names[l_atom] : cdk.atom_names[i_atom];
              const int param_idx = vk.dihe_param_idx[pos];
              if (nfound == 0) {
                sprintf(buffer_data, "       Atom Names            Indices         Amplitude "
                        "   Phase    N\n");
                di_body.push_back(buffer_data);
                sprintf(buffer_data, "  ---- ---- ---- ----  ---- ---- ---- ----  ---------- "
                        "---------- --\n");
                di_body.push_back(buffer_data);
              }
              const int i_index = (forward) ? i_atom + 1 : l_atom + 1;
              const int j_index = (forward) ? j_atom + 1 : k_atom + 1;
              const int k_index = (forward) ? k_atom + 1 : j_atom + 1;
              const int l_index = (forward) ? l_atom + 1 : i_atom + 1;
              sprintf(buffer_data, "  %c%c%c%c %c%c%c%c %c%c%c%c %c%c%c%c  %4d %4d %4d %4d  "
                      "%10.4lf %10.4lf %2d\n", i_name.x, i_name.y, i_name.z, i_name.w, j_name.x,
                      j_name.y, j_name.z, j_name.w, k_name.x, k_name.y, k_name.z, k_name.w,
                      l_name.x, l_name.y, l_name.z, l_name.w, i_index, j_index, k_index, l_index,
                      vk.dihe_amp[param_idx], vk.dihe_phi[param_idx] * 180.0 / symbols::pi,
                      static_cast<int>(std::round(vk.dihe_freq[param_idx])));
              di_body.push_back(buffer_data);
              nfound++;
            }
          }
        }
        break;
      case StateVariable::UREY_BRADLEY:
        {
          const char4 i_tref = atom_types[0];
          const char4 k_tref = atom_types[1];
          for (int pos = 0; pos < vk.nubrd; pos++) {
            const int i_atom = vk.ubrd_i_atoms[pos];
            const int k_atom = vk.ubrd_k_atoms[pos];
            const char4 i_type = cdk.atom_types[i_atom];
            const char4 k_type = cdk.atom_types[k_atom];
            if ((i_type == i_tref && k_type == k_tref) || (i_type == k_tref && k_type == i_tref)) {
              const bool forward = (i_type == i_tref);
              const char4 i_name = (forward) ? cdk.atom_names[i_atom] : cdk.atom_names[k_atom];
              const char4 k_name = (forward) ? cdk.atom_names[k_atom] : cdk.atom_names[i_atom];
              const int param_idx = vk.ubrd_param_idx[pos];
              if (nfound == 0) {
                sprintf(buffer_data, " Atom Names   Indices     Eq. L0    Stiffness\n");
                di_body.push_back(buffer_data);
                sprintf(buffer_data, "  ---- ----  ---- ----  ---------- ----------\n");
                di_body.push_back(buffer_data);
              }
              const int i_index = (forward) ? i_atom + 1 : k_atom + 1;
              const int k_index = (forward) ? k_atom + 1 : i_atom + 1;
              sprintf(buffer_data, "  %c%c%c%c %c%c%c%c  %4d %4d  %10.4lf %10.4lf\n", i_name.x,
                      i_name.y, i_name.z, i_name.w, k_name.x, k_name.y, k_name.z, k_name.w,
                      i_index, k_index, vk.ubrd_leq[param_idx], vk.ubrd_keq[param_idx]);
              di_body.push_back(buffer_data);
              nfound++;
            }
          }
        }
        break;
      case StateVariable::CHARMM_IMPROPER:
        {
          const char4 i_tref = atom_types[0];
          const char4 j_tref = atom_types[1];
          const char4 k_tref = atom_types[2];
          const char4 l_tref = atom_types[3];
          for (int pos = 0; pos < vk.ndihe; pos++) {
            const int i_atom = vk.cimp_i_atoms[pos];
            const int j_atom = vk.cimp_j_atoms[pos];
            const int k_atom = vk.cimp_k_atoms[pos];
            const int l_atom = vk.cimp_l_atoms[pos];
            const char4 i_type = cdk.atom_types[i_atom];
            const char4 j_type = cdk.atom_types[j_atom];
            const char4 k_type = cdk.atom_types[k_atom];
            const char4 l_type = cdk.atom_types[l_atom];
            if ((i_type == i_tref && j_type == j_tref && k_type == k_tref && l_type == l_tref) ||
                (i_type == l_tref && j_type == k_tref && k_type == j_tref && l_type == i_tref)) {
              const bool forward = (i_type == i_tref);
              const char4 i_name = (forward) ? cdk.atom_names[i_atom] : cdk.atom_names[l_atom];
              const char4 j_name = (forward) ? cdk.atom_names[j_atom] : cdk.atom_names[k_atom];
              const char4 k_name = (forward) ? cdk.atom_names[k_atom] : cdk.atom_names[j_atom];
              const char4 l_name = (forward) ? cdk.atom_names[l_atom] : cdk.atom_names[i_atom];
              const int param_idx = vk.cimp_param_idx[pos];
              if (nfound == 0) {
                sprintf(buffer_data, "       Atom Names            Indices         Stiffness "
                        "   Phase  \n");
                di_body.push_back(buffer_data);
                sprintf(buffer_data, "  ---- ---- ---- ----  ---- ---- ---- ----  ---------- "
                        "----------\n");
                di_body.push_back(buffer_data);
              }
              const int i_index = (forward) ? i_atom + 1 : l_atom + 1;
              const int j_index = (forward) ? j_atom + 1 : k_atom + 1;
              const int k_index = (forward) ? k_atom + 1 : j_atom + 1;
              const int l_index = (forward) ? l_atom + 1 : i_atom + 1;
              sprintf(buffer_data, "  %c%c%c%c %c%c%c%c %c%c%c%c %c%c%c%c  %4d %4d %4d %4d  "
                      "%10.4lf %10.4lf\n", i_name.x, i_name.y, i_name.z, i_name.w, j_name.x,
                      j_name.y, j_name.z, j_name.w, k_name.x, k_name.y, k_name.z, k_name.w,
                      l_name.x, l_name.y, l_name.z, l_name.w, i_index, j_index, k_index, l_index,
                      vk.cimp_keq[param_idx], vk.cimp_phi[param_idx] * 180.0 / symbols::pi);
              di_body.push_back(buffer_data);
              nfound++;
            }
          }
        }
        break;
      case StateVariable::CMAP:
        {
          const char4 i_tref = atom_types[0];
          const char4 j_tref = atom_types[1];
          const char4 k_tref = atom_types[2];
          const char4 l_tref = atom_types[3];
          const char4 m_tref = atom_types[4];
          for (int pos = 0; pos < vk.ndihe; pos++) {
            const int i_atom = vk.cmap_i_atoms[pos];
            const int j_atom = vk.cmap_j_atoms[pos];
            const int k_atom = vk.cmap_k_atoms[pos];
            const int l_atom = vk.cmap_l_atoms[pos];
            const int m_atom = vk.cmap_m_atoms[pos];
            const char4 i_type = cdk.atom_types[i_atom];
            const char4 j_type = cdk.atom_types[j_atom];
            const char4 k_type = cdk.atom_types[k_atom];
            const char4 l_type = cdk.atom_types[l_atom];
            const char4 m_type = cdk.atom_types[m_atom];
            if ((i_type == i_tref && j_type == j_tref && k_type == k_tref && l_type == l_tref &&
                 m_type == m_tref) ||
                (i_type == m_tref && j_type == l_tref && k_type == k_tref && l_type == j_tref &&
                 m_type == i_tref)) {
              const bool forward = (i_type == i_tref);
              const char4 i_name = (forward) ? cdk.atom_names[i_atom] : cdk.atom_names[m_atom];
              const char4 j_name = (forward) ? cdk.atom_names[j_atom] : cdk.atom_names[l_atom];
              const char4 k_name = (forward) ? cdk.atom_names[k_atom] : cdk.atom_names[k_atom];
              const char4 l_name = (forward) ? cdk.atom_names[l_atom] : cdk.atom_names[j_atom];
              const char4 m_name = (forward) ? cdk.atom_names[m_atom] : cdk.atom_names[i_atom];
              const int param_idx = vk.cimp_param_idx[pos];
              if (nfound == 0) {
                sprintf(buffer_data, "          Atom Names                 Indices          "
                        " Map\n");
                di_body.push_back(buffer_data);
                sprintf(buffer_data, "  ---- ---- ---- ---- ----  ---- ---- ---- ---- ----  "
                        "----\n");
                di_body.push_back(buffer_data);
              }
              const int i_index = (forward) ? i_atom + 1 : m_atom + 1;
              const int j_index = (forward) ? j_atom + 1 : l_atom + 1;
              const int k_index = (forward) ? k_atom + 1 : k_atom + 1;
              const int l_index = (forward) ? l_atom + 1 : j_atom + 1;
              const int m_index = (forward) ? m_atom + 1 : i_atom + 1;
              sprintf(buffer_data, "  %c%c%c%c %c%c%c%c %c%c%c%c %c%c%c%c %c%c%c%c  %4d %4d %4d "
                      "%4d %4d  %4d\n", i_name.x, i_name.y, i_name.z, i_name.w, j_name.x, j_name.y,
                      j_name.z, j_name.w, k_name.x, k_name.y, k_name.z, k_name.w, l_name.x,
                      l_name.y, l_name.z, l_name.w, m_name.x, m_name.y, m_name.z, m_name.w,
                      i_index, j_index, k_index, l_index, m_index, vk.cmap_surf_idx[param_idx]);
              di_body.push_back(buffer_data);
              nfound++;
            }
          }
        }
        break;
      case StateVariable::RESTRAINT:
      case StateVariable::VDW:
      case StateVariable::VDW_ONE_FOUR:
      case StateVariable::ELECTROSTATIC:
      case StateVariable::ELECTROSTATIC_ONE_FOUR:
      case StateVariable::GENERALIZED_BORN:
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
        break;
      }
    }
    break;
  case DataRequestKind::STRING:

    // The text is the user-supplied message
    di_body.push_back(ask.getMessage());
    break;
  }
  data_items.emplace_back(ask, di_body);
  compareExternalRegistryNumbers(data_items.back().getExternalRegistryNumber());
}

//-------------------------------------------------------------------------------------------------
void MdlMol::writeMdl(std::ofstream *foutp, const MdlMolVersion vformat) const {
  const TextFile result(writeMdl(vformat), TextOrigin::RAM);
  writeFrame(foutp, result);
}

//-------------------------------------------------------------------------------------------------
void MdlMol::writeMdl(const std::string &fname, const MdlMolVersion vformat,
                      const PrintSituation expectation) const {
  const std::string activity = (data_item_count > 0) ?
    "Open an output file for writing an MDL MOL format structure." :
    "Open an SDF archive for writing MDL MOL format output with additional data items.";
  std::ofstream foutp = openOutputFile(fname, expectation, activity);
  writeMdl(&foutp, vformat);
  foutp.close();
}

//-------------------------------------------------------------------------------------------------
std::string MdlMol::writeMdl(const MdlMolVersion vformat) const {

  // Build the result based on the MDL MOL leading three lines, which are common to both V2000 and
  // V3000 formats.
  std::string result(title + '\n' + software_details + '\n' + general_comment + '\n');
  std::string buffer(512, ' ');
  char* buffer_data = buffer.data();
  switch (vformat) {
  case MdlMolVersion::V2000:

    // Write out the counts line.
    sprintf(buffer.data(), "%3d%3d%3d   %3d%3d            999 V2000\n", atom_count, bond_count,
            list_count, static_cast<int>(chirality), stext_entry_count);
    result.append(buffer_data, 40);
    
    // Write out the atom block.
    for (int i = 0; i < atom_count; i++) {
      sprintf(buffer_data, "%10.4lf%10.4lf%10.4lf %c%c%c%2d%3d%3d%3d%3d%3d%3d  0  0%3d%3d%3d\n",
              coordinates[i].x, coordinates[i].y, coordinates[i].z, atomic_symbols[i].x,
              atomic_symbols[i].y, atomic_symbols[i].z, getIsotopicShiftCode(i),
              getFormalChargeCode(i), static_cast<int>(parities[i]), getImplicitHydrogenCode(i),
              static_cast<int>(stereo_considerations[i]), static_cast<int>(valence_connections[i]),
              static_cast<int>(hydrogenation_protocol[i]), atom_atom_mapping_count[i],
              static_cast<int>(orientation_stability[i]),
              static_cast<int>(exact_change_enforced[i]));
      result.append(buffer_data, 70);
    }

    // Write out the bond block.
    for (int i = 0; i < bond_count; i++) {

      // Add 1 to the bond atom indices to get back into the file format.  This adjustment is
      // automated for the properties.
      sprintf(buffer.data(), "%3d%3d%3d%3d  0%3d%3d\n", bonds[i].getFirstAtom() + 1,
              bonds[i].getSecondAtom() + 1, static_cast<int>(bonds[i].getOrder()),
              static_cast<int>(bonds[i].getStereochemistry()),
              static_cast<int>(bonds[i].getRingStatus()),
              static_cast<int>(bonds[i].getReactivity()));
      result.append(buffer_data, 22);
    }

    // Write out the atom list block, if appropriate.
    if (property_element_lists == false) {
      for (int i = 0; i < list_count; i++) {
        const int n_entry = element_lists[i].getEntryCount();
        sprintf(buffer_data, "%3d %c    %d", element_lists[i].getAttachmentPoint(),
                element_lists[i].getExclusionCode(), n_entry);
        int nchar = 10;
        for (int j = 0; j < n_entry; j++) {
          sprintf(&buffer_data[nchar], " %3d", element_lists[i].getEntry(j));
          nchar += 4;
        }
        buffer_data[nchar] = '\n';
        nchar++;
        buffer_data[nchar] = '\0';
        result.append(buffer_data, nchar);
      }
    }

    // Write out the S-text block

    // Write out the properties block.
    for (int i = 0; i < properties_count; i++) {
      result.append(properties[i].getMdlText());
    }

    // Write out the terminal line (not considered a property, even in the V2000 format)
    result.append("M  END\n");
    
    break;
  case MdlMolVersion::V3000:
  case MdlMolVersion::UNKNOWN:

    // Make the default writing method the V3000 format.
    break;
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
void MdlMol::writeDataItems(std::ofstream *foutp, const int mol_index) const {
  const TextFile result(writeDataItems(mol_index), TextOrigin::RAM);
  writeFrame(foutp, result);
}

//-------------------------------------------------------------------------------------------------
void MdlMol::writeDataItems(const std::string &fname, const PrintSituation expectation,
                            const int mol_index) const {
  const std::string activity("Open an SDF archive for writing MDL MOL format output with "
                             "additional data items.");
  std::ofstream foutp = openOutputFile(fname, expectation, activity);
  writeDataItems(&foutp, mol_index);
  foutp.close();
}

//-------------------------------------------------------------------------------------------------
std::string MdlMol::writeDataItems(const int mol_index) const {
  for (int i = 0; i < data_item_count; i++) {

    // Compose the header line
    std::string header_line("> ");
    bool required_id = false;
    if (data_items[i].placeInternalRegnoInHeader()) {
      header_line += std::to_string(mol_index);
    }
    if (data_items[i].placeExternalRegnoInHeader()) {
      header_line += "(" + data_items[i].getExternalRegistryNumber() + ")";
    }
    if (data_items[i].placeTitleInHeader()) {
      header_line += "<" + data_items[i].getItemName() + ">";
      required_id = true;
    }
    if (data_items[i].placeMaccsIIFieldInHeader()) {
      header_line += "DT" + std::to_string(data_items[i].getMaccsFieldNumber());
      required_id = true;
    }
    if (data_items[i].noteArchivesInHeader()) {
      header_line += " FROM ARCHIVES";
    }

    // CHECK
    printf("Header = %s\n", header_line.c_str()); 
    // END CHECK
  }
}

//-------------------------------------------------------------------------------------------------
void MdlMol::allocate() {

  // Atom property fields are resized and then set as part of a loop in the parent MdlMol
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
int MdlMol::getIsotopicShiftCode(const int atom_index) const {
  if (property_isotopes) {
    return 0;
  }
  else {
    return isotopic_shifts[atom_index];
  }
  __builtin_unreachable();
}
 
//-------------------------------------------------------------------------------------------------
void MdlMol::interpretFormalCharge(const int charge_in, const int atom_index) {
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
          "entry.  Title of entry: \"" + title + "\".", "MdlMol", "interpretFormalCharge");    
  }
}

//-------------------------------------------------------------------------------------------------
int MdlMol::getFormalChargeCode(const int atom_index) const {
  if (property_formal_charges) {
    return 0;
  }
  switch (formal_charges[atom_index]) {
  case 0:
    return (radicals[atom_index] == RadicalState::DOUBLET) ? 4 : 0;
  case 1:
    return 3;
  case 2:
    return 2;
  case 3:
    return 1;
  case -1:
    return 5;
  case -2:
    return 6;
  case -3:
    return 7;
  default:

    // Formal charges outside the V2000 range must be handled by properties, and any such case
    // implies that all formal charges are handled this way.
    rtErr("A formal charge that cannot be expressed in the atoms block of a V2000 format MDL MOL "
          "file exists, but there are no properties to account for it.", "MdlMol",
          "getFormalChargeCode");
  }
  __builtin_unreachable();
}
 
//-------------------------------------------------------------------------------------------------
MolObjAtomStereo MdlMol::interpretStereoParity(const int setting_in) {
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
          "Title of entry: \"" + title + "\".", "MdlMol", "interpretStereoParity");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
void MdlMol::interpretImplicitHydrogenContent(const int nh_in, const int atom_index) {
  if (nh_in > 5 || nh_in < 0) {
    rtErr("An implicit hydrogen content of " + std::to_string(nh_in) + " would imply " +
          std::to_string(nh_in - 1) + " hydrogens can be inferred around an atom in MDL MOL "
          "entry \"" + title + "\".", "MdlMol", "interpretImplicitHydrogenContent");
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
int MdlMol::getImplicitHydrogenCode(const int atom_index) const {
  switch (hydrogenation_protocol[atom_index]) {
  case HydrogenAssignment::VALENCE_SHELL:
    return (implicit_hydrogens[atom_index] == 0) ? 0 : implicit_hydrogens[atom_index] + 1;
  case HydrogenAssignment::DO_NOT_HYDROGENATE:
    return 1;
  }
  __builtin_unreachable();
}
 
//-------------------------------------------------------------------------------------------------
bool MdlMol::interpretBooleanValue(const int value_in, const std::string &desc) {
  if (value_in != 0 && value_in != 1) {
    rtErr("A directive of " + std::to_string(value_in) + " is invalid when " + desc + ".  Title "
          "of entry: \"" + title + "\".", "MdlMol", "interpretBooleanValue");
  }
  return (value_in == 1);
}

//-------------------------------------------------------------------------------------------------
int MdlMol::interpretValenceNumber(const int count_in) {
  if (count_in == 15) {
    return 0;
  }
  else if (count_in >= 0 && count_in < 15) {
    return count_in;
  }
  else {
    rtErr("An atom cannot have " + std::to_string(count_in) + " valence connections, as is the "
          "case for one atom in MDL MOL entry \"" + title + "\".", "MdlMol",
          "interpretValenceNumber");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
StereoRetention MdlMol::interpretStereoStability(const int code_in) {
  switch (code_in) {
  case 0:
    return StereoRetention::NOT_APPLIED;
  case 1:
    return StereoRetention::INVERTED;
  case 2:
    return StereoRetention::RETAINED;
  default:
    rtErr("A stereochemistry retention code of " + std::to_string(code_in) + " is invalid in MDL "
          "MOL entry \"" + title + "\".", "MdlMol", "interpretStereoStability");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
void MdlMol::updateV2kAtomAttributes() {

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
      property_formal_charges = true;
    }
    if (properties[i].getCode() == char4({ 'R', 'A', 'D', 'M' })) {
      for (int j = 0; j < atom_count; j++) {
        formal_charges[j] = 0;
        radicals[j] = RadicalState::NONE;
      }
      property_radicals = true;
    }
    if (properties[i].getCode() == char4({ 'I', 'S', 'O', 'M' })) {
      for (int j = 0; j < atom_count; j++) {
        isotopic_shifts[j] = 0;
      }
      property_isotopes = true;
    }
    if (properties[i].getCode() == char4({ 'A', 'L', 'S', 'M' })) {
      property_element_lists = true;
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
void MdlMol::hydrogenate() {

}

//-------------------------------------------------------------------------------------------------
void MdlMol::compareExternalRegistryNumbers(const std::string &regno_in) {
  if (regno_in.size() > 0LLU) {
    if (external_regno.size() == 0LLU) {
      external_regno = regno_in;
    }
    else {
      switch (policy) {
      case ExceptionResponse::DIE:
        rtErr("The external registry number associated with a data item to add (" + regno_in +
              ") is inconsistent with the external registry number already associated with the "
              "MDL MOL entry (" + external_regno + ").", "MdlMol",
              "compareExternalRegistryNumbers");
      case ExceptionResponse::WARN:
        rtWarn("The external registry number associated with a data item to add (" + regno_in +
               ") is inconsistent with the external registry number already associated with the "
               "MDL MOL entry (" + external_regno + ").  The existing registry number will take "
               "precedence.", "MdlMol", "compareExternalRegistryNumbers");
        break;
      case ExceptionResponse::SILENT:
        break;
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------
void MdlMol::checkAtomCount(int ext_atom_count) {
  if (ext_atom_count != atom_count) {
    rtErr("The number of atoms coming in (" + std::to_string(ext_atom_count) + ") is not "
          "consistent with the number of atoms in the molecule (" + std::to_string(atom_count) +
          ").", "MdlMol", "checkAtomCount");
  }
}

//-------------------------------------------------------------------------------------------------
std::vector<MolObjBond> operator+(const std::vector<MolObjBond> &lhs,
				  const std::vector<MolObjBond> &rhs) {
  std::vector<MolObjBond> result(lhs.begin(), lhs.end());
  result.insert(result.end(), rhs.begin(), rhs.end());
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<MdlMol> readStructureDataFile(const TextFile &tf, const int low_frame_limit,
                                          const int high_frame_limit,
                                          const CaseSensitivity capitalization,
                                          const ExceptionResponse policy) {
  std::vector<MdlMol> result;

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
std::vector<MdlMol> readStructureDataFile(const TextFile &tf, const CaseSensitivity capitalization,
                                          const ExceptionResponse policy) {
  std::vector<MdlMol> result;

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
std::vector<MdlMol> readStructureDataFile(const std::string &file_name,
                                          const CaseSensitivity capitalization,
                                          const ExceptionResponse policy) {
  const TextFile tf(file_name);
  return readStructureDataFile(tf, capitalization, policy);
}

} // namespace structure
} // namespace stormm
