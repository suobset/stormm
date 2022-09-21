#include "copyright.h"
#include "mdl_mol_format.h"

namespace stormm {
namespace structure {

using parse::readIntegerValue;
using parse::readRealValue;
using parse::TextFileReader;

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
MdlMolObj::MdlMolObj(const TextFile &tf, const int line_start, const int line_end_in):
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
    if (tfr.line_limits[next_delim_loc + 1] - tfr.line_limits[next_delim_loc] >= 4 &&
        lptr[0] == '$' && lptr[1] == '$' && lptr[2] == '$' && lptr[3] == '$') {
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
    const char* counts_line_ptr = &tfr.text[tfr.line_limits[line_start + 3]];
    atom_count   = readIntegerValue(counts_line_ptr, 0, 3);
    bond_count   = readIntegerValue(counts_line_ptr, 3, 3);
    list_count = readIntegerValue(counts_line_ptr, 6, 3);
    const int chir_num = readIntegerValue(counts_line_ptr, 12, 3);
    if (chir_num == 0) {
      chirality = MolObjChirality::ACHIRAL;
    }
    else if (chir_num == 1) {
      chirality = MolObjChirality::CHIRAL;
    }
    else {
      rtErr("Invalid chirality setting detected at line " + std::to_string(line_start + 3) +
            " in .sdf or MDL MOL file " + tf.getFileName() + ".", "MdlMolObj");
    }
    stext_entry_count = readIntegerValue(counts_line_ptr, 15, 3);
    properties_count  = readIntegerValue(counts_line_ptr, 30, 3);

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
      atomic_symbols[iatm] = tf.extractChar4(i, 31, 3);
      isotopic_shifts[iatm] = readIntegerValue(atom_line_ptr, 34, 2);
      if (isotopic_shifts[iatm] > 4 || isotopic_shifts[iatm] < -3) {
        rtErr("A V2000 MOL format entry should not describe an isotopic shift outside the range "
              "[-3, 4].  Shift found: " + std::to_string(isotopic_shifts[iatm]) +
              ".  Title of entry: \"" + title + "\".", "MdlMolObj");
      }
      
      // Standard Template Library vector<bool> works differently from other vectors.  Set its
      // contents in a different manner.
      bool dblt_flag;
      formal_charges[iatm] = interpretFormalCharge(readIntegerValue(atom_line_ptr, 36, 3),
                                                   &dblt_flag);
      doublet_radicals[iatm] = dblt_flag;
      parities[iatm] = interpretStereoParity(readIntegerValue(atom_line_ptr, 39, 3));
      implicit_hydrogens[iatm] =
        interpretImplicitHydrogenContent(readIntegerValue(atom_line_ptr, 42, 3));
      stereo_considerations[iatm] =
        interpretBooleanValue(readIntegerValue(atom_line_ptr, 45, 3),
                              "interpreting stereochemical considerations");
      valence_connections[iatm] = interpretValenceNumber(readIntegerValue(atom_line_ptr, 48, 3));
      if (readIntegerValue(atom_line_ptr, 51, 3) == 1 && implicit_hydrogens[iatm] > 0) {
        rtErr("The H0 designation, indicating that implicit hydrogens are not allowed on atom " +
              std::to_string(iatm + 1) + " of MDL MOL entry \"" +  title + "\", is present but "
              "the number of implicit hydrogens has also been indicated as " +
              std::to_string(implicit_hydrogens[iatm]) + ".", "MdlMolObj");
      }
      atom_atom_mapping_count[iatm] = readIntegerValue(atom_line_ptr, 60, 3);
      orientation_stability[iatm] =
        interpretStereoStability(readIntegerValue(atom_line_ptr, 63, 3));
      exact_change_enforced[iatm] = interpretBooleanValue(readIntegerValue(atom_line_ptr, 66, 3),
                                                          "interpreting exact change flag");
    }
  }

  // Read the bonds block
  if (version == 2000) {
    
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
void MdlMolObj::allocate() {
  coordinates.resize(atom_count);
  atomic_symbols.resize(atom_count);
  atomic_numbers.resize(atom_count);
  formal_charges.resize(atom_count);
  isotopic_shifts.resize(atom_count);
  parities.resize(atom_count);
  implicit_hydrogens.resize(atom_count);
  stereo_considerations.resize(atom_count);
  valence_connections.resize(atom_count);
  atom_atom_mapping_count.resize(atom_count);
  exact_change_enforced.resize(atom_count);
  orientation_stability.resize(atom_count);
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
  if (nh_in >= 5 || nh_in < 1) {
    rtErr("An implicit hydrogen content of " + std::to_string(nh_in) + " would imply " +
          std::to_string(nh_in - 1) + " hydrogens can be inferred around an atom in MDL MOL "
          "entry " + title + "\".", "MdlMolObj", "interpretImplicitHydrogenContent");
  }
  return nh_in - 1;
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
}

//-------------------------------------------------------------------------------------------------
StereoChemicalRetention MdlMolObj::interpretStereoStability(const int code_in) {
  switch (code_in) {
  case 0:
    return StereoChemicalRetention::NOT_APPLIED;
  case 1:
    return StereoChemicalRetention::INVERTED;
  case 2:
    return StereoChemicalRetention::RETAINED;
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
  
} // namespace structure
} // namespace stormm
