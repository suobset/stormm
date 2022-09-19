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
    atomic_numbers{}, formal_charges{}, isotopic_shifts{}, bonds{}, properties{}, title{""},
    software_details{""}, general_comment{""}
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
  }
    
  // Read the atoms block

  // Read the bonds block

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
