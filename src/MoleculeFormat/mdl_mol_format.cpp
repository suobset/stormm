#include "copyright.h"
#include "mdl_mol_format.h"

namespace stormm {
namespace structure {

//-------------------------------------------------------------------------------------------------
MdlMolObj::MdlMolObj():
  atom_count{0},
  bond_count{0},
  sgroup_count{0},
  constraint_count{0},
  chirality{MolObjChirality::ACHIRAL},
  registry_number{-1},
  coordinates{},
  atomic_symbols{},
  atomic_numbers{},
  formal_charges{},
  isotopic_shifts{},
  bonds{},
  properties{},
  title{""},
  software_details{""},
  general_comment{""}
{}

//-------------------------------------------------------------------------------------------------
MdlMolObj::MdlMolObj(const TextFile &tf, const int line_start, const int line_end_in):
  MdlMolObj()
{
  // Default line end of -1 indicates reading to the end of the file
  const int line_end = (line_end_in < 0) ? tf.getLineCount() : line_end_in;
  
  // Begin by reading the molecule name (title), generating software details, and any general
  // comment (always three and only three distinct lines, even if left blank)

  // Read the counts line
  
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

} // namespace structure
} // namespace stormm
