#include "Math/rounding.h"
#include "restraint_apparatus.h"

namespace omni {
namespace restraints {

using card::HybridKind;
using math::roundUp;
  
//-------------------------------------------------------------------------------------------------
RestraintApparatus::RestraintApparatus(const std::vector<BoundedRestraint> &rbasis) :
  restraint_count{static_cast<int>(rbasis.size())},
  position_count{0}, distance_count{0}, angle_count{0}, dihedral_count{0},
  rposn_i_atoms{HybridKind::POINTER, "rst_position_i"},
  rbond_i_atoms{HybridKind::POINTER, "rst_bond_i"},
  rbond_j_atoms{HybridKind::POINTER, "rst_bond_i"},
  rangl_i_atoms{HybridKind::POINTER, "rst_angl_i"},
  rangl_j_atoms{HybridKind::POINTER, "rst_angl_j"},
  rangl_k_atoms{HybridKind::POINTER, "rst_angl_k"},
  rdihe_i_atoms{HybridKind::POINTER, "rst_dihe_i"},
  rdihe_j_atoms{HybridKind::POINTER, "rst_dihe_j"},
  rdihe_k_atoms{HybridKind::POINTER, "rst_dihe_k"},
  rdihe_l_atoms{HybridKind::POINTER, "rst_dihe_l"},
  int_data{HybridKind::ARRAY, "rst_int_data"}
{
  // Count distance, angle, and dihedral restraints
  for (int i = 0; i < restraint_count; i++) {
    const int rord = rbasis[i].getOrder();
    position_count += (rord == 1);
    distance_count += (rord == 2);
    angle_count    += (rord == 3);
    dihedral_count += (rord == 4);
  }

  // Collect distance, angle, and dihedral restraint atoms
  std::vector<int> tmp_rposn_i(position_count);
  std::vector<int> tmp_rbond_i(distance_count);
  std::vector<int> tmp_rbond_j(distance_count);
  std::vector<int> tmp_rangl_i(angle_count);
  std::vector<int> tmp_rangl_j(angle_count);
  std::vector<int> tmp_rangl_k(angle_count);
  std::vector<int> tmp_rdihe_i(dihedral_count);
  std::vector<int> tmp_rdihe_j(dihedral_count);
  std::vector<int> tmp_rdihe_k(dihedral_count);
  std::vector<int> tmp_rdihe_l(dihedral_count);
  int nposnr = 0;
  int nbondr = 0;
  int nanglr = 0;
  int ndiher = 0;
  for (int i = 0; i < restraint_count; i++) {
    const int rord = rbasis[i].getOrder();
    if (rord == 1) {
      tmp_rposn_i[nposnr] = rbasis[i].getAtomIndex(1);
      nposnr++;
    }
    else if (rord == 2) {
      tmp_rbond_i[nbondr] = rbasis[i].getAtomIndex(1);
      tmp_rbond_j[nbondr] = rbasis[i].getAtomIndex(2);
      nbondr++;
    }
    else if (rord == 3) {
      tmp_rangl_i[nanglr] = rbasis[i].getAtomIndex(1);
      tmp_rangl_j[nanglr] = rbasis[i].getAtomIndex(2);
      tmp_rangl_k[nanglr] = rbasis[i].getAtomIndex(3);
      nanglr++;
    }
    else if (rord == 4) {
      tmp_rdihe_i[ndiher] = rbasis[i].getAtomIndex(1);
      tmp_rdihe_j[ndiher] = rbasis[i].getAtomIndex(2);
      tmp_rdihe_k[ndiher] = rbasis[i].getAtomIndex(3);
      tmp_rdihe_l[ndiher] = rbasis[i].getAtomIndex(4);
      ndiher++;
    }
  }
  
  // Allocate the necessary space and set pointers
  const int nelem = (2 * roundUp(distance_count, warp_size_int)) +
                    (3 * roundUp(angle_count, warp_size_int)) +
                    (4 * roundUp(dihedral_count, warp_size_int));
  int_data.resize(nelem);
  size_t ic = 0LLU;
  ic = rbond_i_atoms.putHost(&int_data, tmp_rbond_i, ic, warp_size_zu);
  ic = rbond_j_atoms.putHost(&int_data, tmp_rbond_j, ic, warp_size_zu);
  ic = rangl_i_atoms.putHost(&int_data, tmp_rangl_i, ic, warp_size_zu);
  ic = rangl_j_atoms.putHost(&int_data, tmp_rangl_j, ic, warp_size_zu);
  ic = rangl_k_atoms.putHost(&int_data, tmp_rangl_k, ic, warp_size_zu);
  ic = rdihe_i_atoms.putHost(&int_data, tmp_rdihe_i, ic, warp_size_zu);
  ic = rdihe_j_atoms.putHost(&int_data, tmp_rdihe_j, ic, warp_size_zu);
  ic = rdihe_k_atoms.putHost(&int_data, tmp_rdihe_k, ic, warp_size_zu);
  ic = rdihe_l_atoms.putHost(&int_data, tmp_rdihe_l, ic, warp_size_zu);

  //
}

} // namespace restraints
} // namespace omni
