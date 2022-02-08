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
  time_based_restraints{false},
  rposn_atoms{HybridKind::POINTER, "rst_position_i"},
  rbond_i_atoms{HybridKind::POINTER, "rst_bond_i"},
  rbond_j_atoms{HybridKind::POINTER, "rst_bond_i"},
  rangl_i_atoms{HybridKind::POINTER, "rst_angl_i"},
  rangl_j_atoms{HybridKind::POINTER, "rst_angl_j"},
  rangl_k_atoms{HybridKind::POINTER, "rst_angl_k"},
  rdihe_i_atoms{HybridKind::POINTER, "rst_dihe_i"},
  rdihe_j_atoms{HybridKind::POINTER, "rst_dihe_j"},
  rdihe_k_atoms{HybridKind::POINTER, "rst_dihe_k"},
  rdihe_l_atoms{HybridKind::POINTER, "rst_dihe_l"},
  int_data{HybridKind::ARRAY, "rst_int_data"},
  rposn_init_keq{HybridKind::POINTER, "rposn_init_keq"},
  rposn_final_keq{HybridKind::POINTER, "rposn_final_keq"},
  rposn_init_r{HybridKind::POINTER, "rposn_init_r"},
  rposn_final_r{HybridKind::POINTER, "rposn_final_r"},
  rposn_init_xy{HybridKind::POINTER, "rposn_init_xy"},
  rposn_init_z{HybridKind::POINTER, "rposn_init_z"},
  rposn_final_xy{HybridKind::POINTER, "rposn_final_xy"},
  rposn_final_z{HybridKind::POINTER, "rposn_final_z"},
  rbond_init_keq{HybridKind::POINTER, "rbond_init_keq"},
  rbond_final_keq{HybridKind::POINTER, "rbond_final_keq"},
  rbond_init_r{HybridKind::POINTER, "rbond_init_r"},
  rbond_final_r{HybridKind::POINTER, "rbond_final_r"},
  rangl_init_keq{HybridKind::POINTER, "rangl_init_keq"},
  rangl_final_keq{HybridKind::POINTER, "rangl_final_keq"},
  rangl_init_r{HybridKind::POINTER, "rangl_init_r"},
  rangl_final_r{HybridKind::POINTER, "rangl_final_r"},
  rdihe_init_keq{HybridKind::POINTER, "rdihe_init_keq"},
  rdihe_final_keq{HybridKind::POINTER, "rdihe_final_keq"},
  rdihe_init_r{HybridKind::POINTER, "rdihe_init_r"},
  rdihe_final_r{HybridKind::POINTER, "rdihe_final_r"},
  double_data{HybridKind::ARRAY, "rst_double_data"},
  double2_data{HybridKind::ARRAY, "rst_double2_data"},
  double4_data{HybridKind::ARRAY, "rst_double4_data"},
  sp_rposn_init_keq{HybridKind::POINTER, "sp_rposn_init_keq"},
  sp_rposn_final_keq{HybridKind::POINTER, "sp_rposn_final_keq"},
  sp_rposn_init_r{HybridKind::POINTER, "sp_rposn_init_r"},
  sp_rposn_final_r{HybridKind::POINTER, "sp_rposn_final_r"},
  sp_rposn_init_xy{HybridKind::POINTER, "sp_rposn_init_xy"},
  sp_rposn_init_z{HybridKind::POINTER, "sp_rposn_init_z"},
  sp_rposn_final_xy{HybridKind::POINTER, "sp_rposn_final_xy"},
  sp_rposn_final_z{HybridKind::POINTER, "sp_rposn_final_z"},
  sp_rbond_init_keq{HybridKind::POINTER, "sp_rbond_init_keq"},
  sp_rbond_final_keq{HybridKind::POINTER, "sp_rbond_final_keq"},
  sp_rbond_init_r{HybridKind::POINTER, "sp_rbond_init_r"},
  sp_rbond_final_r{HybridKind::POINTER, "sp_rbond_final_r"},
  sp_rangl_init_keq{HybridKind::POINTER, "sp_rangl_init_keq"},
  sp_rangl_final_keq{HybridKind::POINTER, "sp_rangl_final_keq"},
  sp_rangl_init_r{HybridKind::POINTER, "sp_rangl_init_r"},
  sp_rangl_final_r{HybridKind::POINTER, "sp_rangl_final_r"},
  sp_rdihe_init_keq{HybridKind::POINTER, "sp_rdihe_init_keq"},
  sp_rdihe_final_keq{HybridKind::POINTER, "sp_rdihe_final_keq"},
  sp_rdihe_init_r{HybridKind::POINTER, "sp_rdihe_init_r"},
  sp_rdihe_final_r{HybridKind::POINTER, "sp_rdihe_final_r"},
  float_data{HybridKind::ARRAY, "rst_float_data"},
  float2_data{HybridKind::ARRAY, "rst_float2_data"},
  float4_data{HybridKind::ARRAY, "rst_float4_data"}
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
  std::vector<int> tmp_rposn_init_step(position_count);
  std::vector<int> tmp_rposn_final_step(position_count);
  std::vector<int> tmp_rbond_init_step(distance_count);
  std::vector<int> tmp_rbond_final_step(distance_count);
  std::vector<int> tmp_rangl_init_step(angle_count);
  std::vector<int> tmp_rangl_final_step(angle_count);
  std::vector<int> tmp_rdihe_init_step(dihedral_count);
  std::vector<int> tmp_rdihe_final_step(dihedral_count);
  std::vector<double> tmp_rposn_init_z(position_count);
  std::vector<double> tmp_rposn_final_z(position_count);
  std::vector<double2> tmp_rposn_init_keq(position_count);
  std::vector<double2> tmp_rposn_final_keq(position_count);
  std::vector<double2> tmp_rposn_init_xy(position_count);
  std::vector<double2> tmp_rposn_final_xy(position_count);
  std::vector<double2> tmp_rbond_init_keq(distance_count);
  std::vector<double2> tmp_rbond_final_keq(distance_count);
  std::vector<double2> tmp_rangl_init_keq(angle_count);
  std::vector<double2> tmp_rangl_final_keq(angle_count);
  std::vector<double2> tmp_rdihe_init_keq(dihedral_count);
  std::vector<double2> tmp_rdihe_final_keq(dihedral_count);
  int nposnr = 0;
  int nbondr = 0;
  int nanglr = 0;
  int ndiher = 0;
  for (int i = 0; i < restraint_count; i++) {
    const int rord = rbasis[i].getOrder();
    if (rord == 1) {
      tmp_rposn_i[nposnr]          = rbasis[i].getAtomIndex(1);
      tmp_rposn_init_step[nposnr]  = rbasis[i].getInitialStep();
      tmp_rposn_final_step[nposnr] = rbasis[i].getFinalStep();
      tmp_rposn_init_keq[nposnr]   = rbasis[i].getInitialStiffness();
      tmp_rposn_final_keq[nposnr]  = rbasis[i].getFinalStiffness();
      nposnr++;
    }
    else if (rord == 2) {
      tmp_rbond_i[nbondr]          = rbasis[i].getAtomIndex(1);
      tmp_rbond_j[nbondr]          = rbasis[i].getAtomIndex(2);
      tmp_rbond_init_step[nanglr]  = rbasis[i].getInitialStep();
      tmp_rbond_final_step[nanglr] = rbasis[i].getFinalStep();
      tmp_rbond_init_keq[nbondr]   = rbasis[i].getInitialStiffness();
      tmp_rbond_final_keq[nbondr]  = rbasis[i].getFinalStiffness();
      nbondr++;
    }
    else if (rord == 3) {
      tmp_rangl_i[nanglr]          = rbasis[i].getAtomIndex(1);
      tmp_rangl_j[nanglr]          = rbasis[i].getAtomIndex(2);
      tmp_rangl_k[nanglr]          = rbasis[i].getAtomIndex(3);
      tmp_rangl_init_step[nanglr]  = rbasis[i].getInitialStep();
      tmp_rangl_final_step[nanglr] = rbasis[i].getFinalStep();
      tmp_rangl_init_keq[nanglr]   = rbasis[i].getInitialStiffness();
      tmp_rangl_final_keq[nanglr]  = rbasis[i].getFinalStiffness();
      nanglr++;
    }
    else if (rord == 4) {
      tmp_rdihe_i[ndiher]          = rbasis[i].getAtomIndex(1);
      tmp_rdihe_j[ndiher]          = rbasis[i].getAtomIndex(2);
      tmp_rdihe_k[ndiher]          = rbasis[i].getAtomIndex(3);
      tmp_rdihe_l[ndiher]          = rbasis[i].getAtomIndex(4);
      tmp_rdihe_init_step[ndiher]  = rbasis[i].getInitialStep();
      tmp_rdihe_final_step[ndiher] = rbasis[i].getFinalStep();
      tmp_rdihe_init_keq[ndiher]   = rbasis[i].getInitialStiffness();
      tmp_rdihe_final_keq[ndiher]  = rbasis[i].getFinalStiffness();
      ndiher++;
    }
  }
  
  // Allocate the necessary space and set pointers
  const int padded_posn_count = roundUp(position_count, warp_size_int);
  const int padded_bond_count = roundUp(distance_count, warp_size_int);
  const int padded_angl_count = roundUp(angle_count, warp_size_int);
  const int padded_dihe_count = roundUp(dihedral_count, warp_size_int);
  const int nint_elem = (3 * padded_posn_count) + (4 * padded_bond_count) +
                        (5 * padded_angl_count) + (6 * padded_dihe_count);
  int_data.resize(nint_elem);
  size_t ic = 0LLU;
  ic = rposn_atoms.putHost(&int_data, tmp_rposn_i, ic, warp_size_zu);
  ic = rbond_i_atoms.putHost(&int_data, tmp_rbond_i, ic, warp_size_zu);
  ic = rbond_j_atoms.putHost(&int_data, tmp_rbond_j, ic, warp_size_zu);
  ic = rangl_i_atoms.putHost(&int_data, tmp_rangl_i, ic, warp_size_zu);
  ic = rangl_j_atoms.putHost(&int_data, tmp_rangl_j, ic, warp_size_zu);
  ic = rangl_k_atoms.putHost(&int_data, tmp_rangl_k, ic, warp_size_zu);
  ic = rdihe_i_atoms.putHost(&int_data, tmp_rdihe_i, ic, warp_size_zu);
  ic = rdihe_j_atoms.putHost(&int_data, tmp_rdihe_j, ic, warp_size_zu);
  ic = rdihe_k_atoms.putHost(&int_data, tmp_rdihe_k, ic, warp_size_zu);
  ic = rdihe_l_atoms.putHost(&int_data, tmp_rdihe_l, ic, warp_size_zu);
  ic = rposn_init_step.putHost(&int_data, tmp_rposn_init_step, ic, warp_size_zu);
  ic = rposn_final_step.putHost(&int_data, tmp_rposn_final_step, ic, warp_size_zu);
  ic = rbond_init_step.putHost(&int_data, tmp_rbond_init_step, ic, warp_size_zu);
  ic = rbond_final_step.putHost(&int_data, tmp_rbond_final_step, ic, warp_size_zu);
  ic = rangl_init_step.putHost(&int_data, tmp_rangl_init_step, ic, warp_size_zu);
  ic = rangl_final_step.putHost(&int_data, tmp_rangl_final_step, ic, warp_size_zu);
  ic = rdihe_init_step.putHost(&int_data, tmp_rdihe_init_step, ic, warp_size_zu);
  ic = rdihe_final_step.putHost(&int_data, tmp_rdihe_final_step, ic, warp_size_zu);
  const int ndbl_elem = 2 * padded_posn_count;
  double_data.resize(ndbl_elem);
  size_t dc = 0LLU;
  dc = rposn_init_z.putHost(&double_data, tmp_rposn_init_z, ic, warp_size_zu);
  dc = rposn_final_z.putHost(&double_data, tmp_rposn_final_z, ic, warp_size_zu);
  const int ndbl2_elem = (4 * padded_posn_count) + (2 * padded_bond_count) +
                         (2 * padded_angl_count) + (4 * padded_dihe_count);
  size_t d2c = 0LLU;
  d2c = rposn_init_keq.putHost(&double2_data, tmp_rposn_init_keq, d2c, warp_size_zu);
  d2c = rposn_final_keq.putHost(&double2_data, tmp_rposn_final_keq, d2c, warp_size_zu);
  d2c = rposn_init_xy.putHost(&double2_data, tmp_rposn_init_xy, d2c, warp_size_zu);
  d2c = rposn_final_xy.putHost(&double2_data, tmp_rposn_final_xy, d2c, warp_size_zu);
  d2c = rbond_init_keq.putHost(&double2_data, tmp_rbond_init_keq, d2c, warp_size_zu);
  d2c = rbond_final_keq.putHost(&double2_data, tmp_rbond_final_keq, d2c, warp_size_zu);
  d2c = rangl_init_keq.putHost(&double2_data, tmp_rangl_init_keq, d2c, warp_size_zu);
  d2c = rangl_final_keq.putHost(&double2_data, tmp_rangl_final_keq, d2c, warp_size_zu);
  d2c = rdihe_init_keq.putHost(&double2_data, tmp_rdihe_init_keq, d2c, warp_size_zu);
  d2c = rdihe_final_keq.putHost(&double2_data, tmp_rdihe_final_keq, d2c, warp_size_zu);
}

} // namespace restraints
} // namespace omni
