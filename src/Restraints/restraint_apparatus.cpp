#include "Math/rounding.h"
#include "Reporting/error_format.h"
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
    rposn_init_step{HybridKind::POINTER, "rst_posn_step0"},
    rposn_final_step{HybridKind::POINTER, "rst_posn_stepF"},
    rbond_init_step{HybridKind::POINTER, "rst_bond_step0"},
    rbond_final_step{HybridKind::POINTER, "rst_bond_stepF"},
    rangl_init_step{HybridKind::POINTER, "rst_angl_step0"},
    rangl_final_step{HybridKind::POINTER, "rst_angl_stepF"},
    rdihe_init_step{HybridKind::POINTER, "rst_dihe_step0"},
    rdihe_final_step{HybridKind::POINTER, "rst_dihe_stepF"},
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
    float4_data{HybridKind::ARRAY, "rst_float4_data"},
    ag_pointer{(restraint_count > 0) ? rbasis[0].getTopologyPointer() : nullptr}
{
  // Do not let a restraint apparatus be created with an empty vector of BoundedRestraints
  if (restraint_count == 0) {
    rtErr("A restraint apparatus must be initialized with a topology pointer.  An empty vector "
          "is therefore unacceptable, as there is no individual restraint to suggest a topology.  "
          "Instead, use the overload taking the topology pointer as its sole argument.",
          "RestraintApparatus");
  }

  // Set the topology pointer and check for consistency
  for (int i = 1; i < restraint_count; i++) {
    if (rbasis[i].getTopologyPointer() != ag_pointer) {
      rtErr("All restraints must point to the same topology.", "RestraintApparatus");
    }
  }
  
  // Count distance, angle, and dihedral restraints
  for (int i = 0; i < restraint_count; i++) {
    const int rord = rbasis[i].getOrder();
    position_count += (rord == 1);
    distance_count += (rord == 2);
    angle_count    += (rord == 3);
    dihedral_count += (rord == 4);
  }

  // Allocate Hybrid space and set pointers
  allocate();
  
  // Collect all restraints' atoms and parameters in host-side vectors for upload en masse
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
  std::vector<double4> tmp_rposn_init_r(position_count);
  std::vector<double4> tmp_rposn_final_r(position_count);
  std::vector<double4> tmp_rbond_init_r(distance_count);
  std::vector<double4> tmp_rbond_final_r(distance_count);
  std::vector<double4> tmp_rangl_init_r(angle_count);
  std::vector<double4> tmp_rangl_final_r(angle_count);
  std::vector<double4> tmp_rdihe_init_r(dihedral_count);
  std::vector<double4> tmp_rdihe_final_r(dihedral_count);
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
      double3 refcrd               = rbasis[i].getInitialTargetSite();
      tmp_rposn_init_xy[nposnr]    = {refcrd.x, refcrd.y};
      tmp_rposn_init_z[nposnr]     = refcrd.z;
      refcrd                       = rbasis[i].getFinalTargetSite();
      tmp_rposn_final_xy[nposnr]   = {refcrd.x, refcrd.y};
      tmp_rposn_final_z[nposnr]    = refcrd.z;
      tmp_rposn_init_r[nposnr]     = rbasis[i].getInitialDisplacements();
      tmp_rposn_final_r[nposnr]    = rbasis[i].getFinalDisplacements();
      nposnr++;
    }
    else if (rord == 2) {
      tmp_rbond_i[nbondr]          = rbasis[i].getAtomIndex(1);
      tmp_rbond_j[nbondr]          = rbasis[i].getAtomIndex(2);
      tmp_rbond_init_step[nanglr]  = rbasis[i].getInitialStep();
      tmp_rbond_final_step[nanglr] = rbasis[i].getFinalStep();
      tmp_rbond_init_keq[nbondr]   = rbasis[i].getInitialStiffness();
      tmp_rbond_final_keq[nbondr]  = rbasis[i].getFinalStiffness();
      tmp_rbond_init_r[nposnr]     = rbasis[i].getInitialDisplacements();
      tmp_rbond_final_r[nposnr]    = rbasis[i].getFinalDisplacements();
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
      tmp_rangl_init_r[nposnr]     = rbasis[i].getInitialDisplacements();
      tmp_rangl_final_r[nposnr]    = rbasis[i].getFinalDisplacements();
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
      tmp_rdihe_init_r[nposnr]     = rbasis[i].getInitialDisplacements();
      tmp_rdihe_final_r[nposnr]    = rbasis[i].getFinalDisplacements();
      ndiher++;
    }
  }
  
  // Load the integer data into Hybrid objects
  rposn_atoms.putHost(tmp_rposn_i);
  rbond_i_atoms.putHost(tmp_rbond_i);
  rbond_j_atoms.putHost(tmp_rbond_j);
  rangl_i_atoms.putHost(tmp_rangl_i);
  rangl_j_atoms.putHost(tmp_rangl_j);
  rangl_k_atoms.putHost(tmp_rangl_k);
  rdihe_i_atoms.putHost(tmp_rdihe_i);
  rdihe_j_atoms.putHost(tmp_rdihe_j);
  rdihe_k_atoms.putHost(tmp_rdihe_k);
  rdihe_l_atoms.putHost(tmp_rdihe_l);
  rposn_init_step.putHost(tmp_rposn_init_step);
  rposn_final_step.putHost(tmp_rposn_final_step);
  rbond_init_step.putHost(tmp_rbond_init_step);
  rbond_final_step.putHost(tmp_rbond_final_step);
  rangl_init_step.putHost(tmp_rangl_init_step);
  rangl_final_step.putHost(tmp_rangl_final_step);
  rdihe_init_step.putHost(tmp_rdihe_init_step);
  rdihe_final_step.putHost(tmp_rdihe_final_step);

  // Load double-precision (including double2 and double4) data into Hybrid objects
  rposn_init_z.putHost(tmp_rposn_init_z);
  rposn_final_z.putHost(tmp_rposn_final_z);
  rposn_init_keq.putHost(tmp_rposn_init_keq);
  rposn_final_keq.putHost(tmp_rposn_final_keq);
  rposn_init_xy.putHost(tmp_rposn_init_xy);
  rposn_final_xy.putHost(tmp_rposn_final_xy);
  rbond_init_keq.putHost(tmp_rbond_init_keq);
  rbond_final_keq.putHost(tmp_rbond_final_keq);
  rangl_init_keq.putHost(tmp_rangl_init_keq);
  rangl_final_keq.putHost(tmp_rangl_final_keq);
  rdihe_init_keq.putHost(tmp_rdihe_init_keq);
  rdihe_final_keq.putHost(tmp_rdihe_final_keq);
  rposn_init_r.putHost(tmp_rposn_init_r);
  rposn_final_r.putHost(tmp_rposn_final_r);
  rbond_init_r.putHost(tmp_rbond_init_r);
  rbond_final_r.putHost(tmp_rbond_final_r);
  rangl_init_r.putHost(tmp_rangl_init_r);
  rangl_final_r.putHost(tmp_rangl_final_r);
  rdihe_init_r.putHost(tmp_rdihe_init_r);
  rdihe_final_r.putHost(tmp_rdihe_final_r);

  // Copy double-precision data into single-precision format
  std::vector<float> sp_tmp_rposn_init_z(tmp_rposn_init_z.begin(), tmp_rposn_init_z.end());
  std::vector<float> sp_tmp_rposn_final_z(tmp_rposn_final_z.begin(), tmp_rposn_final_z.end());
  std::vector<float2> sp_tmp_rposn_init_keq(position_count);
  std::vector<float2> sp_tmp_rposn_final_keq(position_count);
  std::vector<float2> sp_tmp_rposn_init_xy(position_count);
  std::vector<float2> sp_tmp_rposn_final_xy(position_count);
  std::vector<float2> sp_tmp_rbond_init_keq(distance_count);
  std::vector<float2> sp_tmp_rbond_final_keq(distance_count);
  std::vector<float2> sp_tmp_rangl_init_keq(angle_count);
  std::vector<float2> sp_tmp_rangl_final_keq(angle_count);
  std::vector<float2> sp_tmp_rdihe_init_keq(dihedral_count);
  std::vector<float2> sp_tmp_rdihe_final_keq(dihedral_count);
  std::vector<float4> sp_tmp_rposn_init_r(position_count);
  std::vector<float4> sp_tmp_rposn_final_r(position_count);
  std::vector<float4> sp_tmp_rbond_init_r(distance_count);
  std::vector<float4> sp_tmp_rbond_final_r(distance_count);
  std::vector<float4> sp_tmp_rangl_init_r(angle_count);
  std::vector<float4> sp_tmp_rangl_final_r(angle_count);
  std::vector<float4> sp_tmp_rdihe_init_r(dihedral_count);
  std::vector<float4> sp_tmp_rdihe_final_r(dihedral_count);
  for (int i = 0; i < position_count; i++) {
    sp_tmp_rposn_init_keq[i]  = vtConv2f(tmp_rposn_init_keq[i]);
    sp_tmp_rposn_final_keq[i] = vtConv2f(tmp_rposn_final_keq[i]);
    sp_tmp_rposn_init_xy[i]   = vtConv2f(tmp_rposn_init_xy[i]);
    sp_tmp_rposn_final_xy[i]  = vtConv2f(tmp_rposn_final_xy[i]);
    sp_tmp_rposn_init_r[i]    = vtConv4f(tmp_rposn_init_r[i]);
    sp_tmp_rposn_final_r[i]   = vtConv4f(tmp_rposn_final_r[i]);
  }
  for (int i = 0; i < distance_count; i++) {
    sp_tmp_rbond_init_keq[i]  = vtConv2f(tmp_rbond_init_keq[i]);
    sp_tmp_rbond_final_keq[i] = vtConv2f(tmp_rbond_final_keq[i]);
    sp_tmp_rbond_init_r[i]    = vtConv4f(tmp_rbond_init_r[i]);
    sp_tmp_rbond_final_r[i]   = vtConv4f(tmp_rbond_final_r[i]);    
  }
  for (int i = 0; i < angle_count; i++) {
    sp_tmp_rangl_init_keq[i]  = vtConv2f(tmp_rangl_init_keq[i]);
    sp_tmp_rangl_final_keq[i] = vtConv2f(tmp_rangl_final_keq[i]);
    sp_tmp_rangl_init_r[i]    = vtConv4f(tmp_rangl_init_r[i]);
    sp_tmp_rangl_final_r[i]   = vtConv4f(tmp_rangl_final_r[i]);    
  }
  for (int i = 0; i < dihedral_count; i++) {
    sp_tmp_rdihe_init_keq[i]  = vtConv2f(tmp_rdihe_init_keq[i]);
    sp_tmp_rdihe_final_keq[i] = vtConv2f(tmp_rdihe_final_keq[i]);
    sp_tmp_rdihe_init_r[i]    = vtConv4f(tmp_rdihe_init_r[i]);
    sp_tmp_rdihe_final_r[i]   = vtConv4f(tmp_rdihe_final_r[i]);    
  }

  // Load single-precision data
  sp_rposn_init_z.putHost(sp_tmp_rposn_init_z);
  sp_rposn_final_z.putHost(sp_tmp_rposn_final_z);
  sp_rposn_init_keq.putHost(sp_tmp_rposn_init_keq);
  sp_rposn_final_keq.putHost(sp_tmp_rposn_final_keq);
  sp_rposn_init_xy.putHost(sp_tmp_rposn_init_xy);
  sp_rposn_final_xy.putHost(sp_tmp_rposn_final_xy);
  sp_rbond_init_keq.putHost(sp_tmp_rbond_init_keq);
  sp_rbond_final_keq.putHost(sp_tmp_rbond_final_keq);
  sp_rangl_init_keq.putHost(sp_tmp_rangl_init_keq);
  sp_rangl_final_keq.putHost(sp_tmp_rangl_final_keq);
  sp_rdihe_init_keq.putHost(sp_tmp_rdihe_init_keq);
  sp_rdihe_final_keq.putHost(sp_tmp_rdihe_final_keq);
  sp_rposn_init_r.putHost(sp_tmp_rposn_init_r);
  sp_rposn_final_r.putHost(sp_tmp_rposn_final_r);
  sp_rbond_init_r.putHost(sp_tmp_rbond_init_r);
  sp_rbond_final_r.putHost(sp_tmp_rbond_final_r);
  sp_rangl_init_r.putHost(sp_tmp_rangl_init_r);
  sp_rangl_final_r.putHost(sp_tmp_rangl_final_r);
  sp_rdihe_init_r.putHost(sp_tmp_rdihe_init_r);
  sp_rdihe_final_r.putHost(sp_tmp_rdihe_final_r);
}

//-------------------------------------------------------------------------------------------------
RestraintApparatus::RestraintApparatus(const RestraintApparatus &original) :
    restraint_count{original.restraint_count},
    position_count{original.position_count},
    distance_count{original.distance_count},
    angle_count{original.angle_count},
    dihedral_count{original.dihedral_count},
    time_based_restraints{original.time_based_restraints},
    rposn_atoms{original.rposn_atoms},
    rbond_i_atoms{original.rbond_i_atoms},
    rbond_j_atoms{original.rbond_i_atoms},
    rangl_i_atoms{original.rangl_i_atoms},
    rangl_j_atoms{original.rangl_i_atoms},
    rangl_k_atoms{original.rangl_k_atoms},
    rdihe_i_atoms{original.rdihe_i_atoms},
    rdihe_j_atoms{original.rdihe_i_atoms},
    rdihe_k_atoms{original.rdihe_k_atoms},
    rdihe_l_atoms{original.rdihe_l_atoms},
    rposn_init_step{original.rposn_init_step},
    rposn_final_step{original.rposn_final_step},
    rbond_init_step{original.rbond_init_step},
    rbond_final_step{original.rbond_final_step},
    rangl_init_step{original.rangl_init_step},
    rangl_final_step{original.rangl_final_step},
    rdihe_init_step{original.rdihe_init_step},
    rdihe_final_step{original.rdihe_final_step},
    int_data{original.int_data},
    rposn_init_keq{original.rposn_init_keq},
    rposn_final_keq{original.rposn_final_keq},
    rposn_init_r{original.rposn_init_r},
    rposn_final_r{original.rposn_final_r},
    rposn_init_xy{original.rposn_init_xy},
    rposn_init_z{original.rposn_init_z},
    rposn_final_xy{original.rposn_final_xy},
    rposn_final_z{original.rposn_final_z},
    rbond_init_keq{original.rbond_init_keq},
    rbond_final_keq{original.rbond_final_keq},
    rbond_init_r{original.rbond_init_r},
    rbond_final_r{original.rbond_final_r},
    rangl_init_keq{original.rangl_init_keq},
    rangl_final_keq{original.rangl_final_keq},
    rangl_init_r{original.rangl_init_r},
    rangl_final_r{original.rangl_final_r},
    rdihe_init_keq{original.rdihe_init_keq},
    rdihe_final_keq{original.rdihe_final_keq},
    rdihe_init_r{original.rdihe_init_r},
    rdihe_final_r{original.rdihe_final_r},
    double_data{original.double_data},
    double2_data{original.double2_data},
    double4_data{original.double4_data},
    sp_rposn_init_keq{original.sp_rposn_init_keq},
    sp_rposn_final_keq{original.sp_rposn_final_keq},
    sp_rposn_init_r{original.sp_rposn_init_r},
    sp_rposn_final_r{original.sp_rposn_final_r},
    sp_rposn_init_xy{original.sp_rposn_init_xy},
    sp_rposn_init_z{original.sp_rposn_init_z},
    sp_rposn_final_xy{original.sp_rposn_final_xy},
    sp_rposn_final_z{original.sp_rposn_final_z},
    sp_rbond_init_keq{original.sp_rbond_init_keq},
    sp_rbond_final_keq{original.sp_rbond_final_keq},
    sp_rbond_init_r{original.sp_rbond_init_r},
    sp_rbond_final_r{original.sp_rbond_final_r},
    sp_rangl_init_keq{original.sp_rangl_init_keq},
    sp_rangl_final_keq{original.sp_rangl_final_keq},
    sp_rangl_init_r{original.sp_rangl_init_r},
    sp_rangl_final_r{original.sp_rangl_final_r},
    sp_rdihe_init_keq{original.sp_rdihe_init_keq},
    sp_rdihe_final_keq{original.sp_rdihe_final_keq},
    sp_rdihe_init_r{original.sp_rdihe_init_r},
    sp_rdihe_final_r{original.sp_rdihe_final_r},
    float_data{original.float_data},
    float2_data{original.float2_data},
    float4_data{original.float4_data},
    ag_pointer{original.ag_pointer}
{
  // Repair pointers to ensure correct access to data already copied over inside the ARRAY-kind
  // Hybrid objects.
  allocate();
}

//-------------------------------------------------------------------------------------------------
RestraintApparatus::RestraintApparatus(RestraintApparatus &&original) :
    restraint_count{original.restraint_count},
    position_count{original.position_count},
    distance_count{original.distance_count},
    angle_count{original.angle_count},
    dihedral_count{original.dihedral_count},
    time_based_restraints{original.time_based_restraints},
    rposn_atoms{std::move(original.rposn_atoms)},
    rbond_i_atoms{std::move(original.rbond_i_atoms)},
    rbond_j_atoms{std::move(original.rbond_i_atoms)},
    rangl_i_atoms{std::move(original.rangl_i_atoms)},
    rangl_j_atoms{std::move(original.rangl_i_atoms)},
    rangl_k_atoms{std::move(original.rangl_k_atoms)},
    rdihe_i_atoms{std::move(original.rdihe_i_atoms)},
    rdihe_j_atoms{std::move(original.rdihe_i_atoms)},
    rdihe_k_atoms{std::move(original.rdihe_k_atoms)},
    rdihe_l_atoms{std::move(original.rdihe_l_atoms)},
    rposn_init_step{std::move(original.rposn_init_step)},
    rposn_final_step{std::move(original.rposn_final_step)},
    rbond_init_step{std::move(original.rbond_init_step)},
    rbond_final_step{std::move(original.rbond_final_step)},
    rangl_init_step{std::move(original.rangl_init_step)},
    rangl_final_step{std::move(original.rangl_final_step)},
    rdihe_init_step{std::move(original.rdihe_init_step)},
    rdihe_final_step{std::move(original.rdihe_final_step)},
    int_data{std::move(original.int_data)},
    rposn_init_keq{std::move(original.rposn_init_keq)},
    rposn_final_keq{std::move(original.rposn_final_keq)},
    rposn_init_r{std::move(original.rposn_init_r)},
    rposn_final_r{std::move(original.rposn_final_r)},
    rposn_init_xy{std::move(original.rposn_init_xy)},
    rposn_init_z{std::move(original.rposn_init_z)},
    rposn_final_xy{std::move(original.rposn_final_xy)},
    rposn_final_z{std::move(original.rposn_final_z)},
    rbond_init_keq{std::move(original.rbond_init_keq)},
    rbond_final_keq{std::move(original.rbond_final_keq)},
    rbond_init_r{std::move(original.rbond_init_r)},
    rbond_final_r{std::move(original.rbond_final_r)},
    rangl_init_keq{std::move(original.rangl_init_keq)},
    rangl_final_keq{std::move(original.rangl_final_keq)},
    rangl_init_r{std::move(original.rangl_init_r)},
    rangl_final_r{std::move(original.rangl_final_r)},
    rdihe_init_keq{std::move(original.rdihe_init_keq)},
    rdihe_final_keq{std::move(original.rdihe_final_keq)},
    rdihe_init_r{std::move(original.rdihe_init_r)},
    rdihe_final_r{std::move(original.rdihe_final_r)},
    double_data{std::move(original.double_data)},
    double2_data{std::move(original.double2_data)},
    double4_data{std::move(original.double4_data)},
    sp_rposn_init_keq{std::move(original.sp_rposn_init_keq)},
    sp_rposn_final_keq{std::move(original.sp_rposn_final_keq)},
    sp_rposn_init_r{std::move(original.sp_rposn_init_r)},
    sp_rposn_final_r{std::move(original.sp_rposn_final_r)},
    sp_rposn_init_xy{std::move(original.sp_rposn_init_xy)},
    sp_rposn_init_z{std::move(original.sp_rposn_init_z)},
    sp_rposn_final_xy{std::move(original.sp_rposn_final_xy)},
    sp_rposn_final_z{std::move(original.sp_rposn_final_z)},
    sp_rbond_init_keq{std::move(original.sp_rbond_init_keq)},
    sp_rbond_final_keq{std::move(original.sp_rbond_final_keq)},
    sp_rbond_init_r{std::move(original.sp_rbond_init_r)},
    sp_rbond_final_r{std::move(original.sp_rbond_final_r)},
    sp_rangl_init_keq{std::move(original.sp_rangl_init_keq)},
    sp_rangl_final_keq{std::move(original.sp_rangl_final_keq)},
    sp_rangl_init_r{std::move(original.sp_rangl_init_r)},
    sp_rangl_final_r{std::move(original.sp_rangl_final_r)},
    sp_rdihe_init_keq{std::move(original.sp_rdihe_init_keq)},
    sp_rdihe_final_keq{std::move(original.sp_rdihe_final_keq)},
    sp_rdihe_init_r{std::move(original.sp_rdihe_init_r)},
    sp_rdihe_final_r{std::move(original.sp_rdihe_final_r)},
    float_data{std::move(original.float_data)},
    float2_data{std::move(original.float2_data)},
    float4_data{std::move(original.float4_data)},
    ag_pointer{original.ag_pointer}
{}

//-------------------------------------------------------------------------------------------------
void RestraintApparatus::allocate() {

  // Allocate the necessary space and set pointers
  const int padded_posn_count = roundUp(position_count, warp_size_int);
  const int padded_bond_count = roundUp(distance_count, warp_size_int);
  const int padded_angl_count = roundUp(angle_count, warp_size_int);
  const int padded_dihe_count = roundUp(dihedral_count, warp_size_int);
  const int nint_elem = (3 * padded_posn_count) + (4 * padded_bond_count) +
                        (5 * padded_angl_count) + (6 * padded_dihe_count);
  const int ndbl_elem = 2 * padded_posn_count;
  const int ndbl2_elem = (4 * padded_posn_count) + (2 * padded_bond_count) +
                         (2 * padded_angl_count) + (2 * padded_dihe_count);
  const int ndbl4_elem = (2 * padded_posn_count) + (2 * padded_bond_count) +
                         (2 * padded_angl_count) + (2 * padded_dihe_count);

  // Resize all ARRAY-kind Hybrids to accept POINTER-kind Hybrids
  int_data.resize(nint_elem);
  double_data.resize(ndbl_elem);
  double2_data.resize(ndbl2_elem);
  double4_data.resize(ndbl4_elem);
  float_data.resize(ndbl_elem);
  float2_data.resize(ndbl2_elem);
  float4_data.resize(ndbl4_elem);

  // Set pointers to integer data
  size_t ic = 0LLU;
  rposn_atoms.setPointer(&int_data, ic, position_count);
  ic += padded_posn_count;
  rbond_i_atoms.setPointer(&int_data, ic, distance_count);
  ic += padded_bond_count;
  rbond_j_atoms.setPointer(&int_data, ic, distance_count);
  ic += padded_bond_count;
  rangl_i_atoms.setPointer(&int_data, ic, angle_count);
  ic += padded_angl_count;
  rangl_j_atoms.setPointer(&int_data, ic, angle_count);
  ic += padded_angl_count;
  rangl_k_atoms.setPointer(&int_data, ic, angle_count);
  ic += padded_angl_count;
  rdihe_i_atoms.setPointer(&int_data, ic, dihedral_count);
  ic += padded_dihe_count;
  rdihe_j_atoms.setPointer(&int_data, ic, dihedral_count);
  ic += padded_dihe_count;
  rdihe_k_atoms.setPointer(&int_data, ic, dihedral_count);
  ic += padded_dihe_count;
  rdihe_l_atoms.setPointer(&int_data, ic, dihedral_count);
  ic += padded_dihe_count;
  rposn_init_step.setPointer(&int_data, ic, position_count);
  ic += padded_posn_count;
  rposn_final_step.setPointer(&int_data, ic, position_count);
  ic += padded_posn_count;
  rbond_init_step.setPointer(&int_data, ic, distance_count);
  ic += padded_bond_count;
  rbond_final_step.setPointer(&int_data, ic, distance_count);
  ic += padded_bond_count;
  rangl_init_step.setPointer(&int_data, ic, angle_count);
  ic += padded_angl_count;
  rangl_final_step.setPointer(&int_data, ic, angle_count);
  ic += padded_angl_count;
  rdihe_init_step.setPointer(&int_data, ic, dihedral_count);
  ic += padded_dihe_count;
  rdihe_final_step.setPointer(&int_data, ic, dihedral_count);

  // Set pointers to double-precision data
  size_t dc = 0LLU;
  rposn_init_z.setPointer(&double_data, ic, position_count);
  dc += padded_posn_count;
  rposn_final_z.setPointer(&double_data, ic, position_count);

  // Set pointers to double2 data
  size_t d2c = 0LLU;
  rposn_init_keq.setPointer(&double2_data, d2c, position_count);
  d2c += padded_posn_count;
  rposn_final_keq.setPointer(&double2_data, d2c, position_count);
  d2c += padded_posn_count;
  rposn_init_xy.setPointer(&double2_data, d2c, position_count);
  d2c += padded_posn_count;
  rposn_final_xy.setPointer(&double2_data, d2c, position_count);
  d2c += padded_posn_count;
  rbond_init_keq.setPointer(&double2_data, d2c, distance_count);
  d2c += padded_bond_count;
  rbond_final_keq.setPointer(&double2_data, d2c, distance_count);
  d2c += padded_bond_count;
  rangl_init_keq.setPointer(&double2_data, d2c, angle_count);
  d2c += padded_angl_count;
  rangl_final_keq.setPointer(&double2_data, d2c, angle_count);
  d2c += padded_angl_count;
  rdihe_init_keq.setPointer(&double2_data, d2c, dihedral_count);
  d2c += padded_dihe_count;
  rdihe_final_keq.setPointer(&double2_data, d2c, dihedral_count);

  // Set pointers to double4 data
  size_t d4c = 0LLU;
  rposn_init_r.setPointer(&double4_data, d4c, position_count);
  d4c += padded_posn_count;
  rposn_final_r.setPointer(&double4_data, d4c, position_count);
  d4c += padded_posn_count;
  rbond_init_r.setPointer(&double4_data, d4c, distance_count);
  d4c += padded_bond_count;
  rbond_final_r.setPointer(&double4_data, d4c, distance_count);
  d4c += padded_bond_count;
  rangl_init_r.setPointer(&double4_data, d4c, angle_count);
  d4c += padded_angl_count;
  rangl_final_r.setPointer(&double4_data, d4c, angle_count);
  d4c += padded_angl_count;
  rdihe_init_r.setPointer(&double4_data, d4c, dihedral_count);
  d4c += padded_dihe_count;
  rdihe_final_r.setPointer(&double4_data, d4c, dihedral_count);

  // Set pointers to single-precision data
  size_t fc = 0LLU;
  sp_rposn_init_z.setPointer(&float_data, ic, position_count);
  fc += padded_posn_count;
  sp_rposn_final_z.setPointer(&float_data, ic, position_count);

  // Set pointers to float2 data
  size_t f2c = 0LLU;
  sp_rposn_init_keq.setPointer(&float2_data, f2c, position_count);
  f2c += padded_posn_count;
  sp_rposn_final_keq.setPointer(&float2_data, f2c, position_count);
  f2c += padded_posn_count;
  sp_rposn_init_xy.setPointer(&float2_data, f2c, position_count);
  f2c += padded_posn_count;
  sp_rposn_final_xy.setPointer(&float2_data, f2c, position_count);
  f2c += padded_posn_count;
  sp_rbond_init_keq.setPointer(&float2_data, f2c, distance_count);
  f2c += padded_bond_count;
  sp_rbond_final_keq.setPointer(&float2_data, f2c, distance_count);
  f2c += padded_bond_count;
  sp_rangl_init_keq.setPointer(&float2_data, f2c, angle_count);
  f2c += padded_angl_count;
  sp_rangl_final_keq.setPointer(&float2_data, f2c, angle_count);
  f2c += padded_angl_count;
  sp_rdihe_init_keq.setPointer(&float2_data, f2c, dihedral_count);
  f2c += padded_dihe_count;
  sp_rdihe_final_keq.setPointer(&float2_data, f2c, dihedral_count);

  // Set pointers to float4 data
  size_t f4c = 0LLU;
  sp_rposn_init_r.setPointer(&float4_data, f4c, position_count);
  f4c += padded_posn_count;
  sp_rposn_final_r.setPointer(&float4_data, f4c, position_count);
  f4c += padded_posn_count;
  sp_rbond_init_r.setPointer(&float4_data, f4c, distance_count);
  f4c += padded_bond_count;
  sp_rbond_final_r.setPointer(&float4_data, f4c, distance_count);
  f4c += padded_bond_count;
  sp_rangl_init_r.setPointer(&float4_data, f4c, angle_count);
  f4c += padded_angl_count;
  sp_rangl_final_r.setPointer(&float4_data, f4c, angle_count);
  f4c += padded_angl_count;
  sp_rdihe_init_r.setPointer(&float4_data, f4c, dihedral_count);
  f4c += padded_dihe_count;
  sp_rdihe_final_r.setPointer(&float4_data, f4c, dihedral_count);  
}

} // namespace restraints
} // namespace omni
