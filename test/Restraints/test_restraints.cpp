#include <vector>
#include "../../src/Constants/scaling.h"
#include "../../src/DataTypes/omni_vector_types.h"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/Parsing/parse.h"
#include "../../src/Parsing/polynumeric.h"
#include "../../src/Potential/scorecard.h"
#include "../../src/Potential/valence_potential.h"
#include "../../src/Random/random.h"
#include "../../src/Reporting/error_format.h"
#include "../../src/Restraints/bounded_restraint.h"
#include "../../src/Restraints/restraint_apparatus.h"
#include "../../src/Restraints/restraint_builder.h"
#include "../../src/Structure/local_arrangement.h"
#include "../../src/Topology/atomgraph.h"
#include "../../src/Trajectory/coordinateframe.h"
#include "../../src/Trajectory/phasespace.h"
#include "../../src/Trajectory/trajectory_enumerators.h"
#include "../../src/UnitTesting/unit_test.h"

using omni::constants::ExceptionResponse;
using omni::constants::tiny;
using omni::data_types::char4;
using omni::data_types::double2;
using omni::data_types::double4;
using omni::diskutil::osSeparator;
using omni::diskutil::DrivePathType;
using omni::diskutil::getDrivePathType;
using omni::energy::ScoreCard;
using omni::energy::EvaluateForce;
using omni::energy::evaluateRestraints;
using omni::errors::rtWarn;
using omni::parse::separateText;
using omni::parse::NumberFormat;
using omni::parse::polyNumericVector;
using omni::random::Xoshiro256ppGenerator;
using omni::structure::distance;
using omni::structure::angle;
using omni::structure::dihedral_angle;
using omni::topology::AtomGraph;
using omni::trajectory::CoordinateFrameReader;
using omni::trajectory::CoordinateFrameWriter;
using omni::trajectory::PhaseSpace;
using omni::trajectory::PhaseSpaceWriter;
using omni::trajectory::TrajectoryKind;
using namespace omni::restraints;
using namespace omni::testing;

//-------------------------------------------------------------------------------------------------
// Perform finite difference calculations of the energy and force contributions for a set of
// restraints.  This function will assemble a vector of interlaced forces on all atoms based on
// changes in the overall energy computed by a finite-difference approximation.
//
// Arguments:
//   ra:          The set of applicable restraints.  Note that no topology is required, as the
//                restraints constitute an energy function all their own, equivalent to an
//                AtomGraph albeit lacking some fundamental details of the system like atom and
//                residue counts.  Also, the restraint apparatus contains a pointer to the
//                relevant topology.
//   cfw:         System coordinates (must be modifiable)
//   fd_delta:    Finite difference move, in Angstroms
//-------------------------------------------------------------------------------------------------
std::vector<double> finiteDifferenceForces(const RestraintApparatus &ra,
                                           CoordinateFrameWriter *cfw,
                                           const double fd_delta = 0.0000001) {
  const RestraintApparatusDpReader rar = ra.dpData();
  const CoordinateFrameReader cfr(*cfw);
  ScoreCard sc(1);
  std::vector<double> result(cfr.natom * 3, 0.0);

  // Map out which atoms are touched by any restraints
  std::vector<bool> atom_touched(cfr.natom, false);
  for (int i = 0; i < rar.nposn; i++) {
    atom_touched[rar.rposn_atoms[i]] = true;
  }
  for (int i = 0; i < rar.nbond; i++) {
    atom_touched[rar.rbond_i_atoms[i]] = true;
    atom_touched[rar.rbond_j_atoms[i]] = true;
  }
  for (int i = 0; i < rar.nangl; i++) {
    atom_touched[rar.rangl_i_atoms[i]] = true;
    atom_touched[rar.rangl_j_atoms[i]] = true;
    atom_touched[rar.rangl_k_atoms[i]] = true;
  }
  for (int i = 0; i < rar.ndihe; i++) {
    atom_touched[rar.rdihe_i_atoms[i]] = true;
    atom_touched[rar.rdihe_j_atoms[i]] = true;
    atom_touched[rar.rdihe_k_atoms[i]] = true;
    atom_touched[rar.rdihe_l_atoms[i]] = true;
  }

  // Loop over all atoms
  const double e0 = evaluateRestraints(rar, cfr, &sc);
  for (int i = 0; i < cfr.natom; i++) {
    if (atom_touched[i] == false) {
      continue;
    }
    cfw->xcrd[i] += fd_delta;
    result[(3 * i)    ] = (e0 - evaluateRestraints(rar, cfr, &sc)) / fd_delta;
    cfw->xcrd[i] -= fd_delta;
    cfw->ycrd[i] += fd_delta;
    result[(3 * i) + 1] = (e0 - evaluateRestraints(rar, cfr, &sc)) / fd_delta;
    cfw->ycrd[i] -= fd_delta;
    cfw->zcrd[i] += fd_delta;
    result[(3 * i) + 2] = (e0 - evaluateRestraints(rar, cfr, &sc)) / fd_delta;
    cfw->zcrd[i] -= fd_delta;
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
// Break a series of restraints into its numerical elements.
//
// Arguments:
//   rst_in:      Input vector of bounded (flat-bottom, harmonic bimodal restraints
//   atom_idx:    Vector of atom positions (returned)
//   k_values:    Vector of k2(0), k3(0), k2(1), k3(1), ... for each restraint in the input list
///               (filled and returned)
//   r_values:    Vector of r1(0), r2(0), r3(0), r4(0), r1(1), r2(1), ... for each restraint in
//                the input list  (filled and returned)
//-------------------------------------------------------------------------------------------------
void digestRestraintList(const std::vector<BoundedRestraint> &rst_in, std::vector<int> *atom_idx,
                         std::vector<double> *k_values, std::vector<double> *r_values)
{
  for (size_t i = 0; i < rst_in.size(); i++) {
    const int orderi = rst_in[i].getOrder();
    for (int j = 0; j < orderi; j++) {
      atom_idx->push_back(rst_in[i].getAtomIndex(j + 1));
    }
    const double2 k23i = rst_in[i].getInitialStiffness();
    const double2 k23f = rst_in[i].getFinalStiffness();
    const double4 ri = rst_in[i].getInitialDisplacements();
    const double4 rf = rst_in[i].getFinalDisplacements();
    k_values->push_back(k23i.x);
    k_values->push_back(k23i.y);
    k_values->push_back(k23f.x);
    k_values->push_back(k23f.y);
    r_values->push_back(ri.x);
    r_values->push_back(ri.y);
    r_values->push_back(ri.z);
    r_values->push_back(ri.w);
    r_values->push_back(rf.x);
    r_values->push_back(rf.y);
    r_values->push_back(rf.z);
    r_values->push_back(rf.w);
  }
}

//-------------------------------------------------------------------------------------------------
int main(int argc, char* argv[]) {

  // Some baseline initialization
  TestEnvironment oe(argc, argv);

  // Section 1
  section("Test the restraint builders");

  // Section 2
  section("Test the RestraintApparatus object");

  // Section 3
  section("Test restraint potential functions");

  // The restraint builders are automated tools for applying restraints to a structure based on
  // its features, to guide motion during dynamics or energy minimizations.
  section(1);
  const char osc = osSeparator();
  const std::string topology_base = oe.getOmniSourcePath() + osc + "test" + osc + "Namelists" +
                                    osc + "topol";
  const std::string coordinate_base = oe.getOmniSourcePath() + osc + "test" + osc + "Namelists" +
                                      osc + "coord";
  const std::string gly_tyr_top_name = topology_base + osc + "gly_tyr.top";
  const std::string gly_tyr_crd_name = coordinate_base + osc + "gly_tyr.inpcrd";
  const std::string gly_lys_top_name = topology_base + osc + "gly_lys.top";
  const std::string gly_lys_crd_name = coordinate_base + osc + "gly_lys.inpcrd";
  const bool input_exists = (getDrivePathType(gly_tyr_top_name) == DrivePathType::FILE &&
                             getDrivePathType(gly_tyr_crd_name) == DrivePathType::FILE &&
                             getDrivePathType(gly_lys_top_name) == DrivePathType::FILE &&
                             getDrivePathType(gly_lys_crd_name) == DrivePathType::FILE);
  const TestPriority do_tests = (input_exists) ? TestPriority::CRITICAL : TestPriority::ABORT;
  const std::string snp_name = oe.getOmniSourcePath() + osc + "test" + osc + "Restraints" + osc +
                               "rst_output.m";
  const bool snp_exists = (getDrivePathType(snp_name) == DrivePathType::FILE);
  const TestPriority do_snps = (snp_exists) ? TestPriority::CRITICAL : TestPriority::ABORT;
  AtomGraph gy_ag(gly_tyr_top_name, ExceptionResponse::SILENT);
  PhaseSpace gy_ps(gly_tyr_crd_name);
  AtomGraph gk_ag(gly_lys_top_name, ExceptionResponse::SILENT);
  PhaseSpace gk_ps(gly_lys_crd_name);
  const CoordinateFrameReader gy_cfr(gy_ps);
  const CoordinateFrameReader gk_cfr(gk_ps);
  CoordinateFrameWriter gy_cfw(&gy_ps);
  CoordinateFrameWriter gk_cfw(&gk_ps);
  const ChemicalFeatures gy_chemfe(&gy_ag, gy_cfr);
  const ChemicalFeatures gk_chemfe(&gk_ag, gk_cfr);
  const std::vector<BoundedRestraint> gy_pos_rstr =
    applyPositionalRestraints(&gy_ag, gy_cfr, ":* & @CA,N,C,O", 1.0);
  const std::vector<BoundedRestraint> gk_pos_rstr =
    applyPositionalRestraints(&gk_ag, gk_cfr, ":* & @CA,N,C,O,NZ", 2.5);
  const std::vector<BoundedRestraint> gk_hb_rstr =
    applyHydrogenBondPreventors(&gk_ag, gk_cfr, 64.0);
  const std::vector<BoundedRestraint> gy_hold_rstr =
    applyHoldingRestraints(&gy_ag, gy_cfr, ":TYR & ! @CA,N,C,O,CB & ! @H*", 2.5);
  const std::vector<int> gy_pos_atom_ans = {  4,  5,  6,  8, 11, 12, 13, 15, 32, 33, 34 };
  const std::vector<int> gk_pos_atom_ans = {  4,  5,  6,  8, 11, 12, 13, 15, 29, 33, 34, 35 };
  const std::vector<int> gk_hb_atom_ans = {
     6,  5,  6, 12,  6, 13,  6, 34,  6, 35, 13,  5, 13,  6, 13, 12, 13, 34, 13, 35, 29,  5, 29,
     6, 29, 12, 29, 13, 29, 34, 29, 35, 35,  5, 35,  6, 35, 12, 35, 13, 35, 34 };
  std::vector<double> k_elem = { 0.0, 1.0, 0.0, 1.0 };
  std::vector<double> r_elem = { 0.0, 0.0, 0.0, 16.0, 0.0, 0.0, 0.0, 16.0 };
  std::vector<double> gy_pos_k_ans, gk_pos_k_ans, gy_pos_r_ans, gk_pos_r_ans;
  std::vector<double> gk_hb_k_ans, gk_hb_r_ans;
  for (int i = 0; i < 11; i++) {
    gy_pos_k_ans.insert(gy_pos_k_ans.end(), k_elem.begin(), k_elem.end());
    gy_pos_r_ans.insert(gy_pos_r_ans.end(), r_elem.begin(), r_elem.end());
  }
  k_elem[1] = 2.5;
  k_elem[3] = 2.5;
  for (int i = 0; i < 12; i++) {
    gk_pos_k_ans.insert(gk_pos_k_ans.end(), k_elem.begin(), k_elem.end());
    gk_pos_r_ans.insert(gk_pos_r_ans.end(), r_elem.begin(), r_elem.end());
  }
  k_elem[0] =   64.0;
  k_elem[1] =    0.0;
  k_elem[2] =   64.0;
  k_elem[3] =    0.0;
  r_elem[1] =    3.1;
  r_elem[2] = 1000.0;
  r_elem[3] = 1100.0;
  r_elem[5] =    3.1;
  r_elem[6] = 1000.0;
  r_elem[7] = 1100.0;
  for (int i = 0; i < 21; i++) {
    gk_hb_k_ans.insert(gk_hb_k_ans.end(), k_elem.begin(), k_elem.end());
    gk_hb_r_ans.insert(gk_hb_r_ans.end(), r_elem.begin(), r_elem.end());
  }
  std::vector<int> gy_pos_atoms, gk_pos_atoms, gk_hb_atoms, gy_hold_atoms;
  std::vector<double> gy_pos_k, gk_pos_k, gy_pos_r, gk_pos_r, gk_hb_k, gk_hb_r;
  std::vector<double> gy_hold_k, gy_hold_r;
  digestRestraintList(gy_pos_rstr, &gy_pos_atoms, &gy_pos_k, &gy_pos_r);
  digestRestraintList(gk_pos_rstr, &gk_pos_atoms, &gk_pos_k, &gk_pos_r);
  digestRestraintList(gk_hb_rstr, &gk_hb_atoms, &gk_hb_k, &gk_hb_r);
  digestRestraintList(gy_hold_rstr, &gy_hold_atoms, &gy_hold_k, &gy_hold_r);
  check(gy_pos_atoms, RelationalOperator::EQUAL, gy_pos_atom_ans, "Positional restraints for the "
        "Gly-Tyr system were not applied to the expected atoms.", do_tests);
  check(gy_pos_k, RelationalOperator::EQUAL, gy_pos_k_ans, "Positional restraints for the "
        "Gly-Tyr system do not incorporate the expected stiffnesses.", do_tests);
  check(gy_pos_r, RelationalOperator::EQUAL, gy_pos_r_ans, "Positional restraints for the "
        "Gly-Tyr system do not incorporate the expected displacements.", do_tests);
  check(gk_pos_atoms, RelationalOperator::EQUAL, gk_pos_atom_ans, "Positional restraints for the "
        "Gly-Lys system were not applied to the expected atoms.", do_tests);
  check(gk_pos_k, RelationalOperator::EQUAL, gk_pos_k_ans, "Positional restraints for the "
        "Gly-Lys system do not incorporate the expected stiffnesses.", do_tests);
  check(gk_pos_r, RelationalOperator::EQUAL, gk_pos_r_ans, "Positional restraints for the "
        "Gly-Lys system do not incorporate the expected displacements.", do_tests);
  check(gk_hb_atoms, RelationalOperator::EQUAL, gk_hb_atom_ans, "Hydrogen bond preventors for the "
        "Gly-Lys system were not applied to the expected atoms.", do_tests);
  check(gk_hb_k, RelationalOperator::EQUAL, gk_hb_k_ans, "Hydrogen bond preventors for the "
        "Gly-Lys system do not incorporate the expected stiffnesses.", do_tests);
  check(gk_hb_r, RelationalOperator::EQUAL, gk_hb_r_ans, "Hydrogen bond preventors for the "
        "Gly-Lys system do not incorporate the expected displacements.", do_tests);
  snapshot(snp_name, polyNumericVector(gy_hold_atoms), "gy_hold_indices", NumberFormat::INTEGER,
           "Holding restraints for the Gly-Tyr system were not applied to the expected atoms.",
           oe.takeSnapshot(), 1.0e-4, 1.0e-6, PrintSituation::OVERWRITE, do_snps);
  snapshot(snp_name, polyNumericVector(gy_hold_k), "gy_hold_kval", NumberFormat::STANDARD_REAL,
           "Holding restraints for the Gly-Tyr system do not incorporate the expected "
           "stiffnesses.", oe.takeSnapshot(), 1.0e-4, 1.0e-6, PrintSituation::APPEND, do_snps);
  snapshot(snp_name, polyNumericVector(gy_hold_r), "gy_hold_rval", NumberFormat::STANDARD_REAL,
           "Holding restraints for the Gly-Tyr system do not incorporate the expected "
           "displacements.", oe.takeSnapshot(), 1.0e-4, 1.0e-6, PrintSituation::APPEND, do_snps);
  CHECK_THROWS_SOFT(BoundedRestraint trick_br(":GLY & @CA", ":GLY & @C", ":LYS & @PH",
                                              ":LYS & @CA", &gk_ag, &gk_chemfe, gk_cfr, 0, 0,
                                              1.0, 1.0, 0.5, 1.5, 1.6, 2.0, 1.0, 1.0, 0.5, 1.5,
                                              1.6, 2.0), "A four-point dihedral angle restraint "
                    "was created based on a nonsensical atom mask.", do_tests);
  CHECK_THROWS_SOFT(BoundedRestraint trick_br(":GLY & @CA,C", ":GLY & @C", ":LYS & @N",
                                              ":LYS & @CA", &gk_ag, &gk_chemfe, gk_cfr, 0, 0,
                                              1.0, 1.0, 0.5, 1.5, 1.6, 2.0, 1.0, 1.0, 0.5, 1.5,
                                              1.6, 2.0), "A four-point dihedral angle restraint "
                    "was created based on an atom mask specifying multiple atoms.", do_tests);
  CHECK_THROWS_SOFT(BoundedRestraint trick_br(-1, -1, -1, -1, &gk_ag, 0, 0, 1.0, 1.0, 0.5, 1.5,
                                              1.6, 2.0, 1.0, 1.0, 0.5, 1.5, 1.6, 2.0),
                    "A restraint was created without a single valid atom.", do_tests);
  CHECK_THROWS_SOFT(BoundedRestraint trick_br(0, 102, 3, 5, &gk_ag, 0, 0, 1.0, 1.0, 0.5, 1.5,
                                              1.6, 2.0, 1.0, 1.0, 0.5, 1.5, 1.6, 2.0),
                    "A restraint was created with an invalid atom index.", do_tests);
  CHECK_THROWS_SOFT(BoundedRestraint trick_br(0, 1, 3, 5, &gk_ag, 0, 0, 1.0, 1.0, 0.5, 1.5,
                                              1.1, 2.0, 1.0, 1.0, 0.5, 1.5, 1.6, 2.0),
                    "A restraint was created with zig-zagging displacment bounds.", do_tests);
  CHECK_THROWS_SOFT(gy_pos_rstr[2].getAtomIndex(2), "An invalid participating atom index (2) was "
                    "produced for a positional restraint.", do_tests);
  CHECK_THROWS_SOFT(gy_hold_rstr[3].getAtomIndex(5), "An invalid participating atom index (5) was "
                    "produced for a dihedral angle restraint.", do_tests);

  // Test the sets of restraints added by each method
  section(2);
  RestraintApparatus gy_ra(gy_pos_rstr);
  check(gy_ra.getPositionalRestraintCount(), RelationalOperator::EQUAL, gy_pos_rstr.size(),
        "The RestraintApparatus did not incorporate the correct number of positional restraints.",
        do_tests);
  check(gy_ra.getTimeDependence() == false, "The RestraintApparatus reports time dependence in "
        "its restraints when none exists.", do_tests);
  const std::vector<BoundedRestraint> reproduce_gy_pos_rstr = gy_ra.getRestraintList();
  std::vector<int> rp_gy_pos_atoms;
  std::vector<double> rp_gy_pos_k, rp_gy_pos_r;
  digestRestraintList(reproduce_gy_pos_rstr, &rp_gy_pos_atoms, &rp_gy_pos_k, &rp_gy_pos_r);
  check(rp_gy_pos_atoms, RelationalOperator::EQUAL, gy_pos_atoms, "Positional restraints for the "
        "Gly-Tyr system were not applied to the correct atoms in the RestraintApparaus.",
        do_tests);
  check(rp_gy_pos_k, RelationalOperator::EQUAL, gy_pos_k, "Positional restraints for the "
        "Gly-Tyr system do not have the correct stiffnesses in the RestraintApparaus.", do_tests);
  check(rp_gy_pos_r, RelationalOperator::EQUAL, gy_pos_r, "Positional restraints for the "
        "Gly-Tyr system do not have the correct displacements in the RestraintApparaus.",
        do_tests);
  check(gy_ra.getTopologyPointer() == &gy_ag, "The topology pointer returned by a "
        "RestraintApparatus is not a pointer to the original topology.", do_tests);
  CHECK_THROWS(RestraintApparatus(std::vector<BoundedRestraint>()), "A restraint "
               "apparatus was created without a valid topology pointer.");
  std::vector<BoundedRestraint> mixed_system_rstr(gy_pos_rstr.begin(), gy_pos_rstr.end());
  mixed_system_rstr.insert(mixed_system_rstr.end(), gk_hb_rstr.begin(), gk_hb_rstr.end());
  CHECK_THROWS(RestraintApparatus gy_test_ra(mixed_system_rstr), "A restraint apparatus was "
               "created based on restraints with inconsistent topologies.");

  // Test the restraint functions with finite difference approximations.  Begin with simple
  // harmonic wells that, while they may have distinct stiffnesses on either side of the minimum,
  // have no flat bottom region and no practical limit on the extent of each harmonic region.
  section(3);
  std::vector<BoundedRestraint> posn_restraints;
  const int natom_gk = gk_ag.getAtomCount();
  for (int i = 0; i < natom_gk; i++) {
    posn_restraints.emplace_back(i, &gk_ag, gk_cfr, 1.0, 1.0, 0.0, 0.0, 0.0, 50.0);
  }
  RestraintApparatus gk_posn_ra(posn_restraints);
  std::vector<BoundedRestraint> bond_restraints;
  for (int i = 0; i < natom_gk - 4; i += 2) {

    // Compute the original particle : particle distance and use that as the equilibrium
    const double leq = distance(i, i + 3, gk_cfr);
    bond_restraints.emplace_back(i, i + 3, &gk_ag, 0.8, 2.6, 0.0, leq, leq, leq + 500.0);
  }
  RestraintApparatus gk_bond_ra(bond_restraints);
  std::vector<BoundedRestraint> angl_restraints;
  for (int i = 0; i < natom_gk - 4; i += 3) {
    const double teq = angle(i, i + 2, i + 3, gk_cfr);
    angl_restraints.emplace_back(i, i + 2, i + 3, &gk_ag, 1.4, 3.7, teq - 3.0, teq, teq,
                                 teq + 3.0);
  }
  RestraintApparatus gk_angl_ra(angl_restraints);
  std::vector<BoundedRestraint> dihe_restraints;
  for (int i = 0; i < natom_gk - 4; i += 4) {
    const double phi_eq = dihedral_angle(i, i + 1, i + 2, i + 3, gk_cfr);
    dihe_restraints.emplace_back(i, i + 1, i + 2, i + 3, &gk_ag, 1.4, 3.7, phi_eq - 3.0, phi_eq,
                                 phi_eq, phi_eq + 3.0);
  }
  RestraintApparatus gk_dihe_ra(dihe_restraints);

  // Initialize a random number generator and create new restraint collections with linear
  // regions likely to be accessible when the molecular configuration becomes perturbed.
  Xoshiro256ppGenerator xsr(9183053);
  std::vector<BoundedRestraint> posn_scramble_restraints;
  for (int i = 0; i < natom_gk; i++) {
    const double r2 = 0.25 * xsr.uniformRandomNumber();
    posn_scramble_restraints.emplace_back(i, &gk_ag, gk_cfr, 1.4, 3.7, 0.0, r2, r2 + 0.5, 50.0);
  }
  RestraintApparatus gk_posn_scrm_ra(posn_scramble_restraints);
  std::vector<BoundedRestraint> bond_scramble_restraints;
  for (int i = 0; i < natom_gk - 4; i += 2) {
    const double leq = distance(i, i + 3, gk_cfr) + 0.25 * xsr.uniformRandomNumber();
    const double k2 = 1.7 * (0.5 - xsr.uniformRandomNumber());
    const double k3 = 1.4 * (0.5 - xsr.uniformRandomNumber());
    bond_scramble_restraints.emplace_back(i, i + 3, &gk_ag, k2, k3, leq - 0.25, leq, leq + 0.5,
                                          leq + 0.85);
  }
  RestraintApparatus gk_bond_scrm_ra(bond_scramble_restraints);
  std::vector<BoundedRestraint> angl_scramble_restraints;
  for (int i = 0; i < natom_gk - 4; i += 2) {
    const double teq = angle(i, i + 2, i + 3, gk_cfr);
    const double k2 = 1.9 * (0.5 - xsr.uniformRandomNumber());
    const double k3 = 1.8 * (0.5 - xsr.uniformRandomNumber());
    angl_scramble_restraints.emplace_back(i, i + 2, i + 3, &gk_ag, k2, k3, teq - 0.15, teq,
                                          teq + 0.15, teq + 5.0);
  }
  RestraintApparatus gk_angl_scrm_ra(angl_scramble_restraints);
  std::vector<BoundedRestraint> dihe_scramble_restraints;
  for (int i = 0; i < natom_gk - 4; i += 2) {
    const double phi_eq = dihedral_angle(i, i + 1, i + 2, i + 3, gk_cfr);
    const double k2 = 1.1 * (0.5 - xsr.uniformRandomNumber());
    const double k3 = 1.3 * (0.5 - xsr.uniformRandomNumber());
    dihe_scramble_restraints.emplace_back(i, i + 1, i + 2, i + 3, &gk_ag, k2, k3, phi_eq - 0.75,
                                          phi_eq - 0.15, phi_eq + 0.10, phi_eq + 0.55);
  }
  RestraintApparatus gk_dihe_scrm_ra(dihe_scramble_restraints);
  
  // Perturb the Gly-Lys system so that the test restraints, all of which are at equilibrium in the
  // original system configuration, restraints will exert forces on the particles they affect.
  for (int i = 0; i < gk_cfw.natom; i++) {
    gk_cfw.xcrd[i] += xsr.gaussianRandomNumber();
    gk_cfw.ycrd[i] += xsr.gaussianRandomNumber();
    gk_cfw.zcrd[i] += xsr.gaussianRandomNumber();
  }
  ScoreCard sc(1);
  evaluateRestraints(gk_posn_ra, &gk_ps, &sc, EvaluateForce::YES);
  std::vector<double> analytic_frc = gk_ps.getInterlacedCoordinates(TrajectoryKind::FORCES);
  std::vector<double> fd_frc = finiteDifferenceForces(gk_posn_ra, &gk_cfw);
  check(analytic_frc, RelationalOperator::EQUAL, Approx(fd_frc).margin(2.0e-6), "Forces computed "
        "by analytic positional restraints do not agree with finite-difference approximations in "
        "the Gly-Lys system.", do_tests);
  gk_ps.initializeForces();
  evaluateRestraints(gk_bond_ra, &gk_ps, &sc, EvaluateForce::YES);
  analytic_frc = gk_ps.getInterlacedCoordinates(TrajectoryKind::FORCES);
  fd_frc = finiteDifferenceForces(gk_bond_ra, &gk_cfw);  
  check(analytic_frc, RelationalOperator::EQUAL, Approx(fd_frc).margin(2.0e-6), "Forces computed "
        "by analytic distance restraints do not agree with finite-difference approximations in "
        "the Gly-Lys system.", do_tests);
  gk_ps.initializeForces();
  evaluateRestraints(gk_angl_ra, &gk_ps, &sc, EvaluateForce::YES);
  analytic_frc = gk_ps.getInterlacedCoordinates(TrajectoryKind::FORCES);
  fd_frc = finiteDifferenceForces(gk_angl_ra, &gk_cfw);
  check(analytic_frc, RelationalOperator::EQUAL, Approx(fd_frc).margin(2.0e-6), "Forces computed "
        "by analytic angle restraints do not agree with finite-difference approximations in "
        "the Gly-Lys system.", do_tests);
  gk_ps.initializeForces();
  evaluateRestraints(gk_dihe_ra, &gk_ps, &sc, EvaluateForce::YES);
  analytic_frc = gk_ps.getInterlacedCoordinates(TrajectoryKind::FORCES);
  fd_frc = finiteDifferenceForces(gk_dihe_ra, &gk_cfw, 0.00000001);
  check(analytic_frc, RelationalOperator::EQUAL, Approx(fd_frc).margin(1.0e-4), "Forces computed "
        "by analytic dihedral angle restraints do not agree with finite-difference approximations "
        "in the Gly-Lys system.", do_tests);
  gk_ps.initializeForces();
  evaluateRestraints(gk_posn_scrm_ra, &gk_ps, &sc, EvaluateForce::YES);
  analytic_frc = gk_ps.getInterlacedCoordinates(TrajectoryKind::FORCES);
  fd_frc = finiteDifferenceForces(gk_posn_scrm_ra, &gk_cfw, 0.0000001);
  check(analytic_frc, RelationalOperator::EQUAL, Approx(fd_frc).margin(2.0e-6), "Forces computed "
        "by analytic position restraints bearing all of the permutations of the bimodal, "
        "flat-bottom well format do not agree with finite-difference approximations in the "
        "Gly-Lys system.", do_tests);
  gk_ps.initializeForces();
  evaluateRestraints(gk_bond_scrm_ra, &gk_ps, &sc, EvaluateForce::YES);
  analytic_frc = gk_ps.getInterlacedCoordinates(TrajectoryKind::FORCES);
  fd_frc = finiteDifferenceForces(gk_bond_scrm_ra, &gk_cfw, 0.0000001);
  check(analytic_frc, RelationalOperator::EQUAL, Approx(fd_frc).margin(2.0e-6), "Forces computed "
        "by analytic distance restraints bearing all of the permutations of the bimodal, "
        "flat-bottom well format do not agree with finite-difference approximations in the "
        "Gly-Lys system.", do_tests);
  gk_ps.initializeForces();
  evaluateRestraints(gk_angl_scrm_ra, &gk_ps, &sc, EvaluateForce::YES);
  analytic_frc = gk_ps.getInterlacedCoordinates(TrajectoryKind::FORCES);
  fd_frc = finiteDifferenceForces(gk_angl_scrm_ra, &gk_cfw, 0.0000001);
  check(analytic_frc, RelationalOperator::EQUAL, Approx(fd_frc).margin(2.0e-6), "Forces computed "
        "by analytic angle restraints bearing all of the permutations of the bimodal, "
        "flat-bottom well format do not agree with finite-difference approximations in the "
        "Gly-Lys system.", do_tests);
  gk_ps.initializeForces();
  evaluateRestraints(gk_dihe_scrm_ra, &gk_ps, &sc, EvaluateForce::YES);
  analytic_frc = gk_ps.getInterlacedCoordinates(TrajectoryKind::FORCES);
  fd_frc = finiteDifferenceForces(gk_dihe_scrm_ra, &gk_cfw, 0.0000001);
  check(analytic_frc, RelationalOperator::EQUAL, Approx(fd_frc).margin(2.0e-6), "Forces computed "
        "by analytic dihedral angle restraints bearing all of the permutations of the bimodal, "
        "flat-bottom well format do not agree with finite-difference approximations in the "
        "Gly-Lys system.", do_tests);

  // Perturb the Gly-Tyr system so that the test restraints, all of which are at equilibrium in the
  // original system configuration, restraints will exert forces on the particles they affect.
  for (int i = 0; i < gy_cfw.natom; i++) {
    gy_cfw.xcrd[i] += xsr.gaussianRandomNumber();
    gy_cfw.ycrd[i] += xsr.gaussianRandomNumber();
    gy_cfw.zcrd[i] += xsr.gaussianRandomNumber();
  }
  evaluateRestraints(gy_ra, &gy_ps, &sc, EvaluateForce::YES);
  analytic_frc = gy_ps.getInterlacedCoordinates(TrajectoryKind::FORCES);
  fd_frc = finiteDifferenceForces(gy_ra, &gy_cfw);
  check(analytic_frc, RelationalOperator::EQUAL, Approx(fd_frc).margin(2.0e-6), "Forces computed "
        "by analytic positional restraints do not agree with finite-difference approximations in "
        "the Gly-Tyr system.", do_tests);
  
  // Add more restraints and try again
  section(2);
  gy_ra.addRestraints(gy_hold_rstr);
  check(gy_ra.getTotalRestraintCount(), RelationalOperator::EQUAL, gy_pos_rstr.size() +
        gy_hold_rstr.size(), "The RestraintApparatus did not report the correct number of "
        "total restraints after incorporating holding restraints for the Gly-Tyr system.",
        do_tests);
  section(3);
  gy_ps.initializeForces();
  evaluateRestraints(gy_ra, &gy_ps, &sc, EvaluateForce::YES);
  analytic_frc = gy_ps.getInterlacedCoordinates(TrajectoryKind::FORCES);
  fd_frc = finiteDifferenceForces(gy_ra, &gy_cfw, 0.00000001);
  check(analytic_frc, RelationalOperator::EQUAL, Approx(fd_frc).margin(2.0e-5), "Forces computed "
        "by analytic positional and holding restraints do not agree with finite-difference "
        "approximations.", do_tests);
  
  // Summary evaluation
  printTestSummary(oe.getVerbosity());

  return 0;
}
