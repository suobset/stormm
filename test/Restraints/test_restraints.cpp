#include <vector>
#include "../../src/Constants/scaling.h"
#include "../../src/DataTypes/omni_vector_types.h"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/Parsing/parse.h"
#include "../../src/Parsing/polynumeric.h"
#include "../../src/Reporting/error_format.h"
#include "../../src/Restraints/bounded_restraint.h"
#include "../../src/Restraints/restraint_apparatus.h"
#include "../../src/Restraints/restraint_builder.h"
#include "../../src/Topology/atomgraph.h"
#include "../../src/Trajectory/coordinateframe.h"
#include "../../src/Trajectory/phasespace.h"
#include "../../src/UnitTesting/unit_test.h"

using omni::constants::ExceptionResponse;
using omni::constants::tiny;
using omni::data_types::char4;
using omni::data_types::double2;
using omni::data_types::double4;
using omni::diskutil::osSeparator;
using omni::diskutil::DrivePathType;
using omni::diskutil::getDrivePathType;
using omni::errors::rtWarn;
using omni::parse::separateText;
using omni::parse::NumberFormat;
using omni::parse::polyNumericVector;
using omni::topology::AtomGraph;
using omni::trajectory::CoordinateFrameReader;
using omni::trajectory::PhaseSpace;
using namespace omni::restraints;
using namespace omni::testing;

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
  section("Test restraint potential functions");

  // Section 3
  section("Test the RestraintApparatus object");

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
  const bool input_exists = (getDrivePathType(gly_tyr_top_name) == DrivePathType::FILE,
                             getDrivePathType(gly_tyr_crd_name) == DrivePathType::FILE,
                             getDrivePathType(gly_lys_crd_name) == DrivePathType::FILE,
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
  
  // Summary evaluation
  printTestSummary(oe.getVerbosity());

  return 0;
}
