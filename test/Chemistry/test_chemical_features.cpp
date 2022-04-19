#include "../../src/Chemistry/chemical_features.h"
#include "../../src/Chemistry/chemistry_enumerators.h"
#include "../../src/Constants/behavior.h"
#include "../../src/Constants/symbol_values.h"
#include "../../src/DataTypes/omni_vector_types.h"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/Parsing/polynumeric.h"
#include "../../src/Reporting/error_format.h"
#include "../../src/Topology/atomgraph.h"
#include "../../src/Topology/atomgraph_abstracts.h"
#include "../../src/Trajectory/phasespace.h"
#include "../../src/Trajectory/trajectory_enumerators.h"
#include "../../src/UnitTesting/stopwatch.h"
#include "../../src/UnitTesting/unit_test.h"

using omni::chemistry::ChiralOrientation;
using omni::constants::ExceptionResponse;
using omni::symbols::amber_ancient_bioq;
using omni::data_types::int2;
using omni::data_types::char4;
using omni::diskutil::DrivePathType;
using omni::diskutil::getDrivePathType;
using omni::diskutil::osSeparator;
using omni::errors::rtWarn;
using omni::parse::polyNumericVector;
using omni::parse::NumberFormat;
using omni::parse::operator!=;
using omni::parse::operator==;
using omni::topology::AtomGraph;
using omni::topology::ChemicalDetailsKit;
using omni::trajectory::CoordinateFileKind;
using omni::trajectory::PhaseSpace;

using namespace omni::chemistry;
using namespace omni::testing;

int main(int argc, char* argv[]) {

  // Lay out the testing environment
  TestEnvironment oe(argc, argv);
  StopWatch timer;

  // Test the existence of topology and coordinate files
  const char osc = osSeparator();
  const std::string base_chem_name = oe.getOmniSourcePath() + osc + "test" + osc + "Chemistry";
  const std::string base_nml_name  = oe.getOmniSourcePath() + osc + "test" + osc + "Namelists";
  const std::string base_nmlp_name = base_nml_name + osc + "topol";
  const std::string base_nmlc_name = base_nml_name + osc + "coord";
  const std::string base_top_name  = oe.getOmniSourcePath() + osc + "test" + osc + "Topology";
  const std::string base_crd_name  = oe.getOmniSourcePath() + osc + "test" + osc + "Trajectory";

  // Ligands with an eight-membered ring that will avoid anti-aromaticity by puckering
  const std::string mol1_top_name  = base_chem_name + osc + "lig1_c8h8.top";
  const std::string mol1_crd_name  = base_chem_name + osc + "lig1_c8h8.inpcrd";
  const std::string mol2_top_name  = base_chem_name + osc + "lig2_c8h8.top";
  const std::string mol2_crd_name  = base_chem_name + osc + "lig2_c8h8.inpcrd";

  // Additional ligands, with and without water
  const std::string mol3_top_name    = base_chem_name + osc + "lig3_1fsj.top";
  const std::string mol3_crd_name    = base_chem_name + osc + "lig3_1fsj.inpcrd";
  const std::string drug_top_name    = base_top_name + osc + "drug_example.top";
  const std::string drug_crd_name    = base_crd_name + osc + "drug_example.inpcrd";
  const std::string drug_vs_top_name = base_top_name + osc + "drug_example_vs.top";
  const std::string drug_vs_crd_name = base_crd_name + osc + "drug_example_vs.inpcrd";

  // Small peptides and proteins
  const std::string ala_top_name     = base_nmlp_name + osc + "ala.top";
  const std::string ala_crd_name     = base_nmlc_name + osc + "ala.inpcrd";
  const std::string gly_top_name     = base_nmlp_name + osc + "gly_gly.top";
  const std::string gly_crd_name     = base_nmlc_name + osc + "gly_gly.inpcrd";
  const std::string phe_top_name     = base_nmlp_name + osc + "phe.top";
  const std::string phe_crd_name     = base_nmlc_name + osc + "phe.inpcrd";
  const std::string trpcage_top_name = base_top_name + osc + "trpcage_in_water.top";
  const std::string trpcage_crd_name = base_crd_name + osc + "trpcage_in_water.inpcrd";
  const std::string ubiquit_top_name = base_top_name + osc + "ubiquitin.top";
  const std::string ubiquit_crd_name = base_crd_name + osc + "ubiquitin.inpcrd";

  // Water boxes
  const std::string tip5p_top_name = base_top_name + osc + "tip5p.top";
  const std::string tip5p_crd_name = base_crd_name + osc + "tip5p.rst";
  
  // Check the existence of all files
  const std::vector<std::string> top_files = { mol1_top_name, mol2_top_name, mol3_top_name,
                                               drug_top_name, drug_vs_top_name, ala_top_name,
                                               gly_top_name, phe_top_name, trpcage_top_name,
                                               ubiquit_top_name, tip5p_top_name };
  const std::vector<std::string> crd_files = { mol1_crd_name, mol2_crd_name, mol3_crd_name,
                                               drug_crd_name, drug_vs_crd_name, ala_crd_name,
                                               gly_crd_name, phe_crd_name, trpcage_crd_name,
                                               ubiquit_crd_name, tip5p_crd_name };
  bool files_exist = true;
  const size_t nsys = top_files.size();
  for (size_t i = 0; i < nsys; i++) {
    files_exist = (files_exist && getDrivePathType(top_files[i]) == DrivePathType::FILE);
    files_exist = (files_exist && getDrivePathType(crd_files[i]) == DrivePathType::FILE);
  }
  const TestPriority do_tests = (files_exist) ? TestPriority::CRITICAL : TestPriority::ABORT;
  if (files_exist == false) {
    rtWarn("Some required files were not found.  Check the ${OMNI_SOURCE} environment variable, "
           "currently set to " + oe.getOmniSourcePath() + ", for validity.  Subsequent tests will "
           "be skipped.", "test_chemical_features");
  }
  std::vector<AtomGraph> sys_ag(nsys);
  std::vector<PhaseSpace> sys_ps(nsys);
  if (files_exist) {
    for (size_t i = 0; i < nsys; i++) {
      sys_ag[i].buildFromPrmtop(top_files[i], ExceptionResponse::SILENT, amber_ancient_bioq,
                                1.2, 2.0, 0.01);
      sys_ps[i].buildFromFile(crd_files[i], CoordinateFileKind::UNKNOWN);
    }
  }
  std::vector<int> first_mol_size(nsys);
  std::vector<ChemicalFeatures> sys_chem;
  const MapRotatableGroups mapg_yes = MapRotatableGroups::YES;
  const MapRotatableGroups mapg_no  = MapRotatableGroups::NO;
  if (files_exist) {
    sys_chem.reserve(nsys);
    for (size_t i = 0; i < nsys; i++) {
      const int2 first_mol_lims = sys_ag[i].getMoleculeLimits(0);
      first_mol_size[i] = first_mol_lims.y - first_mol_lims.x;
      sys_chem.emplace_back(&sys_ag[i], CoordinateFrameReader(sys_ps[i]),
                            (first_mol_size[i] < 120) ? mapg_yes : mapg_no, 300.0, &timer);
    }
  }
  else {
    sys_chem.resize(nsys);  
  }
  std::vector<int> ring_counts(nsys);
  std::vector<int> fused_ring_counts(nsys);
  std::vector<int> mutable_ring_counts(nsys);
  std::vector<int> aromatic_group_counts(nsys);  
  std::vector<int> polar_h_counts(nsys);
  std::vector<int> hbond_donor_counts(nsys);
  std::vector<int> hbond_acceptor_counts(nsys);
  std::vector<int> chiral_center_counts(nsys);
  std::vector<int> rotatable_bond_counts(nsys);
  std::vector<int> cis_trans_bond_counts(nsys);
  for (size_t i = 0; i < nsys; i++) {
    ring_counts[i] = sys_chem[i].getRingCount();
    fused_ring_counts[i] = sys_chem[i].getFusedRingCount();
    mutable_ring_counts[i] = sys_chem[i].getRingCount();
    aromatic_group_counts[i] = sys_chem[i].getAromaticGroupCount();
    polar_h_counts[i] = sys_chem[i].getPolarHydrogenCount();
    hbond_donor_counts[i] = sys_chem[i].getHydrogenBondDonorCount();
    hbond_acceptor_counts[i] = sys_chem[i].getHydrogenBondAcceptorCount();
    chiral_center_counts[i] = sys_chem[i].getChiralCenterCount();
    rotatable_bond_counts[i] = sys_chem[i].getRotatableBondCount();
    cis_trans_bond_counts[i] = sys_chem[i].getCisTransBondCount();
  }
  std::vector<int> ring_cnt_ans           = {    1,    1,    4, 1226, 1226,    0,    0,    1,
                                              1567,  486,  216 };
  std::vector<int> fused_ring_cnt_ans     = {    0,    0,    0,    0,    0,    0,    0,    0,
                                                 1,    0,    0 };
  std::vector<int> mutable_ring_cnt_ans   = {    1,    1,    4, 1226, 1226,    0,    0,    1,
                                              1567,  486,  216 };
  std::vector<int> aromatic_group_cnt_ans = {    1,    1,    3,    1,    1,    0,    0,    1,
                                                 2,    3,    0 };
  std::vector<int> polar_h_cnt_ans        = {    0,    0,    6, 2450, 2450,    2,    3,    2,
                                              3155, 1088,  432 };
  std::vector<int> hbond_donor_cnt_ans    = {    0,    0,    4, 1226, 1226,    2,    3,    2,
                                              1587,  580,  216 };
  std::vector<int> hbond_acceptor_cnt_ans = {    0,    1,    7, 1228, 1228,    4,    6,    4,
                                              1614,  681,  216 };
  std::vector<int> chiral_center_cnt_ans  = {    0,    0,    0,    1,    1,    1,    0,    1,
                                                18,   82,    0 };
  std::vector<int> rotatable_bond_cnt_ans = {    1,    3,    9,    8,    8,    4,    7,    6,
                                                 0,    0,    0 };
  std::vector<int> cis_trans_bond_cnt_ans = {    0,    0,    0,    0,    0,    0,    0,    0,
                                                 0,    0,    0 };
  std::vector<int> trp_cage_lchir = sys_chem[8].getChiralCenters(ChiralOrientation::SINISTER);
  std::vector<int> trp_cage_lchir_ans;
  const ChemicalDetailsKit trp_cdk = sys_ag[8].getChemicalDetailsKit();
  for (int i = 0; i < trp_cdk.natom; i++) {
    if ((trp_cdk.atom_names[i] == char4({'C', 'A', ' ', ' '}) &&
         trp_cdk.res_names[sys_ag[8].getResidueIndex(i)] != char4({'G', 'L', 'Y', ' '})) ||
        (trp_cdk.atom_names[i] == char4({'C', 'B', ' ', ' '}) &&
         trp_cdk.res_names[sys_ag[8].getResidueIndex(i)] == char4({'I', 'L', 'E', ' '}))) {
      trp_cage_lchir_ans.push_back(i);
    }
  }
  check(ring_counts, RelationalOperator::EQUAL, ring_cnt_ans, "Overall counts of ring systems do "
        "not meet expectations.", do_tests);
  check(fused_ring_counts, RelationalOperator::EQUAL, fused_ring_cnt_ans, "Counts of fused ring "
        "systems do not meet expectations.", do_tests);
  check(mutable_ring_counts, RelationalOperator::EQUAL, mutable_ring_cnt_ans, "Counts of mutable "
        "ring systems do not meet expectations.", do_tests);
  check(aromatic_group_counts, RelationalOperator::EQUAL, aromatic_group_cnt_ans, "Counts of "
        "aromatic ring systems do not meet expectations.", do_tests);
  check(polar_h_counts, RelationalOperator::EQUAL, polar_h_cnt_ans, "Counts of polar hydrogens "
        "do not meet expectations.", do_tests);
  check(hbond_donor_counts, RelationalOperator::EQUAL, hbond_donor_cnt_ans, "Counts of hydrogen "
        "bond donors do not meet expectations.", do_tests);
  check(hbond_acceptor_counts, RelationalOperator::EQUAL, hbond_acceptor_cnt_ans, "Counts of "
        "hydrogen bond acceptors do not meet expectations.", do_tests);
  check(chiral_center_counts, RelationalOperator::EQUAL, chiral_center_cnt_ans, "Counts of chiral "
        "centers in various systems do not agree.", do_tests);
  check(rotatable_bond_counts, RelationalOperator::EQUAL, rotatable_bond_cnt_ans, "Counts of "
        "rotatable bonds do not meet expectations.", do_tests);
  check(cis_trans_bond_counts, RelationalOperator::EQUAL, cis_trans_bond_cnt_ans, "Counts of "
        "cis-trans bonds do not meet expectations.", do_tests);
  check(trp_cage_lchir_ans, RelationalOperator::EQUAL, trp_cage_lchir, "Chiral center indices for "
        "L-chiral centers (should be amino acid CA atoms, excluding glycine, plus the isoleucine "
        "CB atom) do not meet expectations.", do_tests);
  const std::string fc_name = base_chem_name + osc + "formal_charges.m";
  const std::string bo_name = base_chem_name + osc + "bond_orders.m";
  const std::string ro_name = base_chem_name + osc + "rotating_groups.m";
  const std::string ch_name = base_chem_name + osc + "chiral_atoms.m";
  const bool snps_exist = (getDrivePathType(fc_name) == DrivePathType::FILE &&
                           getDrivePathType(bo_name) == DrivePathType::FILE &&
                           getDrivePathType(ro_name) == DrivePathType::FILE &&
                           getDrivePathType(ch_name) == DrivePathType::FILE);
  const TestPriority do_snps = (snps_exist) ? TestPriority::CRITICAL : TestPriority::ABORT;
  if (snps_exist == false) {
    rtWarn("Snapshot files " + fc_name + ", " + bo_name + ", and " + ro_name + " must be "
           "accessible in order to check formal charge and bond order calculations, "
           "respectively.  Check the ${OMNI_SOURCE} environment variable for validity.  "
           "Subsequent tests will be skipped.", "test_chemical_features");
  }
  bool ch_unwritten = true;
  for (size_t i = 0; i < nsys; i++) {
    snapshot(fc_name, polyNumericVector(sys_chem[i].getFormalCharges()), std::string("fc_") +
             std::to_string(i), 1.0e-6, "Formal charges computed for the system described by " +
             sys_ag[i].getFileName() + " do not meet expectations.", oe.takeSnapshot(), 1.0e-8,
             NumberFormat::STANDARD_REAL,
             (i == 0LLU) ? PrintSituation::OVERWRITE : PrintSituation::APPEND, do_snps);
    snapshot(bo_name, polyNumericVector(sys_chem[i].getBondOrders()), std::string("bo_") +
             std::to_string(i), 1.0e-6, "Bond orders computed for the system described by " +
             sys_ag[i].getFileName() + " do not meet expectations.", oe.takeSnapshot(), 1.0e-8,
             NumberFormat::STANDARD_REAL,
             (i == 0LLU) ? PrintSituation::OVERWRITE : PrintSituation::APPEND, do_snps);
    if (sys_chem[i].getChiralCenterCount() > 0) {
      snapshot(ch_name, polyNumericVector(sys_chem[i].getChiralCenters()), std::string("chcen_") +
               std::to_string(i), 1.0e-6, "Chiral centers detected for the system described by " +
               sys_ag[i].getFileName() + " do not meet expectations.", oe.takeSnapshot(), 1.0e-8,
               NumberFormat::INTEGER,
               (ch_unwritten) ? PrintSituation::OVERWRITE : PrintSituation::APPEND, do_snps);
      ch_unwritten = false;
    }
    check(sum<double>(sys_chem[i].getFormalCharges()), RelationalOperator::EQUAL,
          Approx(sum<double>(sys_ag[i].getPartialCharge<double>())).margin(1.0e-4), "The sum of "
          "formal charges computed for " + sys_ag[i].getFileName() + " does not match the sum of "
          "partial charges given in the topology.", do_tests);
  }
  
  // Check the rotatable bond groups, and the inversion groups generated for smaller structures
  bool ro_unwritten = true;
  for (size_t i = 0; i < nsys; i++) {
    const int nrotor = sys_chem[i].getRotatableBondCount();
    if (nrotor > 0) {
      const std::vector<RotatorGroup> all_rotors = sys_chem[i].getRotatableBondGroups();
      std::vector<int> rotor_ids;
      for (int j = 0; j < nrotor; j++) {
        rotor_ids.push_back(all_rotors[j].root_atom);
        rotor_ids.push_back(all_rotors[j].pivot_atom);
        for (size_t k = 0; k < all_rotors[j].rotatable_atoms.size(); k++) {
          rotor_ids.push_back(all_rotors[j].rotatable_atoms[k]);
        }
      }
      snapshot(ro_name, polyNumericVector(rotor_ids), std::string("rotators_") + std::to_string(i),
               1.0e-6, "Rotatable atom IDs obtained for the system described by " +
               sys_ag[i].getFileName() + " do not meet expectations.", oe.takeSnapshot(), 1.0e-8,
               NumberFormat::INTEGER,
               (ro_unwritten) ? PrintSituation::OVERWRITE : PrintSituation::APPEND, do_snps);
      ro_unwritten = false;
    }
    const int nchiral = sys_chem[i].getChiralCenterCount();
    if (nchiral > 0) {
      const std::vector<RotatorGroup> all_invertors = sys_chem[i].getChiralInversionGroups();
      std::vector<int> invertor_ids;
      for (int j = 0; j < nchiral; j++) {
        invertor_ids.push_back(all_invertors[j].root_atom);
        invertor_ids.push_back(all_invertors[j].pivot_atom);
        for (size_t k = 0; k < all_invertors[j].rotatable_atoms.size(); k++) {
          invertor_ids.push_back(all_invertors[j].rotatable_atoms[k]);
        }
      }
      snapshot(ro_name, polyNumericVector(invertor_ids), std::string("invertors_") +
               std::to_string(i), 1.0e-6, "Chiral inversion group atom IDs obtained for the "
               "system described by " + sys_ag[i].getFileName() + " do not meet expectations.",
               oe.takeSnapshot(), 1.0e-8, NumberFormat::INTEGER,
               (ro_unwritten) ? PrintSituation::OVERWRITE : PrintSituation::APPEND, do_snps);
      ro_unwritten = false;
    }    
  }

  // Ring detection in some polycyclic systems
  const std::vector<std::string> ring_top_names = {
    base_chem_name + osc + "anthracene_like.top",
    base_chem_name + osc + "morphine_like.top",
    base_chem_name + osc + "pyrene_like.top"
  };
  const std::vector<std::string> ring_crd_names = {
    base_chem_name + osc + "anthracene_like.inpcrd",
    base_chem_name + osc + "morphine_like.inpcrd",
    base_chem_name + osc + "pyrene_like.inpcrd"
  };
  const size_t nring_mols = ring_top_names.size();
  bool rings_exist = true;
  for (size_t i = 0; i < nring_mols; i++) {
    rings_exist = (rings_exist && getDrivePathType(ring_top_names[i]) == DrivePathType::FILE &&
                   getDrivePathType(ring_top_names[i]) == DrivePathType::FILE);
  }
  if (rings_exist == false) {
    rtWarn("Topologies and coordinates for polycyclic test molecules, i.e. " + ring_top_names[0] +
           " and " + ring_crd_names[0] + ", were not found.  Check the ${OMNI_SOURCE} environment "
           "variable for validity.  Subsequent tests will be skipped.", "test_chemical_features");
  }
  std::vector<AtomGraph> ring_ag(nring_mols);
  std::vector<PhaseSpace> ring_ps(nring_mols);
  for (size_t i = 0; i < nring_mols; i++) {
    ring_ag[i].buildFromPrmtop(ring_top_names[i], ExceptionResponse::SILENT, amber_ancient_bioq,
                               1.2, 2.0, 0.01);
    ring_ps[i].buildFromFile(ring_crd_names[i], CoordinateFileKind::UNKNOWN);    
  }
  std::vector<ChemicalFeatures> ring_chem;
  if (rings_exist) {
    ring_chem.reserve(nring_mols);
    for (size_t i = 0; i < nring_mols; i++) {
      const int2 first_mol_lims = ring_ag[i].getMoleculeLimits(0);
      first_mol_size[i] = first_mol_lims.y - first_mol_lims.x;
      ring_chem.emplace_back(&ring_ag[i], CoordinateFrameReader(ring_ps[i]), mapg_yes, 300.0,
                             &timer);
    }
  }
  else {
    ring_chem.resize(nring_mols);
  }
  const TestPriority do_rings = (rings_exist) ? TestPriority::CRITICAL : TestPriority::ABORT;
  std::vector<int> polycyclic_ring_counts(nring_mols);
  std::vector<int> polycyclic_rotator_counts(nring_mols);
  std::vector<int> polycyclic_cistrans_counts(nring_mols);
  for (size_t i = 0; i < nring_mols; i++) {
    polycyclic_ring_counts[i] = ring_chem[i].getRingCount();
    polycyclic_rotator_counts[i] = ring_chem[i].getRotatableBondCount();
    polycyclic_cistrans_counts[i] = ring_chem[i].getCisTransBondCount();
  }
  const std::vector<int> polycyclic_ring_counts_answer = { 4, 5, 8 };
  const std::vector<int> polycyclic_rotator_counts_answer = { 0, 0, 0 };
  const std::vector<int> polycyclic_cistrans_counts_answer = { 0, 0, 0 };
  check(polycyclic_ring_counts, RelationalOperator::EQUAL, polycyclic_ring_counts_answer,
        "The number of rings detected in polycyclic systems does not meet expectations.",
        do_rings);
  check(polycyclic_rotator_counts, RelationalOperator::EQUAL, polycyclic_rotator_counts_answer,
        "The number of rotatable bonds detected in polycyclic systems does not meet "
        "expectations. (There should be no rotatable bonds.)", do_rings);
  check(polycyclic_cistrans_counts, RelationalOperator::EQUAL, polycyclic_cistrans_counts_answer,
        "The number of cis-trans isomerization sites detected in polycyclic systems does not meet "
        "expectations. (There should be no bonds capable of cis-trans isomerization.)", do_rings);
  
  // Summary evaluation
  if (oe.getDisplayTimingsOrder()) {
    timer.printResults();
  }
  printTestSummary(oe.getVerbosity());
}
