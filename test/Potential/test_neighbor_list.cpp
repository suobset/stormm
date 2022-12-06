#include "copyright.h"
#include "../../src/Constants/behavior.h"
#include "../../src/DataTypes/stormm_vector_types.h"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/Parsing/parse.h"
#include "../../src/Potential/forward_exclusionmask.h"
#include "../../src/Potential/valence_potential.h"
#include "../../src/Potential/scorecard.h"
#include "../../src/Reporting/error_format.h"
#include "../../src/Reporting/summary_file.h"
#include "../../src/Topology/atomgraph.h"
#include "../../src/Topology/atomgraph_abstracts.h"
#include "../../src/Trajectory/phasespace.h"
#include "../../src/UnitTesting/unit_test.h"

#ifndef STORMM_USE_HPC
using stormm::int2;
#endif
using stormm::constants::ExceptionResponse;
using stormm::diskutil::DrivePathType;
using stormm::diskutil::getBaseName;
using stormm::diskutil::getDrivePathType;
using stormm::diskutil::osSeparator;
using stormm::errors::rtWarn;
using stormm::parse::NumberFormat;
using stormm::parse::polyNumericVector;
using stormm::review::stormmSplash;
using stormm::topology::AtomGraph;
using stormm::topology::NonbondedKit;
using stormm::trajectory::CoordinateFileKind;
using stormm::trajectory::PhaseSpace;
using stormm::trajectory::TrajectoryKind;
using namespace stormm::energy;
using namespace stormm::testing;

//-------------------------------------------------------------------------------------------------
// Loop over all exclusions in a given array based on the supplied bounds.  Check for non-reflexive
// exclusions and also for self-exclusions (neither of which should exist).  Modify formal
// arguments if either condition is found.  Return a vector of booleans to indicate that exclusions
// listed in the topology were found within a compact "forward" exclusion mask.
//
// Arguments
//   natom:                     Number of atoms in the system
//   bounds:                    Bounds array for exclusions in excl_list
//   excl_list:                 Exclusions list, rather like EXCLUDED_ATOMS_LIST in an Amber prmtop
//                              but with double-counting for all exclusions
//   femask:                    A very compact expression of the system's masks
//   self_exclusion_detected:   Indicator of self exclusions detected in error (returned)
//   non_reflexive_exclusions:  Indicator of non-reflexive elements of excl_list (returned)
//-------------------------------------------------------------------------------------------------
std::vector<bool> checkForwardMaskExcl(const int natom, const int* bounds, const int* excl_list,
                                       const ForwardExclusionMask &femask,
                                       bool *self_exclusion_detected,
                                       std::vector<int2> *non_reflexive_exclusions) {
  std::vector<bool> result(bounds[natom], false);
  for (int i = 0; i < natom; i++) {
    for (int j = bounds[i]; j < bounds[i + 1]; j++) {
      const int j_atom_idx = excl_list[j];
      *self_exclusion_detected = (*self_exclusion_detected || j_atom_idx == i);
      bool mirror_found = false;
      for (int k = bounds[j_atom_idx]; k < bounds[j_atom_idx + 1]; k++) {
        mirror_found = (mirror_found || excl_list[k] == i);
      }
      if (mirror_found == false) {
        non_reflexive_exclusions->push_back({i, j_atom_idx});
      }
      result[j] = femask.testExclusion(i, j_atom_idx);
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
// main
//-------------------------------------------------------------------------------------------------
int main(const int argc, const char* argv[]) {

  // Some baseline initialization
  TestEnvironment oe(argc, argv);
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmSplash();
  }

  // Section 1
  section("Forward exclusion list analysis");

  // Locate topologies and coordinate files
  const char osc = osSeparator();
  const std::string base_top_name = oe.getStormmSourcePath() + osc + "test" + osc + "Topology";
  const std::string base_crd_name = oe.getStormmSourcePath() + osc + "test" + osc + "Trajectory";
  const std::string base_ptl_name = oe.getStormmSourcePath() + osc + "test" + osc + "Potential";
  const std::string trpi_top_name = base_top_name + osc + "trpcage.top";
  const std::string trpi_crd_name = base_crd_name + osc + "trpcage.inpcrd";
  const std::string trpw_top_name = base_top_name + osc + "trpcage_in_water.top";
  const std::string trpw_crd_name = base_crd_name + osc + "trpcage_in_water.inpcrd";
  const std::string dhfr_top_name = base_top_name + osc + "dhfr_cmap.top";
  const std::string dhfr_crd_name = base_crd_name + osc + "dhfr_cmap.inpcrd";
  const std::string alad_top_name = base_top_name + osc + "ala_dipeptide.top";
  const std::string alad_crd_name = base_crd_name + osc + "ala_dipeptide.inpcrd";
  const std::string ahec_top_name = base_top_name + osc + "ahector.top";
  const std::string ahec_crd_name = base_crd_name + osc + "ahector.inpcrd";
  const bool systems_exist = (getDrivePathType(trpi_top_name) == DrivePathType::FILE &&
                              getDrivePathType(trpi_crd_name) == DrivePathType::FILE &&
                              getDrivePathType(trpw_top_name) == DrivePathType::FILE &&
                              getDrivePathType(trpw_crd_name) == DrivePathType::FILE &&
                              getDrivePathType(dhfr_top_name) == DrivePathType::FILE &&
                              getDrivePathType(dhfr_crd_name) == DrivePathType::FILE &&
                              getDrivePathType(alad_top_name) == DrivePathType::FILE &&
                              getDrivePathType(alad_crd_name) == DrivePathType::FILE &&
                              getDrivePathType(ahec_top_name) == DrivePathType::FILE &&
                              getDrivePathType(ahec_crd_name) == DrivePathType::FILE);
  const TestPriority do_tests = (systems_exist) ? TestPriority::CRITICAL : TestPriority::ABORT;
  if (systems_exist == false) {
    rtWarn("Files for the Trp-cage miniprotein, the DHFR globular protein (with CHARMM potential "
           "details), alanine dipeptide, and scorpion toxin were not found.  These files should "
           "befound in the ${STORMM_SOURCE}/test/Topology and ${STORMM_SOURCE}/test/Trajectory "
           "directories.  Check the $STORMM_SOURCE environment variable.  A number of tests will "
           "be skipped.", "test_neighbor_list");
  }
  AtomGraph trpi_ag, trpw_ag, dhfr_ag, alad_ag, ahec_ag;
  PhaseSpace trpi_ps, trpw_ps, dhfr_ps, alad_ps, ahec_ps;
  if (systems_exist) {
    trpi_ag.buildFromPrmtop(trpi_top_name, ExceptionResponse::SILENT);
    trpi_ps.buildFromFile(trpi_crd_name, CoordinateFileKind::AMBER_INPCRD);
    trpw_ag.buildFromPrmtop(trpw_top_name, ExceptionResponse::SILENT);
    trpw_ps.buildFromFile(trpw_crd_name, CoordinateFileKind::AMBER_INPCRD);
    dhfr_ag.buildFromPrmtop(dhfr_top_name, ExceptionResponse::SILENT);
    dhfr_ps.buildFromFile(dhfr_crd_name, CoordinateFileKind::AMBER_INPCRD);
    alad_ag.buildFromPrmtop(alad_top_name, ExceptionResponse::SILENT);
    alad_ps.buildFromFile(alad_crd_name, CoordinateFileKind::AMBER_INPCRD);
    ahec_ag.buildFromPrmtop(ahec_top_name, ExceptionResponse::SILENT);
    ahec_ps.buildFromFile(ahec_crd_name, CoordinateFileKind::AMBER_INPCRD);
  }

  // Set up the forward exclusion lists and check each topology's exclusions
  section(1);
  const std::vector<AtomGraph*> all_topologies = { &trpi_ag, &trpw_ag, &dhfr_ag, &alad_ag,
                                                   &ahec_ag };
  const std::vector<ForwardExclusionMask> all_femasks = { ForwardExclusionMask(&trpi_ag),
                                                          ForwardExclusionMask(&trpw_ag),
                                                          ForwardExclusionMask(&dhfr_ag),
                                                          ForwardExclusionMask(&alad_ag),
                                                          ForwardExclusionMask(&ahec_ag) };
  const std::vector<int> primary_mask_count_answer = { 99, 99, 140, 11, 145 };
  const std::vector<int> extended_mask_count_answer = { 0, 0, 0, 0, 20 };
  const std::vector<int> secondary_mask_count_answer = { 0, 0, 0, 0, 0 };
  std::vector<int> primary_mask_counts(all_topologies.size());
  std::vector<int> extended_mask_counts(all_topologies.size());
  std::vector<int> secondary_mask_counts(all_topologies.size());

  // Scan the exclusion lists in each topology, making sure that the forward masks cover all
  // of them.  Also check the topology exclusions themselves, in one additional context, to
  // make sure that they contain no self exclusions and that every exclusion is reflexive.
  for (size_t i = 0; i < all_topologies.size(); i++) {
    bool self_exclusion_detected = false;
    std::vector<int2> non_reflexive_exclusions;
    const NonbondedKit<double> nbk = all_topologies[i]->getDoublePrecisionNonbondedKit();
    const std::vector<bool> nb11x_satisfied =
      checkForwardMaskExcl(nbk.natom, nbk.nb11_bounds, nbk.nb11x, all_femasks[i],
                           &self_exclusion_detected, &non_reflexive_exclusions);
    const std::vector<bool> nb12x_satisfied =
      checkForwardMaskExcl(nbk.natom, nbk.nb12_bounds, nbk.nb12x, all_femasks[i],
                           &self_exclusion_detected, &non_reflexive_exclusions);
    const std::vector<bool> nb13x_satisfied =
      checkForwardMaskExcl(nbk.natom, nbk.nb13_bounds, nbk.nb13x, all_femasks[i],
                           &self_exclusion_detected, &non_reflexive_exclusions);
    const std::vector<bool> nb14x_satisfied =
      checkForwardMaskExcl(nbk.natom, nbk.nb14_bounds, nbk.nb14x, all_femasks[i],
                           &self_exclusion_detected, &non_reflexive_exclusions);
    check(self_exclusion_detected == false, "Self exclusions were detected in topology " +
          all_topologies[i]->getFileName() + ".", do_tests);
    check(non_reflexive_exclusions.size() == 0LLU, "Non-reflexive exclusions were detected in "
          "topology " + all_topologies[i]->getFileName() + ".", do_tests);
    if (nb11x_satisfied.size() > 0LLU) {
      check(nb11x_satisfied, RelationalOperator::EQUAL,
            std::vector<bool>(nbk.nb11_bounds[nbk.natom], true), "Some 1:1 exclusions in "
            "topology " + all_topologies[i]->getFileName() + " were not satisfied by the forward "
            "mask.", do_tests);
    }
    if (nb12x_satisfied.size() > 0LLU) {
      check(nb12x_satisfied, RelationalOperator::EQUAL,
            std::vector<bool>(nbk.nb12_bounds[nbk.natom], true), "Some 1:2 exclusions in "
            "topology " + all_topologies[i]->getFileName() + " were not satisfied by the forward "
            "mask.", do_tests);
    }
    if (nb13x_satisfied.size() > 0LLU) {
      check(nb13x_satisfied, RelationalOperator::EQUAL,
            std::vector<bool>(nbk.nb13_bounds[nbk.natom], true), "Some 1:3 exclusions in "
            "topology " + all_topologies[i]->getFileName() + " were not satisfied by the forward "
            "mask.", do_tests);
    }
    if (nb14x_satisfied.size() > 0LLU) {
      check(nb14x_satisfied, RelationalOperator::EQUAL,
            std::vector<bool>(nbk.nb14_bounds[nbk.natom], true), "Some 1:4 exclusions in "
            "topology " + all_topologies[i]->getFileName() + " were not satisfied by the forward "
            "mask.", do_tests);
    }
    primary_mask_counts[i]   = all_femasks[i].getPrimaryMaskCount();
    extended_mask_counts[i]  = all_femasks[i].getExtendedMaskCount();
    secondary_mask_counts[i] = all_femasks[i].getSecondaryMaskCount();
    check((nbk.nb11_bounds[nbk.natom] + nbk.nb12_bounds[nbk.natom] + nbk.nb13_bounds[nbk.natom] +
           nbk.nb14_bounds[nbk.natom]) / 2, RelationalOperator::EQUAL,
          all_femasks[i].getTotalExclusionsCount(), "The total number of exclusions expressed in "
          "topology " + getBaseName(all_topologies[i]->getFileName()) + " is not refelcted in the "
          "corresponding ForwardExclusionsMask.", do_tests);
  }
  check(primary_mask_counts, RelationalOperator::EQUAL, primary_mask_count_answer, "Unique "
        "primary exclusion mask counts in each system do not match.", do_tests);
  check(extended_mask_counts, RelationalOperator::EQUAL, extended_mask_count_answer, "Unique "
        "extended exclusion mask counts in each system do not match.", do_tests);
  check(secondary_mask_counts, RelationalOperator::EQUAL, secondary_mask_count_answer, "Unique "
        "secondary exclusion mask counts in each system do not match.", do_tests);

  // Print results
  printTestSummary(oe.getVerbosity());

  return countGlobalTestFailures();
}
