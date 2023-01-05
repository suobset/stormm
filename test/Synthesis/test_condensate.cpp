#include "../../src/DataTypes/stormm_vector_types.h"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/Numerics/split_fixed_precision.h"
#include "../../src/Random/random.h"
#include "../../src/Reporting/summary_file.h"
#include "../../src/Structure/rmsd.h"
#include "../../src/Structure/rmsd_plan.h"
#include "../../src/Synthesis/condensate.h"
#include "../../src/Synthesis/phasespace_synthesis.h"
#include "../../src/Synthesis/synthesis_enumerators.h"
#include "../../src/Topology/atomgraph.h"
#include "../../src/Trajectory/phasespace.h"
#include "../../src/UnitTesting/test_system_manager.h"
#include "../../src/UnitTesting/unit_test.h"

using namespace stormm::data_types;
using namespace stormm::diskutil;
using namespace stormm::numerics;
using namespace stormm::random;
using namespace stormm::review;
using namespace stormm::structure;
using namespace stormm::synthesis;
using namespace stormm::testing;

//-------------------------------------------------------------------------------------------------
// main
//-------------------------------------------------------------------------------------------------
int main(const int argc, const char* argv[]) {

  // Baseline initialization
  const TestEnvironment oe(argc, argv);
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmSplash();
  }

  // Section 1: test the veracity of the Condensate object
  section("Verify Condensate contents");

  // Section 2: test RMSD calculations with the Condensate
  section("Test RMSD calcuations");
  
  // Read all systems
  const char osc = osSeparator();
  const std::string base_top_dir = oe.getStormmSourcePath() + osc + "test" + osc + "Topology";
  const std::string base_crd_dir = oe.getStormmSourcePath() + osc + "test" + osc + "Trajectory";
  TestSystemManager tsm(base_top_dir, "top",
                        { "symmetry_L1", "stereo_L1", "symmetry_L1_vs",
                          "stereo_L1_vs", "bromobenzene", "bromobenzene_vs" }, base_crd_dir,
                        "inpcrd", { "symmetry_L1", "stereo_L1", "symmetry_L1_vs",
                                    "stereo_L1_vs", "bromobenzene", "bromobenzene_vs" });

  // Form a synthesis of the coordinates.  Perturb each structure slightly.
  Xoroshiro128pGenerator xrs(517389207);
  std::vector<AtomGraph*> ag_list;
  std::vector<PhaseSpace> ps_list;
  const int n_structures = 100;
  ag_list.reserve(n_structures);
  ps_list.reserve(n_structures);
  for (int i = 0; i < n_structures; i++) {
    const int ag_id = static_cast<double>(tsm.getSystemCount()) * xrs.uniformRandomNumber();
    ag_list.push_back(tsm.getTopologyPointer(ag_id));
    ps_list.push_back(tsm.exportPhaseSpace(ag_id));
    PhaseSpaceWriter psw = ps_list.back().data();
    for (int j = 0; j < psw.natom; j++) {
      psw.xcrd[j] += 0.5 * (xrs.uniformRandomNumber() - 0.5);
      psw.ycrd[j] += 0.5 * (xrs.uniformRandomNumber() - 0.5);
      psw.zcrd[j] += 0.5 * (xrs.uniformRandomNumber() - 0.5);
    }
  }
  PhaseSpaceSynthesis poly_ps(ps_list, ag_list, 48);

  // Create a condensate from the synthesis.
  Condensate cdns_dbl(poly_ps, CondensationLevel::DOUBLE);
  Condensate cdns_flt(poly_ps, CondensationLevel::FLOAT);
  section(1);
  const int n_test_pt = 100;
  std::vector<double> xdbl_dev(n_test_pt), ydbl_dev(n_test_pt), zdbl_dev(n_test_pt);
  std::vector<double> xflt_dev(n_test_pt), yflt_dev(n_test_pt), zflt_dev(n_test_pt);
  CondensateWriter cdw_dbl = cdns_dbl.data();
  CondensateWriter cdw_flt = cdns_flt.data();
  PsSynthesisWriter poly_psw = poly_ps.data();
  for (int i = 0; i < n_test_pt; i++) {
    const int sys_idx = poly_psw.system_count * xrs.uniformRandomNumber();
    const int atom_idx = static_cast<int>(static_cast<double>(poly_psw.atom_counts[sys_idx]) *
                                          xrs.uniformRandomNumber()) +
                         poly_psw.atom_starts[sys_idx];
    xdbl_dev[i] = (hostInt95ToDouble(poly_psw.xcrd[atom_idx], poly_psw.xcrd_ovrf[atom_idx]) *
                   poly_psw.inv_gpos_scale) - cdw_dbl.xcrd[atom_idx];
    ydbl_dev[i] = (hostInt95ToDouble(poly_psw.ycrd[atom_idx], poly_psw.ycrd_ovrf[atom_idx]) *
                   poly_psw.inv_gpos_scale) - cdw_dbl.ycrd[atom_idx];
    zdbl_dev[i] = (hostInt95ToDouble(poly_psw.zcrd[atom_idx], poly_psw.zcrd_ovrf[atom_idx]) *
                   poly_psw.inv_gpos_scale) - cdw_dbl.zcrd[atom_idx];
    xflt_dev[i] = (hostInt95ToDouble(poly_psw.xcrd[atom_idx], poly_psw.xcrd_ovrf[atom_idx]) *
                   poly_psw.inv_gpos_scale) - cdw_flt.xcrd_sp[atom_idx];
    yflt_dev[i] = (hostInt95ToDouble(poly_psw.ycrd[atom_idx], poly_psw.ycrd_ovrf[atom_idx]) *
                   poly_psw.inv_gpos_scale) - cdw_flt.ycrd_sp[atom_idx];
    zflt_dev[i] = (hostInt95ToDouble(poly_psw.zcrd[atom_idx], poly_psw.zcrd_ovrf[atom_idx]) *
                   poly_psw.inv_gpos_scale) - cdw_flt.zcrd_sp[atom_idx];
  }
  check(xdbl_dev, RelationalOperator::EQUAL, std::vector<double>(n_test_pt, 0.0), "Cartesian X "
        "coordinates kept in the Condensate differ from the original PhaseSpaceSynthesis, even "
        "when recorded in double-precision.", tsm.getTestingStatus());
  check(ydbl_dev, RelationalOperator::EQUAL, std::vector<double>(n_test_pt, 0.0), "Cartesian Y "
        "coordinates kept in the Condensate differ from the original PhaseSpaceSynthesis, even "
        "when recorded in double-precision.", tsm.getTestingStatus());
  check(zdbl_dev, RelationalOperator::EQUAL, std::vector<double>(n_test_pt, 0.0), "Cartesian Z "
        "coordinates kept in the Condensate differ from the original PhaseSpaceSynthesis, even "
        "when recorded in double-precision.", tsm.getTestingStatus());
  check(xflt_dev, RelationalOperator::EQUAL,
        Approx(std::vector<double>(n_test_pt, 0.0)).margin(6.0e-6), "Cartesian X coordinates kept "
        "in the Condensate differ from the original PhaseSpaceSynthesis, even when recorded in "
        "double-precision.", tsm.getTestingStatus());
  check(yflt_dev, RelationalOperator::EQUAL,
        Approx(std::vector<double>(n_test_pt, 0.0)).margin(6.0e-6), "Cartesian Y coordinates kept "
        "in the Condensate differ from the original PhaseSpaceSynthesis, even when recorded in "
        "double-precision.", tsm.getTestingStatus());
  check(zflt_dev, RelationalOperator::EQUAL,
        Approx(std::vector<double>(n_test_pt, 0.0)).margin(6.0e-6), "Cartesian Z coordinates kept "
        "in the Condensate differ from the original PhaseSpaceSynthesis, even when recorded in "
        "double-precision.", tsm.getTestingStatus());
  check(cdw_dbl.natr_insr, RelationalOperator::EQUAL, tsm.getSystemCount(), "The number of "
        "all-to-reference instructions does not agree with the number of systems for a CPU-bound "
        "Condensate.", tsm.getTestingStatus());
  check(cdw_dbl.nata_insr, RelationalOperator::EQUAL, tsm.getSystemCount(), "The number of "
        "all-to-all instructions does not agree with the number of systems for a CPU-bound "
        "Condensate.", tsm.getTestingStatus());
  const std::vector<int> insr_counts_ans = { 18, 14, 19, 19, 14, 16 };
  std::vector<int> insr_counts(tsm.getSystemCount());
  for (int i = 0; i < cdw_dbl.natr_insr; i++) {
    insr_counts[i] = (cdw_dbl.ata_insr[i].w >> 16);
  }
  check(insr_counts, RelationalOperator::EQUAL, insr_counts_ans, "Calculation counts in each "
        "all-to-all work unit do not meet expectations.", tsm.getTestingStatus());

  // Try an RMSD calculation with the Condensate guiding the process on the CPU.  Compare it to
  // an RMSD calculation using the synthesis directly.
  RMSDPlan rplan(poly_ps);
  Hybrid<double> atr_rmsd_a(rplan.getReferenceRMSDSize());
  Hybrid<double> atr_rmsd_b(rplan.getReferenceRMSDSize());
  Hybrid<double> atr_rmsd_c(rplan.getReferenceRMSDSize());
  Hybrid<double> ata_rmsd_a(rplan.getRMSDMatrixSize());
  Hybrid<double> ata_rmsd_b(rplan.getRMSDMatrixSize());
  Hybrid<double> ata_rmsd_c(rplan.getRMSDMatrixSize());
  Hybrid<int> examples(poly_ps.getUniqueTopologyExampleIndices());
  rmsd(rplan, poly_ps, examples, &atr_rmsd_a);
  rmsd(rplan, poly_ps, &ata_rmsd_a);
  rmsd(rplan, poly_ps, cdns_dbl, examples, &atr_rmsd_b);
  rmsd(rplan, poly_ps, cdns_dbl, &ata_rmsd_b);
  rmsd(rplan, poly_ps, cdns_flt, examples, &atr_rmsd_c);
  rmsd(rplan, poly_ps, cdns_flt, &ata_rmsd_c);
  const std::vector<double> atr_rmsd_av = atr_rmsd_a.readHost();
  const std::vector<double> atr_rmsd_bv = atr_rmsd_b.readHost();
  const std::vector<double> atr_rmsd_cv = atr_rmsd_c.readHost();
  const std::vector<double> ata_rmsd_av = ata_rmsd_a.readHost();
  const std::vector<double> ata_rmsd_bv = ata_rmsd_b.readHost();
  const std::vector<double> ata_rmsd_cv = ata_rmsd_c.readHost();
  section(2);
  check(atr_rmsd_bv, RelationalOperator::EQUAL, atr_rmsd_av, "RMSD (all to reference) "
        "calculations guided by a Condensate object's instructions on the CPU do not match RMSD "
        "calculations performed on the original PhaseSpaceSynthesis.", tsm.getTestingStatus());
  check(ata_rmsd_bv, RelationalOperator::EQUAL, ata_rmsd_av, "RMSD (all to all) calculations "
        "guided by a Condensate object's instructions on the CPU do not match RMSD calculations "
        "performed on the original PhaseSpaceSynthesis.", tsm.getTestingStatus());
  check(ata_rmsd_cv, RelationalOperator::EQUAL, Approx(ata_rmsd_av).margin(1.5e-5), "RMSD (all to "
        "all) calculations guided by a Condensate object's instructions on the CPU do not match "
        "RMSD calculations performed on the original PhaseSpaceSynthesis when the Condensate is "
        "cast in 32 bit floating-point precision.", tsm.getTestingStatus());
  
  // Print results
  printTestSummary(oe.getVerbosity());

  return countGlobalTestFailures();
}
