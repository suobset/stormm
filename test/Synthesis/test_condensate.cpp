#include "../../src/Constants/behavior.h"
#include "../../src/Constants/scaling.h"
#include "../../src/DataTypes/common_types.h"
#include "../../src/DataTypes/stormm_vector_types.h"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/Math/matrix_ops.h"
#include "../../src/Math/vector_ops.h"
#include "../../src/Numerics/split_fixed_precision.h"
#include "../../src/Parsing/parse.h"
#include "../../src/Random/random.h"
#include "../../src/Reporting/summary_file.h"
#include "../../src/Structure/rmsd.h"
#include "../../src/Structure/rmsd_plan.h"
#include "../../src/Synthesis/condensate.h"
#include "../../src/Synthesis/phasespace_synthesis.h"
#include "../../src/Synthesis/synthesis_enumerators.h"
#include "../../src/Topology/atomgraph.h"
#include "../../src/Trajectory/coordinate_copy.h"
#include "../../src/Trajectory/trajectory_enumerators.h"
#include "../../src/Trajectory/phasespace.h"
#include "../../src/UnitTesting/test_system_manager.h"
#include "../../src/UnitTesting/unit_test.h"

using namespace stormm::constants;
using namespace stormm::data_types;
using namespace stormm::diskutil;
using namespace stormm::math;
using namespace stormm::numerics;
using namespace stormm::parse;
using namespace stormm::random;
using namespace stormm::review;
using namespace stormm::structure;
using namespace stormm::synthesis;
using namespace stormm::testing;
using namespace stormm::trajectory;

//-------------------------------------------------------------------------------------------------
// Quantify the difference in Cartesian X, Y, and Z coordinates.  Various overloads serve different
// coordinate formats.
//
// Arguments:
//   xcrd:      Experimental Cartesian X coordinates
//   ycrd:      Experimental Cartesian Y coordinates
//   zcrd:      Experimental Cartesian Z coordinates
//   xref:      Reference Cartesian X coordinates
//   yref:      Reference Cartesian Y coordinates
//   zref:      Reference Cartesian Z coordinates
//   natom:     The number of atoms to analyze
//   copy_msg:  Message to display, describing the type of copying that led to the experimental
//              coordinates, in the event that any inconsistencies are found
//   do_test:   Indication of whether the test is feasible, based on prior success in reading
//              critical coordiante files
//   overall:   The expected overall difference that should be found
//   bounds:    The expected overall difference that should be found
//-------------------------------------------------------------------------------------------------
template <typename Tcrd, typename Tref>
void diffCoordinates(const Tcrd* xcrd, const Tcrd* ycrd, const Tcrd* zcrd, const int* xcrd_ovrf,
                     const int* ycrd_ovrf, const int* zcrd_ovrf, const double crd_scale,
                     const Tref* xref, const Tref* yref, const Tref* zref, const int* xref_ovrf,
                     const int* yref_ovrf, const int* zref_ovrf, const double ref_scale,
                     const int natom, const double* umat, const double* umat_ref,
                     const std::string &copy_msg, const TestPriority do_test,
                     const double overall = tiny, const double bounds = tiny) {
  bool xoff = false;
  bool yoff = false;
  bool zoff = false;
  double rmsd = 0.0;
  std::vector<double> x_errors, y_errors, z_errors;
  for (int i = 0; i < natom; i++) {
    const double tx_crd = (xcrd_ovrf == nullptr) ?
                          static_cast<double>(xcrd[i]) / crd_scale :
                          hostInt95ToDouble(xcrd[i], xcrd_ovrf[i]) / crd_scale;
    const double ty_crd = (ycrd_ovrf == nullptr) ?
                          static_cast<double>(ycrd[i]) / crd_scale :
                          hostInt95ToDouble(ycrd[i], ycrd_ovrf[i]) / crd_scale;
    const double tz_crd = (zcrd_ovrf == nullptr) ?
                          static_cast<double>(zcrd[i]) / crd_scale :
                          hostInt95ToDouble(zcrd[i], zcrd_ovrf[i]) / crd_scale;
    const double tx_ref = (xref_ovrf == nullptr) ?
                          static_cast<double>(xref[i]) / ref_scale :
                          hostInt95ToDouble(xref[i], xref_ovrf[i]) / ref_scale;
    const double ty_ref = (yref_ovrf == nullptr) ?
                          static_cast<double>(yref[i]) / ref_scale :
                          hostInt95ToDouble(yref[i], yref_ovrf[i]) / ref_scale;
    const double tz_ref = (zref_ovrf == nullptr) ?
                          static_cast<double>(zref[i]) / ref_scale :
                          hostInt95ToDouble(zref[i], zref_ovrf[i]) / ref_scale;
    xoff = (xoff || fabs(tx_crd - tx_ref) > bounds);
    yoff = (yoff || fabs(ty_crd - ty_ref) > bounds);
    zoff = (zoff || fabs(tz_crd - tz_ref) > bounds);
    if (fabs(tx_crd - tx_ref) > bounds) {
      xoff = true;
      x_errors.push_back(fabs(tx_crd - tx_ref));
    }
    if (fabs(ty_crd - ty_ref) > bounds) {
      yoff = true;
      y_errors.push_back(fabs(ty_crd - ty_ref));
    }
    if (fabs(tz_crd - tz_ref) > bounds) {
      zoff = true;
      z_errors.push_back(fabs(tz_crd - tz_ref));
    }
    const double dx = tx_crd - tx_ref;
    const double dy = ty_crd - ty_ref;
    const double dz = tz_crd - tz_ref;
    rmsd += (dx * dx) + (dy * dy) + (dz * dz);
  }
  if (umat != nullptr && umat_ref != nullptr) {
    std::vector<double> umat_stl(9), umat_ref_stl(9);
    for (int i = 0; i < 9; i++) {
      umat_stl[i] = umat[i];
      umat_ref_stl[i] = umat_ref[i];
    }
    check(umat_stl, RelationalOperator::EQUAL, umat_ref_stl, "Transformation matrices for two "
          "coordinate sets disagree after copying " + copy_msg + ".", do_test);
  }
  rmsd = sqrt(rmsd / static_cast<double>(natom));
  std::string fail_msg;
  if (xoff) {
    fail_msg = "X (" + realToString(mean(x_errors), 11, 4, NumberFormat::SCIENTIFIC) + " +/- " +
               realToString(variance(x_errors, VarianceMethod::STANDARD_DEVIATION), 11, 4,
                            NumberFormat::SCIENTIFIC) + ")";
  }
  if (yoff) {
    fail_msg += (xoff) ? ", Y" : "Y";
    fail_msg += " (" + realToString(mean(y_errors), 11, 4, NumberFormat::SCIENTIFIC) + " +/- " +
                realToString(variance(y_errors, VarianceMethod::STANDARD_DEVIATION), 11, 4,
                             NumberFormat::SCIENTIFIC) + ")";
  }
  if (zoff) {
    fail_msg += (xoff || yoff) ? ", Z" : "Z";
    fail_msg += " (" + realToString(mean(z_errors), 11, 4, NumberFormat::SCIENTIFIC) + " +/- " +
                realToString(variance(z_errors, VarianceMethod::STANDARD_DEVIATION), 11, 4,
                             NumberFormat::SCIENTIFIC) + ")";
  }
  check((xoff || yoff || zoff) == false, "Coordinates exceeded the specified deviations after "
        "copying " + copy_msg + ".  Deviations occur in " + fail_msg + ".", do_test);
  if (overall > 1.5 * tiny) {
    check(rmsd, RelationalOperator::EQUAL, Approx(overall, ComparisonType::RELATIVE, 0.01),
          "The root mean squared coordinate deviation obtained after copying " + copy_msg +
          " does not meet expectations.", do_test);
  }
  else {
    check(rmsd, RelationalOperator::EQUAL, Approx(0.0).margin(overall), "The root mean squared "
          "coordinate deviation obtained after copying " + copy_msg + " does not meet "
          "expectations.", do_test);
  }
}

template <typename Tcrd, typename Tref>
void diffCoordinates(const Tcrd* xcrd, const Tcrd* ycrd, const Tcrd* zcrd, const int* xcrd_ovrf,
                     const int* ycrd_ovrf, const int* zcrd_ovrf, const double crd_scale,
                     const Tref* xref, const Tref* yref, const Tref* zref, const int* xref_ovrf,
                     const int* yref_ovrf, const int* zref_ovrf, const double ref_scale,
                     const int natom, const std::string &copy_msg, const TestPriority do_test,
                     const double overall = tiny, const double bounds = tiny) {
  diffCoordinates(xcrd, ycrd, zcrd, xcrd_ovrf, ycrd_ovrf, zcrd_ovrf, crd_scale, xref, yref, zref,
                  xref_ovrf, yref_ovrf, zref_ovrf, ref_scale, natom, nullptr, nullptr, copy_msg,
                  do_test, overall, bounds);
}

template <typename Tcrd, typename Tref>
void diffCoordinates(const Tcrd* xcrd, const Tcrd* ycrd, const Tcrd* zcrd, const double crd_scale,
                     const Tref* xref, const Tref* yref, const Tref* zref, const double ref_scale,
                     const int natom, const std::string &copy_msg, const TestPriority do_test,
                     const double overall = 0.0, const double bounds = tiny) {
  diffCoordinates(xcrd, ycrd, zcrd, nullptr, nullptr, nullptr, crd_scale, xref, yref, zref,
                  nullptr, nullptr, nullptr, ref_scale, natom, nullptr, nullptr, copy_msg, do_test,
                  overall, bounds);
}

template <typename Tcrd, typename Tref>
void diffCoordinates(const Tcrd* xcrd, const Tcrd* ycrd, const Tcrd* zcrd, const Tref* xref,
                     const Tref* yref, const Tref* zref, const int natom,
                     const std::string &copy_msg, const TestPriority do_test,
                     const double overall = 0.0, const double bounds = tiny) {
  diffCoordinates(xcrd, ycrd, zcrd, nullptr, nullptr, nullptr, 1.0, xref, yref, zref, nullptr,
                  nullptr, nullptr, 1.0, natom, nullptr, nullptr, copy_msg, do_test, overall,
                  bounds);
}

template <typename Tcrd, typename Tref>
void diffCoordinates(const Tcrd* xcrd, const Tcrd* ycrd, const Tcrd* zcrd, const double crd_scale,
                     const Tref* xref, const Tref* yref, const Tref* zref, const double ref_scale,
                     const int natom, const double* umat, const double* umat_ref,
                     const std::string &copy_msg, const TestPriority do_test,
                     const double overall = 0.0, const double bounds = tiny) {
  diffCoordinates(xcrd, ycrd, zcrd, nullptr, nullptr, nullptr, crd_scale, xref, yref, zref,
                  nullptr, nullptr, nullptr, ref_scale, natom, umat, umat_ref,
                  copy_msg, do_test, overall, bounds);
}

template <typename Tcrd, typename Tref>
void diffCoordinates(const Tcrd* xcrd, const Tcrd* ycrd, const Tcrd* zcrd, const Tref* xref,
                     const Tref* yref, const Tref* zref, const int natom, const double* umat,
                     const double* umat_ref, const std::string &copy_msg,
                     const TestPriority do_test, const double overall = 0.0,
                     const double bounds = tiny) {
  diffCoordinates(xcrd, ycrd, zcrd, nullptr, nullptr, nullptr, 1.0, xref, yref, zref, nullptr,
                  nullptr, nullptr, 1.0, natom, umat, umat_ref, copy_msg, do_test, overall,
                  bounds);
}


// CHECK
void printPhases(const PhaseSpaceReader &recv_r) {
  for (int i = 0; i < 4; i++) {
    printf("  %12.8lf  %12.8lf  %12.8lf    ", recv_r.xcrd[i], recv_r.xvel[i], recv_r.xfrc[i]);
    printf("  %12.8lf  %12.8lf  %12.8lf\n", recv_r.xalt[i], recv_r.vxalt[i], recv_r.fxalt[i]);
  }
  printf("\n");
}

void printPhases(const PsSynthesisReader &poly_psr, const int sysno) {
  for (int i = 0; i < 4; i++) {
    const int j = poly_psr.atom_starts[sysno] + i;
    printf("  %12.8lf  %12.8lf  %12.8lf    ",
           hostInt95ToDouble(poly_psr.xcrd[j], poly_psr.xcrd_ovrf[j]) * poly_psr.inv_gpos_scale,
           hostInt95ToDouble(poly_psr.xvel[j], poly_psr.xvel_ovrf[j]) * poly_psr.inv_vel_scale,
           hostInt95ToDouble(poly_psr.xfrc[j], poly_psr.xfrc_ovrf[j]) * poly_psr.inv_frc_scale);
    printf("  %12.8lf  %12.8lf  %12.8lf\n",
           hostInt95ToDouble(poly_psr.xalt[j], poly_psr.xalt_ovrf[j]) * poly_psr.inv_gpos_scale,
           hostInt95ToDouble(poly_psr.vxalt[j], poly_psr.vxalt_ovrf[j]) * poly_psr.inv_vel_scale,
           hostInt95ToDouble(poly_psr.fxalt[j], poly_psr.fxalt_ovrf[j]) * poly_psr.inv_frc_scale);
  }
  printf("\n");
}
// END CHECK

//-------------------------------------------------------------------------------------------------
// Replicate a series of integers.
//
// Arguments:
//   length:  The length of the original series { 0, 1, 2, ..., length - 1 }
//   nrep:    The number of times to replicate the series
//-------------------------------------------------------------------------------------------------
std::vector<int> replicateSeries(const int length, const int nrep) {
  std::vector<int> orig(length);
  for (int i = 0; i < length; i++) {
    orig[i] = i;
  }
  return tileVector(orig, nrep);
}

//-------------------------------------------------------------------------------------------------
// Create arrays of various single-system coordinate objects from the entries in a test system
// manager.
//
// Arguments:
//   tsm:     The basis test system manager, containing the coordinates to use as templates
//   tsm_cf:  Array of coordinate frames for each system, built and returned
//   tsm_ps:  Array of PhaseSpace objects for each system, built and returned
//   tsm_cs:  Array of coordinate series for each system, built and returned.  After the first
//            frame, subsequent frames' atoms will be randomly perturbed.
//-------------------------------------------------------------------------------------------------
template <typename T>
PhaseSpaceSynthesis spawnMutableCoordinateObjects(const TestSystemManager &tsm,
                                                  std::vector<CoordinateFrame> *tsm_cf,
                                                  std::vector<CoordinateFrameReader> *tsm_cfr,
                                                  std::vector<PhaseSpace> *tsm_ps,
                                                  std::vector<PhaseSpaceReader> *tsm_psr,
                                                  std::vector<CoordinateSeries<T>> *tsm_cs,
                                                  std::vector<CoordinateSeriesReader<T>> *tsm_csr,
                                                  const int random_seed = 57983) {
  Xoroshiro128pGenerator xrs(random_seed);
  tsm_cf->reserve(tsm.getSystemCount());
  tsm_ps->reserve(tsm.getSystemCount());
  tsm_cs->reserve(tsm.getSystemCount());
  const int small_mol_frames = 4;
  const int nbits = (isFloatingPointScalarType<T>()) ? 0 : 54;
  const double random_factor = pow(2.0, nbits - 4); 
  for (int i = 0; i < tsm.getSystemCount(); i++) {
    tsm_cf->emplace_back(tsm.exportCoordinateFrame(i));
    tsm_ps->emplace_back(tsm.exportPhaseSpace(i));
    tsm_cs->emplace_back(tsm.exportCoordinateFrame(i), small_mol_frames, nbits);
  }
  for (int i = 0; i < tsm.getSystemCount(); i++) {
    CoordinateSeriesWriter<T> icsw = tsm_cs->data()[i].data();
    for (int j = 1; j < small_mol_frames; j++) {
      int ij_start = j * roundUp(icsw.natom, warp_size_int);
      for (int k = ij_start; k < ij_start + icsw.natom; k++) {
        icsw.xcrd[k] += xrs.gaussianRandomNumber() * random_factor;
        icsw.ycrd[k] += xrs.gaussianRandomNumber() * random_factor;
        icsw.zcrd[k] += xrs.gaussianRandomNumber() * random_factor;
      }
      if (icsw.unit_cell != UnitCellType::NONE) {
        const int bdim_offset = j * roundUp(6, warp_size_int);
        const int xfrm_offset = j * roundUp(9, warp_size_int);
        for (int k = 0; k < 6; k++) {
          icsw.boxdim[bdim_offset + k] += (0.01 * xrs.uniformRandomNumber());
        }
        computeBoxTransform(&icsw.boxdim[bdim_offset], &icsw.umat[xfrm_offset],
                            &icsw.invu[xfrm_offset]);
      }
    }
  }

  // Create a vector of readers for the original, double-precision coordinates
  tsm_cfr->reserve(tsm.getSystemCount());
  tsm_psr->reserve(tsm.getSystemCount());
  tsm_csr->reserve(tsm.getSystemCount());
  for (int i = 0; i < tsm.getSystemCount(); i++) {
    tsm_cfr->emplace_back(tsm_cf->data()[i].data());
    tsm_psr->emplace_back(tsm_ps->data()[i].data());
    tsm_csr->emplace_back(tsm_cs->data()[i].data());
  }
  
  // Return a synthesis of everything, times two
  const std::vector<int> replicator = replicateSeries(tsm.getSystemCount(), 2);
  return PhaseSpaceSynthesis(*tsm_ps, tsm.getTopologyPointer(), replicator);
}

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
  Condensate cdns_dbl(poly_ps, PrecisionModel::DOUBLE);
  Condensate cdns_flt(poly_ps, PrecisionModel::SINGLE);
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

  // With the most complex coordinate objects now established, check various coordinate copying
  // overloads.
  std::vector<CoordinateFrame> tsm_cf;
  std::vector<PhaseSpace> tsm_ps;
  std::vector<CoordinateSeries<llint>> tsm_cs;
  std::vector<CoordinateFrameReader> tsm_cfr;
  std::vector<PhaseSpaceReader> tsm_psr;
  std::vector<CoordinateSeriesReader<llint>> tsm_csr;
  PhaseSpaceSynthesis tsm_psyn = spawnMutableCoordinateObjects(tsm, &tsm_cf, &tsm_cfr, &tsm_ps,
                                                               &tsm_psr, &tsm_cs, &tsm_csr);
  Condensate tsm_cdns(tsm_psyn);

  // Test coordinate copying into a PhaseSpace object
  PhaseSpace recv(tsm_cf[0].getAtomCount());
  const PhaseSpaceReader recv_r = recv.data();
  coordCopy(&recv, tsm_cf[0]);
  coordCopy(&recv, TrajectoryKind::VELOCITIES, tsm_cs[0].exportFrame(1));
  coordCopy(&recv, TrajectoryKind::FORCES, tsm_cs[0].exportFrame(2));
  coordCopy(&recv, TrajectoryKind::POSITIONS, CoordinateCycle::ALTERNATE, tsm_cf[0]);
  coordCopy(&recv, TrajectoryKind::FORCES, CoordinateCycle::ALTERNATE, tsm_cdns,
            tsm.getSystemCount());

  // The X, Y, and Z coordinates of the PhaseSpace object recv's PRIMARY coordinates should match
  // those of the test system manager's first coordinate frame, exactly.
  diffCoordinates<double, double>(recv_r.xcrd, recv_r.ycrd, recv_r.zcrd, tsm_cfr[0].xcrd,
                                  tsm_cfr[0].ycrd, tsm_cfr[0].zcrd, recv_r.natom,
                                  "CoordinateFrame => PhaseSpace(POSITIONS, PRIMARY)",
                                  tsm.getTestingStatus());

  // The PhaseSpace object recv's PRIMARY velocities and forces should match those of the second
  // and third frames of the series created for the same system, to a precision of about one part
  // in sixteen quadrillion (just over 64-bit floating point precision).
  const int second_frame_offset =     roundUp(tsm_csr[0].natom, warp_size_int);
  const int third_frame_offset  = 2 * roundUp(tsm_csr[0].natom, warp_size_int);
  diffCoordinates<double, llint>(recv_r.xvel, recv_r.yvel, recv_r.zvel, 1.0,
                                 &tsm_csr[0].xcrd[second_frame_offset],
                                 &tsm_csr[0].ycrd[second_frame_offset],
                                 &tsm_csr[0].zcrd[second_frame_offset], tsm_csr[0].gpos_scale,
                                 recv_r.natom,
                                 "CoordinateSeries(1) => PhaseSpace(VELOCITIES, PRIMARY)",
                                 tsm.getTestingStatus());
  diffCoordinates<double, llint>(recv_r.xfrc, recv_r.yfrc, recv_r.zfrc, 1.0,
                                 &tsm_csr[0].xcrd[third_frame_offset],
                                 &tsm_csr[0].ycrd[third_frame_offset],
                                 &tsm_csr[0].zcrd[third_frame_offset], tsm_csr[0].gpos_scale,
                                 recv_r.natom,
                                 "CoordinateSeries(2) => PhaseSpace(FORCES, PRIMARY)",
                                 tsm.getTestingStatus());
  const CondensateWriter tsm_cdnsw = tsm_cdns.data();
  const int second_set_offset = tsm_cdnsw.atom_starts[tsm.getSystemCount()];
  diffCoordinates<double, float>(recv_r.fxalt, recv_r.fyalt, recv_r.fzalt,
                                 &tsm_cdnsw.xcrd_sp[second_set_offset],
                                 &tsm_cdnsw.ycrd_sp[second_set_offset],
                                 &tsm_cdnsw.zcrd_sp[second_set_offset], recv_r.natom,
                                 "Condensate(" + std::to_string(tsm_cdnsw.system_count) +
                                 ") => PhaseSpace(FORCES, ALTERNATE)", tsm.getTestingStatus());

  // Bring synthesis objects to the fore, filling out more aspects of the recv PhaseSpace object.
  CoordinateFrame buffer_a_cf(tsm_cf[0].getAtomCount()), buffer_b_cf(tsm_cf[0].getAtomCount());
  coordCopy(&buffer_a_cf, tsm_psyn, tsm.getSystemCount());
  coordCopy(&tsm_psyn, 0, TrajectoryKind::FORCES, CoordinateCycle::ALTERNATE, buffer_a_cf);
  coordCopy(&buffer_b_cf, tsm_psyn, 0, TrajectoryKind::FORCES, CoordinateCycle::ALTERNATE);
  coordCopy(&recv, TrajectoryKind::FORCES, CoordinateCycle::ALTERNATE, buffer_a_cf);
  diffCoordinates<double, double>(recv_r.fxalt, recv_r.fyalt, recv_r.fzalt, tsm_cfr[0].xcrd,
                                  tsm_cfr[0].ycrd, tsm_cfr[0].zcrd, recv_r.natom,
                                  "CoordinateFrame => PhaseSpace(FORCES, ALTERNATE)",
                                  tsm.getTestingStatus(), 3.477e-9, 5.0e-9);
  CHECK_THROWS(coordCopy(&recv, TrajectoryKind::VELOCITIES, CoordinateCycle::ALTERNATE, tsm_cdns,
                         1), "Systems with mismatching numbers of atoms were submitted to a copy "
               "operation between PhaseSpace and Condensate objects.");
  coordCopy(&recv, TrajectoryKind::VELOCITIES, CoordinateCycle::ALTERNATE, tsm_cdns, 0);
  diffCoordinates<double, double>(recv_r.vxalt, recv_r.vyalt, recv_r.vzalt, tsm_cfr[0].xcrd,
                                  tsm_cfr[0].ycrd, tsm_cfr[0].zcrd, recv_r.natom,
                                  "CoordinateFrame => PhaseSpace(FORCES, ALTERNATE)",
                                  tsm.getTestingStatus(), 1.437e-7, 2.7e-7);

  // Push coordinates into a new PhaseSpaceSynthesis with higher precision
  const std::vector<int> replicator = replicateSeries(tsm.getSystemCount(), 2);
  PhaseSpaceSynthesis test_psyn(tsm_ps, replicator, tsm.getTopologyPointer(), replicator, 62, 24,
                                61, 60);
  PsSynthesisWriter test_psynw = test_psyn.data();
  coordCopy(&test_psyn, 1, TrajectoryKind::VELOCITIES, tsm_cf[1]);
  const int sys1_start = test_psynw.atom_starts[1];
  diffCoordinates(&test_psynw.xvel[sys1_start], &test_psynw.yvel[sys1_start],
                  &test_psynw.zvel[sys1_start], &test_psynw.xvel_ovrf[sys1_start],
                  &test_psynw.yvel_ovrf[sys1_start], &test_psynw.zvel_ovrf[sys1_start],
                  test_psynw.vel_scale, tsm_cfr[1].xcrd, tsm_cfr[1].ycrd, tsm_cfr[1].zcrd,
                  nullptr, nullptr, nullptr, 1.0, tsm_cfr[1].natom,
                  "CoordinateFrame => PhaseSpaceSynthesis(VELOCITIES, PRIMARY)",
                  tsm.getTestingStatus());
  coordCopy(&test_psyn, 1, TrajectoryKind::FORCES, CoordinateCycle::ALTERNATE, tsm_cs[1], 3);
  const int fourth_frame_offset = 3 * roundUp(tsm_csr[1].natom, warp_size_int);
  diffCoordinates(&test_psynw.fxalt[sys1_start], &test_psynw.fyalt[sys1_start],
                  &test_psynw.fzalt[sys1_start], &test_psynw.fxalt_ovrf[sys1_start],
                  &test_psynw.fyalt_ovrf[sys1_start], &test_psynw.fzalt_ovrf[sys1_start],
                  test_psynw.frc_scale, &tsm_csr[1].xcrd[fourth_frame_offset],
                  &tsm_csr[1].ycrd[fourth_frame_offset], &tsm_csr[1].zcrd[fourth_frame_offset],
                  nullptr, nullptr, nullptr, tsm_csr[1].gpos_scale, tsm_cfr[1].natom,
                  "CoordinateFrame => PhaseSpaceSynthesis(FORCES, ALTERNATE)",
                  tsm.getTestingStatus());
  coordCopy(&test_psyn, 2, TrajectoryKind::FORCES, CoordinateCycle::ALTERNATE, tsm_cdns, 2);
  const int sys2_start = test_psynw.atom_starts[2];
  diffCoordinates(&test_psynw.fxalt[sys2_start], &test_psynw.fyalt[sys2_start],
                  &test_psynw.fzalt[sys2_start], &test_psynw.fxalt_ovrf[sys2_start],
                  &test_psynw.fyalt_ovrf[sys2_start], &test_psynw.fzalt_ovrf[sys2_start],
                  test_psynw.frc_scale, &tsm_cdnsw.xcrd_sp[sys2_start],
                  &tsm_cdnsw.ycrd_sp[sys2_start], &tsm_cdnsw.zcrd_sp[sys2_start], nullptr, nullptr,
                  nullptr, 1.0, tsm_cdnsw.atom_counts[2],
                  "Condensate => PhaseSpaceSynthesis(FORCES, ALTERNATE)", tsm.getTestingStatus());
  coordCopy(&test_psyn, tsm.getSystemCount(), recv);
  CoordinateFrameWriter buffer_a_cfw = buffer_a_cf.data();
  diffCoordinates(&test_psynw.fxalt[second_set_offset], &test_psynw.fyalt[second_set_offset],
                  &test_psynw.fzalt[second_set_offset], &test_psynw.fxalt_ovrf[second_set_offset],
                  &test_psynw.fyalt_ovrf[second_set_offset],
                  &test_psynw.fzalt_ovrf[second_set_offset], test_psynw.frc_scale,
                  buffer_a_cfw.xcrd, buffer_a_cfw.ycrd, buffer_a_cfw.zcrd, nullptr, nullptr,
                  nullptr, 1.0, buffer_a_cfw.natom,
                  "Condensate => PhaseSpaceSynthesis(FORCES, ALTERNATE)",tsm.getTestingStatus());

  // Create more test systems with periodic boundary conditions.  Test that the transformation
  // matrices are transferred properly, while also checking other combinations of origin and
  // destination for coordCopy() overloads.
  const std::vector<std::string> pbc_fi_names = { "symmetry_C1_in_water", "symmetry_C2_in_water",
                                                  "symmetry_C3_in_water", "symmetry_C4_in_water" };
  TestSystemManager tsm_pbc(base_top_dir, "top", pbc_fi_names, base_crd_dir, "inpcrd",
                            pbc_fi_names);
  std::vector<CoordinateFrame> tsm_pbc_cf;
  std::vector<PhaseSpace> tsm_pbc_ps;
  std::vector<CoordinateSeries<double>> tsm_pbc_cs;
  std::vector<CoordinateFrameReader> tsm_pbc_cfr;
  std::vector<PhaseSpaceReader> tsm_pbc_psr;
  std::vector<CoordinateSeriesReader<double>> tsm_pbc_csr;
  PhaseSpaceSynthesis tsm_pbc_psyn = spawnMutableCoordinateObjects(tsm_pbc, &tsm_pbc_cf,
                                                                   &tsm_pbc_cfr, &tsm_pbc_ps,
                                                                   &tsm_pbc_psr, &tsm_pbc_cs,
                                                                   &tsm_pbc_csr);
  Condensate tsm_pbc_cdns(tsm_pbc_cs[2], PrecisionModel::DOUBLE);
  PhaseSpace recv_pbc(tsm_pbc_cfr[2].natom);
  PhaseSpaceReader recv_pbc_r = recv_pbc.data();
  check(tsm_pbc_cdns.getBasis() == CondensateSource::SERIES, "A condensate based on a coordinate "
        "series does not properly report its basis.", tsm_pbc.getTestingStatus());
  check(tsm_pbc_cdns.ownsCoordinates() == false, "A condensate based on a coordinate series of "
        "double-precision reals should rely on the source to store its data, but does not.",
        tsm_pbc.getTestingStatus());
  coordCopy(&recv_pbc, TrajectoryKind::POSITIONS, CoordinateCycle::ALTERNATE, tsm_pbc_cdns, 3);
  coordCopy(&recv_pbc, TrajectoryKind::FORCES, CoordinateCycle::ALTERNATE, tsm_pbc_cdns, 2);
  const CoordinateFrame sm2_f3 = tsm_pbc_cs[2].exportFrame(3);
  const CoordinateFrameReader sm2_f3r = sm2_f3.data();
  diffCoordinates(recv_pbc_r.xalt, recv_pbc_r.yalt, recv_pbc_r.zalt, sm2_f3r.xcrd, sm2_f3r.ycrd,
                  sm2_f3r.zcrd, sm2_f3r.natom, recv_pbc_r.umat_alt, sm2_f3r.umat,
                  "CoordinateSeries => PhaseSpace(POSITIONS, ALTERNATE)",
                  tsm_pbc.getTestingStatus());
  const CoordinateFrame sm2_f2 = tsm_pbc_cs[2].exportFrame(2);
  const CoordinateFrameReader sm2_f2r = sm2_f2.data();
  diffCoordinates(recv_pbc_r.fxalt, recv_pbc_r.fyalt, recv_pbc_r.fzalt, sm2_f2r.xcrd, sm2_f2r.ycrd,
                  sm2_f2r.zcrd, sm2_f2r.natom,
                  "CoordinateSeries => PhaseSpace(FORCES, ALTERNATE)",
                  tsm_pbc.getTestingStatus());

  // Print results
  printTestSummary(oe.getVerbosity());

  return countGlobalTestFailures();
}
