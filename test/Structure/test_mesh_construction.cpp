#include <string>
#include <vector>
#include "copyright.h"
#include "../../src/DataTypes/common_types.h"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/Numerics/split_fixed_precision.h"
#include "../../src/Reporting/error_format.h"
#include "../../src/Reporting/summary_file.h"
#include "../../src/Structure/background_mesh.h"
#include "../../src/Topology/atomgraph.h"
#include "../../src/Topology/atomgraph_enumerators.h"
#include "../../src/UnitTesting/approx.h"
#include "../../src/UnitTesting/stopwatch.h"
#include "../../src/UnitTesting/test_system_manager.h"
#include "../../src/UnitTesting/unit_test.h"

using stormm::data_types::ullint;
using stormm::numerics::int95ToDouble;
using stormm::review::stormmSplash;
using namespace stormm::diskutil;
using namespace stormm::structure;
using namespace stormm::testing;
using namespace stormm::topology;

//-------------------------------------------------------------------------------------------------
// main
//-------------------------------------------------------------------------------------------------
int main(const int argc, const char* argv[]) {

  // Some baseline initialization
  TestEnvironment oe(argc, argv);
  StopWatch timer;
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmSplash();
  }

  // Section 1
  section("Establish various mesh objects");

  // Create various meshes.  Test constructor overloads and traps.
  section(1);
  char osc = osSeparator();
  const std::string crd_base_name = oe.getStormmSourcePath() + osc + "test" + osc + "Trajectory";
  const std::string top_base_name = oe.getStormmSourcePath() + osc + "test" + osc + "Topology";
  TestSystemManager tsm(top_base_name, "top", { "trpcage" },
                        crd_base_name, "inpcrd", { "trpcage" });
  AtomGraph *trpi_ag = tsm.getTopologyPointer(0);
  trpi_ag->modifyAtomMobility(0, trpi_ag->getAtomCount(), MobilitySetting::OFF);
  const std::vector<bool> trpi_mobile = trpi_ag->getAtomMobility();
  const CoordinateFrame trpi_cf = tsm.exportCoordinateFrame(0);
  BackgroundMesh<ullint> bgm_a(GridDetail::OCCLUSION, NonbondedPotential::CLASH);
  bgm_a.setTopologyPointer(tsm.getTopologyPointer(0));
  bgm_a.setCoordinatePointer(&trpi_cf);
  bgm_a.setMeshParameters(10.0, 1.0);
  bgm_a.setProbeRadius(1.4);
  bgm_a.computeField();
  BackgroundMeshWriter<double, ullint> bgmw_a = bgm_a.dpData();
  const std::vector<double> trpi_mesh_dims = { static_cast<double>(bgmw_a.dims.na),
                                               static_cast<double>(bgmw_a.dims.nb),
                                               static_cast<double>(bgmw_a.dims.nc),
                                               bgmw_a.dims.invu[0], bgmw_a.dims.invu[1],
                                               bgmw_a.dims.invu[2], bgmw_a.dims.invu[3],
                                               bgmw_a.dims.invu[4], bgmw_a.dims.invu[5],
                                               bgmw_a.dims.invu[6], bgmw_a.dims.invu[7],
                                               bgmw_a.dims.invu[8] };
  const std::vector<double> trpi_mesh_dims_ans = { 45.0, 43.0, 37.0,  1.0,  0.0,  0.0,
                                                    0.0,  1.0,  0.0,  0.0,  0.0,  1.0 };
  check(trpi_mesh_dims, RelationalOperator::EQUAL, trpi_mesh_dims_ans, "The mesh dimensions "
        "computed for Trp-cage do not meet expectations.", tsm.getTestingStatus());
  BackgroundMesh<double> bgm_b(GridDetail::NONBONDED_FIELD, NonbondedPotential::ELECTROSTATIC,
                               tsm.getTopologyPointer(0), &trpi_cf, 2.5, 1.5);
  BackgroundMeshWriter<double, double> bgmw_b = bgm_b.dpData();
  BackgroundMesh<double> bgm_c(GridDetail::NONBONDED_FIELD, NonbondedPotential::VAN_DER_WAALS,
                               3.15061, 0.1521, VdwCombiningRule::LORENTZ_BERTHELOT,
                               tsm.getTopologyPointer(0), &trpi_cf, 2.5, 1.5);
  BackgroundMeshWriter<double, double> bgmw_c = bgm_c.dpData();
  std::vector<double> gpt_x(bgmw_b.dims.na), gpt_y(bgmw_b.dims.nb), gpt_z(bgmw_b.dims.nb);
  int95ToDouble(gpt_x.data(), bgmw_b.avec_abs_x, bgmw_b.avec_abs_x_ovrf, bgmw_b.dims.na,
                bgmw_b.dims.inv_scale);
  int95ToDouble(gpt_y.data(), bgmw_b.bvec_y, bgmw_b.bvec_y_ovrf, bgmw_b.dims.nb,
                bgmw_b.dims.inv_scale);
  int95ToDouble(gpt_z.data(), bgmw_b.cvec_z, bgmw_b.cvec_z_ovrf, bgmw_b.dims.nc,
                bgmw_b.dims.inv_scale);
  addScalarToVector(&gpt_y, int95ToDouble(bgmw_b.dims.orig_y) * bgmw_b.dims.inv_scale);
  addScalarToVector(&gpt_z, int95ToDouble(bgmw_b.dims.orig_z) * bgmw_b.dims.inv_scale);
  std::vector<double> lj_check, lj_check_ans, qq_check, qq_check_ans;
  const std::vector<double> trpi_sig = trpi_ag->getLennardJonesSigma<double>();
  const std::vector<double> trpi_eps = trpi_ag->getLennardJonesEpsilon<double>();
  const NonbondedKit<double> trpi_nbk = trpi_ag->getDoublePrecisionNonbondedKit();
  const std::vector<bool> trpi_frozen = trpi_ag->getAtomMobility();
  const CoordinateFrameReader trpi_cfr = trpi_cf.data();
  for (int i = 0; i < bgmw_b.dims.na; i++) {
    for (int j = 0; j < bgmw_b.dims.nb; j++) {
      for (int k = 0; k < bgmw_b.dims.nc; k++) {
        double u_vdw  = 0.0;
        double u_elec = 0.0;
        for (int m = 0; m < trpi_cfr.natom; m++) {
          if (trpi_mobile[m]) {
            continue;
          }
          const double dx = gpt_x[i] - trpi_cfr.xcrd[m];
          const double dy = gpt_y[j] - trpi_cfr.ycrd[m];
          const double dz = gpt_z[k] - trpi_cfr.zcrd[m];
          const double r2 = (dx * dx) + (dy * dy) + (dz * dz);
          const double invr2 = 1.0 / r2;
          const double invr  = sqrt(invr2);
          const double invr6 = invr2 * invr2 * invr2;
          const int m_idx = trpi_nbk.lj_idx[m];
          const double m_q = trpi_nbk.charge[m];
          const double m_sig = 0.5 * (trpi_sig[m_idx] + 3.15061);
          const double m_eps = sqrt(trpi_eps[m_idx] * 0.1521);
          const double m_sig6 = pow(m_sig, 6.0);
          u_vdw += 4.0 * m_eps * ((m_sig6 * m_sig6 * invr6) - m_sig6) * invr6;
          u_elec += m_q * invr;
        }
        if (u_vdw < 2.0) {
          lj_check_ans.push_back(u_vdw);
          qq_check_ans.push_back(u_elec * trpi_nbk.coulomb_constant);
          const size_t mesh_coef_idx = 64 * ((((k * bgmw_b.dims.nb) + j) * bgmw_b.dims.na) + i);
          lj_check.push_back(bgmw_c.coeffs[mesh_coef_idx]);
          qq_check.push_back(bgmw_b.coeffs[mesh_coef_idx]);
        }
      }
    }
  }
  check(lj_check, RelationalOperator::EQUAL, lj_check_ans, "Lennard-Jones energies computed for "
        "mesh grid points do not match the constant coefficients for the corresponding mesh "
        "elements.", tsm.getTestingStatus());
  check(qq_check, RelationalOperator::EQUAL, qq_check_ans, "Electrostatic energies computed for "
        "mesh grid points do not match the constant coefficients for the corresponding mesh "
        "elements.", tsm.getTestingStatus());
  
  // CHECK
  for (int i = 0; i < bgmw_b.dims.na; i++) {
    for (int j = 0; j < bgmw_b.dims.nb; j++) {
      printf("|");
      for (int k = 0; k < bgmw_b.dims.nc; k++) {
        const size_t coeff_idx = 64 * ((((k * bgmw_c.dims.nb) + j) * bgmw_c.dims.na) + i);
        const double vdw_j = bgmw_c.coeffs[coeff_idx];
        const double vdw_e = bgmw_b.coeffs[coeff_idx] * 0.1;
        char u_code  = ' ';
        char u_code2 = ' ';
        if (vdw_j > 2.0) {
          u_code  = '[';
          u_code2 = ']';
        }
        else {
          if (vdw_e < -2.5) {
            u_code  = '#';
            u_code2 = '#';
          }
          else if (vdw_e < -2.0) {
            u_code  = '#';
            u_code2 = '=';
          }
          else if (vdw_e < -1.5) {
            u_code  = '=';
            u_code2 = '=';
          }
          else if (vdw_e < -1.0) {
            u_code  = '=';
            u_code2 = '-';
          }
          else if (vdw_e < -0.5) {
            u_code  = '-';
            u_code2 = '-';
          }
          else if (vdw_e < 0.0) {
            u_code  = ' ';
            u_code2 = ' ';
          }
          else if (vdw_e < 0.5) {
            u_code  = '.';
            u_code2 = '.';
          }
          else if (vdw_e < 1.0) {
            u_code  = '.';
            u_code2 = '+';
          }
          else if (vdw_e < 1.5) {
            u_code  = '+';
            u_code2 = '+';
          }
          else if (vdw_e < 2.0) {
            u_code  = '+';
            u_code2 = 'X';
          }
          else if (vdw_e < 2.5) {
            u_code  = 'X';
            u_code2 = 'X';
          }
          else {
            u_code  = '@';
            u_code2 = '@';
          }
        }
        printf("%c%c", u_code, u_code2);
      }
      printf("|\n");
    }
    printf("\n");
  }
  // END CHECK

  // Summary evaluation
  if (oe.getDisplayTimingsOrder()) {
    timer.assignTime(0);
    timer.printResults();
  }
  printTestSummary(oe.getVerbosity());

  return countGlobalTestFailures();
}
