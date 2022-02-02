#include "../../src/Constants/behavior.h"
#include "../../src/DataTypes/common_types.h"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/Potential/energy_enumerators.h"
#include "../../src/Potential/scorecard.h"
#include "../../src/Potential/valence_potential.h"
#include "../../src/Random/random.h"
#include "../../src/Reporting/error_format.h"
#include "../../src/Topology/atomgraph.h"
#include "../../src/Topology/atomgraph_abstracts.h"
#include "../../src/Trajectory/phasespace.h"
#include "../../src/Trajectory/trajectory_enumerators.h"
#include "../../src/UnitTesting/unit_test.h"

using omni::constants::ExceptionResponse;
using omni::data_types::llint;
using omni::data_types::double3;
using omni::data_types::float3;
using omni::data_types::llint3;
using omni::diskutil::DrivePathType;
using omni::diskutil::getDrivePathType;
using omni::diskutil::osSeparator;
using omni::energy::EvaluateForce;
using omni::energy::ScoreCard;
using omni::random::Xoshiro256ppGenerator;
using omni::topology::AtomGraph;
using omni::topology::ValenceKit;
using omni::trajectory::CoordinateFileKind;
using omni::trajectory::PhaseSpace;
using omni::trajectory::TrajectoryKind;
using namespace omni::testing;

//-------------------------------------------------------------------------------------------------
// Compute the force between two particles given a set of harmonic bond parameters.
//
//
//-------------------------------------------------------------------------------------------------
template <typename T, typename T3>
T3 bond_force(T3 crd1, T3 crd2, T equil, T stiff) {
}

//-------------------------------------------------------------------------------------------------
// main
//-------------------------------------------------------------------------------------------------
int main(int argc, char* argv[]) {

  // Baseline variables
  TestEnvironment oe(argc, argv);

  // Perform arithmetic operations similar to a single bond force computation for a range of
  // distances and stiffness constants that emulates what would be encountered in a typical
  // simulation.  Compute the quantities with coordinates and parameter constants in double
  // precision, then in single precision, then in split single precision with fixed-precision
  // intermediaries and accumulation.
  const int nsample = 256;
  const int ndistance = 1200;
  const int npts = nsample * ndistance;
  const double distance_discretization = 0.001;
  Xoshiro256ppGenerator xsr_rng(90384011);
  std::vector<double3> crd1(npts), crd2(npts), frc1(npts);
  std::vector<float3> f_crd1(npts), f_crd2(npts), f_frc(npts);
  std::vector<llint3> lli_crd1(npts), lli_crd2(npts), lli_frc(npts);
  std::vector<double> equil(npts), stiff(npts);
  std::vector<float> f_equil(npts), f_stiff(npts);
  std::vector<llint> lli_equil(npts), lli_stiff(npts);
  for (int i = 0; i < npts; i++) {

    // Select random coordinates
    crd1[i].x = xsr_rng.gaussianRandomNumber();
    crd1[i].y = xsr_rng.gaussianRandomNumber();
    crd1[i].z = xsr_rng.gaussianRandomNumber();
    crd2[i].x = xsr_rng.gaussianRandomNumber();
    crd2[i].y = xsr_rng.gaussianRandomNumber();
    crd2[i].z = xsr_rng.gaussianRandomNumber();

    // Extend or contract the distance along the same direction to meet a particular value
    const double target_distance = 0.8 +
                                   ((static_cast<double>(i) + xsr_rng.uniformRandomNumber()) *
                                    distance_discretization);
    const double pdx = crd2[i].x - crd1[i].x;
    const double pdy = crd2[i].y - crd1[i].y;
    const double pdz = crd2[i].z - crd1[i].z;
    const double factor = target_distance / sqrt((pdx * pdx) + (pdy * pdy) + (pdz * pdz));
    crd2[i].x = crd1[i].x + (factor * pdx);
    crd2[i].y = crd1[i].y + (factor * pdy);
    crd2[i].z = crd1[i].z + (factor * pdz);

    // The locations are now established.  Recompute displacements and total distance to keep
    // things in familiar nomenclature.
    const double dx = crd2[i].x - crd1[i].x;
    const double dy = crd2[i].y - crd1[i].y;
    const double dz = crd2[i].z - crd1[i].z;
    const double r = sqrt((dx * dx) + (dy * dy) + (dz * dz));

    // Select a random equilibrium length, based on the target length, and a stiffness constant
    stiff[i] = 350.0 + (250.0 * xsr_rng.uniformRandomNumber());
    equil[i] = target_distance + ((60.0 / stiff[i]) * xsr_rng.gaussianRandomNumber());
    const double dl = target_distance - equil[i];
    const double fmag = 2.0 * stiff[i] * dl / r;
    frc1[i].x = fmag * dx;
    frc1[i].y = fmag * dy;
    frc1[i].z = fmag * dz;

    // CHECK
    if (i < 10) {
      printf("  %12.7lf %12.7lf %12.7lf\n", frc1[i].x, frc1[i].y, frc1[i].z);
    }
    // END CHECK
    
    // Perform the computation with all coordinates reduced to floating-point numbers
    f_crd1[i].x = crd1[i].x;
    f_crd1[i].y = crd1[i].y;
    f_crd1[i].z = crd1[i].z;
    f_crd2[i].x = crd2[i].x;
    f_crd2[i].y = crd2[i].y;
    f_crd2[i].z = crd2[i].z;
    f_equil[i] = equil[i];
    f_stiff[i] = stiff[i];    
  }

  // Read topology
  const char osc = osSeparator();
  const std::string topology_home = oe.getOmniSourcePath() + osc + "test" + osc + "Topology";
  const std::string trpcage_top = topology_home + osc + "trpcage.top";
  const std::string coordinate_home = oe.getOmniSourcePath() + osc + "test" + osc + "Trajectory";
  const std::string trpcage_crd = coordinate_home + osc + "trpcage.inpcrd";
  AtomGraph trpcage_ag;
  PhaseSpace trpcage_ps;
  ScoreCard all_systems_sc(1);
  const int trpcage_idx = 0;
  if (getDrivePathType(trpcage_top) == DrivePathType::FILE) {
    trpcage_ag.buildFromPrmtop(trpcage_top, ExceptionResponse::SILENT);
  }
  else {
    rtErr("The Trp-cage topology " + trpcage_top + " was not found.  The benchmark cannot run.",
          "SplitForceAccumulation", "split_valence");
  }
  if (getDrivePathType(trpcage_crd) == DrivePathType::FILE) {
    trpcage_ps.buildFromFile(trpcage_crd, CoordinateFileKind::AMBER_INPCRD);
  }
  else {
    rtErr("The Trp-cage topology " + trpcage_crd + " was not found.  The benchmark cannot run.",
          "SplitForceAccumulation", "split_valence");
  }

  // Evaluate the forces due to bond and bond angle terms in double precision computations
  const TrajectoryKind tkind = TrajectoryKind::FORCES;
  trpcage_ps.initializeForces();
  const double trpcage_bond_e = evaluateBondTerms(trpcage_ag, &trpcage_ps, &all_systems_sc,
                                                  EvaluateForce::YES, trpcage_idx);
  const std::vector<double> trpcage_bond_frc = trpcage_ps.getInterlacedCoordinates(tkind);
  trpcage_ps.initializeForces();
  const double trpcage_angl_e = evaluateAngleTerms(trpcage_ag, &trpcage_ps, &all_systems_sc,
                                                   EvaluateForce::YES, trpcage_idx);
  const std::vector<double> trpcage_angl_frc = trpcage_ps.getInterlacedCoordinates(tkind);
  trpcage_ps.initializeForces();

  // Compute the forces in single precision.  Accumulate in fixed precision.
  
  // CHECK
  
  printf("BondAngle = [\n");
  for (int i = 0; i < trpcage_ag.getAtomCount(); i += 25) {
    printf("  %12.7lf %12.7lf %12.7lf %12.7lf %12.7lf %12.7lf\n", trpcage_bond_frc[3 * i],
           trpcage_bond_frc[(3 * i) + 1], trpcage_bond_frc[(3 * i) + 2], trpcage_angl_frc[3 * i],
           trpcage_angl_frc[(3 * i) + 1], trpcage_angl_frc[(3 * i) + 2]);
  }
  printf("];\n");
  // END CHECK
}
