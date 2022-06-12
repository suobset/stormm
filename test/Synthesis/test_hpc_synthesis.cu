// -*-c++-*-
#include <cuda.h>
#include <cuda_runtime.h>
#include <cusolverDn.h>
#include <nvml.h>
#include "../../src/Accelerator/hpc_config.cuh"
#include "../../src/Constants/scaling.h"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/Math/rounding.h"
#include "../../src/Namelists/nml_files.h"
#include "../../src/Parsing/textfile.h"
#include "../../src/Potential/cacheresource.h"
#include "../../src/Potential/hpc_valence_potential.cuh"
#include "../../src/Reporting/error_format.h"
#include "../../src/Synthesis/phasespace_synthesis.h"
#include "../../src/Synthesis/systemcache.h"
#include "../../src/Synthesis/atomgraph_synthesis.h"
#include "../../src/Synthesis/synthesis_abstracts.h"
#include "../../src/Synthesis/valence_workunit.h"
#include "../../src/Topology/atomgraph.h"
#include "../../src/Topology/atomgraph_abstracts.h"
#include "../../src/Trajectory/phasespace.h"
#include "../../src/UnitTesting/unit_test.h"
#include "../../src/UnitTesting/stopwatch.h"

using namespace omni::card;
using namespace omni::constants;
using namespace omni::errors;
using namespace omni::diskutil;
using namespace omni::energy;
using namespace omni::math;
using namespace omni::parse;
using namespace omni::synthesis;
using namespace omni::testing;
using namespace omni::topology;
using namespace omni::trajectory;

//-------------------------------------------------------------------------------------------------
// main
//-------------------------------------------------------------------------------------------------
int main(const int argc, const char* argv[]) {

  // Some baseline initialization
  TestEnvironment oe(argc, argv);
  StopWatch timer;
  HpcConfig gpu_config(ExceptionResponse::WARN);
  std::vector<int> my_gpus = gpu_config.getGpuDevice(1);
  GpuDetails gpu = gpu_config.getGpuInfo(my_gpus[0]);
  int nblocks = gpu.getSMPCount();
  int nthreads = gpu.getMaxThreadsPerBlock();
  int max_vwu_atoms = maximum_valence_work_unit_atoms;
  if (gpu.getArchMajor() == 6 && gpu.getArchMinor() == 1) {
    nblocks *= 2;
    nthreads /= 2;
    max_vwu_atoms /= 2;
  }
  
  // Section 1
  section("Coordinate compilation and staging");
  
  // Section 2
  section("Topology compilation and staging");

  // Collect coordinates and topologies
  const char osc = osSeparator();
  std::string buffer("&files\n  -p ");
  buffer += oe.getOmniSourcePath() + osc + "test" + osc + "Namelists" + osc + "topol" + osc +
            ".*.top\n  -c ";
  buffer += oe.getOmniSourcePath() + osc + "test" + osc + "Namelists" + osc + "coord" + osc +
            ".*.inpcrd\n&end\n";
  const TextFile tf(buffer, TextOrigin::RAM);
  int start_line = 0;
  FilesControls fcon(tf, &start_line);
  const SystemCache sysc(fcon, ExceptionResponse::SILENT, MapRotatableGroups::NO, &timer);

  // Form the syntheses of topologies and coordinates
  section(1);
  const int nsys = sysc.getSystemCount();
  const TestPriority do_tests = (nsys > 0) ? TestPriority::CRITICAL : TestPriority::ABORT;
  if (nsys == 0) {
    rtWarn("No system topology and coordinate pairs were detected.  Subsequent tests will be "
           "skipped.", "test_hpc_synthesis");
  }
  std::vector<int> topology_indices(nsys, 0);
  for (int i = 0; i < nsys; i++) {
    topology_indices[i] = i;
  }
  AtomGraphSynthesis poly_ag(sysc.getSystemTopologyPointerCC(), topology_indices,
                             ExceptionResponse::WARN, &timer);
  PhaseSpaceSynthesis poly_ps(sysc);
  check(poly_ag.getSystemCount(), RelationalOperator::EQUAL, poly_ps.getSystemCount(),
        "PhaseSpaceSynthesis and AtomGraphSynthesis objects formed from the same SystemCache have "
        "different numbers of systems inside of them.", do_tests);

  // Upload the compiled systems and check the results
  PsSynthesisWriter psw = poly_ps.data();
  poly_ag.upload();
  poly_ps.upload();
  std::vector<double> gpu_charges = poly_ag.getPartialCharges<double>(HybridTargetLevel::DEVICE);
  int padded_atom_count = 0;
  for (int i = 0; i < nsys; i++) {
    if (i == nsys - 1) {
      padded_atom_count += poly_ag.getSystemTopologyPointer(i)->getAtomCount();
    }
    else {
      padded_atom_count += roundUp(poly_ag.getSystemTopologyPointer(i)->getAtomCount(),
                                   warp_size_int);
    }
  }
  std::vector<double> rbt_charges(padded_atom_count, 0.0);
  int atom_offset = 0;
  for (int i = 0; i < nsys; i++) {
    const AtomGraph *iag_ptr = poly_ag.getSystemTopologyPointer(i);
    std::vector<double> ichg = iag_ptr->getPartialCharge<double>();
    const int natom = iag_ptr->getAtomCount();
    for (int j = 0; j < natom; j++) {
      rbt_charges[atom_offset + j] = ichg[j];
    }
    const int padded_natom = roundUp(natom, warp_size_int);
    if (i < nsys - 1) {
      for (int j = natom; j < padded_natom; j++) {
        gpu_charges[atom_offset + j] = 0.0;
      }
    }
    atom_offset += padded_natom;
  }
  check(gpu_charges, RelationalOperator::EQUAL, rbt_charges, "Charges pulled from the GPU in an "
        "AtomGraphSynthesis object do not meet expectations.", do_tests);

  // Allocate resources for various kernels
  CacheResource tb_space(nblocks, max_vwu_atoms);
  MolecularMechanicsControls mmctrl;
  mmctrl.primeWorkUnitCounters(gpu, poly_ag);
  ScoreCard sc(nsys, 1, 32);
  
  // Launch the valence evaluation kernel
  launchValenceSp(poly_ag, &mmctrl, &poly_ps, &sc, &tb_space, EvaluateForce::YES,
                  EvaluateEnergy::YES, VwuGoal::ACCUMULATE, ForceAccumulationMethod::SPLIT, gpu);
  
  // Summary evaluation
  if (oe.getDisplayTimingsOrder()) {
    timer.assignTime(0);
    timer.printResults();
  }
  printTestSummary(oe.getVerbosity());

  return 0;
}
