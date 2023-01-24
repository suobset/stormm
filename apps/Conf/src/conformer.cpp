#include "../../../src/Accelerator/hpc_config.h"
#include "../../../src/Chemistry/chemistry_enumerators.h"
#include "../../../src/FileManagement/file_listing.h"
#include "../../../src/FileManagement/file_util.h"
#include "../../../src/Math/series_ops.h"
#include "../../../src/Math/summation.h"
#include "../../../src/Math/vector_ops.h"
#ifdef STORMM_USE_HPC
#  include "../../../src/Accelerator/kernel_manager.h"
#  include "../../../src/Accelerator/gpu_details.h"
#  include "../../../src/MolecularMechanics/hpc_minimization.h"
#else
#  include "../../../src/MolecularMechanics/minimization.h"
#endif
#include "../../../src/MoleculeFormat/mdlmol.h"
#include "../../../src/Namelists/nml_minimize.h"
#include "../../../src/Namelists/nml_precision.h"
#include "../../../src/Namelists/nml_random.h"
#include "../../../src/Namelists/nml_solvent.h"
#include "../../../src/Namelists/user_settings.h"
#include "../../../src/Potential/scorecard.h"
#include "../../../src/Potential/static_exclusionmask.h"
#include "../../../src/Random/random.h"
#include "../../../src/Reporting/error_format.h"
#include "../../../src/Reporting/help_messages.h"
#include "../../../src/Reporting/reporting_enumerators.h"
#include "../../../src/Structure/clash_detection.h"
#include "../../../src/Synthesis/atomgraph_synthesis.h"
#include "../../../src/Synthesis/condensate.h"
#include "../../../src/Synthesis/phasespace_synthesis.h"
#include "../../../src/Synthesis/static_mask_synthesis.h"
#include "../../../src/Synthesis/systemcache.h"
#include "../../../src/Topology/atomgraph_enumerators.h"
#include "../../../src/Trajectory/coordinate_copy.h"
#include "../../../src/UnitTesting/stopwatch.h"
#include "../../../src/UnitTesting/unit_test.h"
#include "analysis.h"
#include "setup.h"

using stormm::random::Xoshiro256ppGenerator;
using stormm::testing::StopWatch;
using namespace stormm::card;
using namespace stormm::chemistry;
using namespace stormm::diskutil;
using namespace stormm::display;
using namespace stormm::energy;
using namespace stormm::errors;
using namespace stormm::math;
using namespace stormm::mm;
using namespace stormm::namelist;
using namespace stormm::structure;
using namespace stormm::synthesis;
using namespace stormm::trajectory;
using namespace stormm::topology;
using namespace conf_app::setup;
using namespace conf_app::analysis;

//-------------------------------------------------------------------------------------------------
// Display a general help message for this program.
//-------------------------------------------------------------------------------------------------
void displayGeneralHelpMessage(const std::vector<std::string> &all_namelists) {
  const std::string base_msg =
    terminalFormat("A program for probing conformations of chemical structures and returning "
                   "those with the lowest energy.\n\n", "conformer");
  std::string list_of_nml;
  for (size_t i = 0; i < all_namelists.size(); i++) {
    list_of_nml += "  - " + all_namelists[i] + "\n";
  }
  const std::string nml_msg =
    terminalFormat("Applicable namelists (re-run with one of these terms as the command-line "
                   "argument, IN QUOTES, i.e. \"&files\" or '&files', for further details):\n" +
                   list_of_nml, nullptr, nullptr, 0, 2, 2);
  printf("%s", base_msg.c_str());
  printf("%s", nml_msg.c_str());
  printf("\n");
}

//-------------------------------------------------------------------------------------------------
// main
//-------------------------------------------------------------------------------------------------
int main(int argc, const char* argv[]) {

  // Check for a help message
  const std::vector<std::string> my_namelists = {
    "&files", "&conformer", "&restraint", "&solvent", "&random", "&minimize", "&report",
    "&precision"
  };
  if (detectHelpSignal(argc, argv)) {
    displayGeneralHelpMessage(my_namelists);
    return 0;
  }
  if (displayNamelistHelp(argc, argv, my_namelists)) { 
    return 0;
  }

  // Wall time tracking
  StopWatch master_timer("Master timings for conformer.omni");
  master_timer.addCategory(tm_input_parsing);
  master_timer.addCategory(tm_feature_detection);
  master_timer.addCategory(tm_coordinate_expansion);
  master_timer.addCategory("Energy minimization");

  // Read information from the command line and initialize the UserSettings object
  UserSettings ui(argc, argv, AppName::CONFORMER);

  // Get details of the GPU to use
#ifdef STORMM_USE_HPC
  const HpcConfig gpu_config(ExceptionResponse::WARN);
  const std::vector<int> my_gpus = gpu_config.getGpuDevice(1);
  const GpuDetails gpu = gpu_config.getGpuInfo(my_gpus[0]);
#endif
  
  // Boot up the master random number generator
  const RandomControls rngcon = ui.getRandomNamelistInfo();
  Xoshiro256ppGenerator xrs(rngcon.getRandomSeed(), rngcon.getWarmupCycleCount());
  
  // Read topologies and coordinate files
  std::vector<MdlMol> sdf_recovery;
  SystemCache sc(ui.getFilesNamelistInfo(), &sdf_recovery, ui.getExceptionBehavior(),
                 MapRotatableGroups::NO, ui.getPrintingPolicy(), &master_timer);
  master_timer.assignTime(1);
  
  // Prepare the workspace for multiple rounds of refinement.  Establish the conditions for
  // generative modeling.  Create non-bonded exclusions masks for each topology in CPU (and, if
  // applicable, GPU) applications to facilitate clash detection routines and later energy
  // calculations.
  setGenerativeConditions(ui, &sc, sdf_recovery, &master_timer);
  const std::vector<int> replica_counts = calculateReplicaCounts(ui.getConformerNamelistInfo(),
                                                                 sc, &master_timer);  
  PhaseSpaceSynthesis sandbox = buildReplicaWorkspace(replica_counts, sc, &master_timer);
  const std::vector<AtomGraph*> unique_topologies = sandbox.getUniqueTopologies();
  const std::vector<int> conformer_topology_indices = sandbox.getUniqueTopologyIndices();
  const std::vector<AtomGraph*> replica_topologies = sandbox.getSystemTopologyPointer();
  const std::vector<RestraintApparatus*> replica_restraints =
    buildReplicaRestraints(replica_counts, sc, &master_timer);
  if (replica_topologies.size() != replica_restraints.size()) {
    rtErr("Replica topology and restraint apparatus pointer vectors must contain the same number "
          "of elements (" + std::to_string(replica_topologies.size()) + " / " +
          std::to_string(replica_restraints.size()) + ").", "main");
  }
  AtomGraphSynthesis sandbox_ag(unique_topologies, replica_restraints, conformer_topology_indices, 
                                incrementingSeries<int>(0, replica_restraints.size(), 1));
#ifdef STORMM_USE_HPC
  StaticExclusionMaskSynthesis sandbox_sems(sandbox_ag.getTopologyPointers(),
                                            sandbox_ag.getTopologyIndices());
#endif
  // Parse the rotatable bonds, cis-trans isomers, and chiral centers in each system to prepare
  // a much larger list of all the possible conformers that each system might be able to access.
  expandConformers(&sandbox, replica_counts, ui, sc, &xrs, &master_timer);
  
  // Complete the work unit construction.  Upload data to the GPU, if applicable.
#ifdef STORMM_USE_HPC
  InitializationTask init_task;
  switch (sandbox_ag.getImplicitSolventModel()) {
  case ImplicitSolventModel::NONE:
    init_task = InitializationTask::GENERAL_MINIMIZATION;
    break;
  case ImplicitSolventModel::HCT_GB:
  case ImplicitSolventModel::OBC_GB:
  case ImplicitSolventModel::OBC_GB_II:
  case ImplicitSolventModel::NECK_GB:
  case ImplicitSolventModel::NECK_GB_II:
    init_task = InitializationTask::GB_MINIMIZATION;
    break;
  }
  sandbox_ag.loadNonbondedWorkUnits(sandbox_sems, init_task, 0, gpu);
  sandbox_ag.upload();
  sandbox_sems.upload();
  sandbox.upload();
#endif
  master_timer.assignTime(3);

  // Loop over all systems and perform energy minimizations.
  const PrecisionControls preccon = ui.getPrecisionNamelistInfo();
  std::vector<StaticExclusionMask> sandbox_masks;
  sandbox_masks.reserve(unique_topologies.size());
  for (int i = 0; i < unique_topologies.size(); i++) {
    sandbox_masks.emplace_back(unique_topologies[i]);
  }
#ifdef STORMM_USE_HPC
  ScoreCard emin = launchMinimization(sandbox_ag, sandbox_sems, &sandbox,
                                      ui.getMinimizeNamelistInfo(), gpu,
                                      preccon.getValenceMethod(), preccon.getEnergyScalingBits());
#else
  ScoreCard emin(sandbox.getSystemCount(), ui.getMinimizeNamelistInfo().getTotalCycles() + 1,
                 preccon.getEnergyScalingBits());
  const MinimizeControls mincon = ui.getMinimizeNamelistInfo();
  const int nconf = sandbox.getSystemCount();
  for (int i = 0; i < nconf; i++) {

    // Pull the coordinates out of the synthesis, back into a PhaseSpace object which will keep
    // them in double-precision for the entire minimization calculation, then re-import the
    // modified coordinates to the synthesis to maintain consistency with the HPC protocol.
    // Subsequent analyses can focus on the synthesis.
    PhaseSpace psi = sandbox.exportSystem(i);
    const AtomGraph *agi = sandbox_ag.getSystemTopologyPointer(i);
    ScoreCard emin_i = minimize(&psi, agi, sandbox_masks[conformer_topology_indices[i]], mincon);
    emin.import(emin_i, i, 0);
    coordCopy(&sandbox, i, psi);
  }
#endif
  master_timer.assignTime(4);

  // Download the results from the HPC device, if applicable
#ifdef STORMM_USE_HPC
  sandbox.download();
  emin.download();
  Condensate sandbox_snapshot(sandbox, PrecisionModel::SINGLE, gpu);
#else
  Condensate sandbox_snapshot(sandbox);
#endif

  // For each rotatable bond, compute the value obtained in each conformer.  This will produce a
  // map of the viable minima for various conformations.

  // Get a vector of the best conformations by their indices in the synthesis.
  const std::vector<int> best_confs = filterMinimizedStructures(sandbox, sandbox_masks, emin,
                                                                ui.getConformerNamelistInfo());

  // Print the best conformations for each topological system.
  printResults(sandbox, best_confs, emin, sc, sdf_recovery, ui.getFilesNamelistInfo());

  // Print timings results
  master_timer.printResults();

  return 0;
}
