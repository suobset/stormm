#include <vector>
#include "../../../src/Constants/behavior.h"
#include "../../../src/Chemistry/chemistry_enumerators.h"
#include "../../../src/FileManagement/file_enumerators.h"
#include "../../../src/MolecularMechanics/minimization.h"
#include "../../../src/MoleculeFormat/mdlmol.h"
#include "../../../src/MoleculeFormat/mdlmol_refinement.h"
#include "../../../src/MoleculeFormat/molecule_format_enumerators.h"
#include "../../../src/Namelists/nml_minimize.h"
#include "../../../src/Namelists/user_settings.h"
#include "../../../src/Reporting/error_format.h"
#include "../../../src/Reporting/help_messages.h"
#include "../../../src/Synthesis/phasespace_synthesis.h"
#include "../../../src/Synthesis/systemcache.h"
#include "../../../src/Topology/atomgraph_enumerators.h"
#include "../../../src/Trajectory/trajectory_enumerators.h"
#include "../../../src/Trajectory/phasespace.h"
#include "../../../src/UnitTesting/stopwatch.h"
#include "../../../src/UnitTesting/unit_test.h"

using namespace stormm::constants;
using namespace stormm::chemistry;
using namespace stormm::diskutil;
using namespace stormm::display;
using namespace stormm::energy;
using namespace stormm::errors;
using namespace stormm::mm;
using namespace stormm::namelist;
using namespace stormm::restraints;
using namespace stormm::structure;
using namespace stormm::synthesis;
using namespace stormm::testing;
using namespace stormm::topology;
using namespace stormm::trajectory;

//-------------------------------------------------------------------------------------------------
// main
//-------------------------------------------------------------------------------------------------
int main(int argc, const char* argv[]) {

  // Wall time tracking
  StopWatch timer("Timings for ffrefine.stormm");

  // Parse the command line
  CommandLineParser clip("ffrefine.stormm", "A program for performing energy minimizations with "
                         "experimental force fields.");
  clip.addStandardApplicationInputs();
  const std::vector<std::string> my_namelists = { "&files", "&minimize", "&restraint", "&report" };
  clip.addControlBlocks(my_namelists);
  if (displayNamelistHelp(argc, argv, my_namelists) && clip.doesProgramExitOnHelp()) {
    return 0;
  }
  clip.parseUserInput(argc, argv);

  // Read information from the command line and initialize the UserSettings object
  UserSettings ui(clip, { "-pe", "-ce" });
  
  // Read topologies and coordinate files.  Assemble critical details about each system.
  SystemCache sysc(ui.getFilesNamelistInfo(), ui.getRestraintNamelistInfo(),
                   ui.getDynamicsNamelistInfo(), ui.getExceptionBehavior(),
                   MapRotatableGroups::YES, ui.getPrintingPolicy(), &timer);
  const int system_count = sysc.getSystemCount();
  std::vector<MdlMol> sdf_recovery = sysc.getStructureDataEntry();
  customizeDataItems(&sdf_recovery, sysc, ui.getReportNamelistInfo());

  // Perform minimizations as requested.
  const MinimizeControls mincon = ui.getMinimizeNamelistInfo();
  const int mini_timings = timer.addCategory("Minimization");
  std::vector<ScoreCard> all_mme;
  if (ui.getMinimizePresence()) {
    all_mme.reserve(system_count);
    for (int i = 0; i < system_count; i++) {
      PhaseSpace *ps = sysc.getCoordinatePointer(i);
      const RestraintApparatus& ra = sysc.getRestraints(i);
      switch(ps->getUnitCellType()) {
      case UnitCellType::NONE:
        all_mme.emplace_back(minimize(ps, sysc.getSystemTopology(i), ra,
                                      sysc.getSystemStaticMask(i), mincon));
        timer.assignTime(mini_timings);
      case UnitCellType::ORTHORHOMBIC:
      case UnitCellType::TRICLINIC:
        break;
      }
    }
  }
  
  // Print restart files from energy minimization
  if (mincon.getCheckpointProduction()) {
    for (int i = 0; i < system_count; i++) {
      const PhaseSpace ps = sysc.getCoordinates(i);
      switch (sysc.getCheckpointKind(i)) {
      case CoordinateFileKind::AMBER_CRD:
      case CoordinateFileKind::AMBER_INPCRD:
      case CoordinateFileKind::AMBER_ASCII_RST:
      case CoordinateFileKind::AMBER_NETCDF:
      case CoordinateFileKind::AMBER_NETCDF_RST:
        ps.exportToFile(sysc.getCheckpointName(i), 0.0, TrajectoryKind::POSITIONS,
                        sysc.getCheckpointKind(i),
                        sysc.getPrintingProtocol(CoordinateFileRole::CHECKPOINT, i));
        break;
      case CoordinateFileKind::SDF:
        sdf_recovery[i].impartCoordinates(ps);
        updateDataItemReadouts(&sdf_recovery[i], sysc, all_mme[i]);
        sdf_recovery[i].writeMdl(sysc.getCheckpointName(i), MdlMolVersion::V2000,
                                 sysc.getPrintingProtocol(CoordinateFileRole::CHECKPOINT, i));
        sdf_recovery[i].writeDataItems(sysc.getCheckpointName(i), PrintSituation::APPEND);
        break;        
      case CoordinateFileKind::PDB:

        // CHECK
        rtErr("This functionality is not yet implemented.", "main");
        // END CHECK
        
        break;
      case CoordinateFileKind::UNKNOWN:
        break;
      }
    }
  }

  // Collate the topologies in preparation to operate on the new coordinate population
  timer.assignTime(0);
  timer.printResults();
  
  return 0;
}
