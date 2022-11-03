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
// Display a general helpmessage for this program
//-------------------------------------------------------------------------------------------------
void displayGeneralHelpMessage() {
  const std::string base_msg =
    terminalFormat("A program for performing energy minimizations on many structures and then "
                   "editing parameters for cyclical refinement of force constants based on their "
                   "structural consequences.\n\n", "ffrefine");
  const std::string nml_msg =
    terminalFormat("Applicable namelists (re-run with one of these terms as the command-line "
                   "argument, IN QUOTES, i.e. \"&files\" or '&files', for further details):\n"
                   "  - &files\n  - &restraint\n  - &minimize\n  - &report\n", nullptr, nullptr,
                   0, 2, 2);
  printf("%s", base_msg.c_str());
  printf("%s", nml_msg.c_str());
  printf("\n");
}

//-------------------------------------------------------------------------------------------------
// main
//-------------------------------------------------------------------------------------------------
int main(int argc, const char* argv[]) {

  // Wall time tracking
  StopWatch timer("Timings for ffrefine.omni");

  // Check for a help message
  if (detectHelpSignal(argc, argv)) {
    displayGeneralHelpMessage();
    return 0;
  }
  if (displayNamelistHelp(argc, argv, { "&files", "&restraint", "&minimize", "&report" })) {
    return 0;
  }
  
  // Read information from the command line and initialize the UserSettings object
  UserSettings ui(argc, argv, AppName::FFREFINE);
  
  // Read topologies and coordinate files.  Assemble critical details about each system.
  std::vector<MdlMol> sdf_recovery;
  SystemCache sysc(ui.getFilesNamelistInfo(), ui.getRestraintNamelistInfo(),
                   &sdf_recovery, ui.getExceptionBehavior(), MapRotatableGroups::YES,
                   ui.getPrintingPolicy(), &timer);
  customizeDataItems(&sdf_recovery, sysc, ui.getReportNamelistInfo());

  // Perform minimizations as requested.
  const int system_count = sysc.getSystemCount();
  const MinimizeControls mincon = ui.getMinimizeNamelistInfo();
  const int mini_timings = timer.addCategory("Minimization");
  std::vector<ScoreCard> all_mme;
  if (ui.getMinimizePresence()) {
    all_mme.reserve(system_count);
    for (int i = 0; i < system_count; i++) {
      PhaseSpace *ps = sysc.getCoordinatePointer(i);
      const RestraintApparatus& ra = sysc.getRestraintReference(i);
      switch(ps->getUnitCellType()) {
      case UnitCellType::NONE:
        all_mme.emplace_back(minimize(ps, sysc.getSystemTopologyReference(i), ra,
                                      sysc.getSystemStaticMaskReference(i), mincon));
        timer.assignTime(mini_timings);
      case UnitCellType::ORTHORHOMBIC:
      case UnitCellType::TRICLINIC:
        break;
      }
      
      // CHECK
#if 0
      printf("Energy progression in system %2d\n", i);
      printf("STEP    BOND    ANGLE   DIHEDRAL IMPROPER  ELEC 1-4  VDW 1-4    ELEC     VDW   "
             "  RESTRAINT    TOTAL  \n");
      printf("----  -------- -------- -------- --------  -------- --------  -------- --------"
             "  ---------  ---------\n");
      const std::vector<double> bond_nrg = all_mme[i].reportEnergyHistory(StateVariable::BOND, 0);
      const std::vector<double> angl_nrg = all_mme[i].reportEnergyHistory(StateVariable::ANGLE, 0);
      const std::vector<double> dihe_nrg =
        all_mme[i].reportEnergyHistory(StateVariable::PROPER_DIHEDRAL, 0);
      const std::vector<double> impr_nrg =
        all_mme[i].reportEnergyHistory(StateVariable::IMPROPER_DIHEDRAL, 0);
      const std::vector<double> qq14_nrg =
        all_mme[i].reportEnergyHistory(StateVariable::ELECTROSTATIC_ONE_FOUR, 0);
      const std::vector<double> lj14_nrg =
        all_mme[i].reportEnergyHistory(StateVariable::VDW_ONE_FOUR, 0);
      const std::vector<double> qqnb_nrg =
        all_mme[i].reportEnergyHistory(StateVariable::ELECTROSTATIC, 0);
      const std::vector<double> ljnb_nrg =
          all_mme[i].reportEnergyHistory(StateVariable::VDW, 0);
      const std::vector<double> rstr_nrg =
        all_mme[i].reportEnergyHistory(StateVariable::RESTRAINT, 0);
      all_mme[i].computePotentialEnergy();      
      all_mme[i].computeTotalEnergy();
      const std::vector<double> totl_nrg =
        all_mme[i].reportEnergyHistory(StateVariable::POTENTIAL_ENERGY, 0);        
      for (int j = 0; j < all_mme[i].getSampleSize(); j += 50) {
        printf("%4d  %8.4lf %8.4lf %8.4lf %8.4lf  %8.4lf %8.4lf  %8.4lf %8.4lf  %9.4lf  %9.4lf\n",
               j, bond_nrg[j], angl_nrg[j], dihe_nrg[j], impr_nrg[j], qq14_nrg[j], lj14_nrg[j],
               qqnb_nrg[j], ljnb_nrg[j], rstr_nrg[j], totl_nrg[j]);
      }
      const AtomGraph *iag_ptr = sysc.getSystemTopologyPointer(i);
      const ImplicitSolventModel i_ism = iag_ptr->getImplicitSolventModel();
      if (i >= 0) {
        printf("Label %2d = %s (%s)\n", i, sysc.getSystemLabel(i).c_str(),
               getImplicitSolventModelName(i_ism).c_str());
        PhaseSpaceWriter psw = ps->data();
        for (int i = 0; i < psw.natom; i++) {
          printf("  %9.4lf %9.4lf %9.4lf\n", psw.xcrd[i], psw.ycrd[i], psw.zcrd[i]);
        }
        printf("\n");
      }
#endif
      // END CHECK
    }
  }
  
  // Print restart files from energy minimization
  if (mincon.getCheckpointProduction()) {
    for (int i = 0; i < system_count; i++) {
      const PhaseSpace ps = sysc.getCoordinateReference(i);
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
