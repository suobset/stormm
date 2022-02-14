#include <cstdlib>
#include "../../../src/Namelists/nml_minimize.h"
#include "../../../src/Namelists/nml_random.h"
#include "../../../src/Namelists/nml_solvent.h"
#include "../../../src/Parsing/parse.h"
#include "../../../src/Parsing/polynumeric.h"
#include "../../../src/Reporting/error_format.h"
#include "user_settings.h"

namespace conf_app {
namespace user_input {

using omni::errors::rtErr;
using omni::errors::rtWarn;
using omni::namelist::extractImplicitSolventModel;
using omni::parse::NumberFormat;
using omni::parse::verifyNumberFormat;

//-------------------------------------------------------------------------------------------------
MoleculeSystem::MoleculeSystem() :
    topology_file_name{}, coordinate_file_name{}, frame_start{0}, frame_end{0}, replica_count{0},
    coordinate_kind{CoordinateFileKind::UNKNOWN}
{}

//-------------------------------------------------------------------------------------------------
MoleculeSystem::MoleculeSystem(const std::string &topology_file_in,
                               const std::string &coordinate_file_in, const int frame_start_in,
                               const int frame_end_in, const int replica_count_in,
                               const CoordinateFileKind coordinate_kind_in) :
    topology_file_name{topology_file_in},
    coordinate_file_name{coordinate_file_in},
    frame_start{frame_start_in},
    frame_end{(replica_count_in > 1) ? frame_start : frame_end_in},
    replica_count{replica_count_in},
    coordinate_kind{coordinate_kind_in}
{}

//-------------------------------------------------------------------------------------------------
UserSettings::UserSettings(const int argc, const char* argv[]) :
  input_file{std::string("cgen.in")},
  topology_file_names{},
  coordinate_file_names{},
  systems{},
  common_core_mask{std::string("")},
  report_file{std::string("cgen.out")},
  all_free_trajectory_frames{false},
  system_count{0},
  igseed{omni::namelist::default_random_seed},
  random_stream_count{omni::namelist::default_random_streams},
  minimization_steps{omni::namelist::default_minimize_maxcyc},
  steepest_descent_steps{omni::namelist::default_minimize_ncyc},
  initial_move_length{omni::namelist::default_minimize_dx0},
  convergence_criterion{omni::namelist::default_minimize_drms},
  igb{extractImplicitSolventModel(omni::namelist::default_solvent_igb)},
  born_radii_cutoff{omni::namelist::default_solvent_rgbmax},
  dielectric{omni::namelist::default_solvent_extdiel},
  salt_concentration{omni::namelist::default_solvent_saltcon},
  topology_cache{},
  initial_coordinates_cache{},
  topology_indices{}
{
  // Detect command line arguments, and note that their presence overrides similar directives
  // in the input deck.
  bool cli_igseed = false;
  for (int i = 1; i < argc; i++) {
    if (i < argc - 1 && strcmp(argv[i], "-igseed") == 0) {
      if (verifyNumberFormat(argv[i + 1], NumberFormat::INTEGER) == false) {
        rtErr("The random seed " + std::string(argv[i + 1]) + " is not valid.", "UserSettings");
      }
      igseed = atoi(argv[i + 1]);
    }
  }

  // Process the input file

}

//-------------------------------------------------------------------------------------------------
int UserSettings::getSystemCount() const {
  return system_count;
}

} // namespace user_input
} // namespace conf_app
