#include "command.h"

namespace conf_app {
namespace user_input {

//-------------------------------------------------------------------------------------------------
MoleculeSystem::MoleculeSystem() :
  topology_file_name{}, coordinate_file_name{}, frame_start{0}, frame_end{0},
  coordinate_kind{CoordinateFileKind::UNKNOWN}
{}

} // namespace user_input
} // namespace conf_app
