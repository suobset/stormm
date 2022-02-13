#include "user_settings.h"

namespace conf_app {
namespace user_input {

//-------------------------------------------------------------------------------------------------
MoleculeSystem::MoleculeSystem() :
    topology_file_name{}, coordinate_file_name{}, frame_start{0}, frame_end{0},
    coordinate_kind{CoordinateFileKind::UNKNOWN}
{}

//-------------------------------------------------------------------------------------------------
MoleculeSystem::MoleculeSystem(const std::string &topology_file_in,
                               const std::string &coordinate_file_in, const int frame_start_in,
                               const int frame_end_in,
                               const CoordinateFileKind coordinate_kind_in) :
    topology_file_name{topology_file_in},
    coordinate_file_name{coordinate_file_in},
    frame_start{frame_start_in},
    frame_end{frame_end_in},
    coordinate_kind{coordinate_kind_in}
{}

} // namespace user_input
} // namespace conf_app
