#include "copyright.h"
#include "write_annotated_frame.h"

namespace stormm {
namespace trajectory {

//-------------------------------------------------------------------------------------------------
void writeFrame(std::ofstream *foutp, const std::string &filename, const TextFile &tf,
                const PhaseSpace &ps, const int total_width, const int fraction,
                const std::vector<int> &crd_lines, int x_start, int y_start, int z_start) {
}

//-------------------------------------------------------------------------------------------------
void writeFrame(std::ofstream *foutp, const std::string &filename, const TextFile &tf,
                const CoordinateFrame &cf, const int total_width, const int fraction,
                const std::vector<std::vector<int>> &crd_lines, int x_start, int y_start,
                int z_start) {
}
  
} // namespace trajectory
} // namespace stormm
