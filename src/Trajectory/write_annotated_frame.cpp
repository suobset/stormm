#include "copyright.h"
#include "write_annotated_frame.h"

namespace stormm {
namespace trajectory {

//-------------------------------------------------------------------------------------------------
void writeFrame(std::ofstream *foutp, const std::string &filename, const TextFile &tf,
                const double* xcrd, const double* ycrd, const double* zcrd, const int total_width,
                const int fraction, const std::vector<int> &crd_lines, const int x_start,
                const int y_start, const int z_start) {

}

//-------------------------------------------------------------------------------------------------
void writeFrame(std::ofstream *foutp, const std::string &filename, const TextFile &tf,
                const PhaseSpaceReader &psr, const int total_width, const int fraction,
                const std::vector<int> &crd_lines, int x_start, int y_start, int z_start) {
  writeFrame(foutp, filename, tf, psr.xcrd, psr.ycrd, psr.zcrd, total_width, fraction, crd_lines,
             x_start, y_start, z_start);
}

//-------------------------------------------------------------------------------------------------
void writeFrame(std::ofstream *foutp, const std::string &filename, const TextFile &tf,
                const CoordinateFrameReader &cfr, const int total_width, const int fraction,
                const std::vector<std::vector<int>> &crd_lines, const int x_start,
                const int y_start, const int z_start) {
  writeFrame(foutp, filename, tf, cfr.xcrd, cfr.ycrd, cfr.zcrd, total_width, fraction, crd_lines,
             x_start, y_start, z_start);
}

//-------------------------------------------------------------------------------------------------
void writeFrame(std::ofstream *foutp, const std::string &filename, const TextFile &tf,
                const PhaseSpace &ps, const int total_width, const int fraction,
                const std::vector<int> &crd_lines, int x_start, int y_start, int z_start) {
  const PhaseSpaceReader psr = ps.data();
  writeFrame(foutp, filename, tf, psr, total_width, fraction, crd_lines, x_start, y_start,
             z_start);
}

//-------------------------------------------------------------------------------------------------
void writeFrame(std::ofstream *foutp, const std::string &filename, const TextFile &tf,
                const CoordinateFrame &cf, const int total_width, const int fraction,
                const std::vector<std::vector<int>> &crd_lines, const int x_start,
                const int y_start, const int z_start) {
  const CoordinateFrameReader cfr = cf.data();
  writeFrame(foutp, filename, tf, psr, total_width, fraction, crd_lines, x_start, y_start,
             z_start);
}

//-------------------------------------------------------------------------------------------------
void writeFrame(std::ofstream *foutp, const std::string &filename, const TextFile &tf,
                const PhaseSpaceSynthesis &poly_ps, const int total_width, const int fraction,
                const std::vector<std::vector<int>> &crd_lines, const int x_start,
                const int y_start, const int z_start) {
  const PsSynthesisReader poly_psr = poly_ps.data();
  for (int i = 0; i < poly_psr.system_count; i++) {
    if (
  }
}

} // namespace trajectory
} // namespace stormm
