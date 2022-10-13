#include "copyright.h"
#include "write_annotated_frame.h"

namespace stormm {
namespace trajectory {

//-------------------------------------------------------------------------------------------------
void writeFrame(std::ofstream *foutp, const std::string &filename, const TextFile &tf,
                const double* xcrd, const double* ycrd, const double* zcrd,
                const CoordinateSwapPlan &excision) {

}

//-------------------------------------------------------------------------------------------------
void writeFrame(std::ofstream *foutp, const std::string &filename, const TextFile &tf,
                const PhaseSpaceReader &psr, const CoordinateSwapPlan &excision) {
  writeFrame(foutp, filename, tf, psr.xcrd, psr.ycrd, psr.zcrd, excision);
}

//-------------------------------------------------------------------------------------------------
void writeFrame(std::ofstream *foutp, const std::string &filename, const TextFile &tf,
                const CoordinateFrameReader &cfr, const CoordinateSwapPlan &excision) {
  writeFrame(foutp, filename, tf, cfr.xcrd, cfr.ycrd, cfr.zcrd, excision);
}

//-------------------------------------------------------------------------------------------------
void writeFrame(std::ofstream *foutp, const std::string &filename, const TextFile &tf,
                const PhaseSpace &ps, const CoordinateSwapPlan &excision) {
  const PhaseSpaceReader psr = ps.data();
  writeFrame(foutp, filename, tf, psr, excision);
}

//-------------------------------------------------------------------------------------------------
void writeFrame(std::ofstream *foutp, const std::string &filename, const TextFile &tf,
                const CoordinateFrame &cf, const CoordinateSwapPlan &excision) {
  const CoordinateFrameReader cfr = cf.data();
  writeFrame(foutp, filename, tf, cfr, excision);
}

//-------------------------------------------------------------------------------------------------
void writeFrame(std::ofstream *foutp, const std::string &filename,
                const std::vector<TextFile> &tf_list, const PsSynthesisReader &poly_psr,
                const std::vector<CoordinateSwapPlan> &excision_list,
                const std::vector<int> &plan_indices) {
  for (int i = 0; i < poly_psr.system_count; i++) {
    const int natom  = poly_psr.atom_counts[i];
    const int offset = poly_psr.atom_starts[i];
    std::vector<double> xcrd(natom), ycrd(natom), zcrd(natom);
    for (int i = 0; i < natom; i++) {
      xcrd[i] = poly_psr.xcrd[offset + i];
      ycrd[i] = poly_psr.ycrd[offset + i];
      zcrd[i] = poly_psr.zcrd[offset + i];
    }
    writeFrame(foutp, filename, tf_list[i], xcrd.data(), ycrd.data(), zcrd.data(),
               excision_list[plan_indices[i]]);
  }
}

//-------------------------------------------------------------------------------------------------
void writeFrame(std::ofstream *foutp, const std::string &filename,
                const std::vector<TextFile> &tf_list, const PhaseSpaceSynthesis &poly_ps,
                const std::vector<CoordinateSwapPlan> &excision_list,
                const std::vector<int> &plan_indices) {
  const PsSynthesisReader poly_psr = poly_ps.data();
  writeFrame(foutp, filename, tf_list, poly_psr, excision_list, plan_indices);
}

} // namespace trajectory
} // namespace stormm
