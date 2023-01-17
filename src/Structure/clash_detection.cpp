#include "copyright.h"
#include "clash_detection.h"

namespace stormm {
namespace structure {

//-------------------------------------------------------------------------------------------------
bool detectClash(const CoordinateFrameReader &cfr, const ValenceKit<double> &vk,
                 const NonbondedKit<double> &nbk, const StaticExclusionMaskReader &maskr,
                 const double elec_limit, const double vdw_ratio) {
  return detectClash<double, double>(cfr.xcrd, cfr.ycrd, cfr.zcrd, vk, nbk, maskr, elec_limit,
                                     vdw_ratio);
}

//-------------------------------------------------------------------------------------------------
bool detectClash(const CoordinateFrame &cf, const AtomGraph &ag,
                 const StaticExclusionMask &mask, const double elec_limit,
                 const double vdw_ratio) {
  return detectClash(cf.data(), ag.getDoublePrecisionValenceKit(),
                     ag.getDoublePrecisionNonbondedKit(), mask.data(), elec_limit, vdw_ratio);
}

//-------------------------------------------------------------------------------------------------
bool detectClash(const CoordinateFrame *cf, const AtomGraph *ag,
                 const StaticExclusionMask &mask, const double elec_limit,
                 const double vdw_ratio) {
  return detectClash(cf->data(), ag->getDoublePrecisionValenceKit(),
                     ag->getDoublePrecisionNonbondedKit(), mask.data(), elec_limit, vdw_ratio);
}

//-------------------------------------------------------------------------------------------------
bool detectClash(const PhaseSpaceReader &psr, const ValenceKit<double> &vk,
                 const NonbondedKit<double> &nbk, const StaticExclusionMaskReader &maskr,
                 const double elec_limit, const double vdw_ratio) {
  return detectClash<double, double>(psr.xcrd, psr.ycrd, psr.zcrd, vk, nbk, maskr, elec_limit,
                                     vdw_ratio);
}

//-------------------------------------------------------------------------------------------------
bool detectClash(const PhaseSpace &ps, const AtomGraph &ag, const StaticExclusionMask &mask,
                 const double elec_limit, const double vdw_ratio) {
  return detectClash(ps.data(), ag.getDoublePrecisionValenceKit(),
                     ag.getDoublePrecisionNonbondedKit(), mask.data(), elec_limit, vdw_ratio);
}

//-------------------------------------------------------------------------------------------------
bool detectClash(const PhaseSpace *ps, const AtomGraph *ag, const StaticExclusionMask &mask,
                 const double elec_limit, const double vdw_ratio) {
  return detectClash(ps->data(), ag->getDoublePrecisionValenceKit(),
                     ag->getDoublePrecisionNonbondedKit(), mask.data(), elec_limit, vdw_ratio);
}

} // namespace structure
} // namespace stormm
