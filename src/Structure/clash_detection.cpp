#include "copyright.h"
#include "clash_detection.h"

namespace stormm {
namespace structure {

//-------------------------------------------------------------------------------------------------
bool detectVanDerWaalsClash(const CoordinateFrameReader &cfr, const NonbondedKit<double> &nbk,
                            const StaticExclusionMaskReader &maskr, const double ratio) {
  return detectVanDerWaalsClash<double, double>(cfr.xcrd, cfr.ycrd, cfr.zcrd, nbk, maskr, ratio);
}

//-------------------------------------------------------------------------------------------------
bool detectVanDerWaalsClash(const CoordinateFrame &cf, const AtomGraph &ag,
                            const StaticExclusionMask &mask, const double ratio) {
  return detectVanDerWaalsClash(cf.data(), ag.getDoublePrecisionNonbondedKit(), mask.data(),
                                ratio);
}

//-------------------------------------------------------------------------------------------------
bool detectVanDerWaalsClash(const CoordinateFrame *cf, const AtomGraph *ag,
                            const StaticExclusionMask &mask, const double ratio) {
  return detectVanDerWaalsClash(cf->data(), ag->getDoublePrecisionNonbondedKit(), mask.data(),
                                ratio);
}

//-------------------------------------------------------------------------------------------------
bool detectVanDerWaalsClash(const PhaseSpaceReader &psr, const NonbondedKit<double> &nbk,
                            const StaticExclusionMaskReader &maskr, const double ratio) {
  return detectVanDerWaalsClash<double, double>(psr.xcrd, psr.ycrd, psr.zcrd, nbk, maskr, ratio);
}

//-------------------------------------------------------------------------------------------------
bool detectVanDerWaalsClash(const PhaseSpace &ps, const AtomGraph &ag,
                            const StaticExclusionMask &mask, const double ratio) {
  return detectVanDerWaalsClash(ps.data(), ag.getDoublePrecisionNonbondedKit(), mask.data(),
                                ratio);
}

//-------------------------------------------------------------------------------------------------
bool detectVanDerWaalsClash(const PhaseSpace *ps, const AtomGraph *ag,
                            const StaticExclusionMask &mask, const double ratio) {
  return detectVanDerWaalsClash(ps->data(), ag->getDoublePrecisionNonbondedKit(), mask.data(),
                                ratio);
}

} // namespace structure
} // namespace stormm
