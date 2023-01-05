// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace structure {

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
bool detectVanDerWaalsClash(const Tcoord* xcrd, const Tcoord* ycrd, const Tcoord* zcrd,
                            const NonbondedKit<Tcalc> &nbk, const StaticExclusionMaskReader &maskr,
                            const Tcalc	ratio, const Tcalc inv_scale) {

  // If this point has been reached, no clash was detected.
  return false;
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
bool detectVanDerWaalsClash(const CoordinateSeriesReader<Tcoord> &csr, const size_t frame,
                            const NonbondedKit<Tcalc> &nbk, const StaticExclusionMaskReader &maskr,
                            Tcalc ratio) {
  const size_t padded_atoms = roundUp(csr.natom, warp_size_int);
  return detectVanDerWaalsClash(&csr.xcrd[frame * padded_atoms], &csr.ycrd[frame * padded_atoms],
                                &csr.zcrd[frame * padded_atoms], nbk, maskr, ratio,
                                csr.inv_gpos_scale);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
bool detectVanDerWaalsClash(const CoordinateSeries<Tcoord> *cs, const int frame,
                            const AtomGraph *ag, const StaticExclusionMask &mask, Tcalc ratio) {
  const size_t ct = std::type_index(typeid(Tcalc)).hash_code();
  if (ct == double_type_index) {
    return detectVanDerWaalsClash<Tcoord, double>(cs->data(), frame,
                                                  ag->getDoublePrecisionNonbondedKit(),
                                                  mask.data(), ratio);
  }
  else if (ct == float_type_index) {
    return detectVanDerWaalsClash<Tcoord, float>(cs->data(), frame,
                                                 ag->getSinglePrecisionNonbondedKit(),
                                                 mask.data(), ratio);
  }
  else {
    rtErr("The clash detection must be performed in double or float.  " +
          getStormmScalarTypeName<Tcalc> + " is invalid.", "detectVanDerWaalsClash");
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
bool detectVanDerWaalsClash(const CoordinateSeries<Tcoord> &cs, const int frame,
                            const AtomGraph &ag, const StaticExclusionMask &mask, Tcalc ratio) {
  return detectVanDerWaalsClash<Tcoord, Tcalc>(cs.getSelfPointer(), frame, ag.getSelfPointer(),
                                               mask, ratio);
}

} // namespace structure
} // namespace stormm
