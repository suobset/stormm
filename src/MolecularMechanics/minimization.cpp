#include "minimization.h"

namespace stormm {
namespace mm {

using energy::EvaluateForce;
using math::invertSquareMatrix;
using math::matrixVectorMultiply;
using structure::placeVirtualSites;
using structure::transmitVirtualSiteForces;
  
//-------------------------------------------------------------------------------------------------
ScoreCard minimize(PhaseSpace *ps, const AtomGraph &ag, const RestraintApparatus &ra,
                   const StaticExclusionMask &se, const MinimizeControls &mincon,
                   const int nrg_scale_bits) {
  const NonbondedKit<double> nbk = ag.getDoublePrecisionNonbondedKit();
  const ValenceKit<double> vk = ag.getDoublePrecisionValenceKit();
  const VirtualSiteKit<double> vsk = ag.getDoublePrecisionVirtualSiteKit();
  PhaseSpaceWriter psw = ps->data();
  return minimize<double, double,
                  double, double2, double4>(psw.xcrd, psw.ycrd, psw.zcrd, psw.xfrc, psw.yfrc,
                                            psw.zfrc, psw.xvel, psw.yvel, psw.zvel, psw.xprv,
                                            psw.yprv, psw.zprv, vk, nbk,
                                            ra.getDoublePrecisionAbstract(), vsk, se.data(),
                                            mincon, nrg_scale_bits);
}

//-------------------------------------------------------------------------------------------------
ScoreCard minimize(PhaseSpace *ps, const AtomGraph *ag, const StaticExclusionMask &se,
                   const MinimizeControls &mincon, const int nrg_scale_bits) {
  const NonbondedKit<double> nbk = ag->getDoublePrecisionNonbondedKit();
  const ValenceKit<double> vk = ag->getDoublePrecisionValenceKit();
  const VirtualSiteKit<double> vsk = ag->getDoublePrecisionVirtualSiteKit();
  const RestraintApparatus ra(ag);
  PhaseSpaceWriter psw = ps->data();
  return minimize<double, double,
                  double, double2, double4>(psw.xcrd, psw.ycrd, psw.zcrd, psw.xfrc, psw.yfrc,
                                            psw.zfrc, psw.xvel, psw.yvel, psw.zvel, psw.xprv,
                                            psw.yprv, psw.zprv, vk, nbk,
                                            ra.getDoublePrecisionAbstract(), vsk, se.data(),
                                            mincon, nrg_scale_bits);
}

//-------------------------------------------------------------------------------------------------
ScoreCard minimize(PhaseSpaceWriter psw, const ValenceKit<double> &vk,
                   const NonbondedKit<double> &nbk,
                   const RestraintKit<double, double2, double4> &rar,
                   const VirtualSiteKit<double> &vsk, const StaticExclusionMaskReader &ser,
                   const MinimizeControls &mincon, int nrg_scale_bits) {
  return minimize<double, double,
                  double, double2, double4>(psw.xcrd, psw.ycrd, psw.zcrd, psw.xfrc, psw.yfrc,
                                            psw.zfrc, psw.xvel, psw.yvel, psw.zvel, psw.xprv,
                                            psw.yprv, psw.zprv, vk, nbk, rar, vsk, ser, mincon,
                                            nrg_scale_bits);
}
  
} // namespace mm
} // namespace stormm
