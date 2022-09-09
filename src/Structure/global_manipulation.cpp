#include "copyright.h"
#include "global_manipulation.h"
#include "structure_utils.h"

namespace stormm {
namespace structure {

using topology::UnitCellType;
  
//-------------------------------------------------------------------------------------------------
void rotateCoordinates(CoordinateFrameWriter *cfw, const VirtualSiteKit<double> &vsk,
                       const double alpha, const double beta, const double gamma,
                       const int lower_limit, const int upper_limit) {
  coordinateBoundsCheck(lower_limit, upper_limit, cfw->natom, "rotateCoordinates");
  const int actual_upper_limit = (upper_limit <= lower_limit) ? cfw->natom : upper_limit;
  rotateCoordinates<double, double>(cfw->xcrd, cfw->ycrd, cfw->zcrd, alpha, beta, gamma,
                                    lower_limit, actual_upper_limit);
  placeVirtualSites<double, double>(cfw->xcrd, cfw->ycrd, cfw->zcrd, nullptr, nullptr,
                                    UnitCellType::NONE, vsk);
}


//-------------------------------------------------------------------------------------------------
void rotateCoordinates(CoordinateFrame *cf, const AtomGraph &ag, const double alpha,
                       const double beta, const double gamma, const int lower_limit,
                       const int upper_limit) {
  CoordinateFrameWriter cfw = cf->data();
  coordinateBoundsCheck(lower_limit, upper_limit, cfw.natom, "rotateCoordinates");
  const int actual_upper_limit = (upper_limit <= lower_limit) ? cfw.natom : upper_limit;
  rotateCoordinates<double, double>(cfw.xcrd, cfw.ycrd, cfw.zcrd, alpha, beta, gamma, lower_limit,
                                    actual_upper_limit);
  placeVirtualSites<double, double>(cfw.xcrd, cfw.ycrd, cfw.zcrd, nullptr, nullptr,
                                    UnitCellType::NONE, ag.getDoublePrecisionVirtualSiteKit());
}

//-------------------------------------------------------------------------------------------------
void rotateCoordinates(PhaseSpaceWriter *psw, const VirtualSiteKit<double> &vsk,
                       const double alpha, const double beta, const double gamma,
                       const int lower_limit, const int upper_limit) {
  coordinateBoundsCheck(lower_limit, upper_limit, psw->natom, "rotateCoordinates");
  const int actual_upper_limit = (upper_limit <= lower_limit) ? psw->natom : upper_limit;
  rotateCoordinates<double, double>(psw->xcrd, psw->ycrd, psw->zcrd, alpha, beta, gamma,
                                    lower_limit, actual_upper_limit);
  placeVirtualSites<double, double>(psw->xcrd, psw->ycrd, psw->zcrd, nullptr, nullptr,
                                    UnitCellType::NONE, vsk);
}

//-------------------------------------------------------------------------------------------------
void rotateCoordinates(PhaseSpace *ps, const AtomGraph &ag, const double alpha, const double beta,
                       const double gamma, const int lower_limit, const int upper_limit) {
  PhaseSpaceWriter psw = ps->data();
  coordinateBoundsCheck(lower_limit, upper_limit, psw.natom, "rotateCoordinates");
  const int actual_upper_limit = (upper_limit <= lower_limit) ? psw.natom : upper_limit;
  rotateCoordinates<double, double>(psw.xcrd, psw.ycrd, psw.zcrd, alpha, beta, gamma, lower_limit,
                                    actual_upper_limit);
  placeVirtualSites<double, double>(psw.xcrd, psw.ycrd, psw.zcrd, nullptr, nullptr,
                                    UnitCellType::NONE, ag.getDoublePrecisionVirtualSiteKit());
}

//-------------------------------------------------------------------------------------------------
void translateCoordinates(CoordinateFrameWriter *cfw, const VirtualSiteKit<double> &vsk,
                          const double xmove, const double ymove, const double zmove,
                          const int lower_limit, const int upper_limit) {
  coordinateBoundsCheck(lower_limit, upper_limit, cfw->natom, "translateCoordinates");
  const int actual_upper_limit = (upper_limit <= lower_limit) ? cfw->natom : upper_limit;
  translateCoordinates<double, double>(cfw->xcrd, cfw->ycrd, cfw->zcrd, xmove, ymove, zmove,
                                       lower_limit, actual_upper_limit);
  placeVirtualSites<double, double>(cfw->xcrd, cfw->ycrd, cfw->zcrd, nullptr, nullptr,
                                    UnitCellType::NONE, vsk);
}


//-------------------------------------------------------------------------------------------------
void translateCoordinates(CoordinateFrame *cf, const AtomGraph &ag, const double xmove,
                          const double ymove, const double zmove, const int lower_limit,
                          const int upper_limit) {
  CoordinateFrameWriter cfw = cf->data();
  coordinateBoundsCheck(lower_limit, upper_limit, cfw.natom, "translateCoordinates");
  const int actual_upper_limit = (upper_limit <= lower_limit) ? cfw.natom : upper_limit;
  translateCoordinates<double, double>(cfw.xcrd, cfw.ycrd, cfw.zcrd, xmove, ymove, zmove,
                                       lower_limit, actual_upper_limit);
  placeVirtualSites<double, double>(cfw.xcrd, cfw.ycrd, cfw.zcrd, nullptr, nullptr,
                                    UnitCellType::NONE, ag.getDoublePrecisionVirtualSiteKit());
}

//-------------------------------------------------------------------------------------------------
void translateCoordinates(PhaseSpaceWriter *psw, const VirtualSiteKit<double> &vsk,
                          const double xmove, const double ymove, const double zmove,
                          const int lower_limit, const int upper_limit) {
  coordinateBoundsCheck(lower_limit, upper_limit, psw->natom, "translateCoordinates");
  const int actual_upper_limit = (upper_limit <= lower_limit) ? psw->natom : upper_limit;
  translateCoordinates<double, double>(psw->xcrd, psw->ycrd, psw->zcrd, xmove, ymove, zmove,
                                       lower_limit, actual_upper_limit);
  placeVirtualSites<double, double>(psw->xcrd, psw->ycrd, psw->zcrd, nullptr, nullptr,
                                    UnitCellType::NONE, vsk);
}

//-------------------------------------------------------------------------------------------------
void translateCoordinates(PhaseSpace *ps, const AtomGraph &ag, const double xmove,
                          const double ymove, const double zmove, const int lower_limit,
                          const int upper_limit) {
  PhaseSpaceWriter psw = ps->data();
  coordinateBoundsCheck(lower_limit, upper_limit, psw.natom, "translateCoordinates");
  const int actual_upper_limit = (upper_limit <= lower_limit) ? psw.natom : upper_limit;
  translateCoordinates<double, double>(psw.xcrd, psw.ycrd, psw.zcrd, xmove, ymove, zmove,
                                       lower_limit, actual_upper_limit);
  placeVirtualSites<double, double>(psw.xcrd, psw.ycrd, psw.zcrd, nullptr, nullptr,
                                    UnitCellType::NONE, ag.getDoublePrecisionVirtualSiteKit());
}

} // namespace structure
} // namespace stormm
