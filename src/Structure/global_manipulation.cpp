#include "global_manipulation.h"
#include "structure_utils.h"
#include "virtual_site_handling.h"

namespace omni {
namespace structure {

using topology::UnitCellType;
  
//-------------------------------------------------------------------------------------------------
std::vector<double> beardRotationMatrix(const double om_x, const double om_y, const double om_z) {
  std::vector<double> result(9);

  // Convenient quantities
  const double om2_x = om_x*om_x;
  const double om2_y = om_y*om_y;
  const double om2_z = om_z*om_z;
  const double om2 = om2_x + om2_y + om2_z;
  const double om = sqrt(om2);
  const double inv_om = 1.0/om;
  const double inv_om2 = inv_om*inv_om;
  const double cos_om = cos(om);
  const double sin_om = sin(om);

  // Compute rotation matrix
  result[0] = (((om2_y + om2_z) * cos_om) + om2_x) * inv_om2;
  result[3] = ((om_x * om_y * inv_om2) * (1.0 - cos_om)) + ((om_z * inv_om) * sin_om);
  result[6] = ((om_x * om_z * inv_om2) * (1.0 - cos_om)) - ((om_y * inv_om) * sin_om);
  result[1] = ((om_x * om_y * inv_om2) * (1.0 - cos_om)) - ((om_z * inv_om) * sin_om);
  result[4] = (((om2_x + om2_z) * cos_om) + om2_y) * inv_om2;
  result[7] = ((om_y * om_z * inv_om2) * (1.0 - cos_om)) + ((om_x * inv_om) * sin_om);
  result[2] = ((om_x * om_z * inv_om2) * (1.0 - cos_om)) + ((om_y * inv_om) * sin_om);
  result[5] = ((om_y * om_z * inv_om2) * (1.0 - cos_om)) - ((om_x * inv_om) * sin_om);
  result[8] = (((om2_x + om2_y) * cos_om) + om2_z) * inv_om2;
  return result;
}

//-------------------------------------------------------------------------------------------------
void rotateCoordinates(double* xcrd, double* ycrd, double* zcrd, const double alpha,
                       const double beta, const double gamma, const int lower_limit,
                       const int upper_limit, const VirtualSiteKit<double> *vsk) {
  const std::vector<double> rmat = beardRotationMatrix(alpha, beta, gamma);
  for (int i = lower_limit; i < upper_limit; i++) {
    const double xc = xcrd[i];
    const double yc = ycrd[i];
    const double zc = zcrd[i];
    xcrd[i] = (rmat[0] * xc) + (rmat[3] * yc) + (rmat[6] * zc);
    ycrd[i] = (rmat[1] * xc) + (rmat[4] * yc) + (rmat[7] * zc);
    zcrd[i] = (rmat[2] * xc) + (rmat[5] * yc) + (rmat[8] * zc);
  }

  // Reposition virtual sites, if information is present to do so.  This assumes no unit cell
  // with periodic boundary conditions, as is approriate for coordinates undergoing global
  // rotation.  Dereferencing the pointer is not costly so long as the abstract contains only
  // constants and pointers.
  if (vsk != nullptr) {
    placeVirtualSites(xcrd, ycrd, zcrd, nullptr, nullptr, UnitCellType::NONE, *vsk);
  }
}

//-------------------------------------------------------------------------------------------------
void rotateCoordinates(CoordinateFrameWriter *cfw, const VirtualSiteKit<double> &vsk,
                       const double alpha, const double beta, const double gamma,
                       const int lower_limit, const int upper_limit) {
  coordinateBoundsCheck(lower_limit, upper_limit, cfw->natom, "rotateCoordinates");
  const int actual_upper_limit = (upper_limit <= lower_limit) ? cfw->natom : upper_limit;
  rotateCoordinates(cfw->xcrd, cfw->ycrd, cfw->zcrd, alpha, beta, gamma, lower_limit,
                    actual_upper_limit);
  placeVirtualSites(cfw->xcrd, cfw->ycrd, cfw->zcrd, nullptr, nullptr, UnitCellType::NONE, vsk);
}


//-------------------------------------------------------------------------------------------------
void rotateCoordinates(CoordinateFrame *cf, const AtomGraph &ag, const double alpha,
                       const double beta, const double gamma, const int lower_limit,
                       const int upper_limit) {
  CoordinateFrameWriter cfw = cf->data();
  coordinateBoundsCheck(lower_limit, upper_limit, cfw.natom, "rotateCoordinates");
  const int actual_upper_limit = (upper_limit <= lower_limit) ? cfw.natom : upper_limit;
  rotateCoordinates(cfw.xcrd, cfw.ycrd, cfw.zcrd, alpha, beta, gamma, lower_limit,
                    actual_upper_limit);
  placeVirtualSites(cfw.xcrd, cfw.ycrd, cfw.zcrd, nullptr, nullptr, UnitCellType::NONE,
                    ag.getDoublePrecisionVirtualSiteKit());
}

//-------------------------------------------------------------------------------------------------
void rotateCoordinates(PhaseSpaceWriter *psw, const VirtualSiteKit<double> &vsk,
                       const double alpha, const double beta, const double gamma,
                       const int lower_limit, const int upper_limit) {
  coordinateBoundsCheck(lower_limit, upper_limit, psw->natom, "rotateCoordinates");
  const int actual_upper_limit = (upper_limit <= lower_limit) ? psw->natom : upper_limit;
  rotateCoordinates(psw->xcrd, psw->ycrd, psw->zcrd, alpha, beta, gamma, lower_limit,
                    actual_upper_limit);
  placeVirtualSites(psw->xcrd, psw->ycrd, psw->zcrd, nullptr, nullptr, UnitCellType::NONE, vsk);
}

//-------------------------------------------------------------------------------------------------
void rotateCoordinates(PhaseSpace *ps, const AtomGraph &ag, const double alpha, const double beta,
                       const double gamma, const int lower_limit, const int upper_limit) {
  PhaseSpaceWriter psw = ps->data();
  coordinateBoundsCheck(lower_limit, upper_limit, psw.natom, "rotateCoordinates");
  const int actual_upper_limit = (upper_limit <= lower_limit) ? psw.natom : upper_limit;
  rotateCoordinates(psw.xcrd, psw.ycrd, psw.zcrd, alpha, beta, gamma, lower_limit,
                    actual_upper_limit);
  placeVirtualSites(psw.xcrd, psw.ycrd, psw.zcrd, nullptr, nullptr, UnitCellType::NONE,
                    ag.getDoublePrecisionVirtualSiteKit());
}

//-------------------------------------------------------------------------------------------------
void translateCoordinates(double* xcrd, double* ycrd, double* zcrd, const double xmove,
                          const double ymove, const double zmove, const int lower_limit,
                          const int upper_limit, const VirtualSiteKit<double> *vsk) {
  for (int i = lower_limit; i < upper_limit; i++) {
    xcrd[i] += xmove;
    ycrd[i] += ymove;
    zcrd[i] += zmove;
  }

  // Reposition virtual sites, if information is present to do so.  This assumes no unit cell
  // with periodic boundary conditions, as is approriate for coordinates undergoing global
  // rotation.  Dereferencing the pointer is not costly so long as the abstract contains only
  // constants and pointers.
  if (vsk != nullptr) {
    placeVirtualSites(xcrd, ycrd, zcrd, nullptr, nullptr, UnitCellType::NONE, *vsk);
  }
}

//-------------------------------------------------------------------------------------------------
void translateCoordinates(CoordinateFrameWriter *cfw, const VirtualSiteKit<double> &vsk,
                          const double xmove, const double ymove, const double zmove,
                          const int lower_limit, const int upper_limit) {
  coordinateBoundsCheck(lower_limit, upper_limit, cfw->natom, "translateCoordinates");
  const int actual_upper_limit = (upper_limit <= lower_limit) ? cfw->natom : upper_limit;
  translateCoordinates(cfw->xcrd, cfw->ycrd, cfw->zcrd, xmove, ymove, zmove, lower_limit,
                    actual_upper_limit);
  placeVirtualSites(cfw->xcrd, cfw->ycrd, cfw->zcrd, nullptr, nullptr, UnitCellType::NONE, vsk);
}


//-------------------------------------------------------------------------------------------------
void translateCoordinates(CoordinateFrame *cf, const AtomGraph &ag, const double xmove,
                          const double ymove, const double zmove, const int lower_limit,
                          const int upper_limit) {
  CoordinateFrameWriter cfw = cf->data();
  coordinateBoundsCheck(lower_limit, upper_limit, cfw.natom, "translateCoordinates");
  const int actual_upper_limit = (upper_limit <= lower_limit) ? cfw.natom : upper_limit;
  translateCoordinates(cfw.xcrd, cfw.ycrd, cfw.zcrd, xmove, ymove, zmove, lower_limit,
                    actual_upper_limit);
  placeVirtualSites(cfw.xcrd, cfw.ycrd, cfw.zcrd, nullptr, nullptr, UnitCellType::NONE,
                    ag.getDoublePrecisionVirtualSiteKit());
}

//-------------------------------------------------------------------------------------------------
void translateCoordinates(PhaseSpaceWriter *psw, const VirtualSiteKit<double> &vsk,
                          const double xmove, const double ymove, const double zmove,
                          const int lower_limit, const int upper_limit) {
  coordinateBoundsCheck(lower_limit, upper_limit, psw->natom, "translateCoordinates");
  const int actual_upper_limit = (upper_limit <= lower_limit) ? psw->natom : upper_limit;
  translateCoordinates(psw->xcrd, psw->ycrd, psw->zcrd, xmove, ymove, zmove, lower_limit,
                    actual_upper_limit);
  placeVirtualSites(psw->xcrd, psw->ycrd, psw->zcrd, nullptr, nullptr, UnitCellType::NONE, vsk);
}

//-------------------------------------------------------------------------------------------------
void translateCoordinates(PhaseSpace *ps, const AtomGraph &ag, const double xmove,
                          const double ymove, const double zmove, const int lower_limit,
                          const int upper_limit) {
  PhaseSpaceWriter psw = ps->data();
  coordinateBoundsCheck(lower_limit, upper_limit, psw.natom, "translateCoordinates");
  const int actual_upper_limit = (upper_limit <= lower_limit) ? psw.natom : upper_limit;
  translateCoordinates(psw.xcrd, psw.ycrd, psw.zcrd, xmove, ymove, zmove, lower_limit,
                    actual_upper_limit);
  placeVirtualSites(psw.xcrd, psw.ycrd, psw.zcrd, nullptr, nullptr, UnitCellType::NONE,
                    ag.getDoublePrecisionVirtualSiteKit());
}

} // namespace structure
} // namespace omni
