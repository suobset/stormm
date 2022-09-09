#include <cmath>
#include "copyright.h"
#include "Math/vector_ops.h"
#include "Reporting/error_format.h"
#include "Topology/atomgraph_enumerators.h"
#include "Trajectory/coordinateframe.h"
#include "local_arrangement.h"

namespace stormm {
namespace structure {

using math::crossProduct;
using trajectory::CoordinateFrameWriter;
  
//-------------------------------------------------------------------------------------------------
void imageCoordinates(PhaseSpace *ps, const ImagingMethod style) {
  CoordinateFrameWriter cfw(ps);
  imageCoordinates(cfw.xcrd, cfw.ycrd, cfw.zcrd, cfw.natom, cfw.umat, cfw.invu, cfw.unit_cell,
                   style);
}

//-------------------------------------------------------------------------------------------------
void imageCoordinates(CoordinateFrame *cf, const ImagingMethod style) {
  CoordinateFrameWriter cfw = cf->data();
  imageCoordinates(cfw.xcrd, cfw.ycrd, cfw.zcrd, cfw.natom, cfw.umat, cfw.invu, cfw.unit_cell,
                   style);
}

//-------------------------------------------------------------------------------------------------
double distance(const int atom_i, const int atom_j, const CoordinateFrameReader &cfr) {
  return distance(atom_i, atom_j, cfr.xcrd, cfr.ycrd, cfr.zcrd, cfr.umat, cfr.invu, cfr.unit_cell);
}

//-------------------------------------------------------------------------------------------------
double distance(const int atom_i, const int atom_j, const CoordinateFrame &cf) {
  return distance(atom_i, atom_j, cf.data());
}

//-------------------------------------------------------------------------------------------------
double distance(const int atom_i, const int atom_j, const PhaseSpace &ps) {
  return distance(atom_i, atom_j, CoordinateFrameReader(ps));
}

//-------------------------------------------------------------------------------------------------
double angle(const int atom_i, const int atom_j, const int atom_k,
             const CoordinateFrameReader &cfr) {
  return angle(atom_i, atom_j, atom_k, cfr.xcrd, cfr.ycrd, cfr.zcrd, cfr.umat, cfr.invu,
               cfr.unit_cell);
}

//-------------------------------------------------------------------------------------------------
double angle(const int atom_i, const int atom_j, const int atom_k, const CoordinateFrame &cf) {
  return angle(atom_i, atom_j, atom_k, cf.data());
}

//-------------------------------------------------------------------------------------------------
double angle(const int atom_i, const int atom_j, const int atom_k, const PhaseSpace &ps) {
  return angle(atom_i, atom_j, atom_k, CoordinateFrameReader(ps));
}

//-------------------------------------------------------------------------------------------------
double dihedral_angle(int atom_i, int atom_j, int atom_k, int atom_l,
                      const CoordinateFrameReader &cfr) {
  return dihedral_angle(atom_i, atom_j, atom_k, atom_l, cfr.xcrd, cfr.ycrd, cfr.zcrd, cfr.umat,
                        cfr.invu, cfr.unit_cell);
}

//-------------------------------------------------------------------------------------------------
double dihedral_angle(const int atom_i, const int atom_j, const int atom_k, const int atom_l,
                      const CoordinateFrame &cf) {
  return dihedral_angle(atom_i, atom_j, atom_k, atom_l, cf.data());
}

//-------------------------------------------------------------------------------------------------
double dihedral_angle(const int atom_i, const int atom_j, const int atom_k, const int atom_l,
                      const PhaseSpace &ps) {
  return dihedral_angle(atom_i, atom_j, atom_k, atom_l, CoordinateFrameReader(ps));
}

} // namespace structure
} // namespace stormm
