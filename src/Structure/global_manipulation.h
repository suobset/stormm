// -*-c++-*-
#ifndef STORMM_GLOBAL_MANIPULATION_H
#define STORMM_GLOBAL_MANIPULATION_H

#include <vector>
#include "copyright.h"
#include "DataTypes/common_types.h"
#include "Topology/atomgraph.h"
#include "Topology/atomgraph_abstracts.h"
#include "Trajectory/coordinateframe.h"
#include "Trajectory/phasespace.h"
#include "virtual_site_handling.h"

namespace stormm {
namespace structure {

using data_types::isSignedIntegralScalarType;
using topology::AtomGraph;
using topology::VirtualSiteKit;
using trajectory::CoordinateFrame;
using trajectory::CoordinateFrameWriter;
using trajectory::PhaseSpace;
using trajectory::PhaseSpaceWriter;

/// \brief Compute an unbiased rotation matrix gien orthogonal rotations about each of three lab
///        frame axes.  The treatment is taken as a replacement to extrinsic Euler angle
///        rotations.  Its derivation can be found in:
///
///        Daniel A. Beard and Tamar Schlick. "Unbiased Rotational Moves for Rigid-Body Dynamics."
///        Biophysical Journal 85:2973-2976.  2003.
///
/// \param om_x  Rotation about the lab frame X axis (also can be the Euler angle alpha)
/// \param om_y  Rotation about the lab frame Y axis (also can be the Euler angle beta)
/// \param om_z  Rotation about the lab frame Z axis (also can be the Euler angle gamma)
template <typename T>
std::vector<T> beardRotationMatrix(const T om_x, const T om_y, const T om_z);

/// \brief Rotate coordinates or a subset of coordinates.  All overloaded forms of this function
///        can accept a VirtualSite kit abstract that will direct the repositioning of virtual
///        sites, in case one or more frames were subject to partial rotations.
///
/// Overloaded:
///   - Accept raw pointers with templated types
///   - Accept a pointer to a PhaseSpace or CoordinateFrame object, or abstracts thereof
///
/// \param xcrd         Cartesian X coordinates of the particles
/// \param ycrd         Cartesian Y coordinates of the particles
/// \param zcrd         Cartesian Z coordinates of the particles
/// \param cfw          Coordinates of all particles
/// \param cf           Coordinates of all particles
/// \param psw          Coordinates of all particles
/// \param ps           Coordinates of all particles
/// \param alpha        The first Euler rotation angle, about the lab frame x axis
/// \param beta         The second Euler rotation angle, about the lab frame y axis
/// \param gamma        The third Euler rotation angle, about the lab frame z axis
/// \param lower_limit  The lower limit of atoms to rotate
/// \param upper_limit  The upper limit of atoms to rotate
/// \param vsk          Topology abstract governing virtual site placement
/// \param ag           System topology (contains information on virtual site placement)
/// \{
template <typename Tcoord, typename Tcalc>
void rotateCoordinates(Tcoord* xcrd, Tcoord* ycrd, Tcoord* zcrd, Tcalc alpha, Tcalc beta,
                       Tcalc gamma, int lower_limit, int upper_limit,
                       const VirtualSiteKit<Tcalc> *vsk = nullptr,
                       Tcalc globalpos_scale_factor = 1.0);

void rotateCoordinates(CoordinateFrameWriter *cfw, const VirtualSiteKit<double> &vsk, double alpha,
                       double beta, double gamma, int lower_limit = 0, int upper_limit = 0);

void rotateCoordinates(CoordinateFrame *cf, const AtomGraph &ag, double alpha, double beta,
                       double gamma, int lower_limit = 0, int upper_limit = 0);

void rotateCoordinates(PhaseSpaceWriter *psw, const VirtualSiteKit<double> &vsk, double alpha,
                       double beta, double gamma, int lower_limit = 0, int upper_limit = 0);

void rotateCoordinates(PhaseSpace *ps, const AtomGraph &ag, double alpha, double beta,
                       double gamma, int lower_limit = 0, int upper_limit = 0);
/// \}

/// \brief Translate coordinates or a subset of coordinates.  All overloaded forms of this function
///        can accept a VirtualSite kit abstract that will direct the repositioning of virtual
///        sites, in case one or more frames were subject to partial rotations.
///
/// Overloaded:
///   - Accept raw pointers with templated types
///   - Accept a pointer to a PhaseSpace or CoordinateFrame object, or abstracts thereof
///
/// \param xcrd         Cartesian X coordinates of the particles
/// \param ycrd         Cartesian Y coordinates of the particles
/// \param zcrd         Cartesian Z coordinates of the particles
/// \param cfw          Coordinates of all particles
/// \param cf           Coordinates of all particles
/// \param psw          Coordinates of all particles
/// \param ps           Coordinates of all particles
/// \param xmove        Translation to impart along the Cartesian x axis
/// \param ymove        Translation to impart along the Cartesian y axis
/// \param zmove        Translation to impart along the Cartesian z axis
/// \param lower_limit  The lower limit of atoms to rotate
/// \param upper_limit  The upper limit of atoms to rotate
/// \param vsk          Topology abstract governing virtual site placement
/// \param ag           System topology (contains information on virtual site placement)
/// \{
template <typename Tcoord, typename Tcalc>
void translateCoordinates(Tcoord* xcrd, Tcoord* ycrd, Tcoord* zcrd, Tcalc xmove, Tcalc ymove,
                          Tcalc zmove, int lower_limit, int upper_limit,
                          const VirtualSiteKit<Tcalc> *vsk = nullptr,
                          Tcalc globalpos_scale_factor = 1.0);

void translateCoordinates(CoordinateFrameWriter *cfw, const VirtualSiteKit<double> &vsk,
                          double xmove, double ymove, double zmove, int lower_limit = 0,
                          int upper_limit = 0);

void translateCoordinates(CoordinateFrame *cf, const AtomGraph &ag, double xmove, double ymove,
                          double zmove, int lower_limit = 0, int upper_limit = 0);

void translateCoordinates(PhaseSpaceWriter *psw, const VirtualSiteKit<double> &vsk, double xmove,
                          double ymove, double zmove, int lower_limit = 0, int upper_limit = 0);

void translateCoordinates(PhaseSpace *ps, const AtomGraph &ag, double xmove, double ymove,
                          double zmove, int lower_limit = 0, int upper_limit = 0);
/// \}

} // namespace structure
} // namespace stormm

#include "global_manipulation.tpp"

#endif
