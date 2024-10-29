// -*-c++-*-
#ifndef STORMM_GLOBAL_MANIPULATION_H
#define STORMM_GLOBAL_MANIPULATION_H

#include <vector>
#include "copyright.h"
#include "Constants/behavior.h"
#include "DataTypes/common_types.h"
#include "Reporting/error_format.h"
#include "Synthesis/condensate.h"
#include "Synthesis/phasespace_synthesis.h"
#include "Topology/atomgraph.h"
#include "Topology/atomgraph_abstracts.h"
#include "Trajectory/coordinateframe.h"
#include "Trajectory/coordinate_series.h"
#include "Trajectory/phasespace.h"
#include "Trajectory/trajectory_enumerators.h"
#include "virtual_site_handling.h"

namespace stormm {
namespace structure {

using constants::PrecisionModel;
using data_types::isSignedIntegralScalarType;
using data_types::isSignedIntegralHpcVectorType;
using synthesis::Condensate;
using synthesis::CondensateWriter;
using synthesis::PhaseSpaceSynthesis;
using synthesis::PsSynthesisWriter;
using topology::AtomGraph;
using topology::VirtualSiteKit;
using trajectory::CoordinateFrame;
using trajectory::CoordinateFrameWriter;
using trajectory::CoordinateSeries;
using trajectory::CoordinateSeriesWriter;
using trajectory::PhaseSpace;
using trajectory::PhaseSpaceWriter;
using trajectory::TrajectoryKind;

/// \brief Compute an unbiased rotation matrix gien orthogonal rotations about each of three lab
///        frame axes.  The treatment is taken as a replacement to extrinsic Euler angle
///        rotations.  Its derivation can be found in:
///
///        Daniel A. Beard and Tamar Schlick. "Unbiased Rotational Moves for Rigid-Body Dynamics."
///        Biophysical Journal 85:2973-2976.  2003.
///
/// \param om_x  Rotation about the lab frame X axis (also can be the Euler angle alpha), valid
///              range [ -pi, pi )
/// \param om_y  Rotation about the lab frame Y axis (also can be the Euler angle beta), valid
///              range [ 0, pi )
/// \param om_z  Rotation about the lab frame Z axis (also can be the Euler angle gamma), valid
///              range [ -pi, pi )
template <typename T>
std::vector<T> beardRotationMatrix(const T om_x, const T om_y, const T om_z);

/// \brief Rotate coordinates or a subset of coordinates.  Some overloaded forms of this function
///        can accept a VirtualSite kit abstract that will direct the repositioning of virtual
///        sites, in case one or more frames were subject to partial rotations.  More general forms
///        of the function work with rotation matrices or triads of unit vectors describing the
///        new coordinate axes.
///
/// Overloaded:
///   - Accept raw pointers with templated types
///   - Accept a pointer to a PhaseSpace or CoordinateFrame object, or abstracts thereof
///   - Rotate coordinates of a point expressed as three scalar numbers
///   - Rotate coordinates of many points expressed as a triad of vectors: C-style arrays, Standard
///     template library vectors, or Hybrid objects
///   - Rotate a three-tuple of coordinates
///   - Provide the rotation as a triad of vectors or a 3 x 3 rotation matrix representing the
///     unit cell axes of the new coordinate system
///
/// \param x       First component of the particle's or particles' coordinates
/// \param y       Second component of the particle's or particles' coordinates
/// \param z       Third component of the particle's or particles' coordinates
/// \param n       The number of coordinates, trusted length of x, y, and z when rotating multiple
///                coordinates at once
/// \param c       Tuple of coordinates with the X, Y, and Z components in the "x", "y", and "z"
///                members, respectively
/// \param axis_a  Unit vector describing the first dimension of the new coordinate system
/// \param axis_b  Unit vector describing the second dimension of the new coordinate system
/// \param axis_c  Unit vector describing the third dimension of the new coordinate system
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

template <typename T, typename Tcalc>
void rotateCoordinates(T *x, T *y, T *z, const Tcalc* axis_a, const Tcalc* axis_b,
                       const Tcalc* axis_c, Tcalc globalpos_scale_factor = 1.0);

template <typename T, typename Tcalc>
void rotateCoordinates(T *x, T *y, T *z, const Tcalc* umat, Tcalc globalpos_scale_factor = 1.0);

template <typename T, typename Tcalc>
void rotateCoordinates(T* x, T* y, T* z, size_t n, const Tcalc* axis_a, const Tcalc* axis_b,
                       const Tcalc* axis_c, Tcalc globalpos_scale_factor = 1.0);

template <typename T, typename Tcalc>
void rotateCoordinates(T* x, T* y, T* z, size_t n, const Tcalc* umat,
                       Tcalc globalpos_scale_factor = 1.0);

template <typename T, typename Tcalc>
void rotateCoordinates(std::vector<T> *x, std::vector<T> *y, std::vector<T> *z,
                       const std::vector<Tcalc> &axis_a, const std::vector<Tcalc> &axis_b,
                       const std::vector<Tcalc> &axis_c, Tcalc globalpos_scale_factor = 1.0);

template <typename T, typename Tcalc>
void rotateCoordinates(std::vector<T> *x, std::vector<T> *y, std::vector<T> *z,
                       const std::vector<Tcalc> &umat, Tcalc globalpos_scale_factor = 1.0);

template <typename T, typename Tcalc>
void rotateCoordinates(Hybrid<T> *x, Hybrid<T> *y, Hybrid<T> *z,
                       const Hybrid<Tcalc> &axis_a, const Hybrid<Tcalc> &axis_b,
                       const Hybrid<Tcalc> &axis_c, Tcalc globalpos_scale_factor = 1.0);

template <typename T, typename Tcalc>
void rotateCoordinates(Hybrid<T> *x, Hybrid<T> *y, Hybrid<T> *z, const Hybrid<Tcalc> &umat,
                       Tcalc globalpos_scale_factor = 1.0);

template <typename T3, typename Tcalc>
void rotateCoordinates(T3 *c, const T3 axis_a, const T3 axis_b, const T3 axis_c,
                       Tcalc globalpos_scale_factor = 1.0);

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

template <typename Tcoord, typename Tcalc>
void rotateCoordinates(CoordinateSeriesWriter<Tcoord> *csw, size_t frame_index,
                       const VirtualSiteKit<Tcalc> &vsk, Tcalc alpha, Tcalc beta, Tcalc gamma,
                       int lower_limit = 0, int upper_limit = 0);

template <typename Tcoord, typename Tcalc>
void rotateCoordinates(CoordinateSeries<Tcoord> *cs, size_t frame_index, const AtomGraph &ag,
                       Tcalc alpha, Tcalc beta, Tcalc gamma, int lower_limit = 0,
                       int upper_limit = 0);
                       
template <typename Tcalc>
void rotateCoordinates(PhaseSpaceSynthesis *poly_ps, int system_index, Tcalc alpha, Tcalc beta,
                       Tcalc gamma, int lower_limit = 0, int upper_limit = 0);

template <typename Tcalc>
void rotateCoordinates(CondensateWriter *cdnsw, int system_index, const VirtualSiteKit<Tcalc> &vsk,
                       Tcalc alpha, Tcalc beta, Tcalc gamma, int lower_limit = 0,
                       int upper_limit = 0);

void rotateCoordinates(Condensate *cdns, int system_index, double alpha, double beta,
                       double gamma, int lower_limit = 0, int upper_limit = 0);

void rotateCoordinates(Condensate *cdns, int system_index, const AtomGraph &ag, double alpha,
                       double beta, double gamma, int lower_limit = 0, int upper_limit = 0);
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

template <typename Tcoord, typename Tcalc>
void translateCoordinates(CoordinateSeriesWriter<Tcoord> *csw, Tcalc xmove, Tcalc ymove,
                          Tcalc zmove, int lower_limit = 0, int upper_limit = 0,
                          const VirtualSiteKit<Tcalc> *vsk = nullptr);

template <typename Tcoord, typename Tcalc>
void translateCoordinates(CoordinateSeries<Tcoord> *cs, const AtomGraph &ag, Tcalc xmove,
                          Tcalc ymove, Tcalc zmove, int lower_limit = 0, int upper_limit = 0);

template <typename Tcalc>
void translateCoordinates(PsSynthesisWriter *poly_psw, int system_index,
                          const VirtualSiteKit<Tcalc> &vsk, Tcalc xmove, Tcalc ymove,
                          Tcalc zmove, int lower_limit = 0, int upper_limit = 0);

template <typename Tcalc>
void translateCoordinates(PhaseSpaceSynthesis *poly_ps, int system_index, const AtomGraph &ag,
                          Tcalc xmove, Tcalc ymove, Tcalc zmove, int lower_limit = 0,
                          int upper_limit = 0);

template <typename Tcalc>
void translateCoordinates(CondensateWriter *cdnsw, int system_index,
                          const VirtualSiteKit<Tcalc> &vsk, Tcalc xmove, Tcalc ymove,
                          Tcalc zmove, int lower_limit = 0, int upper_limit = 0);

void translateCoordinates(Condensate *cdns, int system_index, double xmove, double ymove,
                          double zmove, int lower_limit = 0, int upper_limit = 0);

void translateCoordinates(Condensate *cdns, int system_index, const AtomGraph &ag,
                          double xmove, double ymove, double zmove, int lower_limit = 0,
                          int upper_limit = 0);
/// \}

} // namespace structure
} // namespace stormm

#include "global_manipulation.tpp"

#endif
