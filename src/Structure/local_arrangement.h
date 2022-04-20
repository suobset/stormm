// -*-c++-*-
#ifndef OMNI_LOCAL_ARRANGEMENT_H
#define OMNI_LOCAL_ARRANGEMENT_H

#include <cmath>
#include "Math/vector_ops.h"
#include "Reporting/error_format.h"
#include "Topology/atomgraph_enumerators.h"
#include "Trajectory/coordinateframe.h"
#include "Accelerator/hybrid.h"
#include "Trajectory/coordinateframe.h"
#include "Trajectory/phasespace.h"
#include "structure_enumerators.h"

namespace omni {
namespace structure {

using card::Hybrid;
using topology::UnitCellType;
using trajectory::CoordinateFrame;
using trajectory::CoordinateFrameReader;
using trajectory::PhaseSpace;
using trajectory::PhaseSpaceReader;

/// \brief Image a single value to a range.
///        
/// \param x      The value to re-image
/// \param range  Defines the interval over which to image the value
/// \param style  Determines whether to image the value to the interval [ 0, range ) or
///               [ -0.5 * range, 0.5 * range ).
double imageValue(double x, double range, ImagingMethod style);

/// \brief Image a coordinate (or distance between coordinates) to fit within the interval
///        [-0.5, 0.5], taking into account box specifications.
///
/// Overloaded:
///   - Operate on a single (x, y, z) tuple
///   - Operate on three C-style arrays of real numbers x, y, and z
///   - Operate on three std::vectors of real numbers x, y, and z
///   - Operate on three Hybrid arrays of real numbers x, y, and z
///
/// \param x          Cartesian x coordinate of one or more particles
/// \param y          Cartesian y coordinate of one or more particles
/// \param z          Cartesian z coordinate of one or more particles
/// \param length     Number of particles or coordinate tuples to re-image
/// \param umat       Transformation matrix to go into fractional coordinates
/// \param invu       Transformation matrix to go from fraction coordinates back to real space
/// \param unit_cell  Shape of the unit cell
/// \{
void imageCoordinates(double *x, double *y, double *z, const double* umat, const double* invu,
                      UnitCellType unit_cell, ImagingMethod style);

void imageCoordinates(double* x, double* y, double* z, const int length, const double* umat,
                      const double* invu, UnitCellType unit_cell, ImagingMethod style);

void imageCoordinates(std::vector<double> *x, std::vector<double> *y, std::vector<double> *z,
                      const double* umat, const double* invu, UnitCellType unit_cell,
                      ImagingMethod style);

void imageCoordinates(Hybrid<double> *x, Hybrid<double> *y, Hybrid<double> *z,
                      const double* umat, const double* invu, UnitCellType unit_cell,
                      ImagingMethod style);

void imageCoordinates(PhaseSpace *ps, ImagingMethod style);

void imageCoordinates(CoordinateFrame *cf, ImagingMethod style);
/// \}

/// \brief Compute the distance between two points in a coordinate set.  All forms of this
///        function work with double-precision coordinates on the host.
///
/// Overloaded:
///   - Accept two atom indices and raw pointers to the coordinates and box specifications
///   - Accept two atom indices and a CoordinateFrameReader object
///   - Accept two atom indices and a CoordinateFrame object
///   - Accept two atom indices and a PhaseSpace object
/// \{
double distance(int atom_i, int atom_j, const double* xcrd, const double* ycrd, const double* zcrd,
                const double* umat, const double* invu, UnitCellType unit_cell);
double distance(int atom_i, int atom_j, const CoordinateFrameReader &cfr);
double distance(int atom_i, int atom_j, const CoordinateFrame &cf);
double distance(int atom_i, int atom_j, const PhaseSpace &ps);
/// \}

/// \brief Compute the angle between three points in a coordinate set.  All forms of this
///        function work with double-precision coordinates on the host.  The result is returned in
///        radians.
///
/// Overloaded:
///   - Accept three atom indices and raw pointers to the coordinates and box specifications
///   - Accept three atom indices and a CoordinateFrameReader object
///   - Accept three atom indices and a CoordinateFrame object
///   - Accept three atom indices and a PhaseSpace object
/// \{
double angle(int atom_i, int atom_j, int atom_k, const double* xcrd, const double* ycrd,
             const double* zcrd, const double* umat, const double* invu, UnitCellType unit_cell);
double angle(int atom_i, int atom_j, int atom_k, const CoordinateFrameReader &cfr);
double angle(int atom_i, int atom_j, int atom_k, const CoordinateFrame &cf);
double angle(int atom_i, int atom_j, int atom_k, const PhaseSpace &ps);
/// \}

/// \brief Compute the dihedral angle between three points in a coordinate set.  All forms of this
///        function work with double-precision coordinates on the host.  The result is returned in
///        radians.
///
/// Overloaded:
///   - Accept four atom indices and raw pointers to the coordinates and box specifications
///   - Accept four atom indices and a CoordinateFrameReader object
///   - Accept four atom indices and a CoordinateFrame object
///   - Accept four atom indices and a PhaseSpace object
/// \{
double dihedral_angle(int atom_i, int atom_j, int atom_k, int atom_l, const double* xcrd,
                      const double* ycrd, const double* zcrd, const double* umat,
                      const double* invu, UnitCellType unit_cell);
double dihedral_angle(int atom_i, int atom_j, int atom_k, int atom_l,
                      const CoordinateFrameReader &cfr);
double dihedral_angle(int atom_i, int atom_j, int atom_k, int atom_l, const CoordinateFrame &cf);
double dihedral_angle(int atom_i, int atom_j, int atom_k, int atom_l, const PhaseSpace &ps);
/// \}

} // namespace geometry
} // namespace omni

#endif
