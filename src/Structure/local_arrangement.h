// -*-c++-*-
#ifndef OMNI_LOCAL_ARRANGEMENT_H
#define OMNI_LOCAL_ARRANGEMENT_H

#include <cmath>
#include "DataTypes/common_types.h"
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
using data_types::isSignedIntegralScalarType;
using math::crossProduct;
using topology::UnitCellType;
using trajectory::CoordinateFrame;
using trajectory::CoordinateFrameReader;
using trajectory::CoordinateFrameWriter;
using trajectory::PhaseSpace;
using trajectory::PhaseSpaceReader;

/// \brief Image a single value to a range.
///
/// \param x                  The value to re-image
/// \param range              Defines the interval over which to image the value
/// \param style              Determines whether to image the value to the interval [ 0, range )
///                           or [ -0.5 * range, 0.5 * range ).
/// \param gpos_scale_factor  Conversion factor for fixed precision coordinates
template <typename Tcoord, typename Tcalc>
Tcoord imageValue(Tcoord x, Tcalc range, ImagingMethod style, Tcalc gpos_scale_factor = 1.0);

/// \brief Image a coordinate (or distance between coordinates) to fit within the interval
///        [-0.5, 0.5], taking into account box specifications.
///
/// Overloaded:
///   - Operate on a single (x, y, z) tuple
///   - Operate on three C-style arrays of real numbers x, y, and z
///   - Operate on three std::vectors of real numbers x, y, and z
///   - Operate on three Hybrid arrays of real numbers x, y, and z
///
/// \param x                  Cartesian x coordinate of one or more particles
/// \param y                  Cartesian y coordinate of one or more particles
/// \param z                  Cartesian z coordinate of one or more particles
/// \param length             Number of particles or coordinate tuples to re-image
/// \param umat               Transformation matrix to go into fractional coordinates
/// \param invu               Transformation matrix from fractional coordinates back to real space
/// \param unit_cell          Shape of the unit cell
/// \param gpos_scale_factor  Conversion factor for fixed precision coordinates
/// \{
template <typename Tcoord, typename Tcalc>
void imageCoordinates(Tcoord *x, Tcoord *y, Tcoord *z, const Tcalc* umat, const Tcalc* invu,
                      UnitCellType unit_cell, ImagingMethod style, Tcalc gpos_scale_factor = 1.0);

template <typename Tcoord, typename Tcalc>
void imageCoordinates(Tcoord* x, Tcoord* y, Tcoord* z, const int length, const Tcalc* umat,
                      const Tcalc* invu, UnitCellType unit_cell, ImagingMethod style,
                      Tcalc gpos_scale_factor = 1.0);

template <typename Tcoord, typename Tcalc>
void imageCoordinates(std::vector<Tcoord> *x, std::vector<Tcoord> *y, std::vector<Tcoord> *z,
                      const Tcalc* umat, const Tcalc* invu, UnitCellType unit_cell,
                      ImagingMethod style, Tcalc gpos_scale_factor = 1.0);

template <typename Tcoord, typename Tcalc>
void imageCoordinates(Hybrid<Tcoord> *x, Hybrid<Tcoord> *y, Hybrid<Tcoord> *z,
                      const Tcalc* umat, const Tcalc* invu, UnitCellType unit_cell,
                      ImagingMethod style, Tcalc gpos_scale_factor = 1.0);

void imageCoordinates(PhaseSpace *ps, ImagingMethod style);

void imageCoordinates(CoordinateFrame *cf, ImagingMethod style);
/// \}

/// \brief Compute the distance between two points in a coordinate set.  All forms of this
///        function work with double-precision coordinates on the host.
///
/// Overloaded:
///   - Accept two atom indices and raw pointers to templated coordinate arrays (float, double, or
///     fixed precision signed integers) and templated box specifications (float or double, which
///     will determine the numerical precision of the internal calculations)
///   - Accept two atom indices and any of the coordinate storage objects
///
/// \param atom_i             Topological index of the first atom in the system
/// \param atom_j             Topological index of the second atom in the system
/// \param xcrd               Cartesian x coordinate of all particles in the system
/// \param ycrd               Cartesian y coordinate of all particles in the system
/// \param zcrd               Cartesian z coordinate of all particles in the system
/// \param umat               Transformation matrix to go into fractional coordinates
/// \param invu               Transformation matrix from fractional coordinates back to real space
/// \param unit_cell          Shape of the unit cell
/// \param gpos_scale_factor  Conversion factor for fixed precision coordinates
/// \{
template <typename Tcoord, typename Tcalc>
Tcalc distance(int atom_i, int atom_j, const Tcoord* xcrd, const Tcoord* ycrd, const Tcoord* zcrd,
               const Tcalc* umat, const Tcalc* invu, UnitCellType unit_cell,
               Tcalc gpos_scale_factor = 1.0);

double distance(int atom_i, int atom_j, const CoordinateFrameReader &cfr);

double distance(int atom_i, int atom_j, const CoordinateFrame &cf);

double distance(int atom_i, int atom_j, const PhaseSpace &ps);
/// \}

/// \brief Compute the angle between three points in a coordinate set.  All forms of this
///        function work with double-precision coordinates on the host.  The result is returned in
///        radians.
///
/// Overloaded:
///   - Accept three atom indices and raw pointers to templated coordinate arrays (float, double,
///     or fixed precision signed integers) and templated box specifications (float or double,
///     which will determine the numerical precision of the internal calculations)
///   - Accept three atom indices and any of the coordinate storage objects
///
/// Parameters for this function follow from descriptions in distance, above
/// \{
template <typename Tcoord, typename Tcalc>
Tcalc angle(int atom_i, int atom_j, int atom_k, const Tcoord* xcrd, const Tcoord* ycrd,
            const Tcoord* zcrd, const Tcalc* umat, const Tcalc* invu, UnitCellType unit_cell,
            Tcalc gpos_scale_factor = 1.0);

double angle(int atom_i, int atom_j, int atom_k, const CoordinateFrameReader &cfr);

double angle(int atom_i, int atom_j, int atom_k, const CoordinateFrame &cf);

double angle(int atom_i, int atom_j, int atom_k, const PhaseSpace &ps);
/// \}

/// \brief Compute the dihedral angle between three points in a coordinate set.  All forms of this
///        function work with double-precision coordinates on the host.  The result is returned in
///        radians.
///
/// Overloaded:
///   - Accept four atom indices and raw pointers to templated coordinate arrays (float, double,
///     or fixed precision signed integers) and templated box specifications (float or double,
///     which will determine the numerical precision of the internal calculations)
///   - Accept four atom indices and any of the coordinate storage objects
///
/// Parameters for this function follow from descriptions in distance, above
/// \{
template <typename Tcoord, typename Tcalc>
Tcalc dihedral_angle(int atom_i, int atom_j, int atom_k, int atom_l, const Tcoord* xcrd,
                     const Tcoord* ycrd, const Tcoord* zcrd, const Tcalc* umat,
                     const Tcalc* invu, UnitCellType unit_cell, Tcalc gpos_scale_factor = 1.0);

double dihedral_angle(int atom_i, int atom_j, int atom_k, int atom_l,
                      const CoordinateFrameReader &cfr);

double dihedral_angle(int atom_i, int atom_j, int atom_k, int atom_l, const CoordinateFrame &cf);

double dihedral_angle(int atom_i, int atom_j, int atom_k, int atom_l, const PhaseSpace &ps);
/// \}

} // namespace structure
} // namespace omni

#include "local_arrangement.tpp"

#endif
