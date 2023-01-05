// -*-c++-*-
#ifndef STORMM_CLASH_DETECTION_H
#define STORMM_CLASH_DETECTION_H

#include "copyright.h"
#include "DataTypes/common_types.h"
#include "Math/rounding.h"
#include "Potential/static_exclusionmask.h"
#include "Reporting/error_format.h"
#include "Topology/atomgraph.h"
#include "Topology/atomgraph_abstracts.h"
#include "Trajectory/coordinateframe.h"
#include "Trajectory/coordinate_series.h"
#include "Trajectory/phasespace.h"

namespace stormm {
namespace structure {

using data_types::getStormmScalarTypeName;
using energy::StaticExclusionMask;
using energy::StaticExclusionMaskReader;
using topology::AtomGraph;
using topology::NonbondedKit;
using trajectory::CoordinateFrame;
using trajectory::CoordinateFrameReader;
using trajectory::CoordinateSeries;
using trajectory::CoordinateSeriesReader;
using trajectory::PhaseSpace;
using trajectory::PhaseSpaceReader;
  
/// \brief Detect a van-der Waals clash between particles based on a minimum required ratio against
///        any given pair's non-bonded sigma ratio.  Return TRUE if a clash is found, FALSE if not.
///
/// Overloaded:
///   - Provide various coordinate objects, by const pointer or by const reference, or their
///     read-only abstracts
///   - Provide the original topology and exclusion mask, or their abstracts
///
/// \param xcrd       Cartesian X coordinates of all particles
/// \param ycrd       Cartesian Y coordinates of all particles
/// \param zcrd       Cartesian Z coordinates of all particles
/// \param cf         The system coordinates
/// \param ps         The system coordinates (the current point in the time cycle will be queried)
/// \param cs         The system coordinates (a frame number is required)
/// \param frame      Frame index to examine, if the coordinate object provided is a
///                   CoordinateSeries or a Condensate
/// \param cfr        Read-only abstract of the coordinates
/// \param psr        Read-only abstract of the coordinates
/// \param ag         System topology, containing van-der Waals parameters and each atom's type
///                   index
/// \param nbk        Non-bonded abstract for the system topology
/// \param mask       Static exclusion mask for the system, indicating exclusions in an all-to-all
///                   interaction context
/// \param maskr      Read-only abstract for the exclusion mask
/// \param ratio      The minimum required ratio of the distance between any pair of particles and
///                   mutual (non-bonded) sigma parameters
/// \param inv_scale  Inverse scaling factor for converting coordinates into Angstrom units
/// \{
template <typename Tcoord, typename Tcalc>
bool detectVanDerWaalsClash(const Tcoord* xcrd, const Tcoord* ycrd, const Tcoord* zcrd,
                            const NonbondedKit<Tcalc> &nbk, const StaticExclusionMaskReader &maskr,
                            Tcalc ratio = 0.5, Tcalc inv_scale = 1.0);

bool detectVanDerWaalsClash(const CoordinateFrameReader &cfr, const NonbondedKit<double> &nbk,
                            const StaticExclusionMaskReader &maskr, double ratio = 0.5);

bool detectVanDerWaalsClash(const CoordinateFrame &cf, const AtomGraph &ag,
                            const StaticExclusionMask &mask, double ratio = 0.5);

bool detectVanDerWaalsClash(const CoordinateFrame *cf, const AtomGraph *ag,
                            const StaticExclusionMask &mask, double ratio = 0.5);

bool detectVanDerWaalsClash(const PhaseSpaceReader &psr, const NonbondedKit<double> &nbk,
                            const StaticExclusionMaskReader &maskr, double ratio = 0.5);

bool detectVanDerWaalsClash(const PhaseSpace &ps, const AtomGraph &ag,
                            const StaticExclusionMask &mask, double ratio = 0.5);

bool detectVanDerWaalsClash(const PhaseSpace *ps, const AtomGraph *ag,
                            const StaticExclusionMask &mask, double ratio = 0.5);

template <typename Tcoord, typename Tcalc>
bool detectVanDerWaalsClash(const CoordinateSeriesReader<Tcoord> &csr, int frame,
                            const NonbondedKit<double> &nbk,
                            const StaticExclusionMaskReader &maskr, double ratio = 0.5);

template <typename Tcoord, typename Tcalc>
bool detectVanDerWaalsClash(const CoordinateSeries<Tcoord> *cs, int frame, const AtomGraph *ag,
                            const StaticExclusionMask &mask, Tcalc ratio = 0.5);

template <typename Tcoord, typename Tcalc>
bool detectVanDerWaalsClash(const CoordinateSeries<Tcoord> &cs, int frame, const AtomGraph &ag,
                            const StaticExclusionMask &mask, Tcalc ratio = 0.5);
/// \}

} // namespace structure
} // namespace stormm

#include "clash_detection.tpp"

#endif
