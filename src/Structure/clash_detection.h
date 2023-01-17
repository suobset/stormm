// -*-c++-*-
#ifndef STORMM_CLASH_DETECTION_H
#define STORMM_CLASH_DETECTION_H

#include "copyright.h"
#include "Constants/scaling.h"
#include "DataTypes/common_types.h"
#include "Math/rounding.h"
#include "Namelists/nml_minimize.h"
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
using data_types::isFloatingPointScalarType;
using energy::StaticExclusionMask;
using energy::StaticExclusionMaskReader;
using energy::supertile_length;
using energy::tile_length;
using energy::tile_lengths_per_supertile;
using math::roundUp;
using namelist::default_minimize_clash_ratio;
using namelist::default_minimize_clash_r0;
using topology::AtomGraph;
using topology::NonbondedKit;
using topology::ValenceKit;
using trajectory::CoordinateFrame;
using trajectory::CoordinateFrameReader;
using trajectory::CoordinateSeries;
using trajectory::CoordinateSeriesReader;
using trajectory::PhaseSpace;
using trajectory::PhaseSpaceReader;

/// \brief The minimum system size needed in order to engage a neighbor list based search.
constexpr int clash_direct_calculation_size_limit = 32;

/// \brief Compute the maximum distance at which two particles can be considered to participate in
///        a van-der Waals clash.
///
/// \param nbk         Non-bonded abstract for the system topology
/// \param elec_limit  The limiting distance at which two particles (with presumed electrostatic
///                    properties) will be considered to clash
/// \param vdw_ratio   The minimum required ratio of the distance between any pair of particles and
///                    mutual (non-bonded) sigma parameters
template <typename Tcalc>
Tcalc maxClashingDistance(const NonbondedKit<Tcalc> &nbk, const Tcalc elec_limit,
                          const Tcalc vdw_ratio);

/// \brief Implement a trivial test to see whether the (non-imaged) Cartesian ranges of two
///        cached sets of atom coordinates might overlap.  Return TRUE if there is a possibility,
///        FALSE if not.
///
/// \param cachi_xcrd  Cartesian X coordinates of the first set of particles
/// \param cachj_xcrd  Cartesian X coordinates of the second set of particles
/// \param cachi_ycrd  Cartesian Y coordinates of the first set of particles
/// \param cachj_ycrd  Cartesian Y coordinates of the second set of particles
/// \param cachi_zcrd  Cartesian Z coordinates of the first set of particles
/// \param cachj_zcrd  Cartesian Z coordinates of the second set of particles
/// \param ni_atoms    Number of atoms in the first set
/// \param nj_atoms    Number of atoms in the second set
/// \param max_clash   Maximum distance at which two particles might be deemed to clash, by either
///                    a raw distance or a van-der Waals sigma ratio criterion
template <typename Tcoord, typename Tcalc>
bool trivialClashCheck(const std::vector<Tcalc> &cachi_xcrd, const std::vector<Tcalc> &cachj_xcrd,
                       const std::vector<Tcalc> &cachi_ycrd, const std::vector<Tcalc> &cachj_ycrd,
                       const std::vector<Tcalc> &cachi_zcrd, const std::vector<Tcalc> &cachj_zcrd,
                       int ni_atoms, int nj_atoms, Tcalc max_clash);
  
/// \brief Perform a direct comparison of particle positions in a molecule to determine whether
///        there are any clashes.
///
/// \param xcrd        Cartesian X coordinates of all particles
/// \param ycrd        Cartesian Y coordinates of all particles
/// \param zcrd        Cartesian Z coordinates of all particles
/// \param nbk         Non-bonded abstract for the system topology
/// \param maskr       Read-only abstract for the exclusion mask
/// \param vdw_ratio   The minimum required ratio of the distance between any pair of particles and
///                    mutual (non-bonded) sigma parameters
/// \param elec_limit  The limiting distance at which two particles (with presumed electrostatic
///                    properties) will be considered to clash
/// \param inv_scale   Inverse scaling factor for converting coordinates into Angstrom units
template <typename Tcoord, typename Tcalc>
bool directClashTesting(const Tcoord* xcrd, const Tcoord* ycrd, const Tcoord* zcrd,
                        const ValenceKit<Tcalc> &vk, const NonbondedKit<Tcalc> &nbk,
                        const StaticExclusionMaskReader &maskr,
                        Tcalc elec_limit = default_minimize_clash_r0,
                        Tcalc vdw_ratio = default_minimize_clash_ratio, Tcalc inv_scale = 1.0);

/// \brief Detect a van-der Waals clash between particles based on a minimum required ratio against
///        any given pair's non-bonded sigma ratio.  Return TRUE if a clash is found, FALSE if not.
///
/// Overloaded:
///   - Provide various coordinate objects, by const pointer or by const reference, or their
///     read-only abstracts
///   - Provide the original topology and exclusion mask, or their abstracts
///
/// Parameter descriptions in this routine follow from directClashTesting() above, with the
/// additions of:
///
/// \param cf          The system coordinates
/// \param ps          The system coordinates (the current point in the time cycle will be queried)
/// \param cs          The system coordinates (a frame number is required)
/// \param frame       Frame index to examine, if the coordinate object provided is a
///                    CoordinateSeries or a Condensate
/// \param cfr         Read-only abstract of the coordinates
/// \param psr         Read-only abstract of the coordinates
/// \param ag          System topology, containing van-der Waals parameters and each atom's type
///                    index
/// \param mask        Static exclusion mask for the system, indicating exclusions in an all-to-all
///                    interaction context
/// \{
template <typename Tcoord, typename Tcalc>
bool detectClash(const Tcoord* xcrd, const Tcoord* ycrd, const Tcoord* zcrd,
                 const ValenceKit<Tcalc> &vk, const NonbondedKit<Tcalc> &nbk,
                 const StaticExclusionMaskReader &maskr,
                 Tcalc elec_limit = default_minimize_clash_r0,
                 Tcalc vdw_ratio = default_minimize_clash_ratio, Tcalc inv_scale = 1.0);

bool detectClash(const CoordinateFrameReader &cfr, const ValenceKit<double> &vk,
                 const NonbondedKit<double> &nbk, const StaticExclusionMaskReader &maskr,
                 double elec_limit = default_minimize_clash_r0,
                 double vdw_ratio = default_minimize_clash_ratio);

bool detectClash(const CoordinateFrame &cf, const AtomGraph &ag, const StaticExclusionMask &mask,
                 double elec_limit = default_minimize_clash_r0,
                 double vdw_ratio = default_minimize_clash_ratio);

bool detectClash(const CoordinateFrame *cf, const AtomGraph *ag, const StaticExclusionMask &mask,
                 double elec_limit = default_minimize_clash_r0,
                 double vdw_ratio = default_minimize_clash_ratio);

bool detectClash(const PhaseSpaceReader &psr, const ValenceKit<double> &vk,
                 const NonbondedKit<double> &nbk, const StaticExclusionMaskReader &maskr,
                 double elec_limit = default_minimize_clash_r0,
                 double vdw_ratio = default_minimize_clash_ratio);

bool detectClash(const PhaseSpace &ps, const AtomGraph &ag, const StaticExclusionMask &mask,
                 double elec_limit = default_minimize_clash_r0,
                 double vdw_ratio = default_minimize_clash_ratio);

bool detectClash(const PhaseSpace *ps, const AtomGraph *ag, const StaticExclusionMask &mask,
                 double elec_limit = default_minimize_clash_r0,
                 double vdw_ratio = default_minimize_clash_ratio);

template <typename Tcoord, typename Tcalc>
bool detectClash(const CoordinateSeriesReader<Tcoord> &csr, int frame, const ValenceKit<Tcalc> &vk,
                 const NonbondedKit<double> &nbk, const StaticExclusionMaskReader &maskr,
                 Tcalc elec_limit = default_minimize_clash_r0,
                 Tcalc vdw_ratio = default_minimize_clash_ratio, Tcalc inv_scale = 1.0);

template <typename Tcoord, typename Tcalc>
bool detectClash(const CoordinateSeries<Tcoord> *cs, int frame, const AtomGraph *ag,
                 const StaticExclusionMask &mask, Tcalc elec_limit = default_minimize_clash_r0,
                 Tcalc vdw_ratio = default_minimize_clash_ratio, Tcalc inv_scale = 1.0);

template <typename Tcoord, typename Tcalc>
bool detectClash(const CoordinateSeries<Tcoord> &cs, int frame, const AtomGraph &ag,
                 const StaticExclusionMask &mask, Tcalc elec_limit = default_minimize_clash_r0,
                 Tcalc vdw_ratio = default_minimize_clash_ratio, Tcalc inv_scale = 1.0);
/// \}

} // namespace structure
} // namespace stormm

#include "clash_detection.tpp"

#endif
