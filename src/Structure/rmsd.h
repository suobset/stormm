// -*-c++-*-
#ifndef OMNI_STRUCTURE_RMSD_H
#define OMNI_STRUCTURE_RMSD_H

#include "Topology/atomgraph.h"
#include "Topology/atomgraph_abstracts.h"
#include "Trajectory/coordinateframe.h"
#include "Trajectory/phasespace.h"
#include "structure_enumerators.h"

namespace omni {
namespace structure {

using topology::AtomGraph;
using topology::ChemicalDetailsKit;
using trajectory::CoordinateFrame;
using trajectory::CoordinateFrameReader;
using trajectory::PhaseSpace;
using trajectory::PhaseSpaceReader;
 
/// \brief Compute the positional RMSD between two sets of coordinates.
///
/// Overloaded:
///   - Use const pointers to coordinates and mass
///   - Accept PhaseSpace objects and the topology as const pointers or references, or take
///     the relevant abstracts
///   - Accept CoordinateFrame objects and the topology as const pointers or references, or take
///     the relevant abstracts
///   - Accept CoordinateSeries objects and the topology as const pointers or references, or
///     take the relevant abstracts, and return a Standard Template Library vector of results
///     for all frames of one series to all frames of the other
///
/// \param xcrd_a       Cartesian X coordinates of the first frame
/// \param xcrd_b       Cartesian X coordinates of the second frame
/// \param ycrd_a       Cartesian Y coordinates of the first frame
/// \param ycrd_b       Cartesian Y coordinates of the second frame
/// \param zcrd_a       Cartesian Z coordinates of the first frame
/// \param zcrd_b       Cartesian Z coordinates of the second frame
/// \param masses       Masses of all particles
/// \param method       Alignment and mass-weighting scheme for computing RMSD
/// \param psr_a        Coordinate abstract for the first frame
/// \param psr_b        Coordinate abstract for the second frame
/// \param cdk          Chemical details from the topology describing both frames (for masses)
/// \param cfr_a        Coordinate abstract for the first frame
/// \param cfr_b        Coordinate abstract for the second frame
/// \param ps_a         Coordinates for the first frame
/// \param ps_b         Coordinates for the second frame
/// \param cf_a         Coordinates for the first frame
/// \param cf_b         Coordinates for the second frame
/// \param ag           Topology for both frames (contains masses)
/// \param lower_limit  Lower limit of the atoms to cover in RMSD calculation
/// \param upper_limit  Upper limit of the atoms to cover in RMSD calculation
/// \{
double rmsd(const double* xcrd_a, const double* ycrd_a, const double* zcrd_a, const double* xcrd_b,
            const double* ycrd_b, const double* zcrd_b, const double* masses, RmsdMethod method,
            int lower_limit, int upper_limit);

double rmsd(const PhaseSpaceReader &psr_a, const PhaseSpaceReader &psr_b,
            const ChemicalDetailsKit &cdk, RmsdMethod method, int lower_limit, int upper_limit);

double rmsd(const PhaseSpace &ps_a, const PhaseSpace &ps_b, const AtomGraph &ag, RmsdMethod method,
            int lower_limit, int upper_limit);

double rmsd(const CoordinateFrameReader &cfr_a, const CoordinateFrameReader &cfr_b,
            const ChemicalDetailsKit &cdk, RmsdMethod method, int lower_limit, int upper_limit);

double rmsd(const CoordinateFrame &cf_a, const CoordinateFrame &cf_b, const AtomGraph &ag,
            RmsdMethod method, int lower_limit, int upper_limit);
/// \}

} // namespace structure
} // namespace omni

#endif

