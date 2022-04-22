// -*-c++-*-
#ifndef OMNI_MM_EVALUATION_H
#define OMNI_MM_EVALUATION_H

#include "Potential/energy_enumerators.h"
#include "Potential/scorecard.h"
#include "Potential/static_exclusionmask.h"
#include "Restraints/restraint_apparatus.h"
#include "Topology/atomgraph.h"
#include "Topology/atomgraph_abstracts.h"
#include "Topology/atomgraph_enumerators.h"
#include "Trajectory/coordinateframe.h"
#include "Trajectory/phasespace.h"

namespace omni {
namespace mm {

using energy::EvaluateForce;
using energy::ScoreCard;
using energy::StaticExclusionMask;
using restraints::RestraintApparatusDpReader;
using topology::AtomGraph;
using topology::NonbondedKit;
using topology::UnitCellType;
using topology::ValenceKit;
using topology::VirtualSiteKit;
using trajectory::CoordinateFrame;
using trajectory::CoordinateFrameReader;
using trajectory::PhaseSpace;
using trajectory::PhaseSpaceWriter;
  
/// \brief Evaluate molecular mechanics valence energies and forces in a single system.  Energy
///        results are collected in fixed precision in the ScoreCard object fed to this routine.
///
/// \param xcrd        Cartesian X positions of all particles
/// \param ycrd        Cartesian Y positions of all particles
/// \param zcrd        Cartesian Z positions of all particles
/// \param umat        Transformation matrix taking coordinates into fractional (unit cell) space
/// \param invu        Transformation matrix to take fractional coordinates back into real space
/// \param unit_cell   The type of simulation cell, to guide re-imaging considerations
/// \param xfrc        Forces acting on all particles in the Cartesian X direction (accumulated
///                    and returned)
/// \param yfrc        Forces acting on all particles in the Cartesian Y direction
/// \param zfrc        Forces acting on all particles in the Cartesian Z direction
/// \param sc          Object to hold the resulting energies (also keeps running sums)
/// \param vk          Valence parameters abstract from the original topology
/// \param nbk         Non-bonded parameter abstract from the original topology
/// \param eval_force  Directive to have forces evaluated or not
void evalValeMM(double* xcrd, double* ycrd, double* zcrd, double* umat, double* invu,
                UnitCellType unit_cell, double* xfrc, double* yfrc, double* zfrc, ScoreCard *sc,
                const ValenceKit<double> &vk, const NonbondedKit<double> &nbk,
                EvaluateForce eval_force);

/// \brief Evaluate molecular mechanics energies and forces due to valence interactions and
///        NMR restraints in a single system.  Energy results are collected in fixed precision in
///        the ScoreCard object fed to this routine.  Descriptions of input variables follow from
///        evalValeMM() above, with the additional inclusion of:
///
/// \param rar   Restraint apparatus abstract guiding NMR restraint potentials, supplemental to
///              the system's topology
/// \param step  The step number, which may affect restraint activation
void evalValeRestMM(double* xcrd, double* ycrd, double* zcrd, double* umat, double* invu,
                    UnitCellType unit_cell, double* xfrc, double* yfrc, double* zfrc,
                    ScoreCard *sc, const ValenceKit<double> &vk, const NonbondedKit<double> &nbk,
                    const RestraintApparatusDpReader &rar, EvaluateForce eval_force, int step);

/// \brief Evaluate molecular mechanics energies and forces due to valence interactions, NMR
///        restraints, and non-bonded interactions in a single system.  Energy results are
///        collected in fixed precision in the ScoreCard object fed to this routine.  Descriptions
///        of input variables follow from evalValeMM() and evalValeRestMM() above.
void evalNonbValeRestMM(double* xcrd, double* ycrd, double* zcrd, double* umat, double* invu,
                        UnitCellType unit_cell, double* xfrc, double* yfrc, double* zfrc,
                        ScoreCard *sc, const ValenceKit<double> &vk,
                        const NonbondedKit<double> &nbk, const StaticExclusionMask &se,
                        const RestraintApparatusDpReader &rar, EvaluateForce eval_force, int step);

} // namespace mm
} // namespace omni

#endif
