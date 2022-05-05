// -*-c++-*-
#ifndef OMNI_MM_EVALUATION_H
#define OMNI_MM_EVALUATION_H

#include "Potential/energy_enumerators.h"
#include "Potential/nonbonded_potential.h"
#include "Potential/scorecard.h"
#include "Potential/static_exclusionmask.h"
#include "Potential/valence_potential.h"
#include "Restraints/restraint_apparatus.h"
#include "Topology/atomgraph.h"
#include "Topology/atomgraph_abstracts.h"
#include "Topology/atomgraph_enumerators.h"
#include "Trajectory/coordinateframe.h"
#include "Trajectory/phasespace.h"

namespace omni {
namespace mm {

using energy::EvaluateForce;
using energy::evaluateBondTerms;
using energy::evaluateAngleTerms;
using energy::evaluateDihedralTerms;
using energy::evaluateUreyBradleyTerms;
using energy::evaluateCharmmImproperTerms;
using energy::evaluateCmapTerms;
using energy::evaluateAttenuated14Terms;
using energy::ScoreCard;
using energy::StaticExclusionMask;
using energy::StaticExclusionMaskReader;
using restraints::RestraintApparatus;
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
/// Overloaded:
///   - Take raw pointers to coordinates, box transformations, and forces, along with abstracts
///     to the required parameters
///   - Take a PhaseSpace object and the relevant topology, or abstracts thereof
///
/// \param xcrd          Cartesian X positions of all particles
/// \param ycrd          Cartesian Y positions of all particles
/// \param zcrd          Cartesian Z positions of all particles
/// \param umat          Transformation matrix taking coordinates into fractional (unit cell) space
/// \param invu          Transformation matrix to take fractional coordinates back into real space
/// \param unit_cell     The type of simulation cell, to guide re-imaging considerations
/// \param xfrc          Forces acting on all particles in the Cartesian X direction (accumulated
///                      and returned)
/// \param yfrc          Forces acting on all particles in the Cartesian Y direction
/// \param zfrc          Forces acting on all particles in the Cartesian Z direction
/// \param ps            PhaseSpace object with all coordinates, box information, and forces
/// \param psw           PhaseSpace object abstract with coordinates, box information, and forces
/// \param sc            Object to hold the resulting energies (also keeps running sums)
/// \param vk            Valence parameters abstract from the original topology
/// \param nbk           Non-bonded parameter abstract from the original topology
/// \param ag            System topology from which relevant parameter abstracts can be obtained
/// \param eval_force    Directive to have forces evaluated or not
/// \param system_index  Index of the system within the score card
/// \{
template <typename Tcoord, typename Tforce, typename Tcalc>
void evalValeMM(Tcoord* xcrd, Tcoord* ycrd, Tcoord* zcrd, double* umat, double* invu,
                UnitCellType unit_cell, Tforce* xfrc, Tforce* yfrc, Tforce* zfrc, ScoreCard *sc,
                const ValenceKit<Tcalc> &vk, const NonbondedKit<Tcalc> &nbk,
                EvaluateForce eval_force, int system_index = 0, Tcalc inv_gpos_factor = 1.0,
                Tcalc force_factor = 1.0);

void evalValeMM(PhaseSpaceWriter psw, ScoreCard *sc, const ValenceKit<double> &vk,
                const NonbondedKit<double> &nbk, EvaluateForce eval_force, int system_index = 0);

void evalValeMM(PhaseSpace *ps, ScoreCard *sc, const AtomGraph &ag, EvaluateForce eval_force,
                int system_index = 0);

void evalValeMM(PhaseSpace *ps, ScoreCard *sc, const AtomGraph *ag, EvaluateForce eval_force,
                int system_index = 0);
/// \}

/// \brief Evaluate molecular mechanics energies and forces due to valence interactions and
///        NMR restraints in a single system.  Energy results are collected in fixed precision in
///        the ScoreCard object fed to this routine.
///
/// Overloaded:
///   - Take raw pointers to coordinates, box transformations, and forces, along with abstracts
///     to the required parameters
///   - Take a PhaseSpace object and the relevant topology plus restraint collection, or abstracts
///     thereof
///
/// Descriptions of input variables follow from evalValeMM() above, with the addition of:
///
/// \param rar   Restraint apparatus abstract guiding NMR restraint potentials, supplemental to
///              the system's topology
/// \param ra    Restraint apparatus from which an appropriate abstract can be derived
/// \param step  The step number, which may affect restraint activation
/// \{
void evalValeRestMM(double* xcrd, double* ycrd, double* zcrd, double* umat, double* invu,
                    UnitCellType unit_cell, double* xfrc, double* yfrc, double* zfrc,
                    ScoreCard *sc, const ValenceKit<double> &vk, const NonbondedKit<double> &nbk,
                    const RestraintApparatusDpReader &rar, EvaluateForce eval_force,
                    int system_index = 0, int step = 0);

void evalValeRestMM(PhaseSpaceWriter psw, ScoreCard *sc, const ValenceKit<double> &vk,
                    const NonbondedKit<double> &nbk, const RestraintApparatusDpReader &rar,
                    EvaluateForce eval_force, int system_index = 0, int step = 0);

void evalValeRestMM(PhaseSpace *ps, ScoreCard *sc, const AtomGraph &ag,
                    const RestraintApparatus &ra, EvaluateForce eval_force, int system_index = 0,
                    int step = 0);

void evalValeRestMM(PhaseSpace *ps, ScoreCard *sc, const AtomGraph *ag,
                    const RestraintApparatus &ra, EvaluateForce eval_force, int system_index = 0,
                    int step = 0);
/// \}
  
/// \brief Evaluate the molecular mechanics energies and forces due to valence and non-bonded
///        pair interactions.  Energy results are collected in fixed precision in the ScoreCard
///        object fed to this routine.
///
/// Overloaded:
///   - Take raw pointers to coordinates, box transformations, and forces, along with abstracts
///     to the required parameters
///   - Take a PhaseSpace object and the relevant topology plus restraint collection, or abstracts
///     thereof
///
/// Descriptions of input variables follow from evalValeMM() and evalValeRestMM() above.
/// \{
void evalNonbValeMM(double* xcrd, double* ycrd, double* zcrd, double* umat, double* invu,
                    UnitCellType unit_cell, double* xfrc, double* yfrc, double* zfrc,
                    ScoreCard *sc, const ValenceKit<double> &vk,
                    const NonbondedKit<double> &nbk, const StaticExclusionMaskReader &ser,
                    EvaluateForce eval_force, int system_index = 0);

void evalNonbValeMM(PhaseSpaceWriter psw, ScoreCard *sc, const ValenceKit<double> &vk,
                    const NonbondedKit<double> &nbk, const StaticExclusionMaskReader &ser,
                    EvaluateForce eval_force, int system_index = 0);

void evalNonbValeMM(PhaseSpace *ps, ScoreCard *sc, const AtomGraph &ag,
                    const StaticExclusionMask &se, EvaluateForce eval_force, int system_index = 0);

void evalNonbValeMM(PhaseSpace *ps, ScoreCard *sc, const AtomGraph *ag,
                    const StaticExclusionMask &se, EvaluateForce eval_force, int system_index = 0);
/// \}
  
/// \brief Evaluate molecular mechanics energies and forces due to valence interactions, NMR
///        restraints, and non-bonded interactions in a single system.  Energy results are
///        collected in fixed precision in the ScoreCard object fed to this routine.
///
/// Overloaded:
///   - Take raw pointers to coordinates, box transformations, and forces, along with abstracts
///     to the required parameters
///   - Take a PhaseSpace object and the relevant topology plus restraint collection, or abstracts
///     thereof
///
/// Descriptions of input variables follow from evalValeMM() and evalValeRestMM() above, with the
/// addition of:
///
/// \param ser  Abstract for a map of non-bonded exclusions for all pairs of atoms in the system
/// \param se   Permanent map of non-bonded exclusions for all pairs of atoms in the system
/// \{
void evalNonbValeRestMM(double* xcrd, double* ycrd, double* zcrd, double* umat, double* invu,
                        UnitCellType unit_cell, double* xfrc, double* yfrc, double* zfrc,
                        ScoreCard *sc, const ValenceKit<double> &vk,
                        const NonbondedKit<double> &nbk, const StaticExclusionMaskReader &ser,
                        const RestraintApparatusDpReader &rar, EvaluateForce eval_force,
                        int system_index = 0, int step = 0);

void evalNonbValeRestMM(PhaseSpaceWriter psw, ScoreCard *sc, const ValenceKit<double> &vk,
                        const NonbondedKit<double> &nbk, const StaticExclusionMaskReader &ser,
                        const RestraintApparatusDpReader &rar, EvaluateForce eval_force,
                        int system_index = 0, int step = 0);

void evalNonbValeRestMM(PhaseSpace *ps, ScoreCard *sc, const AtomGraph &ag,
                        const StaticExclusionMask &se, const RestraintApparatus &ra,
                        EvaluateForce eval_force, int system_index = 0, int step = 0);

void evalNonbValeRestMM(PhaseSpace *ps, ScoreCard *sc, const AtomGraph *ag,
                        const StaticExclusionMask &se, const RestraintApparatus &ra,
                        EvaluateForce eval_force, int system_index = 0, int step = 0);
/// \}

} // namespace mm
} // namespace omni

#include "mm_evaluation.tpp"

#endif
