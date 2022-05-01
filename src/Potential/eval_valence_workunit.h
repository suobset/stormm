// -*-c++-*-
#ifndef OMNI_EVAL_VALENCE_WORKUNIT_H
#define OMNI_EVAL_VALENCE_WORKUNIT_H

#include "Restraints/restraint_apparatus.h"
#include "Synthesis/valence_workunit.h"
#include "Topology/atomgraph.h"
#include "Topology/atomgraph_abstracts.h"
#include "Trajectory/phasespace.h"
#include "Trajectory/coordinateframe.h"
#include "energy_enumerators.h"
#include "scorecard.h"

namespace omni {
namespace energy {

using restraints::RestraintApparatus;
using restraints::RestraintApparatusDpReader;
using synthesis::ValenceWorkUnit;
using synthesis::VwuTask;
using synthesis::VwuGoal;
using topology::AtomGraph;
using topology::UnitCellType;
using topology::NonbondedKit;
using topology::ValenceKit;
using topology::VirtualSiteKit;
using trajectory::CoordinateFrame;
using trajectory::PhaseSpace;

/// \brief Carry out the local evaluations need for each valence work unit using locally cached
///        data.  This routine makes no assumptions about const-ness of the coordinates to permit
///        maximum flexibility.
void localVwuEvaluation(const ValenceKit<double> vk, const VirtualSiteKit<double> vsk,
                        const NonbondedKit<double> nbk, const RestraintApparatusDpReader rar,
                        const double* sh_charges, const int* sh_lj_idx, double* sh_xcrd,
                        double* sh_ycrd, double* sh_zcrd, const double* umat, const double* invu,
                        UnitCellType unit_cell, double* sh_xfrc, double* sh_yfrc,
                        double* sh_zfrc, ScoreCard *ecard, const int sysid,
                        const ValenceWorkUnit &my_vwu, EvaluateForce eval_force,
                        VwuTask activity, VwuGoal purpose, int step_number);

/// \brief Initialize the appropriate energies in preparation for a loop over valence work units.
///
/// \param ecard     The energy tracking object
/// \param activity  The activity that the work unit will be performing (correspondence with the
///                  energetic state variables is handled internally)
/// \param sysid     Index of the system of interest, the one to initialize
void evalVwuInitEnergy(ScoreCard *ecard, const VwuTask activity, const int sysid);

/// \brief Evaluate force and energy-related tasks in a list of valence work units given a system
///        or synthesis of systems with one or more topologies and coordinate sets.
///
/// Overloaded:
///   - Evaluate for a single topology and coordinate set
///   - Evaluate for a synthesis of topologies and coordinate sets
///   - Accept pointers or references to the above objects, or in the case of a single system
///     offer the option of passing critical topology abstracts by value with pointers to C-style
///     arrays for coordinates and forces
///   - Evaluate forces, or energies alone
/// \{
void evalValenceWorkUnits(const ValenceKit<double> vk, const VirtualSiteKit<double> vsk,
                          const NonbondedKit<double> nbk, const RestraintApparatusDpReader rar,
                          double* xcrd, double* ycrd, double* zcrd, const double* umat,
                          const double* invu, UnitCellType unit_cell, double* xfrc, double* yfrc,
                          double* zfrc, ScoreCard *ecard, int sysid,
                          const std::vector<ValenceWorkUnit> &vwu_list,
                          EvaluateForce eval_force = EvaluateForce::NO,
                          VwuTask activity = VwuTask::ALL_TASKS,
                          VwuGoal goal = VwuGoal::ACCUMULATE_FORCES, int step_number = 0);

void evalValenceWorkUnits(const AtomGraph *ag, PhaseSpace *ps, const RestraintApparatus *ra,
                          ScoreCard *ecard, int sysid,
                          const std::vector<ValenceWorkUnit> &vwu_list,
                          EvaluateForce eval_force = EvaluateForce::NO,
                          VwuTask activity = VwuTask::ALL_TASKS,
                          VwuGoal goal = VwuGoal::ACCUMULATE_FORCES, int step_number = 0);

void evalValenceWorkUnits(const ValenceKit<double> vk, const VirtualSiteKit<double> vsk,
                          const NonbondedKit<double> nbk, const RestraintApparatusDpReader rar,
                          const double* xcrd, const double* ycrd, const double* zcrd,
                          const double* umat, const double* invu, UnitCellType unit_cell,
                          double* xfrc, double* yfrc, double* zfrc, ScoreCard *ecard,
                          const int sysid, const std::vector<ValenceWorkUnit> &vwu_list,
                          EvaluateForce eval_force, VwuTask activity, int step_number);

void evalValenceWorkUnits(const AtomGraph &ag, const PhaseSpace &ps, const RestraintApparatus &ra,
                          ScoreCard *ecard, int sysid,
                          const std::vector<ValenceWorkUnit> &vwu_list,
                          VwuTask activity = VwuTask::ALL_TASKS, int step_number = 0);
/// \}

} // namespace energy
} // namespace omni

#endif
