// -*-c++-*-
#ifndef OMNI_EVAL_VALENCE_WORKUNIT_H
#define OMNI_EVAL_VALENCE_WORKUNIT_H

#include "Synthesis/valence_workunit.h"
#include "Topology/atomgraph.h"
#include "Topology/atomgraph_abstracts.h"
#include "Trajectory/phasespace.h"
#include "Trajectory/coordinateframe.h"
#include "energy_enumerators.h"
#include "scorecard.h"

namespace omni {
namespace energy {

using synthesis::ValenceWorkUnit;
using synthesis::VwuTask;
using topology::AtomGraph;
using topology::UnitCellType;
using topology::NonbondedKit;
using topology::ValenceKit;
using topology::VirtualSiteKit;
using trajectory::CoordinateFrame;
using trajectory::PhaseSpace;
  
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
                          const NonbondedKit<double> nbk, const double* xcrd, const double* ycrd,
                          const double* zcrd, const double* umat, const double* invu,
                          UnitCellType unit_cell, double* xfrc, double* yfrc, double* zfrc,
                          ScoreCard *ecard, int sysid,
                          const std::vector<ValenceWorkUnit> &vwu_list,
                          EvaluateForce eval_force = EvaluateForce::NO,
                          VwuTask activity = VwuTask::ALL_TASKS);

void evalValenceWorkUnits(const AtomGraph *ag, PhaseSpace *ps, ScoreCard *ecard, int sysid,
                          const std::vector<ValenceWorkUnit> &vwu_list,
                          EvaluateForce eval_force = EvaluateForce::NO,
                          VwuTask activity = VwuTask::ALL_TASKS);

void evalValenceWorkUnits(const AtomGraph &ag, const PhaseSpace &ps, ScoreCard *ecard, int sysid,
                          const std::vector<ValenceWorkUnit> &vwu_list,
                          VwuTask activity = VwuTask::ALL_TASKS);
/// \}

} // namespace energy
} // namespace omni

#endif
