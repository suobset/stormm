// -*-c++-*-
#ifndef OMNI_NONBONDED_POTENTIAL_H
#define OMNI_NONBONDED_POTENTIAL_H

#include "Constants/generalized_born.h"
#include "Topology/atomgraph.h"
#include "Trajectory/phasespace.h"
#include "scorecard.h"
#include "static_exclusionmask.h"
#include "energy_enumerators.h"

namespace omni {
namespace energy {

using topology::AtomGraph;
using topology::ImplicitSolventKit;
using trajectory::PhaseSpace;
using namespace generalized_born_defaults;

/// \brief Evaluate the non-bonded electrostatic energy using an all-to-all approach with no
///        imaging considerations (isolated boundary conditions).
///
/// \param ag            System topology
/// \param se            Exclusion mask providing bits for all atom pairs
/// \param ps            Coordinates and force accumulators (modified by this function)
/// \param ecard         Energy components and other state variables (volume, temperature, etc.)
///                      (modified by this function)
/// \param eval_force    Flag to have forces also evaluated
/// \param system_index  Index of the system to which this energy contributes
/// \{
double2 evaluateNonbondedEnergy(const AtomGraph &ag, const StaticExclusionMask &se, PhaseSpace *ps,
                                ScoreCard *ecard,
                                EvaluateForce eval_elec_force = EvaluateForce::NO,
                                EvaluateForce eval_vdw_force  = EvaluateForce::NO,
                                int system_index = 0);
/// \}

/// \brief Evaluate the non-bonded Generalized Born energy of a system of particles (no cutoff is
///        applied in computation of the radii or evaluation of the energy and forces)
///
/// \param ag            System topology
/// \param ngb_tables    "Neck" Generalized Born tables from some pre-loaded cache of constants
/// \param ps            Coordinates and force accumulators (modified by this function)
/// \param ecard         Energy components and other state variables (volume, temperature, etc.)
///                      (modified by this function)
/// \param eval_force    Flag to have forces also evaluated
/// \param system_index  Index of the system to which this energy contributes
double evaluateGeneralizedBornEnergy(const AtomGraph &ag,
                                     const NeckGeneralizedBornTable &ngb_tables, PhaseSpace *ps,
                                     ScoreCard *ecard, const EvaluateForce eval_force,
                                     const int system_index);

} // namespace energy
} // namespace omni

#endif
