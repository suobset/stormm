// -*-c++-*-
#ifndef OMNI_NONBONDED_POTENTIAL_H
#define OMNI_NONBONDED_POTENTIAL_H

#include "Constants/generalized_born.h"
#include "Topology/atomgraph.h"
#include "Trajectory/coordinateframe.h"
#include "Trajectory/phasespace.h"
#include "scorecard.h"
#include "static_exclusionmask.h"
#include "energy_enumerators.h"

namespace omni {
namespace energy {

using topology::AtomGraph;
using topology::ImplicitSolventKit;
using topology::NonbondedKit;
using topology::UnitCellType;
using trajectory::CoordinateFrame;
using trajectory::CoordinateFrameReader;
using trajectory::CoordinateFrameWriter;
using trajectory::PhaseSpace;
using trajectory::PhaseSpaceWriter;
using namespace generalized_born_defaults;

/// \brief Evaluate the non-bonded electrostatic energy using an all-to-all approach with no
///        imaging considerations (isolated boundary conditions).
///
/// Overloaded:
///   - Operate on raw pointers for the coordinates, box transformations, and forces
///   - Pass the original topology by pointer or by reference
///   - Pass the non-bonded parameter kit by value
///   - Evaluate based on a PhaseSpace object, with the option to compute forces
///   - Evaluate energy only, based on a CoordinateFrame abstract
///
/// \param ag            System topology
/// \param nbk           Non-bonded parameters abstract taken from the original topology
/// \param se            Exclusion mask providing bits for all atom pairs
/// \param ps            Coordinates and force accumulators (modified by this function)
/// \param psw           Coordinates and force accumulators (modified by this function)
/// \param cfr           Coordinates of all particles, plus box dimensions (if needed)
/// \param cfw           Coordinates of all particles, plus box dimensions (if needed)
/// \param xcrd          Cartesian X coordinates of all particles
/// \param ycrd          Cartesian Y coordinates of all particles
/// \param zcrd          Cartesian Z coordinates of all particles
/// \param xfrc          Cartesian X forces acting on all particles
/// \param yfrc          Cartesian X forces acting on all particles
/// \param zfrc          Cartesian X forces acting on all particles
/// \param umat          Box space transformation matrix
/// \param invu          Inverse transformation matrix, fractional coordinates back to real space
/// \param unit_cell     The unit cell type, i.e. triclinic
/// \param ecard         Energy components and other state variables (volume, temperature, etc.)
///                      (modified by this function)
/// \param eval_force    Flag to have forces also evaluated
/// \param system_index  Index of the system to which this energy contributes
/// \{
double2 evaluateNonbondedEnergy(const NonbondedKit<double> nbk, const StaticExclusionMask &se,
                                const double* xcrd, const double* ycrd, const double* zcrd,
                                const double* umat, const double* invu,
                                UnitCellType unit_cell, double* xfrc, double* yfrc,
                                double* zfrc, ScoreCard *ecard,
                                EvaluateForce eval_elec_force = EvaluateForce::NO,
                                EvaluateForce eval_vdw_force = EvaluateForce::NO,
                                int system_index = 0);

double2 evaluateNonbondedEnergy(const NonbondedKit<double> nbk, const StaticExclusionMask &se,
                                PhaseSpaceWriter psw, ScoreCard *ecard,
                                EvaluateForce eval_elec_force = EvaluateForce::NO,
                                EvaluateForce eval_vdw_force  = EvaluateForce::NO,
                                int system_index = 0);

double2 evaluateNonbondedEnergy(const AtomGraph &ag, const StaticExclusionMask &se, PhaseSpace *ps,
                                ScoreCard *ecard,
                                EvaluateForce eval_elec_force = EvaluateForce::NO,
                                EvaluateForce eval_vdw_force  = EvaluateForce::NO,
                                int system_index = 0);

double2 evaluateNonbondedEnergy(const AtomGraph *ag, const StaticExclusionMask &se, PhaseSpace *ps,
                                ScoreCard *ecard,
                                EvaluateForce eval_elec_force = EvaluateForce::NO,
                                EvaluateForce eval_vdw_force  = EvaluateForce::NO,
                                int system_index = 0);

double2 evaluateNonbondedEnergy(const NonbondedKit<double> nbk, const StaticExclusionMask &se,
                                CoordinateFrameReader cfr, ScoreCard *ecard, int system_index = 0);

double2 evaluateNonbondedEnergy(const NonbondedKit<double> nbk, const StaticExclusionMask &se,
                                const CoordinateFrameWriter &cfw, ScoreCard *ecard,
                                int system_index = 0);

double2 evaluateNonbondedEnergy(const AtomGraph &ag, const StaticExclusionMask &se,
                                const CoordinateFrameReader &cfr, ScoreCard *ecard,
                                int system_index = 0);

double2 evaluateNonbondedEnergy(const AtomGraph *ag, const StaticExclusionMask &se,
                                const CoordinateFrameReader &cfr, ScoreCard *ecard,
                                int system_index = 0);
/// \}

/// \brief Evaluate the non-bonded Generalized Born energy of a system of particles (no cutoff is
///        applied in computation of the radii or evaluation of the energy and forces)
///
/// Overloaded:
///   - Operate on raw pointers for the coordinates, box transformations, and forces
///   - Pass the original topology by pointer or by reference,
///   - Pass the non-bonded parameter kit by value
///   - Evaluate based on a PhaseSpace object, with the option to compute forces
///   - Evaluate energy only, based on a CoordinateFrame abstract
///
/// \param ag            System topology
/// \param nbk           Non-bonded parameters abstract taken from the original topology
/// \param ngb_tables    "Neck" Generalized Born tables from some pre-loaded cache of constants
/// \param ps            Coordinates and force accumulators (modified by this function)
/// \param psw           Coordinates and force accumulators (modified by this function)
/// \param cfr           Coordinates of all particles, plus box dimensions (if needed)
/// \param cfw           Coordinates of all particles, plus box dimensions (if needed)
/// \param xcrd          Cartesian X coordinates of all particles
/// \param ycrd          Cartesian Y coordinates of all particles
/// \param zcrd          Cartesian Z coordinates of all particles
/// \param xfrc          Cartesian X forces acting on all particles
/// \param yfrc          Cartesian X forces acting on all particles
/// \param zfrc          Cartesian X forces acting on all particles
/// \param umat          Box space transformation matrix
/// \param invu          Inverse transformation matrix, fractional coordinates back to real space
/// \param unit_cell     The unit cell type, i.e. triclinic
/// \param ecard         Energy components and other state variables (volume, temperature, etc.)
///                      (modified by this function)
/// \param eval_force    Flag to have forces also evaluated
/// \param system_index  Index of the system to which this energy contributes
/// \{
double evaluateGeneralizedBornEnergy(const NonbondedKit<double> nbk,
                                     const ImplicitSolventKit<double> isk,
                                     const NeckGeneralizedBornTable &ngb_tables,
                                     const double* xcrd, const double* ycrd, const double* zcrd,
                                     const double* umat, const double* invu,
                                     UnitCellType unit_cell, double* xfrc, double* yfrc,
                                     double* zfrc, ScoreCard *ecard,
                                     EvaluateForce eval_force = EvaluateForce::NO,
                                     int system_index = 0);

double evaluateGeneralizedBornEnergy(const NonbondedKit<double> nbk,
                                     const ImplicitSolventKit<double> isk,
                                     const NeckGeneralizedBornTable &ngb_tables,
                                     PhaseSpaceWriter psw, ScoreCard *ecard,
                                     EvaluateForce eval_force = EvaluateForce::NO,
                                     int system_index = 0);

double evaluateGeneralizedBornEnergy(const AtomGraph &ag,
                                     const NeckGeneralizedBornTable &ngb_tables, PhaseSpace *ps,
                                     ScoreCard *ecard,
                                     EvaluateForce eval_force = EvaluateForce::NO,
                                     int system_index = 0);

double evaluateGeneralizedBornEnergy(const AtomGraph *ag,
                                     const NeckGeneralizedBornTable &ngb_tables, PhaseSpace *ps,
                                     ScoreCard *ecard,
                                     EvaluateForce eval_force = EvaluateForce::NO,
                                     int system_index = 0);

double evaluateGeneralizedBornEnergy(const NonbondedKit<double> nbk,
                                     const ImplicitSolventKit<double> isk,
                                     const NeckGeneralizedBornTable &ngb_tables,
                                     const CoordinateFrameReader cfr, ScoreCard *ecard,
                                     int system_index = 0);

double evaluateGeneralizedBornEnergy(const NonbondedKit<double> nbk,
                                     const ImplicitSolventKit<double> isk,
                                     const NeckGeneralizedBornTable &ngb_tables,
                                     const CoordinateFrameWriter &cfw, ScoreCard *ecard,
                                     int system_index = 0);

double evaluateGeneralizedBornEnergy(const AtomGraph &ag,
                                     const NeckGeneralizedBornTable &ngb_tables,
                                     const CoordinateFrameReader &cfr, ScoreCard *ecard,
                                     int system_index = 0);

double evaluateGeneralizedBornEnergy(const AtomGraph *ag,
                                     const NeckGeneralizedBornTable &ngb_tables,
                                     const CoordinateFrameReader &cfr, ScoreCard *ecard,
                                     int system_index = 0);
/// \}

} // namespace energy
} // namespace omni

#endif
