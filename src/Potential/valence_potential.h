// -*-c++-*-
#ifndef OMNI_VALENCE_POTENTIAL_H
#define OMNI_VALENCE_POTENTIAL_H

#include "Restraints/restraint_apparatus.h"
#include "Topology/atomgraph.h"
#include "Topology/atomgraph_abstracts.h"
#include "Topology/atomgraph_enumerators.h"
#include "Trajectory/coordinateframe.h"
#include "Trajectory/phasespace.h"
#include "energy_enumerators.h"
#include "scorecard.h"

namespace omni {
namespace energy {

using restraints::RestraintApparatus;
using restraints::RestraintApparatusDpReader;
using topology::AtomGraph;
using topology::ValenceKit;
using topology::NonbondedKit;
using topology::UnitCellType;
using trajectory::CoordinateFrameReader;
using trajectory::CoordinateFrameWriter;
using trajectory::PhaseSpace;
using trajectory::PhaseSpaceWriter;

/// \brief A value which approaches one from below.  This is used in dihedral and similar
///        computations to detect when a value is nearing 1.0 and might generate a singularity
///        in some denominator.  It is nearly the closest value to 1.0 that can be represented in
///        an IEEE-754 format 32-bit floating point number (the format can go a couple of bits'
///        worth of precision closer, but this is 1 part in about a million).
/// \{
constexpr double asymptotic_to_one_lf = 0.99999904632568359375;
constexpr float  asymptotic_to_one_f  = (float)asymptotic_to_one_lf;
/// \}

/// \brief A value which captures 1 / (1 - asymptotic_to_one), to put a cap on the value of such
///        fractions 1 / (1 - x) as x -> 1.
/// \{
constexpr double inverse_one_minus_asymptote_lf = 1048576.0;
constexpr float  inverse_one_minus_asymptote_f = (float)1048576.0;
/// \}

/// \brief Evaluate the bond energy contributions based on a topology and a coordinate set.
///        These simple routines can serve as a check on much more complex routines involving
///        streamlined data structures and GPU execution.
///
/// Overloaded:
///   - Evaluate based on raw pointers to coordinates, box transformations, and forces
///   - Evaluate based on a PhaseSpace object, with the option to compute and store forces
///   - Evaluate energy only based on a CoordinateFrame abstract
///   - Pass a topology by pointer, by reference, or just the ValenceKit abstract by value
///
/// \param ag            System topology
/// \param vk            Valence parameters abstract from the system topology
/// \param ps            Coordinates, box size, and force accumulators (modified by this function)
/// \param psw           Coordinates, box size, and force accumulators (modified by this function)
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
double evaluateBondTerms(const ValenceKit<double> vk, const double* xcrd, const double* ycrd,
                         const double* zcrd, const double* umat, const double* invu,
                         UnitCellType unit_cell, double* xfrc, double* yfrc, double* zfrc,
                         ScoreCard *ecard, EvaluateForce eval_force = EvaluateForce::NO,
                         int system_index = 0);
  
double evaluateBondTerms(const ValenceKit<double> vk, PhaseSpaceWriter psw, ScoreCard *ecard,
                         EvaluateForce eval_force = EvaluateForce::NO, int system_index = 0);

double evaluateBondTerms(const AtomGraph &ag, PhaseSpace *ps, ScoreCard *ecard,
                         EvaluateForce eval_force = EvaluateForce::NO, int system_index = 0);

double evaluateBondTerms(const AtomGraph *ag, PhaseSpace *ps, ScoreCard *ecard,
                         EvaluateForce eval_force = EvaluateForce::NO, int system_index = 0);

double evaluateBondTerms(const ValenceKit<double> vk, const CoordinateFrameReader cfr,
                         ScoreCard *ecard, int system_index = 0);

double evaluateBondTerms(const ValenceKit<double> vk, const CoordinateFrameWriter &cfw,
                         ScoreCard *ecard, int system_index = 0);

double evaluateBondTerms(const AtomGraph &ag, const CoordinateFrameReader &cfr,
                         ScoreCard *ecard, int system_index = 0);

double evaluateBondTerms(const AtomGraph *ag, const CoordinateFrameReader &cfr,
                         ScoreCard *ecard, int system_index = 0);
/// \}
  
/// \brief Evaluate the angle bending energy contributions with a simple routine based on a
///        topology and a coordinate set.  These simple routines can serve as a check on much more
///        complex routines involving streamlined data structures and GPU execution.
///
/// Overloaded:
///   - Evaluate based on raw pointers to coordinates, box transformations, and forces
///   - Evaluate based on a PhaseSpace object, with the option to compute and store forces
///   - Evaluate energy only based on a CoordinateFrame abstract
///   - Pass a topology by pointer, by reference, or just the ValenceKit abstract by value
///
/// \param ag            System topology
/// \param vk            Valence parameters abstract from the system topology
/// \param ps            Coordinates, box size, and force accumulators (modified by this function)
/// \param psw           Coordinates, box size, and force accumulators (modified by this function)
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
double evaluateAngleTerms(const ValenceKit<double> vk, const double* xcrd, const double* ycrd,
                          const double* zcrd, const double* umat, const double* invu,
                          UnitCellType unit_cell, double* xfrc, double* yfrc, double* zfrc,
                          ScoreCard *ecard, EvaluateForce eval_force = EvaluateForce::NO,
                          int system_index = 0);

double evaluateAngleTerms(const ValenceKit<double> vk, PhaseSpaceWriter psw, ScoreCard *ecard,
                          EvaluateForce eval_force = EvaluateForce::NO, int system_index = 0);
                          
double evaluateAngleTerms(const AtomGraph &ag, PhaseSpace *ps, ScoreCard *ecard,
                          EvaluateForce eval_force = EvaluateForce::NO, int system_index = 0);

double evaluateAngleTerms(const AtomGraph *ag, PhaseSpace *ps, ScoreCard *ecard,
                          EvaluateForce eval_force = EvaluateForce::NO, int system_index = 0);

double evaluateAngleTerms(const ValenceKit<double> vk, const CoordinateFrameReader cfr,
                          ScoreCard *ecard, int system_index = 0);

double evaluateAngleTerms(const ValenceKit<double> vk, const CoordinateFrameWriter &cfw,
                          ScoreCard *ecard, int system_index = 0);

double evaluateAngleTerms(const AtomGraph &ag, const CoordinateFrameReader &cfr,
                          ScoreCard *ecard, int system_index = 0);

double evaluateAngleTerms(const AtomGraph *ag, const CoordinateFrameReader &cfr,
                          ScoreCard *ecard, int system_index = 0);
/// \}

/// \brief Evaluate the proper and improper dihedral energy contributions with a simple routine
///        based on a topology and a PhaseSpace object to store forces in double precision.  The
///        contributions of proper and improper dihedrals are stored in the x and y components
///        of the resulting double-precision tuple, respectively.  These results can be compared
///        to the fixed-precision accumulators in the energy tracking object.
///
/// Overloaded:
///   - Evaluate based on raw pointers to coordinates, box transformations, and forces
///   - Evaluate based on a PhaseSpace object, with the option to compute and store forces
///   - Evaluate energy only based on a CoordinateFrame abstract
///   - Pass a topology by pointer, by reference, or just the ValenceKit abstract by value
///
/// \param ag            System topology
/// \param vk            Valence parameters abstract from the system topology
/// \param ps            Coordinates, box size, and force accumulators (modified by this function)
/// \param psw           Coordinates, box size, and force accumulators (modified by this function)
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
double2 evaluateDihedralTerms(const ValenceKit<double> vk, const double* xcrd, const double* ycrd,
                              const double* zcrd, const double* umat, const double* invu,
                              UnitCellType unit_cell, double* xfrc, double* yfrc, double* zfrc,
                              ScoreCard *ecard, EvaluateForce eval_force = EvaluateForce::NO,
                              int system_index = 0);

double2 evaluateDihedralTerms(const ValenceKit<double> vk, PhaseSpaceWriter psw, ScoreCard *ecard,
                              EvaluateForce eval_force = EvaluateForce::NO, int system_index = 0);

double2 evaluateDihedralTerms(const AtomGraph &ag, PhaseSpace *ps, ScoreCard *ecard,
                              EvaluateForce eval_force = EvaluateForce::NO, int system_index = 0);

double2 evaluateDihedralTerms(const AtomGraph *ag, PhaseSpace *ps, ScoreCard *ecard,
                              EvaluateForce eval_force = EvaluateForce::NO, int system_index = 0);

double2 evaluateDihedralTerms(const ValenceKit<double> vk, const CoordinateFrameReader cfr,
                              ScoreCard *ecard, int system_index = 0);

double2 evaluateDihedralTerms(const ValenceKit<double> vk, const CoordinateFrameWriter &cfw,
                              ScoreCard *ecard, int system_index = 0);

double2 evaluateDihedralTerms(const AtomGraph &ag, const CoordinateFrameReader &cfr,
                              ScoreCard *ecard, int system_index = 0);

double2 evaluateDihedralTerms(const AtomGraph *ag, const CoordinateFrameReader &cfr,
                              ScoreCard *ecard, int system_index = 0);
/// \}
  
/// \brief Evaluate Urey-Bradley harmonic angle interactions with a simple routine.  This looks
///        almost exactly like the bond computations but it kept in this separate routine for
///        simplicity of bookkeeping.  The resulting energy accumulated in double-precision can be
///        compared to the result obtained with fixed-precision accumulation.
///
/// Overloaded:
///   - Evaluate based on raw pointers to coordinates, box transformations, and forces
///   - Evaluate based on a PhaseSpace object, with the option to compute and store forces
///   - Evaluate energy only based on a CoordinateFrame abstract
///   - Pass a topology by pointer, by reference, or just the ValenceKit abstract by value
///
/// \param ag            System topology
/// \param vk            Valence parameters abstract from the system topology
/// \param ps            Coordinates, box size, and force accumulators (modified by this function)
/// \param psw           Coordinates, box size, and force accumulators (modified by this function)
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
double evaluateUreyBradleyTerms(const ValenceKit<double> vk, const double* xcrd,
                                const double* ycrd, const double* zcrd, const double* umat,
                                const double* invu, const UnitCellType unit_cell, double* xfrc,
                                double* yfrc, double* zfrc, ScoreCard *ecard,
                                EvaluateForce eval_force = EvaluateForce::NO,
                                int system_index = 0);

double evaluateUreyBradleyTerms(const ValenceKit<double> vk, PhaseSpaceWriter psw,
                                ScoreCard *ecard, EvaluateForce eval_force = EvaluateForce::NO,
                                int system_index = 0);

double evaluateUreyBradleyTerms(const AtomGraph &ag, PhaseSpace *ps, ScoreCard *ecard,
                                EvaluateForce eval_force = EvaluateForce::NO,
                                int system_index = 0);

double evaluateUreyBradleyTerms(const AtomGraph *ag, PhaseSpace *ps, ScoreCard *ecard,
                                EvaluateForce eval_force = EvaluateForce::NO,
                                int system_index = 0);

double evaluateUreyBradleyTerms(const ValenceKit<double> vk, const CoordinateFrameReader cfr,
                                ScoreCard *ecard, int system_index = 0);

double evaluateUreyBradleyTerms(const ValenceKit<double> vk, const CoordinateFrameWriter &cfw,
                                ScoreCard *ecard, int system_index = 0);

double evaluateUreyBradleyTerms(const AtomGraph &ag, const CoordinateFrameReader &cfr,
                                ScoreCard *ecard, int system_index = 0);

double evaluateUreyBradleyTerms(const AtomGraph *ag, const CoordinateFrameReader &cfr,
                                ScoreCard *ecard, int system_index = 0);
/// \}
  
/// \brief Evaluate CHARMM harmonic improper dihedral terms with a simple routine.  This
///        contributes to a separate improper energy accumulator, CHARMM_IMPROPER, in the
///        ScoreCard object.  The double-precision output of this function can be compared to the
///        result in fixed-precision accumulation.
///
/// Overloaded:
///   - Evaluate based on raw pointers to coordinates, box transformations, and forces
///   - Evaluate based on a PhaseSpace object, with the option to compute and store forces
///   - Evaluate energy only based on a CoordinateFrame abstract
///   - Pass a topology by pointer, by reference, or just the ValenceKit abstract by value
///
/// \param ag            System topology
/// \param vk            Valence parameters abstract from the system topology
/// \param ps            Coordinates, box size, and force accumulators (modified by this function)
/// \param psw           Coordinates, box size, and force accumulators (modified by this function)
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
double evaluateCharmmImproperTerms(const ValenceKit<double> vk, const double* xcrd,
                                   const double* ycrd, const double* zcrd, const double* umat,
                                   const double* invu, UnitCellType unit_cell, double* xfrc,
                                   double* yfrc, double* zfrc, ScoreCard *ecard,
                                   EvaluateForce eval_force = EvaluateForce::NO,
                                   int system_index = 0);

double evaluateCharmmImproperTerms(const ValenceKit<double> vk, PhaseSpaceWriter psw,
                                   ScoreCard *ecard, EvaluateForce eval_force = EvaluateForce::NO,
                                   int system_index = 0);

double evaluateCharmmImproperTerms(const AtomGraph &ag, PhaseSpace *ps, ScoreCard *ecard,
                                   EvaluateForce eval_force = EvaluateForce::NO,
                                   int system_index = 0);

double evaluateCharmmImproperTerms(const AtomGraph *ag, PhaseSpace *ps, ScoreCard *ecard,
                                   EvaluateForce eval_force = EvaluateForce::NO,
                                   int system_index = 0);

double evaluateCharmmImproperTerms(const ValenceKit<double> vk, const CoordinateFrameReader cfr,
                                   ScoreCard *ecard, int system_index = 0);

double evaluateCharmmImproperTerms(const ValenceKit<double> vk, const CoordinateFrameWriter &cfw,
                                   ScoreCard *ecard, int system_index = 0);

double evaluateCharmmImproperTerms(const AtomGraph &ag, const CoordinateFrameReader &cfr,
                                   ScoreCard *ecard, int system_index = 0);

double evaluateCharmmImproperTerms(const AtomGraph *ag, const CoordinateFrameReader &cfr,
                                   ScoreCard *ecard, int system_index = 0);
/// \}

/// \brief Evaluate CHARMM CMAP two-dimensional cubic spline potentials.  As with all other
///        functions in this library, the results of double-precision accumulation are returned
///        for comparison to fixed-precision results.
///
/// Overloaded:
///   - Evaluate based on raw pointers to coordinates, box transformations, and forces
///   - Evaluate based on a PhaseSpace object, with the option to compute and store forces
///   - Evaluate energy only based on a CoordinateFrame abstract
///   - Pass a topology by pointer, by reference, or just the ValenceKit abstract by value
///
/// \param ag            System topology
/// \param vk            Valence parameters abstract from the system topology
/// \param ps            Coordinates, box size, and force accumulators (modified by this function)
/// \param psw           Coordinates, box size, and force accumulators (modified by this function)
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
double evaluateCmapTerms(const ValenceKit<double> vk, const double* xcrd, const double* ycrd,
                         const double* zcrd, const double* umat, const double* invu,
                         UnitCellType unit_cell, double* xfrc, double* yfrc, double* zfrc,
                         ScoreCard *ecard, EvaluateForce eval_force = EvaluateForce::NO,
                         int system_index = 0);

double evaluateCmapTerms(const ValenceKit<double> vk, PhaseSpaceWriter psw, ScoreCard *ecard,
                         EvaluateForce eval_force = EvaluateForce::NO, int system_index = 0);

double evaluateCmapTerms(const AtomGraph &ag, PhaseSpace *ps, ScoreCard *ecard,
                         EvaluateForce eval_force = EvaluateForce::NO, int system_index = 0);

double evaluateCmapTerms(const AtomGraph *ag, PhaseSpace *ps, ScoreCard *ecard,
                         EvaluateForce eval_force = EvaluateForce::NO, int system_index = 0);

double evaluateCmapTerms(const ValenceKit<double> vk, const CoordinateFrameReader cfr,
                         ScoreCard *ecard, int system_index = 0);

double evaluateCmapTerms(const ValenceKit<double> vk, const CoordinateFrameWriter &cfw,
                         ScoreCard *ecard, int system_index = 0);

double evaluateCmapTerms(const AtomGraph &ag, const CoordinateFrameReader &cfr, ScoreCard *ecard,
                         int system_index = 0);

double evaluateCmapTerms(const AtomGraph *ag, const CoordinateFrameReader &cfr, ScoreCard *ecard,
                         int system_index = 0);
/// \}

/// \brief Evaluate 1:4 non-bonded pair interactions.  This requires a suprising amount of
///        bookkeeping to make it performant, but the result is straightforward and this reference
///        routine will work from that setup.
///
/// Overloaded:
///   - Evaluate based on raw pointers to coordinates, box transformations, and forces
///   - Evaluate based on a PhaseSpace object, with the option to compute and store forces
///   - Evaluate energy only based on a CoordinateFrame abstract
///   - Pass a topology by pointer, by reference, or just the ValenceKit abstract by value
///
/// \param ag              System topology
/// \param vk              Valence parameters abstract from the system topology
/// \param nbk             Non-bonded parameters abstract from the system topology
/// \param ps              Coordinates, box size, and force accumulators (modified by this
///                        function)
/// \param psw             Coordinates, box size, and force accumulators (modified by this
///                        function)
/// \param cfr             Coordinates of all particles, plus box dimensions (if needed)
/// \param cfw             Coordinates of all particles, plus box dimensions (if needed)
/// \param xcrd             Cartesian X coordinates of all particles
/// \param ycrd             Cartesian Y coordinates of all particles
/// \param zcrd             Cartesian Z coordinates of all particles
/// \param xfrc             Cartesian X forces acting on all particles
/// \param yfrc             Cartesian X forces acting on all particles
/// \param zfrc             Cartesian X forces acting on all particles
/// \param umat             Box space transformation matrix
/// \param invu             Inverse transformation, fractional coordinates back to real space
/// \param unit_cell        The unit cell type, i.e. triclinic
/// \param ecard            Energy components and other state variables (volume, temperature, etc.)
///                         (modified by this function)
/// \param eval_elec_force  Flag to have electrostatic forces evaluated
/// \param eval_vdw_force   Flag to have van-der Waals (Lennard-Jones) forces evaluated
/// \param system_index     Index of the system to which this energy contributes
/// \{
double2 evaluateAttenuated14Terms(const ValenceKit<double> vk, const NonbondedKit<double> nbk,
                                  const double* xcrd, const double* ycrd, const double* zcrd,
                                  const double* umat, const double* invu,
                                  UnitCellType unit_cell, double* xfrc, double* yfrc,
                                  double* zfrc, ScoreCard *ecard,
                                  EvaluateForce eval_elec_force = EvaluateForce::NO,
                                  EvaluateForce eval_vdw_force = EvaluateForce::NO,
                                  int system_index = 0);

double2 evaluateAttenuated14Terms(const ValenceKit<double> vk, const NonbondedKit<double> nbk,
                                  PhaseSpaceWriter psw, ScoreCard *ecard,
                                  EvaluateForce eval_elec_force = EvaluateForce::NO,
                                  EvaluateForce eval_vdw_force = EvaluateForce::NO,
                                  int system_index = 0);

double2 evaluateAttenuated14Terms(const AtomGraph &ag, PhaseSpace *ps, ScoreCard *ecard,
                                  EvaluateForce eval_elec_force = EvaluateForce::NO,
                                  EvaluateForce eval_vdw_force = EvaluateForce::NO,
                                  int system_index = 0);

double2 evaluateAttenuated14Terms(const AtomGraph *ag, PhaseSpace *ps, ScoreCard *ecard,
                                  EvaluateForce eval_elec_force = EvaluateForce::NO,
                                  EvaluateForce eval_vdw_force = EvaluateForce::NO,
                                  int system_index = 0);

double2 evaluateAttenuated14Terms(const ValenceKit<double> vk, const NonbondedKit<double> nbk,
                                  const CoordinateFrameReader cfr, ScoreCard *ecard,
                                  int system_index = 0);

double2 evaluateAttenuated14Terms(const ValenceKit<double> vk, const NonbondedKit<double> nbk,
                                  const CoordinateFrameWriter &cfw, ScoreCard *ecard,
                                  int system_index = 0);

double2 evaluateAttenuated14Terms(const AtomGraph &ag, const CoordinateFrameReader &cfr,
                                  ScoreCard *ecard, int system_index = 0);

double2 evaluateAttenuated14Terms(const AtomGraph *ag, const CoordinateFrameReader &cfr,
                                  ScoreCard *ecard, int system_index = 0);
/// \}

/// \brief Compute the mixture of end-point values that will determine the actual strength and
///        displacement settings of a flat-bottom bimodal harmonic restraint.  The flag about a
///        RestraintApparatus having time-dependent restraints is mostly for convenience, a way to
///        tell whether there is any time-dependent restraint in the collection at all.  The initial
///        and final settings of the steps for each restraint encode whether there is actual time
///        dependence in the result.
///
/// \param step_number  The current step number of the simulation (may include energy minimization
///                     step counts)
/// \param init_step    The initial step at which the restraint engages
/// \param final_step   The final step at which the restraint becomes mature
double2 computeRestraintMixture(int step_number, int init_step, int final_step);

/// \brief Compute critical elements of the restraining potential: its difference from the target
///        value that determines some harmonic stiffness penalty, the harmonic penalty stiffness,
///        and the energy contribution.
///
/// \param init_k   Initial stiffness parameters
/// \param final_k  Final stiffness parameters
/// \param init_r   Initial displacement parameters
/// \param final_r  Final displacement parameters
/// \param mixwt    Pre-calculated mixing factor for combining initial and final parameters
/// \param dr       The measured value of the restraint coordinate among its participating atoms
double3 restraintDelta(const double2 init_k, const double2 final_k, const double4 init_r,
                       const double4 final_r, const double2 mixwt, double dr);

/// \brief Evaluate flat-bottom bimodal harmonic restraints of the form used in Amber's sander and
///        pmemd programs for NMR annealing calculations.  Time dependence of the restraints is
///        recognized, although there is not the same diversity of time evolution functions in
///        OMNI.
///
/// Overloaded:
///   - Evaluate based on a PhaseSpace object, with the option to compute and store forces
///   - Evaluate energy only based on a CoordinateFrame abstract
///   - Pass the restraint collection by pointer, by reference, or the reader abstract by value
///
/// \param ra            Restraint apparatus applicable to the topology that describes the system
/// \param rar           Double-precision reader for the restraint apparatus 
/// \param ps            Coordinates, box size, and force accumulators (modified by this function)
/// \param psw           Coordinates, box size, and force accumulators (modified by this function)
/// \param cfr           Coordinates of all particles, plus box dimensions (if needed)
/// \param cfw           Coordinates of all particles, plus box dimensions (if needed)
/// \param ecard         Energy components and other state variables (volume, temperature, etc.)
///                      (modified by this function)
/// \param eval_force    Flag to have forces also evaluated
/// \param system_index  Index of the system to which this energy contributes
/// \param step_number   The step number at which the energy is being evaluated (may determine the
///                      restraint parameters by mixing their endpoint values)
/// \{
double evaluateRestraints(const RestraintApparatusDpReader rar, const double* xcrd,
                          const double* ycrd, const double* zcrd, const double* umat,
                          const double* invu, const UnitCellType unit_cell, double* xfrc,
                          double* yfrc, double* zfrc, ScoreCard *ecard,
                          EvaluateForce eval_force = EvaluateForce::NO, int system_index = 0,
                          int step_number = 0);

double evaluateRestraints(const RestraintApparatusDpReader rar, PhaseSpaceWriter psw,
                          ScoreCard *ecard, EvaluateForce eval_force = EvaluateForce::NO,
                          int system_index = 0, int step_number = 0);

double evaluateRestraints(const RestraintApparatus &ra, PhaseSpace *ps, ScoreCard *ecard,
                          EvaluateForce eval_force = EvaluateForce::NO, int system_index = 0,
                          int step_number = 0);

double evaluateRestraints(const RestraintApparatus *ra, PhaseSpace *ps, ScoreCard *ecard,
                          EvaluateForce eval_force = EvaluateForce::NO, int system_index = 0,
                          int step_number = 0);

double evaluateRestraints(const RestraintApparatusDpReader rar, const CoordinateFrameReader cfr,
                          ScoreCard *ecard, int system_index = 0, int step_number = 0);

double evaluateRestraints(const RestraintApparatusDpReader rar, const CoordinateFrameWriter &cfw,
                          ScoreCard *ecard, int system_index = 0, int step_number = 0);

double evaluateRestraints(const RestraintApparatus &ra, const CoordinateFrameReader cfr,
                          ScoreCard *ecard, int system_index = 0, int step_number = 0);

double evaluateRestraints(const RestraintApparatus *ra, const CoordinateFrameReader cfr,
                          ScoreCard *ecard, int system_index = 0, int step_number = 0);
/// \}
  
} // namespace energy
} // namespace omni

#endif
