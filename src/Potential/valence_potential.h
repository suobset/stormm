// -*-c++-*-
#ifndef OMNI_VALENCE_POTENTIAL_H
#define OMNI_VALENCE_POTENTIAL_H

#include "Trajectory/coordinateframe.h"
#include "Trajectory/phasespace.h"
#include "Topology/atomgraph.h"
#include "Topology/atomgraph_abstracts.h"
#include "energy_enumerators.h"
#include "scorecard.h"

namespace omni {
namespace energy {

using topology::AtomGraph;
using topology::ValenceKit;
using topology::NonbondedKit;
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
/// \param ecard         Energy components and other state variables (volume, temperature, etc.)
///                      (modified by this function)
/// \param eval_force    Flag to have forces also evaluated
/// \param system_index  Index of the system to which this energy contributes
/// \{
double evaluateBondTerms(const ValenceKit<double> vk, PhaseSpaceWriter psw, ScoreCard *ecard,
                         EvaluateForce eval_force = EvaluateForce::NO, int system_index = 0);

double evaluateBondTerms(const AtomGraph &ag, PhaseSpace *ps, ScoreCard *ecard,
                         EvaluateForce eval_force = EvaluateForce::NO, int system_index = 0);

double evaluateBondTerms(const AtomGraph *ag, PhaseSpace *ps, ScoreCard *ecard,
                         EvaluateForce eval_force = EvaluateForce::NO, int system_index = 0);

double evaluateBondTerms(const ValenceKit<double> vk, const CoordinateFrameReader cfr,
                         ScoreCard *ecard, const int system_index);

double evaluateBondTerms(const ValenceKit<double> vk, const CoordinateFrameWriter &cfw,
                         ScoreCard *ecard, const int system_index);
/// \}
  
/// \brief Evaluate the angle bending energy contributions with a simple routine based on a
///        topology and a coordinate set.  These simple routines can serve as a check on much more
///        complex routines involving streamlined data structures and GPU execution.
///
/// Overloaded:
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
/// \param ecard         Energy components and other state variables (volume, temperature, etc.)
///                      (modified by this function)
/// \param eval_force    Flag to have forces also evaluated
/// \param system_index  Index of the system to which this energy contributes
/// \{
double evaluateAngleTerms(const ValenceKit<double> vk, PhaseSpaceWriter psw, ScoreCard *ecard,
                          EvaluateForce eval_force = EvaluateForce::NO, int system_index = 0);
                          
double evaluateAngleTerms(const AtomGraph &ag, PhaseSpace *ps, ScoreCard *ecard,
                          EvaluateForce eval_force = EvaluateForce::NO, int system_index = 0);

double evaluateAngleTerms(const AtomGraph *ag, PhaseSpace *ps, ScoreCard *ecard,
                          EvaluateForce eval_force = EvaluateForce::NO, int system_index = 0);

double evaluateAngleTerms(const ValenceKit<double> vk, const CoordinateFrameReader cfr,
                          ScoreCard *ecard, const int system_index);

double evaluateAngleTerms(const ValenceKit<double> vk, const CoordinateFrameWriter &cfw,
                          ScoreCard *ecard, const int system_index);
/// \}

/// \brief Evaluate the proper and improper dihedral energy contributions with a simple routine
///        based on a topology and a PhaseSpace object to store forces in double precision.  The
///        contributions of proper and improper dihedrals are stored in the x and y components
///        of the resulting double-precision tuple, respectively.  These results can be compared
///        to the fixed-precision accumulators in the energy tracking object.
///
/// Overloaded:
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
/// \param ecard         Energy components and other state variables (volume, temperature, etc.)
///                      (modified by this function)
/// \param eval_force    Flag to have forces also evaluated
/// \param system_index  Index of the system to which this energy contributes
/// \{
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
/// \}
  
/// \brief Evaluate Urey-Bradley harmonic angle interactions with a simple routine.  This looks
///        almost exactly like the bond computations but it kept in this separate routine for
///        simplicity of bookkeeping.  The resulting energy accumulated in double-precision can be
///        compared to the result obtained with fixed-precision accumulation.
///
/// Overloaded:
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
/// \param ecard         Energy components and other state variables (volume, temperature, etc.)
///                      (modified by this function)
/// \param eval_force    Flag to have forces also evaluated
/// \param system_index  Index of the system to which this energy contributes
/// \{
double evaluateUreyBradleyTerms(const ValenceKit<double> vk, PhaseSpaceWriter psw,
                                ScoreCard *ecard, const EvaluateForce eval_force,
                                const int system_index);

double evaluateUreyBradleyTerms(const AtomGraph &ag, PhaseSpace *ps, ScoreCard *ecard,
                                const EvaluateForce eval_force, const int system_index);

double evaluateUreyBradleyTerms(const AtomGraph *ag, PhaseSpace *ps, ScoreCard *ecard,
                                const EvaluateForce eval_force, const int system_index);

double evaluateUreyBradleyTerms(const ValenceKit<double> vk, const CoordinateFrameReader cfr,
                                ScoreCard *ecard, const int system_index);

double evaluateUreyBradleyTerms(const ValenceKit<double> vk, const CoordinateFrameWriter &cfw,
                                ScoreCard *ecard, const int system_index);
/// \}
  
/// \brief Evaluate CHARMM harmonic improper dihedral terms with a simple routine.  This
///        contributes to a separate improper energy accumulator, CHARMM_IMPROPER, in the
///        ScoreCard object.  The double-precision output of this function can be compared to the
///        result in fixed-precision accumulation.
///
/// Overloaded:
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
/// \param ecard         Energy components and other state variables (volume, temperature, etc.)
///                      (modified by this function)
/// \param eval_force    Flag to have forces also evaluated
/// \param system_index  Index of the system to which this energy contributes
/// \{
double evaluateCharmmImproperTerms(const ValenceKit<double> vk, PhaseSpaceWriter psw,
                                   ScoreCard *ecard, const EvaluateForce eval_force,
                                   const int system_index);

double evaluateCharmmImproperTerms(const AtomGraph &ag, PhaseSpace *ps, ScoreCard *ecard,
                                   const EvaluateForce eval_force, const int system_index);

double evaluateCharmmImproperTerms(const AtomGraph *ag, PhaseSpace *ps, ScoreCard *ecard,
                                   const EvaluateForce eval_force, const int system_index);

double evaluateCharmmImproperTerms(const ValenceKit<double> vk, const CoordinateFrameReader cfr,
                                   ScoreCard *ecard, const int system_index);

double evaluateCharmmImproperTerms(const ValenceKit<double> vk, const CoordinateFrameWriter &cfw,
                                   ScoreCard *ecard, const int system_index);
/// \}

/// \brief Evaluate CHARMM CMAP two-dimensional cubic spline potentials.  As with all other
///        functions in this library, the results of double-precision accumulation are returned
///        for comparison to fixed-precision results.
///
/// Overloaded:
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
/// \param ecard         Energy components and other state variables (volume, temperature, etc.)
///                      (modified by this function)
/// \param eval_force    Flag to have forces also evaluated
/// \param system_index  Index of the system to which this energy contributes
/// \{
double evaluateCmapTerms(const ValenceKit<double> vk, PhaseSpaceWriter psw, ScoreCard *ecard,
                         const EvaluateForce eval_force, const int system_index);

double evaluateCmapTerms(const AtomGraph &ag, PhaseSpace *ps, ScoreCard *ecard,
                         const EvaluateForce eval_force, const int system_index);

double evaluateCmapTerms(const AtomGraph *ag, PhaseSpace *ps, ScoreCard *ecard,
                         const EvaluateForce eval_force, const int system_index);

double evaluateCmapTerms(const ValenceKit<double> vk, const CoordinateFrameReader cfr,
                         ScoreCard *ecard, const int system_index);

double evaluateCmapTerms(const ValenceKit<double> vk, const CoordinateFrameWriter &cfw,
                         ScoreCard *ecard, const int system_index);
/// \}

/// \brief Evaluate 1:4 non-bonded pair interactions.  This requires a suprising amount of
///        bookkeeping to make it performant, but the result is straightforward and this reference
///        routine will work from that setup.
///
/// Overloaded:
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
/// \param ecard            Energy components and other state variables (volume, temperature, etc.)
///                         (modified by this function)
/// \param eval_elec_force  Flag to have electrostatic forces evaluated
/// \param eval_vdw_force   Flag to have van-der Waals (Lennard-Jones) forces evaluated
/// \param system_index     Index of the system to which this energy contributes
/// \{
double2 evaluateAttenuated14Terms(const ValenceKit<double> vk, const NonbondedKit<double> nbk,
                                  PhaseSpaceWriter psw, ScoreCard *ecard,
                                  const EvaluateForce eval_elec_force,
                                  const EvaluateForce eval_vdw_force, const int system_index);

double2 evaluateAttenuated14Terms(const AtomGraph &ag, PhaseSpace *ps, ScoreCard *ecard,
                                  const EvaluateForce eval_elec_force,
                                  const EvaluateForce eval_vdw_force, const int system_index);

double2 evaluateAttenuated14Terms(const AtomGraph *ag, PhaseSpace *ps, ScoreCard *ecard,
                                  const EvaluateForce eval_elec_force,
                                  const EvaluateForce eval_vdw_force, const int system_index);

double2 evaluateAttenuated14Terms(const ValenceKit<double> vk, const NonbondedKit<double> nbk,
                                  const CoordinateFrameReader cfr,
                                  ScoreCard *ecard, const int system_index);

double2 evaluateAttenuated14Terms(const ValenceKit<double> vk, const NonbondedKit<double> nbk,
                                  const CoordinateFrameWriter &cfw,
                                  ScoreCard *ecard, const int system_index);
/// \}

} // namespace energy
} // namespace omni

#endif
