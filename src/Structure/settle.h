// -*-c++-*-
#ifndef STORMM_SETTLE_H
#define STORMM_SETTLE_H

#include <typeinfo>
#include <typeindex>
#include <unistd.h>
#include "copyright.h"
#include "Constants/behavior.h"
#include "Constants/fixed_precision.h"
#include "DataTypes/common_types.h"
#include "DataTypes/stormm_vector_types.h"
#include "Math/vector_ops.h"
#include "Synthesis/atomgraph_synthesis.h"
#include "Synthesis/phasespace_synthesis.h"
#include "Synthesis/synthesis_abstracts.h"
#include "Synthesis/synthesis_enumerators.h"
#include "Topology/atomgraph.h"
#include "Topology/atomgraph_abstracts.h"
#include "Trajectory/phasespace.h"
#include "global_manipulation.h"

namespace stormm {
namespace structure {

using constants::PrecisionModel;
using data_types::isSignedIntegralScalarType;
using numerics::globalpos_scale_nonoverflow_bits;
using numerics::velocity_scale_nonoverflow_bits;
using stmath::crossProduct;
using stmath::normalize;
using synthesis::AtomGraphSynthesis;
using synthesis::maximum_valence_work_unit_atoms;
using synthesis::PhaseSpaceSynthesis;
using synthesis::SyValenceKit;
using synthesis::SyAtomUpdateKit;
using synthesis::vwu_abstract_length;
using synthesis::VwuAbstractMap;
using topology::AtomGraph;
using topology::ConstraintKit;
using trajectory::PhaseSpace;
using trajectory::PhaseSpaceWriter;
  
/// \brief Apply velocity restraints to a particular group of atoms.
///
/// \param oxy_idx      Index of the heavy (oxygen) atom in the arrays xcrd, ycrd, ..., zvel
/// \param hd1_idx      Index of the first light (hydrogen) atom in the arrays xcrd, ..., zvel
/// \param hd2_idx      Index of the second light (hydrogen) atom in the arrays xcrd, ..., zvel
/// \param m_o          Mass of the heavy atom
/// \param m_h          Mass of each of the light atoms
/// \param m_oh         Combined mass of the heavy atom and one light atom
/// \param xcrd         Cartesian X locations of the particles
/// \param ycrd         Cartesian Y locations of the particles
/// \param zcrd         Cartesian Z locations of the particles
/// \param xvel         Cartesian X velcoties of the particles (updated and returned)
/// \param yvel         Cartesian Y velcoties of the particles (updated and returned)
/// \param zvel         Cartesian Z velcoties of the particles (updated and returned)
/// \param xcrd_ovrf    Overflow bits for Cartesian X locations of the particles
/// \param ycrd_ovrf    Overflow bits for Cartesian Y locations of the particles
/// \param zcrd_ovrf    Overflow bits for Cartesian Z locations of the particles
/// \param xvel_ovrf    Overflow bits for Cartesian X velocities of the particles
/// \param yvel_ovrf    Overflow bits for Cartesian Y velocities of the particles
/// \param zvel_ovrf    Overflow bits for Cartesian Z velocities of the particles
/// \param gpos_scale   Scaling factor for fixed-precision particle positions
/// \param vel_scale    Scaling factor for fixed-precision particle velocties
template <typename Tcoord, typename Tcalc>
void settleGroupVelocity(int oxy_idx, int hd1_idx, int hd2_idx, Tcalc m_o, Tcalc m_h, Tcalc m_oh,
                         const Tcoord* xcrd, const Tcoord* ycrd, const Tcoord* zcrd, Tcoord* xvel,
                         Tcoord* yvel, Tcoord* zvel, const int* xcrd_ovrf = nullptr,
                         const int* ycrd_ovrf = nullptr, const int* zcrd_ovrf = nullptr,
                         int* xvel_ovrf = nullptr, int* yvel_ovrf = nullptr,
                         int* zvel_ovrf = nullptr, Tcalc gpos_scale = 1.0, Tcalc vel_scale = 1.0);

/// \brief Apply velocity constraints for rigid waters and perhaps other groups conforming to the
///        rigid symmetry required by the SETTLE approach.  Descriptions of input parameters follow
///        from settleGroupVelocity(), above, in addition to:
///
/// \param psw       Contains current coordinates and developing velocities of the system
/// \param poly_ps   Contains current coordinates and developing velocities of multiple systems
/// \param ag        Topology for the system of interest, from which constraint parameters can be
///                  extracted
/// \param poly_ag   Topology synthesis for multiple systems, from which constraint parameters can
///                  be extracted
/// \param cnst      Constraint parameters, including SETTLE group geometry and mass information
/// \param poly_vk   Valence abstract for the topology synthesis, contains valence work units with
///                  sets of imported atoms and constraint instruction codes
/// \param poly_auk  Constraint parameters, including SETTLE group geometry and mass information
///                  for a synthesis of systems
/// \{
template <typename Tcoord, typename Tcalc>
void settleVelocities(const Tcoord* xcrd, const Tcoord* ycrd, const Tcoord* zcrd, Tcoord* xvel,
                      Tcoord* yvel, Tcoord* zvel, const ConstraintKit<Tcalc> &cnst,
                      const int* xcrd_ovrf = nullptr, const int* ycrd_ovrf = nullptr,
                      const int* zcrd_ovrf = nullptr, int* xvel_ovrf = nullptr,
                      int* yvel_ovrf = nullptr, int* zvel_ovrf = nullptr, Tcalc gpos_scale = 1.0,
                      Tcalc vel_scale = 1.0);

template <typename Tcalc>
void settleVelocities(PhaseSpaceWriter *psw, const ConstraintKit<Tcalc> &cnst);

template <typename Tcalc, typename Tcalc2, typename Tcalc4>
void settleVelocities(PsSynthesisWriter *poly_psw, const SyValenceKit<Tcalc> &poly_vk,
                      const SyAtomUpdateKit<Tcalc, Tcalc2, Tcalc4> &poly_auk);

void settleVelocities(PhaseSpace *ps, const AtomGraph *ag, PrecisionModel prec);

void settleVelocities(PhaseSpace *ps, const AtomGraph &ag, PrecisionModel prec);

void settleVelocities(PhaseSpaceSynthesis *poly_ps, const AtomGraphSynthesis *poly_ag,
                      PrecisionModel prec = PrecisionModel::SINGLE);  

void settleVelocities(PhaseSpaceSynthesis *poly_ps, const AtomGraphSynthesis &poly_ag,
                      PrecisionModel prec = PrecisionModel::SINGLE);  
/// \}

/// \brief Utility function to encapsulate particle position and velocity updates at the end of
///        the SETTLE position constraint application, for the case of int95_t coordinates and one
///        of three Cartesian dimensions.
///
/// \param atm_idx        The index of the atom in question
/// \param p_to_v_factor  Scaling factor incorporating the inverse time step, inverse global
///                       position scaling factor, and the velocity scaling factor (this scales
///                       positions down to real-valued numbers in units of Angstroms, then scales
///                       that result up to a fixed-precision representation of the velocity)
/// \param r_adj          The displacement of the atom from the origin of the local frame of
///                       reference in which constraint work was done.  This displacement is
///                       expected to be provided such that the rotation back to the original
///                       orientation of the SETTLE group (e.g. water) has already been applied
///                       and with the appropriate fixed-precision position scaling folded in.
/// \param r_ref          Reference location of the original particle relative to the coordinate
///                       origin in which constraints were applied.  This value is expected to
///                       have the fixed-precision position scaling folded in.
/// \param r_ref_ovrf     Overflow bits for r_ref
/// \param r_dev          Array of developing particle positions, modified and returned
/// \param r_dev_ovrf     Overflow bits for positions in r_dev
/// \param v_dev          Array of developing particle velocities, modified and returned
/// \param v_dev_ovrf     Overflow bits for velocities in v_dev
template <typename Tcoord, typename Tcalc>
void updateSettlePVC(int atm_idx, Tcalc p_to_v_factor, Tcalc r_adj, Tcoord r_ref, int r_ref_ovrf,
                     Tcoord* r_dev, int* r_dev_ovrf, Tcoord* v_dev, int* v_dev_ovrf);
  
/// \brief Correct the geometry of a group of three atoms using the SETTLE analytic technique.  The
///        atoms must be such that there is one heavy atom and two light atoms each connected to
///        the heavy atom by a consistent angle and the same bond length.  Descriptions of input
///        parameters follow from settleGroupVelocity(), above, in addition to:
///
/// \param frac_heavy  Fraction of mass, out of the entire system, concentrated in the heavy atom
/// \param frac_light  Fraction of mass, out of the entire system, concentrated in the light atom
/// \param ra          Internal distance measurement relating to the displacement of the heavy
///                    atom from the center of mass
/// \param rb          Internal distance measurement relating to the Cartesian X displacement of
///                    each light atom from the center of mass in the primed coordinate system
/// \param rc          Internal distance measurement relating to the Cartesian Y displacement of
///                    each light atom from the center of mass in the primed coordinate system
/// \param xalt        Developing Cartesian X coordinates of particles in the SETTLE group
/// \param yalt        Developing Cartesian Y coordinates of particles in the SETTLE group
/// \param zalt        Developing Cartesian Z coordinates of particles in the SETTLE group
/// \param vxalt       Developing Cartesian X velocities of particles in the SETTLE group
/// \param vyalt       Developing Cartesian Y velocities of particles in the SETTLE group
/// \param vzalt       Developing Cartesian Z velocities of particles in the SETTLE group
/// \param dt          The time step (for making the velocity adjustment)
template <typename Tcoord, typename Tcalc, typename Tcalc3>
void settleGroupPosition(int oxy_idx, int hd1_idx, int hd2_idx, Tcalc frac_heavy, Tcalc frac_light,
                         Tcalc ra, Tcalc rb, Tcalc rc, const Tcoord* xcrd, const Tcoord* ycrd,
                         const Tcoord* zcrd, Tcoord* xalt, Tcoord* yalt, Tcoord* zalt,
                         Tcoord* vxalt,	Tcoord*	vyalt, Tcoord* vzalt, Tcalc dt,
                         const int* xcrd_ovrf = nullptr, const int* ycrd_ovrf = nullptr,
                         const int* zcrd_ovrf = nullptr, int* xalt_ovrf = nullptr,
                         int* yalt_ovrf = nullptr, int* zalt_ovrf = nullptr,
                         int* vxalt_ovrf = nullptr, int* vyalt_ovrf = nullptr,
                         int* vzalt_ovrf = nullptr, Tcalc gpos_scale = 1.0, Tcalc vel_scale = 1.0);
  
/// \brief Apply positional constraints for rigid waters and perhaps other groups conforming to the
///        rigid symmetry required by the SETTLE approach.  Descriptions of input parameters follow
///        from settleGroupPosition() and settleGroupVelocities(), above.
/// \{
template <typename Tcoord, typename Tcalc, typename Tcalc3>
void settlePositions(const Tcoord* xcrd, const Tcoord* ycrd, const Tcoord* zcrd, Tcoord* xalt,
                     Tcoord* yalt, Tcoord* zalt, Tcoord* vxalt, Tcoord* vyalt, Tcoord* vzalt,
                     const ConstraintKit<Tcalc> &cnst, Tcalc dt, const int* xcrd_ovrf = nullptr,
                     const int* ycrd_ovrf = nullptr, const int* zcrd_ovrf = nullptr,
                     int* xalt_ovrf = nullptr, int* yalt_ovrf = nullptr, int* zalt_ovrf = nullptr,
                     int* vxalt_ovrf = nullptr, int* vyalt_ovrf = nullptr,
                     int* vzalt_ovrf = nullptr, Tcalc gpos_scale = 1.0, Tcalc vel_scale = 1.0);

template <typename Tcalc, typename Tcalc3>
void settlePositions(PhaseSpaceWriter *psw, const ConstraintKit<Tcalc> &cnst, Tcalc dt);

void settlePositions(PhaseSpace *ps, const AtomGraph *ag, double dt,
                     PrecisionModel prec = PrecisionModel::DOUBLE);

void settlePositions(PhaseSpace *ps, const AtomGraph &ag, double dt,
                     PrecisionModel prec = PrecisionModel::DOUBLE);

template <typename Tcalc, typename Tcalc2, typename Tcalc3, typename Tcalc4>
void settlePositions(PsSynthesisWriter *poly_psw, const SyValenceKit<Tcalc> &poly_vk,
                     const SyAtomUpdateKit<Tcalc, Tcalc2, Tcalc4> &poly_auk, Tcalc dt);

void settlePositions(PhaseSpaceSynthesis *poly_ps, const AtomGraphSynthesis *poly_ag, double dt,
                     const PrecisionModel prec = PrecisionModel::SINGLE);

void settlePositions(PhaseSpaceSynthesis *poly_ps, const AtomGraphSynthesis &poly_ag, double dt,
                     const PrecisionModel prec = PrecisionModel::SINGLE);
/// \}

/// \brief Set positions for a system in the ideal geometry based on the plane of atoms in any
///        SETTLE groups.  This is not an energy-conserving move, rather a utility for generating
///        a plausible reference state.
///
/// Overloaded:
///   - Provide the original coordinates and topology (the function will operate on positions for
///     the current stage of the object's time cycle)
///   - Provide the abstracts of the coordinate object and constraint parameters
///
/// \param ps    Contains positions of all atoms in the system
/// \param psw   Contains positions of all atoms in the system.  This should be taken at the
///              current stage in the original object's time cycle.
/// \param ag    The original topology
/// \param cnst  Constraint parameters for SETTLE groups
/// \param prec  The precision model in which to evaluate the positions
/// \{
template <typename T>
void idealizeSettleReference(PhaseSpaceWriter *psw, const ConstraintKit<T> &cnst);

void idealizeSettleReference(PhaseSpace *ps, const AtomGraph *ag,
                             PrecisionModel prec = PrecisionModel::DOUBLE);

void idealizeSettleReference(PhaseSpace *ps, const AtomGraph &ag,
                             PrecisionModel prec = PrecisionModel::DOUBLE);
/// \}

} // namespace structure
} // namespace stormm

#include "settle.tpp"

#endif
