// -*-c++-*-
#ifndef OMNI_EVAL_SYNTHESIS_H
#define OMNI_EVAL_SYNTHESIS_H

#include "DataTypes/common_types.h"
#include "Math/vector_ops.h"
#include "Potential/scorecard.h"
#include "Potential/valence_potential.h"
#include "Potential/eval_valence_workunit.h"
#include "Synthesis/phasespace_synthesis.h"
#include "Synthesis/synthesis_abstracts.h"
#include "Synthesis/synthesis_enumerators.h"
#include "Synthesis/valence_workunit.h"

namespace omni {
namespace energy {

using energy::ScoreCard;
using energy::evalHarmonicStretch;
using energy::commitVwuEnergies;
using math::readBitFromMask;
using synthesis::maximum_valence_work_unit_atoms;
using synthesis::PsSynthesisWriter;
using synthesis::SyValenceKit;
using synthesis::vwu_abstract_length;
using synthesis::VwuAbstractMap;
using synthesis::VwuGoal;
using synthesis::VwuTask;

/// \brief Carry out the instructions in a single valence work unit, as presented in the topology
///        synthesis.  This routine is called by various incarnations of the evalSyValenceEnergy()
///        function.  While any of the local position, velocity, or force data passed into it
///        could be modified, the calling functions may impose const-ness on the global positions,
///        velocities, and forces of various particles.
///
/// \param syvk         Consensus tables of valence parameters and instructions
/// \param sh_xcrd      
/// \param sh_ycrd
/// \param sh_zcrd
/// \param sh_xvel
/// \param sh_yvel
/// \param sh_zvel
/// \param sh_xfrc
/// \param sh_yfrc
/// \param sh_zfrc
/// \param ecard
/// \param vwu_idx
/// \param eval_force
/// \param activity
/// \param purpose
/// \param step_number
template <typename Tcalc>
void synthesisVwuEvaluation(const SyValenceKit<Tcalc> syvk, llint* sh_xcrd, llint* sh_ycrd,
                            llint* sh_zcrd, llint* sh_xvel, llint* sh_yvel, llint* sh_zvel,
                            llint* sh_xfrc, llint* sh_yfrc, llint* sh_zfrc, double inv_gpos_scale,
                            double force_scale, ScoreCard *ecard, int vwu_idx,
                            EvaluateForce eval_force, VwuTask activity, VwuGoal purpose,
                            int step_number);

/// \brief Evaluate all work units in an AtomGraphSynthesis (synthesis of topologies).  The
///        
template <typename Tcalc>
void evalSyValenceEnergy(const SyValenceKit<Tcalc> syvk, PsSynthesisWriter psyw, ScoreCard *ecard,
                         EvaluateForce eval_force, VwuTask activity, VwuGoal purpose,
                         int step_number);

} // namespace energy
} // namespace omni

#include "eval_synthesis.tpp"

#endif
