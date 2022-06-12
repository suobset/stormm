// -*-c++-*-
#ifndef OMNI_VALENCE_POTENTIAL_CUH
#define OMNI_VALENCE_POTENTIAL_CUH

#include "Accelerator/gpu_details.h"
#include "Constants/fixed_precision.h"
#include "MolecularMechanics/mm_controls.h"
#include "Potential/cacheresource.h"
#include "Potential/energy_enumerators.h"
#include "Potential/scorecard.h"
#include "Synthesis/atomgraph_synthesis.h"
#include "Synthesis/phasespace_synthesis.h"
#include "Synthesis/synthesis_abstracts.h"
#include "Synthesis/synthesis_enumerators.h"

#define EXCL_GMEM_OFSET (blockIdx.x * maximum_valence_work_unit_atoms)

namespace omni {
namespace energy {

using card::GpuDetails;
using mm::MMControlKit;
using mm::MolecularMechanicsControls;
using numerics::ForceAccumulationMethod;
using synthesis::AtomGraphSynthesis;
using synthesis::PhaseSpaceSynthesis;
using synthesis::PsSynthesisWriter;
using synthesis::SyValenceKit;
using synthesis::VwuGoal;
  
/// \brief Evaluate valence work units and move atoms.
///
/// Overloaded:
///   - Perform work in single or double precision
///   - Compute forces, energy, or both
///   - Move particles, or instead accumulate forces, energies, or both
///
/// \param poly_vk  Valence parameters based on consensus tables from a topology synthesis
/// \param poly_ps  Coordinates, velocities, and forces of all systems
/// \param gmem_x   Exclusive space in global memory arrays reserved for each thread block, to be
///                 brought into free L1 cache
/// \{
void launchValenceDp(const SyValenceKit<double> &poly_vk, MMControlKit<double> *ctrl,
                     PsSynthesisWriter *poly_psw, ScoreCardWriter *scw,
                     CacheResourceKit<double> *gmem_r, EvaluateForce eval_force,
                     EvaluateEnergy eval_energy, VwuGoal purpose, const GpuDetails &gpu);

void launchValenceDp(const AtomGraphSynthesis &poly_ag, MolecularMechanicsControls *mmctrl,
                     PhaseSpaceSynthesis *poly_ps, ScoreCard *sc, CacheResource *tb_space,
                     EvaluateForce eval_force, EvaluateEnergy eval_energy, VwuGoal purpose,
                     const GpuDetails &gpu);

void launchValenceSp(const SyValenceKit<float> &poly_vk, MMControlKit<float> *ctrl,
                     PsSynthesisWriter *poly_psw, ScoreCardWriter *scw,
                     CacheResourceKit<float> *gmem_r, EvaluateForce eval_force,
                     EvaluateEnergy eval_energy, VwuGoal purpose,
                     ForceAccumulationMethod force_sum, const GpuDetails &gpu);

void launchValenceSp(const AtomGraphSynthesis &poly_ag, MolecularMechanicsControls *mmctrl,
                     PhaseSpaceSynthesis *poly_ps, ScoreCard *sc, CacheResource *tb_space,
                     EvaluateForce eval_force, EvaluateEnergy eval_energy, VwuGoal purpose,
                     ForceAccumulationMethod force_sum, const GpuDetails &gpu);
/// \}
  
} // namespace energy
} // namespace omni

#endif
