// -*-c++-*-
#ifndef OMNI_NONBONDED_POTENTIAL_CUH
#define OMNI_NONBONDED_POTENTIAL_CUH

#include "Accelerator/gpu_details.h"
#include "MolecularMechanics/mm_controls.h"
#include "Synthesis/phasespace_synthesis.h"
#include "Synthesis/static_mask_synthesis.h"
#include "Synthesis/synthesis_abstracts.h"
#include "cacheresource.h"
#include "energy_enumerators.h"
#include "scorecard.h"

namespace omni {
namespace energy {

using card::GpuDetails;
using card::KernelManager;
using mm::MMControlKit;
using mm::MolecularMechanicsControls;
using numerics::ForceAccumulationMethod;
using synthesis::AtomGraphSynthesis;
using synthesis::PhaseSpaceSynthesis;
using synthesis::PsSynthesisWriter;
using synthesis::SeMaskSynthesisReader;
using synthesis::StaticExclusionMaskSynthesis;
using synthesis::SyNonbondedKit;

/// \brief Set the __shared__ memory configuration for various nonbonded interaction kernels
void nonbondedKernelSetup();

/// \brief Obtain information on launch bounds and block-specific requirements for each version of
///        the non-bonded interactions kernel.  Deposit the results in a developing object that
///        will later record launch grid dimensions for managing the kernels.
///
/// \param wisdom  Object to store the kernel specifications obtained
void queryNonbondedKernelRequirements(KernelManager *wisdom);

/// \brief Evaluate nonbonded work units based on tile groups.  All of the kernels launched by
///        these functions will compute forces, energies, or both, and if forces are computed the
///        results will be dumped back into global accumulators.  None of the kernel launchers
///        will move particles.
/// \{
void launchNonbondedTileGroupsDp(const SyNonbondedKit<double> &poly_nbk,
                                 const SeMaskSynthesisReader &poly_ser, MMControlKit<double> *ctrl,
                                 PsSynthesisWriter *poly_psw, ScoreCardWriter *scw,
                                 CacheResourceKit<double> *gmem_r, const EvaluateForce eval_force,
                                 const EvaluateEnergy eval_energy, const GpuDetails &gpu);

void launchNonbondedTileGroupsDp(const AtomGraphSynthesis &poly_ag,
                                 const StaticExclusionMaskSynthesis &poly_se,
                                 MolecularMechanicsControls *mmctrl, PhaseSpaceSynthesis *poly_ps,
                                 ScoreCard *sc, CacheResource *tb_space,
                                 const EvaluateForce eval_force, const EvaluateEnergy eval_energy,
                                 const GpuDetails &gpu);

void launchNonbondedTileGroupsSp(const SyNonbondedKit<float> &poly_nbk,
                                 const SeMaskSynthesisReader &poly_ser, MMControlKit<float> *ctrl,
                                 PsSynthesisWriter *poly_psw, ScoreCardWriter *scw,
                                 CacheResourceKit<float> *gmem_r, const EvaluateForce eval_force,
                                 const EvaluateEnergy eval_energy,
                                 const ForceAccumulationMethod force_sum, const GpuDetails &gpu);

void launchNonbondedTileGroupsSp(const AtomGraphSynthesis &poly_ag,
                                 const StaticExclusionMaskSynthesis &poly_se,
                                 MolecularMechanicsControls *mmctrl, PhaseSpaceSynthesis *poly_ps,
                                 ScoreCard *sc, CacheResource *tb_space,
                                 const EvaluateForce eval_force, const EvaluateEnergy eval_energy,
                                 const ForceAccumulationMethod force_sum, const GpuDetails &gpu);
/// \}

} // namespace energy
} // namespace omni

#endif
