// -*-c++-*-
#ifndef OMNI_VALENCE_POTENTIAL_CUH
#define OMNI_VALENCE_POTENTIAL_CUH

#include "Accelerator/gpu_details.h"
#include "Accelerator/kernel_manager.h"
#include "Constants/fixed_precision.h"
#include "DataTypes/omni_vector_types.h"
#include "MolecularMechanics/mm_controls.h"
#include "Synthesis/atomgraph_synthesis.h"
#include "Synthesis/phasespace_synthesis.h"
#include "Synthesis/synthesis_abstracts.h"
#include "Synthesis/synthesis_enumerators.h"
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
using synthesis::SyValenceKit;
using synthesis::SyRestraintKit;
using synthesis::VwuGoal;

/// \brief Set the __shared__ memory configuration for various valence interaction kernels.
void valenceKernelSetup();

/// \brief Test various subdivisions of a thread block to see whether the workload can be better
///        distributed.  This is perilous if the kernel allocates significant static __shared__
///        memory resources per block, but can be advantageous if the systems are small and there
///        are many to compute.  This assumes that the kernel is compiled with a minimum of one
///        block per multiprocessor in its launch bounds, and the highest possible thread count
///        for that one block.  The function returns a tuple containing the best identified
///        block multiplier and the best corresponding block size, i.e. "split a kernel of 896
///        maximum threads into 7 blocks of 128."
///
/// \param max_threads  The maximum number of threads for which the kernel is compiled
/// \param smp_count    The number of streaming multiprocessors on the device
/// \param vwu_size     Maximum size of each valence work unit
/// \param vwu_count    The total number of valence work units that the kernel will process
int2 testValenceKernelSubdivision(const int max_threads, const int smp_count, const int vwu_size,
                                  const int vwu_count);

/// \brief Obtain information on launch bounds and block-specific requirements for each version of
///        the valence interactions kernel.  Deposit the results in a developing object that will
///        later record launch grid dimensions for managing the kernels.
///
/// \param wisdom   Object to store the kernel specifications obtained
/// \param poly_ag  A collection of topologies describing the workload
void queryValenceKernelRequirements(KernelManager *wisdom, const AtomGraphSynthesis &poly_ag);

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
void launchValenceDp(const SyValenceKit<double> &poly_vk,
                     const SyRestraintKit<double, double2, double4> &poly_rk,
                     MMControlKit<double> *ctrl, PsSynthesisWriter *poly_psw, ScoreCardWriter *scw,
                     CacheResourceKit<double> *gmem_r, EvaluateForce eval_force,
                     EvaluateEnergy eval_energy, VwuGoal purpose, const KernelManager &launcher);

void launchValenceDp(const AtomGraphSynthesis &poly_ag, MolecularMechanicsControls *mmctrl,
                     PhaseSpaceSynthesis *poly_ps, ScoreCard *sc, CacheResource *tb_space,
                     EvaluateForce eval_force, EvaluateEnergy eval_energy, VwuGoal purpose,
                     const KernelManager &launcher);

void launchValenceSp(const SyValenceKit<float> &poly_vk,
                     const SyRestraintKit<float, float2, float4> &poly_rk,
                     MMControlKit<float> *ctrl, PsSynthesisWriter *poly_psw, ScoreCardWriter *scw,
                     CacheResourceKit<float> *gmem_r, EvaluateForce eval_force,
                     EvaluateEnergy eval_energy, VwuGoal purpose,
                     ForceAccumulationMethod force_sum, const KernelManager &launcher);
  
void launchValenceSp(const AtomGraphSynthesis &poly_ag, MolecularMechanicsControls *mmctrl,
                     PhaseSpaceSynthesis *poly_ps, ScoreCard *sc, CacheResource *tb_space,
                     EvaluateForce eval_force, EvaluateEnergy eval_energy, VwuGoal purpose,
                     ForceAccumulationMethod force_sum, const KernelManager &launcher);
/// \}
  
} // namespace energy
} // namespace omni

#endif
