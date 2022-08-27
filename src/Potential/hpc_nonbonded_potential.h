// -*-c++-*-
#ifndef STORMM_NONBONDED_POTENTIAL_CUH
#define STORMM_NONBONDED_POTENTIAL_CUH

#include "copyright.h"
#include "Accelerator/gpu_details.h"
#include "Accelerator/kernel_manager.h"
#include "Constants/behavior.h"
#include "MolecularMechanics/mm_controls.h"
#include "Potential/energy_enumerators.h"
#include "Synthesis/implicit_solvent_workspace.h"
#include "Synthesis/phasespace_synthesis.h"
#include "Synthesis/static_mask_synthesis.h"
#include "Synthesis/synthesis_abstracts.h"
#include "Synthesis/synthesis_enumerators.h"
#include "cacheresource.h"
#include "energy_enumerators.h"
#include "scorecard.h"

namespace stormm {
namespace energy {

using card::GpuDetails;
using card::KernelManager;
using constants::PrecisionModel;
using mm::MMControlKit;
using mm::MolecularMechanicsControls;
using numerics::AccumulationMethod;
using synthesis::AtomGraphSynthesis;
using synthesis::ImplicitSolventWorkspace;
using synthesis::ISWorkspaceKit;
using synthesis::NbwuKind;
using synthesis::PhaseSpaceSynthesis;
using synthesis::PsSynthesisWriter;
using synthesis::SeMaskSynthesisReader;
using synthesis::StaticExclusionMaskSynthesis;
using synthesis::SyNonbondedKit;
using synthesis::VwuGoal;

/// \brief Set the __shared__ memory configuration for various nonbonded interaction kernels
void nonbondedKernelSetup();

/// \brief Obtain information on launch bounds and block-specific requirements for each version of
///        the non-bonded interactions kernel.  Deposit the results in a developing object that
///        will later record launch grid dimensions for managing the kernels.
///
/// \param prec      The desired precision model for the kernel
/// \param kind      The type of non-bonded work units to evaluate
/// \param eval_frc  Whether the kernel shall evaluate forces
/// \param eval_nrg  Whether the kernel shall evaluate energies
/// \param acc_meth  Accumulation method for forces
cudaFuncAttributes queryNonbondedKernelRequirements(PrecisionModel prec, NbwuKind kind,
                                                    EvaluateForce eval_frc,
                                                    EvaluateEnergy eval_nrg,
                                                    AccumulationMethod acc_meth);

/// \brief Obtain information on launch bounds and block-specific requirements for each version of
///        the Born radii computation kernel.  Deposit the results in a developing object that
///        will later record launch grid dimensions for managing the kernels.  This kernel's
///        parameter descriptions follow from queryNonbondedKernelRequirements() above.
cudaFuncAttributes queryBornRadiiKernelRequirements(PrecisionModel prec, NbwuKind kind,
                                                    AccumulationMethod acc_meth);

/// \brief Obtain information on launch bounds and block-specific requirements for each version of
///        the Born radii derivative computation kernel.  Deposit the results in a developing
///        object that will later record launch grid dimensions for managing the kernels.  This
///        kernel's parameter descriptions follow from queryNonbondedKernelRequirements() above.
cudaFuncAttributes queryBornDerivativeKernelRequirements(PrecisionModel prec, NbwuKind kind,
                                                         AccumulationMethod acc_meth);

/// \brief Evaluate Generalized Born radii for the current state.  Prepare for Born radii
///        derivative calculations.  The accumulation of the radii and their derivatives is
///        handled in fixed-precision arithmetic, using the Cartesian X and Y force components of
///        the thread block workspace as accumulators.
///
/// Overloaded:
///   - Accept abstracts that imply the precision level of kernel to use, and launch parameters
///     assumed to fit that kernel.
///   - Accept the original objects and make the necessary abstracts before launching (this is
///     useful for testing purposes)
///
/// \param kind         The type non-bonded work to perform, indicating the kernel to launch
/// \param poly_nbk     Non-bonded parameters of all systems
/// \param poly_ag      Compiled topologies of all systems
/// \param ngb_kit      Neck Generalized Born tables for certain GB models, translated to a
///                     specific precision model
/// \param ngb_tab      Neck Generalized Born tables for certain GB models
/// \param ctrl         Abstract for molecular mechanics progress counters and run bounds
/// \param mmctrl       Progress counters and run bounds
/// \param poly_psw     Abstract for coordinates and forces of all systems
/// \param poly_ps      Coordinate and force compilation for all systems
/// \param gmem_r       Abstract for thread block specific resources
/// \param tb_space     Cache resources for the kernel launch
/// \param iswk         Abstract for the implicit solvent accumulator workspace
/// \param isw          Implicit solvent accumulators for connecting each of the Generalized Born
///                     calculation stages
/// \param acc_meth     Accumulation method for Born radii
/// \param bt           Block and thread counts for the kernel, given the precision and force or
///                     energy computation requirements (an "abstract" of the KernelManager object)
/// \param launcher     Contains launch parameters for all kernels in STORMM
/// \{
void launchBornRadiiCalculation(NbwuKind kind, const SyNonbondedKit<double, double2> &poly_nbk,
                                MMControlKit<double> *ctrl, PsSynthesisWriter *poly_psw,
                                CacheResourceKit<double> *gmem_r, ISWorkspaceKit<double> *iswk,
                                const int2 bt);

void launchBornRadiiCalculation(NbwuKind kind, const SyNonbondedKit<float, float2> &poly_nbk,
                                MMControlKit<float> *ctrl, PsSynthesisWriter *poly_psw,
                                CacheResourceKit<float> *gmem_r, ISWorkspaceKit<float> *iswk,
                                AccumulationMethod acc_meth, const int2 bt);

void launchBornRadiiCalculation(PrecisionModel prec, const AtomGraphSynthesis &poly_ag,
                                MolecularMechanicsControls *mmctrl, PhaseSpaceSynthesis *poly_ps,
                                CacheResource *tb_space, ImplicitSolventWorkspace *isw,
                                const KernelManager &launcher,
                                AccumulationMethod acc_meth = AccumulationMethod::AUTOMATIC);
/// \}

/// \brief Evaluate nonbonded work units based on the strategy determined in the topology
///        synthesis.  All of the kernels launched by these functions will compute forces,
///        energies, or both, and if forces are computed the results will be dumped back into
///        global accumulators.  None of the kernel launcher will move particles.
///
/// Overloaded:
///   - Accept abstracts and launch parameters for the specific kernel (this is slightly more
///     performant by allowing re-use of the abstracts and parameters if they are created at the
///     start of some loop)
///   - Accept the original objects and make the necessary abstracts before launching (this is
///     useful for testing purposes)
///
/// \param kind         The type non-bonded work to perform, indicating the kernel to launch
/// \param poly_nbk     Non-bonded parameters of all systems
/// \param poly_ag      Compiled topologies of all systems
/// \param poly_ser     Abstract for the exclusion masks of all systems
/// \param poly_se      Exclusion masks for the synthesis of topologies
/// \param ctrl         Abstract for molecular mechanics progress counters and run bounds
/// \param mmctrl       Progress counters and run bounds
/// \param poly_psw     Abstract for coordinates and forces of all systems
/// \param poly_ps      Coordinate and force compilation for all systems
/// \param scw          Abstract for the energy tracking on all systems
/// \param sc           Energy tracking for all systems (this must be pre-sized to accommodate the
///                     entire synthesis)
/// \param gmem_r       Abstract for thread block specific resources
/// \param tb_space     Cache resources for the kernel launch
/// \param eval_force   Whether to evaluate forces
/// \param eval_energy  Whether to evaluate energies
/// \param acc_meth     Accumulation method for forces
/// \param bt           Block and thread counts for the kernel, given the precision and force or
///                     energy computation requirements (an "abstract" of the KernelManager object)
/// \param launcher     Contains launch parameters for all kernels in STORMM
/// \{
void launchNonbonded(NbwuKind kind, const SyNonbondedKit<double, double2> &poly_nbk,
                     const SeMaskSynthesisReader &poly_ser, MMControlKit<double> *ctrl,
                     PsSynthesisWriter *poly_psw, ScoreCardWriter *scw,
                     CacheResourceKit<double> *gmem_r, const EvaluateForce eval_force,
                     const EvaluateEnergy eval_energy, const int2 bt);

void launchNonbonded(NbwuKind kind, const SyNonbondedKit<float, float2> &poly_nbk,
                     const SeMaskSynthesisReader &poly_ser, MMControlKit<float> *ctrl,
                     PsSynthesisWriter *poly_psw, ScoreCardWriter *scw,
                     CacheResourceKit<float> *gmem_r, EvaluateForce eval_force,
                     EvaluateEnergy eval_energy, AccumulationMethod force_sum, const int2 bt);

void launchNonbonded(PrecisionModel prec, const AtomGraphSynthesis &poly_ag,
                     const StaticExclusionMaskSynthesis &poly_se,
                     MolecularMechanicsControls *mmctrl, PhaseSpaceSynthesis *poly_ps,
                     ScoreCard *sc, CacheResource *tb_space, EvaluateForce eval_force,
                     EvaluateEnergy eval_energy, AccumulationMethod force_sum,
                     const KernelManager &launcher);

void launchNonbonded(PrecisionModel prec, const AtomGraphSynthesis &poly_ag,
                     const StaticExclusionMaskSynthesis &poly_se,
                     MolecularMechanicsControls *mmctrl, PhaseSpaceSynthesis *poly_ps,
                     ScoreCard *sc, CacheResource *tb_space, EvaluateForce eval_force,
                     EvaluateEnergy eval_energy, const KernelManager &launcher);
/// \}

} // namespace energy
} // namespace stormm

#endif
