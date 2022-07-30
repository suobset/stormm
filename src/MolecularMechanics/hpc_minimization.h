// -*-c++-*-
#ifndef OMNI_HPC_MINIMIZATION_H
#define OMNI_HPC_MINIMIZATION_H

#include "Accelerator/gpu_details.h"
#include "Accelerator/kernel_manager.h"
#include "Constants/behavior.h"
#include "Constants/fixed_precision.h"
#include "MolecularMechanics/line_minimization.h"
#include "MolecularMechanics/mm_controls.h"
#include "Potential/cacheresource.h"
#include "Potential/scorecard.h"
#include "Math/reduction_bridge.h"
#include "Synthesis/atomgraph_synthesis.h"
#include "Synthesis/phasespace_synthesis.h"
#include "Synthesis/static_mask_synthesis.h"

namespace omni {
namespace mm {

using card::GpuDetails;
using card::KernelManager;
using constants::PrecisionModel;
using energy::CacheResource;
using energy::ScoreCard;
using math::ReductionBridge;
using mm::MolecularMechanicsControls;
using mm::LineMinimization;
using numerics::ForceAccumulationMethod;
using synthesis::AtomGraphSynthesis;
using synthesis::PhaseSpaceSynthesis;
using synthesis::StaticExclusionMaskSynthesis;

/// \brief Set the __shared__ memory access size for the conjugate gradient particle advancement
///        kernels.
void minimizationKernelSetup();

/// \brief Get the launch parameters for conjugate gradient particle advancement kernels.
///
/// \param prec  The precision model of the coordinate object (the particle advancement is always
///              performed in double-precision real numbers)
cudaFuncAttributes queryMinimizationKernelRequirements(const PrecisionModel prec);

/// \brief Run energy minimizations of all structures in a synthesis based on user input or other
///        data contained in a molecular mechanics control object.
///
/// \param prec         The precision to use in general arithmetic.  Also activates an extended
///                     mode of fixed-precision data representation in the coordinates object.
/// \param poly_ag      Compilation of all systems' topologies (a "synthesis" of systems)
/// \param poly_se      Synthesis of exclusion masks, tailored to the topology synthesis
/// \param poly_ps      Compilation of coordinates, forces, and other useful arrays for minimizing
///                     each structure
/// \param mmctrl       Control data for the minimization run, plus progress counters for GPU
///                     kernel management
/// \param sc           Energy tracking object allocated to cover all systems in the synthesis
/// \param vale_cache   Pre-allocated scratch space for valence kernel thread blocks
/// \param nonb_cache   Pre-allocated scratch space for non-bonded kernel thread blocks
/// \param rbg          Pre-allocated storage space for reduction work units
/// \param line_record  Stores line minimization progress, energies and move lengths
/// \param acc_meth     Choice of force accumulation method
/// \param gpu          Details of the GPU in use
/// \param launcher     Holds parameters for all kernels, will be accessed for the correct launch
///                     configurations of the kernels relevant to the workload at hand
void launchMinimization(const PrecisionModel prec, const AtomGraphSynthesis &poly_ag,
                        const StaticExclusionMaskSynthesis &poly_se,
                        PhaseSpaceSynthesis *poly_ps, MolecularMechanicsControls *mmctrl,
                        ScoreCard *sc, CacheResource *vale_cache, CacheResource *nonb_cache,
                        ReductionBridge *rbg, LineMinimization *line_record,
                        const ForceAccumulationMethod acc_meth, const GpuDetails &gpu,
                        const KernelManager &launcher);

} // namespace mm
} // namespace omni

#endif
