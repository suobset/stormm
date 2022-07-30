// -*-c++-*-
#ifndef OMNI_HPC_MINIMIZATION_H
#define OMNI_HPC_MINIMIZATION_H

#include "Accelerator/gpu_details.h"
#include "Accelerator/kernel_manager.h"
#include "Constants/behavior.h"
#include "Constants/fixed_precision.h"
#include "MolecularMechanics/line_minimization.h"
#include "MolecularMechanics/mm_controls.h"
#include "Potential/cache_resource.h"
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
using Synthesis::StaticExclusionMaskSynthesis;

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
/// \param prec
/// \param poly_ag
/// \param poly_se
/// \param poly_ps
/// \param mmctrl
/// \param sc
/// \param vale_cache
/// \param nonb_cache
/// \param rbg
/// \param line_record
/// \param acc_meth
/// \param launcher
void launchMinimization(const PrecisionModel prec, const AtomGraphSynthesis &poly_ag,
                        const StaticExclusionMaskSynthesis &poly_se,
                        PhaseSpaceSynthesis *poly_ps, MolecularMechanicsControls *mmctrl,
                        ScoreCard *sc, CacheResource *vale_cache, CacheResource *nonb_cache,
                        ReductionBridge *rbg, LineMinimization *line_record,
                        const ForceAccumulationMethod acc_meth, const KernelManager &launcher);

} // namespace mm
} // namespace omni

#endif
