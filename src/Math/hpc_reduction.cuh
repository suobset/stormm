// -*-c++-*-
#ifndef OMNI_HPC_REDUCTION_H
#define OMNI_HPC_REDUCTION_H

#include "Accelerator/kernel_manager.h"
#include "MolecularMechanics/mm_controls.h"
#include "Synthesis/atomgraph_synthesis.h"
#include "Synthesis/phasespace_synthesis.h"
#include "reduction_abstracts.h"
#include "reduction_bridge.h"
#include "reduction_enumerators.h"

namespace omni {
namespace math {

using card::KernelManager;
using constants::PrecisionModel;
using mm::MolecularMechanicsControls;
using mm::MMControlKit;
using synthesis::AtomGraphSynthesis;
using synthesis::PhaseSpaceSynthesis;

/// \brief Obtain the kernel function attributes for one of the reduction kernels.
///
/// \param prec     The precision level at which coordinates and forces are stored--accumulations
///                 by various reduction kernels are always handled in double precision real
///                 numbers, but specifying SINGLE or DOUBLE here indicates whether the kernel
///                 will access overflow arrays for the extended fixed-precision data format.
/// \param purpose  The procedure, e.g. conjugate gradient move vector determination, to be
///                 accomplished by reductions across all systems
/// \param process  Strategy for performing the reduction, one of two stages or all together
cudaFuncAttributes queryReductionKernelRequirements(const PrecisionModel prec,
                                                    const ReductionGoal purpose,
                                                    const ReductionStage process);

/// \brief 
void launchConjugateGradientDp(const ReductionKit &redk, ConjGradSubstrate *cgsbs,
                               MMControlKit<double> *ctrl, const KernelManager &launcher);

/// \brief
void launchConjugateGradientDp(const AtomGraphSynthesis poly_ag, PhaseSpaceSynthesis *poly_ps,
                               ReductionBridge *rbg, MolecularMechanicsControls *mmctrl,
                               const KernelManager &launcher);

/// \brief
void launchConjugateGradientSp(const ReductionKit &redk, ConjGradSubstrate *cgsbs,
                               MMControlKit<float> *ctrl, const KernelManager &launcher);

/// \brief
void launchConjugateGradientSp(const AtomGraphSynthesis poly_ag, PhaseSpaceSynthesis *poly_ps,
                               ReductionBridge *rbg, MolecularMechanicsControls *mmctrl,
                               const KernelManager &launcher);
  
} // synthesis
} // omni

#endif
