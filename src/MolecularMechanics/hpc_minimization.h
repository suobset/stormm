// -*-c++-*-
#ifndef STORMM_HPC_MINIMIZATION_H
#define STORMM_HPC_MINIMIZATION_H

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

namespace stormm {
namespace mm {

using card::GpuDetails;
using card::KernelManager;
using constants::PrecisionModel;
using energy::CacheResource;
using energy::ScoreCard;
using energy::ScoreCardWriter;
using math::ReductionBridge;
using math::ReductionKit;
using numerics::ForceAccumulationMethod;
using synthesis::AtomGraphSynthesis;
using synthesis::PhaseSpaceSynthesis;
using synthesis::PsSynthesisWriter;
using synthesis::StaticExclusionMaskSynthesis;

/// \brief Set the __shared__ memory access size for the conjugate gradient particle advancement
///        kernels.
void minimizationKernelSetup();

/// \brief Get the launch parameters for conjugate gradient particle advancement kernels.
///
/// \param prec  The precision model of the coordinate object (the particle advancement is always
///              performed in double-precision real numbers)
cudaFuncAttributes queryMinimizationKernelRequirements(const PrecisionModel prec);

/// \brief Launch the line minimization kernel in a standalone capacity.  This function expedites
///        debugging and provides abstraction of the line movement function for other applications.
///
/// \param prec         The precision model in use by the coordinate set (indicates whether the
///                     extended fixed-precision number format is in use)
/// \param poly_psw     Coordinates for each system and forces acting on each system
/// \param redk         Reduction work units to inform the limits of each thread block's activity
/// \param scw          Read-only access to the current energy of each system
/// \param lmw          Writeable abstract for a line minimization manager
/// \param move_number  The number of the move--values of 0 to 3 will lead to different behavior
/// \param redu_lp      Launch parameters for the kernel (an abstract of the KernelManager object)
/// \{
void launchLineAdvance(PrecisionModel prec, PsSynthesisWriter *poly_psw, const ReductionKit &redk,
                       const ScoreCardWriter &scw, LinMinWriter *lmw, int move_number,
                       const int2 redu_lp);

void launchLineAdvance(PrecisionModel prec, PhaseSpaceSynthesis *poly_ps,
                       const AtomGraphSynthesis &poly_ag, ScoreCard *sc,
                       LineMinimization *line_record, int move_number,
                       const KernelManager &launcher);
/// \}

/// \brief Run energy minimizations of all structures in a synthesis based on user input or other
///        data contained in a molecular mechanics control object.
///
/// \param prec           The precision to use in general arithmetic.  Also activates an extended
///                       mode of fixed-precision data representation in the coordinates object.
/// \param poly_ag        Compilation of all systems' topologies (a "synthesis" of systems)
/// \param poly_se        Synthesis of exclusion masks, tailored to the topology synthesis
/// \param poly_ps        Compilation of coordinates, forces, and other useful arrays for
///                       minimizing each structure
/// \param mmctrl_fe      Control data for the minimization run, plus progress counters for GPU
///                       kernel management in force + energy computations
/// \param mmctrl_xe      Control data for the minimization run, plus progress counters for GPU
///                       kernel management in energy-only computations
/// \param sc             Energy tracking object allocated to cover all systems in the synthesis
/// \param vale_fe_cache  Pre-allocated scratch space for valence kernel thread blocks in force +
///                       energy computations
/// \param nonb_fe_cache  Pre-allocated scratch space for non-bonded kernel thread blocks in
///                       energy-only computations
/// \param vale_xe_cache  Pre-allocated scratch space for valence kernel thread blocks in force +
///                       energy computations
/// \param nonb_xe_cache  Pre-allocated scratch space for non-bonded kernel thread blocks in
///                       energy-only computations
/// \param rbg            Pre-allocated storage space for reduction work units
/// \param line_record    Stores line minimization progress, energies and move lengths
/// \param acc_meth       Choice of force accumulation method
/// \param gpu            Details of the GPU in use
/// \param launcher       Holds parameters for all kernels, will be accessed for the correct
///                       launch configurations of the kernels relevant to the workload at hand
void launchMinimization(const PrecisionModel prec, const AtomGraphSynthesis &poly_ag,
                        const StaticExclusionMaskSynthesis &poly_se,
                        PhaseSpaceSynthesis *poly_ps, MolecularMechanicsControls *mmctrl_fe,
                        MolecularMechanicsControls *mmctrl_xe, ScoreCard *sc,
                        CacheResource *vale_fe_cache, CacheResource *nonb_fe_cache,
                        CacheResource *vale_xe_cache, CacheResource *nonb_xe_cache,
                        ReductionBridge *rbg, LineMinimization *line_record,
                        const ForceAccumulationMethod acc_meth, const GpuDetails &gpu,
                        const KernelManager &launcher);

} // namespace mm
} // namespace stormm

#endif
