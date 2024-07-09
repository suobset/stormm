// -*-c++-*-
#ifndef STORMM_HPC_PME_POTENTIAL_H

#include "copyright.h"
#include "Constants/behavior.h"
#include "MolecularMechanics/mm_controls.h"
#include "Synthesis/synthesis_abstracts.h"
#include "energy_enumerators.h"
#include "local_exclusionmask.h"
#include "pme_potential.h"
#include "ppitable.h"
#include "scorecard.h"
#include "tile_manager.h"

namespace stormm {
namespace energy {

using constants::PrecisionModel;
using mm::MMControlKit;
using mm::MolecularMechanicsControls;
using synthesis::SyNonbondedKit;

/// \brief Return critical attributes of a selected particle-particle pair interactions kernel
///        based on selected features.
///
/// \param prec
/// \param eval_frc
/// \param eval_nrg
/// \param neighbor_list
/// \param has_tiny_box        
/// \param collision_handling
cudaFuncAttributes queryPmePairsKernelRequirements(PrecisionModel prec, EvaluateForce eval_frc,
                                                   EvaluateEnergy eval_nrg,
                                                   NeighborListKind neighbor_list,
                                                   TinyBoxPresence has_tiny_box,
                                                   ClashResponse collision_handling);

/// \brief Launch the appropriate kernel to evaluate particle-particle pair interactions in a
///        neighbor list for periodic simulations.
///
/// Overloaded:
///   - Supply abstracts, differentiating single- and double-precision calculations
///   - Supply the original objects, differentiating a unified neighbor list from separated
///     neighbor lists
///
/// \param poly_nbk
/// \param lemr
/// \param tlpn
/// \param nrg_tab
/// \param scw
/// \param cgw
/// \param cgw_qq
/// \param cgw_lj
/// \param ctrl
/// \param eval_frc
/// \param eval_nrg
/// \param has_tiny_box
/// \param bt
/// \param clash_distance
/// \param clash_ratio
/// \{
void launchPmePairs(const SyNonbondedKit<double, double2> &poly_nbk,
                    const LocalExclusionMaskReader &lemr, const TilePlan &tlpn,
                    const PPIKit<double, double4> &nrg_tab, ScoreCardWriter *scw,
                    CellGridWriter<double, llint, double, double4> *cgw,
                    MMControlKit<double> *ctrl, EvaluateForce eval_frc, EvaluateEnergy eval_nrg,
                    TinyBoxPresence has_tiny_box, const int2 bt, double clash_distance,
                    double clash_ratio);

void launchPmePairs(const SyNonbondedKit<double, double2> &poly_nbk,
                    const LocalExclusionMaskReader &lemr, const TilePlan &tlpn,
                    const PPIKit<double, double4> &nrg_tab, ScoreCardWriter *scw,
                    CellGridWriter<double, llint, double, double4> *cgw_qq,
                    CellGridWriter<double, llint, double, double4> *cgw_lj,
                    MMControlKit<double> *ctrl, EvaluateForce eval_frc, EvaluateEnergy eval_nrg,
                    TinyBoxPresence has_tiny_box, const int2 bt, double clash_distance,
                    double clash_ratio);

void launchPmePairs(const SyNonbondedKit<float, float2> &poly_nbk,
                    const LocalExclusionMaskReader &lemr, const TilePlan &tlpn,
                    const PPIKit<float, float4> &nrg_tab, ScoreCardWriter *scw,
                    CellGridWriter<float, int, float, float4> *cgw,
                    MMControlKit<float> *ctrl, EvaluateForce eval_frc, EvaluateEnergy eval_nrg,
                    TinyBoxPresence has_tiny_box, const int2 bt, double clash_distance,
                    double clash_ratio);

void launchPmePairs(const SyNonbondedKit<float, float2> &poly_nbk,
                    const LocalExclusionMaskReader &lemr, const TilePlan &tlpn,
                    const PPIKit<float, float4> &nrg_tab, ScoreCardWriter *scw,
                    CellGridWriter<float, int, float, float4> *cgw_qq,
                    CellGridWriter<float, int, float, float4> *cgw_lj,
                    MMControlKit<float> *ctrl, EvaluateForce eval_frc, EvaluateEnergy eval_nrg,
                    TinyBoxPresence has_tiny_box, const int2 bt, double clash_distance,
                    double clash_ratio);
/// \}

} // namespace energy
} // namespace stormm

#endif
