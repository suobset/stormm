// -*-c++-*-
#include "copyright.h"
#include "Constants/behavior.h"
#include "Constants/hpc_bounds.h"
#include "Constants/fixed_precision.h"
#include "DataTypes/common_types.h"
#include "Math/vector_ops.h"
#include "Parsing/parse.h"
#include "Potential/energy_enumerators.h"
#ifdef STORMM_USE_HPC
#include "Math/hpc_reduction.h"
#include "Math/reduction_workunit.h"
#include "Potential/cellgrid.h"
#include "Potential/hpc_map_density.h"
#include "Potential/hpc_nonbonded_potential.h"
#include "Potential/hpc_pme_potential.h"
#include "Potential/hpc_valence_potential.h"
#include "Structure/hpc_rmsd.h"
#include "Structure/hpc_virtual_site_handling.h"
#include "Trajectory/hpc_integration.h"
#endif
#include "Synthesis/valence_workunit.h"
#include "core_kernel_manager.h"

namespace stormm {
namespace card {

using constants::ExceptionResponse;
#ifdef STORMM_USE_HPC
using energy::queryBornRadiiKernelRequirements;
using energy::queryBornDerivativeKernelRequirements;
using energy::queryGeneralQMapKernelRequirements;
using energy::queryMigrationKernelRequirements;
using energy::queryNonbondedKernelRequirements;
using energy::queryPMEPairsKernelRequirements;
using energy::queryShrAccQMapKernelRequirements;
using energy::queryValenceKernelRequirements;
using stmath::optReductionKernelSubdivision;
using stmath::queryReductionKernelRequirements;
using structure::queryRMSDKernelRequirements;
using structure::queryVirtualSiteKernelRequirements;
using trajectory::queryIntegrationKernelRequirements;
#endif
using numerics::globalpos_scale_nonoverflow_bits;
using stmath::findBin;
using synthesis::maximum_valence_work_unit_atoms;
using synthesis::half_valence_work_unit_atoms;
using synthesis::quarter_valence_work_unit_atoms;
using synthesis::minimum_valence_work_unit_atoms;
using parse::CaseSensitivity;
using parse::strcmpCased;

//-------------------------------------------------------------------------------------------------
CoreKlManager::CoreKlManager(const GpuDetails &gpu_in, const AtomGraphSynthesis *poly_ag) :
    KernelManager(gpu_in),
    valence_kernel_width{poly_ag->getValenceThreadBlockSize()},
    nonbond_block_multiplier_dp{nonbondedBlockMultiplier(gpu_in, poly_ag->getUnitCellType(),
                                                         PrecisionModel::DOUBLE,
                                                         poly_ag->getImplicitSolventModel())},
    nonbond_block_multiplier_sp{nonbondedBlockMultiplier(gpu_in, poly_ag->getUnitCellType(),
                                                         PrecisionModel::SINGLE,
                                                         poly_ag->getImplicitSolventModel())},
    gbradii_block_multiplier_dp{gbRadiiBlockMultiplier(gpu_in, PrecisionModel::DOUBLE)},
    gbradii_block_multiplier_sp{gbRadiiBlockMultiplier(gpu_in, PrecisionModel::SINGLE)},
    gbderiv_block_multiplier_dp{gbDerivativeBlockMultiplier(gpu_in, PrecisionModel::DOUBLE)},
    gbderiv_block_multiplier_sp{gbDerivativeBlockMultiplier(gpu_in, PrecisionModel::SINGLE)},
    gen_qmap_block_multiplier_dp{densityMappingBlockMultiplier(gpu_in, PrecisionModel::DOUBLE,
                                                               QMapMethod::GENERAL_PURPOSE)},
    gen_qmap_block_multiplier_sp{densityMappingBlockMultiplier(gpu_in, PrecisionModel::SINGLE,
                                                               QMapMethod::GENERAL_PURPOSE)},
    sac_qmap_block_multiplier_dp{densityMappingBlockMultiplier(gpu_in, PrecisionModel::DOUBLE,
                                                               QMapMethod::ACC_SHARED)},
    sac_qmap_block_multiplier_sp{densityMappingBlockMultiplier(gpu_in, PrecisionModel::SINGLE,
                                                               QMapMethod::ACC_SHARED)},
    reduction_block_multiplier{reductionBlockMultiplier()},
    virtual_site_kernel_width{poly_ag->getValenceThreadBlockSize()},
    rmsd_block_multiplier_dp{rmsdBlockMultiplier(PrecisionModel::DOUBLE)},
    rmsd_block_multiplier_sp{rmsdBlockMultiplier(PrecisionModel::SINGLE)}
{
#ifdef STORMM_USE_HPC
  const std::vector<ClashResponse> clash_policy = { ClashResponse::NONE, ClashResponse::FORGIVE };
  const std::vector<std::string> clash_ext = { "", "NonClash" };
  const std::string valence_kwidth_ext = getEnumerationName(valence_kernel_width);
  for (int i = 0; i < clash_policy.size(); i++) {
    const std::string i_ext = clash_ext[i] + valence_kwidth_ext;

    // Valence kernel entries
    catalogValenceKernel(PrecisionModel::DOUBLE, EvaluateForce::NO, EvaluateEnergy::YES,
                         AccumulationMethod::SPLIT, VwuGoal::ACCUMULATE, clash_policy[i],
                         valence_kernel_width, "kdsValenceEnergyAccumulation" + i_ext);
    catalogValenceKernel(PrecisionModel::DOUBLE, EvaluateForce::YES, EvaluateEnergy::NO,
                         AccumulationMethod::SPLIT, VwuGoal::MOVE_PARTICLES, clash_policy[i],
                         valence_kernel_width, "kdsValenceAtomUpdate" + i_ext);
    catalogValenceKernel(PrecisionModel::DOUBLE, EvaluateForce::YES, EvaluateEnergy::NO,
                         AccumulationMethod::SPLIT, VwuGoal::ACCUMULATE, clash_policy[i],
                         valence_kernel_width, "kdsValenceForceAccumulation" + i_ext);
    catalogValenceKernel(PrecisionModel::DOUBLE, EvaluateForce::YES, EvaluateEnergy::YES,
                         AccumulationMethod::SPLIT, VwuGoal::MOVE_PARTICLES, clash_policy[i],
                         valence_kernel_width, "kdsValenceEnergyAtomUpdate" + i_ext);
    catalogValenceKernel(PrecisionModel::DOUBLE, EvaluateForce::YES, EvaluateEnergy::YES,
                         AccumulationMethod::SPLIT, VwuGoal::ACCUMULATE, clash_policy[i],
                         valence_kernel_width, "kdsValenceForceEnergyAccumulation" + i_ext);
    catalogValenceKernel(PrecisionModel::SINGLE, EvaluateForce::NO, EvaluateEnergy::YES,
                         AccumulationMethod::SPLIT, VwuGoal::ACCUMULATE, clash_policy[i],
                         valence_kernel_width, "kfsValenceEnergyAccumulation" + i_ext);
    catalogValenceKernel(PrecisionModel::SINGLE, EvaluateForce::YES, EvaluateEnergy::NO,
                         AccumulationMethod::SPLIT, VwuGoal::MOVE_PARTICLES, clash_policy[i],
                         valence_kernel_width, "kfsValenceAtomUpdate" + i_ext);
    catalogValenceKernel(PrecisionModel::SINGLE, EvaluateForce::YES, EvaluateEnergy::NO,
                         AccumulationMethod::WHOLE, VwuGoal::MOVE_PARTICLES, clash_policy[i],
                         valence_kernel_width, "kfValenceAtomUpdate" + i_ext);
    catalogValenceKernel(PrecisionModel::SINGLE, EvaluateForce::YES, EvaluateEnergy::NO,
                         AccumulationMethod::SPLIT, VwuGoal::ACCUMULATE, clash_policy[i],
                         valence_kernel_width, "kfsValenceForceAccumulation" + i_ext);
    catalogValenceKernel(PrecisionModel::SINGLE, EvaluateForce::YES, EvaluateEnergy::NO,
                         AccumulationMethod::WHOLE, VwuGoal::ACCUMULATE, clash_policy[i],
                         valence_kernel_width, "kfValenceForceAccumulation" + i_ext);
    catalogValenceKernel(PrecisionModel::SINGLE, EvaluateForce::YES, EvaluateEnergy::YES,
                         AccumulationMethod::WHOLE, VwuGoal::MOVE_PARTICLES, clash_policy[i],
                         valence_kernel_width, "kfValenceEnergyAtomUpdate" + i_ext);
    catalogValenceKernel(PrecisionModel::SINGLE, EvaluateForce::YES, EvaluateEnergy::YES,
                         AccumulationMethod::SPLIT, VwuGoal::MOVE_PARTICLES, clash_policy[i],
                         valence_kernel_width, "kfsValenceEnergyAtomUpdate" + i_ext);
    catalogValenceKernel(PrecisionModel::SINGLE, EvaluateForce::YES, EvaluateEnergy::YES,
                         AccumulationMethod::WHOLE, VwuGoal::ACCUMULATE, clash_policy[i],
                         valence_kernel_width, "kfValenceForceEnergyAccumulation" + i_ext);
    catalogValenceKernel(PrecisionModel::SINGLE, EvaluateForce::YES, EvaluateEnergy::YES,
                         AccumulationMethod::SPLIT, VwuGoal::ACCUMULATE, clash_policy[i],
                         valence_kernel_width, "kfsValenceForceEnergyAccumulation" + i_ext);
#if 0
    catalogValenceKernel(
#endif
    // Non-bonded kernel entries
    const std::vector<ImplicitSolventModel> is_models = { ImplicitSolventModel::NONE,
                                                          ImplicitSolventModel::HCT_GB,
                                                          ImplicitSolventModel::NECK_GB };
    const std::vector<std::string> is_model_names = { "Vacuum", "GB", "GBNeck" };
    for (int j = 0; j < 3; j++) {
      switch (poly_ag->getUnitCellType()) {
      case UnitCellType::NONE:
        catalogNonbondedKernel(PrecisionModel::DOUBLE, NbwuKind::TILE_GROUPS, EvaluateForce::NO,
                               EvaluateEnergy::YES, AccumulationMethod::SPLIT, is_models[j],
                               clash_policy[i],
                               "ktgd" + is_model_names[j] + "Energy" + clash_ext[i]);
        catalogNonbondedKernel(PrecisionModel::DOUBLE, NbwuKind::TILE_GROUPS, EvaluateForce::YES,
                               EvaluateEnergy::NO, AccumulationMethod::SPLIT, is_models[j],
                               clash_policy[i],
                               "ktgds" + is_model_names[j] + "Force" + clash_ext[i]);
        catalogNonbondedKernel(PrecisionModel::DOUBLE, NbwuKind::TILE_GROUPS, EvaluateForce::YES,
                               EvaluateEnergy::YES, AccumulationMethod::SPLIT, is_models[j],
                               clash_policy[i],
                               "ktgds" + is_model_names[j] + "ForceEnergy" + clash_ext[i]);
        catalogNonbondedKernel(PrecisionModel::SINGLE, NbwuKind::TILE_GROUPS, EvaluateForce::NO,
                               EvaluateEnergy::YES, AccumulationMethod::SPLIT, is_models[j],
                               clash_policy[i],
                               "ktgf" + is_model_names[j] + "Energy" + clash_ext[i]);
        catalogNonbondedKernel(PrecisionModel::SINGLE, NbwuKind::TILE_GROUPS, EvaluateForce::YES,
                               EvaluateEnergy::NO, AccumulationMethod::SPLIT, is_models[j],
                               clash_policy[i],
                               "ktgfs" + is_model_names[j] + "Force" + clash_ext[i]);
        catalogNonbondedKernel(PrecisionModel::SINGLE, NbwuKind::TILE_GROUPS, EvaluateForce::YES,
                               EvaluateEnergy::NO, AccumulationMethod::WHOLE, is_models[j],
                               clash_policy[i],
                               "ktgfs" + is_model_names[j] + "Force" + clash_ext[i]);
        catalogNonbondedKernel(PrecisionModel::SINGLE, NbwuKind::TILE_GROUPS, EvaluateForce::YES,
                               EvaluateEnergy::YES, AccumulationMethod::SPLIT, is_models[j],
                               clash_policy[i],
                               "ktgfs" + is_model_names[j] + "ForceEnergy" + clash_ext[i]);
        catalogNonbondedKernel(PrecisionModel::SINGLE, NbwuKind::TILE_GROUPS, EvaluateForce::YES,
                               EvaluateEnergy::YES, AccumulationMethod::WHOLE, is_models[j],
                               clash_policy[i],
                               "ktgfs" + is_model_names[j] + "ForceEnergy" + clash_ext[i]);

        // Generalized Born radii and radial derivative kernel entries
        if (is_models[j] != ImplicitSolventModel::NONE && clash_policy[i] == ClashResponse::NONE) {
          catalogBornRadiiKernel(PrecisionModel::DOUBLE, NbwuKind::TILE_GROUPS,
                                 AccumulationMethod::SPLIT, is_models[j],
                                 "ktgdsCalculate" + is_model_names[j] + "Radii");
          catalogBornRadiiKernel(PrecisionModel::SINGLE, NbwuKind::TILE_GROUPS,
                                 AccumulationMethod::SPLIT, is_models[j],
                                 "ktgfsCalculate" + is_model_names[j] + "Radii");
          catalogBornRadiiKernel(PrecisionModel::SINGLE, NbwuKind::TILE_GROUPS,
                                 AccumulationMethod::WHOLE, is_models[j],
                                 "ktgfCalculate" + is_model_names[j] + "Radii");
          catalogBornDerivativeKernel(PrecisionModel::DOUBLE, NbwuKind::TILE_GROUPS,
                                      AccumulationMethod::SPLIT, is_models[j],
                                      "ktgdsCalculate" + is_model_names[j] + "Derivatives");
          catalogBornDerivativeKernel(PrecisionModel::SINGLE, NbwuKind::TILE_GROUPS,
                                      AccumulationMethod::SPLIT, is_models[j],
                                      "ktgfsCalculate" + is_model_names[j] + "Derivatives");
          catalogBornDerivativeKernel(PrecisionModel::SINGLE, NbwuKind::TILE_GROUPS,
                                      AccumulationMethod::WHOLE, is_models[j],
                                      "ktgfCalculate" + is_model_names[j] + "Derivatives");
        }
        break;
      case UnitCellType::ORTHORHOMBIC:
      case UnitCellType::TRICLINIC:
        break;
      }
    }
  }

  // Stand-alone integration kernel entries
  const std::vector<IntegrationStage> time_step_parts = { IntegrationStage::VELOCITY_ADVANCE,
                                                          IntegrationStage::VELOCITY_CONSTRAINT,
                                                          IntegrationStage::POSITION_ADVANCE,
                                                          IntegrationStage::GEOMETRY_CONSTRAINT };
  const std::vector<std::string> time_step_abbrevs = { "VelAdv", "VelCnst", "PosAdv", "GeomCnst" };
  for (size_t i = 0; i < time_step_parts.size(); i++) {
    catalogIntegrationKernel(PrecisionModel::DOUBLE, AccumulationMethod::SPLIT,
                             valence_kernel_width, time_step_parts[i],
                             "kdsIntegration" + time_step_abbrevs[i] + valence_kwidth_ext);
    catalogIntegrationKernel(PrecisionModel::SINGLE, AccumulationMethod::SPLIT,
                             valence_kernel_width, time_step_parts[i],
                             "kfsIntegration" + time_step_abbrevs[i] + valence_kwidth_ext);
    catalogIntegrationKernel(PrecisionModel::SINGLE, AccumulationMethod::WHOLE,
                             valence_kernel_width, time_step_parts[i],
                             "kfIntegration" + time_step_abbrevs[i] + valence_kwidth_ext);
  }
  
  // PME density mapping (spreading) kernel entries
  const std::vector<PrecisionModel> all_prec = { PrecisionModel::DOUBLE, PrecisionModel::SINGLE };
  const std::vector<bool> use_short_format = { true, false };
  const std::vector<std::string> short_formats = { "", "sf" };
  for (int order = 4; order <= 6; order++) {
    for (size_t i = 0; i < all_prec.size(); i++) {
      std::string prec_ordr;
      switch (all_prec[i]) {
      case PrecisionModel::DOUBLE:
        prec_ordr = "d";
        break;
      case PrecisionModel::SINGLE:
        prec_ordr = "f";
        break;
      }
      prec_ordr += std::to_string(order);
      catalogGeneralQMapKernel(all_prec[i], int_type_index, order,
                               "ksi" + prec_ordr + "MapDensity");
      catalogGeneralQMapKernel(all_prec[i], llint_type_index, order,
                               "kli" + prec_ordr + "MapDensity");
      catalogGeneralQMapKernel(all_prec[i], float_type_index, order,
                               "ksr" + prec_ordr + "MapDensity");
      catalogGeneralQMapKernel(all_prec[i], double_type_index, order,
                               "klr" + prec_ordr + "MapDensity");
      for (size_t j = 0; j < all_prec.size(); j++) {
        std::string aprec_ordr = prec_ordr;
        switch (all_prec[j]) {
        case PrecisionModel::DOUBLE:
          aprec_ordr += "d";
          break;
        case PrecisionModel::SINGLE:
          aprec_ordr += "s";
          break;
        }
        for (int k = 0; k < 2; k++) {
          std::string bprec_ordr = aprec_ordr + short_formats[k];
          catalogShrAccQMapKernel(all_prec[i], all_prec[j], use_short_format[k], int_type_index,
                                  order, "kSAsi" + bprec_ordr + "MapDensity");
          catalogShrAccQMapKernel(all_prec[i], all_prec[j], use_short_format[k], llint_type_index,
                                  order, "kSAli" + bprec_ordr + "MapDensity");
          catalogShrAccQMapKernel(all_prec[i], all_prec[j], use_short_format[k], float_type_index,
                                  order, "kSAsr" + bprec_ordr + "MapDensity");
          catalogShrAccQMapKernel(all_prec[i], all_prec[j], use_short_format[k], double_type_index,
                                  order, "kSAlr" + bprec_ordr + "MapDensity");
        }
      }
    }
  }

  // Reduction kernel entries
  const int reduction_div = optReductionKernelSubdivision(poly_ag->getSystemAtomCounts(), gpu_in);
  catalogReductionKernel(PrecisionModel::DOUBLE, ReductionGoal::CONJUGATE_GRADIENT,
                         ReductionStage::GATHER, reduction_div, "kdgtConjGrad");
  catalogReductionKernel(PrecisionModel::DOUBLE, ReductionGoal::CONJUGATE_GRADIENT,
                         ReductionStage::SCATTER, reduction_div, "kdscConjGrad");
  catalogReductionKernel(PrecisionModel::DOUBLE, ReductionGoal::CONJUGATE_GRADIENT,
                         ReductionStage::RESCALE, reduction_div, "kdrsConjGrad");
  catalogReductionKernel(PrecisionModel::DOUBLE, ReductionGoal::CONJUGATE_GRADIENT,
                         ReductionStage::ALL_REDUCE, reduction_div, "kdrdConjGrad");
  catalogReductionKernel(PrecisionModel::SINGLE, ReductionGoal::CONJUGATE_GRADIENT,
                         ReductionStage::GATHER, reduction_div, "kfgtConjGrad");
  catalogReductionKernel(PrecisionModel::SINGLE, ReductionGoal::CONJUGATE_GRADIENT,
                         ReductionStage::SCATTER, reduction_div, "kfscConjGrad");
  catalogReductionKernel(PrecisionModel::SINGLE, ReductionGoal::CONJUGATE_GRADIENT,
                         ReductionStage::RESCALE, reduction_div, "kfrsConjGrad");
  catalogReductionKernel(PrecisionModel::SINGLE, ReductionGoal::CONJUGATE_GRADIENT,
                         ReductionStage::ALL_REDUCE, reduction_div, "kfrdConjGrad");

  // Virtual site kernel entries
  const std::string vs_kwidth_ext = getEnumerationName(virtual_site_kernel_width);
  catalogVirtualSiteKernel(PrecisionModel::DOUBLE, VirtualSiteActivity::PLACEMENT,
                           virtual_site_kernel_width, "kdPlaceVirtualSites" + vs_kwidth_ext);
  catalogVirtualSiteKernel(PrecisionModel::SINGLE, VirtualSiteActivity::PLACEMENT,
                           virtual_site_kernel_width, "kfPlaceVirtualSites" + vs_kwidth_ext);
  catalogVirtualSiteKernel(PrecisionModel::DOUBLE, VirtualSiteActivity::TRANSMIT_FORCES,
                           virtual_site_kernel_width, "kdTransmitVSiteForces" + vs_kwidth_ext);
  catalogVirtualSiteKernel(PrecisionModel::SINGLE, VirtualSiteActivity::TRANSMIT_FORCES,
                           virtual_site_kernel_width, "kfTransmitVSiteForces" + vs_kwidth_ext);

  // RMSD kernel entries
  catalogRMSDKernel(PrecisionModel::DOUBLE, RMSDTask::REFERENCE, "kdComputeRMSDToReference");
  catalogRMSDKernel(PrecisionModel::SINGLE, RMSDTask::REFERENCE, "kfComputeRMSDToReference");
  catalogRMSDKernel(PrecisionModel::DOUBLE, RMSDTask::MATRIX, "kdComputeRMSDMatrix");
  catalogRMSDKernel(PrecisionModel::SINGLE, RMSDTask::MATRIX, "kfComputeRMSDMatrix");

  // PME pair interactions
  const std::vector<NeighborListKind> all_neighbor_lists = { NeighborListKind::DUAL,
                                                             NeighborListKind::MONO };
  const std::vector<TinyBoxPresence> all_box_wrappings = { TinyBoxPresence::YES,
                                                           TinyBoxPresence::NO };
  std::string kname("kxx");
  kname.reserve(96);
  for (size_t i = 0; i < all_prec.size(); i++) {
    switch (all_prec[i]) {
    case PrecisionModel::DOUBLE:
      kname[1] = 'd';
      break;
    case PrecisionModel::SINGLE:
      kname[1] = 'f';
      break;
    }
    for (size_t j = 0; j < all_prec.size(); j++) {
      switch (all_prec[j]) {
      case PrecisionModel::DOUBLE:
        kname[2] = 'd';
        break;
      case PrecisionModel::SINGLE:
        kname[2] = 'f';
        break;
      }
      for (size_t k = 0; k < all_neighbor_lists.size(); k++) {
        std::string kextn;
        switch (all_neighbor_lists[k]) {
        case NeighborListKind::MONO:
          break;
        case NeighborListKind::DUAL:
          kextn += "Dual";
          break;
        }
        const int kextn_ngbr_len = kextn.size();
        for (size_t m = 0; m < all_box_wrappings.size(); m++) {
          kextn.resize(kextn_ngbr_len);
          switch (all_box_wrappings[m]) {
          case TinyBoxPresence::NO:
            break;
          case TinyBoxPresence::YES:
            kextn += "Tiny";
            break;
          }
          const int kextn_wrap_len = kextn.size();
          for (size_t n = 0; n < clash_policy.size(); n++) {
            kextn.resize(kextn_wrap_len);
            switch (clash_policy[n]) {
            case ClashResponse::NONE:
              break;
            case ClashResponse::FORGIVE:
              kextn += "NonClash";
              break;
            }
            const int kextn_bump_len = kextn.size();
            kname.resize(3);
            kname += "TowerPlate";
            kname += "FE" + kextn;
            catalogPMEPairsKernel(all_prec[i], all_prec[j], all_neighbor_lists[k],
                                  all_box_wrappings[m], EvaluateForce::YES, EvaluateEnergy::YES,
                                  clash_policy[n], kname);
            kname.resize(13);
            kname += "FX" + kextn;
            catalogPMEPairsKernel(all_prec[i], all_prec[j], all_neighbor_lists[k],
                                  all_box_wrappings[m], EvaluateForce::YES, EvaluateEnergy::NO,
                                  clash_policy[n], kname);
            kname.resize(13);
            kname += "XE" + kextn;
            catalogPMEPairsKernel(all_prec[i], all_prec[j], all_neighbor_lists[k],
                                  all_box_wrappings[m], EvaluateForce::NO, EvaluateEnergy::YES,
                                  clash_policy[n], kname);
          }
        }
      }
    }
  }

  // Particle migration kernels
  for (size_t i = 0; i < all_prec.size(); i++) {
    for (size_t j = 0; j < all_neighbor_lists.size(); j++) {
      for (int k = 1; k <= 2; k++) {
        std::string kname = "k";
        switch (all_prec[i]) {
        case PrecisionModel::DOUBLE:
          kname += "d";
          break;
        case PrecisionModel::SINGLE:
          kname += "f";
          break;
        }
        kname += "Migration";
        kname += (k == 1) ? "One" : "Two";
        switch (all_neighbor_lists[j]) {
        case NeighborListKind::DUAL:
          kname += "Dual";
          break;
        case NeighborListKind::MONO:
          break;
        }
        catalogMigrationKernel(all_prec[i], all_neighbor_lists[j], k,
                               globalpos_scale_nonoverflow_bits - 1, kname);
        if (k == 1) {
          catalogMigrationKernel(all_prec[i], all_neighbor_lists[j], k,
                                 globalpos_scale_nonoverflow_bits + 1, kname);
        }
      }
    }
  }
#endif
}

//-------------------------------------------------------------------------------------------------
CoreKlManager::CoreKlManager(const GpuDetails &gpu_in, const AtomGraphSynthesis &poly_ag) :
    CoreKlManager(gpu_in, poly_ag.getSelfPointer())
{}

//-------------------------------------------------------------------------------------------------
void CoreKlManager::catalogValenceKernel(const PrecisionModel prec, const EvaluateForce eval_force,
                                         const EvaluateEnergy eval_nrg,
                                         const AccumulationMethod acc_meth,
                                         const VwuGoal purpose,
                                         const ClashResponse collision_handling,
                                         const ValenceKernelSize kwidth,
                                         const std::string &kernel_name) {
  const std::string k_key = valenceKernelKey(prec, eval_force, eval_nrg, acc_meth, purpose,
                                             collision_handling, kwidth);
  std::map<std::string, KernelFormat>::iterator it = k_dictionary.find(k_key);
  if (it != k_dictionary.end()) {
    rtErr("Valence kernel identifier " + k_key + " already exists in the kernel map.",
          "CoreKlManager", "catalogValenceKernel");
  }
#ifdef STORMM_USE_HPC
#  ifdef STORMM_USE_CUDA
  const cudaFuncAttributes attr = queryValenceKernelRequirements(prec, eval_force, eval_nrg,
                                                                 acc_meth, purpose,
                                                                 collision_handling, kwidth);

  // Infer the maximum number of blocks per streaming multiprocessor based on the kernel width and
  // the thread count of the largest possible kernel.
  int block_mult;
  switch (prec) {
  case PrecisionModel::DOUBLE:
    block_mult = 1;
    break;
  case PrecisionModel::SINGLE:
    {
      const cudaFuncAttributes attr_xl = queryValenceKernelRequirements(prec, eval_force, eval_nrg,
                                                                        acc_meth, purpose,
                                                                        collision_handling,
                                                                        ValenceKernelSize::XL);
      block_mult = attr_xl.maxThreadsPerBlock / attr.maxThreadsPerBlock;
    }
    break;
  }
  k_dictionary[k_key] = KernelFormat(attr, 2 * block_mult, 1, gpu, kernel_name);
#  endif
#else
  k_dictionary[k_key] = KernelFormat();
#endif
}

//-------------------------------------------------------------------------------------------------
void CoreKlManager::catalogValenceKernel(const PrecisionModel prec,
                                         const PrecisionModel neighbor_prec,
                                         const NeighborListKind neighbor_layout,
                                         const EvaluateEnergy eval_nrg,
                                         const AccumulationMethod acc_meth,
                                         const ClashResponse collision_handling,
                                         const std::string &kernel_name) {
  const std::string k_key = valenceKernelKey(prec, neighbor_prec, neighbor_layout, eval_nrg,
                                             acc_meth, collision_handling);
  std::map<std::string, KernelFormat>::iterator it = k_dictionary.find(k_key);
  if (it != k_dictionary.end()) {
    rtErr("Valence kernel identifier " + k_key + " already exists in the kernel map.",
          "CoreKlManager", "catalogValenceKernel");
  }
#ifdef STORMM_USE_HPC
#  ifdef STORMM_USE_CUDA
  const cudaFuncAttributes attr = queryValenceKernelRequirements(prec, EvaluateForce::YES,
                                                                 eval_nrg, acc_meth,
                                                                 VwuGoal::MOVE_PARTICLES,
                                                                 collision_handling,
                                                                 ValenceKernelSize::XL, true,
                                                                 neighbor_prec, neighbor_layout);
  k_dictionary[k_key] = KernelFormat(attr, 2, 1, gpu, kernel_name);
#  endif
#else
  k_dictionary[k_key] = KernelFormat();
#endif
}

//-------------------------------------------------------------------------------------------------
void CoreKlManager::catalogIntegrationKernel(const PrecisionModel prec,
                                             const AccumulationMethod acc_meth,
                                             const ValenceKernelSize kwidth,
                                             const IntegrationStage process,
                                             const std::string &kernel_name) {
  const std::string k_key = integrationKernelKey(prec, acc_meth, kwidth, process);
  std::map<std::string, KernelFormat>::iterator it = k_dictionary.find(k_key);
  if (it != k_dictionary.end()) {
    rtErr("Integration kernel identifier " + k_key + " already exists in the kernel map.",
          "CoreKlManager", "catalogIntegrationKernel");
  }
#ifdef STORMM_USE_HPC
#  ifdef STORMM_USE_CUDA
  const cudaFuncAttributes attr = queryIntegrationKernelRequirements(prec, acc_meth, kwidth,
                                                                     process);

  // As above with the valence kernels, infer the block count based on the thread count of the
  // largest possible kernel.
  int block_mult;
  switch (prec) {
  case PrecisionModel::DOUBLE:
    block_mult = 1;
    break;
  case PrecisionModel::SINGLE:
    {
      const cudaFuncAttributes attr_xl = queryIntegrationKernelRequirements(prec, acc_meth,
                                                                            ValenceKernelSize::XL,
                                                                            process);
      block_mult = attr_xl.maxThreadsPerBlock / attr.maxThreadsPerBlock;
    }
    break;
  }
  k_dictionary[k_key] = KernelFormat(attr, 2 * block_mult, 1, gpu, kernel_name);
#  endif
#else
  k_dictionary[k_key] = KernelFormat();
#endif
}

//-------------------------------------------------------------------------------------------------
void CoreKlManager::catalogNonbondedKernel(const PrecisionModel prec, const NbwuKind kind,
                                           const EvaluateForce eval_force,
                                           const EvaluateEnergy eval_nrg,
                                           const AccumulationMethod acc_meth,
                                           const ImplicitSolventModel igb,
                                           const ClashResponse collision_handling,
                                           const std::string &kernel_name) {
  const std::string k_key = nonbondedKernelKey(prec, kind, eval_force, eval_nrg, acc_meth, igb,
                                               collision_handling);
  std::map<std::string, KernelFormat>::iterator it = k_dictionary.find(k_key);
  if (it != k_dictionary.end()) {
    rtErr("Non-bonded kernel identifier " + k_key + " already exists in the kernel map.",
          "CoreKlManager", "catalogNonbondedKernel");
  }
#ifdef STORMM_USE_HPC
#  ifdef STORMM_USE_CUDA
  const cudaFuncAttributes attr = queryNonbondedKernelRequirements(prec, kind, eval_force,
                                                                   eval_nrg, acc_meth, igb,
                                                                   collision_handling);
  int nonbond_block_multiplier;
  switch (prec) {
  case PrecisionModel::DOUBLE:
    nonbond_block_multiplier = nonbond_block_multiplier_dp;
     break;
  case PrecisionModel::SINGLE:
    nonbond_block_multiplier = nonbond_block_multiplier_sp;
    break;
  }
  k_dictionary[k_key] = KernelFormat(attr, nonbond_block_multiplier, 1, gpu, kernel_name);
#  endif
#else
  k_dictionary[k_key] = KernelFormat();
#endif
}

//-------------------------------------------------------------------------------------------------
void CoreKlManager::catalogBornRadiiKernel(const PrecisionModel prec, const NbwuKind kind,
                                           const AccumulationMethod acc_meth,
                                           const ImplicitSolventModel igb,
                                           const std::string &kernel_name) {
  const std::string k_key = bornRadiiKernelKey(prec, kind, acc_meth, igb);
  std::map<std::string, KernelFormat>::iterator it = k_dictionary.find(k_key);
  if (it != k_dictionary.end()) {
    rtErr("Born radii kernel identifier " + k_key + " already exists in the kernel map.",
          "CoreKlManager", "catalogBornRadiiKernel");
  }
#ifdef STORMM_USE_HPC
#  ifdef STORMM_USE_CUDA
  const cudaFuncAttributes attr = queryBornRadiiKernelRequirements(prec, kind, acc_meth, igb);
  int gbradii_block_multiplier;
  switch (prec) {
  case PrecisionModel::DOUBLE:
    gbradii_block_multiplier = gbradii_block_multiplier_dp;
    break;
  case PrecisionModel::SINGLE:
    gbradii_block_multiplier = gbradii_block_multiplier_sp;
    break;
  }
  k_dictionary[k_key] = KernelFormat(attr, gbradii_block_multiplier, 1, gpu, kernel_name);
#  endif
#else
  k_dictionary[k_key] = KernelFormat();
#endif
}

//-------------------------------------------------------------------------------------------------
void CoreKlManager::catalogBornDerivativeKernel(const PrecisionModel prec, const NbwuKind kind,
                                                const AccumulationMethod acc_meth,
                                                const ImplicitSolventModel igb,
                                                const std::string &kernel_name) {
  const std::string k_key = bornDerivativeKernelKey(prec, kind, acc_meth, igb);
  std::map<std::string, KernelFormat>::iterator it = k_dictionary.find(k_key);
  if (it != k_dictionary.end()) {
    rtErr("Born radii derivative kernel identifier " + k_key + " already exists in the kernel "
          "map.", "CoreKlManager", "catalogBornRadiiKernel");
  }
#ifdef STORMM_USE_HPC
#  ifdef STORMM_USE_CUDA
  const cudaFuncAttributes attr = queryBornDerivativeKernelRequirements(prec, kind, acc_meth, igb);
  int gbderiv_block_multiplier;
  switch (prec) {
  case PrecisionModel::DOUBLE:
    gbderiv_block_multiplier = gbderiv_block_multiplier_dp;
    break;
  case PrecisionModel::SINGLE:
    gbderiv_block_multiplier = gbderiv_block_multiplier_sp;
    break;
  }
  k_dictionary[k_key] = KernelFormat(attr, gbderiv_block_multiplier, 1, gpu, kernel_name);
#  endif
#else
  k_dictionary[k_key] = KernelFormat();
#endif
}

//-------------------------------------------------------------------------------------------------
void CoreKlManager::catalogGeneralQMapKernel(const PrecisionModel prec, const size_t cg_tmat,
                                             const int order, const std::string &kernel_name) {
  const std::string k_key = generalQMapKernelKey(prec, cg_tmat, order);
  std::map<std::string, KernelFormat>::iterator it = k_dictionary.find(k_key);
  if (it != k_dictionary.end()) {
    rtErr("Particle density mapping kernel identifier " + k_key + " already exists in the kernel "
          "map.", "CoreKlManager", "catalogGeneralQMapKernel");
  }
#ifdef STORMM_USE_HPC
#  ifdef STORMM_USE_CUDA
  const cudaFuncAttributes attr = queryGeneralQMapKernelRequirements(prec, cg_tmat, order);
  int gen_qmap_block_multiplier;
  switch (prec) {
  case PrecisionModel::DOUBLE:
    gen_qmap_block_multiplier = gen_qmap_block_multiplier_dp[order];
    break;
  case PrecisionModel::SINGLE:
    gen_qmap_block_multiplier = gen_qmap_block_multiplier_sp[order];
    break;
  }
  k_dictionary[k_key] = KernelFormat(attr, gen_qmap_block_multiplier, 1, gpu, kernel_name);
#  endif
#else
  k_dictionary[k_key] = KernelFormat();
#endif
}

//-------------------------------------------------------------------------------------------------
void CoreKlManager::catalogShrAccQMapKernel(const PrecisionModel calc_prec,
                                            const PrecisionModel acc_prec,
                                            const bool overflow, const size_t cg_tmat,
                                            const int order, const std::string &kernel_name) {
  const std::string k_key = shrAccQMapKernelKey(calc_prec, acc_prec, overflow, cg_tmat, order);
  std::map<std::string, KernelFormat>::iterator it = k_dictionary.find(k_key);
  if (it != k_dictionary.end()) {
    rtErr("Particle density mapping kernel identifier " + k_key + " already exists in the kernel "
          "map.", "CoreKlManager", "catalogShrAccQMapKernel");
  }
#ifdef STORMM_USE_HPC
#  ifdef STORMM_USE_CUDA
  const cudaFuncAttributes attr = queryShrAccQMapKernelRequirements(calc_prec, acc_prec, overflow,
                                                                    cg_tmat, order);
  int sac_qmap_block_multiplier;
  switch (calc_prec) {
  case PrecisionModel::DOUBLE:
    sac_qmap_block_multiplier = sac_qmap_block_multiplier_dp[order];
    break;
  case PrecisionModel::SINGLE:
    sac_qmap_block_multiplier = sac_qmap_block_multiplier_sp[order];
    break;
  }
  k_dictionary[k_key] = KernelFormat(attr, sac_qmap_block_multiplier, 1, gpu, kernel_name);
#  endif
#else
  k_dictionary[k_key] = KernelFormat();
#endif
}

//-------------------------------------------------------------------------------------------------
void CoreKlManager::catalogReductionKernel(const PrecisionModel prec, const ReductionGoal purpose,
                                           const ReductionStage process, const int subdivision,
                                           const std::string &kernel_name) {
  const std::string k_key = reductionKernelKey(prec, purpose, process);
  std::map<std::string, KernelFormat>::iterator it = k_dictionary.find(k_key);
  if (it != k_dictionary.end()) {
    rtErr("Reduction kernel identifier " + k_key + " already exists in the kernel map.",
          "CoreKlManager", "catalogReductionKernel");
  }
#ifdef STORMM_USE_HPC
#  ifdef STORMM_USE_CUDA
  const cudaFuncAttributes attr = queryReductionKernelRequirements(prec, purpose, process);
  k_dictionary[k_key] = KernelFormat(attr, reduction_block_multiplier, subdivision, gpu,
                                     kernel_name);
#  endif
#else
  k_dictionary[k_key] = KernelFormat();
#endif
}

//-------------------------------------------------------------------------------------------------
void CoreKlManager::catalogVirtualSiteKernel(const PrecisionModel prec,
                                             const VirtualSiteActivity purpose,
                                             const ValenceKernelSize kwidth,
                                             const std::string &kernel_name) {
  const std::string k_key = virtualSiteKernelKey(prec, purpose, kwidth);
  std::map<std::string, KernelFormat>::iterator it = k_dictionary.find(k_key);
  if (it != k_dictionary.end()) {
    rtErr("Virtual site handling kernel identifier " + k_key + " already exists in the kernel "
          "map.", "CoreKlManager", "catalogVirtualSiteKernel");
  }
#ifdef STORMM_USE_HPC
#  ifdef STORMM_USE_CUDA
  const cudaFuncAttributes attr = queryVirtualSiteKernelRequirements(prec, purpose, kwidth);
  int virtual_site_block_multiplier;
  switch (prec) {
  case PrecisionModel::DOUBLE:
    virtual_site_block_multiplier = 1;
    break;
  case PrecisionModel::SINGLE:
    {
      const cudaFuncAttributes attr_xl = queryVirtualSiteKernelRequirements(prec, purpose,
                                                                            ValenceKernelSize::XL);
      virtual_site_block_multiplier = attr_xl.maxThreadsPerBlock / attr.maxThreadsPerBlock;
    }
    break;
  }
  k_dictionary[k_key] = KernelFormat(attr, 2 * virtual_site_block_multiplier, 1, gpu, kernel_name);
#  endif
#else
  k_dictionary[k_key] = KernelFormat();
#endif
}

//-------------------------------------------------------------------------------------------------
void CoreKlManager::catalogRMSDKernel(const PrecisionModel prec, const RMSDTask order,
                                      const std::string &kernel_name) {
  const std::string k_key = rmsdKernelKey(prec, order);
  std::map<std::string, KernelFormat>::iterator it = k_dictionary.find(k_key);
  if (it != k_dictionary.end()) {
    rtErr("RMSD calculation kernel identifier " + k_key + " already exists in the kernel map.",
          "CoreKlManager", "catalogRMSDKernel");
  }
#ifdef STORMM_USE_HPC
#  ifdef STORMM_USE_CUDA
  const cudaFuncAttributes attr = queryRMSDKernelRequirements(prec, order);
  int rmsd_block_multiplier;
  switch (prec) {
  case PrecisionModel::DOUBLE:
    rmsd_block_multiplier = rmsd_block_multiplier_dp;
    break;
  case PrecisionModel::SINGLE:
    rmsd_block_multiplier = rmsd_block_multiplier_sp;
    break;
  }
  k_dictionary[k_key] = KernelFormat(attr, rmsd_block_multiplier, 1, gpu, kernel_name);
#  endif
#else
  k_dictionary[k_key] = KernelFormat();
#endif
}

//-------------------------------------------------------------------------------------------------
void CoreKlManager::catalogPMEPairsKernel(const PrecisionModel coord_prec,
                                          const PrecisionModel calc_prec,
                                          const NeighborListKind grid_configuration,
                                          const TinyBoxPresence has_tiny_box,
                                          const EvaluateForce eval_frc,
                                          const EvaluateEnergy eval_nrg,
                                          const ClashResponse mitigation,
                                          const std::string &kernel_name) {
  const std::string k_key = pmePairsKernelKey(coord_prec, calc_prec, grid_configuration,
                                              has_tiny_box, eval_frc, eval_nrg, mitigation);
  std::map<std::string, KernelFormat>::iterator it = k_dictionary.find(k_key);
  if (it != k_dictionary.end()) {
    rtErr("PME particle-particle pair interactions kernel identifier " + k_key + " already exists "
          "in the kernel map.", "CoreKlManager", "catalogPMEPairsKernel");
  }
#ifdef STORMM_USE_HPC
#  ifdef STORMM_USE_CUDA
  const cudaFuncAttributes attr = queryPMEPairsKernelRequirements(coord_prec, calc_prec,
                                                                  grid_configuration, eval_frc,
                                                                  eval_nrg, has_tiny_box,
                                                                  mitigation);
  k_dictionary[k_key] = KernelFormat(attr, pmePairsBlockMultiplier(gpu, coord_prec, calc_prec), 1,
                                     gpu, kernel_name);
#  endif
#else
  k_dictionary[k_key] = KernelFormat();
#endif
}

//-------------------------------------------------------------------------------------------------
void CoreKlManager::catalogMigrationKernel(const PrecisionModel coord_prec,
                                           const NeighborListKind grid_configuration,
                                           const int stage_number, const int gpos_bits,
                                           const std::string &kernel_name) {
  const std::string k_key = migrationKernelKey(coord_prec, grid_configuration, stage_number,
                                               gpos_bits);
  std::map<std::string, KernelFormat>::iterator it = k_dictionary.find(k_key);
  if (it != k_dictionary.end()) {
    rtErr("PME particle migration kernel identifier " + k_key + " already exists in the kernel "
          "map.", "CoreKlManager", "catalogMigrationKernel");
  }
#ifdef STORMM_USE_HPC
#  ifdef STORMM_USE_CUDA
  const cudaFuncAttributes attr = queryMigrationKernelRequirements(coord_prec,
                                                                   grid_configuration,
                                                                   stage_number, gpos_bits);
  k_dictionary[k_key] = KernelFormat(attr, migrationBlockMultiplier(), 1, gpu, kernel_name);
#  endif
#else
  k_dictionary[k_key] = KernelFormat();
#endif
}

//-------------------------------------------------------------------------------------------------
int2 CoreKlManager::getValenceKernelDims(const PrecisionModel prec, const EvaluateForce eval_force,
                                         const EvaluateEnergy eval_nrg,
                                         const AccumulationMethod acc_meth,
                                         const VwuGoal purpose,
                                         const ClashResponse collision_handling) const {
  const std::string k_key = valenceKernelKey(prec, eval_force, eval_nrg, acc_meth, purpose,
                                             collision_handling, valence_kernel_width);
  if (k_dictionary.find(k_key) == k_dictionary.end()) {
    rtErr("Valence kernel identifier " + k_key + " was not found in the kernel map.",
          "CoreKlManager", "getValenceKernelDims");
  }
  return k_dictionary.at(k_key).getLaunchParameters();
}

//-------------------------------------------------------------------------------------------------
int2 CoreKlManager::getPmeValenceKernelDims(const PrecisionModel prec,
                                            const PrecisionModel neighbor_prec,
                                            const NeighborListKind neighbor_layout,
                                            const EvaluateEnergy eval_nrg,
                                            const AccumulationMethod acc_meth,
                                            const ClashResponse collision_handling) const {
  const std::string k_key = valenceKernelKey(prec, EvaluateForce::YES, eval_nrg, acc_meth,
                                             VwuGoal::MOVE_PARTICLES, collision_handling,
                                             ValenceKernelSize::XL);
  if (k_dictionary.find(k_key) == k_dictionary.end()) {
    rtErr("Valence kernel identifier " + k_key + " was not found in the kernel map.",
          "CoreKlManager", "getValenceKernelDims");
  }
  return k_dictionary.at(k_key).getLaunchParameters();
}

//-------------------------------------------------------------------------------------------------
int2 CoreKlManager::getIntegrationKernelDims(const PrecisionModel prec,
                                             const AccumulationMethod acc_meth,
                                             const IntegrationStage process) const {
  const std::string k_key = integrationKernelKey(prec, acc_meth, valence_kernel_width, process);
  if (k_dictionary.find(k_key) == k_dictionary.end()) {
    rtErr("Integration kernel identifier " + k_key + " was not found in the kernel map.",
          "CoreKlManager", "getIntegrationKernelDims");
  }
  return k_dictionary.at(k_key).getLaunchParameters();
}

//-------------------------------------------------------------------------------------------------
int2 CoreKlManager::getNonbondedKernelDims(const PrecisionModel prec, const NbwuKind kind,
                                           const EvaluateForce eval_force,
                                           const EvaluateEnergy eval_nrg,
                                           const AccumulationMethod acc_meth,
                                           const ImplicitSolventModel igb,
                                           const ClashResponse collision_handling) const {
  const std::string k_key = nonbondedKernelKey(prec, kind, eval_force, eval_nrg, acc_meth, igb,
                                               collision_handling);
  if (k_dictionary.find(k_key) == k_dictionary.end()) {
    rtErr("Non-bonded kernel identifier " + k_key + " was not found in the kernel map.",
          "CoreKlManager", "getNonbondedKernelDims");
  }
  return k_dictionary.at(k_key).getLaunchParameters();
}

//-------------------------------------------------------------------------------------------------
int2 CoreKlManager::getBornRadiiKernelDims(const PrecisionModel prec, const NbwuKind kind,
                                           const AccumulationMethod acc_meth,
                                           const ImplicitSolventModel igb) const {
  if (igb == ImplicitSolventModel::NONE) {
    return { 0, 0 };
  }
  const std::string k_key = bornRadiiKernelKey(prec, kind, acc_meth, igb);
  if (k_dictionary.find(k_key) == k_dictionary.end()) {
    rtErr("Born radii computation kernel identifier " + k_key + " was not found in the kernel "
          "map.", "CoreKlManager", "getBornRadiiKernelDims");
  }
  return k_dictionary.at(k_key).getLaunchParameters();
}

//-------------------------------------------------------------------------------------------------
int2 CoreKlManager::getBornDerivativeKernelDims(const PrecisionModel prec, const NbwuKind kind,
                                                const AccumulationMethod acc_meth,
                                                const ImplicitSolventModel igb) const {
  if (igb == ImplicitSolventModel::NONE) {
    return { 0, 0 };
  }
  const std::string k_key = bornDerivativeKernelKey(prec, kind, acc_meth, igb);
  if (k_dictionary.find(k_key) == k_dictionary.end()) {
    rtErr("Born radii derivative computation kernel identifier " + k_key + " was not found in "
          "the kernel map.", "CoreKlManager", "getBornRadiiKernelDims");
  }
  return k_dictionary.at(k_key).getLaunchParameters();
}

//-------------------------------------------------------------------------------------------------
int2 CoreKlManager::getDensityMappingKernelDims(const QMapMethod approach,
                                                const PrecisionModel calc_prec,
                                                const PrecisionModel acc_prec, const bool overflow,
                                                const size_t cg_tmat, const int order) const {
  std::string k_key;
  switch (approach) {
  case QMapMethod::ACC_SHARED:
    k_key = shrAccQMapKernelKey(calc_prec, acc_prec, overflow, cg_tmat, order);
    break;
  case QMapMethod::GENERAL_PURPOSE:
    k_key = generalQMapKernelKey(calc_prec, cg_tmat, order);
    break;
  case QMapMethod::AUTOMATIC:
    rtErr("A density mapping strategy must be decided before requesting kernel launch parameters.",
          "CoreKlManager", "getDensityMappingKernelDims");
  }
  if (k_dictionary.find(k_key) == k_dictionary.end()) {
    rtErr("Particle density mapping kernel identifier " + k_key + " was not found in "
          "the kernel map.", "CoreKlManager", "getDensityMappingKernelDims");
  }
  return k_dictionary.at(k_key).getLaunchParameters();
}

//-------------------------------------------------------------------------------------------------
int2 CoreKlManager::getDensityMappingKernelDims(const QMapMethod approach,
                                                const PrecisionModel prec,
                                                const size_t cg_tmat, const int order) const {
  switch (approach) {
  case QMapMethod::GENERAL_PURPOSE:
  case QMapMethod::AUTOMATIC:

    // The approach must be GENERAL_PURPOSE, or the overloaded form that is called will raise an
    // exception. Mark the need for overflow bits as TRUE.  This is not a factor in selecting
    // kernels for the naive mapping approach, although the kernels do make a decision internally
    // based on the bit count.
    return getDensityMappingKernelDims(approach, prec, prec, true, cg_tmat, order);
  case QMapMethod::ACC_SHARED:
    rtErr("Both the calculation and accumulation precisions are required inputs for kernels that "
          "carry out " + getEnumerationName(approach) + ".", "CoreKlManager",
          "getDensityMappingKernelDims");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
int2 CoreKlManager::getReductionKernelDims(const PrecisionModel prec, const ReductionGoal purpose,
                                           const ReductionStage process) const {
  const std::string k_key = reductionKernelKey(prec, purpose, process);
  if (k_dictionary.find(k_key) == k_dictionary.end()) {
    rtErr("Reduction kernel identifier " + k_key + " was not found in the kernel map.",
          "CoreKlManager", "getReductionKernelDims");
  }
  return k_dictionary.at(k_key).getLaunchParameters();
}

//-------------------------------------------------------------------------------------------------
int2 CoreKlManager::getVirtualSiteKernelDims(const PrecisionModel prec,
                                             const VirtualSiteActivity purpose) const {
  const std::string k_key = virtualSiteKernelKey(prec, purpose, virtual_site_kernel_width);
  if (k_dictionary.find(k_key) == k_dictionary.end()) {
    rtErr("Virtual site handling kernel identifier " + k_key + " was not found in the kernel map.",
          "CoreKlManager", "getVirtualSiteKernelDims");
  }
  return k_dictionary.at(k_key).getLaunchParameters();
}

//-------------------------------------------------------------------------------------------------
int2 CoreKlManager::getRMSDKernelDims(const PrecisionModel prec, const RMSDTask order) const {
  const std::string k_key = rmsdKernelKey(prec, order);
  if (k_dictionary.find(k_key) == k_dictionary.end()) {
    rtErr("RMSD calculation kernel identifier " + k_key + " was not found in the kernel map.",
          "CoreKlManager", "getRMSDKernelDims");
  }
  return k_dictionary.at(k_key).getLaunchParameters();
}

//-------------------------------------------------------------------------------------------------
int2 CoreKlManager::getPMEPairsKernelDims(const PrecisionModel coord_prec,
                                          const PrecisionModel calc_prec,
                                          const NeighborListKind grid_configuration,
                                          const TinyBoxPresence has_tiny_box,
                                          const EvaluateForce eval_frc,
                                          const EvaluateEnergy eval_nrg,
                                          const ClashResponse mitigation) const {
  const std::string k_key = pmePairsKernelKey(coord_prec, calc_prec, grid_configuration,
                                              has_tiny_box, eval_frc, eval_nrg, mitigation);
  if (k_dictionary.find(k_key) == k_dictionary.end()) {
    rtErr("PME particle-particle pair interactions kernel identifier " + k_key + " was not found "
          "in the kernel map.", "CoreKlManager", "getPMEPairsKernelDims");
  }
  return k_dictionary.at(k_key).getLaunchParameters();
}

//-------------------------------------------------------------------------------------------------
int2 CoreKlManager::getMigrationKernelDims(const PrecisionModel coord_prec,
                                           const NeighborListKind grid_configuration,
                                           const int stage_number, const int gpos_bits,
                                           const int chain_count) const {
  const std::string k_key = migrationKernelKey(coord_prec, grid_configuration, stage_number,
                                               gpos_bits);
  if (k_dictionary.find(k_key) == k_dictionary.end()) {
    rtErr("PME particle migration kernel identifier " + k_key + " was not found in the kernel "
          "map.", "CoreKlManager", "getMigrationKernelDims");
  }

  // While the topology itself may offer sufficient clues about the workflow for modifying other
  // kernel launches, the chain count in CellGrid objects, which the developer will not be held
  // responsible for having available when the kernel manager is constructed, informs the
  // subdivision of the migration kernels.  Apply that here.
  int2 result = k_dictionary.at(k_key).getLaunchParameters();
  const int nsmp = gpu.getSMPCount();
  if (chain_count == 0 || chain_count > nsmp * 2) {
    result.x *= 2;
    result.y /= 2;
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
int nonbondedBlockMultiplier(const GpuDetails &gpu, const UnitCellType unit_cell,
                             const PrecisionModel prec, const ImplicitSolventModel igb) {
#ifdef STORMM_USE_HPC
#  ifdef STORMM_USE_CUDA
  switch (unit_cell) {
  case UnitCellType::NONE:
    switch (igb) {
    case ImplicitSolventModel::NONE:
      switch (prec) {
      case PrecisionModel::DOUBLE:
        return 3;
      case PrecisionModel::SINGLE:
        return (gpu.getArchMajor() == 7 && gpu.getArchMinor() >= 5) ? 4 : 5;
      }
      break;
    case ImplicitSolventModel::HCT_GB:
    case ImplicitSolventModel::OBC_GB:
    case ImplicitSolventModel::OBC_GB_II:
    case ImplicitSolventModel::NECK_GB:
    case ImplicitSolventModel::NECK_GB_II:
      switch (prec) {
      case PrecisionModel::DOUBLE:
        return 3;
      case PrecisionModel::SINGLE:
        return (gpu.getArchMajor() == 7 && gpu.getArchMinor() >= 5) ? 4 : 5;
      }
      break;
    }
    break;
  case UnitCellType::ORTHORHOMBIC:
  case UnitCellType::TRICLINIC:
    switch (prec) {
    case PrecisionModel::DOUBLE:
      return 2;
    case PrecisionModel::SINGLE:
      return (gpu.getArchMajor() == 7 && gpu.getArchMinor() >= 5) ? 2 : 3;
    }
  }
#  else
  // Other vendors are not known to make GPUs that have special requirements
  switch (unit_cell) {
  case UnitCellType::NONE:
    switch (igb) {
    case ImplicitSolventModel::NONE:
      switch (prec) {
      case PrecisionModel::DOUBLE:
        return 3;
      case PrecisionModel::SINGLE:
        return 5;
      }
      break;
    case ImplicitSolventModel::HCT_GB:
    case ImplicitSolventModel::OBC_GB:
    case ImplicitSolventModel::OBC_GB_II:
    case ImplicitSolventModel::NECK_GB:
    case ImplicitSolventModel::NECK_GB_II:
      switch (prec) {
      case PrecisionModel::DOUBLE:
        return 3;
      case PrecisionModel::SINGLE:
        return 5;
      }
      break;
    }
    break;
  case UnitCellType::ORTHORHOMBIC:
  case UnitCellType::TRICLINIC:
    return 3;
  }
#  endif
  __builtin_unreachable();
#else
  return 1;
#endif  
}

//-------------------------------------------------------------------------------------------------
int gbRadiiBlockMultiplier(const GpuDetails &gpu, const PrecisionModel prec) {
#ifdef STORMM_USE_HPC
  switch (prec) {
  case PrecisionModel::DOUBLE:
    return 3;
  case PrecisionModel::SINGLE:
    return (gpu.getArchMajor() == 7 && gpu.getArchMinor() >= 5) ? 4 : 5;
  }
  __builtin_unreachable();
#else
  return 1;
#endif  
}

//-------------------------------------------------------------------------------------------------
int gbDerivativeBlockMultiplier(const GpuDetails &gpu, const PrecisionModel prec) {
#ifdef STORMM_USE_HPC
  switch (prec) {
  case PrecisionModel::DOUBLE:
    return 3;
  case PrecisionModel::SINGLE:
    return (gpu.getArchMajor() == 7 && gpu.getArchMinor() >= 5) ? 4 : 5;
  }
  __builtin_unreachable();
#else
  return 1;
#endif  
}

//-------------------------------------------------------------------------------------------------
std::vector<int> densityMappingBlockMultiplier(const GpuDetails &gpu, const PrecisionModel prec,
                                               const QMapMethod approach) {
#ifdef STORMM_USE_HPC
  switch (approach) {
  case QMapMethod::ACC_SHARED:
    switch (prec) {
    case PrecisionModel::DOUBLE:
      return std::vector<int>(9, 4);
    case PrecisionModel::SINGLE:
      if (gpu.getArchMajor() == 7 && gpu.getArchMinor() >= 5) {
        return std::vector<int>(9, 4);
      }
      else {
        return std::vector<int>(9, 4);
      }
      break;
    }
    break;
  case QMapMethod::GENERAL_PURPOSE:
    switch (prec) {
    case PrecisionModel::DOUBLE:
      if (gpu.getArchMajor() == 7 && gpu.getArchMinor() >= 5) {
        return { 0, 4, 4, 4, 4, 4, 4, 4, 4 };
      }
      else {
        return { 0, 6, 6, 6, 6, 5, 4, 4, 4 };
      }
      break;
    case PrecisionModel::SINGLE:
      if (gpu.getArchMajor() == 7 && gpu.getArchMinor() >= 5) {
        return std::vector<int>(9, 4);
      }
      else {
        return { 0, 6, 6, 6, 6, 5, 4, 4, 4 };
      }
      break;
    }
    break;
  case QMapMethod::AUTOMATIC:
    rtErr("Specify a particular method for density accumulation.",
          "densityMappingBlockMultiplier");
  }
  __builtin_unreachable();
#else
  return std::vector<int>(9, 1);
#endif
}

//-------------------------------------------------------------------------------------------------
int reductionBlockMultiplier() {
#ifdef STORMM_USE_HPC
  return 4;
#else
  return 1;
#endif
}

//-------------------------------------------------------------------------------------------------
int virtualSiteBlockMultiplier(const PrecisionModel prec) {
#ifdef STORMM_USE_HPC
  switch (prec) {
  case PrecisionModel::DOUBLE:
    return 3;
  case PrecisionModel::SINGLE:
    return 4;
  }
  __builtin_unreachable();
#else
  return 1;
#endif
}

//-------------------------------------------------------------------------------------------------
int rmsdBlockMultiplier(const PrecisionModel prec) {
#ifdef STORMM_USE_HPC
  return 4;
#else
  return 1;
#endif
}

//-------------------------------------------------------------------------------------------------
int pmePairsBlockMultiplier(const GpuDetails &gpu, const PrecisionModel coord_prec,
                            const PrecisionModel calc_prec) {
#ifdef STORMM_USE_HPC
  return 2;
#else
  return 1;
#endif
}

//-------------------------------------------------------------------------------------------------
int migrationBlockMultiplier() {
#ifdef STORMM_USE_HPC
  return 2;
#else
  return 1;
#endif
}

//-------------------------------------------------------------------------------------------------
std::string valenceKernelWidthExtension(const PrecisionModel prec,
                                        const ValenceKernelSize kwidth) {

  // Only append a size modifier if the kernel operates in single-precision.  All kernels that
  // carry out instructions found in valence work units with double-precision arithmetic are
  // sized to accommodate the largest possible work unit atom capacity, "one size fits all."
  switch (prec) {
  case PrecisionModel::DOUBLE:
    return std::string("");
  case PrecisionModel::SINGLE:
    switch (kwidth) {
    case ValenceKernelSize::XL:
      return std::string("_xl");
    case ValenceKernelSize::LG:
      return std::string("_lg");
    case ValenceKernelSize::MD:
      return std::string("_md");
    case ValenceKernelSize::SM:
      return std::string("_sm");
    }
    break;
  }
  __builtin_unreachable();
}
  
//-------------------------------------------------------------------------------------------------
std::string valenceKernelKey(const PrecisionModel prec, const EvaluateForce eval_force,
                             const EvaluateEnergy eval_nrg, const AccumulationMethod acc_meth,
                             const VwuGoal purpose, const ClashResponse collision_handling,
                             const ValenceKernelSize kwidth) {
  std::string k_key("vale_");
  switch (prec) {
  case PrecisionModel::DOUBLE:
    k_key += "d";
    break;
  case PrecisionModel::SINGLE:
    k_key += "f";
    break;
  }
  switch (eval_force) {
  case EvaluateForce::YES:
    k_key += "f";
    break;
  case EvaluateForce::NO:
    k_key += "x";
    break;
  }
  switch (eval_nrg) {
  case EvaluateEnergy::YES:
    k_key += "e";
    break;
  case EvaluateEnergy::NO:
    k_key += "x";
    break;
  }
  if (eval_force == EvaluateForce::YES) {
    switch (acc_meth) {
    case AccumulationMethod::SPLIT:
      k_key += "s";
      break;
    case AccumulationMethod::WHOLE:
      k_key += "w";
      break;
    case AccumulationMethod::AUTOMATIC:
      break;
    }
    switch (purpose) {
    case VwuGoal::ACCUMULATE:
      k_key += "a";
      break;
    case VwuGoal::MOVE_PARTICLES:
      k_key += "m";
      break;
    }
  }
  switch (collision_handling) {
  case ClashResponse::NONE:
    k_key += "_cl";
    break;
  case ClashResponse::FORGIVE:
    k_key += "_nc";
    break;
  }
  k_key += valenceKernelWidthExtension(prec, kwidth);
  return k_key;
}

//-------------------------------------------------------------------------------------------------
std::string valenceKernelKey(const PrecisionModel prec, const PrecisionModel neighbor_prec,
                             const NeighborListKind neighbor_layout, const EvaluateEnergy eval_nrg,
                             const AccumulationMethod acc_meth,
                             const ClashResponse collision_handling) {
  std::string k_key = valenceKernelKey(prec, EvaluateForce::YES, eval_nrg, acc_meth,
                                       VwuGoal::MOVE_PARTICLES, collision_handling,
                                       ValenceKernelSize::XL);
  k_key += "_pme";
  switch (neighbor_layout) {
  case NeighborListKind::DUAL:
    k_key += "dual_";
    break;
  case NeighborListKind::MONO:
    k_key += "_";
    break;
  }
  switch (neighbor_prec) {
  case PrecisionModel::DOUBLE:
    k_key += "d";
    break;
  case PrecisionModel::SINGLE:
    k_key += "f";
    break;
  }
  return k_key;
}

//-------------------------------------------------------------------------------------------------
std::string integrationKernelKey(PrecisionModel prec, AccumulationMethod acc_meth, 
                                 ValenceKernelSize kwidth, IntegrationStage process) {
  std::string k_key("intg_");
  switch (prec) {
  case PrecisionModel::DOUBLE:
    k_key += "d";
    break;
  case PrecisionModel::SINGLE:
    k_key += "f";
    break;
  }
  switch (acc_meth) {
  case AccumulationMethod::SPLIT:
    k_key += "s";
    break;
  case AccumulationMethod::WHOLE:
    k_key += "w";
    break;
  case AccumulationMethod::AUTOMATIC:
    break;
  }
  k_key += valenceKernelWidthExtension(prec, kwidth);
  switch (process) {
  case IntegrationStage::VELOCITY_ADVANCE:
    k_key += "_va";
    break;
  case IntegrationStage::VELOCITY_CONSTRAINT:
    k_key += "_vc";
    break;
  case IntegrationStage::POSITION_ADVANCE:
    k_key += "_pa";
    break;
  case IntegrationStage::GEOMETRY_CONSTRAINT:
    k_key += "_gc";
    break;
  }
  return k_key;
}

//-------------------------------------------------------------------------------------------------
std::string nonbondedKernelKey(const PrecisionModel prec, const NbwuKind kind,
                               const EvaluateForce eval_force, const EvaluateEnergy eval_nrg,
                               const AccumulationMethod acc_meth, const ImplicitSolventModel igb,
                               const ClashResponse collision_handling) {
  std::string k_key("nonb_");
  switch (prec) {
  case PrecisionModel::DOUBLE:
    k_key += "d";
    break;
  case PrecisionModel::SINGLE:
    k_key += "f";
    break;
  }
  switch (kind) {
  case NbwuKind::TILE_GROUPS:
    k_key += "tg";
    break;
  case NbwuKind::SUPERTILES:
    k_key += "st";
    break;
  case NbwuKind::HONEYCOMB:
    k_key += "hc";
    break;
  case NbwuKind::UNKNOWN:
    rtWarn("No kernels are available for an " + getEnumerationName(kind) + " non-bonded work unit "
           "layout.  This is likely the result of a failure to use the loadNonbondedWorkUnits() "
           "function in the AtomGraphSynthesis with which this kernel manager is associated.",
           "nonbondedKernelKey");
    break;
  }
  switch (igb) {
  case ImplicitSolventModel::NONE:
    k_key += "_vac_";
    break;
  case ImplicitSolventModel::HCT_GB:
  case ImplicitSolventModel::OBC_GB:
  case ImplicitSolventModel::OBC_GB_II:
    k_key += "_gbs_";
    break;
  case ImplicitSolventModel::NECK_GB:
  case ImplicitSolventModel::NECK_GB_II:
    k_key += "_gbn_";
    break;
  }
  switch (eval_force) {
  case EvaluateForce::YES:
    k_key += "f";
    break;
  case EvaluateForce::NO:
    k_key += "x";
    break;
  }
  switch (eval_nrg) {
  case EvaluateEnergy::YES:
    k_key += "e";
    break;
  case EvaluateEnergy::NO:
    k_key += "x";
    break;
  }
  if (eval_force == EvaluateForce::YES) {
    switch (acc_meth) {
    case AccumulationMethod::SPLIT:
      k_key += "s";
      break;
    case AccumulationMethod::WHOLE:
      k_key += "w";
      break;
    case AccumulationMethod::AUTOMATIC:
      break;
    }
  }
  switch (collision_handling) {
  case ClashResponse::NONE:
    k_key += "_cl";
    break;
  case ClashResponse::FORGIVE:
    k_key += "_nc";
    break;
  }
  return k_key;
}

//-------------------------------------------------------------------------------------------------
std::string appendBornKernelKey(const PrecisionModel prec, const NbwuKind kind,
                                const AccumulationMethod acc_meth,
                                const ImplicitSolventModel igb) {
  std::string app_key("");
  switch (prec) {
  case PrecisionModel::DOUBLE:
    app_key += "d";
    break;
  case PrecisionModel::SINGLE:
    app_key += "f";
    break;
  }
  switch (kind) {
  case NbwuKind::TILE_GROUPS:
    app_key += "tg";
    break;
  case NbwuKind::SUPERTILES:
    app_key += "st";
    break;
  case NbwuKind::HONEYCOMB:
    app_key += "hc";
    break;
  case NbwuKind::UNKNOWN:
    break;
  }
  switch (acc_meth) {
  case AccumulationMethod::SPLIT:
    app_key += "s";
    break;
  case AccumulationMethod::WHOLE:
    app_key += "w";
    break;
  case AccumulationMethod::AUTOMATIC:
    break;
  }
  switch (igb) {
  case ImplicitSolventModel::NONE:
    app_key += "_vac";
    break;
  case ImplicitSolventModel::HCT_GB:
  case ImplicitSolventModel::OBC_GB:
  case ImplicitSolventModel::OBC_GB_II:
    app_key += "_gbs";
    break;
  case ImplicitSolventModel::NECK_GB:
  case ImplicitSolventModel::NECK_GB_II:
    app_key += "_gbn";
    break;
  }
  return app_key;
}

//-------------------------------------------------------------------------------------------------
std::string bornRadiiKernelKey(const PrecisionModel prec, const NbwuKind kind,
                               const AccumulationMethod acc_meth, const ImplicitSolventModel igb) {
  std::string k_key("gbrd_");
  return k_key + appendBornKernelKey(prec, kind, acc_meth, igb);
}

//-------------------------------------------------------------------------------------------------
std::string bornDerivativeKernelKey(const PrecisionModel prec, const NbwuKind kind,
                                    const AccumulationMethod acc_meth,
                                    const ImplicitSolventModel igb) {
  std::string k_key("gbdv_");
  return k_key + appendBornKernelKey(prec, kind, acc_meth, igb);
}

//-------------------------------------------------------------------------------------------------
std::string appendQMapKernelKey(const PrecisionModel prec, const size_t cg_tmat,
                                const int order) {
  std::string result;
  if (cg_tmat == int_type_index) {
    result += "si_";
  }
  else if (cg_tmat == llint_type_index) {
    result += "li_";
  }
  else if (cg_tmat == float_type_index) {
    result += "sr_";
  }
  else if (cg_tmat == double_type_index) {
    result += "lr_";
  }
  switch (prec) {
  case PrecisionModel::DOUBLE:
    result += "d_";
    break;
  case PrecisionModel::SINGLE:
    result += "f_";
    break;
  }
  result += std::to_string(order);
  return result;

}
  
//-------------------------------------------------------------------------------------------------
std::string generalQMapKernelKey(const PrecisionModel prec, const size_t cg_tmat,
                                 const int order) {
  return "qmap_" + appendQMapKernelKey(prec, cg_tmat, order);
}

//-------------------------------------------------------------------------------------------------
std::string shrAccQMapKernelKey(const PrecisionModel calc_prec, const PrecisionModel acc_prec,
                                const bool overflow, const size_t cg_tmat, const int order) {
  std::string base_key("qmap_sacc_");
  switch (acc_prec) {
  case PrecisionModel::DOUBLE:
    base_key += (overflow) ? "d_xf_" : "d_sf_";
    break;
  case PrecisionModel::SINGLE:
    base_key += (overflow) ? "s_xf_" : "s_sf_";
    break;
  }
  return base_key + appendQMapKernelKey(calc_prec, cg_tmat, order);
}

//-------------------------------------------------------------------------------------------------
std::string reductionKernelKey(const PrecisionModel prec, const ReductionGoal purpose,
                               const ReductionStage process) {
  std::string k_key("redc_");
  switch (prec) {
  case PrecisionModel::DOUBLE:
    k_key += "d";
    break;
  case PrecisionModel::SINGLE:
    k_key += "f";
    break;
  }
  switch (purpose) {
  case ReductionGoal::NORMALIZE:
  case ReductionGoal::CENTER_ON_ZERO:
    break;
  case ReductionGoal::CONJUGATE_GRADIENT:
    k_key += "_cg";
    break;
  }
  switch (process) {
  case ReductionStage::GATHER:
    k_key += "_gt";
    break;
  case ReductionStage::SCATTER:
    k_key += "_sc";
    break;
  case ReductionStage::RESCALE:
    k_key += "_rs";
    break;
  case ReductionStage::ALL_REDUCE:
    k_key += "_rd";
    break;
  }
  return k_key;
}

//-------------------------------------------------------------------------------------------------
std::string virtualSiteKernelKey(const PrecisionModel prec, const VirtualSiteActivity purpose,
                                 const ValenceKernelSize kwidth) {
  std::string k_key("vste_");
  switch (prec) {
  case PrecisionModel::DOUBLE:
    k_key += "d";
    break;
  case PrecisionModel::SINGLE:
    k_key += "f";
    break;
  }
  switch (purpose) {
  case VirtualSiteActivity::PLACEMENT:
    k_key += "_pl";
    break;
  case VirtualSiteActivity::TRANSMIT_FORCES:
    k_key += "_xm";
    break;
  }
  k_key += valenceKernelWidthExtension(prec, kwidth);
  return k_key;
}

//-------------------------------------------------------------------------------------------------
std::string rmsdKernelKey(const PrecisionModel prec, const RMSDTask order) {
  std::string k_key("rmsd_");
  switch (prec) {
  case PrecisionModel::DOUBLE:
    k_key += "d";
    break;
  case PrecisionModel::SINGLE:
    k_key += "f";
    break;
  }
  switch (order) {
  case RMSDTask::REFERENCE:
    k_key += "_r";
    break;
  case RMSDTask::MATRIX:
    k_key += "_m";
    break;
  }
  return k_key;
}

//-------------------------------------------------------------------------------------------------
std::string pmePairsKernelKey(const PrecisionModel coord_prec, const PrecisionModel calc_prec,
                              const NeighborListKind grid_configuration,
                              const TinyBoxPresence has_tiny_box, const EvaluateForce eval_frc,
                              const EvaluateEnergy eval_nrg, const ClashResponse mitigation) {
  std::string k_key("pme_tp_");
  switch (coord_prec) {
  case PrecisionModel::DOUBLE:
    k_key += "d";
    break;
  case PrecisionModel::SINGLE:
    k_key += "f";
    break;
  }
  switch (calc_prec) {
  case PrecisionModel::DOUBLE:
    k_key += "d";
    break;
  case PrecisionModel::SINGLE:
    k_key += "f";
    break;
  }
  switch (grid_configuration) {
  case NeighborListKind::DUAL:
    k_key += "pr";
    break;
  case NeighborListKind::MONO:
    k_key += "un";
    break;
  }
  switch (eval_frc) {
  case EvaluateForce::YES:
    switch (eval_nrg) {
    case EvaluateEnergy::YES:
      k_key += "fe";
      break;
    case EvaluateEnergy::NO:
      k_key += "f";
      break;
    }
    break;
  case EvaluateForce::NO:
    switch (eval_nrg) {
    case EvaluateEnergy::YES:
      k_key += "e";
      break;
    case EvaluateEnergy::NO:

      // This erroneous case will be trapped elsewhere
      break;
    }
    break;
  }
  switch (has_tiny_box) {
  case TinyBoxPresence::YES:
    k_key += "t";
    break;
  case TinyBoxPresence::NO:
    k_key += "s";
    break;
  }
  switch (mitigation) {
  case ClashResponse::FORGIVE:
    k_key += "nc";
    break;
  case ClashResponse::NONE:
    k_key += "cl";
    break;
  }
  return k_key;
}

//-------------------------------------------------------------------------------------------------
std::string migrationKernelKey(const PrecisionModel coord_prec,
                               const NeighborListKind grid_configuration,
                               const int stage_number, const int gpos_bits) {
  std::string k_key("migr_");
  switch (coord_prec) {
  case PrecisionModel::DOUBLE:
    k_key += "d";
    break;
  case PrecisionModel::SINGLE:
    k_key += "f";
    break;
  }
  switch (grid_configuration) {
  case NeighborListKind::DUAL:
    k_key += "pr";
    break;
  case NeighborListKind::MONO:
    k_key += "un";
    break;
  }
  if (stage_number < 1 || stage_number > 2) {
    rtErr("The accepted stages of particle migration are 1 and 2.  " +
          std::to_string(stage_number) + " is invalid.", "migrationKernelKey");
  }
  k_key += std::to_string(stage_number);
  if (gpos_bits > 0 && gpos_bits <= globalpos_scale_nonoverflow_bits) {
    k_key += "_lr";
  }
  else {
    k_key += "_fn";
  }
  return k_key;
}

} // namespace card
} // namespace stormm
