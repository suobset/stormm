// -*-c++-*-
#include "copyright.h"
#include "Accelerator/hybrid.h"
#include "Constants/hpc_bounds.h"
#include "DataTypes/stormm_vector_types.h"
#include "Synthesis/synthesis_enumerators.h"
#include "Synthesis/valence_workunit.h"
#include "Topology/atomgraph_enumerators.h"
#include "hpc_virtual_site_handling.h"

namespace stormm {
namespace structure {

using card::HybridTargetLevel;
using numerics::max_llint_accumulation;
using numerics::max_int_accumulation_f;
using synthesis::maximum_valence_work_unit_atoms;
using synthesis::vwu_abstract_length;
using synthesis::VwuAbstractMap;
using topology::VirtualSiteKind;
  
#include "Math/rounding.cui"
#include "Math/vector_formulas.cui"
#include "Numerics/accumulation.cui"
  
#define VSITE_STANDALONE
#define TCALC  float
#  define VSITE_STANDALONE_BLOCK_MULTIPLIER 4
#  define TCALC2 float2
#  define TCALC3 float3
#  define TCALC4 float4
#  define TCALC_IS_SINGLE
#  define COS_FUNC cosf
#  define SIN_FUNC sinf
#  define SQRT_FUNC sqrtf
#  define LLCONV_FUNC __float2ll_rn
#    define KERNEL_NAME kfPlaceVirtualSites
#      include "virtual_site_placement.cui"
#    undef KERNEL_NAME
#    define KERNEL_NAME kfTransmitVSiteForces
#      include "virtual_site_transmission.cui"
#    undef KERNEL_NAME
#  undef TCALC2
#  undef TCALC3
#  undef TCALC4
#  undef TCALC_IS_SINGLE
#  undef COS_FUNC
#  undef SIN_FUNC
#  undef SQRT_FUNC
#  undef LLCONV_FUNC
#  undef VSITE_STANDALONE_BLOCK_MULTIPLIER
#undef TCALC

#define TCALC  double
#  define VSITE_STANDALONE_BLOCK_MULTIPLIER 3
#  define TCALC2 double2
#  define TCALC3 double3
#  define TCALC4 double4
#  define COS_FUNC cos
#  define SIN_FUNC sin
#  define SQRT_FUNC sqrt
#    define KERNEL_NAME kdPlaceVirtualSites
#      include "virtual_site_placement.cui"
#    undef KERNEL_NAME
#    define KERNEL_NAME kdTransmitVSiteForces
#      include "virtual_site_transmission.cui"
#    undef KERNEL_NAME
#  undef TCALC2
#  undef TCALC3
#  undef TCALC4
#  undef COS_FUNC
#  undef SIN_FUNC
#  undef SQRT_FUNC
#  undef VSITE_STANDALONE_BLOCK_MULTIPLIER
#undef TCALC
#undef VSITE_STANDALONE

//-------------------------------------------------------------------------------------------------
extern cudaFuncAttributes queryVirtualSiteKernelRequirements(const PrecisionModel prec,
                                                             const VirtualSiteActivity purpose) {
  cudaFuncAttributes result;
  switch (prec) {
  case PrecisionModel::DOUBLE:
    switch (purpose) {
    case VirtualSiteActivity::PLACEMENT:
      if (cudaFuncGetAttributes(&result, kdPlaceVirtualSites) != cudaSuccess) {
        rtErr("Error obtaining attributes for kernel kdVirtualSitePlacement.",
              "queryVirtualSiteKernelRequirements");
      }
      break;
    case VirtualSiteActivity::TRANSMIT_FORCES:
      if (cudaFuncGetAttributes(&result, kdTransmitVSiteForces) != cudaSuccess) {
        rtErr("Error obtaining attributes for kernel kdTransmitVSiteForces.",
              "queryVirtualSiteKernelRequirements");
      }
      break;
    }
    break;
  case PrecisionModel::SINGLE:
    switch (purpose) {
    case VirtualSiteActivity::PLACEMENT:
      if (cudaFuncGetAttributes(&result, kfPlaceVirtualSites) != cudaSuccess) {
        rtErr("Error obtaining attributes for kernel kfVirtualSitePlacement.",
              "queryVirtualSiteKernelRequirements");
      }
      break;
    case VirtualSiteActivity::TRANSMIT_FORCES:
      if (cudaFuncGetAttributes(&result, kfTransmitVSiteForces) != cudaSuccess) {
        rtErr("Error obtaining attributes for kernel kfTransmitVSiteForces.",
              "queryVirtualSiteKernelRequirements");
      }
      break;
    }
    break;
  }
  return result;
}
  
//-------------------------------------------------------------------------------------------------
void launchVirtualSitePlacement(PsSynthesisWriter *poly_psw, CacheResourceKit<double> *gmem_r,
                                const SyValenceKit<double> &poly_vk,
                                const SyAtomUpdateKit<double, double2, double4> &poly_auk,
                                const int2 bt) {
  kdPlaceVirtualSites<<<bt.x, bt.y>>>(*poly_psw, poly_vk, poly_auk, *gmem_r);
}

//-------------------------------------------------------------------------------------------------
void launchVirtualSitePlacement(PsSynthesisWriter *poly_psw, CacheResourceKit<float> *gmem_r,
                                const SyValenceKit<float> &poly_vk,
                                const SyAtomUpdateKit<float, float2, float4> &poly_auk,
                                const int2 bt) {
  kfPlaceVirtualSites<<<bt.x, bt.y>>>(*poly_psw, poly_vk, poly_auk, *gmem_r);
}

//-------------------------------------------------------------------------------------------------
void launchTransmitVSiteForces(PsSynthesisWriter *poly_psw, CacheResourceKit<double> *gmem_r,
                               const SyValenceKit<double> &poly_vk,
                               const SyAtomUpdateKit<double, double2, double4> &poly_auk,
                               const int2 bt) {
  kdTransmitVSiteForces<<<bt.x, bt.y>>>(*poly_psw, poly_vk, poly_auk, *gmem_r);
}

//-------------------------------------------------------------------------------------------------
void launchTransmitVSiteForces(PsSynthesisWriter *poly_psw, CacheResourceKit<float> *gmem_r,
                               const SyValenceKit<float> &poly_vk,
                               const SyAtomUpdateKit<float, float2, float4> &poly_auk,
                               const int2 bt) {
  kfTransmitVSiteForces<<<bt.x, bt.y>>>(*poly_psw, poly_vk, poly_auk, *gmem_r);
}

//-------------------------------------------------------------------------------------------------
void launchVirtualSiteHandling(const PrecisionModel prec, const VirtualSiteActivity purpose,
                               PhaseSpaceSynthesis *poly_ps, CacheResource *tb_space,
                               const AtomGraphSynthesis &poly_ag, const CoreKlManager &launcher) {

  // Bail out if there are no virtual sites in the synthesis
  if (poly_ag.getVirtualSiteCount() == 0) {
    return;
  }
  const HybridTargetLevel devc_tier = HybridTargetLevel::DEVICE;
  PsSynthesisWriter poly_psw = poly_ps->data(devc_tier);
  const int2 lp = launcher.getVirtualSiteKernelDims(prec, VirtualSiteActivity::PLACEMENT);
  switch (prec) {
  case PrecisionModel::DOUBLE:
    {
      CacheResourceKit<double> gmem_r = tb_space->dpData(devc_tier);
      const SyValenceKit<double> poly_vk = poly_ag.getDoublePrecisionValenceKit(devc_tier);
      const SyAtomUpdateKit<double, double2, double4> poly_auk =
        poly_ag.getDoublePrecisionAtomUpdateKit(devc_tier);
      switch (purpose) {
      case VirtualSiteActivity::PLACEMENT:
        launchVirtualSitePlacement(&poly_psw, &gmem_r, poly_vk, poly_auk, lp);
        break;
      case VirtualSiteActivity::TRANSMIT_FORCES:
        launchTransmitVSiteForces(&poly_psw, &gmem_r, poly_vk, poly_auk, lp);
        break;
      }
    }
    break;
  case PrecisionModel::SINGLE:
    {
      CacheResourceKit<float> gmem_r = tb_space->spData(devc_tier);
      const SyValenceKit<float> poly_vk = poly_ag.getSinglePrecisionValenceKit(devc_tier);
      const SyAtomUpdateKit<float, float2, float4> poly_auk =
        poly_ag.getSinglePrecisionAtomUpdateKit(devc_tier);
      switch (purpose) {
      case VirtualSiteActivity::PLACEMENT:
        launchVirtualSitePlacement(&poly_psw, &gmem_r, poly_vk, poly_auk, lp);
        break;
      case VirtualSiteActivity::TRANSMIT_FORCES:
        launchTransmitVSiteForces(&poly_psw, &gmem_r, poly_vk, poly_auk, lp);
        break;
      }
    }
    break;
  }
}

} // namespace structure
} // namespace stormm
