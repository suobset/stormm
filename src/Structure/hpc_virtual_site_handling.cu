// -*-c++-*-
#include "copyright.h"
#include "Constants/hpc_bounds.h"
#include "DataTypes/stormm_vector_types.h"
#include "Synthesis/synthesis_enumerators.h"
#include "Synthesis/valence_workunit.h"
#include "Topology/atomgraph_enumerators.h"
#include "hpc_virtual_site_handling.h"

namespace stormm {
namespace structure {

using data_types::int95_t;
using numerics::max_llint_accumulation;
using numerics::max_int_accumulation_f;
using synthesis::vwu_abstract_length;
using synthesis::VwuAbstractMap;
using topology::VirtualSiteKind;
  
#include "Math/rounding.cui"
#include "Math/vector_formulas.cui"
#include "Numerics/accumulation.cui"
  
#define VSITE_STANDALONE
#define TCALC  float
#define TCALC2 float2
#define TCALC3 float3
#define TCALC4 float4
#define TCALC_IS_SINGLE
#define COS_FUNC cosf
#define SIN_FUNC sinf
#define SQRT_FUNC sqrtf
#define LLCONV_FUNC __float2ll_rn
#  define KERNEL_NAME kfPlaceVirtualSites
#    include "virtual_site_placement.cui"
#  undef KERNEL_NAME
#undef TCALC
#undef TCALC2
#undef TCALC3
#undef TCALC4
#undef TCALC_IS_SINGLE
#undef COS_FUNC
#undef SIN_FUNC
#undef SQRT_FUNC
#undef LLCONV_FUNC

#define TCALC  double
#define TCALC2 double2
#define TCALC3 double3
#define TCALC4 double4
#define COS_FUNC cos
#define SIN_FUNC sin
#define SQRT_FUNC sqrt
#  define KERNEL_NAME kdPlaceVirtualSites
#    include "virtual_site_placement.cui"
#  undef KERNEL_NAME
#undef TCALC
#undef TCALC2
#undef TCALC3
#undef TCALC4
#undef COS_FUNC
#undef SIN_FUNC
#undef SQRT_FUNC
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
              "queryValenceKernelRequirements");
      }
      break;
    case VirtualSiteActivity::TRANSMIT_FORCES:
      break;
    }
    break;
  case PrecisionModel::SINGLE:
    switch (purpose) {
    case VirtualSiteActivity::PLACEMENT:
      if (cudaFuncGetAttributes(&result, kfPlaceVirtualSites) != cudaSuccess) {
        rtErr("Error obtaining attributes for kernel kfVirtualSitePlacement.",
              "queryValenceKernelRequirements");
      }
      break;
    case VirtualSiteActivity::TRANSMIT_FORCES:
      break;
    }
    break;
  }
  return result;
}
  
//-------------------------------------------------------------------------------------------------
void launchVirtualSitePlacement(PsSynthesisWriter *poly_psw, CacheResourceKit<double> *gmem_r,
                                const SyValenceKit<double> &poly_vk,
                                const SyAtomUpdateKit<double2, double4> &poly_auk, const int2 bt) {
  kdPlaceVirtualSites<<<bt.x, bt.y>>>(*poly_psw, poly_vk, poly_auk, *gmem_r);
}

//-------------------------------------------------------------------------------------------------
void launchVirtualSitePlacement(PsSynthesisWriter *poly_psw, CacheResourceKit<float> *gmem_r,
                                const SyValenceKit<float> &poly_vk,
                                const SyAtomUpdateKit<float2, float4> &poly_auk, const int2 bt) {
  kfPlaceVirtualSites<<<bt.x, bt.y>>>(*poly_psw, poly_vk, poly_auk, *gmem_r);
}

//-------------------------------------------------------------------------------------------------
void launchVirtualSitePlacement(PrecisionModel prec, PhaseSpaceSynthesis *poly_ps,
                                CacheResource *tb_space, const AtomGraphSynthesis &poly_ag,
                                const KernelManager &launcher) {
  PsSynthesisWriter poly_psw = poly_ps->data();
  const int2 lp = launcher.getVirtualSiteKernelDims(prec, VirtualSiteActivity::PLACEMENT);
  switch (prec) {
  case PrecisionModel::DOUBLE:
    {
      CacheResourceKit<double> gmem_r = tb_space->dpData();
      const SyValenceKit<double> poly_vk = poly_ag.getDoublePrecisionValenceKit();
      const SyAtomUpdateKit<double2, double4> poly_auk = poly_ag.getDoublePrecisionAtomUpdateKit();
      launchVirtualSitePlacement(&poly_psw, &gmem_r, poly_vk, poly_auk, lp);
    }
    break;
  case PrecisionModel::SINGLE:
    {
      CacheResourceKit<float> gmem_r = tb_space->spData();
      const SyValenceKit<float> poly_vk = poly_ag.getSinglePrecisionValenceKit();
      const SyAtomUpdateKit<float2, float4> poly_auk = poly_ag.getSinglePrecisionAtomUpdateKit();
      launchVirtualSitePlacement(&poly_psw, &gmem_r, poly_vk, poly_auk, lp);
    }
    break;
  }
}

} // namespace structure
} // namespace stormm
