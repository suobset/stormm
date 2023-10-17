#include "copyright.h"
#include "map_density.h"

namespace stormm {
namespace energy {

//-------------------------------------------------------------------------------------------------
void matchThemes(const NonbondedTheme pm_theme, const NonbondedTheme cg_theme) {
  bool problem = false;
  switch (pm_theme) {
  case NonbondedTheme::ELECTROSTATIC:
    switch (cg_theme) {
    case NonbondedTheme::ELECTROSTATIC:
    case NonbondedTheme::ALL:
      break;
    case NonbondedTheme::VAN_DER_WAALS:
      problem = true;
    }
    break;
  case NonbondedTheme::VAN_DER_WAALS:
    switch (cg_theme) {
    case NonbondedTheme::ELECTROSTATIC:
      problem = true;
    case NonbondedTheme::VAN_DER_WAALS:
    case NonbondedTheme::ALL:
      break;
    }
    break;
  case NonbondedTheme::ALL:
    rtErr("A particle-mesh interaction grid cannot carry more than one kind of non-bonded "
          "potential.", "matchThemes");
    break;
  }
}

#ifdef STORMM_USE_HPC
//-------------------------------------------------------------------------------------------------
void mapDensity(PMIGridWriter *pm_wrt, PMIGridAccumulator *pm_acc, MMControlKit<double> *ctrl,
                const CellGridReader<void, void, void, void> &v_cgr, const size_t cg_tmat,
                const SyNonbondedKit<double, double2> &synbk, const int block_count, const int2 lp,
                const QMapMethod approach, PMIGrid *pm) {
  const HybridTargetLevel devc_tier = HybridTargetLevel::DEVICE;
  switch (approach) {
  case QMapMethod::ACC_SHARED:
    launchShrAccDensityKernel(pm_wrt, ctrl, v_cgr, cg_tmat, synbk, lp);
    break;
  case QMapMethod::GENERAL_PURPOSE:
    launchPMIGridInitialization(pm_acc, block_count);
    launchGenPrpDensityKernel(pm_acc, v_cgr, cg_tmat, synbk, lp);
    launchPMIGridRealConversion(pm_wrt, *pm_acc, block_count);    
    break;
  case QMapMethod::AUTOMATIC:
    break;
  }

  // Mark that the particle-mesh interaction grids are now presented in real format
  pm->setRealDataFormat();
}

//-------------------------------------------------------------------------------------------------
void mapDensity(PMIGridWriter *pm_wrt, PMIGridAccumulator *pm_acc, MMControlKit<float> *ctrl,
                const CellGridReader<void, void, void, void> &v_cgr, const size_t cg_tmat,
                const SyNonbondedKit<float, float2> &synbk, const int block_count, const int2 lp,
                const QMapMethod approach, PMIGrid *pm) {
  const HybridTargetLevel devc_tier = HybridTargetLevel::DEVICE;
  switch (approach) {
  case QMapMethod::ACC_SHARED:
    launchShrAccDensityKernel(pm_wrt, ctrl, v_cgr, cg_tmat, synbk, lp);
    break;
  case QMapMethod::GENERAL_PURPOSE:
    launchPMIGridInitialization(pm_acc, block_count);
    launchGenPrpDensityKernel(pm_acc, v_cgr, cg_tmat, synbk, lp);
    launchPMIGridRealConversion(pm_wrt, *pm_acc, block_count);    
    break;
  case QMapMethod::AUTOMATIC:
    break;
  }

  // Mark that the particle-mesh interaction grids are now presented in real format
  pm->setRealDataFormat();
}
#endif

//-------------------------------------------------------------------------------------------------
void mapDensity(PMIGrid *pm, const AtomGraphSynthesis *poly_ag) {

  // Extract the cell grid pointer and unroll its templating
  const size_t cg_tmat  = pm->getCellGridMatrixTypeID();
  const size_t cg_tacc  = pm->getCellGridAccumulatorTypeID();
  const size_t cg_tcalc = pm->getCellGridCalculationTypeID();
  if (cg_tmat == double_type_index) {

    // The type of the cell dimension matrices implies the format of the coordinate data (a
    // four-tuple of the matrix dimension data type).
    unrollMapDensityCall<double, double4>(pm, cg_tacc, cg_tcalc, poly_ag);
  }
  else if (cg_tmat == float_type_index) {
    unrollMapDensityCall<float, float4>(pm, cg_tacc, cg_tcalc, poly_ag);
  }
  else if (cg_tmat == llint_type_index) {
    unrollMapDensityCall<llint, llint4>(pm, cg_tacc, cg_tcalc, poly_ag);
  }
  else if (cg_tmat == int_type_index) {
    unrollMapDensityCall<int, int4>(pm, cg_tacc, cg_tcalc, poly_ag);
  }
}

//-------------------------------------------------------------------------------------------------
void mapDensity(PMIGrid *pm, const AtomGraphSynthesis &poly_ag) {
  mapDensity(pm, poly_ag.getSelfPointer());
}

//-------------------------------------------------------------------------------------------------
std::vector<double> mapDensity(const CoordinateFrame *cf, const AtomGraph *ag,
                               const NonbondedTheme theme, const int grid_dim_a,
                               const int grid_dim_b, const int grid_dim_c, const int order) {
  return mapDensity<double>(cf->data(), ag->getDoublePrecisionNonbondedKit(), theme, grid_dim_a,
                            grid_dim_b, grid_dim_c, order);
}

//-------------------------------------------------------------------------------------------------
std::vector<double> mapDensity(const CoordinateFrame &cf, const AtomGraph &ag,
                               const NonbondedTheme theme, const int grid_dim_a,
                               const int grid_dim_b, const int grid_dim_c, const int order) {
  return mapDensity<double>(cf.data(), ag.getDoublePrecisionNonbondedKit(), theme, grid_dim_a,
                            grid_dim_b, grid_dim_c, order);
}

} // namespace energy
} // namespace stormm
