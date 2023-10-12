#include "copyright.h"
#include "mapping_resource.h"

namespace stormm {
namespace energy {

//-------------------------------------------------------------------------------------------------
MappingResourceKit::MappingResourceKit(const int max_blocks_in, const int max_grid_points_in,
                                       const PrecisionModel mode_in, const QMapMethod approach_in,
                                       int* overflow_in, double* dstash_in, float* fstash_in) : 
    max_blocks{max_blocks_in}, max_grid_points{max_grid_points_in}, mode{mode_in},
    approach{approach_in}, overflow{overflow_in}, dstash{dstash_in}, fstash{fstash_in}
{}

//-------------------------------------------------------------------------------------------------
MappingResource::MappingResource(const PMIGrid *pm, const CoreKlManager &launcher,
                                 const ExceptionResponse policy) :
    block_limit{0}, grid_points_per_block{0}, mode{pm->getMode()},
    approach{pm->getWorkUnitConfiguration()},
    overflow{HybridKind::ARRAY, "maprsrc_overflow"},
    double_buffer{HybridKind::ARRAY, "maprsrc_dstash"},
    float_buffer{HybridKind::ARRAY, "maprsrc_fstash"}
{
  const size_t cg_tcalc = pm->getCellGridCalculationTypeID();
  
  // Determine the number of blocks for the kernel that will be launched.  Allocate arrays as
  // needed.
  int2 lp;
  if (cg_tcalc == double_type_index) {
    lp = launcher.getDensityMappingKernelDims(approach, PrecisionModel::DOUBLE, mode,
                                              pm->getCellGridMatrixTypeID(),
                                              pm->getInterpolationOrder());
  }
  else if (cg_tcalc == float_type_index) {
    lp = launcher.getDensityMappingKernelDims(approach, PrecisionModel::SINGLE, mode,
                                              pm->getCellGridMatrixTypeID(),
                                              pm->getInterpolationOrder());
  }
  else {
    rtErr("Only " + getStormmScalarTypeName<double>() + " and " +
          getStormmScalarTypeName<float>() + " types are allowed in the GPU density mapping "
          "kernels.", "MappingResource");
  }
  switch (approach) {
  case QMapMethod::GENERAL_PURPOSE:
  case QMapMethod::AUTOMATIC:

    // It may be an error to create a mapping resource for kernels that will not require any such
    // thing.  However, allow the object to exist in its present, un-allocated state for developers
    // that wish to always have an object even when using general-purpose density mapping kernels.
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("A mapping resource should not be created for a set of particle-mesh interaction "
            "grids set to use the general-purpose mapping kernel, or an unspecified approach.",
            "MappingResource");
    case ExceptionResponse::WARN:
      rtWarn("A mapping resource is unecessary for a set of particle-mesh interaction grids set "
             "to use the general-purpose mapping kernel, or an unspecified approach.",
            "MappingResource");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
    break;
  case QMapMethod::ACC_REGISTER:
    block_limit = lp.x;
    grid_points_per_block = pm->getLargestWorkUnitGridPoints();
    switch (mode) {
    case PrecisionModel::DOUBLE:
      double_buffer.resize(block_limit * grid_points_per_block);
      break;
    case PrecisionModel::SINGLE:
      float_buffer.resize(block_limit * grid_points_per_block);
      break;
    }
    break;
  case QMapMethod::ACC_SHARED:
    block_limit = lp.x;
    grid_points_per_block = pm->getLargestWorkUnitGridPoints();
    overflow.resize(block_limit * grid_points_per_block);
    break;
  }
}

//-------------------------------------------------------------------------------------------------
MappingResource::MappingResource(const PMIGrid &pm, const CoreKlManager &launcher,
                                 const ExceptionResponse policy) :
  MappingResource(pm.getSelfPointer(), launcher, policy)
{}

//-------------------------------------------------------------------------------------------------
int MappingResource::getBlockCount() const {
  return block_limit;
}
  
//-------------------------------------------------------------------------------------------------
int MappingResource::getGridPointsPerBlock() const {
  return grid_points_per_block;
}

//-------------------------------------------------------------------------------------------------
PrecisionModel MappingResource::getMode() const {
  return mode;
}

//-------------------------------------------------------------------------------------------------
QMapMethod MappingResource::getApproach() const {
  return approach;
}

//-------------------------------------------------------------------------------------------------
MappingResourceKit MappingResource::data(const HybridTargetLevel tier) {
  return MappingResourceKit(block_limit, grid_points_per_block, mode, approach,
                            overflow.data(tier), double_buffer.data(tier),
                            float_buffer.data(tier));
}

} // namespace energy
} // namespace stormm
