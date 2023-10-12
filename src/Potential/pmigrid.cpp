#ifdef STORMM_USE_HPC
#  ifdef STORMM_USE_CUDA
#include <cuda_runtime.h>
#  endif
#endif
#include "copyright.h"
#include "Math/vector_ops.h"
#include "Numerics/split_fixed_precision.h"
#include "Synthesis/brickwork.h"
#include "pmigrid.h"

namespace stormm {
namespace energy {

using numerics::hostInt95ToDouble;
using numerics::hostDoubleToInt95;
using numerics::hostSplitFPSum;
using synthesis::Brickwork;
using stmath::prefixSumInPlace;
using stmath::PrefixSumType;

//-------------------------------------------------------------------------------------------------
PMIGridWriter::PMIGridWriter(const NonbondedTheme theme_in, const PrecisionModel mode_in,
                             const int fp_bits_in, const int nsys_in, const int order_in,
                             const int wu_count_in, const uint4* dims_in, double* ddata_in,
                             float* fdata_in, const uint* work_units_in) :
    theme{theme_in}, mode{mode_in}, shacc_fp_scale{static_cast<float>(pow(2.0, fp_bits_in))},
    nsys{nsys_in}, order{order_in}, wu_count{wu_count_in}, dims{dims_in}, ddata{ddata_in},
    fdata{fdata_in}, work_units{work_units_in}
{}

//-------------------------------------------------------------------------------------------------
PMIGridReader::PMIGridReader(const NonbondedTheme theme_in, const PrecisionModel mode_in,
                             const int fp_bits_in, const int nsys_in, const int order_in,
                             const uint4* dims_in, const double* ddata_in, const float* fdata_in) :
    theme{theme_in}, mode{mode_in}, shacc_fp_scale{static_cast<float>(pow(2.0, fp_bits_in))},
    nsys{nsys_in}, order{order_in}, dims{dims_in}, ddata{ddata_in}, fdata{fdata_in}
{}

//-------------------------------------------------------------------------------------------------
PMIGridReader::PMIGridReader(const PMIGridWriter &w) :
    theme{w.theme}, mode{w.mode}, shacc_fp_scale{w.shacc_fp_scale}, nsys{w.nsys}, order{w.order},
    dims{w.dims}, ddata{w.ddata}, fdata{w.fdata}
{}

//-------------------------------------------------------------------------------------------------
PMIGridReader::PMIGridReader(const PMIGridWriter *w) :
    theme{w->theme}, mode{w->mode}, shacc_fp_scale{w->shacc_fp_scale}, nsys{w->nsys},
    order{w->order}, dims{w->dims}, ddata{w->ddata}, fdata{w->fdata}
{}

//-------------------------------------------------------------------------------------------------
PMIGridAccumulator::PMIGridAccumulator(const NonbondedTheme theme_in, const PrecisionModel mode_in,
                                       const int fp_bits_in, const int nsys_in, const int order_in,
                                       const int wu_count_in, const uint4* dims_in,
                                       double* ddata_in, float* fdata_in, int* overflow_in,
                                       const uint* work_units_in) :
    theme{theme_in}, mode{mode_in}, fp_bits{fp_bits_in},
    fp_scale{static_cast<float>(pow(2.0, fp_bits_in))},
    nsys{nsys_in}, order{order_in},
    order_squared{order_in * order_in},
    order_cubed{order_in * order_in * order_in},
    wu_count{wu_count_in},
    dims{dims_in},
    lldata{reinterpret_cast<llint*>(ddata_in)},
    idata{reinterpret_cast<int*>(fdata_in)},
    overflow{overflow_in},
    work_units{work_units_in}
{}

//-------------------------------------------------------------------------------------------------
PMIGridFPReader::PMIGridFPReader(const NonbondedTheme theme_in, const PrecisionModel mode_in,
                                 const int fp_bits_in, const int nsys_in, const int order_in,
                                 const uint4* dims_in, const double* ddata_in,
                                 const float* fdata_in, const int* overflow_in) :
    theme{theme_in}, mode{mode_in}, fp_bits{fp_bits_in},
    fp_scale{static_cast<float>(pow(2.0, fp_bits_in))},
    nsys{nsys_in}, order{order_in}, dims{dims_in},
    lldata{reinterpret_cast<const llint*>(ddata_in)},
    idata{reinterpret_cast<const int*>(fdata_in)},
    overflow{overflow_in}
{}

//-------------------------------------------------------------------------------------------------
PMIGridFPReader::PMIGridFPReader(const PMIGridAccumulator &w) :
    theme{w.theme}, mode{w.mode}, fp_bits{w.fp_bits}, fp_scale{w.fp_scale}, nsys{w.nsys},
    order{w.order}, dims{w.dims}, lldata{w.lldata}, idata{w.idata}, overflow{w.overflow}
{}

//-------------------------------------------------------------------------------------------------
PMIGridFPReader::PMIGridFPReader(const PMIGridAccumulator *w) :
    theme{w->theme}, mode{w->mode}, fp_bits{w->fp_bits}, fp_scale{w->fp_scale}, nsys{w->nsys},
    order{w->order}, dims{w->dims}, lldata{w->lldata}, idata{w->idata}, overflow{w->overflow}
{}

//-------------------------------------------------------------------------------------------------
NonbondedTheme PMIGrid::getTheme() const {
  return theme;
}

//-------------------------------------------------------------------------------------------------
PrecisionModel PMIGrid::getMode() const {
  return mode;
}

//-------------------------------------------------------------------------------------------------
bool PMIGrid::fixedPrecisionEnabled() const {
  return (fp_accumulation_bits > 0);
}

//-------------------------------------------------------------------------------------------------
QMapMethod PMIGrid::getRecommendedMappingMethod() const {
  return recommendation;
}

//-------------------------------------------------------------------------------------------------
QMapMethod PMIGrid::getWorkUnitConfiguration() const {
  return work_unit_configuration;
}

//-------------------------------------------------------------------------------------------------
bool PMIGrid::dataIsReal() const {
  return data_is_real;
}

//-------------------------------------------------------------------------------------------------
int PMIGrid::getFixedPrecisionBits() const {
  return fp_accumulation_bits;
}

//-------------------------------------------------------------------------------------------------
int PMIGrid::getSharedFixedPrecisionBits() const {
  return shared_fp_accumulation_bits;
}

//-------------------------------------------------------------------------------------------------
int PMIGrid::getSystemCount() const {
  return system_count;
}

//-------------------------------------------------------------------------------------------------
int PMIGrid::getInterpolationOrder() const {
  return b_spline_order;
}

//-------------------------------------------------------------------------------------------------
size_t PMIGrid::getTotalCapacity() const {
  return capacity;
}

//-------------------------------------------------------------------------------------------------
uint4 PMIGrid::getGridDimensions(const int system_index) const {
  return grid_dimensions.readHost(system_index);
}

//-------------------------------------------------------------------------------------------------
int PMIGrid::getGridDimensions(const int system_index, const UnitCellAxis uc_axis) const {
  const uint4 gdim = grid_dimensions.readHost(system_index);
  switch (uc_axis) {
  case UnitCellAxis::A:
    return gdim.x;
  case UnitCellAxis::B:
    return gdim.y;
  case UnitCellAxis::C:
    return gdim.z;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
int PMIGrid::getGridDimensions(const int system_index, const CartesianDimension uc_axis) const {
  const uint4 gdim = grid_dimensions.readHost(system_index);
  switch (uc_axis) {
  case CartesianDimension::X:
    return gdim.x;
  case CartesianDimension::Y:
    return gdim.y;
  case CartesianDimension::Z:
    return gdim.z;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
int PMIGrid::getWorkUnitCount() const {
  return work_unit_count;
}

//-------------------------------------------------------------------------------------------------
PMIGridWriter PMIGrid::data(const HybridTargetLevel tier) {
  return PMIGridWriter(theme, mode, shared_fp_accumulation_bits, system_count, b_spline_order,
                       work_unit_count, grid_dimensions.data(tier), dgrid_stack.data(tier),
                       fgrid_stack.data(tier), work_units.data(tier));
}

//-------------------------------------------------------------------------------------------------
const PMIGridReader PMIGrid::data(const HybridTargetLevel tier) const {

  // Raise an exception if a real-valued reader is requested for fixed-precision integer data.
  if (data_is_real == false) {
    rtErr("Data is currently represented as fixed-precision integer values with " +
          std::to_string(fp_accumulation_bits) + " bits of detail after the decimal.  Convert to "
          "real-valued data before trying to read it as floating-point values.", "PMIGrid",
          "data");
  }
  return PMIGridReader(theme, mode, shared_fp_accumulation_bits, system_count, b_spline_order,
                       grid_dimensions.data(tier), dgrid_stack.data(tier), fgrid_stack.data(tier));

}

//-------------------------------------------------------------------------------------------------
PMIGridAccumulator PMIGrid::fpData(const HybridTargetLevel tier) {

  // Raise an exception if a fixed-precision writer is requested for an object that is not yet
  // prepared to support such accumulation.
  if (overflow_stack.size() != capacity) {
    rtErr("Overflow accumulators must be allocated in order to accumulate density in "
          "fixed-precision.", "PMIGrid", "fpData");
  }
  return PMIGridAccumulator(theme, mode, fp_accumulation_bits, system_count, b_spline_order,
                            work_unit_count, grid_dimensions.data(tier), dgrid_stack.data(tier),
                            fgrid_stack.data(tier), overflow_stack.data(tier),
                            work_units.data(tier));
}

//-------------------------------------------------------------------------------------------------
const PMIGridFPReader PMIGrid::fpData(const HybridTargetLevel tier) const {

  // Raise an exception if a real-valued reader is requested for fixed-precision integer data.
  if (data_is_real) {
    rtErr("Data is currently represented as real-valued, floating point numbers.  Re-initialize "
          "and accumulate the data as fixed-precision integers before attempting to view it in "
          "such a format.  Current fixed-precision detail bit count: " +
          std::to_string(fp_accumulation_bits) + ".", "PMIGrid", "fpData");
  }
  return PMIGridFPReader(theme, mode, fp_accumulation_bits, system_count, b_spline_order,
                         grid_dimensions.data(tier), dgrid_stack.data(tier),
                         fgrid_stack.data(tier), overflow_stack.data(tier));
}

//-------------------------------------------------------------------------------------------------
const CellGridReader<void, void, void, void>
PMIGrid::getTemplateFreeCellGridReader(const HybridTargetLevel tier) const {
  if (cg_tmat == int_type_index) {
    return unrollTemplateFreeCGReader<int, int4>(tier);
  }
  else if (cg_tmat == llint_type_index) {
    return unrollTemplateFreeCGReader<llint, llint4>(tier);    
  }
  else if (cg_tmat == float_type_index) {
    return unrollTemplateFreeCGReader<float, float4>(tier);
  }
  else if (cg_tmat == double_type_index) {
    return unrollTemplateFreeCGReader<double, double4>(tier);
  }
  else {
    rtErr("The valid types for the CellGrid's coordinate representation are int, llint, float, "
          "and double.", "PMIGrid", "getTemplateFreeCellGridReader");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
size_t PMIGrid::getCellGridMatrixTypeID() const {
  return cg_tmat;
}

//-------------------------------------------------------------------------------------------------
size_t PMIGrid::getCellGridAccumulatorTypeID() const {
  return cg_tacc;
}

//-------------------------------------------------------------------------------------------------
size_t PMIGrid::getCellGridCalculationTypeID() const {
  return cg_tcalc;
}

//-------------------------------------------------------------------------------------------------
size_t PMIGrid::getCellGridCoordinateTypeID() const {
  return cg_tcrd;
}

//-------------------------------------------------------------------------------------------------
double PMIGrid::getTotalOnGrid(const int system_index) const {
  validateSystemIndex(system_index);
  const uint4 gdims = grid_dimensions.readHost(system_index);
  const size_t illim = gdims.w;
  const size_t ihlim = illim + static_cast<size_t>(gdims.x * gdims.y * gdims.z);
  double result = 0.0;
  if (data_is_real) {
    const double* ddata_ptr = dgrid_stack.data();
    const float* fdata_ptr = fgrid_stack.data();
    const double scale_up = pow(2.0, 72.0);
    const double scale_down = 1.0 / scale_up;
    int95_t iresult = { 0LL, 0 };
    switch (mode) {
    case PrecisionModel::DOUBLE:
      for (size_t i = illim; i < ihlim; i++) {
        const int95_t i_val = hostDoubleToInt95(ddata_ptr[i] * scale_up);
        iresult = hostSplitFPSum(i_val, iresult);
      }
      break;
    case PrecisionModel::SINGLE:
      for (size_t i = illim; i < ihlim; i++) {
        const int95_t i_val = hostDoubleToInt95(static_cast<double>(fdata_ptr[i]) * scale_up);
        iresult = hostSplitFPSum(i_val, iresult);
      }
      break;
    }
    result = hostInt95ToDouble(iresult) * scale_down;
  }
  else {
    const double inv_scl = pow(2.0, -fp_accumulation_bits);
    const llint* lldata_ptr = reinterpret_cast<const llint*>(dgrid_stack.data());
    const int* idata_ptr = reinterpret_cast<const int*>(fgrid_stack.data());
    const int* ovrf_ptr = overflow_stack.data();
    switch (mode) {
    case PrecisionModel::DOUBLE:
      {
        int95_t iresult = { 0LL, 0 };
        for (size_t i = illim; i < ihlim; i++) {
          iresult = hostSplitFPSum(iresult, lldata_ptr[i], ovrf_ptr[i]);
        }
        result = hostInt95ToDouble(iresult) * inv_scl;
      }
      break;
    case PrecisionModel::SINGLE:
      {
        int2 iresult = { 0, 0 };
        for (size_t i = illim; i < ihlim; i++) {
          iresult = hostSplitFPSum(iresult, idata_ptr[i], ovrf_ptr[i]);
        }
        result = hostInt63ToDouble(iresult) * inv_scl;
      }
      break;
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<double> PMIGrid::getGrid(const int system_index) const {
  validateSystemIndex(system_index);
  const uint4 gdims = grid_dimensions.readHost(system_index);
  const size_t llim = gdims.w;
  const size_t hlim = gdims.w + (gdims.x * gdims.y * gdims.z);
  std::vector<double> result(hlim - llim);
  switch (mode) {
  case PrecisionModel::DOUBLE:
    if (data_is_real) {
      const double* ddata_ptr = dgrid_stack.data();
      for (size_t i = llim; i < hlim; i++) {
        result[i - llim] = ddata_ptr[i];
      }
    }
    else {
      const PMIGridFPReader fpr = fpData();
      const double inv_scale = 1.0 / fpr.fp_scale;
      for (size_t i = llim; i < hlim; i++) {
        result[i - llim] = hostInt95ToDouble(fpr.lldata[i], fpr.overflow[i]) * inv_scale;
      }
    }
    break;
  case PrecisionModel::SINGLE:
    if (data_is_real) {
      const float* fdata_ptr = fgrid_stack.data();
      for (size_t i = llim; i < hlim; i++) {
        result[i - llim] = fdata_ptr[i];
      }
    }
    else {
      const PMIGridFPReader fpr = fpData();
      const double inv_scale = 1.0 / fpr.fp_scale;
      for (size_t i = llim; i < hlim; i++) {
        result[i - llim] = hostInt63ToDouble(fpr.idata[i], fpr.overflow[i]) * inv_scale;
      }
    }
    break;
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
const PhaseSpaceSynthesis* PMIGrid::getCoordinateSynthesisPointer() const {
  return poly_ps_pointer;
}

//-------------------------------------------------------------------------------------------------
const PMIGrid* PMIGrid::getSelfPointer() const {
  return this;
}

#ifdef STORMM_USE_HPC
//-------------------------------------------------------------------------------------------------
void PMIGrid::upload() {
  grid_dimensions.upload();
  dgrid_stack.upload();
  fgrid_stack.upload();
  overflow_stack.upload();
  work_units.upload();
}

//-------------------------------------------------------------------------------------------------
void PMIGrid::download() {
  grid_dimensions.download();
  dgrid_stack.download();
  fgrid_stack.download();
  overflow_stack.download();
  work_units.download();
}
#endif

//-------------------------------------------------------------------------------------------------
void PMIGrid::setMode(const PrecisionModel mode_in) {

  // If the mode changes, free one of the Hybrid data arrays and allocate the other as needed.
  // Set pointers to reinterpret (re-use) the newly allocated data array for fixed-precision
  // integer accumulation as appropriate.
  switch (mode_in) {
  case PrecisionModel::DOUBLE:
    fgrid_stack.resize(0);
    fgrid_stack.shrinkToFit();
    dgrid_stack.resize(capacity);
    break;
  case PrecisionModel::SINGLE:
    dgrid_stack.resize(0);
    dgrid_stack.shrinkToFit();
    fgrid_stack.resize(capacity);
    break;
  }
  mode = mode_in;    
}

//-------------------------------------------------------------------------------------------------
void PMIGrid::setRealDataFormat(const bool real_in) {
  data_is_real = real_in;
}

//-------------------------------------------------------------------------------------------------
void PMIGrid::prepareFixedPrecisionModel(const int fp_accumulation_bits_in,
                                         const int shared_fp_bits_in) {
  if (fp_accumulation_bits_in < 0) {
    rtErr("A fixed precision accumulation in " + std::to_string(fp_accumulation_bits_in) +
          " is invalid.", "PMIGrid", "prepareFixedPrecisionModel");
  }
  fp_accumulation_bits = fp_accumulation_bits_in;
  data_is_real = (fp_accumulation_bits == 0);
  if (fp_accumulation_bits > 0) {
    validateFixedPrecisionBits(fp_accumulation_bits);
    if (overflow_stack.size() == 0) {
      overflow_stack.resize(capacity);
    }
    shared_fp_accumulation_bits = fp_accumulation_bits;
  }
  else {
    if (overflow_stack.size() > 0) {
      overflow_stack.resize(0);
      overflow_stack.shrinkToFit();
    }

    // Set the __shared__ accumulation fixed precision bits according to user input, or default
    // values if necessary.  There must always be a significant fixed precision model for
    // __shared__ accumulation in case one of the kernels is called.
    if (shared_fp_bits_in > 0) {
      shared_fp_accumulation_bits = shared_fp_bits_in;
      validateFixedPrecisionBits(shared_fp_accumulation_bits);
    }
    else {
      switch (theme) {
      case NonbondedTheme::ELECTROSTATIC:
        switch (mode) {
        case PrecisionModel::DOUBLE:
          shared_fp_accumulation_bits = default_qspread_fp_bits_dp;
          break;
        case PrecisionModel::SINGLE:
          shared_fp_accumulation_bits = default_qspread_fp_bits_sp;
          break;
        }
        break;
      case NonbondedTheme::VAN_DER_WAALS:
        switch (mode) {
        case PrecisionModel::DOUBLE:
          shared_fp_accumulation_bits = default_ljspread_fp_bits_dp;
          break;
        case PrecisionModel::SINGLE:
          shared_fp_accumulation_bits = default_ljspread_fp_bits_dp;
          break;
        }
        break;
      case NonbondedTheme::ALL:
        break;
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------
void PMIGrid::setRecommendedMappingMethod(const GpuDetails &gpu) {
  const CellGridReader<void, void, void, void> cgr_vc = getTemplateFreeCellGridReader();
  if (b_spline_order <= 6 && b_spline_order <= cgr_vc.mesh_ticks + 1) {
    recommendation = QMapMethod::ACC_REGISTER;
  }
  else {

    // Check that the grid points of at least four spatial decomposition cells will fit in the
    // allotted __shared__memory space.
    const int cell_gp = cgr_vc.mesh_ticks * cgr_vc.mesh_ticks * cgr_vc.mesh_ticks;
    switch (mode) {
    case PrecisionModel::DOUBLE:
      if (shared_acc_buffer_size / (cell_gp * (sizeof(int) + sizeof(llint))) >= 4) {
        recommendation = QMapMethod::ACC_SHARED;
      }
      else {
        recommendation = QMapMethod::GENERAL_PURPOSE;
      }
      break;
    case PrecisionModel::SINGLE:
      if (shared_acc_buffer_size / (cell_gp * (2 * sizeof(int))) >= 4) {
        recommendation = QMapMethod::ACC_SHARED;
      }
      else {
        recommendation = QMapMethod::GENERAL_PURPOSE;
      }
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
void PMIGrid::prepareWorkUnits(const QMapMethod approach, const GpuDetails &gpu) {
  shared_acc_buffer_size = findSharedBufferSize(gpu);

  // Obtain a read-only, void-casted abstract of the attached CellGrid.
  const CellGridReader<void, void, void, void> cgr_vc = getTemplateFreeCellGridReader();

  // Compute the brick layout for the associated CellGrid object.
  std::vector<int3> ext_dims(system_count);
  for (int i = 0; i < system_count; i++) {
    const ullint iv = cgr_vc.system_cell_grids[i];
    ext_dims[i] = { static_cast<int>((iv >> 28) & 0xfff), static_cast<int>((iv >> 40) & 0xfff),
                    static_cast<int>(iv >> 52) };
  }

  // Mark the work unit configuration
  work_unit_configuration = approach;
  Brickwork bw;
  int halo;
  switch (approach) {
  case QMapMethod::GENERAL_PURPOSE:
    return;
  case QMapMethod::AUTOMATIC:
    setRecommendedMappingMethod(gpu);
    prepareWorkUnits(recommendation, gpu);
    return;
  case QMapMethod::ACC_REGISTER:
    {
      // Check that other parameters are suitable for the register accumulation kernel.  If not,
      // recursively call the work unit preparation with a feasible method.
      if (b_spline_order >= 6 || b_spline_order > cgr_vc.mesh_ticks + 1) {
        recommendation = QMapMethod::GENERAL_PURPOSE;
        prepareWorkUnits(recommendation, gpu);
        return;
      }
      switch (b_spline_order) {
      case 4:
        bw = Brickwork(ext_dims, max4s_atom_bearing_region_adim, max4s_atom_bearing_cross_section,
                       1, 0, max4s_grid_mapping_volume, gpu.getSMPCount());
        break;
      case 5:
        bw = Brickwork(ext_dims, max5s_atom_bearing_region_adim, max5s_atom_bearing_cross_section,
                       1, 0, max5s_grid_mapping_volume, gpu.getSMPCount());
        break;
      case 6:
        bw = Brickwork(ext_dims, max6s_atom_bearing_region_adim, max6s_atom_bearing_cross_section,
                       1, 0, max6s_grid_mapping_volume, gpu.getSMPCount());
        break;
      }
      halo = 1;
    }
    break;
  case QMapMethod::ACC_SHARED:
    {
      // Compute the largest cross-sectional area based on the allowed size of the __shared__
      // memory buffer, the precision of the accumulation, and the density of the grid.  The
      // maximum cross section is first computed in terms of the grid-mapping region and then
      // updated with the necessary density mapping halo. 
      const int cell_gp = cgr_vc.mesh_ticks * cgr_vc.mesh_ticks * cgr_vc.mesh_ticks;
      int max_buffered_cells;
      switch (mode) {
      case PrecisionModel::DOUBLE:
        max_buffered_cells = shared_acc_buffer_size / (cell_gp * (sizeof(llint) + sizeof(int)));
        break;
      case PrecisionModel::SINGLE:
        max_buffered_cells = shared_acc_buffer_size / (cell_gp * (2 * sizeof(int)));
        break;
      }
      halo = computeMappingHalo();
      
      // If the volume of a single cell would exceed the allowed limits of the __shared__ memory
      // per block, recursively call the function to use the general-purpose mapping kernel.
      if (max_buffered_cells == 0) {
        setRecommendedMappingMethod(gpu);
        prepareWorkUnits(QMapMethod::GENERAL_PURPOSE, gpu);
        return;
      }
      int max_xsection = std::min((halo + 3) * (halo + 3),
                                  (max_buffered_cells + halo) * (halo + 1));

      // If the maximum cross section permits no cells, a new mapping method must be determined.
      if (max_xsection == 0) {
        recommendation = QMapMethod::GENERAL_PURPOSE;
        prepareWorkUnits(recommendation, gpu);
        return;
      }
      
      // Compute the brickwork for this set of systems.
      bw = Brickwork(ext_dims, max_shared_acc_atom_bearing_region_adim, max_xsection, halo, 0,
                     max_buffered_cells, gpu.getSMPCount());
    }
    break;
  }

  // Cases that do not invlove work unit production will have returned before reaching this point.
  work_unit_count = bw.getBrickCount();
  work_units.resize(density_mapping_wu_size * work_unit_count);
  std::vector<uint> tmp_wu;
  for (int i = 0; i < work_unit_count; i++) {
    const int3 i_orig = bw.getBrickOrigin(i);
    const int3 i_dims = bw.getBrickLengths(i);
    addWorkUnit(&tmp_wu, bw.getSystemMembership(i), i_orig.x, i_orig.y, i_orig.z, i_dims.x,
                i_dims.y, i_dims.z, cgr_vc, halo);
  }
  work_units.putHost(tmp_wu);
}

//-------------------------------------------------------------------------------------------------
void PMIGrid::prepareWorkUnits(const GpuDetails &gpu) {
  
  // Force evaluation of the recommended mapping method if it has not already been done
  if (recommendation == QMapMethod::AUTOMATIC) {
    setRecommendedMappingMethod(gpu);
  }
  prepareWorkUnits(recommendation, gpu);
}

//-------------------------------------------------------------------------------------------------
int PMIGrid::getLargestWorkUnitGridPoints() const {
  const CellGridReader<void, void, void, void> cgr_vc = getTemplateFreeCellGridReader();

  // The largest work unit will be at the front of the list, per ordering in prepareWorkUnits().
  // However, the size of the work unit refers to the atom-bearing region, which may not guarantee
  // that the work unit with the largest volume of atoms has the largest volume of grid points to
  // map.
  const int gp_per_cell = cgr_vc.mesh_ticks * cgr_vc.mesh_ticks * cgr_vc.mesh_ticks;
  const uint* wu_ptr = work_units.data();
  uint largest_gm_region = 0;
  switch (work_unit_configuration) {
  case QMapMethod::GENERAL_PURPOSE:
  case QMapMethod::AUTOMATIC:
    return 0;
  case QMapMethod::ACC_REGISTER:
    for (int i = 0; i < work_unit_count; i++) {
      const int ofs = i * density_mapping_wu_size;
      largest_gm_region = std::max(largest_gm_region,
                                   (wu_ptr[ofs + 20] - 1) * (wu_ptr[ofs + 21] - 1) *
                                   (wu_ptr[ofs + 22] - 1));
    }
    break;
  case QMapMethod::ACC_SHARED:
    for (int i = 0; i < work_unit_count; i++) {
      const int ofs = i * density_mapping_wu_size;
      largest_gm_region = std::max(largest_gm_region,
                                   wu_ptr[ofs + 10] * wu_ptr[ofs + 11] * wu_ptr[ofs + 12]);
    }
    break;
  }
  return largest_gm_region * gp_per_cell;    
}

//-------------------------------------------------------------------------------------------------
int PMIGrid::computeMappingHalo() const {
  const CellGridReader<void, void, void, void> cgr_vc = getTemplateFreeCellGridReader();
  int result = 1;
  int cg_ext = cgr_vc.mesh_ticks + 1;
  while (cg_ext < b_spline_order) {
    cg_ext += cgr_vc.mesh_ticks;
    result++;
  }
  return result;
}
    
//-------------------------------------------------------------------------------------------------
void PMIGrid::initialize(const HybridTargetLevel tier, const GpuDetails &gpu) {
  switch (tier) {
  case HybridTargetLevel::HOST:
    if (fp_accumulation_bits > 0) {
      PMIGridAccumulator pm_acc = fpData(tier);
      for (int pos = 0; pos < pm_acc.nsys; pos++) {
        const uint4 pdims = pm_acc.dims[pos];
        const uint ilim = pdims.w + (pdims.x * pdims.y * pdims.z);
        switch (mode) {
        case PrecisionModel::DOUBLE:
          for (uint i = pdims.w; i < ilim; i++) {
            pm_acc.lldata[i] = 0LL;
            pm_acc.overflow[i] = 0;
          }
          break;
        case PrecisionModel::SINGLE:
          for (uint i = pdims.w; i < ilim; i++) {
            pm_acc.idata[i] = 0LL;
            pm_acc.overflow[i] = 0;
          }
          break;
        }
      }
    }
    else {
      PMIGridWriter pm_wrt = data(tier);
      for (int pos = 0; pos < pm_wrt.nsys; pos++) {
        const uint4 pdims = pm_wrt.dims[pos];
        const uint ilim = pdims.w + (pdims.x * pdims.y * pdims.z);
        switch (mode) {
        case PrecisionModel::DOUBLE:
          for (uint i = pdims.w; i < ilim; i++) {
            pm_wrt.ddata[i] = 0.0;
          }
          break;
        case PrecisionModel::SINGLE:
          for (uint i = pdims.w; i < ilim; i++) {
            pm_wrt.fdata[i] = 0.0f;
          }
          break;
        }
      }
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    if (fp_accumulation_bits > 0) {
      PMIGridAccumulator pm_acc = fpData(tier);
      launchPMIGridInitialization(&pm_acc, gpu);
    }
    break;
#endif
  }

  // Note whether data has been initialized for a fixed-precision or real-valued representation.
  data_is_real = (fp_accumulation_bits == 0);
}

//-------------------------------------------------------------------------------------------------
void PMIGrid::convertToReal(const HybridTargetLevel tier, const GpuDetails &gpu) {

  // Return immediately if this object is already working in real-valued numbers.
  if (fp_accumulation_bits == 0) {
    return;
  }

  // Raise an exception if the data is already real-valued.
  if (data_is_real) {
    rtErr("Data in the object is already converted to real values.", "PMIGrid", "convertToReal");
  }

  // Perform the conversion in host or device memory.
  PMIGridAccumulator pm_acc = fpData(tier);
  PMIGridWriter pm_wrt = data(tier);
  switch (tier) {
  case HybridTargetLevel::HOST:
    {
      const double conv_scale = pow(2.0, -fp_accumulation_bits);
      const float conv_scalef = conv_scale;
      for (int i = 0; i < system_count; i++) {
        const uint4 gdims = pm_acc.dims[i];
        const uint jlim = gdims.w + (gdims.x * gdims.y * gdims.z);
        switch (mode) {
        case PrecisionModel::DOUBLE:
          for (uint j = gdims.w; j < jlim; j++) {
            pm_wrt.ddata[j] = hostInt95ToDouble(pm_acc.lldata[j], pm_acc.overflow[j]) * conv_scale;
          }
          break;
        case PrecisionModel::SINGLE:
          for (uint j = gdims.w; j < jlim; j++) {
            pm_wrt.fdata[j] = hostInt63ToFloat(pm_acc.idata[j], pm_acc.overflow[j]) * conv_scalef;
          }
          break;
        }
      }
    }
    break;
#ifdef STORMM_USE_HPC    
  case HybridTargetLevel::DEVICE:
    {
      launchPMIGridRealConversion(&pm_wrt, pm_acc, gpu);
#  ifdef STORMM_USE_CUDA
      if (cudaDeviceSynchronize() != cudaSuccess) {
        rtErr("Error in device synchronization after launching conversion to real-valued "
              "representation.", "PMIGrid", "convertToReal");
      }
#  endif
    }
    break;
#endif
  }

  // Note that the data is present in real-valued quantities.
  data_is_real = true;
}

//-------------------------------------------------------------------------------------------------
void PMIGrid::validateFixedPrecisionBits(const int fp_bits) const {
  int min_problem = -1;
  int max_problem = -1;
  switch (theme) {
  case NonbondedTheme::ELECTROSTATIC:
    if (fp_bits < minimum_qspread_fp_bits) {
      min_problem = minimum_qspread_fp_bits;
    }
    switch (mode) {
    case PrecisionModel::DOUBLE:
      if (fp_bits > maximum_qspread_fp_bits_dp) {
        max_problem = maximum_qspread_fp_bits_dp;
      }
      break;
    case PrecisionModel::SINGLE:
      if (fp_bits > maximum_qspread_fp_bits_sp) {
        max_problem = maximum_qspread_fp_bits_sp;
      }
      break;
    }
    break;
  case NonbondedTheme::VAN_DER_WAALS:
    if (fp_bits < minimum_ljspread_fp_bits) {
      min_problem = minimum_ljspread_fp_bits;
    }
    switch (mode) {
    case PrecisionModel::DOUBLE:
      if (fp_bits > maximum_ljspread_fp_bits_dp) {
        max_problem = maximum_ljspread_fp_bits_dp;
      }
      break;
    case PrecisionModel::SINGLE:
      if (fp_bits > maximum_ljspread_fp_bits_sp) {
        max_problem = maximum_ljspread_fp_bits_sp;
      }
      break;
    }
    break;
  case NonbondedTheme::ALL:
    break;
  }
  if (min_problem > 0) {
    rtErr("A minimum of " + std::to_string(min_problem) + " bits is required to ensure stable "
          "representation of the density on a " + getEnumerationName(theme) + " particle-mesh "
          "interaction grid.", "PMIGrid", "validateFixedPrecisionBits");
  }
  if (max_problem > 0) {
    rtErr("A maximum of " + std::to_string(max_problem) + " bits is allowed to guard against "
          "overflow of the density representation on a " + getEnumerationName(theme) +
          " particle-mesh interaction grid.", "PMIGrid", "validateFixedPrecisionBits");
  }
}

//-------------------------------------------------------------------------------------------------
void PMIGrid::validateSystemIndex(const int system_index) const {
  if (system_index < 0 || system_index >= poly_ps_pointer->getSystemCount()) {
    rtErr("System index " + std::to_string(system_index) + " is invalid for an associated "
          "synthesis of " + std::to_string(poly_ps_pointer->getSystemCount()) + " systems.",
          "PMIGrid", "validateSystemIndex");
  }
}

//-------------------------------------------------------------------------------------------------
int PMIGrid::findSharedBufferSize(const GpuDetails &gpu) const {
#ifdef STORMM_USE_HPC
#  ifdef STORMM_USE_CUDA
  if (gpu.getArchMajor() >= 7) {
    return 24000;
  }
  else {
    return 15360;
  }
#  else
  return 10800;
#  endif
#else
  return 24000;
#endif
}

//-------------------------------------------------------------------------------------------------
void PMIGrid::addWorkUnit(std::vector<uint> *result, const int sysid, const int pos_a,
                          const int pos_b, const int pos_c, const int gmap_a, const int gmap_b,
                          const int gmap_c, const CellGridReader<void, void, void, void> &cgr_vc,
                          const int halo) {
  const ullint cg_loc = cgr_vc.system_cell_grids[sysid];
  const int system_init_cell = (cg_loc & 0xfffffffLLU);
  const int system_cell_adim = ((cg_loc >> 28) & 0xffLLU);
  const int system_cell_bdim = ((cg_loc >> 40) & 0xffLLU);
  const int system_cell_cdim = ((cg_loc >> 52) & 0xffLLU);
  
  // Tailor the work units to the chosen mapping approach
  switch (work_unit_configuration) {
  case QMapMethod::ACC_REGISTER:
    {
      // Push the new work unit to the back, padding any of the unused 32 indices with zeros.  The
      // register-based accumulation comes with an implicit assumption that all atoms that map to
      // any given grid point can be found in the same spatial decomposition cell, or an adjacent
      // cell down / back / to the left of it.
      int ncs = 0;
      int reim_i = pos_a - 1;
      reim_i += ((reim_i < 0) - (reim_i >= system_cell_adim)) * system_cell_adim;
      for (int j = pos_b - 1; j < pos_b + gmap_b; j++) {
        const int reim_j = j + (((j < 0) - (j >= system_cell_bdim)) * system_cell_bdim);
        for (int k = pos_c - 1; k < pos_c + gmap_c; k++) {
         const int reim_k = k - (((k < 0) - (k >= system_cell_cdim)) * system_cell_cdim);
           const int bc_cell_init = system_init_cell +
                                   (((reim_k * system_cell_bdim) + reim_j) * system_cell_adim);
          result->push_back((bc_cell_init + reim_i) | ((gmap_a + 1) << 28));
          ncs++;
        }
      }
      for (int i = ncs; i < 16; i++) {
        result->push_back(0);
      }
      result->push_back(sysid);
      result->push_back(reim_i);
      result->push_back(pos_b - 1 + ((pos_b - 1 < 0) * system_cell_bdim));
      result->push_back(pos_c - 1 + ((pos_c - 1 < 0) * system_cell_cdim));
      result->push_back(gmap_a + 1);
      result->push_back(gmap_b + 1);
      result->push_back(gmap_c + 1);
      const uint4 pmig_dims = grid_dimensions.readHost(sysid);
      result->push_back(pmig_dims.x);
      result->push_back(pmig_dims.y);
      result->push_back(pmig_dims.z);
      result->push_back(pmig_dims.w);
      int warp_strides_per_cell;  
      warp_strides_per_cell = ((cgr_vc.mesh_ticks * cgr_vc.mesh_ticks * cgr_vc.mesh_ticks) +
                               warp_size_int - 1) / warp_size_int;
      result->push_back((gmap_a * gmap_b * gmap_c * warp_strides_per_cell) |
                        (warp_strides_per_cell << 16));

      // Indices 28-31 are filled in with zeros
      for (int i = 0; i < 4; i++) {
        result->push_back(0);
      }
    }
    break;
  case QMapMethod::ACC_SHARED:
    {
      int reim_i = pos_a - halo;
      reim_i += ((reim_i < 0) - (reim_i >= system_cell_adim)) * system_cell_adim;
      int reim_j = pos_b - halo;
      reim_j += ((reim_j < 0) - (reim_j >= system_cell_bdim)) * system_cell_bdim;
      int reim_k = pos_c - halo;
      reim_k += ((reim_k < 0) - (reim_k >= system_cell_cdim)) * system_cell_cdim;
      result->push_back(sysid);
      result->push_back(reim_i);
      result->push_back(reim_j);
      result->push_back(reim_k);
      result->push_back(gmap_a + halo);
      result->push_back(gmap_b + halo);
      result->push_back(gmap_c + halo);
      result->push_back(pos_a);
      result->push_back(pos_b);
      result->push_back(pos_c);
      result->push_back(gmap_a);
      result->push_back(gmap_b);
      result->push_back(gmap_c);
      result->push_back(system_cell_adim);
      result->push_back(system_cell_bdim);
      result->push_back(system_cell_cdim);
      result->push_back((gmap_b + halo) * (gmap_c + halo));
      result->push_back(cgr_vc.mesh_ticks * cgr_vc.mesh_ticks * cgr_vc.mesh_ticks *
                        gmap_a * gmap_b * gmap_c);
      const uint4 pmig_dims = grid_dimensions.readHost(sysid);
      result->push_back(pmig_dims.x);
      result->push_back(pmig_dims.y);
      result->push_back(pmig_dims.z);
      result->push_back(pmig_dims.w);
      const ullint cg_dims = cgr_vc.system_cell_grids[sysid];
      const uint cell_start_idx = (cg_dims & 0xfffffffLLU);
      const uint cell_na = ((cg_dims >> 28) & 0xfffLLU);
      const uint cell_nb = ((cg_dims >> 40) & 0xfffLLU);
      const uint cell_nc = (cg_dims >> 52);
      result->push_back(cell_na);
      result->push_back(cell_nb);
      result->push_back(cell_nc);
      result->push_back(cell_start_idx);
      result->push_back(cgr_vc.system_chain_bounds[sysid]);

      // Indices 27-31 are filled with zeros
      for (int i = 0; i < 5; i++) {
        result->push_back(0);
      }
    }
    break;
  case QMapMethod::GENERAL_PURPOSE:
  case QMapMethod::AUTOMATIC:

    // These cases will never be reached and imply no work units
    break;
  }
}

} // namespace energy
} // namespace stormm
