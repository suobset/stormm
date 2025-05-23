#include "copyright.h"
#include "Constants/fixed_precision.h"
#include "Reporting/error_format.h"
#include "coordinate_util.h"

namespace stormm {
namespace trajectory {

using numerics::globalpos_scale_nonoverflow_bits;
using numerics::velocity_scale_nonoverflow_bits;
using numerics::force_scale_nonoverflow_bits;
  
//-------------------------------------------------------------------------------------------------
std::string nameCoordinateType(const size_t ct_coords) {
  if (ct_coords == cf_type_index) {
    return std::string("CoordinateFrame");
  }
  else if (ct_coords == ps_type_index) {
    return std::string("PhaseSpace");
  }
  else if (ct_coords == cs_dtype_index) {
    return std::string("CoordinateSeries<float64_t>");
  }  
  else if (ct_coords == cs_ftype_index) {
    return std::string("CoordinateSeries<float32_t>");
  }  
  else if (ct_coords == cs_stype_index) {
    return std::string("CoordinateSeries<int16_t>");
  }  
  else if (ct_coords == cs_itype_index) {
    return std::string("CoordinateSeries<int32_t>");
  }  
  else if (ct_coords == cs_ltype_index) {
    return std::string("CoordinateSeries<int64_t>");
  }  
  else if (ct_coords == poly_ps_type_index) {
    return std::string("PhaseSpaceSynthesis");
  }  
  else if (ct_coords == cdns_type_index) {
    return std::string("Condensate");
  }
  else {
    return std::string("Unrecognized coordinate object type.", "nameCoordinateType");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
bool overflowRequired(const TrajectoryKind kind, const int scaling_bits, const int* x_ovrf,
                      const int* y_ovrf, const int* z_ovrf) {
  bool result;
  switch (kind)	{
  case TrajectoryKind::POSITIONS:
    result = (scaling_bits > globalpos_scale_nonoverflow_bits);
    break;
  case TrajectoryKind::VELOCITIES:
    result = (scaling_bits > velocity_scale_nonoverflow_bits);
    break;
  case TrajectoryKind::FORCES:
    result = (scaling_bits > force_scale_nonoverflow_bits);
    break;
  }
  if (result && (x_ovrf == nullptr || y_ovrf == nullptr || z_ovrf == nullptr)) {
    rtErr("Overflow bits may be required for " + getEnumerationName(kind) + " a fixed-precision "
          "representation with " + std::to_string(scaling_bits) + " after the decimal, but "
          "arrays of overflow bits were not presented.", "overflowRequired");
  }
  return result;
}
  
} // namespace trajectory
} // namespace stormm
