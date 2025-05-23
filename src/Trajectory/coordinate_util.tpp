// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace trajectory {

//-------------------------------------------------------------------------------------------------
template <typename Tdest, typename Torig>
void checkCopyValidity(const Tdest *destination, const Torig &origin,
                       const HybridTargetLevel destination_tier,
                       const HybridTargetLevel origin_tier) {
  const std::string dest_obj = nameCoordinateType(std::type_index(typeid(Tdest)).hash_code());
  const std::string orig_obj = nameCoordinateType(std::type_index(typeid(Torig)).hash_code());
  switch (destination_tier) {
  case HybridTargetLevel::HOST:
    confirmCpuMemory(destination->getFormat(), "The destination object (" + dest_obj +
                     ") does not have memory allocated on the CPU host to receive coordinates.",
                     "coordCopy");
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    confirmGpuMemory(destination->getFormat(), "The destination object (" + dest_obj +
                     ") does not have memory allocated on the GPU device to receive coordinates.",
                     "coordCopy");
    break;
#endif
  }
  switch (origin_tier) {
  case HybridTargetLevel::HOST:
    confirmCpuMemory(origin.getFormat(), "The origin object (" + orig_obj + ") does not have "
                     "memory allocated on the CPU host to provide coordinates.", "coordCopy");
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    confirmGpuMemory(origin.getFormat(), "The origin object (" + orig_obj + ") does not have "
                     "memory allocated on the GPU device to provide coordinates.", "coordCopy");
    break;
#endif
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T> UnitCellType determineUnitCellType(const T* dims) {
  if (fabs(dims[3]) < 1.0e-6 || fabs(dims[4]) < 1.0e-6 || fabs(dims[5]) < 1.0e-6) {
    return UnitCellType::NONE;
  }
  else if (fabs(dims[3] - (0.5 * symbols::pi)) < 1.0e-6 ||
           fabs(dims[4] - (0.5 * symbols::pi)) < 1.0e-6 ||
           fabs(dims[5] - (0.5 * symbols::pi)) < 1.0e-6) {
    return UnitCellType::ORTHORHOMBIC;
  }
  else {
   return UnitCellType::TRICLINIC;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename T>
UnitCellType determineUnitCellType(const std::vector<T> &dims, const int offset) {
  if (offset + 6 > dims.size()) {
    rtErr("An offset of " + std::to_string(offset) + " is invalid for obtaining box dimensions "
          "from a Standard Template Library vector of size " + std::to_string(dims.size()) + ".",
          "determineUnitCellType");
  }
  return determineUnitCellType(&dims.data()[offset]);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
UnitCellType determineUnitCellType(const Hybrid<T> &dims, const int offset) {
  if (offset + 6 > dims.size()) {
    rtErr("An offset of " + std::to_string(offset) + " is invalid for obtaining box dimensions "
          "from a Hybrid object of size " + std::to_string(dims.size()) + ".",
          "determineUnitCellType");
  }
  return determineUnitCellType(&dims.data()[offset]);
}

} // namespace trajectory
} // namespace stormm
