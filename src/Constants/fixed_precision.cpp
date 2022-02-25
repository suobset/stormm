#include "Reporting/error_format.h"
#include "fixed_precision.h"

namespace omni {
namespace numerics {

//-------------------------------------------------------------------------------------------------
std::string fixedPrecisionRangeErrorMessage(const int choice, const int min_val,
                                            const int max_val) {
  if (choice < min_val) {
    return std::string(" is too small, and would result in an inaccurate simulation.");
  }
  if (choice > max_val) {
    return std::string(" is too large, and might overflow the 64-bit integer format.");
  }
  return std::string("");
}

//-------------------------------------------------------------------------------------------------
void checkGlobalPositionBits(const int choice) {
  if (choice < min_globalpos_scale_bits || choice > max_globalpos_scale_bits) {
    rtErr("Global position scaling is allowed in a fixed precision range of " +
          std::to_string(min_globalpos_scale_bits) + " to " +
          std::to_string(max_globalpos_scale_bits) + ".  A value of " +
          std::to_string(choice) +
          fixedPrecisionRangeErrorMessage(choice, min_globalpos_scale_bits,
                                          max_globalpos_scale_bits), "checkGlobalPositionBits");
  }
}

//-------------------------------------------------------------------------------------------------
void checkLocalPositionBits(const int choice) {
  if (choice < min_localpos_scale_bits || choice > max_localpos_scale_bits) {
    rtErr("Local position scaling is allowed in a fixed precision range of " +
          std::to_string(min_localpos_scale_bits) + " to " +
          std::to_string(max_localpos_scale_bits) + " bits.  A value of " +
          std::to_string(choice) + fixedPrecisionRangeErrorMessage(choice, min_localpos_scale_bits,
                                                                   max_localpos_scale_bits),
          "checkLocalPositionBits");
  }
}

//-------------------------------------------------------------------------------------------------
void checkVelocityBits(const int choice) {
  if (choice < min_velocity_scale_bits || choice > max_velocity_scale_bits) {
    rtErr("Velocity scaling is allowed in a fixed precision range of " +
          std::to_string(min_velocity_scale_bits) + " to " +
          std::to_string(max_velocity_scale_bits) + " bits.  A value of " +
          std::to_string(choice) + fixedPrecisionRangeErrorMessage(choice, min_velocity_scale_bits,
                                                                   max_velocity_scale_bits),
          "checkVelocityBits");
  }
}

//-------------------------------------------------------------------------------------------------
void checkForceBits(const int choice) {
  if (choice < min_force_scale_bits || choice > max_force_scale_bits) {
    rtErr("Force accumulation is allowed in a fixed precision range of " +
          std::to_string(min_force_scale_bits) + " to " +
          std::to_string(max_force_scale_bits) + " bits.  A value of " + std::to_string(choice) +
          fixedPrecisionRangeErrorMessage(choice, min_force_scale_bits, max_force_scale_bits),
          "checkForceBits");
  }
}

} // namespace numerics
} // namespace omni
