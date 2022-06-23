#include "Parsing/parse.h"
#include "Reporting/error_format.h"
#include "fixed_precision.h"

namespace omni {
namespace numerics {

using parse::strcmpCased;
  
//-------------------------------------------------------------------------------------------------
PrecisionLevel translatePrecisionLevel(const std::string &choice, const ExceptionResponse policy) {
  if (strcmpCased(choice, std::string("single"))) {
    return PrecisionLevel::SINGLE;
  }
  else if (strcmpCased(choice, std::string("single_plus"))) {
    return PrecisionLevel::SINGLE_PLUS;
  }
  else if (strcmpCased(choice, std::string("double"))) {
    return PrecisionLevel::DOUBLE;
  }
  else {
    rtErr("Invalid request for precision level " + choice + ".", "translatePrecisionLevel");
  }
  __builtin_unreachable();
}
  
//-------------------------------------------------------------------------------------------------
std::string getPrecisionLevelName(const PrecisionLevel plevel) {
  switch (plevel) {
  case PrecisionLevel::SINGLE:
    return std::string("SINGLE");
  case PrecisionLevel::SINGLE_PLUS:
    return std::string("SINGLE_PLUS");
  case PrecisionLevel::DOUBLE:
    return std::string("DOUBLE");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
ForceAccumulationMethod translateForceAccumulationMethod(const std::string &choice,
                                                         const ExceptionResponse policy) {
  if (strcmpCased(choice, std::string("split"))) {
    return ForceAccumulationMethod::SPLIT;
  }
  else if (strcmpCased(choice, std::string("whole"))) {    
    return ForceAccumulationMethod::WHOLE;
  }
  else if (strcmpCased(choice, std::string("automatic"))) {    
    return ForceAccumulationMethod::AUTOMATIC;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getForceAccumulationMethodName(const ForceAccumulationMethod method) {
  switch (method) {
  case ForceAccumulationMethod::SPLIT:
    return std::string("SPLIT");
  case ForceAccumulationMethod::WHOLE:
    return std::string("WHOLE");
  case ForceAccumulationMethod::AUTOMATIC:
    return std::string("AUTOMATIC");
  }
  __builtin_unreachable();
}

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

//-------------------------------------------------------------------------------------------------
void checkEnergyBits(const int choice) {
  if (choice < min_energy_scale_bits || choice > max_energy_scale_bits) {
    rtErr("Energy accumulation is allowed in a fixed precision range of " +
          std::to_string(min_energy_scale_bits) + " to " +
          std::to_string(max_energy_scale_bits) + " bits.  A value of " + std::to_string(choice) +
          fixedPrecisionRangeErrorMessage(choice, min_energy_scale_bits, max_energy_scale_bits),
          "checkEnergyBits");
  }
}

//-------------------------------------------------------------------------------------------------
void checkChargeMeshBits(const int choice, const PrecisionLevel pmodel) {
  if (choice < min_charge_mesh_scale_bits) {
    rtErr("Charge mesh accumulation must take place with at least " +
          std::to_string(min_charge_mesh_scale_bits) + " bits.  A values of " +
          std::to_string(choice) +
          fixedPrecisionRangeErrorMessage(choice, min_charge_mesh_scale_bits, 64),
          "checkChargeMeshBits");
  }
  switch (pmodel) {
  case PrecisionLevel::SINGLE:
  case PrecisionLevel::SINGLE_PLUS:
    if (choice > 29) {
      rtErr("Charge mesh accumulation in a " + getPrecisionLevelName(pmodel) + " precision model "
            "will take place with a 32-bit signed integer accumulation grid.  A value of " +
            std::to_string(choice) + fixedPrecisionRangeErrorMessage(choice, 8, 29),
            "checkChargeMeshBits");
    }
    break;
  case PrecisionLevel::DOUBLE:
    if (choice > max_charge_mesh_scale_bits) {
      rtErr("Charge mesh accumulation in a " + getPrecisionLevelName(pmodel) + " precision model "
            "will take place with a 64-bit signed integer accumulation grid.  A value of " +
            std::to_string(choice) +
            fixedPrecisionRangeErrorMessage(choice, 8, max_charge_mesh_scale_bits),
            "checkChargeMeshBits");
    }
    break;
  }
}

} // namespace numerics
} // namespace omni
