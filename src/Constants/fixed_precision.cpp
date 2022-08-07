#include "Parsing/parse.h"
#include "Reporting/error_format.h"
#include "fixed_precision.h"

namespace stormm {
namespace numerics {

using constants::getPrecisionModelName;
using parse::strcmpCased;

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
void checkChargeMeshBits(const int choice, const PrecisionModel pmodel) {
  if (choice < min_charge_mesh_scale_bits) {
    rtErr("Charge mesh accumulation must take place with at least " +
          std::to_string(min_charge_mesh_scale_bits) + " bits.  A values of " +
          std::to_string(choice) +
          fixedPrecisionRangeErrorMessage(choice, min_charge_mesh_scale_bits, 64),
          "checkChargeMeshBits");
  }
  switch (pmodel) {
  case PrecisionModel::SINGLE:
    if (choice > 31) {
      rtErr("Charge mesh accumulation in a " + getPrecisionModelName(pmodel) + " precision model "
            "will take place with a 32-bit signed integer accumulation grid.  A value of " +
            std::to_string(choice) + fixedPrecisionRangeErrorMessage(choice, 8, 31),
            "checkChargeMeshBits");
    }
    break;
  case PrecisionModel::DOUBLE:
    if (choice > max_charge_mesh_scale_bits) {
      rtErr("Charge mesh accumulation in a " + getPrecisionModelName(pmodel) + " precision model "
            "will take place with a 64-bit signed integer accumulation grid.  A value of " +
            std::to_string(choice) +
            fixedPrecisionRangeErrorMessage(choice, 8, max_charge_mesh_scale_bits),
            "checkChargeMeshBits");
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void splitRealConversion(const float fval, int *primary, int *overflow) {
  int ival;
  if (fabsf(fval) >= max_int_accumulation_f) {
    const int spillover = fval / max_int_accumulation_f;
    *primary = fval - (static_cast<float>(spillover) * max_int_accumulation_f);
    *overflow = spillover;
  }
  else {
    *primary = fval;
    *overflow = 0;
  }
}

//-------------------------------------------------------------------------------------------------
void splitRealConversion(const double dval, llint *primary, int *overflow) {
  if (fabs(dval) >= max_llint_accumulation) {
    const int spillover = dval / max_llint_accumulation;
    *primary = dval - (static_cast<double>(spillover) * max_llint_accumulation);
    *overflow = spillover;
  }
  else {
    *primary = dval;
    *overflow = 0;
  }
}

//-------------------------------------------------------------------------------------------------
void splitRealAccumulation(const float fval, int *primary, int *overflow) {
  int ival;
  if (fabsf(fval) >= max_int_accumulation_f) {
    const int spillover = fval / max_int_accumulation_f;
    ival = fval - (static_cast<float>(spillover) * max_int_accumulation_f);
    *overflow += spillover;
  }
  else {
    ival = fval;
  }
  const int prim_old = *primary;
  *primary += ival;
  const int prim_old_plus_ival = prim_old + ival;
  if ((prim_old ^ prim_old_plus_ival) < 0 && (prim_old ^ ival) >= 0) {
    *overflow += (1 - (2 * (ival < 0))) * 2;
  }
}

//-------------------------------------------------------------------------------------------------
void splitRealAccumulation(const double dval, llint *primary, int *overflow) {
  llint ival;
  if (fabs(dval) >= max_llint_accumulation) {
    const int spillover = dval / max_llint_accumulation;
    ival = dval - (static_cast<double>(spillover) * max_llint_accumulation);
    *overflow += spillover;
  }
  else {
    ival = dval;
  }
  const llint prim_old = *primary;
  *primary += ival;
  const llint prim_old_plus_ival = prim_old + ival;
  if ((prim_old ^ prim_old_plus_ival) < 0LL && (prim_old ^ ival) >= 0LL) {
    *overflow += (1 - (2 * (ival < 0LL))) * 2;
  }
}

//-------------------------------------------------------------------------------------------------
ForceAccumulationMethod chooseForceAccumulationMethod(const int frc_bits) {
  if (frc_bits <= 24) {
    return ForceAccumulationMethod::SPLIT;
  }
  else {
    return ForceAccumulationMethod::WHOLE;
  }
  __builtin_unreachable();
}
  
} // namespace numerics
} // namespace stormm
