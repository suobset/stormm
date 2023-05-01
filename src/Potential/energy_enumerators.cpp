#include "copyright.h"
#include "Constants/behavior.h"
#include "Parsing/parse.h"
#include "energy_enumerators.h"

namespace stormm {
namespace energy {

using constants::CaseSensitivity;
using parse::strcmpCased;

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const EvaluateForce input) {
  switch (input) {
  case EvaluateForce::NO:
    return std::string("NO");
  case EvaluateForce::YES:
    return std::string("YES");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const EvaluateEnergy input) {
  switch (input) {
  case EvaluateEnergy::NO:
    return std::string("NO");
  case EvaluateEnergy::YES:
    return std::string("YES");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const EvaluateVirial input) {
  switch (input) {
  case EvaluateVirial::NO:
    return std::string("NO");
  case EvaluateVirial::YES:
    return std::string("YES");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const DihedralStyle input) {
  switch (input) {
  case DihedralStyle::COSINE:
    return std::string("COSINE");
  case DihedralStyle::HARMONIC:
    return std::string("HARMONIC");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const StateVariable input) {
  switch (input) {
  case StateVariable::BOND:
    return std::string("BOND");
  case StateVariable::ANGLE:
    return std::string("ANGLE");
  case StateVariable::PROPER_DIHEDRAL:
    return std::string("PROPER_DIHEDRAL");
  case StateVariable::IMPROPER_DIHEDRAL:
    return std::string("IMPROPER_DIHEDRAL");
  case StateVariable::UREY_BRADLEY:
    return std::string("UREY_BRADLEY");
  case StateVariable::CHARMM_IMPROPER:
    return std::string("CHARMM_IMPROPER");
  case StateVariable::CMAP:
    return std::string("CMAP");
  case StateVariable::VDW:
    return std::string("VDW");
  case StateVariable::VDW_ONE_FOUR:
    return std::string("VDW_ONE_FOUR");
  case StateVariable::ELECTROSTATIC:
    return std::string("ELECTROSTATIC");
  case StateVariable::ELECTROSTATIC_ONE_FOUR:
    return std::string("ELECTROSTATIC_ONE_FOUR");
  case StateVariable::GENERALIZED_BORN:
    return std::string("GENERALIZED_BORN");
  case StateVariable::RESTRAINT:
    return std::string("RESTRAINT");
  case StateVariable::KINETIC:
    return std::string("KINETIC");
  case StateVariable::PRESSURE:
    return std::string("PRESSURE");
  case StateVariable::VIRIAL_11:
    return std::string("VIRIAL_11");
  case StateVariable::VIRIAL_12:
    return std::string("VIRIAL_12");
  case StateVariable::VIRIAL_22:
    return std::string("VIRIAL_22");
  case StateVariable::VIRIAL_13:
    return std::string("VIRIAL_13");
  case StateVariable::VIRIAL_23:
    return std::string("VIRIAL_23");
  case StateVariable::VIRIAL_33:
    return std::string("VIRIAL_33");
  case StateVariable::VOLUME:
    return std::string("VOLUME");
  case StateVariable::TEMPERATURE_ALL:
    return std::string("TEMPERATURE_ALL");
  case StateVariable::TEMPERATURE_PROTEIN:
    return std::string("TEMPERATURE_PROTEIN");
  case StateVariable::TEMPERATURE_LIGAND:
    return std::string("TEMPERATURE_LIGAND");
  case StateVariable::TEMPERATURE_SOLVENT:
    return std::string("TEMPERATURE_SOLVENT");
  case StateVariable::DU_DLAMBDA:
    return std::string("DU_DLAMBDA");
  case StateVariable::POTENTIAL_ENERGY:
    return std::string("POTENTIAL_ENERGY");
  case StateVariable::TOTAL_ENERGY:
    return std::string("TOTAL_ENERGY");
  case StateVariable::ALL_STATES:
    return std::string("ALL_STATES");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const NonbondedPotential input) {
  switch (input) {
  case NonbondedPotential::ELECTROSTATIC:
    return std::string("ELECTROSTATIC");
  case NonbondedPotential::VAN_DER_WAALS:
    return std::string("VAN_DER_WAALS");
  case NonbondedPotential::CLASH:
    return std::string("CLASH");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const VdwCombiningRule input) {
  switch (input) {
  case VdwCombiningRule::LORENTZ_BERTHELOT:
    return std::string("LORENTZ_BERTHELOT");
  case VdwCombiningRule::GEOMETRIC:
    return std::string("GEOMETRIC");
  case VdwCombiningRule::NBFIX:
    return std::string("NBFIX");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const ClashResponse input) {
  switch (input) {
  case ClashResponse::NONE:
    return std::string("NONE");
  case ClashResponse::FORGIVE:
    return std::string("FORGIVE");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const EnergySample input) {
  switch (input) {
  case EnergySample::TIME_SERIES:
    return std::string("TIME_SERIES");
  case EnergySample::FINAL:
    return std::string("FINAL");
  case EnergySample::TIME_AVERAGE:
    return std::string("TIME_AVERAGE");
  }
  __builtin_unreachable();
}
  
} // namespace energy
} // namespace stormm
