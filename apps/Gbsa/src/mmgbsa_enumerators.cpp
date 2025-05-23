#include "copyright.h"
#include "../../../src/Constants/behavior.h"
#include "../../../src/Parsing/parse.h"
#include "../../../src/Reporting/error_format.h"
#include "mmgbsa_enumerators.h"

namespace mmgbsa {

using stormm::errors::rtErr;
using stormm::parse::strcmpCased;
using stormm::constants::CaseSensitivity;
  
//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const MMGBSAProblemKind input) {
  switch (input) {
  case MMGBSAProblemKind::LIGAND:
    return std::string("LIGAND");
  case MMGBSAProblemKind::RECEPTOR:
    return std::string("RECEPTOR");
  case MMGBSAProblemKind::COMPLEX:
    return std::string("COMPLEX");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const WeightingScheme input) {
  switch (input) {
  case WeightingScheme::BOLTZMANN:
    return std::string("BOLTZMANN");
  case WeightingScheme::FLAT:
    return std::string("FLAT");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const EnergyReport input) {
  switch (input) {
  case EnergyReport::COMPLETE:
    return std::string("COMPLETE");
  case EnergyReport::COMPOUND_AVERAGES:
    return std::string("COMPOUND_AVERAGES");
  case EnergyReport::BEST_RESULTS:
    return std::string("BEST_RESULTS");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const PrintedPoseFormat input) {
  switch (input) {
  case PrintedPoseFormat::PDB:
    return std::string("PDB");
  case PrintedPoseFormat::AMBER:
    return std::string("AMBER");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
WeightingScheme translateWeightingScheme(const std::string &input) {
  if (strcmpCased(input, "boltzmann", CaseSensitivity::NO)) {
    return WeightingScheme::BOLTZMANN;
  }
  else if (strcmpCased(input, "flat", CaseSensitivity::NO) ||
           strcmpCased(input, "even", CaseSensitivity::NO) ||
           strcmpCased(input, "equal", CaseSensitivity::NO) ||
           strcmpCased(input, "none", CaseSensitivity::NO)) {
    return WeightingScheme::FLAT;
  }
  else {
    rtErr("The input " + input + " cooresponds to no known weighting scheme.",
          "translateWeightingScheme");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
EnergyReport translateEnergyReport(const std::string &input) {
  if (strcmpCased(input, "complete", CaseSensitivity::NO) ||
      strcmpCased(input, "all", CaseSensitivity::NO)) {
    return EnergyReport::COMPLETE;
  }
  else if (strcmpCased(input, "averages", CaseSensitivity::NO)) {
    return EnergyReport::COMPOUND_AVERAGES;
  }
  else if (strcmpCased(input, "best", CaseSensitivity::NO)) {
    return EnergyReport::BEST_RESULTS;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
PrintedPoseFormat translatePrintedPoseFormat(const std::string &input) {
  if (strcmpCased(input, "pdb", CaseSensitivity::NO) ||
      strcmpCased(input, "ProteinDataBank", CaseSensitivity::NO)) {
    return PrintedPoseFormat::PDB;
  }
  else if (strcmpCased(input, "amber", CaseSensitivity::NO)) {
    return PrintedPoseFormat::AMBER;
  }
  __builtin_unreachable();
}

} // namespace mmgbsa
