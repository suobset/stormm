#include <fstream>
#include <iostream>
#include <string>
#include "copyright.h"
#include "Constants/behavior.h"
#include "Parsing/parse.h"
#include "Reporting/error_format.h"
#include "trajectory_enumerators.h"

namespace stormm {
namespace trajectory {

using constants::CaseSensitivity;
using parse::strcmpCased;

//-------------------------------------------------------------------------------------------------
std::string getCoordinateFileKindDescription(const CoordinateFileKind cfkind) {
  switch (cfkind) {
  case CoordinateFileKind::AMBER_CRD:
    return std::string("Amber ascii trajectory");
  case CoordinateFileKind::AMBER_INPCRD:
    return std::string("Amber input coordinates");
  case CoordinateFileKind::AMBER_ASCII_RST:
    return std::string("Amber ascii restart");
  case CoordinateFileKind::AMBER_NETCDF:
    return std::string("Amber NetCDF binary trajectory");
  case CoordinateFileKind::AMBER_NETCDF_RST:
    return std::string("Amber NetCDF binary restart");
  case CoordinateFileKind::SDF:
    return std::string("MDL MOL / SDF");
  case CoordinateFileKind::UNKNOWN:
    return std::string("Unknown coordinate file format");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const CoordinateFileKind cfkind) {
  switch (cfkind) {
  case CoordinateFileKind::AMBER_CRD:
    return std::string("AMBER_CRD");
  case CoordinateFileKind::AMBER_INPCRD:
    return std::string("AMBER_INPCRD");
  case CoordinateFileKind::AMBER_ASCII_RST:
    return std::string("AMBER_ASCII_RST");
  case CoordinateFileKind::AMBER_NETCDF:
    return std::string("AMBER_NETCDF");
  case CoordinateFileKind::AMBER_NETCDF_RST:
    return std::string("AMBER_NETCDF_RST");
  case CoordinateFileKind::SDF:
    return std::string("SDF");
  case CoordinateFileKind::UNKNOWN:
    return std::string("UNKNOWN");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const AncdfVariable key) {
  switch (key) {
  case AncdfVariable::NCFRAME:
    return std::string("frame");
  case AncdfVariable::NCSPATIAL:
    return std::string("spatial");
  case AncdfVariable::NCATOM:
    return std::string("atom");
  case AncdfVariable::NCCELL_SPATIAL:
    return std::string("cell_spatial");
  case AncdfVariable::NCCELL_LENGTHS:
    return std::string("cell_lengths");
  case AncdfVariable::NCCELL_ANGULAR:
    return std::string("cell_angular");
  case AncdfVariable::NCCELL_ANGLES:
    return std::string("cell_angles");
  case AncdfVariable::NCCOORDS:
    return std::string("coordinates");
  case AncdfVariable::NCVELO:
    return std::string("velocities");
  case AncdfVariable::NCTEMPERATURE:
    return std::string("temp0");
  case AncdfVariable::NCTIME:
    return std::string("time");
  case AncdfVariable::NCLABEL:
    return std::string("label");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const CoordinateFileRole cpkind) {
  switch (cpkind) {
  case CoordinateFileRole::INITIATE:
    return std::string("INITIATE");
  case CoordinateFileRole::TRAJECTORY:
    return std::string("TRAJECTORY");
  case CoordinateFileRole::CHECKPOINT:
    return std::string("CHECKPOINT");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const TrajectoryFusion protocol) {
  switch (protocol) {
  case TrajectoryFusion::ON:
    return std::string("ON");
  case TrajectoryFusion::OFF:
    return std::string("OFF");
  case TrajectoryFusion::AUTO:
    return std::string("AUTO");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const TrajectoryKind input) {
  switch (input) {
  case TrajectoryKind::POSITIONS:
    return std::string("POSITIONS");
  case TrajectoryKind::VELOCITIES:
    return std::string("VELOCITIES");
  case TrajectoryKind::FORCES:
    return std::string("FORCES");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const CoordinateCycle orientation) {
  switch (orientation) {
  case CoordinateCycle::PRIMARY:
    return std::string("PRIMARY");
  case CoordinateCycle::ALTERNATE:
    return std::string("ALTERNATE");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
CoordinateFileKind translateCoordinateFileKind(const std::string &name_in) {
  if (strcmpCased(name_in, "AMBER_CRD", CaseSensitivity::NO)) {
    return CoordinateFileKind::AMBER_CRD;
  }
  else if (strcmpCased(name_in, "AMBER_INPCRD", CaseSensitivity::NO)) {
    return CoordinateFileKind::AMBER_INPCRD;
  }
  else if (strcmpCased(name_in, "AMBER_ASCII_RST", CaseSensitivity::NO)) {
    return CoordinateFileKind::AMBER_ASCII_RST;
  }
  else if (strcmpCased(name_in, "AMBER_NETCDF", CaseSensitivity::NO)) {
    return CoordinateFileKind::AMBER_NETCDF;
  }
  else if (strcmpCased(name_in, "AMBER_NETCDF_RST", CaseSensitivity::NO)) {
    return CoordinateFileKind::AMBER_NETCDF_RST;
  }
  else if (strcmpCased(name_in, "SDF", CaseSensitivity::NO) ||
           strcmpCased(name_in, "MDL_MOL", CaseSensitivity::NO)) {
    return CoordinateFileKind::SDF;
  }
  else if (strcmpCased(name_in, "UNKNOWN", CaseSensitivity::NO)) {
    return CoordinateFileKind::UNKNOWN;
  }
  else {
    rtErr("Unrecognized coordinate file enumeration " + name_in + ".",
          "translateCoordinateFileKind");
  }
  __builtin_unreachable();
}

} // namespace trajectory
} // namespace stormm
