#include "trajectory_enumerators.h"

namespace omni {
namespace trajectory {

//-------------------------------------------------------------------------------------------------
std::string getCoordinateFileKindName(const CoordinateFileKind cfkind) {
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
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getAncdfVariableName(const AncdfVariable key) {
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

} // namespace trajectory
} // namespace omni
