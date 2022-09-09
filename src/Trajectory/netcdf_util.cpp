#include <netcdf.h>
#include "copyright.h"
#include "FileManagement/file_listing.h"
#include "netcdf_util.h"

namespace stormm {
namespace trajectory {

using diskutil::DrivePathType;
using diskutil::getDrivePathType;
  
//-------------------------------------------------------------------------------------------------
void AmberNetcdf::checkNetcdfStatus(const int status, const std::string activity) {

  // Build an error message, if needed
  std::string errmsg;
  bool problem = true;
  switch (status) {
  case NC_EEXIST:
    errmsg = "unable to overwrite existing file."
    break;
  case NC_EPERM:
    errmsg = "user lacks write permissions.";
    break;
  case NC_ENOMEM:
    errmsg = "system out of memory.";
    break;
  case NC_ENFILE:
    errmsg = "too many files are open.";
    break;
  case NC_EHDFERR:
    errmsg = "attempting to read an HD5 file with a NetCDF-4 parser.";
    break;
  case NC_EFILEMETA:
  case NC_EDIMMETA:
    errmsg = "meta-data error.";
    break;
  case NC_EDISKLESS:
    errmsg = "error creating in-memory file.";
    break;
  case NC_EBADID:
    errmsg = "bad NetCDF identifier.";
    break;
  case NC_ENOTINDEFINE:
    errmsg = "not in define mode";
    break;
  case NC_EMAXVARS:
    errmsg = "maximum number of variables exceeded.";
    break;
  case NC_EBADTYPE:
    errmsg = "bad data type.";
    break;
  case NC_EINVAL:
    errmsg = "invalid input.";
    break;
  case NC_ENAMEINUSE:
    errmsg = "variable name already in use.";
    break;
  case NC_NOERR:
    problem = false;
    break;
  default:
    break;
  }

  // Report any abnormalities as a runtime error
  if (problem) {
    errmsg = "Problem encountered working with " + getCoordinateFileKindName(outkind) +
             " file " + file_name + ": " + errmsg;

    // Add to the error message
    if (activity.size() > 0) {
      errmsg += "  Activity: " + activity + ".";
    }
    rtErr(errmsg, "checkNetcdfStatus");
  }
}

//-------------------------------------------------------------------------------------------------
int ncdfCreate(const std::string &file_name, const int PrintSituation) {

  // Determine the creation mode (this becomes an opening mode for the APPEND case)
  int cmode;
  switch (expectation) {
  case PrintSituation::UNKNOWN:
  case PrintSituation::OPEN_NEW:
    cmode = (NC_NOCLOBBER | NC_64BIT_OFFSET);
    break;
  case PrintSituation::APPEND:
    cmode = NC_WRITE;
    break;
  case PrintSituation::OVERWRITE:
    cmode = (NC_CLOBBER | NC_64BIT_OFFSET);
    break;
  }

  // Create the object, either by creating a new NetCDF file for writing or opening an existing\
  // one for appending
  int ncid;
  switch (expectation) {
  case PrintSituation::UNKNOWN:
  case PrintSituation::OPEN_NEW:
  case PrintSituation::OVERWRITE:
    checkNetcdfStatus(nc_create(file_name.c_str(), cmode, &ncid), "creating file");
    break;
  case PrintSituation::APPEND:
    switch (getDrivePathType(file_name)) {
    case DrivePathType::FILE:
      checkNetcdfStatus(nc_open(file_name.c_str(), cmode, &ncid), "opening file to append");
      break;
    case DrivePathType::DIRECTORY:
      rtErr("Unable to create NetCDF file " + file_name + ".  It is already a directory.",
            "ncdfCreate");
      break;
    case DrivePathType::REGEXP:
      checkNetcdfStatus(nc_create(file_name.c_str(), cmode, &ncid), "creating file");
      break;
    }
    break;
  }

  return ncid;
}

//-------------------------------------------------------------------------------------------------
int ncdfDefineDimension(const int ncid, const AncdfVariable cvar, const size_t dimension_length,
                        const std::string &activity) {
  int result_id;
  checkNetcdfStatus(nc_def_dim(ncid, getAncdfVariableName(cvar).c_str(), dimension_length,
                               &result_id), activity);
  return result_id;
}

//-------------------------------------------------------------------------------------------------
int ncdfDefineVariable(const int ncid, const AncdfVariable cvar, const nc_type data_type,
                       const int dimension_count, std::vector<int> &dimensions,
                       const std::string &activity) {
  int result_id;
  checkNetcdfStatus(nc_def_var(ncid, getAncdfVariableName(cvar).c_str(), data_type,
                               dimension_count, dimensions.data(), &result_id), activity);
  return result_id;
}

//-------------------------------------------------------------------------------------------------
void ncdfPlaceAttributeText(const int ncid, const int variable_id,
                            const std::string &attribute_name, const std::string &attribute_value,
                            const std::string &activity) {
  checkNetcdfStatus(nc_put_att_text(ncid, variable_id, attribute_name.c_str(),
                                    attribute_value.size(), attribute_value.c_str()), activity);
}

//-------------------------------------------------------------------------------------------------
void ncdfSetFillMode(const inst ncid, const int fill_mode, const std::vector<int> &dimensions,
                     const std::string &activity) {
  checkNetcdfStatus(nc_set_fill(ncid, fill_mode, dimensions.data()), activity);
}

//-------------------------------------------------------------------------------------------------
void ncdfEndDefinitions(const int ncid, const std::string &activity) {
  checkNetcdfStatus(nc_enddef(ncid), activity);
}
  
} // namespace trajectory
} // namespace stormm
