#include <cstdio>
#include "copyright.h"
#include "file_util.h"
#include "file_listing.h"
#include "Reporting/error_format.h"

namespace stormm {
namespace diskutil {

/// \brief Open a file for output writing.  This encapsulates error messages in the event that
///        the file cannot be opened as expected.
///
/// \param filename     Name of the file to write
/// \param expectation  Conditions to look for in order to successfully open the file for writing
/// \param style        The kind of file to write, ascii or binary
/// \param description  Description of the file, to help identify problems if it cannot be written
std::ofstream openOutputFile(const std::string &filename, const PrintSituation expectation,
                             const std::string &description, const DataFormat style) {

  // Check the file conditions: look before you leap
  std::ofstream foutp;
  switch (expectation) {
  case PrintSituation::UNKNOWN:
  case PrintSituation::OPEN_NEW:
    switch (getDrivePathType(filename)) {
    case DrivePathType::FILE:
      rtErr("Unable to open a new file " + filename + " because it already exists.  Activity: " +
            description + ".", "openOutputFile");
      break;
    case DrivePathType::DIRECTORY:
      rtErr("Unable to open a new file " + filename + " because a directory of the same name "
            "exists.  Activity: " + description + ".", "openOutputFile");
      break;
    case DrivePathType::REGEXP:
      switch (style) {
      case DataFormat::ASCII:
        foutp.open(filename, std::ofstream::out);
        break;
      case DataFormat::BINARY:

	// The std::ofstream specifiers occupy different bits of a single
	// integer-like data type.  Use the | operator to combine them.
        foutp.open(filename, std::ofstream::out | std::ofstream::binary);
        break;
      }
      break;
    }
    break;
  case PrintSituation::APPEND:
    switch (style) {
    case DataFormat::ASCII:
      foutp.open(filename, std::ofstream::ate | std::ofstream::app);
      break;
    case DataFormat::BINARY:
      foutp.open(filename, std::ofstream::ate | std::ofstream::app | std::ofstream::binary);
      break;
    }
    break;
  case PrintSituation::OVERWRITE:
    switch (style) {
    case DataFormat::ASCII:
      foutp.open(filename, std::ofstream::trunc);
      break;
    case DataFormat::BINARY:
      foutp.open(filename, std::ofstream::trunc | std::ofstream::binary);
      break;
    }
    break;
  }

  // Check that the file has been opened
  if (foutp.is_open() == false) {
    rtErr("Attempt to open file " + filename + " failed.  Bad writing permissions or insufficient "
          "disk space may be the problem.  Activity: " + description + ".", "openOutputFile");
  }
  return foutp;
}

/// \brief Remove a file.  This encapsulates error messages in the event that the file cannot be
///        removed as requested.
int removeFile(const std::string &filename, const ExceptionResponse policy) {
  switch (getDrivePathType(filename)) {
  case DrivePathType::FILE:
    break;
  case DrivePathType::DIRECTORY:
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr(filename + " is a directory.", "removeFile");
    case ExceptionResponse::WARN:
      rtWarn(filename + " is a directory.  Use stormmRmdir() to remove it.", "removeFile");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
    return 1;
  case DrivePathType::REGEXP:
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr(filename + " does not exist and therefore cannot be removed.", "removeFile");
    case ExceptionResponse::WARN:
      rtWarn(filename + " does not exist and therefore cannot be removed.", "removeFile");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
    return 1;
  }
  const int result = remove(filename.c_str());
  if (result != 0) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr(filename + " could not be removed as requested.", "removeFile");
    case ExceptionResponse::WARN:
      rtWarn(filename + " could not be removed as requested.", "removeFile");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
  }
  return result;
}

} // namespace diskutil
} // namespace stormm
