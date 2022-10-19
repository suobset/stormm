// -*-c++-*-
#ifndef STORMM_FILE_UTIL_H
#define STORMM_FILE_UTIL_H

#include <fstream>
#include <iostream>
#include <string>
#include "copyright.h"
#include "Constants/behavior.h"
#include "file_enumerators.h"

namespace stormm {
namespace diskutil {

using constants::ExceptionResponse;

/// \brief Open a file for output writing.  This encapsulates error messages in the event that
///        the file cannot be opened as expected.
///
/// \param filename     Name of the file to write
/// \param expectation  Conditions to look for in order to successfully open the file for writing
/// \param style        The kind of file to write, ascii or binary
/// \param description  Description of the file, to help identify problems if it cannot be written
std::ofstream openOutputFile(const std::string &filename,
                             PrintSituation expectation = PrintSituation::OPEN_NEW,
                             const std::string &description = std::string(""),
                             DataFormat style = DataFormat::ASCII);

/// \brief Remove a file.  This encapsulates error messages in the event that the file cannot be
///        removed as requested.
///
/// \param filename  Name of the file to remove
/// \param policy    Indicates what to do if the file cannot be removed for some reason
int removeFile(const std::string &filename, ExceptionResponse policy = ExceptionResponse::WARN);

} // namespace diskutil
} // namespace stormm

#endif
