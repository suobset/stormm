#include <fstream>
#include <iostream>
#include <string>
#include "Parsing/parse.h"
#include "Parsing/polynumeric.h"
#include "Parsing/textfile.h"
#include "Reporting/error_format.h"
#include "trajectory_enumerators.h"

namespace omni {
namespace trajectory {

using parse::NumberFormat;
using parse::separateText;
using parse::TextFile;
using parse::TextFileReader;
using parse::TextOrigin;
using parse::verifyNumberFormat;

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
  case CoordinateFileKind::UNKNOWN:
    return std::string("Unknown coordinate file format");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
CoordinateFileKind detectCoordinateFileKind(const std::string &file_name,
                                            const std::string &caller) {

  // Read the first 2kB of the file, if that much exists.  Look for chracters that might indicate
  // it is a binary file.  If no such characters are found, try to parse it as one of the known
  // ASCII formats.
  std::ifstream finp;
  finp.open(file_name.c_str());
  if (finp.is_open() == false) {
    if (caller.size() == 0) {
      rtErr(file_name + " was not found.", "detectCoordinateFileKind");
    }
    else {
      rtErr(file_name + " was not found when called from " + caller + ".",
            "detectCoordinateFileKind");
    }
  }
  const int maxchar = 2048;
  std::vector<char> buffer(maxchar);
  int pos = 0;
  char c;
  bool is_binary = false;
  while (pos < maxchar && finp.get(c) && is_binary == false) {
    is_binary = (is_binary || static_cast<int>(c) < 0);
    buffer[pos] = c;
    pos++;
  }
  finp.close();
  if (is_binary) {

    // TBD: Determine the binary format.
    return CoordinateFileKind::UNKNOWN;
  }
  const size_t nchar = pos;
  std::string first_part(nchar, ' ');
  for (size_t i = 0; i < nchar; i++) {
    first_part[i] = buffer[i];
  }
  const TextFile tf(file_name, TextOrigin::RAM, first_part, "detectCoordinateFileKind");
  const TextFileReader tfr = tf.data();

  // Test for an Amber coordinate file
  if (tfr.line_count < 3) {
    return CoordinateFileKind::UNKNOWN;
  }
  const std::vector<std::string> line_two = separateText(&tfr.text[tfr.line_limits[1]],
                                                         tfr.line_limits[2] - tfr.line_limits[1]);
  if (line_two.size() == 0) {
    return CoordinateFileKind::UNKNOWN;
  }
  if (line_two.size() == 2 &&
      (verifyNumberFormat(line_two[1].c_str(), NumberFormat::SCIENTIFIC) ||
       verifyNumberFormat(line_two[1].c_str(), NumberFormat::STANDARD_REAL))) {
    if (tfr.line_limits[3] - tfr.line_limits[2] >= 36 &&
        verifyNumberFormat(tfr.text, NumberFormat::STANDARD_REAL, tfr.line_limits[2], 12) &&
        verifyNumberFormat(tfr.text, NumberFormat::STANDARD_REAL, tfr.line_limits[2] + 12, 12) &&
        verifyNumberFormat(tfr.text, NumberFormat::STANDARD_REAL, tfr.line_limits[2] + 24, 12)) {
      return CoordinateFileKind::AMBER_ASCII_RST;
    }
    else {
      return CoordinateFileKind::UNKNOWN;
    }
  }
  else if (line_two.size() == 1 &&
           verifyNumberFormat(line_two[0].c_str(), NumberFormat::INTEGER)) {
    if (tfr.line_limits[3] - tfr.line_limits[2] >= 36 &&
        verifyNumberFormat(tfr.text, NumberFormat::STANDARD_REAL, tfr.line_limits[2], 12) &&
        verifyNumberFormat(tfr.text, NumberFormat::STANDARD_REAL, tfr.line_limits[2] + 12, 12) &&
        verifyNumberFormat(tfr.text, NumberFormat::STANDARD_REAL, tfr.line_limits[2] + 24, 12)) {
      return CoordinateFileKind::AMBER_INPCRD;
    }
    else {
      return CoordinateFileKind::UNKNOWN;
    }
  }
  else if (line_two.size() >= 3 &&
           verifyNumberFormat(line_two[0].c_str(), NumberFormat::STANDARD_REAL) &&
           verifyNumberFormat(line_two[1].c_str(), NumberFormat::STANDARD_REAL) &&
           verifyNumberFormat(line_two[2].c_str(), NumberFormat::STANDARD_REAL)) {
    if (tfr.line_limits[2] - tfr.line_limits[1] >= 24 &&
        verifyNumberFormat(tfr.text, NumberFormat::STANDARD_REAL, tfr.line_limits[1], 8) &&
        verifyNumberFormat(tfr.text, NumberFormat::STANDARD_REAL, tfr.line_limits[1] + 8, 8) &&
        verifyNumberFormat(tfr.text, NumberFormat::STANDARD_REAL, tfr.line_limits[1] + 16, 8)) {
      return CoordinateFileKind::AMBER_CRD;
    }
    else {
      return CoordinateFileKind::UNKNOWN;
    }
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
