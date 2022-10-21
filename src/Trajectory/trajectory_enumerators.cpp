#include <fstream>
#include <iostream>
#include <string>
#include "copyright.h"
#include "Parsing/ascii_numbers.h"
#include "Parsing/parse.h"
#include "Parsing/polynumeric.h"
#include "Reporting/error_format.h"
#include "trajectory_enumerators.h"

namespace stormm {
namespace trajectory {

using parse::NumberFormat;
using parse::readIntegerValue;
using parse::separateText;
using parse::strcmpCased;
using parse::TextFileReader;
using parse::TextOrigin;
using parse::verifyNumberFormat;

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
DataFormat getTrajectoryFormat(CoordinateFileKind cfkind) {
  switch (cfkind) {
  case CoordinateFileKind::AMBER_CRD:
  case CoordinateFileKind::AMBER_INPCRD:
  case CoordinateFileKind::AMBER_ASCII_RST:
  case CoordinateFileKind::SDF:
    return DataFormat::ASCII;
  case CoordinateFileKind::AMBER_NETCDF:
  case CoordinateFileKind::AMBER_NETCDF_RST:
    return DataFormat::BINARY;
  case CoordinateFileKind::UNKNOWN:
    rtWarn("Unable to determine the nature of an UNKNOWN kind trajectory file.  Reporting BINARY.",
           "getTrajectoryFormat");
    return DataFormat::BINARY;    
  }
  __builtin_unreachable();
}
  
//-------------------------------------------------------------------------------------------------
CoordinateFileKind translateCoordinateFileKind(const std::string &name_in) {
  if (strcmpCased(name_in, "AMBER_CRD")) {
    return CoordinateFileKind::AMBER_CRD;
  }
  else if (strcmpCased(name_in, "AMBER_INPCRD")) {
    return CoordinateFileKind::AMBER_INPCRD;
  }
  else if (strcmpCased(name_in, "AMBER_ASCII_RST")) {
    return CoordinateFileKind::AMBER_ASCII_RST;
  }
  else if (strcmpCased(name_in, "AMBER_NETCDF")) {
    return CoordinateFileKind::AMBER_NETCDF;
  }
  else if (strcmpCased(name_in, "AMBER_NETCDF_RST")) {
    return CoordinateFileKind::AMBER_NETCDF_RST;
  }
  else if (strcmpCased(name_in, "SDF") || strcmpCased(name_in, "MDL_MOL")) {
    return CoordinateFileKind::SDF;
  }
  else if (strcmpCased(name_in, "UNKNOWN")) {
    return CoordinateFileKind::UNKNOWN;
  }
  else {
    rtErr("Unrecognized coordinate file enumeration " + name_in + ".",
          "translateCoordinateFileKind");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
CoordinateFileKind detectCoordinateFileKind(const TextFile &tf) {
  const TextFileReader tfr = tf.data();

  // Test for an Amber coordinate file
  if (tfr.line_count < 3) {
    return CoordinateFileKind::UNKNOWN;
  }
  const std::vector<std::string> line_two = separateText(&tfr.text[tfr.line_limits[1]],
                                                         tfr.line_limits[2] - tfr.line_limits[1]);
  if (line_two.size() == 2 && verifyNumberFormat(line_two[0].c_str(), NumberFormat::INTEGER) &&
      stoi(line_two[0]) > 0 &&
      (verifyNumberFormat(line_two[1].c_str(), NumberFormat::SCIENTIFIC) ||
       verifyNumberFormat(line_two[1].c_str(), NumberFormat::STANDARD_REAL))) {

    // Check the third line to confirm
    if (tfr.line_limits[3] - tfr.line_limits[2] >= 36 &&
        verifyNumberFormat(tfr.text, NumberFormat::STANDARD_REAL, tfr.line_limits[2], 12) &&
        verifyNumberFormat(tfr.text, NumberFormat::STANDARD_REAL, tfr.line_limits[2] + 12, 12) &&
        verifyNumberFormat(tfr.text, NumberFormat::STANDARD_REAL, tfr.line_limits[2] + 24, 12)) {

      // Some input files will still have a time stamp.  A check for the presence of actual
      // velocities will be needed later, as this routine may have only been fed a stub from the
      // head of the file.
      return CoordinateFileKind::AMBER_ASCII_RST;
    }
  }
  else if (line_two.size() == 1 &&
           verifyNumberFormat(line_two[0].c_str(), NumberFormat::INTEGER) &&
           stoi(line_two[0]) > 0) {

    // Check the third line to confirm
    if (tfr.line_limits[3] - tfr.line_limits[2] >= 36 &&
        verifyNumberFormat(tfr.text, NumberFormat::STANDARD_REAL, tfr.line_limits[2], 12) &&
        verifyNumberFormat(tfr.text, NumberFormat::STANDARD_REAL, tfr.line_limits[2] + 12, 12) &&
        verifyNumberFormat(tfr.text, NumberFormat::STANDARD_REAL, tfr.line_limits[2] + 24, 12)) {
      return CoordinateFileKind::AMBER_INPCRD;
    }
  }
  else if (tfr.line_limits[2] - tfr.line_limits[1] >= 24 &&
           verifyNumberFormat(tfr.text, NumberFormat::STANDARD_REAL, tfr.line_limits[1], 8) &&
           verifyNumberFormat(tfr.text, NumberFormat::STANDARD_REAL, tfr.line_limits[1] + 8, 8) &&
           verifyNumberFormat(tfr.text, NumberFormat::STANDARD_REAL, tfr.line_limits[1] + 16, 8)) {
    return CoordinateFileKind::AMBER_CRD;
  }
  
  // If the fourth line looks like an MDL MOL-format counts line, return that result after checking
  // the first atom line.
  if (tfr.line_count > 4 &&
      verifyContents(tfr, 3, 0, 3, NumberFormat::INTEGER) &&
      verifyContents(tfr, 3, 3, 3, NumberFormat::INTEGER) &&
      readIntegerValue(&tfr.text[tfr.line_limits[3]], 0, 3) >= 0) {
    if (verifyNumberFormat(tfr.text, NumberFormat::STANDARD_REAL, tfr.line_limits[4], 10) &&
        verifyNumberFormat(tfr.text, NumberFormat::STANDARD_REAL, tfr.line_limits[4] + 10, 10) &&
        verifyNumberFormat(tfr.text, NumberFormat::STANDARD_REAL, tfr.line_limits[4] + 20, 10)) {
      return CoordinateFileKind::SDF;
    }
  }

  // If the phrase "V2000" or "V3000" can be detected anywhere in the fourth through tenth lines,
  // return that result as well.
  const int nsearch_lines = std::min(10, tfr.line_count);
  for (int i = 3; i < nsearch_lines; i++) {
    const int jlim = tfr.line_limits[i + 1] - 4;
    for (int j = tfr.line_limits[i]; j < jlim; j++) {
      if (tfr.text[j] == 'V' && (tfr.text[j + 1] == '2' || tfr.text[j + 1] == '3') &&
          tfr.text[j + 2] == '0' && tfr.text[j + 3] == '0' && tfr.text[j + 4] == '0') {
        return CoordinateFileKind::SDF;
      }
    }
  }
  
  // If the file type has not yet been ascertained, return UNKNOWN
  return CoordinateFileKind::UNKNOWN;
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
  return detectCoordinateFileKind(tf);
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

} // namespace trajectory
} // namespace stormm
