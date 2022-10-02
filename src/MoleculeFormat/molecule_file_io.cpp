#include "copyright.h"
#include "molecule_file_io.h"

namespace stormm {
namespace structure {

//-------------------------------------------------------------------------------------------------
int getMdlFormatEnd(const TextFileReader &tfr, const int line_start, const int line_end_in) {
  const int line_end = (line_end_in < 0) ? tfr.line_count : line_end_in;
  int result = line_start;
  bool found = false;
  while (found == false && result < line_end) {
    const char* lptr = &tfr.text[tfr.line_limits[result]];
    if (tfr.line_limits[result + 1] - tfr.line_limits[result] >= 6 &&
        lptr[0] == 'M' && lptr[1] == ' ' && lptr[2] == ' ' && lptr[3] == 'E' && lptr[4] == 'N' &&
        lptr[5] == 'D') {
      found = true;
    }
    else {
      result++;
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
int getMdlFormatEnd(const TextFile &tf, const int line_start, const int line_end_in) {
  return getMdlFormatEnd(tf.data(), line_start, line_end_in);
}

//-------------------------------------------------------------------------------------------------
int findMolObjVersion(const char* text, const int nchar) {
  for (int i = 0; i < nchar; i++) {
    if (text[i] == 'V' || text[i] == 'v') {
      if (i < nchar - 4 && text[i + 2] == '0' && text[i + 3] == '0' && text[i + 4] == '0') {
        if (text[i + 1] == '2') {
          return 2000;
        }
        else if (text[i + 1] == '3') {
          return 3000;
        }
      }
    }
  }
  return 2000;
}

//-------------------------------------------------------------------------------------------------
std::vector<int2> findSdfMolEntryLimits(const TextFile &tf) {
  const TextFileReader tfr = tf.data();
  int nsection = 0;
  int last_delimiter_line = -1;
  std::vector<int2> result;
  for (int i = 0; i < tfr.line_count; i++) {
    if (tfr.line_lengths[i] >= 4 &&
        tfr.text[tfr.line_limits[i]    ] == '$' && tfr.text[tfr.line_limits[i] + 1] == '$' &&
        tfr.text[tfr.line_limits[i] + 2] == '$' && tfr.text[tfr.line_limits[i] + 3] == '$') {

      // A MOL entry must have at least four lines, so any text preceding a $$$$ marker must be at
      // least four lines long to count as a valid entry worth parsing
      if ((last_delimiter_line >=  0 && i - last_delimiter_line >= 4) ||
          (last_delimiter_line == -1 && i >= 4)) {

        // Estimate the number of frames based on the length of the first
        if (nsection == 0) {
          const int nframe_est = (tfr.line_count / (i - last_delimiter_line)) + 1;
          result.reserve(nframe_est);
        }
        result.push_back({ last_delimiter_line + 1, i });
        nsection++;
      }
      last_delimiter_line = i;
    }
  }
  if (tfr.line_count - last_delimiter_line >= 4) {
    result.push_back({ last_delimiter_line, tfr.line_count });
  }
  result.shrink_to_fit();
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<double3> extractSdfCoordinates(const TextFile &tf, const int line_start_in,
                                           const int line_end_in) {
  std::vector<double3> result;

  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<double3> extractSdfCoordinates(const TextFile &tf, const int frame_number) {
  const std::vector<int2> limits = findSdfMolEntryLimits(tf);
  if (frame_number < 0 || frame_number >= static_cast<int>(limits.size())) {
    rtErr("The requested frame index (" + std::to_string(frame_number) + ") does not exist in an "
          "SD file with " + std::to_string(limits.size()) + " MDL MOL entries.",
          "extractSdfCoordinates");
  }
  return extractSdfCoordinates(tf, limits[frame_number].x, limits[frame_number].y);
}

} // namespace structure
} // namespace stormm
