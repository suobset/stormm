#include "copyright.h"
#include "Parsing/parse.h"
#include "residue_dictionary.h"

namespace stormm {
namespace chemistry {

using parse::uppercase;

//-------------------------------------------------------------------------------------------------
bool isAminoAcid(const char rx, const char ry, const char rz) {
  return ((rx == 'A' && ry == 'L' && rz == 'A') || (rx == 'A' && ry == 'R' && rz == 'G') ||
          (rx == 'A' && ry == 'S' && rz == 'N') || (rx == 'A' && ry == 'S' && rz == 'P') ||
          (rx == 'A' && ry == 'S' && rz == 'H') || (rx == 'C' && ry == 'Y' && rz == 'S') ||
          (rx == 'G' && ry == 'L' && rz == 'H') || (rx == 'G' && ry == 'L' && rz == 'N') ||
          (rx == 'G' && ry == 'L' && rz == 'U') || (rx == 'G' && ry == 'L' && rz == 'Y') ||
          (rx == 'H' && ry == 'I' && rz == 'D') || (rx == 'H' && ry == 'I' && rz == 'E') ||
          (rx == 'H' && ry == 'I' && rz == 'P') || (rx == 'H' && ry == 'I' && rz == 'S') ||
          (rx == 'I' && ry == 'L' && rz == 'E') || (rx == 'L' && ry == 'E' && rz == 'U') ||
          (rx == 'L' && ry == 'Y' && rz == 'N') || (rx == 'L' && ry == 'Y' && rz == 'S') ||
          (rx == 'M' && ry == 'E' && rz == 'T') || (rx == 'P' && ry == 'H' && rz == 'R') ||
          (rx == 'P' && ry == 'R' && rz == 'O') || (rx == 'S' && ry == 'E' && rz == 'R') ||
          (rx == 'T' && ry == 'H' && rz == 'R') || (rx == 'T' && ry == 'R' && rz == 'P') ||
          (rx == 'T' && ry == 'Y' && rz == 'R') || (rx == 'V' && ry == 'A' && rz == 'L'));
}
  
//-------------------------------------------------------------------------------------------------
bool isAminoAcid(const char4 residue_name) {
  char4 rn = uppercase(residue_name);
  return ((isAminoAcid(rn.x, rn.y, rn.z) && (rn.w == ' ' || rn.w == 'C' || rn.w == 'N')) ||
          ((rn.x == 'C' || rn.x == 'N' || rn.x == ' ') && isAminoAcid(rn.y, rn.z, rn.w)));
}

//-------------------------------------------------------------------------------------------------
bool isNucleicAcid(const char rx, const char ry, const char rz) {
  return ((rx == 'A' && ry == 'D' && rz == 'E') || (rx == 'C' && ry == 'Y' && rz == 'T') ||
          (rx == 'G' && ry == 'U' && rz == 'A') || (rx == 'T' && ry == 'H' && rz == 'Y') ||
          (rx == 'U' && ry == 'R' && rz == 'A'));
}
  
//-------------------------------------------------------------------------------------------------
bool isNucleicAcid(const char4 residue_name) {
  char4 rn = uppercase(residue_name);
  return ((isNucleicAcid(rn.x, rn.y, rn.z) && rn.w == ' ') ||
          (rn.x == ' ' && isNucleicAcid(rn.y, rn.z, rn.w)));
}

} // namespace chemistry
} // namespace stormm
