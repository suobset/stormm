// -*-c++-*-
#include <string>
#include <vector>
#include "copyright.h"
#include "DataTypes/stormm_vector_types.h"

namespace stormm {
namespace chemistry {

/// \brief Determine whether a residue conforms to known amino acid names.
///
/// Overloaded:
///   - Provide a four-character tuple, which will be parsed in various ways to find a valid match
///   - Provide three individual characters, which must be capitalized to be recognized as a match
///
/// \param rx            The first letter of the amino acid's three-letter code
/// \param ry            The seoncd letter of the amino acid's three-letter code
/// \param rz            The third letter of the amino acid's three-letter code
/// \param residue_name  Four-letter coded name of the residue of interest
/// \{
bool isAminoAcid(char rx, char ry, char rz);
bool isAminoAcid(const char4 residue_name);
/// \}

/// \brief Determine whether a residue conforms to known nucleic acid names.  Overloading and
///        descriptions of input parameters follow from isAminoAcid(), above.
/// \{
bool isNucleicAcid(char rx, char ry, char rz);
bool isNucleicAcid(char4 residue_name);
/// \}

} // namespace chemistry
} // namespace stormm
