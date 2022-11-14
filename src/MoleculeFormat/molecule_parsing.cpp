#include "copyright.h"
#include "molecule_parsing.h"

namespace stormm {
namespace structure {

//-------------------------------------------------------------------------------------------------
AtomMask maskFromSdfDataItem(const std::string &item_name, const MdlMol &molecule,
                             const AtomGraph *ag, const ChemicalFeatures *chemfe,
                             const CoordinateFrame &cf) {
  std::vector<std::string> data_lines = molecule.getDataItemContent(item_name);

  // If there is no data item, return an empty mask.
  if (data_lines.size() == 0LLU) {
    return AtomMask(ag, chemfe);
  }

  // Concatenate the strings, sans carriage returns if they have not already been removed.
  std::string all_data;
  for (size_t i = 0LLU; i < data_lines.size(); i++) {
    const size_t nchar = data_lines[i].size();
    for (size_t j = 0LLU; i < nchar; i++) {
      if (data_lines[i][j] == '\n') {
        data_lines[i][j] = ' ';
      }
    }
    all_data.append(data_lines[i]);

    // Add an extra space between data lines.  This should not disrupt any valid atom masks, and
    // will prevent lists of integers or atom names from having the final entry of one line be
    // concatenated as one number or word in front of the first entry on the next line.
    all_data += ' ';
  }
  
  // Try to construct an AtomMask.
  try {
    AtomMask result(all_data, ag, chemfe, cf.data());
    return result;
  }
  catch (std::runtime_error) {

    // If the data item content could not be processed as an atom mask, try parsing it as as list
    // of integers, which would be assumed to indicate atom indices (starting at 1, and thus
    // decremented for C / C++ array indexing).  The text is assumed not to contain comments or
    // quoted strings.  Commas will be counted as delimiters.
    const std::vector<std::string> words = separateText(all_data, {}, {}, { "\"", "\"", "'", "'" },
                                                        { "," });
    const int word_count = words.size();
    bool all_integer = true;
    for (int i = 0; i < word_count; i++) {
      all_integer = (all_integer && verifyNumberFormat(words.c_str(), NumberFormat::INTEGER));
    }
    AtomMask result(ag);
    if (all_integer == true) {
      for (
    }
  }
}

//-------------------------------------------------------------------------------------------------
AtomMask maskFromSdfDataItem(const std::string &item_name, const MdlMol &molecule,
                             const AtomGraph *ag) {
  const CoordinateFrame cf = molecule.exportCoordinateFrame();
  const ChemicalFeatures chemfe(ag, cf.data());
  maskFromSdfDataItem(item_name, molecule, ag, chemfe, cf);
}

} // namespace structure
} // namespace stormm
