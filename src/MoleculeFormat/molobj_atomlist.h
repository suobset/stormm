// -*-c++-*-
#ifndef STORMM_MOLOBJ_ATOMLIST_H
#define STORMM_MOLOBJ_ATOMLIST_H

#include <string>
#include <vector>
#include "copyright.h"
#include "Parsing/textfile.h"

namespace stormm {
namespace structure {

using parse::TextFile;

/// \brief An atom list entry (this object can be assembled either from one of the deprecated
///        V2000 format lines after the bonds block, or from one of the "M  ALS" properties)
class MolObjAtomList {
public:

  /// \brief The constructor can take all member variables (and all come with default values to
  ///        let this form of the constructor serve as the blank object constructor), or a pointer
  ///        to the line of a text file from which the information shall come.
  ///
  /// \param tf           Text of the original .sdf or .mol file, read into RAM
  /// \param line_number  Number of the line on which to read the data
  /// \param title        The title of the structure, if known, for error tracing purposes
  /// \{
  MolObjAtomList(const std::vector<int> &atomic_numbers_in = {}, bool exclusions_in = false,
                 int atom_attachment_in = 0);

  MolObjAtomList(const TextFile &tf, int line_number, const std::string &title = std::string(""));
  /// \}

private:
  int entry_count;                  ///< The number of atomic elements (identified by Z-numbers)
                                    ///<   in this list
  std::vector<int> atomic_numbers;  ///< Atomic (Z-) numbers of atoms that are to be excluded or
                                    ///<   included by processing this list
  bool exclusions;                  ///< Indicate whether this list covers atomic numbers which
                                    ///<   are to be excluded (TRUE) or included (FALSE)
  int atom_attachment;              ///< Attachment point of the list, an index of an atom in the
                                    ///<   molecule itself
};

} // namespace structure
} // namespace stormm

#endif
