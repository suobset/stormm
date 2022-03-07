// -*-c++-*-
#ifndef OMNI_TOPOLOGY_BOUNDS_CHECKS_H
#define OMNI_TOPOLOGY_BOUNDS_CHECKS_H

namespace omni {
namespace topology {

/// \brief Check an atom index against the available content, reporting an error if the bounds are
///        exceeded.
///
/// \param index          Index of the atom of interest
/// \param max_atoms      The maximum number of available atoms
/// \param class_caller   Name of the calling function or class, to be included in error messages
/// \param method_caller  Name of a class member function calling this routine, if applicable
void atomValidityCheck(int index, int max_atoms, const char* class_caller = nullptr,
                       const char* method_caller = nullptr);

/// \brief Check an atom range against the available content, reporting an error if the bounds are
///        exceeded.
///
/// \param low_index      Atom index at the left end of the range of interest
/// \param high_index     The right end of the range of interest (this can be the maximum number of
///                       atoms in the system--the group is assumed to be [low_index, high_index)
/// \param max_atoms      The maximum number of available atoms
/// \param class_caller   Name of the calling function or class, to be included in error messages
/// \param method_caller  Name of a class member function calling this routine, if applicable
void atomValidityCheck(int low_index, int high_index, int max_atoms,
                       const char* class_caller = nullptr, const char* method_caller = nullptr);

/// \brief Check a residue index against the available content, reporting an error if the bounds
///        are exceeded.
///
/// \param index          Index of the residue of interest
/// \param max_atoms      The maximum number of available residues
/// \param class_caller   Name of the calling function or class, to be included in error messages
/// \param method_caller  Name of a class member function calling this routine, if applicable
void residueValidityCheck(int index, int max_residues, const char* class_caller = nullptr,
                          const char* method_caller = nullptr);

/// \brief Check an residue range against the available content, reporting an error if the bounds
///        are exceeded.
///
/// \param low_index      Residue index at the left end of the range of interest
/// \param high_index     The right end of the range of interest (this can be the systems' maximum
///                       number of residues--the group is assumed to be [low_index, high_index)
/// \param max_residues   The maximum number of available residues
/// \param class_caller   Name of the calling function or class, to be included in error messages
/// \param method_caller  Name of a class member function calling this routine, if applicable
void residueValidityCheck(int low_index, int high_index, int max_residues,
                          const char* class_caller = nullptr, const char* method_caller = nullptr);
                          

/// \brief Check a molecule index against the available content, reporting an error if the bounds
///        are exceeded.
///
/// \param index          Index of the molecule of interest
/// \param max_atoms      The maximum number of available molecules
/// \param class_caller   Name of the calling function or class, to be included in error messages
/// \param method_caller  Name of a class member function calling this routine, if applicable
void moleculeValidityCheck(int index, int max_molecules, const char* class_caller = nullptr,
                           const char* method_caller = nullptr);

/// \brief Check an molecule range against the available content, reporting an error if the bounds
///        are exceeded.
///
/// \param low_index      Molecule index at the left end of the range of interest
/// \param high_index     The right end of the range of interest (this can be the system's maximum
///                       number of molecules--the group is assumed to be [low_index, high_index)
/// \param max_molecules  The maximum number of available molecules
/// \param class_caller   Name of the calling function or class, to be included in error messages
/// \param method_caller  Name of a class member function calling this routine, if applicable
void moleculeValidityCheck(int low_index, int high_index, int max_molecules,
                           const char* class_caller = nullptr,
                           const char* method_caller = nullptr);

} // namespace topology
} // namespace omni

#endif
