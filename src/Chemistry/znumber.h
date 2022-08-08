// -*-c++-*-
#ifndef STORMM_ZNUMBER_H
#define STORMM_ZNUMBER_H

#include <vector>
#include "copyright.h"
#include "DataTypes/stormm_vector_types.h"

namespace stormm {
namespace chemistry {

/// \brief Obtain atomic Z-numbers for a list of atoms, assuming that the masses follow the
///        weighted natural abundances 
///
/// \param masses  Masses of the atoms in question
std::vector<int> massToZNumber(const std::vector<double> &masses);

/// \brief Convert a series of atomic numbers to element symbols from the periodic table
///
/// Overloaded:
///   - Obtain the symbol for a single Z-number
///   - Obtain symbols for a vector of Z-numbers
///
/// \param atomic_numbers  Z-numbers of the atoms in question (can in clude 0, for a virtual site)
/// \{
char2 zNumberToSymbol(const int atomic_number);
std::vector<char2> zNumberToSymbol(const std::vector<int> &atomic_numbers);
/// \}

/// \brief Convert a series of atomic symbols to elemental Z-numbers
///
/// \param atomic_symbols  Periodic table symbols of the elements in question
std::vector<int> symbolToZNumber(const std::vector<char2> &atomic_symbols);

} // namespace chemistry
} // namespace stormm

#endif

