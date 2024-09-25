// -*-c++-*-
#ifndef STORMM_HOSE_POPC_H
#define STORMM_HOSE_POPC_H

#include "copyright.h"
#include "DataTypes/common_types.h"

namespace stormm {
namespace numerics {

/// \brief Count the number of active bits in a short unsigned integer.
///
/// \param x  The integer to analyze
int hostPopcs(ushort x);

/// \brief Count the number of active bits in a long unsigned integer.
///
/// \param x  The integer to analyze
int hostPopc(uint x);

/// \brief Count the number of active bits in a long long unsigned integer.
///
/// \param x  The integer to analyze
int hostPopcll(ullint x);

/// \brief Find the position of the first bit set to one in a short unsigned integer.
///
/// \param x  The integer to analyze
int hostFfss(ushort x);

/// \brief Find the position of the first bit set to one in a long unsigned integer.
///
/// \param x  The integer to analyze
int hostFfs(uint x);

/// \brief Find the position of the first bit set to one in a long long unsigned integer.
///
/// \param x  The integer to analyze
int hostFfsll(ullint x);

} // namespace numerics
} // namespace stormm

#endif
