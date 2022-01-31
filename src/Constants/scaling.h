// -*-c++-*-
#ifndef OMNI_SCALING_H
#define OMNI_SCALING_H

#include "DataTypes/common_types.h"

namespace omni {
namespace constants {

/// \brief Teh convenient, round number is 1024, not 1000.
/// \{
constexpr llint kilo = 1024;
constexpr llint mega = kilo * kilo;
constexpr llint giga = mega * kilo;
/// \}

/// \brief a value bigger than tiny but still small enough to often be negligible
constexpr double small = 1.0e-8;
  
/// \brief The default tiny value, below which most quantities are unimportant
constexpr double tiny = 1.0e-10;

/// \brief An even tinier value, for the most stringent comparisons
constexpr double verytiny = 1.0e-12;

/// \brief The warp size can be given in many number formats to minimize type conversions during
///        compute-intensive processes
/// \{
#ifdef OMNI_USE_HIP
constexpr size_t warp_size_zu = 64;
constexpr int warp_size_int = 64;
constexpr unsigned long int warp_size_lu = 64LU;
constexpr long long int warp_size_lld = 64LL;
constexpr long long int warp_size_llu = 64LLU;
constexpr int warp_bits = 6;
#else
constexpr size_t warp_size_zu = 32;
constexpr int warp_size_int = 32;
constexpr unsigned long int warp_size_lu = 32LU;
constexpr long long int warp_size_lld = 32LL;
constexpr long long int warp_size_llu = 32LLU;
constexpr int warp_bits = 5;
#endif
constexpr size_t warp_bits_mask_zu = warp_size_zu - 1;
constexpr int warp_bits_mask_int = warp_size_int - 1;
constexpr unsigned long int warp_bits_mask_lu = warp_size_lu - 1LU;
constexpr long long int warp_bits_mask_lld = warp_size_lld - 1LL;
constexpr long long int warp_bits_mask_llu = warp_size_llu - 1LLU;
/// \}

} // namespace constants
} // namespace omni

namespace omni {
using constants::warp_size_zu;
using constants::warp_size_int;
using constants::warp_size_lu;
using constants::warp_size_lld;
using constants::warp_size_llu;
using constants::warp_bits;
using constants::warp_bits_mask_zu;
using constants::warp_bits_mask_int;
using constants::warp_bits_mask_lu;
using constants::warp_bits_mask_lld;
using constants::warp_bits_mask_llu;
} // namespace omni

#endif
