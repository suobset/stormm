#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include "Constants/symbol_values.h"
#include "Parsing/polynumeric.h"
#include "Reporting/error_format.h"
#include "random.h"

namespace omni {
namespace random {

using parse::PolyNumeric;

//-------------------------------------------------------------------------------------------------
Ran2Generator::Ran2Generator(const int igseed) :
  iunstable_a{-1},
  iunstable_b{123456789},
  state_sample{0},
  state_vector{0}
{
  // Initialize by calling the uniform generator
  uniformRandomNumber();

  // Compute additional random numbers based on a pseudo-random seed.  In order to make this
  // generator diverge from others based on different seeds, we must first run through about
  // sixteen cycles of the uniform generator in order to run through correlated results left
  // over after the common initialization.  The entire state vector needs to change in order
  // to produce a decorrelated random number.
  iunstable_a = igseed;
  for (int i = 0; i < 15; i++) {
    uniformRandomNumber();
  }
}

//-------------------------------------------------------------------------------------------------
double Ran2Generator::uniformRandomNumber() {

  // Local variable to avoid dereferencing the self pointer
  int lcl_iunstbl_a  = iunstable_a;
  int lcl_iunstbl_b  = iunstable_b;
  int lcl_state_sample = state_sample;

  // Populate the state vector
  int unstbl_quotient;
  if (lcl_iunstbl_a <= 0) {
    lcl_iunstbl_a = (-lcl_iunstbl_a < 1) ? 1 : -lcl_iunstbl_a;
    lcl_iunstbl_b = lcl_iunstbl_a;
    for (int j = ntab + 7; j >= 0 ; j--) {
      unstbl_quotient = lcl_iunstbl_a / iq1;
      lcl_iunstbl_a = ia1 * (lcl_iunstbl_a - unstbl_quotient * iq1) - (unstbl_quotient * ir1);
      lcl_iunstbl_a += (lcl_iunstbl_a < 0) * im1;
      state_vector[j] = lcl_iunstbl_a;
    }
    lcl_state_sample = state_vector[0];
  }
  unstbl_quotient = lcl_iunstbl_a / iq1;
  lcl_iunstbl_a = ia1 * (lcl_iunstbl_a - unstbl_quotient * iq1) - unstbl_quotient*ir1;
  lcl_iunstbl_a += (lcl_iunstbl_a < 0) * im1;
  unstbl_quotient = lcl_iunstbl_b / iq2;
  lcl_iunstbl_b = ia2 * (lcl_iunstbl_b - unstbl_quotient * iq2) - unstbl_quotient * ir2;
  lcl_iunstbl_b += (lcl_iunstbl_b < 0) * im2;
  int swap_index = lcl_state_sample / ndiv;
  lcl_state_sample = state_vector[swap_index] - lcl_iunstbl_b;
  state_vector[swap_index] = lcl_iunstbl_a;

  // Generate the alternative number so as to return the lesser of this or rand2_max
  lcl_state_sample += (lcl_state_sample < 1) * im1_minus_one;
  double temp = am * lcl_state_sample;

  // Copy lcl values back to modifiable struct member variable
  iunstable_a = lcl_iunstbl_a;
  iunstable_b = lcl_iunstbl_b;
  state_sample = lcl_state_sample;

  return std::min(ran2_max, temp);
}

//-------------------------------------------------------------------------------------------------
double Ran2Generator::gaussianRandomNumber() {

  using symbols::twopi;

  const double x1 = std::sqrt(-2.0 * std::log(uniformRandomNumber()));
  const double x2 = std::sin(twopi * uniformRandomNumber());
  return x1 * x2;
}

//-------------------------------------------------------------------------------------------------
Xoroshiro128pGenerator::Xoroshiro128pGenerator(const int igseed, const int niter) :
    state{seed128(igseed)}
{
  // Trap cases where ullint is not 64 bit (this should never be an issue)
  if (sizeof(ullint) != 8) {
    rtErr("The size of an unsigned long long integer (ullint) must be 64 bits in order to "
          "instantiate the xoroshiro128+ random number generator.", "Xoroshiro128pGenerator");
  }

  // Trap cases where double is not 64 bit or float is not 32 bit (this should never be an issue)
  if (sizeof(double) != 8) {
    rtErr("The size of a double-precision real number (double) must by 64 bits in order to "
          "instantiate the xoroshiro128+ random number generator.", "Xoroshiro128pGenerator");
  }
  if (sizeof(float) != 4) {
    rtErr("The size of a single-precision real number (float) must by 32 bits in order to "
          "instantiate the xoroshiro128+ random number generator.", "Xoroshiro128pGenerator");
  }

  // Run some iterations to get the bit occupancy to roughly half zeros and half ones
  for (int i = 0; i < niter; i++) {
    next();
  }
}

//-------------------------------------------------------------------------------------------------
Xoroshiro128pGenerator::Xoroshiro128pGenerator(const ullint2 state_in) :
  state{state_in}
{}
    
//-------------------------------------------------------------------------------------------------
double Xoroshiro128pGenerator::uniformRandomNumber() {
  const ullint rndbits = next();
  PolyNumeric work;
  work.ulli = (((rndbits >> 12) & 0xfffffffffffffLLU) | 0x3ff0000000000000LLU);
  return work.d - 1.0;
}

//-------------------------------------------------------------------------------------------------
double Xoroshiro128pGenerator::gaussianRandomNumber() {

  using symbols::twopi;

  const double x1 = std::sqrt(-2.0 * std::log(uniformRandomNumber()));
  const double x2 = std::sin(twopi * uniformRandomNumber());
  return x1 * x2;
}

//-------------------------------------------------------------------------------------------------
void Xoroshiro128pGenerator::jump() {
  const ullint2 stride = { xrs128p_jump_i, xrs128p_jump_ii };
  fastForward(stride);
}

//-------------------------------------------------------------------------------------------------
void Xoroshiro128pGenerator::longJump() {
  const ullint2 stride = { xrs128p_longjump_i, xrs128p_longjump_ii };
  fastForward(stride);
}

//-------------------------------------------------------------------------------------------------
ullint2 Xoroshiro128pGenerator::revealState() const {
  return state;
}

//-------------------------------------------------------------------------------------------------
ullint2 Xoroshiro128pGenerator::seed128(const int igseed) {

  // It is within the C and C++ standards to set an unsigned integer equal to a signed integer,
  // even a negative value.  It changes the interpretation, no more.
  const ullint uiseed = igseed;
  const int nbits = sizeof(uint) * 8;
  const int nfill = sizeof(ullint) * 8;
  const int nspc  = nfill / nbits;
  int offset = 0;
  ullint2 result = { 0LLU, 0LLU };
  for (int i = 0; i < nbits; i++) {
    result.x |= (((uiseed >> i) & 0x1LLU) << ((nspc * i) + offset));
    offset++;
    offset *= (offset != nspc);
    result.y |= (((uiseed >> i) & 0x1LLU) << ((nspc * i) + offset));
    offset++;
    offset *= (offset != nspc);
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
ullint Xoroshiro128pGenerator::next() {
  const ullint s0 = state.x;
  ullint       s1 = state.y;
  const ullint result = s0 + s1;
  s1 ^= s0;
  state.x = (((s0 << 24) | (s0 >> (64 - 24))) ^ s1 ^ (s1 << 16));
  state.y =  ((s1 << 37) | (s1 >> (64 - 37)));
  return result;
}

//-------------------------------------------------------------------------------------------------
void Xoroshiro128pGenerator::fastForward(const ullint2 stride) {
  ullint s0 = 0LLU;
  ullint s1 = 0LLU;
  for (int b = 0; b < 64; b++) {
    if (stride.x & (0x1LLU << b)) {
      s0 ^= state.x;
      s1 ^= state.y;
    }
    next();
  }
  for (int b = 0; b < 64; b++) {
    if (stride.y & (0x1LLU << b)) {
      s0 ^= state.x;
      s1 ^= state.y;
    }
    next();
  }
  state.x = s0;
  state.y = s1;
}

//-------------------------------------------------------------------------------------------------
Xoshiro256ppGenerator::Xoshiro256ppGenerator(const int igseed, const int niter) :
    state{seed256(igseed)}
{
  // Trap cases where ullint is not 64 bit (this should never be an issue)
  if (sizeof(ullint) != 8) {
    rtErr("The size of an unsigned long long integer (ullint) must be 64 bits in order to "
          "instantiate the xoshiro256++ random number generator.", "Xoshiro256ppGenerator");
  }

  // Trap cases where double is not 64 bit or float is not 32 bit (this should never be an issue)
  if (sizeof(double) != 8) {
    rtErr("The size of a double-precision real number (double) must by 64 bits in order to "
          "instantiate the xoroshiro256++ random number generator.", "Xoshiro256ppGenerator");
  }
  if (sizeof(float) != 4) {
    rtErr("The size of a single-precision real number (float) must by 32 bits in order to "
          "instantiate the xoroshiro256++ random number generator.", "Xoshiro256ppGenerator");
  }

  // Run some iterations to get the bit occupancy to roughly half zeros and half ones
  for (int i = 0; i < niter; i++) {
    next();
  }
}

//-------------------------------------------------------------------------------------------------
Xoshiro256ppGenerator::Xoshiro256ppGenerator(const ullint4 state_in) :
    state{state_in}
{}

//-------------------------------------------------------------------------------------------------
double Xoshiro256ppGenerator::uniformRandomNumber() {
  const ullint rndbits = next();
  PolyNumeric work;
  work.ulli = (((rndbits >> 12) & 0xfffffffffffffLLU) | 0x3ff0000000000000LLU);
  return work.d - 1.0;
}

//-------------------------------------------------------------------------------------------------
double Xoshiro256ppGenerator::gaussianRandomNumber() {

  using symbols::twopi;

  const double x1 = std::sqrt(-2.0 * std::log(uniformRandomNumber()));
  const double x2 = std::sin(twopi * uniformRandomNumber());
  return x1 * x2;
}

//-------------------------------------------------------------------------------------------------
void Xoshiro256ppGenerator::jump() {
  const ullint4 stride = { xrs256pp_jump_i,   xrs256pp_jump_ii,
                           xrs256pp_jump_iii, xrs256pp_jump_iv };
  fastForward(stride);
}

//-------------------------------------------------------------------------------------------------
void Xoshiro256ppGenerator::longJump() {
  const ullint4 stride = { xrs256pp_longjump_i,   xrs256pp_longjump_ii,
                           xrs256pp_longjump_iii, xrs256pp_longjump_iv };
  fastForward(stride);
}

//-------------------------------------------------------------------------------------------------
ullint4 Xoshiro256ppGenerator::revealState() const {
  return state;
}

//-------------------------------------------------------------------------------------------------
ullint4 Xoshiro256ppGenerator::seed256(const int igseed) {

  // It is within the C and C++ standards to set an unsigned integer equal to a signed integer,
  // even a negative value.  It changes the interpretation, no more.
  const ullint uiseed = igseed;
  const ullint ujseed = (igseed ^ -1);
  const int nbits = sizeof(uint) * 8;
  const int nfill = sizeof(ullint) * 8;
  const int nspc  = nfill / nbits;
  ullint4 result = { 0LLU, 0LLU, 0LLU, 0LLU };
  for (int i = 0; i < nbits; i++) {
    result.x |= (((uiseed >> i) & 0x1LLU) << ((nspc * i)));
    result.y |= (((uiseed >> i) & 0x1LLU) << ((nspc * i) + 1));
    result.z |= (((ujseed >> i) & 0x1LLU) << ((nspc * i)));
    result.w |= (((ujseed >> i) & 0x1LLU) << ((nspc * i) + 1));
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
ullint Xoshiro256ppGenerator::next() {
  const ullint sxsw = state.x + state.w;
  const ullint result = state.x + ((sxsw << 23) | (sxsw >> (64 - 23)));
  const ullint t = (state.y << 17);
  state.z ^= state.x;
  state.w ^= state.y;
  state.y ^= state.z;
  state.x ^= state.w;
  state.z ^= t;
  state.w = ((state.w << 45) | (state.w >> (64 - 45)));
  return result;
}

//-------------------------------------------------------------------------------------------------
void Xoshiro256ppGenerator::fastForward(ullint4 stride) {
  ullint s0 = 0LLU;
  ullint s1 = 0LLU;
  ullint s2 = 0LLU;
  ullint s3 = 0LLU;
  for (int b = 0; b < 64; b++) {
    if (stride.x & (0x1LLU << b)) {
      s0 ^= state.x;
      s1 ^= state.y;
      s2 ^= state.z;
      s3 ^= state.w;
    }
    next();
  }
  for (int b = 0; b < 64; b++) {
    if (stride.y & (0x1LLU << b)) {
      s0 ^= state.x;
      s1 ^= state.y;
      s2 ^= state.z;
      s3 ^= state.w;
    }
    next();
  }
  for (int b = 0; b < 64; b++) {
    if (stride.z & (0x1LLU << b)) {
      s0 ^= state.x;
      s1 ^= state.y;
      s2 ^= state.z;
      s3 ^= state.w;
    }
    next();
  }
  for (int b = 0; b < 64; b++) {
    if (stride.w & (0x1LLU << b)) {
      s0 ^= state.x;
      s1 ^= state.y;
      s2 ^= state.z;
      s3 ^= state.w;
    }
    next();
  }
  state.x = s0;
  state.y = s1;
  state.z = s2;
  state.w = s3;
}

} // namespace random
} // namespace omni
