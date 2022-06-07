// -*-c++-*-
#ifndef OMNI_RANDOM_H
#define OMNI_RANDOM_H

#include "Accelerator/hybrid.h"
#include "Constants/scaling.h"
#include "DataTypes/common_types.h"
#include "DataTypes/omni_vector_types.h"
#include "Math/rounding.h"

namespace omni {
namespace random {

using card::Hybrid;
using card::HybridTargetLevel;
using constants::giga_zu;
  
/// \brief Constants for the "ran2" pseudo-random number generator
/// \{
constexpr int im1 = 2147483563;
constexpr int im2 = 2147483399;
constexpr double am = 1.0 / im1;
constexpr int im1_minus_one = im1 - 1;
constexpr int ia1 = 40014;
constexpr int ia2 = 40692;
constexpr int iq1 = 53668;
constexpr int iq2 = 52774;
constexpr int ir1 = 12211;
constexpr int ir2 = 3791;
constexpr int ntab = 32;
constexpr int ndiv = 1 + im1_minus_one / ntab;
constexpr double double_increment = 2.2205e-16;
constexpr double ran2_max = 1.0 - double_increment;
/// \}
  
/// \brief Stores the state of a Ran2 pseudo-random number generator.  Member functions produce
///        random numbers along various distributions, as required.  While it is not as performant
///        to have a member function, these random number generators are intended for convenience
///        and unit testing purposes.  Developers that wish to use higher-performance random number
///        generators should use the Xoroshiro128Generator (below) or link other libraries for the
///        C++ layer (i.e. Intel MKL) or libraries that run on the GPU layer (i.e. cuRAND). 
class Ran2Generator {
public:

  /// \brief Initialize the unit testing on-board pseudo-random number generator.
  ///
  /// \param igseed  The pseudo-randome number seed
  Ran2Generator(int igseed = 827493);  

  /// \brief Return a single random number distributed over a uniform distribution.  This is an
  ///        internal generator based on an integer vector state which will provide reproducible
  ///        results across platforms.
  double uniformRandomNumber();

  /// \brief Return a normally distributed random number.  This works off of the uniform random
  ///        number generator and will thus advance the state vector of the RandomNumberGenerator
  ///        that produces the result.
  double gaussianRandomNumber();

private:
  int iunstable_a;
  int iunstable_b;
  int state_sample;
  int state_vector[ntab + 8];
};

/// \brief Constants for Xoroshiro128+ generator.
/// \{
constexpr ullint xrs128p_jump_i      = 0xdf900294d8f554a5LLU;
constexpr ullint xrs128p_jump_ii     = 0x170865df4b3201fcLLU;
constexpr ullint xrs128p_longjump_i  = 0xd2a98b26625eee7bLLU;
constexpr ullint xrs128p_longjump_ii = 0xdddf9b1090aa7ac1LLU;
/// \}

/// \brief Constants for Xoshiro256++ generator.
/// \{
constexpr ullint xrs256pp_jump_i       = 0x180ec6d33cfd0abaLLU;
constexpr ullint xrs256pp_jump_ii      = 0xd5a61266f0c9392cLLU;
constexpr ullint xrs256pp_jump_iii     = 0xa9582618e03fc9aaLLU;
constexpr ullint xrs256pp_jump_iv      = 0x39abdc4529b1661cLLU;
constexpr ullint xrs256pp_longjump_i   = 0x76e15d3efefdcbbfLLU;
constexpr ullint xrs256pp_longjump_ii  = 0xc5004e441c522fb3LLU;
constexpr ullint xrs256pp_longjump_iii = 0x77710069854ee241LLU;
constexpr ullint xrs256pp_longjump_iv  = 0x39109bb02acbe635LLU;
/// \}

/// \brief The "Xoroshiro128+" random number generator.  It's decent, but not recommended for
///        situations where a quarter million or more streams are producing random number sequences
///        in unison.  Those streams can end up containing correlated patterns.  This generator has
///        been shown to fail BigCrush.
class Xoroshiro128pGenerator {
public:

  /// \brief The constructor can start or restart the generator.
  ///
  /// Overloaded:
  ///   - Initialize the generator upon construction.
  ///   - Take the given state.
  ///
  /// \param igseed    The pseudo-randome number seed
  /// \param niter     Number of iterations to run through while equilibrating the generator
  /// \param state_in  State to accept
  Xoroshiro128pGenerator(int igseed = 827493, int niter = 25);
  Xoroshiro128pGenerator(const ullint2 state_in);
  
  /// \brief Return a single random number distributed over a uniform distribution.  This is an
  ///        internal generator based on an integer vector state which will provide reproducible
  ///        results across platforms.
  double uniformRandomNumber();

  /// \brief Return a normally distributed random number.  This works off of the uniform random
  ///        number generator and will thus advance the state vector of the RandomNumberGenerator
  ///        that produces the result.
  double gaussianRandomNumber();

  /// \brief Jump forward 2^64 iterations in the sequence.
  void jump();
  
  /// \brief Jump forward 2^96 iterations in the sequence.
  void longJump();

  /// \brief Reveal the current state of the generator.
  ullint2 revealState() const;

  /// \brief Set the current state of the generator.
  void setState();
  
private:
  ullint2 state;  ///< 128-bit state vector for the generator
  
  /// \brief Transform a 32-bit int into 128 well-spaced bits for the seed
  ///
  /// \param igseed  Random seed passed down from the constructor
  ullint2 seed128(int igseed);
  
  /// \brief Iterate to the next random number
  ullint next();

  /// \brief Jump forward in the sequence according to the jump or long-jump perturbations
  ///
  /// \param stride  The perturbation to use.  This will be either xrs128p_jump_(i,ii) or
  ///                xrs128p_longjump_(i,ii).
  void fastForward(const ullint2 stride);
};

/// \brief An series of "Xorshiro128+" generators, with state vectors for all of them and the
///        means for seeding the series based on long jumps from a single state vector.
template <typename T> class Xoroshiro128pSeries {
public:

  /// \brief The constructor can work off of a simple random initialization seed integer or a
  ///        specific state vector for the first generator in the series.
  ///
  /// \{
  Xoroshiro128pSeries(size_t generators_in = 1LLU, size_t depth_in = 2LLU, int igseed_in = 827493,
                      int niter = 25, size_t bank_limit = 4LLU * giga);

  Xoroshiro128pSeries(const ullint2 state_in, int generators_in = 0);
  /// \}
  
private:
  size_t generators;       ///< The quantity of (pseudo-) random number generators in the series
  size_t depth;            ///< The depth of each generator's random numbers in the memory bank
  int igseed;              ///< Integer used to seed this generator series, for records-keeping
                           ///<   purposes
  Hybrid<ullint2> states;  ///< Thestate vectors for all random number generators in the series
  Hybrid<T> bank;          ///< Bank of random numbers pre-computed and saved for later use.  The
                           ///<   numbers are stored only in the specifictype format of the series,
                           ///<   float or double, to save space and avoid ambiguities as to when
                           ///<   the numbers should be refreshed if multiple levels of precision
                           ///<   in the random numbers are required.
};
  
/// \brief The "Xoshiro256++" random number generator.  While not cryptographically useful, it is
///        a rock-solid random number generator for both floating-point and 64-bit integer results.
class Xoshiro256ppGenerator {
public:

  /// \brief The constructor can start or restart the generator.
  ///
  /// Overloaded:
  ///   - Initialize the generator upon construction.
  ///   - Take the given state.
  ///
  /// \param igseed    The pseudo-randome number seed
  /// \param niter     Number of iterations to run through while equilibrating the generator
  /// \param state_in  State to accept
  Xoshiro256ppGenerator(int igseed = 827493, int niter = 50);
  Xoshiro256ppGenerator(const ullint4 state_in);

  /// \brief Return a single random number distributed over a uniform distribution.  This is an
  ///        internal generator based on an integer vector state which will provide reproducible
  ///        results across platforms.
  double uniformRandomNumber();

  /// \brief Return a normally distributed random number.  This works off of the uniform random
  ///        number generator and will thus advance the state vector of the RandomNumberGenerator
  ///        that produces the result.
  double gaussianRandomNumber();

  /// \brief Jump forward 2^64 iterations in the sequence.
  void jump();
  
  /// \brief Jump forward 2^96 iterations in the sequence.
  void longJump();

  /// \brief Reveal the current state of the generator.
  ullint4 revealState() const;

  /// \brief Reveal the random bit string.
  ullint revealBitString() const;
  
private:
  ullint4 state;  ///< 128-bit state vector for the generator
  
  /// \brief Transform a 32-bit int into 128 well-spaced bits for the seed
  ///
  /// \param igseed  Random seed passed down from the constructor
  ullint4 seed256(int igseed);
  
  /// \brief Iterate to the next random number
  ullint next();

  /// \brief Jump forward in the sequence according to the jump or long-jump perturbations
  ///
  /// \param stride  The perturbation to use.  This will be either xrs256pp_jump_(i,ii,iii,iv) or
  ///                xrs256pp_longjump_(i,ii,iii,iv).
  void fastForward(const ullint4 stride);
};
  
} // namespace random
} // namespace omni

#include "random.tpp"

#endif
