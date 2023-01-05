// -*-c++-*-
#ifndef STORMM_RANDOM_H
#define STORMM_RANDOM_H

#include <vector>
#include "copyright.h"
#include "Accelerator/hybrid.h"
#include "Constants/scaling.h"
#include "Constants/hpc_bounds.h"
#include "DataTypes/common_types.h"
#include "DataTypes/stormm_vector_types.h"
#include "Math/rounding.h"

namespace stormm {
namespace random {

using card::Hybrid;
using card::HybridKind;
using card::HybridTargetLevel;
using constants::giga_zu;
using data_types::isFloatingPointScalarType;
using math::roundUp;

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

/// \brief Enumerate the types of random numbers that can be generated.
enum class RandomNumberKind {
  UNIFORM,   ///< Uniform random number distribution
  GAUSSIAN   ///< Normal distribution of random numbers
};

/// \brief List the random number generator types that can power a RandomNumberMill object.
enum class RandomAlgorithm {
  XOROSHIRO_128P,  ///< Xoroshiro128+ generator (see below, fails BigCrush and not advised for
                   ///<   powering mills with > 1024 generator streams)
  XOSHIRO_256PP    ///< Xoshiro256++ generator (see below--high quality generator)
};

/// \brief Define the order in which various random number generators will fill a matrix.
enum class RngFillMode {
  COLUMNS,  ///< Fill the column-major matrix one column at a time
  ROWS      ///< Fill the column-major matrix in row-major order, one row at a time
};
  
/// \brief The default random seed for STORMM
constexpr int default_random_seed = 827493;

/// \brief Numbers of "scrub cycles" expected to sufficiently randomize non-zero seeds in any of
///        the various XOR-shift generators.
/// \{
constexpr int default_xoroshiro128p_scrub = 25;
constexpr int default_xoshiro256pp_scrub  = 50;
/// \}
  
/// \brief The maximum number of long jumps to take with any of the scrambled linear pseudo-random
///        number generators.  The CPU will take the long jumps, while individual GPU threads make
///        standard jumps ahead in the pseudo-random sequence.  The CPU is 10-30 times faster
///        than a single GPU thread, so this number could in fact scale with the number of state
///        vectors being requested, but the optimization would be on the order of milliseconds
///        and making this a constant makes it easier to reproduce the sequence at a given point.
constexpr int max_xo_long_jumps = 1024;

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
  Ran2Generator(int igseed = default_random_seed);  

  /// \brief Return a single random number distributed over a uniform distribution.  This is an
  ///        internal generator based on an integer vector state which will provide reproducible
  ///        results across platforms.
  ///
  /// Overloaded:
  ///   - Issue a single random number
  ///   - Issue multiple random numbers as a std::vector<double>
  ///   - Issue random numbers to fill a Standard Template Library vector of arbitrary type (these
  ///     random numbers will still be issued as double-precision, using a procedure that draws on
  ///     53 bits of the 64-bit string)
  ///
  /// \param count     The quantity of random numbers to produce
  /// \param rows      The total rows of a matrix to fill up with random numbers
  /// \param columns   The total columns of a matrix to fill up with random numbers
  /// \param xv        Pre-allocated array to fill with random numbers
  /// \param fp_scale  Fixed precision scaling factor to apply when filling arrays of integer type
  /// \{
  double uniformRandomNumber();
  std::vector<double> uniformRandomNumber(size_t count);
  std::vector<double> uniformRandomNumber(size_t rows, size_t columns,
                                          RngFillMode mode = RngFillMode::COLUMNS);

  template <typename T>
  void uniformRandomNumber(std::vector<T> *xv, double fp_scale = 1.0);
  /// \}

  /// \brief Return a normally distributed random number.  This works off of the uniform random
  ///        number generator and will thus advance the state vector of the RandomNumberGenerator
  ///        that produces the result.  This functon is overloaded in a manner similar to the
  ///        uniformRandomNumber() member function (see above) and accepts input parameters with
  ///        similar descriptions.
  /// \{
  double gaussianRandomNumber();
  std::vector<double> gaussianRandomNumber(size_t count);
  std::vector<double> gaussianRandomNumber(size_t rows, size_t columns,
                                           RngFillMode mode = RngFillMode::COLUMNS);

  template <typename T>
  void gaussianRandomNumber(std::vector<T> *xv, double fp_scale = 1.0);
  /// \}

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
  /// \{
  Xoroshiro128pGenerator(int igseed = default_random_seed, int niter = 25);
  Xoroshiro128pGenerator(const ullint2 state_in);
  /// \}
  
  /// \brief Return a single random number distributed over a uniform distribution.  This is an
  ///        internal generator based on an integer vector state which will provide reproducible
  ///        results across platforms.
  ///
  /// Overloaded:
  ///   - Issue a single random number
  ///   - Issue multiple random numbers as a std::vector<double>
  ///   - Issue random numbers to fill a Standard Template Library vector of arbitrary type (these
  ///     random numbers will still be issued as double-precision, using a procedure that draws on
  ///     53 bits of the 64-bit string)
  ///
  /// \param count     The quantity of random numbers to produce
  /// \param rows      The total rows of a matrix to fill up with random numbers
  /// \param columns   The total columns of a matrix to fill up with random numbers
  /// \param xv        Pre-allocated array to fill with random numbers
  /// \param fp_scale  Fixed precision scaling factor to apply when filling arrays of integer type
  /// \{
  double uniformRandomNumber();
  std::vector<double> uniformRandomNumber(size_t count);
  std::vector<double> uniformRandomNumber(size_t rows, size_t columns,
                                          RngFillMode mode = RngFillMode::COLUMNS);

  template <typename T>
  void uniformRandomNumber(std::vector<T> *xv, double fp_scale = 1.0);
  /// \}

  /// \brief Return a single-precision random number obtained from a uniform distribution.
  ///        Various overloads as well as descriptions for their parameters follow from the
  ///        uniformRandomNumber() member function, above.
  /// \{
  float spUniformRandomNumber();
  std::vector<float> spUniformRandomNumber(size_t count);
  std::vector<float> spUniformRandomNumber(size_t rows, size_t columns,
                                           RngFillMode mode = RngFillMode::COLUMNS);

  template <typename T>
  void spUniformRandomNumber(std::vector<T> *xv, float fp_scale = 1.0);
  /// \}

  /// \brief Return a normally distributed random number.  This works off of the uniform random
  ///        number generator and will thus advance the state vector of the RandomNumberGenerator
  ///        that produces the result.  This functon is overloaded in a manner similar to the
  ///        uniformRandomNumber() member function (see above) and accepts input parameters with
  ///        similar descriptions.
  /// \{
  double gaussianRandomNumber();
  std::vector<double> gaussianRandomNumber(size_t count);
  std::vector<double> gaussianRandomNumber(size_t rows, size_t columns,
                                           RngFillMode mode = RngFillMode::COLUMNS);

  template <typename T>
  void gaussianRandomNumber(std::vector<T> *xv, double fp_scale = 1.0);
  /// \}

  /// \brief Return a single-precision random number obtained from a normal distribution.  This
  ///        functon is overloaded in a manner similar to the uniformRandomNumber() member function
  ///        (see above) and accepts input parameters with similar descriptions.
  /// \{
  float spGaussianRandomNumber();
  std::vector<float> spGaussianRandomNumber(size_t count);
  std::vector<float> spGaussianRandomNumber(size_t rows, size_t columns,
                                            RngFillMode mode = RngFillMode::COLUMNS);

  template <typename T>
  void spGaussianRandomNumber(std::vector<T> *xv, float fp_scale = 1.0);
  /// \}

  /// \brief Jump forward 2^64 iterations in the sequence.
  void jump();
  
  /// \brief Jump forward 2^96 iterations in the sequence.
  void longJump();

  /// \brief Reveal the current state of the generator.
  ullint2 revealState() const;

  /// \brief Reveal the random bit string.
  ullint revealBitString() const;

  /// \brief Set the current state of the generator.
  void setState(const ullint2 state_in);
  
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
  Xoshiro256ppGenerator(int igseed = default_random_seed, int niter = 50);
  Xoshiro256ppGenerator(const ullint4 state_in);

  /// \brief Return a single random number distributed over a uniform distribution.  This is an
  ///        internal generator based on an integer vector state which will provide reproducible
  ///        results across platforms.
  ///
  /// Overloaded:
  ///   - Issue a single random number
  ///   - Issue multiple random numbers as a std::vector<double>
  ///   - Issue random numbers to fill a Standard Template Library vector of arbitrary type (these
  ///     random numbers will still be issued as double-precision, using a procedure that draws on
  ///     53 bits of the 64-bit string)
  ///
  /// \param count     The quantity of random numbers to produce
  /// \param rows      The total rows of a matrix to fill up with random numbers
  /// \param columns   The total columns of a matrix to fill up with random numbers
  /// \param xv        Pre-allocated array to fill with random numbers
  /// \param fp_scale  Fixed precision scaling factor to apply when filling arrays of integer type
  /// \{
  double uniformRandomNumber();
  std::vector<double> uniformRandomNumber(size_t count);
  std::vector<double> uniformRandomNumber(size_t rows, size_t columns,
                                          RngFillMode mode = RngFillMode::COLUMNS);

  template <typename T>
  void uniformRandomNumber(std::vector<T> *xv, double fp_scale = 1.0);
  /// \}

  /// \brief Return a single-precision random number obtained from a uniform distribution.
  ///        Various overloads as well as descriptions for their parameters follow from the
  ///        uniformRandomNumber() member function, above.
  /// \{
  float spUniformRandomNumber();
  std::vector<float> spUniformRandomNumber(size_t count);
  std::vector<float> spUniformRandomNumber(size_t rows, size_t columns,
                                           RngFillMode mode = RngFillMode::COLUMNS);

  template <typename T>
  void spUniformRandomNumber(std::vector<T> *xv, float fp_scale = 1.0);
  /// \}

  /// \brief Return a normally distributed random number.  This works off of the uniform random
  ///        number generator and will thus advance the state vector of the RandomNumberGenerator
  ///        that produces the result.  This functon is overloaded in a manner similar to the
  ///        uniformRandomNumber() member function (see above) and accepts input parameters with
  ///        similar descriptions.
  /// \{
  double gaussianRandomNumber();
  std::vector<double> gaussianRandomNumber(size_t count);
  std::vector<double> gaussianRandomNumber(size_t rows, size_t columns,
                                           RngFillMode mode = RngFillMode::COLUMNS);

  template <typename T>
  void gaussianRandomNumber(std::vector<T> *xv, double fp_scale = 1.0);
  /// \}

  /// \brief Return a single-precision random number obtained from a normal distribution.  This
  ///        functon is overloaded in a manner similar to the uniformRandomNumber() member function
  ///        (see above) and accepts input parameters with similar descriptions.
  /// \{
  float spGaussianRandomNumber();
  std::vector<float> spGaussianRandomNumber(size_t count);
  std::vector<float> spGaussianRandomNumber(size_t rows, size_t columns,
                                            RngFillMode mode = RngFillMode::COLUMNS);

  template <typename T>
  void spGaussianRandomNumber(std::vector<T> *xv, float fp_scale = 1.0);
  /// \}

  /// \brief Jump forward 2^64 iterations in the sequence.
  void jump();
  
  /// \brief Jump forward 2^96 iterations in the sequence.
  void longJump();

  /// \brief Reveal the current state of the generator.
  ullint4 revealState() const;

  /// \brief Reveal the random bit string.
  ullint revealBitString() const;

  /// \brief Set the current state of the generator.
  void setState(const ullint4 state_in);
  
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

/// \brief An series of "Xorshiro128+" generators, with state vectors for all of them and the
///        means for seeding the series based on long jumps from a single state vector.
template <typename T> class RandomNumberMill {
public:

  /// \brief The constructor can work off of a simple random initialization seed integer or a
  ///        specific state vector for the first generator in the series.
  ///
  /// \param generators_in   The count of generators in the series
  /// \param depth_in        Quantity of random numbers from each generator to store in the bank
  /// \param init_kind       Style in which to initialize the random numbers in the bank,
  ///                        i.e. from a uniform or a normal distribution
  /// \param igseed_in       The seed for the first generator in the series
  /// \param niter           The number of iterations to use in initializing each generator
  /// \param bank_limit      The maximum length of the random number cache
  /// \param state_in        The state to apply to generator zero, thus determining the initial
  ///                        states of all other generators
  /// \{
  RandomNumberMill(const ullint2 state_in, size_t generators_in = 1LLU, size_t depth_in = 2LLU,
                   RandomNumberKind init_kind = RandomNumberKind::GAUSSIAN,
                   size_t bank_limit = constants::giga_zu);

  RandomNumberMill(const ullint4 state_in, size_t generators_in = 1LLU, size_t depth_in = 2LLU,
                   RandomNumberKind init_kind = RandomNumberKind::GAUSSIAN,
                   size_t bank_limit = constants::giga_zu);

  RandomNumberMill(size_t generators_in = 1LLU, size_t depth_in = 2LLU,
                   RandomAlgorithm style_in = RandomAlgorithm::XOSHIRO_256PP,
                   RandomNumberKind init_kind = RandomNumberKind::GAUSSIAN,
                   int igseed_in = default_random_seed, int niter = default_xoshiro256pp_scrub,
                   size_t bank_limit = constants::giga_zu);
  /// \}

  /// \brief Get the number of (forward-jumped) generators.
  size_t getGeneratorCount() const;

  /// \brief Get the depth of the bank for each generator.
  size_t getDepth() const;
  
  /// \brief Get the number of generators for which to refesh all banked values at one time.
  size_t getRefreshStride() const;

  /// \brief Get a random number out of the bank.
  ///
  /// \param generator_index  The generator series from which to obtain the value
  /// \param layer_index      Layer from which to obtain the value
  T getBankValue(size_t generator_index, size_t layer_index) const;
  
  /// \brief Populate a portion of this object's bank with random numbers from each of the
  ///        respective generators.
  ///
  /// \param first_gen  Index of the first generator to draw random numbers from
  /// \param last_gen   Index of the generator before which to stop drawing new random numbers
  void uniformRandomNumbers(size_t first_gen, size_t last_gen);

  /// \brief Populate a portion of this object's bank with normally distributed random numbers
  ///        from each of the respective generators.
  ///
  /// \param first_gen  Index of the first generator to draw random numbers from
  /// \param last_gen   Index of the generator before which to stop drawing new random numbers
  void gaussianRandomNumbers(size_t first_gen, size_t last_gen);

private:
  RandomAlgorithm style;     ///< The kind of random number generator that this mill uses
  size_t generators;         ///< The quantity of (pseudo-) random number generators in the series
  size_t depth;              ///< The depth of each generator's random numbers in the memory bank
  size_t refresh_stride;     ///< The number of segments in which the generators series can be
                             ///<   refreshed, calling a subset of the generators to recalculate
                             ///<   the banked random values for that generator / lane.
  Hybrid<ullint2> state_xy;  ///< State vectors for all random number generators in the series,
                             ///<   sufficient for a Xoroshiro128+ generator or any other requiring
                             ///<   a 128-bit state
  Hybrid<ullint2> state_zw;  ///< State vectors for all random number generators in the series,
                             ///<   extended for Xoshiro256++ and any other generator requiring a
                             ///<   256-bit state vector
  Hybrid<T> bank;            ///< Bank of random numbers pre-computed and saved for later use.  The
                             ///<   numbers are stored only in the specifictype format of the
                             ///<   series, float or double, to save space and avoid ambiguities
                             ///<   as to when the numbers should be refreshed if multiple levels
                             ///<   of precision in the random numbers are required.

  /// \brief Check the inputs for validity and safety, to avoid excessive allocations.
  ///
  /// \param bank_limit  The maximum length of the random number cache
  void checkDimensions(size_t bank_limit);

  /// \brief Prime the bank with its first complement of random numbers.
  ///
  /// \param init_kind  Style in which to initialize the random numbers in the bank, i.e. from a
  ///                   uniform or a normal distribution
  void initializeBank(RandomNumberKind init_kind);
};

/// \brief Initialize an array of Xoroshiro128p generators using CPU-based code.
///
/// Overloaded:
///   - Operate on a C-style array
///   - Operate on a Standard Template Library Vector of RNG state vectors
///   - Operate on a Hybrid object (HOST only--see the overloaded form in hpc_random.h for GPU
///     functionality)
///
/// \brief state_vector  The array of RNG state vectors
/// \brief n_generators  Trusted length of the state_vector array, the number of RN generators
/// \brief igseed        Seed to apply in creating the first random number generator, from which
///                      all others are long jumps plus short jumps 
/// \param scrub_cycles  Number of times to run the random number generator and discard results in
///                      order to raise the quality of random numbers coming out of it
/// \{
void initXoroshiro128pArray(ullint2* state_vector, int n_generators, int igseed,
                            int scrub_cycles = default_xoroshiro128p_scrub);
void initXoroshiro128pArray(std::vector<ullint2> *state_vector, int igseed,
                            int scrub_cycles = default_xoroshiro128p_scrub);
void initXoroshiro128pArray(Hybrid<ullint2> *state_vector, int igseed,
                            int scrub_cycles = default_xoroshiro128p_scrub);
/// \}

/// \brief Initialize an array of Xoroshiro128p generators using CPU-based code.
///
/// Overloaded:
///   - Operate on a C-style array
///   - Operate on a Standard Template Library Vector of RNG state vectors
///   - Operate on a Hybrid object (HOST only--see the overloaded form in hpc_random.h for GPU
///     functionality)
///
/// \brief state_xy      
/// \brief n_generators  Trusted length of the state_vector array, the number of RN generators
/// \brief igseed        Seed to apply in creating the first random number generator, from which
///                      all others are long jumps plus short jumps
/// \param scrub_cycles  Number of times to run the random number generator and discard results in
///                      order to raise the quality of random numbers coming out of it
/// \{
void initXoshiro256ppArray(ullint2* state_xy, ullint2* state_zw, int n_generators, int igseed,
                           int scrub_cycles = default_xoshiro256pp_scrub);
void initXoshiro256ppArray(std::vector<ullint2> *state_xy, std::vector<ullint2> *state_zw,
                           int igseed, int scrub_cycles = default_xoshiro256pp_scrub);
void initXoshiro256ppArray(Hybrid<ullint2> *state_xy, Hybrid<ullint2> *state_zw, int igseed,
                           int scrub_cycles = default_xoshiro256pp_scrub);
/// \}

/// \brief Fill a cache of random numbers using a series (an array) of state vectors (generators).
///        The second of two ullint2 state vectors may be supplied as nullptr if only 128-bit
///        states are in use.
///
/// Overloaded:
///   - Work with 128-bit or 256-bit generator states.
///   - Operate on C-style arrays, Standard Template Library Vectors, or Hybrid objects
///   - An alternative form in the hpc_random library makes use of the eponymous GPU kernels.
///
/// \param state_xy     First halves of each 256-bit generator state vector, or the array of
///                     128-bit state vectors (state_zw should be the null pointer in this case)
/// \param state_zw     Second havles of each 256-bit generator state vector, if required
/// \param cache        Array of real-valued random number results to fill out (the template
///                     parameter indicates the data type of this array)
/// \param length       Trusted length of state_xy and state_zw
/// \param depth        Quantity of random numbers for each generator to produce during checkout
///                     from the state vector arrays.  The (warp size-padded) length parameter
///                     times depth gives the size of the cache array.
/// \param method       Method for generating random numbers (some XOR-shift technique)
/// \param product      Shape of the random distribution over which to take results
/// \param index_start  Starting index at which to begin drawing upon random number generators
/// \param index_end    Upper bound of random number generators from which to draw results
/// \{
template <typename T>
void fillRandomCache(ullint2* state_xy, ullint2* state_zw, T* cache, size_t length, size_t depth,
                     RandomAlgorithm method, RandomNumberKind product, size_t index_start,
                     size_t index_end);

template <typename T>
void fillRandomCache(std::vector<ullint2> *state_xy, std::vector<ullint2> *state_zw,
                     std::vector<T> *cache, size_t length, size_t depth, RandomAlgorithm method,
                     RandomNumberKind product, size_t index_start, size_t index_end);

template <typename T>
void fillRandomCache(Hybrid<ullint2> *state_xy, Hybrid<ullint2> *state_zw, std::vector<T> *cache,
                     size_t length, size_t depth, RandomAlgorithm method, RandomNumberKind product,
                     size_t index_start, size_t index_end);
/// \}

} // namespace random
} // namespace stormm

#include "random.tpp"

#endif
