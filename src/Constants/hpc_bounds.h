// -*-c++-*-
#ifndef OMNI_HPC_BOUNDS_H
#define OMNI_HPC_BOUNDS_H

namespace omni {
namespace constants {

/// \brief Enumerator for general block sizes
enum class BlockSize {
  TINY,    ///< Tiny block of 128 threads
  SMALL,   ///< Small block of 256 threads
  MEDIUM,  ///< Medium-sized block of 512 threads
  LARGE    ///< Largest block of 1024 threads
};

/// \brief Increment for block size tests used in determining the optimal size for kernel launches
constexpr int block_size_test_increment = 64;
  
/// \brief HPC block sizes, based on a maximum number of 1024 threads per block (common to both
///        CUDA and HIP).  The tiny block size represents a minimum size for any one thread block,
///        and the large block size represents a maximum.  Blocks requiring thread counts other
///        than those presented here are possible, within those limits.
/// \{
constexpr int tiny_block_size   = 128;
constexpr int small_block_size  = 256;
constexpr int medium_block_size = 512;
constexpr int large_block_size  = 1024;
/// \}

} // namespace constants
} // namespace omni

// Share the major block bounds
namespace omni {
  using constants::BlockSize;
  using constants::large_block_size;
  using constants::medium_block_size;
  using constants::small_block_size;
  using constants::tiny_block_size;
}

#endif
