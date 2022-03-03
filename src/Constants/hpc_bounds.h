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

/// \brief HPC block sizes, based on a maximum number of 1024 threads per block (common to both
///        CUDA and HIP)
/// \{
constexpr int tiny_block_size   = 128;
constexpr int small_block_size  = 256;
constexpr int medium_block_size = 512;
constexpr int large_block_size  = 1024;
/// \}

}
}

// Share the major block bounds
namespace omni {
  using constants::BlockSize;
  using constants::large_block_size;
  using constants::medium_block_size;
  using constants::small_block_size;
  using constants::tiny_block_size;
}

#endif
