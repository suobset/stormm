// -*-c++-*-
#ifndef OMNI_HPC_BOUNDS_H
#define OMNI_HPC_BOUNDS_H

namespace omni {
namespace constants {

/// \brief HPC block sizes, based on a maximum number of 1024 threads per block (common to both
///        CUDA and HIP)
/// \{
constexpr int large_block_size  = 1024;
constexpr int medium_block_size = 512;
constexpr int small_block_size  = 256;
constexpr int tiny_block_size   = 128;
/// \}

}
}

// Share the major block bounds
namespace omni {
  using constants::large_block_size;
  using constants::medium_block_size;
  using constants::small_block_size;
  using constants::tiny_block_size;
}

#endif
