// -*-c++-*-
#ifndef OMNI_HPC_STATUS
#define OMNI_HPC_STATUS

#include <vector>
#include <string>
#ifdef OMNI_USE_HPC
#  ifdef OMNI_USE_CUDA
#    include <cuda_runtime.h>
#  endif
#endif
#include "Constants/scaling.h"

namespace omni {
namespace card {

/// \brief Detect the amount of GPU memory occupied by a given process to determine whether it is
///        a significant process.  Occupying a megabyte of GPU RAM is one way for a process to
///        qualify as significant.
/// \{
constexpr long long int significant_gpu_memory = constants::mega;
/// \}
  
/// \brief Pertinent aspects of one particular GPU.  Condensing the data for each GPU in this
///        manner helps to ensure that one cache line will obtain all statistics for a single GPU.
class GpuDetails {
public:

  /// \brief Constructors include a blank constructor (which automatically labels the GPU as
  ///        unavailable) and constructors based on cudaDeviceProp or hipDeviceProp.
  ///
  /// \param devprop    A CUDA device properties object reported by cudaGetDeviceProperties()
  /// \param dev_index  Index of the GPU in a longer list produced by the CUDA runtime library
  /// \{
  GpuDetails();
#ifdef OMNI_USE_HPC
#  ifdef OMNI_USE_CUDA
  GpuDetails(const cudaDeviceProp &devprop, int dev_index);
#  endif
#endif
  /// \}

  /// \brief Get the availability of the GPU
  bool getAvailability() const;

  /// \brief Get whether the architecture of the GPU is supported
  bool getGpuSupported() const;

  /// \brief Get the major architectural number of the GPU
  int getArchMajor() const;

  /// \brief Get the minor architectural number of the GPU
  int getArchMinor() const;

  /// \brief Get the number of streaming multiprocessors on the GPU
  int getSMPCount() const;

  /// \brief Get the amount of RAM on the GPU, in megabytes (assuming that all of the device RAM
  ///        is available for any given process)
  int getCardRam() const;

  /// \brief Get the maximum number of threads per block supported by this GPU
  int getMaxThreadsPerBlock() const;
  
  /// \brief Get the maximum number of threads per streaming multiprocessor on this GPU
  int getMaxThreadsPerSMP() const;

  /// \brief Get the maximum number of blocks per streaming multiprocessor on this GPU
  int getMaxBlocksPerSMP() const;

  /// \brief Get the maximum amount of L1 __shared__ memory per block on this GPU
  int getMaxSharedPerBlock() const;
  
  /// \brief Get the available L1 __shared__ memory per streaming multiprocessor on this GPU
  int getMaxSharedPerSMP() const;

  /// \brief Get the maximum number of registers available per block on this GPU
  int getRegistersPerBlock() const;
  
  /// \brief Get the maximum number of registers available per streaming multiprocessor on this GPU
  int getRegistersPerSMP() const;

  /// \brief Get the name of the GPU
  std::string getCardName() const;

private:
  bool available;            ///< Flag to indicate whether a GPU is available for the program's use
  bool supported;            ///< Flag to indicate whether OMNI supports this GPU
  int arch_major;            ///< Major architecture numbers for each GPU
  int arch_minor;            ///< Minor architecture numbers for each GPU
  int smp_count;             ///< Number of streaming multiprocessors in each GPU
  int card_ram;              ///< The amount of RAM available on each GPU, in megabtyes
  int max_threads_per_block; ///< The maximum number of threads per thread block
  int max_threads_per_smp;   ///< Number of threads one streaming multiprocessor (SMP) can handle
  int max_blocks_per_smp;    ///< Maximum number of blocks permissible on one SMP
  int max_shared_per_block;  ///< Maximum shared memory available per block (bytes)
  int max_shared_per_smp;    ///< Maximum shared memory available per SMP (bytes)
  int registers_per_block;   ///< Number of registers available for each thread block
  int registers_per_smp;     ///< Size of the register file on each SMP
  std::string card_name;     ///< Name of the card according to the server
};

} // namespace card
} // namespace omni

/// \brief ***Global*** GPU descriptor that describes no valid GPU.  This is the equivalent of
///        nullptr for the GpuDetails object, and if passed to various functions that might launch
///        a CUDA kernel will trigger the corresponding CPU process instead.
extern omni::card::GpuDetails null_gpu;

#endif
