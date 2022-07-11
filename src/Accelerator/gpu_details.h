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
#include "Constants/behavior.h"
#include "Constants/scaling.h"
#include "Constants/fixed_precision.h"
#include "Potential/energy_enumerators.h"
#include "Synthesis/synthesis_enumerators.h"
#include "DataTypes/omni_vector_types.h"

namespace omni {
namespace card {

using constants::PrecisionModel;
#ifndef OMNI_USE_HPC
using data_types::int2;
#endif
using energy::EvaluateForce;
using energy::EvaluateEnergy;
using numerics::ForceAccumulationMethod;
using numerics::PrecisionLevel;
using synthesis::VwuGoal;

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

/// \brief A class to guide the implementation of GPU kernels, with selected thread counts per
///        block and block counts per launch grid for a specific GPU based on the workload.  This
///        object is constructed in a stepwise manner, with each kind of work unit contributing
///        new launch specifications.
class KernelManager {
public:

  /// \brief The constructor will fill in values as if this were a single-threaded CPU "launch."
  ///        Filling the attributes of this object here, while it might seem convenient, would
  ///        require including each of the headers for the various HPC units.  Including this
  ///        object in their headers would, in turn, lead to circular dependencies.  Member
  ///        variables of this object will therefore be filled by the selectLaunchParameters()
  ///        function wrapped function in a separate library that includes all of the appropriate
  ///        HPC units.
  KernelManager();

  /// \brief Get the block and thread counts for the valence kernel.
  ///
  /// \param prec        The type of floating point numbers in which the kernel shall work
  /// \param eval_force  Indication of whether the kernel will evaluate forces on atoms
  /// \param eval_nrg    Indication of whether to evaluate the energy of the system as a whole
  /// \param acc_meth    The force accumulation method (SPLIT or WHOLE, AUTOMATIC will produce an
  ///                    error in this context)
  int2 getValenceKernelDims(PrecisionModel prec, EvaluateForce eval_force, EvaluateEnergy eval_nrg,
                            ForceAccumulationMethod acc_meth) const;

  /// \brief Get the block and thread counts for the non-bonded kernel.  Parameters descriptions
  ///        for this function follow from getValenceKernelDims() above.
  int2 getNonbondedKernelDims(PrecisionModel prec, EvaluateForce eval_force,
                              EvaluateEnergy eval_nrg, ForceAccumulationMethod acc_meth) const;

  /// \brief Get the block and thread counts for a reduction kernel.
  int2 getReductionKernelDims() const;

  /// \brief Set the register, maximum block size, and threads counts for one of the valence
  ///        kernels.  This function complements getValenceKernelDims(), although the other
  ///        function reports numbers based on this functions input information and some further
  ///        analysis.
  ///
  /// \param thread_limit    Maximum number of threads any one block of the kernel can launch with
  /// \param register_count  The number of registers per thread in this kernel
  /// \param shared_usage    The maximum amount of __shared__ memory that the kernel might
  ///                        allocate, per block, in bytes
  void setValenceKernelAttributes(int thread_limit, int register_count, int shared_usage);

  /// \brief Set the register, maximum block size, and threads counts for one of the non-bonded
  ///        kernels.  Parameter descriptions for this function follow from
  ///        setValenceKernelAttributes() above.
  void setNonbondedKernelAttributes(int thread_limit, int register_count, int shared_usage);

  /// \brief Set the register, maximum block size, and threads counts for one of the reduction
  ///        kernels.  Parameter descriptions for this function follow from
  ///        setValenceKernelAttributes() above.
  void setReductionKernelAttributes(int thread_limit, int register_count, int shared_usage);
  
private:

  // In each of the following int4 tuples, the selected number of threads per block is given in
  // the x member, the number of blocks in the grid is given by the y member, the maximum number
  // of threads per block and the number of registers per thread are packed into the z member's
  // low and high 16 bits, respectively, and the maximum amount of __shared__ memory (statically
  // allocated plus dynamically allocatable) is stored in the w member.  Each piece of information
  // is accessible using the appropriate private member functions.
  int4 valence_kernel_de_dims;    ///< Optimal dimensions for the valence kernel launch grid with
                                  ///<   double-precision, energy-only mode
  int4 valence_kernel_dfs_dims;   ///< Optimal dimensions for the valence kernel launch grid with
                                  ///<   double-precision, force-only, split accumulation
  int4 valence_kernel_dfes_dims;  ///< Optimal dimensions for the valence kernel launch grid with
                                  ///<   double-precision, force and energy computations, split
                                  ///<   accumulation
  int4 valence_kernel_fe_dims;    ///< Optimal dimensions for the valence kernel launch grid with
                                  ///<   single-precision, energy-only mode
  int4 valence_kernel_ffs_dims;   ///< Optimal dimensions for the valence kernel launch grid with
                                  ///<   single-precision, force-only, split accumulation
  int4 valence_kernel_ffes_dims;  ///< Optimal dimensions for the valence kernel launch grid with
                                  ///<   single-precision, force and energy computations, split
                                  ///<   accumulation
  int4 valence_kernel_ffw_dims;   ///< Optimal dimensions for the valence kernel launch grid with
                                  ///<   single-precision, force-only, whole accumulation
  int4 valence_kernel_ffew_dims;  ///< Optimal dimensions for the valence kernel launch grid with
                                  ///<   single-precision, force and energy computations, whole
                                  ///<   accumulation
  int4 nonbond_kernel_de_dims;    ///< Optimal dimensions for the non-bonded kernel launch grid
                                  ///<   with double-precision, energy computations only
  int4 nonbond_kernel_dfs_dims;   ///< Optimal dimensions for the non-bonded kernel launch grid
                                  ///<   with double-precision, force computations only
  int4 nonbond_kernel_dfes_dims;  ///< Optimal dimensions for the non-bonded kernel launch grid
                                  ///<   with double-precision, force and energy computations
  int4 nonbond_kernel_fe_dims;    ///< Optimal dimensions for the non-bonded kernel launch grid
                                  ///<   with single-precision, energy computations only
  int4 nonbond_kernel_ffs_dims;   ///< Optimal dimensions for the non-bonded kernel launch grid
                                  ///<   with single-precision, force computations only, split
                                  ///<   force accumulation
  int4 nonbond_kernel_ffes_dims;  ///< Optimal dimensions for the non-bonded kernel launch grid
                                  ///<   with single-precision, force and energy computations,
                                  ///<   split force accumulation
  int4 nonbond_kernel_ffw_dims;   ///< Optimal dimensions for the non-bonded kernel launch grid
                                  ///<   with single-precision, force computations only, whole
                                  ///<   64-bit int accumulation
  int4 nonbond_kernel_ffew_dims;  ///< Optimal dimensions for the non-bonded kernel launch grid
                                  ///<   with single-precision, force and energy computations,
                                  ///<   whole 64-bit integer force accumulation
  int4 reduction_kernel_dims;     ///< Optimal dimensions for the reduction kernel launch grid

  /// \brief Get the selected number of threads per block for launching a given kernel based on
  ///        the available GPU.
  ///
  /// \param kval  Tuple encoding launch parameters for the kernel of interest
  int getSelectedBlockDim(const int4 kval) const;

  /// \brief Get the selected number of blocks in the launch grid for a given kernel based on the
  ///        available GPU.
  ///
  /// \param kval  Tuple encoding launch parameters for the kernel of interest
  int getSelectedGridDim(const int4 kval) const;

  /// \brief Get the maximum number of threads per block that the compilation of a given kernel
  ///        can tolerate.
  ///
  /// \param kval  Tuple encoding launch parameters for the kernel of interest
  int getMaximumThreadsPerBlock(const int4 kval) const;

  /// \brief Get the number of registers required by each thread of a kernel, as it was compiled.
  ///
  /// \param kval  Tuple encoding launch parameters for the kernel of interest
  int getRegistersPerThread(const int4 kval) const;

  /// \brief Get the largest possible amount of __shared__ memory needed by a function.  This may
  ///        overestimate the true usage, but should be a consideration when selecting the number
  ///        of blocks to put on each streaming multiprocessor.
  int getMaximumSharedMemoryPerBlock(const int4 kval) const;

  /// \brief Choose the grid dimension and thread count for a kernel, based on its known register
  ///        and __shared__ memory usage, maximum threads per block, and the specs of the GPU.
  ///
  /// \param registers_per_thread     Register requirements for each thread in the kernel (this
  ///                                 number will be updated according to a table hard-coded inside
  ///                                 the function to pad the register requirements as necessary
  ///                                 for a realistic kernel thread usage)
  /// \param kernel_max_threads       Maximumthreads per block for which the kernel is compiled
  /// \param shared_memory_per_block  The maximum amount of __shared__ memory that each block might
  ///                                 need
  /// \param gpu                      Specific attributes of the GPU chosen at runtime
  int4 setLaunchDims(int registers_per_thread, int kernel_max_threads,
                     int shared_memory_per_block, const GpuDetails &gpu);
};
  
} // namespace card
} // namespace omni

#ifndef OMNI_USE_HPC
/// \brief ***Global*** GPU descriptor that describes no valid GPU.  This is the equivalent of
///        nullptr for the GpuDetails object, and if passed to various functions that might launch
///        a CUDA kernel will trigger the corresponding CPU process instead.  An equivalent
///        expression of this occurs in hpc_config.cuh if OMNI_USE_HPC is defined.
extern omni::card::GpuDetails null_gpu;
#endif

#endif
