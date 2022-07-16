// -*-c++-*-
#ifndef OMNI_KERNEL_MANAGER_H
#define OMNI_KERNEL_MANAGER_H

#include <map>
#include <string>
#include "Constants/behavior.h"
#include "Constants/fixed_precision.h"
#include "DataTypes/omni_vector_types.h"
#include "Potential/energy_enumerators.h"
#include "Synthesis/synthesis_enumerators.h"
#include "gpu_details.h"

namespace omni {
namespace card {

using constants::PrecisionModel;
#ifndef OMNI_USE_HPC
using data_types::int2;
#endif
using energy::EvaluateForce;
using energy::EvaluateEnergy;
using numerics::ForceAccumulationMethod;
using synthesis::NbwuKind;
using synthesis::ReductionStage;
using synthesis::VwuGoal;
  
/// \brief Encapsulate the operations to store and retrieve information about a kernel's format.
class KernelFormat {
public:

  /// \brief The constructor lays out a blank object.  Kernel formats are found by querying their
  ///        specifications with HPC library functions.
  KernelFormat();

  /// \brief Get the optimal block and grid sizes for kernel launches with the present GPU.
  int2 getLaunchParameters() const;

  /// \brief Get the register usage of the kernel.
  int getRegisterUsage() const;

  /// \brief Get the maximum thread count for a single block in the kernel launch.
  int getBlockSizeLimit() const;

  /// \brief Get the amount of __shared__ memory needed by any one block.
  int getSharedMemoryRequirement() const;

  /// \brief Get the name of this kernel
  const std::string& getKernelName() const;
  
  /// \brief Set the register usage, block size limit, and shared memory usage for a kernel.
  ///
  /// Overloaded:
  ///   - Provide explicit instructions on whether to consider breaking up the blocks into smaller
  ///     units
  ///   - Assume that the largest possible block size is always to be used
  ///
  /// \param lb_max_threads_per_block  Maximum threads per block, as stated in the launch bounds
  /// \param lb_min_blocks_per_smp     Minimum blocks per multiprocessor, from the launch bounds
  /// \param register_usage_in         Input register usage
  /// \param shared_usage_in           Input __shared__ memory usage
  /// \param block_subdivision         Preferred block multiplicity (this will compound the input
  ///                                  minimum number of blocks per multiprocessor)
  /// \param gpu                  Details of the available GPU (likely passed in from a
  ///                             KernelManager struct containing many KernelFormat objects)
  /// \param kernel_name_in       The name of the kernel, for reporting purposes later (optional)
  /// \{
  void build(int lb_max_threads_per_block, int lb_min_blocks_per_smp, int register_usage_in,
             int shared_usage_in, int block_subdivision, const GpuDetails &gpu,
             const std::string &kernel_name_in = std::string(""));

  void build(int lb_max_threads_per_block, int lb_min_blocks_per_smp, int register_usage_in,
             int shared_usage_in, const GpuDetails &gpu,
             const std::string &kernel_name_in = std::string(""));
  /// \}
  
private:
  int block_dimension;         ///< Computed optimal block dimension to use in kernel launches
  int grid_dimension;          ///< Computed optimal grid size to use in kernel launches
  int register_usage;          ///< The number of registers needed by each thread of the kernel
                               ///<   as it is compiled for the current executable
  int block_size_limit;        ///< The largest block size usable by the kernel launch (exceeding
                               ///<   this will cause the kernel launch to fail)
  int shared_usage;            ///< The maximum amount of __shared__ memory needed by each block
  std::string kernel_name;     ///< Name of the kernel, for reporting purposes
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
  ///
  /// \param gpu_in                  Details of the GPU in use (this is relevant, as it will be
  ///                                used to interpret the layout of any kernels)
  KernelManager(const GpuDetails &gpu);

  /// \brief Get the block and thread counts for the valence kernel.
  ///
  /// \param prec        The type of floating point numbers in which the kernel shall work
  /// \param eval_force  Indication of whether the kernel will evaluate forces on atoms
  /// \param eval_nrg    Indication of whether to evaluate the energy of the system as a whole
  /// \param acc_meth    The force accumulation method (SPLIT or WHOLE, AUTOMATIC will produce an
  ///                    error in this context)
  int2 getValenceKernelDims(PrecisionModel prec, EvaluateForce eval_force, EvaluateEnergy eval_nrg,
                            ForceAccumulationMethod acc_meth, VwuGoal purpose) const;

  /// \brief Get the block and thread counts for the non-bonded kernel.  Parameters descriptions
  ///        for this function follow from getValenceKernelDims() above, with the addition of:
  ///
  /// \param kind  The type of non-bonded work unit: tile groups, supertiles, or honeycomb
  ///              being relevant
  int2 getNonbondedKernelDims(PrecisionModel prec, NbwuKind kind, EvaluateForce eval_force,
                              EvaluateEnergy eval_nrg, ForceAccumulationMethod acc_meth) const;

  /// \brief Get the block and thread counts for a reduction kernel.
  ///
  /// \param process  The reduction step to perform
  int2 getReductionKernelDims(ReductionStage process) const;

  /// \brief Get the GPU information for the active GPU.
  const GpuDetails& getGpu() const;
  
  /// \brief Set the register, maximum block size, and threads counts for one of the valence
  ///        kernels.  This function complements getValenceKernelDims(), although the other
  ///        function reports numbers based on this functions input information and some further
  ///        analysis.
  ///
  /// \param prec            The type of floating point numbers in which the kernel shall work
  /// \param eval_force      Indication of whether the kernel will evaluate forces on atoms
  /// \param eval_nrg        Indication of whether to evaluate the energy of the system as a whole
  /// \param acc_meth        The force accumulation method (SPLIT or WHOLE, AUTOMATIC will produce
  ///                        an error in this context)
  /// \param thread_limit    Maximum number of threads any one block of the kernel can launch with
  /// \param register_count  The number of registers per thread in this kernel
  /// \param shared_usage    The maximum amount of __shared__ memory that the kernel might
  ///                        allocate, per block, in bytes
  /// \param block_mult      The requested multiplicity of blocks on each streaming multiprocessor
  void catalogValenceKernel(PrecisionModel prec, EvaluateForce eval_force, EvaluateEnergy eval_nrg,
                            ForceAccumulationMethod acc_meth, VwuGoal purpose,
                            int lb_max_threads_per_block, int lb_max_blocks_per_smp,
                            int register_count, int shared_usage, int block_subdivision = 1);

  /// \brief Set the register, maximum block size, and threads counts for one of the non-bonded
  ///        kernels.  Parameter descriptions for this function follow from
  ///        setValenceKernelAttributes() above, with the addition of:
  ///
  /// \param kind  The type of non-bonded work unit: tile groups, supertiles, or honeycomb
  ///              being relevant
  void catalogNonbondedKernel(PrecisionModel prec, NbwuKind kind, EvaluateForce eval_force,
                              EvaluateEnergy eval_nrg, ForceAccumulationMethod acc_meth,
                              int lb_max_threads_per_block, int lb_max_blocks_per_smp,
                              int register_count, int shared_usage, int block_subdivision = 1);

  /// \brief Set the register, maximum block size, and threads counts for one of the reduction
  ///        kernels.  Parameter descriptions for this function follow from
  ///        setValenceKernelAttributes() above, with the addition of:
  ///
  /// \param process  How far to take the reduction operation
  void setReductionKernelAttributes(ReductionStage process, int thread_limit, int register_count,
                                    int shared_usage, int block_mult = 0);

  /// \brief Print out the kernel launch parameters found for this workload.
  ///
  /// \param kernel_name  Name of the kernel for which to print the parameters (if blank, no
  ///                     kernels' parameters will be printed, and if "ALL" (case-insensitive),
  ///                     all kernels' parameters will be printed)
  void printLaunchParameters(const std::string &kernel_name = std::string("")) const;
  
private:

  /// The details of the GPU in use are simply copied into this object.
  GpuDetails gpu;
  
  /// Store the resource requirements and selected launch parameters for a variety of kernels.
  /// is accessible using the appropriate private member functions.  A descriptive name for the
  /// varaiable matches the kernel, followed by a letter code with the following meanings:
  /// \{
  KernelFormat valence_kernel_de_dims;
  KernelFormat valence_kernel_dfsm_dims;
  KernelFormat valence_kernel_dfsa_dims;
  KernelFormat valence_kernel_dfesm_dims;
  KernelFormat valence_kernel_dfesa_dims;
  KernelFormat valence_kernel_fe_dims;
  KernelFormat valence_kernel_ffsm_dims;
  KernelFormat valence_kernel_ffwm_dims;
  KernelFormat valence_kernel_ffsa_dims;
  KernelFormat valence_kernel_ffwa_dims;
  KernelFormat valence_kernel_ffewm_dims;
  KernelFormat valence_kernel_ffesm_dims;
  KernelFormat valence_kernel_ffewa_dims;
  KernelFormat valence_kernel_ffesa_dims;
  KernelFormat nonbond_kernel_de_dims;
  KernelFormat nonbond_kernel_dfs_dims;
  KernelFormat nonbond_kernel_dfes_dims;
  KernelFormat nonbond_kernel_fe_dims;
  KernelFormat nonbond_kernel_ffs_dims;
  KernelFormat nonbond_kernel_ffw_dims;
  KernelFormat nonbond_kernel_ffes_dims;
  KernelFormat nonbond_kernel_ffew_dims;
  KernelFormat reduction_kernel_gt_dims;
  KernelFormat reduction_kernel_sc_dims;
  KernelFormat reduction_kernel_ar_dims;
  /// \}

  std::map<std::string, KernelFormat> k_dictionary;
};

/// \brief Obtain a unique string identifier for one of the valence kernels.  Each identifier
///        begins with "vale_" and is then appended with letter codes for different aspects
///        according to the following system:
///        - { d, f }      Perform calculations in double (d) or float (f) arithmetic
///        - { e, f, fe }  Compute energies (e), forces (f), or both (ef)
///        - { s, w }      Accumulate forces in split integers (s) or whole integers (w)
///        - { m, a }      Move atoms (m) or accumulate forces or energies (a)
///
/// \param prec        The type of floating point numbers in which the kernel shall work
/// \param eval_force  Indication of whether the kernel will evaluate forces on atoms
/// \param eval_nrg    Indication of whether to evaluate the energy of the system as a whole
/// \param acc_meth    The force accumulation method (SPLIT or WHOLE, AUTOMATIC will produce an
///                    error in this context)
/// \param purpose     The intended action to take with computed forces
std::string valenceKernelKey(PrecisionModel prec, EvaluateForce eval_force,
                             EvaluateEnergy eval_nrg, ForceAccumulationMethod acc_meth,
                             VwuGoal purpose);

/// \brief Obtain a unique string identifier for one of the non-bonded kernels.  Each identifier
///        begins with "nonb_" and is then appended with letter codes for different aspects
///        according to the following system:
///        - { d, f }        Perform calculations in double (d) or float (f) arithmetic
///        - { tg, st, hc }  Use a "tile groups" or "supertiles" strategy for breaking down
///                          systems with isolated boundary conditions, or a "honeycomb" strategy
///                          for breaking down systems with periodic boundary conditions.
///        - { e, f, fe }    Compute energies (e), forces (f), or both (ef)
///        - { s, w }        Accumulate forces in split integers (s) or whole integers (w)
///
/// \param prec        The type of floating point numbers in which the kernel shall work
/// \param kind        The type of non-bonded work unit to evaluate
/// \param eval_force  Indication of whether the kernel will evaluate forces on atoms
/// \param eval_nrg    Indication of whether to evaluate the energy of the system as a whole
/// \param acc_meth    The force accumulation method (SPLIT or WHOLE, AUTOMATIC will produce an
///                    error in this context)
std::string nonbondedKernelKey(PrecisionModel prec, NbwuKind kind, EvaluateForce eval_force,
                               EvaluateEnergy eval_nrg, ForceAccumulationMethod acc_meth);

} // namespace card
} // namespace omni

#endif
