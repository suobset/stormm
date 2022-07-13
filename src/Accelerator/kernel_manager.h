// -*-c++-*-
#ifndef OMNI_KERNEL_MANAGER_H
#define OMNI_KERNEL_MANAGER_H

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

  /// \brief Set the register usage, block size limit, and shared memory usage for a kernel.
  ///
  /// Overloaded:
  ///   - Provide explicit instructions on whether to consider breaking up the blocks into smaller
  ///     units
  ///   - Assume that the largest possible block size is always to be used
  ///
  /// \param register_usage_in    Input register usage
  /// \param block_size_limit_in  Input block size limit
  /// \param shared_usage_in      Input __shared__ memory usage
  /// \param block_dimension_in   Preferred block size (this will replace the maximum size, if
  ///                             provided)
  /// \param gpu                  Details of the available GPU (likely passed in from a
  ///                             KernelManager struct containing many KernelFormat objects)
  /// \{
  void build(int register_usage_in, int block_size_limit_in, int shared_usage_in,
             int block_dimension_in, const GpuDetails &gpu);

  void build(int register_usage_in, int block_size_limit_in, int shared_usage_in,
             const GpuDetails &gpu);
  /// \}
  
private:
  int block_dimension;         ///< Computed optimal block dimension to use in kernel launches
  int grid_dimension;          ///< Computed optimal grid size to use in kernel launches
  int register_usage;          ///< The number of registers needed by each thread of the kernel
                               ///<   as it is compiled for the current executable
  int block_size_limit;        ///< The largest block size usable by the kernel launch (exceeding
                               ///<   this will cause the kernel launch to fail)
  int shared_usage;            ///< The maximum amount of __shared__ memory needed by each block
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
  /// \param gpu_in  Details of the GPU in use (this is relevant, as it will be used to interpret
  ///                the layout of any kernels)
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
  ///        for this function follow from getValenceKernelDims() above.
  int2 getNonbondedKernelDims(PrecisionModel prec, EvaluateForce eval_force,
                              EvaluateEnergy eval_nrg, ForceAccumulationMethod acc_meth) const;

  /// \brief Get the block and thread counts for a reduction kernel.
  ///
  /// \param process  The reduction step to perform
  int2 getReductionKernelDims(ReductionStage process) const;

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
  void setValenceKernelAttributes(PrecisionModel prec, EvaluateForce eval_force,
                                  EvaluateEnergy eval_nrg, ForceAccumulationMethod acc_meth,
                                  VwuGoal purpose, int thread_limit, int register_count,
                                  int shared_usage);

  /// \brief Set the register, maximum block size, and threads counts for one of the non-bonded
  ///        kernels.  Parameter descriptions for this function follow from
  ///        setValenceKernelAttributes() above, with the addition of:
  ///
  /// \param kind  The type of non-bonded work unit: tile groups, supertiles, or honeycomb
  ///              being relevant
  void setNonbondedKernelAttributes(PrecisionModel prec, EvaluateForce eval_force,
                                    EvaluateEnergy eval_nrg, ForceAccumulationMethod acc_meth,
                                    NbwuKind kind, int thread_limit, int register_count,
                                    int shared_usage);

  /// \brief Set the register, maximum block size, and threads counts for one of the reduction
  ///        kernels.  Parameter descriptions for this function follow from
  ///        setValenceKernelAttributes() above, with the addition of:
  ///
  /// \param process  How far to take the reduction operation
  void setReductionKernelAttributes(ReductionStage process, int thread_limit, int register_count,
                                    int shared_usage);
  
private:

  /// The details of the GPU in use are simply copied into this object.
  GpuDetails gpu;
  
  /// In each of the following int4 tuples, the selected number of threads per block is given in
  /// the x member, the number of blocks in the grid is given by the y member, the maximum number
  /// of threads per block and the number of registers per thread are packed into the z member's
  /// low and high 16 bits, respectively, and the maximum amount of __shared__ memory (statically
  /// allocated plus dynamically allocatable) is stored in the w member.  Each piece of information
  /// is accessible using the appropriate private member functions.  A descriptive name for the
  /// varaiable matches the kernel, followed by a letter code with the following meanings:
  ///   - { d, f }      The kernel performs calculations in double (d) or float (f) arithmetic
  ///   - { e, f, fe }  The kernel computes energies (e), forces (f), or both (ef)
  ///   - { s, w }      The kernel accumulates forces in split integers (s) or whole integers (w)
  ///   - { m, a }      The kernel moves atoms (m) or accumulates forces or energies (a)
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
};

} // namespace card
} // namespace omni

#endif
