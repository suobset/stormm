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
  int4 valence_kernel_de_dims;
  int4 valence_kernel_dfsm_dims;
  int4 valence_kernel_dfsa_dims;
  int4 valence_kernel_dfesm_dims;
  int4 valence_kernel_dfesa_dims;
  int4 valence_kernel_fe_dims;
  int4 valence_kernel_ffsm_dims;
  int4 valence_kernel_ffwm_dims;
  int4 valence_kernel_ffsa_dims;
  int4 valence_kernel_ffwa_dims;
  int4 valence_kernel_ffewm_dims;
  int4 valence_kernel_ffesm_dims;
  int4 valence_kernel_ffewa_dims;
  int4 valence_kernel_ffesa_dims;
  int4 nonbond_kernel_de_dims;
  int4 nonbond_kernel_dfs_dims;
  int4 nonbond_kernel_dfes_dims;
  int4 nonbond_kernel_fe_dims;
  int4 nonbond_kernel_ffs_dims;
  int4 nonbond_kernel_ffw_dims;
  int4 nonbond_kernel_ffes_dims;
  int4 nonbond_kernel_ffew_dims;
  int4 reduction_kernel_gt_dims;
  int4 reduction_kernel_sc_dims;
  int4 reduction_kernel_ar_dims;
  /// \}
  
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
  /// \param kernel_layout   Information on kernel requirements, including register usage (this
  ///                        number will be updated according to a table hard-coded inside the
  ///                        function to pad the register requirements as necessary), maximum
  ///                        threads per block for which the kernel is compiled, and the maximum
  ///                        amount of __shared__ memory that each block might need
  /// \param gpu             Specific attributes of the GPU chosen at runtime
  int2 setLaunchDims(int4 kernel_layout, const GpuDetails &gpu);
};

} // namespace card
} // namespace omni

#endif
