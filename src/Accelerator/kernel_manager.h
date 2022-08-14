// -*-c++-*-
#ifndef STORMM_KERNEL_MANAGER_H
#define STORMM_KERNEL_MANAGER_H

#include <map>
#include <string>
#ifdef STORMM_USE_HPC
#  ifdef STORMM_USE_CUDA
#    include <cuda_runtime.h>
#  endif
#endif
#include "copyright.h"
#include "Constants/behavior.h"
#include "Constants/fixed_precision.h"
#include "DataTypes/stormm_vector_types.h"
#include "Math/reduction.h"
#include "Potential/energy_enumerators.h"
#include "Structure/structure_enumerators.h"
#include "Synthesis/atomgraph_synthesis.h"
#include "Synthesis/synthesis_enumerators.h"
#include "Topology/atomgraph_enumerators.h"
#include "gpu_details.h"

namespace stormm {
namespace card {

using constants::PrecisionModel;
#ifndef STORMM_USE_HPC
using data_types::int2;
#endif
using energy::EvaluateForce;
using energy::EvaluateEnergy;
using numerics::ForceAccumulationMethod;
using structure::VirtualSiteActivity;
using synthesis::AtomGraphSynthesis;
using synthesis::NbwuKind;
using math::ReductionGoal;
using math::ReductionStage;
using synthesis::VwuGoal;
using topology::UnitCellType;

/// \brief Encapsulate the operations to store and retrieve information about a kernel's format.
class KernelFormat {
public:

  /// \brief The constructor takes launch bounds and other information that can be plucked from a
  ///        cudaFuncAttributes object.
  ///
  /// Overloaded:
  ///   - Construct a blank object
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
  /// \param attr                      Result of a CUDA runtime query to get kernel specifications
  /// \param gpu                       Details of the available GPU (likely passed in from a
  ///                                  KernelManager struct containing many KernelFormat objects)
  /// \param kernel_name_in            Name of the kernel, for reporting purposes later (optional)
  /// \{
  KernelFormat();
  
  KernelFormat(int lb_max_threads_per_block, int lb_min_blocks_per_smp, int register_usage_in,
               int shared_usage_in, int block_subdivision, const GpuDetails &gpu,
               const std::string &kernel_name_in = std::string(""));

  KernelFormat(int lb_max_threads_per_block, int lb_min_blocks_per_smp, int register_usage_in,
               int shared_usage_in, const GpuDetails &gpu,
               const std::string &kernel_name_in = std::string(""));

#ifdef STORMM_USE_HPC
#  ifdef STORMM_USE_CUDA
  KernelFormat(const cudaFuncAttributes &attr, int lb_min_blocks_per_smp, int block_subdivision,
               const GpuDetails &gpu, const std::string &kernel_name_in = std::string(""));
#  endif
#endif
  /// \}

  /// \brief Take the default copy and move constructors as well as assignment operators.
  /// \{
  KernelFormat(const KernelFormat &original) = default;
  KernelFormat(KernelFormat &&original) = default;
  KernelFormat& operator=(const KernelFormat &other) = default;
  KernelFormat& operator=(KernelFormat &&other) = default;
  /// \}
  
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
  
private:
  int block_size_limit;        ///< The largest block size usable by the kernel launch (exceeding
                               ///<   this will cause the kernel launch to fail)
  int shared_usage;            ///< The maximum amount of __shared__ memory needed by each block
  int block_dimension;         ///< Computed optimal block dimension to use in kernel launches
  int grid_dimension;          ///< Computed optimal grid size to use in kernel launches
  int register_usage;          ///< The number of registers needed by each thread of the kernel
                               ///<   as it is compiled for the current executable
  std::string kernel_name;     ///< Name of the kernel, for reporting purposes
};

/// \brief A class to guide the implementation of GPU kernels, with selected thread counts per
///        block and block counts per launch grid for a specific GPU based on the workload.  This
///        object is constructed in a stepwise manner, with each kind of work unit contributing
///        new launch specifications.
class KernelManager {
public:

  /// \brief The constructor will fill in values as if this were a single-threaded CPU "launch."
  ///
  /// \param gpu_in   Details of the GPU in use (this is relevant, as it will be used to interpret
  ///                 the layout of any kernels)
  /// \param poly_ag  Topologies for all systems, offering details of the workload
  KernelManager(const GpuDetails &gpu_in, const AtomGraphSynthesis &poly_ag);

  /// \brief Get the architecture-specific block multiplier.  This will run a minimum number of
  ///        blocks per streaming multiprocessor on some cards, specifically NVIDIA's GTX 1080-Ti,
  ///        when large blocks cannot use more than 32k registers in all.
  int getArchBlockMultiplier() const;
  
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
  /// \param prec     The type of floating point numbers represented by the kernel's substrate
  ///                 (all reductions happenin double-precision real numbers)
  /// \param purpose  The procedure requested, for which the appropriate kernel shall be found 
  /// \param process  The reduction step to perform
  int2 getReductionKernelDims(PrecisionModel prec, ReductionGoal purpose,
                              ReductionStage process) const;

  /// \brief Get the block and thread counts for a virtual site placement or force transmission
  ///        kernel.
  ///
  /// \param prec     Implies a level of detail in the fixed-precision coordinate representation,
  ///                 and specifies the precision in which to perform the placement or force
  ///                 transmission operations
  /// \param purpose  The process to perform with standalone virtual site kernels: place virtual
  ///                 sites or transmit forces from them to frame atoms with mass
  int2 getVirtualSiteKernelDims(PrecisionModel prec, VirtualSiteActivity purpose) const;

  /// \brief Get the GPU information for the active GPU.
  const GpuDetails& getGpu() const;
  
  /// \brief Print out the kernel launch parameters found for this workload.
  ///
  /// \param kernel_name  Name of the kernel for which to print the parameters (if blank, no
  ///                     kernels' parameters will be printed, and if "ALL" (case-insensitive),
  ///                     all kernels' parameters will be printed)
  void printLaunchParameters(const std::string &kernel_name = std::string("")) const;
  
private:

  /// The details of the GPU in use are simply copied into this object.
  GpuDetails gpu;

  /// The workload-specific block multiplier for valence kernels.  No provision is needed for
  /// NVIDIA GTX 1080-Ti, as the blocks will have at least multiplicity 2 on each streaming
  /// multiprocessor, but the workload may lead to a choice of smaller blocks for higher
  /// throughput.  The maximum number of threads per streaming multiprocessor caps at 1024 and
  /// the block multiplicity caps at 4, so no provisions are needed for Turing architectures,
  /// either. 
  int valence_block_multiplier;

  /// The architecture-specific block multiplier for non-bonded kernels, essentially a provision
  /// for NVIDIA Turing cards.  The non-bonded kernels are not subject to the same constraints and
  /// block multiplicity as the valence kernels, as the non-bonded kernels will always run at
  /// least two blocks per streaming multiprocessor.  This value does depend on the unit cell type
  /// of the topology synthesis at hand.
  int nonbond_block_multiplier;

  /// The workload-specific block multiplier for reduction kernels.  Like the valence kernels, the
  /// thread count per streaming multiprocessor will not go above 1024 (this time out of bandwidth
  /// limitations), but the block multiplicity (which starts at 4) could be increased.
  int reduction_block_multiplier;

  /// The architecture-specific block multiplier for virtual site handling kernels.
  int virtual_site_block_multiplier;
  
  /// Store the resource requirements and selected launch parameters for a variety of kernels.
  /// Keys are determined according to the free functions further on in this library.
  std::map<std::string, KernelFormat> k_dictionary;

  /// \brief Set the register, maximum block size, and thread counts for one of the valence
  ///        kernels.  This function complements getValenceKernelDims(), although the other
  ///        function reports numbers based on this functions input information and some further
  ///        analysis.
  ///
  /// \param prec         The type of floating point numbers in which the kernel shall work
  /// \param eval_force   Indication of whether the kernel will evaluate forces on atoms
  /// \param eval_nrg     Indication of whether to evaluate the energy of the system as a whole
  /// \param acc_meth     The force accumulation method (SPLIT or WHOLE, AUTOMATIC will produce
  ///                     an error in this context)
  /// \param subdivision  Number of times that the basic valence kernel should be subdivided
  /// \param kernel_name  [Optional] Name of the kernel in the actual code
  void catalogValenceKernel(PrecisionModel prec, EvaluateForce eval_force, EvaluateEnergy eval_nrg,
                            ForceAccumulationMethod acc_meth, VwuGoal purpose, int subdivision,
                            const std::string &kernel_name = std::string(""));

  /// \brief Set the register, maximum block size, and thread counts for one of the non-bonded
  ///        kernels.  Parameter descriptions for this function follow from
  ///        setValenceKernelAttributes() above, with the addition of:
  ///
  /// \param kind         The type of non-bonded work unit: tile groups, supertiles, or honeycomb
  ///                     being relevant
  void catalogNonbondedKernel(PrecisionModel prec, NbwuKind kind, EvaluateForce eval_force,
                              EvaluateEnergy eval_nrg, ForceAccumulationMethod acc_meth,
                              const std::string &kernel_name = std::string(""));

  /// \brief Set the register, maximum block size, and thread counts for one of the reduction
  ///        kernels.  Parameter descriptions for this function follow from
  ///        setValenceKernelAttributes() above, with the addition of:
  ///
  /// \param purpose      Reason for doing the reduction, i.e. conjugate gradient transformation
  /// \param process      How far to take the reduction operation
  void catalogReductionKernel(PrecisionModel prec, ReductionGoal purpose, ReductionStage process,
                              int subdivision, const std::string &kernel_name = std::string(""));

  /// \brief Set the register, maximum block size, and thread counts for one of the virtual site
  ///        placement kernels.  The first two parameters follow from getVirtualSiteKernelDims(),
  ///        above.
  ///
  /// \param prec         The type of floating point numbers in which the kernel shall work
  /// \param purpose      The process to perform with standalone virtual site kernels
  /// \param subdivision  Number of times that the basic virtual site kernel should be subdivided
  /// \param kernel_name  [Optional] Name of the kernel in the actual code
  void catalogVirtualSiteKernel(PrecisionModel prec, VirtualSiteActivity purpose, int subdivision,
                                const std::string &kernel_name = std::string(""));
};

/// \brief Obtain the workload-specific block multiplier for valence interaction kernels.
int valenceBlockMultiplier();

/// \brief Obtain the architecture-specific block multiplier for non-bonded interaction kernels.
///
/// \param gpu        Details of the GPU that will perform the calculations
/// \param unit_cell  The unit cell type of the systems to evaluate
int nonbondedBlockMultiplier(const GpuDetails &gpu, UnitCellType unit_cell);

/// \brief Obtain the workload-specific block multiplier for reduction kernels.
int reductionBlockMultiplier();

/// \brief Obtain the workload-specific block multiplier for virtual site handling kernels.  The
///        typical kernel will run on 256 threads.
int virtualSiteBlockMultiplier();
  
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

/// \brief Obtain a unique string identifier for one of the reduction kernels.  Each identifier
///        begins with "redc_" and is then appended with letter codes for different aspects
///        according to the following system:
///        - { d, f }        Perform calculations in double (d) or float (f) arithmetic
///        - { cg }          Calculate the conjugate gradient vector
///        - { gt, sc, rd }  Perform a gathering, scattering, or combined all-reduce operation
///
/// \param prec     The type of floating point numbers in which the kernel shall work
/// \param purpose  The reason for doing the reduction--the purpose of each kernel is unique
/// \param process  The reduction stage to perform
std::string reductionKernelKey(PrecisionModel prec, ReductionGoal purpose, ReductionStage process);

/// \brief Obtain a unique string identifier for one of the virtaul site handling kernels.  Each
///        identifier begins with "vste_" and is then appended with letter codes for different
///        activities according to the following system:
///        - { d, f }    Perform calculations in double (d) or float (f) arithmetic
///        - { pl, xm }  Place particles or transmit forces to atoms with mass
///
/// \param prec     The type of floating point numbers in which the kernel shall work
/// \param purpose  The process to perform with standalone virtual site kernels
std::string virtualSiteKernelKey(PrecisionModel prec, VirtualSiteActivity process);
  
} // namespace card
} // namespace stormm

#endif
