// -*-c++-*-
#ifndef OMNI_VALENCE_POTENTIAL_CUH
#define OMNI_VALENCE_POTENTIAL_CUH

#define EXCL_GMEM_OFSET (blockIdx.x * maximum_valence_work_unit_atoms)

namespace omni {
namespace energy {

/// \brief Convert a floating point number into a fixed-precision representation with two integers,
///        using atomic additions to both accumulators.
///
/// \param fval                The value to convert (pre-inflated by whatever scaling factor for
///                            the number of fixed-precision bits after the decimal)
/// \param pos                 Position of the atom in the various accumulator arrays
/// \param sh_primary          Array of primary accumulators (low 32 bits), most likely resident
///                            in __shared__ L1 memory
/// \param sh_overflow_active  Array of flags to indicate that at least one of each consecutive set
///                            of 32 atoms' forces has received an overflow accumulation.  This is
///                            stored in __shared__ memory as its size is relatively small.
/// \param overflow            Array of overflow accumulators (high 31 bits) resident in global
///                            memory
__device__ __forceinline__ void splitForceContribution(const float fval, const int pos,
                                                       int* sh_primary, int* sh_overflow_active,
                                                       int* overflow);
  
/// \brief Evaluate valence work units and move atoms.
///
/// Overloaded:
///   - Perform work in single or double precision
///   - Compute forces, energy, or both
///   - Move particles
///
/// \param poly_vk  Valence parameters based on consensus tables from a topology synthesis
/// \param poly_ps  Coordinates, velocities, and forces of all systems
/// \param gmem_x   Exclusive space in global memory arrays reserved for each thread block, to be
///                 brought into free L1 cache
/// \{
__global__ void __launch_bounds__(large_block_size, 1)
kfValenceForceAccumulation(const SyValenceKit<float> poly_vk, PsSynthesisWriter poly_ps,
                           CacheResource<float> gmem_r);

__global__ void __launch_bounds__(large_block_size, 1)
kdValenceForceAccumulation(const SyValenceKit<double> poly_vk, PsSynthesisWriter poly_ps,
                           CacheResource<double> gmem_r);
/// \}
  
} // namespace energy
} // namespace omni

#endif
