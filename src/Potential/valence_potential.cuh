// -*-c++-*-
#ifndef OMNI_VALENCE_POTENTIAL_CUH
#define OMNI_VALENCE_POTENTIAL_CUH

#define EXCL_GMEM_OFSET (blockIdx.x * maximum_valence_work_unit_atoms)

namespace omni {
namespace energy {
    
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
