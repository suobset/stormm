#include "Math/rounding.h"
#include "cacheresource.h"

namespace omni {
namespace energy {

using card::HybridKind;
using math::roundUp;
  
//-------------------------------------------------------------------------------------------------
CacheResource::CacheResource(const int block_limit_in, const int atom_limit_in) :
    block_limit{block_limit_in},
    atom_limit{roundUp(atom_limit_in, warp_size_int)},
    x_coordinates{HybridKind::POINTER, "cache_xcrd"},
    y_coordinates{HybridKind::POINTER, "cache_ycrd"},
    z_coordinates{HybridKind::POINTER, "cache_zcrd"},
    x_velocities{HybridKind::POINTER, "cache_xvel"},
    y_velocities{HybridKind::POINTER, "cache_yvel"},
    z_velocities{HybridKind::POINTER, "cache_zvel"},
    x_force_overflow{HybridKind::POINTER, "cache_xfov"},
    y_force_overflow{HybridKind::POINTER, "cache_yfov"},
    z_force_overflow{HybridKind::POINTER, "cache_zfov"},
    charges{static_cast<size_t>(block_limit * atom_limit), "cache_charges"},
    sp_charges{static_cast<size_t>(block_limit * atom_limit), "sp_cache_charges"},
    lennard_jones_indices{HybridKind::POINTER, "cache_ljidx"},
    int_data{static_cast<size_t>(4 * block_limit * atom_limit), "cache_int_data"},
    llint_data{static_cast<size_t>(6 * block_limit * atom_limit), "cache_llint_data"}
{
  if (block_limit <= 0 || atom_limit <= 0) {
    rtErr("Device resources private to individual blocks cannot be allocated for " +
          std::to_string(block_limit) + " blocks with " + std::to_string(atom_limit) + " atoms.",
          "CacheResource");
  }
  const int per_item = block_limit * atom_limit;
  x_coordinates.setPointer(&llint_data,            0, per_item);
  y_coordinates.setPointer(&llint_data,     per_item, per_item);
  z_coordinates.setPointer(&llint_data, 2 * per_item, per_item);
  x_velocities.setPointer(&llint_data,  3 * per_item, per_item);
  y_velocities.setPointer(&llint_data,  4 * per_item, per_item);
  z_velocities.setPointer(&llint_data,  5 * per_item, per_item);
  x_force_overflow.setPointer(&int_data,                 0, per_item);
  y_force_overflow.setPointer(&int_data,          per_item, per_item);
  z_force_overflow.setPointer(&int_data,      2 * per_item, per_item);
  lennard_jones_indices.setPointer(&int_data, 3 * per_item, per_item);
}

//-------------------------------------------------------------------------------------------------
CacheResourceKit<double>
CacheResource::getDoublePrecisionCacheResourceKit(const HybridTargetLevel tier) {
  return CacheResourceKit<double>(block_limit, atom_limit, x_coordinates.data(tier),
                                  y_coordinates.data(tier), z_coordinates.data(tier),
                                  x_velocities.data(tier), y_velocities.data(tier),
                                  z_velocities.data(tier), x_force_overflow.data(tier),
                                  y_force_overflow.data(tier), z_force_overflow.data(tier),
                                  charges.data(tier), lennard_jones_indices.data(tier));
}

//-------------------------------------------------------------------------------------------------
CacheResourceKit<float>
CacheResource::getSinglePrecisionCacheResourceKit(const HybridTargetLevel tier) {
  return CacheResourceKit<float>(block_limit, atom_limit, x_coordinates.data(tier),
                                 y_coordinates.data(tier), z_coordinates.data(tier),
                                 x_velocities.data(tier), y_velocities.data(tier),
                                 z_velocities.data(tier), x_force_overflow.data(tier),
                                 y_force_overflow.data(tier), z_force_overflow.data(tier),
                                 sp_charges.data(tier), lennard_jones_indices.data(tier));
}

} // namespace energy
} // namespace omni
