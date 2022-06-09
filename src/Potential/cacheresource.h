// -*-c++-*-
#ifndef OMNI_CACHERESOURCE_H
#define OMNI_CACHERESOURCE_H

#include "Accelerator/hybrid.h"
#include "DataTypes/common_types.h"

namespace omni {
namespace energy {

using card::Hybrid;
using card::HybridTargetLevel;

/// \brief Abstract for the CacheResource object, accessible as a C-style struct and suitable for
///        passing to GPU kernels as a kernel argument.
template <typename T> struct CacheResourceKit {
  llint* xcrd;           // Cartesian X coordinates of locally cached particles
  llint* ycrd;           // Cartesian Y coordinates of locally cached particles
  llint* zcrd;           // Cartesian Z coordinates of locally cached particles
  llint* xvel;           // Cartesian X velocities of locally cached particles
  llint* yvel;           // Cartesian Y velocities of locally cached particles
  llint* zvel;           // Cartesian Z velocities of locally cached particles
  llint* xfrc_overflow;  // Cartesian X force overflow buffers
  llint* yfrc_overflow;  // Cartesian Y force overflow buffers
  llint* zfrc_overflow;  // Cartesian Z force overflow buffers
  T* charges;            // Charge parameters for locally cached particles
  int* lj_idx;           // Lennard-Jones indices of locally cached particles
};

/// \brief An object to hold temporary data for a particular work unit (whether bonded or
///        non-bonded), resident in GMEM but private to a particular thread block.  This object
///        must be allocated in such a way as to be ready to hold the private workspaces of any
///        thread blocks that will make use of it.  The object is arranged to support operations
///        based on PhaseSpaceSynthesis respresentations of systems' coordinates and forces.
class CacheResource {
public:

  /// \brief The constructor does not take a GPU description, but instead maximum numbers of
  ///        blocks and atoms per block that might be required, allocating space as appropriate.
  CacheResource(const int max_block_count, const int max_atom_count);

  /// \brief Get a set of pointers to this object with single-precision representations for the
  ///        charges.
  CacheResourceKit<float>
  getSinglePrecisionCacheResourceKit(HybridTargetLevel tier = HybridTargetLevel::HOST);

  /// \brief Get a set of pointers to this object with double-precision representations for the
  ///        charges.
  CacheResourceKit<double>
  getDoublePrecisionCacheResourceKit(HybridTargetLevel tier = HybridTargetLevel::HOST);

private:
  Hybrid<llint> x_coordinates;        ///< Cartesian X coordinates of locally cached particles
  Hybrid<llint> y_coordinates;        ///< Cartesian Y coordinates of locally cached particles
  Hybrid<llint> z_coordinates;        ///< Cartesian Z coordinates of locally cached particles
  Hybrid<llint> x_velocities;         ///< Cartesian X velocities of locally cached particles
  Hybrid<llint> y_velocities;         ///< Cartesian Y velocities of locally cached particles
  Hybrid<llint> z_velocities;         ///< Cartesian Z velocities of locally cached particles
  Hybrid<llint> x_force_overflow;     ///< Cartesian X force overflow buffers
  Hybrid<llint> y_force_overflow;     ///< Cartesian Y force overflow buffers
  Hybrid<llint> z_force_overflow;     ///< Cartesian Z force overflow buffers
  Hybrid<double> charges;             ///< Charge parameters for locally cached particles
  Hybrid<float> sp_charges;           ///< Charge parameters for locally cached particles (single
                                      ///<   precision)
  Hybrid<int> lennard_jones_indices;  ///< Lennard-Jones indices of locally cached particles  
};

} // namespace energy
} // namespace omni

#endif
