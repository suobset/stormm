// -*-c++-*-
#include "copyright.h"
#include "coordinate_graft.h"
#include "hpc_coordinate_graft.cuh"

namespace stormm {
namespace trajectory {

//-------------------------------------------------------------------------------------------------
void launchGraftCoordXYZ(TrajectoryKind kind, void* xdest, void* ydest, void *zdest,
                         const double dest_scale, const void* xorig, const void* yorig,
                         const void* zorig, const double orig_scale, const int2* atom_list,
                         const int nxfer, const size_t dest_type_index,
                         const size_t orig_type_index) {
  if (dest_type_index == double_type_index) {
#if 0
    double* d_xdest = reinterpret_cast<double*>(xdest);
    double* d_ydest = reinterpret_cast<double*>(ydest);
    double* d_zdest = reinterpret_cast<double*>(zdest);
#endif
  }
  else if (dest_type_index == float_type_index) {

  }
  else if (dest_type_index == int_type_index) {
    
  }
  else if (dest_type_index == llint_type_index) {
    
  }
  else if (dest_type_index == short_type_index) {
    
  }
}

} // namespace trajectory
} // namespace stormm
