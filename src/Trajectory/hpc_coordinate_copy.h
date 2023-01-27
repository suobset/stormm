// -*-c++-*-
#ifndef HPC_COORDINATE_COPY_H
#define HPC_COORDINATE_COPY_H

#include "copyright.h"
#include "Accelerator/gpu_details.h"
#include "Synthesis/condensate.h"
#include "Synthesis/phasespace_synthesis.h"
#include "coordinateframe.h"
#include "phasespace.h"

namespace stormm {
namespace trajectory {

using card::GpuDetails;
  
/// \brief Extract the coordinates of one object on the host or device to another on the host or
///        device.  This will launch a templated CUDA kernel to accomplish the copying.  This
///        functon trusts that the coordinate arrays it is being fed will be accessible to the
///        device.
///
/// \param xdest        Cartesian X coordinates for all particles in the destination system(s)
/// \param ydest        Cartesian X coordinates for all particles in the destination system(s)
/// \param zdest        Cartesian X coordinates for all particles in the destination system(s)
/// \param dest_scale   Scaling factor for the destination coordinates
/// \param dest_offset  Offset at which to begin writing coordinates for each atom to the
///                     destination
/// \param ct_dest      Specifier code for the destination data type
/// \param xorig        Cartesian X coordinates for all particles in the original system(s)
/// \param yorig        Cartesian Y coordinates for all particles in the original system(s)
/// \param zorig        Cartesian Z coordinates for all particles in the original system(s)
/// \param orig_scale   Scaling factor for the origin coordinates
/// \param orig_offset  Offset at which to begin reading coordinates for each atom from the origin
/// \param ct_orig      Specifier code for the origin data type
/// \param natom        The number of atoms' coordinates to copy
/// \param gpu          Details of the GPU to use
void launchCopyCoordinateXYZ(void* xdest, void* ydest, void* zdest, double dest_scale,
                             size_t dest_offset, size_t ct_dest, const void* xorig,
                             const void* yorig, const void* zorig, double orig_scale,
                             size_t orig_offset, size_t ct_orig, int natom, const GpuDetails &gpu);

} // namespace trajectory
} // namespace stormm

#endif
