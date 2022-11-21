// -*-c++-*-
#ifndef STORMM_HPC_PHASE_SPACE_SYNTHESIS_CUH
#define STORMM_HPC_PHASE_SPACE_SYNTHESIS_CUH

#include <vector>
#ifdef STORMM_USE_CUDA
#include <cuda.h>
#include <cuda_runtime.h>
#endif
#include "copyright.h"
#include "Constants/hpc_bounds.h"
#include "phasespace_synthesis.h"

namespace stormm {
namespace synthesis {

using constants::large_block_size;

/// \brief Transfer a subset of system coordinates (and perhaps box dimensions) data from specific
///        systems within a PhaseSpaceSynthesis object betwee the host and the accelerator device.
///        PhaseSpaceSynthesis is a central object within STORMM and can be huge (multiple
///        gigabytes worth of coordinate, velocity, and force data.  As such, it deserves a
///        dedicated kernel for managing data transfer to and from the host.
///
/// \param destination  Collection of pointers to PhaseSpaceSynthesis data (this must be visible
///                     on the device, but could be device-only memory or host-mapped memory)
/// \param source       Collection of pointers to PhaseSpaceSynthesis data (this must be visible
///                     on the device, but could be device-only memory or host-mapped memory)
/// \param low_index    Lower bound of systems to upload
/// \param high_index   Upper bound of systems to upload (the range is [low_index, high_index))
__global__ void __launch_bounds__(large_block_size, 1)
kSystemTransfer(PsSynthesisWriter destination, PsSynthesisWriter source, int low_index,
                int high_index, const TrajectoryKind material);

/// \brief Initialize forces for one or all systems of a PhaseSpaceSynthesis on the GPU.
///
/// \param psyw   Writeable abstract for the PhaseSpaceSynthesis, containing system limits and
///               pointers to all coordinates and forces
/// \param index  Index of the system to initialize; if negative, all systems will be initialized.
__global__ void __launch_bounds__(large_block_size, 1)
kPsyInitializeForces(PsSynthesisWriter psyw, const int index);

/// \brief Initialize critical buffers in the phase space (specifically, prior coordinates and
///        velocities) which would otherwise not be used in energy minimization calculations.
///
/// \param psyw   Writeable abstract for the PhaseSpaceSynthesis
__global__ void __launch_bounds__(large_block_size, 1)
kPsyPrimeConjugateGradient(PsSynthesisWriter psyw);

/// \brief Import the Cartesian X, Y, and Z components of poositions, velocities, or forces of one
///        system within the synthesis.
///
/// \param x_recv             Long-long integer component of the fixed-precision arrays that the
///                           imported Cartesian X data shall replace
/// \param x_recv_ovrf        Integer overflow component of the fixed-precision arrays that the
///                           imported Cartesian X data shall replace
/// \param y_recv             Long-long integer component of Y fixed-precision arrays
/// \param y_recv_ovrf        Integer overflow component of Y fixed-precision arrays
/// \param z_recv             Long-long integer component of Z fixed-precision arrays
/// \param z_recv_ovrf        Integer overflow component of Z fixed-precision arrays
/// \param box_xform          Box space transformation matrix series within the synthesis (the
///                           matrices have units of inverse Angstroms)
/// \param inverse_xform      Inverse transformation matrix series within the synthesis (the
///                           matrices have units of Angstroms)
/// \param box_dimension      Box dimension series within the synthesis (units of Angstroms)
/// \param box_vectors        Primary bits for the fixed precision box vectors matrix
/// \param box_vector_ovrf    Overflow bits for the fixed precision box vectors matrix
/// \param atom_starts        Starting positions for the atoms of each system in the synthesis
/// \param atom_counts        Counts of atoms in each system of the synthesis
/// \param x_import           Input Cartesian X coordinates (positions, velocities, or forces)
/// \param y_import           Input Cartesian Y coordinates
/// \param z_import           Input Cartesian Z coordinates
/// \param box_xform_in       Transformation matrix to take coordinates into fractional space  
///                           (for positions only--provide nullptr for velocities or forces)
/// \param inverse_xform_in   Transformation matrix to take coordinates back to real space.  The
///                           units of elements in this matrix are Angstroms.
/// \param box_dimensions_in  Dimensions of the box, in internal units of Angstroms
/// \param system_index       Index of the system within this synthesis that the imported
///                           coordinates shall replace
/// \param coversion_factor   Scaling factor to take the input X, Y, and Z data into the apprpriate
///                           fixed-precision scale for incorporation into the synthesis
template <typename T>
__global__ void __launch_bounds__(large_block_size, 1)
kPsyImportSystemData(llint* x_recv, int* x_recv_ovrf, llint* y_recv, int* y_recv_ovrf,
                     llint* z_recv, int* z_recv_ovrf, double* box_xform, double* inverse_xform,
                     double* box_dimensions, llint* box_vectors, int* box_vector_ovrf,
                     const int* atom_starts, const int* atom_counts, const T* x_import,
                     const T* y_import, const T* z_import, const double* box_xform_in,
                     const double* inverse_xform_in, const double* box_dimensions_in,
                     const int system_index, const double conversion_factor) {

  // Transfer the box information, if it exists
  switch (kind) {
  case TrajectoryKind::POSITIONS:
    if (threadIdx.x < 9) {
      box_xform[box_offset + threadIdx.x] = box_xform_in[threadIdx.x];
      inverse_xform[box_offset + threadIdx.x] = inverse_xform_in[threadIdx.x];
    }
    if (threadIdx.x >= warp_size_int && threadIdx.x < warp_size_int + 9) {
      const int95_t fpbv = doubleToInt95(inverse_xform_in[threadIdx.x] * globalpos_scale);
      box_vectors[box_offset + threadIdx.x] = fpbv.x;
      box_vector_ovrf[box_offset + threadIdx.x] = fpbv.y;
    }
    if (threadIdx.x >= twice_warp_size_int && threadIdx.x < twice_warp_size_int + 6) {
      box_dim_ptr[dim_offset + threadIdx.x] = box_dimensions_in[threadIdx.x];
    }
    break;
  case TrajectoryKind::VELOCITIES:
  case TrajectoryKind::FORCES:
    break;
  }

  // Transfer atomic data
  const int pos_start = psyw.atom_starts[system_index];
  const int pos_end   = pos_start + psyw.atom_counts[system_index];
  const int stride = pos_end - pos_start;
  const int padded_stride = devcRoundUp(stride, warp_size_int);
  int pos = threadIdx.x;
  while (pos < padded_stride) {
    if (pos < stride) {
      const int95_t fpx = doubleToInt95((double)(x_import[pos]) * conversion_factor);
      const size_t ip = pos + pos_start;
      x_recv[ip]      = fpx.x;
      x_recv_ovrf[ip] = fpx.y;
    }
    pos += blockDim.x * gridDim.x;
  }
  while (pos < 2 * padded_stride) {
    const int rel_pos = pos - padded_stride;
    if (rel_pos < padded_stride) {
      const int95_t fpy = doubleToInt95((double)(y_import[rel_pos]) * conversion_factor);
      const size_t ip = rel_pos + pos_start;
      y_recv[ip]      = fpy.x;
      y_recv_ovrf[ip] = fpy.y;
    }
    pos += blockDim.x * gridDim.x;
  }
  while (pos < 3 * padded_stride) {
    const int rel_pos = pos - (2 * padded_stride);
    if (rel_pos < padded_stride) {
      const int95_t fpz = doubleToInt95((double)(z_import[rel_pos]) * conversion_factor);
      const size_t ip = rel_pos + pos_start;
      z_recv[ip]      = fpz.x;
      z_recv_ovrf[ip] = fpz.y;
    }
    pos += blockDim.x * gridDim.x;
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void launchPsyImportSystemData(llint* x_recv, int* x_recv_ovrf, llint* y_recv, int* y_recv_ovrf,
                               llint* z_recv, int* z_recv_ovrf, double* box_xform,
                               double* inverse_xform, double* box_dimensions, llint* box_vectors,
                               int* box_vector_ovrf, const int* atom_starts,
                               const int* atom_counts, const T* x_import, const T* y_import,
                               const T* z_import, const double* box_xform_in,
                               const double* inverse_xform_in, const double* box_dimensions_in,
                               const int system_index, const double conversion_factor,
                               const GpuDetails &gpu) {
  kPsyImportSystemData<<<large_block_size,
                         gpu.getSMPCount()>>>(x_recv, x_recv_ovrf, y_recv, y_recv_ovrf, z_recv,
                                              z_recv_ovrf, box_xform, inverse_xform,
                                              box_dimensions, box_vectors, box_vector_ovrf,
                                              atom_starts, atom_counts, x_import, y_import,
                                              z_import, box_xform_in, inverse_xform_in,
                                              box_dimensions_in, system_index, conversion_factor);
}

} // namespace synthesis
} // namespace stormm

#endif
