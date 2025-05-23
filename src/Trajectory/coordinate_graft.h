// -*-c++-*-
#ifndef STORMM_COORDINATE_GRAFT_H
#define STORMM_COORDINATE_GRAFT_H

#include <vector>
#include "copyright.h"
#include "Accelerator/card_utilities.h"
#include "Accelerator/gpu_details.h"
#include "Accelerator/gpu_enumerators.h"
#include "Numerics/split_fixed_precision.h"
#include "Synthesis/phasespace_synthesis.h"
#include "Synthesis/condensate.h"
#include "coordinateframe.h"
#include "coordinate_util.h"
#include "phasespace.h"
#include "coordinate_series.h"
#include "trajectory_enumerators.h"

namespace stormm {
namespace trajectory {

using card::GpuDetails;
using card::HpcKernelSync;
using card::Hybrid;
using card::HybridTargetLevel;
using numerics::hostChangeFPBits;
using numerics::hostDoubleToInt95;
using numerics::hostInt95ToDouble;
using synthesis::Condensate;
using synthesis::CondensateReader;
using synthesis::CondensateWriter;
using synthesis::PhaseSpaceSynthesis;
using synthesis::PsSynthesisReader;
using synthesis::PsSynthesisWriter;
  
#ifdef STORMM_USE_HPC
/// \brief Launch the appropriate kernel to perform coordinate grafting using a GPU.  This makes
///        use of templated unrolling based on the destination type, similar to the coordinate
///        copying functions (see hpc_coordinate_copy.cu).
///
/// \param xdest            Cartesian X coordinates from the destination structure
/// \param ydest            Cartesian Y coordinates from the destination structure
/// \param zdest            Cartesian Z coordinates from the destination structure
/// \param dest_scale       The scaling factor of coordinates in the origin, if the origin is
///                         represented in fixed-precision
/// \param xorig            Cartesian X coordinates from the original structure
/// \param yorig            Cartesian Y coordinates from the original structure
/// \param zorig            Cartesian Z coordinates from the original structure
/// \param orig_scale       The scaling factor of coordinates in the origin, if the origin is
///                         represented in fixed-precision
/// \param atom_list        The list of atoms to transfer, with destination atom indices given in
///                         the "x" member of each tuple and origin atom indices given in the "y"
///                         member
/// \param nxfer            The total number of atoms with coordinates to transfer
/// \param dest_type_index  Type index of the destination coordinate set (see the runtime-derived
///                         values in src/DataTypes/common_types.h)
/// \param orig_type_index  Type index of the original coordinate set
void launchGraftCoordXYZ(TrajectoryKind kind, void* xdest, void* ydest, void *zdest,
                         double dest_scale, const void* xorig, const void* yorig,
                         const void* zorig, double orig_scale, const int2* atom_list,
                         const int nxfer, size_t dest_type_index, size_t orig_type_index);
#endif
  
/// \brief The low-level implementation to graft coordinates (whether positions, velocities, or
///        forces) with Cartesian X, Y, and Z components from selected atoms of one object onto
///        selected atoms of another.  This is a general implementation with some arguments set to
///        defaults of null pointers to avoid the need for overloading.
///
/// Overloaded:
///   - A templated form to accept any standard pointers to integer or real data
///   - A special case for int95_t split fixed-precision data.  This avoids compiler warnings that
///     would occur if the templated version were required to handle the overflow bit arrays and
///     typecast real numbers to long long integers in int95_t initialization.
///
/// \param xdest       Array of Cartesian X coordinates in the destination structure
/// \param ydest       Array of Cartesian Y coordinates in the destination structure
/// \param zdest       Array of Cartesian Z coordinates in the destination structure
/// \param dest_scale  The scaling factor by which coordinates in the destination structure are
///                    represented, if in fixed precision
/// \param xorig       Array of Cartesian X coordinates in the origin structure
/// \param yorig       Array of Cartesian Y coordinates in the origin structure
/// \param zorig       Array of Cartesian Z coordinates in the origin structure
/// \param orig_scale  The scaling factor by which coordinates in the origin structure are 
///                    represented, if in fixed precision
/// \param atom_list   List of atom index correspondence pairs, tuples giving the detinatio index
///                    in the "x" member and the origin index in the "y" member
/// \param nxfer       The trusted length of atom_list
/// \param xdest_ovrf  Array of elements with overflow bits for xdest.  The presence of this array
///                    as a non-null pointer implies that coordinates are presented in split
///                    fixed-precison integers.
/// \param ydest_ovrf  Array of elements with overflow bits for ydest
/// \param zdest_ovrf  Array of elements with overflow bits for zdest
/// \param xorig_ovrf  Array of elements with overflow bits for xorig
/// \param yorig_ovrf  Array of elements with overflow bits for yorig
/// \param zorig_ovrf  Array of elements with overflow bits for zorig
/// \{
template <typename Tdest, typename Torig>
void graftCoordinateXYZ(TrajectoryKind kind, Tdest* xdest, Tdest* ydest, Tdest* zdest,
                        const double dest_scale, const Torig* xorig, const Torig* yorig,
                        const Torig* zorig, const double orig_scale, const int2* atom_list,
                        const int nxfer);

void graftCoordinateXYZ(const TrajectoryKind kind, llint* xdest, llint* ydest, llint* zdest,
                        const double dest_scale, const llint* xorig, const llint* yorig,
                        const llint* zorig, const double orig_scale, const int2* atom_list,
                        const int nxfer, int* xdest_ovrf, int* ydest_ovrf, int* zdest_ovrf,
                        const int* xorig_ovrf, const int* yorig_ovrf, const int* zorig_ovrf);
/// \}

/// \brief Graft (some or all) coordinates from one object onto another.  The coordinates might be
///        positions, velocities, or forces.  The direction of the graft is given by the object
///        that comes first (the destination object).
///
/// \param dest        The destination structure
/// \param frame_dest  Frame of the destination object, valid in cases when the destination is a
///                    CoordinateSeries
/// \param orig        The original structure
/// \param frame_orig  Frame of the original object, valid in cases when the origin is a
///                    CoordinateSeries
/// \param atom_list   The list of atoms to transfer.  Destination atom indices are indicated by
///                    the "x" members of each tuple, origin atom indices by the "y" members of
///                    each tuple.
/// \param dest_tier   Indicate whether to load coordinates into the destination on the CPU host
///                    or GPU device memory
/// \param orig_tier   Indicate whether to take coordinates from memory in the origin object on the
///                    CPU host or GPU device
/// \param gpu         Details of the GPU that may carry out the transfer
/// \param sync        The synchronization protocol between the CPU and GPU
/// \{
void coordGraft(CoordinateFrameWriter *dest, const CoordinateFrameReader &orig,
                const int2* atom_list, int nxfer,
                HybridTargetLevel dest_tier = HybridTargetLevel::HOST,
                HybridTargetLevel orig_tier = HybridTargetLevel::HOST,
                const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

void coordGraft(CoordinateFrame *dest, const CoordinateFrame &orig, const Hybrid<int2> &atom_list,
                HybridTargetLevel dest_tier = HybridTargetLevel::HOST,
                HybridTargetLevel orig_tier = HybridTargetLevel::HOST,
                const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

void coordGraft(CoordinateFrameWriter *dest, const PhaseSpaceReader &orig, const int2* atom_list,
                int nxfer, TrajectoryKind kind,
                HybridTargetLevel dest_tier = HybridTargetLevel::HOST,
                HybridTargetLevel orig_tier = HybridTargetLevel::HOST,
                const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

void coordGraft(CoordinateFrame *dest, const PhaseSpace &orig, const Hybrid<int2> &atom_list,
                TrajectoryKind kind, CoordinateCycle orientation,
                HybridTargetLevel dest_tier = HybridTargetLevel::HOST,
                HybridTargetLevel orig_tier = HybridTargetLevel::HOST,
                const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

void coordGraft(CoordinateFrame *dest, const PhaseSpace &orig, const Hybrid<int2> &atom_list,
                TrajectoryKind kind, HybridTargetLevel dest_tier = HybridTargetLevel::HOST,
                HybridTargetLevel orig_tier = HybridTargetLevel::HOST,
                const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

void coordGraft(CoordinateFrameWriter *dest, const PhaseSpaceReader &orig, const int2* atom_list,
                int nxfer, HybridTargetLevel dest_tier = HybridTargetLevel::HOST,
                HybridTargetLevel orig_tier = HybridTargetLevel::HOST,
                const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

void coordGraft(CoordinateFrame *dest, const PhaseSpace &orig, const Hybrid<int2> &atom_list,
                HybridTargetLevel dest_tier = HybridTargetLevel::HOST,
                HybridTargetLevel orig_tier = HybridTargetLevel::HOST,
                const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

template <typename Torig>
void coordGraft(CoordinateFrameWriter *dest, const CoordinateSeriesReader<Torig> &orig,
                const int2* atom_list, int nxfer,
                HybridTargetLevel dest_tier = HybridTargetLevel::HOST,
                HybridTargetLevel orig_tier = HybridTargetLevel::HOST,
                const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

template <typename Torig>
void coordGraft(CoordinateFrame *dest, const CoordinateSeries<Torig> &orig,
                const Hybrid<int2> &atom_list,
                HybridTargetLevel dest_tier = HybridTargetLevel::HOST,
                HybridTargetLevel orig_tier = HybridTargetLevel::HOST,
                const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

void coordGraft(CoordinateFrameWriter *dest, const PsSynthesisReader &orig,
                const int2* atom_list, int nxfer, TrajectoryKind kind,
                HybridTargetLevel dest_tier = HybridTargetLevel::HOST,
                HybridTargetLevel orig_tier = HybridTargetLevel::HOST,
                const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

void coordGraft(CoordinateFrame *dest, const PhaseSpaceSynthesis &orig,
                const Hybrid<int2> &atom_list, TrajectoryKind kind, CoordinateCycle orientation,
                HybridTargetLevel dest_tier = HybridTargetLevel::HOST,
                HybridTargetLevel orig_tier = HybridTargetLevel::HOST,
                const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

void coordGraft(CoordinateFrameWriter *dest, const PsSynthesisReader &orig, const int2* atom_list,
                int nxfer, HybridTargetLevel dest_tier = HybridTargetLevel::HOST,
                HybridTargetLevel orig_tier = HybridTargetLevel::HOST,
                const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

void coordGraft(CoordinateFrame *dest, const PhaseSpaceSynthesis &orig,
                const Hybrid<int2> &atom_list, TrajectoryKind kind,
                HybridTargetLevel dest_tier = HybridTargetLevel::HOST,
                HybridTargetLevel orig_tier = HybridTargetLevel::HOST,
                const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

void coordGraft(CoordinateFrame *dest, const PhaseSpaceSynthesis &orig,
                const Hybrid<int2> &atom_list,
                HybridTargetLevel dest_tier = HybridTargetLevel::HOST,
                HybridTargetLevel orig_tier = HybridTargetLevel::HOST,
                const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

void coordGraft(CoordinateFrameWriter *dest, const CondensateReader &orig, const int2* atom_list,
                int nxfer, HybridTargetLevel dest_tier = HybridTargetLevel::HOST,
                HybridTargetLevel orig_tier = HybridTargetLevel::HOST,
                const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

void coordGraft(CoordinateFrame *dest, const Condensate &orig, const Hybrid<int2> &atom_list,
                HybridTargetLevel dest_tier = HybridTargetLevel::HOST,
                HybridTargetLevel orig_tier = HybridTargetLevel::HOST,
                const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);
/// \}

} // namespace trajectory
} // namespace stormm

#include "coordinate_graft.tpp"

#endif
