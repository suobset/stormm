// -*-c++-*-
#ifndef STORMM_COORDINATE_COPY_H
#define STORMM_COORDINATE_COPY_H

#include "Accelerator/card_utilities.h"
#include "Accelerator/gpu_details.h"
#include "Accelerator/gpu_enumerators.h"
#include "Constants/behavior.h"
#include "DataTypes/common_types.h"
#include "DataTypes/stormm_vector_types.h"
#include "Math/matrix_ops.h"
#include "Math/rounding.h"
#include "Numerics/split_fixed_precision.h"
#include "Synthesis/condensate.h"
#include "Synthesis/phasespace_synthesis.h"
#include "Synthesis/synthesis_enumerators.h"
#include "coordinateframe.h"
#include "coordinate_series.h"
#include "phasespace.h"
#include "trajectory_enumerators.h"
#ifdef STORMM_USE_HPC
#  include "hpc_coordinate_copy.h"
#endif

namespace stormm {
namespace trajectory {

using card::GpuDetails;
using card::HpcKernelSync;
using card::HybridTargetLevel;
using constants::PrecisionModel;
using data_types::isFloatingPointScalarType;
using data_types::isSignedIntegralScalarType;
using data_types::getStormmScalarTypeName;
using stmath::extractBoxDimensions;
using stmath::roundUp;
using numerics::hostDoubleToInt95;
using numerics::hostChangeFPBits;
using synthesis::Condensate;
using synthesis::CondensateReader;
using synthesis::StructureSource;
using synthesis::CondensateWriter;
using synthesis::PhaseSpaceSynthesis;
using synthesis::PsSynthesisReader;
using synthesis::PsSynthesisWriter;

/// \brief Validate the number of atoms in two different coordinate objects and report an error
///        if they disagree, appropriate to call attention to one of the coordCopy() overloads.
///
/// \param destination_atoms  The number of atoms expected by the destinaton object
/// \param origin_atoms       The number of atoms available in the origin object
void coordCopyValidateAtomCounts(int destination_atoms, int origin_atoms);

/// \brief Validate the number of frames in a coordinate series against the one requested.  Produce
///        an error indicating the source of the problem in one of the coordCopy() overloads.
///
/// \param frame_index  The frame requested (to write into, or read from)
/// \param frame_count  The frame count in the series
void coordCopyValidateFrameIndex(int frame_index, int frame_count);

/// \brief Validate the number of systems in a synthesis, against the index requested.  Produce an
///        error indicating the source of the problem in one of the coordCopy() overloads.
///
/// \param system_index  The system requested (to write into, or read from)
/// \param system_count  The system count in the synthesis
void coordCopyValidateSystemIndex(int system_index, int system_count);

/// \brief Copy the box information from one object to another.
///
/// Overloaded:
///   - Accept all box dimension arrays individually
///   - Accept Condensate abstracts, implying that the six-membered box dimensions array is not
///     needed or present
///
/// \param dest_boxdim  Destination box dimensions
/// \param dest_umat    Destination unit cell transformation matrix (takes coordinates into box
///                     space)
/// \param dest_invu    Destination inverse transformation matrix (takes coordinates back into
///                     real space)
/// \param index_dest   Index of the system or frame in the destination coordinate object
/// \param orig_boxdim  Original object box dimensions
/// \param orig_umat    Original unit cell transformation matrix (takes coordinates into box
///                     space)
/// \param orig_invu    Original inverse transformation matrix (takes coordinates back into real
///                     space)
/// \param index_orig   Index of the system or frame in the original coordinate object
/// \param boxvecs      Fixed-precision box vectors used by PhaseSpaceSynthesis objects
/// \param boxvec_ovrf  Overflow bits for boxvecs
/// \{
void copyBoxInformation(double* dest_boxdim, double* dest_umat, double* dest_invu,
                        int index_dest, const double* orig_boxdim, const double* orig_umat,
                        const double* orig_invu, int index_orig, llint* boxvecs = nullptr,
                        int* boxvec_ovrf = nullptr);

void copyBoxInformation(double* dest_boxdim, double* dest_umat, double* dest_invu,
                        int index_dest, const double* orig_boxdim, const double* orig_umat,
                        const double* orig_invu, llint* boxvecs = nullptr,
                        int* boxvec_ovrf = nullptr);

void copyBoxInformation(double* dest_boxdim, double* dest_umat, double* dest_invu,
                        const double* orig_boxdim, const double* orig_umat,
                        const double* orig_invu, int index_orig = 0);
/// \}

/// \brief Copy one set of XYZ coordinates, trusting that both share the same number of atoms.
///
/// Overloaded:
///   - Copy real-valued data with no conversion (types will be checked)
///   - Copy real-valued data to what could be fixed-precision format (with a conversion on the
///     destination)
///   - Copy what could be fixed-precision data to real-valued data (with a conversion on the
///     origin)
///   - Copy what could be fixed-precision data to other fixed-precision data (with conversions on
///     both ends)
///   - Copy 95-bit integer data to any other format, including 95-bit integer data in a different
///     scaling factor
///
/// \param xdest        Destination object Cartesian X coordinates
/// \param ydest        Destination object Cartesian Y coordinates
/// \param zdest        Destination object Cartesian Z coordinates
/// \param xorig        Origin object Cartesian X coordinates
/// \param yorig        Origin object Cartesian Y coordinates
/// \param zorig        Origin object Cartesian Z coordinates
/// \param dest_scale   Scaling factor applicable to coordinates in the destination arrays
/// \param orig_scale   Scaling factor applicable to coordinates in the origin arrays
/// \param dest_offset  Offset for atoms in each of the destination arrays
/// \param orig_offset  Offset for atoms in each of the origin arrays
/// \param natom        Trusted number of atoms
/// \param dest_tier    
/// \{
template <typename Tdest, typename Torig>
void copyCoordinateXYZ(Tdest* xdest, Tdest* ydest, Tdest* zdest, const Torig* xorig,
                       const Torig* yorig, const Torig* zorig, int natom);

template <typename Tdest, typename Torig>
void copyCoordinateXYZ(Tdest* xdest, Tdest* ydest, Tdest* zdest, double dest_scale,
                       const Torig* xorig, const Torig* yorig, const Torig* zorig,
                       double orig_scale, int natom);

template <typename Tdest, typename Torig>
void copyCoordinateXYZ(Tdest* xdest, Tdest* ydest, Tdest* zdest, const Torig* xorig,
                       const Torig* yorig, const Torig* zorig, double orig_scale, int natom);

template <typename Tdest, typename Torig>
void copyCoordinateXYZ(Tdest* xdest, Tdest* ydest, Tdest* zdest, double dest_scale,
                       const Torig* xorig, const Torig* yorig, const Torig* zorig,
                       double orig_scale, int natom);

template <typename Tdest>
void copyCoordinateXYZ(Tdest* xdest, Tdest* ydest, Tdest* zdest,
                       const llint* xorig, const int* xorig_ovrf, const llint* yorig,
                       const int* yorig_ovrf, const llint* zorig, const int* zorig_ovrf,
                       double orig_scale, int natom);

template <typename Tdest>
void copyCoordinateXYZ(Tdest* xdest, Tdest* ydest, Tdest* zdest, double dest_scale,
                       const llint* xorig, const int* xorig_ovrf, const llint* yorig,
                       const int* yorig_ovrf, const llint* zorig, const int* zorig_ovrf,
                       double orig_scale, int natom);

template <typename Torig>
void copyCoordinateXYZ(llint* xdest, int* xdest_ovrf, llint* ydest, int* ydest_ovrf, llint* zdest,
                       int* zdest_ovrf, double dest_scale, const Torig* xorig, const Torig* yorig,
                       const Torig* zorig, int natom);

template <typename Torig>
void copyCoordinateXYZ(llint* xdest, int* xdest_ovrf, llint* ydest, int* ydest_ovrf, llint* zdest,
                       int* zdest_ovrf, double dest_scale, const Torig* xorig, const Torig* yorig,
                       const Torig* zorig, double orig_scale, int natom);

void copyCoordinateXYZ(llint* xdest, int* xdest_ovrf, llint* ydest, int* ydest_ovrf, llint* zdest,
                       int* zdest_ovrf, double dest_scale, const llint* xorig,
                       const int* xorig_ovrf, const llint* yorig, const int* yorig_ovrf,
                       const llint* zorig, const int* zorig_ovrf, double orig_scale, int natom);
/// \}

/// \brief Make a deep copy of coordinates from one object to another.  The oder of arguments fed
///        to this function is inspired by the C standard library strcpy() function.  Here, the
///        goal is to copy strings of coordinates.  Checks will be made to ensure that the two
///        systems have the same number of atoms.
///
/// Overloaded:
///   - Copy into a CoordinateFrame, a PhaseSpace object, a specific frame of a CoordinateSeries,
///     a particular system of a PhaseSpaceSynthesis object, or a particular system of a
///     Condensate object.
///   - Copy from any of the above.
///   - Provide pointers and const references to either object, or abstracts for each of them.
///
/// \param destination      Pre-allocated object in which to put the copied coordinates (modified
///                         and returned)
/// \param kind             The type of coordinates (positions, velocities or forces) to assign in
///                         the destination, if the origin might be ambiguous
/// \param orientation      Point in the time cycle to copy the coordinates into, if applicable
/// \param origin           The source of coordinates to be copied
/// \param frame_dest       Frame number of some coordinate series into which coordinates shall be
///                         placed
/// \param frame_orig       Frame number of some coordinate series from which coordinates shall be
///                         copied
/// \param index_dest       System index of some synthesis object into which coordinates shall be
///                         placed
/// \param index_orig       System index of some synthesis object from which coordinates shall be
///                         copied
/// \param dest_atom_start  Starting atom index in the destination array (required if only
///                         abstracts are provided for multi-structure objects, as the abstracts
///                         may not be valid on the host to retrieve this information)
/// \param orig_atom_start  Starting atom index in the origin array
/// \param natom            The number of atoms in the system to copy (required for multi-structure
///                         objects when only abstracts are provided, as the abstracts may not be
///                         valid on the host to retrieve this information)
/// \{
void coordCopy(CoordinateFrameWriter *destination, const CoordinateFrameReader &origin,
               HybridTargetLevel destination_tier = HybridTargetLevel::HOST,
               HybridTargetLevel origin_tier = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

void coordCopy(CoordinateFrame *destination, const CoordinateFrame &origin,
               HybridTargetLevel destination_tier = HybridTargetLevel::HOST,
               HybridTargetLevel origin_tier = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

void coordCopy(CoordinateFrameWriter *destination, const PhaseSpaceReader &origin,
               TrajectoryKind kind, HybridTargetLevel destination_tier = HybridTargetLevel::HOST,
               HybridTargetLevel origin_tier = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

void coordCopy(CoordinateFrameWriter *destination, const PhaseSpaceReader &origin,
               HybridTargetLevel destination_tier = HybridTargetLevel::HOST,
               HybridTargetLevel origin_tier = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

void coordCopy(CoordinateFrame *destination, const PhaseSpace &origin, TrajectoryKind kind,
               CoordinateCycle orientation,
               HybridTargetLevel destination_tier = HybridTargetLevel::HOST,
               HybridTargetLevel origin_tier = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

void coordCopy(CoordinateFrame *destination, const PhaseSpace &origin, TrajectoryKind kind,
               HybridTargetLevel destination_tier = HybridTargetLevel::HOST,
               HybridTargetLevel origin_tier = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

void coordCopy(CoordinateFrame *destination, const PhaseSpace &origin,
               HybridTargetLevel destination_tier = HybridTargetLevel::HOST,
               HybridTargetLevel origin_tier = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

template <typename T>
void coordCopy(CoordinateFrame *destination, const CoordinateSeries<T> &origin,
               int frame_orig, HybridTargetLevel destination_tier = HybridTargetLevel::HOST,
               HybridTargetLevel origin_tier = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

void coordCopy(CoordinateFrameWriter *destination, const PsSynthesisReader &origin,
               int orig_atom_start, int index_orig, TrajectoryKind kind,
               HybridTargetLevel destination_tier = HybridTargetLevel::HOST,
               HybridTargetLevel origin_tier = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

void coordCopy(CoordinateFrame *destination, const PhaseSpaceSynthesis &origin, int index_orig,
               TrajectoryKind kind, CoordinateCycle orientation,
               HybridTargetLevel destination_tier = HybridTargetLevel::HOST,
               HybridTargetLevel origin_tier = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

void coordCopy(CoordinateFrame *destination, const PhaseSpaceSynthesis &origin, int index_orig,
               TrajectoryKind kind,
               HybridTargetLevel destination_tier = HybridTargetLevel::HOST,
               HybridTargetLevel origin_tier = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

void coordCopy(CoordinateFrame *destination, const PhaseSpaceSynthesis &origin, int index_orig,
               HybridTargetLevel destination_tier = HybridTargetLevel::HOST,
               HybridTargetLevel origin_tier = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

void coordCopy(CoordinateFrameWriter *destination, const CondensateReader &origin,
               int orig_atom_start, int index_orig,
               HybridTargetLevel destination_tier = HybridTargetLevel::HOST,
               HybridTargetLevel origin_tier = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

void coordCopy(CoordinateFrame *destination, const Condensate &origin, int index_orig,
               HybridTargetLevel destination_tier = HybridTargetLevel::HOST,
               HybridTargetLevel origin_tier = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);
  
void coordCopy(PhaseSpaceWriter *destination, TrajectoryKind kind,
               const CoordinateFrameReader &origin,
               HybridTargetLevel destination_tier = HybridTargetLevel::HOST,
               HybridTargetLevel origin_tier = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

void coordCopy(PhaseSpaceWriter *destination, const CoordinateFrameReader &origin,
               HybridTargetLevel destination_tier = HybridTargetLevel::HOST,
               HybridTargetLevel origin_tier = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

void coordCopy(PhaseSpace *destination, TrajectoryKind kind, CoordinateCycle orientation,
               const CoordinateFrame &origin,
               HybridTargetLevel destination_tier = HybridTargetLevel::HOST,
               HybridTargetLevel origin_tier = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

void coordCopy(PhaseSpace *destination, TrajectoryKind kind, const CoordinateFrame &origin,
               HybridTargetLevel destination_tier = HybridTargetLevel::HOST,
               HybridTargetLevel origin_tier = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

void coordCopy(PhaseSpace *destination, const CoordinateFrame &origin,
               HybridTargetLevel destination_tier = HybridTargetLevel::HOST,
               HybridTargetLevel origin_tier = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

void coordCopy(PhaseSpaceWriter *destination, const PhaseSpaceReader &origin,
               HybridTargetLevel destination_tier = HybridTargetLevel::HOST,
               HybridTargetLevel origin_tier = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

void coordCopy(PhaseSpace *destination, const PhaseSpace &origin,
               HybridTargetLevel destination_tier = HybridTargetLevel::HOST,
               HybridTargetLevel origin_tier = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

template <typename T>
void coordCopy(PhaseSpace *destination, TrajectoryKind kind, CoordinateCycle orientation,
               const CoordinateSeries<T> &origin, int frame_orig,
               HybridTargetLevel destination_tier = HybridTargetLevel::HOST,
               HybridTargetLevel origin_tier = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

template <typename T>
void coordCopy(PhaseSpace *destination, TrajectoryKind kind, const CoordinateSeries<T> &origin,
               int frame_orig, HybridTargetLevel destination_tier = HybridTargetLevel::HOST,
               HybridTargetLevel origin_tier = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

template <typename T>
void coordCopy(PhaseSpace *destination, const CoordinateSeries<T> &origin, int frame_orig,
               HybridTargetLevel destination_tier = HybridTargetLevel::HOST,
               HybridTargetLevel origin_tier = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

void coordCopy(PhaseSpaceWriter *destination, const PsSynthesisReader &origin, int orig_atom_start,
               int index_orig, HybridTargetLevel destination_tier = HybridTargetLevel::HOST,
               HybridTargetLevel origin_tier = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

void coordCopy(PhaseSpace *destination, const PhaseSpaceSynthesis &origin, int index_orig,
               HybridTargetLevel destination_tier = HybridTargetLevel::HOST,
               HybridTargetLevel origin_tier = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);
  
void coordCopy(PhaseSpaceWriter *destination, TrajectoryKind kind, const CondensateReader &origin,
               int orig_atom_start, int index_orig,
               HybridTargetLevel destination_tier = HybridTargetLevel::HOST,
               HybridTargetLevel origin_tier = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

void coordCopy(PhaseSpaceWriter *destination, const CondensateReader &origin, int orig_atom_start,
               int index_orig, HybridTargetLevel destination_tier = HybridTargetLevel::HOST,
               HybridTargetLevel origin_tier = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

void coordCopy(PhaseSpace *destination, TrajectoryKind kind, CoordinateCycle orientation,
               const Condensate &origin, int index_orig,
               HybridTargetLevel destination_tier = HybridTargetLevel::HOST,
               HybridTargetLevel origin_tier = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

void coordCopy(PhaseSpace *destination, TrajectoryKind kind, const Condensate &origin,
               int index_orig, HybridTargetLevel destination_tier = HybridTargetLevel::HOST,
               HybridTargetLevel origin_tier = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

void coordCopy(PhaseSpace *destination, const Condensate &origin, int index_orig,
               HybridTargetLevel destination_tier = HybridTargetLevel::HOST,
               HybridTargetLevel origin_tier = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

template <typename T>
void coordCopy(CoordinateSeries<T> *destination, int frame_dest,
               const CoordinateFrameReader &origin,
               HybridTargetLevel destination_tier = HybridTargetLevel::HOST,
               HybridTargetLevel origin_tier = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

template <typename T>
void coordCopy(CoordinateSeries<T> *destination, int frame_dest, const CoordinateFrame &origin,
               HybridTargetLevel destination_tier = HybridTargetLevel::HOST,
               HybridTargetLevel origin_tier = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

template <typename T>
void coordCopy(CoordinateSeries<T> *destination, int frame_dest, const PhaseSpaceReader &origin,
               const TrajectoryKind kind,
               HybridTargetLevel destination_tier = HybridTargetLevel::HOST,
               HybridTargetLevel origin_tier = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

template <typename T>
void coordCopy(CoordinateSeries<T> *destination, int frame_dest, const PhaseSpaceReader &origin,
               HybridTargetLevel destination_tier, HybridTargetLevel origin_tier,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

template <typename T>
void coordCopy(CoordinateSeries<T> *destination, int frame_dest, const PhaseSpace &origin,
               TrajectoryKind kind, CoordinateCycle orientation,
               HybridTargetLevel destination_tier = HybridTargetLevel::HOST,
               HybridTargetLevel origin_tier = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

template <typename T>
void coordCopy(CoordinateSeries<T> *destination, int frame_dest, const PhaseSpace &origin,
               TrajectoryKind kind, HybridTargetLevel destination_tier = HybridTargetLevel::HOST,
               HybridTargetLevel origin_tier = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

template <typename T>
void coordCopy(CoordinateSeries<T> *destination, int frame_dest, const PhaseSpace &origin,
               HybridTargetLevel destination_tier = HybridTargetLevel::HOST,
               HybridTargetLevel origin_tier = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

template <typename Tdest, typename Torig>
void coordCopy(CoordinateSeriesWriter<Tdest> *destination, size_t frame_dest,
               const CoordinateSeriesReader<Torig> &origin, size_t frame_orig,
               HybridTargetLevel destination_tier = HybridTargetLevel::HOST,
               HybridTargetLevel origin_tier = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

template <typename Tdest, typename Torig>
void coordCopy(CoordinateSeries<Tdest> *destination, size_t frame_dest,
               const CoordinateSeries<Torig> &origin, size_t frame_orig,
               HybridTargetLevel destination_tier = HybridTargetLevel::HOST,
               HybridTargetLevel origin_tier = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);
  
template <typename T>
void coordCopy(CoordinateSeriesWriter<T> *destination, size_t frame_dest,
               const PsSynthesisReader &origin, int orig_atom_start, int index_orig,
               TrajectoryKind kind = TrajectoryKind::POSITIONS,
               HybridTargetLevel destination_tier = HybridTargetLevel::HOST,
               HybridTargetLevel origin_tier = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

template <typename T>
void coordCopy(CoordinateSeries<T> *destination, int frame_dest,
               const PhaseSpaceSynthesis &origin, int index_orig,
               TrajectoryKind kind = TrajectoryKind::POSITIONS,
	       CoordinateCycle orientation = CoordinateCycle::PRIMARY,
               HybridTargetLevel destination_tier = HybridTargetLevel::HOST,
               HybridTargetLevel origin_tier = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

template <typename T>
void coordCopy(CoordinateSeriesWriter<T> *destination, size_t frame_dest,
               const CondensateReader &origin, size_t orig_atom_start, size_t index_orig,
               HybridTargetLevel destination_tier = HybridTargetLevel::HOST,
               HybridTargetLevel origin_tier = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

template <typename T>
void coordCopy(CoordinateSeries<T> *destination, int frame_dest, const Condensate &origin,
               int index_orig, HybridTargetLevel destination_tier = HybridTargetLevel::HOST,
               HybridTargetLevel origin_tier = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

void coordCopy(PhaseSpaceSynthesis *destination, int index_dest, TrajectoryKind kind,
               CoordinateCycle orientation, const CoordinateFrameReader &origin,
               HybridTargetLevel destination_tier = HybridTargetLevel::HOST,
               HybridTargetLevel origin_tier = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

void coordCopy(PhaseSpaceSynthesis *destination, int index_dest, TrajectoryKind kind,
               const CoordinateFrameReader &origin,
               HybridTargetLevel destination_tier = HybridTargetLevel::HOST,
               HybridTargetLevel origin_tier = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

void coordCopy(PhaseSpaceSynthesis *destination, int index_dest,
               const CoordinateFrameReader &origin,
               HybridTargetLevel destination_tier = HybridTargetLevel::HOST,
               HybridTargetLevel origin_tier = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

void coordCopy(PhaseSpaceSynthesis *destination, int index_dest, TrajectoryKind kind,
               CoordinateCycle orientation, const CoordinateFrame &origin,
               HybridTargetLevel destination_tier = HybridTargetLevel::HOST,
               HybridTargetLevel origin_tier = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

void coordCopy(PhaseSpaceSynthesis *destination, int index_dest, TrajectoryKind kind,
               const CoordinateFrame &origin,
               HybridTargetLevel destination_tier = HybridTargetLevel::HOST,
               HybridTargetLevel origin_tier = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

void coordCopy(PhaseSpaceSynthesis *destination, int index_dest, const CoordinateFrame &origin,
               HybridTargetLevel destination_tier = HybridTargetLevel::HOST,
               HybridTargetLevel origin_tier = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

void coordCopy(PhaseSpaceSynthesis *destination, int index_dest, const PhaseSpaceReader &origin,
               HybridTargetLevel destination_tier = HybridTargetLevel::HOST,
               HybridTargetLevel origin_tier = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);
               
void coordCopy(PhaseSpaceSynthesis *destination, int index_dest, const PhaseSpace &origin,
               HybridTargetLevel destination_tier = HybridTargetLevel::HOST,
               HybridTargetLevel origin_tier = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

template <typename T>
void coordCopy(PsSynthesisWriter *destination, int dest_atom_start, int index_dest,
               const TrajectoryKind kind, const CoordinateSeriesReader<T> &origin, int frame_orig,
               HybridTargetLevel destination_tier = HybridTargetLevel::HOST,
               HybridTargetLevel origin_tier = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

template <typename T>
void coordCopy(PsSynthesisWriter *destination, int dest_atom_start, int index_dest,
               const CoordinateSeriesReader<T> &origin, int frame_orig,
               HybridTargetLevel destination_tier, HybridTargetLevel origin_tier,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

template <typename T>
void coordCopy(PhaseSpaceSynthesis *destination, int index_dest, TrajectoryKind kind,
               CoordinateCycle orientation, const CoordinateSeries<T> &origin, int frame_orig,
               HybridTargetLevel destination_tier = HybridTargetLevel::HOST,
               HybridTargetLevel origin_tier = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

template <typename T>
void coordCopy(PhaseSpaceSynthesis *destination, int index_dest, TrajectoryKind kind,
               const CoordinateSeries<T> &origin, int frame_orig,
               HybridTargetLevel destination_tier = HybridTargetLevel::HOST,
               HybridTargetLevel origin_tier = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

template <typename T>
void coordCopy(PhaseSpaceSynthesis *destination, int index_dest,
               const CoordinateSeries<T> &origin, int frame_orig,
               HybridTargetLevel destination_tier = HybridTargetLevel::HOST,
               HybridTargetLevel origin_tier = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

void coordCopy(PsSynthesisWriter *destination, int dest_atom_start, int index_dest,
               TrajectoryKind kind, const CondensateReader &origin, int orig_atom_start,
               int index_orig, int natom,
               HybridTargetLevel destination_tier = HybridTargetLevel::HOST,
               HybridTargetLevel origin_tier = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

void coordCopy(PsSynthesisWriter *destination, int dest_atom_start, int index_dest,
               const CondensateReader &origin, int orig_atom_start, int index_orig, int natom,
               HybridTargetLevel destination_tier = HybridTargetLevel::HOST,
               HybridTargetLevel origin_tier = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

void coordCopy(PhaseSpaceSynthesis *destination, int index_dest, TrajectoryKind kind,
               CoordinateCycle orientation, const Condensate &origin, int index_orig,
               HybridTargetLevel destination_tier = HybridTargetLevel::HOST,
               HybridTargetLevel origin_tier = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

void coordCopy(PhaseSpaceSynthesis *destination, int index_dest, TrajectoryKind kind,
               const Condensate &origin, int index_orig,
               HybridTargetLevel destination_tier = HybridTargetLevel::HOST,
               HybridTargetLevel origin_tier = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

void coordCopy(PhaseSpaceSynthesis *destination, int index_dest, const Condensate &origin,
               int index_orig, HybridTargetLevel destination_tier = HybridTargetLevel::HOST,
               HybridTargetLevel origin_tier = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

void coordCopy(PsSynthesisWriter *destination, int dest_atom_start, int index_dest,
               const PsSynthesisReader &origin, int orig_atom_start, int index_orig, int natom,
               HybridTargetLevel destination_tier = HybridTargetLevel::HOST,
               HybridTargetLevel origin_tier = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);
  
void coordCopy(PhaseSpaceSynthesis *destination, int index_dest,
               const PhaseSpaceSynthesis &origin, int index_orig,
               HybridTargetLevel destination_tier = HybridTargetLevel::HOST,
               HybridTargetLevel origin_tier = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

void coordCopy(CondensateWriter *destination, int dest_atom_start, int index_dest,
               const CoordinateFrameReader &origin,
               HybridTargetLevel destination_tier = HybridTargetLevel::HOST,
               HybridTargetLevel origin_tier = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

void coordCopy(Condensate *destination, int index_dest, const CoordinateFrame &origin,
               HybridTargetLevel destination_tier = HybridTargetLevel::HOST,
               HybridTargetLevel origin_tier = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

void coordCopy(CondensateWriter *destination, int dest_atom_start, int index_dest,
               const PhaseSpaceReader &origin, TrajectoryKind kind,
               HybridTargetLevel destination_tier = HybridTargetLevel::HOST,
               HybridTargetLevel origin_tier = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

void coordCopy(CondensateWriter *destination, int dest_atom_start, int index_dest,
               const PhaseSpaceReader &origin, HybridTargetLevel destination_tier,
               HybridTargetLevel origin_tier, const GpuDetails &gpu = null_gpu,
               HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

void coordCopy(Condensate *destination, int index_dest, const PhaseSpace &origin,
               TrajectoryKind kind, CoordinateCycle orientation,
               HybridTargetLevel destination_tier = HybridTargetLevel::HOST,
               HybridTargetLevel origin_tier = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

void coordCopy(Condensate *destination, int index_dest, const PhaseSpace &origin,
               TrajectoryKind kind, HybridTargetLevel destination_tier = HybridTargetLevel::HOST,
               HybridTargetLevel origin_tier = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);
  
void coordCopy(Condensate *destination, int index_dest, const PhaseSpace &origin,
               HybridTargetLevel destination_tier = HybridTargetLevel::HOST,
               HybridTargetLevel origin_tier = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);
  
template <typename T>
void coordCopy(CondensateWriter *destination, size_t dest_atom_start, int index_dest,
               const CoordinateSeriesReader<T> &origin, size_t frame_orig,
               HybridTargetLevel destination_tier = HybridTargetLevel::HOST,
               HybridTargetLevel origin_tier = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

template <typename T>
void coordCopy(Condensate *destination, int index_dest, const CoordinateSeries<T> &origin,
               size_t frame_orig, HybridTargetLevel destination_tier = HybridTargetLevel::HOST,
               HybridTargetLevel origin_tier = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);


void coordCopy(CondensateWriter *destination, int dest_atom_start, int index_dest,
               const PsSynthesisReader &origin, int orig_atom_start, int index_orig, int natom,
               TrajectoryKind kind, HybridTargetLevel destination_tier = HybridTargetLevel::HOST,
               HybridTargetLevel origin_tier = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

void coordCopy(CondensateWriter *destination, int dest_atom_start, int index_dest,
               const PsSynthesisReader &origin, int orig_atom_start, int index_orig, int natom,
               HybridTargetLevel destination_tier, HybridTargetLevel origin_tier,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

void coordCopy(Condensate *destination, int index_dest, const PhaseSpaceSynthesis &origin,
               int index_orig, TrajectoryKind kind, CoordinateCycle orientation,
               HybridTargetLevel destination_tier = HybridTargetLevel::HOST,
               HybridTargetLevel origin_tier = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

void coordCopy(Condensate *destination, int index_dest, const PhaseSpaceSynthesis &origin,
               int index_orig, TrajectoryKind kind,
               HybridTargetLevel destination_tier = HybridTargetLevel::HOST,
               HybridTargetLevel origin_tier = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

void coordCopy(Condensate *destination, int index_dest, const PhaseSpaceSynthesis &origin,
               int index_orig, HybridTargetLevel destination_tier = HybridTargetLevel::HOST,
               HybridTargetLevel origin_tier = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

void coordCopy(CondensateWriter *destination, int dest_atom_start, int index_dest,
               const CondensateReader &origin, int orig_atom_start, int index_orig, int natom,
               HybridTargetLevel destination_tier = HybridTargetLevel::HOST,
               HybridTargetLevel origin_tier = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

void coordCopy(Condensate *destination, int index_dest, const Condensate &origin, int index_orig,
               HybridTargetLevel destination_tier = HybridTargetLevel::HOST,
               HybridTargetLevel origin_tier = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);
/// \}

} // namespace trajectory
} // namespace stormm

#include "coordinate_copy.tpp"

#endif
