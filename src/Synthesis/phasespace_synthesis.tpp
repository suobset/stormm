// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace synthesis {

//-------------------------------------------------------------------------------------------------
template <typename T>
void PhaseSpaceSynthesis::import(const T* x_import, const T* y_import, const T* z_import,
                                 const double* box_xform_in, const double* inverse_xform_in,
                                 const double* box_dimensions_in, const int system_index,
                                 const double inverse_scaling_factor, const TrajectoryKind kind,
                                 const CoordinateCycle orientation, const HybridTargetLevel tier) {
  llint *x_recv, *y_recv, *z_recv;
  int *x_recv_ovrf, *y_recv_ovrf, *z_recv_ovrf;
  double *box_xform_ptr, *inverse_xform_ptr, *box_dim_ptr;
  double conv_factor;
  switch (kind) {
  case TrajectoryKind::POSITIONS:
    conv_factor = inverse_scaling_factor * globalpos_scale;
    switch (orientation) {
    case CoordinateCycle::PAST:
      x_recv      = x_prior_coordinates.data(tier);
      y_recv      = y_prior_coordinates.data(tier);
      z_recv      = z_prior_coordinates.data(tier);
      x_recv_ovrf = x_prior_coord_overflow.data(tier);
      y_recv_ovrf = y_prior_coord_overflow.data(tier);
      z_recv_ovrf = z_prior_coord_overflow.data(tier);
      break;
    case CoordinateCycle::PRESENT:
      x_recv      = x_coordinates.data(tier);
      y_recv      = y_coordinates.data(tier);
      z_recv      = z_coordinates.data(tier);
      x_recv_ovrf = x_coordinate_overflow.data(tier);
      y_recv_ovrf = y_coordinate_overflow.data(tier);
      z_recv_ovrf = z_coordinate_overflow.data(tier);
      box_xform_ptr = box_space_transforms.data(tier);
      inverse_xform_ptr = inverse_transforms.data(tier);
      box_dim_ptr = box_dimensions.data(tier);
      break;
    case CoordinateCycle::FUTURE:
      x_recv      = x_future_coordinates.data(tier);
      y_recv      = y_future_coordinates.data(tier);
      z_recv      = z_future_coordinates.data(tier);
      x_recv_ovrf = x_future_coord_overflow.data(tier);
      y_recv_ovrf = y_future_coord_overflow.data(tier);
      z_recv_ovrf = z_future_coord_overflow.data(tier);
      break;
    }
    break;
  case TrajectoryKind::VELOCITIES:
    conv_factor = inverse_scaling_factor * velocity_scale;
    switch (orientation) {
    case CoordinateCycle::PAST:
      x_recv      = x_prior_velocities.data(tier);
      y_recv      = y_prior_velocities.data(tier);
      z_recv      = z_prior_velocities.data(tier);
      x_recv_ovrf = x_prior_velocity_overflow.data(tier);
      y_recv_ovrf = y_prior_velocity_overflow.data(tier);
      z_recv_ovrf = z_prior_velocity_overflow.data(tier);
      break;
    case CoordinateCycle::PRESENT:
      x_recv      = x_velocities.data(tier);
      y_recv      = y_velocities.data(tier);
      z_recv      = z_velocities.data(tier);
      x_recv_ovrf = x_velocity_overflow.data(tier);
      y_recv_ovrf = y_velocity_overflow.data(tier);
      z_recv_ovrf = z_velocity_overflow.data(tier);
      break;
    case CoordinateCycle::FUTURE:
      x_recv      = x_future_velocities.data(tier);
      y_recv      = y_future_velocities.data(tier);
      z_recv      = z_future_velocities.data(tier);
      x_recv_ovrf = x_future_velocity_overflow.data(tier);
      y_recv_ovrf = y_future_velocity_overflow.data(tier);
      z_recv_ovrf = z_future_velocity_overflow.data(tier);
      break;
    }
    break;
  case TrajectoryKind::FORCES:
    conv_factor = inverse_scaling_factor * force_scale;
    switch (orientation) {
    case CoordinateCycle::PAST:
      x_recv      = x_prior_forces.data(tier);
      y_recv      = y_prior_forces.data(tier);
      z_recv      = z_prior_forces.data(tier);
      x_recv_ovrf = x_prior_force_overflow.data(tier);
      y_recv_ovrf = y_prior_force_overflow.data(tier);
      z_recv_ovrf = z_prior_force_overflow.data(tier);
      break;
    case CoordinateCycle::PRESENT:
      x_recv      = x_forces.data(tier);
      y_recv      = y_forces.data(tier);
      z_recv      = z_forces.data(tier);
      x_recv_ovrf = x_force_overflow.data(tier);
      y_recv_ovrf = y_force_overflow.data(tier);
      z_recv_ovrf = z_force_overflow.data(tier);
      break;
    case CoordinateCycle::FUTURE:
      x_recv      = x_future_forces.data(tier);
      y_recv      = y_future_forces.data(tier);
      z_recv      = z_future_forces.data(tier);
      x_recv_ovrf = x_future_force_overflow.data(tier);
      y_recv_ovrf = y_future_force_overflow.data(tier);
      z_recv_ovrf = z_future_force_overflow.data(tier);
      break;
    }
    break;
  }
  switch (tier) {
  case HybridTargetLevel::HOST:
    {
      const int pos_start   = atom_starts.readHost(system_index);
      const int pos_end     = pos_start + atom_counts.readHost(system_index);
      const int box_offset  = roundUp(9, warp_size_int) * system_index;
      const int dim_offset  = roundUp(6, warp_size_int) * system_index;
      switch (kind) {
      case TrajectoryKind::POSITIONS:
        for (size_t i = 0; i < 9LLU; i++) {
          box_xform_ptr[box_offset + i] = box_xform_in[i];
          inverse_xform_ptr[box_offset + i] = inverse_xform_in[i];
          const int95_t fpbv = doubleToInt95(inverse_xform_in[i] * globalpos_scale);
          box_vectors.putHost(fpbv.x, box_offset + i);
          box_vector_overflow.putHost(fpbv.y, box_offset + i);
        }
        for (size_t i = 0; i < 6LLU; i++) {
          box_dim_ptr[dim_offset + i] = box_dimensions_in[i];
        }
        break;
      case TrajectoryKind::VELOCITIES:
      case TrajectoryKind::FORCES:
        break;
      }
      for (int i = pos_start; i < pos_end; i++) {
        const size_t ip = i - pos_start;
        const int95_t fpx = doubleToInt95(static_cast<double>(x_import[ip]) * conv_factor);
        const int95_t fpy = doubleToInt95(static_cast<double>(y_import[ip]) * conv_factor);
        const int95_t fpz = doubleToInt95(static_cast<double>(z_import[ip]) * conv_factor);
        x_recv[i]      = fpx.x;
        x_recv_ovrf[i] = fpx.y;
        y_recv[i]      = fpy.x;
        y_recv_ovrf[i] = fpy.y;
        z_recv[i]      = fpz.x;
        z_recv_ovrf[i] = fpz.y;
      }
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    
    break;
#endif
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void PhaseSpaceSynthesis::import(const CoordinateSeriesReader<T> &csr, const int frame_index,
                                 const int system_index, const TrajectoryKind kind,
                                 const CoordinateCycle orientation, const HybridTargetLevel tier) {
  const size_t xfrm_offset = frame_index * roundUp(9, warp_size_int);
  const size_t bdim_offset = frame_index * roundUp(6, warp_size_int);
  const size_t atom_offset = static_cast<size_t>(frame_index) *
                             static_cast<size_t>(roundUp(csr.natom, warp_size_int));
  import(&csr.xcrd[atom_offset], &csr.ycrd[atom_offset], &csr.zcrd[atom_offset],
         &csr.umat[xfrm_offset], &csr.invu[xfrm_offset], &csr.boxdim[bdim_offset], system_index,
         csr.inv_gpos_scale, kind, orientation, tier);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void PhaseSpaceSynthesis::import(const CoordinateSeriesWriter<T> &csw, const int frame_index,
                                 const int system_index, const TrajectoryKind kind,
                                 const CoordinateCycle orientation, const HybridTargetLevel tier) {
  import(csw, frame_index, system_index, kind, orientation, tier);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void PhaseSpaceSynthesis::import(const CoordinateSeries<T> &cs, const int frame_index,
                                 const int system_index, const TrajectoryKind kind,
                                 const CoordinateCycle orientation, const HybridTargetLevel tier) {
  import(cs.data(), frame_index, system_index, kind, orientation, tier);
}
  
} // namespace synthesis
} // namespace stormm
