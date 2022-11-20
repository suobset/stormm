// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace synthesis {

//-------------------------------------------------------------------------------------------------
template <typename T>
void PhaseSpaceSynthesis::import(const T* x_import, const T* y_import, const T* z_import,
                                 const double* box_transform, const double* inverse_transform,
                                 const double* box_dimensions, const int system_index,
                                 const double inverse_scaling_factor, const TrajectoryKind kind,
                                 const CoordinateCycle orientation, const HybridTargetLevel tier) {
  switch (tier) {
  case HybridTargetLevel::HOST:
    {
      llint *x_recv, *y_recv, *z_recv, *x_recv_ovrf, *y_recv_ovrf, *z_recv_ovrf;
      const int pos_start   = atom_starts.readHost(system_index);
      const int pos_end     = pos_start + atom_counts.readHost(system_index);
      const int box_offset  = roundUp(9, warp_size_int) * system_index;
      const int dim_offset  = roundUp(6, warp_size_int) * system_index;
      double conv_factor;
      switch (trajkind) {
      case TrajectoryKind::POSITIONS:
        conv_factor = inverse_scaling_factor * globalpos_scale;
        switch (orientation) {
        case CoordinateCycle::PAST:
          x_recv      = x_prior_coordinates.data();
          y_recv      = y_prior_coordinates.data();
          z_recv      = z_prior_coordinates.data();
          x_recv_ovrf = x_prior_coord_overflow.data();
          y_recv_ovrf = y_prior_coord_overflow.data();
          z_recv_ovrf = z_prior_coord_overflow.data();
          break;
        case CoordinateCycle::PRESENT:
          x_recv      = x_coordinates.data();
          y_recv      = y_coordinates.data();
          z_recv      = z_coordinates.data();
          x_recv_ovrf = x_coordinate_overflow.data();
          y_recv_ovrf = y_coordinate_overflow.data();
          z_recv_ovrf = z_coordinate_overflow.data();
          for (int i = 0; i < 9; i++) {
            box_space_transforms.putHost(box_transform[i], box_offset + i);
            inverse_transforms.putHost(inverse_transform[i], box_offset + i);
            box_dimensions.putHost(inverse_transform[i], dim_offset + i);
            const int95_t fpbv = doubleToint95(inverse_transform[i] * globalpos_scale);
            box_vectors.putHost(fpbv.x, box_offset + i);
            box_vectors_overflow.putHost(fpbv.y, box_offset + i);
          }
          break;
        case CoordinateCycle::FUTURE:
          x_recv      = x_future_coordinates.data();
          y_recv      = y_future_coordinates.data();
          z_recv      = z_future_coordinates.data();
          x_recv_ovrf = x_future_coord_overflow.data();
          y_recv_ovrf = y_future_coord_overflow.data();
          z_recv_ovrf = z_future_coord_overflow.data();
          break;
        }
        break;
      case TrajectoryKind::VELOCITIES:
        conv_factor = inverse_scaling_factor * velocity_scale;
        switch (orientation) {
        case CoordinateCycle::PAST:
          x_recv      = x_prior_velocities.data();
          y_recv      = y_prior_velocities.data();
          z_recv      = z_prior_velocities.data();
          x_recv_ovrf = x_prior_velocity_overflow.data();
          y_recv_ovrf = y_prior_velocity_overflow.data();
          z_recv_ovrf = z_prior_velocity_overflow.data();
          break;
        case CoordinateCycle::PRESENT:
          x_recv      = x_velocities.data();
          y_recv      = y_velocities.data();
          z_recv      = z_velocities.data();
          x_recv_ovrf = x_velocity_overflow.data();
          y_recv_ovrf = y_velocity_overflow.data();
          z_recv_ovrf = z_velocity_overflow.data();
          break;
        case CoordinateCycle::FUTURE:
          x_recv      = x_future_velocities.data();
          y_recv      = y_future_velocities.data();
          z_recv      = z_future_velocities.data();
          x_recv_ovrf = x_future_velocity_overflow.data();
          y_recv_ovrf = y_future_velocity_overflow.data();
          z_recv_ovrf = z_future_velocity_overflow.data();
          break;
        }
        break;
      case TrajectoryKind::FORCES:
        conv_factor = inverse_scaling_factor * force_scale;
        switch (orientation) {
        case CoordinateCycle::PAST:
          x_recv      = x_prior_forces.data();
          y_recv      = y_prior_forces.data();
          z_recv      = z_prior_forces.data();
          x_recv_ovrf = x_prior_force_overflow.data();
          y_recv_ovrf = y_prior_force_overflow.data();
          z_recv_ovrf = z_prior_force_overflow.data();
          break;
        case CoordinateCycle::PRESENT:
          x_recv      = x_forces.data();
          y_recv      = y_forces.data();
          z_recv      = z_forces.data();
          x_recv_ovrf = x_force_overflow.data();
          y_recv_ovrf = y_force_overflow.data();
          z_recv_ovrf = z_force_overflow.data();
          break;
        case CoordinateCycle::FUTURE:
          x_recv      = x_future_forces.data();
          y_recv      = y_future_forces.data();
          z_recv      = z_future_forces.data();
          x_recv_ovrf = x_future_force_overflow.data();
          y_recv_ovrf = y_future_force_overflow.data();
          z_recv_ovrf = z_future_force_overflow.data();
          break;
        }
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
void PhaseSpaceSynthesis::import(const CoordinateSeries<T> &cs, const int frame_index,
                                 const int system_index, const TrajectoryKind kind,
                                 const CoordinateCycle orientation, const HybridTargetLevel tier) {

}
  
} // namespace synthesis
} // namespace stormm
