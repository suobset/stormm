#include "Constants/fixed_precision.h"
#include "Math/rounding.h"
#include "Math/matrix_ops.h"
#include "phasespace_synthesis.h"

namespace omni {
namespace trajectory {

using cuda::HybridKind;
using math::roundUp;
using math::invertSquareMatrix;
using numerics::global_position_scale_lf;
using numerics::global_position_scale_bits;
using numerics::global_force_scale_lf;
using topology::UnitCellType;

//-------------------------------------------------------------------------------------------------
PsSynthesisReader::PsSynthesisReader(const int system_count_in, const UnitCellType unit_cell_in,
                                     const ThermostatKind heat_bath_kind_in,
                                     const BarostatKind piston_kind_in, const double time_step_in,
                                     const int* atom_starts_in, const int* atom_counts_in,
                                     const llint* boxvecs_in, const double* umat_in,
                                     const double* invu_in, const double* boxdims_in,
                                     const float* sp_umat_in, const float* sp_invu_in,
                                     const float* sp_boxdims_in, const longlong4* xyz_qlj_in,
                                     const double* xvel_in, const double* yvel_in,
                                     const double* zvel_in, const llint* xfrc_in,
                                     const llint* yfrc_in, const llint* zfrc_in) :
    system_count{system_count_in}, unit_cell{unit_cell_in}, heat_bath_kind{heat_bath_kind_in},
    piston_kind{piston_kind_in}, time_step{time_step_in}, atom_starts{atom_starts_in},
    atom_counts{atom_counts_in}, boxvecs{boxvecs_in}, umat{umat_in}, invu{invu_in},
    boxdims{boxdims_in}, sp_umat{sp_umat_in}, sp_invu{sp_invu_in}, sp_boxdims{sp_boxdims_in},
    xyz_qlj{xyz_qlj_in}, xvel{xvel_in}, yvel{yvel_in}, zvel{zvel_in}, xfrc{xfrc_in}, yfrc{yfrc_in},
    zfrc{zfrc_in}
{}

//-------------------------------------------------------------------------------------------------
PsSynthesisWriter::PsSynthesisWriter(const int system_count_in, const UnitCellType unit_cell_in,
                                     const ThermostatKind heat_bath_kind_in,
                                     const BarostatKind piston_kind_in, const double time_step_in,
                                     const int* atom_starts_in, const int* atom_counts_in,
                                     const llint* boxvecs_in, const double* umat_in,
                                     const double* invu_in, const double* boxdims_in,
                                     const float* sp_umat_in, const float* sp_invu_in,
                                     const float* sp_boxdims_in, longlong4* xyz_qlj_in,
                                     double* xvel_in, double* yvel_in, double* zvel_in,
                                     llint* xfrc_in, llint* yfrc_in, llint* zfrc_in) :
    system_count{system_count_in}, unit_cell{unit_cell_in}, heat_bath_kind{heat_bath_kind_in},
    piston_kind{piston_kind_in}, time_step{time_step_in}, atom_starts{atom_starts_in},
    atom_counts{atom_counts_in}, boxvecs{boxvecs_in}, umat{umat_in}, invu{invu_in},
    boxdims{boxdims_in}, sp_umat{sp_umat_in}, sp_invu{sp_invu_in}, sp_boxdims{sp_boxdims_in},
    xyz_qlj{xyz_qlj_in}, xvel{xvel_in}, yvel{yvel_in}, zvel{zvel_in}, xfrc{xfrc_in}, yfrc{yfrc_in},
    zfrc{zfrc_in}
{}

//-------------------------------------------------------------------------------------------------
PhaseSpaceSynthesis::PhaseSpaceSynthesis(const std::vector<PhaseSpace> &ps_list,
                                         const double time_step_in,
                                         const std::vector<AtomGraph*> &ag_list,
                                         const std::vector<Thermostat> &heat_baths_in,
                                         const std::vector<Barostat> &pistons_in) :
    system_count{static_cast<int>(ps_list.size())},
    unit_cell{UnitCellType::NONE},
    heat_bath_kind{ThermostatKind::NONE},
    piston_kind{BarostatKind::NONE},
    time_step{time_step_in},
    atom_starts{ps_list.size(), "labframe_starts"},
    atom_counts{ps_list.size(), "labframe_counts"},
    xyz_qlj{Hybrid<longlong4>(HybridKind::ARRAY, "labframe_xyz_idqlj")},
    x_velocities{Hybrid<double>(HybridKind::POINTER, "labframe_vx")},
    y_velocities{Hybrid<double>(HybridKind::POINTER, "labframe_vy")},
    z_velocities{Hybrid<double>(HybridKind::POINTER, "labframe_vz")},
    x_forces{Hybrid<llint>(HybridKind::POINTER, "labframe_fx")},
    y_forces{Hybrid<llint>(HybridKind::POINTER, "labframe_fy")},
    z_forces{Hybrid<llint>(HybridKind::POINTER, "labframe_fz")},
    box_vectors{Hybrid<llint>(HybridKind::POINTER, "labframe_boxdims")},
    box_space_transforms{Hybrid<double>(HybridKind::POINTER, "labframe_umat")},
    inverse_transforms{Hybrid<double>(HybridKind::POINTER, "labframe_invu")},
    box_dimensions{Hybrid<double>(HybridKind::POINTER, "labframe_dims")},
    sp_box_space_transforms{Hybrid<float>(HybridKind::POINTER, "labframe_umat")},
    sp_inverse_transforms{Hybrid<float>(HybridKind::POINTER, "labframe_invu")},
    sp_box_dimensions{Hybrid<float>(HybridKind::POINTER, "labframe_dims")},
    llint_data{Hybrid<llint>(HybridKind::ARRAY, "labframe_llint")},
    double_data{Hybrid<double>(HybridKind::ARRAY, "labframe_double")},
    float_data{Hybrid<float>(HybridKind::ARRAY, "labframe_float")},
    topologies{std::vector<AtomGraph*>(ag_list.begin(), ag_list.end())}
{
  // Check validity of input arrays
  if (system_count == 0) {
    rtErr("At least one PhaseSpace object must be provided.", "PhaseSpaceSynthesis");
  }
  else if (ag_list.size() != system_count) {
    rtErr("One topology pointer must be provided for each coordinate system (currently " +
          std::to_string(ag_list.size()) + " topology pointers and " +
          std::to_string(system_count) + " PhaseSpace objects.", "PhaseSpaceSynthesis");
  }
  else if (heat_baths_in.size() != system_count) {
    rtErr("One thermostat must be provided for each system (currently " +
          std::to_string(system_count) + " systems and " + std::to_string(heat_baths_in.size()) +
          " thermostats.", "PhaseSpaceSynthesis");
  }
  else if (pistons_in.size() != system_count) {
    rtErr("One barostat must be provided for each system (currently " +
          std::to_string(system_count) + " systems and " + std::to_string(pistons_in.size()) +
          " barostats.", "PhaseSpaceSynthesis");
  }

  // Check validity of thermostats: all thermostats must be of the same type, but temperatures
  // and details of the coupling strength may vary between systems
  heat_bath_kind = heat_baths_in[0].kind;
  for (int i = 1; i < system_count; i++) {
    if (heat_baths_in[i].kind != heat_baths_in[0].kind) {
      rtErr("All system thermostats must be of the same class.  System " + std::to_string(i) +
            " uses a " + getThermostatName(heat_baths_in[i].kind) + " thermostat whereas the "
            "first system locks the method to " + getThermostatName(heat_baths_in[0].kind) + ".",
            "PhaseSpaceSynthesis");
    }
  }

  // Check validity of barostats: all barostats must be of the same type, but pressures and
  // details of the rescaling may vary between systems
  piston_kind = pistons_in[0].kind;
  for (int i = 1; i < system_count; i++) {
    if (pistons_in[i].kind != pistons_in[0].kind) {
      rtErr("All system barostats must be of the same class.  System " + std::to_string(i) +
            " uses a " + getBarostatName(pistons_in[i].kind) + " barostat whereas the first "
            "system locks the method to " + getBarostatName(pistons_in[0].kind) + ".",
            "PhaseSpaceSynthesis");
    }
  }

  // Check that coordinates match topologies.  Set atom starts and counts in the process.
  int acc_limit = 0;
  for (int i = 0; i < system_count; i++) {
    const int natom = ps_list[i].getAtomCount();
    atom_counts.putHost(natom, i);
    if (natom != ag_list[i]->getAtomCount()) {
      rtErr("Input topology and coordinate sets disagree on atom counts (" +
            std::to_string(ag_list[i]->getAtomCount()) + " vs. " + std::to_string(natom) + ").",
            "PhaseSpaceSynthesis");
    }
    atom_starts.putHost(acc_limit, i);
    acc_limit += roundUp(natom, warp_size_int);
  }

  // Allocate data and set internal pointers
  const int atom_stride = acc_limit;
  const int mtrx_stride = roundUp(9, warp_size_int);
  const int xfrm_stride = system_count * mtrx_stride;
  xyz_qlj.resize(atom_stride);
  llint_data.resize((3 * atom_stride) + xfrm_stride);
  x_forces.setPointer(&llint_data,                  0, atom_stride);
  y_forces.setPointer(&llint_data,        atom_stride, atom_stride);
  z_forces.setPointer(&llint_data,    2 * atom_stride, atom_stride);
  box_vectors.setPointer(&llint_data, 3 * atom_stride, xfrm_stride);
  double_data.resize(3 * xfrm_stride);
  box_space_transforms.setPointer(&double_data,               0, xfrm_stride);
  inverse_transforms.setPointer(&double_data,       xfrm_stride, xfrm_stride);
  box_dimensions.setPointer(&double_data,       2 * xfrm_stride, xfrm_stride);
  float_data.resize(3 * xfrm_stride);
  sp_box_space_transforms.setPointer(&float_data,               0, xfrm_stride);
  sp_inverse_transforms.setPointer(&float_data,       xfrm_stride, xfrm_stride);
  sp_box_dimensions.setPointer(&float_data,       2 * xfrm_stride, xfrm_stride);

  // Loop over all systems
  for (int i = 0; i < system_count; i++) {

    // Get a reader for the PhaseSpace object's host-side data
    PhaseSpaceReader psr = ps_list[i].data();

    // Assign coordinates, converting from double precision to fixed precision format
    std::vector<int> lj_indices = ag_list[i]->getLennardJonesIndex();
    std::vector<int> charge_indices = ag_list[i]->getChargeIndex();
    double* xvel_ptr = x_velocities.data();
    double* yvel_ptr = y_velocities.data();
    double* zvel_ptr = z_velocities.data();
    llint* xfrc_ptr = x_forces.data();
    llint* yfrc_ptr = y_forces.data();
    llint* zfrc_ptr = z_forces.data();
    const int asi = atom_starts.readHost(i);
    for (int j = 0; j < psr.natom; j++) {
      llint4 atompt;

      // The x, y, and z components hold x, y, and z positions in a "lab frame" unwrapped form.
      atompt.x = static_cast<llint>(psr.xcrd[j] * global_position_scale_lf);
      atompt.y = static_cast<llint>(psr.ycrd[j] * global_position_scale_lf);
      atompt.z = static_cast<llint>(psr.zcrd[j] * global_position_scale_lf);
      const llint qi = charge_indices[j];
      const llint qlj = lj_indices[j];

      // The w component holds a bit-pack representation of the charge index (up to 8,388,608
      // unique charge types are possible), the Lennard-Jones atom type index (up to 512 unique
      // Lennard-Jones types are possible) and the atom index (up to 2,147,483,648 atoms can be
      // indexed, but the memory will break for other reasons long before that limit is reached).
      atompt.w = ((qi << 41) | (qlj << 32) | (asi + j));
      const size_t asi_j = asi + j;
      xyz_qlj.putHost(atompt, asi_j);
      xvel_ptr[asi_j] = psr.xvel[j];
      yvel_ptr[asi_j] = psr.yvel[j];
      zvel_ptr[asi_j] = psr.zvel[j];
      xfrc_ptr[asi_j] = static_cast<llint>(psr.xfrc[j] * global_force_scale_lf);
      yfrc_ptr[asi_j] = static_cast<llint>(psr.yfrc[j] * global_force_scale_lf);
      zfrc_ptr[asi_j] = static_cast<llint>(psr.zfrc[j] * global_force_scale_lf);
    }

    // Handle the box space transformation.  The transformation matrices of the labframe will
    // become slightly desynchronized from the original PhaseSpace object, but this is in the
    // interest of making a system which can be represented in fixed precision, pixelated
    // coordinates with box vectors which are likewise multiples of the positional discretization.
    for (int j = 0; j < 9; j++) {
      const llint lli_invu = psr.invu[j] * global_position_scale_lf;
      const double d_invu = static_cast<double>(lli_invu >> global_position_scale_bits);
      inverse_transforms.putHost(d_invu, i * mtrx_stride + j);
      sp_inverse_transforms.putHost(d_invu, i * mtrx_stride + j);
    }
    const double* invu_ptr = &inverse_transforms.data()[i * mtrx_stride];
    double* umat_ptr = &box_space_transforms.data()[i * mtrx_stride];
    invertSquareMatrix(invu_ptr, umat_ptr, 3);
    for (int j = 0; j < 6; j++) {
      box_dimensions.putHost(psr.boxdim[j], i * roundUp(6, warp_size_int) + j);
      sp_box_dimensions.putHost(psr.boxdim[j], i * roundUp(6, warp_size_int) + j);
    }
  }
  for (int i = 0; i < system_count * mtrx_stride; i++) {
    sp_box_space_transforms.putHost(box_space_transforms.readHost(i), i);
  }
}

#ifdef OMNI_USE_HPC
//-------------------------------------------------------------------------------------------------
void PhaseSpaceSynthesis::upload() {
  atom_starts.upload();
  atom_counts.upload();
  xyz_qlj.upload();
  llint_data.upload();
  double_data.upload();
  float_data.upload();
}

//-------------------------------------------------------------------------------------------------
void PhaseSpaceSynthesis::download() {
  atom_starts.download();
  atom_counts.download();
  xyz_qlj.download();
  llint_data.download();
  double_data.download();
  float_data.download();
}
#endif

//-------------------------------------------------------------------------------------------------
void PhaseSpaceSynthesis::velocityVerletCoordinateUpdate() {

}

//-------------------------------------------------------------------------------------------------
void PhaseSpaceSynthesis::velocityVerletVelocityUpdate() {

}

//-------------------------------------------------------------------------------------------------
void PhaseSpaceSynthesis::assignMaxwellVelocities() {

}

//-------------------------------------------------------------------------------------------------
void PhaseSpaceSynthesis::berendsenThermocoupling() {

}

//-------------------------------------------------------------------------------------------------
const PsSynthesisReader PhaseSpaceSynthesis::data(HybridTargetLevel tier) const {
  return PsSynthesisReader(system_count, unit_cell, heat_bath_kind, piston_kind, time_step,
                           atom_starts.data(tier), atom_counts.data(tier), box_vectors.data(tier),
                           box_space_transforms.data(tier), inverse_transforms.data(tier),
                           box_dimensions.data(tier), sp_box_space_transforms.data(tier),
                           sp_inverse_transforms.data(tier), sp_box_dimensions.data(tier),
                           xyz_qlj.data(tier), x_velocities.data(tier), y_velocities.data(tier),
                           z_velocities.data(tier), x_forces.data(tier), y_forces.data(tier),
                           z_forces.data(tier));
}

//-------------------------------------------------------------------------------------------------
PsSynthesisWriter PhaseSpaceSynthesis::data(HybridTargetLevel tier) {
  return PsSynthesisWriter(system_count, unit_cell, heat_bath_kind, piston_kind, time_step,
                           atom_starts.data(tier), atom_counts.data(tier), box_vectors.data(tier),
                           box_space_transforms.data(tier), inverse_transforms.data(tier),
                           box_dimensions.data(tier), sp_box_space_transforms.data(tier),
                           sp_inverse_transforms.data(tier), sp_box_dimensions.data(tier),
                           xyz_qlj.data(tier), x_velocities.data(tier), y_velocities.data(tier),
                           z_velocities.data(tier), x_forces.data(tier), y_forces.data(tier),
                           z_forces.data(tier));
}

} // namespace trajectory
} // namespace omni
