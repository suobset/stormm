#include <cmath>
#ifdef STORMM_USE_HPC
#  ifdef STORMM_USE_CUDA
#  include <cuda_runtime.h>
#  endif
#endif
#include "copyright.h"
#include "Constants/hpc_bounds.h"
#include "Constants/fixed_precision.h"
#include "FileManagement/file_listing.h"
#include "Math/rounding.h"
#include "Math/matrix_ops.h"
#include "Math/vector_ops.h"
#include "Trajectory/write_frame.h"
#include "phasespace_synthesis.h"
#ifdef STORMM_USE_HPC
#  include "hpc_phasespace_synthesis.h"
#endif

namespace stormm {
namespace synthesis {

using card::HybridKind;
using diskutil::splitPath;
using diskutil::DataFormat;
using diskutil::DrivePathType;
using diskutil::getDrivePathType;
using math::roundUp;
using math::invertSquareMatrix;
using math::tileVector;
using numerics::checkGlobalPositionBits;
using numerics::checkLocalPositionBits;
using numerics::checkVelocityBits;
using numerics::checkForceBits;
using numerics::force_scale_nonoverflow_bits;
using numerics::globalpos_scale_nonoverflow_bits;
using numerics::max_llint_accumulation;
using numerics::splitRealConversion;
using numerics::transform_scale_lf;
using numerics::velocity_scale_nonoverflow_bits;
using topology::UnitCellType;
using trajectory::PhaseSpaceWriter;
using trajectory::PhaseSpaceReader;
using trajectory::writeFrame;
  
//-------------------------------------------------------------------------------------------------
PsSynthesisWriter::PsSynthesisWriter(const int system_count_in, const UnitCellType unit_cell_in,
                                     const ThermostatKind heat_bath_kind_in,
                                     const BarostatKind piston_kind_in, const double time_step_in,
                                     const int* atom_starts_in, const int* atom_counts_in,
                                     const double gpos_scale_in, const double lpos_scale_in,
                                     const double vel_scale_in, const double frc_scale_in,
                                     const int gpos_bits_in, const int lpos_bits_in,
                                     const int vel_bits_in, const int frc_bits_in,
                                     llint* boxvecs_in, double* umat_in, double* invu_in,
                                     double* boxdims_in, float* sp_umat_in, float* sp_invu_in,
                                     float* sp_boxdims_in, llint* xcrd_in, llint* ycrd_in,
                                     llint* zcrd_in, int* xcrd_ovrf_in, int* ycrd_ovrf_in,
                                     int* zcrd_ovrf_in, llint* xvel_in, llint* yvel_in,
                                     llint* zvel_in, int* xvel_ovrf_in, int* yvel_ovrf_in,
                                     int* zvel_ovrf_in, llint* xfrc_in, llint* yfrc_in,
                                     llint* zfrc_in, int* xfrc_ovrf_in, int* yfrc_ovrf_in,
                                     int* zfrc_ovrf_in, llint* xprv_in, llint* yprv_in,
                                     llint* zprv_in, int* xprv_ovrf_in, int* yprv_ovrf_in,
                                     int* zprv_ovrf_in, llint* xnxt_in, llint* ynxt_in,
                                     llint* znxt_in, int* xnxt_ovrf_in, int* ynxt_ovrf_in,
                                     int* znxt_ovrf_in) :
    system_count{system_count_in}, unit_cell{unit_cell_in}, heat_bath_kind{heat_bath_kind_in},
    piston_kind{piston_kind_in}, time_step{time_step_in}, atom_starts{atom_starts_in},
    atom_counts{atom_counts_in}, gpos_scale{gpos_scale_in}, lpos_scale{lpos_scale_in},
    vel_scale{vel_scale_in}, frc_scale{frc_scale_in},
    inv_gpos_scale{1.0 / gpos_scale_in},
    inv_lpos_scale{1.0 / lpos_scale_in},
    inv_vel_scale{1.0 / vel_scale_in},
    inv_frc_scale{1.0 / frc_scale_in},
    gpos_scale_f{static_cast<float>(gpos_scale)},
    lpos_scale_f{static_cast<float>(lpos_scale)},
    vel_scale_f{static_cast<float>(vel_scale)},
    frc_scale_f{static_cast<float>(frc_scale)},
    inv_gpos_scale_f{static_cast<float>(inv_gpos_scale)},
    inv_lpos_scale_f{static_cast<float>(inv_lpos_scale)},
    inv_vel_scale_f{static_cast<float>(inv_vel_scale)},
    inv_frc_scale_f{static_cast<float>(inv_frc_scale)},
    gpos_bits{gpos_bits_in}, lpos_bits{lpos_bits_in}, vel_bits{vel_bits_in}, frc_bits{frc_bits_in},
    boxvecs{boxvecs_in}, umat{umat_in}, invu{invu_in}, boxdims{boxdims_in}, sp_umat{sp_umat_in},
    sp_invu{sp_invu_in}, sp_boxdims{sp_boxdims_in}, xcrd{xcrd_in}, ycrd{ycrd_in}, zcrd{zcrd_in},
    xcrd_ovrf{xcrd_ovrf_in}, ycrd_ovrf{ycrd_ovrf_in}, zcrd_ovrf{zcrd_ovrf_in}, 
    xvel{xvel_in}, yvel{yvel_in}, zvel{zvel_in}, xvel_ovrf{xvel_ovrf_in}, yvel_ovrf{yvel_ovrf_in},
    zvel_ovrf{zvel_ovrf_in}, xfrc{xfrc_in}, yfrc{yfrc_in}, zfrc{zfrc_in}, xfrc_ovrf{xfrc_ovrf_in},
    yfrc_ovrf{yfrc_ovrf_in}, zfrc_ovrf{zfrc_ovrf_in}, xprv{xprv_in}, yprv{yprv_in}, zprv{zprv_in},
    xprv_ovrf{xprv_ovrf_in}, yprv_ovrf{yprv_ovrf_in}, zprv_ovrf{zprv_ovrf_in}, xnxt{xnxt_in},
    ynxt{ynxt_in}, znxt{znxt_in}, xnxt_ovrf{xnxt_ovrf_in}, ynxt_ovrf{ynxt_ovrf_in},
    znxt_ovrf{znxt_ovrf_in}
{}

//-------------------------------------------------------------------------------------------------
PsSynthesisReader::PsSynthesisReader(const int system_count_in, const UnitCellType unit_cell_in,
                                     const ThermostatKind heat_bath_kind_in,
                                     const BarostatKind piston_kind_in, const double time_step_in,
                                     const int* atom_starts_in, const int* atom_counts_in,
                                     const double gpos_scale_in, const double lpos_scale_in,
                                     const double vel_scale_in, const double frc_scale_in,
                                     const int gpos_bits_in, const int lpos_bits_in,
                                     const int vel_bits_in, const int frc_bits_in,
                                     const llint* boxvecs_in, const double* umat_in,
                                     const double* invu_in, const double* boxdims_in,
                                     const float* sp_umat_in, const float* sp_invu_in,
                                     const float* sp_boxdims_in, const llint* xcrd_in,
                                     const llint* ycrd_in, const llint* zcrd_in,
                                     const int* xcrd_ovrf_in, const int* ycrd_ovrf_in,
                                     const int* zcrd_ovrf_in, const llint* xvel_in,
                                     const llint* yvel_in, const llint* zvel_in,
                                     const int* xvel_ovrf_in, const int* yvel_ovrf_in,
                                     const int* zvel_ovrf_in, const llint* xfrc_in,
                                     const llint* yfrc_in, const llint* zfrc_in,
                                     const int* xfrc_ovrf_in, const int* yfrc_ovrf_in,
                                     const int* zfrc_ovrf_in, const llint* xprv_in,
                                     const llint* yprv_in, const llint* zprv_in,
                                     const int* xprv_ovrf_in, const int* yprv_ovrf_in,
                                     const int* zprv_ovrf_in, const llint* xnxt_in,
                                     const llint* ynxt_in, const llint* znxt_in,
                                     const int* xnxt_ovrf_in, const int* ynxt_ovrf_in,
                                     const int* znxt_ovrf_in) :
    system_count{system_count_in}, unit_cell{unit_cell_in}, heat_bath_kind{heat_bath_kind_in},
    piston_kind{piston_kind_in}, time_step{time_step_in}, atom_starts{atom_starts_in},
    atom_counts{atom_counts_in}, gpos_scale{gpos_scale_in}, lpos_scale{lpos_scale_in},
    vel_scale{vel_scale_in}, frc_scale{frc_scale_in},
    inv_gpos_scale{1.0 / gpos_scale_in},
    inv_lpos_scale{1.0 / lpos_scale_in},
    inv_vel_scale{1.0 / vel_scale_in},
    inv_frc_scale{1.0 / frc_scale_in},
    gpos_scale_f{static_cast<float>(gpos_scale)},
    lpos_scale_f{static_cast<float>(lpos_scale)},
    vel_scale_f{static_cast<float>(vel_scale)},
    frc_scale_f{static_cast<float>(frc_scale)},
    inv_gpos_scale_f{static_cast<float>(inv_gpos_scale)},
    inv_lpos_scale_f{static_cast<float>(inv_lpos_scale)},
    inv_vel_scale_f{static_cast<float>(inv_vel_scale)},
    inv_frc_scale_f{static_cast<float>(inv_frc_scale)},
    gpos_bits{gpos_bits_in}, lpos_bits{lpos_bits_in}, vel_bits{vel_bits_in}, frc_bits{frc_bits_in},
    boxvecs{boxvecs_in}, umat{umat_in}, invu{invu_in}, boxdims{boxdims_in}, sp_umat{sp_umat_in},
    sp_invu{sp_invu_in}, sp_boxdims{sp_boxdims_in}, xcrd{xcrd_in}, ycrd{ycrd_in}, zcrd{zcrd_in},
    xcrd_ovrf{xcrd_ovrf_in}, ycrd_ovrf{ycrd_ovrf_in}, zcrd_ovrf{zcrd_ovrf_in}, 
    xvel{xvel_in}, yvel{yvel_in}, zvel{zvel_in}, xvel_ovrf{xvel_ovrf_in}, yvel_ovrf{yvel_ovrf_in},
    zvel_ovrf{zvel_ovrf_in}, xfrc{xfrc_in}, yfrc{yfrc_in}, zfrc{zfrc_in}, xfrc_ovrf{xfrc_ovrf_in},
    yfrc_ovrf{yfrc_ovrf_in}, zfrc_ovrf{zfrc_ovrf_in}, xprv{xprv_in}, yprv{yprv_in}, zprv{zprv_in},
    xprv_ovrf{xprv_ovrf_in}, yprv_ovrf{yprv_ovrf_in}, zprv_ovrf{zprv_ovrf_in}, xnxt{xnxt_in},
    ynxt{ynxt_in}, znxt{znxt_in}, xnxt_ovrf{xnxt_ovrf_in}, ynxt_ovrf{ynxt_ovrf_in},
    znxt_ovrf{znxt_ovrf_in}
{}

//-------------------------------------------------------------------------------------------------
PsSynthesisReader::PsSynthesisReader(const PsSynthesisWriter &psyw) :
    system_count{psyw.system_count},
    unit_cell{psyw.unit_cell},
    heat_bath_kind{psyw.heat_bath_kind},
    piston_kind{psyw.piston_kind},
    time_step{psyw.time_step},
    atom_starts{psyw.atom_starts},
    atom_counts{psyw.atom_counts},
    gpos_scale{psyw.gpos_scale},
    lpos_scale{psyw.lpos_scale},
    vel_scale{psyw.vel_scale},
    frc_scale{psyw.frc_scale},
    inv_gpos_scale{1.0 / gpos_scale},
    inv_lpos_scale{1.0 / lpos_scale},
    inv_vel_scale{1.0 / vel_scale},
    inv_frc_scale{1.0 / frc_scale},
    gpos_scale_f{static_cast<float>(gpos_scale)},
    lpos_scale_f{static_cast<float>(lpos_scale)},
    vel_scale_f{static_cast<float>(vel_scale)},
    frc_scale_f{static_cast<float>(frc_scale)},
    inv_gpos_scale_f{static_cast<float>(inv_gpos_scale)},
    inv_lpos_scale_f{static_cast<float>(inv_lpos_scale)},
    inv_vel_scale_f{static_cast<float>(inv_vel_scale)},
    inv_frc_scale_f{static_cast<float>(inv_frc_scale)},
    gpos_bits{psyw.gpos_bits},
    lpos_bits{psyw.lpos_bits},
    vel_bits{psyw.vel_bits},
    frc_bits{psyw.frc_bits},
    boxvecs{psyw.boxvecs},
    umat{psyw.umat},
    invu{psyw.invu},
    boxdims{psyw.boxdims},
    sp_umat{psyw.sp_umat},
    sp_invu{psyw.sp_invu},
    sp_boxdims{psyw.sp_boxdims},
    xcrd{psyw.xcrd}, ycrd{psyw.ycrd}, zcrd{psyw.zcrd},
    xcrd_ovrf{psyw.xcrd_ovrf}, ycrd_ovrf{psyw.ycrd_ovrf}, zcrd_ovrf{psyw.zcrd_ovrf},
    xvel{psyw.xvel}, yvel{psyw.yvel}, zvel{psyw.zvel},
    xvel_ovrf{psyw.xvel_ovrf}, yvel_ovrf{psyw.yvel_ovrf}, zvel_ovrf{psyw.zvel_ovrf},
    xfrc{psyw.xfrc}, yfrc{psyw.yfrc}, zfrc{psyw.zfrc},
    xfrc_ovrf{psyw.xfrc_ovrf}, yfrc_ovrf{psyw.yfrc_ovrf}, zfrc_ovrf{psyw.zfrc_ovrf},
    xprv{psyw.xprv}, yprv{psyw.yprv}, zprv{psyw.zprv},
    xprv_ovrf{psyw.xprv_ovrf}, yprv_ovrf{psyw.yprv_ovrf}, zprv_ovrf{psyw.zprv_ovrf},
    xnxt{psyw.xnxt}, ynxt{psyw.ynxt}, znxt{psyw.znxt},
    xnxt_ovrf{psyw.xnxt_ovrf}, ynxt_ovrf{psyw.ynxt_ovrf}, znxt_ovrf{psyw.znxt_ovrf}
{}
  
//-------------------------------------------------------------------------------------------------
PhaseSpaceSynthesis::PhaseSpaceSynthesis(const std::vector<PhaseSpace> &ps_list,
                                         const std::vector<AtomGraph*> &ag_list,
                                         const std::vector<Thermostat> &heat_baths_in,
                                         const std::vector<Barostat> &pistons_in,
                                         const double time_step_in,
                                         const int globalpos_scale_bits_in,
                                         const int localpos_scale_bits_in,
                                         const int velocity_scale_bits_in,
                                         const int force_scale_bits_in) :
    system_count{static_cast<int>(ps_list.size())},
    unit_cell{UnitCellType::NONE},
    cycle_position{CoordinateCycle::PRESENT},
    heat_bath_kind{ThermostatKind::NONE},
    piston_kind{BarostatKind::NONE},
    time_step{time_step_in},
    globalpos_scale{pow(2.0, globalpos_scale_bits_in)},
    localpos_scale{pow(2.0, localpos_scale_bits_in)},
    velocity_scale{pow(2.0, velocity_scale_bits_in)},
    force_scale{pow(2.0, force_scale_bits_in)},
    inverse_globalpos_scale{1.0 / globalpos_scale},
    inverse_localpos_scale{1.0 / localpos_scale},
    inverse_velocity_scale{1.0 / velocity_scale},
    inverse_force_scale{1.0 / force_scale},
    globalpos_scale_bits{globalpos_scale_bits_in},
    localpos_scale_bits{localpos_scale_bits_in},
    velocity_scale_bits{velocity_scale_bits_in},
    force_scale_bits{force_scale_bits_in},
    atom_starts{HybridKind::POINTER, "labframe_starts"},
    atom_counts{HybridKind::POINTER, "labframe_counts"},
    x_coordinates{HybridKind::POINTER, "labframe_xpos"},
    y_coordinates{HybridKind::POINTER, "labframe_ypos"},
    z_coordinates{HybridKind::POINTER, "labframe_zpos"},
    x_coordinate_overflow{HybridKind::POINTER, "labframe_xpos_ovrf"},
    y_coordinate_overflow{HybridKind::POINTER, "labframe_ypos_ovrf"},
    z_coordinate_overflow{HybridKind::POINTER, "labframe_zpos_ovrf"},
    x_velocities{HybridKind::POINTER, "labframe_vx"},
    y_velocities{HybridKind::POINTER, "labframe_vy"},
    z_velocities{HybridKind::POINTER, "labframe_vz"},
    x_velocity_overflow{HybridKind::POINTER, "labframe_vx_ovrf"},
    y_velocity_overflow{HybridKind::POINTER, "labframe_vy_ovrf"},
    z_velocity_overflow{HybridKind::POINTER, "labframe_vz_ovrf"},
    x_forces{HybridKind::POINTER, "labframe_fx"},
    y_forces{HybridKind::POINTER, "labframe_fy"},
    z_forces{HybridKind::POINTER, "labframe_fz"},
    x_force_overflow{HybridKind::POINTER, "labframe_fx_ovrf"},
    y_force_overflow{HybridKind::POINTER, "labframe_fy_ovrf"},
    z_force_overflow{HybridKind::POINTER, "labframe_fz_ovrf"},
    x_prior_coordinates{HybridKind::POINTER, "labframe_x_prev"},
    y_prior_coordinates{HybridKind::POINTER, "labframe_y_prev"},
    z_prior_coordinates{HybridKind::POINTER, "labframe_z_prev"},
    x_prior_coord_overflow{HybridKind::POINTER, "labframe_xprv_ovrf"},
    y_prior_coord_overflow{HybridKind::POINTER, "labframe_yprv_ovrf"},
    z_prior_coord_overflow{HybridKind::POINTER, "labframe_zprv_ovrf"},
    x_future_coordinates{HybridKind::POINTER, "labframe_x_next"},
    y_future_coordinates{HybridKind::POINTER, "labframe_y_next"},
    z_future_coordinates{HybridKind::POINTER, "labframe_z_next"},
    x_future_coord_overflow{HybridKind::POINTER, "labframe_xnxt_ovrf"},
    y_future_coord_overflow{HybridKind::POINTER, "labframe_ynxt_ovrf"},
    z_future_coord_overflow{HybridKind::POINTER, "labframe_znxt_ovrf"},
    box_vectors{HybridKind::POINTER, "labframe_boxdims"},
    box_space_transforms{HybridKind::POINTER, "labframe_umat"},
    inverse_transforms{HybridKind::POINTER, "labframe_invu"},
    box_dimensions{HybridKind::POINTER, "labframe_dims"},
    sp_box_space_transforms{HybridKind::POINTER, "labframe_umat"},
    sp_inverse_transforms{HybridKind::POINTER, "labframe_invu"},
    sp_box_dimensions{HybridKind::POINTER, "labframe_dims"},
    int_data{HybridKind::ARRAY, "labframe_int"},
    llint_data{HybridKind::ARRAY, "labframe_llint"},
    double_data{HybridKind::ARRAY, "labframe_double"},
    float_data{HybridKind::ARRAY, "labframe_float"},
    topologies{ag_list.begin(), ag_list.end()}
{
  // Check validity of input
  checkGlobalPositionBits(globalpos_scale_bits);
  checkLocalPositionBits(localpos_scale_bits);
  checkVelocityBits(velocity_scale_bits);
  checkForceBits(force_scale_bits);
  const int nbaths = heat_baths_in.size();
  const int npumps = pistons_in.size();
  if (system_count == 0) {
    rtErr("At least one PhaseSpace object must be provided.", "PhaseSpaceSynthesis");
  }
  else if (ag_list.size() != system_count) {
    rtErr("One topology pointer must be provided for each coordinate system (currently " +
          std::to_string(ag_list.size()) + " topology pointers and " +
          std::to_string(system_count) + " PhaseSpace objects.", "PhaseSpaceSynthesis");
  }
  else if (nbaths != 1 && nbaths != system_count) {
    rtErr("One thermostat must be provided, or one thermostat for each system (currently " +
          std::to_string(system_count) + " systems and " + std::to_string(nbaths) +
          " thermostats).", "PhaseSpaceSynthesis");
  }
  else if (npumps != 1 && npumps != system_count) {
    rtErr("One barostat must be provided, or one barostat for each system (currently " +
          std::to_string(system_count) + " systems and " + std::to_string(npumps) +
          " barostats).", "PhaseSpaceSynthesis");
  }

  // Check validity of thermostats: all thermostats must be of the same type, but temperatures
  // and details of the coupling strength may vary between systems
  heat_bath_kind = heat_baths_in[0].getKind();
  for (int i = 1; i < nbaths; i++) {
    if (heat_baths_in[i].getKind() != heat_baths_in[0].getKind()) {
      rtErr("All system thermostats must be of the same class.  System " + std::to_string(i) +
            " uses a " + getThermostatName(heat_baths_in[i].getKind()) + " thermostat whereas the "
            "first system locks the method to " + getThermostatName(heat_baths_in[0].getKind()) +
            ".", "PhaseSpaceSynthesis");
    }
  }
  
  // Check validity of barostats: all barostats must be of the same type, but pressures and
  // details of the rescaling may vary between systems
  piston_kind = pistons_in[0].getKind();
  for (int i = 1; i < npumps; i++) {
    if (pistons_in[i].getKind() != pistons_in[0].getKind()) {
      rtErr("All system barostats must be of the same class.  System " + std::to_string(i) +
            " uses a " + getBarostatName(pistons_in[i].getKind()) + " barostat whereas the first "
            "system locks the method to " + getBarostatName(pistons_in[0].getKind()) + ".",
            "PhaseSpaceSynthesis");
    }
  }

  // Allocate data and set internal pointers
  int atom_stride = 0;
  for (int i = 0; i < system_count; i++) {
    atom_stride += roundUp(ps_list[i].getAtomCount(), warp_size_int);
  }
  allocate(atom_stride);

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

  // Loop over all systems
  const int mtrx_stride = roundUp(9, warp_size_int);
  const int dim_stride = roundUp(6, warp_size_int);
  for (int i = 0; i < system_count; i++) {

    // Get a reader for the PhaseSpace object's host-side data
    PhaseSpaceReader psr = ps_list[i].data();

    // Assign coordinates, converting from double precision to fixed precision format
    llint* xpos_ptr = x_coordinates.data();
    llint* ypos_ptr = y_coordinates.data();
    llint* zpos_ptr = z_coordinates.data();
    int* xpos_ovrf_ptr  = x_coordinate_overflow.data();
    int* ypos_ovrf_ptr  = y_coordinate_overflow.data();
    int* zpos_ovrf_ptr  = z_coordinate_overflow.data();
    llint* xvel_ptr = x_velocities.data();
    llint* yvel_ptr = y_velocities.data();
    llint* zvel_ptr = z_velocities.data();
    int* xvel_ovrf_ptr  = x_velocity_overflow.data();
    int* yvel_ovrf_ptr  = y_velocity_overflow.data();
    int* zvel_ovrf_ptr  = z_velocity_overflow.data();
    llint* xfrc_ptr = x_forces.data();
    llint* yfrc_ptr = y_forces.data();
    llint* zfrc_ptr = z_forces.data();
    int* xfrc_ovrf_ptr  = x_force_overflow.data();
    int* yfrc_ovrf_ptr  = y_force_overflow.data();
    int* zfrc_ovrf_ptr  = z_force_overflow.data();
    const int asi = atom_starts.readHost(i);
    for (int j = 0; j < psr.natom; j++) {

      // The x, y, and z components hold x, y, and z positions in a "lab frame" unwrapped form.
      // In this conversion, lower-precision forms of each coordinate (positions, velocities, or
      // forces), those that most simulations will run with, are rounded to the nearest integer
      // values rather than truncated rounding towards zero.  Very high precision formats, those
      // with more than 53-54 bits after the decimal, are likely to need no rounding as the
      // fixed precision format will already be more precise than the double-precision floating
      // point number's mantissa.  However, there is a "no man's land" in this implementation,
      // between the non-overflow bit counts for each coordinate (36-44 bits) and the point at
      // which the fixed-precision representation will capture all of the information in nearly
      // all numbers, where rounding towards zero will occur.
      const size_t asi_j = asi + j;
      if (globalpos_scale_bits <= globalpos_scale_nonoverflow_bits) {
        xpos_ptr[asi_j] = llround(psr.xcrd[j] * globalpos_scale);
        ypos_ptr[asi_j] = llround(psr.ycrd[j] * globalpos_scale);
        zpos_ptr[asi_j] = llround(psr.zcrd[j] * globalpos_scale);
      }
      else {
        splitRealConversion(psr.xcrd[j] * globalpos_scale, &xpos_ptr[asi_j],
                            &xpos_ovrf_ptr[asi_j]);
        splitRealConversion(psr.ycrd[j] * globalpos_scale, &ypos_ptr[asi_j],
                            &ypos_ovrf_ptr[asi_j]);
        splitRealConversion(psr.zcrd[j] * globalpos_scale, &zpos_ptr[asi_j],
                            &zpos_ovrf_ptr[asi_j]);
      }
      if (velocity_scale_bits <= velocity_scale_nonoverflow_bits) {
        xvel_ptr[asi_j] = llround(psr.xvel[j] * velocity_scale);
        yvel_ptr[asi_j] = llround(psr.yvel[j] * velocity_scale);
        zvel_ptr[asi_j] = llround(psr.zvel[j] * velocity_scale);
      }
      else {
        splitRealConversion(psr.xvel[j] * globalpos_scale, &xvel_ptr[asi_j],
                            &xvel_ovrf_ptr[asi_j]);
        splitRealConversion(psr.yvel[j] * globalpos_scale, &yvel_ptr[asi_j],
                            &yvel_ovrf_ptr[asi_j]);
        splitRealConversion(psr.zvel[j] * globalpos_scale, &zvel_ptr[asi_j],
                            &zvel_ovrf_ptr[asi_j]);
      }
      if (force_scale_bits <= force_scale_nonoverflow_bits) {
        xfrc_ptr[asi_j] = llround(psr.xfrc[j] * force_scale);
        yfrc_ptr[asi_j] = llround(psr.yfrc[j] * force_scale);
        zfrc_ptr[asi_j] = llround(psr.zfrc[j] * force_scale);
      }
      else {
        splitRealConversion(psr.xfrc[j] * globalpos_scale, &xfrc_ptr[asi_j],
                            &xfrc_ovrf_ptr[asi_j]);
        splitRealConversion(psr.yfrc[j] * globalpos_scale, &yfrc_ptr[asi_j],
                            &yfrc_ovrf_ptr[asi_j]);
        splitRealConversion(psr.zfrc[j] * globalpos_scale, &zfrc_ptr[asi_j],
                            &zfrc_ovrf_ptr[asi_j]);
      }
    }

    // Handle the box space transformation.  The transformation matrices of the labframe will
    // become slightly desynchronized from the original PhaseSpace object, but this is in the
    // interest of making a system which can be represented in fixed precision, pixelated
    // coordinates with box vectors which are likewise multiples of the positional discretization.
    for (int j = 0; j < 9; j++) {
      const llint lli_invu = llround(psr.invu[j] * transform_scale_lf);
      const double d_invu = static_cast<double>(lli_invu);
      inverse_transforms.putHost(d_invu, (i * mtrx_stride) + j);
      sp_inverse_transforms.putHost(d_invu, (i * mtrx_stride) + j);
    }
    const double* invu_ptr = &inverse_transforms.data()[i * mtrx_stride];
    double* umat_ptr = &box_space_transforms.data()[i * mtrx_stride];
    invertSquareMatrix(invu_ptr, umat_ptr, 3);
    for (int j = 0; j < 6; j++) {
      box_dimensions.putHost(psr.boxdim[j], (i * dim_stride) + j);
      sp_box_dimensions.putHost(psr.boxdim[j], (i * dim_stride) + j);
    }
  }
  for (int i = 0; i < system_count * mtrx_stride; i++) {
    sp_box_space_transforms.putHost(box_space_transforms.readHost(i), i);
  }
}

//-------------------------------------------------------------------------------------------------
PhaseSpaceSynthesis::PhaseSpaceSynthesis(const std::vector<PhaseSpace> &ps_list,
                                         const std::vector<AtomGraph*> &ag_list,
                                         const std::vector<int> &index_key,
                                         const std::vector<Thermostat> &heat_baths_in,
                                         const std::vector<Barostat> &pistons_in,
                                         const double time_step_in,
                                         const int globalpos_scale_bits_in,
                                         const int localpos_scale_bits_in,
                                         const int velocity_scale_bits_in,
                                         const int force_scale_bits_in) :
    PhaseSpaceSynthesis(tileVector(ps_list, index_key), tileVector(ag_list, index_key),
                        heat_baths_in, pistons_in, time_step_in, globalpos_scale_bits_in,
                        localpos_scale_bits_in, velocity_scale_bits_in, force_scale_bits_in)
{}

//-------------------------------------------------------------------------------------------------
PhaseSpaceSynthesis::PhaseSpaceSynthesis(const std::vector<PhaseSpace> &ps_list,
                                         const std::vector<int> &ps_index_key,
                                         const std::vector<AtomGraph*> &ag_list,
                                         const std::vector<int> &ag_index_key,
                                         const std::vector<Thermostat> &heat_baths_in,
                                         const std::vector<Barostat> &pistons_in,
                                         const double time_step_in,
                                         const int globalpos_scale_bits_in,
                                         const int localpos_scale_bits_in,
                                         const int velocity_scale_bits_in,
                                         const int force_scale_bits_in) :
    PhaseSpaceSynthesis(tileVector(ps_list, ps_index_key), tileVector(ag_list, ag_index_key),
                        heat_baths_in, pistons_in, time_step_in, globalpos_scale_bits_in,
                        localpos_scale_bits_in, velocity_scale_bits_in, force_scale_bits_in)
{}

//-------------------------------------------------------------------------------------------------
PhaseSpaceSynthesis::PhaseSpaceSynthesis(const SystemCache &sysc,
                                         const std::vector<Thermostat> &heat_baths_in,
                                         const std::vector<Barostat> &pistons_in,
                                         const double time_step_in,
                                         const int globalpos_scale_bits_in,
                                         const int localpos_scale_bits_in,
                                         const int velocity_scale_bits_in,
                                         const int force_scale_bits_in) :
    PhaseSpaceSynthesis(sysc.getCoordinateReference(), sysc.getSystemTopologyPointerCC(),
                        heat_baths_in, pistons_in, time_step_in, globalpos_scale_bits_in,
                        localpos_scale_bits_in, velocity_scale_bits_in, force_scale_bits_in)
{}

//-------------------------------------------------------------------------------------------------
PhaseSpaceSynthesis::PhaseSpaceSynthesis(const std::vector<PhaseSpace> &ps_list,
                                         const std::vector<AtomGraph*> &ag_list,
                                         const int globalpos_scale_bits_in,
                                         const int localpos_scale_bits_in,
                                         const int velocity_scale_bits_in,
                                         const int force_scale_bits_in) :
    PhaseSpaceSynthesis(ps_list, ag_list, { Thermostat() }, { Barostat() }, 1.0,
                        globalpos_scale_bits_in, localpos_scale_bits_in, velocity_scale_bits_in,
                        force_scale_bits_in)
{}

//-------------------------------------------------------------------------------------------------
PhaseSpaceSynthesis::PhaseSpaceSynthesis(const std::vector<PhaseSpace> &ps_list,
                                         const std::vector<AtomGraph*> &ag_list,
                                         const std::vector<int> &index_key,
                                         const int globalpos_scale_bits_in,
                                         const int localpos_scale_bits_in,
                                         const int velocity_scale_bits_in,
                                         const int force_scale_bits_in) :
    PhaseSpaceSynthesis(tileVector(ps_list, index_key), tileVector(ag_list, index_key),
                        { Thermostat() }, { Barostat() }, 1.0, globalpos_scale_bits_in,
                        localpos_scale_bits_in, velocity_scale_bits_in, force_scale_bits_in)
{}

//-------------------------------------------------------------------------------------------------
PhaseSpaceSynthesis::PhaseSpaceSynthesis(const std::vector<PhaseSpace> &ps_list,
                                         const std::vector<int> &ps_index_key,
                                         const std::vector<AtomGraph*> &ag_list,
                                         const std::vector<int> &ag_index_key,
                                         const int globalpos_scale_bits_in,
                                         const int localpos_scale_bits_in,
                                         const int velocity_scale_bits_in,
                                         const int force_scale_bits_in) :
    PhaseSpaceSynthesis(tileVector(ps_list, ps_index_key), tileVector(ag_list, ag_index_key),
                        { Thermostat() }, { Barostat() }, 1.0, globalpos_scale_bits_in,
                        localpos_scale_bits_in, velocity_scale_bits_in, force_scale_bits_in)
{}

//-------------------------------------------------------------------------------------------------
PhaseSpaceSynthesis::PhaseSpaceSynthesis(const SystemCache &sysc,
                                         const int globalpos_scale_bits_in,
                                         const int localpos_scale_bits_in,
                                         const int velocity_scale_bits_in,
                                         const int force_scale_bits_in) :
    PhaseSpaceSynthesis(sysc.getCoordinateReference(), sysc.getSystemTopologyPointerCC(),
                        { Thermostat() }, { Barostat() }, 1.0, globalpos_scale_bits_in,
                        localpos_scale_bits_in, velocity_scale_bits_in, force_scale_bits_in)
{}

//-------------------------------------------------------------------------------------------------
PhaseSpaceSynthesis::PhaseSpaceSynthesis(const PhaseSpaceSynthesis &original) :
    system_count{original.system_count},
    unit_cell{original.unit_cell},
    cycle_position{original.cycle_position},
    heat_bath_kind{original.heat_bath_kind},
    piston_kind{original.piston_kind},
    time_step{original.time_step},
    globalpos_scale{original.globalpos_scale},
    localpos_scale{original.localpos_scale},
    velocity_scale{original.velocity_scale},
    force_scale{original.force_scale},
    inverse_globalpos_scale{original.inverse_globalpos_scale},
    inverse_localpos_scale{original.inverse_localpos_scale},
    inverse_velocity_scale{original.inverse_velocity_scale},
    inverse_force_scale{original.inverse_force_scale},
    globalpos_scale_bits{original.globalpos_scale_bits},
    localpos_scale_bits{original.localpos_scale_bits},
    velocity_scale_bits{original.velocity_scale_bits},
    force_scale_bits{original.force_scale_bits},
    atom_starts{original.atom_starts},
    atom_counts{original.atom_counts},
    x_coordinates{original.x_coordinates},
    y_coordinates{original.y_coordinates},
    z_coordinates{original.z_coordinates},
    x_coordinate_overflow{original.x_coordinate_overflow},
    y_coordinate_overflow{original.y_coordinate_overflow},
    z_coordinate_overflow{original.z_coordinate_overflow},
    x_velocities{original.x_velocities},
    y_velocities{original.y_velocities},
    z_velocities{original.z_velocities},
    x_velocity_overflow{original.x_velocity_overflow},
    y_velocity_overflow{original.y_velocity_overflow},
    z_velocity_overflow{original.z_velocity_overflow},
    x_forces{original.x_forces},
    y_forces{original.y_forces},
    z_forces{original.z_forces},
    x_force_overflow{original.x_force_overflow},
    y_force_overflow{original.y_force_overflow},
    z_force_overflow{original.z_force_overflow},
    x_prior_coordinates{original.x_prior_coordinates},
    y_prior_coordinates{original.y_prior_coordinates},
    z_prior_coordinates{original.z_prior_coordinates},
    x_prior_coord_overflow{original.x_prior_coord_overflow},
    y_prior_coord_overflow{original.y_prior_coord_overflow},
    z_prior_coord_overflow{original.z_prior_coord_overflow},
    x_future_coordinates{original.x_future_coordinates},
    y_future_coordinates{original.y_future_coordinates},
    z_future_coordinates{original.z_future_coordinates},
    x_future_coord_overflow{original.x_future_coord_overflow},
    y_future_coord_overflow{original.y_future_coord_overflow},
    z_future_coord_overflow{original.z_future_coord_overflow},
    box_vectors{original.box_vectors},
    box_space_transforms{original.box_space_transforms},
    inverse_transforms{original.inverse_transforms},
    box_dimensions{original.box_dimensions},
    sp_box_space_transforms{original.sp_box_space_transforms},
    sp_inverse_transforms{original.sp_inverse_transforms},
    sp_box_dimensions{original.sp_box_dimensions},
    int_data{original.int_data},
    llint_data{original.llint_data},
    double_data{original.double_data},
    float_data{original.float_data},
    topologies{original.topologies}
{
  // The allocate function again handle pointer repair, just like in the PhaseSpace object.
  // Sum the atom stride based on the AtomGraph pointers, as the PhaseSpace objects that created
  // the original must have been in agreement in order for it to exist in the first place.  The
  // Hybrid objects will not be resized as they already have the proper sizes.
  int atom_stride = 0;
  for (int i = 0; i < system_count; i++) {
    atom_stride += roundUp(topologies[i]->getAtomCount(), warp_size_int);
  }
  allocate(atom_stride);
}

//-------------------------------------------------------------------------------------------------
PhaseSpaceSynthesis::PhaseSpaceSynthesis(PhaseSpaceSynthesis &&original) :
    system_count{original.system_count},
    unit_cell{original.unit_cell},
    cycle_position{original.cycle_position},
    heat_bath_kind{original.heat_bath_kind},
    piston_kind{original.piston_kind},
    time_step{original.time_step},
    globalpos_scale{original.globalpos_scale},
    localpos_scale{original.localpos_scale},
    velocity_scale{original.velocity_scale},
    force_scale{original.force_scale},
    inverse_globalpos_scale{original.inverse_globalpos_scale},
    inverse_localpos_scale{original.inverse_localpos_scale},
    inverse_velocity_scale{original.inverse_velocity_scale},
    inverse_force_scale{original.inverse_force_scale},
    globalpos_scale_bits{original.globalpos_scale_bits},
    localpos_scale_bits{original.localpos_scale_bits},
    velocity_scale_bits{original.velocity_scale_bits},
    force_scale_bits{original.force_scale_bits},
    atom_starts{std::move(original.atom_starts)},
    atom_counts{std::move(original.atom_counts)},
    x_coordinates{std::move(original.x_coordinates)},
    y_coordinates{std::move(original.y_coordinates)},
    z_coordinates{std::move(original.z_coordinates)},
    x_coordinate_overflow{std::move(original.x_coordinate_overflow)},
    y_coordinate_overflow{std::move(original.y_coordinate_overflow)},
    z_coordinate_overflow{std::move(original.z_coordinate_overflow)},
    x_velocities{std::move(original.x_velocities)},
    y_velocities{std::move(original.y_velocities)},
    z_velocities{std::move(original.z_velocities)},
    x_velocity_overflow{std::move(original.x_velocity_overflow)},
    y_velocity_overflow{std::move(original.y_velocity_overflow)},
    z_velocity_overflow{std::move(original.z_velocity_overflow)},
    x_forces{std::move(original.x_forces)},
    y_forces{std::move(original.y_forces)},
    z_forces{std::move(original.z_forces)},
    x_force_overflow{std::move(original.x_force_overflow)},
    y_force_overflow{std::move(original.y_force_overflow)},
    z_force_overflow{std::move(original.z_force_overflow)},
    x_prior_coordinates{std::move(original.x_prior_coordinates)},
    y_prior_coordinates{std::move(original.y_prior_coordinates)},
    z_prior_coordinates{std::move(original.z_prior_coordinates)},
    x_prior_coord_overflow{std::move(original.x_prior_coord_overflow)},
    y_prior_coord_overflow{std::move(original.y_prior_coord_overflow)},
    z_prior_coord_overflow{std::move(original.z_prior_coord_overflow)},
    x_future_coordinates{std::move(original.x_future_coordinates)},
    y_future_coordinates{std::move(original.y_future_coordinates)},
    z_future_coordinates{std::move(original.z_future_coordinates)},
    x_future_coord_overflow{std::move(original.x_future_coord_overflow)},
    y_future_coord_overflow{std::move(original.y_future_coord_overflow)},
    z_future_coord_overflow{std::move(original.z_future_coord_overflow)},
    box_vectors{std::move(original.box_vectors)},
    box_space_transforms{std::move(original.box_space_transforms)},
    inverse_transforms{std::move(original.inverse_transforms)},
    box_dimensions{std::move(original.box_dimensions)},
    sp_box_space_transforms{std::move(original.sp_box_space_transforms)},
    sp_inverse_transforms{std::move(original.sp_inverse_transforms)},
    sp_box_dimensions{std::move(original.sp_box_dimensions)},
    int_data{std::move(original.int_data)},
    llint_data{std::move(original.llint_data)},
    double_data{std::move(original.double_data)},
    float_data{std::move(original.float_data)}
{}

//-------------------------------------------------------------------------------------------------
const PsSynthesisReader PhaseSpaceSynthesis::data(const HybridTargetLevel tier) const {
  return data(cycle_position, tier);
}

//-------------------------------------------------------------------------------------------------
const PsSynthesisReader PhaseSpaceSynthesis::data(const CoordinateCycle orientation,
                                                  const HybridTargetLevel tier) const {
  switch (orientation) {
  case CoordinateCycle::PAST:
    return PsSynthesisReader(system_count, unit_cell, heat_bath_kind, piston_kind, time_step,
                             atom_starts.data(tier), atom_counts.data(tier), globalpos_scale,
                             localpos_scale, velocity_scale, force_scale, globalpos_scale_bits,
                             localpos_scale_bits, velocity_scale_bits, force_scale_bits,
                             box_vectors.data(tier), box_space_transforms.data(tier),
                             inverse_transforms.data(tier), box_dimensions.data(tier),
                             sp_box_space_transforms.data(tier), sp_inverse_transforms.data(tier),
                             sp_box_dimensions.data(tier), x_prior_coordinates.data(tier),
                             y_prior_coordinates.data(tier), z_prior_coordinates.data(tier),
                             x_prior_coord_overflow.data(tier), y_prior_coord_overflow.data(tier),
                             z_prior_coord_overflow.data(tier), x_velocities.data(tier),
                             y_velocities.data(tier), z_velocities.data(tier),
                             x_velocity_overflow.data(tier), y_velocity_overflow.data(tier),
                             z_velocity_overflow.data(tier), x_forces.data(tier),
                             y_forces.data(tier), z_forces.data(tier), x_force_overflow.data(tier),
                             y_force_overflow.data(tier), z_force_overflow.data(tier),
                             x_future_coordinates.data(tier), y_future_coordinates.data(tier),
                             z_future_coordinates.data(tier), x_future_coord_overflow.data(tier),
                             y_future_coord_overflow.data(tier),
                             z_future_coord_overflow.data(tier),
                             x_coordinates.data(tier), y_coordinates.data(tier),
                             z_coordinates.data(tier), x_coordinate_overflow.data(tier),
                             y_coordinate_overflow.data(tier), z_coordinate_overflow.data(tier));
  case CoordinateCycle::PRESENT:
    return PsSynthesisReader(system_count, unit_cell, heat_bath_kind, piston_kind, time_step,
                             atom_starts.data(tier), atom_counts.data(tier), globalpos_scale,
                             localpos_scale, velocity_scale, force_scale, globalpos_scale_bits,
                             localpos_scale_bits, velocity_scale_bits, force_scale_bits,
                             box_vectors.data(tier), box_space_transforms.data(tier),
                             inverse_transforms.data(tier), box_dimensions.data(tier),
                             sp_box_space_transforms.data(tier), sp_inverse_transforms.data(tier),
                             sp_box_dimensions.data(tier), x_coordinates.data(tier),
                             y_coordinates.data(tier), z_coordinates.data(tier),
                             x_coordinate_overflow.data(tier), y_coordinate_overflow.data(tier),
                             z_coordinate_overflow.data(tier), x_velocities.data(tier),
                             y_velocities.data(tier), z_velocities.data(tier),
                             x_velocity_overflow.data(tier), y_velocity_overflow.data(tier),
                             z_velocity_overflow.data(tier), x_forces.data(tier),
                             y_forces.data(tier), z_forces.data(tier), x_force_overflow.data(tier),
                             y_force_overflow.data(tier), z_force_overflow.data(tier),
                             x_prior_coordinates.data(tier), y_prior_coordinates.data(tier),
                             z_prior_coordinates.data(tier), x_prior_coord_overflow.data(tier),
                             y_prior_coord_overflow.data(tier), z_prior_coord_overflow.data(tier),
                             x_future_coordinates.data(tier), y_future_coordinates.data(tier),
                             z_future_coordinates.data(tier), x_future_coord_overflow.data(tier),
                             y_future_coord_overflow.data(tier),
                             z_future_coord_overflow.data(tier));
  case CoordinateCycle::FUTURE:
    return PsSynthesisReader(system_count, unit_cell, heat_bath_kind, piston_kind, time_step,
                             atom_starts.data(tier), atom_counts.data(tier), globalpos_scale,
                             localpos_scale, velocity_scale, force_scale, globalpos_scale_bits,
                             localpos_scale_bits, velocity_scale_bits, force_scale_bits,
                             box_vectors.data(tier), box_space_transforms.data(tier),
                             inverse_transforms.data(tier), box_dimensions.data(tier),
                             sp_box_space_transforms.data(tier), sp_inverse_transforms.data(tier),
                             sp_box_dimensions.data(tier), x_future_coordinates.data(tier),
                             y_future_coordinates.data(tier), z_future_coordinates.data(tier),
                             x_future_coord_overflow.data(tier),
                             y_future_coord_overflow.data(tier),
                             z_future_coord_overflow.data(tier), x_velocities.data(tier),
                             y_velocities.data(tier), z_velocities.data(tier),
                             x_velocity_overflow.data(tier), y_velocity_overflow.data(tier),
                             z_velocity_overflow.data(tier), x_forces.data(tier),
                             y_forces.data(tier), z_forces.data(tier), x_force_overflow.data(tier),
                             y_force_overflow.data(tier), z_force_overflow.data(tier),
                             x_coordinates.data(tier), y_coordinates.data(tier),
                             z_coordinates.data(tier), x_coordinate_overflow.data(tier),
                             y_coordinate_overflow.data(tier), z_coordinate_overflow.data(tier),
                             x_prior_coordinates.data(tier), y_prior_coordinates.data(tier),
                             z_prior_coordinates.data(tier), x_prior_coord_overflow.data(tier),
                             y_prior_coord_overflow.data(tier), z_prior_coord_overflow.data(tier));
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
PsSynthesisWriter PhaseSpaceSynthesis::data(const HybridTargetLevel tier) {
  return data(cycle_position, tier);
}

//-------------------------------------------------------------------------------------------------
PsSynthesisWriter PhaseSpaceSynthesis::data(const CoordinateCycle orientation,
                                            const HybridTargetLevel tier) {
  switch (orientation) {
  case CoordinateCycle::PAST:
    return PsSynthesisWriter(system_count, unit_cell, heat_bath_kind, piston_kind, time_step,
                             atom_starts.data(tier), atom_counts.data(tier), globalpos_scale,
                             localpos_scale, velocity_scale, force_scale, globalpos_scale_bits,
                             localpos_scale_bits, velocity_scale_bits, force_scale_bits,
                             box_vectors.data(tier), box_space_transforms.data(tier),
                             inverse_transforms.data(tier), box_dimensions.data(tier),
                             sp_box_space_transforms.data(tier), sp_inverse_transforms.data(tier),
                             sp_box_dimensions.data(tier), x_prior_coordinates.data(tier),
                             y_prior_coordinates.data(tier), z_prior_coordinates.data(tier),
                             x_prior_coord_overflow.data(tier), y_prior_coord_overflow.data(tier),
                             z_prior_coord_overflow.data(tier), x_velocities.data(tier),
                             y_velocities.data(tier), z_velocities.data(tier),
                             x_velocity_overflow.data(tier), y_velocity_overflow.data(tier),
                             z_velocity_overflow.data(tier), x_forces.data(tier),
                             y_forces.data(tier), z_forces.data(tier), x_force_overflow.data(tier),
                             y_force_overflow.data(tier), z_force_overflow.data(tier),
                             x_future_coordinates.data(tier), y_future_coordinates.data(tier),
                             z_future_coordinates.data(tier), x_future_coord_overflow.data(tier),
                             y_future_coord_overflow.data(tier),
                             z_future_coord_overflow.data(tier),
                             x_coordinates.data(tier), y_coordinates.data(tier),
                             z_coordinates.data(tier), x_coordinate_overflow.data(tier),
                             y_coordinate_overflow.data(tier), z_coordinate_overflow.data(tier));
  case CoordinateCycle::PRESENT:
    return PsSynthesisWriter(system_count, unit_cell, heat_bath_kind, piston_kind, time_step,
                             atom_starts.data(tier), atom_counts.data(tier), globalpos_scale,
                             localpos_scale, velocity_scale, force_scale, globalpos_scale_bits,
                             localpos_scale_bits, velocity_scale_bits, force_scale_bits,
                             box_vectors.data(tier), box_space_transforms.data(tier),
                             inverse_transforms.data(tier), box_dimensions.data(tier),
                             sp_box_space_transforms.data(tier), sp_inverse_transforms.data(tier),
                             sp_box_dimensions.data(tier), x_coordinates.data(tier),
                             y_coordinates.data(tier), z_coordinates.data(tier),
                             x_coordinate_overflow.data(tier), y_coordinate_overflow.data(tier),
                             z_coordinate_overflow.data(tier), x_velocities.data(tier),
                             y_velocities.data(tier), z_velocities.data(tier),
                             x_velocity_overflow.data(tier), y_velocity_overflow.data(tier),
                             z_velocity_overflow.data(tier), x_forces.data(tier),
                             y_forces.data(tier), z_forces.data(tier), x_force_overflow.data(tier),
                             y_force_overflow.data(tier), z_force_overflow.data(tier),
                             x_prior_coordinates.data(tier), y_prior_coordinates.data(tier),
                             z_prior_coordinates.data(tier), x_prior_coord_overflow.data(tier),
                             y_prior_coord_overflow.data(tier), z_prior_coord_overflow.data(tier),
                             x_future_coordinates.data(tier), y_future_coordinates.data(tier),
                             z_future_coordinates.data(tier), x_future_coord_overflow.data(tier),
                             y_future_coord_overflow.data(tier),
                             z_future_coord_overflow.data(tier));
  case CoordinateCycle::FUTURE:
    return PsSynthesisWriter(system_count, unit_cell, heat_bath_kind, piston_kind, time_step,
                             atom_starts.data(tier), atom_counts.data(tier), globalpos_scale,
                             localpos_scale, velocity_scale, force_scale, globalpos_scale_bits,
                             localpos_scale_bits, velocity_scale_bits, force_scale_bits,
                             box_vectors.data(tier), box_space_transforms.data(tier),
                             inverse_transforms.data(tier), box_dimensions.data(tier),
                             sp_box_space_transforms.data(tier), sp_inverse_transforms.data(tier),
                             sp_box_dimensions.data(tier), x_future_coordinates.data(tier),
                             y_future_coordinates.data(tier), z_future_coordinates.data(tier),
                             x_future_coord_overflow.data(tier),
                             y_future_coord_overflow.data(tier),
                             z_future_coord_overflow.data(tier), x_velocities.data(tier),
                             y_velocities.data(tier), z_velocities.data(tier),
                             x_velocity_overflow.data(tier), y_velocity_overflow.data(tier),
                             z_velocity_overflow.data(tier), x_forces.data(tier),
                             y_forces.data(tier), z_forces.data(tier), x_force_overflow.data(tier),
                             y_force_overflow.data(tier), z_force_overflow.data(tier),
                             x_coordinates.data(tier), y_coordinates.data(tier),
                             z_coordinates.data(tier), x_coordinate_overflow.data(tier),
                             y_coordinate_overflow.data(tier), z_coordinate_overflow.data(tier),
                             x_prior_coordinates.data(tier), y_prior_coordinates.data(tier),
                             z_prior_coordinates.data(tier), x_prior_coord_overflow.data(tier),
                             y_prior_coord_overflow.data(tier), z_prior_coord_overflow.data(tier));
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
int PhaseSpaceSynthesis::getSystemCount() const {
  return system_count;
}

//-------------------------------------------------------------------------------------------------
UnitCellType PhaseSpaceSynthesis::getUnitCellType() const {
  return unit_cell;
}

//-------------------------------------------------------------------------------------------------
CoordinateCycle PhaseSpaceSynthesis::getCyclePosition() const {
  return cycle_position;
}

//-------------------------------------------------------------------------------------------------
int PhaseSpaceSynthesis::getGlobalPositionBits() const {
  return globalpos_scale_bits;
}

//-------------------------------------------------------------------------------------------------
int PhaseSpaceSynthesis::getLocalPositionBits() const {
  return localpos_scale_bits;
}

//-------------------------------------------------------------------------------------------------
int PhaseSpaceSynthesis::getVelocityBits() const {
  return velocity_scale_bits;
}

//-------------------------------------------------------------------------------------------------
int PhaseSpaceSynthesis::getForceAccumulationBits() const {
  return force_scale_bits;
}

#ifdef STORMM_USE_HPC
//-------------------------------------------------------------------------------------------------
PsSynthesisWriter PhaseSpaceSynthesis::deviceViewToHostData() {
  return deviceViewToHostData(cycle_position);
}

//-------------------------------------------------------------------------------------------------
PsSynthesisWriter PhaseSpaceSynthesis::deviceViewToHostData(const CoordinateCycle orientation) {

  // Some of this could be expedited by retracing the manner in which some of the ARRAY-kind
  // Hybrids were allocated to infer the values of subsequent pointers, but
  // cudaHostGetDevicePointer is not a high-latency or laborious call.  Reduce complexity in the
  // central Hybrid object by generating the pointers every time this special case occurs.
  int*    devc_atom_starts;
  int*    devc_atom_counts;
  llint*  devc_boxvecs;
  double* devc_umat;
  double* devc_invu;
  double* devc_boxdims;
  float*  devc_sp_umat;
  float*  devc_sp_invu;
  float*  devc_sp_boxdims;
  llint*  devc_xcrd;
  llint*  devc_ycrd;
  llint*  devc_zcrd;
  int*    devc_xcrd_ovrf;
  int*    devc_ycrd_ovrf;
  int*    devc_zcrd_ovrf;
  llint*  devc_xvel;
  llint*  devc_yvel;
  llint*  devc_zvel;
  int*    devc_xvel_ovrf;
  int*    devc_yvel_ovrf;
  int*    devc_zvel_ovrf;
  llint*  devc_xfrc;
  llint*  devc_yfrc;
  llint*  devc_zfrc;
  int*    devc_xfrc_ovrf;
  int*    devc_yfrc_ovrf;
  int*    devc_zfrc_ovrf;
  llint*  devc_xprv;
  llint*  devc_yprv;
  llint*  devc_zprv;
  int*    devc_xprv_ovrf;
  int*    devc_yprv_ovrf;
  int*    devc_zprv_ovrf;
  llint*  devc_xnxt;
  llint*  devc_ynxt;
  llint*  devc_znxt;
  int*    devc_xnxt_ovrf;
  int*    devc_ynxt_ovrf;
  int*    devc_znxt_ovrf;
  bool problem = false;
#  ifdef STORMM_USE_CUDA
  problem = (problem || cudaHostGetDevicePointer((void **)&devc_atom_starts,
                                                 (void *)atom_starts.data(), 0) != cudaSuccess);
  problem = (problem || cudaHostGetDevicePointer((void **)&devc_atom_counts,
                                                 (void *)atom_counts.data(), 0) != cudaSuccess);
  problem = (problem || cudaHostGetDevicePointer((void **)&devc_boxvecs,
                                                 (void *)box_vectors.data(), 0) != cudaSuccess);
  problem = (problem ||
             cudaHostGetDevicePointer((void **)&devc_umat,
                                      (void *)box_space_transforms.data(), 0) != cudaSuccess);
  problem = (problem ||  
             cudaHostGetDevicePointer((void **)&devc_invu,
                                      (void *)inverse_transforms.data(), 0) != cudaSuccess);
  problem = (problem || cudaHostGetDevicePointer((void **)&devc_boxdims,
                                                 (void *)box_dimensions.data(), 0) != cudaSuccess);
  problem = (problem ||
             cudaHostGetDevicePointer((void **)&devc_sp_umat,
                                      (void *)sp_box_space_transforms.data(), 0) != cudaSuccess);
  problem = (problem ||  
             cudaHostGetDevicePointer((void **)&devc_sp_invu,
                                      (void *)sp_inverse_transforms.data(), 0) != cudaSuccess);
  problem = (problem ||
             cudaHostGetDevicePointer((void **)&devc_sp_boxdims,
                                      (void *)sp_box_dimensions.data(), 0) != cudaSuccess);
  problem = (problem || cudaHostGetDevicePointer((void **)&devc_xcrd,
                                                 (void *)x_coordinates.data(), 0) != cudaSuccess);
  problem = (problem || cudaHostGetDevicePointer((void **)&devc_ycrd,
                                                 (void *)y_coordinates.data(), 0) != cudaSuccess);
  problem = (problem || cudaHostGetDevicePointer((void **)&devc_zcrd,
                                                 (void *)z_coordinates.data(), 0) != cudaSuccess);
  problem = (problem ||
             cudaHostGetDevicePointer((void **)&devc_xcrd_ovrf,
                                      (void *)x_coordinate_overflow.data(), 0) != cudaSuccess);
  problem = (problem ||
             cudaHostGetDevicePointer((void **)&devc_ycrd_ovrf,
                                      (void *)y_coordinate_overflow.data(), 0) != cudaSuccess);
  problem = (problem ||
             cudaHostGetDevicePointer((void **)&devc_zcrd_ovrf,
                                      (void *)z_coordinate_overflow.data(), 0) != cudaSuccess);
  problem = (problem || cudaHostGetDevicePointer((void **)&devc_xvel,
                                                 (void *)x_velocities.data(), 0) != cudaSuccess);
  problem = (problem || cudaHostGetDevicePointer((void **)&devc_yvel,
                                                 (void *)y_velocities.data(), 0) != cudaSuccess);
  problem = (problem || cudaHostGetDevicePointer((void **)&devc_zvel,
                                                 (void *)z_velocities.data(), 0) != cudaSuccess);
  problem = (problem ||
             cudaHostGetDevicePointer((void **)&devc_xvel_ovrf,
                                      (void *)x_velocity_overflow.data(), 0) != cudaSuccess);
  problem = (problem ||
             cudaHostGetDevicePointer((void **)&devc_yvel_ovrf,
                                      (void *)y_velocity_overflow.data(), 0) != cudaSuccess);
  problem = (problem ||
             cudaHostGetDevicePointer((void **)&devc_zvel_ovrf,
                                      (void *)z_velocity_overflow.data(), 0) != cudaSuccess);
  problem = (problem || cudaHostGetDevicePointer((void **)&devc_xfrc,
                                                 (void *)x_forces.data(), 0) != cudaSuccess);
  problem = (problem || cudaHostGetDevicePointer((void **)&devc_yfrc,
                                                 (void *)y_forces.data(), 0) != cudaSuccess);
  problem = (problem || cudaHostGetDevicePointer((void **)&devc_zfrc,
                                                 (void *)z_forces.data(), 0) != cudaSuccess);
  problem = (problem ||
             cudaHostGetDevicePointer((void **)&devc_xfrc_ovrf,
                                      (void *)x_force_overflow.data(), 0) != cudaSuccess);
  problem = (problem ||
             cudaHostGetDevicePointer((void **)&devc_yfrc_ovrf,
                                      (void *)y_force_overflow.data(), 0) != cudaSuccess);
  problem = (problem ||
             cudaHostGetDevicePointer((void **)&devc_zfrc_ovrf,
                                      (void *)z_force_overflow.data(), 0) != cudaSuccess);
  problem = (problem ||
             cudaHostGetDevicePointer((void **)&devc_xprv,
                                      (void *)x_prior_coordinates.data(), 0) != cudaSuccess);
  problem = (problem ||
             cudaHostGetDevicePointer((void **)&devc_yprv,
                                      (void *)y_prior_coordinates.data(), 0) != cudaSuccess);
  problem = (problem ||
             cudaHostGetDevicePointer((void **)&devc_zprv,
                                      (void *)z_prior_coordinates.data(), 0) != cudaSuccess);
  problem = (problem ||
             cudaHostGetDevicePointer((void **)&devc_xprv_ovrf,
                                      (void *)x_prior_coord_overflow.data(), 0) != cudaSuccess);
  problem = (problem ||
             cudaHostGetDevicePointer((void **)&devc_yprv_ovrf,
                                      (void *)y_prior_coord_overflow.data(), 0) != cudaSuccess);
  problem = (problem ||
             cudaHostGetDevicePointer((void **)&devc_zprv_ovrf,
                                      (void *)z_prior_coord_overflow.data(), 0) != cudaSuccess);
  problem = (problem ||
             cudaHostGetDevicePointer((void **)&devc_xnxt,
                                      (void *)x_prior_coordinates.data(), 0) != cudaSuccess);
  problem = (problem ||
             cudaHostGetDevicePointer((void **)&devc_ynxt,
                                      (void *)y_prior_coordinates.data(), 0) != cudaSuccess);
  problem = (problem ||
             cudaHostGetDevicePointer((void **)&devc_znxt,
                                      (void *)z_prior_coordinates.data(), 0) != cudaSuccess);
  problem = (problem ||
             cudaHostGetDevicePointer((void **)&devc_xnxt_ovrf,
                                      (void *)x_prior_coord_overflow.data(), 0) != cudaSuccess);
  problem = (problem ||
             cudaHostGetDevicePointer((void **)&devc_ynxt_ovrf,
                                      (void *)y_prior_coord_overflow.data(), 0) != cudaSuccess);
  problem = (problem ||
             cudaHostGetDevicePointer((void **)&devc_znxt_ovrf,
                                      (void *)z_prior_coord_overflow.data(), 0) != cudaSuccess);
#  endif
  if (problem) {
    rtErr("Unable to get device-mapped pointers to host memory.  Types of memory used in each "
          "array: " + std::string(int_data.getLabel().name) + " " + int_data.getLabel().format +
          ", " + std::string(llint_data.getLabel().name) + " " + llint_data.getLabel().format +
          ", " + std::string(double_data.getLabel().name) + " " + double_data.getLabel().format +
          ", " + std::string(float_data.getLabel().name) + " " + float_data.getLabel().format +
          ".", "PhaseSpaceSynthesis", "deviceViewToHostData");
  }
  switch (orientation) {
  case CoordinateCycle::PAST:
    return PsSynthesisWriter(system_count, unit_cell, heat_bath_kind, piston_kind, time_step,
                             devc_atom_starts, devc_atom_counts, globalpos_scale,
                             localpos_scale, velocity_scale, force_scale, globalpos_scale_bits,
                             localpos_scale_bits, velocity_scale_bits, force_scale_bits,
                             devc_boxvecs, devc_umat, devc_invu, devc_boxdims, devc_sp_umat,
                             devc_sp_invu, devc_sp_boxdims, devc_xprv, devc_yprv, devc_zprv,
                             devc_xprv_ovrf, devc_yprv_ovrf, devc_zprv_ovrf, devc_xvel, devc_yvel,
                             devc_zvel, devc_xvel_ovrf, devc_yvel_ovrf, devc_zvel_ovrf, devc_xfrc,
                             devc_yfrc, devc_zfrc, devc_xfrc_ovrf, devc_yfrc_ovrf, devc_zfrc_ovrf,
                             devc_xnxt, devc_ynxt, devc_znxt, devc_xnxt_ovrf, devc_ynxt_ovrf,
                             devc_znxt_ovrf, devc_xcrd, devc_ycrd, devc_zcrd, devc_xcrd_ovrf,
                             devc_ycrd_ovrf, devc_zcrd_ovrf);
  case CoordinateCycle::PRESENT:
    return PsSynthesisWriter(system_count, unit_cell, heat_bath_kind, piston_kind, time_step,
                             devc_atom_starts, devc_atom_counts, globalpos_scale,
                             localpos_scale, velocity_scale, force_scale, globalpos_scale_bits,
                             localpos_scale_bits, velocity_scale_bits, force_scale_bits,
                             devc_boxvecs, devc_umat, devc_invu, devc_boxdims, devc_sp_umat,
                             devc_sp_invu, devc_sp_boxdims, devc_xcrd, devc_ycrd, devc_zcrd,
                             devc_xcrd_ovrf, devc_ycrd_ovrf, devc_zcrd_ovrf, devc_xvel, devc_yvel,
                             devc_zvel, devc_xvel_ovrf, devc_yvel_ovrf, devc_zvel_ovrf, devc_xfrc,
                             devc_yfrc, devc_zfrc, devc_xfrc_ovrf, devc_yfrc_ovrf, devc_zfrc_ovrf,
                             devc_xprv, devc_yprv, devc_zprv, devc_xprv_ovrf, devc_yprv_ovrf,
                             devc_zprv_ovrf, devc_xnxt, devc_ynxt, devc_znxt, devc_xnxt_ovrf,
                             devc_ynxt_ovrf, devc_znxt_ovrf);
  case CoordinateCycle::FUTURE:
    return PsSynthesisWriter(system_count, unit_cell, heat_bath_kind, piston_kind, time_step,
                             devc_atom_starts, devc_atom_counts, globalpos_scale,
                             localpos_scale, velocity_scale, force_scale, globalpos_scale_bits,
                             localpos_scale_bits, velocity_scale_bits, force_scale_bits,
                             devc_boxvecs, devc_umat, devc_invu, devc_boxdims, devc_sp_umat,
                             devc_sp_invu, devc_sp_boxdims, devc_xnxt, devc_ynxt, devc_znxt,
                             devc_xnxt_ovrf, devc_ynxt_ovrf, devc_znxt_ovrf, devc_xvel, devc_yvel,
                             devc_zvel, devc_xvel_ovrf, devc_yvel_ovrf, devc_zvel_ovrf, devc_xfrc,
                             devc_yfrc, devc_zfrc, devc_xfrc_ovrf, devc_yfrc_ovrf, devc_zfrc_ovrf,
                             devc_xcrd, devc_ycrd, devc_zcrd, devc_xcrd_ovrf, devc_ycrd_ovrf,
                             devc_zcrd_ovrf, devc_xprv, devc_yprv, devc_zprv, devc_xprv_ovrf,
                             devc_yprv_ovrf, devc_zprv_ovrf);
  }
}

//-------------------------------------------------------------------------------------------------
void PhaseSpaceSynthesis::upload() {
  atom_starts.upload();
  atom_counts.upload();
  llint_data.upload();
  int_data.upload();
  double_data.upload();
  float_data.upload();
}

//-------------------------------------------------------------------------------------------------
void PhaseSpaceSynthesis::upload(const TrajectoryKind kind) {
  const int stride = roundUp(atom_starts.readHost(system_count - 1) +
                             atom_counts.readHost(system_count - 1), warp_size_int);
  switch (kind) {
  case TrajectoryKind::POSITIONS:
    x_coordinates.upload(0, 9 * stride);
    if (globalpos_scale_bits > globalpos_scale_nonoverflow_bits) {
      x_coordinate_overflow.upload(0, 9 * stride);
    }
    double_data.upload();
    float_data.upload();
    box_vectors.upload();
    break;
  case TrajectoryKind::VELOCITIES:
    x_velocities.upload(0, 3 * stride);
    if (velocity_scale_bits > velocity_scale_nonoverflow_bits) {
      x_velocity_overflow.upload(0, 3 * stride);
    }
    break;
  case TrajectoryKind::FORCES:
    x_forces.upload(0, 3 * stride);
    if (force_scale_bits > force_scale_nonoverflow_bits) {
      x_force_overflow.upload(0, 3 * stride);
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void PhaseSpaceSynthesis::upload(const TrajectoryKind kind, const int system_index,
                                 const GpuDetails &gpu) {
  upload(kind, system_index, system_index + 1, gpu);
}

//-------------------------------------------------------------------------------------------------
void PhaseSpaceSynthesis::upload(const TrajectoryKind kind, const int system_lower_bound,
                                 const int system_upper_bound, const GpuDetails &gpu) {
  PsSynthesisWriter devc_view = data(HybridTargetLevel::DEVICE);
  PsSynthesisWriter host_view = deviceViewToHostData();
  systemTransfer(&devc_view, &host_view, kind, system_lower_bound, system_upper_bound, gpu);
}

//-------------------------------------------------------------------------------------------------
void PhaseSpaceSynthesis::download() {
  atom_starts.download();
  atom_counts.download();
  llint_data.download();
  int_data.download();
  double_data.download();
  float_data.download();
}

//-------------------------------------------------------------------------------------------------
void PhaseSpaceSynthesis::download(const TrajectoryKind kind) {
  const int stride = roundUp(atom_starts.readHost(system_count - 1) +
                             atom_counts.readHost(system_count - 1), warp_size_int);
  switch (kind) {
  case TrajectoryKind::POSITIONS:
    x_coordinates.download(0, 9 * stride);
    if (globalpos_scale_bits > globalpos_scale_nonoverflow_bits) {
      x_coordinate_overflow.download(0, 9 * stride);
    }
    double_data.download();
    float_data.download();
    box_vectors.download();
    break;
  case TrajectoryKind::VELOCITIES:
    x_velocities.download(0, 3 * stride);
    if (velocity_scale_bits > velocity_scale_nonoverflow_bits) {
      x_velocity_overflow.download(0, 3 * stride);
    }
    break;
  case TrajectoryKind::FORCES:
    x_forces.download(0, 3 * stride);
    if (force_scale_bits > force_scale_nonoverflow_bits) {
      x_force_overflow.download(0, 3 * stride);
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void PhaseSpaceSynthesis::download(const TrajectoryKind kind, const int system_index,
                                   const GpuDetails &gpu) {
  download(kind, system_index, system_index + 1, gpu);
}

//-------------------------------------------------------------------------------------------------
void PhaseSpaceSynthesis::download(const TrajectoryKind kind, const int system_lower_bound,
                                   const int system_upper_bound, const GpuDetails &gpu) {
  PsSynthesisWriter devc_view = data(HybridTargetLevel::DEVICE);
  PsSynthesisWriter host_view = deviceViewToHostData();
  systemTransfer(&host_view, &devc_view, kind, system_lower_bound, system_upper_bound, gpu);
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
void PhaseSpaceSynthesis::extractSystem(PhaseSpace *ps, const int index,
                                        const HybridTargetLevel tier) const {
  PhaseSpaceWriter psw = ps->data();
  if (atom_counts.readHost(index) != psw.natom) {
    rtErr("A PhaseSpace object sized for " + std::to_string(psw.natom) + " atoms is not prepared "
          "to accept a system of " + std::to_string(atom_counts.readHost(index)) + " atoms from "
          "this synthesis.", "PhaseSpaceSynthesis", "extractPhaseSpace");
  }  
#ifdef STORMM_USE_HPC
  std::vector<llint> xcrd_buff, ycrd_buff, zcrd_buff, xvel_buff, yvel_buff, zvel_buff;
  std::vector<llint> xfrc_buff, yfrc_buff, zfrc_buff, box_buff, inv_buff, dim_buff;
#endif
  int atom_offset = atom_starts.readHost(index);
  int mtrx_offset = index * roundUp(9, warp_size_int);
  int bdim_offset = index * roundUp(6, warp_size_int);
  const double crd_deflation = inverse_globalpos_scale;
  const double vel_deflation = inverse_velocity_scale;
  const double frc_deflation = inverse_force_scale;
  const double crd_ovrf_deflation = crd_deflation * max_llint_accumulation;
  const double vel_ovrf_deflation = vel_deflation * max_llint_accumulation;
  const double frc_ovrf_deflation = frc_deflation * max_llint_accumulation;
  switch (tier) {
  case HybridTargetLevel::HOST:
    {
      // The phase space object should be extracted in a manner that orients its "present" towards
      // the current coordinates from the phase space synthesis, wherever the larger synthesis
      // sits in its own cycle.
      const PsSynthesisReader tpsr = this->data(cycle_position, tier);
      for (int i = 0; i < psw.natom; i++) {
        psw.xcrd[i] = static_cast<double>(tpsr.xcrd[atom_offset + i]) * crd_deflation;
        psw.ycrd[i] = static_cast<double>(tpsr.ycrd[atom_offset + i]) * crd_deflation;
        psw.zcrd[i] = static_cast<double>(tpsr.zcrd[atom_offset + i]) * crd_deflation;
        psw.xvel[i] = static_cast<double>(tpsr.xvel[atom_offset + i]) * vel_deflation;
        psw.yvel[i] = static_cast<double>(tpsr.yvel[atom_offset + i]) * vel_deflation;
        psw.zvel[i] = static_cast<double>(tpsr.zvel[atom_offset + i]) * vel_deflation;
        psw.xfrc[i] = static_cast<double>(tpsr.xfrc[atom_offset + i]) * frc_deflation;
        psw.yfrc[i] = static_cast<double>(tpsr.yfrc[atom_offset + i]) * frc_deflation;
        psw.zfrc[i] = static_cast<double>(tpsr.zfrc[atom_offset + i]) * frc_deflation;
        psw.xprv[i] = static_cast<double>(tpsr.xprv[atom_offset + i]) * crd_deflation;
        psw.yprv[i] = static_cast<double>(tpsr.yprv[atom_offset + i]) * crd_deflation;
        psw.zprv[i] = static_cast<double>(tpsr.zprv[atom_offset + i]) * crd_deflation;
        psw.xnxt[i] = static_cast<double>(tpsr.xnxt[atom_offset + i]) * crd_deflation;
        psw.ynxt[i] = static_cast<double>(tpsr.ynxt[atom_offset + i]) * crd_deflation;
        psw.znxt[i] = static_cast<double>(tpsr.znxt[atom_offset + i]) * crd_deflation;
      }
      if (globalpos_scale_bits > globalpos_scale_nonoverflow_bits) {
        for (int i = 0; i < psw.natom; i++) {
          psw.xcrd[i] += static_cast<double>(tpsr.xcrd_ovrf[atom_offset + i]) * crd_ovrf_deflation;
          psw.ycrd[i] += static_cast<double>(tpsr.ycrd_ovrf[atom_offset + i]) * crd_ovrf_deflation;
          psw.zcrd[i] += static_cast<double>(tpsr.zcrd_ovrf[atom_offset + i]) * crd_ovrf_deflation;
          psw.xprv[i] += static_cast<double>(tpsr.xprv_ovrf[atom_offset + i]) * crd_ovrf_deflation;
          psw.yprv[i] += static_cast<double>(tpsr.yprv_ovrf[atom_offset + i]) * crd_ovrf_deflation;
          psw.zprv[i] += static_cast<double>(tpsr.zprv_ovrf[atom_offset + i]) * crd_ovrf_deflation;
          psw.xnxt[i] += static_cast<double>(tpsr.xnxt_ovrf[atom_offset + i]) * crd_ovrf_deflation;
          psw.ynxt[i] += static_cast<double>(tpsr.ynxt_ovrf[atom_offset + i]) * crd_ovrf_deflation;
          psw.znxt[i] += static_cast<double>(tpsr.znxt_ovrf[atom_offset + i]) * crd_ovrf_deflation;
        }
      }
      if (velocity_scale_bits > velocity_scale_nonoverflow_bits) {
        for (int i = 0; i < psw.natom; i++) {
          psw.xvel[i] += static_cast<double>(tpsr.xvel_ovrf[atom_offset + i]) * vel_ovrf_deflation;
          psw.yvel[i] += static_cast<double>(tpsr.yvel_ovrf[atom_offset + i]) * vel_ovrf_deflation;
          psw.zvel[i] += static_cast<double>(tpsr.zvel_ovrf[atom_offset + i]) * vel_ovrf_deflation;
        }
      }
      if (force_scale_bits > force_scale_nonoverflow_bits) {
        const int* xfrc_ovrf_ptr = x_force_overflow.data();
        const int* yfrc_ovrf_ptr = y_force_overflow.data();
        const int* zfrc_ovrf_ptr = z_force_overflow.data();
        for (int i = 0; i < psw.natom; i++) {
          psw.xfrc[i] += static_cast<double>(tpsr.xfrc_ovrf[atom_offset + i]) * frc_ovrf_deflation;
          psw.yfrc[i] += static_cast<double>(tpsr.yfrc_ovrf[atom_offset + i]) * frc_ovrf_deflation;
          psw.zfrc[i] += static_cast<double>(tpsr.zfrc_ovrf[atom_offset + i]) * frc_ovrf_deflation;
        }
      }
      for (int i = 0; i < 9; i++) {
        psw.umat[i] = tpsr.umat[mtrx_offset + i];
        psw.invu[i] = tpsr.invu[mtrx_offset + i];
      }
      for (int i = 0; i < 6; i++) {
        psw.boxdim[i] = tpsr.boxdims[bdim_offset + i];
      }
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    {
      const std::vector<llint> xcrd_buff = x_coordinates.readDevice(atom_offset, psw.natom);
      const std::vector<llint> ycrd_buff = y_coordinates.readDevice(atom_offset, psw.natom);
      const std::vector<llint> zcrd_buff = z_coordinates.readDevice(atom_offset, psw.natom);
      const std::vector<llint> xvel_buff = x_velocities.readDevice(atom_offset, psw.natom);
      const std::vector<llint> yvel_buff = y_velocities.readDevice(atom_offset, psw.natom);
      const std::vector<llint> zvel_buff = z_velocities.readDevice(atom_offset, psw.natom);
      const std::vector<llint> xfrc_buff = x_forces.readDevice(atom_offset, psw.natom);
      const std::vector<llint> yfrc_buff = y_forces.readDevice(atom_offset, psw.natom);
      const std::vector<llint> zfrc_buff = z_forces.readDevice(atom_offset, psw.natom);
      const std::vector<llint> xprv_buff = x_prior_coordinates.readDevice(atom_offset, psw.natom);
      const std::vector<llint> yprv_buff = y_prior_coordinates.readDevice(atom_offset, psw.natom);
      const std::vector<llint> zprv_buff = z_prior_coordinates.readDevice(atom_offset, psw.natom);
      const std::vector<llint> xnxt_buff = x_future_coordinates.readDevice(atom_offset, psw.natom);
      const std::vector<llint> ynxt_buff = y_future_coordinates.readDevice(atom_offset, psw.natom);
      const std::vector<llint> znxt_buff = z_future_coordinates.readDevice(atom_offset, psw.natom);
      const std::vector<double> box_buff = box_space_transforms.readDevice(mtrx_offset, 9);
      const std::vector<double> inv_buff = inverse_transforms.readDevice(mtrx_offset, 9);
      const std::vector<double> dim_buff = box_dimensions.readDevice(bdim_offset, 9);
      for (int i = 0; i < psw.natom; i++) {
        psw.xcrd[i] = static_cast<double>(xcrd_buff[i]) * crd_deflation;
        psw.ycrd[i] = static_cast<double>(ycrd_buff[i]) * crd_deflation;
        psw.zcrd[i] = static_cast<double>(zcrd_buff[i]) * crd_deflation;
        psw.xvel[i] = static_cast<double>(xvel_buff[i]) * vel_deflation;
        psw.yvel[i] = static_cast<double>(yvel_buff[i]) * vel_deflation;
        psw.zvel[i] = static_cast<double>(zvel_buff[i]) * vel_deflation;
        psw.xfrc[i] = static_cast<double>(xfrc_buff[i]) * frc_deflation;
        psw.yfrc[i] = static_cast<double>(yfrc_buff[i]) * frc_deflation;
        psw.zfrc[i] = static_cast<double>(zfrc_buff[i]) * frc_deflation;
        psw.xprv[i] = static_cast<double>(xprv_buff[i]) * crd_deflation;
        psw.yprv[i] = static_cast<double>(yprv_buff[i]) * crd_deflation;
        psw.zprv[i] = static_cast<double>(zprv_buff[i]) * crd_deflation;
        psw.xnxt[i] = static_cast<double>(xnxt_buff[i]) * crd_deflation;
        psw.ynxt[i] = static_cast<double>(ynxt_buff[i]) * crd_deflation;
        psw.znxt[i] = static_cast<double>(znxt_buff[i]) * crd_deflation;
      }
      if (globalpos_scale_bits > globalpos_scale_nonoverflow_bits) {
        const std::vector<int> xcrd_ovrf_buff = x_coordinate_overflow.readDevice(atom_offset,
                                                                                 psw.natom);
        const std::vector<int> ycrd_ovrf_buff = y_coordinate_overflow.readDevice(atom_offset,
                                                                                 psw.natom);
        const std::vector<int> zcrd_ovrf_buff = z_coordinate_overflow.readDevice(atom_offset,
                                                                                 psw.natom);
        const std::vector<int> xprv_ovrf_buff = x_prior_coord_overflow.readDevice(atom_offset,
                                                                                  psw.natom);
        const std::vector<int> yprv_ovrf_buff = y_prior_coord_overflow.readDevice(atom_offset,
                                                                                  psw.natom);
        const std::vector<int> zprv_ovrf_buff = z_prior_coord_overflow.readDevice(atom_offset,
                                                                                  psw.natom);
        for (int i = 0; i < psw.natom; i++) {
          psw.xcrd[i] += static_cast<double>(xcrd_ovrf_buff[i]) * crd_ovrf_deflation;
          psw.ycrd[i] += static_cast<double>(ycrd_ovrf_buff[i]) * crd_ovrf_deflation;
          psw.zcrd[i] += static_cast<double>(zcrd_ovrf_buff[i]) * crd_ovrf_deflation;
          psw.xprv[i] += static_cast<double>(xprv_ovrf_buff[i]) * crd_ovrf_deflation;
          psw.yprv[i] += static_cast<double>(yprv_ovrf_buff[i]) * crd_ovrf_deflation;
          psw.zprv[i] += static_cast<double>(zprv_ovrf_buff[i]) * crd_ovrf_deflation;
        }
      }
      if (velocity_scale_bits > velocity_scale_nonoverflow_bits) {
        const std::vector<int> xvel_ovrf_buff = x_velocity_overflow.readDevice(atom_offset,
                                                                               psw.natom);
        const std::vector<int> yvel_ovrf_buff = y_velocity_overflow.readDevice(atom_offset,
                                                                               psw.natom);
        const std::vector<int> zvel_ovrf_buff = z_velocity_overflow.readDevice(atom_offset,
                                                                               psw.natom);
        for (int i = 0; i < psw.natom; i++) {
          psw.xvel[i] += static_cast<double>(xvel_ovrf_buff[i]) * vel_ovrf_deflation;
          psw.yvel[i] += static_cast<double>(yvel_ovrf_buff[i]) * vel_ovrf_deflation;
          psw.zvel[i] += static_cast<double>(zvel_ovrf_buff[i]) * vel_ovrf_deflation;
        }
      }
      if (force_scale_bits > force_scale_nonoverflow_bits) {
        const std::vector<int> xfrc_ovrf_buff = x_force_overflow.readDevice(atom_offset,
                                                                            psw.natom);
        const std::vector<int> yfrc_ovrf_buff = y_force_overflow.readDevice(atom_offset,
                                                                            psw.natom);
        const std::vector<int> zfrc_ovrf_buff = z_force_overflow.readDevice(atom_offset,
                                                                            psw.natom);
        for (int i = 0; i < psw.natom; i++) {
          psw.xfrc[i] += static_cast<double>(xfrc_ovrf_buff[i]) * frc_ovrf_deflation;
          psw.yfrc[i] += static_cast<double>(yfrc_ovrf_buff[i]) * frc_ovrf_deflation;
          psw.zfrc[i] += static_cast<double>(zfrc_ovrf_buff[i]) * frc_ovrf_deflation;
        }
      }
      for (int i = 0; i < 9; i++) {
        psw.umat[i] = box_buff[i];
        psw.invu[i] = inv_buff[i];
      }
      for (int i = 0; i < 6; i++) {
        psw.boxdim[i] = dim_buff[i];
      }
    }
    break;
#endif
  }
}

//-------------------------------------------------------------------------------------------------
void PhaseSpaceSynthesis::extractCoordinates(PhaseSpace *ps, const int index,
                                             const TrajectoryKind trajkind,
                                             const HybridTargetLevel tier) const {
  PhaseSpaceWriter psw = ps->data();
  if (atom_counts.readHost(index) != psw.natom) {
    rtErr("A PhaseSpace object sized for " + std::to_string(psw.natom) + " atoms is not prepared "
          "to accept a system of " + std::to_string(atom_counts.readHost(index)) + " atoms from "
          "this synthesis.", "PhaseSpaceSynthesis", "extractPhaseSpace");
  }
  const int atom_offset = atom_starts.readHost(index);
  const int mtrx_offset = index * roundUp(9, warp_size_int);
  const int bdim_offset = index * roundUp(6, warp_size_int);

  // Download the relevant data, if necessary.
  const PsSynthesisReader tpsr = this->data(cycle_position, tier);
  switch (tier) {
  case HybridTargetLevel::HOST:
    switch (trajkind) {
    case TrajectoryKind::POSITIONS:
      {
        // The phase space object should be extracted in a manner that orients its "present"
        // towards the current coordinates from the phase space synthesis, wherever the larger
        // synthesis sits in its own cycle.
        const double crd_deflation = inverse_globalpos_scale;
        for (int i = 0; i < psw.natom; i++) {
          psw.xcrd[i] = static_cast<double>(tpsr.xcrd[atom_offset + i]) * crd_deflation;
          psw.ycrd[i] = static_cast<double>(tpsr.ycrd[atom_offset + i]) * crd_deflation;
          psw.zcrd[i] = static_cast<double>(tpsr.zcrd[atom_offset + i]) * crd_deflation;
          psw.xprv[i] = static_cast<double>(tpsr.xprv[atom_offset + i]) * crd_deflation;
          psw.yprv[i] = static_cast<double>(tpsr.yprv[atom_offset + i]) * crd_deflation;
          psw.zprv[i] = static_cast<double>(tpsr.zprv[atom_offset + i]) * crd_deflation;
        }
        if (globalpos_scale_bits > globalpos_scale_nonoverflow_bits) {
          const double crd_ovrf_deflation = crd_deflation * max_llint_accumulation;
          for (int i = 0; i < psw.natom; i++) {
            const size_t aoi_zu = atom_offset + i; 
            psw.xcrd[i] += static_cast<double>(tpsr.xcrd_ovrf[aoi_zu]) * crd_ovrf_deflation;
            psw.ycrd[i] += static_cast<double>(tpsr.ycrd_ovrf[aoi_zu]) * crd_ovrf_deflation;
            psw.zcrd[i] += static_cast<double>(tpsr.zcrd_ovrf[aoi_zu]) * crd_ovrf_deflation;
            psw.xprv[i] += static_cast<double>(tpsr.xprv_ovrf[aoi_zu]) * crd_ovrf_deflation;
            psw.yprv[i] += static_cast<double>(tpsr.yprv_ovrf[aoi_zu]) * crd_ovrf_deflation;
            psw.zprv[i] += static_cast<double>(tpsr.zprv_ovrf[aoi_zu]) * crd_ovrf_deflation;
            psw.xnxt[i] += static_cast<double>(tpsr.xnxt_ovrf[aoi_zu]) * crd_ovrf_deflation;
            psw.ynxt[i] += static_cast<double>(tpsr.ynxt_ovrf[aoi_zu]) * crd_ovrf_deflation;
            psw.znxt[i] += static_cast<double>(tpsr.znxt_ovrf[aoi_zu]) * crd_ovrf_deflation;
          }
        }
        for (int i = 0; i < 9; i++) {
          psw.umat[i] = tpsr.umat[mtrx_offset + i];
          psw.invu[i] = tpsr.invu[mtrx_offset + i];
        }
        for (int i = 0; i < 6; i++) {
          psw.boxdim[i] = tpsr.boxdims[bdim_offset + i];        
        }
      }
      break;
    case TrajectoryKind::VELOCITIES:
      {
        const double vel_deflation = inverse_velocity_scale;
        for (int i = 0; i < psw.natom; i++) {
          psw.xvel[i] = static_cast<double>(tpsr.xvel[atom_offset + i]) * vel_deflation;
          psw.yvel[i] = static_cast<double>(tpsr.yvel[atom_offset + i]) * vel_deflation;
          psw.zvel[i] = static_cast<double>(tpsr.zvel[atom_offset + i]) * vel_deflation;
        }
        if (velocity_scale_bits > velocity_scale_nonoverflow_bits) {
          const double vel_ovrf_deflation = vel_deflation * max_llint_accumulation;
          for (int i = 0; i < psw.natom; i++) {
            const size_t aoi_zu = atom_offset + i;
            psw.xvel[i] += static_cast<double>(tpsr.xvel_ovrf[aoi_zu]) * vel_ovrf_deflation;
            psw.yvel[i] += static_cast<double>(tpsr.yvel_ovrf[aoi_zu]) * vel_ovrf_deflation;
            psw.zvel[i] += static_cast<double>(tpsr.zvel_ovrf[aoi_zu]) * vel_ovrf_deflation;
          }
        }
      }
      break;
    case TrajectoryKind::FORCES:
      {
        const double frc_deflation = inverse_force_scale;
        for (int i = 0; i < psw.natom; i++) {
          psw.xfrc[i] = static_cast<double>(tpsr.xfrc[atom_offset + i]) * frc_deflation;
          psw.yfrc[i] = static_cast<double>(tpsr.yfrc[atom_offset + i]) * frc_deflation;
          psw.zfrc[i] = static_cast<double>(tpsr.zfrc[atom_offset + i]) * frc_deflation;
        }
        if (force_scale_bits > force_scale_nonoverflow_bits) {
          const double frc_ovrf_deflation = frc_deflation * max_llint_accumulation;
          for (int i = 0; i < psw.natom; i++) {
            const size_t aoi_zu = atom_offset + i;
            psw.xfrc[i] += static_cast<double>(tpsr.xfrc_ovrf[aoi_zu]) * frc_ovrf_deflation;
            psw.yfrc[i] += static_cast<double>(tpsr.yfrc_ovrf[aoi_zu]) * frc_ovrf_deflation;
            psw.zfrc[i] += static_cast<double>(tpsr.zfrc_ovrf[aoi_zu]) * frc_ovrf_deflation;
          }
        }
      }
      break;
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    switch (trajkind) {
    case TrajectoryKind::POSITIONS:
      {
        const std::vector<llint> xcrd_buff = x_coordinates.readDevice(atom_offset, psw.natom);
        const std::vector<llint> ycrd_buff = y_coordinates.readDevice(atom_offset, psw.natom);
        const std::vector<llint> zcrd_buff = z_coordinates.readDevice(atom_offset, psw.natom);
        const std::vector<llint> xprv_buff = x_prior_coordinates.readDevice(atom_offset,
                                                                            psw.natom);
        const std::vector<llint> yprv_buff = y_prior_coordinates.readDevice(atom_offset,
                                                                            psw.natom);
        const std::vector<llint> zprv_buff = z_prior_coordinates.readDevice(atom_offset,
                                                                            psw.natom);
        const std::vector<llint> xnxt_buff = x_future_coordinates.readDevice(atom_offset,
                                                                             psw.natom);
        const std::vector<llint> ynxt_buff = y_future_coordinates.readDevice(atom_offset,
                                                                             psw.natom);
        const std::vector<llint> znxt_buff = z_future_coordinates.readDevice(atom_offset,
                                                                             psw.natom);
        const std::vector<double> box_buff = box_space_transforms.readDevice(mtrx_offset, 9);
        const std::vector<double> inv_buff = inverse_transforms.readDevice(mtrx_offset, 9);
        const std::vector<double> dim_buff = box_dimensions.readDevice(bdim_offset, 6);
        const double crd_deflation = inverse_globalpos_scale;
        for (int i = 0; i < psw.natom; i++) {
          psw.xcrd[i] = static_cast<double>(xcrd_buff[i]) * crd_deflation;
          psw.ycrd[i] = static_cast<double>(ycrd_buff[i]) * crd_deflation;
          psw.zcrd[i] = static_cast<double>(zcrd_buff[i]) * crd_deflation;
          psw.xprv[i] = static_cast<double>(xprv_buff[i]) * crd_deflation;
          psw.yprv[i] = static_cast<double>(yprv_buff[i]) * crd_deflation;
          psw.zprv[i] = static_cast<double>(zprv_buff[i]) * crd_deflation;
          psw.xnxt[i] = static_cast<double>(xnxt_buff[i]) * crd_deflation;
          psw.ynxt[i] = static_cast<double>(ynxt_buff[i]) * crd_deflation;
          psw.znxt[i] = static_cast<double>(znxt_buff[i]) * crd_deflation;
        }
        if (globalpos_scale_bits > globalpos_scale_nonoverflow_bits) {
          const std::vector<int> xcrd_ovrf_buff = x_coordinate_overflow.readDevice(atom_offset,
                                                                                   psw.natom);
          const std::vector<int> ycrd_ovrf_buff = y_coordinate_overflow.readDevice(atom_offset,
                                                                                   psw.natom);
          const std::vector<int> zcrd_ovrf_buff = z_coordinate_overflow.readDevice(atom_offset,
                                                                                   psw.natom);
          const std::vector<int> xprv_ovrf_buff = x_prior_coord_overflow.readDevice(atom_offset,
                                                                                    psw.natom);
          const std::vector<int> yprv_ovrf_buff = y_prior_coord_overflow.readDevice(atom_offset,
                                                                                    psw.natom);
          const std::vector<int> zprv_ovrf_buff = z_prior_coord_overflow.readDevice(atom_offset,
                                                                                    psw.natom);
          const std::vector<int> xnxt_ovrf_buff = x_future_coord_overflow.readDevice(atom_offset,
                                                                                     psw.natom);
          const std::vector<int> ynxt_ovrf_buff = y_future_coord_overflow.readDevice(atom_offset,
                                                                                     psw.natom);
          const std::vector<int> znxt_ovrf_buff = z_future_coord_overflow.readDevice(atom_offset,
                                                                                     psw.natom);
          const double crd_ovrf_deflation = crd_deflation * max_llint_accumulation;
          for (int i = 0; i < psw.natom; i++) {
            psw.xcrd[i] += static_cast<double>(xcrd_ovrf_buff[i]) * crd_ovrf_deflation;
            psw.ycrd[i] += static_cast<double>(ycrd_ovrf_buff[i]) * crd_ovrf_deflation;
            psw.zcrd[i] += static_cast<double>(zcrd_ovrf_buff[i]) * crd_ovrf_deflation;
            psw.xprv[i] += static_cast<double>(xprv_ovrf_buff[i]) * crd_ovrf_deflation;
            psw.yprv[i] += static_cast<double>(yprv_ovrf_buff[i]) * crd_ovrf_deflation;
            psw.zprv[i] += static_cast<double>(zprv_ovrf_buff[i]) * crd_ovrf_deflation;
            psw.xnxt[i] += static_cast<double>(xnxt_ovrf_buff[i]) * crd_ovrf_deflation;
            psw.ynxt[i] += static_cast<double>(ynxt_ovrf_buff[i]) * crd_ovrf_deflation;
            psw.znxt[i] += static_cast<double>(znxt_ovrf_buff[i]) * crd_ovrf_deflation;
          }
        }
        for (int i = 0; i < 9; i++) {
          psw.umat[i] = box_buff[i];
          psw.invu[i] = inv_buff[i];
        }
        for (int i = 0; i < 6; i++) {
          psw.boxdim[i] = dim_buff[i];
        }
      }
      break;
    case TrajectoryKind::VELOCITIES:
      {
        const std::vector<llint> xvel_buff = x_velocities.readDevice(atom_offset, psw.natom);
        const std::vector<llint> yvel_buff = y_velocities.readDevice(atom_offset, psw.natom);
        const std::vector<llint> zvel_buff = z_velocities.readDevice(atom_offset, psw.natom);
        const double vel_deflation = inverse_velocity_scale;
        for (int i = 0; i < psw.natom; i++) {
          psw.xvel[i] = static_cast<double>(xvel_buff[i]) * vel_deflation;
          psw.yvel[i] = static_cast<double>(yvel_buff[i]) * vel_deflation;
          psw.zvel[i] = static_cast<double>(zvel_buff[i]) * vel_deflation;
        }
        if (velocity_scale_bits > velocity_scale_nonoverflow_bits) {
          const std::vector<int> xvel_ovrf_buff = x_velocity_overflow.readDevice(atom_offset,
                                                                                 psw.natom);
          const std::vector<int> yvel_ovrf_buff = y_velocity_overflow.readDevice(atom_offset,
                                                                                 psw.natom);
          const std::vector<int> zvel_ovrf_buff = z_velocity_overflow.readDevice(atom_offset,
                                                                                 psw.natom);
          const double vel_ovrf_deflation = vel_deflation * max_llint_accumulation;
          for (int i = 0; i < psw.natom; i++) {
            psw.xvel[i] += static_cast<double>(xvel_ovrf_buff[i]) * vel_ovrf_deflation;
            psw.yvel[i] += static_cast<double>(yvel_ovrf_buff[i]) * vel_ovrf_deflation;
            psw.zvel[i] += static_cast<double>(zvel_ovrf_buff[i]) * vel_ovrf_deflation;
          }
        }
      }
      break;
    case TrajectoryKind::FORCES:
      {
        const std::vector<llint> xfrc_buff = x_forces.readDevice(atom_offset, psw.natom);
        const std::vector<llint> yfrc_buff = y_forces.readDevice(atom_offset, psw.natom);
        const std::vector<llint> zfrc_buff = z_forces.readDevice(atom_offset, psw.natom);
        const double frc_deflation = inverse_force_scale;
        for (int i = 0; i < psw.natom; i++) {
          psw.xfrc[i] = static_cast<double>(xfrc_buff[i]) * frc_deflation;
          psw.yfrc[i] = static_cast<double>(yfrc_buff[i]) * frc_deflation;
          psw.zfrc[i] = static_cast<double>(zfrc_buff[i]) * frc_deflation;
        }
        if (force_scale_bits > force_scale_nonoverflow_bits) {
          const std::vector<int> xfrc_ovrf_buff = x_force_overflow.readDevice(atom_offset,
                                                                              psw.natom);
          const std::vector<int> yfrc_ovrf_buff = y_force_overflow.readDevice(atom_offset,
                                                                              psw.natom);
          const std::vector<int> zfrc_ovrf_buff = z_force_overflow.readDevice(atom_offset,
                                                                              psw.natom);
          const double frc_ovrf_deflation = frc_deflation * max_llint_accumulation;
          for (int i = 0; i < psw.natom; i++) {
            psw.xfrc[i] += static_cast<double>(xfrc_ovrf_buff[i]) * frc_ovrf_deflation;
            psw.yfrc[i] += static_cast<double>(yfrc_ovrf_buff[i]) * frc_ovrf_deflation;
            psw.zfrc[i] += static_cast<double>(zfrc_ovrf_buff[i]) * frc_ovrf_deflation;
          }
        }
      }
      break;
    }
    break;
#endif
  }
}

//-------------------------------------------------------------------------------------------------
PhaseSpace PhaseSpaceSynthesis::exportSystem(const int index, const HybridTargetLevel tier) const {
  if (index < 0 || index >= system_count) {
    rtErr("Index " + std::to_string(index) + " is invalid for a collection of " +
          std::to_string(system_count) + " systems.", "PhaseSpaceSynthesis", "exportSystem");
  }
  PhaseSpace result(atom_counts.readHost(index), unit_cell);
  extractSystem(&result, index, tier);
  return result;
}
  
//-------------------------------------------------------------------------------------------------
CoordinateFrame PhaseSpaceSynthesis::exportCoordinates(const int index,
                                                       const TrajectoryKind trajkind,
                                                       const HybridTargetLevel tier) const {
  if (index < 0 || index >= system_count) {
    rtErr("Index " + std::to_string(index) + " is invalid for a collection of " +
          std::to_string(system_count) + " systems.", "PhaseSpaceSynthesis", "exportSystem");
  }
  CoordinateFrame result(atom_counts.readHost(index));
  const int astart = atom_starts.readHost(index);
  CoordinateFrameWriter rsw = result.data();
  std::vector<llint> xbuffer, ybuffer, zbuffer;
  switch (tier) {
  case HybridTargetLevel::HOST:
    switch (trajkind) {
    case TrajectoryKind::POSITIONS:
      switch (cycle_position) {
      case CoordinateCycle::PAST:
        xbuffer = x_prior_coordinates.readHost(astart, rsw.natom);
        ybuffer = y_prior_coordinates.readHost(astart, rsw.natom);
        zbuffer = z_prior_coordinates.readHost(astart, rsw.natom);
        break;
      case CoordinateCycle::PRESENT:
        xbuffer = x_coordinates.readHost(astart, rsw.natom);
        ybuffer = y_coordinates.readHost(astart, rsw.natom);
        zbuffer = z_coordinates.readHost(astart, rsw.natom);
        break;
      case CoordinateCycle::FUTURE:
        xbuffer = x_future_coordinates.readHost(astart, rsw.natom);
        ybuffer = y_future_coordinates.readHost(astart, rsw.natom);
        zbuffer = z_future_coordinates.readHost(astart, rsw.natom);
        break;
      }
      break;
    case TrajectoryKind::VELOCITIES:
      xbuffer = x_velocities.readHost(astart, rsw.natom);
      ybuffer = y_velocities.readHost(astart, rsw.natom);
      zbuffer = z_velocities.readHost(astart, rsw.natom);
      break;
    case TrajectoryKind::FORCES:
      xbuffer = x_forces.readHost(astart, rsw.natom);
      ybuffer = y_forces.readHost(astart, rsw.natom);
      zbuffer = z_forces.readHost(astart, rsw.natom);
      break;
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    switch (trajkind) {
    case TrajectoryKind::POSITIONS:
      switch (cycle_position) {
      case CoordinateCycle::PAST:
        xbuffer = x_prior_coordinates.readDevice(astart, rsw.natom);
        ybuffer = y_prior_coordinates.readDevice(astart, rsw.natom);
        zbuffer = z_prior_coordinates.readDevice(astart, rsw.natom);
        break;
      case CoordinateCycle::PRESENT:
        xbuffer = x_coordinates.readDevice(astart, rsw.natom);
        ybuffer = y_coordinates.readDevice(astart, rsw.natom);
        zbuffer = z_coordinates.readDevice(astart, rsw.natom);
        break;
      case CoordinateCycle::FUTURE:
        xbuffer = x_future_coordinates.readDevice(astart, rsw.natom);
        ybuffer = y_future_coordinates.readDevice(astart, rsw.natom);
        zbuffer = z_future_coordinates.readDevice(astart, rsw.natom);
        break;
      }
      break;
    case TrajectoryKind::VELOCITIES:
      xbuffer = x_velocities.readDevice(astart, rsw.natom);
      ybuffer = y_velocities.readDevice(astart, rsw.natom);
      zbuffer = z_velocities.readDevice(astart, rsw.natom);
      break;
    case TrajectoryKind::FORCES:
      xbuffer = x_forces.readDevice(astart, rsw.natom);
      ybuffer = y_forces.readDevice(astart, rsw.natom);
      zbuffer = z_forces.readDevice(astart, rsw.natom);
      break;
    }
    break;
#endif
  }
  double deflation;
  switch (trajkind) {
  case TrajectoryKind::POSITIONS:
    deflation = inverse_globalpos_scale;
    break;
  case TrajectoryKind::VELOCITIES:
    deflation = inverse_velocity_scale;
    break;
  case TrajectoryKind::FORCES:
    deflation = inverse_force_scale;
    break;
  }
  for (int i = 0; i < rsw.natom; i++) {
    rsw.xcrd[i] = static_cast<double>(xbuffer[i]) * deflation;
    rsw.ycrd[i] = static_cast<double>(ybuffer[i]) * deflation;
    rsw.zcrd[i] = static_cast<double>(zbuffer[i]) * deflation;
  }
  if (globalpos_scale_bits > globalpos_scale_nonoverflow_bits) {
    std::vector<int> x_ovrf_buffer, y_ovrf_buffer, z_ovrf_buffer;
    switch (tier) {
    case HybridTargetLevel::HOST:
      switch (trajkind) {
      case TrajectoryKind::POSITIONS:
        switch (cycle_position) {
        case CoordinateCycle::PAST:
          x_ovrf_buffer = x_prior_coord_overflow.readHost(astart, rsw.natom);
          y_ovrf_buffer = y_prior_coord_overflow.readHost(astart, rsw.natom);
          z_ovrf_buffer = z_prior_coord_overflow.readHost(astart, rsw.natom);
          break;
        case CoordinateCycle::PRESENT:
          x_ovrf_buffer = x_coordinate_overflow.readHost(astart, rsw.natom);
          y_ovrf_buffer = y_coordinate_overflow.readHost(astart, rsw.natom);
          z_ovrf_buffer = z_coordinate_overflow.readHost(astart, rsw.natom);
          break;
        case CoordinateCycle::FUTURE:
          x_ovrf_buffer = x_future_coord_overflow.readHost(astart, rsw.natom);
          y_ovrf_buffer = y_future_coord_overflow.readHost(astart, rsw.natom);
          z_ovrf_buffer = z_future_coord_overflow.readHost(astart, rsw.natom);
          break;
        }
        break;
      case TrajectoryKind::VELOCITIES:
        x_ovrf_buffer = x_velocity_overflow.readHost(astart, rsw.natom);
        y_ovrf_buffer = y_velocity_overflow.readHost(astart, rsw.natom);
        z_ovrf_buffer = z_velocity_overflow.readHost(astart, rsw.natom);
        break;
      case TrajectoryKind::FORCES:
        x_ovrf_buffer = x_force_overflow.readHost(astart, rsw.natom);
        y_ovrf_buffer = y_force_overflow.readHost(astart, rsw.natom);
        z_ovrf_buffer = z_force_overflow.readHost(astart, rsw.natom);
        break;
      }
      break;
#ifdef STORMM_USE_HPC
    case HybridTargetLevel::DEVICE:
      switch (trajkind) {
      case TrajectoryKind::POSITIONS:
        switch (cycle_position) {
        case CoordinateCycle::PAST:
          x_ovrf_buffer = x_prior_coord_overflow.readDevice(astart, rsw.natom);
          y_ovrf_buffer = y_prior_coord_overflow.readDevice(astart, rsw.natom);
          z_ovrf_buffer = z_prior_coord_overflow.readDevice(astart, rsw.natom);
          break;
        case CoordinateCycle::PRESENT:
          x_ovrf_buffer = x_coordinate_overflow.readDevice(astart, rsw.natom);
          y_ovrf_buffer = y_coordinate_overflow.readDevice(astart, rsw.natom);
          z_ovrf_buffer = z_coordinate_overflow.readDevice(astart, rsw.natom);
          break;
        case CoordinateCycle::FUTURE:
          x_ovrf_buffer = x_future_coord_overflow.readDevice(astart, rsw.natom);
          y_ovrf_buffer = y_future_coord_overflow.readDevice(astart, rsw.natom);
          z_ovrf_buffer = z_future_coord_overflow.readDevice(astart, rsw.natom);
          break;
        }
        break;
      case TrajectoryKind::VELOCITIES:
        x_ovrf_buffer = x_velocity_overflow.readDevice(astart, rsw.natom);
        y_ovrf_buffer = y_velocity_overflow.readDevice(astart, rsw.natom);
        z_ovrf_buffer = z_velocity_overflow.readDevice(astart, rsw.natom);
        break;
      case TrajectoryKind::FORCES:
        x_ovrf_buffer = x_force_overflow.readDevice(astart, rsw.natom);
        y_ovrf_buffer = y_force_overflow.readDevice(astart, rsw.natom);
        z_ovrf_buffer = z_force_overflow.readDevice(astart, rsw.natom);
        break;
      }
      break;
#endif
    }
    const double ovrf_deflation = deflation * max_llint_accumulation;
    for (int i = 0; i < rsw.natom; i++) {
      rsw.xcrd[i] += static_cast<double>(x_ovrf_buffer[i]) * ovrf_deflation;
      rsw.ycrd[i] += static_cast<double>(y_ovrf_buffer[i]) * ovrf_deflation;
      rsw.zcrd[i] += static_cast<double>(z_ovrf_buffer[i]) * ovrf_deflation;
    }  
  }
  return result;
}
  
//-------------------------------------------------------------------------------------------------
#ifdef STORMM_USE_HPC
void PhaseSpaceSynthesis::initializeForces(const GpuDetails &gpu, const HybridTargetLevel tier,
                                           const int index)
#else
void PhaseSpaceSynthesis::initializeForces(const int index)
#endif
{
  if (index >= system_count) {
    rtErr("Index " + std::to_string(index) + " is invalid for a collection of " +
          std::to_string(system_count) + " systems.", "PhaseSpaceSynthesis", "initializeForces");
  }
#ifdef STORMM_USE_HPC
  switch (tier) {
  case HybridTargetLevel::HOST:
    {
#endif
      llint* xptr = x_forces.data();
      llint* yptr = y_forces.data();
      llint* zptr = z_forces.data();
      int* x_ovrf_ptr = x_force_overflow.data();
      int* y_ovrf_ptr = y_force_overflow.data();
      int* z_ovrf_ptr = z_force_overflow.data();
      if (index < 0) {
        for (int i = 0; i < system_count; i++) {
          const int jmin = atom_starts.readHost(i);
          const int jmax = jmin + atom_counts.readHost(i);
          for (int j = jmin; j < jmax; j++) {
            xptr[j] = 0LL;
            yptr[j] = 0LL;
            zptr[j] = 0LL;
          }
          if (force_scale_bits > force_scale_nonoverflow_bits) {
            for (int j = jmin; j < jmax; j++) {
              x_ovrf_ptr[j] = 0;
              y_ovrf_ptr[j] = 0;
              z_ovrf_ptr[j] = 0;
            }
          }
        }
      }
      else {
        const int imin = atom_starts.readHost(index);
        const int imax = imin + atom_counts.readHost(index);
        for (int i = imin; i < imax; i++) {
          xptr[i] = 0LL;
          yptr[i] = 0LL;
          zptr[i] = 0LL;
        }
        if (force_scale_bits > force_scale_nonoverflow_bits) {
          for (int i = imin; i < imax; i++) {
            x_ovrf_ptr[i] = 0;
            y_ovrf_ptr[i] = 0;
            z_ovrf_ptr[i] = 0;
          }
        }
      }
#ifdef STORMM_USE_HPC
    }
    break;
  case HybridTargetLevel::DEVICE:
    {
      PsSynthesisWriter dptr = data(HybridTargetLevel::DEVICE);
      psyInitializeForces(&dptr, index, gpu);
    }
    break;
  }
#endif
}

//-------------------------------------------------------------------------------------------------
#ifdef STORMM_USE_HPC
void PhaseSpaceSynthesis::primeConjugateGradientCalculation(const GpuDetails &gpu,
                                                            const HybridTargetLevel tier)
#else
void PhaseSpaceSynthesis::primeConjugateGradientCalculation()
#endif
{
#ifdef STORMM_USE_HPC
  switch (tier) {
  case HybridTargetLevel::HOST:
    {
#endif
      llint* xfrc_ptr = x_forces.data();
      llint* yfrc_ptr = y_forces.data();
      llint* zfrc_ptr = z_forces.data();
      int* xfrc_ovrf_ptr = x_force_overflow.data();
      int* yfrc_ovrf_ptr = y_force_overflow.data();
      int* zfrc_ovrf_ptr = z_force_overflow.data();
      llint* xprv_ptr = x_forces.data();
      llint* yprv_ptr = y_forces.data();
      llint* zprv_ptr = z_forces.data();
      int* xprv_ovrf_ptr = x_force_overflow.data();
      int* yprv_ovrf_ptr = y_force_overflow.data();
      int* zprv_ovrf_ptr = z_force_overflow.data();
      llint* xvel_ptr = x_velocities.data();
      llint* yvel_ptr = y_velocities.data();
      llint* zvel_ptr = z_velocities.data();
      int* xvel_ovrf_ptr = x_velocity_overflow.data();
      int* yvel_ovrf_ptr = y_velocity_overflow.data();
      int* zvel_ovrf_ptr = z_velocity_overflow.data();
      for (int i = 0; i < system_count; i++) {
        const int jmin = atom_starts.readHost(i);
        const int jmax = jmin + atom_counts.readHost(i);
        for (int j = jmin; j < jmax; j++) {
          xprv_ptr[j] = xfrc_ptr[j];
          yprv_ptr[j] = yfrc_ptr[j];
          zprv_ptr[j] = zfrc_ptr[j];
          xvel_ptr[j] = 0LL;
          yvel_ptr[j] = 0LL;
          zvel_ptr[j] = 0LL;
          if (force_scale_bits > force_scale_nonoverflow_bits) {
            xprv_ovrf_ptr[j] = xfrc_ovrf_ptr[j];
            yprv_ovrf_ptr[j] = yfrc_ovrf_ptr[j];
            zprv_ovrf_ptr[j] = zfrc_ovrf_ptr[j];
            xvel_ovrf_ptr[j] = 0;
            yvel_ovrf_ptr[j] = 0;
            zvel_ovrf_ptr[j] = 0;
          }
        }
      }      
#ifdef STORMM_USE_HPC
    }
    break;
  case HybridTargetLevel::DEVICE:
    {
      PsSynthesisWriter dptr = data(HybridTargetLevel::DEVICE);
      psyPrimeConjugateGradient(&dptr, gpu);
    }
    break;
  }
#endif
}

//-------------------------------------------------------------------------------------------------
#ifdef STORMM_USE_HPC
void PhaseSpaceSynthesis::printTrajectory(const std::vector<int> &system_indices,
                                          const std::string &file_name, const double current_time,
                                          const CoordinateFileKind output_kind,
                                          const PrintSituation expectation, const GpuDetails &gpu)
#else
void PhaseSpaceSynthesis::printTrajectory(const std::vector<int> &system_indices,
                                          const std::string &file_name, const double current_time,
                                          const CoordinateFileKind output_kind,
                                          const PrintSituation expectation)
#endif
{
  // Bail out if there are no frames to print
  const size_t nframe = system_indices.size();
  if (nframe == 0LLU) {
    return;
  }
#ifdef STORMM_USE_HPC
  int low_frame = nframe;
  int high_frame = 0;
  for (size_t i = 0; i < nframe; i++) {
    low_frame = std::min(system_indices[i], low_frame);
    high_frame = std::max(system_indices[i], high_frame);
  }
  high_frame++;
  switch (output_kind) {
  case CoordinateFileKind::AMBER_CRD:
  case CoordinateFileKind::AMBER_INPCRD:
  case CoordinateFileKind::AMBER_NETCDF:
    download(TrajectoryKind::POSITIONS, low_frame, high_frame, gpu);
    break;
  case CoordinateFileKind::AMBER_ASCII_RST:
  case CoordinateFileKind::AMBER_NETCDF_RST:
    download(TrajectoryKind::POSITIONS, low_frame, high_frame, gpu);
    download(TrajectoryKind::VELOCITIES, low_frame, high_frame, gpu);
    break;
  case CoordinateFileKind::UNKNOWN:
    break;
  }
#endif
  
  // Open and write, or append the file with the named frames by transferring the necessary data
  // to temporary arrays (the fixed-precision conversion to real numbers is needed, and the
  // pointers in the writer object are convenient).  Trajectories will write all frames to a
  // single file.  Restart file formats will write individual frames to separate files.
  std::ofstream foutp;
  char buffer[128];

  // If the request if to print a trajectory file, check that all frames have the same sizes.
  PrintSituation actual_expectation = expectation;
  switch (output_kind) {
  case CoordinateFileKind::AMBER_CRD:
    {
      if (expectation == PrintSituation::UNKNOWN) {
        actual_expectation = PrintSituation::APPEND;
      }
      const int natom = atom_counts.readHost(system_indices[0]);
      for (size_t i = 1; i < nframe; i++) {
        if (atom_counts.readHost(system_indices[i]) != natom) {
          rtErr("Frames of different sizes cannot be printed to a single trajectory.",
                "PhaseSpaceSynthesis", "printTrajectory");
        }
      }
      const DataFormat style = (output_kind == CoordinateFileKind::AMBER_CRD) ? DataFormat::ASCII :
                                                                                DataFormat::BINARY;
      const bool fi_exists = (getDrivePathType(file_name) == DrivePathType::FILE);
      foutp = openOutputFile(file_name, actual_expectation, "Open an output trajectory for "
                             "writing PhaseSpaceSynthesis contents.", style);

      // Initialize the trajectory if this is the first time printing it
      if (fi_exists == false || actual_expectation == PrintSituation::OVERWRITE ||
          actual_expectation == PrintSituation::OPEN_NEW) {
        sprintf(buffer, "Generated by STORMM\n");
        foutp.write(buffer, strlen(buffer));
      }
    }
    break;
  case CoordinateFileKind::AMBER_NETCDF:
    break;
  case CoordinateFileKind::AMBER_INPCRD:
  case CoordinateFileKind::AMBER_ASCII_RST:
  case CoordinateFileKind::AMBER_NETCDF_RST:
    if (expectation == PrintSituation::APPEND) {
      rtErr("Input coordinate and restart files cannot be appended.  They can only be written if "
            "no file is present, overwritten if permitted by the user / developer, or generate an "
            "error if the file exists and overwriting is not enabled.", "PhaseSpaceSynthesis",
            "printTrajectory");
    }
    else if (expectation == PrintSituation::UNKNOWN) {
      actual_expectation = PrintSituation::OPEN_NEW;
    }
    break;
  case CoordinateFileKind::UNKNOWN:
    rtErr("Printing request for unknown file type.", "PhaseSpaceSynthesis", "printTrajectory");
  }

  // Proceed frame by frame
  std::vector<double> tmp_xcrd, tmp_ycrd, tmp_zcrd, tmp_xvel, tmp_yvel, tmp_zvel;
  const PsSynthesisReader tpsr = this->data(cycle_position);
  for (size_t i = 0; i < nframe; i++) {
    const int fr_start = atom_starts.readHost(i);
    const int fr_end   = fr_start + atom_counts.readHost(i);
    const int frame_atom_count = fr_end - fr_start;

    // Transfer the particle positions.  Resize the holding arrays at each frame: if the frames
    // are for a trajectory and therefore all the same size, the resizing will take no effort.
    // But, if restart files were requested, the holding arrays may need to be of different sizes.
    tmp_xcrd.resize(frame_atom_count);
    tmp_ycrd.resize(frame_atom_count);
    tmp_zcrd.resize(frame_atom_count);
    for (int j = fr_start; j < fr_end; j++) {
      tmp_xcrd[j - fr_start] = static_cast<double>(tpsr.xcrd[j]) * inverse_globalpos_scale;
      tmp_ycrd[j - fr_start] = static_cast<double>(tpsr.ycrd[j]) * inverse_globalpos_scale;
      tmp_zcrd[j - fr_start] = static_cast<double>(tpsr.zcrd[j]) * inverse_globalpos_scale;
    }
    if (globalpos_scale_bits > globalpos_scale_nonoverflow_bits) {
      const double inv_ovrf_scale = inverse_globalpos_scale * max_llint_accumulation;
      for (int j = fr_start; j < fr_end; j++) {
        tmp_xcrd[j - fr_start] += static_cast<double>(tpsr.xcrd_ovrf[j]) * inv_ovrf_scale;
        tmp_ycrd[j - fr_start] += static_cast<double>(tpsr.ycrd_ovrf[j]) * inv_ovrf_scale;
        tmp_zcrd[j - fr_start] += static_cast<double>(tpsr.zcrd_ovrf[j]) * inv_ovrf_scale;
      }
    }

    // Transfer the particle velocities, if necessary
    switch (output_kind) {
    case CoordinateFileKind::AMBER_CRD:
    case CoordinateFileKind::AMBER_NETCDF:
    case CoordinateFileKind::AMBER_INPCRD:
    case CoordinateFileKind::UNKNOWN:
      break;
    case CoordinateFileKind::AMBER_ASCII_RST:
    case CoordinateFileKind::AMBER_NETCDF_RST:
      {
        tmp_xvel.resize(frame_atom_count);
        tmp_yvel.resize(frame_atom_count);
        tmp_zvel.resize(frame_atom_count);
        for (int j = fr_start; j < fr_end; j++) {
          tmp_xvel[j - fr_start] = static_cast<double>(tpsr.xvel[j]) * inverse_velocity_scale;
          tmp_yvel[j - fr_start] = static_cast<double>(tpsr.yvel[j]) * inverse_velocity_scale;
          tmp_zvel[j - fr_start] = static_cast<double>(tpsr.zvel[j]) * inverse_velocity_scale;
        }
        if (globalpos_scale_bits > globalpos_scale_nonoverflow_bits) {
          const double inv_ovrf_scale = inverse_velocity_scale * max_llint_accumulation;
          for (int j = fr_start; j < fr_end; j++) {
            tmp_xvel[j - fr_start] += static_cast<double>(tpsr.xvel_ovrf[j]) * inv_ovrf_scale;
            tmp_yvel[j - fr_start] += static_cast<double>(tpsr.yvel_ovrf[j]) * inv_ovrf_scale;
            tmp_zvel[j - fr_start] += static_cast<double>(tpsr.zvel_ovrf[j]) * inv_ovrf_scale;
          }
        }
      }
      break;
    }
    std::vector<double> tmp_boxdims(6);
    switch (unit_cell) {
    case UnitCellType::NONE:
      break;
    case UnitCellType::ORTHORHOMBIC:
    case UnitCellType::TRICLINIC:
      {
        const int dim_stride = roundUp(6, warp_size_int);
        for (int j = 0; j < 6; j++) {
          tmp_boxdims[i] = tpsr.boxdims[(i * dim_stride) + j];
        }
      }
    }
    
    // Open and write an individual restart file, or continue writing to a trajectory file
    switch (output_kind) {
    case CoordinateFileKind::AMBER_CRD:
      writeFrame(&foutp, file_name, output_kind, tmp_xcrd, tmp_ycrd, tmp_zcrd, tmp_xvel, tmp_yvel,
                 tmp_zvel, unit_cell, tmp_boxdims);
      break;
    case CoordinateFileKind::AMBER_NETCDF:
      break;
    case CoordinateFileKind::AMBER_INPCRD:
    case CoordinateFileKind::AMBER_ASCII_RST:
      {
        std::string before_dot, after_dot;
        splitPath(file_name, &before_dot, &after_dot);
        std::string aug_file_name = before_dot + "_" + std::to_string(i + 1) + "." + after_dot;
        const DataFormat style = (output_kind == CoordinateFileKind::AMBER_INPCRD ||
                                  output_kind == CoordinateFileKind::AMBER_ASCII_RST) ?
                                 DataFormat::ASCII : DataFormat::BINARY;
        foutp = openOutputFile(aug_file_name, actual_expectation, "Open an output trajectory for "
                               "writing frame " + std::to_string(i + 1) + " of a "
                               "PhaseSpaceSynthesis object's contents.", style);
        sprintf(buffer, "Generated by STORMM\n");
        sprintf(&buffer[strlen(buffer)], "%8d %15.7e\n", frame_atom_count, current_time);
        foutp.write(buffer, strlen(buffer));
        writeFrame(&foutp, file_name, output_kind, tmp_xcrd, tmp_ycrd, tmp_zcrd, tmp_xvel,
                   tmp_yvel, tmp_zvel, unit_cell, tmp_boxdims);
        foutp.close();
      }
      break;
    case CoordinateFileKind::AMBER_NETCDF_RST:
      break;
    case CoordinateFileKind::UNKNOWN:
      break;
    }
  }

  // Close trajectory files spanning the entire synthesis
  switch (output_kind) {
  case CoordinateFileKind::AMBER_CRD:
  case CoordinateFileKind::AMBER_NETCDF:
    foutp.close();
    break;
  case CoordinateFileKind::AMBER_INPCRD:
  case CoordinateFileKind::AMBER_ASCII_RST:
  case CoordinateFileKind::AMBER_NETCDF_RST:
  case CoordinateFileKind::UNKNOWN:
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void PhaseSpaceSynthesis::allocate(const size_t atom_stride) {
  const size_t system_stride = roundUp<size_t>(system_count, warp_size_zu);
  int_data.resize((2LLU * system_stride) + (15LLU * atom_stride));
  atom_counts.setPointer(&int_data, 0LLU, system_count);
  atom_starts.setPointer(&int_data, system_stride, system_count);
  const size_t xfrm_stride = system_count * roundUp<size_t>(9, warp_size_zu);
  const size_t twosys = 2LLU * system_stride;
  llint_data.resize((15LLU * atom_stride) + xfrm_stride);
  x_coordinates.setPointer(&llint_data,                       0LLU, atom_stride);
  y_coordinates.setPointer(&llint_data,                atom_stride, atom_stride);
  z_coordinates.setPointer(&llint_data,         2LLU * atom_stride, atom_stride);
  x_prior_coordinates.setPointer(&llint_data,   3LLU * atom_stride, atom_stride);
  y_prior_coordinates.setPointer(&llint_data,   4LLU * atom_stride, atom_stride);
  z_prior_coordinates.setPointer(&llint_data,   5LLU * atom_stride, atom_stride);
  x_future_coordinates.setPointer(&llint_data,  6LLU * atom_stride, atom_stride);
  y_future_coordinates.setPointer(&llint_data,  7LLU * atom_stride, atom_stride);
  z_future_coordinates.setPointer(&llint_data,  8LLU * atom_stride, atom_stride);
  x_velocities.setPointer(&llint_data,          9LLU * atom_stride, atom_stride);
  y_velocities.setPointer(&llint_data,         10LLU * atom_stride, atom_stride);
  z_velocities.setPointer(&llint_data,         11LLU * atom_stride, atom_stride);
  x_forces.setPointer(&llint_data,             12LLU * atom_stride, atom_stride);
  y_forces.setPointer(&llint_data,             13LLU * atom_stride, atom_stride);
  z_forces.setPointer(&llint_data,             14LLU * atom_stride, atom_stride);
  box_vectors.setPointer(&llint_data,          15LLU * atom_stride, xfrm_stride);
  x_coordinate_overflow.setPointer(&int_data,                           twosys, atom_stride);
  y_coordinate_overflow.setPointer(&int_data,            atom_stride  + twosys, atom_stride);
  z_coordinate_overflow.setPointer(&int_data,   ( 2LLU * atom_stride) + twosys, atom_stride);
  x_prior_coord_overflow.setPointer(&int_data,  ( 3LLU * atom_stride) + twosys, atom_stride);
  y_prior_coord_overflow.setPointer(&int_data,  ( 4LLU * atom_stride) + twosys, atom_stride);
  z_prior_coord_overflow.setPointer(&int_data,  ( 5LLU * atom_stride) + twosys, atom_stride);
  x_future_coord_overflow.setPointer(&int_data, ( 6LLU * atom_stride) + twosys, atom_stride);
  y_future_coord_overflow.setPointer(&int_data, ( 7LLU * atom_stride) + twosys, atom_stride);
  z_future_coord_overflow.setPointer(&int_data, ( 8LLU * atom_stride) + twosys, atom_stride);
  x_velocity_overflow.setPointer(&int_data,     ( 9LLU * atom_stride) + twosys, atom_stride);
  y_velocity_overflow.setPointer(&int_data,     (10LLU * atom_stride) + twosys, atom_stride);
  z_velocity_overflow.setPointer(&int_data,     (11LLU * atom_stride) + twosys, atom_stride);
  x_force_overflow.setPointer(&int_data,        (12LLU * atom_stride) + twosys, atom_stride);
  y_force_overflow.setPointer(&int_data,        (13LLU * atom_stride) + twosys, atom_stride);
  z_force_overflow.setPointer(&int_data,        (14LLU * atom_stride) + twosys, atom_stride);
  double_data.resize(3LLU * xfrm_stride);
  box_space_transforms.setPointer(&double_data,               0LLU, xfrm_stride);
  inverse_transforms.setPointer(&double_data,          xfrm_stride, xfrm_stride);
  box_dimensions.setPointer(&double_data,       2LLU * xfrm_stride, xfrm_stride);
  float_data.resize(3LLU * xfrm_stride);
  sp_box_space_transforms.setPointer(&float_data,             0LLU, xfrm_stride);
  sp_inverse_transforms.setPointer(&float_data,        xfrm_stride, xfrm_stride);
  sp_box_dimensions.setPointer(&float_data,     2LLU * xfrm_stride, xfrm_stride);
}

} // namespace trajectory
} // namespace stormm
