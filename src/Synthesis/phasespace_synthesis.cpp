#include <cmath>
#ifdef STORMM_USE_HPC
#  ifdef STORMM_USE_CUDA
#  include <cuda_runtime.h>
#  endif
#endif
#include "copyright.h"
#include "Constants/hpc_bounds.h"
#include "FileManagement/file_listing.h"
#include "Math/matrix_ops.h"
#include "Math/rounding.h"
#include "Math/summation.h"
#include "Math/vector_ops.h"
#include "Numerics/split_fixed_precision.h"
#include "Trajectory/write_frame.h"
#include "phasespace_synthesis.h"
#ifdef STORMM_USE_HPC
#  include "hpc_phasespace_synthesis.h"
#endif

namespace stormm {
namespace synthesis {

using card::HybridKind;
using card::getEnumerationName;
using diskutil::splitPath;
using diskutil::DataFormat;
using diskutil::DrivePathType;
using diskutil::getDrivePathType;
using math::invertSquareMatrix;
using math::prefixSumInPlace;
using math::PrefixSumType;
using math::roundUp;
using math::sum;
using math::tileVector;
using numerics::checkGlobalPositionBits;
using numerics::checkLocalPositionBits;
using numerics::checkVelocityBits;
using numerics::checkForceBits;
using numerics::force_scale_nonoverflow_bits;
using numerics::globalpos_scale_nonoverflow_bits;
using numerics::max_llint_accumulation;
using numerics::velocity_scale_nonoverflow_bits;
using topology::UnitCellType;
using trajectory::PhaseSpaceWriter;
using trajectory::PhaseSpaceReader;
using trajectory::writeFrame;
  
//-------------------------------------------------------------------------------------------------
PsSynthesisWriter::PsSynthesisWriter(const int system_count_in, const int unique_topology_count_in,
                                     const UnitCellType unit_cell_in,
                                     const ThermostatKind heat_bath_kind_in,
                                     const BarostatKind piston_kind_in, const double time_step_in,
                                     const int* atom_starts_in, const int* atom_counts_in,
                                     const int* common_ag_list_in, const int* common_ag_bounds_in,
                                     const int* unique_ag_idx_in, const int* replica_idx_in,
                                     const double gpos_scale_in, const double lpos_scale_in,
                                     const double vel_scale_in, const double frc_scale_in,
                                     const int gpos_bits_in, const int lpos_bits_in,
                                     const int vel_bits_in, const int frc_bits_in,
                                     llint* boxvecs_in, int* boxvec_ovrf_in, double* umat_in,
                                     double* invu_in, double* boxdims_in, llint* alt_boxvecs_in,
                                     int* alt_boxvec_ovrf_in, double* umat_alt_in,
                                     double* invu_alt_in, double* alt_boxdims_in, llint* xcrd_in,
                                     llint* ycrd_in, llint* zcrd_in, int* xcrd_ovrf_in,
                                     int* ycrd_ovrf_in, int* zcrd_ovrf_in, llint* xvel_in,
                                     llint* yvel_in, llint* zvel_in, int* xvel_ovrf_in,
                                     int* yvel_ovrf_in, int* zvel_ovrf_in, llint* xfrc_in,
                                     llint* yfrc_in, llint* zfrc_in, int* xfrc_ovrf_in,
                                     int* yfrc_ovrf_in, int* zfrc_ovrf_in, llint* xalt_in,
                                     llint* yalt_in, llint* zalt_in, int* xalt_ovrf_in,
                                     int* yalt_ovrf_in, int* zalt_ovrf_in, llint* vxalt_in,
                                     llint* vyalt_in, llint* vzalt_in, int* vxalt_ovrf_in,
                                     int* vyalt_ovrf_in, int* vzalt_ovrf_in, llint* fxalt_in,
                                     llint* fyalt_in, llint* fzalt_in, int* fxalt_ovrf_in,
                                     int* fyalt_ovrf_in, int* fzalt_ovrf_in) :
    system_count{system_count_in}, unique_topology_count{unique_topology_count_in},
    unit_cell{unit_cell_in}, heat_bath_kind{heat_bath_kind_in}, piston_kind{piston_kind_in},
    time_step{time_step_in}, atom_starts{atom_starts_in}, atom_counts{atom_counts_in},
    common_ag_list{common_ag_list_in}, common_ag_bounds{common_ag_bounds_in},
    unique_ag_idx{unique_ag_idx_in}, replica_idx{replica_idx_in}, gpos_scale{gpos_scale_in},
    lpos_scale{lpos_scale_in}, vel_scale{vel_scale_in}, frc_scale{frc_scale_in},
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
    boxvecs{boxvecs_in}, boxvec_ovrf{boxvec_ovrf_in}, umat{umat_in}, invu{invu_in},
    boxdims{boxdims_in}, alt_boxvecs{alt_boxvecs_in}, alt_boxvec_ovrf{alt_boxvec_ovrf_in},
    umat_alt{umat_alt_in}, invu_alt{invu_alt_in}, alt_boxdims{alt_boxdims_in}, xcrd{xcrd_in},
    ycrd{ycrd_in}, zcrd{zcrd_in}, xcrd_ovrf{xcrd_ovrf_in}, ycrd_ovrf{ycrd_ovrf_in},
    zcrd_ovrf{zcrd_ovrf_in}, xvel{xvel_in}, yvel{yvel_in}, zvel{zvel_in}, xvel_ovrf{xvel_ovrf_in},
    yvel_ovrf{yvel_ovrf_in}, zvel_ovrf{zvel_ovrf_in}, xfrc{xfrc_in}, yfrc{yfrc_in}, zfrc{zfrc_in},
    xfrc_ovrf{xfrc_ovrf_in}, yfrc_ovrf{yfrc_ovrf_in}, zfrc_ovrf{zfrc_ovrf_in}, xalt{xalt_in},
    yalt{yalt_in}, zalt{zalt_in}, xalt_ovrf{xalt_ovrf_in}, yalt_ovrf{yalt_ovrf_in},
    zalt_ovrf{zalt_ovrf_in}, vxalt{vxalt_in}, vyalt{vyalt_in}, vzalt{vzalt_in},
    vxalt_ovrf{vxalt_ovrf_in}, vyalt_ovrf{vyalt_ovrf_in}, vzalt_ovrf{vzalt_ovrf_in},
    fxalt{fxalt_in}, fyalt{fyalt_in}, fzalt{fzalt_in}, fxalt_ovrf{fxalt_ovrf_in},
    fyalt_ovrf{fyalt_ovrf_in}, fzalt_ovrf{fzalt_ovrf_in}
{}

//-------------------------------------------------------------------------------------------------
PsSynthesisReader::PsSynthesisReader(const int system_count_in, const int unique_topology_count_in,
                                     const UnitCellType unit_cell_in,
                                     const ThermostatKind heat_bath_kind_in,
                                     const BarostatKind piston_kind_in, const double time_step_in,
                                     const int* atom_starts_in, const int* atom_counts_in,
                                     const int* common_ag_list_in, const int* common_ag_bounds_in,
                                     const int*	unique_ag_idx_in, const int* replica_idx_in,
                                     const double gpos_scale_in, const double lpos_scale_in,
                                     const double vel_scale_in, const double frc_scale_in,
                                     const int gpos_bits_in, const int lpos_bits_in,
                                     const int vel_bits_in, const int frc_bits_in,
                                     const llint* boxvecs_in, const int* boxvec_ovrf_in,
                                     const double* umat_in, const double* invu_in,
                                     const double* boxdims_in, const llint* alt_boxvecs_in,
                                     const int* alt_boxvec_ovrf_in, const double* umat_alt_in,
                                     const double* invu_alt_in, const double* alt_boxdims_in,
                                     const llint* xcrd_in, const llint* ycrd_in,
                                     const llint* zcrd_in, const int* xcrd_ovrf_in,
                                     const int* ycrd_ovrf_in, const int* zcrd_ovrf_in,
                                     const llint* xvel_in, const llint* yvel_in,
                                     const llint* zvel_in, const int* xvel_ovrf_in,
                                     const int* yvel_ovrf_in, const int* zvel_ovrf_in,
                                     const llint* xfrc_in, const llint* yfrc_in,
                                     const llint* zfrc_in, const int* xfrc_ovrf_in,
                                     const int* yfrc_ovrf_in, const int* zfrc_ovrf_in,
                                     const llint* xalt_in, const llint* yalt_in,
                                     const llint* zalt_in, const int* xalt_ovrf_in,
                                     const int* yalt_ovrf_in, const int* zalt_ovrf_in,
                                     const llint* vxalt_in, const llint* vyalt_in,
                                     const llint* vzalt_in, const int* vxalt_ovrf_in,
                                     const int* vyalt_ovrf_in, const int* vzalt_ovrf_in,
                                     const llint* fxalt_in, const llint* fyalt_in,
                                     const llint* fzalt_in, const int* fxalt_ovrf_in,
                                     const int* fyalt_ovrf_in, const int* fzalt_ovrf_in) :
    system_count{system_count_in}, unique_topology_count{unique_topology_count_in},
    unit_cell{unit_cell_in}, heat_bath_kind{heat_bath_kind_in}, piston_kind{piston_kind_in},
    time_step{time_step_in}, atom_starts{atom_starts_in}, atom_counts{atom_counts_in},
    common_ag_list{common_ag_list_in}, common_ag_bounds{common_ag_bounds_in},
    unique_ag_idx{unique_ag_idx_in}, replica_idx{replica_idx_in}, gpos_scale{gpos_scale_in},
    lpos_scale{lpos_scale_in}, vel_scale{vel_scale_in}, frc_scale{frc_scale_in},
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
    boxvecs{boxvecs_in}, boxvec_ovrf{boxvec_ovrf_in}, umat{umat_in}, invu{invu_in},
    boxdims{boxdims_in}, alt_boxvecs{boxvecs_in}, alt_boxvec_ovrf{alt_boxvec_ovrf_in},
    umat_alt{umat_alt_in}, invu_alt{invu_alt_in}, alt_boxdims{alt_boxdims_in}, xcrd{xcrd_in},
    ycrd{ycrd_in}, zcrd{zcrd_in}, xcrd_ovrf{xcrd_ovrf_in}, ycrd_ovrf{ycrd_ovrf_in},
    zcrd_ovrf{zcrd_ovrf_in}, xvel{xvel_in}, yvel{yvel_in}, zvel{zvel_in}, xvel_ovrf{xvel_ovrf_in},
    yvel_ovrf{yvel_ovrf_in}, zvel_ovrf{zvel_ovrf_in}, xfrc{xfrc_in}, yfrc{yfrc_in}, zfrc{zfrc_in},
    xfrc_ovrf{xfrc_ovrf_in}, yfrc_ovrf{yfrc_ovrf_in}, zfrc_ovrf{zfrc_ovrf_in}, xalt{xalt_in},
    yalt{yalt_in}, zalt{zalt_in}, xalt_ovrf{xalt_ovrf_in}, yalt_ovrf{yalt_ovrf_in},
    zalt_ovrf{zalt_ovrf_in}, vxalt{vxalt_in}, vyalt{vyalt_in}, vzalt{vzalt_in},
    vxalt_ovrf{vxalt_ovrf_in}, vyalt_ovrf{vyalt_ovrf_in}, vzalt_ovrf{vzalt_ovrf_in},
    fxalt{fxalt_in}, fyalt{fyalt_in}, fzalt{fzalt_in}, fxalt_ovrf{fxalt_ovrf_in},
    fyalt_ovrf{fyalt_ovrf_in}, fzalt_ovrf{fzalt_ovrf_in}
{}

//-------------------------------------------------------------------------------------------------
PsSynthesisReader::PsSynthesisReader(const PsSynthesisWriter &psyw) :
    system_count{psyw.system_count},
    unique_topology_count{psyw.unique_topology_count},
    unit_cell{psyw.unit_cell},
    heat_bath_kind{psyw.heat_bath_kind},
    piston_kind{psyw.piston_kind},
    time_step{psyw.time_step},
    atom_starts{psyw.atom_starts},
    atom_counts{psyw.atom_counts},
    common_ag_list{psyw.common_ag_list},
    common_ag_bounds{psyw.common_ag_bounds},
    unique_ag_idx{psyw.unique_ag_idx},
    replica_idx{psyw.replica_idx},
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
    boxvec_ovrf{psyw.boxvec_ovrf},
    umat{psyw.umat},
    invu{psyw.invu},
    boxdims{psyw.boxdims},
    alt_boxvecs{psyw.alt_boxvecs},
    alt_boxvec_ovrf{psyw.alt_boxvec_ovrf},
    umat_alt{psyw.umat_alt},
    invu_alt{psyw.invu_alt},
    alt_boxdims{psyw.alt_boxdims},
    xcrd{psyw.xcrd}, ycrd{psyw.ycrd}, zcrd{psyw.zcrd},
    xcrd_ovrf{psyw.xcrd_ovrf}, ycrd_ovrf{psyw.ycrd_ovrf}, zcrd_ovrf{psyw.zcrd_ovrf},
    xvel{psyw.xvel}, yvel{psyw.yvel}, zvel{psyw.zvel},
    xvel_ovrf{psyw.xvel_ovrf}, yvel_ovrf{psyw.yvel_ovrf}, zvel_ovrf{psyw.zvel_ovrf},
    xfrc{psyw.xfrc}, yfrc{psyw.yfrc}, zfrc{psyw.zfrc},
    xfrc_ovrf{psyw.xfrc_ovrf}, yfrc_ovrf{psyw.yfrc_ovrf}, zfrc_ovrf{psyw.zfrc_ovrf},
    xalt{psyw.xalt}, yalt{psyw.yalt}, zalt{psyw.zalt},
    xalt_ovrf{psyw.xalt_ovrf}, yalt_ovrf{psyw.yalt_ovrf}, zalt_ovrf{psyw.zalt_ovrf},
    vxalt{psyw.vxalt}, vyalt{psyw.vyalt}, vzalt{psyw.vzalt},
    vxalt_ovrf{psyw.vxalt_ovrf}, vyalt_ovrf{psyw.vyalt_ovrf}, vzalt_ovrf{psyw.vzalt_ovrf},
    fxalt{psyw.fxalt}, fyalt{psyw.fyalt}, fzalt{psyw.fzalt},
    fxalt_ovrf{psyw.fxalt_ovrf}, fyalt_ovrf{psyw.fyalt_ovrf}, fzalt_ovrf{psyw.fzalt_ovrf}
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
    unique_topology_count{0},
    unit_cell{UnitCellType::NONE},
    cycle_position{CoordinateCycle::PRIMARY},
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
    shared_topology_instances{HybridKind::POINTER, "labframe_instances"},
    shared_topology_instance_bounds{HybridKind::POINTER, "labframe_instance_bnds"},
    unique_topology_reference{HybridKind::POINTER, "labframe_unique_idx"},
    shared_topology_instance_index{HybridKind::POINTER, "labframe_shared_idx"},
    x_coordinates{HybridKind::POINTER, "labframe_xpos"},
    y_coordinates{HybridKind::POINTER, "labframe_ypos"},
    z_coordinates{HybridKind::POINTER, "labframe_zpos"},
    x_coordinate_overflow{HybridKind::POINTER, "labframe_xpos_ovrf"},
    y_coordinate_overflow{HybridKind::POINTER, "labframe_ypos_ovrf"},
    z_coordinate_overflow{HybridKind::POINTER, "labframe_zpos_ovrf"},
    x_alt_coordinates{HybridKind::POINTER, "labframe_x_prev"},
    y_alt_coordinates{HybridKind::POINTER, "labframe_y_prev"},
    z_alt_coordinates{HybridKind::POINTER, "labframe_z_prev"},
    x_alt_coord_overflow{HybridKind::POINTER, "labframe_xalt_ovrf"},
    y_alt_coord_overflow{HybridKind::POINTER, "labframe_yalt_ovrf"},
    z_alt_coord_overflow{HybridKind::POINTER, "labframe_zalt_ovrf"},
    x_velocities{HybridKind::POINTER, "labframe_vx"},
    y_velocities{HybridKind::POINTER, "labframe_vy"},
    z_velocities{HybridKind::POINTER, "labframe_vz"},
    x_velocity_overflow{HybridKind::POINTER, "labframe_vx_ovrf"},
    y_velocity_overflow{HybridKind::POINTER, "labframe_vy_ovrf"},
    z_velocity_overflow{HybridKind::POINTER, "labframe_vz_ovrf"},
    x_alt_velocities{HybridKind::POINTER, "labframe_vx_prev"},
    y_alt_velocities{HybridKind::POINTER, "labframe_vy_prev"},
    z_alt_velocities{HybridKind::POINTER, "labframe_vz_prev"},
    x_alt_velocity_overflow{HybridKind::POINTER, "labframe_valtx_ovrf"},
    y_alt_velocity_overflow{HybridKind::POINTER, "labframe_valty_ovrf"},
    z_alt_velocity_overflow{HybridKind::POINTER, "labframe_valtz_ovrf"},
    x_forces{HybridKind::POINTER, "labframe_fx"},
    y_forces{HybridKind::POINTER, "labframe_fy"},
    z_forces{HybridKind::POINTER, "labframe_fz"},
    x_force_overflow{HybridKind::POINTER, "labframe_fx_ovrf"},
    y_force_overflow{HybridKind::POINTER, "labframe_fy_ovrf"},
    z_force_overflow{HybridKind::POINTER, "labframe_fz_ovrf"},
    x_alt_forces{HybridKind::POINTER, "labframe_fx_prev"},
    y_alt_forces{HybridKind::POINTER, "labframe_fy_prev"},
    z_alt_forces{HybridKind::POINTER, "labframe_fz_prev"},
    x_alt_force_overflow{HybridKind::POINTER, "labframe_faltx_ovrf"},
    y_alt_force_overflow{HybridKind::POINTER, "labframe_falty_ovrf"},
    z_alt_force_overflow{HybridKind::POINTER, "labframe_faltz_ovrf"},
    box_vectors{HybridKind::POINTER, "labframe_boxvecs"},
    box_vector_overflow{HybridKind::POINTER, "labframe_boxvec_ovrf"},
    box_space_transforms{HybridKind::POINTER, "labframe_umat"},
    inverse_transforms{HybridKind::POINTER, "labframe_invu"},
    box_dimensions{HybridKind::POINTER, "labframe_dims"},
    alt_box_vectors{HybridKind::POINTER, "labframe_boxvecs"},
    alt_box_vector_overflow{HybridKind::POINTER, "labframe_boxvec_ovrf"},
    alt_box_transforms{HybridKind::POINTER, "labframe_umat"},
    alt_inverse_transforms{HybridKind::POINTER, "labframe_invu"},
    alt_box_dimensions{HybridKind::POINTER, "labframe_dims"},
    int_data{HybridKind::ARRAY, "labframe_int"},
    llint_data{HybridKind::ARRAY, "labframe_llint"},
    double_data{HybridKind::ARRAY, "labframe_double"},
    topologies{ag_list},
    unique_topologies{}
{
  // Check validity of input
  if (ag_list.size() != ps_list.size()) {
    rtErr("The number of input topologies (" + std::to_string(ag_list.size()) + ") must match the "
          "number of systems (" + std::to_string(ps_list.size()) + ").", "PhaseSpaceSynthesis");
  }
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

  // Find the unique topologies
  std::vector<bool> unique(system_count, true);
  for (int i = 0; i < system_count; i++) {
    if (unique[i]) {
      for (int j = i + 1; j < system_count; j++) {
        unique[j] = (unique[j] && (topologies[j] != topologies[i]));
      }
    }
  }
  unique_topology_count = sum<int>(unique);
  unique_topologies.reserve(unique_topology_count);
  for (int i = 0; i < system_count; i++) {
    if (unique[i]) {
      unique_topologies.push_back(const_cast<AtomGraph*>(topologies[i]));
    }
  }  
  
  // Allocate data and set internal pointers
  int atom_stride = 0;
  for (int i = 0; i < system_count; i++) {
    atom_stride += roundUp(ps_list[i].getAtomCount(), warp_size_int);
  }
  allocate(atom_stride);

  // Survey all systems and list all examples using each unique topology.
  int* stib_ptr = shared_topology_instance_bounds.data();
  for (int i = 0; i < system_count; i++) {
    const AtomGraph* iag_ptr = topologies[i];
    for (int j = 0; j < unique_topology_count; j++) {
      if (iag_ptr == unique_topologies[j]) {
        stib_ptr[j] += 1;
      }
    }
  }
  prefixSumInPlace(&shared_topology_instance_bounds, PrefixSumType::EXCLUSIVE,
                   "PhaseSpaceSynthesis");
  std::vector<int> stib_counters = shared_topology_instance_bounds.readHost();
  if (stib_counters.back() != system_count) {
    rtErr("Counts of systems linked to each unique topology are incorrect.",
          "PhaseSpaceSynthesis");
  }
  int* sti_ptr  = shared_topology_instances.data();
  int* replica_ptr = shared_topology_instance_index.data();
  int* utr_ptr = unique_topology_reference.data();
  for (int i = 0; i < system_count; i++) {
    const AtomGraph* iag_ptr = topologies[i];
    for (int j = 0; j < unique_topology_count; j++) {
      if (iag_ptr == unique_topologies[j]) {
        sti_ptr[stib_counters[j]] = i;
        replica_ptr[i] = stib_counters[j] - stib_ptr[j];
        utr_ptr[i] = j;
        stib_counters[j] += 1;
      }
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

  // Loop over all systems
  const int mtrx_stride = roundUp(9, warp_size_int);
  const int dim_stride = roundUp(6, warp_size_int);
  llint* xpos_ptr = x_coordinates.data();
  llint* ypos_ptr = y_coordinates.data();
  llint* zpos_ptr = z_coordinates.data();
  int* xpos_ovrf_ptr  = x_coordinate_overflow.data();
  int* ypos_ovrf_ptr  = y_coordinate_overflow.data();
  int* zpos_ovrf_ptr  = z_coordinate_overflow.data();
  llint* xpos_alt_ptr = x_alt_coordinates.data();
  llint* ypos_alt_ptr = y_alt_coordinates.data();
  llint* zpos_alt_ptr = z_alt_coordinates.data();
  int* xpos_alt_ovrf_ptr  = x_alt_coord_overflow.data();
  int* ypos_alt_ovrf_ptr  = y_alt_coord_overflow.data();
  int* zpos_alt_ovrf_ptr  = z_alt_coord_overflow.data();
  llint* xvel_ptr = x_velocities.data();
  llint* yvel_ptr = y_velocities.data();
  llint* zvel_ptr = z_velocities.data();
  int* xvel_ovrf_ptr  = x_velocity_overflow.data();
  int* yvel_ovrf_ptr  = y_velocity_overflow.data();
  int* zvel_ovrf_ptr  = z_velocity_overflow.data();
  llint* xvel_alt_ptr = x_alt_velocities.data();
  llint* yvel_alt_ptr = y_alt_velocities.data();
  llint* zvel_alt_ptr = z_alt_velocities.data();
  int* xvel_alt_ovrf_ptr  = x_alt_velocity_overflow.data();
  int* yvel_alt_ovrf_ptr  = y_alt_velocity_overflow.data();
  int* zvel_alt_ovrf_ptr  = z_alt_velocity_overflow.data();
  llint* xfrc_ptr = x_forces.data();
  llint* yfrc_ptr = y_forces.data();
  llint* zfrc_ptr = z_forces.data();
  int* xfrc_ovrf_ptr  = x_force_overflow.data();
  int* yfrc_ovrf_ptr  = y_force_overflow.data();
  int* zfrc_ovrf_ptr  = z_force_overflow.data();
  llint* xfrc_alt_ptr = x_alt_forces.data();
  llint* yfrc_alt_ptr = y_alt_forces.data();
  llint* zfrc_alt_ptr = z_alt_forces.data();
  int* xfrc_alt_ovrf_ptr  = x_alt_force_overflow.data();
  int* yfrc_alt_ovrf_ptr  = y_alt_force_overflow.data();
  int* zfrc_alt_ovrf_ptr  = z_alt_force_overflow.data();
  for (int i = 0; i < system_count; i++) {

    // Get a reader for the PhaseSpace object's host-side data
    PhaseSpaceReader psr = ps_list[i].data();

    // Assign coordinates, converting from double precision to fixed precision format
    const int asi = atom_starts.readHost(i);

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
    hostDoubleToInt95(psr.xcrd, psr.ycrd, psr.zcrd, &xpos_ptr[asi], &xpos_ovrf_ptr[asi],
                      &ypos_ptr[asi], &ypos_ovrf_ptr[asi], &zpos_ptr[asi], &zpos_ovrf_ptr[asi],
                      psr.natom, globalpos_scale);
    hostDoubleToInt95(psr.xalt, psr.yalt, psr.zalt, &xpos_alt_ptr[asi], &xpos_alt_ovrf_ptr[asi],
                      &ypos_alt_ptr[asi], &ypos_alt_ovrf_ptr[asi], &zpos_alt_ptr[asi],
                      &zpos_alt_ovrf_ptr[asi], psr.natom, globalpos_scale);
    hostDoubleToInt95(psr.xvel, psr.yvel, psr.zvel, &xvel_ptr[asi], &xvel_ovrf_ptr[asi],
                      &yvel_ptr[asi], &yvel_ovrf_ptr[asi], &zvel_ptr[asi], &zvel_ovrf_ptr[asi],
                      psr.natom, velocity_scale);
    hostDoubleToInt95(psr.vxalt, psr.vyalt, psr.vzalt, &xvel_alt_ptr[asi], &xvel_alt_ovrf_ptr[asi],
                      &yvel_alt_ptr[asi], &yvel_alt_ovrf_ptr[asi], &zvel_alt_ptr[asi],
                      &zvel_alt_ovrf_ptr[asi], psr.natom, velocity_scale);
    hostDoubleToInt95(psr.xfrc, psr.yfrc, psr.zfrc, &xfrc_ptr[asi], &xfrc_ovrf_ptr[asi],
                      &yfrc_ptr[asi], &yfrc_ovrf_ptr[asi], &zfrc_ptr[asi], &zfrc_ovrf_ptr[asi],
                      psr.natom, force_scale);
    hostDoubleToInt95(psr.fxalt, psr.fyalt, psr.fzalt, &xfrc_alt_ptr[asi], &xfrc_alt_ovrf_ptr[asi],
                      &yfrc_alt_ptr[asi], &yfrc_alt_ovrf_ptr[asi], &zfrc_alt_ptr[asi],
                      &zfrc_alt_ovrf_ptr[asi], psr.natom, force_scale);

    // Handle the box space transformation.  The transformation matrices of the labframe will
    // become slightly desynchronized from the original PhaseSpace object, but this is in the
    // interest of making a system which can be represented in fixed precision, pixelated
    // coordinates with box vectors which are likewise multiples of the positional discretization.
    llint* bv_ptr = box_vectors.data();
    int* bv_ovrf_ptr = box_vector_overflow.data();
    hostDoubleToInt95(psr.invu, &bv_ptr[i * mtrx_stride], &bv_ovrf_ptr[i * mtrx_stride], 9,
                      globalpos_scale);
    for (int j = 0; j < 9; j++) {
      const size_t imj = (i * mtrx_stride) + j;
      const double d_invu = hostInt95ToDouble(bv_ptr[imj], bv_ovrf_ptr[imj]);
      inverse_transforms.putHost(d_invu, imj);
      alt_inverse_transforms.putHost(d_invu, imj);
    }
    const double* invu_ptr = inverse_transforms.data();
    double* umat_ptr = box_space_transforms.data();
    invertSquareMatrix(&invu_ptr[i * mtrx_stride], &umat_ptr[i * mtrx_stride], 3);
    for (int j = 0; j < 6; j++) {
      box_dimensions.putHost(psr.boxdim[j], (i * dim_stride) + j);
      alt_box_dimensions.putHost(psr.boxdim[j], (i * dim_stride) + j);
    }

    // Copy the box transformations into the ALTERNATE space
    llint* alt_bv_ptr = alt_box_vectors.data();
    int* alt_bv_ovrf_ptr = alt_box_vector_overflow.data();
    double* alt_umat_ptr = alt_box_transforms.data();
    for (int j = 0; j < 9; j++) {
      const size_t imj = (i * mtrx_stride) + j;
      alt_bv_ptr[imj]      = bv_ptr[imj];
      alt_bv_ovrf_ptr[imj] = bv_ovrf_ptr[imj];
      alt_umat_ptr[imj]    = umat_ptr[imj];
    }
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
    unique_topology_count{original.unique_topology_count},
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
    shared_topology_instances{original.shared_topology_instances},
    shared_topology_instance_bounds{original.shared_topology_instance_bounds},
    unique_topology_reference{original.unique_topology_reference},
    shared_topology_instance_index{original.shared_topology_instance_index},
    x_coordinates{original.x_coordinates},
    y_coordinates{original.y_coordinates},
    z_coordinates{original.z_coordinates},
    x_coordinate_overflow{original.x_coordinate_overflow},
    y_coordinate_overflow{original.y_coordinate_overflow},
    z_coordinate_overflow{original.z_coordinate_overflow},
    x_alt_coordinates{original.x_alt_coordinates},
    y_alt_coordinates{original.y_alt_coordinates},
    z_alt_coordinates{original.z_alt_coordinates},
    x_alt_coord_overflow{original.x_alt_coord_overflow},
    y_alt_coord_overflow{original.y_alt_coord_overflow},
    z_alt_coord_overflow{original.z_alt_coord_overflow},
    x_velocities{original.x_velocities},
    y_velocities{original.y_velocities},
    z_velocities{original.z_velocities},
    x_velocity_overflow{original.x_velocity_overflow},
    y_velocity_overflow{original.y_velocity_overflow},
    z_velocity_overflow{original.z_velocity_overflow},
    x_alt_velocities{original.x_alt_velocities},
    y_alt_velocities{original.y_alt_velocities},
    z_alt_velocities{original.z_alt_velocities},
    x_alt_velocity_overflow{original.x_alt_velocity_overflow},
    y_alt_velocity_overflow{original.y_alt_velocity_overflow},
    z_alt_velocity_overflow{original.z_alt_velocity_overflow},
    x_forces{original.x_forces},
    y_forces{original.y_forces},
    z_forces{original.z_forces},
    x_force_overflow{original.x_force_overflow},
    y_force_overflow{original.y_force_overflow},
    z_force_overflow{original.z_force_overflow},
    x_alt_forces{original.x_alt_forces},
    y_alt_forces{original.y_alt_forces},
    z_alt_forces{original.z_alt_forces},
    x_alt_force_overflow{original.x_alt_force_overflow},
    y_alt_force_overflow{original.y_alt_force_overflow},
    z_alt_force_overflow{original.z_alt_force_overflow},
    box_vectors{original.box_vectors},
    box_vector_overflow{original.box_vector_overflow},
    box_space_transforms{original.box_space_transforms},
    inverse_transforms{original.inverse_transforms},
    box_dimensions{original.box_dimensions},
    alt_box_vectors{original.alt_box_vectors},
    alt_box_vector_overflow{original.alt_box_vector_overflow},
    alt_box_transforms{original.alt_box_transforms},
    alt_inverse_transforms{original.alt_inverse_transforms},
    alt_box_dimensions{original.alt_box_dimensions},
    int_data{original.int_data},
    llint_data{original.llint_data},
    double_data{original.double_data},
    topologies{original.topologies},
    unique_topologies{original.unique_topologies}
{
  // The allocate function again handles pointer repair, just like in the PhaseSpace object.
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
    unique_topology_count{original.unique_topology_count},
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
    shared_topology_instances{std::move(original.shared_topology_instances)},
    shared_topology_instance_bounds{std::move(original.shared_topology_instance_bounds)},
    unique_topology_reference{std::move(original.unique_topology_reference)},
    shared_topology_instance_index{std::move(original.shared_topology_instance_index)},
    x_coordinates{std::move(original.x_coordinates)},
    y_coordinates{std::move(original.y_coordinates)},
    z_coordinates{std::move(original.z_coordinates)},
    x_coordinate_overflow{std::move(original.x_coordinate_overflow)},
    y_coordinate_overflow{std::move(original.y_coordinate_overflow)},
    z_coordinate_overflow{std::move(original.z_coordinate_overflow)},
    x_alt_coordinates{std::move(original.x_alt_coordinates)},
    y_alt_coordinates{std::move(original.y_alt_coordinates)},
    z_alt_coordinates{std::move(original.z_alt_coordinates)},
    x_alt_coord_overflow{std::move(original.x_alt_coord_overflow)},
    y_alt_coord_overflow{std::move(original.y_alt_coord_overflow)},
    z_alt_coord_overflow{std::move(original.z_alt_coord_overflow)},
    x_velocities{std::move(original.x_velocities)},
    y_velocities{std::move(original.y_velocities)},
    z_velocities{std::move(original.z_velocities)},
    x_velocity_overflow{std::move(original.x_velocity_overflow)},
    y_velocity_overflow{std::move(original.y_velocity_overflow)},
    z_velocity_overflow{std::move(original.z_velocity_overflow)},
    x_alt_velocities{std::move(original.x_alt_velocities)},
    y_alt_velocities{std::move(original.y_alt_velocities)},
    z_alt_velocities{std::move(original.z_alt_velocities)},
    x_alt_velocity_overflow{std::move(original.x_alt_velocity_overflow)},
    y_alt_velocity_overflow{std::move(original.y_alt_velocity_overflow)},
    z_alt_velocity_overflow{std::move(original.z_alt_velocity_overflow)},
    x_forces{std::move(original.x_forces)},
    y_forces{std::move(original.y_forces)},
    z_forces{std::move(original.z_forces)},
    x_force_overflow{std::move(original.x_force_overflow)},
    y_force_overflow{std::move(original.y_force_overflow)},
    z_force_overflow{std::move(original.z_force_overflow)},
    x_alt_forces{std::move(original.x_alt_forces)},
    y_alt_forces{std::move(original.y_alt_forces)},
    z_alt_forces{std::move(original.z_alt_forces)},
    x_alt_force_overflow{std::move(original.x_alt_force_overflow)},
    y_alt_force_overflow{std::move(original.y_alt_force_overflow)},
    z_alt_force_overflow{std::move(original.z_alt_force_overflow)},
    box_vectors{std::move(original.box_vectors)},
    box_vector_overflow{std::move(original.box_vector_overflow)},
    box_space_transforms{std::move(original.box_space_transforms)},
    inverse_transforms{std::move(original.inverse_transforms)},
    box_dimensions{std::move(original.box_dimensions)},
    alt_box_vectors{std::move(original.alt_box_vectors)},
    alt_box_vector_overflow{std::move(original.alt_box_vector_overflow)},
    alt_box_transforms{std::move(original.alt_box_transforms)},
    alt_inverse_transforms{std::move(original.alt_inverse_transforms)},
    alt_box_dimensions{std::move(original.alt_box_dimensions)},
    int_data{std::move(original.int_data)},
    llint_data{std::move(original.llint_data)},
    double_data{std::move(original.double_data)},
    topologies{std::move(original.topologies)},
    unique_topologies{std::move(original.unique_topologies)}
{}

//-------------------------------------------------------------------------------------------------
const PhaseSpaceSynthesis* PhaseSpaceSynthesis::getSelfPointer() const {
  return this;
}

//-------------------------------------------------------------------------------------------------
const PsSynthesisReader PhaseSpaceSynthesis::data(const HybridTargetLevel tier) const {
  return data(cycle_position, tier);
}

//-------------------------------------------------------------------------------------------------
const PsSynthesisReader PhaseSpaceSynthesis::data(const CoordinateCycle orientation,
                                                  const HybridTargetLevel tier) const {
  switch (orientation) {
  case CoordinateCycle::PRIMARY:
    return PsSynthesisReader(system_count, unique_topology_count, unit_cell, heat_bath_kind,
                             piston_kind, time_step, atom_starts.data(tier),
                             atom_counts.data(tier), shared_topology_instances.data(tier),
                             shared_topology_instance_bounds.data(tier),
                             unique_topology_reference.data(tier),
                             shared_topology_instance_index.data(tier), globalpos_scale,
                             localpos_scale, velocity_scale, force_scale, globalpos_scale_bits,
                             localpos_scale_bits, velocity_scale_bits, force_scale_bits,
                             box_vectors.data(tier), box_vector_overflow.data(tier),
                             box_space_transforms.data(tier), inverse_transforms.data(tier),
                             box_dimensions.data(tier), alt_box_vectors.data(tier),
                             alt_box_vector_overflow.data(tier), alt_box_transforms.data(tier),
                             alt_inverse_transforms.data(tier), alt_box_dimensions.data(tier),
                             x_coordinates.data(tier), y_coordinates.data(tier),
                             z_coordinates.data(tier), x_coordinate_overflow.data(tier),
                             y_coordinate_overflow.data(tier), z_coordinate_overflow.data(tier),
                             x_velocities.data(tier), y_velocities.data(tier),
                             z_velocities.data(tier), x_velocity_overflow.data(tier),
                             y_velocity_overflow.data(tier), z_velocity_overflow.data(tier),
                             x_forces.data(tier), y_forces.data(tier), z_forces.data(tier),
                             x_force_overflow.data(tier), y_force_overflow.data(tier),
                             z_force_overflow.data(tier), x_alt_coordinates.data(tier),
                             y_alt_coordinates.data(tier), z_alt_coordinates.data(tier),
                             x_alt_coord_overflow.data(tier), y_alt_coord_overflow.data(tier),
                             z_alt_coord_overflow.data(tier), x_alt_velocities.data(tier),
                             y_alt_velocities.data(tier), z_alt_velocities.data(tier),
                             x_alt_velocity_overflow.data(tier),
                             y_alt_velocity_overflow.data(tier),
                             z_alt_velocity_overflow.data(tier), x_alt_forces.data(tier),
                             y_alt_forces.data(tier), z_alt_forces.data(tier),
                             x_alt_force_overflow.data(tier),
                             y_alt_force_overflow.data(tier),
                             z_alt_force_overflow.data(tier));
  case CoordinateCycle::ALTERNATE:
    return PsSynthesisReader(system_count, unique_topology_count, unit_cell, heat_bath_kind,
                             piston_kind, time_step, atom_starts.data(tier),
                             atom_counts.data(tier), shared_topology_instances.data(tier),
                             shared_topology_instance_bounds.data(tier),
                             unique_topology_reference.data(tier),
                             shared_topology_instance_index.data(tier), globalpos_scale,
                             localpos_scale, velocity_scale, force_scale, globalpos_scale_bits,
                             localpos_scale_bits, velocity_scale_bits, force_scale_bits,
                             alt_box_vectors.data(tier), alt_box_vector_overflow.data(tier),
                             alt_box_transforms.data(tier), alt_inverse_transforms.data(tier),
                             alt_box_dimensions.data(tier), box_vectors.data(tier),
                             box_vector_overflow.data(tier), box_space_transforms.data(tier),
                             inverse_transforms.data(tier), box_dimensions.data(tier),
                             x_alt_coordinates.data(tier), y_alt_coordinates.data(tier),
                             z_alt_coordinates.data(tier), x_alt_coord_overflow.data(tier),
                             y_alt_coord_overflow.data(tier), z_alt_coord_overflow.data(tier),
                             x_alt_velocities.data(tier), y_alt_velocities.data(tier),
                             z_alt_velocities.data(tier), x_alt_velocity_overflow.data(tier),
                             y_alt_velocity_overflow.data(tier),
                             z_alt_velocity_overflow.data(tier), x_alt_forces.data(tier),
                             y_alt_forces.data(tier), z_alt_forces.data(tier),
                             x_alt_force_overflow.data(tier), y_alt_force_overflow.data(tier),
                             z_alt_force_overflow.data(tier), x_coordinates.data(tier),
                             y_coordinates.data(tier), z_coordinates.data(tier),
                             x_coordinate_overflow.data(tier), y_coordinate_overflow.data(tier),
                             z_coordinate_overflow.data(tier), x_velocities.data(tier),
                             y_velocities.data(tier), z_velocities.data(tier),
                             x_velocity_overflow.data(tier), y_velocity_overflow.data(tier),
                             z_velocity_overflow.data(tier), x_forces.data(tier),
                             y_forces.data(tier), z_forces.data(tier), x_force_overflow.data(tier),
                             y_force_overflow.data(tier), z_force_overflow.data(tier));
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
  case CoordinateCycle::PRIMARY:
    return PsSynthesisWriter(system_count, unique_topology_count, unit_cell, heat_bath_kind,
                             piston_kind, time_step, atom_starts.data(tier),
                             atom_counts.data(tier), shared_topology_instances.data(tier),
                             shared_topology_instance_bounds.data(tier),
                             unique_topology_reference.data(tier),
                             shared_topology_instance_index.data(tier), globalpos_scale,
                             localpos_scale, velocity_scale, force_scale, globalpos_scale_bits,
                             localpos_scale_bits, velocity_scale_bits, force_scale_bits,
                             box_vectors.data(tier), box_vector_overflow.data(tier),
                             box_space_transforms.data(tier), inverse_transforms.data(tier),
                             box_dimensions.data(tier), alt_box_vectors.data(tier),
                             alt_box_vector_overflow.data(tier), alt_box_transforms.data(tier),
                             alt_inverse_transforms.data(tier), alt_box_dimensions.data(tier),
                             x_coordinates.data(tier), y_coordinates.data(tier),
                             z_coordinates.data(tier), x_coordinate_overflow.data(tier),
                             y_coordinate_overflow.data(tier), z_coordinate_overflow.data(tier),
                             x_velocities.data(tier), y_velocities.data(tier),
                             z_velocities.data(tier), x_velocity_overflow.data(tier),
                             y_velocity_overflow.data(tier), z_velocity_overflow.data(tier),
                             x_forces.data(tier), y_forces.data(tier), z_forces.data(tier),
                             x_force_overflow.data(tier), y_force_overflow.data(tier),
                             z_force_overflow.data(tier), x_alt_coordinates.data(tier),
                             y_alt_coordinates.data(tier), z_alt_coordinates.data(tier),
                             x_alt_coord_overflow.data(tier), y_alt_coord_overflow.data(tier),
                             z_alt_coord_overflow.data(tier), x_alt_velocities.data(tier),
                             y_alt_velocities.data(tier), z_alt_velocities.data(tier),
                             x_alt_velocity_overflow.data(tier),
                             y_alt_velocity_overflow.data(tier),
                             z_alt_velocity_overflow.data(tier), x_alt_forces.data(tier),
                             y_alt_forces.data(tier), z_alt_forces.data(tier),
                             x_alt_force_overflow.data(tier),
                             y_alt_force_overflow.data(tier),
                             z_alt_force_overflow.data(tier));
  case CoordinateCycle::ALTERNATE:
    return PsSynthesisWriter(system_count, unique_topology_count, unit_cell, heat_bath_kind,
                             piston_kind, time_step, atom_starts.data(tier),
                             atom_counts.data(tier), shared_topology_instances.data(tier),
                             shared_topology_instance_bounds.data(tier),
                             unique_topology_reference.data(tier),
                             shared_topology_instance_index.data(tier), globalpos_scale,
                             localpos_scale, velocity_scale, force_scale, globalpos_scale_bits,
                             localpos_scale_bits, velocity_scale_bits, force_scale_bits,
                             alt_box_vectors.data(tier), alt_box_vector_overflow.data(tier),
                             alt_box_transforms.data(tier), alt_inverse_transforms.data(tier),
                             alt_box_dimensions.data(tier), box_vectors.data(tier),
                             box_vector_overflow.data(tier), box_space_transforms.data(tier),
                             inverse_transforms.data(tier), box_dimensions.data(tier),
                             x_alt_coordinates.data(tier), y_alt_coordinates.data(tier),
                             z_alt_coordinates.data(tier), x_alt_coord_overflow.data(tier),
                             y_alt_coord_overflow.data(tier), z_alt_coord_overflow.data(tier),
                             x_alt_velocities.data(tier), y_alt_velocities.data(tier),
                             z_alt_velocities.data(tier), x_alt_velocity_overflow.data(tier),
                             y_alt_velocity_overflow.data(tier),
                             z_alt_velocity_overflow.data(tier), x_alt_forces.data(tier),
                             y_alt_forces.data(tier), z_alt_forces.data(tier),
                             x_alt_force_overflow.data(tier), y_alt_force_overflow.data(tier),
                             z_alt_force_overflow.data(tier), x_coordinates.data(tier),
                             y_coordinates.data(tier), z_coordinates.data(tier),
                             x_coordinate_overflow.data(tier), y_coordinate_overflow.data(tier),
                             z_coordinate_overflow.data(tier), x_velocities.data(tier),
                             y_velocities.data(tier), z_velocities.data(tier),
                             x_velocity_overflow.data(tier), y_velocity_overflow.data(tier),
                             z_velocity_overflow.data(tier), x_forces.data(tier),
                             y_forces.data(tier), z_forces.data(tier), x_force_overflow.data(tier),
                             y_force_overflow.data(tier), z_force_overflow.data(tier));
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
int PhaseSpaceSynthesis::getSystemCount() const {
  return system_count;
}

//-------------------------------------------------------------------------------------------------
int PhaseSpaceSynthesis::getUniqueTopologyCount() const {
  return unique_topology_count;
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

//-------------------------------------------------------------------------------------------------
const AtomGraph* PhaseSpaceSynthesis::getSystemTopologyPointer(const int system_index) const {
  if (system_index < 0 || system_index >= system_count) {
    rtErr("System index " + std::to_string(system_index) + " is invalid for a synthesis with " +
          std::to_string(system_count) + " coordinate sets.", "PhaseSpaceSynthesis",
          "getSystemTopologyPointer");
  }
  return topologies[system_index];
}

//-------------------------------------------------------------------------------------------------
const std::vector<AtomGraph*>& PhaseSpaceSynthesis::getSystemTopologyPointer() const {
  return topologies;
}

//-------------------------------------------------------------------------------------------------
const std::vector<AtomGraph*>& PhaseSpaceSynthesis::getUniqueTopologies() const {
  return unique_topologies;
}

//-------------------------------------------------------------------------------------------------
int PhaseSpaceSynthesis::getUniqueTopologyIndex(const int system_index) const {
  return unique_topology_reference.readHost(system_index);
}

//-------------------------------------------------------------------------------------------------
std::vector<int> PhaseSpaceSynthesis::getUniqueTopologyIndices() const {
  return unique_topology_reference.readHost();
}

//-------------------------------------------------------------------------------------------------
std::vector<int> PhaseSpaceSynthesis::getUniqueTopologyExampleIndices() const {
  std::vector<int> result(unique_topology_count);
  const int* sti_ptr = shared_topology_instances.data();
  for (int i = 0; i < unique_topology_count; i++) {
    result[i] = sti_ptr[shared_topology_instance_bounds.readHost(i)];
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
int PhaseSpaceSynthesis::getTopologyInstanceCount(const int topology_index) const {
  if (topology_index < 0 || topology_index >= unique_topology_count) {
    rtErr("Topology index " + std::to_string(topology_index) + " is invalid for a collection of "
          "systems based on " + std::to_string(unique_topology_count) + " unique topologies.",
          "PhaseSpaceSynthesis", "getSystemIndicesByTopology");
  }
  return shared_topology_instance_bounds.readHost(topology_index + 1) -
         shared_topology_instance_bounds.readHost(topology_index);
}

//-------------------------------------------------------------------------------------------------
std::vector<int> PhaseSpaceSynthesis::getSystemIndicesByTopology(const int topology_index) const {
  if (topology_index < 0 || topology_index >= unique_topology_count) {
    rtErr("Topology index " + std::to_string(topology_index) + " is invalid for a collection of "
          "systems based on " + std::to_string(unique_topology_count) + " unique topologies.",
          "PhaseSpaceSynthesis", "getSystemIndicesByTopology");
  }
  const int llim = shared_topology_instance_bounds.readHost(topology_index);
  const int hlim = shared_topology_instance_bounds.readHost(topology_index + 1);
  const int* sti_ptr = shared_topology_instances.data();
  std::vector<int> result(hlim - llim);
  for (int i = llim; i < hlim; i++) {
    result[i - llim] = sti_ptr[i];
  }
  return result;
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
  int *devc_atom_starts, *devc_atom_counts;
  int *devc_shared_topology_instances, *devc_shared_topology_instance_bounds;
  int *devc_unique_topology_reference, *devc_shared_topology_instance_index;
  llint *devc_boxvecs, *devc_alt_boxvecs;
  int *devc_boxvec_ovrf, *devc_alt_boxvec_ovrf;
  double *devc_umat, *devc_invu, *devc_boxdims, *devc_umat_alt, *devc_invu_alt, *devc_alt_boxdims;
  llint *devc_xcrd,  *devc_ycrd,  *devc_zcrd,  *devc_xvel,  *devc_yvel,  *devc_zvel;
  llint *devc_xfrc,  *devc_yfrc,  *devc_zfrc,  *devc_xalt,  *devc_yalt,  *devc_zalt;
  llint *devc_vxalt, *devc_vyalt, *devc_vzalt, *devc_fxalt, *devc_fyalt, *devc_fzalt;
  int *devc_xcrd_ovrf,  *devc_ycrd_ovrf,  *devc_zcrd_ovrf;
  int *devc_xvel_ovrf,  *devc_yvel_ovrf,  *devc_zvel_ovrf;
  int *devc_xfrc_ovrf,  *devc_yfrc_ovrf,  *devc_zfrc_ovrf;
  int *devc_xalt_ovrf,  *devc_yalt_ovrf,  *devc_zalt_ovrf;
  int *devc_vxalt_ovrf, *devc_vyalt_ovrf, *devc_vzalt_ovrf;
  int *devc_fxalt_ovrf, *devc_fyalt_ovrf, *devc_fzalt_ovrf;

  bool problem = false;
#  ifdef STORMM_USE_CUDA
  // Obtain pointers to system indexing and unit cell dimensions stored on the host.  The pointers
  // are valid on the device.
  problem = (problem || cudaHostGetDevicePointer((void **)&devc_atom_starts,
                                                 (void *)atom_starts.data(), 0) != cudaSuccess);
  problem = (problem || cudaHostGetDevicePointer((void **)&devc_atom_counts,
                                                 (void *)atom_counts.data(), 0) != cudaSuccess);
  problem = (problem ||
             cudaHostGetDevicePointer((void **)&devc_shared_topology_instances,
                                      (void *)shared_topology_instances.data(), 0) != cudaSuccess);
  problem = (problem ||
             cudaHostGetDevicePointer((void **)&devc_shared_topology_instance_bounds,
                                      (void *)shared_topology_instance_bounds.data(), 0) !=
             cudaSuccess);
  problem = (problem || cudaHostGetDevicePointer((void **)&devc_boxvecs,
                                                 (void *)box_vectors.data(), 0) != cudaSuccess);
  problem = (problem ||
             cudaHostGetDevicePointer((void **)&devc_boxvec_ovrf,
                                      (void *)box_vector_overflow.data(), 0) != cudaSuccess);
  problem = (problem ||
             cudaHostGetDevicePointer((void **)&devc_umat,
                                      (void *)box_space_transforms.data(), 0) != cudaSuccess);
  problem = (problem ||  
             cudaHostGetDevicePointer((void **)&devc_invu,
                                      (void *)inverse_transforms.data(), 0) != cudaSuccess);
  problem = (problem || cudaHostGetDevicePointer((void **)&devc_boxdims,
                                                 (void *)box_dimensions.data(), 0) != cudaSuccess);
  problem = (problem ||
             cudaHostGetDevicePointer((void **)&devc_alt_boxvecs,
                                      (void *)alt_box_vectors.data(), 0) != cudaSuccess);
  problem = (problem ||
             cudaHostGetDevicePointer((void **)&devc_alt_boxvec_ovrf,
                                      (void *)alt_box_vector_overflow.data(), 0) != cudaSuccess);
  problem = (problem ||
             cudaHostGetDevicePointer((void **)&devc_umat_alt,
                                      (void *)alt_box_transforms.data(), 0) != cudaSuccess);
  problem = (problem ||  
             cudaHostGetDevicePointer((void **)&devc_invu_alt,
                                      (void *)alt_inverse_transforms.data(), 0) != cudaSuccess);
  problem = (problem ||
             cudaHostGetDevicePointer((void **)&devc_alt_boxdims,
                                      (void *)alt_box_dimensions.data(), 0) != cudaSuccess);
  problem = (problem ||
             cudaHostGetDevicePointer((void **)&devc_unique_topology_reference,
                                      (void *)unique_topology_reference.data(), 0) != cudaSuccess);
  problem = (problem ||
             cudaHostGetDevicePointer((void **)&devc_shared_topology_instance_index,
                                      (void *)shared_topology_instance_index.data(), 0) !=
             cudaSuccess);
  
  // Obtain pointers to the host-data for atomic positions.
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

  // Obtain pointers to the host-data for particle velocities.
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

  // Obtain pointers to the host-data for forces acting on all particles.
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

  // Obtain pointers to the host data for prior positions.
  problem = (problem ||
             cudaHostGetDevicePointer((void **)&devc_xalt,
                                      (void *)x_alt_coordinates.data(), 0) != cudaSuccess);
  problem = (problem ||
             cudaHostGetDevicePointer((void **)&devc_yalt,
                                      (void *)y_alt_coordinates.data(), 0) != cudaSuccess);
  problem = (problem ||
             cudaHostGetDevicePointer((void **)&devc_zalt,
                                      (void *)z_alt_coordinates.data(), 0) != cudaSuccess);
  problem = (problem ||
             cudaHostGetDevicePointer((void **)&devc_xalt_ovrf,
                                      (void *)x_alt_coord_overflow.data(), 0) != cudaSuccess);
  problem = (problem ||
             cudaHostGetDevicePointer((void **)&devc_yalt_ovrf,
                                      (void *)y_alt_coord_overflow.data(), 0) != cudaSuccess);
  problem = (problem ||
             cudaHostGetDevicePointer((void **)&devc_zalt_ovrf,
                                      (void *)z_alt_coord_overflow.data(), 0) != cudaSuccess);

  // Obtain pointers to the host data for prior velocities.
  problem = (problem ||
             cudaHostGetDevicePointer((void **)&devc_vxalt,
                                      (void *)x_alt_velocities.data(), 0) != cudaSuccess);
  problem = (problem ||
             cudaHostGetDevicePointer((void **)&devc_vyalt,
                                      (void *)y_alt_velocities.data(), 0) != cudaSuccess);
  problem = (problem ||
             cudaHostGetDevicePointer((void **)&devc_vzalt,
                                      (void *)z_alt_velocities.data(), 0) != cudaSuccess);
  problem = (problem ||
             cudaHostGetDevicePointer((void **)&devc_vxalt_ovrf,
                                      (void *)x_alt_velocity_overflow.data(), 0) != cudaSuccess);
  problem = (problem ||
             cudaHostGetDevicePointer((void **)&devc_vyalt_ovrf,
                                      (void *)y_alt_velocity_overflow.data(), 0) != cudaSuccess);
  problem = (problem ||
             cudaHostGetDevicePointer((void **)&devc_vzalt_ovrf,
                                      (void *)z_alt_velocity_overflow.data(), 0) != cudaSuccess);

  // Obtain pointers to the host data for prior forces.
  problem = (problem ||
             cudaHostGetDevicePointer((void **)&devc_fxalt,
                                      (void *)x_alt_forces.data(), 0) != cudaSuccess);
  problem = (problem ||
             cudaHostGetDevicePointer((void **)&devc_fyalt,
                                      (void *)y_alt_forces.data(), 0) != cudaSuccess);
  problem = (problem ||
             cudaHostGetDevicePointer((void **)&devc_fzalt,
                                      (void *)z_alt_forces.data(), 0) != cudaSuccess);
  problem = (problem ||
             cudaHostGetDevicePointer((void **)&devc_fxalt_ovrf,
                                      (void *)x_alt_force_overflow.data(), 0) != cudaSuccess);
  problem = (problem ||
             cudaHostGetDevicePointer((void **)&devc_fyalt_ovrf,
                                      (void *)y_alt_force_overflow.data(), 0) != cudaSuccess);
  problem = (problem ||
             cudaHostGetDevicePointer((void **)&devc_fzalt_ovrf,
                                      (void *)z_alt_force_overflow.data(), 0) != cudaSuccess);
#  endif // STORMM_USE_CUDA
  if (problem) {
    rtErr("Unable to get device-mapped pointers to host memory.  Types of memory used in each "
          "array: " + std::string(int_data.getLabel().name) + " " + int_data.getLabel().format +
          ", " + std::string(llint_data.getLabel().name) + " " + llint_data.getLabel().format +
          ", " + std::string(double_data.getLabel().name) + " " + double_data.getLabel().format +
          ".", "PhaseSpaceSynthesis", "deviceViewToHostData");
  }
  switch (orientation) {
  case CoordinateCycle::PRIMARY:
    return PsSynthesisWriter(system_count, unique_topology_count, unit_cell, heat_bath_kind,
                             piston_kind, time_step, devc_atom_starts, devc_atom_counts,
                             devc_shared_topology_instances, devc_shared_topology_instance_bounds,
                             devc_unique_topology_reference, devc_shared_topology_instance_index,
                             globalpos_scale, localpos_scale, velocity_scale, force_scale,
                             globalpos_scale_bits, localpos_scale_bits, velocity_scale_bits,
                             force_scale_bits, devc_boxvecs, devc_boxvec_ovrf, devc_umat,
                             devc_invu, devc_boxdims, devc_alt_boxvecs, devc_alt_boxvec_ovrf,
                             devc_umat_alt, devc_invu_alt, devc_alt_boxdims, devc_xcrd, devc_ycrd,
                             devc_zcrd, devc_xcrd_ovrf, devc_ycrd_ovrf, devc_zcrd_ovrf, devc_xvel,
                             devc_yvel, devc_zvel, devc_xvel_ovrf, devc_yvel_ovrf, devc_zvel_ovrf,
                             devc_xfrc, devc_yfrc, devc_zfrc, devc_xfrc_ovrf, devc_yfrc_ovrf,
                             devc_zfrc_ovrf, devc_xalt, devc_yalt, devc_zalt, devc_xalt_ovrf,
                             devc_yalt_ovrf, devc_zalt_ovrf, devc_vxalt, devc_vyalt, devc_vzalt,
                             devc_vxalt_ovrf, devc_vyalt_ovrf, devc_vzalt_ovrf, devc_fxalt,
                             devc_fyalt, devc_fzalt, devc_fxalt_ovrf, devc_fyalt_ovrf,
                             devc_fzalt_ovrf);
  case CoordinateCycle::ALTERNATE:
    return PsSynthesisWriter(system_count, unique_topology_count, unit_cell, heat_bath_kind,
                             piston_kind, time_step, devc_atom_starts, devc_atom_counts,
                             devc_shared_topology_instances, devc_shared_topology_instance_bounds,
                             devc_unique_topology_reference, devc_shared_topology_instance_index,
                             globalpos_scale, localpos_scale, velocity_scale, force_scale,
                             globalpos_scale_bits, localpos_scale_bits, velocity_scale_bits,
                             force_scale_bits, devc_alt_boxvecs, devc_alt_boxvec_ovrf,
                             devc_umat_alt, devc_invu_alt, devc_alt_boxdims, devc_boxvecs,
                             devc_boxvec_ovrf, devc_umat, devc_invu, devc_boxdims, devc_xalt,
                             devc_yalt, devc_zalt, devc_xalt_ovrf, devc_yalt_ovrf, devc_zalt_ovrf,
                             devc_vxalt, devc_vyalt, devc_vzalt, devc_vxalt_ovrf, devc_vyalt_ovrf,
                             devc_vzalt_ovrf, devc_fxalt, devc_fyalt, devc_fzalt, devc_fxalt_ovrf,
                             devc_fyalt_ovrf, devc_fzalt_ovrf, devc_xcrd, devc_ycrd, devc_zcrd,
                             devc_xcrd_ovrf, devc_ycrd_ovrf, devc_zcrd_ovrf, devc_xvel, devc_yvel,
                             devc_zvel, devc_xvel_ovrf, devc_yvel_ovrf, devc_zvel_ovrf, devc_xfrc,
                             devc_yfrc, devc_zfrc, devc_xfrc_ovrf, devc_yfrc_ovrf, devc_zfrc_ovrf);
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
void PhaseSpaceSynthesis::upload() {
  atom_starts.upload();
  atom_counts.upload();
  llint_data.upload();
  int_data.upload();
  double_data.upload();
}

//-------------------------------------------------------------------------------------------------
void PhaseSpaceSynthesis::upload(const TrajectoryKind kind) {
  const size_t atom_stride = roundUp<size_t>(atom_starts.readHost(system_count - 1) +
                                             atom_counts.readHost(system_count - 1), warp_size_zu);
  const size_t xfrm_stride = system_count * roundUp<size_t>(9, warp_size_zu);
  switch (kind) {
  case TrajectoryKind::POSITIONS:
    llint_data.upload(0, (6LLU * atom_stride) + (2LLU * xfrm_stride));
    if (globalpos_scale_bits > globalpos_scale_nonoverflow_bits) {
      int_data.upload(0, (6LLU * atom_stride) + (2LLU * xfrm_stride));
    }
    double_data.upload();
    break;
  case TrajectoryKind::VELOCITIES:
    llint_data.upload(( 6LLU * atom_stride) + (2LLU * xfrm_stride), 6LLU * atom_stride);
    if (velocity_scale_bits > velocity_scale_nonoverflow_bits) {
      int_data.upload(( 6LLU * atom_stride) + (2LLU * xfrm_stride), 6LLU * atom_stride);
    }
    break;
  case TrajectoryKind::FORCES:
    llint_data.upload((12LLU * atom_stride) + (2LLU * xfrm_stride), 6LLU * atom_stride);
    if (force_scale_bits > force_scale_nonoverflow_bits) {
      int_data.upload((12LLU * atom_stride) + (2LLU * xfrm_stride), 6LLU * atom_stride);
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
}

//-------------------------------------------------------------------------------------------------
void PhaseSpaceSynthesis::download(const TrajectoryKind kind) {
  const size_t atom_stride = roundUp<size_t>(atom_starts.readHost(system_count - 1) +
                                             atom_counts.readHost(system_count - 1), warp_size_zu);
  const size_t xfrm_stride = system_count * roundUp<size_t>(9, warp_size_zu);
  switch (kind) {
  case TrajectoryKind::POSITIONS:
    llint_data.download(0, (6LLU * atom_stride) + (2LLU * xfrm_stride));
    if (globalpos_scale_bits > globalpos_scale_nonoverflow_bits) {
      int_data.download(0, (6LLU * atom_stride) + (2LLU * xfrm_stride));
    }
    double_data.download();
    box_vectors.download();
    break;
  case TrajectoryKind::VELOCITIES:
    llint_data.download(( 6LLU * atom_stride) + (2LLU * xfrm_stride), 6LLU * atom_stride);
    if (velocity_scale_bits > velocity_scale_nonoverflow_bits) {
      int_data.download(( 6LLU * atom_stride) + (2LLU * xfrm_stride), 6LLU * atom_stride);
    }
    break;
  case TrajectoryKind::FORCES:
    llint_data.download((12LLU * atom_stride) + (2LLU * xfrm_stride), 6LLU * atom_stride);
    if (force_scale_bits > force_scale_nonoverflow_bits) {
      int_data.download((12LLU * atom_stride) + (2LLU * xfrm_stride), 6LLU * atom_stride);
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
                                        const HybridTargetLevel origin,
                                        const HybridTargetLevel destination,
                                        const GpuDetails &gpu) const {
  PhaseSpaceWriter psw = ps->data();
  if (atom_counts.readHost(index) != psw.natom) {
    rtErr("A PhaseSpace object sized for " + std::to_string(psw.natom) + " atoms is not prepared "
          "to accept a system of " + std::to_string(atom_counts.readHost(index)) + " atoms from "
          "this synthesis.", "PhaseSpaceSynthesis", "extractPhaseSpace");
  }
  int atom_offset = atom_starts.readHost(index);
  int mtrx_offset = index * roundUp(9, warp_size_int);
  int bdim_offset = index * roundUp(6, warp_size_int);
  const double crd_deflation = inverse_globalpos_scale;
  const double vel_deflation = inverse_velocity_scale;
  const double frc_deflation = inverse_force_scale;
  switch (origin) {
  case HybridTargetLevel::HOST:
    switch (destination) {
    case HybridTargetLevel::HOST:
      {
        // The phase space object should be extracted in a manner that orients its "present"
        // towards the current coordinates from the phase space synthesis, wherever the larger
        // synthesis sits in its own cycle.
        const PsSynthesisReader tpsr = this->data(cycle_position, origin);
        hostInt95ToDouble(psw.xcrd, psw.ycrd, psw.zcrd, &tpsr.xcrd[atom_offset],
                          &tpsr.xcrd_ovrf[atom_offset], &tpsr.ycrd[atom_offset],
                          &tpsr.ycrd_ovrf[atom_offset], &tpsr.zcrd[atom_offset],
                          &tpsr.zcrd_ovrf[atom_offset], psw.natom, crd_deflation);
        hostInt95ToDouble(psw.xalt, psw.yalt, psw.zalt, &tpsr.xalt[atom_offset],
                          &tpsr.xalt_ovrf[atom_offset], &tpsr.yalt[atom_offset],
                          &tpsr.yalt_ovrf[atom_offset], &tpsr.zalt[atom_offset],
                          &tpsr.zalt_ovrf[atom_offset], psw.natom, crd_deflation);
        hostInt95ToDouble(psw.xvel, psw.yvel, psw.zvel, &tpsr.xvel[atom_offset],
                          &tpsr.xvel_ovrf[atom_offset], &tpsr.yvel[atom_offset],
                          &tpsr.yvel_ovrf[atom_offset], &tpsr.zvel[atom_offset],
                          &tpsr.zvel_ovrf[atom_offset], psw.natom, vel_deflation);
        hostInt95ToDouble(psw.vxalt, psw.vyalt, psw.vzalt, &tpsr.vxalt[atom_offset],
                          &tpsr.vxalt_ovrf[atom_offset], &tpsr.vyalt[atom_offset],
                          &tpsr.vyalt_ovrf[atom_offset], &tpsr.vzalt[atom_offset],
                          &tpsr.vzalt_ovrf[atom_offset], psw.natom, vel_deflation);
        hostInt95ToDouble(psw.xfrc, psw.yfrc, psw.zfrc, &tpsr.xfrc[atom_offset],
                          &tpsr.xfrc_ovrf[atom_offset], &tpsr.yfrc[atom_offset],
                          &tpsr.yfrc_ovrf[atom_offset], &tpsr.zfrc[atom_offset],
                          &tpsr.zfrc_ovrf[atom_offset], psw.natom, frc_deflation);
        hostInt95ToDouble(psw.fxalt, psw.fyalt, psw.fzalt, &tpsr.fxalt[atom_offset],
                          &tpsr.fxalt_ovrf[atom_offset], &tpsr.fyalt[atom_offset],
                          &tpsr.fyalt_ovrf[atom_offset], &tpsr.fzalt[atom_offset],
                          &tpsr.fzalt_ovrf[atom_offset], psw.natom, frc_deflation);
        for (int i = 0; i < 9; i++) {
          psw.umat[i]     = tpsr.umat[mtrx_offset + i];
          psw.invu[i]     = tpsr.invu[mtrx_offset + i];
          psw.umat_alt[i] = tpsr.umat_alt[mtrx_offset + i];
          psw.invu_alt[i] = tpsr.invu_alt[mtrx_offset + i];
        }
        for (int i = 0; i < 6; i++) {
          psw.boxdim[i]     = tpsr.boxdims[bdim_offset + i];
          psw.boxdim_alt[i] = tpsr.alt_boxdims[bdim_offset + i];
        }
      }
      break;
#ifdef STORMM_USE_HPC
    case HybridTargetLevel::DEVICE:
      {
        PhaseSpaceWriter psw = ps->data(destination);
        extractSystem(&psw, index, gpu, origin);
      }
      break;
#endif
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    {
      PhaseSpaceWriter psw = ps->data(destination);
      extractSystem(&psw, index, gpu, origin);
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
                                                       const CoordinateCycle orientation,
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
  std::vector<double> umat_buffer, invu_buffer, bdim_buffer;
  const size_t xfrm_offset = static_cast<size_t>(index) * roundUp<size_t>(9, warp_size_int);
  const size_t bdim_offset = static_cast<size_t>(index) * roundUp<size_t>(6, warp_size_int);
  switch (tier) {
  case HybridTargetLevel::HOST:
    switch (trajkind) {
    case TrajectoryKind::POSITIONS:
      switch (orientation) {
      case CoordinateCycle::ALTERNATE:
        xbuffer = x_alt_coordinates.readHost(astart, rsw.natom);
        ybuffer = y_alt_coordinates.readHost(astart, rsw.natom);
        zbuffer = z_alt_coordinates.readHost(astart, rsw.natom);
        umat_buffer = alt_box_transforms.readHost(xfrm_offset, 9);
        invu_buffer = alt_inverse_transforms.readHost(xfrm_offset, 9);
        bdim_buffer = alt_box_dimensions.readHost(bdim_offset, 6);
        break;
      case CoordinateCycle::PRIMARY:
        xbuffer = x_coordinates.readHost(astart, rsw.natom);
        ybuffer = y_coordinates.readHost(astart, rsw.natom);
        zbuffer = z_coordinates.readHost(astart, rsw.natom);
        umat_buffer = box_space_transforms.readHost(xfrm_offset, 9);
        invu_buffer = inverse_transforms.readHost(xfrm_offset, 9);
        bdim_buffer = box_dimensions.readHost(bdim_offset, 6);
        break;
      }
      break;
    case TrajectoryKind::VELOCITIES:
      switch (orientation) {
      case CoordinateCycle::ALTERNATE:
        xbuffer = x_alt_velocities.readHost(astart, rsw.natom);
        ybuffer = y_alt_velocities.readHost(astart, rsw.natom);
        zbuffer = z_alt_velocities.readHost(astart, rsw.natom);
        break;
      case CoordinateCycle::PRIMARY:
        xbuffer = x_velocities.readHost(astart, rsw.natom);
        ybuffer = y_velocities.readHost(astart, rsw.natom);
        zbuffer = z_velocities.readHost(astart, rsw.natom);
        break;
      }
      break;
    case TrajectoryKind::FORCES:
      switch (orientation) {
      case CoordinateCycle::ALTERNATE:
        xbuffer = x_alt_forces.readHost(astart, rsw.natom);
        ybuffer = y_alt_forces.readHost(astart, rsw.natom);
        zbuffer = z_alt_forces.readHost(astart, rsw.natom);
        break;
      case CoordinateCycle::PRIMARY:
        xbuffer = x_forces.readHost(astart, rsw.natom);
        ybuffer = y_forces.readHost(astart, rsw.natom);
        zbuffer = z_forces.readHost(astart, rsw.natom);
        break;
      }
      break;
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    switch (trajkind) {
    case TrajectoryKind::POSITIONS:
      switch (orientation) {
      case CoordinateCycle::ALTERNATE:
        xbuffer = x_alt_coordinates.readDevice(astart, rsw.natom);
        ybuffer = y_alt_coordinates.readDevice(astart, rsw.natom);
        zbuffer = z_alt_coordinates.readDevice(astart, rsw.natom);
        umat_buffer = alt_box_transforms.readDevice(xfrm_offset, 9);
        invu_buffer = alt_inverse_transforms.readDevice(xfrm_offset, 9);
        bdim_buffer = alt_box_dimensions.readDevice(bdim_offset, 6);
        break;
      case CoordinateCycle::PRIMARY:
        xbuffer = x_coordinates.readDevice(astart, rsw.natom);
        ybuffer = y_coordinates.readDevice(astart, rsw.natom);
        zbuffer = z_coordinates.readDevice(astart, rsw.natom);
        umat_buffer = box_space_transforms.readDevice(xfrm_offset, 9);
        invu_buffer = inverse_transforms.readDevice(xfrm_offset, 9);
        bdim_buffer = box_dimensions.readDevice(bdim_offset, 6);
        break;
      }
      break;
    case TrajectoryKind::VELOCITIES:
      switch (orientation) {
      case CoordinateCycle::ALTERNATE:
        xbuffer = x_alt_velocities.readDevice(astart, rsw.natom);
        ybuffer = y_alt_velocities.readDevice(astart, rsw.natom);
        zbuffer = z_alt_velocities.readDevice(astart, rsw.natom);
        break;
      case CoordinateCycle::PRIMARY:
        xbuffer = x_velocities.readDevice(astart, rsw.natom);
        ybuffer = y_velocities.readDevice(astart, rsw.natom);
        zbuffer = z_velocities.readDevice(astart, rsw.natom);
        break;
      }
      break;
    case TrajectoryKind::FORCES:
      switch (orientation) {
      case CoordinateCycle::ALTERNATE:
        xbuffer = x_alt_forces.readDevice(astart, rsw.natom);
        ybuffer = y_alt_forces.readDevice(astart, rsw.natom);
        zbuffer = z_alt_forces.readDevice(astart, rsw.natom);
        break;
      case CoordinateCycle::PRIMARY:
        xbuffer = x_forces.readDevice(astart, rsw.natom);
        ybuffer = y_forces.readDevice(astart, rsw.natom);
        zbuffer = z_forces.readDevice(astart, rsw.natom);
        break;
      }
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
  if (trajkind == TrajectoryKind::POSITIONS) {
    for (int i = 0; i < 9; i++) {
      rsw.umat[i] = umat_buffer[i];
      rsw.invu[i] = invu_buffer[i];
    }
    for (int i = 0; i < 6; i++) {
      rsw.boxdim[i] = bdim_buffer[i];
    }
  }
  if (globalpos_scale_bits > globalpos_scale_nonoverflow_bits) {
    std::vector<int> x_ovrf_buffer, y_ovrf_buffer, z_ovrf_buffer;
    switch (tier) {
    case HybridTargetLevel::HOST:
      switch (trajkind) {
      case TrajectoryKind::POSITIONS:
        switch (orientation) {
        case CoordinateCycle::ALTERNATE:
          x_ovrf_buffer = x_alt_coord_overflow.readHost(astart, rsw.natom);
          y_ovrf_buffer = y_alt_coord_overflow.readHost(astart, rsw.natom);
          z_ovrf_buffer = z_alt_coord_overflow.readHost(astart, rsw.natom);
          break;
        case CoordinateCycle::PRIMARY:
          x_ovrf_buffer = x_coordinate_overflow.readHost(astart, rsw.natom);
          y_ovrf_buffer = y_coordinate_overflow.readHost(astart, rsw.natom);
          z_ovrf_buffer = z_coordinate_overflow.readHost(astart, rsw.natom);
          break;
        }
        break;
      case TrajectoryKind::VELOCITIES:
        switch (orientation) {
        case CoordinateCycle::ALTERNATE:
          x_ovrf_buffer = x_alt_velocity_overflow.readHost(astart, rsw.natom);
          y_ovrf_buffer = y_alt_velocity_overflow.readHost(astart, rsw.natom);
          z_ovrf_buffer = z_alt_velocity_overflow.readHost(astart, rsw.natom);
          break;
        case CoordinateCycle::PRIMARY:
          x_ovrf_buffer = x_velocity_overflow.readHost(astart, rsw.natom);
          y_ovrf_buffer = y_velocity_overflow.readHost(astart, rsw.natom);
          z_ovrf_buffer = z_velocity_overflow.readHost(astart, rsw.natom);
          break;
        }
        break;
      case TrajectoryKind::FORCES:
        switch (orientation) {
        case CoordinateCycle::ALTERNATE:
          x_ovrf_buffer = x_alt_force_overflow.readHost(astart, rsw.natom);
          y_ovrf_buffer = y_alt_force_overflow.readHost(astart, rsw.natom);
          z_ovrf_buffer = z_alt_force_overflow.readHost(astart, rsw.natom);
          break;
        case CoordinateCycle::PRIMARY:
          x_ovrf_buffer = x_force_overflow.readHost(astart, rsw.natom);
          y_ovrf_buffer = y_force_overflow.readHost(astart, rsw.natom);
          z_ovrf_buffer = z_force_overflow.readHost(astart, rsw.natom);
          break;
        }
        break;
      }
      break;
#ifdef STORMM_USE_HPC
    case HybridTargetLevel::DEVICE:
      switch (trajkind) {
      case TrajectoryKind::POSITIONS:
        switch (orientation) {
        case CoordinateCycle::ALTERNATE:
          x_ovrf_buffer = x_alt_coord_overflow.readDevice(astart, rsw.natom);
          y_ovrf_buffer = y_alt_coord_overflow.readDevice(astart, rsw.natom);
          z_ovrf_buffer = z_alt_coord_overflow.readDevice(astart, rsw.natom);
          break;
        case CoordinateCycle::PRIMARY:
          x_ovrf_buffer = x_coordinate_overflow.readDevice(astart, rsw.natom);
          y_ovrf_buffer = y_coordinate_overflow.readDevice(astart, rsw.natom);
          z_ovrf_buffer = z_coordinate_overflow.readDevice(astart, rsw.natom);
          break;
        }
        break;
      case TrajectoryKind::VELOCITIES:
        switch (orientation) {
        case CoordinateCycle::ALTERNATE:
          x_ovrf_buffer = x_alt_velocity_overflow.readDevice(astart, rsw.natom);
          y_ovrf_buffer = y_alt_velocity_overflow.readDevice(astart, rsw.natom);
          z_ovrf_buffer = z_alt_velocity_overflow.readDevice(astart, rsw.natom);
          break;
        case CoordinateCycle::PRIMARY:
          x_ovrf_buffer = x_velocity_overflow.readDevice(astart, rsw.natom);
          y_ovrf_buffer = y_velocity_overflow.readDevice(astart, rsw.natom);
          z_ovrf_buffer = z_velocity_overflow.readDevice(astart, rsw.natom);
          break;
        }
        break;
      case TrajectoryKind::FORCES:
        switch (orientation) {
        case CoordinateCycle::ALTERNATE:
          x_ovrf_buffer = x_alt_force_overflow.readDevice(astart, rsw.natom);
          y_ovrf_buffer = y_alt_force_overflow.readDevice(astart, rsw.natom);
          z_ovrf_buffer = z_alt_force_overflow.readDevice(astart, rsw.natom);
          break;
        case CoordinateCycle::PRIMARY:
          x_ovrf_buffer = x_force_overflow.readDevice(astart, rsw.natom);
          y_ovrf_buffer = y_force_overflow.readDevice(astart, rsw.natom);
          z_ovrf_buffer = z_force_overflow.readDevice(astart, rsw.natom);
          break;
        }
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
CoordinateFrame PhaseSpaceSynthesis::exportCoordinates(const int index,
                                                       const TrajectoryKind trajkind,
                                                       const HybridTargetLevel tier) const {
  return exportCoordinates(index, cycle_position, trajkind, tier);
}

//-------------------------------------------------------------------------------------------------
#ifdef STORMM_USE_HPC
void PhaseSpaceSynthesis::initializeForces(const CoordinateCycle orientation,
                                           const GpuDetails &gpu,
                                           const HybridTargetLevel tier, const int index)
#else
void PhaseSpaceSynthesis::initializeForces(const CoordinateCycle orientation, const int index)
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
      llint* xptr;
      llint* yptr;
      llint* zptr;
      int* x_ovrf_ptr;
      int* y_ovrf_ptr;
      int* z_ovrf_ptr;
      switch (orientation) {
      case CoordinateCycle::ALTERNATE:
        xptr = x_alt_forces.data();
        yptr = y_alt_forces.data();
        zptr = z_alt_forces.data();
        x_ovrf_ptr = x_alt_force_overflow.data();
        y_ovrf_ptr = y_alt_force_overflow.data();
        z_ovrf_ptr = z_alt_force_overflow.data();
        break;
      case CoordinateCycle::PRIMARY:
        xptr = x_forces.data();
        yptr = y_forces.data();
        zptr = z_forces.data();
        x_ovrf_ptr = x_force_overflow.data();
        y_ovrf_ptr = y_force_overflow.data();
        z_ovrf_ptr = z_force_overflow.data();
        break;
      }
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
      PsSynthesisWriter dptr = data(orientation, HybridTargetLevel::DEVICE);
      psyInitializeForces(&dptr, index, gpu);
    }
    break;
  }
#endif
}

//-------------------------------------------------------------------------------------------------
#ifdef STORMM_USE_HPC
void PhaseSpaceSynthesis::initializeForces(const GpuDetails &gpu, const HybridTargetLevel tier,
                                           const int index) {
  initializeForces(cycle_position, gpu, tier, index);
}
#else
void PhaseSpaceSynthesis::initializeForces(const int index) {
  initializeForces(cycle_position, index);
}
#endif

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
      llint* xalt_ptr = x_forces.data();
      llint* yalt_ptr = y_forces.data();
      llint* zalt_ptr = z_forces.data();
      int* xalt_ovrf_ptr = x_force_overflow.data();
      int* yalt_ovrf_ptr = y_force_overflow.data();
      int* zalt_ovrf_ptr = z_force_overflow.data();
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
          xalt_ptr[j] = xfrc_ptr[j];
          yalt_ptr[j] = yfrc_ptr[j];
          zalt_ptr[j] = zfrc_ptr[j];
          xvel_ptr[j] = 0LL;
          yvel_ptr[j] = 0LL;
          zvel_ptr[j] = 0LL;
          if (force_scale_bits > force_scale_nonoverflow_bits) {
            xalt_ovrf_ptr[j] = xfrc_ovrf_ptr[j];
            yalt_ovrf_ptr[j] = yfrc_ovrf_ptr[j];
            zalt_ovrf_ptr[j] = zfrc_ovrf_ptr[j];
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
void PhaseSpaceSynthesis::printTrajectory(const std::vector<int> &system_indices,
                                          const std::string &file_name, const double current_time,
                                          const CoordinateFileKind output_kind,
                                          const PrintSituation expectation) const {

  // Bail out if there are no frames to print
  const size_t nframe = system_indices.size();
  if (nframe == 0LLU) {
    return;
  }
  
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
      const bool fi_exists = (getDrivePathType(file_name) == DrivePathType::FILE);
      foutp = openOutputFile(file_name, actual_expectation, "Open an output trajectory for "
                             "writing PhaseSpaceSynthesis contents.", DataFormat::ASCII);

      // Initialize the trajectory if this is the first time printing it
      if (fi_exists == false || actual_expectation == PrintSituation::OVERWRITE ||
          actual_expectation == PrintSituation::OPEN_NEW) {
        sprintf(buffer, "Generated by STORMM\n");
        foutp.write(buffer, strlen(buffer));
      }
    }
    break;
  case CoordinateFileKind::AMBER_NETCDF:
    if (expectation == PrintSituation::UNKNOWN) {
      actual_expectation = PrintSituation::OPEN_NEW;
    }
    foutp = openOutputFile(file_name, actual_expectation, "Open an SD file archive for writing "
                           "PhaseSpaceSynthesis contents.", DataFormat::BINARY);
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
  case CoordinateFileKind::SDF:
    rtErr("Some information needed to print an SD file is not present in a this object alone.  "
          "The program must use the writeFrame() method and pass this object to one of its "
          "overloads.", "PhaseSpaceSynthesis", "printTrajectory");
    break;
  case CoordinateFileKind::UNKNOWN:
    rtErr("Printing request for unknown file type.", "PhaseSpaceSynthesis", "printTrajectory");
  }

  // Proceed frame by frame
  std::vector<double> tmp_xcrd, tmp_ycrd, tmp_zcrd, tmp_xvel, tmp_yvel, tmp_zvel;
  const PsSynthesisReader tpsr = this->data(cycle_position);
  for (size_t i = 0; i < nframe; i++) {
    const int fr_start = atom_starts.readHost(system_indices[i]);
    const int fr_end   = fr_start + atom_counts.readHost(system_indices[i]);
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
    case CoordinateFileKind::SDF:
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
    case CoordinateFileKind::SDF:
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
    break;
  case CoordinateFileKind::SDF:
  case CoordinateFileKind::UNKNOWN:
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void PhaseSpaceSynthesis::import(const PhaseSpaceReader &psr, const int system_index,
                                 const CoordinateCycle orientation,
                                 const HybridTargetLevel tier) {
  CoordinateCycle focus = orientation;
  import(psr.xcrd, psr.ycrd, psr.zcrd, psr.umat, psr.invu, psr.boxdim, system_index, 1.0,
         TrajectoryKind::POSITIONS, focus, tier);
  import(psr.xvel, psr.yvel, psr.zvel, nullptr, nullptr, nullptr, system_index, 1.0,
         TrajectoryKind::VELOCITIES, focus, tier);
  import(psr.xfrc, psr.yfrc, psr.zfrc, nullptr, nullptr, nullptr, system_index, 1.0,
         TrajectoryKind::FORCES, focus, tier);
  switch (orientation) {
  case CoordinateCycle::ALTERNATE:
    focus = CoordinateCycle::PRIMARY;
    break;
  case CoordinateCycle::PRIMARY:
    focus = CoordinateCycle::ALTERNATE;
    break;
  }
  import(psr.xalt, psr.yalt, psr.zalt, psr.umat_alt, psr.invu_alt, psr.boxdim_alt, system_index,
         1.0, TrajectoryKind::POSITIONS, focus, tier);
  import(psr.vxalt, psr.vyalt, psr.vzalt, nullptr, nullptr, nullptr, system_index, 1.0,
         TrajectoryKind::VELOCITIES, focus, tier);
  import(psr.fxalt, psr.fyalt, psr.fzalt, nullptr, nullptr, nullptr, system_index, 1.0,
         TrajectoryKind::FORCES, focus, tier);
}

//-------------------------------------------------------------------------------------------------
void PhaseSpaceSynthesis::import(const PhaseSpaceWriter &psw, const int system_index,
                                 const CoordinateCycle orientation,
                                 const HybridTargetLevel tier) {
  import(psw, system_index, orientation, tier);
}

//-------------------------------------------------------------------------------------------------
void PhaseSpaceSynthesis::import(const PhaseSpace &ps, const int system_index,
                                 const CoordinateCycle orientation,
                                 const HybridTargetLevel tier) {
  import(ps.data(tier), system_index, orientation, tier);
}

//-------------------------------------------------------------------------------------------------
void PhaseSpaceSynthesis::import(const CoordinateFrameReader &cfr, const int system_index,
                                 const TrajectoryKind kind, const CoordinateCycle orientation,
                                 const HybridTargetLevel tier) {
  import(cfr.xcrd, cfr.ycrd, cfr.zcrd, cfr.umat, cfr.invu, cfr.boxdim, system_index, 1.0, kind,
         orientation, tier);
}

//-------------------------------------------------------------------------------------------------
void PhaseSpaceSynthesis::import(const CoordinateFrameWriter &cfw, const int system_index,
                                 const TrajectoryKind kind, const CoordinateCycle orientation,
                                 const HybridTargetLevel tier) {
  import(cfw, system_index, kind, orientation, tier);
}

//-------------------------------------------------------------------------------------------------
void PhaseSpaceSynthesis::import(const CoordinateFrame &cf, const int system_index,
                                 const TrajectoryKind kind, const CoordinateCycle orientation,
                                 const HybridTargetLevel tier) {
  import(cf.data(tier), system_index, kind, orientation, tier);
}

//-------------------------------------------------------------------------------------------------
void PhaseSpaceSynthesis::allocate(const size_t atom_stride) {
  const size_t system_stride = roundUp<size_t>(system_count, warp_size_zu);
  const size_t xfrm_stride = system_count * roundUp<size_t>(9, warp_size_zu);
  const size_t dbl_xfrm_stride = 2LLU * xfrm_stride;
  int_data.resize((6LLU * system_stride) + (18LLU * atom_stride) + dbl_xfrm_stride +
                  roundUp(static_cast<size_t>(unique_topology_count + 1), warp_size_zu));
  atom_counts.setPointer(&int_data, 0LLU, system_count);
  atom_starts.setPointer(&int_data, system_stride, system_count);
  const size_t twosys   = (2LLU * system_stride);
  const size_t twosysx  = twosys + xfrm_stride;
  const size_t twosysdx = twosys + (2LLU * xfrm_stride);
  llint_data.resize((18LLU * atom_stride) + dbl_xfrm_stride);
  x_coordinates.setPointer(&llint_data,                                       0LLU, atom_stride);
  y_coordinates.setPointer(&llint_data,                                atom_stride, atom_stride);
  z_coordinates.setPointer(&llint_data,                         2LLU * atom_stride, atom_stride);
  x_alt_coordinates.setPointer(&llint_data,                     3LLU * atom_stride, atom_stride);
  y_alt_coordinates.setPointer(&llint_data,                     4LLU * atom_stride, atom_stride);
  z_alt_coordinates.setPointer(&llint_data,                     5LLU * atom_stride, atom_stride);
  box_vectors.setPointer(&llint_data,                           6LLU * atom_stride, xfrm_stride);  
  alt_box_vectors.setPointer(&llint_data,       xfrm_stride + ( 6LLU * atom_stride), xfrm_stride);
  x_velocities.setPointer(&llint_data,      dbl_xfrm_stride + ( 6LLU * atom_stride), atom_stride);
  y_velocities.setPointer(&llint_data,      dbl_xfrm_stride + ( 7LLU * atom_stride), atom_stride);
  z_velocities.setPointer(&llint_data,      dbl_xfrm_stride + ( 8LLU * atom_stride), atom_stride);
  x_alt_velocities.setPointer(&llint_data,  dbl_xfrm_stride + ( 9LLU * atom_stride), atom_stride);
  y_alt_velocities.setPointer(&llint_data,  dbl_xfrm_stride + (10LLU * atom_stride), atom_stride);
  z_alt_velocities.setPointer(&llint_data,  dbl_xfrm_stride + (11LLU * atom_stride), atom_stride);
  x_forces.setPointer(&llint_data,          dbl_xfrm_stride + (12LLU * atom_stride), atom_stride);
  y_forces.setPointer(&llint_data,          dbl_xfrm_stride + (13LLU * atom_stride), atom_stride);
  z_forces.setPointer(&llint_data,          dbl_xfrm_stride + (14LLU * atom_stride), atom_stride);
  x_alt_forces.setPointer(&llint_data,      dbl_xfrm_stride + (15LLU * atom_stride), atom_stride);
  y_alt_forces.setPointer(&llint_data,      dbl_xfrm_stride + (16LLU * atom_stride), atom_stride);
  z_alt_forces.setPointer(&llint_data,      dbl_xfrm_stride + (17LLU * atom_stride), atom_stride);
  x_coordinate_overflow.setPointer(&int_data,                             twosys, atom_stride);
  y_coordinate_overflow.setPointer(&int_data,            atom_stride  +   twosys, atom_stride);
  z_coordinate_overflow.setPointer(&int_data,   ( 2LLU * atom_stride) +   twosys, atom_stride);
  x_alt_coord_overflow.setPointer(&int_data,    ( 3LLU * atom_stride) +   twosys, atom_stride);
  y_alt_coord_overflow.setPointer(&int_data,    ( 4LLU * atom_stride) +   twosys, atom_stride);
  z_alt_coord_overflow.setPointer(&int_data,    ( 5LLU * atom_stride) +   twosys, atom_stride);
  box_vector_overflow.setPointer(&int_data,     ( 6LLU * atom_stride) +   twosys, xfrm_stride);
  alt_box_vector_overflow.setPointer(&int_data, ( 6LLU * atom_stride) +  twosysx, xfrm_stride);
  x_velocity_overflow.setPointer(&int_data,     ( 6LLU * atom_stride) + twosysdx, atom_stride);
  y_velocity_overflow.setPointer(&int_data,     ( 7LLU * atom_stride) + twosysdx, atom_stride);
  z_velocity_overflow.setPointer(&int_data,     ( 8LLU * atom_stride) + twosysdx, atom_stride);
  x_alt_velocity_overflow.setPointer(&int_data, ( 9LLU * atom_stride) + twosysdx, atom_stride);
  y_alt_velocity_overflow.setPointer(&int_data, (10LLU * atom_stride) + twosysdx, atom_stride);
  z_alt_velocity_overflow.setPointer(&int_data, (11LLU * atom_stride) + twosysdx, atom_stride);
  x_force_overflow.setPointer(&int_data,        (12LLU * atom_stride) + twosysdx, atom_stride);
  y_force_overflow.setPointer(&int_data,        (13LLU * atom_stride) + twosysdx, atom_stride);
  z_force_overflow.setPointer(&int_data,        (14LLU * atom_stride) + twosysdx, atom_stride);
  x_alt_force_overflow.setPointer(&int_data,    (15LLU * atom_stride) + twosysdx, atom_stride);
  y_alt_force_overflow.setPointer(&int_data,    (16LLU * atom_stride) + twosysdx, atom_stride);
  z_alt_force_overflow.setPointer(&int_data,    (17LLU * atom_stride) + twosysdx, atom_stride);
  double_data.resize(6LLU * xfrm_stride);
  box_space_transforms.setPointer(&double_data,                   0LLU, xfrm_stride);
  inverse_transforms.setPointer(&double_data,              xfrm_stride, xfrm_stride);
  box_dimensions.setPointer(&double_data,           2LLU * xfrm_stride, xfrm_stride);
  alt_box_transforms.setPointer(&double_data,       3LLU * xfrm_stride, xfrm_stride);
  alt_inverse_transforms.setPointer(&double_data,   4LLU * xfrm_stride, xfrm_stride);
  alt_box_dimensions.setPointer(&double_data,       5LLU * xfrm_stride, xfrm_stride);

  // Target the arrays detailing unique system instances to the back of the int_data array, as
  // the alignment of int_data and llint_data then works better for fixed-precision coordinates.
  // This is strictly organizational and has no effect on performance.
  shared_topology_instances.setPointer(&int_data, (18LLU * atom_stride) + twosysdx, system_count);
  shared_topology_instance_bounds.setPointer(&int_data, (18LLU * atom_stride) + twosysdx +
                                             system_stride, unique_topology_count + 1);
  unique_topology_reference.setPointer(&int_data, (18LLU * atom_stride) + twosysdx +
                                       (2LLU * system_stride), system_count);
  shared_topology_instance_index.setPointer(&int_data, (18LLU * atom_stride) + twosysdx +
                                       (3LLU * system_stride), system_count);
}

} // namespace trajectory
} // namespace stormm
