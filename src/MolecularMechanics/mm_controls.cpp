#include "Constants/scaling.h"
#include "Topology/atomgraph_enumerators.h"
#include "mm_controls.h"

namespace omni {
namespace mm {

using card::HybridKind;
using topology::UnitCellType;
  
//-------------------------------------------------------------------------------------------------
MolecularMechanicsControls::MolecularMechanicsControls(const double time_step_in,
                                                       const double rattle_tol_in,
                                                       const double initial_step_in) :
  step_number{0}, time_step{time_step_in}, rattle_tol{rattle_tol_in},
  initial_step{initial_step_in},
  vwu_progress{HybridKind::POINTER, "mm_vwu_counters"},
  nbwu_progress{HybridKind::POINTER, "mm_nbwu_counters"},
  pmewu_progress{HybridKind::POINTER, "mm_pmewu_counters"},
  gather_wu_progress{HybridKind::POINTER, "mm_gtwu_counters"},
  scatter_wu_progress{HybridKind::POINTER, "mm_scwu_counters"},
  all_reduce_wu_progress{HybridKind::POINTER, "mm_rdwu_counters"},
  progress_data{12 * warp_size_int, "work_unit_prog_data"}
{
  vwu_progress.setPointer(&progress_data,                             0, 2 * warp_size_int);
  nbwu_progress.setPointer(&progress_data,            2 * warp_size_int, 2 * warp_size_int);
  pmewu_progress.setPointer(&progress_data,           4 * warp_size_int, 2 * warp_size_int);
  gather_wu_progress.setPointer(&progress_data,       6 * warp_size_int, 2 * warp_size_int);
  scatter_wu_progress.setPointer(&progress_data,      8 * warp_size_int, 2 * warp_size_int);
  all_reduce_wu_progress.setPointer(&progress_data,  10 * warp_size_int, 2 * warp_size_int);
}

//-------------------------------------------------------------------------------------------------
MolecularMechanicsControls::MolecularMechanicsControls(const DynamicsControls &user_input) :
    MolecularMechanicsControls()
{
  time_step = user_input.getTimeStep();  
  rattle_tol = user_input.getRattleTolerance();
}

//-------------------------------------------------------------------------------------------------
MolecularMechanicsControls::MolecularMechanicsControls(const MinimizeControls &user_input) :
    MolecularMechanicsControls()
{
  initial_step = user_input.getInitialStep();
}

//-------------------------------------------------------------------------------------------------
MolecularMechanicsControls::
MolecularMechanicsControls(const MolecularMechanicsControls &original) :
  step_number{original.step_number},
  time_step{original.time_step},
  rattle_tol{original.rattle_tol},
  vwu_progress{original.vwu_progress},
  nbwu_progress{original.nbwu_progress},
  pmewu_progress{original.pmewu_progress},
  gather_wu_progress{original.gather_wu_progress},
  scatter_wu_progress{original.scatter_wu_progress},
  all_reduce_wu_progress{original.all_reduce_wu_progress},
  progress_data{original.progress_data}
{
  vwu_progress.swapTarget(&progress_data);
  nbwu_progress.swapTarget(&progress_data);
  pmewu_progress.swapTarget(&progress_data);
  gather_wu_progress.swapTarget(&progress_data);
  scatter_wu_progress.swapTarget(&progress_data);
  all_reduce_wu_progress.swapTarget(&progress_data);
}

//-------------------------------------------------------------------------------------------------
MolecularMechanicsControls&
MolecularMechanicsControls::operator=(const MolecularMechanicsControls &other) {

  // Guard against self assignment
  if (this == &other) {
    return *this;
  }
  step_number = other.step_number;
  time_step = other.time_step;
  rattle_tol = other.rattle_tol;
  vwu_progress = other.vwu_progress;
  nbwu_progress = other.nbwu_progress;
  pmewu_progress = other.pmewu_progress;
  gather_wu_progress = other.gather_wu_progress;
  scatter_wu_progress = other.scatter_wu_progress;
  all_reduce_wu_progress = other.all_reduce_wu_progress;
  progress_data = other.progress_data;

  // Repair pointers and return the result
  vwu_progress.swapTarget(&progress_data);
  nbwu_progress.swapTarget(&progress_data);
  pmewu_progress.swapTarget(&progress_data);
  gather_wu_progress.swapTarget(&progress_data);
  scatter_wu_progress.swapTarget(&progress_data);
  all_reduce_wu_progress.swapTarget(&progress_data);
  return *this;
}

//-------------------------------------------------------------------------------------------------
MolecularMechanicsControls::MolecularMechanicsControls(MolecularMechanicsControls &&original) :
  step_number{original.step_number},
  time_step{original.time_step},
  rattle_tol{original.rattle_tol},
  vwu_progress{std::move(original.vwu_progress)},
  nbwu_progress{std::move(original.nbwu_progress)},
  pmewu_progress{std::move(original.pmewu_progress)},
  gather_wu_progress{std::move(original.gather_wu_progress)},
  scatter_wu_progress{std::move(original.scatter_wu_progress)},
  all_reduce_wu_progress{std::move(original.all_reduce_wu_progress)},
  progress_data{std::move(original.progress_data)}
{}

//-------------------------------------------------------------------------------------------------
MolecularMechanicsControls&
MolecularMechanicsControls::operator=(MolecularMechanicsControls &&other) {

  // Guard against self assignment
  if (this == &other) {
    return *this;
  }
  step_number = other.step_number;
  time_step = other.time_step;
  rattle_tol = other.rattle_tol;
  vwu_progress = std::move(other.vwu_progress);
  nbwu_progress = std::move(other.nbwu_progress);
  pmewu_progress = std::move(other.pmewu_progress);
  gather_wu_progress = std::move(other.gather_wu_progress);
  scatter_wu_progress = std::move(other.scatter_wu_progress);
  all_reduce_wu_progress = std::move(other.all_reduce_wu_progress);
  progress_data = std::move(other.progress_data);
  return *this;
}

//-------------------------------------------------------------------------------------------------
MMControlKit<double> MolecularMechanicsControls::dpData(const HybridTargetLevel tier) {
  return MMControlKit<double>(step_number, time_step, rattle_tol, initial_step,
                              vwu_progress.data(tier), nbwu_progress.data(tier),
                              pmewu_progress.data(tier), gather_wu_progress.data(tier),
                              scatter_wu_progress.data(tier), all_reduce_wu_progress.data(tier));
}

//-------------------------------------------------------------------------------------------------
MMControlKit<float> MolecularMechanicsControls::spData(const HybridTargetLevel tier) {
  return MMControlKit<float>(step_number, time_step, rattle_tol, initial_step,
                             vwu_progress.data(tier), nbwu_progress.data(tier),
                             pmewu_progress.data(tier), gather_wu_progress.data(tier),
                             scatter_wu_progress.data(tier), all_reduce_wu_progress.data(tier));
}

//-------------------------------------------------------------------------------------------------
void MolecularMechanicsControls::primeWorkUnitCounters(const GpuDetails &gpu,
                                                       const AtomGraphSynthesis &poly_ag) {
  const int arch_major = gpu.getArchMajor();
  const int arch_minor = gpu.getArchMinor();
  int vwu_block_count = gpu.getSMPCount();
  int pmewu_block_count = gpu.getSMPCount();
  if (arch_major == 6 && arch_minor == 1) {
    vwu_block_count *= 2;
    pmewu_block_count *= 2;
  }
  int nbwu_block_count = gpu.getSMPCount();
  if (poly_ag.getUnitCellType() == UnitCellType::NONE) {

    // Use four non-bonded blocks per SM on Turing, five on all other NVIDIA architectures
    if (arch_major == 7 && arch_minor >= 5) { 
      nbwu_block_count *= 4;
    }
    else {
      nbwu_block_count *= 5;
    }
  }
  else {

    // Use two non-bonded blocks per SM on Turing, three on all other NVIDIA architectures
    if (arch_major == 7 && arch_minor >= 5) { 
      nbwu_block_count *= 2;
    }
    else {
      nbwu_block_count *= 3;
    }
  }
  int gtwu_block_count = gpu.getSMPCount();
  int scwu_block_count = gpu.getSMPCount();
  int rdwu_block_count = gpu.getSMPCount();
  for (int i = 0; i < twice_warp_size_int; i++) {
    vwu_progress.putHost(vwu_block_count, i);
    nbwu_progress.putHost(nbwu_block_count, i);
    pmewu_progress.putHost(pmewu_block_count, i);
    gather_wu_progress.putHost(gtwu_block_count, i);
    scatter_wu_progress.putHost(scwu_block_count, i);
    all_reduce_wu_progress.putHost(rdwu_block_count, i);
  }
#ifdef OMNI_USE_HPC
  upload();
#endif
}

//-------------------------------------------------------------------------------------------------
void MolecularMechanicsControls::incrementStep() {
  step_number += 1;
}
  
#ifdef OMNI_USE_HPC
//-------------------------------------------------------------------------------------------------
void MolecularMechanicsControls::upload() {
  progress_data.upload();
}

//-------------------------------------------------------------------------------------------------
void MolecularMechanicsControls::download() {
  progress_data.download();
}
#endif
  
} // namespace mm
} // namespace omni
