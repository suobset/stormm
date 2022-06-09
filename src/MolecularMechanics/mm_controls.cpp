#include "mm_controls.h"

namespace omni {
namespace mm {

using card::HybridKind;

//-------------------------------------------------------------------------------------------------
MolecularMechanicsControls::MolecularMechanicsControls(const double time_step_in,
                                                       const double rattle_tol_in) :
  step_number{0}, time_step{time_step_in}, rattle_tol{rattle_tol_in},
  vwu_progress{HybridKind::POINTER, "mm_vwu_counters"},
  nbwu_progress{HybridKind::POINTER, "mm_nbwu_counters"},
  pmewu_progress{HybridKind::POINTER, "mm_pmewu_counters"},
  progress_data{warp_size_int * 6, "work_unit_prog_data"}
{
  vwu_progress.setPointer(&progress_data,   0,                 2 * warp_size_int);
  nbwu_progress.setPointer(&progress_data,  2 * warp_size_int, 2 * warp_size_int);
  pmewu_progress.setPointer(&progress_data, 4 * warp_size_int, 2 * warp_size_int);
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
  progress_data{original.progress_data}
{
  vwu_progress.swapTarget(&progress_data);
  nbwu_progress.swapTarget(&progress_data);
  pmewu_progress.swapTarget(&progress_data);
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
  progress_data = other.progress_data;

  // Repair pointers and return the result
  vwu_progress.swapTarget(&progress_data);
  nbwu_progress.swapTarget(&progress_data);
  pmewu_progress.swapTarget(&progress_data);
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
  progress_data = std::move(other.progress_data);
  return *this;
}

//-------------------------------------------------------------------------------------------------
MMControlKit<double> MolecularMechanicsControls::dpData(const HybridTargetLevel tier) {
  return MMControlKit<double>(step_number, time_step, rattle_tol, vwu_progress.data(tier),
                              nbwu_progress.data(tier), pmewu_progress.data(tier));
}

//-------------------------------------------------------------------------------------------------
MMControlKit<float> MolecularMechanicsControls::spData(const HybridTargetLevel tier) {
  return MMControlKit<float>(step_number, time_step, rattle_tol, vwu_progress.data(tier),
                             nbwu_progress.data(tier), pmewu_progress.data(tier));
}


  
} // namespace mm
} // namespace omni
