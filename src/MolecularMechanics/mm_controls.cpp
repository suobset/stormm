#include "Constants/fixed_precision.h"
#include "Constants/scaling.h"
#include "Potential/energy_enumerators.h"
#include "Synthesis/synthesis_enumerators.h"
#include "mm_controls.h"

namespace omni {
namespace mm {

using card::HybridKind;
using energy::EvaluateEnergy;
using energy::EvaluateForce;
using numerics::ForceAccumulationMethod;
using synthesis::VwuGoal;

//-------------------------------------------------------------------------------------------------
MolecularMechanicsControls::MolecularMechanicsControls(const double time_step_in,
                                                       const double rattle_tol_in,
                                                       const double initial_step_in) :
  step_number{0}, sd_cycles{0}, max_cycles{0}, time_step{time_step_in}, rattle_tol{rattle_tol_in},
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
  max_cycles = user_input.getStepCount();
  time_step = user_input.getTimeStep();  
  rattle_tol = user_input.getRattleTolerance();
}

//-------------------------------------------------------------------------------------------------
MolecularMechanicsControls::MolecularMechanicsControls(const MinimizeControls &user_input) :
    MolecularMechanicsControls()
{
  sd_cycles = user_input.getSteepestDescentCycles();
  max_cycles = user_input.getTotalCycles();
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
void MolecularMechanicsControls::primeWorkUnitCounters(const KernelManager &launcher,
                                                       const PrecisionModel prec,
                                                       const AtomGraphSynthesis &poly_ag) {
  const GpuDetails wgpu = launcher.getGpu();
  const int arch_major = wgpu.getArchMajor();
  const int arch_minor = wgpu.getArchMinor();
  const int smp_count = wgpu.getSMPCount();

  // The numbers of blocks that will be launched in each grid are critical for priming the work
  // unit progress counters.  As such, they (the grid launch size) should be consistent across
  // different variants of each kernel for a given precision level, even if the exact numbers of
  // threads per block have to vary based on what each block of that kernel variant is required to
  // do.  Here, it should suffice to query the launch parameters of just one of the blocks.
  const int2 vwu_lp = launcher.getValenceKernelDims(prec, EvaluateForce::YES, EvaluateEnergy::YES,
                                                    ForceAccumulationMethod::SPLIT,
                                                    VwuGoal::ACCUMULATE);
  const int2 nbwu_lp = launcher.getNonbondedKernelDims(prec, poly_ag.getNonbondedWorkType(),
                                                       EvaluateForce::YES, EvaluateEnergy::YES,
                                                       ForceAccumulationMethod::SPLIT);
  int vwu_block_count = vwu_lp.x;
  int nbwu_block_count = nbwu_lp.x;
  int pmewu_block_count = smp_count;
  int gtwu_block_count = smp_count;
  int scwu_block_count = smp_count;
  int rdwu_block_count = smp_count;
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
