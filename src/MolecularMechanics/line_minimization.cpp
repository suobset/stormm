#include "Constants/scaling.h"
#include "Math/rounding.h"
#include "Reporting/error_format.h"
#include "line_minimization.h"

namespace omni {
namespace mm {

using card::HybridKind;
using math::roundUp;

//------------------------------------------------------------------------------------------------
LinMinWriter::LinMinWriter(const int nsys_in, double* move_a_in, double* move_b_in,
                           double* move_c_in, double* nrg_a_in, double* nrg_b_in,
                           double* nrg_c_in, double* nrg_d_in) :
    nsys{nsys_in}, move_a{move_a_in}, move_b{move_b_in}, move_c{move_c_in}, nrg_a{nrg_a_in},
    nrg_b{nrg_b_in}, nrg_c{nrg_c_in}, nrg_d{nrg_d_in}
{}

//------------------------------------------------------------------------------------------------
LinMinReader::LinMinReader(const int nsys_in, const double* move_a_in, const double* move_b_in,
                           const double* move_c_in, const double* nrg_a_in,
                           const double* nrg_b_in, const double* nrg_c_in,
                           const double* nrg_d_in) :
    nsys{nsys_in}, move_a{move_a_in}, move_b{move_b_in}, move_c{move_c_in}, nrg_a{nrg_a_in},
    nrg_b{nrg_b_in}, nrg_c{nrg_c_in}, nrg_d{nrg_d_in}
{}

//------------------------------------------------------------------------------------------------
LineMinimization::LineMinimization(const int system_count_in) :
  system_count{system_count_in},
  move_a{HybridKind::POINTER, "linmin_move_a"},
  move_b{HybridKind::POINTER, "linmin_move_b"},
  move_c{HybridKind::POINTER, "linmin_move_c"},
  energy_a{HybridKind::POINTER, "linmin_nrg_a"},
  energy_b{HybridKind::POINTER, "linmin_nrg_b"},
  energy_c{HybridKind::POINTER, "linmin_nrg_c"},
  energy_d{HybridKind::POINTER, "linmin_nrg_d"},
  storage{static_cast<size_t>(roundUp(system_count, warp_size_int) * 7), "linmin_storage"}
{
  const size_t padded_nsys = roundUp(system_count, warp_size_int);
  move_a.setPointer(&storage,                  0, system_count);
  move_b.setPointer(&storage,        padded_nsys, system_count);
  move_c.setPointer(&storage, 2LLU * padded_nsys, system_count);
  energy_a.setPointer(&storage, 3LLU * padded_nsys, system_count);
  energy_b.setPointer(&storage, 4LLU * padded_nsys, system_count);
  energy_c.setPointer(&storage, 5LLU * padded_nsys, system_count);
  energy_d.setPointer(&storage, 6LLU * padded_nsys, system_count);
}

//------------------------------------------------------------------------------------------------
LineMinimization::LineMinimization(const LineMinimization &original) :
  system_count{original.system_count},
  move_a{original.move_a},
  move_b{original.move_b},
  move_c{original.move_c},
  energy_a{original.energy_a},
  energy_b{original.energy_b},
  energy_c{original.energy_c},
  energy_d{original.energy_d},
  storage{original.storage}
{
  // Repair pointers
  move_a.swapTarget(&storage);
  move_b.swapTarget(&storage);
  move_c.swapTarget(&storage);
  energy_a.swapTarget(&storage);
  energy_b.swapTarget(&storage);
  energy_c.swapTarget(&storage);
  energy_d.swapTarget(&storage);
}

//------------------------------------------------------------------------------------------------
LineMinimization& LineMinimization::operator=(const LineMinimization &other) {
  if (this == &other) {
    return *this;
  }
  system_count = other.system_count;
  move_a = std::move(other.move_a);
  move_b = std::move(other.move_b);
  move_c = std::move(other.move_c);
  energy_a = std::move(other.energy_a);
  energy_b = std::move(other.energy_b);
  energy_c = std::move(other.energy_c);
  energy_d = std::move(other.energy_d);
  return *this;
}

//------------------------------------------------------------------------------------------------
int LineMinimization::getSystemCount() const {
  return system_count;
}

//------------------------------------------------------------------------------------------------
std::vector<double> LineMinimization::getMoveLength(const int move_index,
                                                    const HybridTargetLevel tier) const {
  bool problem = false;
  switch (tier) {
  case HybridTargetLevel::HOST:
    switch (move_index) {
    case 0:
      return move_a.readHost();
    case 1:
      return move_b.readHost();
    case 2:
      return move_c.readHost();
    default:
      problem = true;
      break;
    }
    break;
#ifdef OMNI_USE_HPC
  case HybridTargetLevel::DEVICE:
    switch (move_index) {
    case 0:
      return move_a.readDevice();
    case 1:
      return move_b.readDevice();
    case 2:
      return move_c.readDevice();
    default:
      problem = true;
      break;
    }
    break;
#endif
  }
  if (problem) {    
    rtErr("Up to three moves are taken along any computed gradient during a line minimization, "
          "prior to solving a cubic polynomial to determine the optimal move.  A move index of " +
          std::to_string(move_index) + " is invalid.", "LineMinimization", "getMoveLength");
  }
  __builtin_unreachable();
}

//------------------------------------------------------------------------------------------------
double LineMinimization::getMoveLength(const int move_index, const int system_index,
                                       const HybridTargetLevel tier) const {
  bool problem = false;
  switch (tier) {
  case HybridTargetLevel::HOST:
    switch (move_index) {
    case 0:
      return move_a.readHost(system_index);
    case 1:
      return move_b.readHost(system_index);
    case 2:
      return move_c.readHost(system_index);
    default:
      problem = true;
      break;
    }
    break;
#ifdef OMNI_USE_HPC
  case HybridTargetLevel::DEVICE:
    switch (move_index) {
    case 0:
      return move_a.readDevice(system_index);
    case 1:
      return move_b.readDevice(system_index);
    case 2:
      return move_c.readDevice(system_index);
    default:
      problem = true;
      break;
    }
    break;
#endif
  }
  if (problem) {    
    rtErr("Up to three moves are taken along any computed gradient during a line minimization, "
          "prior to solving a cubic polynomial to determine the optimal move.  A move index of " +
          std::to_string(move_index) + " is invalid.", "LineMinimization", "getMoveLength");
  }
  __builtin_unreachable();
}

//------------------------------------------------------------------------------------------------
std::vector<double> LineMinimization::getEnergy(const int move_index,
                                                const HybridTargetLevel tier) const {
  bool problem = false;
  switch (tier) {
  case HybridTargetLevel::HOST:
    switch (move_index) {
    case 0:
      return energy_a.readHost();
    case 1:
      return energy_b.readHost();
    case 2:
      return energy_c.readHost();
    case 3:
      return energy_d.readHost();
    default:
      problem = true;
      break;
    }
    break;
#ifdef OMNI_USE_HPC
  case HybridTargetLevel::DEVICE:
    switch (move_index) {
    case 0:
      return energy_a.readDevice();
    case 1:
      return energy_b.readDevice();
    case 2:
      return energy_c.readDevice();
    case 3:
      return energy_d.readDevice();
    default:
      problem = true;
      break;
    }
    break;
#endif
  }
  if (problem) {    
    rtErr("Up to three moves are taken along any computed gradient during a line minimization, "
          "prior to solving a cubic polynomial to determine the optimal move.  A move index of " +
          std::to_string(move_index) + " is invalid.", "LineMinimization", "getEnergy");
  }
  __builtin_unreachable();
}

//------------------------------------------------------------------------------------------------
double LineMinimization::getEnergy(const int move_index,
                                   const int system_index, const HybridTargetLevel tier) const {
  bool problem = false;
  switch (tier) {
  case HybridTargetLevel::HOST:
    switch (move_index) {
    case 0:
      return energy_a.readHost(system_index);
    case 1:
      return energy_b.readHost(system_index);
    case 2:
      return energy_c.readHost(system_index);
    case 3:
      return energy_d.readHost(system_index);
    default:
      problem = true;
      break;
    }
    break;
#ifdef OMNI_USE_HPC
  case HybridTargetLevel::DEVICE:
    switch (move_index) {
    case 0:
      return energy_a.readDevice(system_index);
    case 1:
      return energy_b.readDevice(system_index);
    case 2:
      return energy_c.readDevice(system_index);
    case 3:
      return energy_d.readDevice(system_index);
    default:
      problem = true;
      break;
    }
    break;
#endif
  }
  if (problem) {
    rtErr("Up to three moves are taken along any computed gradient during a line minimization, "
          "prior to solving a cubic polynomial to determine the optimal move.  A move index of " +
          std::to_string(move_index) + " is invalid.", "LineMinimization", "getEnergy");
  }
  __builtin_unreachable();
}

//------------------------------------------------------------------------------------------------
LinMinReader LineMinimization::data(const HybridTargetLevel tier) const {
  return LinMinReader(system_count, move_a.data(tier), move_b.data(tier), move_c.data(tier),
                      energy_a.data(tier), energy_b.data(tier), energy_c.data(tier),
                      energy_d.data(tier));
}

//------------------------------------------------------------------------------------------------
LinMinWriter LineMinimization::data(const HybridTargetLevel tier) {
  return LinMinWriter(system_count, move_a.data(tier), move_b.data(tier), move_c.data(tier),
                      energy_a.data(tier), energy_b.data(tier), energy_c.data(tier),
                      energy_d.data(tier));
}

} // namespace mm
} // namespace omni
