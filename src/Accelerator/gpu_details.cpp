// -*-c++-*-
#include <cstdio>
#include <cstdlib>
#include <vector>
#include "Reporting/error_format.h"
#include "gpu_details.h"

#ifndef OMNI_USE_HPC
omni::card::GpuDetails null_gpu;
#endif

namespace omni {
namespace card {

//-------------------------------------------------------------------------------------------------
GpuDetails::GpuDetails() :
    available{false}, supported{false}, arch_major{0}, arch_minor{0}, smp_count{0}, card_ram{0},
    max_threads_per_block{0}, max_threads_per_smp{0}, max_blocks_per_smp{0},
    max_shared_per_block{0}, max_shared_per_smp{0}, registers_per_block{0}, registers_per_smp{0},
    card_name{std::string("blank_gpu")}
{}

//-------------------------------------------------------------------------------------------------
bool GpuDetails::getAvailability() const {
  return available;
}

//-------------------------------------------------------------------------------------------------
bool GpuDetails::getGpuSupported() const {
  return supported;
}

//-------------------------------------------------------------------------------------------------
int GpuDetails::getArchMajor() const {
  return arch_major;
}

//-------------------------------------------------------------------------------------------------
int GpuDetails::getArchMinor() const {
  return arch_minor;
}

//-------------------------------------------------------------------------------------------------
int GpuDetails::getSMPCount() const {
  return smp_count;
}

//-------------------------------------------------------------------------------------------------
int GpuDetails::getCardRam() const {
  return card_ram;
}

//-------------------------------------------------------------------------------------------------
int GpuDetails::getMaxThreadsPerBlock() const {
  return max_threads_per_block;
}

//-------------------------------------------------------------------------------------------------
int GpuDetails::getMaxThreadsPerSMP() const {
  return max_threads_per_smp;
}

//-------------------------------------------------------------------------------------------------
int GpuDetails::getMaxBlocksPerSMP() const {
  return max_blocks_per_smp;
}

//-------------------------------------------------------------------------------------------------
int GpuDetails::getMaxSharedPerBlock() const {
  return max_shared_per_block;
}

//-------------------------------------------------------------------------------------------------
int GpuDetails::getMaxSharedPerSMP() const {
  return max_shared_per_smp;
}

//-------------------------------------------------------------------------------------------------
int GpuDetails::getRegistersPerBlock() const {
  return registers_per_block;
}

//-------------------------------------------------------------------------------------------------
int GpuDetails::getRegistersPerSMP() const {
  return registers_per_smp;
}

//-------------------------------------------------------------------------------------------------
std::string GpuDetails::getCardName() const {
  return card_name;
}

} // namespace card
} // namespace omni
