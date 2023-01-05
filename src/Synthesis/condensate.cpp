#include "copyright.h"
#include "Constants/behavior.h"
#include "Constants/fixed_precision.h"
#include "condensate.h"

namespace stormm {
namespace synthesis {

using card::HybridKind;
using constants::CartesianDimension;
using numerics::max_llint_accumulation;
using numerics::globalpos_scale_nonoverflow_bits;

//-------------------------------------------------------------------------------------------------
CondensateWriter::CondensateWriter(const CondensationLevel mode_in, const int natr_insr_in,
                                   const int nata_insr_in, float* xcrd_sp_in, float* ycrd_sp_in,
                                   float* zcrd_sp_in, double* xcrd_in, double* ycrd_in,
                                   double* zcrd_in, const int4* atr_insr_in,
                                   const int4* ata_insr_in) :
    mode{mode_in}, natr_insr{natr_insr_in}, nata_insr{nata_insr_in}, xcrd_sp{xcrd_sp_in},
    ycrd_sp{ycrd_sp_in}, zcrd_sp{zcrd_sp_in}, xcrd{xcrd_in}, ycrd{ycrd_in}, zcrd{zcrd_in},
    atr_insr{atr_insr_in}, ata_insr{ata_insr_in}
{}

//-------------------------------------------------------------------------------------------------
CondensateReader::CondensateReader(const CondensationLevel mode_in, const int natr_insr_in,
                                   const int nata_insr_in, const float* xcrd_sp_in,
                                   const float* ycrd_sp_in, const float* zcrd_sp_in,
                                   const double* xcrd_in, const double* ycrd_in,
                                   const double* zcrd_in, const int4* atr_insr_in,
                                   const int4* ata_insr_in) :
    mode{mode_in}, natr_insr{natr_insr_in}, nata_insr{nata_insr_in}, xcrd_sp{xcrd_sp_in},
    ycrd_sp{ycrd_sp_in}, zcrd_sp{zcrd_sp_in}, xcrd{xcrd_in}, ycrd{ycrd_in}, zcrd{zcrd_in},
    atr_insr{atr_insr_in}, ata_insr{ata_insr_in}
{}

//-------------------------------------------------------------------------------------------------
CondensateReader::CondensateReader(const CondensateWriter &cdw) :
    mode{cdw.mode}, natr_insr{cdw.natr_insr}, nata_insr{cdw.nata_insr}, xcrd_sp{cdw.xcrd_sp},
    ycrd_sp{cdw.ycrd_sp}, zcrd_sp{cdw.zcrd_sp}, xcrd{cdw.xcrd}, ycrd{cdw.ycrd}, zcrd{cdw.zcrd},
    atr_insr{cdw.atr_insr}, ata_insr{cdw.ata_insr}
{}

//-------------------------------------------------------------------------------------------------
Condensate::Condensate(const PhaseSpaceSynthesis *poly_ps_in, const CondensationLevel mode_in,
                       const GpuDetails &gpu) :
    mode{mode_in}, holds_own_data{false}, csptr_data_type{0}, atr_instruction_count{0},
    ata_instruction_count{0},
    x_coordinates_sp{HybridKind::POINTER, "cdns_xcrd_sp"},
    y_coordinates_sp{HybridKind::POINTER, "cdns_ycrd_sp"},
    z_coordinates_sp{HybridKind::POINTER, "cdns_zcrd_sp"},
    x_coordinates{HybridKind::POINTER, "cdns_xcrd"},
    y_coordinates{HybridKind::POINTER, "cdns_ycrd"},
    z_coordinates{HybridKind::POINTER, "cdns_zcrd"},
    atr_instructions{HybridKind::ARRAY, "cdns_atr_insr"},
    ata_instructions{HybridKind::ARRAY, "cdns_ata_insr"},
    pps_ptr{const_cast<PhaseSpaceSynthesis*>(poly_ps_in)},
    cs_ptr{nullptr},
    float_data{HybridKind::ARRAY, "cdns_floats"},
    double_data{HybridKind::ARRAY, "cdns_doubles"}
{
  rebuild(pps_ptr, mode, gpu);
}

//-------------------------------------------------------------------------------------------------
Condensate::Condensate(const Condensate &original) :
    mode{original.mode},
    holds_own_data{original.holds_own_data},
    csptr_data_type{original.csptr_data_type},
    atr_instruction_count{original.atr_instruction_count},
    ata_instruction_count{original.ata_instruction_count},
    x_coordinates_sp{original.x_coordinates_sp},
    y_coordinates_sp{original.y_coordinates_sp},
    z_coordinates_sp{original.z_coordinates_sp},
    x_coordinates{original.x_coordinates},
    y_coordinates{original.y_coordinates},
    z_coordinates{original.z_coordinates},
    atr_instructions{original.atr_instructions},
    ata_instructions{original.ata_instructions},
    pps_ptr{original.pps_ptr},
    cs_ptr{original.cs_ptr},
    float_data{original.float_data},
    double_data{original.double_data}
{
  repairPointers();
}

//-------------------------------------------------------------------------------------------------
Condensate::Condensate(Condensate &&original) :
    mode{original.mode},
    holds_own_data{original.holds_own_data},
    csptr_data_type{original.csptr_data_type},
    atr_instruction_count{original.atr_instruction_count},
    ata_instruction_count{original.ata_instruction_count},
    x_coordinates_sp{std::move(original.x_coordinates_sp)},
    y_coordinates_sp{std::move(original.y_coordinates_sp)},
    z_coordinates_sp{std::move(original.z_coordinates_sp)},
    x_coordinates{std::move(original.x_coordinates)},
    y_coordinates{std::move(original.y_coordinates)},
    z_coordinates{std::move(original.z_coordinates)},
    atr_instructions{std::move(original.atr_instructions)},
    ata_instructions{std::move(original.ata_instructions)},
    pps_ptr{original.pps_ptr},
    cs_ptr{original.cs_ptr},
    float_data{std::move(original.float_data)},
    double_data{std::move(original.double_data)}
{}

//-------------------------------------------------------------------------------------------------
Condensate& Condensate::operator=(const Condensate &other) {

  // Guard against self-assignment
  if (this == &other) {
    return *this;
  }

  // Typical copying
  mode = other.mode;
  holds_own_data = other.holds_own_data;
  csptr_data_type = other.csptr_data_type;
  atr_instruction_count = other.atr_instruction_count;
  ata_instruction_count = other.ata_instruction_count;
  x_coordinates_sp = other.x_coordinates_sp;
  y_coordinates_sp = other.y_coordinates_sp;
  z_coordinates_sp = other.z_coordinates_sp;
  x_coordinates = other.x_coordinates;
  y_coordinates = other.y_coordinates;
  z_coordinates = other.z_coordinates;
  atr_instructions = other.atr_instructions;
  ata_instructions = other.ata_instructions;
  pps_ptr = other.pps_ptr;
  cs_ptr = other.cs_ptr;
  float_data = other.float_data;
  double_data = other.double_data;

  // Repair pointers and return the result
  repairPointers();
  return *this;
}

//-------------------------------------------------------------------------------------------------
Condensate& Condensate::operator=(Condensate &&other) {

  // Guard against self-assignment
  if (this == &other) {
    return *this;
  }

  // Typical copying
  mode = other.mode;
  holds_own_data = other.holds_own_data;
  csptr_data_type = other.csptr_data_type;
  atr_instruction_count = other.atr_instruction_count;
  ata_instruction_count = other.ata_instruction_count;
  x_coordinates_sp = std::move(other.x_coordinates_sp);
  y_coordinates_sp = std::move(other.y_coordinates_sp);
  z_coordinates_sp = std::move(other.z_coordinates_sp);
  x_coordinates = std::move(other.x_coordinates);
  y_coordinates = std::move(other.y_coordinates);
  z_coordinates = std::move(other.z_coordinates);
  atr_instructions = std::move(other.atr_instructions);
  ata_instructions = std::move(other.ata_instructions);
  pps_ptr = std::move(other.pps_ptr);
  cs_ptr = std::move(other.cs_ptr);
  float_data = std::move(other.float_data);
  double_data = std::move(other.double_data);

  return *this;
}

//-------------------------------------------------------------------------------------------------
Condensate::Condensate(const PhaseSpaceSynthesis &poly_ps_in, const CondensationLevel mode_in,
                       const GpuDetails &gpu) :
    Condensate(poly_ps_in.getSelfPointer(), mode_in, gpu)
{}

//-------------------------------------------------------------------------------------------------
CondensationLevel Condensate::getMode() const {
  return mode;
}

//-------------------------------------------------------------------------------------------------
bool Condensate::ownsCoordinates() const {
  return holds_own_data;
}

//-------------------------------------------------------------------------------------------------
int Condensate::getATRInstructionCount() const {
  return atr_instruction_count;
}

//-------------------------------------------------------------------------------------------------
int Condensate::getATAInstructionCount() const {
  return ata_instruction_count;
}

//-------------------------------------------------------------------------------------------------
const PhaseSpaceSynthesis* Condensate::getSynthesisPointer() const {
  return pps_ptr;
}

//-------------------------------------------------------------------------------------------------
const CondensateReader Condensate::data(const HybridTargetLevel tier) const {
  return CondensateReader(mode, atr_instruction_count, ata_instruction_count,
                          x_coordinates_sp.data(tier), y_coordinates_sp.data(tier),
                          z_coordinates_sp.data(tier), x_coordinates.data(tier),
                          y_coordinates.data(tier), z_coordinates.data(tier),
                          atr_instructions.data(tier), ata_instructions.data(tier));
};  

//-------------------------------------------------------------------------------------------------
CondensateWriter Condensate::data(const HybridTargetLevel tier) {
  return CondensateWriter(mode, atr_instruction_count, ata_instruction_count,
                          x_coordinates_sp.data(tier), y_coordinates_sp.data(tier),
                          z_coordinates_sp.data(tier), x_coordinates.data(tier),
                          y_coordinates.data(tier), z_coordinates.data(tier),
                          atr_instructions.data(tier), ata_instructions.data(tier));
};

#ifdef STORMM_USE_HPC
//-------------------------------------------------------------------------------------------------
void Condensate::upload() {
  atr_instructions.upload();
  ata_instructions.upload();
  float_data.upload();
  double_data.upload();
}

//-------------------------------------------------------------------------------------------------
void Condensate::download() {
  atr_instructions.download();
  ata_instructions.download();
  float_data.download();
  double_data.download();
}
#endif

//-------------------------------------------------------------------------------------------------
void Condensate::rebuild(const PhaseSpaceSynthesis *poly_ps_in, const CondensationLevel mode_in,
                         const GpuDetails &gpu) {
  pps_ptr = const_cast<PhaseSpaceSynthesis*>(poly_ps_in);
  cs_ptr = nullptr;
  const PsSynthesisReader poly_psr = pps_ptr->data();
  const int last_sys = poly_psr.system_count - 1;
  const size_t padded_atoms = poly_psr.atom_starts[last_sys] +
                              roundUp(poly_psr.atom_counts[last_sys], warp_size_int);
  mode = mode_in;
  holds_own_data = true;
  switch (mode) {
  case CondensationLevel::DOUBLE:
    float_data.resize(0);
    float_data.shrinkToFit();
    x_coordinates_sp.setPointer(&float_data, 0, 0);
    y_coordinates_sp.setPointer(&float_data, 0, 0);
    z_coordinates_sp.setPointer(&float_data, 0, 0);
    double_data.resize(3LLU * padded_atoms);
    double_data.shrinkToFit();
    x_coordinates.setPointer(&double_data,                   0, padded_atoms);
    y_coordinates.setPointer(&double_data,        padded_atoms, padded_atoms);
    z_coordinates.setPointer(&double_data, 2LLU * padded_atoms, padded_atoms);
    break;
  case CondensationLevel::FLOAT:
    float_data.resize(3LLU * padded_atoms);
    float_data.shrinkToFit();
    x_coordinates_sp.setPointer(&float_data,                   0, padded_atoms);
    y_coordinates_sp.setPointer(&float_data,        padded_atoms, padded_atoms);
    z_coordinates_sp.setPointer(&float_data, 2LLU * padded_atoms, padded_atoms);
    double_data.resize(0);
    double_data.shrinkToFit();
    x_coordinates.setPointer(&double_data, 0, 0);
    y_coordinates.setPointer(&double_data, 0, 0);
    z_coordinates.setPointer(&double_data, 0, 0);
    break;
  }

  // Map the work units.  Condensates based on CoordinateSeries objects have a similar process,
  // but with only one topology to consider it becomes different enough that that builder has its
  // own code.
  if (gpu != null_gpu) {

    // Estimate the number of comparisons for the all-to-reference and all-to-all cases based on
    // the raw number of systems using each topology.
    llint np_atr = poly_psr.system_count;
    llint np_ata = 0LL;
    for (int i = 0; i < poly_psr.unique_topology_count; i++) {
      const llint nsys = poly_psr.common_ag_bounds[i + 1] - poly_psr.common_ag_bounds[i];
      np_ata += nsys * nsys;
    }
    const llint nsmp = gpu.getSMPCount();
    np_atr = roundUp(np_atr, nsmp) / nsmp;
    np_atr = roundUp(np_atr, static_cast<llint>(medium_block_size / warp_size_int));
    if (np_ata > 32768) {
      np_ata = 16;
    }
    else {
      np_ata = 8;
    }

    // The preferred sizes of each work unit are now determined.  Lay out as many work units as are
    // needed to cover all replicas of all systems.
    llint nwu_atr = 0LL;
    llint nwu_ata = 0LL;
    for (int i = 0; i < poly_psr.unique_topology_count; i++) {
      const llint nsys = poly_psr.common_ag_bounds[i + 1] - poly_psr.common_ag_bounds[i];
      nwu_atr += roundUp(nsys, np_atr) / np_atr;
      const llint ata_strides = roundUp(nsys, np_ata) / np_ata;
      nwu_ata += ata_strides * (ata_strides + 1LL) / 2LL;      
    }
    atr_instruction_count = nwu_atr;
    ata_instruction_count = nwu_ata;
    atr_instructions.resize(atr_instruction_count);
    ata_instructions.resize(ata_instruction_count);
    int4* atr_insr_ptr = atr_instructions.data();
    int4* ata_insr_ptr = ata_instructions.data();
    int atr_insr_idx = 0;
    int ata_insr_idx = 0;
    for (int i = 0; i < poly_psr.unique_topology_count; i++) {
      const int nsys = poly_psr.common_ag_bounds[i + 1] - poly_psr.common_ag_bounds[i];
      const int inp_atr = np_atr;
      const int inp_ata = np_ata;
      const int atr_strides = roundUp(nsys, inp_atr) / np_atr;
      const int ata_strides = roundUp(nsys, inp_ata) / np_ata;
      for (int j = 0; j < atr_strides; j++) {
        const int rdims = std::min(nsys - (j * inp_atr), inp_atr);
        const int4 rtmp = { j * inp_atr, 0, i, rdims };
        atr_insr_ptr[atr_insr_idx] = rtmp;
        atr_insr_idx++;
      }
      for (int j = 0; j < ata_strides; j++) {
        const int xdims = std::min(nsys - (j * inp_ata), inp_ata);
        for (int k = 0; k <= j; k++) {
          const int ydims = std::min(nsys - (k * inp_ata), inp_ata);
          const int4 mtmp = { j * inp_ata, k * inp_ata, i, ((ydims << 16) | xdims) };
          ata_insr_ptr[ata_insr_idx] = mtmp;
          ata_insr_idx++;
        }
      }
    }
  }
  else {

    // If there is no GPU available, all work will be done on the CPU.  Allocate one work unit per
    // system for either case.
    atr_instruction_count = 0;
    ata_instruction_count = 0;
    for (int i = 0; i < poly_psr.unique_topology_count; i++) {
      const int ni_sys = poly_psr.common_ag_bounds[i + 1] - poly_psr.common_ag_bounds[i];
      const int atr_strides = roundUp(ni_sys, 1024) / 1024;
      atr_instruction_count += atr_strides;
      ata_instruction_count += atr_strides * (atr_strides + 1) / 2;
    }
    atr_instructions.resize(atr_instruction_count);
    ata_instructions.resize(ata_instruction_count);
    int4* atr_insr_ptr = atr_instructions.data();
    int4* ata_insr_ptr = ata_instructions.data();
    int ata_insr_idx = 0;
    for (int i = 0; i < poly_psr.unique_topology_count; i++) {
      const int ni_sys = poly_psr.common_ag_bounds[i + 1] - poly_psr.common_ag_bounds[i];
      const int atr_strides = roundUp(ni_sys, 1024) / 1024;
      for (int j = 0; j < atr_strides; j++) {
        const int xdim = std::min(ni_sys - (j * 1024), 1024);
        const int4 rtmp = { j * 1024, 0, i, xdim };
        atr_insr_ptr[i] = rtmp;
        for (int k = 0; k <= j; k++) {
          const int ydim = (k == j) ? xdim : 1024;
          const int4 mtmp = { j * 1024, k * 1024, i, ((ydim << 16) | xdim) };
          ata_insr_ptr[ata_insr_idx] = mtmp;
          ata_insr_idx++;
        }
      }
    }
  }
  update();
}

//-------------------------------------------------------------------------------------------------
void Condensate::rebuild(const PhaseSpaceSynthesis &poly_ps_in, const CondensationLevel mode_in,
                         const GpuDetails &gpu) {
  rebuild(poly_ps_in.getSelfPointer(), mode_in, gpu);
}

//-------------------------------------------------------------------------------------------------
void Condensate::rebuild(const CoordinateSeries<double> *cs_in, const CondensationLevel mode_in,
                         const GpuDetails &gpu) {

  // See the templated rebuild() member function for explanation of the reinterpret_cast() call.
  cs_ptr = reinterpret_cast<CoordinateSeries<int>*>(const_cast<CoordinateSeries<double>*>(cs_in));
  pps_ptr = nullptr;
  csptr_data_type = double_type_index;
  const CoordinateSeriesReader<double> csr = cs_in->data();
  const	size_t padded_atoms = static_cast<size_t>(csr.nframe) *
                              roundUp(static_cast<size_t>(csr.natom), warp_size_zu);
  mode = mode_in;
  const CartesianDimension xdim = CartesianDimension::X;
  const CartesianDimension ydim = CartesianDimension::Y;
  const CartesianDimension zdim = CartesianDimension::Z;
  switch (mode) {
  case CondensationLevel::DOUBLE:
    holds_own_data = false;
    float_data.resize(0);
    float_data.shrinkToFit();
    x_coordinates_sp.setPointer(&float_data, 0, 0);
    y_coordinates_sp.setPointer(&float_data, 0, 0);
    z_coordinates_sp.setPointer(&float_data, 0, 0);
    double_data.resize(0);
    double_data.shrinkToFit();
    x_coordinates.setPointer(const_cast<Hybrid<double>*>(cs_in->getCoordinatePointer(xdim)),
                             0, padded_atoms);
    y_coordinates.setPointer(const_cast<Hybrid<double>*>(cs_in->getCoordinatePointer(ydim)),
                             0, padded_atoms);
    z_coordinates.setPointer(const_cast<Hybrid<double>*>(cs_in->getCoordinatePointer(zdim)),
                             0, padded_atoms);
    break;
  case CondensationLevel::FLOAT:
    holds_own_data = true;
    float_data.resize(3LLU * padded_atoms);
    float_data.shrinkToFit();
    x_coordinates_sp.setPointer(&float_data,                   0, padded_atoms);
    y_coordinates_sp.setPointer(&float_data,        padded_atoms, padded_atoms);
    z_coordinates_sp.setPointer(&float_data, 2LLU * padded_atoms, padded_atoms);
    double_data.resize(0);
    double_data.shrinkToFit();
    x_coordinates.setPointer(&double_data, 0, 0);
    y_coordinates.setPointer(&double_data, 0, 0);
    z_coordinates.setPointer(&double_data, 0, 0);
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void Condensate::rebuild(const CoordinateSeries<float> *cs_in, const CondensationLevel mode_in,
                         const GpuDetails &gpu) {

  // See the templated rebuild() member function for explanation of the reinterpret_cast() call.
  cs_ptr = reinterpret_cast<CoordinateSeries<int>*>(const_cast<CoordinateSeries<float>*>(cs_in));
  pps_ptr = nullptr;
  csptr_data_type = float_type_index;
  const CoordinateSeriesReader<float> csr = cs_in->data();
  const	size_t padded_atoms = static_cast<size_t>(csr.nframe) *
                              roundUp(static_cast<size_t>(csr.natom), warp_size_zu);
  mode = mode_in;
  const CartesianDimension xdim = CartesianDimension::X;
  const CartesianDimension ydim = CartesianDimension::Y;
  const CartesianDimension zdim = CartesianDimension::Z;
  switch (mode) {
  case CondensationLevel::DOUBLE:
    holds_own_data = true;
    float_data.resize(0);
    float_data.shrinkToFit();
    x_coordinates_sp.setPointer(&float_data, 0, 0);
    y_coordinates_sp.setPointer(&float_data, 0, 0);
    z_coordinates_sp.setPointer(&float_data, 0, 0);
    double_data.resize(3LLU * padded_atoms);
    double_data.shrinkToFit();
    x_coordinates.setPointer(&double_data,                   0, padded_atoms);
    y_coordinates.setPointer(&double_data,        padded_atoms, padded_atoms);
    z_coordinates.setPointer(&double_data, 2LLU * padded_atoms, padded_atoms);
    break;
  case CondensationLevel::FLOAT:
    holds_own_data = false;
    float_data.resize(0);
    float_data.shrinkToFit();
    x_coordinates_sp.setPointer(const_cast<Hybrid<float>*>(cs_in->getCoordinatePointer(xdim)),
                                0, padded_atoms);
    y_coordinates_sp.setPointer(const_cast<Hybrid<float>*>(cs_in->getCoordinatePointer(ydim)),
                                0, padded_atoms);
    z_coordinates_sp.setPointer(const_cast<Hybrid<float>*>(cs_in->getCoordinatePointer(zdim)),
                                0, padded_atoms);
    double_data.resize(0);
    double_data.shrinkToFit();
    x_coordinates.setPointer(&double_data, 0, 0);
    y_coordinates.setPointer(&double_data, 0, 0);
    z_coordinates.setPointer(&double_data, 0, 0);
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void Condensate::update(const HybridTargetLevel tier, const GpuDetails &gpu) {

  // Updating is only relevant if the object holds its own, separate copy of the coordinates.
  if (holds_own_data) {
    if (pps_ptr != nullptr) {
      switch (tier) {
      case HybridTargetLevel::HOST:
        {
          const PsSynthesisReader poly_psr = pps_ptr->data(tier);
          CondensateWriter cdw = this->data(tier);
          for (int i = 0; i < poly_psr.system_count; i++) {
            const int llim = poly_psr.atom_starts[i];
            const int hlim = llim + poly_psr.atom_counts[i];
            for (int j = llim; j < hlim; j++) {
              double xij = poly_psr.xcrd[j];
              double yij = poly_psr.ycrd[j];
              double zij = poly_psr.zcrd[j];
              if (poly_psr.gpos_bits >= globalpos_scale_nonoverflow_bits) {
                xij += static_cast<double>(poly_psr.xcrd_ovrf[j]) * max_llint_accumulation;
                yij += static_cast<double>(poly_psr.ycrd_ovrf[j]) * max_llint_accumulation;
                zij += static_cast<double>(poly_psr.zcrd_ovrf[j]) * max_llint_accumulation;
              }
              xij *= poly_psr.inv_gpos_scale;
              yij *= poly_psr.inv_gpos_scale;
              zij *= poly_psr.inv_gpos_scale;
              switch (mode) {
              case CondensationLevel::DOUBLE:
                cdw.xcrd[j] = xij;
                cdw.ycrd[j] = yij;
                cdw.zcrd[j] = zij;
                break;
              case CondensationLevel::FLOAT:
                cdw.xcrd_sp[j] = xij;
                cdw.ycrd_sp[j] = yij;
                cdw.zcrd_sp[j] = zij;
                break;
              }
            }
          }
        }
        break;
#ifdef STORMM_USE_HPC
      case HybridTargetLevel::DEVICE:
        launchCondensateUpdate(gpu);
        break;
#endif
      }
    }
    else if (cs_ptr != nullptr) {

      // If not based on a PhaseSpaceSynthesis, the Condensate must be based on a CoordinateSeries.
      // Assess the type of the CoordinateSeries, making a switch-like apparatus to fire off the
      // templated overload of this function with the correct CoordinateSeries<T> pointer.
      if (csptr_data_type == double_type_index) {
        update(reinterpret_cast<CoordinateSeries<double>*>(cs_ptr), tier, gpu);
      }
      else if (csptr_data_type == float_type_index) {
        update(reinterpret_cast<CoordinateSeries<float>*>(cs_ptr), tier, gpu);
      }
      else if (csptr_data_type == llint_type_index) {
        update(reinterpret_cast<CoordinateSeries<llint>*>(cs_ptr), tier, gpu);
      }
      else if (csptr_data_type == int_type_index) {
        update(reinterpret_cast<CoordinateSeries<int>*>(cs_ptr), tier, gpu);
      }
      else if (csptr_data_type == short_type_index) {
        update(reinterpret_cast<CoordinateSeries<int>*>(cs_ptr), tier, gpu);
      }
      else if (csptr_data_type == char_type_index) {
        update(reinterpret_cast<CoordinateSeries<char>*>(cs_ptr), tier, gpu);
      }
      else {
        rtErr("A CoordinateSeries must be typed as double, float, or some signed integer in "
              "order to submit for analysis.  Check the data type, or call the update() function "
              "by directly supplying a pointer to the original CoordinateSeries.", "Condensate",
              "update");
      }
    }
    else {
      rtErr("There is no current coordinate synthesis or series to base the object upon.",
            "Condensate", "update");
    }
  }
}

//-------------------------------------------------------------------------------------------------
void Condensate::repairPointers() {

  // Repairs only occur if the object holds its own data.  Otherwise, a Condensate that points at
  // the contents of a CoordinateSeries can be copied with no need for pointer manipulation.
  if (holds_own_data) {
    x_coordinates_sp.swapTarget(&float_data);
    y_coordinates_sp.swapTarget(&float_data);
    z_coordinates_sp.swapTarget(&float_data);
    x_coordinates.swapTarget(&double_data);
    y_coordinates.swapTarget(&double_data);
    z_coordinates.swapTarget(&double_data);
  }
}
  
} // namespace synthesis
} // namespace stormm
