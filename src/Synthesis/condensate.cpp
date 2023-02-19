#include "copyright.h"
#include "Constants/fixed_precision.h"
#include "condensate.h"

namespace stormm {
namespace synthesis {

using card::HybridKind;
using numerics::max_llint_accumulation;
using numerics::globalpos_scale_nonoverflow_bits;

//-------------------------------------------------------------------------------------------------
CondensateWriter::CondensateWriter(const PrecisionModel mode_in, const CondensateSource basis_in,
                                   const int system_count_in, const int natr_insr_src_in,
                                   const int natr_insr_top_in, const int natr_insr_lbl_in,
                                   const int nata_insr_src_in, const int nata_insr_top_in,
                                   const int nata_insr_lbl_in, const size_t* atom_starts_in,
                                   const int* atom_counts_in, float* xcrd_sp_in, float* ycrd_sp_in,
                                   float* zcrd_sp_in, double* xcrd_in, double* ycrd_in,
                                   double* zcrd_in, double* umat_in, double* invu_in,
                                   const int4* atr_insr_src_in, const int4* atr_insr_top_in,
                                   const int4* atr_insr_lbl_in, const int* atr_group_src_in,
                                   const int* atr_group_top_in, const int* atr_group_lbl_in,
                                   const int4* ata_insr_src_in, const int4* ata_insr_top_in,
                                   const int4* ata_insr_lbl_in, const int* ata_group_src_in,
                                   const int* ata_group_top_in, const int* ata_group_lbl_in) :
    mode{mode_in}, basis{basis_in}, system_count{system_count_in}, natr_insr_src{natr_insr_src_in},
    natr_insr_top{natr_insr_top_in}, natr_insr_lbl{natr_insr_lbl_in},
    nata_insr_src{nata_insr_src_in}, nata_insr_top{nata_insr_top_in},
    nata_insr_lbl{nata_insr_lbl_in}, atom_starts{atom_starts_in}, atom_counts{atom_counts_in},
    xcrd_sp{xcrd_sp_in}, ycrd_sp{ycrd_sp_in}, zcrd_sp{zcrd_sp_in}, xcrd{xcrd_in}, ycrd{ycrd_in},
    zcrd{zcrd_in}, umat{umat_in}, invu{invu_in}, atr_insr_src{atr_insr_src_in},
    atr_insr_top{atr_insr_top_in}, atr_insr_lbl{atr_insr_lbl_in}, atr_group_src{atr_group_src_in},
    atr_group_top{atr_group_top_in}, atr_group_lbl{atr_group_lbl_in},
    ata_insr_src{ata_insr_src_in}, ata_insr_top{ata_insr_top_in}, ata_insr_lbl{ata_insr_lbl_in},
    ata_group_src{ata_group_src_in}, ata_group_top{ata_group_top_in},
    ata_group_lbl{ata_group_lbl_in}
{}

//-------------------------------------------------------------------------------------------------
CondensateReader::CondensateReader(const PrecisionModel mode_in, const CondensateSource basis_in,
                                   const int system_count_in, const int natr_insr_src_in,
                                   const int natr_insr_top_in, const int natr_insr_lbl_in,
                                   const int nata_insr_src_in, const int nata_insr_top_in,
                                   const int nata_insr_lbl_in, const size_t* atom_starts_in,
                                   const int* atom_counts_in, const float* xcrd_sp_in,
                                   const float* ycrd_sp_in, const float* zcrd_sp_in,
                                   const double* xcrd_in, const double* ycrd_in,
                                   const double* zcrd_in, const double* umat_in,
                                   const double* invu_in, const int4* atr_insr_src_in,
                                   const int4* atr_insr_top_in, const int4* atr_insr_lbl_in,
                                   const int* atr_group_src_in, const int* atr_group_top_in,
                                   const int* atr_group_lbl_in, const int4* ata_insr_src_in,
                                   const int4* ata_insr_top_in, const int4* ata_insr_lbl_in,
                                   const int* ata_group_src_in, const int* ata_group_top_in,
                                   const int* ata_group_lbl_in) :
    mode{mode_in}, basis{basis_in}, system_count{system_count_in}, natr_insr_src{natr_insr_src_in},
    natr_insr_top{natr_insr_top_in}, natr_insr_lbl{natr_insr_lbl_in},
    nata_insr_src{nata_insr_src_in}, nata_insr_top{nata_insr_top_in},
    nata_insr_lbl{nata_insr_lbl_in}, atom_starts{atom_starts_in}, atom_counts{atom_counts_in},
    xcrd_sp{xcrd_sp_in}, ycrd_sp{ycrd_sp_in}, zcrd_sp{zcrd_sp_in}, xcrd{xcrd_in}, ycrd{ycrd_in},
    zcrd{zcrd_in}, umat{umat_in}, invu{invu_in}, atr_insr_src{atr_insr_src_in},
    atr_insr_top{atr_insr_top_in}, atr_insr_lbl{atr_insr_lbl_in}, atr_group_src{atr_group_src_in},
    atr_group_top{atr_group_top_in}, atr_group_lbl{atr_group_lbl_in},
    ata_insr_src{ata_insr_src_in}, ata_insr_top{ata_insr_top_in}, ata_insr_lbl{ata_insr_lbl_in},
    ata_group_src{ata_group_src_in}, ata_group_top{ata_group_top_in},
    ata_group_lbl{ata_group_lbl_in}
{}

//-------------------------------------------------------------------------------------------------
CondensateReader::CondensateReader(const CondensateWriter &cdw) :
    mode{cdw.mode}, basis{cdw.basis}, system_count{cdw.system_count},
    natr_insr_src{cdw.natr_insr_src}, natr_insr_top{cdw.natr_insr_top},
    natr_insr_lbl{cdw.natr_insr_lbl}, nata_insr_src{cdw.nata_insr_src},
    nata_insr_top{cdw.nata_insr_top}, nata_insr_lbl{cdw.nata_insr_lbl},
    atom_starts{cdw.atom_starts}, atom_counts{cdw.atom_counts}, xcrd_sp{cdw.xcrd_sp},
    ycrd_sp{cdw.ycrd_sp}, zcrd_sp{cdw.zcrd_sp}, xcrd{cdw.xcrd}, ycrd{cdw.ycrd}, zcrd{cdw.zcrd},
    umat{cdw.umat}, invu{cdw.invu}, atr_insr_src{cdw.atr_insr_src}, atr_insr_top{cdw.atr_insr_top},
    atr_insr_lbl{cdw.atr_insr_lbl}, atr_group_src{cdw.atr_group_src},
    atr_group_top{cdw.atr_group_top}, atr_group_lbl{cdw.atr_group_lbl},
    ata_insr_src{cdw.ata_insr_src}, ata_insr_top{cdw.ata_insr_top}, ata_insr_lbl{cdw.ata_insr_lbl},
    ata_group_src{cdw.ata_group_src}, ata_group_top{cdw.ata_group_top},
    ata_group_lbl{cdw.ata_group_lbl}
{}

//-------------------------------------------------------------------------------------------------
Condensate::Condensate(const PhaseSpaceSynthesis *poly_ps_in, const PrecisionModel mode_in,
                       const GpuDetails &gpu) :
    mode{mode_in}, basis{CondensateSource::SYNTHESIS},
    system_count{0},
    holds_own_data{false}, csptr_data_type{0},
    atr_instruction_count_src{0}, atr_instruction_count_top{0}, atr_instruction_count_lbl{0},
    ata_instruction_count_src{0}, ata_instruction_count_top{0}, ata_instruction_count_lbl{0},
    atom_starts{HybridKind::ARRAY, "cdns_atom_starts"},
    atom_counts{HybridKind::ARRAY, "cdns_atom_counts"},
    x_coordinates_sp{HybridKind::POINTER, "cdns_xcrd_sp"},
    y_coordinates_sp{HybridKind::POINTER, "cdns_ycrd_sp"},
    z_coordinates_sp{HybridKind::POINTER, "cdns_zcrd_sp"},
    x_coordinates{HybridKind::POINTER, "cdns_xcrd"},
    y_coordinates{HybridKind::POINTER, "cdns_ycrd"},
    z_coordinates{HybridKind::POINTER, "cdns_zcrd"},
    umat{HybridKind::POINTER, "cdns_umat"},
    invu{HybridKind::POINTER, "cdns_invu"},
    atr_instructions_src{HybridKind::ARRAY, "cdns_atr_insr_src"},
    atr_instructions_top{HybridKind::ARRAY, "cdns_atr_insr_top"},
    atr_instructions_lbl{HybridKind::ARRAY, "cdns_atr_insr_lbl"},
    atr_instruction_groups_src{HybridKind::ARRAY, "cdns_atr_group_src"},
    atr_instruction_groups_top{HybridKind::ARRAY, "cdns_atr_group_top"},
    atr_instruction_groups_lbl{HybridKind::ARRAY, "cdns_atr_group_lbl"},
    ata_instructions_src{HybridKind::ARRAY, "cdns_ata_insr_src"},
    ata_instructions_top{HybridKind::ARRAY, "cdns_ata_insr_top"},
    ata_instructions_lbl{HybridKind::ARRAY, "cdns_ata_insr_lbl"},
    ata_instruction_groups_src{HybridKind::ARRAY, "cdns_ata_group_src"},
    ata_instruction_groups_top{HybridKind::ARRAY, "cdns_ata_group_top"},
    ata_instruction_groups_lbl{HybridKind::ARRAY, "cdns_ata_group_lbl"},
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
    basis{original.basis},
    system_count{original.system_count},
    holds_own_data{original.holds_own_data},
    csptr_data_type{original.csptr_data_type},
    atr_instruction_count_src{original.atr_instruction_count_src},
    atr_instruction_count_top{original.atr_instruction_count_top},
    atr_instruction_count_lbl{original.atr_instruction_count_lbl},
    ata_instruction_count_src{original.ata_instruction_count_src},
    ata_instruction_count_top{original.ata_instruction_count_top},
    ata_instruction_count_lbl{original.ata_instruction_count_lbl},
    atom_starts{original.atom_starts},
    atom_counts{original.atom_counts},
    x_coordinates_sp{original.x_coordinates_sp},
    y_coordinates_sp{original.y_coordinates_sp},
    z_coordinates_sp{original.z_coordinates_sp},
    x_coordinates{original.x_coordinates},
    y_coordinates{original.y_coordinates},
    z_coordinates{original.z_coordinates},
    umat{original.umat},
    invu{original.invu},
    atr_instructions_src{original.atr_instructions_src},
    atr_instructions_top{original.atr_instructions_top},
    atr_instructions_lbl{original.atr_instructions_lbl},
    atr_instruction_groups_src{original.atr_instruction_groups_src},
    atr_instruction_groups_top{original.atr_instruction_groups_top},
    atr_instruction_groups_lbl{original.atr_instruction_groups_lbl},
    ata_instructions_src{original.ata_instructions_src},
    ata_instructions_top{original.ata_instructions_top},
    ata_instructions_lbl{original.ata_instructions_lbl},
    ata_instruction_groups_src{original.ata_instruction_groups_src},
    ata_instruction_groups_top{original.ata_instruction_groups_top},
    ata_instruction_groups_lbl{original.ata_instruction_groups_lbl},
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
    basis{original.basis},
    system_count{original.system_count},
    holds_own_data{original.holds_own_data},
    csptr_data_type{original.csptr_data_type},
    atr_instruction_count_src{original.atr_instruction_count_src},
    atr_instruction_count_top{original.atr_instruction_count_top},
    atr_instruction_count_lbl{original.atr_instruction_count_lbl},
    ata_instruction_count_src{original.ata_instruction_count_src},
    ata_instruction_count_top{original.ata_instruction_count_top},
    ata_instruction_count_lbl{original.ata_instruction_count_lbl},
    atom_starts{std::move(original.atom_starts)},
    atom_counts{std::move(original.atom_counts)},
    x_coordinates_sp{std::move(original.x_coordinates_sp)},
    y_coordinates_sp{std::move(original.y_coordinates_sp)},
    z_coordinates_sp{std::move(original.z_coordinates_sp)},
    x_coordinates{std::move(original.x_coordinates)},
    y_coordinates{std::move(original.y_coordinates)},
    z_coordinates{std::move(original.z_coordinates)},
    umat{std::move(original.umat)},
    invu{std::move(original.invu)},
    atr_instructions_src{std::move(original.atr_instructions_src)},
    atr_instructions_top{std::move(original.atr_instructions_top)},
    atr_instructions_lbl{std::move(original.atr_instructions_lbl)},
    atr_instruction_groups_src{std::move(original.atr_instruction_groups_src)},
    atr_instruction_groups_top{std::move(original.atr_instruction_groups_top)},
    atr_instruction_groups_lbl{std::move(original.atr_instruction_groups_lbl)},
    ata_instructions_src{std::move(original.ata_instructions_src)},
    ata_instructions_top{std::move(original.ata_instructions_top)},
    ata_instructions_lbl{std::move(original.ata_instructions_lbl)},
    ata_instruction_groups_src{std::move(original.ata_instruction_groups_src)},
    ata_instruction_groups_top{std::move(original.ata_instruction_groups_top)},
    ata_instruction_groups_lbl{std::move(original.ata_instruction_groups_lbl)},
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
  basis = other.basis;
  system_count = other.system_count;
  holds_own_data = other.holds_own_data;
  csptr_data_type = other.csptr_data_type;
  atr_instruction_count_src = other.atr_instruction_count_src;
  atr_instruction_count_top = other.atr_instruction_count_top;
  atr_instruction_count_lbl = other.atr_instruction_count_lbl;
  ata_instruction_count_src = other.ata_instruction_count_src;
  ata_instruction_count_top = other.ata_instruction_count_top;
  ata_instruction_count_lbl = other.ata_instruction_count_lbl;
  atom_starts = other.atom_starts;
  atom_counts = other.atom_counts;
  x_coordinates_sp = other.x_coordinates_sp;
  y_coordinates_sp = other.y_coordinates_sp;
  z_coordinates_sp = other.z_coordinates_sp;
  x_coordinates = other.x_coordinates;
  y_coordinates = other.y_coordinates;
  z_coordinates = other.z_coordinates;
  umat = other.umat;
  invu = other.invu;
  atr_instructions_src = other.atr_instructions_src;
  atr_instructions_top = other.atr_instructions_top;
  atr_instructions_lbl = other.atr_instructions_lbl;
  atr_instruction_groups_src = other.atr_instruction_groups_src;
  atr_instruction_groups_top = other.atr_instruction_groups_top;
  atr_instruction_groups_lbl = other.atr_instruction_groups_lbl;
  ata_instructions_src = other.ata_instructions_src;
  ata_instructions_top = other.ata_instructions_top;
  ata_instructions_lbl = other.ata_instructions_lbl;
  ata_instruction_groups_src = other.ata_instruction_groups_src;
  ata_instruction_groups_top = other.ata_instruction_groups_top;
  ata_instruction_groups_lbl = other.ata_instruction_groups_lbl;
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
  basis = other.basis;
  system_count = other.system_count;
  holds_own_data = other.holds_own_data;
  csptr_data_type = other.csptr_data_type;
  atr_instruction_count_src = other.atr_instruction_count_src;
  atr_instruction_count_top = other.atr_instruction_count_top;
  atr_instruction_count_lbl = other.atr_instruction_count_lbl;
  ata_instruction_count_src = other.ata_instruction_count_src;
  ata_instruction_count_top = other.ata_instruction_count_top;
  ata_instruction_count_lbl = other.ata_instruction_count_lbl;
  atom_starts = std::move(other.atom_starts);
  atom_counts = std::move(other.atom_counts);
  x_coordinates_sp = std::move(other.x_coordinates_sp);
  y_coordinates_sp = std::move(other.y_coordinates_sp);
  z_coordinates_sp = std::move(other.z_coordinates_sp);
  x_coordinates = std::move(other.x_coordinates);
  y_coordinates = std::move(other.y_coordinates);
  z_coordinates = std::move(other.z_coordinates);
  umat = std::move(other.umat);
  invu = std::move(other.invu);
  atr_instructions_src = std::move(other.atr_instructions_src);
  atr_instructions_top = std::move(other.atr_instructions_top);
  atr_instructions_lbl = std::move(other.atr_instructions_lbl);
  atr_instruction_groups_src = std::move(other.atr_instruction_groups_src);
  atr_instruction_groups_top = std::move(other.atr_instruction_groups_top);
  atr_instruction_groups_lbl = std::move(other.atr_instruction_groups_lbl);
  ata_instructions_src = std::move(other.ata_instructions_src);
  ata_instructions_top = std::move(other.ata_instructions_top);
  ata_instructions_lbl = std::move(other.ata_instructions_lbl);
  ata_instruction_groups_src = std::move(other.ata_instruction_groups_src);
  ata_instruction_groups_top = std::move(other.ata_instruction_groups_top);
  ata_instruction_groups_lbl = std::move(other.ata_instruction_groups_lbl);
  pps_ptr = std::move(other.pps_ptr);
  cs_ptr = std::move(other.cs_ptr);
  float_data = std::move(other.float_data);
  double_data = std::move(other.double_data);
  return *this;
}

//-------------------------------------------------------------------------------------------------
Condensate::Condensate(const PhaseSpaceSynthesis &poly_ps_in, const PrecisionModel mode_in,
                       const GpuDetails &gpu) :
    Condensate(poly_ps_in.getSelfPointer(), mode_in, gpu)
{}

//-------------------------------------------------------------------------------------------------
PrecisionModel Condensate::getMode() const {
  return mode;
}

//-------------------------------------------------------------------------------------------------
CondensateSource Condensate::getBasis() const {
  return basis;
}

//-------------------------------------------------------------------------------------------------
int Condensate::getSystemCount() const {
  return system_count;
}

//-------------------------------------------------------------------------------------------------
bool Condensate::ownsCoordinates() const {
  return holds_own_data;
}

//-------------------------------------------------------------------------------------------------
int Condensate::getATRInstructionCount(const SystemGrouping organization) const {
  switch (organization) {
  case SystemGrouping::SOURCE:
    return atr_instruction_count_src;
  case SystemGrouping::TOPOLOGY:
    return atr_instruction_count_top;
  case SystemGrouping::LABEL:
    return atr_instruction_count_lbl;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
int Condensate::getATAInstructionCount(const SystemGrouping organization) const {
  switch (organization) {
  case SystemGrouping::SOURCE:
    return ata_instruction_count_src;
  case SystemGrouping::TOPOLOGY:
    return ata_instruction_count_top;
  case SystemGrouping::LABEL:
    return ata_instruction_count_lbl;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::vector<int4> Condensate::getATRInstructions(const SystemGrouping organization) const {
  switch (organization) {
  case SystemGrouping::SOURCE:
    return atr_instructions_src.readHost();
  case SystemGrouping::TOPOLOGY:
    return atr_instructions_top.readHost();
  case SystemGrouping::LABEL:
    return atr_instructions_lbl.readHost();
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::vector<int> Condensate::getATRInstructionGroups(const SystemGrouping organization) const {
  switch (organization) {
  case SystemGrouping::SOURCE:
    return atr_instruction_groups_src.readHost();
  case SystemGrouping::TOPOLOGY:
    return atr_instruction_groups_top.readHost();
  case SystemGrouping::LABEL:
    return atr_instruction_groups_lbl.readHost();
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::vector<int4> Condensate::getATAInstructions(const SystemGrouping organization) const {
  switch (organization) {
  case SystemGrouping::SOURCE:
    return ata_instructions_src.readHost();
  case SystemGrouping::TOPOLOGY:
    return ata_instructions_top.readHost();
  case SystemGrouping::LABEL:
    return ata_instructions_lbl.readHost();
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::vector<int> Condensate::getATAInstructionGroups(const SystemGrouping organization) const {
  switch (organization) {
  case SystemGrouping::SOURCE:
    return ata_instruction_groups_src.readHost();
  case SystemGrouping::TOPOLOGY:
    return ata_instruction_groups_top.readHost();
  case SystemGrouping::LABEL:
    return ata_instruction_groups_lbl.readHost();
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
const PhaseSpaceSynthesis* Condensate::getSynthesisPointer() const {
  return pps_ptr;
}

//-------------------------------------------------------------------------------------------------
const CondensateReader Condensate::data(const HybridTargetLevel tier) const {
  return CondensateReader(mode, basis, system_count, atr_instruction_count_src,
                          atr_instruction_count_top, atr_instruction_count_lbl,
                          ata_instruction_count_src, ata_instruction_count_top,
                          ata_instruction_count_lbl, atom_starts.data(tier),
                          atom_counts.data(tier), x_coordinates_sp.data(tier),
                          y_coordinates_sp.data(tier), z_coordinates_sp.data(tier),
                          x_coordinates.data(tier), y_coordinates.data(tier),
                          z_coordinates.data(tier), umat.data(tier), invu.data(tier),
                          atr_instructions_src.data(tier), atr_instructions_top.data(tier),
                          atr_instructions_lbl.data(tier), atr_instruction_groups_src.data(tier),
                          atr_instruction_groups_top.data(tier),
                          atr_instruction_groups_lbl.data(tier), ata_instructions_src.data(tier),
                          ata_instructions_top.data(tier), ata_instructions_lbl.data(tier),
                          ata_instruction_groups_src.data(tier),
                          ata_instruction_groups_top.data(tier),
                          ata_instruction_groups_lbl.data(tier));
};  

//-------------------------------------------------------------------------------------------------
CondensateWriter Condensate::data(const HybridTargetLevel tier) {
  return CondensateWriter(mode, basis, system_count, atr_instruction_count_src,
                          atr_instruction_count_top, atr_instruction_count_lbl,
                          ata_instruction_count_src, ata_instruction_count_top,
                          ata_instruction_count_lbl, atom_starts.data(tier),
                          atom_counts.data(tier), x_coordinates_sp.data(tier),
                          y_coordinates_sp.data(tier), z_coordinates_sp.data(tier),
                          x_coordinates.data(tier), y_coordinates.data(tier),
                          z_coordinates.data(tier), umat.data(tier), invu.data(tier),
                          atr_instructions_src.data(tier), atr_instructions_top.data(tier),
                          atr_instructions_lbl.data(tier), atr_instruction_groups_src.data(tier),
                          atr_instruction_groups_top.data(tier),
                          atr_instruction_groups_lbl.data(tier), ata_instructions_src.data(tier),
                          ata_instructions_top.data(tier), ata_instructions_lbl.data(tier),
                          ata_instruction_groups_src.data(tier),
                          ata_instruction_groups_top.data(tier),
                          ata_instruction_groups_lbl.data(tier));
};

#ifdef STORMM_USE_HPC
//-------------------------------------------------------------------------------------------------
void Condensate::upload() {
  atom_starts.upload();
  atom_counts.upload();
  atr_instructions_src.upload();
  atr_instructions_top.upload();
  atr_instructions_lbl.upload();
  atr_instruction_groups_src.upload();
  atr_instruction_groups_top.upload();
  atr_instruction_groups_lbl.upload();
  ata_instructions_src.upload();
  ata_instructions_top.upload();
  ata_instructions_lbl.upload();
  ata_instruction_groups_src.upload();
  ata_instruction_groups_top.upload();
  ata_instruction_groups_lbl.upload();
  float_data.upload();
  double_data.upload();
}

//-------------------------------------------------------------------------------------------------
void Condensate::download() {
  atom_starts.download();
  atom_counts.download();
  atr_instructions_src.download();
  atr_instructions_top.download();
  atr_instructions_lbl.download();
  atr_instruction_groups_src.download();
  atr_instruction_groups_top.download();
  atr_instruction_groups_lbl.download();
  ata_instructions_src.download();
  ata_instructions_top.download();
  ata_instructions_lbl.download();
  ata_instruction_groups_src.download();
  ata_instruction_groups_top.download();
  ata_instruction_groups_lbl.download();
  float_data.download();
  double_data.download();
}
#endif

//-------------------------------------------------------------------------------------------------
void Condensate::rebuild(const PhaseSpaceSynthesis *poly_ps_in, const PrecisionModel mode_in,
                         const GpuDetails &gpu) {

  // The basis must be set to a synthesis, even if that synthesis is the nullptr.
  basis = CondensateSource::SYNTHESIS;
  
  // Exit if the synthesis is the null pointer.  Only PhaseSpaceSynthesis objects will come in as
  // nullptr
  if (poly_ps_in == nullptr) {
    system_count = 0;
    atom_starts.resize(0);
    atom_counts.resize(0);
    float_data.resize(0);
    float_data.shrinkToFit();
    double_data.resize(0);
    double_data.shrinkToFit();
    x_coordinates.setPointer(&double_data, 0, 0);
    y_coordinates.setPointer(&double_data, 0, 0);
    z_coordinates.setPointer(&double_data, 0, 0);
    x_coordinates_sp.setPointer(&float_data, 0, 0);
    y_coordinates_sp.setPointer(&float_data, 0, 0);
    z_coordinates_sp.setPointer(&float_data, 0, 0);
    umat.setPointer(&double_data, 0, 0);
    invu.setPointer(&double_data, 0, 0);
    atr_instruction_count_src = 0;
    atr_instruction_count_top = 0;
    atr_instruction_count_lbl = 0;
    ata_instruction_count_src = 0;
    ata_instruction_count_top = 0;
    ata_instruction_count_lbl = 0;
    atr_instructions_src.resize(0);
    atr_instructions_top.resize(0);
    atr_instructions_lbl.resize(0);
    atr_instruction_groups_src.resize(0);
    atr_instruction_groups_top.resize(0);
    atr_instruction_groups_lbl.resize(0);
    ata_instructions_src.resize(0);
    ata_instructions_top.resize(0);
    ata_instructions_lbl.resize(0);
    ata_instruction_groups_src.resize(0);
    ata_instruction_groups_top.resize(0);
    ata_instruction_groups_lbl.resize(0);
    return;
  }

  // Build based on the synthesis
  pps_ptr = const_cast<PhaseSpaceSynthesis*>(poly_ps_in);
  cs_ptr = nullptr;
  const PsSynthesisReader poly_psr = pps_ptr->data();
  system_count = poly_psr.system_count;
  atom_starts.resize(system_count);
  atom_counts.resize(system_count);
  size_t* atom_starts_ptr = atom_starts.data();
  int*    atom_counts_ptr = atom_counts.data();
  for (int i = 0; i < system_count; i++) {
    atom_starts_ptr[i] = poly_psr.atom_starts[i];
    atom_counts_ptr[i] = poly_psr.atom_counts[i];
  }
  const int last_sys = poly_psr.system_count - 1;
  const size_t padded_atoms = poly_psr.atom_starts[last_sys] +
                              roundUp(poly_psr.atom_counts[last_sys], warp_size_int);
  const size_t xfrm_spacing = system_count * roundUp<size_t>(9, warp_size_zu);
  mode = mode_in;
  holds_own_data = true;
  switch (mode) {
  case PrecisionModel::DOUBLE:
    float_data.resize(0);
    float_data.shrinkToFit();
    x_coordinates_sp.setPointer(&float_data, 0, 0);
    y_coordinates_sp.setPointer(&float_data, 0, 0);
    z_coordinates_sp.setPointer(&float_data, 0, 0);
    double_data.resize((3LLU * padded_atoms) + (2 * xfrm_spacing));
    double_data.shrinkToFit();
    x_coordinates.setPointer(&double_data,                   0, padded_atoms);
    y_coordinates.setPointer(&double_data,        padded_atoms, padded_atoms);
    z_coordinates.setPointer(&double_data, 2LLU * padded_atoms, padded_atoms);
    umat.setPointer(&double_data,  3LLU * padded_atoms                , xfrm_spacing);
    invu.setPointer(&double_data, (3LLU * padded_atoms) + xfrm_spacing, xfrm_spacing);
    break;
  case PrecisionModel::SINGLE:
    float_data.resize(3LLU * padded_atoms);
    float_data.shrinkToFit();
    x_coordinates_sp.setPointer(&float_data,                   0, padded_atoms);
    y_coordinates_sp.setPointer(&float_data,        padded_atoms, padded_atoms);
    z_coordinates_sp.setPointer(&float_data, 2LLU * padded_atoms, padded_atoms);
    double_data.resize(2 * xfrm_spacing);
    double_data.shrinkToFit();
    x_coordinates.setPointer(&double_data, 0, 0);
    y_coordinates.setPointer(&double_data, 0, 0);
    z_coordinates.setPointer(&double_data, 0, 0);
    umat.setPointer(&double_data, 0, xfrm_spacing);
    invu.setPointer(&double_data, xfrm_spacing, xfrm_spacing);
    break;
  }

  // Map the work units.  Condensates based on CoordinateSeries objects have a similar process,
  // but with many topologies to consider the PhaseSpaceSynthesis is different enough that the
  // builder has its own code.
  computeTopologyWorkUnits(gpu);
  update();
}

//-------------------------------------------------------------------------------------------------
void Condensate::rebuild(const PhaseSpaceSynthesis &poly_ps_in, const PrecisionModel mode_in,
                         const GpuDetails &gpu) {
  rebuild(poly_ps_in.getSelfPointer(), mode_in, gpu);
}

//-------------------------------------------------------------------------------------------------
void Condensate::setWorkUnits(const int* system_list, const int* bounds_list, const int partitions,
                              const SystemGrouping organization, const GpuDetails &gpu) {
  if (pps_ptr == nullptr) {
    rtErr("Work units (other than those for a partition by topology, which is automatic for any "
          "Condensate object) may not be computed without building based on a coordinate "
          "synthesis.", "Condensate", "setWorkUnits");
  }
  const PsSynthesisReader poly_psr = pps_ptr->data();
  switch (organization) {
  case SystemGrouping::SOURCE:
    generateWorkUnits(system_list, poly_psr.unique_ag_idx, bounds_list, partitions,
                      &atr_instructions_src, &atr_instruction_groups_src, &ata_instructions_src,
                      &ata_instruction_groups_src, &atr_instruction_count_src,
                      &ata_instruction_count_src, gpu);
    break;
  case SystemGrouping::TOPOLOGY:
    rtErr("Work units for a partition of the systems based on a topology are automatically "
          "comuted when building from any synthesis.", "Condensate", "setWorkUnits");
  case SystemGrouping::LABEL:
    generateWorkUnits(system_list, poly_psr.unique_ag_idx, bounds_list, partitions,
                      &atr_instructions_lbl, &atr_instruction_groups_lbl, &ata_instructions_lbl,
                      &ata_instruction_groups_lbl, &atr_instruction_count_lbl,
                      &ata_instruction_count_lbl, gpu);
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void Condensate::setWorkUnits(const std::vector<int> &system_list,
                              const std::vector<int> &bounds_list,
                              const SystemGrouping organization, const GpuDetails &gpu) {
  setWorkUnits(system_list.data(), bounds_list.data(), bounds_list.size() - 1, organization, gpu);
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
            const size_t llim = poly_psr.atom_starts[i];
            const size_t hlim = llim + poly_psr.atom_counts[i];
            for (size_t j = llim; j < hlim; j++) {
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
              case PrecisionModel::DOUBLE:
                cdw.xcrd[j] = xij;
                cdw.ycrd[j] = yij;
                cdw.zcrd[j] = zij;
                break;
              case PrecisionModel::SINGLE:
                cdw.xcrd_sp[j] = xij;
                cdw.ycrd_sp[j] = yij;
                cdw.zcrd_sp[j] = zij;
                break;
              }
            }
            const size_t xllim = static_cast<size_t>(i) * roundUp<size_t>(9, warp_size_zu);
            const size_t xhlim = xllim + 9LLU;
            for (size_t j = xllim; j < xhlim; j++) {
              cdw.umat[j] = poly_psr.umat[j];
              cdw.invu[j] = poly_psr.invu[j];
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
    umat.swapTarget(&double_data);
    invu.swapTarget(&double_data);
  }
}

//-------------------------------------------------------------------------------------------------
void Condensate::generateWorkUnits(const int* system_list, const int* topology_index_list,
                                   const int* bounds_list, const int partitions,
                                   Hybrid<int4> *atr_instructions,
                                   Hybrid<int> *atr_instruction_groups,
                                   Hybrid<int4> *ata_instructions,
                                   Hybrid<int> *ata_instruction_groups, int *atr_instruction_count,
                                   int *ata_instruction_count, const GpuDetails &gpu) {
  const int system_count = bounds_list[partitions];
  
  // Compute the number of warps on the GPU and the number of bits per unsigned int
  const llint gpu_nwarp = (gpu == null_gpu) ?
                          large_block_size / warp_size_int :
                          gpu.getSMPCount() * (gpu.getMaxThreadsPerBlock() / warp_size_int);
  const int uint_bits = sizeof(uint) * 8;
  const int half_int_bits = uint_bits / 2;
  const llint np_atr = system_count;
  llint np_ata = 0LL;
  for (int i = 0; i < partitions; i++) {
    const llint ni = bounds_list[i + 1] - bounds_list[i];
    np_ata += ni * ni;
  }
  llint atr_per_warp = ceil(static_cast<double>(np_atr) / static_cast<double>(gpu_nwarp));
  atr_per_warp = std::min(static_cast<llint>(uint_bits), atr_per_warp);
  llint ata_per_warp = ceil(static_cast<double>(np_ata) / static_cast<double>(gpu_nwarp));
  const int ata_stride = std::min(static_cast<llint>(half_int_bits),
                                  static_cast<llint>(ceil(sqrt(ata_per_warp))));;
  ata_per_warp = ata_stride * ata_stride;
  
  // Lay out a temporary array for all-to-one instructions
  std::vector<int4> tmp_atr_instructions;
  std::vector<int> tmp_atr_instruction_groups;
  tmp_atr_instructions.reserve((np_atr / atr_per_warp) + partitions);
  tmp_atr_instruction_groups.reserve(tmp_atr_instructions.size());
  for (int i = 0; i < partitions; i++) {
    const int group_lower_bound = bounds_list[i];
    const int group_upper_bound = bounds_list[i + 1];
    const int top_idx = topology_index_list[system_list[group_lower_bound]];
    for (int j = group_lower_bound; j < group_upper_bound; j += atr_per_warp) {
      const int nrep = std::min(static_cast<int>(atr_per_warp), group_upper_bound - j);
      tmp_atr_instructions.push_back({ j - group_lower_bound, i, top_idx, nrep });
      tmp_atr_instruction_groups.push_back(i);
    }
  }
  *atr_instruction_count = tmp_atr_instructions.size();

  // Lay out a temporary array for all-to-all instructions.
  std::vector<int4> tmp_ata_instructions;
  std::vector<int> tmp_ata_instruction_groups;
  size_t ata_instruction_guess = 0;
  for (int i = 0; i < partitions; i++) {
    const size_t nstride = (bounds_list[i + 1] - bounds_list[i] + ata_stride - 1) / ata_stride;
    ata_instruction_guess += nstride * (nstride + 1) / 2;
  }
  tmp_ata_instructions.reserve(ata_instruction_guess);
  tmp_ata_instruction_groups.reserve(ata_instruction_guess);
  for (int i = 0; i < partitions; i++) {
    const int group_lower_bound = bounds_list[i];
    const int group_upper_bound = bounds_list[i + 1];
    const int top_idx = topology_index_list[system_list[group_lower_bound]];
    for (int j = group_lower_bound; j < group_upper_bound; j += ata_stride) {
      const int njrep = (j + ata_stride < group_upper_bound) ? ata_stride : group_upper_bound - j;
      for (int k = group_lower_bound; k <= j; k += ata_stride) {
        const int nkrep = (k + ata_stride < group_upper_bound) ? ata_stride :
                                                                 group_upper_bound - k;
        tmp_ata_instructions.push_back({ j - group_lower_bound, k - group_lower_bound, top_idx,
                                         (nkrep << half_int_bits) | njrep });
        tmp_ata_instruction_groups.push_back(i);
      }
    }
  }
  *ata_instruction_count = tmp_ata_instructions.size();

  // Load the temporary arrays into the object.
  atr_instructions->resize(*atr_instruction_count);
  ata_instructions->resize(*ata_instruction_count);
  atr_instructions->putHost(tmp_atr_instructions);
  ata_instructions->putHost(tmp_ata_instructions);
  atr_instruction_groups->resize(*atr_instruction_count);
  ata_instruction_groups->resize(*ata_instruction_count);
  atr_instruction_groups->putHost(tmp_atr_instruction_groups);
  ata_instruction_groups->putHost(tmp_ata_instruction_groups);
}

//-------------------------------------------------------------------------------------------------
void Condensate::computeTopologyWorkUnits(const GpuDetails &gpu) {
  const PsSynthesisReader poly_psr = pps_ptr->data();
  generateWorkUnits(poly_psr.common_ag_list, poly_psr.unique_ag_idx, poly_psr.common_ag_bounds,
                    poly_psr.unique_topology_count, &atr_instructions_top,
                    &atr_instruction_groups_top, &ata_instructions_top,
                    &ata_instruction_groups_top, &atr_instruction_count_top,
                    &ata_instruction_count_top, gpu);
}

} // namespace synthesis
} // namespace stormm
