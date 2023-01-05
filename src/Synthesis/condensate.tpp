// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace synthesis {

//-------------------------------------------------------------------------------------------------
template <typename T>
Condensate::Condensate(const CoordinateSeries<T> *cs_in, const CondensationLevel mode_in,
                       const GpuDetails &gpu) :
    Condensate(nullptr, mode_in, gpu)
{
  rebuild(cs_in, mode, gpu);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
Condensate::Condensate(const CoordinateSeries<T> &cs_in, const CondensationLevel mode_in,
                       const GpuDetails &gpu) :
    Condensate(cs_in.getSelfPointer(), mode_in, gpu)
{}

//-------------------------------------------------------------------------------------------------
template <typename T>
const CoordinateSeries<T>* Condensate::getSeriesPointer() const {
  return reinterpret_cast<CoordinateSeries<T>*>(cs_ptr);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void Condensate::rebuild(const CoordinateSeries<T> *cs_in, const CondensationLevel mode_in,
                         const GpuDetails &gpu) {
  
  // Reinterpret the templated CoordinatesSeries<T> pointer as a CoordinateSeries of an arbitrary
  // type.  This prevents the Condensate class as a whole from taking on a template requirement.
  cs_ptr = reinterpret_cast<CoordinateSeries<int>*>(const_cast<CoordinateSeries<T>*>(cs_in));
  pps_ptr = nullptr;
  csptr_data_type = std::type_index(typeid(T)).hash_code();
  const CoordinateSeriesReader<T> csr = cs_in->data();
  const size_t padded_atoms = static_cast<size_t>(csr.nframe) *
                              roundUp(static_cast<size_t>(csr.natom), warp_size_zu);
  mode = mode_in;
  const size_t ct = std::type_index(typeid(T)).hash_code();

  // If the CoordinateSeries is of type float or double, the explicit instantiation will have been
  // taken.
  holds_own_data = false;
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
  update(cs_in);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void Condensate::listWorkInstructions(const CoordinateSeries<T> *cs_in, const GpuDetails &gpu) {
  const CoordinateSeriesReader<T> csr = cs_in->data();
  if (gpu != null_gpu) {

    // Estimate the number of comparisons for the all-to-reference and all-to-all cases based on
    // the raw number of systems using each topology.
    const llint nframe_ll = csr.nframe;
    llint np_atr = nframe_ll;
    llint np_ata = np_atr * np_atr;
    const llint nsmp = gpu.getSMPCount();
    np_atr = roundUp(np_atr, nsmp) / nsmp;
    np_atr = roundUp(np_atr, static_cast<llint>(medium_block_size / warp_size_int));
    if (np_ata > 32768) {
      np_ata = 16;
    }
    else {
      np_ata = 8;
    }
    const int atr_strides = roundUp(nframe_ll, np_atr) / np_atr;
    const int nwu_atr = atr_strides;
    const int ata_strides = roundUp(nframe_ll, np_ata) / np_ata; 
    const llint nwu_ata = static_cast<llint>(ata_strides) *
                          (static_cast<llint>(ata_strides) + 1LL) / 2LL;
    atr_instruction_count = nwu_atr;
    ata_instruction_count = nwu_ata;
    atr_instructions.resize(atr_instruction_count);
    ata_instructions.resize(ata_instruction_count);
    int4* atr_insr_ptr = atr_instructions.data();
    int4* ata_insr_ptr = ata_instructions.data();
    const int inp_atr = np_atr;
    for (int i = 0; i < atr_strides; i++) {
      const int rdims = std::min(csr.nframe - (i * inp_atr), inp_atr);
      const int4 rtmp = { i * inp_atr, 0, 0, rdims };
      atr_insr_ptr[i] = rtmp;
    }
    size_t ata_insr_idx = 0;
    const int inp_ata = np_ata;
    for (int i = 0LL; i < ata_strides; i++) {
      const int xdims = std::min(csr.nframe - (i * np_ata), np_ata);
      for (int j = 0; j <= i; j++) {
        const int ydims = std::min(csr.nframe - (j * inp_ata), inp_ata);
        const int4 mtmp = { i * inp_ata, j * inp_ata, 0, ((ydims << 16) | xdims) };
        ata_insr_ptr[ata_insr_idx] = mtmp;
        ata_insr_idx++;
      }
    }
  }
  else {
    atr_instruction_count = roundUp(csr.nframe, 1024) / 1024;
    ata_instruction_count = atr_instruction_count * (atr_instruction_count + 1) / 2;
    atr_instructions.resize(atr_instruction_count);
    ata_instructions.resize(ata_instruction_count);
    int4* atr_insr_ptr = atr_instructions.data();
    int4* ata_insr_ptr = ata_instructions.data();
    int ata_insr_idx = 0;
    for (int i = 0; i < atr_instruction_count; i++) {
      const int xdim = std::min(csr.nframe - (i * 1024), 1024);
      const int4 rtmp = { i * 1024, 0, i, xdim };
      for (int j = 0; j <= i; j++) {
        const int ydim = (j == i) ? xdim : 1024;
        const int4 mtmp = { i * 1024, j * 1024, 0, ((ydim << 16) | xdim) };
        ata_insr_ptr[ata_insr_idx] = mtmp;
        ata_insr_idx++;
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void Condensate::update(const CoordinateSeries<T> *cs_basis, const HybridTargetLevel tier,
                        const GpuDetails &gpu) {

  // Do not perform any updates if the Condensate already points to the CoordinateSeries
  if (holds_own_data == false) {
    return;
  }

  // Produce an error if the Condensate is not ready to accept the CoordinateSeries in question.
  if (cs_basis != reinterpret_cast<CoordinateSeries<T>*>(cs_ptr)) {
    rtErr("The pointer to the CoordinateSeries supplied does not match the pointer stored "
          "internally.", "Condensate", "Update");
  }

  // Copy the data
  switch (tier) {
  case HybridTargetLevel::HOST:
    {
      const CoordinateSeriesReader<T> csr = cs_basis->data(tier);
      CondensateWriter cdw = this->data(tier);
      const size_t natom_zu = csr.natom;
      const size_t padded_natom = roundUp(natom_zu, warp_size_zu);
      const size_t nframe_zu = csr.nframe;
      for (size_t i = 0; i < nframe_zu; i++) {
        const size_t llim = i * padded_natom;
        const size_t hlim = llim + natom_zu;
        for (size_t j = llim; j < hlim; j++) {
          double xij = csr.xcrd[j] * csr.inv_gpos_scale;
          double yij = csr.ycrd[j] * csr.inv_gpos_scale;
          double zij = csr.zcrd[j] * csr.inv_gpos_scale;
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

//-------------------------------------------------------------------------------------------------
template <typename T>
void Condensate::update(const CoordinateSeries<T> &cs_basis, const HybridTargetLevel tier,
                        const GpuDetails &gpu) {
  update(cs_basis.getSelfPointer(), tier, gpu);
}

} // namespace synthesis
} // namespace stormm
