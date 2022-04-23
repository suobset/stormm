#include <cmath>
#include "DataTypes/common_types.h"
#include "Math/rounding.h"
#include "scorecard.h"

namespace omni {
namespace energy {

using math::roundUp;
using numerics::default_energy_scale_f;
using numerics::default_energy_scale_lf;
using numerics::default_inverse_energy_scale_f;
using numerics::default_inverse_energy_scale_lf;

//-------------------------------------------------------------------------------------------------
ScoreCardReader::ScoreCardReader(const int system_count_in, const int data_stride_in,
                                 const int sampled_step_count_in, const float nrg_scale_f_in,
                                 const double nrg_scale_lf_in, const float inverse_nrg_scale_f_in,
                                 const double inverse_nrg_scale_lf_in,
                                 const llint* instantaneous_accumulators_in,
                                 const double* running_accumulators_in,
                                 const double* squared_accumulators_in,
                                 const llint* time_series_in) :
    system_count{system_count_in}, data_stride{data_stride_in},
    sampled_step_count{sampled_step_count_in}, nrg_scale_f{nrg_scale_f_in},
    nrg_scale_lf{nrg_scale_lf_in}, inverse_nrg_scale_f{inverse_nrg_scale_f_in},
    inverse_nrg_scale_lf{inverse_nrg_scale_lf_in},
    instantaneous_accumulators{instantaneous_accumulators_in},
    running_accumulators{running_accumulators_in},
    squared_accumulators{squared_accumulators_in},
    time_series{time_series_in}
{}

//-------------------------------------------------------------------------------------------------
ScoreCardWriter::ScoreCardWriter(const int system_count_in, const int data_stride_in,
                                 const int sampled_step_count_in, const float nrg_scale_f_in,
                                 const double nrg_scale_lf_in, const float inverse_nrg_scale_f_in,
                                 const double inverse_nrg_scale_lf_in,
                                 llint* instantaneous_accumulators_in,
                                 double* running_accumulators_in, double* squared_accumulators_in,
                                 llint* time_series_in) :
    system_count{system_count_in}, data_stride{data_stride_in},
    sampled_step_count{sampled_step_count_in}, nrg_scale_f{nrg_scale_f_in},
    nrg_scale_lf{nrg_scale_lf_in}, inverse_nrg_scale_f{inverse_nrg_scale_f_in},
    inverse_nrg_scale_lf{inverse_nrg_scale_lf_in},
    instantaneous_accumulators{instantaneous_accumulators_in},
    running_accumulators{running_accumulators_in},
    squared_accumulators{squared_accumulators_in},
    time_series{time_series_in}
{}

//-------------------------------------------------------------------------------------------------
ScoreCard::ScoreCard(const int system_count_in, const int capacity_in,
                     const int nrg_scale_bits_in) :
    system_count{system_count_in},
    data_stride{roundUp(static_cast<int>(StateVariable::ALL_STATES), warp_size_int)},
    sampled_step_count{0},
    sample_capacity{capacity_in},
    nrg_scale_bits{nrg_scale_bits_in},
    nrg_scale_f{0.0},
    nrg_scale_lf{0.0},
    inverse_nrg_scale_f{0.0},
    inverse_nrg_scale_lf{0.0},
    instantaneous_accumulators{static_cast<size_t>(data_stride * system_count), "scorecard_acc"},
    running_accumulators{static_cast<size_t>(data_stride * system_count), "score_running_acc"},
    squared_accumulators{static_cast<size_t>(data_stride * system_count), "score_squared_acc"},
    time_series_accumulators{static_cast<size_t>(data_stride * system_count * sample_capacity),
                             "score_card_ts"}
{
  // Check the fixed precision bit count: there are limits that the user or developer will be
  // allowed, for the sake of energy estimates with some degree of accuracy and remaining within
  // the bounds of the long long integer accumulator format.
  if (nrg_scale_bits < 11) {
    rtErr("Energy and virial accumulation must take place in a precision of at least one part in "
          "2048 of one kcal/mol.  A precision of " + std::to_string(nrg_scale_bits) +
          " bits after the decimal is unacceptably low.", "ScoreCard");
  }
  else if (nrg_scale_bits > 40) {
    rtErr("Energy and virial accumulation can overflow the 64-bit integer accumulators.  Storing "
          "the numbers to a precision of more than one part in 1.0995 trillionths of a kcal/mol, "
          "a precision of " + std::to_string(nrg_scale_bits) + " after the decimal, is "
          "unnecessary and risks producing overflow in large systems.", "ScoreCard");
  }

  // Determine the energy scaling factors
  if (nrg_scale_bits == default_energy_scale_bits) {
    nrg_scale_lf = default_energy_scale_lf;
    nrg_scale_f = default_energy_scale_f;
    inverse_nrg_scale_lf = default_inverse_energy_scale_lf;
    inverse_nrg_scale_f = default_inverse_energy_scale_f;
  }
  else {
    int ib = 0;
    nrg_scale_lf = 1.0;
    while (ib + 10 <= nrg_scale_bits) {
      nrg_scale_lf *= 1024.0;
      ib += 10;
    }
    while (ib < nrg_scale_bits) {
      nrg_scale_lf *= 2.0;
      ib++;
    }
    nrg_scale_f = (float)nrg_scale_lf;
    inverse_nrg_scale_lf = 1.0 / nrg_scale_lf;
    inverse_nrg_scale_f = (float)1.0 / nrg_scale_f;
  }
}

//-------------------------------------------------------------------------------------------------
int ScoreCard::getSystemCount() const {
  return system_count;
}

//-------------------------------------------------------------------------------------------------
int ScoreCard::getSampleSize() const {
  return sampled_step_count;
}

//-------------------------------------------------------------------------------------------------
int ScoreCard::getEnergyScaleBits() const {
  return nrg_scale_bits;
}

#ifdef OMNI_USE_HPC
//-------------------------------------------------------------------------------------------------
void ScoreCard::upload() {
  instantaneous_accumulators.upload();
  running_accumulators.upload();
  squared_accumulators.upload();
  time_series_accumulators.upload();
}

//-------------------------------------------------------------------------------------------------
void ScoreCard::download() {
  instantaneous_accumulators.download();
  running_accumulators.download();
  squared_accumulators.download();
  time_series_accumulators.download();
}
#endif

//-------------------------------------------------------------------------------------------------
const ScoreCardReader ScoreCard::data(const HybridTargetLevel tier) const {
  return ScoreCardReader(system_count, data_stride, sampled_step_count, nrg_scale_f, nrg_scale_lf,
                         inverse_nrg_scale_f, inverse_nrg_scale_lf,
                         instantaneous_accumulators.data(tier), running_accumulators.data(tier),
                         squared_accumulators.data(tier), time_series_accumulators.data(tier));
}

//-------------------------------------------------------------------------------------------------
ScoreCardWriter ScoreCard::data(const HybridTargetLevel tier) {
  return ScoreCardWriter(system_count, data_stride, sampled_step_count, nrg_scale_f, nrg_scale_lf,
                         inverse_nrg_scale_f, inverse_nrg_scale_lf,
                         instantaneous_accumulators.data(tier), running_accumulators.data(tier),
                         squared_accumulators.data(tier), time_series_accumulators.data(tier));
}

//-------------------------------------------------------------------------------------------------
void ScoreCard::reserve(const int new_capacity) {
  if (new_capacity > sample_capacity) {

    // Allocate one space beyond the requested capacity, so that reserve() can be called at the
    // outset of a molecular dynamics simulation for the projected number of sampled steps and
    // not trigger a resize when the final step gets logged.
    sample_capacity = new_capacity + 1;
    time_series_accumulators.resize(static_cast<size_t>(system_count * data_stride) *
                                    static_cast<size_t>(sample_capacity));
  }
}

//-------------------------------------------------------------------------------------------------
void ScoreCard::contribute(const StateVariable var, const llint amount, const int system_index) {
  const size_t slot = static_cast<int>(var) + (data_stride * system_index);
  instantaneous_accumulators.putHost(amount, slot);
  const double dbl_amount = static_cast<double>(amount) * inverse_nrg_scale_lf;
  running_accumulators.putHost(running_accumulators.readHost(slot) + dbl_amount, slot);
  squared_accumulators.putHost(squared_accumulators.readHost(slot) + dbl_amount * dbl_amount,
                               slot);
  const size_t ts_slot = slot + (static_cast<size_t>(data_stride * system_index) *
                                 static_cast<size_t>(sampled_step_count));
  time_series_accumulators.putHost(amount, ts_slot);
}

//-------------------------------------------------------------------------------------------------
void ScoreCard::incrementSampleCount() {
  sampled_step_count += 1;
  if (sampled_step_count == sample_capacity) {
    sample_capacity *= 2;
    time_series_accumulators.resize(static_cast<size_t>(system_count * data_stride) *
                                    static_cast<size_t>(sample_capacity));
  }
}

//-------------------------------------------------------------------------------------------------
std::vector<double> ScoreCard::reportTotalEnergies() {
#ifdef OMNI_USE_HPC
  instantaneous_accumulators.download();
#endif
  const int nvar = static_cast<int>(StateVariable::ALL_STATES);
  const int padded_nvar = roundUp(nvar, warp_size_int);
  std::vector<double> result(system_count);
  const llint* inst_acc_ptr = instantaneous_accumulators.data();
  for (int i = 0; i < system_count; i++) {
    llint lacc = 0LL;
    for (int j = 0; j < nvar; j++) {
      lacc += inst_acc_ptr[(i * padded_nvar) + j];
    }
    result[i] = inverse_nrg_scale_lf * static_cast<double>(lacc);
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
double ScoreCard::reportTotalEnergy(const int system_index) {
#ifdef OMNI_USE_HPC
  instantaneous_accumulators.download();
#endif
  const int nvar = static_cast<int>(StateVariable::ALL_STATES);
  const int padded_nvar = roundUp(nvar, warp_size_int);
  const llint* inst_acc_ptr = instantaneous_accumulators.data();
  llint lacc = 0LL;
  for (int i = 0; i < nvar; i++) {
    lacc += inst_acc_ptr[(system_index * padded_nvar) + i];
  }
  return inverse_nrg_scale_lf * static_cast<double>(lacc);
}

//-------------------------------------------------------------------------------------------------
std::vector<double> ScoreCard::reportInstantaneousStates() {
#ifdef OMNI_USE_HPC
  instantaneous_accumulators.download();
#endif
  const int nvar = static_cast<int>(StateVariable::ALL_STATES);
  const int padded_nvar = roundUp(nvar, warp_size_int);
  std::vector<double> result(nvar * system_count);
  const llint* inst_acc_ptr = instantaneous_accumulators.data();
  for (int i = 0; i < system_count; i++) {
    for (int j = 0; j < nvar; j++) {
      const llint lacc = inst_acc_ptr[(i * padded_nvar) + j];
      result[(i * nvar) + j] = inverse_nrg_scale_lf * static_cast<double>(lacc);
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<double> ScoreCard::reportInstantaneousStates(const int system_index) {
#ifdef OMNI_USE_HPC
  instantaneous_accumulators.download();
#endif
  const int nvar = static_cast<int>(StateVariable::ALL_STATES);
  const int offset = system_index * roundUp(nvar, warp_size_int);
  std::vector<double> result(nvar);
  const llint* inst_acc_ptr = instantaneous_accumulators.data();
  for (int i = 0; i < nvar; i++) {
    result[i] = inverse_nrg_scale_lf * static_cast<double>(inst_acc_ptr[offset + i]);
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<double> ScoreCard::reportInstantaneousStates(const StateVariable aspect) {
#ifdef OMNI_USE_HPC
  instantaneous_accumulators.download();
#endif
  const int aspect_no = static_cast<int>(aspect);
  const int nvar = static_cast<int>(StateVariable::ALL_STATES);
  const int padded_nvar = roundUp(nvar, warp_size_int);
  std::vector<double> result(system_count);
  const llint* inst_acc_ptr = instantaneous_accumulators.data();
  for (int i = 0; i < system_count; i++) {
    const llint lacc = inst_acc_ptr[(i * padded_nvar) + aspect_no];
    result[i] = inverse_nrg_scale_lf * static_cast<double>(lacc);
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
double ScoreCard::reportInstantaneousStates(const StateVariable aspect, const int system_index) {
#ifdef OMNI_USE_HPC
  instantaneous_accumulators.download();
#endif
  const int aspect_no = static_cast<int>(aspect);
  const int nvar = static_cast<int>(StateVariable::ALL_STATES);
  const int padded_nvar = roundUp(nvar, warp_size_int);
  const llint lacc = instantaneous_accumulators.readHost((system_index * padded_nvar) + aspect_no);
  return inverse_nrg_scale_lf * static_cast<double>(lacc);
}

//-------------------------------------------------------------------------------------------------
std::vector<double> ScoreCard::reportAverageStates() {
#ifdef OMNI_USE_HPC
  running_accumulators.download();
#endif
  const int nvar = static_cast<int>(StateVariable::ALL_STATES);
  const int padded_nvar = roundUp(nvar, warp_size_int);
  std::vector<double> result(nvar * system_count);
  const double* run_acc_ptr = running_accumulators.data();
  const double nsamp = static_cast<double>(sampled_step_count);
  for (int i = 0; i < system_count; i++) {
    for (int j = 0; j < nvar; j++) {
      result[(i * nvar) + j] = run_acc_ptr[(i * padded_nvar) + j] / nsamp;
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<double> ScoreCard::reportAverageStates(const int system_index) {
#ifdef OMNI_USE_HPC
  running_accumulators.download();
#endif
  const int nvar = static_cast<int>(StateVariable::ALL_STATES);
  const int offset = system_index * roundUp(nvar, warp_size_int);
  std::vector<double> result(nvar);
  const double* run_acc_ptr = running_accumulators.data();
  const double nsamp = static_cast<double>(sampled_step_count);
  for (int i = 0; i < nvar; i++) {
    result[i] = run_acc_ptr[offset + i] / nsamp;
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<double> ScoreCard::reportAverageStates(const StateVariable aspect) {
#ifdef OMNI_USE_HPC
  running_accumulators.download();
#endif
  const int aspect_no = static_cast<int>(aspect);
  const int nvar = static_cast<int>(StateVariable::ALL_STATES);
  const int padded_nvar = roundUp(nvar, warp_size_int);
  std::vector<double> result(system_count);
  const double* run_acc_ptr = running_accumulators.data();
  const double nsamp = static_cast<double>(sampled_step_count);
  for (int i = 0; i < system_count; i++) {
    result[i] = run_acc_ptr[(i * padded_nvar) + aspect_no] / nsamp;
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
double ScoreCard::reportAverageStates(const StateVariable aspect, const int system_index) {
#ifdef OMNI_USE_HPC
  running_accumulators.download();
#endif
  const int aspect_no = static_cast<int>(aspect);
  const int nvar = static_cast<int>(StateVariable::ALL_STATES);
  const int offset = system_index * roundUp(nvar, warp_size_int);
  const double nsamp = static_cast<double>(sampled_step_count);
  return running_accumulators.readHost(offset + aspect_no) / nsamp;
}

//-------------------------------------------------------------------------------------------------
std::vector<double> ScoreCard::reportVarianceOfStates() {
#ifdef OMNI_USE_HPC
  running_accumulators.download();
  squared_accumulators.download();
#endif
  const int nvar = static_cast<int>(StateVariable::ALL_STATES);
  const int padded_nvar = roundUp(nvar, warp_size_int);
  std::vector<double> result(nvar * system_count, 0.0);
  const double* run_acc_ptr = running_accumulators.data();
  const double* sqr_acc_ptr = squared_accumulators.data();
  const double nsamp = static_cast<double>(sampled_step_count);
  for (int i = 0; i < system_count; i++) {
    for (int j = 0; j < nvar; j++) {
      const double s1 = run_acc_ptr[(i * padded_nvar) + j];
      const double s2 = sqr_acc_ptr[(i * padded_nvar) + j];
      result[(i * nvar) + j] = sqrt((nsamp * s2) - (s1 * s1)) / sqrt(nsamp * (nsamp - 1.0));
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<double> ScoreCard::reportVarianceOfStates(const int system_index) {
#ifdef OMNI_USE_HPC
  running_accumulators.download();
  squared_accumulators.download();
#endif
  const int nvar = static_cast<int>(StateVariable::ALL_STATES);
  const int offset = system_index * roundUp(nvar, warp_size_int);
  std::vector<double> result(nvar);
  const double* run_acc_ptr = running_accumulators.data();
  const double* sqr_acc_ptr = squared_accumulators.data();
  const double nsamp = static_cast<double>(sampled_step_count);
  for (int i = 0; i < nvar; i++) {
    const double s1 = run_acc_ptr[offset + i];
    const double s2 = sqr_acc_ptr[offset + i];
    result[i] = sqrt((nsamp * s2) - (s1 * s1)) / sqrt(nsamp * (nsamp - 1.0));
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<double> ScoreCard::reportVarianceOfStates(const StateVariable aspect) {
#ifdef OMNI_USE_HPC
  running_accumulators.download();
  squared_accumulators.download();
#endif
  const int aspect_no = static_cast<int>(aspect);
  const int nvar = static_cast<int>(StateVariable::ALL_STATES);
  const int padded_nvar = roundUp(nvar, warp_size_int);
  std::vector<double> result(system_count, 0.0);
  const double* run_acc_ptr = running_accumulators.data();
  const double* sqr_acc_ptr = squared_accumulators.data();
  const double nsamp = static_cast<double>(sampled_step_count);
  for (int i = 0; i < system_count; i++) {
    const double s1 = run_acc_ptr[(i * padded_nvar) + aspect_no];
    const double s2 = sqr_acc_ptr[(i * padded_nvar) + aspect_no];
    result[i] = sqrt((nsamp * s2) - (s1 * s1)) / sqrt(nsamp * (nsamp - 1.0));
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
double ScoreCard::reportVarianceOfStates(const StateVariable aspect, const int system_index) {
#ifdef OMNI_USE_HPC
  running_accumulators.download();
  squared_accumulators.download();
#endif
  const int aspect_no = static_cast<int>(aspect);
  const int nvar = static_cast<int>(StateVariable::ALL_STATES);
  const int offset = system_index * roundUp(nvar, warp_size_int);
  const double nsamp = static_cast<double>(sampled_step_count);
  const double s1 = running_accumulators.readHost(offset + aspect_no);
  const double s2 = squared_accumulators.readHost(offset + aspect_no);
  return sqrt((nsamp * s2) - (s1 * s1)) / sqrt(nsamp * (nsamp - 1.0));
}

} // namespace energy
} // namespace omni
