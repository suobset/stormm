// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace topology {

//-------------------------------------------------------------------------------------------------
template <typename T>
ParameterUnion<T>::ParameterUnion(const T* comp_xa, const T* comp_ya, const T* comp_za,
                                  const T* comp_wa, const T* comp_va, const int aparam_count,
                                  const T* comp_xb, const T* comp_yb, const T* comp_zb,
                                  const T* comp_wb, const T* comp_vb, const int bparam_count,
                                  const double match_tol) :
    total_parameters{0},
    component_count{1 + (comp_ya != nullptr) + (comp_za != nullptr) + (comp_wa != nullptr) +
                    (comp_va != nullptr)},
    set_count{1}, union_parameter_x{}, union_parameter_y{}, union_parameter_z{},
    union_parameter_w{}, union_parameter_v{}, set_to_consensus_map{}
{
  // Make the union parameter set the contents of set A
  total_parameters = aparam_count;
  set_to_consensus_map.resize(1);
  set_to_consensus_map[0] = incrementingSeries(0, total_parameters);
  union_parameter_x.resize(total_parameters);
  union_parameter_y.resize(total_parameters);
  union_parameter_z.resize(total_parameters);
  union_parameter_w.resize(total_parameters);
  union_parameter_v.resize(total_parameters);
  for (int i = 0; i < total_parameters; i++) {
    union_parameter_x[i] = comp_xa[i];
    union_parameter_y[i] = (component_count >= 2) ? comp_ya[i] : 0.0;
    union_parameter_z[i] = (component_count >= 3) ? comp_za[i] : 0.0;
    union_parameter_w[i] = (component_count >= 4) ? comp_wa[i] : 0.0;
    union_parameter_v[i] = (component_count == 5) ? comp_va[i] : 0.0;
  }

  // Add the second set, if it exists
  if (comp_xb != nullptr) {
    addSet(comp_xb, comp_yb, comp_zb, comp_wb, comp_vb, bparam_count, match_tol);
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
ParameterUnion<T>::ParameterUnion(const T* comp_xa, const T* comp_ya, const int aparam_count_in,
                                  const T* comp_xb, const T* comp_yb, const int bparam_count_in,
                                  const double match_tol) :
    ParameterUnion(comp_xa, comp_ya, nullptr, nullptr, nullptr, aparam_count_in, comp_xb, comp_yb,
                   nullptr, nullptr, nullptr, bparam_count_in, match_tol)
{}

//-------------------------------------------------------------------------------------------------
template <typename T>
ParameterUnion<T>::ParameterUnion(const T* comp_xa, const int aparam_count_in, const T* comp_xb,
                                  const int bparam_count_in, const double match_tol) :
    ParameterUnion(comp_xa, nullptr, nullptr, nullptr, nullptr, aparam_count_in, comp_xb, nullptr,
                   nullptr, nullptr, nullptr, bparam_count_in, match_tol)
{}

//-------------------------------------------------------------------------------------------------
template <typename T>
ParameterUnion<T>::ParameterUnion(const T* comp_xa, const T* comp_ya, const T* comp_za,
                                  const T* comp_wa, const T* comp_va, const int aparam_count_in) :
    ParameterUnion(comp_xa, comp_ya, comp_za, comp_wa, comp_va, aparam_count_in, nullptr, nullptr,
                   nullptr, nullptr, nullptr, 0, constants::small)
{}

//-------------------------------------------------------------------------------------------------
template <typename T>
ParameterUnion<T>::ParameterUnion(const T* comp_xa, const T* comp_ya, const T* comp_za,
                                  const int aparam_count_in) :
    ParameterUnion(comp_xa, comp_ya, comp_za, nullptr, nullptr, aparam_count_in, nullptr, nullptr,
                   nullptr, nullptr, nullptr, 0, constants::small)
{}

//-------------------------------------------------------------------------------------------------
template <typename T>
ParameterUnion<T>::ParameterUnion(const T* comp_xa, const T* comp_ya, const int aparam_count_in) :
    ParameterUnion(comp_xa, comp_ya, nullptr, nullptr, nullptr, aparam_count_in, nullptr, nullptr,
                   nullptr, nullptr, nullptr, 0, constants::small)
{}

//-------------------------------------------------------------------------------------------------
template <typename T>
ParameterUnion<T>::ParameterUnion(const T* comp_xa, const int aparam_count_in) :
    ParameterUnion(comp_xa, nullptr, nullptr, nullptr, nullptr, aparam_count_in, nullptr, nullptr,
                   nullptr, nullptr, nullptr, 0, constants::small)
{}

//-------------------------------------------------------------------------------------------------
template <typename T>
ParameterUnion<T>::ParameterUnion(const std::vector<T> &comp_xa, const std::vector<T> &comp_ya,
                                  const std::vector<T> &comp_za, const std::vector<T> &comp_wa,
                                  const std::vector<T> &comp_va, const std::vector<T> &comp_xb,
                                  const std::vector<T> &comp_yb, const std::vector<T> &comp_zb,
                                  const std::vector<T> &comp_wb, const std::vector<T> &comp_vb,
                                  const double match_tol) :
    ParameterUnion(comp_xa.data(), comp_ya.data(), comp_za.data(), comp_wa.data(), comp_va.data(),
                   static_cast<int>(minValue<size_t>({ comp_xa.size(), comp_ya.size(),
                                                       comp_za.size(), comp_wa.size(),
                                                       comp_va.size() })),
                   comp_xb.data(), comp_yb.data(), comp_zb.data(), comp_wb.data(), comp_vb.data(),
                   static_cast<int>(minValue<size_t>({ comp_xb.size(), comp_yb.size(),
                                                       comp_zb.size(), comp_wb.size(),
                                                       comp_vb.size() })),
                   match_tol)
{
  // Check that the parameter components were of the same size.  It is an error otherwise.
  if (comp_xa.size() != comp_ya.size() || comp_xa.size() != comp_za.size() ||
      comp_xa.size() != comp_wa.size() || comp_xa.size() != comp_va.size()) {
    rtErr("Different numbers of components were presented for the first parameter set (" +
          std::to_string(comp_xa.size()) + ", " + std::to_string(comp_ya.size()) + ", " +
          std::to_string(comp_za.size()) + ", " + std::to_string(comp_wa.size()) + ", and " +
          std::to_string(comp_va.size()) + ").", "ParameterUnion");
  }
  if (comp_xb.size() != comp_yb.size() || comp_xb.size() != comp_zb.size() ||
      comp_xb.size() != comp_wb.size() || comp_xb.size() != comp_vb.size()) {
    rtErr("Different numbers of components were presented for the second parameter set (" +
          std::to_string(comp_xb.size()) + ", " + std::to_string(comp_yb.size()) + ", " +
          std::to_string(comp_zb.size()) + ", " + std::to_string(comp_wb.size()) + ", and " +
          std::to_string(comp_vb.size()) + ").", "ParameterUnion");
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
ParameterUnion<T>::ParameterUnion(const std::vector<T> &comp_xa, const std::vector<T> &comp_ya,
                                  const std::vector<T> &comp_za, const std::vector<T> &comp_xb,
                                  const std::vector<T> &comp_yb, const std::vector<T> &comp_zb,
                                  const double match_tol) :
    ParameterUnion(comp_xa.data(), comp_ya.data(), comp_za.data(), nullptr, nullptr,
                   static_cast<int>(std::min(comp_xa.size(),
                                             std::min(comp_ya.size(), comp_za.size()))),
                   comp_xb.data(), comp_yb.data(), comp_zb.data(), nullptr, nullptr,
                   static_cast<int>(std::min(comp_xb.size(),
                                             std::min(comp_yb.size(), comp_zb.size()))),
                   match_tol)
{
  // Check that the parameter components were of the same size.  It is an error otherwise.
  if (comp_xa.size() != comp_ya.size() || comp_xa.size() != comp_za.size()) {
    rtErr("Different numbers of components were presented for the first parameter set (" +
          std::to_string(comp_xa.size()) + ", " + std::to_string(comp_ya.size()) + ", and " +
          std::to_string(comp_za.size()) + ").", "ParameterUnion");
  }
  if (comp_xb.size() != comp_yb.size() || comp_xb.size() != comp_zb.size()) {
    rtErr("Different numbers of components were presented for the second parameter set (" +
          std::to_string(comp_xb.size()) + ", " + std::to_string(comp_yb.size()) + ", and " +
          std::to_string(comp_zb.size()) + ").", "ParameterUnion");
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
ParameterUnion<T>::ParameterUnion(const std::vector<T> &comp_xa, const std::vector<T> &comp_ya,
                                  const std::vector<T> &comp_xb, const std::vector<T> &comp_yb,
                                  const double match_tol) :
    ParameterUnion(comp_xa.data(), comp_ya.data(), nullptr, nullptr, nullptr,
                   std::min(comp_xa.size(), comp_ya.size()), comp_xb.data(), comp_yb.data(),
                   nullptr, nullptr, nullptr, std::min(comp_xb.size(), comp_yb.size()), match_tol)
{
  // Check that the parameter components were of the same size.  It is an error otherwise.
  if (comp_xa.size() != comp_ya.size()) {
    rtErr("Different numbers of components were presented for the first parameter set (" +
          std::to_string(comp_xa.size()) + " and " + std::to_string(comp_ya.size()) + ").",
          "ParameterUnion");
  }
  if (comp_xb.size() != comp_yb.size()) {
    rtErr("Different numbers of components were presented for the second parameter set (" +
          std::to_string(comp_xb.size()) + " and " + std::to_string(comp_yb.size()) + ").",
          "ParameterUnion");
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
ParameterUnion<T>::ParameterUnion(const std::vector<T> &comp_xa, const std::vector<T> &comp_xb,
                                  const double match_tol) :
    ParameterUnion(comp_xa.data(), nullptr, nullptr, nullptr, nullptr, comp_xa.size(),
                   comp_xb.data(), nullptr, nullptr, nullptr, nullptr, comp_xb.size(), match_tol)
{}

//-------------------------------------------------------------------------------------------------
template <typename T>
ParameterUnion<T>::ParameterUnion(const std::vector<T> &comp_xa, const std::vector<T> &comp_ya,
                                  const std::vector<T> &comp_za) :
    ParameterUnion(comp_xa.data(), comp_ya.data(), comp_za.data(),
                   std::min(comp_xa.size(), std::min(comp_ya.size(), comp_za.size())),
                   nullptr, nullptr, nullptr, 0, constants::small)
{}

//-------------------------------------------------------------------------------------------------
template <typename T>
ParameterUnion<T>::ParameterUnion(const std::vector<T> &comp_xa, const std::vector<T> &comp_ya) :
    ParameterUnion(comp_xa.data(), comp_ya.data(), nullptr,
                   std::min(comp_xa.size(), comp_ya.size()),
                   nullptr, nullptr, nullptr, 0, constants::small)
{}

//-------------------------------------------------------------------------------------------------
template <typename T>
ParameterUnion<T>::ParameterUnion(const std::vector<T> &comp_xa) :
    ParameterUnion(comp_xa.data(), nullptr, nullptr, comp_xa.size(), nullptr, nullptr, nullptr, 0,
                   constants::small)
{}

//-------------------------------------------------------------------------------------------------
template <typename T> int ParameterUnion<T>::getInputSetParameterCount(const int set_index) const {
  validateSetIndex(set_index, "getInputSetParameterCount");
  return set_to_consensus_map[set_index].size();
}

//-------------------------------------------------------------------------------------------------
template <typename T> int ParameterUnion<T>::getUniqueParameterCount() const {
  return total_parameters;
}

//-------------------------------------------------------------------------------------------------
template <typename T> std::vector<T> ParameterUnion<T>::getUnion(const int comp_idx) const {
  if (comp_idx >= component_count) {
    rtErr("The object contains " + std::to_string(component_count) + " components per parameter.  "
          "Component index " + std::to_string(comp_idx) + " is invalid.", "ParameterUnion",
          "getUnion");
  }
  switch (comp_idx) {
  case 0:
    return union_parameter_x;
  case 1:
    return union_parameter_y;
  case 2:
    return union_parameter_z;
  case 3:
    return union_parameter_w;
  case 4:
    return union_parameter_v;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename T>
const std::vector<int>& ParameterUnion<T>::getCorrespondence(const int set_index) const {
  validateSetIndex(set_index, "getCorrespondence");
  return set_to_consensus_map[set_index];
}

//-------------------------------------------------------------------------------------------------
template <typename T>
int ParameterUnion<T>::getCorrespondence(const int set_index, const int param_index) const {
  validateSetIndex(set_index);
  if (param_index < 0 || param_index >= set_to_consensus_map[set_index].size()) {
    rtErr("Parameter index " + std::to_string(param_index) + " is invalid for an input set of " +
          std::to_string(set_to_consensus_map[set_index].size()) + " input parameters.",
          "ParameterUnion", "getCorrespondence");
  }
  return set_to_consensus_map[set_index][param_index];
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void ParameterUnion<T>::addSet(const T* comp_x, const T* comp_y, const T* comp_z, const T* comp_w,
                               const T* comp_v, const int nparm, const double tol) {

  // Check that the set to be added comports with the existing parameter layout
  if (component_count == 2 && comp_z != nullptr) {
    rtErr("The object's layout is for two parameter components, but a third was provided.",
          "ParameterUnion", "addSet");
  }
  if (component_count == 3 && comp_z == nullptr) {
    rtErr("The object's layout is for three parameter components, but only two were provided.",
          "ParameterUnion", "addSet");
  }
  
  // Count the number of parameters in the new set that can be matched to parameters in union
  int nmatch = 0;
  std::vector<int> correspondence(nparm, -1);
  for (int i = 0; i < nparm; i++) {
    for (int j = 0; j < total_parameters; j++) {
      if (fabs(comp_x[i] - union_parameter_x[j]) <= tol &&
          (component_count < 2 || (fabs(comp_y[i] - union_parameter_y[j]) <= tol)) &&
          (component_count < 3 || (fabs(comp_z[i] - union_parameter_z[j]) <= tol)) &&
          (component_count < 4 || (fabs(comp_w[i] - union_parameter_w[j]) <= tol)) &&
          (component_count < 5 || (fabs(comp_v[i] - union_parameter_v[j]) <= tol))) {
        correspondence[i] = j;
        nmatch++;
        break;
      }
    }
  }
  int next_idx = total_parameters;
  total_parameters += nparm - nmatch;
  const T value_zero = 0.0;
  union_parameter_x.resize(total_parameters, value_zero);
  union_parameter_y.resize(total_parameters, value_zero);
  union_parameter_z.resize(total_parameters, value_zero);
  union_parameter_w.resize(total_parameters, value_zero);
  union_parameter_v.resize(total_parameters, value_zero);
  for (int i = 0; i < nparm; i++) {
    if (correspondence[i] == -1) {
      union_parameter_x[next_idx] = comp_x[i];
      if (component_count >= 2) {
        union_parameter_y[next_idx] = comp_y[i];
      }
      if (component_count >= 3) {
        union_parameter_z[next_idx] = comp_z[i];
      }
      if (component_count >= 4) {
        union_parameter_w[next_idx] = comp_w[i];
      }
      if (component_count >= 5) {
        union_parameter_v[next_idx] = comp_v[i];
      }
      correspondence[i] = next_idx;
      next_idx++;
    }
  }
  set_count += 1;
  set_to_consensus_map.push_back(correspondence);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void ParameterUnion<T>::addSet(const T* comp_x, const T* comp_y, const T* comp_z, const int nparm,
                               const double tol) {
  addSet(comp_x, comp_y, comp_z, nullptr, nullptr, nparm, tol);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void ParameterUnion<T>::addSet(const T* comp_x, const T* comp_y, const int nparm,
                               const double tol) {
  addSet(comp_x, comp_y, nullptr, nullptr, nullptr, nparm, tol);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void ParameterUnion<T>::addSet(const T* comp_x, const int nparm, const double tol) {
  addSet(comp_x, nullptr, nullptr, nullptr, nullptr, nparm, tol);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void ParameterUnion<T>::addSet(const std::vector<T> &comp_x, const std::vector<T> &comp_y,
                               const std::vector<T> &comp_z, const std::vector<T> &comp_w,
                               const std::vector<T> &comp_v, const double tol) {
  addSet(comp_x.data(), comp_y.data(), comp_z.data(), comp_w.data(), comp_v.data(),
         static_cast<int>(minValue<size_t>({ comp_x.size(), comp_y.size(), comp_z.size(),
                                             comp_w.size(), comp_v.size() })), tol);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void ParameterUnion<T>::addSet(const std::vector<T> &comp_x, const std::vector<T> &comp_y,
                               const std::vector<T> &comp_z, const double tol) {
  addSet(comp_x.data(), comp_y.data(), comp_z.data(), nullptr, nullptr,
         static_cast<int>(minValue<size_t>({ comp_x.size(), comp_y.size(), comp_z.size() })), tol);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void ParameterUnion<T>::addSet(const std::vector<T> &comp_x, const std::vector<T> &comp_y,
                               const double tol) {
  addSet(comp_x.data(), comp_y.data(), nullptr, nullptr, nullptr,
         static_cast<int>(std::min(comp_x.size(), comp_y.size())), tol);
}

//-------------------------------------------------------------------------------------------------
template <typename T> void ParameterUnion<T>::validateSetIndex(const int set_index,
                                                               const char* caller) const {
  if (set_index < 0 || set_index >= set_count) {
    rtErr("Set index " + std::to_string(set_index) + " is invalid for a collection of " +
          std::to_string(set_count) + " input sets.", "ParameterUnion", caller);
  }
}

} // namespace topology
} // namespace stormm
