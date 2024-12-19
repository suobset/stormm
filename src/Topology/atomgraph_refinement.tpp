// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace topology {

//-------------------------------------------------------------------------------------------------
template <typename T>
ParameterUnion<T>::ParameterUnion(const T* comp_xa, const T* comp_ya, const T* comp_za,
                                  const size_t aparam_count_in, const T* comp_xb, const T* comp_yb,
                                  const T* comp_zb, const size_t bparam_count_in,
                                  const double match_tol) :
    aparam_count{aparam_count_in}, bparam_count{bparam_count_in}, total_parameters{0},
    union_indices_bmap{std::vector<int>(bparam_count_in, -1)},
    union_parameter_x{},
    union_parameter_y{},
    union_parameter_z{}
{
  // Count the number of parameters in set B that can be matched to parameters in set A
  size_t nmatch = 0;
  for (size_t i = 0; i < bparam_count; i++) {
    for (size_t j = 0; j < aparam_count; j++) {
      if (fabs(comp_xb[i] - comp_xa[j]) <= match_tol &&
          fabs(comp_yb[i] - comp_ya[j]) <= match_tol &&
          (comp_za == nullptr || comp_zb == nullptr ||
           fabs(comp_zb[i] - comp_za[j]) <= match_tol)) {
        union_indices_bmap[i] = j;
        nmatch++;
        break;
      }
    }
  }
  total_parameters = aparam_count + bparam_count - nmatch;
  const T value_zero = 0.0;
  union_parameter_x.resize(total_parameters, value_zero);
  union_parameter_y.resize(total_parameters, value_zero);
  union_parameter_z.resize(total_parameters, value_zero);
  for (size_t i = 0; i < aparam_count; i++) {
    union_parameter_x[i] = comp_xa[i];
    union_parameter_y[i] = comp_ya[i];
    if (comp_za != nullptr) {
      union_parameter_z[i] = comp_za[i];
    }
  }
  size_t curr_idx = aparam_count;
  for (size_t i = 0; i < bparam_count; i++) {
    if (union_indices_bmap[i] == -1) {
      union_parameter_x[curr_idx] = comp_xb[i];
      union_parameter_y[curr_idx] = comp_yb[i];
      if (comp_zb != nullptr) {
        union_parameter_z[curr_idx] = comp_zb[i];
      }
      curr_idx++;
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
ParameterUnion<T>::ParameterUnion(const T* comp_xa, const T* comp_ya, const size_t aparam_count_in,
                                  const T* comp_xb, const T* comp_yb, const size_t bparam_count_in,
                                  const double match_tol) :
    ParameterUnion(comp_xa, comp_ya, nullptr, aparam_count_in, comp_xb, comp_yb, nullptr,
                   bparam_count_in, match_tol)
{}

//-------------------------------------------------------------------------------------------------
template <typename T>
ParameterUnion<T>::ParameterUnion(const std::vector<T> &comp_xa, const std::vector<T> &comp_ya,
                                  const std::vector<T> &comp_za, const std::vector<T> &comp_xb,
                                  const std::vector<T> &comp_yb, const std::vector<T> &comp_zb,
                                  const double match_tol) :
  ParameterUnion(comp_xa.data(), comp_ya.data(), comp_za.data(),
                 std::min(comp_xa.size(), std::min(comp_ya.size(), comp_za.size())),
                 comp_xb.data(), comp_yb.data(), comp_zb.data(),
                 std::min(comp_xb.size(), std::min(comp_yb.size(), comp_zb.size())), match_tol)
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
  ParameterUnion(comp_xa.data(), comp_ya.data(), nullptr, std::min(comp_xa.size(), comp_ya.size()),
                 comp_xb.data(), comp_yb.data(), nullptr, std::min(comp_xb.size(), comp_yb.size()),
                 match_tol)
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
template <typename T> size_t ParameterUnion<T>::getFirstSetParameterCount() const {
  return aparam_count;
}

//-------------------------------------------------------------------------------------------------
template <typename T> size_t ParameterUnion<T>::getSecondSetParameterCount() const {
  return bparam_count;
}

//-------------------------------------------------------------------------------------------------
template <typename T> size_t ParameterUnion<T>::getUniqueParameterCount() const {
  return total_parameters;
}

//-------------------------------------------------------------------------------------------------
template <typename T> std::vector<T> ParameterUnion<T>::getUnion(const int comp_idx) const {
  switch (comp_idx) {
  case 0:
    return union_parameter_x;
  case 1:
    return union_parameter_y;
  case 2:
    return union_parameter_z;
  }
  __builtin_unreachable();
}

} // namespace topology
} // namespace stormm
