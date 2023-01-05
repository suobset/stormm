// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace math {

using data_types::isScalarType;
using data_types::getStormmScalarTypeName;

//-------------------------------------------------------------------------------------------------
template <typename TBase>
void prefixSumInPlace(TBase* vdata, const size_t n_elements, const PrefixSumType style,
                      const char* caller) {
  TBase sum = 0;
  llint llsum = 0;
  switch (style) {
  case PrefixSumType::EXCLUSIVE:
    for (size_t i = 0; i < n_elements; i++) {
      const TBase tmp_sum = vdata[i];
      vdata[i] = sum;
      sum += tmp_sum;
      llsum += static_cast<llint>(tmp_sum);
    }
    break;
  case PrefixSumType::INCLUSIVE:
    for (size_t i = 0; i < n_elements; i++) {
      sum += vdata[i];
      llsum += vdata[i];
      vdata[i] = sum;
    }
    break;
  }
  if (llsum != Approx(sum, ComparisonType::RELATIVE, constants::tiny)) {
    const std::string tbase_name = isScalarType<TBase>() ?
                                   getStormmScalarTypeName<TBase>() :
                                   std::string(std::type_index(typeid(TBase)).name());
    const std::string callfunc = (caller == nullptr) ? "" : ".  Called by " + std::string(caller);
    rtErr("Overflow of numerical format.  Summation of a " + std::to_string(n_elements) +
          "-element series as llint produces " + std::to_string(llsum) +
          ", whereas the array's base type " + tbase_name + " records " + std::to_string(sum) +
          callfunc + ".", "prefixSumInPlace");
  }
}

//-------------------------------------------------------------------------------------------------
template <typename TBase>
void prefixSumInPlace(std::vector<TBase> *v, const PrefixSumType style, const char* caller) {
  prefixSumInPlace(v->data(), v->size(), style, caller);
}

//-------------------------------------------------------------------------------------------------
template <typename TBase>
void prefixSumInPlace(Hybrid<TBase> *v, const PrefixSumType style, const char* caller) {
  prefixSumInPlace(v->data(), v->size(), style, caller);
}

//-------------------------------------------------------------------------------------------------
template <typename TSum, typename TBase>
TSum sum(const TBase* v, const size_t vlen) {
  TSum total = static_cast<TSum>(0);
  for (size_t i = 0; i < vlen; i++) {
    total += static_cast<TSum>(v[i]);
  }
  return total;
}

//-------------------------------------------------------------------------------------------------
template <typename TSum, typename TBase>
TSum sum(const std::vector<TBase> &v) {
  const size_t vlen = v.size();
  TSum total = static_cast<TSum>(0);
  for (size_t i = 0; i < vlen; i++) {
    total += static_cast<TSum>(v[i]);
  }
  return total;
}

//-------------------------------------------------------------------------------------------------
template <typename TSum, typename TBase> TSum sum(const Hybrid<TBase> &hb) {
  const size_t hblen = hb.size();
  TSum total = static_cast<TSum>(0);
  const TBase* hbptr = hb.data();
  for (size_t i = 0; i < hblen; i++) {
    total += static_cast<TSum>(hbptr[i]);
  }
  return total;
}
  
//-------------------------------------------------------------------------------------------------
template <typename TSum, typename TBase>
TSum sumTuple2(const TBase* v, const size_t vlen) {
  TSum total = { 0, 0 };
  for (size_t i = 0; i < vlen; i++) {
    total.x += v[i].x;
    total.y += v[i].y;
  }
  return total;
}

//-------------------------------------------------------------------------------------------------
template <typename TSum, typename TBase>
TSum sumTuple2(const std::vector<TBase> &v) {
  const size_t vlen = v.size();  
  TSum total = { 0, 0 };
  for (size_t i = 0; i < vlen; i++) {
    total.x += v[i].x;
    total.y += v[i].y;
  }
  return total;
}

//-------------------------------------------------------------------------------------------------
template <typename TSum, typename TBase>
TSum sumTuple2(const Hybrid<TBase> &hb, const HybridTargetLevel tier) {
  const size_t hblen = hb.size();
  TSum total = { 0, 0 };
  const TBase* hbptr = hb.data(tier);
  switch (tier) {
  case HybridTargetLevel::HOST:
    for (size_t i = 0; i < hblen; i++) {
      total.x += static_cast<TSum>(hbptr[i].x);
      total.y += static_cast<TSum>(hbptr[i].y);
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    rtErr("Device summation is not yet implemented.", "sumTuple2");
#endif
  }
  return total;
}

//-------------------------------------------------------------------------------------------------
template <typename TSum, typename TBase>
TSum sumTuple3(const TBase* v, const size_t vlen) {
  TSum total = { 0, 0, 0 };
  for (size_t i = 0; i < vlen; i++) {
    total.x += v[i].x;
    total.y += v[i].y;
    total.z += v[i].z;
  }
  return total;
}

//-------------------------------------------------------------------------------------------------
template <typename TSum, typename TBase>
TSum sumTuple3(const std::vector<TBase> &v) {
  const size_t vlen = v.size();  
  TSum total = { 0, 0, 0 };
  for (size_t i = 0; i < vlen; i++) {
    total.x += v[i].x;
    total.y += v[i].y;
    total.z += v[i].z;
  }
  return total;
}

//-------------------------------------------------------------------------------------------------
template <typename TSum, typename TBase>
TSum sumTuple3(const Hybrid<TBase> &hb, const HybridTargetLevel tier) {
  const size_t hblen = hb.size();
  TSum total = { 0, 0, 0 };
  const TBase* hbptr = hb.data(tier);
  switch (tier) {
  case HybridTargetLevel::HOST:
    for (size_t i = 0; i < hblen; i++) {
      total.x += static_cast<TSum>(hbptr[i].x);
      total.y += static_cast<TSum>(hbptr[i].y);
      total.z += static_cast<TSum>(hbptr[i].z);
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    rtErr("Device summation is not yet implemented.", "sumTuple3");
#endif
  }
  return total;
}

//-------------------------------------------------------------------------------------------------
template <typename TSum, typename TBase>
TSum sumTuple4(const TBase* v, const size_t vlen) {
  TSum total = { 0, 0, 0, 0 };
  for (size_t i = 0; i < vlen; i++) {
    total.x += v[i].x;
    total.y += v[i].y;
    total.z += v[i].z;
    total.w += v[i].w;
  }
  return total;
}

//-------------------------------------------------------------------------------------------------
template <typename TSum, typename TBase>
TSum sumTuple4(const std::vector<TBase> &v) {
  const size_t vlen = v.size();  
  TSum total = { 0, 0, 0, 0 };
  for (size_t i = 0; i < vlen; i++) {
    total.x += v[i].x;
    total.y += v[i].y;
    total.z += v[i].z;
    total.w += v[i].w;
  }
  return total;
}

//-------------------------------------------------------------------------------------------------
template <typename TSum, typename TBase>
TSum sumTuple4(const Hybrid<TBase> &hb, const HybridTargetLevel tier) {
  const size_t hblen = hb.size();
  TSum total = { 0, 0, 0, 0 };
  const TBase* hbptr = hb.data(tier);
  switch (tier) {
  case HybridTargetLevel::HOST:
    for (size_t i = 0; i < hblen; i++) {
      total.x += static_cast<TSum>(hbptr[i].x);
      total.y += static_cast<TSum>(hbptr[i].y);
      total.z += static_cast<TSum>(hbptr[i].z);
      total.w += static_cast<TSum>(hbptr[i].w);
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    rtErr("Device summation is not yet implemented.", "sumTuple4");
#endif
  }
  return total;
}

} // namespace math
} // namespace stormm
