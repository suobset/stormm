// -*-c++-*-
namespace omni {
namespace math {

using data_types::isScalarType;
using data_types::getOmniScalarTypeName;

//-------------------------------------------------------------------------------------------------
template <typename TSum, typename TBase>
void prefixSumInPlace(std::vector<TBase> &v, const PrefixSumType style, const char* caller) {
  const size_t n_elements = v.size();
  TBase sum = 0;
  TSum llsum = 0;
  switch (style) {
  case PrefixSumType::EXCLUSIVE:
    for (size_t i = 0; i < n_elements; i++) {
      const TBase tmp_sum = v[i];
      v[i] = sum;
      sum += tmp_sum;
      llsum += static_cast<llint>(tmp_sum);
    }
    break;
  case PrefixSumType::INCLUSIVE:
    for (size_t i = 0; i < n_elements; i++) {
      sum += v[i];
      llsum += v[i];
      v[i] = sum;
    }
    break;
  }
  if (llsum != Approx(sum, ComparisonType::RELATIVE, constants::tiny)) {
    const std::string tsum_name = isScalarType<TSum>() ?
                                  getOmniScalarTypeName<TSum>() :
                                  std::string(std::type_index(typeid(TSum)).name());
    const std::string tbase_name = isScalarType<TBase>() ?
                                   getOmniScalarTypeName<TBase>() :
                                   std::string(std::type_index(typeid(TBase)).name());
    const std::string callfunc = (caller == nullptr) ? "" : ".  Called by " + std::string(caller);
    rtErr("Overflow of numerical format.  Summation of a " + std::to_string(n_elements) +
          "-element series as " + tsum_name + " produces " + std::to_string(llsum) +
          ", whereas the vector's base type " + tbase_name + " records " + std::to_string(sum) +
          callfunc + ".", "prefixSumInPlace");
  }
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
#ifdef OMNI_USE_HPC
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
#ifdef OMNI_USE_HPC
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
#ifdef OMNI_USE_HPC
  case HybridTargetLevel::DEVICE:
    rtErr("Device summation is not yet implemented.", "sumTuple4");
#endif
  }
  return total;
}

} // namespace math
} // namespace omni
