// -*-c++-*-
namespace omni {
namespace topology {

//-------------------------------------------------------------------------------------------------
template <typename T>
void extractBoundedListEntries(std::vector<T> *result, const std::vector<T> &va,
                               const std::vector<int> &va_bounds, const int index) {
  if (index >= static_cast<int>(va_bounds.size()) + 1) {
    rtErr("Index " + std::to_string(index) + " is invalid for a bounds array with " +
          std::to_string(va_bounds.size() - 1LLU) + " elements.", "extractBoundedListEntries");
  }
  const int llim = va_bounds[index];
  const int hlim = va_bounds[index + 1];
  result->resize(hlim - llim);
  T* res_ptr = result->data();
  for (int i = llim; i < hlim; i++) {
    res_ptr[i - llim] = va[i];
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T> std::vector<T> extractBoundedListEntries(const std::vector<T> &va,
                                                               const std::vector<int> &va_bounds,
                                                               const int index) {
  if (index >= static_cast<int>(va_bounds.size()) + 1) {
    rtErr("Index " + std::to_string(index) + " is invalid for a bounds array with " +
          std::to_string(va_bounds.size() - 1LLU) + " elements.", "extractBoundedListEntries");
  }
  const int llim = va_bounds[index];
  const int hlim = va_bounds[index + 1];
  std::vector<T> result(hlim - llim);
  for (int i = llim; i < hlim; i++) {
    result[i - llim] = va[i];
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename T> std::vector<T> extractBoundedListEntries(const std::vector<T> &va,
                                                               const std::vector<int> &va_bounds,
                                                               const std::vector<int> &indices,
                                                               const UniqueValueHandling filter) {
  const int nidx = indices.size();
  int buffer_size = 0LLU;
  for (int pos = 0; pos < nidx; pos++) {
    const int index = indices[pos];
    if (index >= static_cast<int>(va_bounds.size()) + 1) {
      rtErr("Index " + std::to_string(index) + " is invalid for a bounds array with " +
            std::to_string(va_bounds.size() - 1LLU) + " elements.", "extractBoundedListEntries");
    }
    buffer_size += va_bounds[index + 1] - va_bounds[index];
  }
  std::vector<T> result(buffer_size);
  int j = 0;
  for (int pos = 0; pos < nidx; pos++) {
    const int index = indices[pos];
    for (int i = va_bounds[index]; i < va_bounds[index + 1]; i++) {
      result[j] = va[i];
      j++;
    }
  }
  switch (filter) {
  case UniqueValueHandling::UNIQUE_VALUES_ONLY:
    reduceUniqueValues(&result);
    break;
  case UniqueValueHandling::CONFIRM_ALL_COPIES:
    break;
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void extractBoundedListEntries(std::vector<T> *result, const std::vector<T> &va,
                               const std::vector<int> &va_bounds, const std::vector<int> &indices,
                               const UniqueValueHandling filter) {
  const int nidx = indices.size();
  int buffer_size = 0LLU;
  for (int pos = 0; pos < nidx; pos++) {
    const int index = indices[pos];
    if (index >= static_cast<int>(va_bounds.size()) + 1) {
      rtErr("Index " + std::to_string(index) + " is invalid for a bounds array with " +
            std::to_string(va_bounds.size() - 1LLU) + " elements.", "extractBoundedListEntries");
    }
    buffer_size += va_bounds[index + 1] - va_bounds[index];
  }
  result->resize(buffer_size);
  T* res_ptr = result->data();
  int j = 0;
  for (int pos = 0; pos < nidx; pos++) {
    const int index = indices[pos];
    for (int i = va_bounds[index]; i < va_bounds[index + 1]; i++) {
      res_ptr[j] = va[i];
      j++;
    }
  }
  switch (filter) {
  case UniqueValueHandling::UNIQUE_VALUES_ONLY:
    reduceUniqueValues(result);
    break;
  case UniqueValueHandling::CONFIRM_ALL_COPIES:
    break;
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
std::vector<T> getRealParameters(const Hybrid<double> &item, const Hybrid<float> &sp_item,
                                 const HybridTargetLevel tier, const int low_index,
                                 const int high_index, const char* caller, const char* method) {
  const size_t ct = std::type_index(typeid(T)).hash_code();
  if (ct == float_type_index) {
    std::vector<float> tmpv;
    switch (tier) {
    case HybridTargetLevel::HOST:
      tmpv = sp_item.readHost(low_index, high_index - low_index);
      break;
    case HybridTargetLevel::DEVICE:
      tmpv = sp_item.readDevice(low_index, high_index - low_index);
      break;
    }
    return std::vector<T>(tmpv.begin(), tmpv.end());
  }
  else if (ct == double_type_index) {
    std::vector<double> tmpv;
    switch (tier) {
    case HybridTargetLevel::HOST:
      tmpv = item.readHost(low_index, high_index - low_index);
      break;
    case HybridTargetLevel::DEVICE:
      tmpv = item.readDevice(low_index, high_index - low_index);
      break;
    }
    return std::vector<T>(tmpv.begin(), tmpv.end());
  }
  else {
    if (isScalarType<T>()) {
      rtErr("Data is not available in type " + getOmniScalarTypeName<T>() + ".", caller, method);
    }
    else if (isHpcVectorType<T>()) {
      rtErr("Data is not available in type " + getOmniHpcVectorTypeName<T>() + ".", caller,
            method);
    }
    else {
      rtErr("Unrecognized data type " + std::string(std::type_index(typeid(T)).name()) + ".",
            caller, method);
    }
  }
  __builtin_unreachable();
}

} // namespace topology
} // namespace omni
