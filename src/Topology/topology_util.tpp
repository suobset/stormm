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
  result->resize(va_bounds[index + 1] - llim);
  T* res_ptr = result->data();
  for (int i = llim; i < va_bounds[index + 1]; i++) {
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
  std::vector<T> result(va_bounds[index + 1] - llim);
  for (int i = llim; i < va_bounds[index + 1]; i++) {
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

} // namespace topology
} // namespace omni
