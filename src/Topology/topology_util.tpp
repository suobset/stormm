// -*-c++-*-
namespace omni {
namespace topology {

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

} // namespace topology
} // namespace omni
