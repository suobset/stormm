// -*-c++-*-

namespace omni {
namespace math {

//-------------------------------------------------------------------------------------------------
template <typename T> T roundUp(T jagged, T increment) {
  return ((jagged + increment - static_cast<T>(1)) / increment) * increment;
}

} // namespace math
} // namespace omni
