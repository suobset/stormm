// -*-c++-*-
namespace omni {
namespace random {

//-------------------------------------------------------------------------------------------------
template <typename T>
Xoroshiro128pSeries<T>::Xoroshiro128pSeries(const ullint2 state_in, const size_t generators_in,
                                            const size_t depth_in,
                                            const RandomNumberKind init_kind,
                                            const size_t bank_limit) :
    generators{generators_in}, depth{depth_in},
    refresh_stride{(generators + depth - 1LLU) / depth},
    states{generators, "xrs128_states"},
    bank{HybridKind::ARRAY, "xrs128_states"}
{
  // Check that the data type is scalar
  const size_t ct = std::type_index(typeid(T)).hash_code();
  if (isFloatingPointScalarType<T>() == false) {
    rtErr("Random number generator series can only store results in floating point scalar data "
          "types.\n", "Xoroshiro128pSeries");
  }
  const bool t_is_double = (ct == double_type_index);
  
  // Check the bank's dimensions before allocating it, as it is a product of two numbers
  if (generators * depth == 0LLU) {
    rtErr("A random number generator series cannot be constructed with " +
          std::to_string(generators) + " generators and " + std::to_string(depth) + " depth.",
          "Xoroshiro128pSeries");
  }
  if (generators * depth > bank_limit) {
    rtErr("It is not permitted to allocate " + std::to_string(generators * depth) + " random "
          "numbers.  A maximum of " + std::to_string(bank_limit) + " numbers may be stored within "
          "the machine's stated limits.  If no keyword is available to refine this quantity, it "
          "may be necessary to recompile the code with a more permissive built-in limit.",
          "Xoroshiro128pSeries");
  }
  bank.resize(generators * depth);

  // Seed the first random number generator, and from there all the rest.
  Xoroshiro128pGenerator xrs;
  xrs.setState(state_in);
  for (size_t i = 0LLU; i < generators; i++) {
    states.putHost(xrs.revealState(), i);
    xrs.longJump();
  }

  // Initialize the random number bank based on the states of all generators.
  T* bank_ptr = bank.data();
  for (size_t i = 0LLU; i < generators; i++) {
    xrs.setState(states.readHost(i));
    switch (init_kind) {
    case RandomNumberKind::UNIFORM:
      if (t_is_double) {
        for (size_t j = 0LLU; j < depth; j++) {
          bank_ptr[(j * generators) + i] = xrs.uniformRandomNumber();
        }
      }
      else {
        for (size_t j = 0LLU; j < depth; j++) {
          bank_ptr[(j * generators) + i] = xrs.spUniformRandomNumber();
        }
      }
      break;
    case RandomNumberKind::GAUSSIAN:
      if (t_is_double) {
        for (size_t j = 0LLU; j < depth; j++) {
          bank_ptr[(j * generators) + i] = xrs.gaussianRandomNumber();
        }
      }
      else {
        for (size_t j = 0LLU; j < depth; j++) {
          bank_ptr[(j * generators) + i] = xrs.spGaussianRandomNumber();
        }
      }
      break;
    }
    states.putHost(xrs.revealState(), i);
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
Xoroshiro128pSeries<T>::Xoroshiro128pSeries(const size_t generators_in, const size_t depth_in,
                                            const RandomNumberKind init_kind, const int igseed_in,
                                            const int niter, const size_t bank_limit) :
  Xoroshiro128pSeries(Xoroshiro128pGenerator(igseed_in, niter).revealState, generators_in,
                      depth_in, init_kind, bank_limit)
{}

//-------------------------------------------------------------------------------------------------
template <typename T> size_t Xoroshiro128pSeries<T>::getGeneratorCount() const {
  return generators;
}

//-------------------------------------------------------------------------------------------------
template <typename T> size_t Xoroshiro128pSeries<T>::getDepth() const {
  return depth;
}

//-------------------------------------------------------------------------------------------------
template <typename T> size_t Xoroshiro128pSeries<T>::getRefreshStride() const {
  return refresh_stride;
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void Xoroshiro128pSeries<T>::uniformRandomNumbers(const size_t first_gen, const size_t last_gen) {
  Xoroshiro128pGenerator xrs;
  const size_t ct = std::type_index(typeid(T)).hash_code();
  const bool t_is_double = (ct == double_type_index);
  T* bank_ptr = bank.data();
  for (size_t i = first_gen; i < last_gen; i++) {
    xrs.setState(states.readHost(i));
    if (t_is_double) {
      for (size_t j = 0LLU; j < depth; j++) {
        bank_ptr[(j * generators) + i] = xrs.uniformRandomNumber();
      }
    }
    else {
      for (size_t j = 0LLU; j < depth; j++) {
        bank_ptr[(j * generators) + i] = xrs.spUniformRandomNumber();
      }
    }
    states.putHost(xrs.revealState(), i);
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void Xoroshiro128pSeries<T>::gaussianRandomNumbers(const size_t first_gen, const size_t last_gen) {
  Xoroshiro128pGenerator xrs;
  const size_t ct = std::type_index(typeid(T)).hash_code();
  const bool t_is_double = (ct == double_type_index);
  T* bank_ptr = bank.data();
  for (size_t i = first_gen; i < last_gen; i++) {
    xrs.setState(states.readHost(i));
    if (t_is_double) {
      for (size_t j = 0LLU; j < depth; j++) {
        bank_ptr[(j * generators) + i] = xrs.gaussianRandomNumber();
      }
    }
    else {
      for (size_t j = 0LLU; j < depth; j++) {
        bank_ptr[(j * generators) + i] = xrs.spGaussianRandomNumber();
      }
    }
    states.putHost(xrs.revealState(), i);
  }
}

} // namespace random
} // namespace omni
