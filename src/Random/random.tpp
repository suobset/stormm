// -*-c++-*-
namespace omni {
namespace random {

//-------------------------------------------------------------------------------------------------
template <typename T>
Xoroshiro128pSeries<T>::Xoroshiro128pSeries(const size_t generators_in, const size_t depth_in,
                                            const RandomNumberKind init_kind, const int igseed_in,
                                            const int niter, const size_t bank_limit) :
  generators{generators_in}, depth{depth_in}, igseed{igseed_in},
  states{generators, "xrs128_states"},
  bank{HybridKind::ARRAY, "xrs128_states"}
{
  // Check that the depth is a multiple of two if the storage type is float.
  const size_t tcalc_ct = std::type_index(typeid(T)).hash_code();
  
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
  Xoroshiro128pGenerator xrs(igseed, niter);
  states.putHost(xrs.revealState(), 0);
  for (size_t i = 1LLU; i < generators; i++) {
    xrs.longJump();
    states.putHost(xrs.revealState(), i);
  }

  // Initialize the random number bank based on the states of all generators.
  for (size_t i = 0LLU; i < generators; i++) {
    xrs.setState(states.readHost(i));
    for (size_t j = 0LLU; j < depth; j++) {
      xrs.uniformRandomNumber();
    }
  }
}

} // namespace random
} // namespace omni
