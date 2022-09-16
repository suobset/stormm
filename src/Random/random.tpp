// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace random {

//-------------------------------------------------------------------------------------------------
template <typename T>
RandomNumberMill<T>::RandomNumberMill(const ullint2 state_in, const size_t generators_in,
                                      const size_t depth_in, const RandomNumberKind init_kind,
                                      const size_t bank_limit) :
    style{RandomAlgorithm::XOROSHIRO_128P}, generators{generators_in}, depth{depth_in},
    refresh_stride{(generators + depth - 1LLU) / depth},
    state_xy{generators, "state_xy"},
    state_zw{0, "state_zw"},
    bank{HybridKind::ARRAY, "rng_bank"}
{
  checkDimensions(bank_limit);
  bank.resize(generators * depth);

  // Seed the first random number generator, and from there all the rest.
  Xoroshiro128pGenerator xrs;
  xrs.setState(state_in);
  for (size_t i = 0LLU; i < generators; i++) {
    state_xy.putHost(xrs.revealState(), i);
    xrs.longJump();
  }

  // Initialize the random number bank based on the states of all generators.
  initializeBank(init_kind);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
RandomNumberMill<T>::RandomNumberMill(const ullint4 state_in, const size_t generators_in,
                                      const size_t depth_in, const RandomNumberKind init_kind,
                                      const size_t bank_limit) :
    style{RandomAlgorithm::XOSHIRO_256PP}, generators{generators_in}, depth{depth_in},
    refresh_stride{(generators + depth - 1LLU) / depth},
    state_xy{generators, "state_xy"},
    state_zw{generators, "state_zw"},
    bank{HybridKind::ARRAY, "rng_bank"}
{
  checkDimensions(bank_limit);
  bank.resize(generators * depth);

  // Seed the first random number generator, and from there all the rest.
  Xoshiro256ppGenerator xrs;
  xrs.setState(state_in);
  for (size_t i = 0LLU; i < generators; i++) {
    const ullint4 curr_state = xrs.revealState();
    const ullint2 xy_part = { curr_state.x, curr_state.y };
    const ullint2 zw_part = { curr_state.z, curr_state.w };
    state_xy.putHost(xy_part, i);
    state_zw.putHost(zw_part, i);
    xrs.longJump();
  }

  // Initialize the random number bank based on the states of all generators.
  initializeBank(init_kind);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
RandomNumberMill<T>::RandomNumberMill(const size_t generators_in, const size_t depth_in,
                                      const RandomAlgorithm style_in,
                                      const RandomNumberKind init_kind, const int igseed_in,
                                      const int niter, const size_t bank_limit) :
    style{style_in}, generators{generators_in}, depth{depth_in},
    refresh_stride{(generators + depth - 1LLU) / depth},
    state_xy{generators, "state_xy"},
    state_zw{generators, "state_zw"},
    bank{HybridKind::ARRAY, "rng_bank"}
{
  checkDimensions(bank_limit);
  bank.resize(generators * depth);

  // Seed the first random number generator, and from there all the rest.
  switch (style) {
  case RandomAlgorithm::XOROSHIRO_128P:
    {
      Xoroshiro128pGenerator xrs(igseed_in, niter);
      for (size_t i = 0LLU; i < generators; i++) {
        state_xy.putHost(xrs.revealState(), i);
        xrs.longJump();
      }
    }
    break;
  case RandomAlgorithm::XOSHIRO_256PP:
    {
      Xoshiro256ppGenerator xrs(igseed_in, niter);
      for (size_t i = 0LLU; i < generators; i++) {
        const ullint4 curr_state = xrs.revealState();
        const ullint2 xy_part = { curr_state.x, curr_state.y };
        const ullint2 zw_part = { curr_state.z, curr_state.w };
        state_xy.putHost(xy_part, i);
        state_zw.putHost(zw_part, i);
        xrs.longJump();
      }
    }
    break;
  }

  // Initialize the random number bank based on the states of all generators.
  initializeBank(init_kind);

}

//-------------------------------------------------------------------------------------------------
template <typename T> size_t RandomNumberMill<T>::getGeneratorCount() const {
  return generators;
}

//-------------------------------------------------------------------------------------------------
template <typename T> size_t RandomNumberMill<T>::getDepth() const {
  return depth;
}

//-------------------------------------------------------------------------------------------------
template <typename T> size_t RandomNumberMill<T>::getRefreshStride() const {
  return refresh_stride;
}

//-------------------------------------------------------------------------------------------------
template <typename T> T RandomNumberMill<T>::getBankValue(const size_t generator_index,
                                                          const size_t layer_index) const {
  if (generator_index >= generators || layer_index >= depth) {
    rtErr("Value index (" + std::to_string(generator_index) + ", " + std::to_string(layer_index) +
          ") is not valid in a series with " + std::to_string(generators) + " pseudo-random "
          "number generators and " + std::to_string(depth) + " bank depth.", "RandomNumberMill",
          "getBankValue");
  }
  return bank.readHost((layer_index * generators) + generator_index);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void RandomNumberMill<T>::uniformRandomNumbers(const size_t first_gen, const size_t last_gen) {
  Xoroshiro128pGenerator xrs;
  const size_t ct = std::type_index(typeid(T)).hash_code();
  const bool t_is_double = (ct == double_type_index);
  T* bank_ptr = bank.data();
  for (size_t i = first_gen; i < last_gen; i++) {
    xrs.setState(state_xy.readHost(i));
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
    state_xy.putHost(xrs.revealState(), i);
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void RandomNumberMill<T>::gaussianRandomNumbers(const size_t first_gen, const size_t last_gen) {
  Xoroshiro128pGenerator xrs;
  const size_t ct = std::type_index(typeid(T)).hash_code();
  const bool t_is_double = (ct == double_type_index);
  T* bank_ptr = bank.data();
  for (size_t i = first_gen; i < last_gen; i++) {
    xrs.setState(state_xy.readHost(i));
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
    state_xy.putHost(xrs.revealState(), i);
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void RandomNumberMill<T>::checkDimensions(const size_t bank_limit) {

  // Check that the data type is scalar
  if (isFloatingPointScalarType<T>() == false) {
    rtErr("Random number generator series can only store results in floating point scalar data "
          "types.\n", "Xoroshiro128pSeries");
  }

  // Check the bank's dimensions before allocating it, as it is a product of two numbers
  if (generators * depth == 0LLU) {
    rtErr("A random number generator series cannot be constructed with " +
          std::to_string(generators) + " generators and " + std::to_string(depth) + " depth.",
          "Xoroshiro128pSeries");
  }
  if (generators * depth * sizeof(T) > bank_limit) {
    rtErr("It is not permitted to allocate " + std::to_string(generators * depth) + " random "
          "numbers.  A maximum of " + std::to_string(bank_limit) + " numbers may be stored within "
          "the machine's stated limits.  If no keyword is available to refine this quantity, it "
          "may be necessary to recompile the code with a more permissive built-in limit.",
          "Xoroshiro128pSeries");
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void RandomNumberMill<T>::initializeBank(const RandomNumberKind init_kind) {
  T* bank_ptr = bank.data();
  const size_t ct = std::type_index(typeid(T)).hash_code();
  const bool t_is_double = (ct == double_type_index);
  Xoroshiro128pGenerator xrs128p;
  Xoshiro256ppGenerator xrs256pp;
  for (size_t i = 0LLU; i < generators; i++) {
    switch (style) {
    case RandomAlgorithm::XOROSHIRO_128P:
      {
        xrs128p.setState(state_xy.readHost(i));
        switch (init_kind) {
        case RandomNumberKind::UNIFORM:
          if (t_is_double) {
            for (size_t j = 0LLU; j < depth; j++) {
              bank_ptr[(j * generators) + i] = xrs128p.uniformRandomNumber();
            }
          }
          else {
            for (size_t j = 0LLU; j < depth; j++) {
              bank_ptr[(j * generators) + i] = xrs128p.spUniformRandomNumber();
            }
          }
          break;
        case RandomNumberKind::GAUSSIAN:
          if (t_is_double) {
            for (size_t j = 0LLU; j < depth; j++) {
              bank_ptr[(j * generators) + i] = xrs128p.gaussianRandomNumber();
            }
          }
          else {
            for (size_t j = 0LLU; j < depth; j++) {
              bank_ptr[(j * generators) + i] = xrs128p.spGaussianRandomNumber();
            }
          }
          break;
        }
        state_xy.putHost(xrs128p.revealState(), i);
      }
      break;
    case RandomAlgorithm::XOSHIRO_256PP:
      {
        ullint4 init_state;
        ullint2 tmp_state = state_xy.readHost(i);
        init_state.x = tmp_state.x;
        init_state.y = tmp_state.y;
        tmp_state = state_zw.readHost(i);
        init_state.z = tmp_state.x;
        init_state.w = tmp_state.y;
        xrs256pp.setState(init_state);
        switch (init_kind) {
        case RandomNumberKind::UNIFORM:
          if (t_is_double) {
            for (size_t j = 0LLU; j < depth; j++) {
              bank_ptr[(j * generators) + i] = xrs256pp.uniformRandomNumber();
            }
          }
          else {
            for (size_t j = 0LLU; j < depth; j++) {
              bank_ptr[(j * generators) + i] = xrs256pp.spUniformRandomNumber();
            }
          }
          break;
        case RandomNumberKind::GAUSSIAN:
          if (t_is_double) {
            for (size_t j = 0LLU; j < depth; j++) {
              bank_ptr[(j * generators) + i] = xrs256pp.gaussianRandomNumber();
            }
          }
          else {
            for (size_t j = 0LLU; j < depth; j++) {
              bank_ptr[(j * generators) + i] = xrs256pp.spGaussianRandomNumber();
            }
          }
          break;
        }
        init_state = xrs256pp.revealState();
        tmp_state.x = init_state.x;
        tmp_state.y = init_state.y;
        state_xy.putHost(tmp_state, i);
        tmp_state.x = init_state.z;
        tmp_state.y = init_state.w;
        state_zw.putHost(tmp_state, i);
      }
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void fillRandomCache(ullint2* state_xy, ullint2* state_zw, T* cache, const size_t length,
                     const size_t depth, const RandomAlgorithm method,
                     const RandomNumberKind product, const size_t index_start,
                     const size_t index_end) {
  const size_t padded_atom_count = roundUp(length, warp_size_zu);
  Xoroshiro128pGenerator xrs128p;
  Xoshiro256ppGenerator xrs256pp;
  const size_t ct = std::type_index(typeid(T)).hash_code();
  const bool t_is_double = (ct == double_type_index);
  for (int i = index_start; i < index_end; i++) {

    // Set the state of the generator object to the current value in the array of generators.
    // Next, generate a series of random numbers with this seed, then store the resulting state.
    // This is the most economical way to go about the process in terms of overall memory
    // transactions--the state is 256 bits loaded and stored, and each random number is 32 or 64
    // bits (depending on the mode), so best to store at least a handful of random numbers with
    // each checkout of the state.  On the CPU, this process is complicated by the fact that every
    // write of a random number result is going to be a cache miss.  On the GPU, memory coalescence
    // will be ideal.  Store the evolved state vector back in its original place.
    size_t pos = i;
    switch (method) {
    case RandomAlgorithm::XOROSHIRO_128P:
      {
        xrs128p.setState(state_xy[i]);
        if (t_is_double) {
          switch (product) {
          case RandomNumberKind::UNIFORM:
            for (size_t j = 0; j < depth; j++) {
              cache[pos] = xrs128p.uniformRandomNumber();
              pos += padded_atom_count;
            }
            break;
          case RandomNumberKind::GAUSSIAN:
            for (size_t j = 0; j < depth; j++) {
              cache[pos] = xrs128p.gaussianRandomNumber();
              pos += padded_atom_count;
            }
            break;
          }
        }
        else {
          switch (product) {
          case RandomNumberKind::UNIFORM:
            for (size_t j = 0; j < depth; j++) {
              cache[pos] = xrs128p.spUniformRandomNumber();
              pos += padded_atom_count;
            }
            break;
          case RandomNumberKind::GAUSSIAN:
            for (size_t j = 0; j < depth; j++) {
              cache[pos] = xrs128p.spGaussianRandomNumber();
              pos += padded_atom_count;
            }
            break;
          }
        }
        state_xy[i] = xrs128p.revealState();
      }
      break;
    case RandomAlgorithm::XOSHIRO_256PP:
      {
        xrs256pp.setState({ state_xy[i].x, state_xy[i].y, state_zw[i].x, state_zw[i].y });
        if (t_is_double) {
          switch (product) {
          case RandomNumberKind::UNIFORM:
            for (size_t j = 0; j < depth; j++) {
              cache[pos] = xrs256pp.uniformRandomNumber();
              pos += padded_atom_count;
            }
            break;
          case RandomNumberKind::GAUSSIAN:
            for (size_t j = 0; j < depth; j++) {
              cache[pos] = xrs256pp.gaussianRandomNumber();
              pos += padded_atom_count;
            }
            break;
          }
        }
        else {
          switch (product) {
          case RandomNumberKind::UNIFORM:
            for (size_t j = 0; j < depth; j++) {
              cache[pos] = xrs256pp.spUniformRandomNumber();
              pos += padded_atom_count;
            }
            break;
          case RandomNumberKind::GAUSSIAN:
            for (size_t j = 0; j < depth; j++) {
              cache[pos] = xrs256pp.spGaussianRandomNumber();
              pos += padded_atom_count;
            }
            break;
          }
        }
        const ullint4 current_state = xrs256pp.revealState();
        state_xy[i] = { current_state.x, current_state.y };
        state_zw[i] = { current_state.x, current_state.y };
      }
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void fillRandomCache(std::vector<ullint2> *state_xy, std::vector<ullint2> *state_zw,
                     std::vector<T> *cache, const size_t length, const size_t depth,
                     const RandomAlgorithm method, const RandomNumberKind product,
                     const size_t index_start, const size_t index_end) {
  fillRandomCache(state_xy->data(), state_zw->data(), cache->data(), length, depth, method,
                  product, index_start, index_end);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void fillRandomCache(Hybrid<ullint2> *state_xy, Hybrid<ullint2> *state_zw, Hybrid<T> *cache,
                     const size_t length, const size_t depth, const RandomAlgorithm method,
                     const RandomNumberKind product, const size_t index_start,
                     const size_t index_end) {
  fillRandomCache(state_xy->data(), state_zw->data(), cache->data(), length, depth, method,
                  product, index_start, index_end);
}

} // namespace random
} // namespace stormm
