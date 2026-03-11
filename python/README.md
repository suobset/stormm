# STORMM Python Layer (CPU)

This directory provides a Python interface to existing STORMM C++ code. It does not reimplement
STORMM.

## Build

From the repository root:

```bash
cmake -S . -B build -DSTORMM_ENABLE_CUDA=OFF
cmake --build build --target stormm_pybridge -j4
```

This produces `build/libstormm_pybridge.*` and the Python module loads it with `ctypes`.

## Use

Set `PYTHONPATH` to this repository's `python/` directory:

```bash
PYTHONPATH=python python3 -c "from stormm import dynamics; print(dynamics.state_variable_names()[:5])"
```

## End-to-end Demo

```bash
PYTHONPATH=python python3 python/examples/dynamics_end_to_end.py
```

The demo creates, from Python:

- `AtomGraph`
- `PhaseSpace`
- `AtomGraphSynthesis`
- `PhaseSpaceSynthesis`

Then it runs STORMM `dynamics(...)` and reads results back into Python.
