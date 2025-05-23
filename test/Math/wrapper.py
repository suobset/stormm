# pocketfft_wrapper.py
import cffi
import numpy as np
import os

ffi = cffi.FFI()

# Define the C interface we need from pocketfft
ffi.cdef("""
    typedef struct {
        size_t size;
    } shape_t;

    typedef struct {
        size_t size;
    } stride_t;

    void r2c(const shape_t shape, const stride_t stride_in, const stride_t stride_out,
             int axis, bool forward, const double *in, complex double *out, double fct);

    void c2r(const shape_t shape, const stride_t stride_in, const stride_t stride_out,
             int axis, bool forward, const complex double *in, double *out, double fct);

    void c2c(const shape_t shape, const stride_t stride_in, const stride_t stride_out,
             const shape_t axes, bool forward, const complex double *in, 
             complex double *out, double fct);
""")

# Create the wrapper class for PocketFFT
class PocketFFT:
    def __init__(self):
        # Load the PocketFFT library
        self.ffi = ffi
        # Note: You'll need to specify the actual path to your pocketfft header
        self.lib = ffi.dlopen("path/to/pocketfft_hdronly.h")

    def real_fft(self, data):
        """Perform real-to-complex FFT using PocketFFT"""
        shape = np.array([data.shape[0]], dtype=np.size_t)
        stride_in = np.array([sizeof('double')], dtype=np.size_t)
        stride_out = np.array([sizeof('complex double')], dtype=np.size_t)

        # Prepare input and output arrays
        in_data = ffi.from_buffer(data.astype(np.float64))
        out_size = data.shape[0] // 2 + 1
        out_data = np.zeros(out_size, dtype=np.complex128)
        out_ptr = ffi.from_buffer(out_data)

        # Create shape and stride structs
        shape_struct = ffi.new("shape_t*")
        shape_struct.size = data.shape[0]
        stride_in_struct = ffi.new("stride_t*")
        stride_in_struct.size = sizeof('double')
        stride_out_struct = ffi.new("stride_t*")
        stride_out_struct.size = sizeof('complex double')

        # Call PocketFFT r2c
        self.lib.r2c(shape_struct[0], stride_in_struct[0], stride_out_struct[0],
                     0, True, in_data, out_ptr, 1.0)

        return out_data

    def complex_fft(self, data):
        """Perform complex-to-complex FFT using PocketFFT"""
        shape = np.array([data.shape[0]], dtype=np.size_t)
        stride = np.array([sizeof('complex double')], dtype=np.size_t)

        # Prepare input and output arrays
        in_data = ffi.from_buffer(data.astype(np.complex128))
        out_data = np.zeros_like(data, dtype=np.complex128)
        out_ptr = ffi.from_buffer(out_data)

        # Create shape and stride structs
        shape_struct = ffi.new("shape_t*")
        shape_struct.size = data.shape[0]
        stride_struct = ffi.new("stride_t*")
        stride_struct.size = sizeof('complex double')
        axes_struct = ffi.new("shape_t*")
        axes_struct.size = 0

        # Call PocketFFT c2c
        self.lib.c2c(shape_struct[0], stride_struct[0], stride_struct[0],
                     axes_struct[0], True, in_data, out_ptr, 1.0)

        return out_data