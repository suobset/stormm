import numpy as np
from numpy.random import default_rng
import numpy.testing as npt
import cffi
import os

# Initialize CFFI
ffi = cffi.FFI()

# Define the C interface directly
ffi.cdef("""
    void r2c(const size_t shape[], const ptrdiff_t stride_in[], const ptrdiff_t stride_out[],
             size_t axes[], bool forward, const double *in, complex double *out, double fct);

    void c2r(const size_t shape[], const ptrdiff_t stride_in[], const ptrdiff_t stride_out[],
             size_t axes[], bool forward, const complex double *in, double *out, double fct);

    void c2c(const size_t shape[], const ptrdiff_t stride_in[], const ptrdiff_t stride_out[],
             size_t axes[], bool forward, const complex double *in, complex double *out, double fct);
""")

# Load the PocketFFT header-only library
# Note: Update this path to where your pocketfft header is located
lib = ffi.dlopen("path/to/pocketfft_hdronly.h")

class TestFFT:
    def __init__(self, seed=7183529):
        self.rng = default_rng(seed)

    def test_real_fft(self, maxlen=8192, epsilon=5.0e-10):
        """
        Tests the real FFT (RFFT) functionality using PocketFFT directly.
        """
        test_lengths = []
        tl = 1
        while tl < maxlen:
            test_lengths.append(tl)
            if tl < 8:
                tl += 1
            elif tl < maxlen // 8:
                tl = maxlen // 8
            else:
                tl += maxlen // 8

        for length in test_lengths:
            # Generate random data
            grid_values = self.rng.normal(0, 1.0, length)
            original_grid_values = grid_values.copy()

            # Prepare arrays for FFT
            shape = ffi.new("size_t[]", [length])
            stride_in = ffi.new("ptrdiff_t[]", [8])  # sizeof(double)
            stride_out = ffi.new("ptrdiff_t[]", [16])  # sizeof(complex double)
            axes = ffi.new("size_t[]", [0])

            # Create output array for complex data
            complex_size = length // 2 + 1
            complex_data = np.zeros(complex_size, dtype=np.complex128)

            # Get pointers to data
            in_ptr = ffi.cast("double *", ffi.from_buffer(grid_values))
            out_ptr = ffi.cast("complex double *", ffi.from_buffer(complex_data))

            # Forward transform
            lib.r2c(shape, stride_in, stride_out, axes, True, in_ptr, out_ptr, 1.0)

            # Prepare for inverse transform
            fft_data = np.zeros(length, dtype=np.float64)
            back_ptr = ffi.cast("double *", ffi.from_buffer(fft_data))

            # Inverse transform
            lib.c2r(shape, stride_out, stride_in, axes, False, out_ptr, back_ptr, 1.0/length)

            # Compare results
            npt.assert_allclose(
                fft_data, 
                original_grid_values, 
                rtol=epsilon, 
                err_msg=f"RFFT test failed for length {length}"
            )

    def test_complex_fft(self, maxlen=8192, epsilon=5.0e-10):
        """
        Tests the complex FFT (CFFT) functionality using PocketFFT directly.
        """
        test_lengths = []
        tl = 1
        while tl < maxlen:
            test_lengths.append(tl)
            if tl < 4:
                tl += 1
            elif tl < maxlen // 8:
                tl = maxlen // 8
            else:
                tl += maxlen // 8

        for length in test_lengths:
            # Generate complex random data
            complex_data = (self.rng.normal(0, 1.0, length) + 
                          1j * self.rng.normal(0, 1.0, length))
            original_data = complex_data.copy()

            # Prepare arrays for FFT
            shape = ffi.new("size_t[]", [length])
            stride = ffi.new("ptrdiff_t[]", [16])  # sizeof(complex double)
            axes = ffi.new("size_t[]", [0])

            # Create output array
            fft_data = np.zeros(length, dtype=np.complex128)

            # Get pointers to data
            in_ptr = ffi.cast("complex double *", ffi.from_buffer(complex_data))
            out_ptr = ffi.cast("complex double *", ffi.from_buffer(fft_data))

            # Forward transform
            lib.c2c(shape, stride, stride, axes, True, in_ptr, out_ptr, 1.0)

            # Inverse transform
            back_data = np.zeros(length, dtype=np.complex128)
            back_ptr = ffi.cast("complex double *", ffi.from_buffer(back_data))
            lib.c2c(shape, stride, stride, axes, False, out_ptr, back_ptr, 1.0/length)

            # Compare results
            npt.assert_allclose(
                back_data, 
                original_data, 
                rtol=epsilon, 
                err_msg=f"CFFT test failed for length {length}"
            )

    def test_real_3d_fft(self, dims=(8, 8, 8), epsilon=5.0e-10):
        """
        Tests the 3D real FFT functionality using PocketFFT directly.
        """
        if len(dims) != 3:
            raise ValueError(f"Expected 3 dimensions, got {len(dims)}")
        
        if any(d <= 0 for d in dims):
            raise ValueError(f"All dimensions must be positive, got {dims}")

        # Generate 3D random data
        total_size = dims[0] * dims[1] * dims[2]
        grid_values = self.rng.normal(0, 1.0, total_size).reshape(dims)
        original_grid_values = grid_values.copy()

        # Prepare arrays for FFT
        shape = ffi.new("size_t[]", dims)
        stride_in = ffi.new("ptrdiff_t[]", [8, 8 * dims[0], 8 * dims[0] * dims[1]])
        stride_out = ffi.new("ptrdiff_t[]", [16, 16 * (dims[0]//2 + 1), 16 * (dims[0]//2 + 1) * dims[1]])
        axes = ffi.new("size_t[]", [0, 1, 2])

        # Create output array
        complex_shape = (dims[0]//2 + 1, dims[1], dims[2])
        complex_data = np.zeros(complex_shape, dtype=np.complex128)

        # Get pointers to data
        in_ptr = ffi.cast("double *", ffi.from_buffer(grid_values))
        out_ptr = ffi.cast("complex double *", ffi.from_buffer(complex_data))

        # Forward transform
        lib.r2c(shape, stride_in, stride_out, axes, True, in_ptr, out_ptr, 1.0)

        # Inverse transform
        back_data = np.zeros(dims, dtype=np.float64)
        back_ptr = ffi.cast("double *", ffi.from_buffer(back_data))
        lib.c2r(shape, stride_out, stride_in, axes, False, out_ptr, back_ptr, 1.0/total_size)

        # Compare results
        npt.assert_allclose(
            back_data, 
            original_grid_values, 
            rtol=epsilon, 
            err_msg=f"3D RFFT test failed for dimensions {dims}"
        )

def main():
    # Initialize test environment
    fft_tester = TestFFT(seed=7183529)
    
    print("Running a real FFT test using PocketFFT")
    fft_tester.test_real_fft()
    
    print("Running a complex FFT test using PocketFFT")
    fft_tester.test_complex_fft()
    
    print("Running a 3D real FFT test using PocketFFT")
    fft_tester.test_real_3d_fft()

if __name__ == "__main__":
    main()