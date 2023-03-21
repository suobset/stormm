#include "copyright.h"
#include "matrix_ops.h"
#include "tricubic_cell.h"

namespace stormm {
namespace stmath {

//-------------------------------------------------------------------------------------------------
Hybrid<double> getTricubicMatrix() {
  std::vector<double> tmp(4096);
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      for (int k = 0; k < 4 ; k++) {
        const int col_idx = (4 * ((4 * k) + j)) + i;
        for (int xv = 0; xv < 2; xv++) {
          const int xfac = (i == 0) ? 1 : (i == 1) ? xv : (i == 2) ? xv * xv : xv * xv * xv;
          const int dxfac = (i == 0) ? 0 : (i == 1) ? 1 : (i == 2) ? 2 * xv : 3 * xv * xv;
          for (int yv = 0; yv < 2; yv++) {
            const int yfac = (j == 0) ? 1 : (j == 1) ? yv : (j == 2) ? yv * yv : yv * yv * yv;
            const int dyfac = (j == 0) ? 0 : (j == 1) ? 1 : (j == 2) ? 2 * yv : 3 * yv * yv;
            for (int zv = 0; zv < 2; zv++) {
              const int row_idx = (2 * ((2 * zv) + yv)) + xv;
              const int zfac = (k == 0) ? 1 : (k == 1) ? zv : (k == 2) ? zv * zv : zv * zv * zv;
              const int dzfac = (k == 0) ? 0 : (k == 1) ? 1 : (k == 2) ? 2 * zv : 3 * zv * zv;

              // Compute the polynomial coefficient for the function value
              tmp[(col_idx * 64) + row_idx] = xfac * yfac * zfac;

              // Compute the polynomial coefficients for function first derivatives
              tmp[(col_idx * 64) + row_idx +  8] = dxfac * yfac * zfac;
              tmp[(col_idx * 64) + row_idx + 16] = xfac * dyfac * zfac;
              tmp[(col_idx * 64) + row_idx + 24] = xfac * yfac * dzfac;

              // Compute the polynomial coefficients for function cross derivatives
              tmp[(col_idx * 64) + row_idx + 32] = dxfac * dyfac * zfac;
              tmp[(col_idx * 64) + row_idx + 40] = dxfac * yfac * dzfac;
              tmp[(col_idx * 64) + row_idx + 48] = xfac * dyfac * dzfac;

              // Compute the polynomial coefficient for the function triple derivative
              tmp[(col_idx * 64) + row_idx + 56] = dxfac * dyfac * dzfac;
            }
          }
        }
      }
    }
  }

  // Invert the matrix to obtain the tricubic spline coefficients
  std::vector<double> itmp(4096);
  const double* tmp_ptr = tmp.data();
  invertSquareMatrix(tmp_ptr, itmp.data(), 64);
  Hybrid<double> result(4096, "tricubic_eval");
  result.putHost(itmp);

  return result;
}

} // namespace stmath
} // namespace stormm
