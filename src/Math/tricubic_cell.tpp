// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace math {

//-------------------------------------------------------------------------------------------------
template <typename T>
TricubicCell<T>::TricubicCell() :
    coefficients{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
    f{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }, dx{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
    dy{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }, dz{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
    dxx{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }, dyy{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
    dzz{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }, dxy{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
    dxz{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }, dyz{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
    dxyz{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }, origin_x{0.0}, origin_y{0.0}, origin_z{0.0},
    umat{ 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 },
    invu{ 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 }
{}
    
//-------------------------------------------------------------------------------------------------
template <typename T>
TricubicCell<T>::TricubicCell(const std::vector<double> weights_matrix,
                              const std::vector<double> &bounds, const std::vector<T> &f_in,
                              const std::vector<T> &dx_in, const std::vector<T> &dy_in,
                              const std::vector<T> &dz_in, const std::vector<T> &dxy_in,
                              const std::vector<T> &dxz_in, const std::vector<T> &dyz_in,
                              const std::vector<T> &dxyz_in, const std::vector<T> &dxx_in,
                              const std::vector<T> &dyy_in, const std::vector<T> &dzz_in) :
    coefficients{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
    f{ f_in[0], f_in[1], f_in[2], f_in[3], f_in[4], f_in[5], f_in[6], f_in[7] },
    dx{ dx_in[0], dx_in[1], dx_in[2], dx_in[3], dx_in[4], dx_in[5], dx_in[6], dx_in[7] },
    dy{ dy_in[0], dy_in[1], dy_in[2], dy_in[3], dy_in[4], dy_in[5], dy_in[6], dy_in[7] },
    dz{ dz_in[0], dz_in[1], dz_in[2], dz_in[3], dz_in[4], dz_in[5], dz_in[6], dz_in[7] },
    dxx{ dxx_in[0], dxx_in[1], dxx_in[2], dxx_in[3], dxx_in[4], dxx_in[5], dxx_in[6], dxx_in[7] },
    dyy{ dyy_in[0], dyy_in[1], dyy_in[2], dyy_in[3], dyy_in[4], dyy_in[5], dyy_in[6], dyy_in[7] },
    dzz{ dzz_in[0], dzz_in[1], dzz_in[2], dzz_in[3], dzz_in[4], dzz_in[5], dzz_in[6], dzz_in[7] },
    dxy{ dxy_in[0], dxy_in[1], dxy_in[2], dxy_in[3], dxy_in[4], dxy_in[5], dxy_in[6], dxy_in[7] },
    dxz{ dxz_in[0], dxz_in[1], dxz_in[2], dxz_in[3], dxz_in[4], dxz_in[5], dxz_in[6], dxz_in[7] },
    dyz{ dyz_in[0], dyz_in[1], dyz_in[2], dyz_in[3], dyz_in[4], dyz_in[5], dyz_in[6], dyz_in[7] },
    dxyz{ dxyz_in[0], dxyz_in[1], dxyz_in[2], dxyz_in[3], dxyz_in[4], dxyz_in[5], dxyz_in[6],
          dxyz_in[7] }, origin_x{0.0}, origin_y{0.0}, origin_z{0.0},
    umat{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
    invu{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
{
  // Fill out the grid cell's dimensions
  if (bounds.size() > 3) {
    origin_x = bounds[0];
    origin_y = bounds[1];
    origin_z = bounds[2];
  }
  if (bounds.size() == 4LLU) {
    invu[0] = bounds[3];
    invu[4] = bounds[3];
    invu[8] = bounds[3];
  }
  else if (bounds.size() == 6LLU) {
    invu[0] = bounds[3];
    invu[4] = bounds[4];
    invu[8] = bounds[5];
  }
  else if (bounds.size() == 12LLU) {
    for (int i = 0; i < 9; i++) {
      invu[i] = bounds[i + 3];
    }
  }
  else {
    rtErr("Invalid dimesion " + std::to_string(bounds.size()) + " on the bounds array.  There "
          "must be four or six entries.", "TricubicCell");
  }

  // Check the cell lengths
  if (leibnizDeterminant(invu, 3) < 1.0e-8) {
    rtErr("An invalid cell was defined with vectors [ " +
          realToString(invu[0], 7, 4, NumberFormat::STANDARD_REAL) + ", " +
          realToString(invu[1], 7, 4, NumberFormat::STANDARD_REAL) + ", " +
          realToString(invu[2], 7, 4, NumberFormat::STANDARD_REAL) + "], [" +
          realToString(invu[3], 7, 4, NumberFormat::STANDARD_REAL) + ", " +
          realToString(invu[4], 7, 4, NumberFormat::STANDARD_REAL) + ", " +
          realToString(invu[5], 7, 4, NumberFormat::STANDARD_REAL) + "], and [" +
          realToString(invu[6], 7, 4, NumberFormat::STANDARD_REAL) + ", " +
          realToString(invu[7], 7, 4, NumberFormat::STANDARD_REAL) + ", " +
          realToString(invu[8], 7, 4, NumberFormat::STANDARD_REAL) + "].", "TricubicCell");
  }
  else {
    const std::vector<T> invu_copy = { invu[0], invu[1], invu[2], invu[3], invu[4], invu[5],
                                       invu[6], invu[7], invu[8] };
    invertSquareMatrix(invu_copy.data(), umat, 3);
  }
  
  // Fill out the coefficients matrix
  std::vector<double> b(64);
  for (int i = 0; i < 8; i++) {
    b[i     ] = f[i];
    if (bounds.size() <= 6) {
      b[i +  8] = dx_in[i] * invu[0];
      b[i + 16] = dy_in[i] * invu[4];
      b[i + 24] = dz_in[i] * invu[8];
      b[i + 32] = dxy_in[i] * invu[0] * invu[4];
      b[i + 40] = dxz_in[i] * invu[0] * invu[8];
      b[i + 48] = dyz_in[i] * invu[4] * invu[8];
      b[i + 56] = dxyz_in[i] * invu[0] * invu[4] * invu[8];
    }
    else {

      // The general, triclinic case
      b[i +  8] = (dx_in[i] * invu[0]) + (dy_in[i] * invu[3]) + (dz_in[i] * invu[6]);
      b[i + 16] = (dx_in[i] * invu[1]) + (dy_in[i] * invu[4]) + (dz_in[i] * invu[7]);
      b[i + 24] = (dx_in[i] * invu[2]) + (dy_in[i] * invu[5]) + (dz_in[i] * invu[8]);
    }
  }
  std::vector<double> dcoeffs(64);
  matrixVectorMultiply(weights_matrix.data(), b.data(), dcoeffs.data(), 64, 64);
  for (int i = 0; i < 64; i++) {
    coefficients[i] = dcoeffs[i];
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
T TricubicCell<T>::getCoefficient(const int i, const int j, const int k) const {
  if (i > 3 || j > 3 || k > 3 || i < 0 || j < 0 || k < 0) {
    rtErr("A coefficient for x, y, and z powers " + std::to_string(i) + ", " + std::to_string(j) +
          ", and " + std::to_string(k) + " is not acceptable.", "TricubicCell", "getCoefficient");
  }
  return coefficients[(4 * ((4 * k) + j)) + i];
}

//-------------------------------------------------------------------------------------------------
template <typename T>
std::vector<T> TricubicCell<T>::getCoefficients() const {
  std::vector<T> result(64);
  for (size_t i = 0; i < 64LLU; i++) {
    result[i] = coefficients[i];
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void TricubicCell<T>::setCoefficient(const T value, const int i, const int j, const int k) {
  if (i > 3 || j > 3 || k > 3 || i < 0 || j < 0 || k < 0) {
    rtErr("A coefficient for x, y, and z powers " + std::to_string(i) + ", " + std::to_string(j) +
          ", and " + std::to_string(k) + " is not acceptable.", "TricubicCell", "getCoefficient");
  }
  coefficients[(4 * ((4 * k) + j)) + i] = value;
}

//-------------------------------------------------------------------------------------------------
template <typename T>
T TricubicCell<T>::getData(const TricubicBound kind, const int i, const int j, const int k) const {
  switch (kind) {
  case TricubicBound::VALUE:
    return f[(2 * ((2 * i) + j)) + k];
  case TricubicBound::DX:
    return dx[(2 * ((2 * i) + j)) + k];
  case TricubicBound::DY:
    return dy[(2 * ((2 * i) + j)) + k];
  case TricubicBound::DZ:
    return dz[(2 * ((2 * i) + j)) + k];
  case TricubicBound::DXY:
    return dxy[(2 * ((2 * i) + j)) + k];
  case TricubicBound::DXZ:
    return dxz[(2 * ((2 * i) + j)) + k];
  case TricubicBound::DYZ:
    return dyz[(2 * ((2 * i) + j)) + k];
  case TricubicBound::DXYZ:
    return dxyz[(2 * ((2 * i) + j)) + k];
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename T> void TricubicCell<T>::setData(const T value, const TricubicBound kind,
                                                    const int i, const int j, const int k) {
  switch (kind) {
  case TricubicBound::VALUE:
    f[(2 * ((2 * i) + j)) + k] = value;
  case TricubicBound::DX:
    dx[(2 * ((2 * i) + j)) + k] = value;
  case TricubicBound::DY:
    dy[(2 * ((2 * i) + j)) + k] = value;
  case TricubicBound::DZ:
    dz[(2 * ((2 * i) + j)) + k] = value;
  case TricubicBound::DXY:
    dxy[(2 * ((2 * i) + j)) + k] = value;
  case TricubicBound::DXZ:
    dxz[(2 * ((2 * i) + j)) + k] = value;
  case TricubicBound::DYZ:
    dyz[(2 * ((2 * i) + j)) + k] = value;
  case TricubicBound::DXYZ:
    dxyz[(2 * ((2 * i) + j)) + k] = value;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename T> T TricubicCell<T>::getCellOrigin(CartesianDimension dim) const {
  switch (dim) {
  case CartesianDimension::X:
    return origin_x;
  case CartesianDimension::Y:
    return origin_y;
  case CartesianDimension::Z:
    return origin_z;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename T> T TricubicCell<T>::getCellLength(CartesianDimension dim) const {
  switch (dim) {
  case CartesianDimension::X:
    return invu[0];
  case CartesianDimension::Y:
    return invu[4];
  case CartesianDimension::Z:
    return invu[8];
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename T> T TricubicCell<T>::evaluate(const T x, const T y, const T z,
                                                  const TricubicBound kind) const {
  const T xrel = x - origin_x;
  const T yrel = y - origin_y;
  const T zrel = z - origin_z;
  const T x_frac = (umat[0] * xrel) + (umat[3] * yrel) + (umat[6] * zrel);
  const T y_frac = (umat[1] * xrel) + (umat[4] * yrel) + (umat[7] * zrel);
  const T z_frac = (umat[2] * xrel) + (umat[5] * yrel) + (umat[8] * zrel);
  const T zero = 0.0;
  const T one  = 1.0;
  if (x_frac < zero || x_frac >= one || y_frac < zero || y_frac >= one || z_frac < zero ||
      z_frac >= one) {
    rtErr("Point (" + realToString(x, 7, 4, NumberFormat::STANDARD_REAL) + ", " +
          realToString(y, 7, 4, NumberFormat::STANDARD_REAL) + ", " +
          realToString(z, 7, 4, NumberFormat::STANDARD_REAL) + ") is out of bounds for a grid "
          "cell with origin (" + realToString(origin_x, 7, 4, NumberFormat::STANDARD_REAL) + ", " +
          realToString(origin_y, 7, 4, NumberFormat::STANDARD_REAL) + ", " +
          realToString(origin_z, 7, 4, NumberFormat::STANDARD_REAL) + ") and dimensions [ " +
          realToString(invu[0], 7, 4, NumberFormat::STANDARD_REAL) + ", " +
          realToString(invu[1], 7, 4, NumberFormat::STANDARD_REAL) + ", " +
          realToString(invu[2], 7, 4, NumberFormat::STANDARD_REAL) + " ] x [ " +
          realToString(invu[3], 7, 4, NumberFormat::STANDARD_REAL) + ", " +
          realToString(invu[4], 7, 4, NumberFormat::STANDARD_REAL) + ", " +
          realToString(invu[5], 7, 4, NumberFormat::STANDARD_REAL) + " ] x [ " +
          realToString(invu[6], 7, 4, NumberFormat::STANDARD_REAL) + ", " +
          realToString(invu[7], 7, 4, NumberFormat::STANDARD_REAL) + ", " +
          realToString(invu[8], 7, 4, NumberFormat::STANDARD_REAL) + " ].", "TricubicCell",
          "evaluate");
  }
  T result = 0.0;
  const T v_one = 1.0; 
  T xv = v_one;
  switch (kind) {
  case TricubicBound::VALUE:
    for (int i = 0; i < 4; i++) {
      T yv = v_one;
      for (int j = 0; j < 4; j++) {
        T zv = v_one;
        for (int k = 0; k < 4; k++) {
          result += coefficients[(4 * ((4 * k) + j)) + i] * xv * yv * zv;
          zv *= z_frac;
        }
        yv *= y_frac;
      }
      xv *= x_frac;
    }
    break;
  case TricubicBound::DX:
    for (int i = 1; i < 4; i++) {
      T yv = v_one;
      const double difac = i;
      for (int j = 0; j < 4; j++) {
        T zv = v_one;
        for (int k = 0; k < 4; k++) {
          result += coefficients[(4 * ((4 * k) + j)) + i] * difac * xv * yv * zv;
          zv *= z_frac;
        }
        yv *= y_frac;
      }
      xv *= x_frac;
    }
    result *= umat[0];
    break;
  case TricubicBound::DY:
    for (int i = 0; i < 4; i++) {
      T yv = v_one;
      double djfac = 1.0;
      for (int j = 1; j < 4; j++) {
        T zv = v_one;
        for (int k = 0; k < 4; k++) {
          result += coefficients[(4 * ((4 * k) + j)) + i] * xv * djfac * yv * zv;
          zv *= z_frac;
        }
        yv *= y_frac;
        djfac += 1.0;
      }
      xv *= x_frac;
    }
    result *= umat[4];
    break;
  case TricubicBound::DZ:
    for (int i = 0; i < 4; i++) {
      T yv = v_one;
      for (int j = 0; j < 4; j++) {
        T zv = v_one;
        double dkfac = 1.0;
        for (int k = 1; k < 4; k++) {
          result += coefficients[(4 * ((4 * k) + j)) + i] * xv * yv * dkfac * zv;
          zv *= z_frac;
          dkfac += 1.0;
        }
        yv *= y_frac;
      }
      xv *= x_frac;
    }
    result *= umat[8];
    break;
  case TricubicBound::DXY:
    for (int i = 1; i < 4; i++) {
      T yv = v_one;
      const double difac = i;
      double djfac = 1.0;
      for (int j = 1; j < 4; j++) {
        T zv = v_one;
        for (int k = 0; k < 4; k++) {
          result += coefficients[(4 * ((4 * k) + j)) + i] * difac * xv * djfac * yv * zv;
          zv *= z_frac;
        }
        yv *= y_frac;
        djfac += 1.0;
      }
      xv *= x_frac;
    }
    result *= umat[0] * umat[4];
    break;
  case TricubicBound::DXZ:
    for (int i = 1; i < 4; i++) {
      T yv = v_one;
      const double difac = i;
      for (int j = 0; j < 4; j++) {
        T zv = v_one;
        double dkfac = 1.0;
        for (int k = 1; k < 4; k++) {
          result += coefficients[(4 * ((4 * k) + j)) + i] * difac * xv * yv * dkfac * zv;
          zv *= z_frac;
          dkfac += 1.0;
        }
        yv *= y_frac;
      }
      xv *= x_frac;
    }
    result *= umat[0] * umat[8];
    break;
  case TricubicBound::DYZ:
    for (int i = 0; i < 4; i++) {
      T yv = v_one;
      double djfac = 1.0;
      for (int j = 1; j < 4; j++) {
        T zv = v_one;
        double dkfac = 1.0;
        for (int k = 1; k < 4; k++) {
          result += coefficients[(4 * ((4 * k) + j)) + i] * xv * djfac * yv * dkfac * zv;
          zv *= z_frac;
          dkfac += 1.0;
        }
        yv *= y_frac;
        djfac += 1.0;
      }
      xv *= x_frac;
    }
    result *= umat[4] * umat[8];
    break;
  case TricubicBound::DXYZ:
    for (int i = 1; i < 4; i++) {
      T yv = v_one;
      const double difac = i;
      double djfac = 1.0;
      for (int j = 1; j < 4; j++) {
        T zv = v_one;
        double dkfac = 1.0;
        for (int k = 1; k < 4; k++) {
          result += coefficients[(4 * ((4 * k) + j)) + i] * difac * xv * djfac * yv * dkfac * zv;
          zv *= z_frac;
          dkfac += 1.0;
        }
        yv *= y_frac;
        djfac += 1.0;
      }
      xv *= x_frac;
    }
    result *= umat[0] * umat[4] * umat[8];
    break;
  }
  return result;
}

} // namespace math
} // namespace stormm
