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
    dxy{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }, dxz{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
    dyz{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }, dxyz{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
    origin_x{0.0}, origin_y{0.0}, origin_z{0.0}, length_x{1.0}, length_y{1.0}, length_z{1.0}
{}
    
//-------------------------------------------------------------------------------------------------
template <typename T>
TricubicCell<T>::TricubicCell(const std::vector<double> weights_matrix,
                              const std::vector<double> &bounds, const std::vector<T> &f_in,
                              const std::vector<T> &dx_in, const std::vector<T> &dy_in,
                              const std::vector<T> &dz_in, const std::vector<T> &dxy_in,
                              const std::vector<T> &dxz_in, const std::vector<T> &dyz_in,
                              const std::vector<T> &dxyz_in) :
    coefficients{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
    f{ f_in[0], f_in[1], f_in[2], f_in[3], f_in[4], f_in[5], f_in[6], f_in[7] },
    dx{ dx_in[0], dx_in[1], dx_in[2], dx_in[3], dx_in[4], dx_in[5], dx_in[6], dx_in[7] },
    dy{ dy_in[0], dy_in[1], dy_in[2], dy_in[3], dy_in[4], dy_in[5], dy_in[6], dy_in[7] },
    dz{ dz_in[0], dz_in[1], dz_in[2], dz_in[3], dz_in[4], dz_in[5], dz_in[6], dz_in[7] },
    dxy{ dxy_in[0], dxy_in[1], dxy_in[2], dxy_in[3], dxy_in[4], dxy_in[5], dxy_in[6], dxy_in[7] },
    dxz{ dxz_in[0], dxz_in[1], dxz_in[2], dxz_in[3], dxz_in[4], dxz_in[5], dxz_in[6], dxz_in[7] },
    dyz{ dyz_in[0], dyz_in[1], dyz_in[2], dyz_in[3], dyz_in[4], dyz_in[5], dyz_in[6], dyz_in[7] },
    dxyz{ dxyz_in[0], dxyz_in[1], dxyz_in[2], dxyz_in[3], dxyz_in[4], dxyz_in[5], dxyz_in[6],
          dxyz_in[7] }, origin_x{0.0}, origin_y{0.0}, origin_z{0.0}, length_x{1.0}, length_y{1.0},
    length_z{1.0}
{
  // Fill out the coefficients matrix
  std::vector<double> b(64);
  for (int i = 0; i < 8; i++) {
    b[i     ] = f[i];
    b[i +  8] = dx[i];
    b[i + 16] = dy[i];
    b[i + 24] = dz[i];
    b[i + 32] = dxy[i];
    b[i + 40] = dxz[i];
    b[i + 48] = dyz[i];
    b[i + 56] = dxyz[i];
  }
  std::vector<double> dcoeffs(64);
  matrixVectorMultiply(weights_matrix.data(), b.data(), dcoeffs.data(), 64, 64);
  for (int i = 0; i < 64; i++) {
    coefficients[i] = dcoeffs[i];
  }

  // Fill out the grid cell's dimensions
  if (bounds.size() == 4LLU) {
    origin_x = bounds[0];
    origin_y = bounds[1];
    origin_z = bounds[2];
    length_x = bounds[3];
    length_y = bounds[3];
    length_z = bounds[3];
  }
  else if (bounds.size() == 6LLU) {
    origin_x = bounds[0];
    origin_y = bounds[1];
    origin_z = bounds[2];
    length_x = bounds[3];
    length_y = bounds[4];
    length_z = bounds[5];    
  }
  else {
    rtErr("Invalid dimesion " + std::to_string(bounds.size()) + " on the bounds array.  There "
          "must be four or six entries.", "TricubicCell");
  }
  
  // Check the cell lengths
  if (length_x < 0.0 || length_y < 0.0 || length_z < 0.0) {
    rtErr("Invalid cell lengths " + realToString(length_x, 7, 4, NumberFormat::STANDARD_REAL) +
          " x " + realToString(length_y, 7, 4, NumberFormat::STANDARD_REAL) + " x " +
          realToString(length_z, 7, 4, NumberFormat::STANDARD_REAL) + ".", "TricubicCell");
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
    return length_x;
  case CartesianDimension::Y:
    return length_y;
  case CartesianDimension::Z:
    return length_z;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename T> T TricubicCell<T>::evaluate(const T x, const T y, const T z) const {
  if (x < origin_x || y < origin_y || z < origin_z ||
      x >= origin_x + length_x || y >= origin_y + length_y || z >= origin_z + length_z) {
    rtErr("Point (" + realToString(x, 7, 4, NumberFormat::STANDARD_REAL) + ", " +
          realToString(y, 7, 4, NumberFormat::STANDARD_REAL) + ", " +
          realToString(z, 7, 4, NumberFormat::STANDARD_REAL) + ") is out of bounds for a grid "
          "cell with origin (" + realToString(origin_x, 7, 4, NumberFormat::STANDARD_REAL) + ", " +
          realToString(origin_y, 7, 4, NumberFormat::STANDARD_REAL) + ", " +
          realToString(origin_z, 7, 4, NumberFormat::STANDARD_REAL) + ") and dimensions " +
          realToString(length_x, 7, 4, NumberFormat::STANDARD_REAL) + " x " +
          realToString(length_y, 7, 4, NumberFormat::STANDARD_REAL) + " x " +
          realToString(length_z, 7, 4, NumberFormat::STANDARD_REAL) + ".", "TricubicCell",
          "evaluate");
  }
  const T x_frac = (x - origin_x) / length_x;
  const T y_frac = (y - origin_y) / length_y;
  const T z_frac = (z - origin_z) / length_z;
  T result = 0.0;
  const T v_one = 1.0; 
  T xv = v_one;
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
  return result;
}

} // namespace math
} // namespace stormm
