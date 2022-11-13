// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace math {

//-------------------------------------------------------------------------------------------------
template <typename T>
TriCubicCell<T>::TriCubicCell() :
    coefficients{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
    f{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }, dx{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
    dy{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }, dz{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
    dxy{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }, dxz{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
    dyz{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }, dxyz{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
{}
    
//-------------------------------------------------------------------------------------------------
template <typename T>
TriCubicCell<T>::TriCubicCell(const std::vector<double> weights_matrix, const std::vector<T> &f_in,
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
          dxyz_in[7] }
{
  const T zero = 0.0;
  if (sum<T>(f, 8) == zero && sum<T>(dx, 8) == zero && sum<T>(dy, 8) == zero &&
      sum<T>(dz, 8) == zero && sum<T>(dxy, 8) == zero && sum<T>(dxz, 8) == zero &&
      sum<T>(dyz, 8) == zero && sum<T>(dxyz, 8) == zero) {
    return;
  }
  for (int i = 0; i < 8; i++) {

  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
T TriCubicCell<T>::getCoefficient(const int i, const int j, const int k) {
  if (i > 3 || j > 3 || k > 3 || i < 0 || j < 0 || k < 0) {
    rtErr("A coefficient for x, y, and z powers " + std::to_string(i) + ", " + std::to_string(j) +
          ", and " + std::to_string(k) + " is not acceptable.", "TriCubicCell", "getCoefficient");
  }
  return coefficients[(4 * ((4 * k) + j)) + i];
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void TriCubicCell<T>::setCoefficient(const T value, const int i, const int j, const int k) {
  if (i > 3 || j > 3 || k > 3 || i < 0 || j < 0 || k < 0) {
    rtErr("A coefficient for x, y, and z powers " + std::to_string(i) + ", " + std::to_string(j) +
          ", and " + std::to_string(k) + " is not acceptable.", "TriCubicCell", "getCoefficient");
  }
  coefficients[(4 * ((4 * k) + j)) + i] = value;
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void TriCubicCell::getData(const TriCubicBound kind, const int i,  const int j, const int k) {
  switch (kind) {
  case TriCubicKind::VALUE:
    return f[(2 * ((2 * i) + j)) + k];
  case TriCubicKind::DX:
    return dx[(2 * ((2 * i) + j)) + k];
  case TriCubicKind::DY:
    return dy[(2 * ((2 * i) + j)) + k];
  case TriCubicKind::DZ:
    return dz[(2 * ((2 * i) + j)) + k];
  case TriCubicKind::DXY:
    return dxy[(2 * ((2 * i) + j)) + k];
  case TriCubicKind::DXZ:
    return dxz[(2 * ((2 * i) + j)) + k];
  case TriCubicKind::DYZ:
    return dyz[(2 * ((2 * i) + j)) + k];
  case TriCubicKind::DXYZ:
    return dxyz[(2 * ((2 * i) + j)) + k];
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T> void TriCubicCell::setData(const T value, const TriCubicBound kind,
                                                 const int i, const int j, const int k) {
  switch (kind) {
  case TriCubicKind::VALUE:
    f[(2 * ((2 * i) + j)) + k] = value;
  case TriCubicKind::DX:
    dx[(2 * ((2 * i) + j)) + k] = value;
  case TriCubicKind::DY:
    dy[(2 * ((2 * i) + j)) + k] = value;
  case TriCubicKind::DZ:
    dz[(2 * ((2 * i) + j)) + k] = value;
  case TriCubicKind::DXY:
    dxy[(2 * ((2 * i) + j)) + k] = value;
  case TriCubicKind::DXZ:
    dxz[(2 * ((2 * i) + j)) + k] = value;
  case TriCubicKind::DYZ:
    dyz[(2 * ((2 * i) + j)) + k] = value;
  case TriCubicKind::DXYZ:
    dxyz[(2 * ((2 * i) + j)) + k] = value;
  }
}

} // namespace math
} // namespace stormm
