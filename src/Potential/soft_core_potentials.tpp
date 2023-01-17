// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace energy {

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
void quadraticCoreElectrostatics(const Tcalc r, const Tcalc clash_distance, const Tcalc qiqj,
                                 Tcalc *ele_contrib, Tcalc *fmag) {
  const Tcalc value_one = 1.0;
  if (r < clash_distance) {
    const Tcalc value_half = 0.5;
    const Tcalc value_two  = 2.0;
    const Tcalc aparm = -value_half /
                        (clash_distance * clash_distance * (clash_distance + value_one));
    const Tcalc bparm = (value_one / clash_distance) -
                        (aparm * (clash_distance + value_one) * (clash_distance + value_one));
    *ele_contrib += qiqj * ((aparm * (r + value_one) * (r + value_one)) + bparm);
    if (fmag != nullptr && r > static_cast<Tcalc>(constants::tiny)) {
      *fmag += qiqj * ((value_two * aparm) + (value_two * aparm / r));
    }
  }
  else {
    const Tcalc invr = value_one / r;
    const Tcalc invr2 = invr * invr;
    *ele_contrib += qiqj * invr;
    if (fmag != nullptr) {
      *fmag += -(qiqj * invr * invr2);
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
void quarticCoreLennardJones(const Tcalc r, const Tcalc	clash_ratio, const Tcalc lja,
                             const Tcalc ljb, Tcalc *vdw_contrib, Tcalc *fmag) {
  const Tcalc value_zero = 0.0;
  const Tcalc value_nil  = 1.0e-6;
  const Tcalc value_one  = 1.0;
  const bool tcalc_is_double = (std::type_index(typeid(Tcalc)).hash_code() == double_type_index);
  Tcalc ij_sigma;
  if (ljb > value_nil) {
    ij_sigma = (tcalc_is_double) ? sqrt(cbrt(lja / ljb)) : sqrtf(cbrtf(lja / ljb));
  }
  else {
    ij_sigma = value_zero;
  }
  const Tcalc rlimit = clash_ratio * ij_sigma;
  if (r < rlimit) {
    const Tcalc invrlim = value_one / rlimit;
    const Tcalc invrlim2 = invrlim * invrlim;
    const Tcalc invrlim6 = invrlim2 * invrlim2 * invrlim2;
    Tcalc aparm;
    if (tcalc_is_double) {
      aparm = (((6.0 * ljb) - (12.0 * lja * invrlim6)) * invrlim * invrlim6) /
              ((((((4.0 * rlimit) + 12.0) * rlimit) + 12.0) * rlimit) + 4.0);
    }
    else {
      aparm = (((6.0f * ljb) - (12.0f * lja * invrlim6)) * invrlim * invrlim6) /
              ((((((4.0f * rlimit) + 12.0f) * rlimit) + 12.0f) * rlimit) + 4.0f);
    }
    const Tcalc rlimit_plus_one = rlimit + value_one;
    const Tcalc arlimit_po_four = aparm * (rlimit_plus_one * rlimit_plus_one *
                                           rlimit_plus_one * rlimit_plus_one);
    const Tcalc bparm = (((lja * invrlim6) - ljb) * invrlim6) - (arlimit_po_four);
    const Tcalc r_plus_one = r + value_one;
    const Tcalc arpo_three = aparm * r_plus_one * r_plus_one * r_plus_one;
    *vdw_contrib += (arpo_three * r_plus_one) + bparm;
    if (fmag != nullptr && r > static_cast<Tcalc>(constants::tiny)) {
      *fmag += (tcalc_is_double) ? (4.0 * arpo_three / r) : (4.0f * arpo_three / r);
    }
  }
  else {
    const Tcalc invr2 = value_one / (r * r);
    const Tcalc invr4 = invr2 * invr2;
    *vdw_contrib += (lja * invr4 * invr4 * invr4) - (ljb * invr4 * invr2);
    if (fmag != nullptr) {
      if (tcalc_is_double) {
        *fmag += ((6.0 * ljb) - (12.0 * lja * invr4 * invr2)) * invr4 * invr4;
      }
      else {
        *fmag += ((6.0f * ljb) - (12.0f * lja * invr4 * invr2)) * invr4 * invr4;
      }
    }
  }
}
  
} // namespace energy
} // namespace stormm
