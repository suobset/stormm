// -*-c++-*-
#ifndef OMNI_RESTRAINT_APPARATUS_H
#define OMNI_RESTRAINT_APPARATUS_H

#include "bounded_restraint.h"

namespace omni {
namespace restraints {

/// \brief A collection of all restraints pertaining to a specific topology for the purposes of
///        one simulation, energy minimization, or even a single molecular mechanics calculation.
struct RestraintApparatus {

  /// \brief The constructor takes a vector of individual restraints
  RestraintApparatus(const std::vector<BoundedRestraint> &rbasis);

private:
  int restraint_count;
  int distance_count;
  int angle_count;
  int dihedral_count;
  Hybrid<int> rbond_i_atoms;
  Hybrid<int> rbond_j_atoms;
  Hybrid<int> rangl_i_atoms;
  Hybrid<int> rangl_j_atoms;
  Hybrid<int> rangl_k_atoms;
  Hybrid<int> rdihe_i_atoms;
  Hybrid<int> rdihe_j_atoms;
  Hybrid<int> rdihe_k_atoms;
  Hybrid<int> rdihe_l_atoms;
  
  Hybrid<int> int_data;
};

} // namespace restraints
} // namespace omni

#endif
