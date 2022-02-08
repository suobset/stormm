// -*-c++-*-
#ifndef OMNI_RESTRAINT_APPARATUS_H
#define OMNI_RESTRAINT_APPARATUS_H

#include "Accelerator/hybrid.h"
#include "bounded_restraint.h"

namespace omni {
namespace restraints {

using card::Hybrid;

/// \brief A collection of all restraints pertaining to a specific topology for the purposes of
///        one simulation, energy minimization, or even a single molecular mechanics calculation.
struct RestraintApparatus {

  /// \brief The constructor takes a vector of individual restraints
  RestraintApparatus(const std::vector<BoundedRestraint> &rbasis);

private:
  int restraint_count;        ///< Total number of restraints in this apparatus
  int position_count;         ///< Number of positional restraints
  int distance_count;         ///< Number of distance restraints
  int angle_count;            ///< Number of angle (three-point) restraints
  int dihedral_count;         ///< Number of dihedral (four-point) restraints
  Hybrid<int> rposn_i_atoms;  ///< Topological indices of positionally restrained atoms
  Hybrid<int> rbond_i_atoms;  ///< Topological indices of distance restraint i atoms
  Hybrid<int> rbond_j_atoms;  ///< Topological indices of distance restraint j atoms
  Hybrid<int> rangl_i_atoms;  ///< Topological indices of angle restraint i atoms
  Hybrid<int> rangl_j_atoms;  ///< Topological indices of angle restraint j atoms
  Hybrid<int> rangl_k_atoms;  ///< Topological indices of angle restraint k atoms
  Hybrid<int> rdihe_i_atoms;  ///< Topological indices of dihedral (four-point) restraint i atoms
  Hybrid<int> rdihe_j_atoms;  ///< Topological indices of dihedral (four-point) restraint j atoms
  Hybrid<int> rdihe_k_atoms;  ///< Topological indices of dihedral (four-point) restraint k atoms
  Hybrid<int> rdihe_l_atoms;  ///< Topological indices of dihedral (four-point) restraint l atoms
  Hybrid<int> int_data;
};

} // namespace restraints
} // namespace omni

#endif
