// -*-c++-*-
#ifndef OMNI_NMR_RESTRAINT_H
#define OMNI_NMR_RESTRAINT_H

#include <string>
#include "Chemistry/chemical_features.h"
#include "DataTypes/omni_vector_types.h"
#include "Trajectory/coordinateframe.h"
#include "Topology/atomgraph.h"

namespace omni {
namespace restraints {

using chemistry::ChemicalFeatures;
using topology::AtomGraph;
using trajectory::CoordinateFrameReader;

struct BoundedRestraint {

  /// Constructors take either four atom masks (each of whcih must evaluate to exactly one atom)
  /// or four atom numbers (numbers are given for a series 1... atom_count in the topology pointer,
  /// decremented when constructing the mask to index into the actual memory)
  /// \{
  BoundedRestraint(const std::string &mask_i_in, const std::string &mask_j_in, 
                   const std::string &mask_k_in, const std::string &mask_l_in,
                   const AtomGraph *ag_in, const ChemicalFeatures *chemfe,
                   const CoordinateFrameReader &cfr, int init_step_in, int final_step_in,
                   double init_k2_in, double init_k3_in, double init_r1_in, double init_r2_in,
                   double init_r3_in, double init_r4_in, double final_k2_in, double final_k3_in,
                   double final_r1_in, double final_r2_in, double final_r3_in, double final_r4_in,
                   double init_ref_x_crd_in, double init_ref_y_crd_in, double init_ref_z_crd_in,
                   double final_ref_x_crd_in, double final_ref_y_crd_in,
                   double final_ref_z_crd_in);

  BoundedRestraint(int atom_i_in, int atom_j_in, int atom_k_in, int atom_l_in,
                   const AtomGraph *ag_in, int init_step_in, int final_step_in, double init_k2_in,
                   double init_k3_in, double init_r1_in, double init_r2_in, double init_r3_in,
                   double init_r4_in, double final_k2_in, double final_k3_in, double final_r1_in,
                   double final_r2_in, double final_r3_in, double final_r4_in,
                   double init_ref_x_crd_in, double init_ref_y_crd_in, double init_ref_z_crd_in,
                   double final_ref_x_crd_in, double final_ref_y_crd_in,
                   double final_ref_z_crd_in);

  BoundedRestraint(int atom_index, int refr_index, const AtomGraph *ag_in,
                   const CoordinateFrameReader &cfr, double k2_in, double k3_in, double r1_in,
                   double r2_in, double r3_in, double r4_in);
  /// \}

  /// Obtain one of the masks used to specify an atom in this restraint
  ///
  /// \param restrained_atom_number  The 1st, 2nd, 3rd, or 4th atom (specify 1, 2, 3, or 4)
  std::string getAtomMask(int restrained_atom_number) const;

  /// Obtain the topology index of an atom in this restraint
  ///
  /// \param restrained_atom_number  The 1st, 2nd, 3rd, or 4th atom (specify 1, 2, 3, or 4)
  int getAtomIndex(int restrained_atom_number) const;

  /// \brief Get the order of this restraint
  int getOrder() const;

  /// \brief Get the initial step at which to begin applying this restraint
  int getInitialStep() const;

  /// \brief Get the simulation step at which to finish applying this restraint
  int getFinalStep() const;

  /// \brief Get the initial stiffnesses to use when applying this restraint
  double2 getInitialStiffness() const;

  /// \brief Get the final stiffnesses of the restraint in its complete form
  double2 getFinalStiffness() const;

  /// \brief Get the initial displacement parameters to use in applying this restraint
  double4 getInitialDisplacements() const;

  /// \brief Get the final displacement parameters of the restraint in its complete form
  double4 getFinalDisplacements() const;

  /// \brief Get the initial target of a positional restraint
  double3 getInitialTargetSite() const;

  /// \brief Get the final target of a positional restraint
  double3 getFinalTargetSite() const;

  /// Get the topology pointer
  const AtomGraph* getTopologyPointer() const;

private:
  std::string mask_i;   ///< Atom mask specifying atom I (must determine exactly one atom)
  std::string mask_j;   ///< Atom mask specifying atom J (must determine exactly one atom)
  std::string mask_k;   ///< Atom mask specifying atom K (must determine exactly one atom)
  std::string mask_l;   ///< Atom mask specifying atom L (must determine exactly one atom)
  int atom_i;           ///< Index of atom I in the corresponding topology
  int atom_j;           ///< Index of atom J in the corresponding topology
  int atom_k;           ///< Index of atom K in the corresponding topology
  int atom_l;           ///< Index of atom L in the corresponding topology
  int order;            ///< Order of the NMR restraint: 2 = distance restraint between atoms I and
                        ///<   J (implied if K and L are unspecified), 3 = angle restraint between
                        ///<   atoms I, J, and K (implied if L is unspecified), 4 = dihedral
                        ///<   restraint between atoms I, J, K, and L
  int initial_step;     ///< Initial step of the simulation at which to begin applying the
                        ///<   restraint with initial keq and r parameters
  int final_step;       ///< Step of the simulation at which to finish applying the full restraint,
                        ///<   with final kew and r parameters
  double2 initial_keq;  ///< Initial stiffness constants for parabolic restraints between points
                        ///<   r1 and r2 (x member of the tuple) and points r3 and r4 (y member of
                        ///<   the tuple)
  double4 initial_r;    ///< Initial displacement parameters r1 (x), r2 (y), r3 (z), and r4 (w)
  double2 final_keq;    ///< Final stiffness constants
  double4 final_r;      ///< Final displacement parameters
  double3 init_center;  ///< Initial center of the restraint potential (positional restraints only)
  double3 final_center; ///< Final center of the restraint potential (positional restraints only)
  
  /// Pointer to the topology for which this restraint applies
  const AtomGraph *ag_pointer;  
};

} // namespace restraints
} // namespace omni

#endif
