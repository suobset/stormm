// -*-c++-*-
#ifndef OMNI_RESTRAINT_APPARATUS_H
#define OMNI_RESTRAINT_APPARATUS_H

#include "Accelerator/hybrid.h"
#include "Topology/atomgraph.h"
#include "bounded_restraint.h"

namespace omni {
namespace restraints {

using card::Hybrid;

/// \brief A collection of all restraints pertaining to a specific topology for the purposes of
///        one simulation, energy minimization, or even a single molecular mechanics calculation.
struct RestraintApparatus {

  /// \brief The constructor takes a vector of individual restraints
  RestraintApparatus(const std::vector<BoundedRestraint> &rbasis);

  /// \brief The copy constructor works like any other object containing POINTER-kind Hybrids.
  RestraintApparatus(const RestraintApparatus &original);

  /// \brief The move constructor also works like other objects containining POINTER-kind Hybrids.
  RestraintApparatus(RestraintApparatus &&original);
  
  
private:

  // Overall counts and items of interest
  int restraint_count;         ///< Total number of restraints in this apparatus
  int position_count;          ///< Number of positional restraints
  int distance_count;          ///< Number of distance restraints
  int angle_count;             ///< Number of angle (three-point) restraints
  int dihedral_count;          ///< Number of dihedral (four-point) restraints
  bool time_based_restraints;  ///< Flag to indicate whether restraints are at all dependent on
                               ///<   the time step of the simulation (for energy minimizations or
                               ///<   simple molecular mechanics calculations, restraints are
                               ///<   are constant)

  // Integer data 
  Hybrid<int> rposn_atoms;       ///< Topological indices of positionally restrained atoms
  Hybrid<int> rbond_i_atoms;     ///< Topological indices of distance restraint I atoms
  Hybrid<int> rbond_j_atoms;     ///< Topological indices of distance restraint J atoms
  Hybrid<int> rangl_i_atoms;     ///< Topological indices of angle restraint I atoms
  Hybrid<int> rangl_j_atoms;     ///< Topological indices of angle restraint J atoms
  Hybrid<int> rangl_k_atoms;     ///< Topological indices of angle restraint K atoms
  Hybrid<int> rdihe_i_atoms;     ///< Topological indices of dihedral restraint I atoms
  Hybrid<int> rdihe_j_atoms;     ///< Topological indices of dihedral restraint J atoms
  Hybrid<int> rdihe_k_atoms;     ///< Topological indices of dihedral restraint K atoms
  Hybrid<int> rdihe_l_atoms;     ///< Topological indices of dihedral restraint L atoms
  Hybrid<int> rposn_init_step;   ///< Initial step numbers for applying positional restraints
  Hybrid<int> rposn_final_step;  ///< Final step numbers for applying positional restraints
  Hybrid<int> rbond_init_step;   ///< Initial step numbers for applying distance restraints
  Hybrid<int> rbond_final_step;  ///< Final step numbers for applying distance restraints
  Hybrid<int> rangl_init_step;   ///< Initial step numbers for applying angle restraints
  Hybrid<int> rangl_final_step;  ///< Final step numbers for applying angle restraints
  Hybrid<int> rdihe_init_step;   ///< Initial step numbers for applying dihedral restraints
  Hybrid<int> rdihe_final_step;  ///< Final step numbers for applying dihedral restraints
  Hybrid<int> int_data;          ///< Storage space for all integer data in this apparatus

  // Real data in double-precision format
  Hybrid<double2> rposn_init_keq;  ///< Initial stiffnesses for time-dependent positional
                                   ///<   restraints, or the static values of time-independent
                                   ///<   restraints 
  Hybrid<double2> rposn_final_keq; ///< Final stiffnesses for time-dependent positional restraints
                                   ///<   (ignored for time-independent restraints)
  Hybrid<double4> rposn_init_r;    ///< Initial displacements for time-dependent positional
                                   ///<   restraints, or the static values of time-independent
                                   ///<   restraints
  Hybrid<double4> rposn_final_r;   ///< Final displacments for time-dependent positional restraints
                                   ///<   (ignored for time-independent restraints)  
  Hybrid<double2> rposn_init_xy;   ///< Initial X and Y Cartesian coordinates for the target
                                   ///<   location of time-dependent positional restraints, or the
                                   ///<   static values of time-independent restraints
  Hybrid<double> rposn_init_z;     ///< Initial Z Cartesian coordinates for the target
                                   ///<   location of time-dependent positional restraints, or the
                                   ///<   static values of time-independent restraints
  Hybrid<double2> rposn_final_xy;  ///< Final X and Y Cartesian coordinates for the target location
                                   ///<   of time-dependent positional restraints, or the static
                                   ///<   values of time-independent restraints
  Hybrid<double> rposn_final_z;    ///< Final Z Cartesian coordinates for the target location of
                                   ///<   time-dependent positional restraints, or the static
                                   ///<   values of time-independent restraints
  Hybrid<double2> rbond_init_keq;  ///< Initial stiffnesses for time-dependent distance restraints,
                                   ///<   or the static values of time-independent restraints 
  Hybrid<double2> rbond_final_keq; ///< Final stiffnesses for time-dependent distance restraints
                                   ///<   (ignored for time-independent restraints)
  Hybrid<double4> rbond_init_r;    ///< Initial displacements for time-dependent distance
                                   ///<   restraints, or the static values of time-independent
                                   ///<   restraints
  Hybrid<double4> rbond_final_r;   ///< Final displacments for time-dependent distance restraints
                                   ///<   (ignored for time-independent restraints)  
  Hybrid<double2> rangl_init_keq;  ///< Initial stiffnesses for time-dependent angle restraints, or
                                   ///<   the static values of time-independent restraints 
  Hybrid<double2> rangl_final_keq; ///< Final stiffnesses for time-dependent angle restraints
                                   ///<   (ignored for time-independent restraints)
  Hybrid<double4> rangl_init_r;    ///< Initial displacements for time-dependent angle restraints,
                                   ///<   or the static values of time-independent restraints
  Hybrid<double4> rangl_final_r;   ///< Final displacments for time-dependent angle restraints
                                   ///<   (ignored for time-independent restraints)  
  Hybrid<double2> rdihe_init_keq;  ///< Initial stiffnesses for time-dependent dihedral restraints,
                                   ///<   or the static values of time-independent restraints 
  Hybrid<double2> rdihe_final_keq; ///< Final stiffnesses for time-dependent dihedral restraints
                                   ///<   (ignored for time-independent restraints)
  Hybrid<double4> rdihe_init_r;    ///< Initial displacements for time-dependent dihedral
                                   ///<   restraints, or the static values of time-independent
                                   ///<   restraints
  Hybrid<double4> rdihe_final_r;   ///< Final displacments for time-dependent dihedral restraints
                                   ///<   (ignored for time-independent restraints)
  Hybrid<double> double_data;      ///< Storage space for Cartesian Z target coordinates
  Hybrid<double2> double2_data;    ///< Storage space for double-precision stiffness constants and
                                   ///<   X/Y target coordinates
  Hybrid<double4> double4_data;    ///< Storage space for double-precision displacement values

  // Real data in single-precision format
  Hybrid<float2> sp_rposn_init_keq;  ///< Initial stiffnesses for time-dependent positional
                                     ///<   restraints, or the static values of time-independent
                                     ///<   restraints 
  Hybrid<float2> sp_rposn_final_keq; ///< Final stiffnesses for time-dependent positional
                                     ///<   restraints (ignored for time-independent restraints)
  Hybrid<float4> sp_rposn_init_r;    ///< Initial displacements for time-dependent positional
                                     ///<   restraints, or the static values of time-independent
                                     ///<   restraints
  Hybrid<float4> sp_rposn_final_r;   ///< Final displacments for time-dependent positional
                                     ///<   restraints (ignored for time-independent restraints)  
  Hybrid<float2> sp_rposn_init_xy;   ///< Initial X and Y Cartesian coordinates for the target
                                     ///<   location of time-dependent positional restraints, or
                                     ///<   the static values of time-independent restraints
  Hybrid<float> sp_rposn_init_z;     ///< Initial Z Cartesian coordinates for the target
                                     ///<   location of time-dependent positional restraints, or
                                     ///<   the static values of time-independent restraints
  Hybrid<float2> sp_rposn_final_xy;  ///< Final X and Y Cartesian coordinates for the target
                                     ///<   location of time-dependent positional restraints, or
                                     ///<   the static values of time-independent restraints
  Hybrid<float> sp_rposn_final_z;    ///< Final Z Cartesian coordinates for the target location of
                                     ///<   time-dependent positional restraints, or the static
                                     ///<   values of time-independent restraints
  Hybrid<float2> sp_rbond_init_keq;  ///< Initial stiffnesses for time-dependent distance
                                     ///<   restraints, or the static values of time-independent
                                     ///<   restraints 
  Hybrid<float2> sp_rbond_final_keq; ///< Final stiffnesses for time-dependent distance restraints
                                     ///<   (ignored for time-independent restraints)
  Hybrid<float4> sp_rbond_init_r;    ///< Initial displacements for time-dependent distance
                                     ///<   restraints, or the static values of time-independent
                                     ///<   restraints
  Hybrid<float4> sp_rbond_final_r;   ///< Final displacments for time-dependent distance restraints
                                     ///<   (ignored for time-independent restraints)  
  Hybrid<float2> sp_rangl_init_keq;  ///< Initial stiffnesses for time-dependent angle restraints,
                                     ///<   or the static values of time-independent restraints 
  Hybrid<float2> sp_rangl_final_keq; ///< Final stiffnesses for time-dependent angle restraints
                                     ///<   (ignored for time-independent restraints)
  Hybrid<float4> sp_rangl_init_r;    ///< Initial displacements for time-dependent angle
                                     ///<   restraints, or the static values of time-independent
                                     ///<   restraints
  Hybrid<float4> sp_rangl_final_r;   ///< Final displacments for time-dependent angle restraints
                                     ///<   (ignored for time-independent restraints)  
  Hybrid<float2> sp_rdihe_init_keq;  ///< Initial stiffnesses for time-dependent dihedral
                                     ///<   restraints, or the static values of time-independent
                                     ///<   restraints 
  Hybrid<float2> sp_rdihe_final_keq; ///< Final stiffnesses for time-dependent dihedral restraints
                                     ///<   (ignored for time-independent restraints)
  Hybrid<float4> sp_rdihe_init_r;    ///< Initial displacements for time-dependent dihedral
                                     ///<   restraints, or the static values of time-independent
                                     ///<   restraints
  Hybrid<float4> sp_rdihe_final_r;   ///< Final displacments for time-dependent dihedral restraints
                                     ///<   (ignored for time-independent restraints)
  Hybrid<float> float_data;          ///< Storage space for Cartesian Z target coordinates
  Hybrid<float2> float2_data;        ///< Storage space for single-precision stiffness constants
                                     ///<   and X/Y target coordinates
  Hybrid<float4> float4_data;        ///< Storage space for single-precision displacement values

  /// Pointer to the original topology (taken from one of the individual restraints, then used to
  /// check that all individual restraints point to the same topology)
  const AtomGraph *ag_pointer;

  /// \brief Allocate memory and set POINTER-kind Hybrid objects as appropriate.  This function
  ///        encapsulates actions needed by the base constructor as well as copy constructors.
  void allocate();
};

} // namespace restraints
} // namespace omni

#endif
