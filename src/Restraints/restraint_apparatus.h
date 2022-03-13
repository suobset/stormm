// -*-c++-*-
#ifndef OMNI_RESTRAINT_APPARATUS_H
#define OMNI_RESTRAINT_APPARATUS_H

#include "Accelerator/hybrid.h"
#include "Topology/atomgraph.h"
#include "bounded_restraint.h"

namespace omni {
namespace restraints {

using card::Hybrid;
using card::HybridTargetLevel;

/// \brief Double-precision reader abstract for the RestraintApparatus class.  Restraints are very
///        detailed things, as reflected in the tedium of this object.
struct RestraintApparatusDpReader {

  /// \brief The constructor requires a tedious list of arguments.
  RestraintApparatusDpReader(int total_rst_in, int nposn_in, int nbond_in, int nangl_in,
                             int ndihe_in, bool time_dependence_in, const int* rposn_atoms_in,
                             const int* rbond_i_atoms_in, const int* rbond_j_atoms_in,
                             const int* rangl_i_atoms_in, const int* rangl_j_atoms_in,
                             const int* rangl_k_atoms_in, const int* rdihe_i_atoms_in,
                             const int* rdihe_j_atoms_in, const int* rdihe_k_atoms_in,
                             const int* rdihe_l_atoms_in, const int* rposn_init_step_in,
                             const int* rposn_finl_step_in, const int* rbond_init_step_in,
                             const int* rbond_finl_step_in, const int* rangl_init_step_in,
                             const int* rangl_finl_step_in, const int* rdihe_init_step_in,
                             const int* rdihe_finl_step_in, const double2* rposn_init_keq_in,
                             const double2* rposn_finl_keq_in, const double2* rposn_init_xy_in,
                             const double2* rposn_finl_xy_in, const double*  rposn_init_z_in,
                             const double*  rposn_finl_z_in, const double2* rbond_init_keq_in,
                             const double2* rbond_finl_keq_in, const double2* rangl_init_keq_in,
                             const double2* rangl_finl_keq_in, const double2* rdihe_init_keq_in,
                             const double2* rdihe_finl_keq_in, const double4* rposn_init_r_in,
                             const double4* rposn_finl_r_in, const double4* rbond_init_r_in,
                             const double4* rbond_finl_r_in, const double4* rangl_init_r_in,
                             const double4* rangl_finl_r_in, const double4* rdihe_init_r_in,
                             const double4* rdihe_finl_r_in, const AtomGraph *ag_pointer_in);

  /// \brief Take the default copy and move constructors.  The assignment operators will get
  ///        implicitly deleted as this is just a collection of constants.
  /// \{
  RestraintApparatusDpReader(const RestraintApparatusDpReader &original) = default;
  RestraintApparatusDpReader(RestraintApparatusDpReader &&other) = default;
  RestraintApparatusDpReader& operator=(const RestraintApparatusDpReader &original) = default;
  /// \}

  const int total_rst;            ///< Total number of restraints in the apparatus
  const int nposn;                ///< Number of positional restraints
  const int nbond;                ///< Number of distance restraints
  const int nangl;                ///< Number of three-point restraints
  const int ndihe;                ///< Number of four-point dihedral restraints
  const bool time_dependence;     ///< Whether there is a time dependent nature to the restraints
  const int* rposn_atoms;         ///< Atoms which are under positional restraints
  const int* rbond_i_atoms;       ///< Atoms which are under distance restraints (i)
  const int* rbond_j_atoms;       ///< Atoms which are under distance restraints (ii)
  const int* rangl_i_atoms;       ///< Atoms which are under three-point angle restraints (i)
  const int* rangl_j_atoms;       ///< Atoms which are under three-point angle restraints (ii)
  const int* rangl_k_atoms;       ///< Atoms which are under three-point angle restraints (iii)
  const int* rdihe_i_atoms;       ///< Atoms which are under four-point dihedral restraints (i)
  const int* rdihe_j_atoms;       ///< Atoms which are under four-point dihedral restraints (ii)
  const int* rdihe_k_atoms;       ///< Atoms which are under four-point dihedral restraints (iii)
  const int* rdihe_l_atoms;       ///< Atoms which are under four-point dihedral restraints (iv)
  const int* rposn_init_step;     ///< Initial steps for applying positional restraints
  const int* rposn_finl_step;     ///< Final steps for applying positional restraints
  const int* rbond_init_step;     ///< Initial steps for applying distance restraints
  const int* rbond_finl_step;     ///< Final steps for applying distance restraints
  const int* rangl_init_step;     ///< Initial steps for applying three-point angle restraints
  const int* rangl_finl_step;     ///< Final steps for applying three-point angle restraints
  const int* rdihe_init_step;     ///< Initial steps for applying dihedral restraints
  const int* rdihe_finl_step;     ///< Final steps for applying dihedral restraints
  const double2* rposn_init_keq;  ///< Initial harmonic parameters for positional restraints
  const double2* rposn_finl_keq;  ///< Final harmonic parameters for positional restraints
  const double2* rposn_init_xy;   ///< Initial X and Y coordinates for positional restraint targets
  const double2* rposn_finl_xy;   ///< Final X and Y coordinates for positional restraint targets
  const double*  rposn_init_z;    ///< Initial Z coordinates for positional restraint targets
  const double*  rposn_finl_z;    ///< Final Z coordinates for positional restraint targets
  const double2* rbond_init_keq;  ///< Initial harmonic parameters for distance restraints
  const double2* rbond_finl_keq;  ///< Final harmonic parameters for distance restraints
  const double2* rangl_init_keq;  ///< Initial harmonic parameters for angle restraints
  const double2* rangl_finl_keq;  ///< Final harmonic parameters for angle restraints
  const double2* rdihe_init_keq;  ///< Initial harmonic parameters for dihedral restraints
  const double2* rdihe_finl_keq;  ///< Final harmonic parameters for dihedral restraints
  const double4* rposn_init_r;    ///< Initial displacement parameters for positional restraints
  const double4* rposn_finl_r;    ///< Final displacement parameters for positional restraints
  const double4* rbond_init_r;    ///< Initial displacement parameters for distance restraints
  const double4* rbond_finl_r;    ///< Final displacement parameters for distance restraints
  const double4* rangl_init_r;    ///< Initial displacement parameters for angle restraints
  const double4* rangl_finl_r;    ///< Final displacement parameters for angle restraints
  const double4* rdihe_init_r;    ///< Initial displacement parameters for dihedral restraints
  const double4* rdihe_finl_r;    ///< Final displacement parameters for dihedral restraints
  const AtomGraph *ag_pointer;    ///< Pointer to the topology to which this apparatus applies
};

/// \brief Single-precision reader abstract for the RestraintApparatus class.  Restraints are very
///        detailed things, as reflected in the tedium of this object.
struct RestraintApparatusSpReader {

  /// \brief The constructor requires a tedious list of arguments.
  RestraintApparatusSpReader(int total_rst_in, int nposn_in, int nbond_in, int nangl_in,
                             int ndihe_in, bool time_dependence_in, const int* rposn_atoms_in,
                             const int* rbond_i_atoms_in, const int* rbond_j_atoms_in,
                             const int* rangl_i_atoms_in, const int* rangl_j_atoms_in,
                             const int* rangl_k_atoms_in, const int* rdihe_i_atoms_in,
                             const int* rdihe_j_atoms_in, const int* rdihe_k_atoms_in,
                             const int* rdihe_l_atoms_in, const int* rposn_init_step_in,
                             const int* rposn_finl_step_in, const int* rbond_init_step_in,
                             const int* rbond_finl_step_in, const int* rangl_init_step_in,
                             const int* rangl_finl_step_in, const int* rdihe_init_step_in,
                             const int* rdihe_finl_step_in, const float2* rposn_init_keq_in,
                             const float2* rposn_finl_keq_in, const float2* rposn_init_xy_in,
                             const float2* rposn_finl_xy_in, const float*  rposn_init_z_in,
                             const float*  rposn_finl_z_in, const float2* rbond_init_keq_in,
                             const float2* rbond_finl_keq_in, const float2* rangl_init_keq_in,
                             const float2* rangl_finl_keq_in, const float2* rdihe_init_keq_in,
                             const float2* rdihe_finl_keq_in, const float4* rposn_init_r_in,
                             const float4* rposn_finl_r_in, const float4* rbond_init_r_in,
                             const float4* rbond_finl_r_in, const float4* rangl_init_r_in,
                             const float4* rangl_finl_r_in, const float4* rdihe_init_r_in,
                             const float4* rdihe_finl_r_in, const AtomGraph *ag_pointer_in);

  /// \brief Take the default copy and move constructors.  The assignment operators will get
  ///        implicitly deleted as this is just a collection of constants.
  /// \{
  RestraintApparatusSpReader(const RestraintApparatusSpReader &original) = default;
  RestraintApparatusSpReader(RestraintApparatusSpReader &&other) = default;
  RestraintApparatusSpReader& operator=(const RestraintApparatusSpReader &original) = default;
  /// \}

  const int total_rst;           ///< Total number of restraints in the apparatus
  const int nposn;               ///< Number of positional restraints
  const int nbond;               ///< Number of distance restraints
  const int nangl;               ///< Number of three-point restraints
  const int ndihe;               ///< Number of four-point dihedral restraints
  const bool time_dependence;    ///< Whether there is a time dependent nature to the restraints
  const int* rposn_atoms;        ///< Atoms which are under positional restraints
  const int* rbond_i_atoms;      ///< Atoms which are under distance restraints (i)
  const int* rbond_j_atoms;      ///< Atoms which are under distance restraints (ii)
  const int* rangl_i_atoms;      ///< Atoms which are under three-point angle restraints (i)
  const int* rangl_j_atoms;      ///< Atoms which are under three-point angle restraints (ii)
  const int* rangl_k_atoms;      ///< Atoms which are under three-point angle restraints (iii)
  const int* rdihe_i_atoms;      ///< Atoms which are under four-point dihedral restraints (i)
  const int* rdihe_j_atoms;      ///< Atoms which are under four-point dihedral restraints (ii)
  const int* rdihe_k_atoms;      ///< Atoms which are under four-point dihedral restraints (iii)
  const int* rdihe_l_atoms;      ///< Atoms which are under four-point dihedral restraints (iv)
  const int* rposn_init_step;    ///< Initial steps for applying positional restraints
  const int* rposn_finl_step;    ///< Final steps for applying positional restraints
  const int* rbond_init_step;    ///< Initial steps for applying distance restraints
  const int* rbond_finl_step;    ///< Final steps for applying distance restraints
  const int* rangl_init_step;    ///< Initial steps for applying three-point angle restraints
  const int* rangl_finl_step;    ///< Final steps for applying three-point angle restraints
  const int* rdihe_init_step;    ///< Initial steps for applying dihedral restraints
  const int* rdihe_finl_step;    ///< Final steps for applying dihedral restraints
  const float2* rposn_init_keq;  ///< Initial harmonic parameters for positional restraints
  const float2* rposn_finl_keq;  ///< Final harmonic parameters for positional restraints
  const float2* rposn_init_xy;   ///< Initial X and Y coordinates for positional restraint targets
  const float2* rposn_finl_xy;   ///< Final X and Y coordinates for positional restraint targets
  const float*  rposn_init_z;    ///< Initial Z coordinates for positional restraint targets
  const float*  rposn_finl_z;    ///< Final Z coordinates for positional restraint targets
  const float2* rbond_init_keq;  ///< Initial harmonic parameters for distance restraints
  const float2* rbond_finl_keq;  ///< Final harmonic parameters for distance restraints
  const float2* rangl_init_keq;  ///< Initial harmonic parameters for angle restraints
  const float2* rangl_finl_keq;  ///< Final harmonic parameters for angle restraints
  const float2* rdihe_init_keq;  ///< Initial harmonic parameters for dihedral restraints
  const float2* rdihe_finl_keq;  ///< Final harmonic parameters for dihedral restraints
  const float4* rposn_init_r;    ///< Initial displacement parameters for positional restraints
  const float4* rposn_finl_r;    ///< Final displacement parameters for positional restraints
  const float4* rbond_init_r;    ///< Initial displacement parameters for distance restraints
  const float4* rbond_finl_r;    ///< Final displacement parameters for distance restraints
  const float4* rangl_init_r;    ///< Initial displacement parameters for angle restraints
  const float4* rangl_finl_r;    ///< Final displacement parameters for angle restraints
  const float4* rdihe_init_r;    ///< Initial displacement parameters for dihedral restraints
  const float4* rdihe_finl_r;    ///< Final displacement parameters for dihedral restraints
  const AtomGraph *ag_pointer;   ///< Pointer to the topology to which this apparatus applies
};

/// \brief A collection of all restraints pertaining to a specific topology for the purposes of
///        one simulation, energy minimization, or even a single molecular mechanics calculation.
class RestraintApparatus {
public:

  /// \brief The constructor takes a vector of individual restraints
  ///
  /// \param rbasis  A list of restraint objects with which to build the apparatus
  /// \param ag_in   Pointer to the topology for which this apparatus is built (if rbasis is of
  ///                nonzero length, the topology pointer from the first restraint will be
  ///                preferred, otherwise this must be supplied in order to initialize the const
  ///                member variable ag_pointer)
  RestraintApparatus(const std::vector<BoundedRestraint> &rbasis, const AtomGraph *ag_in = nullptr);

  /// \brief The copy constructor works like any other object containing POINTER-kind Hybrids.
  ///
  /// \param original  The object to copy
  RestraintApparatus(const RestraintApparatus &original);

  /// \brief The move constructor also works like other objects containing POINTER-kind Hybrids.
  ///
  /// \param original  The object to move
  RestraintApparatus(RestraintApparatus &&original);

  /// \brief Copy assignment operator
  ///
  /// \param other  The object to move
  RestraintApparatus& operator=(const RestraintApparatus &other);
  
  /// \brief Move assignment operator
  ///
  /// \param other  The object to move
  RestraintApparatus& operator=(RestraintApparatus &&other);

  /// \brief Get the total number of restraints in this apparatus
  int getTotalRestraintCount() const;

  /// \brief Get the number of positional restraints
  int getPositionalRestraintCount() const;
  
  /// \brief Get the number of distance restraints
  int getDistanceRestraintCount() const;
  
  /// \brief Get the number of angle restraints
  int getAngleRestraintCount() const;
  
  /// \brief Get the number of dihedral restraints
  int getDihedralRestraintCount() const;

  /// \brief Get an indication of whether this apparatus uses time-dependent restraints
  bool getTimeDependence() const;

  /// \brief Get a double-precision abstract of this apparatus
  RestraintApparatusDpReader dpData(HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  /// \brief Get a single-precision abstract of this apparatus
  RestraintApparatusSpReader spData(HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  /// \brief Produce a vector of all the restraints in this apparatus, essentially the inverse of
  ///        the constructor.
  std::vector<BoundedRestraint> getRestraintList() const;

  /// \brief Add restraints to the apparatus.
  ///
  /// Overloaded:
  ///   - Add a vector of new restraints
  ///   - Add a solitary new restraint
  ///
  /// \param new_rest  One or more new restraints
  /// \{
  void addRestraints(const std::vector<BoundedRestraint> &new_rest);
  void addRestraint(const BoundedRestraint &new_rest);
  /// \{
  
private:

  // Overall counts and items of interest
  int total_restraint_count;   ///< Total number of restraints in this apparatus
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

  /// \brief Check that the vector of restraints all apply to the same topology.
  ///
  /// \param rbasis  A list of restraint objects with which to build the apparatus
  void checkTopologyPointers(const std::vector<BoundedRestraint> &rbasis);
  
  /// \brief Allocate memory and set POINTER-kind Hybrid objects as appropriate.  This function
  ///        encapsulates actions needed by the base constructor as well as copy constructors.
  ///        Requires prior computation of the numbers of various types of restraints, as set the
  ///        member variables position_count, distance_count, angle_count, and dihedral_count.
  void allocate();

  /// \brief Dissect each restraint in the list provided and populate the RestraintApparatus's
  ///        internal data arrays.  Requires prior computation of the numbers of various types of
  ///        restraints, as set the member variables position_count, distance_count, angle_count,
  ///        and dihedral_count.
  ///
  /// \param rbasis  A list of restraint objects with which to build the apparatus
  void populateInternalArrays(const std::vector<BoundedRestraint> &rbasis);
};

} // namespace restraints
} // namespace omni

#endif
