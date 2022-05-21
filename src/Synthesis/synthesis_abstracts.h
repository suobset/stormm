// -*-c++-*-
#ifndef OMNI_SYNTHESIS_ABSTRACTS_H
#define OMNI_SYNTHESIS_ABSTRACTS_H

#include "DataTypes/common_types.h"
#include "DataTypes/omni_vector_types.h"

namespace omni {
namespace synthesis {

/// \brief Collect the critical valence parameters and indexing information for work unit-based
///        evaluation of the systems in an AtomGraphSynthesis.
template <typename T> struct SyValenceKit {

  /// \brief The constructor takes a straight list of arguments for each member variable.
  explicit SyValenceKit(int nvwu_in, const T coulomb_in, const T* bond_keq_in,
                        const T* bond_leq_in, const T* angl_keq_in, const T* angl_theta_in,
                        const T* dihe_amp_in, const T* dihe_freq_in, const T* dihe_phi_in,
                        const T* attn14_elec_in, const T* attn14_vdw_in, const T* charges_in,
                        const T* lja_14_coeff_in, const T* ljb_14_coeff_in,
                        const T* ljc_14_coeff_in, const int* lj_idx_in, const int* n_lj_types_in,
                        const int* ljabc_offsets_in, const T* ubrd_keq_in, const T* ubrd_leq_in,
                        const T* cimp_keq_in, const T* cimp_phi_in, const int* cmap_dim_in,
                        const T* cmap_patches_in, const int* cmap_patch_bounds_in,
                        const int2* vwu_abstracts_in, const int* vwu_imports_in,
                        const uint2* vwu_manip_in, const uint2* cbnd_insr_in,
                        const uint2* angl_insr_in, const uint2* cdhe_insr_in,
                        const uint* cdhe_ovrt_insr_in, const uint2* cmap_insr_in,
                        const uint* infr14_insr_in, const uint* cbnd_acc_in,
                        const uint* angl_acc_in, const uint* cdhe_acc_in, const uint* cmap_acc_in,
                        const uint* infr14_acc_in);

  // Member variables (all public)
  const int nvwu;                ///< The number of valence work units
  const T coulomb;               ///< Coulomb's constant in kcal-A2/mol-e2
  const T* bond_keq;             ///< Consensus table of harmonic bond stiffness constants
  const T* bond_leq;             ///< Consensus table of harmonic bond equilibrium lengths
  const T* angl_keq;             ///< Consensus table of harmonic angle stiffness constants
  const T* angl_theta;           ///< Consensus table of harmonic angle equilibrium settings
  const T* dihe_amp;             ///< Consensus table of cosine-based dihedral amplitudes
  const T* dihe_freq;            ///< Consensus table of cosine-based dihedral frequencies 
  const T* dihe_phi;             ///< Consensus table of cosine-based dihedral phase angles
  const T* attn14_elec;          ///< Consensus table of 1:4 electrostatic scaling factors
  const T* attn14_vdw;           ///< Consensus table of 1:4 van-der Waals scaling factors
  const T* charges;              ///< Charges on all atoms
  const T* lja_14_coeff;         ///< Lennard-Jones 1:4 interaction A coefficients
  const T* ljb_14_coeff;         ///< Lennard-Jones 1:4 interaction B coefficients
  const T* ljc_14_coeff;         ///< Lennard-Jones 1:4 interaction C coefficients
  const int* lj_idx;             ///< Lennard-Jones type indices for all atoms
  const int* n_lj_types;         ///< Lennard-Jones type counts for all systems
  const int* ljabc_offsets;      ///< Offsets for Lennard-Jones A, B, and C coefficient tables for
                                 ///<   all systems
  const T* ubrd_keq;             ///< Consensus table of Urey-Bradley stiffnesses
  const T* ubrd_leq;             ///< Consensus table of Urey-Bradley I:K separation distances
  const T* cimp_keq;             ///< Consensus table of CHARMM improper dihedral stiffness
                                 ///<   constants
  const T* cimp_phi;             ///< Consensus table of CHARMM improper dihedral phase angles
  const int* cmap_dim;           ///< Dimensions for each member of the consensus table of CMAPs
  const T* cmap_patches;         ///< Consensus table of CMAP pre-computed surface patch data
  const int* cmap_patch_bounds;  ///< Bounds array for cmap_patches
  const int2* vwu_abstracts;     ///< Abstracts for all valence work units
  const int* vwu_imports;        ///< Imported atom tables for all valence work units
  const uint2* vwu_manip;        ///< Manipulation masks (movement in x member, update in y member)
                                 ///<   for all valence work units
  const uint2* cbnd_insr;        ///< Composite bond and Urey-Bradley instructions implemented by
                                 ///<   valence work units
  const uint2* angl_insr;        ///< Harmonic angle instructions implemented by valence work units
  const uint2* cdhe_insr;        ///< Composite cosine-based and CHARMM improper dihedral
                                 ///<   instructions implemented by valence work units
  const uint* cdhe_ovrt_insr;    ///< Instructions for composite cosine-based dihedral overtones
  const uint2* cmap_insr;        ///< CMAP instructions implemented by valence work units
  const uint* infr14_insr;       ///< Inferred 1:4 attenuated interaction instructions implemented
                                 ///<   by valence work units
  const uint* cbnd_acc;          ///< Composite bond energy accumulation masks
  const uint* angl_acc;          ///< Harmonic angle energy accumulation masks
  const uint* cdhe_acc;          ///< Composite dihedral energy accumulation masks
  const uint* cmap_acc;          ///< CMAP energy accumulation masks
  const uint* infr14_acc;        ///< Inferred 1:4 interaction energy accumulation masks
};

/// \brief Collect the critical non-bonded parameters and masking information for work unit-based
///        evaluation of the systems in an AtomGraphSynthesis.
template <typename T> struct SyNonbondedKit {

  /// \brief The constructor takes a straight list of arguments for each member variable.
  explicit SyNonbondedKit();
};
  
} // namespace synthesis
} // namespace omni

#include "synthesis_abstracts.tpp"

#endif

