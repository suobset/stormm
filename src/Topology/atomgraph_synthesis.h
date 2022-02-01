// -*-c++-*-
#ifndef OMNI_ATOMGRAPH_SYNTHESIS_H
#define OMNI_ATOMGRAPH_SYNTHESIS_H

/// \brief A collection of one or more AtomGraph objects, with similar components arranged in
///        contiguous arrays (often padded by the GPU warp size to prevent one system from flowing
///        into another).  Work units covering all systems are laid out in additional, contiguous
///        arrays.  While individual topologies (AtomGraphs) have Hybrid POINTER-kind objects for
///        details such as charges and atom type indices and all of the POINTER-kind objects
///        targeted a single ARRAY-kind object of the correct memory type, the synthesis, which may
///        have many topologies and a great deal more overall information, stores most of its data
///        in a series of ARRAY-kind objects, one for each member variable.
struct AtomGraphSynthesis {

  /// The constructor takes a series of topologies and NMR restraints.  The NMR restraints point
  /// to specific topologies and thereby apply to any coordinate sets that also point to those
  /// topologies.
  AtomGraphSynthesis(const std::vector<AtomGraph*> &topologies_in,
                     const std::vector<int> topology_indices_in);

private:

  // The first member variables pertain to totals across all systems: atoms, potential function
  // terms, and sizes of composite parameter arrays that span all topologies
  int topology_count;       ///< The number of unique topologies represented in this synthesis
  int system_count;         ///< Number of independent coordinate sets present in this synthesis
  int total_atoms;          ///< The total number of atoms, spanning all topologies
  int total_virtual_sites;  ///< The total number of atoms, spanning all topologies
  int total_bond_terms;     ///< The total number of bonds, spanning all topologies
  int total_angl_terms;     ///< The total number of bond angles, spanning all topologies
  int total_dihe_terms;     ///< The total number of dihedrals, spanning all topologies
  int total_ubrd_terms;     ///< The total number of Urey-Bradley angles, spanning all topologies
  int total_cimp_terms;     ///< The total number of CHARMM impropers, spanning all topologies
  int total_cmap_terms;     ///< The total number of CMAPs, spanning all topologies
  int total_atom_types;     ///< Total number of unique atom types, spanning all topologies
  int total_charge_types;   ///< Total number of unique atomic partial charges, spanning all
                            ///<   topologies
  int total_bond_params;    ///< Total number of unique bond parameter sets
  int total_angl_params;    ///< Total number of unique bond angle parameter sets
  int total_dihe_params;    ///< Total number of unique dihedral parameter sets
  int total_ubrd_params;    ///< Total number of unique Urey-Bradley angle parameter sets
  int total_cimp_params;    ///< Total number of unique CHARMM improper parameter sets
  int total_cmap_surfaces;  ///< Total number of unique CMAP surfaces

  // The presence of periodic boundaries or an implicit solvent must be common to all systems, as
  // must the cutoff and other run parameters, although these are contained in out data structres.
  // If an implicit solvent model is in effect, its parameters are taken as common to all systems.
  // Other facets of the the system behavior, such as bond constraints, must also be common to
  // all systems.
  UnitCellType periodic_box_class;    ///< The type of unit cell; either NONE or any combination of
                                      ///< ORTHORHOMBIC and TRICLINIC is acceptable
  ImplicitSolventModel gb_style;      ///< The flavor of Generalized Born (or other implicit
                                      ///<   solvent) to use, i.e. Hawkins / Cramer / Truhlar.
  double dielectric_constant;         ///< Dielectric constant to take for implicit solvent
                                      ///<   calculations
  double salt_concentration;          ///< Salt concentration affecting implicit solvent models
  double coulomb_constant;            ///< Coulomb's constant in units of kcal-A/mol-e^2 (Amber
                                      ///<   differs from other programs in terms of what this is,
                                      ///<   so it can be set here)
  std::string pb_radii_set;           ///< The Poisson-Boltzmann radii set, also used in GB
  ShakeSetting use_bond_constraints;  ///< Toggles use of bond length constraints
  SettleSetting use_settle;           ///< Toggles analytic constraints on rigid water
  char4 water_residue_name;           ///< Name of water residue, compared to residue_names

  /// An array of pointers to the individual topologies that form the basis of this work plan
  /// (this array is topology_count in length and accessible only on the host as it is used to
  /// organize work units, not process the actual energy calculations)
  std::vector<AtomGraph*> topologies;

  /// An array of indices for the source topologies guiding the motion of each system (this array
  /// is system_count in length and accessible only on the host to organize work units)
  Hybrid<int> topology_indices;
  
  // The following arrays are POINTER-kind objects, each directed at a segment of a data array
  // for storing critical topology sizes.  There is one of each of these numbers for every
  // system in this work plan: each of the Hybrid arrays that follows is system_count in length.
  // All of the descriptions below follow from the eponymous member variables in the AtomGraph
  // object and should be taken to mean "... for each topology in the work plan." The pointers
  // target the ARRAY-kind object int_system_data, at intervals of the number of systems in the
  // work plan, with no padding for the warp size.  If these POINTER-kind objects list bounds, the
  // bounds themselves will reflect padding for the HPC warp size.
  Hybrid<int> atom_counts;             ///< Total number of atoms and virtual sites
  Hybrid<int> residue_counts;          ///< Total number of residues, including solvent molecules
  Hybrid<int> molecule_counts;         ///< Total number of molecules in the system
  Hybrid<int> largest_residue_sizes;   ///< Number of atoms in the largest residue
  Hybrid<int> last_solute_residues;    ///< Last residue of the solute, indexed according to the
                                       ///<   overall order in this synthesis rather than the
                                       ///<   original topologies
  Hybrid<int> last_solute_atoms;       ///< Last atom of the solute (a typical solute is a
                                       ///<   biomolecule), indexed according to the overall order
                                       ///<   in this synthesis rather than the original topologies
  Hybrid<int> first_solvent_molecules; ///< First molecule in what is deemed to be solvent, indexed
                                       ///<   according to the overall order in this synthesis
                                       ///<   rather than the original topologies

  // Valence term and off-center particle quantities
  Hybrid<int> urey_bradley_term_counts;      ///< Total number of Urey-Bradley angle stretch terms
  Hybrid<int> charmm_impr_term_counts;       ///< Total number of CHARMM impropers
  Hybrid<int> cmap_term_counts;              ///< Total number of CMAP terms
  Hybrid<int> bond_term_counts;              ///< Total numbers of bonded interactions
  Hybrid<int> angl_term_counts;              ///< Total numbers of bond angle interactions
  Hybrid<int> dihe_term_counts;              ///< Total numbers of dihedral cosine terms
  Hybrid<int> virtual_site_counts;           ///< Number of v-sites / extra points out of all atoms

  // Information relevant to non-bonded calculations
  Hybrid<int> atom_type_counts;        ///< Number of distinct Lennard-Jones types (this is more
                                       ///<   for record-keeping, as all systems in this synthesis
                                       ///<   will work from the same table)
  Hybrid<int> total_exclusion_counts;  ///< Total number of non-bonded exclusions, including 1-4

  // Information relevant to the MD propagation algorithm
  Hybrid<int> rigid_water_counts;       ///< Number of rigid water molecules subject to SETTLE
  Hybrid<int> bond_constraint_counts;   ///< Bonds with lengths constrained by SHAKE or RATTLE
  Hybrid<int> degrees_of_freedom;       ///< Total degrees of freedom, 3N - 6 - constraints in
                                        ///<   each system
  Hybrid<int> nonrigid_particle_counts; ///< A rigid water is one non-rigid particle.  A protein
                                        ///<   with N atoms and no bond length constraints is N
                                        ///<   particles.

  // The following are bounds arrays for the larger data arrays below: more segments of memory,
  // each the number of systems in length.  The contents of these arrays are, with some
  // exceptions, prefix sums over the associated counts arrays above, taking into account zero
  // padding for the stuff they are providing indices into.  A series of systems with atom counts
  // 31, 59, and 78 with an HPC warp size of 32 would generate an atom bounds array stating 0, 32,
  // and 96.  There is no capping value, and the (k+1)th bound minus the kth bounds does not
  // provide information on the exact number associated with the kth system.  Programs should read
  // the exact number of atoms or force field terms per system and the starting bound from one of
  // the arrays below.  The arrays for valence terms spanning all systems provide the global
  // indices of atoms in unified (but not re-ordered) arrays of atom descriptors and also
  // coordinates.  In this sense, these arrays are intermediates between the valence work units
  // that follow and the original AtomGraph objects used to construct this AtomGraphSynthesis.
  Hybrid<int> atom_offsets;           ///< Starting indices for each system's various atomic
                                      ///<   descriptor arrays, which will have eponymous
                                      ///<   names to their AtomGraph counterparts.
  Hybrid<int> residue_offsets;        ///< Starting indices for each system's residue
                                      ///<   descriptor arrays
  Hybrid<int> molecule_offsets;       ///< Starting indices for each system's molecule descriptor
                                      ///<   arrays
  Hybrid<int> ubrd_term_offsets;      ///< Starting indices for Urey-Bradley parameter index
                                      ///<   lists for each system
  Hybrid<int> cimp_term_offsets;      ///< CHARMM improper parameter index lists
  Hybrid<int> cmap_term_offsets;      ///< CMAP atom indexing lists
  Hybrid<int> bond_term_offsets;      ///< Bond stretching term parameter index lists
  Hybrid<int> angl_term_offsets;      ///< Angle bending term parameter index lists
  Hybrid<int> dihe_term_offsets;      ///< Dihedral term parameter index lists
  Hybrid<int> virtual_site_offsets;   ///< Starting indices for each system's virtual site
                                      ///<   descriptor lists
  Hybrid<int> nb_exclusion_offsets;   ///< Starting indices for unified nonbonded exclusion
                                      ///<   lists
  
  /// Data array to hold all system size information, the target of all the above pointers
  Hybrid<int> int_system_data;          

  // Maps of each system and its distinguishable parts.  From here onwards, the AtomGraphSynthesis
  // keeps ARRAY-kind Hybrid Object member variables, with their own data.  As above, many of the
  // names match what is in the AtomGraph object.  The arrays are laid out such that each system's
  // details occur in one contiguous stretch, padded by blank elements up to a multiple of the
  // HPC warp size.  Some of these arrays are indices themselves, but will still need a bounds
  // array, i.e. the point at which to start reading residue limits for the kth system.  Bounds
  // arrays for individual systems within each of these objects are written into the collected
  // integer system data array above, accessed through one of the preceding Hybrid POINTER-kind
  // objects.

  // Atom and residue details
  Hybrid<int> residue_limits;             ///< Atomic indices marking boundaries of individual
                                          ///<   residues
  Hybrid<int> atom_struc_numbers;         ///< Structure atom numbers, such as those taken from a
                                          ///<   PDB file
  Hybrid<int> residue_numbers;            ///< Residue numbers, such as those taken from a PDB file
  Hybrid<int> molecule_limits;            ///< Atomic indices marking the boundaries of molecules
                                          ///<   in the molecule_contents array
  Hybrid<int> atomic_numbers;             ///< Atomic numbers for atoms in all systems
  Hybrid<int> mobile_atoms;               ///< Atom mobility masks for each system
  Hybrid<int> molecule_membership;        ///< Distinct molecules, indexed from 0 for each system,
                                          ///<   indicating molecules to which each atom belongs
  Hybrid<int> molecule_contents;          ///< Indices of atoms making up each molecule, starting
                                          ///<   from the lower bound atom index for each system,
                                          ///<   bounded by molecule_limits
  Hybrid<double> atomic_charges;          ///< Partial charges on all atoms and virtual sites
  Hybrid<double> atomic_masses;           ///< Masses of all atoms
  Hybrid<double> inverse_atomic_masses;   ///< Inverse masses for all atoms
  Hybrid<float> sp_atomic_charges;        ///< Partial charges on all atoms (single precision)
  Hybrid<float> sp_atomic_masses;         ///< Masses of all atoms (single precision)
  Hybrid<float> sp_inverse_atomic_masses; ///< Inverse masses for all atoms (single precision)
  Hybrid<char4> atom_names;               ///< Four letter names of all atoms, i.e. PDB names
  Hybrid<char4> atom_types;               ///< Four letter names of all atom types, i.e. CT
  Hybrid<char4> residue_names;            ///< Four letter names of all residues, i.e. PDB names

  // CHARMM valence term details
  Hybrid<int> urey_bradley_parameter_indices; ///< Index of Urey-Bradley parameters for each term
  Hybrid<int> charmm_impr_parameter_indices;  ///< Index of each CHARMM improper's parameter pair
  Hybrid<int> cmap_surface_indices;           ///< Surface indices for each CMAP term
  Hybrid<int> cmap_surface_dimensions;        ///< Dimensions for every unique CMAP surface
  Hybrid<double> urey_bradley_stiffnesses;    ///< Stiffness constant of each Urey-Bradley stretch
  Hybrid<double> urey_bradley_equilibria;     ///< Equilibrium length of each Urey-Bradley stretch
  Hybrid<double> charmm_impr_stiffnesses;     ///< CHARMM impropers are harmonic, too!
  Hybrid<double> charmm_impr_phase_angles;    ///< The "equilibria" for CHARMM impropers
  Hybrid<double> cmap_surfaces;               ///< Concatenated, column-major format matrices for
                                              ///<   every CMAP surface term
  Hybrid<float> sp_urey_bradley_stiffnesses;  ///< Stiffness constant of each Urey-Bradley stretch,
                                              ///<   in single precision
  Hybrid<float> sp_urey_bradley_equilibria;   ///< Equilibrium length of each Urey-Bradley stretch,
                                              ///<   in single precision
  Hybrid<float> sp_charmm_impr_stiffnesses;   ///< CHARMM impropers are harmonic, too!
  Hybrid<float> sp_charmm_impr_phase_angles;  ///< The "equilibria" for CHARMM impropers
  Hybrid<float> sp_cmap_surfaces;             ///< Concatenated, column-major format matrices for
                                              ///<   every CMAP surface term, in single precision

  // Basic force field valence term details
  Hybrid<int> bond_parameter_indices;    ///< Parameter index of each bond stretching term
  Hybrid<int> angl_parameter_indices;    ///< Parameter index of each angle bending term
  Hybrid<int> dihe_parameter_indices;    ///< Parameter index of each angle bending term
  Hybrid<double> bond_stiffnesses;       ///< Stiffness of each bond stretch, kcal/mol-A^2
  Hybrid<double> bond_equilibria;        ///< Equilibrium lengths of all bonds, A
  Hybrid<double> angl_stiffnesses;       ///< Stiffness of each angle bend, kcal/mol-rad^2
  Hybrid<double> angl_equilibria;        ///< Equilibrium angle for all bending terms, radians
  Hybrid<double> dihe_amplitudes;        ///< Amplitudes of each dihedral cosine term, kcal/mol
  Hybrid<double> dihe_periodicities;     ///< Periodicity of each dihedral / torsion cosine term
  Hybrid<double> dihe_phase_angles;      ///< Phase angle of each dihedral / torsion cosine term
  Hybrid<float> sp_bond_stiffnesses;     ///< Stiffness of each bond stretch (single precision)
  Hybrid<float> sp_bond_equilibria;      ///< Equilibrium lengths of all bonds (single precision)
  Hybrid<float> sp_angl_stiffnesses;     ///< Angle bending stiffnesses (single precision)
  Hybrid<float> sp_angl_equilibria;      ///< Angle bending equilibria (single precision)
  Hybrid<float> sp_dihe_amplitudes;      ///< Amplitudes of torsion cosine terms (single precision)
  Hybrid<float> sp_dihe_periodicities;   ///< Periodicities of torsion terms (single precision)
  Hybrid<float> sp_dihe_phase_angles;    ///< Phase angles of torsion terms (single precision)

  // NMR restraint term details: these function exactly like other parameter sets and are indexed
  // by lists of atoms in the bond work units arrays.  They can be included in the synthesis of
  // AtomGraphs due to their nature as potential terms, whereas the original topologies had to be
  // read from files that did not contain such terms.
  Hybrid<int> nmr_initial_steps;          ///< Initial steps at which to begin applying each NMR
                                          ///<   restraint potential
  Hybrid<int> nmr_final_steps;            ///< Final steps at which to apply the NMR restraint
                                          ///<    potential
  Hybrid<int> nmr_increments;             ///< The number of increments over which to modulate NMR
                                          ///<   restraints
  Hybrid<double2> nmr_k_initial_values;   ///< Stiffnesses of parabolic components of each half of
                                          ///<   the NMR restraint potentials at nmr_initial_steps
  Hybrid<double4> nmr_r_initial_values;   ///< R1 and R2 values for the left-hand side of the NMR
                                          ///<   restraint potentials at nmr_initial_steps
  Hybrid<double2> nmr_k_final_values;     ///< Stiffnesses of the parabolic components of each half
                                          ///<   of the NMR restraint potentials at nmr_final_steps
  Hybrid<double4> nmr_r_final_values;     ///< R1 and R2 values for the left-hand side of the NMR
                                          ///<   restraint potentials at nmr_final_steps
  Hybrid<float2> sp_nmr_k_initial_values; ///< Single precision form of the initial NMR parabolic
                                          ///<   restraint stiffnesses
  Hybrid<float4> sp_nmr_r_initial_values; ///< Single precision form of initial NMR displacement
                                          ///<   demarcations
  Hybrid<float2> sp_nmr_r_final_values;   ///< Single precision form of the final NMR parabolic
                                          ///<   restraint stiffnesses
  Hybrid<float4> sp_nmr_r_final_values;   ///< Single precision form of final NMR displacement
                                          ///<   demarcations

  // Valence work units (VWUs): work units providing instruction sets for the GPU to operate on a
  // continuous, non-rearranged list of atoms and implement all valence terms.  Each VWU pertains
  // to one and only one individual topology from the list and each topology will have at least one
  // VWU, or more depending on what is needed to cover all of its atoms.  The VWU imports a list of
  // atoms, then computes a series of bond, angle, dihedral, and other force field terms based on
  // the cached atoms.  Forces are accumulated on all atoms in __shared__ memory and then
  // contributed back to the global arrays.

  /// A list of atoms that the VWU shall import, with indices into the global array of atoms for
  /// all systems.  Each VWU may import up to half as many atoms as the kernel blocks have threads,
  /// and each VWU takes a stride of that many ints from this array.
  Hybrid<int> atom_imports;

  /// Instruction sets for the bond work units, up to 64 integers for the low and high limits of
  /// up to 32 types of interactions in a given work unit.  Each VWU takes a stride of 64 ints from
  /// this array.  The length of this stride is not a matter of the warp size.
  Hybrid<int> vwu_instruction_sets;

  /// Instructions for bond stretching and Urey-Bradley interactions.  Each uint2 tuple contains
  /// two atom indices in the x member (bits 1-10 and 11-20) and the parameter index of the bond /
  /// Urey-Bradley term in the y member.  A flag in the 21st bit of the x member indicates whether
  /// the energy should contribute to the Urey-Bradley sum (if set to 1) or the harmonic bond sum.
  Hybrid<uint2> bond_instructions;

  /// Instructions for angle bending and three-body NMR restraint interactions.  Each uint2 tuple
  /// contains three atom indices in the x member (bits 1-10, 11-20, and 21-30) and the parameter
  /// index of the angle in the y member.
  Hybrid<uint2> angl_instructions;

  /// Instructions for dihedral and CHARMM improper dihedral interactions.  Each uint2 tuple
  /// contains three atom indices in the x member (bits 1-10, 11-20, and 21-30), plus a flag in
  /// bit 31 to indicate whether the form of the potential should be a cosine term or the CHARMM
  /// harmonic term and another flag in bit 32, relevant only when energy is being computed and
  /// only if the term is not a CHARMM improper, to indicate whether to add the energy to a proper
  /// improper accumulator.  A fourth atom index in the y member (bits 1-10) precedes the index for
  /// 1:4 scaling factors (bits 11-15, up to 32 unique combinations of electrostatic and
  /// Lennard-Jones scaling factors, including zero to indicate that no 1:4 interaction should be
  /// computed).  The final bits (16-32) indiate the dihedral or CHARMM improper dihedral parameter
  /// index (up to 131072 of either, a number larger than nearly any force field, but the limit
  /// applies to the number of unique parameters at work in a particular simulation).
  Hybrid<uint2> dihe_instructions;

  /// Instructions for CMAP interactions.  Each uint2 tuple contains three atom indices in the x
  /// member (bits 1-10, 11-20, and 21-30, plus two more in the y member (bits 1-10 and 11-20).
  /// The CMAP parameter index is given in the final twelve bits of the y member.
  Hybrid<uint2> cmap_instructions;

  /// Instructions for NMR two-body restraints.  Each uint2 tuple contains two atom indices in the
  /// x member (bits 1-10 and 11-20) followed by the index of the slope / intercept series in the
  /// y member.
  Hybrid<uint2> nmr2_instructions;

  /// Instructions for NMR three-body restraints.  Each uint2 tuple contains three atom indices in
  /// the x member (bits 1-10, 11-20, and 21-30) followed by the index of the slope / intercept
  /// series in the y member.
  Hybrid<uint2> nmr3_instructions;

  /// Instructions for NMR four-body restraints.  Each uint2 tuple contains three atom indices in
  /// the x member (bits 1-10, 11-20, and 21-30) followed by a fourth in the y member (bits 1-10)
  /// and finally the index of the slope / intercept series (bits 11-32).
  Hybrid<uint2> nmr4_instructions;
};

#endif
