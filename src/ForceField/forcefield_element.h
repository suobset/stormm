// -*-c++-*-
#ifndef OMNI_FORCEFIELD_ELEMENT_H
#define OMNI_FORCEFIELD_ELEMENT_H

#include <vector>
#include "DataTypes/omni_vector_types.h"
#include "Topology/atomgraph_enumerators.h"
#include "forcefield_enumerators.h"

namespace omni {
namespace modeling {

using topology::VirtualSiteKind;

/// \brief A versatile object for collecting the parameters and scope of applicability of any
///        molecular mechanics force field term.  This relies on enumerations to inform whether
///        the term applies to atom types or atom and residue names and the nature of the term.
class ForceFieldElement {
public:

  /// \brief A variety of constructors can load one or more atoms and their properties.
  ///
  /// Overloaded:
  ///   - Create an empty element with only a ParameterKind (default NONE)
  ///   - Load non-bonded terms describing a specific atom
  ///   - Load valence terms describing a bond, angle, dihedral, Urey-Bradley, or CHARMM improper
  ///   - Load an entire CMAP surface
  ///   - Load a virtual site with frame specifications
  ///
  /// \param kind_in        The kind of force field parameter, obligatory for every constructor
  /// \param atom_i_in      Name or type of atom I in the term
  /// \param atom_j_in      Name or type of atom J in the term
  /// \param atom_k_in      Name or type of atom K in the term
  /// \param atom_l_in      Name or type of atom L in the term
  /// \param atom_m_in      Name or type of atom M in the term
  /// \param resi_i_in      Name or type of residue I in the term
  /// \param resi_j_in      Name or type of residue J in the term
  /// \param resi_k_in      Name or type of residue K in the term
  /// \param resi_l_in      Name or type of residue L in the term
  /// \param resi_m_in      Name or type of residue M in the term
  /// \param prop_a_in      Charge, Lennard-Jones sigma, valence parameter stiffness or amplitude,
  ///                         or virtual site frame dimension 1
  /// \param prop_b_in      Lennard-Jones epsilon, valence parameter equilibrium or phase angle, or
  ///                         virtual site frame dimension 2
  /// \param prop_c_in      Lennard-Jones tertiary parameter (12-6-4 or Buckingham potential),
  ///                         dihedral periodicity, or virtual site frame dimension 3
  /// \param surface_in     Surface values for an entire CMAP
  /// \param frame_type_in  Frame type of the virtual site, if that is what this object contains
  /// \{
  ForceFieldElement(ParameterKind kind_in = ParameterKind::NONE);
  
  ForceFieldElement(ParameterKind kind_in, char4 atom_i_in, double prop_a_in,
                    double prop_b_in = 0.0, double prop_c_in = 0.0);

  ForceFieldElement(ParameterKind kind_in, char4 atom_i_in, char4 atom_j_in, double prop_a_in,
                    double prop_b_in = 0.0, double prop_c_in = 0.0);

  ForceFieldElement(ParameterKind kind_in, char4 atom_i_in, char4 atom_j_in, char4 atom_k_in,
                    double prop_a_in, double prop_b_in = 0.0, double prop_c_in = 0.0);

  ForceFieldElement(ParameterKind kind_in, char4 atom_i_in, char4 atom_j_in, char4 atom_k_in,
                    char4 atom_l_in, double prop_a_in, double prop_b_in = 0.0,
                    double prop_c_in = 0.0);

  ForceFieldElement(ParameterKind kind_in, VirtualSiteKind frame_type_in, char4 atom_i_in,
                    char4 atom_j_in, char4 atom_k_in, char4 residue_i_in, char4 residue_j_in,
                    char4 residue_k_in, double prop_a_in, double prop_b_in = 0.0,
                    double prop_c_in = 0.0);

  ForceFieldElement(ParameterKind kind_in, VirtualSiteKind frame_type_in, char4 atom_i_in,
                    char4 atom_j_in, char4 atom_k_in, char4 atom_l_in, char4 residue_i_in,
                    char4 residue_j_in, char4 residue_k_in, char4 residue_l_in, double prop_a_in,
                    double prop_b_in = 0.0, double prop_c_in = 0.0);

  ForceFieldElement(ParameterKind kind_in, VirtualSiteKind frame_type_in, char4 atom_i_in,
                    char4 atom_j_in, char4 atom_k_in, char4 atom_l_in, char4 atom_m_in,
                    char4 residue_i_in, char4 residue_j_in, char4 residue_k_in, char4 residue_l_in,
                    char4 residue_m_in, double prop_a_in = 0.0, double prop_b_in = 0.0,
                    double prop_c_in = 0.0);

  ForceFieldElement(ParameterKind kind_in, char4 atom_i_in, char4 atom_j_in, char4 atom_k_in,
                    char4 atom_l_in, char4 atom_m_in, char4 residue_i_in, char4 residue_j_in,
                    char4 residue_k_in, char4 residue_l_in, char4 residue_m_in,
                    const std::vector<double> &surface_in);
  /// \}

  /// \brief Get the force field parameter kind
  ParameterKind getKind() const;
  
  /// \brief Get the atom name of the I atom in the term
  char4 getNameOfAtom(char atom_rank = 'I') const;
  
  /// \brief Get the atom type of atom I in the term
  char4 getTypeOfAtom(char atom_rank = 'I') const;
    
  /// \brief Get the resiude name of atom I in the term
  char4 getNameOfResidue(char atom_rank = 'I') const;

  /// \brief Get the charge of an atom with a given atom name and residue nam
  double getCharge() const;

  /// \brief Get the Lennard-Jones sigma parameter of an atom
  double getSigma() const;

  /// \brief Get the Lennard-Jones epsilon parameter of an atom
  double getEpsilon() const;

  /// \brief Get the Lennard-Jones rho parameter of an atom (the third parameter, for 12-6-4
  ///        potentials)
  double getRho() const;

  /// \brief Get the stiffness constant of a bond, angle, Urey-Bradley, or CHARMM improper term.
  double getStiffnessConstant() const;

  /// \brief Get the equilibrium constant of a bond, angle, or Urey-Bradley term.
  double getEquilibriumConstant() const;

  /// \brief Get the amplitude of a cosine-based dihedral term.
  double getAmplitude() const;

  /// \brief Get the phase angle of a cosine-based dihedral or CHARMM improper dihedral term.
  double getPhaseAngle() const;

  /// \brief Get the periodicity of a cosine-based dihedral term
  double getPeriodicity() const;

  /// \brief Get the electrostatic scaling factor for an attenuated 1:4 non-bonded interaction.
  double getElectrostaticScaling() const;
  
  /// \brief Get the van-der Waals scaling factor for an attenuated 1:4 non-bonded interaction.
  double getVanDerWaalsScaling() const;
  
  /// \brief Get the surface associated with a CMAP term.
  std::vector<double> getSurface() const;
  
  /// \brief Get the dimension of a CMAP term's energy surface.  This property is inferred from the
  ///        square root of the size of the input array.
  int getSurfaceDimension() const;  
  
  /// \brief Get the virtual site frame type
  VirtualSiteKind getVirtualSiteFrameType() const;

private:
  ParameterKind kind;    ///< The type of parameter, i.e. BOND or VIRTUAL_SITE_FRAME
  char4 atom_name_i;     ///< Atom or atom type name for atom I of some valence term's arrangement,
                         ///<   or the name of a virtual site atom
  char4 atom_name_j;     ///< Atom or atom type name for atom J of some valence term's arrangement,
                         ///<   or the name of the parent atom in a virtual site frame
  char4 atom_name_k;     ///< Atom or atom type name for atom K of some valence term's arrangement,
                         ///<   or the name of the second frame atom in a virtual site frame
  char4 atom_name_l;     ///< Atom or atom type name for atom L of some valence term's arrangement,
                         ///<   or the name of the third frame atom in a virtual site frame
  char4 atom_name_m;     ///< Atom or atom type name for atom M of some valence term's arrangement,
                         ///<   or the name of the fourth frame atom in a virtual site frame
  char4 residue_name_i;  ///< Residue name for atom I of some valence term arrangement, or a
                         ///<   virtual site atom
  char4 residue_name_j;  ///< Residue name for atom J of some valence term arrangement, or a 
                         ///<   virtual site frame's parent atom
  char4 residue_name_k;  ///< Residue name for atom J of some valence term arrangement, or a 
                         ///<   virtual site frame's second atom
  char4 residue_name_l;  ///< Residue name for atom J of some valence term arrangement, or a 
                         ///<   virtual site frame's third atom
  char4 residue_name_m;  ///< Residue name for atom J of some valence term arrangement, or a 
                         ///<   virtual site frame's fourth atom

  // General-purpose real-values numbers for keeping this parameter's details
  double property_a;  ///< Charge, Lennard-Jones sigma, valence parameter stiffness or amplitude,
                      ///<   or virtual site frame dimension 1
  double property_b;  ///< Lennard-Jones epsilon, valence parameter equilibrium or phase angle, or
                      ///<   virtual site frame dimension 2
  double property_c;  ///< Lennard-Jones tertiary parameter (12-6-4 or Buckingham potential),
                      ///<   dihedral periodicity, or virtual site frame dimension 3

  // Miscellaneous properties
  std::vector<double> surface;  ///< Surface values for an entire CMAP
  int surface_dimension;        ///< Width and length of the square CMAP surface
  VirtualSiteKind frame_type;   ///< Frame type of the virtual site, if that is what this object
                                ///<   contains.
};
  
} // namespace modeling
} // namespace omni

#endif
