#include "forcefield_element.h"

namespace omni {
namespace modeling {

//-------------------------------------------------------------------------------------------------
ForceFieldElement::ForceFieldElement(const ParameterKind kind_in) :
    kind{kind_in}, atom_name_i{' ', ' ', ' ', ' '}, atom_name_j{' ', ' ', ' ', ' '},
    atom_name_k{' ', ' ', ' ', ' '}, atom_name_l{' ', ' ', ' ', ' '},
    atom_name_m{' ', ' ', ' ', ' '}, residue_name_i{' ', ' ', ' ', ' '},
    residue_name_j{' ', ' ', ' ', ' '}, residue_name_k{' ', ' ', ' ', ' '},
    residue_name_l{' ', ' ', ' ', ' '}, residue_name_m{' ', ' ', ' ', ' '}, prop_a{0.0},
    prop_b{0.0}, prop_c{0.0}, surface{}, frame_type{VirtualSiteKind::NONE}
    
//-------------------------------------------------------------------------------------------------
ForceFieldElement::ForceFieldElement(ParameterKind kind_in, char4 atom_i_in,
                                     char4 atom_j_in, double prop_a_in,
                                     double prop_b_in, double prop_c_in) :
  ForceFieldElement(kind_in)
  
} // namespace modeling
} // namespace omni
