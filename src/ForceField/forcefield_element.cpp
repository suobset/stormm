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
{}

//-------------------------------------------------------------------------------------------------
ForceFieldElement::ForceFieldElement(const ParameterKind kind_in, const char4 atom_i_in,
                                     const double prop_a_in, const double prop_b_in,
                                     const double prop_c_in) :
  ForceFieldElement(kind_in)
{
  atom_name_i = atom_i_in;
  prop_a = prop_a_in;
  prop_b = prop_b_in;
  prop_c = prop_c_in;
}

  //-------------------------------------------------------------------------------------------------
ForceFieldElement::ForceFieldElement(const ParameterKind kind_in, const char4 atom_i_in,
                                     const char4 atom_j_in, const double prop_a_in,
                                     const double prop_b_in, const double prop_c_in) :
  ForceFieldElement(kind_in)
{
  atom_name_i = atom_i_in;
  atom_name_j = atom_j_in;
  prop_a = prop_a_in;
  prop_b = prop_b_in;
  prop_c = prop_c_in;
}

//-------------------------------------------------------------------------------------------------
ForceFieldElement::ForceFieldElement(const ParameterKind kind_in, const char4 atom_i_in,
                                     const char4 atom_j_in, const char4 atom_k_in,
                                     const double prop_a_in, const double prop_b_in,
                                     const double prop_c_in) :
  ForceFieldElement(kind_in)
{
  atom_name_i = atom_i_in;
  atom_name_j = atom_j_in;
  atom_name_k = atom_k_in;
  prop_a = prop_a_in;
  prop_b = prop_b_in;
  prop_c = prop_c_in;
}

//-------------------------------------------------------------------------------------------------
ForceFieldElement::ForceFieldElement(const ParameterKind kind_in, const char4 atom_i_in,
                                     const char4 atom_j_in, const char4 atom_k_in,
                                     const char4 atom_l_in, const double prop_a_in,
                                     const double prop_b_in, const double prop_c_in) :
  ForceFieldElement(kind_in)
{
  atom_name_i = atom_i_in;
  atom_name_j = atom_j_in;
  atom_name_k = atom_k_in;
  atom_name_l = atom_l_in;
  prop_a = prop_a_in;
  prop_b = prop_b_in;
  prop_c = prop_c_in;
}

//-------------------------------------------------------------------------------------------------
ForceFieldElement::ForceFieldElement(const ParameterKind kind_in,
                                     const VirtualSiteKind frame_type_in, const char4 atom_i_in,
                                     const char4 atom_j_in, const char4 atom_k_in,
                                     const char4 residue_i_in, const char4 residue_j_in,
                                     const char4 residue_k_in, const double prop_a_in,
                                     const double prop_b_in, const double prop_c_in) :
  ForceFieldElement(kind_in)
{
  atom_name_i = atom_i_in;
  atom_name_j = atom_j_in;
  atom_name_k = atom_k_in;
  residue_name_i = residue_i_in;
  residue_name_j = residue_j_in;
  residue_name_k = residue_k_in;
  prop_a = prop_a_in;
  prop_b = prop_b_in;
  prop_c = prop_c_in;
  frame_type = frame_type_in;
}

//-------------------------------------------------------------------------------------------------
ForceFieldElement::ForceFieldElement(const ParameterKind kind_in,
                                     const VirtualSiteKind frame_type_in, const char4 atom_i_in,
                                     const char4 atom_j_in, const char4 atom_k_in,
                                     const char4 atom_l_in, const char4 residue_i_in,
                                     const char4 residue_j_in, const char4 residue_k_in,
                                     const char4 residue_l_in, const double prop_a_in,
                                     const double prop_b_in, const double prop_c_in) :
  ForceFieldElement(kind_in)
{
  atom_name_i = atom_i_in;
  atom_name_j = atom_j_in;
  atom_name_k = atom_k_in;
  atom_name_l = atom_l_in;
  residue_name_i = residue_i_in;
  residue_name_j = residue_j_in;
  residue_name_k = residue_k_in;
  residue_name_l = residue_l_in;
  prop_a = prop_a_in;
  prop_b = prop_b_in;
  prop_c = prop_c_in;
  frame_type = frame_type_in;
}

//-------------------------------------------------------------------------------------------------
ForceFieldElement::ForceFieldElement(const ParameterKind kind_in,
                                     const VirtualSiteKind frame_type_in, const char4 atom_i_in,
                                     const char4 atom_j_in, const char4 atom_k_in,
                                     const char4 atom_l_in, const char4 atom_m_in,
                                     const char4 residue_i_in, const char4 residue_j_in,
                                     const char4 residue_k_in, const char4 residue_l_in,
                                     const char4 residue_m_in, const double prop_a_in,
                                     const double prop_b_in, const double prop_c_in) :
  ForceFieldElement(kind_in)
{
  atom_name_i = atom_i_in;
  atom_name_j = atom_j_in;
  atom_name_k = atom_k_in;
  atom_name_l = atom_l_in;
  atom_name_m = atom_m_in;
  residue_name_i = residue_i_in;
  residue_name_j = residue_j_in;
  residue_name_k = residue_k_in;
  residue_name_l = residue_l_in;
  residue_name_m = residue_m_in;
  prop_a = prop_a_in;
  prop_b = prop_b_in;
  prop_c = prop_c_in;
  frame_type = frame_type_in;
}

//-------------------------------------------------------------------------------------------------
ForceFieldElement::ForceFieldElement(const ParameterKind kind_in, const char4 atom_i_in,
                                     const char4 atom_j_in, const char4 atom_k_in,
                                     const char4 atom_l_in, const char4 atom_m_in,
                                     const char4 residue_i_in, const char4 residue_j_in,
                                     const char4 residue_k_in, const char4 residue_l_in,
                                     const char4 residue_m_in,
                                     const std::vector<double> &surface_in) :
  ForceFieldElement(kind_in)
{
  atom_name_i = atom_i_in;
  atom_name_j = atom_j_in;
  atom_name_k = atom_k_in;
  atom_name_l = atom_l_in;
  atom_name_m = atom_m_in;
  residue_name_i = residue_i_in;
  residue_name_j = residue_j_in;
  residue_name_k = residue_k_in;
  residue_name_l = residue_l_in;
  residue_name_m = residue_m_in;
  const size_t nsurf = surface_in.size();
  surface.resize(nsurf);
  for (size_t i = 0; i < nsurf; i++) {
    surface[i] = surface_in[i];
  }
}

} // namespace modeling
} // namespace omni
