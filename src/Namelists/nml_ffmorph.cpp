#include "namelist_element.h"
#include "nml_ffmorph.h"

namespace omni {
namespace namelist {

//-------------------------------------------------------------------------------------------------
FFMorphControls::FFMorphControls(const ExceptionResponse policy_in) :
    policy{policy_in},
    harmonic_bonds{}, harmonic_angles{}, cosine_dihedrals{}, urey_bradley_angles{},
    charmm_impropers{}, cmap_surfaces{}, van_der_waals_properties{}, charge_properties{},
    virtual_sites{}
{}

} // namespace namelist 
} // namespace omni
