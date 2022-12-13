#include "rmsd_plan.h"

namespace stormm {
namespace structure {

using card::HybridKind;
  
//-------------------------------------------------------------------------------------------------
RmsdPlanReader::RmsdPlanReader(const int plan_count_in, const RmsdMethod strategy_in,
                               const double mass_fraction_in, const int* alignment_steps_in,
                               const int* core_atoms_in, const int* core_counts_in,
                               const int* core_starts_in, const int* symm_atoms_in,
                               const int4* symm_bounds_in, const int* symm_ranges_in) :
    plan_count{plan_count_in}, strategy{strategy_in}, mass_fraction{mass_fraction_in},
    alignment_steps{alignment_steps_in}, core_atoms{core_atoms_in}, core_counts{core_counts_in},
    core_starts{core_starts_in}, symm_atoms{symm_atoms_in}, symm_bounds{symm_bounds_in},
    symm_ranges{symm_ranges_in}
{}

//-------------------------------------------------------------------------------------------------
RmsdPlan::RmsdPlan(const RmsdMethod strategy_in, const double rmf_in) :
    plan_count{0}, general_strategy{strategy_in}, required_mass_fraction{rmf_in},
    alignment_steps{HybridKind::ARRAY, "rplan_align_steps"},
    asymmetric_core_atoms{HybridKind::ARRAY, "rplan_core_atoms"},
    asymmetric_core_counts{HybridKind::ARRAY, "rplan_core_counts"},
    asymmetric_core_starts{HybridKind::ARRAY, "rplan_core_starts"},
    symmetry_group_atoms{HybridKind::ARRAY, "rplan_symm_atoms"},
    symmetry_group_bounds{HybridKind::ARRAY, "rplan_symm_bounds"},
    symmetry_group_ranges{HybridKind::ARRAY, "rplan_symm_ranges"},
    ag_pointers{}
{}

//-------------------------------------------------------------------------------------------------
RmsdPlan::RmsdPlan(const AtomGraph &ag_in, const CoordinateFrame &cf_in,
                   const RmsdMethod strategy_in, const double rmf_in, const GpuDetails &gpu) :
    RmsdPlan(strategy_in, rmf_in)
{
  // There is one plan that can be defined by the given information.
  plan_count = 1;

  // Get the formal charges, free electron content, chirality, and ring inclusions for the one
  // system.  Each system can, in fact, hold many molecules, but this 
}

} // namespace structure
} // namespace stormm
