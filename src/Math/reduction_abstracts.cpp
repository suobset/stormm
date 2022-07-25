#include "reduction_abstracts.h"

namespace omni {
namespace math {

//-------------------------------------------------------------------------------------------------
ReductionKit::ReductionKit(const int nrdwu_in, const RdwuPerSystem rps_in,
                           const int* rdwu_abstracts_in, const int* atom_counts_in) :
  nrdwu{nrdwu_in}, rps{rps_in}, rdwu_abstracts{rdwu_abstracts_in}, atom_counts{atom_counts_in}
{}

//-------------------------------------------------------------------------------------------------
ReductionKit::ReductionKit(const AtomGraphSynthesis &poly_ag, const HybridTargetLevel tier) :
  nrdwu{poly_ag.getReductionWorkUnitCount()},
  rps{poly_ag.getRdwuPerSystem()},
  rdwu_abstracts{poly_ag.getReductionWorkUnitAbstracts().data(tier)},
  atom_counts{poly_ag.getSystemAtomCounts().data(tier)}
{}

} // namespace math
} // namespace omni
