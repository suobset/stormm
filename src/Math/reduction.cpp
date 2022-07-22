#include "reduction.h"

namespace omni {
namespace math {

//-------------------------------------------------------------------------------------------------
ReductionKit::ReductionKit(const int nrdwu_in, const RdwuPerSystem rps_in,
                           const int* rdwu_abstracts_in, const int* atom_counts_in) :
  nrdwu{nrdwu_in}, rps{rps_in}, rdwu_abstracts{rdwu_abstracts_in}, atom_counts{atom_counts_in}
{}

} // namespace math
} // namespace omni
