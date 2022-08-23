#include "Constants/hpc_bounds.h"
#include "Math/rounding.h"
#include "implicit_solvent_workspace.h"

namespace stormm {
namespace synthesis {

using card::HybridKind;
using math::roundUp;

//-------------------------------------------------------------------------------------------------
ISWorkspaceKit::ISWorkspaceKit(const int fp_bits_in, llint* psi_in, int* psi_overflow_in,
                               llint* sum_deijda_in, int* sum_deijda_overflow_in) :
    fp_bits{fp_bits_in}, psi{psi_in}, psi_overflow{psi_overflow_in}, sum_deijda{sum_deijda_in},
    sum_deijda_overflow{sum_deijda_overflow_in}
{}

//-------------------------------------------------------------------------------------------------
ImplicitSolventWorkspace::ImplicitSolventWorkspace(const Hybrid<int> &atom_starts,
                                                   const Hybrid<int> &atom_counts,
                                                   const int bit_count) :
    fp_bits{bit_count},
    psi{HybridKind::ARRAY, "isw_eff_gbrad"},
    psi_overflow{HybridKind::ARRAY, "isw_eff_gbrad_ovrf"},
    sum_deijda{HybridKind::ARRAY, "isw_sumdeijda"},
    sum_deijda_overflow{HybridKind::ARRAY, "isw_sumdeijda_ovrf"}
{
  const size_t last_sys = atom_starts.size() - 1LLU;
  const int padded_natom = atom_starts.readHost(last_sys) + roundUp(atom_counts.readHost(last_sys),
                                                                    warp_size_int);
  psi.resize(padded_natom);
  psi_overflow.resize(padded_natom);
  sum_deijda.resize(padded_natom);
  sum_deijda_overflow.resize(padded_natom);
}

//-------------------------------------------------------------------------------------------------
ImplicitSolventWorkspace::ImplicitSolventWorkspace(const Hybrid<int> &atom_starts,
                                                   const Hybrid<int> &atom_counts,
                                                   const PrecisionModel prec) :
  ImplicitSolventWorkspace(atom_starts, atom_counts, (prec == PrecisionModel::SINGLE) ? 25 : 55)
{}

//-------------------------------------------------------------------------------------------------
int ImplicitSolventWorkspace::getFixedPrecisionBits() const {
  return fp_bits;
}

//-------------------------------------------------------------------------------------------------
ISWorkspaceKit ImplicitSolventWorkspace::data(const HybridTargetLevel tier) {
  return ISWorkspaceKit(fp_bits, psi.data(tier), psi_overflow.data(tier), sum_deijda.data(tier),
                        sum_deijda_overflow.data(tier));
}

#ifdef STORMM_USE_HPC
//-------------------------------------------------------------------------------------------------
void ImplicitSolventWorkspace::upload() {
  psi.upload();
  psi_overflow.upload();
  sum_deijda.upload();
  sum_deijda_overflow.upload();
}

//-------------------------------------------------------------------------------------------------
void ImplicitSolventWorkspace::download() {
  psi.download();
  psi_overflow.download();
  sum_deijda.download();
  sum_deijda_overflow.download();
}
#endif

} // namespace synthesis
} // namespace stormm
