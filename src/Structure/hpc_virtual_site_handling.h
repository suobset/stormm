// -*-c++-*-
#ifndef STORMM_HPC_VIRTUAL_SITE_HANDLING_H
#define STORMM_HPC_VIRTUAL_SITE_HANDLING_H

#include "Accelerator/kernel_manager.h"
#include "Constants/behavior.h"
#include "Potential/cacheresource.h"
#include "Synthesis/atomgraph_synthesis.h"
#include "Synthesis/phasespace_synthesis.h"
#include "Synthesis/synthesis_abstracts.h"
#include "structure_enumerators.h"

namespace stormm {
namespace structure {

using card::KernelManager;
using constants::PrecisionModel;
using energy::CacheResource;
using energy::CacheResourceKit;
using synthesis::AtomGraphSynthesis;
using synthesis::PhaseSpaceSynthesis;
using synthesis::PsSynthesisWriter;
using synthesis::SyAtomUpdateKit;
using synthesis::SyValenceKit;

/// \brief Get kernel attributes for standalone virtual site handling kernels.
///
/// \param prec     The precision model to use in arithmetic, as well as an indicator of whether
///                 the fixed-precision coordinate representation has extensions to 95 bits
/// \param purpose  Whether to place the virtual sites (after motion of their frame atoms) or to
///                 transmit forces from the virtual sites to their frame atoms
cudaFuncAttributes queryVirtualSiteKernelRequirements(PrecisionModel prec,
                                                      VirtualSiteActivity purpose);
  
/// \brief Place virtual sites based on new atom positions within the current coordinate set.
///
/// Overloaded:
///   - Accept pre-obtained abstracts (the PhaseSpaceSynthesis abstract must point to the correct,
///     current coordinate set
///   - Accept the original coordinate and topology objects (abstracts will be created internally)
///
/// \param prec     The precision model to use
/// \param poly_ps  Coordinates for all atoms in all systems
/// \param poly_ag  Topologies for all systems
/// \{
void launchVirtualSitePlacement(PsSynthesisWriter *poly_psw, CacheResourceKit<double> *gmem_r,
                                const SyValenceKit<double> &poly_vk,
                                const SyAtomUpdateKit<double2, double4> &poly_auk, const int2 bt);

void launchVirtualSitePlacement(PsSynthesisWriter *poly_psw, CacheResourceKit<float> *gmem_r,
                                const SyValenceKit<float> &poly_vk,
                                const SyAtomUpdateKit<float2, float4> &poly_auk, const int2 bt);

void launchVirtualSitePlacement(PrecisionModel prec, PhaseSpaceSynthesis *poly_ps,
                                CacheResource *tb_space, const AtomGraphSynthesis &poly_ag,
                                const KernelManager &launcher);
/// \}

} // namespace structure
} // namespace stormm

#endif
