#include "copyright.h"
#include "Topology/atomgraph_abstracts.h"
#include "nbemulator.h"

namespace stormm {
namespace energy {

using chemistry::ChemicalFeatures;
using chemistry::ChemicalFeaturesReader;
using stmath::foundIn;
using stmath::locateValue;
using stmath::reduceUniqueValues;
using synthesis::SynthesisMapReader;
using topology::NonbondedKit;

//-------------------------------------------------------------------------------------------------
NBEmulator::NBEmulator(const std::string &file_name) :
    stiffness_matrix{0, "nbemulator_amat"},
    row_descriptors{},
    column_descriptors{},
    target_nrg{}
{}

//-------------------------------------------------------------------------------------------------
NBEmulator::NBEmulator(const PhaseSpaceSynthesis &poly_ps, const AtomGraphSynthesis &poly_ag,
                       const SystemCache &sysc, const SynthesisCacheMap &scmap,
                       const EmulatorControls &emulcon) :
    stiffness_matrix{0, "nbemulator_amat"},
    row_descriptors{},
    column_descriptors{},
    target_nrg{}
{
  // Break down the structures according to atomic sources.
  const SynthesisMapReader scmapr = scmap.data();
  const int source_count = emulcon.getSourceCount();
  std::vector<std::vector<int>> source_index_map(poly_ps.getSystemCount());
  for (int i = 0; i < poly_ps.getSystemCount(); i++) {
    source_index_map[i].resize(poly_ps.getAtomCount(i));
    const int cache_idx = scmapr.cache_origins[i];
    const ChemicalFeatures* ichemfe = sysc.getFeaturesPointer(cache_idx);
    const ChemicalFeaturesReader ichemfer = ichemfe->data();

  }
}

} // namespace energy
} // namespace stormm
