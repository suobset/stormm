#include <cmath>
#include "../../../src/Constants/scaling.h"
#include "../../../src/Constants/symbol_values.h"
#include "../../../src/Chemistry/chemistry_enumerators.h"
#include "../../../src/Chemistry/chemical_features.h"
#include "../../../src/Math/rounding.h"
#include "../../../src/Structure/isomerization.h"
#include "../../../src/Structure/rmsd.h"
#include "../../../src/Topology/atomgraph.h"
#include "../../../src/Topology/atomgraph_abstracts.h"
#include "../../../src/Trajectory/coordinate_series.h"
#include "setup.h"

namespace conf_app {
namespace setup {

using stormm::chemistry::ChemicalFeatures;
using stormm::chemistry::ChiralInversionProtocol;
using stormm::chemistry::MapRotatableGroups;
using stormm::chemistry::RotatorGroup;
using stormm::constants::warp_size_int;
using stormm::math::roundUp;
using stormm::namelist::ConformerControls;
using stormm::structure::rmsd;
using stormm::structure::RmsdMethod;
using stormm::structure::rotateAboutBond;
using stormm::topology::AtomGraph;
using stormm::topology::ChemicalDetailsKit;
using stormm::trajectory::PhaseSpace;
using stormm::trajectory::CoordinateFrame;
using stormm::trajectory::CoordinateFrameReader;
using stormm::trajectory::CoordinateSeries;
using stormm::trajectory::CoordinateSeriesWriter;

//-------------------------------------------------------------------------------------------------
PhaseSpaceSynthesis expandConformers(const UserSettings &ui, const SystemCache &sc,
                                     StopWatch *tm) {
  const ConformerControls conf_input = ui.getConformerNamelistInfo();

  // Count the expanded number of systems
  const int ntop = sc.getTopologyCount();
  const int nsys = sc.getSystemCount();
  const std::vector<const AtomGraph*> system_topologies = sc.getSystemTopologyPointer();
  const std::vector<const PhaseSpace*> system_coords    = sc.getCoordinatePointer();
  const std::vector<const AtomGraph*> unique_topologies = sc.getTopologyPointer();
  std::vector<ChemicalFeatures> chemfe_list;
  chemfe_list.reserve(ntop);
  for (int i = 0; i < ntop; i++) {
    const int example_system_idx = sc.getCoordinateExample(i);
    chemfe_list.emplace_back(sc.getTopologyPointer(i),
                             sc.getCoordinateReference(example_system_idx),
                             MapRotatableGroups::YES);
  }
  tm->assignTime(2);
  
  // Loop over all systems, grouping those with the same topology into a coherent group of
  // proto-conformers for coarse-grained sampling of rotatable bonds and chiral centers.
  const int nbond_rotations = conf_input.getRotationSampleCount();
  int nneighborhood = 0;
  for (int i = 0; i < ntop; i++) {
    
    // Determine the combinatorial number of conformers that could populate the coarse-grained
    // search.  This includes bond rotation, chiral inversions, and multiple replicas of the
    // compound in the list of starting coordinates.
    const int ncases = sc.getTopologyCaseCount(i);
    const int nrot_bond = std::min(chemfe_list[i].getRotatableBondCount(),
                                   conf_input.getRotatableBondLimit());
    const int nctx_bond = (conf_input.sampleCisTrans()) ? chemfe_list[i].getCisTransBondCount() :
                                                          0;
    double ln_choices = log(ncases);
    if (nrot_bond > 0 && nbond_rotations > 0) {
      ln_choices += nrot_bond * log(nbond_rotations);
    }
    if (nctx_bond > 0) {
      ln_choices += nctx_bond * log(2.0);
    }
    const int nchiral = (conf_input.sampleChirality()) ? chemfe_list[i].getChiralCenterCount() : 0;
    const std::vector<ChiralInversionProtocol> chiral_protocols;
    if (nchiral > 0) {
      int np2 = 0;
      for (int j = 0; j < nchiral; j++) {
        switch (chiral_protocols[j]) {
        case ChiralInversionProtocol::ROTATE:
        case ChiralInversionProtocol::REFLECT:
          np2++;
          break;
        case ChiralInversionProtocol::DO_NOT_INVERT:
          break;
        }
      }
      ln_choices += np2 * log(2.0);
    }
    int nproto_conf;
    SamplingStrategy strat;
    if (ln_choices < log(conf_input.getSystemTrialCount())) {
      nproto_conf = exp(ln_choices) + 0.99;
      strat = SamplingStrategy::FULL;
    }
    else if (ln_choices < log(static_cast<double>(conf_input.getSystemTrialCount()) * 32.0)) {
      nproto_conf = conf_input.getSystemTrialCount();
      strat = SamplingStrategy::LIMITED;
    }
    else {
      nproto_conf = conf_input.getSystemTrialCount();
      strat = SamplingStrategy::SPARSE;
    }
    
    // Create a series to hold the conformers resulting from the coarse-grained search.
    CoordinateSeries<float> cseries(sc.getTopologyPointer(i)->getAtomCount(), nproto_conf);
    const std::vector<int> top_cases = sc.getTopologicalCases(i);
    cseries.resize(0);
    int fc = 0;
    const int proto_conf_per_case = (nproto_conf + (ncases - 1)) / ncases;
    for (int j = 0; j < ncases; j++) {
      fc += proto_conf_per_case;
      fc = std::min(fc, nproto_conf);
      cseries.resize(fc, sc.getCoordinateReference(top_cases[j]));
    }
        
    // Manipulate this coordinate series using bond rotations as well as chiral inversions.
    // Compute the RMSD matrix and determine a set of diverse conformations.  Cull the results
    // to eliminate ring stabs or severe clashes between tertiary or quaternary atoms.
    fc = 0;
    const std::vector<RotatorGroup> rotatable_groups  = chemfe_list[i].getRotatableBondGroups();
    const std::vector<RotatorGroup> invertible_groups = chemfe_list[i].getChiralInversionGroups();
    CoordinateSeriesWriter<float> cseries_w = cseries.data();
    for (int j = 0; j < ncases; j++) {
      for (int k = 0; k < nrot_bond; k++) {
        for (int m = 0; m < nbond_rotations; m++) {
          const double rval = static_cast<double>(m) * stormm::symbols::twopi /
                              static_cast<double>(nbond_rotations);
          rotateAboutBond(cseries_w, fc, rotatable_groups[k].root_atom,
                          rotatable_groups[k].pivot_atom, rotatable_groups[k].rotatable_atoms,
                          rval);
          fc++;
        }
      }
    }

    // CHECK
#if 0
    ChemicalDetailsKit cdk = sc.getTopologyPointer(i)->getChemicalDetailsKit();
    for (int j = 0; j < ncases * bond_rotation_reps * nchiral_reps; j++) {
      for (int k = 0; k < j; k++) {
        const double this_rmsd = rmsd<float, float>(cseries_w, j, k, cdk, RmsdMethod::ALIGN_MASS,
                                                    0, cdk.natom);
        //printf("  %9.4lf", this_rmsd);
      }
      //printf("\n");
    }
#endif
    // END CHECK
  }
  
  // Create lists of PhaseSpace objects and topology pointers to show how to model each of them
  std::vector<PhaseSpace> ps_list; 
  std::vector<AtomGraph*> ag_list;
#if 0
  int conf_counter = 0;
  for (int i = 0; i < nsys; i++) {
    const int top_idx = sc.getTopologyIndex(i);
    const int nrot_bond = std::min(chemfe_list[top_idx].getRotatableBondCount(),
                                   conf_input.getRotatableBondLimit());
    const int ncopy = pow(nrot_bond, conf_input.getRotationSampleCount());
    for (int j = 0; j < ncopy; j++) {
      ag_list[conf_counter] = const_cast<AtomGraph*>(sc.getTopologyPointer(i));
      ps_list.push_back(PhaseSpace(sc.getCoordinateReference(i)));
      conf_counter++;
    }
  }
#endif
  return PhaseSpaceSynthesis(ps_list, ag_list);
}

} // namespace setup
} // namespace conf_app
