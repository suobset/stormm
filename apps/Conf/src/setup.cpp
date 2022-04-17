#include "../../../src/Constants/symbol_values.h"
#include "../../../src/Chemistry/chemistry_enumerators.h"
#include "../../../src/Chemistry/chemical_features.h"
#include "../../../src/Math/rounding.h"
#include "../../../src/Structure/isomerization.h"
#include "../../../src/Structure/rmsd.h"
#include "../../../src/Topology/atomgraph.h"
#include "../../../src/Topology/atomgraph_abstracts.h"
#include "../../../src/Trajectory/coordinate_series.h"
#include "user_settings.h"
#include "setup.h"

namespace conf_app {
namespace setup {

using user_input::ConformerControls;
using omni::chemistry::ChemicalFeatures;
using omni::chemistry::ChiralInversionProtocol;
using omni::chemistry::MapRotatableGroups;
using omni::chemistry::RotatorGroup;
using omni::math::roundUp;
using omni::structure::rmsd;
using omni::structure::RmsdMethod;
using omni::structure::rotateAboutBond;
using omni::topology::AtomGraph;
using omni::topology::ChemicalDetailsKit;
using omni::trajectory::PhaseSpace;
using omni::trajectory::CoordinateFrame;
using omni::trajectory::CoordinateFrameReader;
using omni::trajectory::CoordinateSeries;
using omni::trajectory::CoordinateSeriesWriter;
  
//-------------------------------------------------------------------------------------------------
PhaseSpaceSynthesis expandConformers(const UserSettings &ui, const SystemCache &sc,
                                     StopWatch *tm) {
  const ConformerControls conf_input = ui.getConformerNamelistInfo();

  // Count the expanded number of systems
  const int ntop = sc.getTopologyCount();
  std::vector<ChemicalFeatures> chemfe_list;
  chemfe_list.reserve(ntop);
  for (int i = 0; i < ntop; i++) {
    const int example_system_idx = sc.getCoordinateExample(i);
    chemfe_list.emplace_back(sc.getTopologyPointer(example_system_idx),
                             sc.getCoordinateReference(example_system_idx),
                             MapRotatableGroups::YES);
  }

  // Loop over all systems, grouping those with the same topology into a coherent group of
  // proto-conformers for coarse-grained sampling of rotatable bonds and chiral centers.
  const int nsys = sc.getSystemCount();
  const int nbond_rotations = conf_input.getRotationSampleCount();
  int nconformer = 0;
  for (int i = 0; i < ntop; i++) {

    // Determine the number of conformers in the coarse-grained search.  This includes bond
    // rotation, chiral inversions, and multiple replicas of the compound in the list of starting
    // coordinates.
    const int nrot_bond = std::min(chemfe_list[i].getRotatableBondCount(),
                                   conf_input.getRotatableBondLimit());
    const int ncases = sc.getTopologyCaseCount(i);
    const int bond_rotation_reps = nrot_bond * nbond_rotations;
    int nproto_conf = ncases * bond_rotation_reps;
    int nchiral_reps = 1;
    if (conf_input.sampleChirality()) {
      const std::vector<ChiralInversionProtocol> chiral_protocols;
      const int nchiral = chemfe_list[i].getChiralCenterCount();
      for (int j = 0; j < nchiral; j++) {
        switch (chiral_protocols[j]) {
        case ChiralInversionProtocol::ROTATE:
        case ChiralInversionProtocol::REFLECT:
          nproto_conf *= 2;
          nchiral_reps *= 2;
          break;
        case ChiralInversionProtocol::DO_NOT_INVERT:
          break;
        }
        if (nproto_conf > conf_input.getSystemTrialCount()) {
          break;
        }
      }
    }
    nproto_conf = std::min(conf_input.getSystemTrialCount(), nproto_conf);

    // Create a series to hold the conformers resulting from the coarse-grained search.
    CoordinateSeries<float> cseries(sc.getCoordinateReference(sc.getCoordinateExample(i)),
                                    nproto_conf);
    const std::vector<int> top_cases = sc.getTopologicalCases(i);
    cseries.resize(0);
    int fc = 0;
    const int proto_conf_per_case = (nproto_conf + (ncases - 1)) / ncases;
    for (int j = 0; j < ncases; j++) {
      fc += proto_conf_per_case;
      fc = std::min(fc, nproto_conf);
      cseries.resize(fc, sc.getCoordinateReference(top_cases[j]));
    }

    // Determine how to sample bond rotations and chiral inversions.  If there are enough slots to
    // sample each permutation, sample them all.  Otherwise, perform a round-robin sampling,
    // explicitly sampling each rotatable bond or chiral inversion in series, with random sampling
    // of other rotatable bonds and chiral centers.
    
    // Manipulate this coordinate series using bond rotations as well as chiral inversions.
    // Compute the RMSD matrix and determine a set of diverse conformations.  Cull the results
    // to eliminate ring stabs or severe clashes between tetiary or quaternary atoms.
    fc = 0;
    const std::vector<RotatorGroup> rotatable_groups  = chemfe_list[i].getRotatableBondGroups();
    const std::vector<RotatorGroup> invertible_groups = chemfe_list[i].getChiralInversionGroups();
    CoordinateSeriesWriter<float> cseries_w = cseries.data();
    for (int j = 0; j < ncases; j++) {
      for (int k = 0; k < nrot_bond; k++) {
        for (int m = 0; m < nbond_rotations; m++) {
          const double rval = static_cast<double>(m) * omni::symbols::twopi /
                              static_cast<double>(nbond_rotations);
          rotateAboutBond(cseries_w, fc, rotatable_groups[k].root_atom,
                          rotatable_groups[k].pivot_atom, rotatable_groups[k].rotatable_atoms,
                          rval);
          fc++;
        }
      }
    }

    // CHECK
    ChemicalDetailsKit cdk = sc.getTopologyPointer(i)->getChemicalDetailsKit();
    for (int j = 0; j < ncases * bond_rotation_reps * nchiral_reps; j++) {
      const CoordinateFrame frm_j = cseries.exportFrame(j);
      const CoordinateFrameReader frm_jr = frm_j.data();
      for (int k = 0; k < j; k++) {
        const CoordinateFrame frm_k = cseries.exportFrame(k);
        const CoordinateFrameReader frm_kr = frm_k.data();
        printf("  %9.4lf", rmsd(frm_jr, frm_kr, cdk, RmsdMethod::ALIGN_MASS));
      }
      printf("\n");
    }
    // END CHECK
  }
  tm->assignTime(2);

  // CHECK
  exit(1);
  // END CHECK
  
  // Create lists of PhaseSpace objects and topology pointers to show how to model each of them
  std::vector<PhaseSpace> ps_list;
  ps_list.reserve(nconformer);
  std::vector<AtomGraph*> ag_list(nconformer);
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
  return PhaseSpaceSynthesis(ps_list, ag_list);
}

} // namespace setup
} // namespace conf_app
