#include "copyright.h"
#include "../../../src/Accelerator/hybrid.h"
#include "../../../src/Analysis/comparison_guide.h"
#include "../../../src/FileManagement/file_enumerators.h"
#include "../../../src/FileManagement/file_listing.h"
#include "../../../src/FileManagement/file_util.h"
#include "../../../src/Math/summation.h"
#include "../../../src/MoleculeFormat/mdlmol_refinement.h"
#include "../../../src/Numerics/split_fixed_precision.h"
#include "../../../src/Parsing/parse.h"
#include "../../../src/Structure/clash_detection.h"
#include "../../../src/Structure/rmsd.h"
#include "../../../src/Structure/rmsd_plan.h"
#include "../../../src/Structure/structure_enumerators.h"
#include "../../../src/Synthesis/condensate.h"
#include "../../../src/Trajectory/trajectory_enumerators.h"
#include "analysis.h"

namespace conf_app {
namespace analysis {

using stormm::analysis::ComparisonGuide;
using stormm::card::Hybrid;
using stormm::diskutil::getBaseName;
using stormm::diskutil::PrintSituation;
using stormm::diskutil::substituteNameExtension;
using stormm::diskutil::splitPath;
using stormm::numerics::hostInt95ToDouble;
using stormm::parse::findStringInVector;
using stormm::stmath::sum;
using stormm::structure::customizeDataItems;
using stormm::structure::detectClash;
using stormm::structure::RMSDPlan;
using stormm::structure::RMSDTask;
using stormm::structure::MdlMolVersion;
using stormm::structure::updateDataItemReadouts;
using stormm::synthesis::CondensateWriter;
using stormm::synthesis::PsSynthesisWriter;
using stormm::synthesis::SynthesisMapReader;
using stormm::synthesis::translateSystemGrouping;
using stormm::topology::AtomGraph;
using stormm::trajectory::CoordinateFileKind;
using stormm::trajectory::TrajectoryFusion;

//-------------------------------------------------------------------------------------------------
std::vector<int> filterMinimizedStructures(const PhaseSpaceSynthesis &poly_ps,
                                           const std::vector<StaticExclusionMask> &poly_ps_masks,
                                           const SystemCache &sc, const SynthesisCacheMap &scmap,
                                           const ScoreCard &emin, const ConformerControls &confcon,
                                           const GpuDetails &gpu) {
  Condensate cdns(poly_ps);
  ComparisonGuide cg(cdns, scmap);
  const std::vector<double> efinal = emin.reportTotalEnergies();
  const PsSynthesisReader poly_psr = poly_ps.data();
  const std::vector<AtomGraph*> unique_ag = poly_ps.getUniqueTopologies();
  std::vector<bool> clashing(poly_psr.system_count, false);
  std::vector<bool> rmsd_eliminated(poly_psr.system_count, false);
  const int nseek = confcon.getFinalStateCount();
  const double rmsd_tol = confcon.getRMSDTolerance();

  // Determine the grouping method
  const SystemGrouping group_method = confcon.getGroupingMethod();
  
  // Seek out clashes and mark those structures as invalid.
  int nbad = 0;
  for (int i = 0; i < poly_psr.system_count; i++) {
    const int iag_idx = poly_psr.unique_ag_idx[i];
    ClashReport clrep(unique_ag[iag_idx]);
    if (detectClash<double>(&cdns, i, unique_ag[iag_idx], &poly_ps_masks[iag_idx],
                            default_minimize_clash_r0, default_minimize_clash_ratio, &clrep)) {
      clashing[i] = true;
      nbad++;
    }
  }
  
  // Find the lowest energy structure, cull those that are within a certain positional rmsd,
  // catalog the result, and repeat until the desired number of unique conformations has been
  // reached or all conformations have been eliminated.
  const RMSDPlan poly_rplan(poly_ps, sc, scmap);
  const SynthesisMapReader scmapr = scmap.data();
  Hybrid<double> rmsd_values(cg.getAllToReferenceOutputSize(), "rmsd_val");
  double* rmsd_val_ptr = rmsd_values.data();
  const int n_result_groups = cg.getPartitionCount(group_method);
  Hybrid<int> reference_frames(n_result_groups, "ref_frame_idx");
  int* ref_frm_ptr = reference_frames.data();
  std::vector<int> successes(n_result_groups * nseek, -1);
  std::vector<int> successes_per_system(n_result_groups, 0);

  // However the results are to be grouped, apply the following cycle: find the best available
  // conformer (lowest energy, not in a clash state, not eliminated by proximity to some previous
  // best conformer selection), eliminate all other conformers within a selected proximity from
  // further consideration.  If no conformers are available within a partition, the search will
  // terminate early.
  for (int rcon = 0; rcon < nseek; rcon++) {
    std::vector<bool> new_conformer_found(n_result_groups, true);
    for (int i = 0; i < n_result_groups; i++) {
      bool init_e_best = true;
      double e_best;
      int best_loc;

      // Loop over all members of the partition of interest.
      int* group_index_ptr;
      int group_llim, group_hlim;
      switch (group_method) {
      case SystemGrouping::SOURCE:
        group_llim = scmapr.csystem_bounds[i];
        group_hlim = scmapr.csystem_bounds[i + 1];
        break;
      case SystemGrouping::TOPOLOGY:
        group_llim = scmapr.ctopol_bounds[i];
        group_hlim = scmapr.ctopol_bounds[i + 1];
        break;
      case SystemGrouping::LABEL:
        group_llim = scmapr.clabel_bounds[i];
        group_hlim = scmapr.clabel_bounds[i + 1];
        break;
      }
      for (int j = group_llim; j < group_hlim; j++) {
        int sys_ij;
        switch (group_method) {
        case SystemGrouping::SOURCE:
          sys_ij = scmapr.csystem_proj[j];
          break;
        case SystemGrouping::TOPOLOGY:
          sys_ij = scmapr.ctopol_proj[j];
          break;
        case SystemGrouping::LABEL:
          sys_ij = scmapr.clabel_proj[j];
          break;
        }
        if (clashing[sys_ij] || rmsd_eliminated[sys_ij]) {
          continue;
        }
        if (init_e_best || e_best > efinal[sys_ij]) {
          e_best = efinal[sys_ij];
          best_loc = sys_ij;
          init_e_best = false;
        }
      }
      if (init_e_best == false) {
        ref_frm_ptr[i] = best_loc;
        const int ni_success = successes_per_system[i]; 
        successes[(i * nseek) + ni_success] = best_loc;
        successes_per_system[i] = ni_success + 1;
      }
      else {
        ref_frm_ptr[i] = group_llim;
        new_conformer_found[i] = false;
      }
    }

    // Compute RMSD values to the lowest-energy structures
    rmsd(cg, poly_rplan, poly_ps, cdns, reference_frames, &rmsd_values, group_method);

    // Eliminate structures with RMSD too close to the best structure found
    for (int i = 0; i < n_result_groups; i++) {
      if (new_conformer_found[i] == false) {
        continue;
      }
      int* group_index_ptr;
      int group_llim, group_hlim;
      switch (group_method) {
      case SystemGrouping::SOURCE:
        group_llim = scmapr.csystem_bounds[i];
        group_hlim = scmapr.csystem_bounds[i + 1];
        break;
      case SystemGrouping::TOPOLOGY:
        group_llim = scmapr.ctopol_bounds[i];
        group_hlim = scmapr.ctopol_bounds[i + 1];
        break;
      case SystemGrouping::LABEL:
        group_llim = scmapr.clabel_bounds[i];
        group_hlim = scmapr.clabel_bounds[i + 1];
        break;
      }
      int irmsd_track = cg.getAllToOneResultOffset(i, group_method);
      for (int j = group_llim; j < group_hlim; j++) {
        if (rmsd_val_ptr[irmsd_track] < rmsd_tol) {
          int sys_ij;
          switch (group_method) {
          case SystemGrouping::SOURCE:
            sys_ij = scmapr.csystem_proj[j];
            break;
          case SystemGrouping::TOPOLOGY:
            sys_ij = scmapr.ctopol_proj[j];
            break;
          case SystemGrouping::LABEL:
            sys_ij = scmapr.clabel_proj[j];
            break;
          }
          rmsd_eliminated[sys_ij] = true;
        }
        irmsd_track++;
      }
    }
  }
  
  // Produce a vector of the conformations remaining
  std::vector<int> result(sum<int>(successes_per_system));
  int nsucc = 0;
  for (int i = 0; i < poly_psr.unique_topology_count; i++) {
    for (int j = 0; j < successes_per_system[i]; j++) {
      result[nsucc] = successes[(i * nseek) + j];
      nsucc++;
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
void printResults(const PhaseSpaceSynthesis &poly_ps, const std::vector<int> &best_confs,
                  const ScoreCard &emin, const SystemCache &sc, const SynthesisCacheMap &scmap,
                  const ConformerControls &confcon, const ReportControls &repcon) {

  // Determine the grouping method
  const SystemGrouping group_method = confcon.getGroupingMethod();
  const int n_result_groups = scmap.getPartitionCount(group_method);
  const PsSynthesisReader poly_psr = poly_ps.data();
  std::vector<std::vector<int>> like_confs(n_result_groups);
  const int nconfs = best_confs.size();
  for (int i = 0; i < nconfs; i++) {
    const int conf_idx = best_confs[i];
    int result_idx;
    switch (group_method) {
    case SystemGrouping::SOURCE:
      result_idx = scmap.getSystemCacheIndex(conf_idx);
      break;
    case SystemGrouping::TOPOLOGY:
      result_idx = poly_psr.unique_ag_idx[conf_idx];
      break;
    case SystemGrouping::LABEL:
      result_idx = sc.getSystemLabelIndex(scmap.getSystemCacheIndex(conf_idx));
      break;
    }
    like_confs[result_idx].push_back(conf_idx);
  }
  const std::vector<AtomGraph*>& unique_topologies = poly_ps.getUniqueTopologies();
  std::vector<std::string> printed_thus_far;
  for (int i = 0; i < n_result_groups; i++) {
    
    // Loop over all conformations sharing this topology.
    for (size_t j = 0; j < like_confs[i].size(); j++) {
      std::string before, after;
      const int cache_idx = scmap.getSystemCacheIndex(like_confs[i][j]);
      const std::string& ij_label = sc.getSystemLabel(cache_idx);
      splitPath(sc.getTrajectoryName(cache_idx), &before, &after);
      const std::string dot_connector = (after.size() > 0) ? "." : "";
      std::string fname;
      switch (sc.getTrajectoryKind(cache_idx)) {
      case CoordinateFileKind::AMBER_CRD:
      case CoordinateFileKind::AMBER_NETCDF:
      case CoordinateFileKind::SDF:
        {
          PrintSituation pr_protocol;
          switch (sc.getTrajectoryFusionProtocol()) {
          case TrajectoryFusion::ON:
          case TrajectoryFusion::AUTO:
            fname = before + dot_connector + after;
            if ((printed_thus_far.size() > 0 && fname == printed_thus_far.back()) ||
                findStringInVector(printed_thus_far, fname) < printed_thus_far.size()) {
              pr_protocol = PrintSituation::APPEND;
            }
            else {
              pr_protocol = sc.getPrintingProtocol();
              printed_thus_far.push_back(fname);
            }
            break;
          case TrajectoryFusion::OFF:
            pr_protocol = sc.getPrintingProtocol();
            if (like_confs[i].size() <= 1) {
              fname = before + dot_connector + after;
            }
            else {
              fname = before + "_" + std::to_string(j) + dot_connector + after;
            }
            break;
          }
          if (sc.getTrajectoryKind(cache_idx) == CoordinateFileKind::SDF) {

            // Build the SD file de novo, or use the existing SD file as a template
            if (sc.getStructureDataEntryPointer(cache_idx)->getAtomCount() == 0) {
              MdlMol tmdl(sc.getFeaturesPointer(cache_idx), poly_psr, like_confs[i][j], 0);
              customizeDataItems(&tmdl, ij_label, *unique_topologies[i],
                                 sc.getRestraints(cache_idx), repcon);
              updateDataItemReadouts(&tmdl, sc, emin, like_confs[i][j]);
              tmdl.writeMdl(fname, MdlMolVersion::V2000, pr_protocol);
              tmdl.writeDataItems(fname, PrintSituation::APPEND);
            }
            else {
              MdlMol tmdl = sc.getStructureDataEntry(cache_idx);
              tmdl.impartCoordinates(poly_ps.exportCoordinates(like_confs[i][j]));
              customizeDataItems(&tmdl, ij_label, *unique_topologies[i],
                                 sc.getRestraints(cache_idx), repcon);
              updateDataItemReadouts(&tmdl, sc, emin, like_confs[i][j]);
              tmdl.writeMdl(fname, MdlMolVersion::V2000, pr_protocol);
              tmdl.writeDataItems(fname, PrintSituation::APPEND);
            }
          }
          else {
            poly_ps.printTrajectory(std::vector<int>(1, like_confs[i][j]), fname, 0.0,
                                    sc.getTrajectoryKind(cache_idx), pr_protocol);
          }
        }
        break;
      case CoordinateFileKind::AMBER_INPCRD:
      case CoordinateFileKind::AMBER_ASCII_RST:
      case CoordinateFileKind::AMBER_NETCDF_RST:
        if (like_confs[i].size() <= 1) {
          fname = before + "." + after;
        }
        else {
          fname = before + "_" + std::to_string(j) + "." + after;
        }
        poly_ps.printTrajectory(std::vector<int>(1, like_confs[i][j]), fname, 0.0,
                                sc.getTrajectoryKind(cache_idx), sc.getPrintingProtocol());
        break;
      case CoordinateFileKind::UNKNOWN:
        break;
      }
    }
  }
}
  
} // namespace analysis
} // namespace conf_app
