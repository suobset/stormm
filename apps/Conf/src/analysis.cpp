#include "copyright.h"
#include "../../../src/Accelerator/hybrid.h"
#include "../../../src/FileManagement/file_enumerators.h"
#include "../../../src/FileManagement/file_listing.h"
#include "../../../src/FileManagement/file_util.h"
#include "../../../src/Math/summation.h"
#include "../../../src/Parsing/parse.h"
#include "../../../src/Structure/clash_detection.h"
#include "../../../src/Structure/rmsd.h"
#include "../../../src/Structure/rmsd_plan.h"
#include "../../../src/Trajectory/trajectory_enumerators.h"
#include "analysis.h"

namespace conf_app {
namespace analysis {

using stormm::card::Hybrid;
using stormm::diskutil::getBaseName;
using stormm::diskutil::PrintSituation;
using stormm::diskutil::substituteNameExtension;
using stormm::diskutil::splitPath;
using stormm::math::sum;
using stormm::parse::findStringInVector;
using stormm::structure::detectClash;
using stormm::structure::RMSDPlan;
using stormm::structure::MdlMolVersion;
using stormm::synthesis::PsSynthesisWriter;
using stormm::topology::AtomGraph;
using stormm::trajectory::CoordinateFileKind;
using stormm::trajectory::TrajectoryFusion;
  
//-------------------------------------------------------------------------------------------------

//-------------------------------------------------------------------------------------------------
std::vector<int> filterMinimizedStructures(const PhaseSpaceSynthesis &poly_ps,
                                           const std::vector<StaticExclusionMask> &poly_ps_masks,
                                           const ScoreCard &emin,
                                           const ConformerControls &confcon) {
  Condensate poly_ps_snapshot(poly_ps);
  const std::vector<double> efinal = emin.reportTotalEnergies();
  const PsSynthesisReader poly_psr = poly_ps.data();
  const std::vector<AtomGraph*> unique_ag = poly_ps.getUniqueTopologies();
  std::vector<bool> clashing(poly_psr.system_count, false);
  std::vector<bool> rmsd_eliminated(poly_psr.system_count, false);
  int nbad = 0;
  const int nseek = confcon.getFinalStateCount();
  const double rmsd_tol = confcon.getRMSDTolerance();

  // Seek out clashes and mark those structures as invalid.
  for (int i = 0; i < poly_psr.unique_topology_count; i++) {
    const int nbad_so_far = nbad;
    for (int j = poly_psr.common_ag_bounds[i]; j < poly_psr.common_ag_bounds[i + 1]; j++) {
      const int sys_ij = poly_psr.common_ag_list[j];
      ClashReport clrep(unique_ag[i]);
      if (detectClash(&poly_ps_snapshot, sys_ij, unique_ag[i], &poly_ps_masks[i],
                      default_minimize_clash_r0, default_minimize_clash_ratio, &clrep)) {
        clashing[sys_ij] = true;
        nbad++;
      }
    }
  }
  
  // Find the lowest energy structure, cull those that are within a certain positional rmsd,
  // catalog the result, and repeat until the desired number of unique conformations has been
  // reached or all conformations have been eliminated.
  const RMSDPlan poly_rplan(poly_ps);
  Hybrid<double> rmsd_values(poly_rplan.getReferenceRMSDSize(), "rmsd_val");
  double* rmsd_val_ptr = rmsd_values.data();
  Hybrid<int> reference_frames(poly_psr.unique_topology_count, "ref_frame_idx");
  int* ref_frm_ptr = reference_frames.data();
  std::vector<int> successes(poly_psr.unique_topology_count * nseek, -1);
  std::vector<int> successes_per_system(poly_psr.unique_topology_count, 0);
  for (int rcon = 0; rcon < nseek; rcon++) {
    for (int i = 0; i < poly_psr.unique_topology_count; i++) {
      bool init_e_best = true;
      double e_best;
      int best_loc;
      for (int j = poly_psr.common_ag_bounds[i]; j < poly_psr.common_ag_bounds[i + 1]; j++) {
        const int sys_ij = poly_psr.common_ag_list[j];
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
        ref_frm_ptr[i] = 0;
      }
    }
    

    // Compute RMSD values to the lowest-energy structures
    rmsd(poly_rplan, poly_ps, poly_ps_snapshot, reference_frames, &rmsd_values);

    // Eliminate structures with RMSD too close to the best structure found
    for (int i = 0; i < poly_psr.unique_topology_count; i++) {
      int irmsd_track = poly_rplan.getReferenceRMSDStart(i);
      for (int j = poly_psr.common_ag_bounds[i]; j < poly_psr.common_ag_bounds[i + 1]; j++) {
        const int sys_ij = poly_psr.common_ag_list[j];
        if (rmsd_val_ptr[irmsd_track] < rmsd_tol) {
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
                  const ScoreCard &emin, const SystemCache &sc,
                  const std::vector<MdlMol> &sdf_recovery, const FilesControls &fcon) {
  const PsSynthesisReader poly_psr = poly_ps.data(); 
  const int ntop = poly_ps.getUniqueTopologyCount();
  std::vector<std::vector<int>> like_confs(ntop);
  const int nconfs = best_confs.size();
  for (int i = 0; i < nconfs; i++) {
    const int conf_idx = best_confs[i];
    const int top_idx = poly_psr.unique_ag_idx[conf_idx];
    like_confs[top_idx].push_back(conf_idx);
  }
  const std::vector<AtomGraph*>& unique_topologies = poly_ps.getUniqueTopologies();
  std::vector<std::string> printed_thus_far;
  for (int i = 0; i < ntop; i++) {

    // Loop over all conformations sharing this topology.
    for (size_t j = 0; j < like_confs[i].size(); j++) {
      std::string before, after;
      const std::string& ij_label = poly_ps.getSystemLabel(like_confs[i][j]);
      const int cache_idx = sc.getFirstMatchingSystemIndex(unique_topologies[i], ij_label);
      splitPath(sc.getTrajectoryName(cache_idx), &before, &after);
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
            fname = before + "." + after;
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
            fname = before + "_" + std::to_string(j - poly_psr.common_ag_bounds[i]) + "." + after;
            break;
          }
          if (sc.getTrajectoryKind(cache_idx) == CoordinateFileKind::SDF) {
            MdlMol tmdl = sdf_recovery[cache_idx];
            tmdl.impartCoordinates(poly_ps.exportCoordinates(like_confs[i][j]));
            tmdl.writeMdl(fname, MdlMolVersion::V2000, pr_protocol);
            tmdl.writeDataItems(fname, PrintSituation::APPEND);
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
        fname = before + "_" + std::to_string(j - poly_psr.common_ag_bounds[i]) + "." + after;
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
