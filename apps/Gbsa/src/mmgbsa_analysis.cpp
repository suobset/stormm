#include <cstdio>
#include "../../../src/copyright.h"
#include "../../../src/Chemistry/chemical_features.h"
#include "../../../src/FileManagement/file_listing.h"
#include "../../../src/Math/summation.h"
#include "../../../src/Parsing/parse.h"
#include "../../../src/Parsing/parsing_enumerators.h"
#include "../../../src/Topology/atomgraph_abstracts.h"
#include "mmgbsa_analysis.h"

namespace mmgbsa {

using stormm::chemistry::ChemicalFeatures;
using stormm::diskutil::extractCommonPaths;
using stormm::diskutil::listCommonPaths;
using stormm::namelist::FilesControls;
using stormm::namelist::ReportControls;
using stormm::parse::JustifyText;
using stormm::parse::minimalRealFormat;
using stormm::review::ReportTable;
using stormm::stmath::sum;
using stormm::topology::AtomGraph;
using stormm::topology::ChemicalDetailsKit;
using stormm::topology::NonbondedKit;
  
//-------------------------------------------------------------------------------------------------
void addLigandTopologySummary(const MMGBSAProblemSet &sandbox, const FilesControls &ficon,
                              const ReportControls &repcon, std::vector<SectionContents> *outp) {
  const int ncomp = sandbox.getComplexCount();
  const int nrcpt = sandbox.getReceptorStates();
  const SystemCache *sysc = sandbox.getCachePointer();
  const int ntop  = sysc->getTopologyCount();

  // Compute the number initial poses for each ligand topology.
  std::vector<int> top_counts(ntop, 0);
  for (int i = 0; i < ncomp; i++) {
    const int top_idx = sandbox.getLigandTopologyCacheIndex(nrcpt + i);
    top_counts[top_idx] += 1;
  }
  
  // Loop over all unique ligand topologies and pull out the rotatable bond counts.
  int nuniq_lig = 0;
  for (int i = 0; i < ntop; i++) {
    nuniq_lig += (top_counts[i] > 0);
  }
  std::vector<std::string> unique_lig_top_names;
  std::vector<int> unique_lig_top_counts, unique_lig_top_atoms, unique_lig_top_rotbonds;
  std::vector<double> unique_lig_top_masses, unique_lig_top_charges;
  unique_lig_top_names.reserve(nuniq_lig);
  unique_lig_top_counts.reserve(nuniq_lig);
  unique_lig_top_atoms.reserve(nuniq_lig);
  unique_lig_top_rotbonds.reserve(nuniq_lig);
  unique_lig_top_masses.reserve(nuniq_lig);
  unique_lig_top_charges.reserve(nuniq_lig);
  for (int i = 0; i < ntop; i++) {
    if (top_counts[i] > 0) {
      const AtomGraph *ag_i = sysc->getTopologyPointer(i);
      const int cache_sys_idx = sysc->getSystemExampleIndex(i);
      const ChemicalFeatures *ft_i = sysc->getFeaturesPointer(cache_sys_idx);
      unique_lig_top_names.push_back(ag_i->getFileName());
      unique_lig_top_counts.push_back(top_counts[i]);
      unique_lig_top_atoms.push_back(ag_i->getAtomCount());
      unique_lig_top_rotbonds.push_back(ft_i->getRotatableBondCount());
      const ChemicalDetailsKit cdk = ag_i->getChemicalDetailsKit();
      double tmass = 0.0;
      for (int i = 0; i < cdk.natom; i++) {
        tmass += cdk.masses[i];
      }
      unique_lig_top_masses[i] = tmass;
      const NonbondedKit<double> nbk = ag_i->getDoublePrecisionNonbondedKit();
      unique_lig_top_charges[i] = sum<double>(nbk.charge, nbk.natom);
    }
  }

  // Create a table of the notable details for each unique ligand topology.
  const std::vector<std::string> common_paths = extractCommonPaths(&unique_lig_top_names, 8, 2);
  const char comment_symbol = commentSymbol(repcon.getOutputSyntax());
  std::vector<std::string> ligand_top_data(nuniq_lig * 6);
  for (int i = 0; i < nuniq_lig; i++) {
    ligand_top_data[i                  ] = std::to_string(unique_lig_top_counts[i]);
    ligand_top_data[i +      nuniq_lig ] = std::to_string(unique_lig_top_atoms[i]);
    ligand_top_data[i + (2 * nuniq_lig)] = std::to_string(unique_lig_top_rotbonds[i]);
    ligand_top_data[i + (3 * nuniq_lig)] = minimalRealFormat(unique_lig_top_masses[i], 0.01, true);
    ligand_top_data[i + (4 * nuniq_lig)] = minimalRealFormat(unique_lig_top_charges[i], 0.1, true);
    ligand_top_data[i + (5 * nuniq_lig)] = std::string(1, comment_symbol) + " " +
                                                  unique_lig_top_names[i];
  }
  const std::vector<std::string> headers = { "Poses", "Atom Count", "Rotatable Bonds",
                                             "Molecular Weight", "Net Charge", "File of Origin" };
  ReportTable lig_summary_tab(ligand_top_data, headers, "ligand_data", repcon.getReportFileWidth(),
                              { JustifyText::RIGHT, JustifyText::RIGHT, JustifyText::RIGHT,
                                JustifyText::LEFT });
  lig_summary_tab.unprotectContent();
  outp->emplace_back("Ligand Summary", ficon.getReportFile());
  outp->back().addNarration("A summary of critical features of each ligand follows.  Molecular "
                            "weights and net charges are presented in atomic units (Daltons and "
                            "proton charges).");
  if (common_paths.size() > 0) {
    outp->back().addNarration("In the table that follows, file names are abbreviated.  The "
                              "following key illustrates the paths:\n" +
                              listCommonPaths(common_paths));
  }
  outp->back().addTable(lig_summary_tab);
}
  
} // namespace mmgbsa
