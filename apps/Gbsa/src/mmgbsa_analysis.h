// -*-c++-*-
#ifndef MMGBSA_ANALYSIS_H
#define MMGBSA_ANALYSIS_H

#include "../../../src/copyright.h"
#include "../../../src/Namelists/nml_files.h"
#include "../../../src/Namelists/nml_report.h"
#include "../../../src/Reporting/section_contents.h"
#include "mmgbsa_problem_set.h"

namespace mmgbsa {

using stormm::namelist::FilesControls;
using stormm::namelist::ReportControls;
using stormm::review::SectionContents;

/// \brief Summarize the chemical features of each ligand topology.
void addLigandTopologySummary(const MMGBSAProblemSet &sandbox, const FilesControls &ficon,
                              const ReportControls &repcon, std::vector<SectionContents> *outp);
  
} // namespace mmgbsa

#endif
