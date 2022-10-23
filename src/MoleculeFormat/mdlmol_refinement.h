// -*-c++-*-
#ifndef STORMM_MDLMOL_REFINEMENT_H
#define STORMM_MDLMOL_REFINEMENT_H

#include <vector>
#include <string>
#include "copyright.h"
#include "Namelists/nml_report.h"
#include "Restraints/restraint_apparatus.h"
#include "Synthesis/systemcache.h"
#include "Topology/atomgraph.h"
#include "mdlmol.h"

namespace stormm {
namespace structure {

using namelist::ReportControls;
using restraints::RestraintApparatus;
using synthesis::SystemCache;
using topology::AtomGraph;

/// \brief Add a series of user-defined data items to one or more MDL MOL entries.
///
/// Overloaded:
///   - Customize the data items of a single MDL MOL object
///   - Customize the data items of a list of MDL MOL objects
///
/// \param mol_entry
/// \param mol_entries
/// \param label
/// \param ag
/// \param ra
/// \param sysc
/// \param repcon
/// \{
void customizeDataItems(MdlMol *mol_entry, const std::string &label, const AtomGraph &ag,
                        const RestraintApparatus &ra, const ReportControls &repcon);

void customizeDataItems(std::vector<MdlMol> *mol_entries, const SystemCache &sysc,
                        const ReportControls &repcon);
/// \}

} // namespace structure
} // namespace stormm

#endif
