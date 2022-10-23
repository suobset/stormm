#include "copyright.h"
#include "mdlmol_refinement.h"

namespace stormm {
namespace structure {

//-------------------------------------------------------------------------------------------------
void customizeDataItems(MdlMol *mol_entry, const std::string &label, const AtomGraph &ag,
                        const RestraintApparatus &ra, const ReportControls &repcon) {
  const int nreq = repcon.getSDFileDataRequestCount();
  const std::vector<MolObjDataRequest>& reqs = repcon.getSDFileDataRequests();
  for (int i = 0; i < nreq; i++) {
    mol_entry->addDataItem(reqs[i], ag, ra);
  }
}

//-------------------------------------------------------------------------------------------------
void customizeDataItems(std::vector<MdlMol> *mol_entries, const SystemCache &sysc,
                        const ReportControls &repcon) {
  const int nmol = mol_entries->size();
  MdlMol* ment_data = mol_entries->data();
  for (int i = 0; i < nmol; i++) {
    customizeDataItems(&ment_data[i], sysc.getSystemLabel(i),
                       sysc.getSystemTopologyReference(i), sysc.getRestraintReference(i),
                       repcon);
  }
}

} // namespace structure
} // namespace stormm
