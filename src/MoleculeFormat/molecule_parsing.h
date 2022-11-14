// -*-c++-*-
#ifndef STORMM_MOLECULE_PARSING_H
#define STORMM_MOLECULE_PARSING_H

#include <string>
#include <vector>
#include "copyright.h"
#include "Chemistry/atommask.h"
#include "Chemistry/chemical_features.h"
#include "MoleculeFormat/mdlmol.h"
#include "Topology/atomgraph.h"
#include "Trajectory/coordinateframe.h"

namespace stormm {
namespace structure {

using chemistry::AtomMask;
using chemistry::ChemicalFeatures;
using structure::MdlMol;
using topology::AtomGraph;
using trajectory::CoordinateFrame;

/// \brief Derive an atom mask (a list of atoms) from information in the data lines of an SD file
///        data item.
///
/// Overloaded:
///   - Accept a pre-made CoordinateFrame object (for AtomMask production) along with the MDL MOL
///     object and pre-determined chemical features along with the topology
///   - Extract the coordinates from the MDL MOL object and the chemical features from the topology
///     as temporary resources
///
/// \param item_name  Name of the data item
/// \param molecule   The MDL MOL entry with appended data items
/// \param ag         Topology of the system of interest, provided for testing whether the text is
///                   a valid atom mask or a list of atom names or numbers
/// \param chemfe     Chemical features detected for the topology of interest
/// \{
AtomMask maskFromSdfDataItem(const std::string &item_name, const MdlMol &molecule,
                             const AtomGraph *ag, const ChemicalFeatures *chemfe,
                             const CoordinateFrame &cf);

AtomMask maskFromSdfDataItem(const std::string &item_name, const MdlMol &molecule,
                             const AtomGraph *ag)
/// \}

} // namespace structure
} // namespace stormm

#endif
