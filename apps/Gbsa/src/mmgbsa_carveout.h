// -*-c++-*-
#ifndef MMGBSA_CARVEOUT_H
#define MMGBSA_CARVEOUT_H

#include "../../../src/copyright.h"
#include "../../../src/Synthesis/systemcache.h"
#include "nml_mmgbsa.h"
#include "mmgbsa_problem_set.h"

namespace mmgbsa {

using stormm::topology::AtomGraph;
using stormm::trajectory::PhaseSpace;

// Given a receptor and a set of ligands, compute the "carevout" atoms of the receptor which lie
// near enough to any one of the ligands to qualify as mobile, or of interest to analysis.
//
// Arguments:
//   gbsacon:  Contains user-provided input on details of MM-GBSA calculations
//   sc:       A collection of all individual receptor and ligand systems needed for MM-GBSA
std::vector<int> receptorCarveOut(const MMGBSAControls &gbsacon, const SystemCache &sc);


} // namespace mmgbsa

#endif
