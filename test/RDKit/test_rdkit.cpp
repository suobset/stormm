#include "copyright.h"

#include "../../src/Reporting/summary_file.h"
#include "../../src/UnitTesting/unit_test.h"

#include <GraphMol/RDKitBase.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

using stormm::errors::rtWarn;
using stormm::review::stormmSplash;
using namespace stormm::testing;

void testCreateRDKitMol() {
  
  std::shared_ptr<RDKit::ROMol> implicit_mol( RDKit::SmilesToMol( "C" ) );

  auto molecule_was_created = ( implicit_mol != nullptr);
  check(molecule_was_created, "RDKit ROMol could not be created from SMILES");

  check(implicit_mol->getNumAtoms(), RelationalOperator::EQUAL, 1, "RDKit ROMol has wrong number of implicit atoms");

  std::shared_ptr<RDKit::ROMol> explicit_mol( RDKit::MolOps::addHs(*implicit_mol) );
  check(explicit_mol->getNumAtoms(), RelationalOperator::EQUAL, 5, "RDKit ROMol has wrong number of explicit atoms");
}

//-------------------------------------------------------------------------------------------------
// main
//-------------------------------------------------------------------------------------------------
int main(const int argc, const char* argv[]) {

  // Some baseline initialization
  TestEnvironment oe(argc, argv);
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmSplash();
  }

  section("Long long int <=> float2 conversion");
  section(1);

  testCreateRDKitMol();

  // Print results
  printTestSummary(oe.getVerbosity());

  return countGlobalTestFailures();
}
