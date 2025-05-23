// -*-c++-*-
#ifndef MMGBSA_TESTING_H
#define MMGBSA_TESTING_H

namespace mmgbsa {

//-------------------------------------------------------------------------------------------------
// Run unit tests of the program-specific code.
//-------------------------------------------------------------------------------------------------
int runUnitTests();
  
//-------------------------------------------------------------------------------------------------
// Run regression tests with some inputs from the main STORMM repository.
//-------------------------------------------------------------------------------------------------
int runRegressionTests();

} // namespace mmgbsa

#endif
