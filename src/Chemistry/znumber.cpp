#include "copyright.h"
#include "UnitTesting/unit_test.h"
#include "znumber.h"
#include "periodic_table.h"

namespace stormm {
namespace chemistry {

using testing::Approx;
using testing::ComparisonType;  
  
//-------------------------------------------------------------------------------------------------
std::vector<int> massToZNumber(const std::vector<double> &masses) {

  // Create a vector of approximate masses for comparisons
  std::vector<Approx> natural_masses;
  for (int i = 0; i < element_maximum_count; i++) {
    natural_masses.push_back(Approx(elemental_masses[i], 0.01));
  }

  // Loop over all atoms in the system
  const int natom = masses.size();
  std::vector<int> result(natom, -1);
  for (int i = 0; i < natom; i++) {
    const double mass_i = masses[i];
    for (int j = 0; j < element_maximum_count; j++) {
      if (natural_masses[j].test(mass_i)) {
	result[i] = j;
	break;
      }
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
char2 zNumberToSymbol(const int atomic_number) {
  if (atomic_number < 0) {
    rtWarn("Atomic number " + std::to_string(atomic_number) + " is invalid and will be "
           "indicated by symbol XX.", "zNumberToSymbol");
  }
  else if (atomic_number > element_maximum_count) {
    rtWarn("Atomic number " + std::to_string(atomic_number) + " is beyond the scope of "
           "stable elements covered by STORMM and will be indicated by symbol XX.",
           "zNumberToSymbol");      
  }
  return elemental_symbols[atomic_number];
}

//-------------------------------------------------------------------------------------------------
std::vector<char2> zNumberToSymbol(const std::vector<int> &atomic_numbers) {
  const int natom = atomic_numbers.size();
  std::vector<char2> symbs(natom);
  for (int i = 0; i < natom; i++) {
    symbs[i] = zNumberToSymbol(atomic_numbers[i]);
  }
  return symbs;
}

//-------------------------------------------------------------------------------------------------
std::vector<int> symbolToZNumber(const std::vector<char2> &atomic_symbols) {
  const int natom = atomic_symbols.size();
  std::vector<int> result(natom, -1);
  for (int i = 0; i < natom; i++) {
    const char2 tmps = atomic_symbols[i];
    for (int j = 0; j < element_maximum_count; j++) {
      if (tmps.x == elemental_symbols[j].x && tmps.y == elemental_symbols[j].y) {
        result[i] = j;
      }
    }
    if (result[i] < 0) {
      std::string tsymbol;
      tsymbol += tmps.x;
      tsymbol += tmps.y;
      rtWarn("No atomic symbol " + tsymbol + " is known.  Atom " + std::to_string(i) +
             "will be represented as a virtual site (symbol VS) with atomic number 0.",
             "symbolToZNumber");
      result[i] = 0;
    }
  }
  return result;
}
  
} // namespace chemistry
} // namespace stormm
