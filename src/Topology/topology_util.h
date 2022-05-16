// -*-c++-*-
#ifndef OMNI_TOPOLOGY_BOUNDS_CHECKS_H
#define OMNI_TOPOLOGY_BOUNDS_CHECKS_H

#include <string>
#include <vector>
#include "Math/sorting.h"
#include "Math/vector_ops.h"
#include "Reporting/error_format.h"
#include "atomgraph_abstracts.h"
#include "topology_limits.h"

namespace omni {
namespace topology {

using math::UniqueValueHandling;
using math::reduceUniqueValues;

/// \brief Extract a series of numbers from a bounded list, based on an index.  Return the result
///        as a std::vector of integers.  If data from multiple bins is requested, offer an option
///        to sort and prune the data for unique values.
///
/// Overloaded:
///   - Accept a single bin index
///   - Accept multiple bin indices
///   - Accept a pre-allocated result vector with one or more bin indices
///
/// \param result     The vector to pack the outcome into (should be pre-allocated, but will be
///                   resized as necessary)
/// \param va         The original list of entries with demarcations given in va_bounds
/// \param va_bounds  Bounds array for va
/// \param index      Index of the bin of interest (checked against the size of va_bounds)
/// \param indices    Vector of indices for the bins of interest (each will be checked against the
///                   size of va_bounds)
/// \param filter     Indication that the result should be reduced to unique values only (applies
///                   to cases with multiple collection bin indices)
/// \{
template <typename T> std::vector<T> extractBoundedListEntries(const std::vector<T> &va,
                                                               const std::vector<int> &va_bounds,
                                                               int index);

template <typename T>
void extractBoundedListEntries(std::vector<T> *result, const std::vector<T> &va,
                               const std::vector<int> &va_bounds, int index);

template <typename T> std::vector<T>
extractBoundedListEntries(const std::vector<T> &va, const std::vector<int> &va_bounds,
                          const std::vector<int> &indices,
                          UniqueValueHandling filter = UniqueValueHandling::UNIQUE_VALUES_ONLY);

template <typename T> void
extractBoundedListEntries(std::vector<T> *result, const std::vector<T> &va,
                          const std::vector<int> &va_bounds, const std::vector<int> &indices,
                          UniqueValueHandling filter = UniqueValueHandling::UNIQUE_VALUES_ONLY);
/// \}
  
/// \brief Make a human-readable list of atoms based on a topology and a vector of indices,
///        providing atom numbers and name as well as residue numbers and names.  This is useful,
///        in particular, when printing error messages to the user.
///
/// \param atom_list  Topological indices for the atoms of interest
/// \param cdk        Chemical details abstract for the topology of interest
std::string writeAtomList(const std::vector<int> &atom_list, const ChemicalDetailsKit &cdk);

} // namespace topology
} // namespace omni

#include "topology_util.tpp"

#endif
