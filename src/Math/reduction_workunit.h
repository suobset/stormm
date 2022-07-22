// -*-c++-*-
#ifndef OMNI_REDUCTION_WORKUNIT_H
#define OMNI_REDUCTION_WORKUNIT_H

#include <cmath>
#include <vector>
#include "Accelerator/gpu_details.h"
#include "DataTypes/common_types.h"

namespace omni {
namespace math {

using card::GpuDetails;

/// \brief The maximum number of slots into which gathering-related kernels can funnel their
///        results.
constexpr int maximum_gathering_results = 1024;

/// \brief A work unit to describe the manner in which groups of atoms in each structure of a
///        synthesis come together to contribute to a single result.  Reduction work units serve
///        one and only one system apiece.  The maximum size of these work units can be specified
///        when they are created, in a manner analogous to valence or non-bonded work units, but
///        is not limited by space in the GPU L1 cache so that, for practical optimizations,
///        gather and scatter operations can be combined into all-reduce operations with a single
///        kernel launch.
class ReductionWorkUnit {
public:

  /// \brief The constructor takes arguments for all members.
  ReductionWorkUnit(int atom_start_in, int atom_end_in, int result_index_in,
                    int dependency_start_in, int dependency_end_in, int system_index_in);

  /// \brief Take all copy and move constructors and assignment operators.
  /// \{
  ReductionWorkUnit(const ReductionWorkUnit &original) = default;
  ReductionWorkUnit(ReductionWorkUnit &&original) = default;
  ReductionWorkUnit& operator=(const ReductionWorkUnit &other) = default;
  ReductionWorkUnit& operator=(ReductionWorkUnit &&other) = default;
  /// \}
  
  /// \brief Get the atom starting index.
  int getAtomStart() const;

  /// \brief Get the upper limit of atoms in this work unit.
  int getAtomEnd() const;

  /// \brief Get the index of whatever result array where this work unit will put its result.
  int getResultIndex() const;

  /// \brief Get the start of dependencies in the result array which pertain to the same system as
  ///        this reduction work unit.  All reduction work units serving the same system will
  ///        contribute their gathering results to contiguous elements of whatever result arrays.
  int getDependencyStart() const;

  /// \brief Get the upper limit of dependencies in the result array which pertain to the same
  ///        system as this reduction work unit.
  int getDependencyEnd() const;

  /// \brief Get the system to which this reduction work unit pertains (each reduction work unit
  ///        will serve one and only one system in a synthesis)
  int getSystemIndex() const;

  /// \brief Produce an abstract containing all of the information, wrapped in a series of eight
  ///        integers.  See the enum class RdwuAbstractMap to understand which element represents
  ///        which value.
  std::vector<int> getAbstract() const;
  
private:
  int atom_start;        ///< Lower limit of atoms in the unified synthesis array which contribute
                         ///<   to this work unit's reduction operation
  int atom_end;          ///< Upper limit of atoms in the unified synthesis array which contribute
                         ///<   to this work unit's reduction operation
  int result_index;      ///< Index in the accumulation arrays to which results should be written
                         ///<   (if a reduction across two kernels is needed at all)
  int dependency_start;  ///< Staring index, in the accumulation arrays, of results from reduction
                         ///<   work units pertaining to the same system
  int dependency_end;    ///< Upper bounding index, in the accumulation arrays, of results from
                         ///<   reduction work units pertaining to the same system
  int system_index;      ///< The system to which the result pertains (this may or may not be
                         ///<   necessary, given the availability of bounds in atom_start and
                         ///<   atom_end)
};

/// \brief Collect pointers to data subject to reduction operations.  Reductions can happen on up
///        to three data sets at once, with one consistent fixed-precision scaling factor (if
///        fixed-precision applies).
template <typename T> struct ReductionSubstrate {

  /// \brief The constructor can take up to three buffers, each with possible overflows to
  ///        accommodate extended fixed-precision formats.
  ///
  /// Overloaded:
  ///   - Prepare to reduce one array of standard data
  ///   - Prepare to reduce three arrays of standard data
  ///   - Accommodate extended fixed-precision formats
  /// \{
  ReductionSubstrate(const T* x_read_in, double* x_buffer_in, T* x_write_in,
                     const int scale_bits_in = 0);

  ReductionSubstrate(const T* x_read_in, const T* y_read_in, const T* z_read_in,
                     double* x_buffer_in, double* y_buffer_in, double* z_buffer_in, T* x_write_in,
                     T* y_write_in, T* z_write_in, const int scale_bits_in = 0);

  ReductionSubstrate(const T* x_read_in, const int* x_read_ovrf_in, double* x_buffer_in,
                     T* x_write_in, int* x_write_ovrf_in, const int scale_bits_in = 0);

  ReductionSubstrate(const T* x_read_in, const int* x_read_ovrf_in, const T* y_read_in,
                     const int* y_read_ovrf_in, const T* z_read_in, const int* z_read_ovrf_in,
                     double* x_buffer_in, double* y_buffer_in, double* z_buffer_in, T* x_write_in,
                     int* x_write_ovrf_in, T* y_write_in, int* y_write_ovrf_in, T* z_write_in,
                     int* z_write_ovrf_in, const int scale_bits_in = 0);
  /// \}

  /// \brief Take the typical copy and move constructors for an abstract with const elements.
  /// \{
  ReductionSubstrate(const ReductionSubstrate &original) = default;
  ReductionSubstrate(ReductionSubstrate &&original) = default;
  /// \}

  const int dim;                ///< The number of dimensions to the data involved in the
                                ///<   reduction, i.e. X/Y/Z coordinate reduction has dimension 3
  const int scale_bits;         ///< The number of bits after the decimal in fixed-precision data
  const double fp_scaling;      ///< Scaling factor for fixed-precision data
  const double inv_fp_scaling;  ///< Inverse of the scaling factor for fixed-precision data
  const T* x_read;              ///< Read-only arrays of reducible data for the 1st dimension
  const T* y_read;              ///< Read-only arrays of reducible data for the 2nd dimension
  const T* z_read;              ///< Read-only arrays of reducible data for the 3rd dimension
  const int* x_read_ovrf;       ///< Overflow arrays for extended fixed-precision format in the
                                ///<   1st dimension of reducible data
  const int* y_read_ovrf;       ///< Overflow arrays for extended fixed-precision format in the
                                ///<   2nd dimension of reducible data
  const int* z_read_ovrf;       ///< Overflow arrays for extended fixed-precision format in the
                                ///<   3rd dimension of reducible data
  double* x_buffer;             ///< Buffer for work unit sums along the 1st dimension
  double* y_buffer;             ///< Buffer for work unit sums along the 2nd dimension
  double* z_buffer;             ///< Buffer for work unit sums along the 3rd dimension
  T* x_write;                   ///< Array accepting results of scattering in the 1st dimension
  T* y_write;                   ///< Array accepting results of scattering in the 2nd dimension
  T* z_write;                   ///< Array accepting results of scattering in the 3rd dimension
  int* x_write_ovrf;            ///< Overflow arrays for scatter results in the 1st dimension,
                                ///<   when extended fixed-precision representations are in effect
  int* y_write_ovrf;            ///< Overflow arrays for scatter results in the 2nd dimension,
                                ///<   when extended fixed-precision representations are in effect
  int* z_write_ovrf;            ///< Overflow arrays for scatter results in the 3rd dimension,
                                ///<   when extended fixed-precision representations are in effect
};

/// \brief Build the reduction (and their components of gathering and scattering) work units for
///        a series of systems of stated sizes.  Only the starting indices and system sizes are
///        needed, and only those are provided (rather than the AtomGraphSynthesis itself, from
///        which these arrays are probably derived), to avoid circular dependencies.
///
/// \param atom_starts     Starting indices for atoms of each system in the collection (the length
///                        of this array will determine the number of systems--provide a vector
///                        with a single entry of zero to get reduction work units for a single
///                        system)
/// \param atom_counts     Atom counts for all systems in the collection
/// \param gpu             Details of the selected GPU to be used for calculations
/// \param launcher        Object to collect wisdom about optimal kernel launch configurations with
///                        the work unit array chosen to suit the collection of systems
/// \param tasks_per_atom  The number of values to reduce across each system.  A center of geometry
///                        computation would have three values (X, Y, and Z), whereas a total
///                        charge summation with a mask would have only one.
std::vector<ReductionWorkUnit> buildReductionWorkUnits(const std::vector<int> &atom_starts,
                                                       const std::vector<int> &atom_counts,
                                                       const GpuDetails &gpu,
                                                       int tasks_per_atom = 1);

} // namespace math
} // namespace omni

#include "reduction_workunit.tpp"

#endif

