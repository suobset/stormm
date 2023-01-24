// -*-c++-*-
#ifndef STORMM_CONDENSATE_H
#define STORMM_CONDENSATE_H

#include "copyright.h"
#include "Accelerator/hybrid.h"
#include "Accelerator/gpu_details.h"
#include "Constants/behavior.h"
#include "Constants/hpc_bounds.h"
#include "DataTypes/common_types.h"
#include "Math/rounding.h"
#include "Synthesis/phasespace_synthesis.h"
#include "Trajectory/coordinate_series.h"
#include "synthesis_enumerators.h"

namespace stormm {
namespace synthesis {

using card::GpuDetails;
using card::Hybrid;
using card::HybridTargetLevel;
using constants::CartesianDimension;
using constants::PrecisionModel;
using math::roundUp;
using trajectory::CoordinateSeries;
using trajectory::CoordinateSeriesReader;
  
/// \brief Writeable abstract for the Condensate class, wherein coordinates (only) can be modified
///        as a consequence of certain analyses.
struct CondensateWriter {

  /// \brief The constructor takes a straight list of all relevant constants and pointers.
  CondensateWriter(PrecisionModel mode_in, CondensateSource basis_in, int system_count_in,
                   int natr_insr_in, int nata_insr_in, const size_t* atom_starts_in,
                   const int* atom_counts_in, float* xcrd_sp, float* ycrd_sp, float* zcrd_sp,
                   double* xcrd, double* ycrd, double* zcrd, double* umat_in, double* invu_in,
                   const int4* atr_insr, const int4* ata_insr);

  /// \brief The usual copy and move constructors aply for an abstract, with copy and move
  ///        assignment operators being implicitly deleted.
  ///
  /// \param original  The original abstract to copy or move
  /// \{
  CondensateWriter(const CondensateWriter &original) = default;
  CondensateWriter(CondensateWriter &&original) = default;
  /// \}

  const PrecisionModel mode;     ///< The compression mode
  const CondensateSource basis;  ///< The original material on which the original object is based
  const int system_count;        ///< The number of systems held by the underlying Condensate
  const int natr_insr;           ///< Number of block instructions needed to span the entire
                                 ///<   coordinate synthesis in all-to-one analyses
  const int nata_insr;           ///< Number of block instructions needed to span the entire
                                 ///<   coordinate synthesis in all-to-all analyses
  const size_t* atom_starts;     ///< Starting indices of each system's atoms in the data arrays
  const int* atom_counts;        ///< Number of atoms in each system
  float* xcrd_sp;                ///< Single-precision Cartesian X coordinates
  float* ycrd_sp;                ///< Single-precision Cartesian Y coordinates
  float* zcrd_sp;                ///< Single-precision Cartesian Z coordinates
  double* xcrd;                  ///< Double-precision Cartesian X coordinates
  double* ycrd;                  ///< Double-precision Cartesian Y coordinates
  double* zcrd;                  ///< Double-precision Cartesian Z coordinates
  double* umat;                  ///< Box transform information for all systems
  double* invu;                  ///< Inverse box transform information for all systems
  const int4* atr_insr;          ///< Block instructions for all-to-one (All To Reference) analyses
  const int4* ata_insr;          ///< Block instructions for _A_ll _T_o _A_ll analyses
};

/// \brief Read-only abstract for the Condensate class.  In most cases, the read-only abstract will
///        be preferred as modifications to coordinates should happen in the original synthesis
///        object with its fixed precision arrays.
struct CondensateReader {

  /// \brief The constructor takes a straight list of all relevant constants and pointers.
  ///        One overloaded form accepts the corresponding writer to convert it to read-only form.
  /// \{
  CondensateReader(PrecisionModel mode_in, CondensateSource basis_in, int system_count_in,
                   int natr_insr_in, int nata_insr_in, const size_t* atom_starts_in,
                   const int* atom_counts_in, const float* xcrd_sp_in, const float* ycrd_sp_in,
                   const float* zcrd_sp_in, const double* xcrd_in, const double* ycrd_in,
                   const double* zcrd_in, const double* umat_in, const double* invu_in,
                   const int4* atr_insr_in, const int4* ata_insr_in);

  CondensateReader(const CondensateWriter &cdw);
  /// \}

  /// \brief The usual copy and move constructors aply for an abstract, with copy and move
  ///        assignment operators being implicitly deleted.
  ///
  /// \param original  The original abstract to copy or move
  /// \{
  CondensateReader(const CondensateReader &original) = default;
  CondensateReader(CondensateReader &&original) = default;
  /// \}

  const PrecisionModel mode;     ///< The compression mode
  const CondensateSource basis;  ///< The original material on which the original object is based
  const int system_count;        ///< The number of systems held by the underlying Condensate
  const int natr_insr;           ///< Number of block instructions needed to span the entire
                                 ///<   coordinate synthesis in all-to-one analyses
  const int nata_insr;           ///< Number of block instructions needed to span the entire
                                 ///<   coordinate synthesis in all-to-all analyses
  const size_t* atom_starts;     ///< Starting indices of each system's atoms in the data arrays
  const int* atom_counts;        ///< Number of atoms in each system
  const float* xcrd_sp;          ///< Single-precision Cartesian X coordinates
  const float* ycrd_sp;          ///< Single-precision Cartesian Y coordinates
  const float* zcrd_sp;          ///< Single-precision Cartesian Z coordinates
  const double* xcrd;            ///< Double-precision Cartesian X coordinates
  const double* ycrd;            ///< Double-precision Cartesian Y coordinates
  const double* zcrd;            ///< Double-precision Cartesian Z coordinates
  const double* umat;            ///< Box transform information for all systems
  const double* invu;            ///< Inverse box transform information for all systems
  const int4* atr_insr;          ///< Block instructions for all-to-one (All To Reference) analyses
  const int4* ata_insr;          ///< Block instructions for _A_ll _T_o _A_ll analyses
};

/// \brief Condense the data format, and possibly offer a reduced representation of coordinates,
///        along with work units detailing which systems should be imported by each block on a
///        given GPU.
class Condensate {
public:

  /// The constructor requires only a complete coordinate synthesis to build the entire object.
  /// The object can be rebuilt post-creation by providing a new PhaseSpaceSynthesis, or
  /// essentially erased by providing a nullptr o conserve memory.
  ///
  /// \param poly_ps_in  The coordinate synthesis
  /// \param mode_in     Compression mode to use
  /// \param gpu         Details of the available GPU
  /// \{
  Condensate(const PhaseSpaceSynthesis *poly_ps_in = nullptr,
             PrecisionModel mode_in = PrecisionModel::SINGLE, const GpuDetails &gpu = null_gpu);

  Condensate(const PhaseSpaceSynthesis &poly_ps_in,
             PrecisionModel mode_in = PrecisionModel::SINGLE, const GpuDetails &gpu = null_gpu);

  template <typename T>
  Condensate(const CoordinateSeries<T> *cs_in,
             PrecisionModel mode_in = PrecisionModel::SINGLE, const GpuDetails &gpu = null_gpu);

  template <typename T>
  Condensate(const CoordinateSeries<T> &cs_in,
             PrecisionModel mode_in = PrecisionModel::SINGLE, const GpuDetails &gpu = null_gpu);
  /// \}

  /// \brief Copy and move constructors as well as copy and move assignment operator overloads
  ///
  /// \param original  The original object to copy or move
  /// \param other     Another object to fill the assignment statement's right-hand side
  /// \{
  Condensate(const Condensate &original);
  Condensate(Condensate &&original);
  Condensate& operator=(const Condensate &original);
  Condensate& operator=(Condensate &&original);
  /// \}
  
  /// \bried Get the compression mode.
  PrecisionModel getMode() const;

  /// \brief Get the basis for the coordinates in the object: synthesis or series
  CondensateSource getBasis() const;

  /// \brief Get the number of systems found in the condensate.
  int getSystemCount() const;
  
  /// \brief Get an indication of whether the Condensate keeps its own copy of the coordinates.
  bool ownsCoordinates() const;
  
  /// \brief Get the number of all-to-one, all-to-reference (ATR) block instructions.
  int getATRInstructionCount() const;

  /// \brief Get the number of all-to-all (ATA) block instructions.
  int getATAInstructionCount() const;
  
  /// \brief Get a const pointer to the original coordinate synthesis.
  const PhaseSpaceSynthesis* getSynthesisPointer() const;

  /// \brief Get a const pointer to the original coordinate series.  The developer must supply the
  ///        template type for the CoordinateSeries in order to re-interpret the arbitrarily stored
  ///        <int> type for the CoordinateSeries pointer.
  template <typename T>
  const CoordinateSeries<T>* getSeriesPointer() const;
  
  /// \brief Get the appropriate abstract based on the const-ness of the object.
  ///
  /// Overloaded:
  ///   - Get a reader for a const Condensate
  ///   - Get a writer for a non-const Condensate (the writer can be quickly converted to a reader)
  ///
  /// \param tier  Obtain pointers at the level of the CPU host or the GPU device
  /// \{
  const CondensateReader data(HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  CondensateWriter data(HybridTargetLevel tier = HybridTargetLevel::HOST);
  /// \}  
  
#ifdef STORMM_USE_HPC
  /// \brief Upload all data from the host to the GPU device.
  void upload();

  /// \brief Download all data from the host to the GPU device.
  void download();

  /// \brief Perform a data transfer of coordinates from the coordinate synthesis upon which this
  ///        object is based, at the level of the GPU device.  This function launches a kernel
  ///        and is called by the update() member function below.
  ///
  /// Overloaded:
  ///   - Pull coordinates from a PhaseSpaceSynthesis object
  ///   - Pull coordinates from a templated CoordinateSeries
  ///
  /// \param cs_basis  A reinterpreted cast of the internally stored CoordinateSeries pointer
  /// \param gpu       Details of the available GPU
  /// \{
  void launchCondensateUpdate(const GpuDetails &gpu);
  /// \}
#endif
  
  /// \brief Rebuild the object using a different PhaseSpaceSynthesis or CoordinateSeries, and
  ///        possibly a new compression mode.
  ///
  /// \param poly_ps_in  The new coordinate synthesis
  /// \param cs_in       The new coordinate series
  /// \param mode_in     Compression mode to use (data for any other mode in use will be cleared)
  /// \param gpu         Details of the available GPU
  /// \{
  void rebuild(const PhaseSpaceSynthesis *poly_ps_in,
               PrecisionModel mode_in = PrecisionModel::SINGLE, const GpuDetails &gpu = null_gpu);

  void rebuild(const PhaseSpaceSynthesis &poly_ps_in,
               PrecisionModel mode_in = PrecisionModel::SINGLE, const GpuDetails &gpu = null_gpu);

  template <typename T>
  void rebuild(const CoordinateSeries<T> *cs_in, const PrecisionModel mode_in,
               const GpuDetails &gpu);
  /// \}

  /// \brief Reload the coordinates based on an updated coordinate synthesis, but without
  ///        re-computing the block instructions or changing the compression level.
  ///
  /// Overloaded::
  ///   - Rely on internally stored pointers to the original coordinate objects
  ///   - Supply a pointer or reference to the CoordinateSeries object, if applicable
  ///
  /// \param cs_basis  The original coordinate series (will be checked for validity)
  /// \param tier      Obtain pointers at the level of the CPU host or the GPU device
  /// \param gpu       Details of the available GPU
  /// \{
  void update(HybridTargetLevel tier = HybridTargetLevel::HOST, const GpuDetails &gpu = null_gpu);
  
  template <typename T>
  void update(const CoordinateSeries<T> *cs_basis,
              HybridTargetLevel tier = HybridTargetLevel::HOST, const GpuDetails &gpu = null_gpu);
  
  template <typename T>
  void update(const CoordinateSeries<T> &cs_basis,
              HybridTargetLevel tier = HybridTargetLevel::HOST, const GpuDetails &gpu = null_gpu);
  /// \}

private:
  PrecisionModel mode;             ///< Mode in which the data from the PhaseSpaceSynthesis object
                                   ///<   is compressed
  CondensateSource basis;          ///< Indicate whether the condensate's coordinates are based on
                                   ///<   PhaseSpaceSynthesis or a CoordinateSeries.
  int system_count;                ///< The number of systems held by the Condensate
  bool holds_own_data;             ///< An indication of whether this object has allocated arrays
                                   ///<   to hold a separate copy of the coordinates it represents.
                                   ///<   A Condensate produced (or rebuilt) based on a
                                   ///<   PhaseSpaceSynthesis object will hold its own data, but a
                                   ///<   Condensate built upon a CoordinateSeries object may or
                                   ///<   may not (if not, its coordinate arrays will be set to
                                   ///<   the appropriate coordinate arrays of the series).
  size_t csptr_data_type;          ///< The original data type of the CoordinateSeries pointer,
                                   ///<   stored here for internal reference.  The pointer itself
                                   ///<   (see cs_ptr, below) is re-cast to an arbitrary <int>
                                   ///<   type to prevent the Condensate object itself from taking
                                   ///<   on a template type requirement.
  int atr_instruction_count;       ///< The number of instructions needed to span the coordinate
                                   ///<   synthesis in all-to-reference (ATR) analyses
  int ata_instruction_count;       ///< The number of instructions needed to span the coordinate
                                   ///<   synthesis in all-to-all (ATA) analyses
  Hybrid<size_t> atom_starts;      ///< Starting points of each system in the coordinate arrays.
                                   ///<   This is a size_t array rather than an integer array due
                                   ///<   the need to accommodate coordinate series, which could
                                   ///<   exceed the practical 2-billion atom limit for syntheses.
  Hybrid<int> atom_counts;         ///< The number of atoms in each system.  While this replicates
                                   ///<   information from the synthesis or series that the object
                                   ///<   is based upon, the amount of information is trivial, and
                                   ///<   it frees the condensate of template considerations when
                                   ///<   trying to determine a map of the contents.
  Hybrid<float> x_coordinates_sp;  ///< Cartesian X coordinates of all particles (single precision)
  Hybrid<float> y_coordinates_sp;  ///< Cartesian Y coordinates of all particles (single precision)
  Hybrid<float> z_coordinates_sp;  ///< Cartesian Z coordinates of all particles (single precision)
  Hybrid<double> x_coordinates;    ///< Cartesian X coordinates of all particles
  Hybrid<double> y_coordinates;    ///< Cartesian Y coordinates of all particles
  Hybrid<double> z_coordinates;    ///< Cartesian Z coordinates of all particles
  Hybrid<double> umat;             ///< Box space transformation matrices
  Hybrid<double> invu;             ///< Inverse box space (back to real space) transformation
                                   ///<   matrices
  Hybrid<int4> atr_instructions;   ///< Block instructions for all-to-reference calculations.
                                   ///<   The reference structures corresponding to each topology
                                   ///<   will be supplied by an auxiliary array of integers giving
                                   ///<   the structures' absolute indices within the synthesis.
                                   ///<   Each instruction gives the starting number of the
                                   ///<   system within the list of similar replicas making use of
                                   ///<   the same topology in its "x" member, leaves the "y"
                                   ///<   member blank, provides the unique topology index of all
                                   ///<   systems to be accessed (also the RMSD plan index) in the
                                   ///<   "z" member, and gives the total number of replicas to be
                                   ///<   processed within the same block instruction in the first
                                   ///<   16 bits of the "w" member.
  Hybrid<int4> ata_instructions;   ///< Block instructions for all-to-all calculations, comparing
                                   ///<   every replica making use of a particular topology to
                                   ///<   every other.  Each instruction gives the starting numbers
                                   ///<   of system indices within the list of similar replicas in
                                   ///<   its "x" and "y" members, provides the unique topology
                                   ///<   (and RMSD plan index) in the "z" member, and indicates
                                   ///<   the number of replicas beyond the starting indices in the
                                   ///<   "x" and "y" members in the "w" member's low and high 16
                                   ///<   bits, respectively.
  PhaseSpaceSynthesis *pps_ptr;    ///< Pointer to the PhaseSpaceSynthesis object upon which this
                                   ///<   Condensate and its block instructions are based.  Each
                                   ///<   Condensate will be based on either a PhaseSpaceSynthesis
                                   ///<   or a CoordinateSeries (see cs_ptr, below), the options
                                   ///<   being exclusive.
  CoordinateSeries<int> *cs_ptr;   ///< Pointer to the CoordinateSeries object upon which this
                                   ///<   Condensate and its block instructions are based.  This
                                   ///<   is an arbitrary type and re-interpreted so that a pointer
                                   ///<   to the original source object may be held without placing
                                   ///<   a template requirement on all Condensate objects in any
                                   ///<   situation.
  Hybrid<float> float_data;        ///< ARRAY-kind Hybrid targeted by the _sp coordinates Hybrids
  Hybrid<double> double_data;      ///< ARRAY-kind Hybrid targeted by the coordinates Hybrids

  /// \brief Repair pointers to the coordinates held by this object in the event that it has been
  ///        copied.
  void repairPointers();

  /// \brief List out the work units for analysis on a CoordinateSeries.  This is templated to work
  ///        with any CoordinateSeries, and a special case of the process used to list work units
  ///        for a PhaseSpaceSynthesis--all frames of the CoordinateSeries are based on a single
  ///        topology.
  ///
  /// \param cs_in  The input coordinate series (passed down all the way from the constructor, to
  ///               avoid the need to recover its true templated type)
  /// \param gpu    Details of the available GPU
  template <typename T>
  void listWorkInstructions(const CoordinateSeries<T> *cs_in, const GpuDetails &gpu = null_gpu);
};

} // namespace synthesis
} // namespace stormm

#include "condensate.tpp"

#endif
