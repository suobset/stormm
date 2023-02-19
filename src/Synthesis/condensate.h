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
#include "Math/series_ops.h"
#include "Trajectory/coordinate_series.h"
#include "phasespace_synthesis.h"
#include "synthesis_enumerators.h"

namespace stormm {
namespace synthesis {

using card::GpuDetails;
using card::Hybrid;
using card::HybridTargetLevel;
using constants::CartesianDimension;
using constants::PrecisionModel;
using math::incrementingSeries;
using math::roundUp;
using trajectory::CoordinateSeries;
using trajectory::CoordinateSeriesReader;
  
/// \brief Writeable abstract for the Condensate class, wherein coordinates (only) can be modified
///        as a consequence of certain analyses.
struct CondensateWriter {

  /// \brief The constructor takes a straight list of all relevant constants and pointers.
  CondensateWriter(PrecisionModel mode_in, CondensateSource basis_in, int system_count_in,
                   int natr_insr_src_in, int natr_insr_top_in, int natr_insr_lbl_in,
                   int nata_insr_src_in, int nata_insr_top_in, int nata_insr_lbl_in,
                   const size_t* atom_starts_in, const int* atom_counts_in, float* xcrd_sp,
                   float* ycrd_sp, float* zcrd_sp, double* xcrd, double* ycrd, double* zcrd,
                   double* umat_in, double* invu_in, const int4* atr_insr_src_in,
                   const int4* atr_insr_top_in, const int4* atr_insr_lbl_in,
                   const int* atr_group_src_in, const int* atr_group_top_in,
                   const int* atr_group_lbl_in, const int4* ata_insr_src_in,
                   const int4* ata_insr_top_in, const int4* ata_insr_lbl_in,
                   const int* ata_group_src_in, const int* ata_group_top_in,
                   const int* ata_group_lbl_in);

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
  const int natr_insr_src;       ///< Number of block instructions spanning the coordinate
                                 ///<   synthesis in all-to-one analyses for groups of systems
                                 ///<   originating in the same systems cache source
  const int natr_insr_top;       ///< Number of block instructions for all-to-one analyses on all
                                 ///<   groups of systems making use of the same topology
  const int natr_insr_lbl;       ///< Number of block instructions for all-to-one analyses on all
                                 ///<   groups of systems making use of a specific cache label  
  const int nata_insr_src;       ///< Number of block instructions spanning the coordinate
                                 ///<   synthesis in all-to-all analyses for groups of systems
                                 ///<   originating in the same systems cache source
  const int nata_insr_top;       ///< Number of block instructions for all-to-all analyses on
                                 ///<   groups of systems making use of the same topology
  const int nata_insr_lbl;       ///< Number of block instructions for all-to-one analyses on all
                                 ///<   groups of systems making use of a specific cache label  
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
  const int4* atr_insr_src;      ///< Block instructions for all-to-one (All To Reference) analyses
                                 ///<   concerning groups originating in each systems cache source
  const int4* atr_insr_top;      ///< Block instructions for all-to-one (All To Reference) analyses
                                 ///<   concerning groups of systems using the same topology
  const int4* atr_insr_lbl;      ///< Block instructions for all-to-one (All To Reference) analyses
                                 ///<   concerning groups of systems using the same cache label
  const int* atr_group_src;      ///< Cache system group indices to put the above instructions in
                                 ///<   context
  const int* atr_group_top;      ///< Topology group indices to put the above instructions in
                                 ///<   context
  const int* atr_group_lbl;      ///< Cache label group indices to put the above instructions in
                                 ///<   context
  const int4* ata_insr_src;      ///< Block instructions for _A_ll _T_o _A_ll analyses concerning
                                 ///<   groups originating in each systems cache source
  const int4* ata_insr_top;      ///< Block instructions for _A_ll _T_o _A_ll analyses concerning
                                 ///<   groups of systems using the same topology
  const int4* ata_insr_lbl;      ///< Block instructions for _A_ll _T_o _A_ll analyses concerning
                                 ///<   groups groups of systems using the same cache label
  const int* ata_group_src;      ///< Cache system group indices to put the above instructions in
                                 ///<   context
  const int* ata_group_top;      ///< Topology group indices to put the above instructions in
                                 ///<   context
  const int* ata_group_lbl;      ///< Cache label group indices to put the above instructions in
                                 ///<   context
};

/// \brief Read-only abstract for the Condensate class.  In most cases, the read-only abstract will
///        be preferred as modifications to coordinates should happen in the original synthesis
///        object with its fixed precision arrays.
struct CondensateReader {

  /// \brief The constructor takes a straight list of all relevant constants and pointers.
  ///        One overloaded form accepts the corresponding writer to convert it to read-only form.
  /// \{
  CondensateReader(PrecisionModel mode_in, CondensateSource basis_in, int system_count_in,
                   int natr_insr_src_in, int natr_insr_top_in, int natr_insr_lbl_in,
                   int nata_insr_src_in, int nata_insr_top_in, int nata_insr_lbl_in,
                   const size_t* atom_starts_in, const int* atom_counts_in, const float* xcrd_sp,
                   const float* ycrd_sp, const float* zcrd_sp, const double* xcrd,
                   const double* ycrd, const double* zcrd, const double* umat_in,
                   const double* invu_in, const int4* atr_insr_src_in,
                   const int4* atr_insr_top_in, const int4* atr_insr_lbl_in,
                   const int* atr_group_src_in, const int* atr_group_top_in,
                   const int* atr_group_lbl_in, const int4* ata_insr_src_in,
                   const int4* ata_insr_top_in, const int4* ata_insr_lbl_in,
                   const int* ata_group_src_in, const int* ata_group_top_in,
                   const int* ata_group_lbl_in);

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
  const int natr_insr_src;       ///< Number of block instructions spanning the coordinate
                                 ///<   synthesis in all-to-one analyses for groups of systems
                                 ///<   originating in the same systems cache source
  const int natr_insr_top;       ///< Number of block instructions for all-to-one analyses on all
                                 ///<   groups of systems making use of the same topology
  const int natr_insr_lbl;       ///< Number of block instructions for all-to-one analyses on all
                                 ///<   groups of systems making use of a specific cache label  
  const int nata_insr_src;       ///< Number of block instructions spanning the coordinate
                                 ///<   synthesis in all-to-all analyses for groups of systems
                                 ///<   originating in the same systems cache source
  const int nata_insr_top;       ///< Number of block instructions for all-to-all analyses on
                                 ///<   groups of systems making use of the same topology
  const int nata_insr_lbl;       ///< Number of block instructions for all-to-one analyses on all
                                 ///<   groups of systems making use of a specific cache label  
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
  const int4* atr_insr_src;      ///< Block instructions for all-to-one (All To Reference) analyses
                                 ///<   concerning groups originating in each systems cache source
  const int4* atr_insr_top;      ///< Block instructions for all-to-one (All To Reference) analyses
                                 ///<   concerning groups of systems using the same topology
  const int4* atr_insr_lbl;      ///< Block instructions for all-to-one (All To Reference) analyses
                                 ///<   concerning groups of systems using the same cache label
  const int* atr_group_src;      ///< Cache system group indices to put the above instructions in
                                 ///<   context
  const int* atr_group_top;      ///< Topology group indices to put the above instructions in
                                 ///<   context
  const int* atr_group_lbl;      ///< Cache label group indices to put the above instructions in
                                 ///<   context
  const int4* ata_insr_src;      ///< Block instructions for _A_ll _T_o _A_ll analyses concerning
                                 ///<   groups originating in each systems cache source
  const int4* ata_insr_top;      ///< Block instructions for _A_ll _T_o _A_ll analyses concerning
                                 ///<   groups of systems using the same topology
  const int4* ata_insr_lbl;      ///< Block instructions for _A_ll _T_o _A_ll analyses concerning
                                 ///<   groups groups of systems using the same cache label
  const int* ata_group_src;      ///< Cache system group indices to put the above instructions in
                                 ///<   context
  const int* ata_group_top;      ///< Topology group indices to put the above instructions in
                                 ///<   context
  const int* ata_group_lbl;      ///< Cache label group indices to put the above instructions in
                                 ///<   context
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
  ///
  /// \param organization  The manner in which systems are grouped
  int getATRInstructionCount(SystemGrouping organization) const;

  /// \brief Get the number of all-to-all (ATA) block instructions.
  ///
  /// \param organization  The manner in which systems are grouped
  int getATAInstructionCount(SystemGrouping organization) const;
  
  /// \brief Get a const pointer to the original coordinate synthesis.
  const PhaseSpaceSynthesis* getSynthesisPointer() const;

  /// \brief Get a const pointer to the original coordinate series.  The developer must supply the
  ///        template type for the CoordinateSeries in order to re-interpret the arbitrarily stored
  ///        <int> type for the CoordinateSeries pointer.
  template <typename T>
  const CoordinateSeries<T>* getSeriesPointer() const;

  /// \brief Get the all-to-one instructions for some grouping of systems.
  ///
  /// \param organization  The manner in which systems are grouped
  std::vector<int4> getATRInstructions(const SystemGrouping organization) const;
  
  /// \brief Get the all-to-one instructions for some grouping of systems.
  ///
  /// \param organization  The manner in which systems are grouped
  std::vector<int> getATRInstructionGroups(const SystemGrouping organization) const;
  
  /// \brief Get the all-to-all instructions for some grouping of systems.
  ///
  /// \param organization  The manner in which systems are grouped
  std::vector<int4> getATAInstructions(const SystemGrouping organization) const;
  
  /// \brief Get the all-to-all instructions for some grouping of systems.
  ///
  /// \param organization  The manner in which systems are grouped
  std::vector<int> getATAInstructionGroups(const SystemGrouping organization) const;
  
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
  void rebuild(const CoordinateSeries<T> *cs_in,
               const PrecisionModel mode_in = PrecisionModel::SINGLE,
               const GpuDetails &gpu = null_gpu);
  /// \}

  /// \brief Add work units for a different partitioning of the systems in the underlying
  ///        synthesis.  This function will check that there is an underlying synthesis (not valid
  ///        for a series), and that the partitioning is not about unique topologies (which is
  ///        already covered when building the object from a series).
  ///
  /// Overloaded:
  ///   - Provide the partitioning scheme in C-style arrays with trusted lengths
  ///   - Provide the partitioning scheme by Standard Template Library Vectors
  ///
  /// \param system_list   Concatenated lists of system indices within the underlying synthesis
  ///                      making up each partition
  /// \param bounds_list   Bounds array for the systems list
  /// \param partitions    The number of partitions, if information is provided in C-style arrays.
  ///                      One less than the trusted length of the bounds_list array.
  /// \param organization  The grouping that the provided lists describe
  /// \param gpu           Details of the available GPU
  /// \{
  void setWorkUnits(const int* system_list, const int* bounds_list, int partitions,
                    SystemGrouping organization, const GpuDetails &gpu);
  
  void setWorkUnits(const std::vector<int> &system_list, const std::vector<int> &bounds_list,
                    SystemGrouping organization, const GpuDetails &gpu);
  /// \}

  /// \brief Reload the coordinates based on an updated coordinate synthesis or series, without
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
  int atr_instruction_count_src;   ///< The number of instructions needed to span the coordinate
                                   ///<   synthesis in all-to-reference (ATR) analyses for systems
                                   ///<   grouped under the same cache source.
  int atr_instruction_count_top;   ///< The number of instructions needed to span the coordinate
                                   ///<   synthesis in all-to-reference (ATR) analyses for systems
                                   ///<   grouped under the same unique topology.
  int atr_instruction_count_lbl;   ///< The number of instructions needed to span the coordinate
                                   ///<   synthesis in all-to-reference (ATR) analyses for systems
                                   ///<   grouped under the same cache label.
  int ata_instruction_count_src;   ///< The number of instructions needed to span the coordinate
                                   ///<   synthesis in all-to-all (ATA) analyses for systems
                                   ///<   grouped under the same cache source.
  int ata_instruction_count_top;   ///< The number of instructions needed to span the coordinate
                                   ///<   synthesis in all-to-all (ATA) analyses for systems
                                   ///<   grouped under the same unique topology.
  int ata_instruction_count_lbl;   ///< The number of instructions needed to span the coordinate
                                   ///<   synthesis in all-to-all (ATA) analyses for systems
                                   ///<   grouped under the same cache label.
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

  /// Block instructions for all-to-reference calculations with various methods of grouping the
  /// systems.  Extensions of each array name correspond to those seen in the SynthesisCacheMap
  /// object (src = source / -sys { ... } keyword entry, top = topology, lbl = cache label).  The
  /// reference structures corresponding to each group will be supplied by an auxiliary array of
  /// integers giving the structures' absolute indices within the synthesis.  Each instruction
  /// gives the index of one structure in the synthesis in its "x" member, the index of the
  /// reference structure in the auxiliary array in its "y" member, provides the unique topology
  /// index of all systems to be accessed (also the RMSD plan index) in the "z" member, and the
  /// number of the structure within its respective group in the "w" member.
  /// \{
  Hybrid<int4> atr_instructions_src;
  Hybrid<int4> atr_instructions_top;
  Hybrid<int4> atr_instructions_lbl;
  /// \}

  /// The instructions above provide offsets within specific system groupings and the indices of
  /// topologies to work from, but the group index is missing.  The following arrays provide group
  /// indices to put each of the instructions above in context.
  /// \{
  Hybrid<int> atr_instruction_groups_src;
  Hybrid<int> atr_instruction_groups_top;
  Hybrid<int> atr_instruction_groups_lbl;
  /// \}

  /// Block instructions for all-to-all calculations with various methods of grouping the systems.
  /// Extensions of each array name correspond to those seen in the SynthesisCacheMap object
  /// (src = source / -sys { ... } keyword entry, top = topology, lbl = cache label).  Each
  /// instruction contains the group indices of two structures in the "x" and "y" members,
  /// which shall be compared according to the topology of index given in the "z" member.  The
  /// low and high 16 bits of the "w" member provide the quantity of additional systems from the
  /// group list that are to participate in the same comparison.  Each instruction can therefore
  /// specify up to 256 structure comparisons.
  /// \{
  Hybrid<int4> ata_instructions_src;
  Hybrid<int4> ata_instructions_top;
  Hybrid<int4> ata_instructions_lbl;
  /// \}

  /// The instructions above provide offsets within specific system groupings and the indices of
  /// topologies to work from, but the group index is missing.  The following arrays provide group
  /// indices to put each of the instructions above in context.
  /// \{
  Hybrid<int> ata_instruction_groups_src;
  Hybrid<int> ata_instruction_groups_top;
  Hybrid<int> ata_instruction_groups_lbl;
  /// \}

  PhaseSpaceSynthesis *pps_ptr;   ///< Pointer to the PhaseSpaceSynthesis object upon which this
                                  ///<   Condensate and its block instructions are based.  Each
                                  ///<   Condensate will be based on either a PhaseSpaceSynthesis
                                  ///<   or a CoordinateSeries (see cs_ptr, below), the options
                                  ///<   being exclusive.
  CoordinateSeries<int> *cs_ptr;  ///< Pointer to the CoordinateSeries object upon which this
                                  ///<   Condensate and its block instructions are based.  This
                                  ///<   is an arbitrary type and re-interpreted so that a pointer
                                  ///<   to the original source object may be held without placing
                                  ///<   a template requirement on all Condensate objects in any
                                  ///<   situation.
  Hybrid<float> float_data;       ///< ARRAY-kind Hybrid targeted by the _sp coordinates Hybrids
  Hybrid<double> double_data;     ///< ARRAY-kind Hybrid targeted by the coordinates Hybrids

  /// \brief Repair pointers to the coordinates held by this object in the event that it has been
  ///        copied.
  void repairPointers();

  /// \brief Generate a set of work units for a particular grouping of the systems in the
  ///        synthesis (or, in a special case, the series).
  ///
  /// \param system_list
  /// \param topology_index_list
  /// \param bounds_list
  /// \param partitions
  /// \param atr_instructions
  /// \param atr_groups
  /// \param ata_instructions
  /// \param ata_groups
  /// \param atr_instruction_count
  /// \param ata_instruction_count
  /// \param gpu
  void generateWorkUnits(const int* system_list, const int* topology_index_list,
                         const int* bounds_list, int partitions, Hybrid<int4> *atr_instructions,
                         Hybrid<int> *atr_instruction_groups, Hybrid<int4> *ata_instructions,
                         Hybrid<int> *ata_instruction_groups, int *atr_instruction_count,
                         int *ata_instruction_count, const GpuDetails &gpu);
  
  /// \brief List out the work units for analysis on a CoordinateSeries.  This is templated to work
  ///        with any CoordinateSeries, and a special case of the process used to list work units
  ///        for a PhaseSpaceSynthesis--all frames of the CoordinateSeries are based on a single
  ///        topology.
  ///
  /// \param cs_in  The input coordinate series (passed down all the way from the constructor, to
  ///               avoid the need to recover its true templated type)
  /// \param gpu    Details of the available GPU
  template <typename T>
  void computeWorkUnits(const CoordinateSeries<T> *cs_in, const GpuDetails &gpu = null_gpu);

  /// \brief Compute work units for topology-based system groupings.  This can be done with only
  ///        the PhaseSpaceSynthesis, and will be the only work units if there is no cache
  ///        provided.
  ///
  /// \param gpu    Details of the available GPU
  void computeTopologyWorkUnits(const GpuDetails &gpu);
};

} // namespace synthesis
} // namespace stormm

#include "condensate.tpp"

#endif
