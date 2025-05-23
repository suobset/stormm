// -*-c++-*-
#ifndef STORMM_PDB_H
#define STORMM_PDB_H

#include <string>
#include <vector>
#include "copyright.h"
#include "Accelerator/hybrid.h"
#include "Accelerator/gpu_enumerators.h"
#include "DataTypes/stormm_vector_types.h"
#include "FileManagement/file_enumerators.h"
#include "Parsing/textfile.h"
#include "Synthesis/condensate.h"
#include "Synthesis/phasespace_synthesis.h"
#include "Topology/atomgraph.h"
#include "Topology/atomgraph_enumerators.h"
#include "Trajectory/coordinateframe.h"
#include "Trajectory/coordinate_series.h"
#include "Trajectory/phasespace.h"
#include "molecule_format_enumerators.h"

namespace stormm {
namespace structure {

using card::default_hpc_format;
using card::HybridFormat;
using diskutil::PrintSituation;
using parse::TextFile;
using synthesis::Condensate;
using synthesis::PhaseSpaceSynthesis;
using topology::AtomGraph;
using topology::UnitCellType;
using trajectory::CoordinateFrame;
using trajectory::CoordinateSeries;
using trajectory::PhaseSpace;

class Pdb {
public:

  /// \brief The constructor can accept a file to read, arrays of data to compile, or a simple
  ///        number of atoms for which to allocate space.
  ///
  /// \param atom_count_in    The anticipated number of atoms in each model entry
  /// \param model_count_in   The anticipated number of models
  /// \param file_name_in     Name of the input PDB file
  /// \param model_index      Indicate a particular model of the PDB to use in constructing the
  ///                         object
  /// \param alternate_prefs  List of symbols that will qualify for including an atom.  White space
  ///                         is an option.
  /// \param xcoordinates_in  Cartesian X coordinates of all particles, perhaps spanning multiple
  ///                         model entries.  If multiple model entries are involved, the inputs
  ///                         for each separate model must be padded by the warp size.
  /// \param ycoordinates_in  Cartesian Y coordinates of all particles
  /// \param zcoordinates_in  Cartesian Z coordinates of all particles
  /// \param atom_names_in    Names of atoms in each entry, common to all entries
  /// \param res_names_in     Names of residues in each entry, listed atom by atom and common to
  ///                         all entries
  /// \param chain_names_in   Names of each chain in the system, common to all entries
  /// \param occupancies_in   The occupancies of each atom, unique to each entry
  /// \param b_factors_in     Temperature factors indicating the spread of positions for each atom
  ///                         about some observed mean.  Unique elements must be provided for each
  ///                         entry.
  /// \param elements_in      IUPAC element symbols for each atom
  /// \param ag               Topology providing atom-specific information for each model entry
  /// \param cf               Coordinates for one system
  /// \param ps               Coordinates for one system
  /// \param cs               Coordinates for multiple copies of one system
  /// \param frames           A list of frames associated with a CoordinateSeries (input cs)
  /// \param frame_number     A specific frame associated with a CoordinateSeries
  /// \{
  Pdb(int atom_count_in = 0, int model_count_in = 1);

  Pdb(const std::string &file_name_in, int model_index = -1,
      const std::vector<char> &alternate_prefs = { ' ', 'A' });

  Pdb(const std::string &file_name_in, const std::vector<int> &model_indices,
      const std::vector<char> &alternate_prefs = { ' ', 'A' });

  Pdb(int atom_count_in, int model_count_in, const std::vector<double> &x_coordinates_in,
      const std::vector<double> &y_coordinates_in, const std::vector<double> &z_coordinates_in,
      const std::vector<char4> &atom_names_in, const std::vector<char4> &residue_names_in,
      const std::vector<char> &chain_names_in, const std::vector<double> &occupancies_in,
      const std::vector<double> &b_factors_in, const std::vector<char2> &elements_in = {},
      const std::vector<int> &atom_numbers_in = {},
      const std::vector<int> &residue_numbers_in = {},
      const std::vector<int> &anisotropy_matrices_in = {},
      const std::vector<PdbAtomKind> &atom_classes_in = {},
      const std::vector<std::string> &header_in = {},
      const std::vector<std::string> &remark_in = {},
      const std::vector<std::string> &extras_in = {});

  Pdb(const AtomGraph *ag, const CoordinateFrame *cf);

  Pdb(const AtomGraph &ag, const CoordinateFrame &cf);

  Pdb(const AtomGraph *ag, const PhaseSpace *ps);  

  Pdb(const AtomGraph &ag, const PhaseSpace &ps);

  template <typename T>
  Pdb(const AtomGraph *ag, const CoordinateSeries<T> *cf, const std::vector<int> &frames = {});

  template <typename T>
  Pdb(const AtomGraph &ag, const CoordinateSeries<T> &cf, const std::vector<int> &frames = {});

  template <typename T> Pdb(const AtomGraph *ag, const CoordinateSeries<T> *cs, int frame_number);

  template <typename T> Pdb(const AtomGraph &ag, const CoordinateSeries<T> &cs, int frame_number);
  /// \}

  /// \brief Get the number of atoms in the structure (the number of atoms per model).
  int getAtomCount() const;

  /// \brief Get the number of models in the object.
  int getModelCount() const;

  /// \brief Get the object's current directive on whether to print formal charges to any output
  ///        files.
  bool getFormalChargePrinting() const;
  
  /// \brief Produce a CoordinateFrame object from a specific model index.
  ///
  /// \param model_index  The model of interest (this will be checked against the actual number of
  ///                     models available)
  /// \param layout       The memory (CPU host and GPU device) layout with which to create the
  ///                     resulting object
  CoordinateFrame exportCoordinateFrame(int model_index,
                                        HybridFormat layout = default_hpc_format) const;

  /// \brief Produce a PhaseSpace object (with zero velocity and force components) based on the
  ///        particle positions in a specific model index.  Descriptions of input parameters follow
  ///        from getCoordinateFrame(), above.
  PhaseSpace exportPhaseSpace(int model_index, HybridFormat layout = default_hpc_format) const;

  /// \brief Produce a CoordinateSeries object based on the particle positions in a sequence of
  ///        models.
  ///
  /// \param model_sequence  The sequence of models to impart to frames of the coordinate series.
  ///                        An empty vector indicates that all models should be imparted, in
  ///                        order.
  template <typename T>
  CoordinateSeries<T> getCoordinateSeries(const std::vector<int> model_sequence = {}) const;

  /// \brief Write the contents of the object to a PDB file.
  ///
  /// \param file_name    Name of the file to write (likely to have a .pdb extension)
  /// \param expectation  Directive as to what to do if the file already exists
  void writeToFile(const std::string &file_name, PrintSituation expectation,
                   const std::vector<int> &model_indices = {}) const;
  
  /// \brief Transcribe information from a topology into the object.  This includes atom names,
  ///        symbols, and formal charges.  Atom classes will be assigned based on recognized amino
  ///        acid names.  Atom symbols for virtual sites will be listed as "VS", and all letters
  ///        in atom symbols will be capitalized.
  ///
  /// Overloaded:
  ///   - Provide the topology by const pointer
  ///   - Provide the topology by const reference
  ///
  /// \param ag  The topology containing information to transcribe
  /// \{
  void loadTopologicalData(const AtomGraph *ag);
  void loadTopologicalData(const AtomGraph &ag);
  /// \}

  /// \brief Load coordinates into the object from various coordinate objects.
  ///
  /// Overloaded:
  ///   - Provide C-style arrays for each of the three Cartesian dimensions
  ///   - Provide one of STORMM's typical coordinate objects, with a system index if necessary
  ///   
  /// \param xcrd_in       Pointer to Cartesian X coordinates, in units of Angstroms.  The order of
  ///                      atoms in these coordinates is expected to match that in the Pdb object.
  /// \param ycrd_in       Pointer to Cartesian Y coordinates, in units of Angstroms
  /// \param zcrd_in       Pointer to Cartesian Z coordinates, in units of Angstroms
  /// \param xcrd_ovrf     Overflow bits for split fixed-precision Cartesian X coordinates
  /// \param ycrd_ovrf     Overflow bits for split fixed-precision Cartesian Y coordinates
  /// \param zcrd_ovrf     Overflow bits for split fixed-precision Cartesian Z coordinates
  /// \param cf            A set of coordinates to draw from
  /// \param ps            A set of coordinates to draw from.  The POSITIONS aspect of this input
  ///                      will be taken, at whatever point in the coordinate cycle the object is
  ///                      currently at.
  /// \param cs            A coordinate series to draw from (one frame, a collection of frames, or
  ///                      all frames in the series can be taken into the object).  If all frames
  ///                      of the series are to be loaded and the series has more frames than the
  ///                      object has models, the loading will stop at the limit of the object.
  /// \param poly_ps       A coordinate synthesis to draw from (a specific system within the
  ///                      synthesis must be specified)
  /// \param cdns          An abridged coordinate synthesis to draw from (a specific system must be
  ///                      specified)
  /// \param system_index  Index of the system to take from a coordinate synthesis
  /// \param frame_index   Index of the frame to take from a coordinate series
  /// \{
  template <typename T>
  void loadCoordinates(const T* xcrd_in, const T* ycrd_in, const T* zcrd_in, int model_index = 0,
                       double gpos_scale = 1.0);

  void loadCoordinates(const llint* xcrd_in, const llint* ycrd_in, const llint* zcrd_in,
                       const int* xcrd_ovrf, const int* ycrd_ovrf, const int* zcrd_ovrf,
                       double gpos_scale = 1.0, int model_index = 0);
  
  void loadCoordinates(const CoordinateFrame *cf, int model_index = 0);

  void loadCoordinates(const CoordinateFrame &cf, int model_index = 0);

  void loadCoordinates(const PhaseSpace *ps, int model_index = 0);

  void loadCoordinates(const PhaseSpace &ps, int model_index = 0);

  template <typename T>
  void loadCoordinates(const CoordinateSeries<T> *cs, const std::vector<int2> &frame_mapping = {});

  template <typename T>
  void loadCoordinates(const CoordinateSeries<T> &cs, const std::vector<int2> &frame_mapping = {});

  void loadCoordinates(const PhaseSpaceSynthesis *poly_ps, int system_index, int model_index = 0);

  void loadCoordinates(const PhaseSpaceSynthesis &poly_ps, int system_index, int model_index = 0);

  void loadCoordinates(const Condensate *cdns, int system_index, int model_index = 0);

  void loadCoordinates(const Condensate &cdns, int system_index, int model_index = 0);
  /// \}

  /// \brief Set whether to print formal charges in any output files.
  ///
  /// \param print_formal_charges_in  The new setting to apply
  void setFormalChargePrinting(bool print_formal_charges_in = true);
  
private:
  std::string file_name;                  ///< Name of the file from which this PDB structure was
                                          ///<   obtained
  int atom_count;                         ///< The number of particles in the structure
  int model_count;                        ///< A Pdb object shares some notions in common with a
                                          ///<   CoordinateSeries object.  It can have an arbitrary
                                          ///<   number of copies of the same system in different
                                          ///<   configurations.
  int padded_atom_count;                  ///< While the Pdb object does not have HPC components
                                          ///<   (Hybrid arrays) it does mimic their layout in that
                                          ///<   each model's data is padded by the warp size.
  int ter_card_count;                     ///< The number of TER cards found in each model entry.
                                          ///<   The TER cards are based on the layout of the first
                                          ///<   model and assumed to be consistent among all
                                          ///<   models.
  UnitCellType unit_cell;                 ///< The type of unit cell.  This applies to all models.
  std::vector<PdbAtomKind> atom_classes;  ///< Indicate the PDB atom classifications
  std::vector<int> atom_numbers;          ///< Store the PDB atom numbers, which may depart from
                                          ///<   their indices in the respective arrays of
                                          ///<   descriptors
  std::vector<char4> atom_names;          ///< Four-character names of atoms, with preceding spaces
                                          ///<   removed.  For printing, white space will be added
                                          ///<   back according to PDB rules.
  std::vector<char4> residue_names;       ///< Four-character names of residues.  The first three
                                          ///<   letters are valid, while the fourth is expected to
                                          ///<   remain blank.
  std::vector<char> chain_names;          ///< Names of chains
  std::vector<int> residue_numbers;       ///< Residue numbers from the PDB file for each atom.
                                          ///<   These correspond to the residue_numbers array in
                                          ///<   the AtomGraph class (src/Topology/atomgraph.h).
  std::vector<double> x_coordinates;      ///< Array of all Cartesian X coordinates for all models
  std::vector<double> y_coordinates;      ///< Array of all Cartesian Y coordinates for all models
  std::vector<double> z_coordinates;      ///< Array of all Cartesian Z coordinates for all models
  std::vector<double> box_dimensions;     ///< And array of box dimensions (A, B, and C vector
                                          ///<   lengths, plus alpha, beta, and gamma box angles)
                                          ///<   for each model entry.  Zero box angles (0.0, 0.0,
                                          ///<   0.0) indicate that no unit cell is available,
                                          ///<   perhaps because the structure describes an NMR
                                          ///<   ensemble.
  std::vector<double> occupancies;        ///< Array of atomic occupancies for each atom, indexed
                                          ///<   alongside the atomic coordinates and padded by the
                                          ///<   warp size if multiple models are present
  std::vector<double> b_factors;          ///< Array of temperature factors
  std::vector<char2> atom_symbols;        ///< Array of atom element names
  std::vector<int> atom_formal_charges;   ///< Formal charges of each atom
  std::vector<int> ter_card_locations;    ///< A list of locations enumerating the atom indices
                                          ///<   (not atom numbers in the atom_numbers array, but
                                          ///<   indices of atoms from [ 0, 1, ..., atom_count )).
                                          ///<   Like other descriptors such as atom_classes and
                                          ///<   atom_numbers, the TER cards apply to all model
                                          ///<   entries.

  // The following are anisotropic temperature factors for each atom, providing the six unique
  // elements of each atom's matrix.
  std::vector<int> anisotropy_xx;  ///< First row and first column
  std::vector<int> anisotropy_xy;  ///< First row and second column
  std::vector<int> anisotropy_xz;  ///< First row and third column
  std::vector<int> anisotropy_yy;  ///< Second row and second column
  std::vector<int> anisotropy_yz;  ///< Second row and third column
  std::vector<int> anisotropy_zz;  ///< Third row and third column

  // The following are output settings which will be called upon when printing the Pdb object's
  // contents to a file.
  bool print_formal_charges;
  
  /// \brief Allocate space for all atom data.  Calls to this function depend on atom_count and
  ///        padded_atom_count being set.
  void allocate();
  
  /// \brief Compute the unit cell transformation matrices for a given model.  This is used when
  ///        exporting coordinate sets.
  ///
  /// \param model_index  The index of the model of interest
  /// \param umat         Pointer to the box space transformation matrix.  It will be allocated,
  ///                     filled, and returned.
  /// \param invu         Pointer to the inverse transformation matrix.  It will be allocated,
  ///                     filled, and returned.
  void computeTransforms(const int model_index, std::vector<double> *umat,
                         std::vector<double> *invu) const;

  /// \brief Count the number of ATOM and HETATM entries in a .pdb file.
  ///
  /// \param tf          The .pdb file contents, taken into RAM
  /// \param line_start  The first line of the range to read
  /// \param line_end    The final line of the range to read.  If set to less than zero the entire
  ///                    file will be considered valid material.
  int countAtomRecords(const TextFile &tf, int line_start = 0, int line_end = -1) const;

  /// \brief Read the ATOM and HETATM records for a particular model from a particular segment of
  ///        a PDB file.  Descriptions of input parameters follow from countAtomRecords(), above,
  ///        in addition to:
  ///
  /// \param atom_start  The mark at which to begin copying atoms into the object's data arrays.
  ///                    Multiple calls of the function may load different models from a single
  ///                    .pdb file into the object.
  void readAtomRecords(const TextFile &tf, int line_start, int line_end, int atom_start,
                       const std::vector<char> &alternate_prefs = { ' ', 'A' });

  /// \brief Seek out the unit cell dimensions, as stated in a CRYST1 record of the PDB file.
  ///
  /// \param tf          ASCII text of the PDB file, transcribed into memory
  /// \param line_start  The starting line of the model to search.  The search will end once a
  ///                    CRYST1 record is found.
  /// \param line_end    The final line of the model to search
  /// \param box_start   The starting index for transcribing box dimensions
  void seekUnitCellDimensions(const TextFile &tf, int line_start, int line_end, int box_start);
  
  /// \brief Validate the model index against the number actually stored in the object.
  ///
  /// \param model_index  The index of interest
  /// \param caller       Name of the calling function (for backtracing purposes)
  void validateModelIndex(int model_index, const char* caller = nullptr) const;

  /// \brief Write one of the stored models to a file.
  ///
  /// \param foutp        Open file stream for output
  /// \param model_index  Index of the stored model to write.  This will be checked for validity.
  void writeModel(std::ofstream *foutp, const int model_index) const;
};

} // namespace structure
} // namespace stormm

#include "pdb.tpp"

#endif
