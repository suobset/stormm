# Verbosity
VB=@

# Source and test directories
SRCDIR=src
TESTDIR=test
BENCHDIR=benchmark
APPDIR=apps
LIBDIR=lib

# STORMM C++ source files
STORMM_CPP_FILES = $(SRCDIR)/Accelerator/hybrid.cpp \
		   $(SRCDIR)/Accelerator/gpu_details.cpp \
		   $(SRCDIR)/Accelerator/kernel_manager.cpp \
		   $(SRCDIR)/Chemistry/atommask.cpp \
		   $(SRCDIR)/Chemistry/atom_equivalence.cpp \
		   $(SRCDIR)/Chemistry/chemical_features.cpp \
		   $(SRCDIR)/Chemistry/indigo.cpp \
		   $(SRCDIR)/Chemistry/znumber.cpp \
		   $(SRCDIR)/Constants/behavior.cpp \
		   $(SRCDIR)/Constants/generalized_born.cpp \
		   $(SRCDIR)/FileManagement/directory_util.cpp \
		   $(SRCDIR)/FileManagement/file_enumerators.cpp \
		   $(SRCDIR)/FileManagement/file_listing.cpp \
		   $(SRCDIR)/FileManagement/file_util.cpp \
		   $(SRCDIR)/ForceField/forcefield_element.cpp \
		   $(SRCDIR)/ForceField/forcefield_enumerators.cpp \
		   $(SRCDIR)/Math/matrix_ops.cpp \
		   $(SRCDIR)/Math/reduction.cpp \
		   $(SRCDIR)/Math/reduction_abstracts.cpp \
		   $(SRCDIR)/Math/reduction_bridge.cpp \
		   $(SRCDIR)/Math/reduction_workunit.cpp \
		   $(SRCDIR)/Math/rounding.cpp \
		   $(SRCDIR)/Math/series_ops.cpp \
		   $(SRCDIR)/Math/statistics.cpp \
		   $(SRCDIR)/Math/tickcounter.cpp \
		   $(SRCDIR)/Math/tricubic_cell.cpp \
		   $(SRCDIR)/Math/vector_ops.cpp \
		   $(SRCDIR)/MolecularMechanics/line_minimization.cpp \
		   $(SRCDIR)/MolecularMechanics/minimization.cpp \
		   $(SRCDIR)/MolecularMechanics/mm_controls.cpp \
		   $(SRCDIR)/MolecularMechanics/mm_evaluation.cpp \
		   $(SRCDIR)/MoleculeFormat/mdlmol.cpp \
		   $(SRCDIR)/MoleculeFormat/mdlmol_atomlist.cpp \
		   $(SRCDIR)/MoleculeFormat/mdlmol_bond.cpp \
		   $(SRCDIR)/MoleculeFormat/mdlmol_dataitem.cpp \
		   $(SRCDIR)/MoleculeFormat/mdlmol_property.cpp \
		   $(SRCDIR)/MoleculeFormat/mdlmol_refinement.cpp \
		   $(SRCDIR)/MoleculeFormat/mdlmol_request.cpp \
		   $(SRCDIR)/MoleculeFormat/molecule_file_io.cpp \
		   $(SRCDIR)/MoleculeFormat/molecule_format_enumerators.cpp \
		   $(SRCDIR)/MoleculeFormat/molecule_parsing.cpp \
		   $(SRCDIR)/Namelists/input.cpp \
		   $(SRCDIR)/Namelists/namelist_element.cpp \
		   $(SRCDIR)/Namelists/namelist_emulator.cpp \
		   $(SRCDIR)/Namelists/namelist_enumerators.cpp \
		   $(SRCDIR)/Namelists/nml_conformer.cpp \
		   $(SRCDIR)/Namelists/nml_dynamics.cpp \
		   $(SRCDIR)/Namelists/nml_ffmorph.cpp \
		   $(SRCDIR)/Namelists/nml_files.cpp \
		   $(SRCDIR)/Namelists/nml_minimize.cpp \
		   $(SRCDIR)/Namelists/nml_precision.cpp \
		   $(SRCDIR)/Namelists/nml_random.cpp \
		   $(SRCDIR)/Namelists/nml_report.cpp \
		   $(SRCDIR)/Namelists/nml_restraint.cpp \
		   $(SRCDIR)/Namelists/nml_solvent.cpp \
		   $(SRCDIR)/Namelists/user_settings.cpp \
		   $(SRCDIR)/Numerics/split_fixed_precision.cpp \
		   $(SRCDIR)/Parsing/ascii_numbers.cpp \
		   $(SRCDIR)/Parsing/citation.cpp \
		   $(SRCDIR)/Parsing/parse.cpp \
		   $(SRCDIR)/Parsing/polynumeric.cpp \
		   $(SRCDIR)/Parsing/tabulation.cpp \
		   $(SRCDIR)/Parsing/textfile.cpp \
		   $(SRCDIR)/Parsing/textguard.cpp \
		   $(SRCDIR)/Potential/cacheresource.cpp \
		   $(SRCDIR)/Potential/energy_enumerators.cpp \
		   $(SRCDIR)/Potential/eval_valence_workunit.cpp \
		   $(SRCDIR)/Potential/forward_exclusionmask.cpp \
		   $(SRCDIR)/Potential/nonbonded_potential.cpp \
		   $(SRCDIR)/Potential/scorecard.cpp \
		   $(SRCDIR)/Potential/static_exclusionmask.cpp \
		   $(SRCDIR)/Potential/valence_potential.cpp \
		   $(SRCDIR)/Random/random.cpp \
		   $(SRCDIR)/Reporting/code_dox.cpp \
		   $(SRCDIR)/Reporting/display.cpp \
		   $(SRCDIR)/Reporting/help_messages.cpp \
		   $(SRCDIR)/Reporting/error_format.cpp \
		   $(SRCDIR)/Reporting/reporting_enumerators.cpp \
		   $(SRCDIR)/Reporting/summary_file.cpp \
		   $(SRCDIR)/Restraints/bounded_restraint.cpp \
		   $(SRCDIR)/Restraints/restraint_apparatus.cpp \
		   $(SRCDIR)/Restraints/restraint_builder.cpp \
		   $(SRCDIR)/Restraints/restraint_enumerators.cpp \
		   $(SRCDIR)/Structure/global_manipulation.cpp \
		   $(SRCDIR)/Structure/isomerization.cpp \
		   $(SRCDIR)/Structure/local_arrangement.cpp \
		   $(SRCDIR)/Structure/mesh_parameters.cpp \
		   $(SRCDIR)/Structure/rmsd.cpp \
		   $(SRCDIR)/Structure/rmsd_plan.cpp \
		   $(SRCDIR)/Structure/structure_ops.cpp \
		   $(SRCDIR)/Structure/structure_utils.cpp \
		   $(SRCDIR)/Structure/virtual_site_handling.cpp \
		   $(SRCDIR)/Synthesis/atomgraph_synthesis.cpp \
		   $(SRCDIR)/Synthesis/implicit_solvent_workspace.cpp \
		   $(SRCDIR)/Synthesis/phasespace_synthesis.cpp \
		   $(SRCDIR)/Synthesis/nonbonded_workunit.cpp \
		   $(SRCDIR)/Synthesis/static_mask_synthesis.cpp \
		   $(SRCDIR)/Synthesis/systemcache.cpp \
		   $(SRCDIR)/Synthesis/valence_workunit.cpp \
		   $(SRCDIR)/Topology/amber_prmtop_util.cpp \
		   $(SRCDIR)/Topology/atomgraph_constructors.cpp \
		   $(SRCDIR)/Topology/atomgraph_detailers.cpp \
		   $(SRCDIR)/Topology/atomgraph_getters.cpp \
		   $(SRCDIR)/Topology/atomgraph_intake.cpp \
		   $(SRCDIR)/Topology/atomgraph_setters.cpp \
		   $(SRCDIR)/Topology/atomgraph_miscellaneous.cpp \
		   $(SRCDIR)/Topology/atomgraph_abstracts.cpp \
		   $(SRCDIR)/Topology/atomgraph_analysis.cpp \
		   $(SRCDIR)/Topology/atomgraph_enumerators.cpp \
		   $(SRCDIR)/Topology/atomgraph_refinement.cpp \
		   $(SRCDIR)/Topology/topology_limits.cpp \
		   $(SRCDIR)/Topology/topology_util.cpp \
		   $(SRCDIR)/Trajectory/amber_ascii.cpp \
		   $(SRCDIR)/Trajectory/barostat.cpp \
		   $(SRCDIR)/Trajectory/coordinateframe.cpp \
		   $(SRCDIR)/Trajectory/coordinate_intake.cpp \
		   $(SRCDIR)/Trajectory/coordinate_swap_plan.cpp \
		   $(SRCDIR)/Trajectory/phasespace.cpp \
		   $(SRCDIR)/Trajectory/thermostat.cpp \
		   $(SRCDIR)/Trajectory/trajectory_enumerators.cpp \
		   $(SRCDIR)/Trajectory/write_annotated_frame.cpp \
		   $(SRCDIR)/Trajectory/write_frame.cpp \
		   $(SRCDIR)/UnitTesting/approx.cpp \
		   $(SRCDIR)/UnitTesting/stopwatch.cpp \
		   $(SRCDIR)/UnitTesting/checklist.cpp \
		   $(SRCDIR)/UnitTesting/file_snapshot.cpp \
		   $(SRCDIR)/UnitTesting/test_environment.cpp \
		   $(SRCDIR)/UnitTesting/test_system_manager.cpp \
		   $(SRCDIR)/UnitTesting/unit_test.cpp \
		   $(SRCDIR)/UnitTesting/unit_test_enumerators.cpp \
		   $(SRCDIR)/UnitTesting/vector_report.cpp

# STORMM C++ header files
STORMM_CPP_HEADERS = $(SRCDIR)/copyright.h \
		     $(SRCDIR)/Accelerator/hybrid.h \
		     $(SRCDIR)/Accelerator/gpu_details.h \
		     $(SRCDIR)/Accelerator/kernel_manager.h \
		     $(SRCDIR)/Chemistry/atommask.h \
		     $(SRCDIR)/Chemistry/atom_equivalence.h \
		     $(SRCDIR)/Chemistry/chemical_features.h \
		     $(SRCDIR)/Chemistry/indigo.h \
		     $(SRCDIR)/Constants/behavior.h \
		     $(SRCDIR)/Constants/fixed_precision.h \
		     $(SRCDIR)/Constants/generalized_born.h \
		     $(SRCDIR)/Constants/scaling.h \
		     $(SRCDIR)/Constants/symbol_values.h \
		     $(SRCDIR)/Chemistry/chemistry_enumerators.h \
		     $(SRCDIR)/Chemistry/periodic_table.h \
		     $(SRCDIR)/Chemistry/znumber.h \
		     $(SRCDIR)/DataTypes/common_types.h \
		     $(SRCDIR)/DataTypes/mixed_types.h \
		     $(SRCDIR)/DataTypes/stormm_vector_types.h \
		     $(SRCDIR)/FileManagement/directory_util.h \
		     $(SRCDIR)/FileManagement/file_enumerators.h \
		     $(SRCDIR)/FileManagement/file_listing.h \
		     $(SRCDIR)/FileManagement/file_util.h \
		     $(SRCDIR)/ForceField/forcefield_element.h \
		     $(SRCDIR)/ForceField/forcefield_enumerators.h \
		     $(SRCDIR)/Math/matrix.h \
		     $(SRCDIR)/Math/matrix_ops.h \
		     $(SRCDIR)/Math/multiplication.h \
		     $(SRCDIR)/Math/reduction.h \
		     $(SRCDIR)/Math/reduction_abstracts.h \
		     $(SRCDIR)/Math/reduction_bridge.h \
		     $(SRCDIR)/Math/reduction_enumerators.h \
		     $(SRCDIR)/Math/reduction_workunit.h \
		     $(SRCDIR)/Math/rounding.h \
		     $(SRCDIR)/Math/series_ops.h \
		     $(SRCDIR)/Math/sorting.h \
		     $(SRCDIR)/Math/statistics.h \
		     $(SRCDIR)/Math/summation.h \
		     $(SRCDIR)/Math/tickcounter.h \
		     $(SRCDIR)/Math/tricubic_cell.h \
		     $(SRCDIR)/Math/vector_ops.h \
		     $(SRCDIR)/MolecularMechanics/line_minimization.h \
		     $(SRCDIR)/MolecularMechanics/minimization.h \
		     $(SRCDIR)/MolecularMechanics/mm_controls.h \
		     $(SRCDIR)/MolecularMechanics/mm_evaluation.h \
		     $(SRCDIR)/MoleculeFormat/mdlmol.h \
		     $(SRCDIR)/MoleculeFormat/mdlmol_atomlist.h \
		     $(SRCDIR)/MoleculeFormat/mdlmol_bond.h \
		     $(SRCDIR)/MoleculeFormat/mdlmol_dataitem.h \
		     $(SRCDIR)/MoleculeFormat/mdlmol_property.h \
		     $(SRCDIR)/MoleculeFormat/mdlmol_refinement.h \
		     $(SRCDIR)/MoleculeFormat/mdlmol_request.h \
		     $(SRCDIR)/MoleculeFormat/molecule_file_io.h \
		     $(SRCDIR)/MoleculeFormat/molecule_format_enumerators.h \
		     $(SRCDIR)/MoleculeFormat/molecule_parsing.h \
		     $(SRCDIR)/Namelists/input.h \
		     $(SRCDIR)/Namelists/namelist_element.h \
		     $(SRCDIR)/Namelists/namelist_emulator.h \
		     $(SRCDIR)/Namelists/namelist_enumerators.h \
		     $(SRCDIR)/Namelists/nml_conformer.h \
		     $(SRCDIR)/Namelists/nml_dynamics.h \
		     $(SRCDIR)/Namelists/nml_ffmorph.h \
		     $(SRCDIR)/Namelists/nml_files.h \
		     $(SRCDIR)/Namelists/nml_minimize.h \
		     $(SRCDIR)/Namelists/nml_precision.h \
		     $(SRCDIR)/Namelists/nml_random.h \
		     $(SRCDIR)/Namelists/nml_report.h \
		     $(SRCDIR)/Namelists/nml_restraint.h \
		     $(SRCDIR)/Namelists/nml_solvent.h \
		     $(SRCDIR)/Namelists/user_settings.h \
		     $(SRCDIR)/Numerics/split_fixed_precision.h \
		     $(SRCDIR)/Parsing/ascii_numbers.h \
		     $(SRCDIR)/Parsing/citation.h \
		     $(SRCDIR)/Parsing/parse.h \
		     $(SRCDIR)/Parsing/polynumeric.h \
		     $(SRCDIR)/Parsing/tabulation.h \
		     $(SRCDIR)/Parsing/textfile.h \
		     $(SRCDIR)/Parsing/textguard.h \
		     $(SRCDIR)/Potential/cacheresource.h \
		     $(SRCDIR)/Potential/energy_abstracts.h \
		     $(SRCDIR)/Potential/energy_enumerators.h \
		     $(SRCDIR)/Potential/eval_valence_workunit.h \
		     $(SRCDIR)/Potential/eval_synthesis.h \
		     $(SRCDIR)/Potential/forward_exclusionmask.h \
		     $(SRCDIR)/Potential/nonbonded_potential.h \
		     $(SRCDIR)/Potential/scorecard.h \
		     $(SRCDIR)/Potential/static_exclusionmask.h \
		     $(SRCDIR)/Potential/valence_potential.h \
		     $(SRCDIR)/Random/random.h \
		     $(SRCDIR)/Reporting/code_dox.h \
		     $(SRCDIR)/Reporting/display.h \
		     $(SRCDIR)/Reporting/help_messages.h \
		     $(SRCDIR)/Reporting/error_format.h \
		     $(SRCDIR)/Reporting/reporting_enumerators.h \
		     $(SRCDIR)/Reporting/summary_file.h \
		     $(SRCDIR)/Restraints/bounded_restraint.h \
		     $(SRCDIR)/Restraints/restraint_apparatus.h \
		     $(SRCDIR)/Restraints/restraint_builder.h \
		     $(SRCDIR)/Restraints/restraint_enumerators.h \
		     $(SRCDIR)/Restraints/restraint_util.h \
		     $(SRCDIR)/Structure/background_mesh.h \
		     $(SRCDIR)/Structure/global_manipulation.h \
		     $(SRCDIR)/Structure/isomerization.h \
		     $(SRCDIR)/Structure/local_arrangement.h \
		     $(SRCDIR)/Structure/mesh_parameters.h \
		     $(SRCDIR)/Structure/rmsd_plan.h \
		     $(SRCDIR)/Structure/structure_enumerators.h \
		     $(SRCDIR)/Structure/structure_ops.h \
		     $(SRCDIR)/Structure/structure_utils.h \
		     $(SRCDIR)/Structure/virtual_site_handling.h \
		     $(SRCDIR)/Synthesis/atomgraph_synthesis.h \
		     $(SRCDIR)/Synthesis/implicit_solvent_workspace.h \
		     $(SRCDIR)/Synthesis/phasespace_synthesis.h \
		     $(SRCDIR)/Synthesis/nonbonded_workunit.h \
		     $(SRCDIR)/Synthesis/static_mask_synthesis.h \
		     $(SRCDIR)/Synthesis/synthesis_abstracts.h \
		     $(SRCDIR)/Synthesis/synthesis_enumerators.h \
		     $(SRCDIR)/Synthesis/systemcache.h \
		     $(SRCDIR)/Synthesis/valence_workunit.h \
		     $(SRCDIR)/Topology/amber_prmtop_util.h \
		     $(SRCDIR)/Topology/atomgraph.h \
		     $(SRCDIR)/Topology/atomgraph_abstracts.h \
		     $(SRCDIR)/Topology/atomgraph_analysis.h \
		     $(SRCDIR)/Topology/atomgraph_enumerators.h \
		     $(SRCDIR)/Topology/atomgraph_intake.h \
		     $(SRCDIR)/Topology/atomgraph_refinement.h \
		     $(SRCDIR)/Topology/topology_limits.h \
		     $(SRCDIR)/Topology/topology_util.h \
		     $(SRCDIR)/Trajectory/amber_ascii.h \
		     $(SRCDIR)/Trajectory/barostat.h \
		     $(SRCDIR)/Trajectory/coordinateframe.h \
		     $(SRCDIR)/Trajectory/coordinate_intake.h \
		     $(SRCDIR)/Trajectory/coordinate_series.h \
		     $(SRCDIR)/Trajectory/coordinate_swap_plan.h \
		     $(SRCDIR)/Trajectory/phasespace.h \
		     $(SRCDIR)/Trajectory/thermostat.h \
		     $(SRCDIR)/Trajectory/trajectory_enumerators.h \
		     $(SRCDIR)/Trajectory/write_annotated_frame.h \
		     $(SRCDIR)/Trajectory/write_frame.h \
		     $(SRCDIR)/UnitTesting/approx.h \
		     $(SRCDIR)/UnitTesting/stopwatch.h \
		     $(SRCDIR)/UnitTesting/checklist.h \
		     $(SRCDIR)/UnitTesting/file_snapshot.h \
		     $(SRCDIR)/UnitTesting/test_environment.h \
		     $(SRCDIR)/UnitTesting/test_system_manager.h \
		     $(SRCDIR)/UnitTesting/unit_test.h \
		     $(SRCDIR)/UnitTesting/unit_test_enumerators.h \
		     $(SRCDIR)/UnitTesting/vector_report.h

# STORMM C++ template source files
STORMM_TPP_FILES = $(SRCDIR)/Accelerator/hybrid.tpp \
		   $(SRCDIR)/Constants/generalized_born.tpp \
		   $(SRCDIR)/DataTypes/common_types.tpp \
		   $(SRCDIR)/DataTypes/stormm_vector_types.tpp \
		   $(SRCDIR)/Math/matrix.tpp \
		   $(SRCDIR)/Math/matrix_ops.tpp \
		   $(SRCDIR)/Math/multiplication.tpp \
		   $(SRCDIR)/Math/reduction.tpp \
		   $(SRCDIR)/Math/reduction_abstracts.tpp \
		   $(SRCDIR)/Math/rounding.tpp \
		   $(SRCDIR)/Math/series_ops.tpp \
		   $(SRCDIR)/Math/summation.tpp \
		   $(SRCDIR)/Math/tickcounter.tpp \
		   $(SRCDIR)/Math/tricubic_cell.tpp \
		   $(SRCDIR)/Math/vector_ops.tpp \
		   $(SRCDIR)/MolecularMechanics/minimization.tpp \
		   $(SRCDIR)/MolecularMechanics/mm_controls.tpp \
		   $(SRCDIR)/MolecularMechanics/mm_evaluation.tpp \
		   $(SRCDIR)/MoleculeFormat/mdlmol.tpp \
		   $(SRCDIR)/Parsing/polynumeric.tpp \
		   $(SRCDIR)/Potential/cacheresource.tpp \
		   $(SRCDIR)/Potential/energy_abstracts.tpp \
		   $(SRCDIR)/Potential/eval_synthesis.tpp \
		   $(SRCDIR)/Potential/nonbonded_potential.tpp \
		   $(SRCDIR)/Potential/scorecard.tpp \
		   $(SRCDIR)/Potential/valence_potential.tpp \
		   $(SRCDIR)/Restraints/restraint_apparatus.tpp \
		   $(SRCDIR)/Restraints/restraint_util.tpp \
		   $(SRCDIR)/Structure/global_manipulation.tpp \
		   $(SRCDIR)/Structure/isomerization.tpp \
		   $(SRCDIR)/Structure/local_arrangement.tpp \
		   $(SRCDIR)/Structure/background_mesh.tpp \
		   $(SRCDIR)/Structure/mesh_parameters.tpp \
		   $(SRCDIR)/Structure/rmsd.tpp \
		   $(SRCDIR)/Structure/virtual_site_handling.tpp \
		   $(SRCDIR)/Synthesis/implicit_solvent_workspace.tpp \
		   $(SRCDIR)/Synthesis/nonbonded_workunit.tpp \
		   $(SRCDIR)/Synthesis/phasespace_synthesis.tpp \
		   $(SRCDIR)/Synthesis/synthesis_abstracts.tpp \
		   $(SRCDIR)/Topology/atomgraph.tpp \
		   $(SRCDIR)/Topology/atomgraph_abstracts.tpp \
		   $(SRCDIR)/Topology/topology_util.tpp \
		   $(SRCDIR)/Trajectory/coordinateframe.tpp \
		   $(SRCDIR)/Trajectory/coordinate_series.tpp \
		   $(SRCDIR)/Trajectory/phasespace.tpp \
		   $(SRCDIR)/Trajectory/thermostat.tpp \
		   $(SRCDIR)/Trajectory/write_annotated_frame.tpp \
		   $(SRCDIR)/UnitTesting/approx.tpp \
		   $(SRCDIR)/UnitTesting/unit_test.tpp

# STORMM C++ object files
STORMM_CPP_OBJS = $(SRCDIR)/Accelerator/hybrid.o \
		  $(SRCDIR)/Accelerator/gpu_details.o \
		  $(SRCDIR)/Accelerator/kernel_manager.o \
		  $(SRCDIR)/Chemistry/atommask.o \
		  $(SRCDIR)/Chemistry/atom_equivalence.o \
		  $(SRCDIR)/Chemistry/chemical_features.o \
		  $(SRCDIR)/Chemistry/indigo.o \
		  $(SRCDIR)/Chemistry/znumber.o \
		  $(SRCDIR)/Constants/behavior.o \
		  $(SRCDIR)/Constants/generalized_born.o \
		  $(SRCDIR)/FileManagement/directory_util.o \
		  $(SRCDIR)/FileManagement/file_enumerators.o \
		  $(SRCDIR)/FileManagement/file_listing.o \
		  $(SRCDIR)/FileManagement/file_util.o \
		  $(SRCDIR)/ForceField/forcefield_element.o \
		  $(SRCDIR)/ForceField/forcefield_enumerators.o \
		  $(SRCDIR)/Math/matrix_ops.o \
		  $(SRCDIR)/Math/reduction.o \
		  $(SRCDIR)/Math/reduction_abstracts.o \
		  $(SRCDIR)/Math/reduction_bridge.o \
		  $(SRCDIR)/Math/reduction_workunit.o \
		  $(SRCDIR)/Math/rounding.o \
		  $(SRCDIR)/Math/series_ops.o \
		  $(SRCDIR)/Math/statistics.o \
		  $(SRCDIR)/Math/tickcounter.o \
		  $(SRCDIR)/Math/tricubic_cell.o \
		  $(SRCDIR)/Math/vector_ops.o \
		  $(SRCDIR)/MolecularMechanics/line_minimization.o \
		  $(SRCDIR)/MolecularMechanics/minimization.o \
		  $(SRCDIR)/MolecularMechanics/mm_controls.o \
		  $(SRCDIR)/MolecularMechanics/mm_evaluation.o \
		  $(SRCDIR)/MoleculeFormat/mdlmol.o \
		  $(SRCDIR)/MoleculeFormat/mdlmol_atomlist.o \
		  $(SRCDIR)/MoleculeFormat/mdlmol_bond.o \
		  $(SRCDIR)/MoleculeFormat/mdlmol_dataitem.o \
		  $(SRCDIR)/MoleculeFormat/mdlmol_property.o \
		  $(SRCDIR)/MoleculeFormat/mdlmol_refinement.o \
		  $(SRCDIR)/MoleculeFormat/mdlmol_request.o \
		  $(SRCDIR)/MoleculeFormat/molecule_file_io.o \
		  $(SRCDIR)/MoleculeFormat/molecule_format_enumerators.o \
		  $(SRCDIR)/MoleculeFormat/molecule_parsing.o \
		  $(SRCDIR)/Namelists/input.o \
		  $(SRCDIR)/Namelists/namelist_element.o \
		  $(SRCDIR)/Namelists/namelist_emulator.o \
		  $(SRCDIR)/Namelists/namelist_enumerators.o \
		  $(SRCDIR)/Namelists/nml_conformer.o \
		  $(SRCDIR)/Namelists/nml_dynamics.o \
		  $(SRCDIR)/Namelists/nml_ffmorph.o \
		  $(SRCDIR)/Namelists/nml_files.o \
		  $(SRCDIR)/Namelists/nml_minimize.o \
		  $(SRCDIR)/Namelists/nml_precision.o \
		  $(SRCDIR)/Namelists/nml_random.o \
		  $(SRCDIR)/Namelists/nml_report.o \
		  $(SRCDIR)/Namelists/nml_restraint.o \
		  $(SRCDIR)/Namelists/nml_solvent.o \
		  $(SRCDIR)/Namelists/user_settings.o \
		  $(SRCDIR)/Numerics/split_fixed_precision.o \
		  $(SRCDIR)/Parsing/ascii_numbers.o \
		  $(SRCDIR)/Parsing/citation.o \
		  $(SRCDIR)/Parsing/parse.o \
		  $(SRCDIR)/Parsing/polynumeric.o \
		  $(SRCDIR)/Parsing/tabulation.o \
		  $(SRCDIR)/Parsing/textfile.o \
		  $(SRCDIR)/Parsing/textguard.o \
		  $(SRCDIR)/Potential/cacheresource.o \
		  $(SRCDIR)/Potential/energy_enumerators.o \
		  $(SRCDIR)/Potential/eval_valence_workunit.o \
		  $(SRCDIR)/Potential/forward_exclusionmask.o \
		  $(SRCDIR)/Potential/nonbonded_potential.o \
		  $(SRCDIR)/Potential/scorecard.o \
		  $(SRCDIR)/Potential/static_exclusionmask.o \
		  $(SRCDIR)/Potential/valence_potential.o \
		  $(SRCDIR)/Random/random.o \
		  $(SRCDIR)/Reporting/code_dox.o \
		  $(SRCDIR)/Reporting/display.o \
		  $(SRCDIR)/Reporting/help_messages.o \
		  $(SRCDIR)/Reporting/error_format.o \
		  $(SRCDIR)/Reporting/reporting_enumerators.o \
		  $(SRCDIR)/Reporting/summary_file.o \
		  $(SRCDIR)/Restraints/bounded_restraint.o \
		  $(SRCDIR)/Restraints/restraint_apparatus.o \
		  $(SRCDIR)/Restraints/restraint_builder.o \
		  $(SRCDIR)/Restraints/restraint_enumerators.o \
		  $(SRCDIR)/Structure/global_manipulation.o \
		  $(SRCDIR)/Structure/isomerization.o \
		  $(SRCDIR)/Structure/local_arrangement.o \
		  $(SRCDIR)/Structure/mesh_parameters.o \
		  $(SRCDIR)/Structure/rmsd.o \
		  $(SRCDIR)/Structure/rmsd_plan.o \
		  $(SRCDIR)/Structure/structure_ops.o \
		  $(SRCDIR)/Structure/structure_utils.o \
		  $(SRCDIR)/Structure/virtual_site_handling.o \
		  $(SRCDIR)/Synthesis/atomgraph_synthesis.o \
		  $(SRCDIR)/Synthesis/implicit_solvent_workspace.o \
		  $(SRCDIR)/Synthesis/phasespace_synthesis.o \
		  $(SRCDIR)/Synthesis/nonbonded_workunit.o \
		  $(SRCDIR)/Synthesis/static_mask_synthesis.o \
		  $(SRCDIR)/Synthesis/systemcache.o \
		  $(SRCDIR)/Synthesis/valence_workunit.o \
		  $(SRCDIR)/Topology/amber_prmtop_util.o \
		  $(SRCDIR)/Topology/atomgraph_abstracts.o \
		  $(SRCDIR)/Topology/atomgraph_analysis.o \
		  $(SRCDIR)/Topology/atomgraph_constructors.o \
		  $(SRCDIR)/Topology/atomgraph_detailers.o \
		  $(SRCDIR)/Topology/atomgraph_enumerators.o \
		  $(SRCDIR)/Topology/atomgraph_getters.o \
		  $(SRCDIR)/Topology/atomgraph_intake.o \
		  $(SRCDIR)/Topology/atomgraph_miscellaneous.o \
		  $(SRCDIR)/Topology/atomgraph_refinement.o \
		  $(SRCDIR)/Topology/atomgraph_setters.o \
		  $(SRCDIR)/Topology/topology_limits.o \
		  $(SRCDIR)/Topology/topology_util.o \
		  $(SRCDIR)/Trajectory/amber_ascii.o \
		  $(SRCDIR)/Trajectory/barostat.o \
		  $(SRCDIR)/Trajectory/coordinateframe.o \
		  $(SRCDIR)/Trajectory/coordinate_intake.o \
		  $(SRCDIR)/Trajectory/coordinate_swap_plan.o \
		  $(SRCDIR)/Trajectory/phasespace.o \
		  $(SRCDIR)/Trajectory/thermostat.o \
		  $(SRCDIR)/Trajectory/trajectory_enumerators.o \
		  $(SRCDIR)/Trajectory/write_annotated_frame.o \
		  $(SRCDIR)/Trajectory/write_frame.o \
		  $(SRCDIR)/UnitTesting/approx.o \
		  $(SRCDIR)/UnitTesting/stopwatch.o \
		  $(SRCDIR)/UnitTesting/checklist.o \
		  $(SRCDIR)/UnitTesting/file_snapshot.o \
		  $(SRCDIR)/UnitTesting/test_environment.o \
		  $(SRCDIR)/UnitTesting/test_system_manager.o \
		  $(SRCDIR)/UnitTesting/unit_test.o \
		  $(SRCDIR)/UnitTesting/unit_test_enumerators.o \
		  $(SRCDIR)/UnitTesting/vector_report.o

# STORMM CUDA source files
STORMM_CUDA_FILES = $(SRCDIR)/Accelerator/hpc_config.cu \
		    $(SRCDIR)/Math/hpc_reduction.cu \
		    $(SRCDIR)/MolecularMechanics/hpc_minimization.cu \
		    $(SRCDIR)/Potential/hpc_nonbonded_potential.cu \
		    $(SRCDIR)/Potential/hpc_valence_potential.cu \
		    $(SRCDIR)/Random/hpc_random.cu \
		    $(SRCDIR)/Structure/hpc_virtual_site_handling.cu \
		    $(SRCDIR)/Synthesis/hpc_implicit_solvent_workspace.cu \
		    $(SRCDIR)/Synthesis/hpc_phasespace_synthesis.cu

# STORMM CUDA header files
STORMM_CUDA_HEADERS = $(SRCDIR)/Constants/hpc_bounds.h \
		      $(SRCDIR)/Accelerator/hpc_config.h \
		      $(SRCDIR)/Accelerator/ptx_macros.h \
		      $(SRCDIR)/Math/hpc_reduction.h \
		      $(SRCDIR)/Math/hpc_summation.cuh \
		      $(SRCDIR)/MolecularMechanics/hpc_minimization.h \
		      $(SRCDIR)/Potential/hpc_nonbonded_potential.h \
		      $(SRCDIR)/Potential/hpc_scorecard.h \
		      $(SRCDIR)/Potential/hpc_valence_potential.h \
		      $(SRCDIR)/Random/hpc_random.h \
		      $(SRCDIR)/Random/hpc_random.cuh \
		      $(SRCDIR)/Structure/hpc_virtual_site_handling.h \
		      $(SRCDIR)/Synthesis/hpc_phasespace_synthesis.h \
		      $(SRCDIR)/Synthesis/hpc_phasespace_synthesis.cuh

# STORMM CUDA included files
STORMM_CUDA_INCLUDED_FILES = $(SRCDIR)/Math/conjugate_gradient.cui \
			     $(SRCDIR)/Math/rounding.cui \
			     $(SRCDIR)/Math/vector_formulas.cui \
			     $(SRCDIR)/MolecularMechanics/line_movement.cui \
			     $(SRCDIR)/Numerics/accumulation.cui \
			     $(SRCDIR)/Potential/nonbonded_potential_tilegroups.cui \
			     $(SRCDIR)/Potential/gbderivative_tilegroups.cui \
			     $(SRCDIR)/Potential/gbradii_tilegroups.cui \
			     $(SRCDIR)/Potential/valence_potential.cui \
			     $(SRCDIR)/Random/xor_shift_rng.cui \
			     $(SRCDIR)/Structure/virtual_site_placement.cui \
			     $(SRCDIR)/Structure/virtual_site_transmission.cui \
			     $(SRCDIR)/Trajectory/thermostat_utilities.cui

# STORMM CUDA object files
STORMM_CUDA_OBJS = $(SRCDIR)/Accelerator/hpc_config.o \
		   $(SRCDIR)/Math/hpc_reduction.o \
		   $(SRCDIR)/MolecularMechanics/hpc_minimization.o \
		   $(SRCDIR)/Potential/hpc_nonbonded_potential.o \
		   $(SRCDIR)/Potential/hpc_scorecard.o \
		   $(SRCDIR)/Potential/hpc_valence_potential.o \
		   $(SRCDIR)/Random/hpc_random.o \
		   $(SRCDIR)/Structure/hpc_virtual_site_handling.o \
		   $(SRCDIR)/Synthesis/hpc_implicit_solvent_workspace.o \
		   $(SRCDIR)/Synthesis/hpc_phasespace_synthesis.o

# Test programs using stormm
STORMM_TEST_PROGS = $(TESTDIR)/bin/test_unit_test \
		    $(TESTDIR)/bin/test_file_management \
		    $(TESTDIR)/bin/test_hybrid \
		    $(TESTDIR)/bin/test_parse \
		    $(TESTDIR)/bin/test_input \
		    $(TESTDIR)/bin/test_namelists \
		    $(TESTDIR)/bin/test_math \
		    $(TESTDIR)/bin/test_numerics \
		    $(TESTDIR)/bin/test_amber_prmtop \
		    $(TESTDIR)/bin/test_topology_mods \
		    $(TESTDIR)/bin/test_amber_coordinates \
		    $(TESTDIR)/bin/test_chemical_features \
		    $(TESTDIR)/bin/test_atommask \
		    $(TESTDIR)/bin/test_local_arrangement \
		    $(TESTDIR)/bin/test_isomerization \
		    $(TESTDIR)/bin/test_synthesis \
		    $(TESTDIR)/bin/test_atomgraph_synthesis \
		    $(TESTDIR)/bin/test_valence_evaluation \
		    $(TESTDIR)/bin/test_nonbonded_evaluation \
		    $(TESTDIR)/bin/test_generalized_born \
		    $(TESTDIR)/bin/test_neighbor_list \
		    $(TESTDIR)/bin/test_restraints \
		    $(TESTDIR)/bin/test_minimization \
		    $(TESTDIR)/bin/test_molecule_formats \
		    $(TESTDIR)/bin/test_mesh_construction

# Test programs using stormm.cuda
STORMM_TEST_CUDA_PROGS = $(TESTDIR)/bin/test_hpc_status \
			 $(TESTDIR)/bin/test_hpc_hybrid \
			 $(TESTDIR)/bin/test_hpc_math \
			 $(TESTDIR)/bin/test_hpc_synthesis \
			 $(TESTDIR)/bin/test_hpc_minimization \
			 $(TESTDIR)/bin/test_hpc_dynamics

# Benchmark programs using stormm
STORMM_BENCH_PROGS = $(BENCHDIR)/bin/valence

# Benchmark programs using stormm.cuda
STORMM_BENCH_CUDA_PROGS = $(BENCHDIR)/bin/accumulate \
			  $(BENCHDIR)/bin/test_nonperiodic_kernels

# Applications using stormm
STORMM_APPS = $(APPDIR)/bin/conformer.stormm \
	      $(APPDIR)/bin/dynamics.stormm \
	      $(APPDIR)/bin/ffrefine.stormm

# Applications using stormm.cuda
STORMM_CUDA_APPS = $(APPDIR)/bin/conformer.stormm.cuda

# Compilation variables
CC=g++
CUCC=nvcc
CUDA_INCLUDES = -I$(SRCDIR) -I${CUDA_HOME}/include
CUDA_LINKS = -L$(SRCDIR) -L${CUDA_HOME}/lib64 -L${CUDA_HOME}/lib64/stubs \
	     -lcurand -lcublas -lcusolver -lcudart -lcudadevrt -lnvidia-ml
CPP_FLAGS = -std=c++17 -fPIC -O0 -g
CUDA_FLAGS = -std=c++17 --compiler-options=-fPIC -O0 -g --ptxas-options="-v"
CUDA_DEFINES = -DSTORMM_USE_HPC -DSTORMM_USE_CUDA
CUDA_ARCHS = -gencode arch=compute_60,code=sm_60 \
	     -gencode arch=compute_61,code=sm_61 \
	     -gencode arch=compute_70,code=sm_70 \
	     -gencode arch=compute_75,code=sm_75 \
	     -gencode arch=compute_80,code=sm_80 \
	     -gencode arch=compute_86,code=sm_86

# Target: compile a C++ source file into a C++ object file
%.o : %.cpp $(STORMM_CPP_HEADERS) $(STORMM_TPP_FILES)
	@echo "[STORMM]  CC $<"
	$(VB)$(CC) $(CPP_FLAGS) $(CPP_DEFINES) -c -o $@ $< $(CPP_INCLUDES)

# Target: compile a CUDA source file into a CUDA object file
%.o : %.cu $(STORMM_CUDA_HEADERS) $(STORMM_CUDA_INCLUDED_FILES)
	@echo "[STORMM]  CUCC $<"
	$(VB)$(CUCC) $(CUDA_FLAGS) $(CUDA_DEFINES) $(CUDA_ARCHS) -c -o $@ $< $(CUDA_INCLUDES)

# Target: the STORMM C++ library (CPU-only executables)
$(LIBDIR)/libstormm.so : CPP_DEFINES =
$(LIBDIR)/libstormm.so : CPP_INCLUDES = -I$(SRCDIR)
$(LIBDIR)/libstormm.so : CPP_LINKS = -L$(SRCDIR)
$(LIBDIR)/libstormm.so : $(STORMM_CPP_OBJS)
	@echo "[STORMM]  Building C++ library..."
	$(VB)$(CC) $(CPP_FLAGS) $(CPP_DEFINES) -shared -o $@ $(STORMM_CPP_OBJS) $(CPP_LINKS)

# Target: the STORMM C++ / CUDA library (CPU / GPU Hybrid executables)
$(LIBDIR)/libstormm_cuda.so : CPP_DEFINES = -DSTORMM_USE_HPC -DSTORMM_USE_CUDA
$(LIBDIR)/libstormm_cuda.so : CPP_INCLUDES = -I$(SRCDIR) -I${CUDA_HOME}/include
$(LIBDIR)/libstormm_cuda.so : CPP_LINKS = -L$(SRCDIR) -L${CUDA_HOME}/lib64 \
					  -L${CUDA_HOME}/lib64/stubs -lcurand -lcublas -lcusolver \
					  -lcudart -lcudadevrt -lnvidia-ml
$(LIBDIR)/libstormm_cuda.so : $(STORMM_CPP_OBJS) $(STORMM_CUDA_OBJS)
	@echo "[STORMM]  Building CUDA library..."
	$(VB)$(CUCC) $(CUDA_FLAGS) $(CUDA_DEFINES) -shared -o $@ $(STORMM_CUDA_OBJS) \
	$(STORMM_CPP_OBJS) $(CUDA_LINKS)

# Target: unit test features program
$(TESTDIR)/bin/test_unit_test : $(LIBDIR)/libstormm.so $(TESTDIR)/UnitTesting/test_unit_test.cpp
	@echo "[STORMM]  Building test_unit_test..."
	$(VB)$(CC) $(CPP_FLAGS) -o $(TESTDIR)/bin/test_unit_test \
	  $(TESTDIR)/UnitTesting/test_unit_test.cpp -L$(LIBDIR) -I$(SRCDIR) -lstormm

# Target: file management features program
$(TESTDIR)/bin/test_file_management : $(LIBDIR)/libstormm.so \
				      $(TESTDIR)/FileManagement/test_file_management.cpp
	@echo "[STORMM]  Building test_file_management..."
	$(VB)$(CC) $(CPP_FLAGS) -o $(TESTDIR)/bin/test_file_management \
	  $(TESTDIR)/FileManagement/test_file_management.cpp -L$(LIBDIR) -I$(SRCDIR) -lstormm

# Target: hybrid object features program (CPU-only)
$(TESTDIR)/bin/test_hybrid : $(LIBDIR)/libstormm.so $(TESTDIR)/Accelerator/test_hybrid.cpp
	@echo "[STORMM]  Building test_hybrid..."
	$(VB)$(CC) $(CPP_FLAGS) -o $(TESTDIR)/bin/test_hybrid \
	  $(TESTDIR)/Accelerator/test_hybrid.cpp -L$(LIBDIR) -I$(SRCDIR) -lstormm

# Target: text parsing features program
$(TESTDIR)/bin/test_parse : $(LIBDIR)/libstormm.so $(TESTDIR)/Parsing/test_parse.cpp
	@echo "[STORMM]  Building test_parse..."
	$(VB)$(CC) $(CPP_FLAGS) -o $(TESTDIR)/bin/test_parse $(TESTDIR)/Parsing/test_parse.cpp \
	  -L$(LIBDIR) -I$(SRCDIR) -lstormm

# Target: namelist input features program
$(TESTDIR)/bin/test_input : $(LIBDIR)/libstormm.so $(TESTDIR)/Parsing/test_input.cpp
	@echo "[STORMM]  Building test_input..."
	$(VB)$(CC) $(CPP_FLAGS) -o $(TESTDIR)/bin/test_input $(TESTDIR)/Parsing/test_input.cpp \
	  -L$(LIBDIR) -I$(SRCDIR) -lstormm

# Target: namelist input features program
$(TESTDIR)/bin/test_namelists : $(LIBDIR)/libstormm.so $(TESTDIR)/Namelists/test_namelists.cpp
	@echo "[STORMM]  Building test_namelists..."
	$(VB)$(CC) $(CPP_FLAGS) -o $(TESTDIR)/bin/test_namelists \
	  $(TESTDIR)/Namelists/test_namelists.cpp -L$(LIBDIR) -I$(SRCDIR) -lstormm

# Target: math operations program
$(TESTDIR)/bin/test_math : $(LIBDIR)/libstormm.so $(TESTDIR)/Math/test_math.cpp
	@echo "[STORMM]  Building test_math..."
	$(VB)$(CC) $(CPP_FLAGS) -o $(TESTDIR)/bin/test_math $(TESTDIR)/Math/test_math.cpp \
	  -L$(LIBDIR) -I$(SRCDIR) -lstormm

# Target: prmtop reading program
$(TESTDIR)/bin/test_amber_prmtop : $(LIBDIR)/libstormm.so \
				   $(TESTDIR)/Topology/test_amber_prmtop.cpp \
				   $(TESTDIR)/Topology/test_amber_prmtop.h \
				   $(TESTDIR)/Topology/test_amber_prmtop.tpp
	@echo "[STORMM]  Building test_amber_prmtop..."
	$(VB)$(CC) $(CPP_FLAGS) -o $(TESTDIR)/bin/test_amber_prmtop \
	  $(TESTDIR)/Topology/test_amber_prmtop.cpp -L$(LIBDIR) -I$(SRCDIR) -lstormm

# Target: prmtop reading program
$(TESTDIR)/bin/test_topology_mods : $(LIBDIR)/libstormm.so \
				    $(TESTDIR)/Topology/test_topology_mods.cpp
	@echo "[STORMM]  Building test_topology_mods..."
	$(VB)$(CC) $(CPP_FLAGS) -o $(TESTDIR)/bin/test_topology_mods \
	  $(TESTDIR)/Topology/test_topology_mods.cpp -L$(LIBDIR) -I$(SRCDIR) -lstormm

# Target: trajectory and restart reading program
$(TESTDIR)/bin/test_amber_coordinates : $(LIBDIR)/libstormm.so \
					$(TESTDIR)/Trajectory/test_amber_coordinates.cpp
	@echo "[STORMM]  Building test_amber_coordinates..."
	$(VB)$(CC) $(CPP_FLAGS) -o $(TESTDIR)/bin/test_amber_coordinates \
	  $(TESTDIR)/Trajectory/test_amber_coordinates.cpp -L$(LIBDIR) -I$(SRCDIR) -lstormm

# Target: chemical feature perception
$(TESTDIR)/bin/test_chemical_features : $(LIBDIR)/libstormm.so \
					$(TESTDIR)/Chemistry/test_chemical_features.cpp
	@echo "[STORMM]  Building test_chemical_features..."
	$(VB)$(CC) $(CPP_FLAGS) -o $(TESTDIR)/bin/test_chemical_features \
	  $(TESTDIR)/Chemistry/test_chemical_features.cpp -L$(LIBDIR) -I$(SRCDIR) -lstormm

# Target: atom mask parsing
$(TESTDIR)/bin/test_atommask : $(LIBDIR)/libstormm.so $(TESTDIR)/Chemistry/test_atommask.cpp
	@echo "[STORMM]  Building test_atommask..."
	$(VB)$(CC) $(CPP_FLAGS) -o $(TESTDIR)/bin/test_atommask \
	  $(TESTDIR)/Chemistry/test_atommask.cpp -L$(LIBDIR) -I$(SRCDIR) -lstormm

# Target: coordinate reimaging and local structure
$(TESTDIR)/bin/test_local_arrangement : $(LIBDIR)/libstormm.so \
				       $(TESTDIR)/Structure/test_local_arrangement.cpp
	@echo "[STORMM]  Building test_local_arrangement..."
	$(VB)$(CC) $(CPP_FLAGS) -o $(TESTDIR)/bin/test_local_arrangement \
	  $(TESTDIR)/Structure/test_local_arrangement.cpp -L$(LIBDIR) -I$(SRCDIR) -lstormm

# Target: isomer generation and "broad brush strokes" manipulation
$(TESTDIR)/bin/test_isomerization : $(LIBDIR)/libstormm.so \
				    $(TESTDIR)/Structure/test_isomerization.cpp
	@echo "[STORMM]  Building test_isomerization..."
	$(VB)$(CC) $(CPP_FLAGS) -o $(TESTDIR)/bin/test_isomerization \
	  $(TESTDIR)/Structure/test_isomerization.cpp -L$(LIBDIR) -I$(SRCDIR) -lstormm

# Target: numerical precision model
$(TESTDIR)/bin/test_numerics : $(LIBDIR)/libstormm.so $(TESTDIR)/Math/test_numerics.cpp
	@echo "[STORMM]  Building test_numerics..."
	$(VB)$(CC) $(CPP_FLAGS) -o $(TESTDIR)/bin/test_numerics \
	  $(TESTDIR)/Math/test_numerics.cpp -L$(LIBDIR) -I$(SRCDIR) -lstormm

# Target: testing the collated coordinate / velocity / force object
$(TESTDIR)/bin/test_synthesis : $(LIBDIR)/libstormm.so \
				$(TESTDIR)/Synthesis/test_synthesis.cpp
	@echo "[STORMM]  Building test_synthesis..."
	$(VB)$(CC) $(CPP_FLAGS) -o $(TESTDIR)/bin/test_synthesis \
	  $(TESTDIR)/Synthesis/test_synthesis.cpp -L$(LIBDIR) -I$(SRCDIR) -lstormm

# Target: prmtop collating object test program
$(TESTDIR)/bin/test_atomgraph_synthesis : $(LIBDIR)/libstormm.so \
					  $(TESTDIR)/Synthesis/test_atomgraph_synthesis.cpp \
					  $(TESTDIR)/Synthesis/assemble_restraints.h \
					  $(TESTDIR)/Synthesis/assemble_restraints.cpp
	@echo "[STORMM]  Building test_atomgraph_synthesis..."
	$(VB)$(CC) $(CPP_FLAGS) -o $(TESTDIR)/bin/test_atomgraph_synthesis \
	  $(TESTDIR)/Synthesis/test_atomgraph_synthesis.cpp \
	  $(TESTDIR)/Synthesis/assemble_restraints.cpp -L$(LIBDIR) -I$(SRCDIR) -lstormm

# Target: valence term evaluation by the most basic CPU routines
$(TESTDIR)/bin/test_valence_evaluation : $(LIBDIR)/libstormm.so \
					 $(TESTDIR)/Potential/test_valence_evaluation.cpp
	@echo "[STORMM]  Building test_valence_evaluation..."
	$(VB)$(CC) $(CPP_FLAGS) -o $(TESTDIR)/bin/test_valence_evaluation \
	  $(TESTDIR)/Potential/test_valence_evaluation.cpp -L$(LIBDIR) -I$(SRCDIR) -lstormm

# Target: non-bonded interaction evaluation by the most basic CPU routines
$(TESTDIR)/bin/test_nonbonded_evaluation : $(LIBDIR)/libstormm.so \
					   $(TESTDIR)/Potential/test_nonbonded_evaluation.cpp
	@echo "[STORMM]  Building test_nonbonded_evaluation..."
	$(VB)$(CC) $(CPP_FLAGS) -o $(TESTDIR)/bin/test_nonbonded_evaluation \
	  $(TESTDIR)/Potential/test_nonbonded_evaluation.cpp -L$(LIBDIR) -I$(SRCDIR) -lstormm

# Target: Generalized Born implicit solvent evaluation by the most basic CPU routines
$(TESTDIR)/bin/test_generalized_born : $(LIBDIR)/libstormm.so \
				       $(TESTDIR)/Potential/test_generalized_born.cpp
	@echo "[STORMM]  Building test_generalized_born..."
	$(VB)$(CC) $(CPP_FLAGS) -o $(TESTDIR)/bin/test_generalized_born \
	  $(TESTDIR)/Potential/test_generalized_born.cpp -L$(LIBDIR) -I$(SRCDIR) -lstormm

# Target: non-bonded neighbor list construction and usage by CPU routines
$(TESTDIR)/bin/test_neighbor_list : $(LIBDIR)/libstormm.so \
				    $(TESTDIR)/Potential/test_neighbor_list.cpp
	@echo "[STORMM]  Building test_neighbor_list..."
	$(VB)$(CC) $(CPP_FLAGS) -o $(TESTDIR)/bin/test_neighbor_list \
	  $(TESTDIR)/Potential/test_neighbor_list.cpp -L$(LIBDIR) -I$(SRCDIR) -lstormm

# Target: restraint construction and execution by CPU routines
$(TESTDIR)/bin/test_restraints : $(LIBDIR)/libstormm.so \
				 $(TESTDIR)/Restraints/test_restraints.cpp
	@echo "[STORMM]  Building test_restraints..."
	$(VB)$(CC) $(CPP_FLAGS) -o $(TESTDIR)/bin/test_restraints \
	  $(TESTDIR)/Restraints/test_restraints.cpp -L$(LIBDIR) -I$(SRCDIR) -lstormm

# Target: molecular mechanics minimization and execution by CPU routines
$(TESTDIR)/bin/test_minimization : $(LIBDIR)/libstormm.so \
				   $(TESTDIR)/MolecularMechanics/test_minimization.cpp
	@echo "[STORMM]  Building test_minimization..."
	$(VB)$(CC) $(CPP_FLAGS) -o $(TESTDIR)/bin/test_minimization \
	  $(TESTDIR)/MolecularMechanics/test_minimization.cpp -L$(LIBDIR) -I$(SRCDIR) -lstormm

# Target: molecular mechanics minimization and execution by CPU routines
$(TESTDIR)/bin/test_molecule_formats : $(LIBDIR)/libstormm.so \
				       $(TESTDIR)/MoleculeFormat/test_molecule_formats.cpp
	@echo "[STORMM]  Building test_molecule_formats..."
	$(VB)$(CC) $(CPP_FLAGS) -o $(TESTDIR)/bin/test_molecule_formats \
	  $(TESTDIR)/MoleculeFormat/test_molecule_formats.cpp -L$(LIBDIR) -I$(SRCDIR) -lstormm

# Target: molecular mechanics minimization and execution by CPU routines
$(TESTDIR)/bin/test_mesh_construction : $(LIBDIR)/libstormm.so \
				        $(TESTDIR)/Structure/test_mesh_construction.cpp
	@echo "[STORMM]  Building test_mesh_construction..."
	$(VB)$(CC) $(CPP_FLAGS) -o $(TESTDIR)/bin/test_mesh_construction \
	  $(TESTDIR)/Structure/test_mesh_construction.cpp -L$(LIBDIR) -I$(SRCDIR) -lstormm

# Target: HPC detection
$(TESTDIR)/bin/test_hpc_status : $(LIBDIR)/libstormm_cuda.so \
				 $(TESTDIR)/Accelerator/test_hpc_status.cpp
	@echo "[STORMM]  Building test_hpc_status..."
	$(VB)$(CUCC) $(CUDA_FLAGS) $(CUDA_DEFINES) $(CUDA_ARCHS) \
	  -o $(TESTDIR)/bin/test_hpc_status $(TESTDIR)/Accelerator/test_hpc_status.cpp \
	  -L$(LIBDIR) -I$(SRCDIR) $(CUDA_LINKS) -lstormm_cuda

# Target: HPC memory functionality
$(TESTDIR)/bin/test_hpc_hybrid : $(LIBDIR)/libstormm_cuda.so \
				 $(TESTDIR)/Accelerator/test_hpc_hybrid.cpp
	@echo "[STORMM]  Building test_hpc_hybrid..."
	$(VB)$(CUCC) $(CUDA_FLAGS) $(CUDA_DEFINES) $(CUDA_ARCHS) \
	  -o $(TESTDIR)/bin/test_hpc_hybrid $(TESTDIR)/Accelerator/test_hpc_hybrid.cpp \
	  -L$(LIBDIR) -I$(SRCDIR) $(CUDA_LINKS) -lstormm_cuda

# Target: HPC basic math kernels
$(TESTDIR)/bin/test_hpc_math : $(LIBDIR)/libstormm_cuda.so $(TESTDIR)/Math/test_hpc_math.cu
	@echo "[STORMM]  Building test_hpc_math..."
	$(VB)$(CUCC) $(CUDA_FLAGS) $(CUDA_DEFINES) $(CUDA_ARCHS) \
	  -o $(TESTDIR)/bin/test_hpc_math $(TESTDIR)/Math/test_hpc_math.cu -L$(LIBDIR) \
	  -I$(SRCDIR) $(CUDA_LINKS) -lstormm_cuda

# Target: HPC molecular system synthesis and associated operations
$(TESTDIR)/bin/test_hpc_synthesis : $(LIBDIR)/libstormm_cuda.so \
				    $(TESTDIR)/Synthesis/test_hpc_synthesis.cu \
				    $(TESTDIR)/Synthesis/assemble_restraints.h \
				    $(TESTDIR)/Synthesis/assemble_restraints.cpp
	@echo "[STORMM]  Building test_hpc_synthesis..."
	$(VB)$(CUCC) $(CUDA_FLAGS) $(CUDA_DEFINES) $(CUDA_ARCHS) \
	  -o $(TESTDIR)/bin/test_hpc_synthesis $(TESTDIR)/Synthesis/test_hpc_synthesis.cu \
	  $(TESTDIR)/Synthesis/assemble_restraints.cpp -L$(LIBDIR) -I$(SRCDIR) $(CUDA_LINKS) \
	  -lstormm_cuda

# Target: HPC molecular system synthesis and associated operations
$(TESTDIR)/bin/test_hpc_minimization : $(LIBDIR)/libstormm_cuda.so \
				       $(TESTDIR)/MolecularMechanics/test_hpc_minimization.cu
	@echo "[STORMM]  Building test_hpc_minimization..."
	$(VB)$(CUCC) $(CUDA_FLAGS) $(CUDA_DEFINES) $(CUDA_ARCHS) \
	  -o $(TESTDIR)/bin/test_hpc_minimization \
	  $(TESTDIR)/MolecularMechanics/test_hpc_minimization.cu -L$(LIBDIR) -I$(SRCDIR) \
	  $(CUDA_LINKS) -lstormm_cuda

# Target: HPC molecular system synthesis and associated operations
$(TESTDIR)/bin/test_hpc_dynamics : $(LIBDIR)/libstormm_cuda.so \
				   $(TESTDIR)/MolecularMechanics/test_hpc_dynamics.cu
	@echo "[STORMM]  Building test_hpc_dynamics..."
	$(VB)$(CUCC) $(CUDA_FLAGS) $(CUDA_DEFINES) $(CUDA_ARCHS) \
	  -o $(TESTDIR)/bin/test_hpc_dynamics \
	  $(TESTDIR)/MolecularMechanics/test_hpc_dynamics.cu -L$(LIBDIR) -I$(SRCDIR) \
	  $(CUDA_LINKS) -lstormm_cuda

# Target: Benchmarking single, double, and fixed-precision computations of valence forces
$(BENCHDIR)/bin/valence : $(LIBDIR)/libstormm.so \
			  $(BENCHDIR)/ForceAccumulation/valence.cpp
	@echo "[STORMM]  Building valence benchmark..."
	$(VB)$(CC) $(CPP_FLAGS) -o $(BENCHDIR)/bin/valence \
	  $(BENCHDIR)/ForceAccumulation/valence.cpp -L$(LIBDIR) -I$(SRCDIR) -lstormm

# Target: Benchmarking split accumulation of fixed-precision values for speed and resource usage
$(BENCHDIR)/bin/accumulate : $(LIBDIR)/libstormm_cuda.so \
			     $(BENCHDIR)/ForceAccumulation/accumulate.cu
	@echo "[STORMM]  Building accumulation benchmark..."
	$(VB)$(CUCC) $(CUDA_FLAGS) $(CUDA_DEFINES) $(CUDA_ARCHS) -o $(BENCHDIR)/bin/accumulate \
	  $(BENCHDIR)/ForceAccumulation/accumulate.cu -L$(LIBDIR) -I$(SRCDIR) $(CUDA_LINKS) \
	  -lstormm_cuda

# Target: Benchmarking kernels needed for dynamics in isolated boundary conditions
$(BENCHDIR)/bin/test_nonperiodic_kernels : $(LIBDIR)/libstormm_cuda.so \
					   $(BENCHDIR)/KernelTesting/test_nonperiodic_kernels.cu
	@echo "[STORMM]  Building non-periodic kernel benchmark..."
	$(VB)$(CUCC) $(CUDA_FLAGS) $(CUDA_DEFINES) $(CUDA_ARCHS) -o \
	  $(BENCHDIR)/bin/test_nonperiodic_kernels \
	  $(BENCHDIR)/KernelTesting/test_nonperiodic_kernels.cu -L$(LIBDIR) -I$(SRCDIR) \
	  $(CUDA_LINKS) -lstormm_cuda

# Target: Conformer generation
$(APPDIR)/bin/conformer.stormm : $(LIBDIR)/libstormm.so $(APPDIR)/Conf/src/conformer.cpp \
				 $(APPDIR)/Conf/src/setup.cpp
	@echo "[STORMM]  Building conformer.stormm..."
	$(VB)$(CC) $(CPP_FLAGS) -o $(APPDIR)/bin/conformer.stormm \
	  $(APPDIR)/Conf/src/conformer.cpp $(APPDIR)/Conf/src/setup.cpp \
	  -L$(LIBDIR) -I$(SRCDIR) -lstormm

# Target: Molecular dynamics
$(APPDIR)/bin/dynamics.stormm : $(LIBDIR)/libstormm.so $(APPDIR)/Dyna/src/dynamics.cpp
	@echo "[STORMM]  Building dynamics.stormm..."
	$(VB)$(CC) $(CPP_FLAGS) -o $(APPDIR)/bin/dynamics.stormm \
	  $(APPDIR)/Dyna/src/dynamics.cpp -L$(LIBDIR) -I$(SRCDIR) -lstormm

# Target: Force field refinement
$(APPDIR)/bin/ffrefine.stormm : $(LIBDIR)/libstormm.so $(APPDIR)/Ffrn/src/ffrefine.cpp
	@echo "[STORMM]  Building ffrefine.stormm..."
	$(VB)$(CC) $(CPP_FLAGS) -o $(APPDIR)/bin/ffrefine.stormm \
	  $(APPDIR)/Ffrn/src/ffrefine.cpp -L$(LIBDIR) -I$(SRCDIR) -lstormm

# Target: Conformer generation with CUDA
$(APPDIR)/bin/conformer.stormm.cuda : $(LIBDIR)/libstormm_cuda.so \
				      $(APPDIR)/Conf/src/conformer.cpp \
				      $(APPDIR)/Conf/src/setup.cpp
	@echo "[STORMM]  Building conformer.stormm.cuda..."
	$(VB)$(CUCC) $(CUDA_FLAGS) $(CUDA_DEFINES) $(CUDA_ARCHS) -o \
	  $(APPDIR)/bin/conformer.stormm.cuda $(APPDIR)/Conf/src/conformer.cpp \
	  $(APPDIR)/Conf/src/setup.cpp -L$(LIBDIR) -I$(SRCDIR) -lstormm_cuda

install : $(LIBDIR)/libstormm.so

clean :
	@echo "[STORMM]  Cleaning CPU libraries"
	$(VB)if [ -e $(SRCDIR)/Accelerator/hybrid.o ] ; then /bin/rm $(SRCDIR)/*/*.o ; fi

cuda : $(LIBDIR)/libstormm_cuda.so

clean.cuda:
	@echo "[STORMM]  Cleaning GPU libraries"
	$(VB)if [ -e $(SRCDIR)/Accelerator/hybrid.o ] ; then /bin/rm $(SRCDIR)/*/*.o ; fi
	$(VB)if [ -e $(SRCDIR)/Accelerator/hpc_config.o ] ; then /bin/rm $(SRCDIR)/*/*.o ; fi

clean.bench:
	@echo "[STORMM]  Cleaning benchmark programs"
	$(VB)if [ -e $(BENCHDIR)/bin/valence ] ; then /bin/rm $(BENCHDIR)/bin/* ; fi

clean.apps:
	@echo "[STORMM]  Cleaning application programs"
	$(VB)if [ -e $(APPDIR)/bin/conformer.stormm ] ; then /bin/rm $(APPDIR)/bin/* ; fi

clean.apps.cuda:
	@echo "[STORMM]  Cleaning CUDA_based application programs"
	$(VB)if [ -e $(APPDIR)/bin/conformer.stormm.cuda ] ; then \
	  /bin/rm $(APPDIR)/bin/*.stormm.cuda ; \
	fi

test.exe : $(STORMM_TEST_PROGS)

test.cuda.exe : $(STORMM_TEST_CUDA_PROGS)

test : $(STORMM_TEST_PROGS)
	for PROG in $(STORMM_TEST_PROGS) ; do \
		echo "[STORMM] Execute $$PROG" ; \
		$$PROG ; \
	done

test.cuda : $(STORMM_TEST_CUDA_PROGS)
	for PROG in $(STORMM_TEST_CUDA_PROGS) ; do \
		echo "[STORMM] Execute $$PROG" ; \
		$$PROG ; \
	done

bench.exe : $(STORMM_BENCH_PROGS)

bench.cuda.exe : $(STORMM_BENCH_CUDA_PROGS)

bench : $(STORMM_BENCH_PROGS)
	for PROG in $(STORMM_BENCH_PROGS) ; do \
		echo "[STORMM] Run benchmark $$PROG" ; \
		$$PROG ; \
	done

bench.cuda : $(STORMM_BENCH_CUDA_PROGS)
	for PROG in $(STORMM_BENCH_CUDA_PROGS) ; do \
		echo "[STORMM] Run benchmark $$PROG" ; \
		$$PROG ; \
	done

apps : $(STORMM_APPS)

apps.cuda : $(STORMM_CUDA_APPS)

yes:  install
