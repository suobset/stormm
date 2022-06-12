# Verbosity
VB=@

# Source and test directories
SRCDIR=src
TESTDIR=test
BENCHDIR=benchmark
APPDIR=apps
LIBDIR=lib

# OMNI C++ source files
OMNI_CPP_FILES = $(SRCDIR)/Accelerator/hybrid.cpp \
		 $(SRCDIR)/Chemistry/atommask.cpp \
		 $(SRCDIR)/Chemistry/chemical_features.cpp \
		 $(SRCDIR)/Chemistry/indigo.cpp \
		 $(SRCDIR)/Chemistry/znumber.cpp \
		 $(SRCDIR)/Constants/fixed_precision.cpp \
		 $(SRCDIR)/Constants/generalized_born.cpp \
		 $(SRCDIR)/FileManagement/directory_util.cpp \
		 $(SRCDIR)/FileManagement/file_listing.cpp \
		 $(SRCDIR)/FileManagement/file_util.cpp \
		 $(SRCDIR)/ForceField/forcefield_element.cpp \
		 $(SRCDIR)/ForceField/forcefield_enumerators.cpp \
		 $(SRCDIR)/Math/matrix_ops.cpp \
		 $(SRCDIR)/Math/rounding.cpp \
		 $(SRCDIR)/Math/series_ops.cpp \
		 $(SRCDIR)/Math/statistics.cpp \
		 $(SRCDIR)/Math/vector_ops.cpp \
		 $(SRCDIR)/MolecularMechanics/minimization.cpp \
		 $(SRCDIR)/MolecularMechanics/mm_controls.cpp \
		 $(SRCDIR)/MolecularMechanics/mm_evaluation.cpp \
		 $(SRCDIR)/MoleculeFormat/mdl_mol_format.cpp \
	         $(SRCDIR)/Namelists/input.cpp \
	         $(SRCDIR)/Namelists/namelist_element.cpp \
	         $(SRCDIR)/Namelists/namelist_emulator.cpp \
	         $(SRCDIR)/Namelists/namelist_enumerations.cpp \
	         $(SRCDIR)/Namelists/nml_conformer.cpp \
	         $(SRCDIR)/Namelists/nml_dynamics.cpp \
	         $(SRCDIR)/Namelists/nml_ffmorph.cpp \
	         $(SRCDIR)/Namelists/nml_files.cpp \
	         $(SRCDIR)/Namelists/nml_minimize.cpp \
	         $(SRCDIR)/Namelists/nml_precision.cpp \
	         $(SRCDIR)/Namelists/nml_random.cpp \
	         $(SRCDIR)/Namelists/nml_restraint.cpp \
	         $(SRCDIR)/Namelists/nml_solvent.cpp \
	         $(SRCDIR)/Namelists/user_settings.cpp \
	         $(SRCDIR)/Parsing/ascii_numbers.cpp \
	         $(SRCDIR)/Parsing/citation.cpp \
		 $(SRCDIR)/Parsing/parse.cpp \
		 $(SRCDIR)/Parsing/polynumeric.cpp \
	         $(SRCDIR)/Parsing/tabulation.cpp \
	         $(SRCDIR)/Parsing/textfile.cpp \
	         $(SRCDIR)/Parsing/textguard.cpp \
	         $(SRCDIR)/Potential/cacheresource.cpp \
	         $(SRCDIR)/Potential/eval_valence_workunit.cpp \
	         $(SRCDIR)/Potential/forward_exclusionmask.cpp \
	         $(SRCDIR)/Potential/nonbonded_potential.cpp \
	         $(SRCDIR)/Potential/scorecard.cpp \
	         $(SRCDIR)/Potential/static_exclusionmask.cpp \
	         $(SRCDIR)/Potential/valence_potential.cpp \
		 $(SRCDIR)/Random/random.cpp \
	         $(SRCDIR)/Reporting/code_dox.cpp \
	         $(SRCDIR)/Reporting/display.cpp \
		 $(SRCDIR)/Reporting/error_format.cpp \
		 $(SRCDIR)/Restraints/bounded_restraint.cpp \
		 $(SRCDIR)/Restraints/restraint_apparatus.cpp \
		 $(SRCDIR)/Restraints/restraint_builder.cpp \
		 $(SRCDIR)/Structure/global_manipulation.cpp \
		 $(SRCDIR)/Structure/isomerization.cpp \
		 $(SRCDIR)/Structure/local_arrangement.cpp \
		 $(SRCDIR)/Structure/rmsd.cpp \
		 $(SRCDIR)/Structure/structure_ops.cpp \
		 $(SRCDIR)/Structure/structure_utils.cpp \
		 $(SRCDIR)/Structure/virtual_site_handling.cpp \
		 $(SRCDIR)/Synthesis/atomgraph_synthesis.cpp \
		 $(SRCDIR)/Synthesis/phasespace_synthesis.cpp \
		 $(SRCDIR)/Synthesis/nonbonded_workunit.cpp \
		 $(SRCDIR)/Synthesis/static_mask_synthesis.cpp \
		 $(SRCDIR)/Synthesis/systemcache.cpp \
		 $(SRCDIR)/Synthesis/valence_workunit.cpp \
	         $(SRCDIR)/Topology/amber_prmtop_util.cpp \
	         $(SRCDIR)/Topology/atomgraph_constructors.cpp \
	         $(SRCDIR)/Topology/atomgraph_detailers.cpp \
	         $(SRCDIR)/Topology/atomgraph_getters.cpp \
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
	         $(SRCDIR)/Trajectory/phasespace.cpp \
	         $(SRCDIR)/Trajectory/thermostat.cpp \
	         $(SRCDIR)/Trajectory/trajectory_enumerators.cpp \
	         $(SRCDIR)/Trajectory/write_frame.cpp \
		 $(SRCDIR)/UnitTesting/approx.cpp \
		 $(SRCDIR)/UnitTesting/stopwatch.cpp \
		 $(SRCDIR)/UnitTesting/checklist.cpp \
		 $(SRCDIR)/UnitTesting/file_snapshot.cpp \
		 $(SRCDIR)/UnitTesting/test_environment.cpp \
		 $(SRCDIR)/UnitTesting/unit_test.cpp \
		 $(SRCDIR)/UnitTesting/vector_report.cpp

# OMNI C++ header files
OMNI_CPP_HEADERS = $(SRCDIR)/Accelerator/hybrid.h \
		   $(SRCDIR)/Chemistry/atommask.h \
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
		   $(SRCDIR)/DataTypes/omni_vector_types.h \
		   $(SRCDIR)/FileManagement/directory_util.h \
		   $(SRCDIR)/FileManagement/file_listing.h \
		   $(SRCDIR)/FileManagement/file_util.h \
		   $(SRCDIR)/ForceField/forcefield_element.h \
		   $(SRCDIR)/ForceField/forcefield_enumerators.h \
		   $(SRCDIR)/Math/matrix.h \
		   $(SRCDIR)/Math/matrix_ops.h \
		   $(SRCDIR)/Math/rounding.h \
		   $(SRCDIR)/Math/series_ops.h \
		   $(SRCDIR)/Math/sorting.h \
	           $(SRCDIR)/Math/statistics.h \
	           $(SRCDIR)/Math/summation.h \
	           $(SRCDIR)/Math/vector_ops.h \
		   $(SRCDIR)/MolecularMechanics/minimization.h \
		   $(SRCDIR)/MolecularMechanics/mm_controls.h \
		   $(SRCDIR)/MolecularMechanics/mm_evaluation.h \
		   $(SRCDIR)/MoleculeFormat/mdl_mol_format.h \
		   $(SRCDIR)/Namelists/input.h \
		   $(SRCDIR)/Namelists/namelist_element.h \
		   $(SRCDIR)/Namelists/namelist_emulator.h \
		   $(SRCDIR)/Namelists/namelist_enumerations.h \
		   $(SRCDIR)/Namelists/nml_conformer.h \
		   $(SRCDIR)/Namelists/nml_dynamics.h \
		   $(SRCDIR)/Namelists/nml_ffmorph.h \
		   $(SRCDIR)/Namelists/nml_files.h \
		   $(SRCDIR)/Namelists/nml_minimize.h \
		   $(SRCDIR)/Namelists/nml_precision.h \
		   $(SRCDIR)/Namelists/nml_random.h \
		   $(SRCDIR)/Namelists/nml_restraint.h \
		   $(SRCDIR)/Namelists/nml_solvent.h \
		   $(SRCDIR)/Namelists/user_settings.h \
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
	           $(SRCDIR)/Reporting/error_format.h \
	           $(SRCDIR)/Restraints/bounded_restraint.h \
	           $(SRCDIR)/Restraints/restraint_apparatus.h \
	           $(SRCDIR)/Restraints/restraint_builder.h \
	           $(SRCDIR)/Restraints/restraint_util.h \
		   $(SRCDIR)/Structure/global_manipulation.h \
		   $(SRCDIR)/Structure/isomerization.h \
		   $(SRCDIR)/Structure/local_arrangement.h \
		   $(SRCDIR)/Structure/structure_enumerators.h \
		   $(SRCDIR)/Structure/rmsd.h \
		   $(SRCDIR)/Structure/structure_ops.h \
		   $(SRCDIR)/Structure/structure_utils.h \
		   $(SRCDIR)/Structure/virtual_site_handling.h \
		   $(SRCDIR)/Synthesis/atomgraph_synthesis.h \
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
		   $(SRCDIR)/Topology/atomgraph_refinement.h \
		   $(SRCDIR)/Topology/topology_limits.h \
		   $(SRCDIR)/Topology/topology_util.h \
		   $(SRCDIR)/Trajectory/amber_ascii.h \
		   $(SRCDIR)/Trajectory/barostat.h \
		   $(SRCDIR)/Trajectory/coordinateframe.h \
		   $(SRCDIR)/Trajectory/coordinate_series.h \
		   $(SRCDIR)/Trajectory/phasespace.h \
		   $(SRCDIR)/Trajectory/thermostat.h \
		   $(SRCDIR)/Trajectory/trajectory_enumerators.h \
		   $(SRCDIR)/Trajectory/write_frame.h \
		   $(SRCDIR)/UnitTesting/approx.h \
		   $(SRCDIR)/UnitTesting/stopwatch.h \
		   $(SRCDIR)/UnitTesting/checklist.h \
		   $(SRCDIR)/UnitTesting/file_snapshot.h \
		   $(SRCDIR)/UnitTesting/test_environment.h \
		   $(SRCDIR)/UnitTesting/unit_test.h \
		   $(SRCDIR)/UnitTesting/unit_test_enumerators.h \
		   $(SRCDIR)/UnitTesting/vector_report.h

# OMNI C++ template source files
OMNI_TPP_FILES = $(SRCDIR)/Accelerator/hybrid.tpp \
		 $(SRCDIR)/Constants/generalized_born.tpp \
		 $(SRCDIR)/DataTypes/common_types.tpp \
		 $(SRCDIR)/DataTypes/omni_vector_types.tpp \
		 $(SRCDIR)/Math/matrix.tpp \
		 $(SRCDIR)/Math/matrix_ops.tpp \
		 $(SRCDIR)/Math/rounding.tpp \
		 $(SRCDIR)/Math/series_ops.tpp \
		 $(SRCDIR)/Math/summation.tpp \
		 $(SRCDIR)/Math/vector_ops.tpp \
		 $(SRCDIR)/MolecularMechanics/minimization.tpp \
		 $(SRCDIR)/MolecularMechanics/mm_controls.tpp \
		 $(SRCDIR)/MolecularMechanics/mm_evaluation.tpp \
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
		 $(SRCDIR)/Structure/rmsd.tpp \
		 $(SRCDIR)/Structure/virtual_site_handling.tpp \
		 $(SRCDIR)/Synthesis/nonbonded_workunit.tpp \
		 $(SRCDIR)/Synthesis/synthesis_abstracts.tpp \
		 $(SRCDIR)/Topology/atomgraph.tpp \
		 $(SRCDIR)/Topology/atomgraph_abstracts.tpp \
		 $(SRCDIR)/Topology/topology_util.tpp \
		 $(SRCDIR)/Trajectory/coordinate_series.tpp \
		 $(SRCDIR)/UnitTesting/approx.tpp \
		 $(SRCDIR)/UnitTesting/unit_test.tpp

# OMNI C++ object files
OMNI_CPP_OBJS = $(SRCDIR)/Accelerator/hybrid.o \
		$(SRCDIR)/Chemistry/atommask.o \
		$(SRCDIR)/Chemistry/chemical_features.o \
		$(SRCDIR)/Chemistry/indigo.o \
		$(SRCDIR)/Chemistry/znumber.o \
		$(SRCDIR)/Constants/fixed_precision.o \
		$(SRCDIR)/Constants/generalized_born.o \
		$(SRCDIR)/FileManagement/directory_util.o \
		$(SRCDIR)/FileManagement/file_listing.o \
		$(SRCDIR)/FileManagement/file_util.o \
		$(SRCDIR)/ForceField/forcefield_element.o \
		$(SRCDIR)/ForceField/forcefield_enumerators.o \
		$(SRCDIR)/Math/matrix_ops.o \
		$(SRCDIR)/Math/rounding.o \
	        $(SRCDIR)/Math/series_ops.o \
	        $(SRCDIR)/Math/statistics.o \
	        $(SRCDIR)/Math/vector_ops.o \
	        $(SRCDIR)/MolecularMechanics/minimization.o \
	        $(SRCDIR)/MolecularMechanics/mm_controls.o \
	        $(SRCDIR)/MolecularMechanics/mm_evaluation.o \
	        $(SRCDIR)/MoleculeFormat/mdl_mol_format.o \
		$(SRCDIR)/Namelists/input.o \
	        $(SRCDIR)/Namelists/namelist_element.o \
	        $(SRCDIR)/Namelists/namelist_emulator.o \
	        $(SRCDIR)/Namelists/namelist_enumerations.o \
	        $(SRCDIR)/Namelists/nml_conformer.o \
	        $(SRCDIR)/Namelists/nml_dynamics.o \
	        $(SRCDIR)/Namelists/nml_ffmorph.o \
	        $(SRCDIR)/Namelists/nml_files.o \
	        $(SRCDIR)/Namelists/nml_minimize.o \
	        $(SRCDIR)/Namelists/nml_precision.o \
	        $(SRCDIR)/Namelists/nml_random.o \
	        $(SRCDIR)/Namelists/nml_restraint.o \
	        $(SRCDIR)/Namelists/nml_solvent.o \
	        $(SRCDIR)/Namelists/user_settings.o \
		$(SRCDIR)/Parsing/ascii_numbers.o \
	        $(SRCDIR)/Parsing/citation.o \
	        $(SRCDIR)/Parsing/parse.o \
	        $(SRCDIR)/Parsing/polynumeric.o \
	        $(SRCDIR)/Parsing/tabulation.o \
	        $(SRCDIR)/Parsing/textfile.o \
	        $(SRCDIR)/Parsing/textguard.o \
	        $(SRCDIR)/Potential/cacheresource.o \
	        $(SRCDIR)/Potential/eval_valence_workunit.o \
	        $(SRCDIR)/Potential/forward_exclusionmask.o \
	        $(SRCDIR)/Potential/nonbonded_potential.o \
	        $(SRCDIR)/Potential/scorecard.o \
	        $(SRCDIR)/Potential/static_exclusionmask.o \
	        $(SRCDIR)/Potential/valence_potential.o \
		$(SRCDIR)/Random/random.o \
		$(SRCDIR)/Reporting/code_dox.o \
		$(SRCDIR)/Reporting/display.o \
	        $(SRCDIR)/Reporting/error_format.o \
	        $(SRCDIR)/Restraints/bounded_restraint.o \
	        $(SRCDIR)/Restraints/restraint_apparatus.o \
	        $(SRCDIR)/Restraints/restraint_builder.o \
		$(SRCDIR)/Structure/global_manipulation.o \
		$(SRCDIR)/Structure/isomerization.o \
		$(SRCDIR)/Structure/local_arrangement.o \
		$(SRCDIR)/Structure/rmsd.o \
		$(SRCDIR)/Structure/structure_ops.o \
		$(SRCDIR)/Structure/structure_utils.o \
		$(SRCDIR)/Structure/virtual_site_handling.o \
		$(SRCDIR)/Synthesis/atomgraph_synthesis.o \
		$(SRCDIR)/Synthesis/phasespace_synthesis.o \
		$(SRCDIR)/Synthesis/nonbonded_workunit.o \
		$(SRCDIR)/Synthesis/static_mask_synthesis.o \
		$(SRCDIR)/Synthesis/systemcache.o \
		$(SRCDIR)/Synthesis/valence_workunit.o \
		$(SRCDIR)/Topology/amber_prmtop_util.o \
		$(SRCDIR)/Topology/atomgraph_constructors.o \
		$(SRCDIR)/Topology/atomgraph_detailers.o \
		$(SRCDIR)/Topology/atomgraph_getters.o \
		$(SRCDIR)/Topology/atomgraph_setters.o \
		$(SRCDIR)/Topology/atomgraph_miscellaneous.o \
		$(SRCDIR)/Topology/atomgraph_abstracts.o \
		$(SRCDIR)/Topology/atomgraph_analysis.o \
		$(SRCDIR)/Topology/atomgraph_enumerators.o \
		$(SRCDIR)/Topology/atomgraph_refinement.o \
		$(SRCDIR)/Topology/topology_limits.o \
		$(SRCDIR)/Topology/topology_util.o \
		$(SRCDIR)/Trajectory/amber_ascii.o \
		$(SRCDIR)/Trajectory/barostat.o \
		$(SRCDIR)/Trajectory/coordinateframe.o \
		$(SRCDIR)/Trajectory/phasespace.o \
		$(SRCDIR)/Trajectory/thermostat.o \
		$(SRCDIR)/Trajectory/trajectory_enumerators.o \
		$(SRCDIR)/Trajectory/write_frame.o \
		$(SRCDIR)/UnitTesting/approx.o \
		$(SRCDIR)/UnitTesting/stopwatch.o \
		$(SRCDIR)/UnitTesting/checklist.o \
		$(SRCDIR)/UnitTesting/file_snapshot.o \
		$(SRCDIR)/UnitTesting/test_environment.o \
		$(SRCDIR)/UnitTesting/unit_test.o \
		$(SRCDIR)/UnitTesting/vector_report.o

# OMNI CUDA source files
OMNI_CUDA_FILES = $(SRCDIR)/Accelerator/hpc_config.cu \
		  $(SRCDIR)/Potential/hpc_valence_potential.cu \
		  $(SRCDIR)/Random/hpc_random.cu \
		  $(SRCDIR)/Synthesis/hpc_phasespace_synthesis.cu

# OMNI CUDA header files
OMNI_CUDA_HEADERS = $(SRCDIR)/Constants/hpc_bounds.h \
		    $(SRCDIR)/Accelerator/gpu_details.h \
		    $(SRCDIR)/Accelerator/hpc_config.cuh \
		    $(SRCDIR)/Accelerator/ptx_macros.h \
		    $(SRCDIR)/Math/hpc_summation.cuh \
		    $(SRCDIR)/Potential/hpc_valence_potential.cuh \
		    $(SRCDIR)/Random/hpc_random.cuh \
		    $(SRCDIR)/Synthesis/hpc_phasespace_synthesis.h \
		    $(SRCDIR)/Synthesis/hpc_phasespace_synthesis.cuh

# OMNI CUDA included files
OMNI_CUDA_INCLUDED_FILES = $(SRCDIR)/Math/rounding.cui \
			   $(SRCDIR)/Potential/accumulation.cui \
			   $(SRCDIR)/Potential/valence_potential.cui

# OMNI CUDA object files
OMNI_CUDA_OBJS = $(SRCDIR)/Accelerator/hpc_config.o \
		 $(SRCDIR)/Potential/hpc_valence_potential.o \
		 $(SRCDIR)/Random/hpc_random.o \
		 $(SRCDIR)/Synthesis/hpc_phasespace_synthesis.o

# Test programs using omni
OMNI_TEST_PROGS = $(TESTDIR)/bin/test_unit_test \
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
	          $(TESTDIR)/bin/test_minimization

# Test programs using omni.cuda
OMNI_TEST_CUDA_PROGS = $(TESTDIR)/bin/test_hpc_status \
		       $(TESTDIR)/bin/test_hpc_hybrid \
		       $(TESTDIR)/bin/test_hpc_math \
		       $(TESTDIR)/bin/test_hpc_synthesis

# Benchmark programs using omni
OMNI_BENCH_PROGS = $(BENCHDIR)/bin/valence

# Benchmark programs using omni.cuda
OMNI_BENCH_CUDA_PROGS = $(BENCHDIR)/bin/accumulate

# Applications using omni
OMNI_APPS = $(APPDIR)/bin/conformer.omni \
	    $(APPDIR)/bin/dynamics.omni \
	    $(APPDIR)/bin/ffrefine.omni

# Compilation variables
CC=g++
CUCC=nvcc
CUDA_INCLUDES = -I$(SRCDIR) -I${CUDA_HOME}/include
CUDA_LINKS = -L$(SRCDIR) -L${CUDA_HOME}/lib64 -L${CUDA_HOME}/lib64/stubs \
	     -lcurand -lcublas -lcusolver -lcudart -lcudadevrt -lnvidia-ml
CPP_FLAGS = -std=c++17 -fPIC -O0 -g
CUDA_FLAGS = -std=c++17 --compiler-options=-fPIC -O0 -g --ptxas-options="-v"
CUDA_DEFINES = -DOMNI_USE_HPC -DOMNI_USE_CUDA
CUDA_ARCHS = -gencode arch=compute_60,code=sm_60 \
             -gencode arch=compute_61,code=sm_61 \
             -gencode arch=compute_70,code=sm_70 \
             -gencode arch=compute_75,code=sm_75 \
             -gencode arch=compute_80,code=sm_80 \
             -gencode arch=compute_80,code=sm_86

# Target: compile a C++ source file into a C++ object file
%.o : %.cpp $(OMNI_CPP_HEADERS) $(OMNI_TPP_FILES)
	@echo "[OMNI]  CC $<"
	$(VB)$(CC) $(CPP_FLAGS) $(CPP_DEFINES) -c -o $@ $< $(CPP_INCLUDES)

# Target: compile a CUDA source file into a CUDA object file
%.o : %.cu $(OMNI_CUDA_HEADERS) $(OMNI_CUDA_INCLUDED_FILES)
	@echo "[OMNI]  CUCC $<"
	$(VB)$(CUCC) $(CUDA_FLAGS) $(CUDA_DEFINES) $(CUDA_ARCHS) -c -o $@ $< $(CUDA_INCLUDES)

# Target: the OMNI C++ library (CPU-only executables)
$(LIBDIR)/libomni.so : CPP_DEFINES =
$(LIBDIR)/libomni.so : CPP_INCLUDES = -I$(SRCDIR)
$(LIBDIR)/libomni.so : CPP_LINKS = -L$(SRCDIR)
$(LIBDIR)/libomni.so : $(OMNI_CPP_OBJS)
	@echo "[OMNI]  Building C++ library..."
	$(VB)$(CC) $(CPP_FLAGS) $(CPP_DEFINES) -shared -o $@ $(OMNI_CPP_OBJS) $(CPP_LINKS)

# Target: the OMNI C++ / CUDA library (CPU / GPU Hybrid executables)
$(LIBDIR)/libomni_cuda.so : CPP_DEFINES = -DOMNI_USE_HPC -DOMNI_USE_CUDA
$(LIBDIR)/libomni_cuda.so : CPP_INCLUDES = -I$(SRCDIR) -I${CUDA_HOME}/include
$(LIBDIR)/libomni_cuda.so : CPP_LINKS = -L$(SRCDIR) -L${CUDA_HOME}/lib64 \
					-L${CUDA_HOME}/lib64/stubs -lcurand -lcublas -lcusolver \
					-lcudart -lcudadevrt -lnvidia-ml
$(LIBDIR)/libomni_cuda.so : $(OMNI_CPP_OBJS) $(OMNI_CUDA_OBJS)
	@echo "[OMNI]  Building CUDA library..."
	$(VB)$(CUCC) $(CUDA_FLAGS) $(CUDA_DEFINES) -shared -o $@ $(OMNI_CUDA_OBJS) \
	$(OMNI_CPP_OBJS) $(CUDA_LINKS)

# Target: unit test features program
$(TESTDIR)/bin/test_unit_test : $(LIBDIR)/libomni.so $(TESTDIR)/UnitTesting/test_unit_test.cpp
	@echo "[OMNI]  Building test_unit_test..."
	$(VB)$(CC) $(CPP_FLAGS) -o $(TESTDIR)/bin/test_unit_test \
	  $(TESTDIR)/UnitTesting/test_unit_test.cpp -L$(LIBDIR) -I$(SRCDIR) -lomni

# Target: file management features program
$(TESTDIR)/bin/test_file_management : $(LIBDIR)/libomni.so \
				      $(TESTDIR)/FileManagement/test_file_management.cpp
	@echo "[OMNI]  Building test_file_management..."
	$(VB)$(CC) $(CPP_FLAGS) -o $(TESTDIR)/bin/test_file_management \
	  $(TESTDIR)/FileManagement/test_file_management.cpp -L$(LIBDIR) -I$(SRCDIR) -lomni

# Target: hybrid object features program (CPU-only)
$(TESTDIR)/bin/test_hybrid : $(LIBDIR)/libomni.so $(TESTDIR)/Accelerator/test_hybrid.cpp
	@echo "[OMNI]  Building test_hybrid..."
	$(VB)$(CC) $(CPP_FLAGS) -o $(TESTDIR)/bin/test_hybrid \
	  $(TESTDIR)/Accelerator/test_hybrid.cpp -L$(LIBDIR) -I$(SRCDIR) -lomni

# Target: text parsing features program
$(TESTDIR)/bin/test_parse : $(LIBDIR)/libomni.so $(TESTDIR)/Parsing/test_parse.cpp
	@echo "[OMNI]  Building test_parse..."
	$(VB)$(CC) $(CPP_FLAGS) -o $(TESTDIR)/bin/test_parse $(TESTDIR)/Parsing/test_parse.cpp \
	  -L$(LIBDIR) -I$(SRCDIR) -lomni

# Target: namelist input features program
$(TESTDIR)/bin/test_input : $(LIBDIR)/libomni.so $(TESTDIR)/Parsing/test_input.cpp
	@echo "[OMNI]  Building test_input..."
	$(VB)$(CC) $(CPP_FLAGS) -o $(TESTDIR)/bin/test_input $(TESTDIR)/Parsing/test_input.cpp \
	  -L$(LIBDIR) -I$(SRCDIR) -lomni

# Target: namelist input features program
$(TESTDIR)/bin/test_namelists : $(LIBDIR)/libomni.so $(TESTDIR)/Namelists/test_namelists.cpp
	@echo "[OMNI]  Building test_namelists..."
	$(VB)$(CC) $(CPP_FLAGS) -o $(TESTDIR)/bin/test_namelists \
	  $(TESTDIR)/Namelists/test_namelists.cpp -L$(LIBDIR) -I$(SRCDIR) -lomni

# Target: math operations program
$(TESTDIR)/bin/test_math : $(LIBDIR)/libomni.so $(TESTDIR)/Math/test_math.cpp
	@echo "[OMNI]  Building test_math..."
	$(VB)$(CC) $(CPP_FLAGS) -o $(TESTDIR)/bin/test_math $(TESTDIR)/Math/test_math.cpp \
	  -L$(LIBDIR) -I$(SRCDIR) -lomni

# Target: prmtop reading program
$(TESTDIR)/bin/test_amber_prmtop : $(LIBDIR)/libomni.so $(TESTDIR)/Topology/test_amber_prmtop.cpp \
				   $(TESTDIR)/Topology/test_amber_prmtop.h \
				   $(TESTDIR)/Topology/test_amber_prmtop.tpp
	@echo "[OMNI]  Building test_amber_prmtop..."
	$(VB)$(CC) $(CPP_FLAGS) -o $(TESTDIR)/bin/test_amber_prmtop \
	  $(TESTDIR)/Topology/test_amber_prmtop.cpp -L$(LIBDIR) -I$(SRCDIR) -lomni

# Target: prmtop reading program
$(TESTDIR)/bin/test_topology_mods : $(LIBDIR)/libomni.so $(TESTDIR)/Topology/test_topology_mods.cpp
	@echo "[OMNI]  Building test_topology_mods..."
	$(VB)$(CC) $(CPP_FLAGS) -o $(TESTDIR)/bin/test_topology_mods \
	  $(TESTDIR)/Topology/test_topology_mods.cpp -L$(LIBDIR) -I$(SRCDIR) -lomni

# Target: trajectory and restart reading program
$(TESTDIR)/bin/test_amber_coordinates : $(LIBDIR)/libomni.so \
					$(TESTDIR)/Trajectory/test_amber_coordinates.cpp
	@echo "[OMNI]  Building test_amber_coordinates..."
	$(VB)$(CC) $(CPP_FLAGS) -o $(TESTDIR)/bin/test_amber_coordinates \
	  $(TESTDIR)/Trajectory/test_amber_coordinates.cpp -L$(LIBDIR) -I$(SRCDIR) -lomni

# Target: chemical feature perception
$(TESTDIR)/bin/test_chemical_features : $(LIBDIR)/libomni.so \
					$(TESTDIR)/Chemistry/test_chemical_features.cpp
	@echo "[OMNI]  Building test_chemical_features..."
	$(VB)$(CC) $(CPP_FLAGS) -o $(TESTDIR)/bin/test_chemical_features \
	  $(TESTDIR)/Chemistry/test_chemical_features.cpp -L$(LIBDIR) -I$(SRCDIR) -lomni

# Target: atom mask parsing
$(TESTDIR)/bin/test_atommask : $(LIBDIR)/libomni.so $(TESTDIR)/Chemistry/test_atommask.cpp
	@echo "[OMNI]  Building test_atommask..."
	$(VB)$(CC) $(CPP_FLAGS) -o $(TESTDIR)/bin/test_atommask \
	  $(TESTDIR)/Chemistry/test_atommask.cpp -L$(LIBDIR) -I$(SRCDIR) -lomni

# Target: coordinate reimaging and local structure
$(TESTDIR)/bin/test_local_arrangement : $(LIBDIR)/libomni.so \
				       $(TESTDIR)/Structure/test_local_arrangement.cpp
	@echo "[OMNI]  Building test_local_arrangement..."
	$(VB)$(CC) $(CPP_FLAGS) -o $(TESTDIR)/bin/test_local_arrangement \
	  $(TESTDIR)/Structure/test_local_arrangement.cpp -L$(LIBDIR) -I$(SRCDIR) -lomni

# Target: isomer generation and "broad brush strokes" manipulation
$(TESTDIR)/bin/test_isomerization : $(LIBDIR)/libomni.so \
				    $(TESTDIR)/Structure/test_isomerization.cpp
	@echo "[OMNI]  Building test_isomerization..."
	$(VB)$(CC) $(CPP_FLAGS) -o $(TESTDIR)/bin/test_isomerization \
	  $(TESTDIR)/Structure/test_isomerization.cpp -L$(LIBDIR) -I$(SRCDIR) -lomni

# Target: numerical precision model
$(TESTDIR)/bin/test_numerics : $(LIBDIR)/libomni.so $(TESTDIR)/Math/test_numerics.cpp
	@echo "[OMNI]  Building test_numerics..."
	$(VB)$(CC) $(CPP_FLAGS) -o $(TESTDIR)/bin/test_numerics \
	  $(TESTDIR)/Math/test_numerics.cpp -L$(LIBDIR) -I$(SRCDIR) -lomni

# Target: testing the collated coordinate / velocity / force object
$(TESTDIR)/bin/test_synthesis : $(LIBDIR)/libomni.so \
				$(TESTDIR)/Synthesis/test_synthesis.cpp
	@echo "[OMNI]  Building test_synthesis..."
	$(VB)$(CC) $(CPP_FLAGS) -o $(TESTDIR)/bin/test_synthesis \
	  $(TESTDIR)/Synthesis/test_synthesis.cpp -L$(LIBDIR) -I$(SRCDIR) -lomni

# Target: prmtop collating object test program
$(TESTDIR)/bin/test_atomgraph_synthesis : $(LIBDIR)/libomni.so \
					  $(TESTDIR)/Synthesis/test_atomgraph_synthesis.cpp
	@echo "[OMNI]  Building test_atomgraph_synthesis..."
	$(VB)$(CC) $(CPP_FLAGS) -o $(TESTDIR)/bin/test_atomgraph_synthesis \
	  $(TESTDIR)/Synthesis/test_atomgraph_synthesis.cpp -L$(LIBDIR) -I$(SRCDIR) -lomni

# Target: valence term evaluation by the most basic CPU routines
$(TESTDIR)/bin/test_valence_evaluation : $(LIBDIR)/libomni.so \
				         $(TESTDIR)/Potential/test_valence_evaluation.cpp
	@echo "[OMNI]  Building test_valence_evaluation..."
	$(VB)$(CC) $(CPP_FLAGS) -o $(TESTDIR)/bin/test_valence_evaluation \
	  $(TESTDIR)/Potential/test_valence_evaluation.cpp -L$(LIBDIR) -I$(SRCDIR) -lomni

# Target: non-bonded interaction evaluation by the most basic CPU routines
$(TESTDIR)/bin/test_nonbonded_evaluation : $(LIBDIR)/libomni.so \
					   $(TESTDIR)/Potential/test_nonbonded_evaluation.cpp
	@echo "[OMNI]  Building test_nonbonded_evaluation..."
	$(VB)$(CC) $(CPP_FLAGS) -o $(TESTDIR)/bin/test_nonbonded_evaluation \
	  $(TESTDIR)/Potential/test_nonbonded_evaluation.cpp -L$(LIBDIR) -I$(SRCDIR) -lomni

# Target: Generalized Born implicit solvent evaluation by the most basic CPU routines
$(TESTDIR)/bin/test_generalized_born : $(LIBDIR)/libomni.so \
				       $(TESTDIR)/Potential/test_generalized_born.cpp
	@echo "[OMNI]  Building test_generalized_born..."
	$(VB)$(CC) $(CPP_FLAGS) -o $(TESTDIR)/bin/test_generalized_born \
	  $(TESTDIR)/Potential/test_generalized_born.cpp -L$(LIBDIR) -I$(SRCDIR) -lomni

# Target: non-bonded neighbor list construction and usage by CPU routines
$(TESTDIR)/bin/test_neighbor_list : $(LIBDIR)/libomni.so \
				    $(TESTDIR)/Potential/test_neighbor_list.cpp
	@echo "[OMNI]  Building test_neighbor_list..."
	$(VB)$(CC) $(CPP_FLAGS) -o $(TESTDIR)/bin/test_neighbor_list \
	  $(TESTDIR)/Potential/test_neighbor_list.cpp -L$(LIBDIR) -I$(SRCDIR) -lomni

# Target: restraint construction and execution by CPU routines
$(TESTDIR)/bin/test_restraints : $(LIBDIR)/libomni.so \
				 $(TESTDIR)/Restraints/test_restraints.cpp
	@echo "[OMNI]  Building test_restraints..."
	$(VB)$(CC) $(CPP_FLAGS) -o $(TESTDIR)/bin/test_restraints \
	  $(TESTDIR)/Restraints/test_restraints.cpp -L$(LIBDIR) -I$(SRCDIR) -lomni

# Target: molecular mechanics minimization and execution by CPU routines
$(TESTDIR)/bin/test_minimization : $(LIBDIR)/libomni.so \
				   $(TESTDIR)/MolecularMechanics/test_minimization.cpp
	@echo "[OMNI]  Building test_minimization..."
	$(VB)$(CC) $(CPP_FLAGS) -o $(TESTDIR)/bin/test_minimization \
	  $(TESTDIR)/MolecularMechanics/test_minimization.cpp -L$(LIBDIR) -I$(SRCDIR) -lomni

# Target: HPC detection
$(TESTDIR)/bin/test_hpc_status : $(LIBDIR)/libomni_cuda.so \
				 $(TESTDIR)/Accelerator/test_hpc_status.cpp
	@echo "[OMNI]  Building test_hpc_status..."
	$(VB)$(CUCC) $(CUDA_FLAGS) $(CUDA_DEFINES) $(CUDA_ARCHS) \
	  -o $(TESTDIR)/bin/test_hpc_status $(TESTDIR)/Accelerator/test_hpc_status.cpp \
	  -L$(LIBDIR) -I$(SRCDIR) $(CUDA_LINKS) -lomni_cuda

# Target: HPC memory functionality
$(TESTDIR)/bin/test_hpc_hybrid : $(LIBDIR)/libomni_cuda.so \
				 $(TESTDIR)/Accelerator/test_hpc_hybrid.cpp
	@echo "[OMNI]  Building test_hpc_hybrid..."
	$(VB)$(CUCC) $(CUDA_FLAGS) $(CUDA_DEFINES) $(CUDA_ARCHS) \
	  -o $(TESTDIR)/bin/test_hpc_hybrid $(TESTDIR)/Accelerator/test_hpc_hybrid.cpp \
	  -L$(LIBDIR) -I$(SRCDIR) $(CUDA_LINKS) -lomni_cuda

# Target: HPC basic math kernels
$(TESTDIR)/bin/test_hpc_math : $(LIBDIR)/libomni_cuda.so $(TESTDIR)/Math/test_hpc_math.cu
	@echo "[OMNI]  Building test_hpc_math..."
	$(VB)$(CUCC) $(CUDA_FLAGS) $(CUDA_DEFINES) $(CUDA_ARCHS) \
	  -o $(TESTDIR)/bin/test_hpc_math $(TESTDIR)/Math/test_hpc_math.cu -L$(LIBDIR) \
	  -I$(SRCDIR) $(CUDA_LINKS) -lomni_cuda

# Target: HPC molecular system synthesis and associated operations
$(TESTDIR)/bin/test_hpc_synthesis : $(LIBDIR)/libomni_cuda.so \
				    $(TESTDIR)/Synthesis/test_hpc_synthesis.cu
	@echo "[OMNI]  Building test_hpc_synthesis..."
	$(VB)$(CUCC) $(CUDA_FLAGS) $(CUDA_DEFINES) $(CUDA_ARCHS) \
	  -o $(TESTDIR)/bin/test_hpc_synthesis $(TESTDIR)/Synthesis/test_hpc_synthesis.cu \
	  -L$(LIBDIR) -I$(SRCDIR) $(CUDA_LINKS) -lomni_cuda

# Target: Benchmarking single, double, and fixed-precision computations of valence forces
$(BENCHDIR)/bin/valence : $(LIBDIR)/libomni.so \
			  $(BENCHDIR)/ForceAccumulation/valence.cpp
	@echo "[OMNI]  Building valence benchmark..."
	$(VB)$(CC) $(CPP_FLAGS) -o $(BENCHDIR)/bin/valence \
	  $(BENCHDIR)/ForceAccumulation/valence.cpp -L$(LIBDIR) -I$(SRCDIR) -lomni

# Target: Benchmarking split accumulation of fixed-precision values for speed and resource usage
$(BENCHDIR)/bin/accumulate : $(LIBDIR)/libomni_cuda.so \
			     $(BENCHDIR)/ForceAccumulation/accumulate.cu
	@echo "[OMNI]  Building accumulation benchmark..."
	$(VB)$(CUCC) $(CUDA_FLAGS) $(CUDA_DEFINES) $(CUDA_ARCHS) -o $(BENCHDIR)/bin/accumulate \
	  $(BENCHDIR)/ForceAccumulation/accumulate.cu -L$(LIBDIR) -I$(SRCDIR) $(CUDA_LINKS) \
	  -lomni_cuda

# Target: Conformer generation
$(APPDIR)/bin/conformer.omni : $(LIBDIR)/libomni.so $(APPDIR)/Conf/src/conformer.cpp \
			       $(APPDIR)/Conf/src/setup.cpp
	@echo "[OMNI]  Building conformer.omni..."
	$(VB)$(CC) $(CPP_FLAGS) -o $(APPDIR)/bin/conformer.omni \
	  $(APPDIR)/Conf/src/conformer.cpp $(APPDIR)/Conf/src/setup.cpp \
	  -L$(LIBDIR) -I$(SRCDIR) -lomni

# Target: Molecular dynamics
$(APPDIR)/bin/dynamics.omni : $(LIBDIR)/libomni.so $(APPDIR)/Dyna/src/dynamics.cpp
	@echo "[OMNI]  Building dynamics.omni..."
	$(VB)$(CC) $(CPP_FLAGS) -o $(APPDIR)/bin/dynamics.omni \
	  $(APPDIR)/Dyna/src/dynamics.cpp -L$(LIBDIR) -I$(SRCDIR) -lomni

# Target: Force field refinement
$(APPDIR)/bin/ffrefine.omni : $(LIBDIR)/libomni.so $(APPDIR)/Ffrn/src/ffrefine.cpp
	@echo "[OMNI]  Building ffrefine.omni..."
	$(VB)$(CC) $(CPP_FLAGS) -o $(APPDIR)/bin/ffrefine.omni \
	  $(APPDIR)/Ffrn/src/ffrefine.cpp -L$(LIBDIR) -I$(SRCDIR) -lomni

install : $(LIBDIR)/libomni.so

clean :
	@echo "[OMNI]  Cleaning CPU libraries"
	$(VB)if [ -e $(SRCDIR)/Accelerator/hybrid.o ] ; then /bin/rm $(SRCDIR)/*/*.o ; fi

cuda : $(LIBDIR)/libomni_cuda.so

clean.cuda:
	@echo "[OMNI]  Cleaning GPU libraries"
	$(VB)if [ -e $(SRCDIR)/Accelerator/hybrid.o ] ; then /bin/rm $(SRCDIR)/*/*.o ; fi
	$(VB)if [ -e $(SRCDIR)/Accelerator/hpc_config.o ] ; then /bin/rm $(SRCDIR)/*/*.o ; fi

clean.bench:
	@echo "[OMNI]  Cleaning benchmark programs"
	$(VB)if [ -e $(BENCHDIR)/bin/valence ] ; then /bin/rm $(BENCHDIR)/bin/* ; fi

clean.apps:
	@echo "[OMNI]  Cleaning application programs"
	$(VB)if [ -e $(APPDIR)/bin/conformer.omni ] ; then /bin/rm $(APPDIR)/bin/* ; fi

test.exe : $(OMNI_TEST_PROGS)

test.cuda.exe : $(OMNI_TEST_CUDA_PROGS)

test : $(OMNI_TEST_PROGS)
	for PROG in $(OMNI_TEST_PROGS) ; do \
		echo "[OMNI] Execute $$PROG" ; \
		$$PROG ; \
	done

test.cuda : $(OMNI_TEST_CUDA_PROGS)
	for PROG in $(OMNI_TEST_CUDA_PROGS) ; do \
		echo "[OMNI] Execute $$PROG" ; \
		$$PROG ; \
	done

bench.exe : $(OMNI_BENCH_PROGS)

bench.cuda.exe : $(OMNI_BENCH_CUDA_PROGS)

bench : $(OMNI_BENCH_PROGS)
	for PROG in $(OMNI_BENCH_PROGS) ; do \
		echo "[OMNI} Run benchmark $$PROG" ; \
		$$PROG ; \
	done

bench.cuda : $(OMNI_BENCH_CUDA_PROGS)
	for PROG in $(OMNI_BENCH_CUDA_PROGS) ; do \
		echo "[OMNI} Run benchmark $$PROG" ; \
		$$PROG ; \
	done

apps : $(OMNI_APPS)

yes:  install
