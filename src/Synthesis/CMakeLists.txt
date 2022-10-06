add_library(atomgraph_synthesis atomgraph_synthesis.cpp)
add_library(implicit_solvent_workspace implicit_solvent_workspace.cpp)
add_library(nonbonded_workunit nonbonded_workunit.cpp)
add_library(phasespace_synthesis phasespace_synthesis.cpp)
add_library(static_mask_synthesis static_mask_synthesis.cpp)
add_library(systemcache systemcache.cpp)
add_library(valence_workunit valence_workunit.cpp)

install(TARGETS atomgraph_synthesis implicit_solvent nonbonded_workunit phasespace_synthesis
        static_mask_synthesis systemcache valence_workunit DESTINATION lib)
install(FILES atomgraph_synthesis.h atomgraph_synthesis.tpp implicit_solvent_workspace.h
        implicit_solvent_workspace.tpp nonbonded_workunit.h nonbonded_workunit.tpp
        phasespace_synthesis.h static_mask_synthesis.h synthesis_abstracts.h
        synthesis_abstracts.tpp synthesis_enumerators.h systemcache.h valence_workunit.h
        DESTINATION include)
