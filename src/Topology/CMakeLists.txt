add_library(amber_prmtop_util amber_prmtop_util.cpp)
add_library(atomgraph atomgraph_constructors.cpp atomgraph_detailers.cpp atomgraph_getters.cpp
            atomgraph_miscellaneous.cpp)
add_library(atomgraph_abstracts atomgraph_abstracts.cpp)
add_library(atomgraph_intake atomgraph_intake.cpp)
add_library(atomgraph_refinement atomgraph_refinement.cpp)
add_library(topology_limits topology_limits.cpp)
add_library(topology_util topology_util.cpp)

install(TARGETS amber_prmtop_util atomgraph atomgraph_abstracts atomgraph_intake
        atomgraph_refinement topology_limits topology_util DESTINATION lib)
install(FILES amber_prmtop_util.h atomgraph_abstracts.h atomgraph_analysis.h
        atomgraph_enumerators.h atomgraph.h atomgraph.tpp atomgraph_intake.h atomgraph_refinement.h
        topology_limits.h topology_util.h DESTINATION include)
