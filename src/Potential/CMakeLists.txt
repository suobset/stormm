add_library(cacheresource cacheresource.cpp)
add_library(eval_valence_workunit eval_valence_workunit.cpp)
add_library(nonbonded_potential nonbonded_potential.cpp)
add_library(static_exclusionmask static_exclusionmask.cpp)
add_library(eval_synthesis eval_synthesis.cpp)
add_library(forward_exclusionmask forward_exclusionmask.cpp)
add_library(scorecard scorecard.cpp)
add_library(valence_potential valence_potential.cpp)

install(TARGETS cacheresource eval_valence_workunit nonbonded_potential static_exclusionmask
        eval_synthesis forward_exclusionmask scorecard valence_potential DESTINATION lib)
install(FILES cacheresource.h energy_abstracts.h energy_enumerators.h eval_synthesis.h
        eval_valence_workunit.h forward_exclusionmask.h nonbonded_potential.h scorecard.h
        static_exclusionmask.h valence_potential.h DESTINATION include)
