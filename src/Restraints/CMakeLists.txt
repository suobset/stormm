add_library(bounded_restraint bounded_restraint.cpp)
add_library(restraint_apparatus restraint_apparatus.cpp)
add_library(restraint_builder restraint_builder.cpp)
add_library(restraint_enumerators restraint_enumerators.cpp)

install(TARGETS bounded_restraint restraint_apparatus restraint_builder restraint_enumerators
        DESTINATION lib)
install(FILES bounded_restraint.h restraint_apparatus.h restraint_apparatus.tpp restraint_builder.h
        restraint_enumerators.h restraint_util.h restraint_util.tpp)
