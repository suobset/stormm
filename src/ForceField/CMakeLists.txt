add_library(forcefield_element forcefield_element.cpp)
add_library(forcefield_enumerators  forcefield_enumerators.cpp)

install(TARGETS forcefield_element forcefield_enumerators DESTINATION lib)
install(FILES forcefield_element.h forcefield_enumerators.h DESTINATION include)
