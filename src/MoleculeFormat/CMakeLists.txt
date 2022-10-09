include_directories(../)

add_library(mdl_mol_format mdl_mol_format.cpp)
add_library(molecule_file_io molecule_file_io.cpp)
add_library(molecule_format_enumerators molecule_format_enumerators.cpp)

install(TARGETS mdl_mol_format molecule_file_io molecule_format_enumerators DESTINATION lib)
install(FILES mdl_mol_format.h molecule_file_io.h molecule_format_enumerators DESTINATION include)
