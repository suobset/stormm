include_directories(../)

add_library(mdl_mol_format mdl_mol_format.cpp)
add_library(molecule_file_io molecule_file_io.cpp)

install(TARGETS mdl_mol_format molecule_file_io DESTINATION lib)
install(FILES mdl_mol_format.h molecule_file_io.h DESTINATION include)
