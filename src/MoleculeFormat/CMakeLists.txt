include_directories(../)

add_library(mdl_mol_obj mdl_mol_obj.cpp)
add_library(molobj_atomlist molobj_atomlist.cpp)
add_library(molobj_bond molobj_bond.cpp)
add_library(molobj_dataitem molobj_dataitem.cpp)
add_library(molobj_property molobj_property.cpp)
add_library(molecule_file_io molecule_file_io.cpp)
add_library(molecule_format_enumerators molecule_format_enumerators.cpp)

install(TARGETS mdl_mol_obj molobj_atomlist molobj_bond molobj_dataitem molobj_property
        molecule_file_io molecule_format_enumerators DESTINATION lib)
install(FILES mdl_mol_obj.h molobj_atomlist.h molobj_bond.h molobj_dataitem.h molobj_property.h
        molecule_file_io.h molecule_format_enumerators DESTINATION include)
