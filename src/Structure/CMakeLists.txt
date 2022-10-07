include_directories(../)

add_library(global_manipulation global_manipulation.cpp)
add_library(isomerization isomerization.cpp)
add_library(local_arrangement local_arrangement.cpp)
add_Library(rmsd rmsd.cpp)
add_library(structure_ops structure_ops.cpp)
add_library(structure_utils structure_utils.cpp)
add_library(virtual_site_handling virtual_site_handling.cpp)

install(TARGETS global_manipulation isomerization local_arrangement rmsd structure_ops
        structure_utils virtual_site_handling DESTINATION lib)
install(FILES global_manipulation.h global_manipulation.tpp isomerization.h local_arrangement.h
        local_arrangement.tpp rmsd.h rmsd.tpp structure_enumerators.h structure_ops.h
        structure_utils.h virtual_site_handling.h virtual_site_handling.tpp DESTINATION include)
