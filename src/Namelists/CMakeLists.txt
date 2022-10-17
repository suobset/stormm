include_directories(../)

add_library(input input.cpp)
add_library(namelist_element namelist_element.cpp)
add_library(namelist_emulator namelist_emulator.cpp)
add_library(nml_conformer nml_conformer.cpp)
add_library(nml_dynamics nml_dynamics.cpp)
add_library(nml_ffmorph nml_ffmorph.cpp)
add_library(nml_files nml_files.cpp)
add_library(nml_minimize nml_minimize.cpp)
add_library(nml_precision nml_precision.cpp)
add_library(nml_random nml_random.cpp)
add_library(nml_restraint nml_restraint.cpp)
add_library(nml_report nml_report.cpp)
add_library(nml_solvent nml_solvent.cpp)
add_library(user_settings user_settings.cpp)

install(TARGETS input namelist_element namelist_emulator nml_conformer nml_dynamics nml_ffmorph
        nml_files nml_minimize nml_precision nml_random nml_report nml_restraint nml_solvent
        user_settings DESTINATION lib)
install(FILES input.h namelist_element.h namelist_emulator.h nml_conformer.h nml_dynamics.h
        nml_ffmorph.h nml_files.h nml_minimize.h nml_precision.h nml_random.h nml_restraint.h
        nml_report.h nml_solvent.h user_settings.h DESTINATION include)
