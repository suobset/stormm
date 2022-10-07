include_directories(../)

add_library(amber_ascii amber_ascii.cpp)
add_library(barostat barostat.cpp)
add_library(coordinateframe coordinateframe.cpp)
add_library(coordinate_intake coordinate_intake.cpp)
add_library(phasespace phasespace.cpp)
add_library(thermostat thermostat.cpp)
add_library(trajectory_enumerators trajectory_enumerators.cpp)
add_library(write_frame write_frame.cpp)

install(TARGETS amber_ascii barostat coordinateframe coordinate_intake phasespace thermostat
        trajectory_enumerators write_frame DESTINATION lib)
install(FILES amber_ascii.h barostat.h coordinateframe.h coordinateframe.tpp coordinate_intake.h
        coordinate_series.h coordinate_series.tpp phasespace.h phasespace.tpp thermostat.h
        thermostat.tpp trajectory_enumerators.h write_frame.h DESTINATION include)

