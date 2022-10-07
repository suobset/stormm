include_directories(../)

add_library(approx approx.cpp)
add_library(checklist checklist.cpp)
add_library(file_snapshot file_snapshot.cpp)
add_library(stopwatch stopwatch.cpp)
add_library(test_environment test_environment.cpp)
add_library(test_system_manager test_system_manager.cpp)
add_library(unit_test unit_test.cpp)
add_library(unit_test_enumerators unit_test_enumerators.cpp)
add_library(vector_report vector_report.cpp)

install(TARGETS approx checklist file_snapshot stopwatch test_environment test_system_manager
        unit_test unit_test_enumerators vector_report DESTINATION lib)
install(FILES approx.h approx.tpp checklist.h file_snapshot.h stopwatch.h test_environment.h
        test_system_manager.h unit_test_enumerators.h unit_test.h unit_test.tpp vector_report.h
        DESTINATION include)
