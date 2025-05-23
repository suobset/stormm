set(APP_NAME "mmgbsa.${STORMM_APP_SUFFIX}")

add_executable(${APP_NAME}
        ${CMAKE_CURRENT_SOURCE_DIR}/src/mmgbsa_analysis.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/mmgbsa_analysis.h
        ${CMAKE_CURRENT_SOURCE_DIR}/src/mmgbsa_carveout.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/mmgbsa_carveout.h
        ${CMAKE_CURRENT_SOURCE_DIR}/src/mmgbsa_problem_set.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/mmgbsa_problem_set.h
        ${CMAKE_CURRENT_SOURCE_DIR}/src/mmgbsa_testing.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/mmgbsa_testing.h
        ${CMAKE_CURRENT_SOURCE_DIR}/src/mmgbsa_enumerators.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/mmgbsa_enumerators.h
        ${CMAKE_CURRENT_SOURCE_DIR}/src/nml_mmgbsa.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/nml_mmgbsa.h)

target_sources(${APP_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src/mmgbsa.cpp)

target_link_libraries(${APP_NAME} ${PROJECT_NAME})

install(TARGETS ${APP_NAME}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
