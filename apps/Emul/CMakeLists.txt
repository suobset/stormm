set(APP_NAME "emulate.${STORMM_APP_SUFFIX}")

add_executable(${APP_NAME}
        ${CMAKE_CURRENT_SOURCE_DIR}/src/emulate_analysis.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/emulate_analysis.h
        ${CMAKE_CURRENT_SOURCE_DIR}/src/emulate_enumerators.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/emulate_enumerators.h)

target_sources(${APP_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src/emulate.cpp)

target_link_libraries(${APP_NAME} ${PROJECT_NAME})

install(TARGETS ${APP_NAME}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
