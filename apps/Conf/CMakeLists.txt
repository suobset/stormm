set(APP_NAME "conformer.${STORMM_APP_SUFFIX}")
set(CUDA_APP_NAME "conformer.${STORMM_APP_SUFFIX}.cuda")

add_executable(${APP_NAME}
        ${CMAKE_CURRENT_SOURCE_DIR}/src/conformer.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/setup.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/setup.h)

add_executable(${CUDA_APP_NAME}
        ${CMAKE_CURRENT_SOURCE_DIR}/src/conformer.cu
        ${CMAKE_CURRENT_SOURCE_DIR}/src/setup.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/setup.h)

target_link_libraries(${APP_NAME} ${PROJECT_NAME})

target_link_libraries(${CUDA_APP_NAME} ${PROJECT_NAME})

install(TARGETS ${APP_NAME}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

install(TARGETS ${CUDA_APP_NAME}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
