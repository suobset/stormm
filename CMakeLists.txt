cmake_minimum_required(VERSION 3.18)

if (CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
    message(FATAL_ERROR "CMake should not be called from the root directory.")
endif()

project(stormm VERSION 0.1 LANGUAGES CXX)

#--- Include any extra utilities / modules required to build ---

option(STORMM_BUILD_APPS "Build the apps" ON)
option(STORMM_BUILD_TESTS "Build the tests" ON)

option(STORMM_ENABLE_CUDA "Enable CUDA support" ON)
option(STORMM_ENABLE_TEST_COVERAGE "Enable test coverage" OFF)

#--- Define any options the user can set ---


#--- Define any global settings / flags ---

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)  # set the -fPIC flag

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -g")

if (STORMM_ENABLE_CUDA)
    message("-- CUDA support will be enabled")

    enable_language(CUDA)

    set(CMAKE_CUDA_STANDARD 17)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)

    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --compiler-options=-O0 -g")

    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -DSTORMM_USE_HPC -DSTORMM_USE_CUDA")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DSTORMM_USE_HPC -DSTORMM_USE_CUDA")

    if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
        set(CMAKE_CUDA_ARCHITECTURES 60 61 70 75 80)
    endif()

    find_package(CUDAToolkit REQUIRED)
    set(STORMM_CUDA_LIBS CUDA::curand CUDA::cublas CUDA::cusolver CUDA::cudart CUDA::nvml)

endif()

if (STORMM_ENABLE_TEST_COVERAGE)
    if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "(Apple)?[Cc]lang")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-instr-generate -fcoverage-mapping")
    elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage -fprofile-arcs -ftest-coverage")
    endif()
endif()

# --- Define any source files ---

## Using GLOB is not considered a 'best practice' as new files aren't always detected by CMakes
## incremental build system, and it may also include any extra source files in the tree that you
## may not want to have included... but is still a very widely used pattern to avoid needed to
## manually list all files associated with the library making it much easier to use and maintain.
file(GLOB_RECURSE STORMM_CPP_FILES CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/src/*.cpp")
file(GLOB_RECURSE STORMM_TPP_FILES CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/src/*.tpp")
file(GLOB_RECURSE STORMM_HPP_FILES CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/src/*.h")

file(GLOB_RECURSE STORMM_CUDA_FILES          CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/src/*.cu")
file(GLOB_RECURSE STORMM_CUDA_HEADERS        CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/src/*.cuh")
file(GLOB_RECURSE STORMM_CUDA_INCLUDED_FILES CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/src/*.cui")

# --- Define the library ---

add_library(${PROJECT_NAME} SHARED ${STORMM_CPP_FILES} ${STORMM_TPP_FILES})

target_include_directories(${PROJECT_NAME}
        PUBLIC
        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
        $<INSTALL_INTERFACE:include/${PROJECT_NAME}>)

if (STORMM_ENABLE_CUDA)
    target_sources(${PROJECT_NAME} PRIVATE ${STORMM_CUDA_FILES} ${STORMM_CUDA_HEADERS} ${STORMM_CUDA_INCLUDED_FILES})
    target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
    target_link_libraries(${PROJECT_NAME} ${STORMM_CUDA_LIBS})
endif()

# --- Create an installable target ---

include(GNUInstallDirs)   # Define 'sensible' defaults for installation directories on UNIX

install(DIRECTORY "${CMAKE_SOURCE_DIR}/src/"
        DESTINATION "include/${PROJECT_NAME}"
        FILES_MATCHING PATTERN "*.h")

install(TARGETS ${PROJECT_NAME}
        EXPORT  ${PROJECT_NAME}Targets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}/${PROJECT_NAME}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}/${PROJECT_NAME}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}/${PROJECT_NAME}
        BUNDLE  DESTINATION ${CMAKE_INSTALL_BINDIR}/${PROJECT_NAME})

#--- Build the tests ---

if (STORMM_BUILD_TESTS)
    include(CTest)
    add_subdirectory(test)
endif()

#--- Build the apps ---

if (STORMM_BUILD_APPS)
    add_subdirectory(apps)
endif()
