# External Packages

To add any exernal packages, STORMM uses the [```ExternalProject``` CMake module](https://cmake.org/cmake/help/latest/module/ExternalProject.html). The ```ExternalProject``` module is used to download, build, and install external packages.

## How ExternalProject CMake Module Works

The ExternalProject module provides the `ExternalProject_Add()` command which handles the complete workflow for external dependencies. Here's how it works:

1. **Download Step**
   - Downloads source code from specified URL (Git, SVN, etc.)
   - Can also use local source directories
   - Verifies downloads using MD5/SHA1 hashes

2. **Configure Step** 
   - Runs CMake configure command on the downloaded source
   - Allows passing custom CMake variables and flags
   - Can specify different generators/toolchains

3. **Build Step**
   - Compiles the external project using specified build tool
   - Supports parallel builds
   - Can customize build commands and targets

4. **Install Step**
   - Installs built artifacts to specified location
   - Controls which components get installed
   - Can set custom install commands

Example usage:

```cmake
ExternalProject_Add(
    my_project
    URL https://example.com/my_project.tar.gz
    CONFIGURE_COMMAND <configure_command>
    BUILD_COMMAND <build_command>
    INSTALL_COMMAND <install_command>
)
```

If an external package has a CMakeLists.txt file, then the ```CONFIGURE_COMMAND``` can be omitted, and we can use the ```BUILD_COMMAND``` to pass the ```CMAKE_ARGS``` to the CMake configure step. The ```INSTALL_COMMAND``` is used to install the built artifacts to the specified location. However, if the external package does not have a CMakeLists.txt file, then the ```CONFIGURE_COMMAND``` must be specified.

If the external package has no CMakeLists.txt or a makefile, then the add_custom_command can be used to build the external package, and add_library can be used to link the external package to the STORMM library.

In the case of STORMM, there are two aspects to incorporating an External Package:

1. The main CMakeLists.txt file in the STORMM repository uses the ExternalProject_Add command to download, build, and install the external package. This main file also uses the add_library command to link the external package to the STORMM library.

2. Each directory in STORMM has it's own CMakeLists.txt file. This file is used to compile and build the code in that directory. A conditional loop is used to link the external package to all the required files in that directory.

## Detailed PocketFFT Integration Explanation

The C version of PocketFFT was integrated into STORMM as an example of handling a package without CMake or makefile support.

### Main CMakeLists.txt File (Case 1: No CMake or makefile support)

#### 1. Initial Setup and Project Definition

```cmake
include(ExternalProject)

ExternalProject_Add(
    PocketFFTProject
    GIT_REPOSITORY https://gitlab.mpcdf.mpg.de/mtr/pocketfft.git
    GIT_TAG master
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
    LOG_DOWNLOAD ON
    SOURCE_DIR ${CMAKE_BINARY_DIR}/third_party/pocketfft
)
```

This section initializes the external project:
- Includes CMake's ExternalProject module
- Defines a new external project named "PocketFFTProject"
- Specifies the git repository to clone from
- Sets empty configure/build/install commands (we'll handle these manually)
- Logs the download process
- Sets where to store the source files

In this case, CONFIGURE_COMMAND, BUILD_COMMAND, and INSTALL_COMMAND are intentionally left empty ("") because PocketFFT is a simple C library without its own build system (no CMake files or Makefile). Instead of using these standard ExternalProject_Add commands:

1. We handle the build process manually using add_custom_command later, which gives us precise control over the compilation process
2. This allows us to:
   - Specify exact compiler flags
   - Handle platform-specific compilation differences
   - Create the static library directly
   - Control the entire build pipeline step-by-step

If we had used the standard commands, they would have tried to run non-existent configure/build scripts. By leaving them empty and using custom commands instead, we can implement the exact build process we need for this simple C library.

### Slight Modification: if PocketFFT had its own CMakeLists.txt file

If PocketFFT had its own CMakeLists.txt file, the ExternalProject_Add command would look like this:

```cmake
ExternalProject_Add(
    PocketFFTProject
    GIT_REPOSITORY https://gitlab.mpcdf.mpg.de/mtr/pocketfft.git
    GIT_TAG master
    CONFIGURE_COMMAND /usr/bin/cmake 
      -DCMAKE_INSTALL_PREFIX=/build/third_party/pocketfft_install
      -DCMAKE_BUILD_TYPE=Release
      -DCMAKE_C_COMPILER=/usr/bin/gcc
      -DBUILD_SHARED_LIBS=OFF
      ../PocketFFTProject
    BUILD_COMMAND /usr/bin/cmake --build . --config Release
    INSTALL_COMMAND /usr/bin/cmake --install . --prefix /build/third_party/pocketfft_install
)
```

The three commands serve different purposes in the build process:

CONFIGURE_COMMAND is the first step that runs CMake's configuration phase. It:
- Sets up build variables and flags
- Checks system dependencies
- Generates build files (like Makefiles)

BUILD_COMMAND is the second step that actually compiles the code:
- Runs after configuration is complete
- Executes the generated build files
- Compiles source files into object files
- Links object files into libraries/executables
- In this example, it builds the Release configuration

INSTALL_COMMAND is the final step that copies built files to their install location:
- Copies headers to include directories
- Copies libraries to lib directories
- Copies executables to bin directories
- Sets up any needed symlinks
- In this example, it installs everything under /build/third_party/pocketfft_install

For a project where this is applicable, this would be all the configuration we need to do on the main CMakeLists.txt file.

#### 2. Library Target Setup

```cmake
add_library(PocketFFT STATIC IMPORTED GLOBAL)

ExternalProject_Get_Property(PocketFFTProject SOURCE_DIR)
set(PocketFFT_SOURCE_DIR ${SOURCE_DIR})

file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
```

This section creates the library framework:
- Creates a placeholder for our library that other targets can link against
- Retrieves the source directory path where PocketFFT will be downloaded
- Ensures the lib directory exists for our compiled library

### 3. Platform-Specific Configuration

```cmake
if(APPLE)
    execute_process(
        COMMAND xcrun --sdk macosx --show-sdk-path
        OUTPUT_VARIABLE SDKROOT
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    set(PLATFORM_FLAGS "-isysroot${SDKROOT}")
else()
    set(PLATFORM_FLAGS "-I/usr/include")
endif()
```

This section handles platform differences:
- Handles platform-specific differences between macOS and Linux
- On macOS: Gets the SDK path and sets appropriate system root flags
- On Linux: Sets system include paths

This section was required as CMake itself does not handle the linking of C/C++ headers into external projects. A manual step was required to link standard headers (in this case, ```math.h```) to PocketFFT. 

### 4. Custom Build Command
```cmake
add_custom_command(
    OUTPUT ${CMAKE_BINARY_DIR}/lib/libPocketFFT.a
    # Compilation step
    COMMAND ${CMAKE_C_COMPILER}
            -std=c99
            -O2
            ${PLATFORM_FLAGS}
            -c ${PocketFFT_SOURCE_DIR}/pocketfft.c
            -o ${CMAKE_BINARY_DIR}/pocketfft.o
    # Library creation step
    COMMAND ${CMAKE_AR} rcs ${CMAKE_BINARY_DIR}/lib/libPocketFFT.a ${CMAKE_BINARY_DIR}/pocketfft.o
    # Linux-specific ranlib step
    COMMAND $<$<NOT:$<PLATFORM_ID:Darwin>>:${CMAKE_RANLIB}> $<$<NOT:$<PLATFORM_ID:Darwin>>:${CMAKE_BINARY_DIR}/lib/libPocketFFT.a>
    DEPENDS PocketFFTProject
    WORKING_DIRECTORY ${PocketFFT_SOURCE_DIR}
    COMMENT "Building PocketFFT library"
    VERBATIM
)
```

The build process consists of three main steps:
1. **Compilation**:
   - Uses C99 standard
   - Applies O2 optimization
   - Creates an object file
2. **Library Creation**:
   - Creates a static library from the object file
3. **Index Generation**:
   - On Linux, runs ranlib to generate an index
   - Dependencies ensure the source is downloaded first

This section runs the terminal commands to compile the PocketFFT library.

### 5. Final Setup and Installation
```cmake
add_custom_target(PocketFFTBuild ALL
    DEPENDS ${CMAKE_BINARY_DIR}/lib/libPocketFFT.a
)

set_target_properties(PocketFFT PROPERTIES
    IMPORTED_LOCATION ${CMAKE_BINARY_DIR}/lib/libPocketFFT.a
)

target_include_directories(PocketFFT
    INTERFACE
        ${PocketFFT_SOURCE_DIR}
)

target_link_libraries(PocketFFT INTERFACE m)

install(FILES ${PocketFFT_SOURCE_DIR}/pocketfft.h DESTINATION include)
install(FILES ${CMAKE_BINARY_DIR}/lib/libPocketFFT.a DESTINATION lib)
```

This final section:
- Creates a build target that ensures the library is built
- Sets up the library's location for linking
- Configures include directories for users of the library
- Links against the math library
- Sets up installation rules for the header and library files

### Directory CMakeLists.txt File

In our case, PocketFFT was integrated into the /test directory. Here's how the test directory's CMakeLists.txt handles PocketFFT integration:

For each test source file in STORMM_TEST_SOURCES, the CMakeLists.txt:

1. Creates an executable from the test source
2. Links the executable to the main STORMM library
3. If PocketFFT is enabled (STORMM_INCLUDE_POCKETFFT is ON):
   - Adds PocketFFTBuild as a dependency to ensure the library is built first
   - Links the test executable to the PocketFFT library

The relevant code from test/CMakeLists.txt:

```cmake
foreach (STORMM_TEST_SOURCE ${STORMM_TEST_SOURCES})

    get_filename_component(STORMM_TEST_NAME ${STORMM_TEST_SOURCE} NAME_WE)

    add_executable(${STORMM_TEST_NAME} ${STORMM_TEST_SOURCE})
    target_link_libraries(${STORMM_TEST_NAME} ${PROJECT_NAME})

    # Add PocketFFT dependency and linking for all test targets
    if(STORMM_INCLUDE_POCKETFFT)
        add_dependencies(${STORMM_TEST_NAME} PocketFFTBuild)
        target_link_libraries(${STORMM_TEST_NAME} PocketFFT)
    endif()
    ...
endforeach()
```

This same pattern needs to be replicated in other STORMM directories that require PocketFFT functionality. Each directory's CMakeLists.txt should:

1. Check if STORMM_INCLUDE_POCKETFFT is enabled
2. Add PocketFFTBuild as a dependency for relevant targets
3. Link PocketFFT to those targets

For example, if the src directory needs PocketFFT, its CMakeLists.txt would include similar conditional linking:

This ensures consistent integration across the entire project while maintaining the conditional nature of the PocketFFT dependency.

---

## Adapting This Approach

This example can serve as a template for integrating other libraries that don't provide CMake support. Key points to consider when adapting:

1. Modify the source location and download method in `ExternalProject_Add`
2. Adjust compilation flags based on the library's requirements
3. Update the build commands to match the library's source files
4. Modify include directories and linked libraries as needed
5. Update installation rules for the library's specific files

The approach demonstrates a robust way to integrate external dependencies while maintaining proper CMake practices and cross-platform compatibility.

