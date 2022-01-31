# Minimum version required
cmake_minimum_required(VERSION 3.1)

# Set the project name
project(OMNI VERSION 0.1)

# Definitions
add_definitions(-std=c++11)

# Specify the C++ 11 standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Subdirectories
add_subdirectory(Chemistry)
add_subdirectory(Constants)
add_subdirectory(Cuda)
add_subdirectory(DataTypes)
add_subdirectory(FileManagement)
add_subdirectory(Math)
add_subdirectory(MoleculeFormat)
add_subdirectory(Parsing)
add_subdirectory(Random)
add_subdirectory(Reporting)
add_subdirectory(Topology)
add_subdirectory(Trajectory)
add_subdirectory(UnitTesting)
