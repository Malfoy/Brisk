cmake_minimum_required(VERSION 3.16)

# Set the project name and version
project(Brisk VERSION 1.0)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_BUILD_TYPE Release)

# Compiler flags
set(CMAKE_CXX_FLAGS "-Ofast -march=native -include cstdint -flto")
set(CMAKE_CXX_FLAGS_DEBUG "-g")

# Path to the kff submodule
set(KFF_DIR ${CMAKE_SOURCE_DIR}/brisk/lib/kff)

# Build kff using ExternalProject_Add
include(ExternalProject)
ExternalProject_Add(kff_project
    SOURCE_DIR ${KFF_DIR}
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/kff_install
    BINARY_DIR ${CMAKE_BINARY_DIR}/kff_build
    BUILD_COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR}/kff_build --target all
    INSTALL_COMMAND ${CMAKE_COMMAND} --install ${CMAKE_BINARY_DIR}/kff_build --prefix ${CMAKE_BINARY_DIR}/kff_install
)

# Add subdirectories for brisk and apps
add_subdirectory(brisk)
# add_subdirectory(apps)

# Add dependencies for kff_project
add_custom_target(kff_dependency ALL DEPENDS kff_project)
