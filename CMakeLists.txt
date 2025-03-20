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


# Add subdirectories for brisk and apps
add_subdirectory(brisk)
add_subdirectory(apps)

