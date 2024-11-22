cmake_minimum_required(VERSION 3.16)

# Set the project name
project(Brisk VERSION 1.0)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_BUILD_TYPE Release)

# Compiler flags
set(CMAKE_CXX_FLAGS "-Ofast -march=native -include cstdint -flto")
set(CMAKE_CXX_FLAGS_DEBUG "-g")

# Add subdirectories
add_subdirectory(brisk)
add_subdirectory(apps)

# Source files for the counter executable
set(SRCS
    counter.cpp  # ONLY source files go here
)

# Add the executable
add_executable(counter ${SRCS})

# Include directories for header files
target_include_directories(counter PRIVATE brisk)

# Link libraries
find_package(ZLIB REQUIRED)
find_package(OpenMP REQUIRED)
target_link_libraries(counter ZLIB::ZLIB OpenMP::OpenMP_CXX brisk kff unordered_dense)
