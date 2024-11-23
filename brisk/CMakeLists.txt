cmake_minimum_required(VERSION 3.16)

project(brisk_lib)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_CXX_FLAGS "-Wall")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -include cstdint")
set(CMAKE_CXX_FLAGS_DEBUG "-g -Ofast")
set(CMAKE_CXX_FLAGS_RELEASE "-Ofast")

# Source files
set(SRCS
    hashing.cpp
    Kmers.cpp
    Decycling.cpp
)

# Header files
set(HEADERS
    Brisk.hpp
    hashing.hpp
    buckets.hpp
    common.hpp
    DenseMenuYo.hpp
    Kmers.hpp
    parameters.hpp
    pow2.hpp
    SuperKmerLight.hpp
    writer.hpp
    Decycling.h
)

include_directories(${CMAKE_CURRENT_BINARY_DIR})
link_directories(${CMAKE_CURRENT_BINARY_DIR})

# Add unordered_dense as a subdirectory
add_subdirectory(lib/kff)
add_subdirectory(lib/unordered_dense)

# Create brisk library
add_library(brisk ${SRCS} ${HEADERS} )

# Include directories for brisk
target_include_directories(brisk
    PUBLIC
    ./
    ${CMAKE_CURRENT_SOURCE_DIR}/lib/unordered_dense/include
)

# Link libraries for brisk
target_link_libraries(brisk
    kff
    unordered_dense
)
