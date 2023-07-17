cmake_minimum_required(VERSION 3.16)

project(brisk_lib)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_CXX_FLAGS "-Wall")
set(CMAKE_CXX_FLAGS_DEBUG "-g -Ofast")
set(CMAKE_CXX_FLAGS_RELEASE "-Ofast")

# Crée des variables avec les fichiers à compiler
set(SRCS
    hashing.cpp
    Kmers.cpp
    Decycling.cpp
    )
    
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

add_subdirectory(lib/kff)
target_link_libraries(kff)

add_subdirectory(lib/unordered_dense)
target_link_libraries(unordered_dense)

add_library(brisk ${SRCS} ${HEADERS})
target_include_directories(brisk
    PUBLIC
    ./
    ${CMAKE_CURRENT_SOURCE_DIR}/lib/unordered_dense/include
)
