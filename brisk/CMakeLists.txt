# Nous voulons un cmake "récent" pour utiliser les dernières fonctionnalités
cmake_minimum_required(VERSION 3.10)

# Notre projet est étiqueté libhello
project(brisk_lib)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_CXX_FLAGS "-Wall")
set(CMAKE_CXX_FLAGS_DEBUG "-g -Ofast")
set(CMAKE_CXX_FLAGS_RELEASE "-Ofast")

# Crée des variables avec les fichiers à compiler
set(SRCS
    Kmers.cpp
    )
    
set(HEADERS
    Brisk.hpp
    buckets.hpp
    common.hpp
    DenseMenuYo.hpp
    Kmers.hpp
    parameters.hpp
    pow2.hpp
    robin_hood.h
    SuperKmerLight.hpp
    writer.hpp
    )

include_directories(${CMAKE_CURRENT_BINARY_DIR})
link_directories(${CMAKE_CURRENT_BINARY_DIR})

add_subdirectory(lib/kff)
target_link_libraries(kff)

add_library(brisk ${SRCS} ${HEADERS})
