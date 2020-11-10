# Nous voulons un cmake "récent" pour utiliser les dernières fonctionnalités
cmake_minimum_required(VERSION 3.10)

# Notre projet est étiqueté libhello
project(brisk_lib)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

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
    )

add_library(brisk_lib ${SRCS} ${HEADERS})
