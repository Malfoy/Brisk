cmake_minimum_required(VERSION 3.16)

# set the project name
project(Brisk)
project(Brisk VERSION 1.0)
# configure_file(kff_io.hpp.in kff_io.hpp)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_BUILD_TYPE Debug)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Ofast -march=native -include cstdint")
set(CMAKE_CXX_FLAGS_DEBUG "-g")

# declare libraries
add_subdirectory(brisk)
link_directories(brisk)
add_subdirectory(apps)

set(SRCS
    counter.cpp
    )
    
set(HEADERS
    strict_fstream.hpp
    CLI11.hpp
    zstr.hpp
    )
# add the executable
add_executable(counter ${SRCS} ${HEADERS})

# link libraries
find_package(ZLIB)
find_package(OpenMP)
target_link_libraries(counter ZLIB::ZLIB OpenMP::OpenMP_CXX brisk kff unordered_dense)
