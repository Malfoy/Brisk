cmake_minimum_required(VERSION 3.19)


# declare libraries

set(SRCS
    colorer.cpp
    )
    
set(HEADERS
    strict_fstream.hpp
    CLI11.hpp
    zstr.hpp
    )

# add the executable
add_executable(colorer ${SRCS} ${HEADERS})

# link libraries
find_package(ZLIB)
find_package(OpenMP)
target_link_libraries(colorer ZLIB::ZLIB OpenMP::OpenMP_CXX brisk kff)
