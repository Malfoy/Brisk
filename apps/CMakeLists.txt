cmake_minimum_required(VERSION 3.16)


# declare libraries

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
target_link_libraries(counter ZLIB::ZLIB OpenMP::OpenMP_CXX brisk kff)
