# See file "BUILD" for instructions

cmake_minimum_required (VERSION 2.6)
project (vecmathlib)

add_executable (example example.cc)
add_executable (loop loop.cc)
add_executable (test test.cc)
add_executable (bench bench.cc)

add_library (instantiations instantiations.cc)

add_subdirectory (pocl)



# # GCC:
# set (CMAKE_CXX_COMPILER "g++")
# set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -g -std=c++11 -march=native -Ofast")

# Clang:
# Note: This fails to link with -O4
set (CMAKE_CXX_COMPILER "clang++-mp-3.3")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -g -std=c++11 -stdlib=libc++ -march=native -O3")

# # Intel:
# set (CMAKE_CXX_COMPILER "icpc")
# set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -g -std=c++11 -fast")
