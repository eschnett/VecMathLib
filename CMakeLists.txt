# See file "BUILD" for instructions

cmake_minimum_required (VERSION 2.6)
project (vecmathlib)

add_executable (example example.cc)
add_executable (loop loop.cc)
add_executable (test test.cc)
add_executable (bench bench.cc)

# GCC:
set (CMAKE_CXX_COMPILER "g++")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -g -std=gnu++11 -march=native -Ofast")

# Clang:
# Note: This fails to link with -O4
# set (CMAKE_CXX_COMPILER "clang++-mp-3.2")
# set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -g -std=gnu++11 -march=native -O3")

# Clang for pocl:
# set (CMAKE_CXX_COMPILER "/Users/eschnett/llvm-3.2/bin/clang++")
# set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -g -std=gnu++11 -march=native -O4")

# Intel:
# set (CMAKE_CXX_COMPILER "icpc")
# set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -g -std=c++11 -fast")
