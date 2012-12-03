# build via ninja: cmake -G Ninja
# build via make: cmake .
# clean: rm CMakeCache.txt

cmake_minimum_required (VERSION 2.6)
project (vecmathlib)

add_executable (example example.cc)
add_executable (test test.cc)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -g -std=gnu++11 -march=native -Ofast")
