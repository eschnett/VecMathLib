# See file "BUILD" for instructions

cmake_minimum_required (VERSION 2.6)

project (vecmathlib)

add_executable (example example.cc)
add_executable (example_float example_float.cc)
add_executable (loop loop.cc)
add_executable (interp interp.cc)
add_executable (selftest selftest.cc)
add_executable (bench bench.cc)

add_library (instantiations instantiations.cc)
