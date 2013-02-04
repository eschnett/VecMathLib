#! /bin/bash

# See file "BUILD" for instructions

rm -f CMakeCache.txt

# Unix Makefiles (build with "make")
cmake -G 'Unix Makefiles'

# Ninja (build with "ninja")
# Note: Ninja is ofteh faster than make, but may not always be available
#cmake -G Ninja
