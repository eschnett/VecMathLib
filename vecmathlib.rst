Vecmathlib
==========

:author: Erik Schnetter <eschnetter@perimeterinstitute.ca>

Vecmathlib provides efficient, accurate, tunable, and most importantly
vectorizable math functions such as sqrt, sin, or atan.

The library is implemented in C++, and intended to be called on SIMD
vectors, e.g. those provided by SSE, AVX, or available in Power7 and
Blue Gene architectures. The same algorithms should also work
efficiently on accelerators such as GPUs.



Implementation
--------------

Vecmathlib consists of three parts:

1. vecmathlib's math algorithms, i.e. the implementations of various
   math functions (e.g. sqrt)
2. SIMD vector classes, wrapping e.g. SSE or AVX vectors
3. a test and benchmarking harness

The algorithms themselves are written in a generic way. They assume an
IEEE floating point layout (consisting of sign bit, exponent, and
mantissa), but work for arbitrary precision and vector sizes. For
example, there is a routine vml_sqrt() that calculates a square root
via an iterative scheme based on Newton's root finding algorithm.
Although not available yet, there can be different implementations
with different performance characteristics for certain math functions.

The SIMD vector classes wrap architecture-specific SIMD capabilities;
for example, there is an implementation of a class realvec<double,4>
based on Intel's AVX instruction set. These classes either provide
math operations and math functions themselves, or implement them via
calls to the generic algorithms. This way, vecmathlib can provide
efficient implementations for all hardware architectures.

It goes without saying that vecmathlib can also be used for scalar
types, e.g. plain float or double, thus providing math functions for
architectures where they are otherwise not available.



Things To Do
------------

Vecmathlib is not finished. Contributions are welcome! There are
several areas where it can be improved:

1. make test harness more systematic, improve coverage
2. research and implement improved algorithms for certain math
   functions
3. implement vector classes for additional hardware architectures
4. make code more portable -- currently works with GCC 4.7, should
   also support pother versions or compilers
5. review C++ class hierarchy, improve design, reduce redundancy
6. measure performance, compare to system library
7. provide "vector" implementation of math functions by calculating
   them element-wise via libc's scalar functions
8. handle inf, nan, negative zero, rounding modes, and everything else
   required for IEEE compliance (if -ffast-math is not used)
