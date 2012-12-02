// -*-C++-*-

#ifndef VECMATHLIB_H
#define VECMATHLIB_H

#include "vec_float.h"
#include "vec_double.h"
#if defined __AVX__             // Intel AVX
#  include "vec_double_avx.h"
#endif

#endif // #ifndef VECMATHLIB_H
