// -*-C++-*-

#ifndef VECMATHLIB_H
#define VECMATHLIB_H

#define VML_DEBUG

#undef VML_HAVE_DENORMALS
#undef VML_HAVE_INF
#undef VML_HAVE_NAN
#undef VML_HAVE_SIGNED_ZERO



#include <cassert>

#ifdef VML_DEBUG
#  define VML_ASSERT(x) assert(x)
#else
#  define VML_ASSERT(x) ((void)0)
#endif

#include "vec_pseudo.h"

#include "vec_float.h"
#if defined __SSE2__            // Intel SSE 2
#  include "vec_float_sse2.h"
#endif
#if defined __AVX__             // Intel AVX
#  include "vec_float_avx.h"
#endif

#include "vec_double.h"
#if defined __SSE2__            // Intel SSE 2
#  include "vec_double_sse2.h"
#endif
#if defined __AVX__             // Intel AVX
#  include "vec_double_avx.h"
#endif

#endif // #ifndef VECMATHLIB_H
