// -*-C++-*-

#ifndef VECMATHLIB_H
#define VECMATHLIB_H

#if defined VML_DEBUG || defined VML_NODEBUG
#  if defined VML_DEBUG && defined VML_NODEBUG
#    error "Only one of VML_DEBUG or VML_NODEBUG may be defined"
#  endif
#else
// default
#  define VML_DEBUG
#endif

#undef VML_HAVE_DENORMALS
#undef VML_HAVE_INF
#undef VML_HAVE_NAN
#define VML_HAVE_SIGNED_ZERO



// This workaround is needed for older libstdc++ versions such as the
// one in Debian 6.0 when compiled with clang++
// http://lists.cs.uiuc.edu/pipermail/cfe-dev/2011-February/013207.html
// The version time stamp used below is the one in Debian 6.0.
#include <cstring>              // pull in __GLIBCXX__
#if defined __GLIBCXX__ && __GLIBCXX__ <= 20101114
namespace std { class type_info; }
#endif



#include <cassert>



#ifdef VML_DEBUG
#  define VML_ASSERT(x) assert(x)
#else
#  define VML_ASSERT(x) ((void)0)
#endif

// Scalarise all vector operations, and use libm's functions (mostly
// useful as fallback)
#include "vec_pseudo.h"

// Use compiler-provided vector types
// #include "vec_builtin.h"

// Scalarise all vector operations; don't use libm, use only
// Vecmathlib's functions (mostly useful for testing Vecmathlib)
#include "vec_test.h"

// Intel SSE 2
#if defined __SSE2__
#  include "vec_float_sse2_scalar.h"
#  include "vec_double_sse2_scalar.h"
#  include "vec_float_sse2.h"
#  include "vec_double_sse2.h"
#endif

// Intel AVX
#if defined __AVX__
#  include "vec_fp8_avx.h"
#  include "vec_fp16_avx.h"
#  include "vec_float_avx.h"
#  include "vec_double_avx.h"
#endif

// IBM Altivec
#if defined __ALTIVEC__
#  include "vec_float_altivec.h"
#endif

#endif // #ifndef VECMATHLIB_H
