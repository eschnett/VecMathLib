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
#undef VML_HAVE_FP_CONTRACT     // can e.g. re-associate
#undef VML_HAVE_INF
#undef VML_HAVE_NAN
#define VML_HAVE_SIGNED_ZERO



// This workaround is needed for older libstdc++ versions such as the
// one in Debian 6.0 when compiled with clang++
// <http://lists.cs.uiuc.edu/pipermail/cfe-dev/2011-February/013207.html>.
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

#if defined __ARM_PCS_VFP       // ARM NEON
// TODO: VFP
#  include "vec_neon_float2.h"
#  include "vec_neon_float4.h"
#endif

#if defined __SSE2__            // Intel SSE 2
#  include "vec_sse_float1.h"
#  include "vec_sse_float4.h"
#  include "vec_sse_double1.h"
#  include "vec_sse_double2.h"
#endif

#if defined __AVX__             // Intel AVX
#  include "vec_avx_fp8_32.h"
#  include "vec_avx_fp16_16.h"
#  include "vec_avx_float8.h"
#  include "vec_avx_double4.h"
#endif

// TODO: MIC

#if defined __ALTIVEC__         // IBM Altivec
#  include "vec_altivec_float4.h"
#endif
#if defined __VSX__             // IBM VSX
#  include "vec_vsx_double2.h"
#endif

// TODO: IBM Blue Gene/P DoubleHummer

#if defined __bgq__ && defined __VECTOR4DOUBLE__ // IBM Blue Gene/Q QPX
// TODO: vec_qpx_float4
#  include "vec_qpx_double4.h"
#endif

#endif // #ifndef VECMATHLIB_H
