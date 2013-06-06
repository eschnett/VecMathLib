// -*-C++-*- Compatibility layer to help instantiante functions to
// create a library that can be called from elsewhere



// Make things go fast (and debugging difficult...)
#define VML_NODEBUG
#include "../vecmathlib.h"

// This must come after including vecmathlib.h
#ifdef VML_AVOID_LIBM
#  define realpseudovec realtestvec
#endif

#include <algorithm>
#include <cstdint>
#include <cstring>

#define pocl_static_assert(b) typedef char _static_assert[(b)?+1:-1]

// Define dummy values
#ifndef cl_khr_fp64
#  undef M_PI
#  define M_PI M_PI_F
#endif

// Define vector types

using std::int32_t;
#define int int32_t
typedef int int2  __attribute__((__ext_vector_type__( 2)));
typedef int int3  __attribute__((__ext_vector_type__( 3)));
typedef int int4  __attribute__((__ext_vector_type__( 4)));
typedef int int8  __attribute__((__ext_vector_type__( 8)));
typedef int int16 __attribute__((__ext_vector_type__(16)));

using std::uint32_t;
#define uint uint32_t
typedef uint uint2  __attribute__((__ext_vector_type__( 2)));
typedef uint uint3  __attribute__((__ext_vector_type__( 3)));
typedef uint uint4  __attribute__((__ext_vector_type__( 4)));
typedef uint uint8  __attribute__((__ext_vector_type__( 8)));
typedef uint uint16 __attribute__((__ext_vector_type__(16)));

#ifdef cles_khr_int64
using std::int64_t;
#define long int64_t
typedef long long2  __attribute__((__ext_vector_type__( 2)));
typedef long long3  __attribute__((__ext_vector_type__( 3)));
typedef long long4  __attribute__((__ext_vector_type__( 4)));
typedef long long8  __attribute__((__ext_vector_type__( 8)));
typedef long long16 __attribute__((__ext_vector_type__(16)));

using std::uint64_t;
#define ulong uint64_t
typedef ulong ulong2  __attribute__((__ext_vector_type__( 2)));
typedef ulong ulong3  __attribute__((__ext_vector_type__( 3)));
typedef ulong ulong4  __attribute__((__ext_vector_type__( 4)));
typedef ulong ulong8  __attribute__((__ext_vector_type__( 8)));
typedef ulong ulong16 __attribute__((__ext_vector_type__(16)));
#endif

typedef float float2  __attribute__((__ext_vector_type__( 2)));
typedef float float3  __attribute__((__ext_vector_type__( 3)));
typedef float float4  __attribute__((__ext_vector_type__( 4)));
typedef float float8  __attribute__((__ext_vector_type__( 8)));
typedef float float16 __attribute__((__ext_vector_type__(16)));

#ifdef cl_khr_fp64
typedef double double2  __attribute__((__ext_vector_type__( 2)));
typedef double double3  __attribute__((__ext_vector_type__( 3)));
typedef double double4  __attribute__((__ext_vector_type__( 4)));
typedef double double8  __attribute__((__ext_vector_type__( 8)));
typedef double double16 __attribute__((__ext_vector_type__(16)));
#endif



#define _cl_fma _cl_std_fma
#define _cl_fmax _cl_std_fmax
#define _cl_fmin _cl_std_fmin



// Generic conversion function
template<typename A, typename B>
static B bitcast(A a)
{
  B b;
  std::memcpy(&b, &a, std::min(sizeof a, sizeof b));
  if (sizeof b > sizeof a) {
    std::memset((char*)&b + sizeof a, 0, sizeof b - sizeof a);
  }
  return b;
}
