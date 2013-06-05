// -*-C++-*-

#ifndef FLOATTYPESS_H
#define FLOATTYPESS_H

//#include <cstdint>

#include <stdint.h>
namespace std {
  typedef unsigned char uint8_t;
  typedef signed char int8_t;
  typedef unsigned short uint16_t;
  typedef short int16_t;
  typedef unsigned int uint32_t;
  typedef int int32_t;
  typedef unsigned long long uint64_t;
  typedef long long int64_t;
}

#include <cstdlib>
#define static_assert(x,y)
#define __builtin_unreachable() (assert(0))

#include <math.h>
namespace std {
  template<typename T> int my_signbit(T x) { return signbit(x); }
#undef signbit
  
  float asinh(float x) { return ::asinh(x); }
  float acosh(float x) { return ::acosh(x); }
  float atanh(float x) { return ::atanh(x); }
  float cbrt(float x) { return ::cbrtf(x); }
  // float ceil(float x) { return ::ceilf(x); }
  float copysign(float x, float y) { return ::copysignf(x, y); }
  float exp2(float x) { return ::exp2f(x); }
  float expm1(float x) { return ::expm1f(x); }
  float fdim(float x, float y) { return ::fdimf(x, y); }
  // float floor(float x) { return ::floorf(x); }
  float fma(float x, float y, float z) { return ::fmaf(x, y, z); }
  float fmax(float x, float y) { return ::fmaxf(x, y); }
  float fmin(float x, float y) { return ::fminf(x, y); }
  float hypot(float x, float y) { return ::hypotf(x, y); }
  int ilogb(float x) { return ::ilogbf(x); }
  float log1p(float x) { return ::log1pf(x); }
  float log2(float x) { return ::log2f(x); }
  float remainder(float x, float y) { return ::remainderf(x, y); }
  float rint(float x) { return ::rintf(x); }
  float round(float x) { return ::roundf(x); }
  bool signbit(float x) { return my_signbit(x); }
  float trunc(float x) { return ::truncf(x); }
  
  double asinh(double x) { return ::asinh(x); }
  double acosh(double x) { return ::acosh(x); }
  double atanh(double x) { return ::atanh(x); }
  double cbrt(double x) { return ::cbrt(x); }
  // double ceil(double x) { return ::ceil(x); }
  double copysign(double x, double y) { return ::copysign(x, y); }
  double exp2(double x) { return ::exp2(x); }
  double expm1(double x) { return ::expm1(x); }
  double fdim(double x, double y) { return ::fdim(x, y); }
  // double floor(double x) { return ::floor(x); }
  double fma(double x, double y, double z) { return ::fma(x, y, z); }
  double fmax(double x, double y) { return ::fmax(x, y); }
  double fmin(double x, double y) { return ::fmin(x, y); }
  double hypot(double x, double y) { return ::hypot(x, y); }
  int ilogb(double x) { return ::ilogb(x); }
  double log1p(double x) { return ::log1p(x); }
  double log2(double x) { return ::log2(x); }
  double remainder(double x, double y) { return ::remainder(x, y); }
  double rint(double x) { return ::rint(x); }
  double round(double x) { return ::round(x); }
  bool signbit(double x) { return my_signbit(x); }
  double trunc(double x) { return ::trunc(x); }
}



namespace vecmathlib {
  
  struct fp8 {
    // 1 bit sign, 4 bits exponent, 3 bits mantissa
    std::uint8_t val;
  };
  
  struct fp16 {
    // 1 bit sign, 5 bits exponent, 10 bits mantissa
    std::uint16_t val;
  };
  
} // namespace vecmathlib

#endif  // #ifndef FLOATTYPES_H
