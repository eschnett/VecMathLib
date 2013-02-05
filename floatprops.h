// -*-C++-*-

#ifndef FLOATPROPS_H
#define FLOATPROPS_H

#include <cmath>
#ifndef __clang__
#  include <cstdint>
#endif
#include <cstring>
#include <iostream>
#include <limits>
#include <sstream>
#include <typeinfo>



#ifdef __clang__

typedef int int32_t;
typedef unsigned int uint32_t;
// typedef long int64_t;
// typedef unsigned long uint64_t;

namespace std {
  float acosh(float a) { return ::acosh(a); }
  float asinh(float a) { return ::asinh(a); }
  float atanh(float a) { return ::atanh(a); }
  float copysign(float a, float b) { return ::copysign(a, b); }
  float exp2(float a) { return ::exp2(a); }
  float expm1(float a) { return ::expm1(a); }
  float fdim(float a, float b) { return ::fdim(a, b); }
  float fma(float a, float b, float c) { return ::fma(a, b, c); }
  float fmax(float a, float b) { return ::fmax(a, b); }
  float fmin(float a, float b) { return ::fmin(a, b); }
  int ilogb(float a) { return ::ilogb(a); }
  bool isfinite(float a) { return true; /*return (isfinite)(a);*/ }
  bool isinf(float a) { return false; /*return (isinf)(a);*/ }
  bool isnan(float a) { return false; /*return (isnan)(a);*/ }
  bool isnormal(float a) { return true; /*return (isnormal)(a);*/ }
  float log1p(float a) { return ::log1p(a); }
  float log2(float a) { return ::log2(a); }
  long lrint(float a) { return ::lrint(a); }
  float remainder(float a, float b) { return ::remainder(a, b); }
  float round(float a) { return ::round(a); }
  float scalbn(float a, int i) { return ::scalbnf(a, i); }
  bool signbit(float a) { return a<0.0f; /*return (signbit)(a);*/ }
  
  double acosh(double a) { return ::acosh(a); }
  double asinh(double a) { return ::asinh(a); }
  double atanh(double a) { return ::atanh(a); }
  double copysign(double a, double b) { return ::copysign(a, b); }
  double exp2(double a) { return ::exp2(a); }
  double expm1(double a) { return ::expm1(a); }
  double fdim(double a, double b) { return ::fdim(a, b); }
  double fma(double a, double b, double c) { return ::fma(a, b, c); }
  double fmax(double a, double b) { return ::fmax(a, b); }
  double fmin(double a, double b) { return ::fmin(a, b); }
  int ilogb(double a) { return ::ilogb(a); }
  bool isfinite(double a) { return true; /*return (isfinite)(a);*/ }
  bool isinf(double a) { return false; /*return (isinf)(a);*/ }
  bool isnan(double a) { return false; /*return (isnan)(a);*/ }
  bool isnormal(double a) { return true; /*return (isnormal)(a);*/ }
  double log1p(double a) { return ::log1p(a); }
  double log2(double a) { return ::log2(a); }
  long lrint(double a) { return ::lrint(a); }
  double remainder(double a, double b) { return ::remainder(a, b); }
  double round(double a) { return ::round(a); }
  double scalbn(double a, int i) { return ::scalbn(a, i); }
  bool signbit(double a) { return a<0.0; /*return (signbit)(a);*/ }
  
  string to_string(int i) { ostringstream buf; buf<<i; return buf.str(); }
}

#endif



namespace vecmathlib {
  
  // A structure describing various properties of a floating point
  // type. Most properties are already described in numeric_limits, so
  // we inherit it.
  template<typename real_t>
  struct floatprops {
    // Some interesting properties are:
    //    digits
    //    epsilon
    //    min_exponent
    //    max_exponent
  };
  
  template<typename int_t>
  struct intprops {
  };
  
  
  
  // Properties of float
  template<>
  struct floatprops<float>: std::numeric_limits<float> {
    typedef float real_t;
    typedef int32_t int_t;
    typedef uint32_t uint_t;
    
    // Ensure the internal representation is what we expect
    static_assert(is_signed, "real_t is not signed");
    static_assert(radix==2, "real_t is not binary");
    
    // Ensure the sizes match
    static_assert(sizeof(real_t) == sizeof(int_t), "int_t has wrong size");
    static_assert(sizeof(real_t) == sizeof(uint_t), "uint_t has wrong size");
    
    // Number of bits in internal representation
    static int const bits = 8 * sizeof(real_t);
    static int const mantissa_bits = digits - 1;
    static int const signbit_bits = 1;
    static int const exponent_bits = bits - mantissa_bits - signbit_bits;
    static int const exponent_offset = 2 - min_exponent;
    static_assert(mantissa_bits + exponent_bits + signbit_bits == bits,
                  "error in bit counts");
    static uint_t const mantissa_mask = (uint_t(1) << mantissa_bits) - 1;
    static uint_t const exponent_mask =
      ((uint_t(1) << exponent_bits) - 1) << mantissa_bits;
    static uint_t const signbit_mask = uint_t(1) << (bits-1);
    static_assert((mantissa_mask & exponent_mask & signbit_mask) == uint_t(0),
                  "error in masks");
    static_assert((mantissa_mask | exponent_mask | signbit_mask) == ~uint_t(0),
                  "error in masks");
    
    // Re-interpret bit patterns
    static inline real_t as_float(int_t x)
    {
      // return *(real_t*)&x;
      // union { int_t i; real_t r; } ir;
      // return ir.i=x, ir.r;
      real_t res;
      std::memcpy(&res, &x, sizeof res);
      return res;
    }
    static inline int_t as_int(real_t x)
    {
      // return *(int_t*)&x;
      // union { real_t r; int_t i; } ri;
      // return ri.r=x, ri.i;
      int_t res;
      std::memcpy(&res, &x, sizeof res);
      return res;
    }
    
    // Convert values
    static inline real_t convert_float(int_t x) { return real_t(x); }
    static inline int_t convert_int(real_t x)
    {
      static_assert(sizeof std::lrint(x) >= sizeof(int_t),
                    "lrint() has wrong return type");
      long res = std::lrint(x);
      if (sizeof std::lrint(x) > sizeof(int_t)) {
        if (res < std::numeric_limits<int_t>::min() ||
            res > std::numeric_limits<int_t>::max())
        {
          return std::numeric_limits<int_t>::min();
        }
      }
      return res;
    }
  };
  
  template<>
  struct intprops<floatprops<float>::int_t> {
    typedef float real_t;
  };
  
  
  
  // Properties of double
  template<>
  struct floatprops<double>: std::numeric_limits<double> {
    typedef double real_t;
    typedef int64_t int_t;
    typedef uint64_t uint_t;
    
    // Ensure the internal representation is what we expect
    static_assert(is_signed, "real_t is not signed");
    static_assert(radix==2, "real_t is not binary");
    
    // Ensure the sizes match
    static_assert(sizeof(real_t) == sizeof(int_t), "int_t has wrong size");
    static_assert(sizeof(real_t) == sizeof(uint_t), "uint_t has wrong size");
    
    // Number of bits in internal representation
    static int const bits = 8 * sizeof(real_t);
    static int const mantissa_bits = digits - 1;
    static int const signbit_bits = 1;
    static int const exponent_bits = bits - mantissa_bits - signbit_bits;
    static int const exponent_offset = 2 - min_exponent;
    static_assert(mantissa_bits + exponent_bits + signbit_bits == bits,
                  "error in bit counts");
    static uint_t const mantissa_mask = (uint_t(1) << mantissa_bits) - 1;
    static uint_t const exponent_mask =
      ((uint_t(1) << exponent_bits) - 1) << mantissa_bits;
    static uint_t const signbit_mask = uint_t(1) << (bits-1);
    static_assert((mantissa_mask & exponent_mask & signbit_mask) == uint_t(0),
                  "error in masks");
    static_assert((mantissa_mask | exponent_mask | signbit_mask) == ~uint_t(0),
                  "error in masks");
    
    // Re-interpret bit patterns
    static inline real_t as_float(int_t x)
    {
      // return *(real_t*)&x;
      // union { int_t i; real_t r; } ir;
      // return ir.i=x, ir.r;
      real_t res;
      std::memcpy(&res, &x, sizeof res);
      return res;
    }
    static inline int_t as_int(real_t x)
    {
      // return *(int_t*)&x;
      // union { real_t r; int_t i; } ri;
      // return ri.r=x, ri.i;
      int_t res;
      std::memcpy(&res, &x, sizeof res);
      return res;
    }
    
    // Convert values
    static inline real_t convert_float(int_t x) { return real_t(x); }
    static inline int_t convert_int(real_t x)
    {
      static_assert(sizeof std::lrint(x) >= sizeof(int_t),
                    "lrint() has wrong return type");
      long res = std::lrint(x);
      if (sizeof std::lrint(x) > sizeof(int_t)) {
        if (res < std::numeric_limits<int_t>::min() ||
            res > std::numeric_limits<int_t>::max())
        {
          return std::numeric_limits<int_t>::min();
        }
      }
      return res;
    }
  };
  
  template<>
  struct intprops<floatprops<double>::int_t> {
    typedef double real_t;
  };
  
  
  
  template<typename int_t>
  inline typename intprops<int_t>::real_t as_float(int_t x)
  {
    typedef typename intprops<int_t>::real_t real_t;
    return floatprops<real_t>::as_float(x);
  }
  
  template<typename real_t>
  inline typename floatprops<real_t>::int_t as_int(real_t x)
  {
    return floatprops<real_t>::as_int(x);
  }
  
} // namespace vecmathlib

#endif  // #ifndef FLOATPROPS_H
