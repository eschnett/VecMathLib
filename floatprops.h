// -*-C++-*-

#ifndef FLOATPROPS_H
#define FLOATPROPS_H

#include <cmath>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <limits>
#include <sstream>
#include <typeinfo>



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
  
  
  
} // namespace vecmathlib

#endif  // #ifndef FLOATPROPS_H
