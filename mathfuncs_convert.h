// -*-C++-*-

#ifndef MATHFUNCS_CONVERT_H
#define MATHFUNCS_CONVERT_H

#include "mathfuncs_base.h"

#include <cmath>



namespace vecmathlib {
  
  template<typename realvec_t>
  realvec_t mathfuncs<realvec_t>::vml_convert_float(intvec_t x)
  {
    // Convert in two passes. Convert as much as possible during the
    // first pass (lobits), so that the second pass (hibits) may be
    // omitted if the high bits are known to be zero.
    int_t lobits = FP::mantissa_bits;
    // int_t hibits = FP::bits - lobits;
    
    // Convert lower bits
    intvec_t xlo = x & IV((U(1) << lobits) - 1);
    // exponent for the equivalent floating point number
    int_t exponent_lo = (FP::exponent_offset + lobits) << FP::mantissa_bits;
    xlo |= exponent_lo;
    // subtract hidden mantissa bit
    realvec_t flo = as_float(xlo) - RV(FP::as_float(exponent_lo));
    
    // Convert upper bits
    // make unsigned by subtracting largest negative number
    // (only do this for the high bits, since they have sufficient
    // precision to handle the overflow)
    x ^= FP::signbit_mask;
    intvec_t xhi = lsr(x, lobits);
    // exponent for the equivalent floating point number
    int_t exponent_hi = (FP::exponent_offset + 2*lobits) << FP::mantissa_bits;
    xhi |= exponent_hi;
    // subtract hidden mantissa bit
    realvec_t fhi = as_float(xhi) - RV(FP::as_float(exponent_hi));
    // add largest negative number again
    fhi -= RV(R(FP::signbit_mask));
    // Ensure that the converted low and high bits are calculated
    // separately, since a real_t doesn't have enough precision to
    // hold all the bits of an int_t
    fhi.barrier();
    
    // Combine results
    return flo + fhi;
  }
  
  
  
  template<typename realvec_t>
  auto mathfuncs<realvec_t>::vml_convert_int(realvec_t x) -> intvec_t
  {
    // Handle zero
    boolvec_t is_zero = x == RV(0.0);
    // Handle overflow
    int_t min_int = FP::signbit_mask;
    int_t max_int = ~FP::signbit_mask;
    boolvec_t is_overflow = x < RV(R(min_int)) || x > RV(R(max_int));
    // Handle negative numbers
    boolvec_t is_negative = signbit(x);
    x = fabs(x);
    
    // Round, by adding a large number that removes the excess
    // precision
    int_t large = U(1) << FP::mantissa_bits;
    x += R(large);
    
    intvec_t exponent = ilogb(x);
    for (int i=0; i<intvec_t::size; ++i) {
      VML_ASSERT(exponent[i] >= FP::mantissa_bits);
    }
    intvec_t ix = as_int(x) & IV(FP::mantissa_mask);
    // add hidden mantissa bit
    ix |= U(1) << FP::mantissa_bits;
    // shift according to exponent
    ix <<= exponent - IV(FP::mantissa_bits);
    
    // Undo the adding above
    ix -= large;
    
    // Handle negative numbers
    ix = ifthen(is_negative, -ix, ix);
    // Handle overflow
    ix = ifthen(is_overflow, IV(min_int), ix);
    // Handle zero
    ix = ifthen(is_zero, IV(I(0)), ix);
    
    return ix;
  }
  
  
  
  template<typename realvec_t>
  realvec_t mathfuncs<realvec_t>::vml_round(realvec_t x)
  {
    realvec_t r = x;
    // Round by adding a large number, destroying all excess precision
    realvec_t offset = copysign(RV(std::scalbn(R(1.0), FP::mantissa_bits)), x);
    r += offset;
    // Ensure the rounding is not optimised away
    r.barrier();
    r -= offset;
    return r;
  }
  
  template<typename realvec_t>
  realvec_t mathfuncs<realvec_t>::vml_ceil(realvec_t x)
  {
    boolvec_t iszero = x == RV(0.0);
    realvec_t offset = RV(0.5) - scalbn(fabs(x), I(-FP::mantissa_bits));
    return ifthen(iszero, x, round(x + offset));
  }
  
  template<typename realvec_t>
  realvec_t mathfuncs<realvec_t>::vml_floor(realvec_t x)
  {
    boolvec_t iszero = x == RV(0.0);
    realvec_t offset = RV(0.5) - scalbn(fabs(x), I(-FP::mantissa_bits));
    return ifthen(iszero, x, round(x - offset));
  }
  
}; // namespace vecmathlib

#endif  // #ifndef MATHFUNCS_CONVERT_H
