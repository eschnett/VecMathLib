// -*-C++-*-

#ifndef MATHFUNCS_FABS_H
#define MATHFUNCS_FABS_H

#include "mathfuncs_base.h"

#include <cassert>
#include <cmath>



namespace vecmathlib {
  
  template<typename realvec_t>
  realvec_t mathfuncs<realvec_t>::vml_copysign(realvec_t x, realvec_t y)
  {
    intvec_t value = as_int(x) & IV(~FP::sign_mask);
    intvec_t sign = as_int(y) & IV(FP::sign_mask);
    return as_float(sign | value);
  }
  
  template<typename realvec_t>
  realvec_t mathfuncs<realvec_t>::vml_fabs(realvec_t x)
  {
    return as_float(as_int(x) & IV(~FP::sign_mask));
  }
  
  template<typename realvec_t>
  auto mathfuncs<realvec_t>::vml_ilogb(realvec_t x) -> intvec_t
  {
    return
      lsr(as_int(x) & IV(FP::exponent_mask), FP::mantissa_bits) -
      IV(FP::exponent_offset);
  }
  
  template<typename realvec_t>
  realvec_t mathfuncs<realvec_t>::vml_scalbn(realvec_t x, intvec_t n)
  {
    return as_float(as_int(x) + (n << FP::mantissa_bits));
    // return x * as_float((n + exponent_offset) << mantissa_bits);
  }
  
  template<typename realvec_t>
  auto mathfuncs<realvec_t>::vml_signbit(realvec_t x) -> boolvec_t
  {
    return convert_bool(as_int(x) & IV(FP::sign_mask));
  }
  
}; // namespace vecmathlib

#endif  // #ifndef MATHFUNCS_FABS_H
