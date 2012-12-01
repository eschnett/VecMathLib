// -*-C++-*-

#ifndef MATHFUNCS_SQRT_H
#define MATHFUNCS_SQRT_H

#include "mathfuncs_base.h"

#include <cassert>
#include <cmath>



namespace vecmathlib {
  
  template<typename realvec_t>
  realvec_t mathfuncs<realvec_t>::vml_sqrt(realvec_t x)
  {
    // Handle special case: zero
    boolvec_t is_zero = x <= RV(0.0);
    x = ifthen(is_zero, RV(1.0), x);
    
    // Initial guess
    assert(all(x > RV(0.0)));
#if 0
    intvec_t ilogb_x = ilogb(x);
    realvec_t r = scalbn(RV(M_SQRT2), ilogb_x >> 1);
    // TODO: divide by M_SQRT2 if ilogb_x % 2 == 1 ?
#else
    real_t correction =
      std::scalbn(R(FP::exponent_offset & 1 ? M_SQRT2 : 1.0),
                  FP::exponent_offset >> 1);
    realvec_t r = lsr(x.as_int(), 1).as_float() * RV(correction);
#endif
    
    // Iterate
    // nmax iterations give an accuracy of 2^nmax binary digits. 5
    // iterations suffice for double precision with its 53 digits.
    int const nmax = 5;
    for (int n=1; n<nmax; ++n) {
      // Step
      assert(all(x > RV(0.0)));
      // Newton method:
      // Solve   f(r) = 0   for   f(r) = x - r^2
      //    r <- r - f(r) / f'(r)
      //    r <- (r + x / r) / 2
      r = RV(0.5) * (r + x / r);
    }
    
    // Handle special case: zero
    r = ifthen(is_zero, RV(0.0), r);
    
    return r;
  }
  
  
  
  template<typename realvec_t>
  realvec_t mathfuncs<realvec_t>::vml_rsqrt(realvec_t x)
  {
    return rcp(sqrt(x));
  }
  
}; // namespace vecmathlib

#endif  // #ifndef MATHFUNCS_SQRT_H
