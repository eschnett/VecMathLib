// -*-C++-*-

#ifndef MATHFUNCS_ASIN_H
#define MATHFUNCS_ASIN_H

#include "mathfuncs_base.h"

#include <cassert>
#include <cmath>



namespace vecmathlib {
  
  template<typename realvec_t>
  realvec_t mathfuncs<realvec_t>::vml_atan(realvec_t x)
  {
    // Handle negative values
    realvec_t x0 = x;
    x = fabs(x);
    
    // Reduce range using 1/x identity
    assert(all(x >= RV(0.0)));
    boolvec_t gt_one = x > RV(1.0);
    x = ifthen(gt_one, rcp(x), x);
    
    // Reduce range again using half-angle formula; see
    // <https://en.wikipedia.org/wiki/Inverse_trigonometric_functions>.
    // This is necessary for good convergence below.
    x = x / (RV(1.0) + sqrt(RV(1.0) + x*x));
    
    // Taylor expansion; see
    // <https://en.wikipedia.org/wiki/Inverse_trigonometric_functions>.
    assert(all(x >= RV(0.0) && x <= RV(0.5)));
    int const nmax = 30;        // ???
    realvec_t y = x / (RV(1.0) + x*x);
    realvec_t x2 = x * y;
    realvec_t r = y;
    for (int n=3; n<nmax; n+=2) {
      y *= RV(R(n-1) / R(n)) * x2;
      r += y;
    }
    
    // Undo second range reduction
    r = RV(2.0) * r;
    
    // Undo range reduction
    r = ifthen(gt_one, RV(M_PI_2) - r, r);
    
    // Handle negative values
    r = copysign(r, x0);
    
    return r;
  }
  
  
  
  template<typename realvec_t>
  realvec_t mathfuncs<realvec_t>::vml_acos(realvec_t x)
  {
    // Handle negative values
    boolvec_t is_negative = signbit(x);
    x = fabs(x);
    
    realvec_t r = RV(2.0) * atan(sqrt(RV(1.0) - x*x) / (RV(1.0) + x));
    
    // Handle negative values
    r = ifthen(is_negative, RV(M_PI) - r, r);
    
    return r;
  }
  
  
  
  template<typename realvec_t>
  realvec_t mathfuncs<realvec_t>::vml_asin(realvec_t x)
  {
    return RV(2.0) * atan(x / (RV(1.0) + sqrt(RV(1.0) - x*x)));
  }
  
}; // namespace vecmathlib

#endif  // #ifndef MATHFUNCS_ASIN_H
