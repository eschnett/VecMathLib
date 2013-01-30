// -*-C++-*-

#ifndef MATHFUNCS_SIN_H
#define MATHFUNCS_SIN_H

#include "mathfuncs_base.h"

#include <cassert>
#include <cmath>



namespace vecmathlib {
  
  template<typename realvec_t>
  realvec_t mathfuncs<realvec_t>::vml_sin(realvec_t x)
  {
    typename realvec_t::scalar_t eps __attribute__((__unused__)) = 0.000001;
    
    // Reduce range: sin(x) = sin(x + 2pi)
    x = remainder(x, RV(2.0*M_PI));
    assert(all(x >= RV(-(1.0+eps)*M_PI) && x <= RV((1.0+eps)*M_PI)));
    
    // Reduce range: sin(x) = -sin(-x)
    realvec_t sign = x;
    x = fabs(x);
    assert(all(x >= RV(0.0) && x <= RV((1.0+eps)*M_PI)));
    
    // Reduce range: sin(x) = sin(pi-x)
    x = ifthen(x > RV(M_PI_2), RV(M_PI) - x, x);
    assert(all(x >= RV(0.0-eps) && x <= RV(M_PI_2)));
    
    // Reduce range: cos(x) = sin(pi/2 - x)
    boolvec_t use_cos = x > RV(0.5*M_PI_2);
    x = ifthen(use_cos, RV(M_PI_2) - x, x);
    
    // Taylor expansion
    // sin(x) = Sum[n=0,nmax,n%2!=0] (-1)^((n-1)/2) x^n / n!
    // cos(x) = Sum[n=0,nmax,n%2==0] (-1)^((n-1)/2) x^n / n!
    realvec_t noffset = ifthen(use_cos, RV(-1.0), RV(0.0));
    int const nmax = 15;
    realvec_t x2 = x*x;
    realvec_t y = ifthen(use_cos, RV(1.0), x);
    realvec_t r = y;
    for (int n=3; n<nmax; n+=2) {
      // sin: (n-1)*n
      // cos: n*(n+1)
      realvec_t rn = RV(FP::convert_float(n)) + noffset;
      y *= - x2 * rcp((rn-RV(1.0))*rn);
      r += y;
    }
    
    // Undo range reduction
    r = copysign(r, sign);
    
    return r;
  }
  
  
  
  template<typename realvec_t>
  realvec_t mathfuncs<realvec_t>::vml_cos(realvec_t x)
  {
    return sin(x + RV(M_PI_2));
  }
  
  template<typename realvec_t>
  realvec_t mathfuncs<realvec_t>::vml_tan(realvec_t x)
  {
    // Reduce range: tan(x) = tan(x + pi)
    x = remainder(x, RV(M_PI));
    assert(all(x >= RV(-M_PI_2) && x <= RV(M_PI_2)));
    
    // Reduce range: tan(x) = -tan(-x)
    realvec_t sign = x;
    x = fabs(x);
    assert(all(x >= RV(0.0) && x <= RV(M_PI_2)));
    
    // Reduce range: cot(x) = -tan(x + pi/2)
    boolvec_t use_cot = x > RV(0.5 * M_PI_2);
    x = ifthen(use_cot, RV(M_PI_2) -x, x);
    assert(all(x >= RV(0.0) && x <= RV(0.5 * M_PI_2)));
    
    // Calculate tan
    realvec_t r = sin(x) / cos(x);
    
    // Undo range reduction
    r = ifthen(use_cot, -rcp(r), r);
    
    // Undo range reducion
    r = copysign(r, sign);
    
    return r;
  }
  
}; // namespace vecmathlib

#endif  // #ifndef MATHFUNCS_SIN_H
