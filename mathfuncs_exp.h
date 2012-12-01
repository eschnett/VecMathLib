// -*-C++-*-

#ifndef MATHFUNCS_EXP_H
#define MATHFUNCS_EXP_H

#include "mathfuncs_base.h"

#include <cassert>
#include <cmath>



namespace vecmathlib {
  
  template<typename realvec_t>
  realvec_t mathfuncs<realvec_t>::vml_exp2(realvec_t x)
  {
    // Rescale
    realvec_t x0 = x;
    realvec_t floor_x = floor(x);
    x -= floor_x;
    intvec_t ifloor_x = convert_int(floor_x);
    
    // Approximate
    assert(all(x >= RV(0.0) && x < RV(1.0)));
    // exp(x) = Sum[n=0,nmax] x^n / n!
    int const nmax = 15;
    x -= RV(0.5);               // shift range to increase accuracy
    x *= RV(M_LN2);
    realvec_t y = RV(M_SQRT2); // x^n / n!, compensated for range shift
    realvec_t r = y;
    for (int n=1; n<nmax; ++n) {
      y *= x * RV(R(1.0) / R(n));
      r += y;
    }
    
    // Undo rescaling
    r = ifthen(x0 < RV(R(FP::min_exponent)), RV(0.0), scalbn(r, ifloor_x));
    
    return r;
  }
  
  
  
  template<typename realvec_t>
  inline
  realvec_t mathfuncs<realvec_t>::vml_exp(realvec_t x)
  {
    return exp2(RV(M_LOG2E) * x);
  }

  template<typename realvec_t>
  inline
  realvec_t mathfuncs<realvec_t>::vml_exp10(realvec_t x)
  {
    return exp(RV(M_LN10) * x);
  }

  template<typename realvec_t>
  inline
  realvec_t mathfuncs<realvec_t>::vml_expm1(realvec_t x)
  {
    return exp(x) - RV(1.0);
  }
  
}; // namespace vecmathlib

#endif  // #ifndef MATHFUNCS_EXP_H
