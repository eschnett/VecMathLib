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
    realvec_t round_x = round(x);
    x -= round_x;
    intvec_t iround_x = convert_int(round_x);
    
    // Approximate
    // exp(x) = Sum[n=0,nmax] x^n / n!
    assert(all(x >= RV(-0.5) && x <= RV(0.5)));
    
    // nmax   max_error
    //    5   4.2e-5
    //    6   2.4e-6
    //    7   1.2e-7
    //   11   2.2e-13
    //   12   6.3e-15
    //   13   1.7e-16
    int const nmax = sizeof(real_t)==4 ? 7 : 11;
    x *= RV(M_LN2);
    realvec_t y = x;            // x^n / n!
    realvec_t r = RV(1.0) + y;
    for (int n=2; n<nmax; ++n) {
      y *= x * RV(R(1.0) / R(n));
      r += y;
    }
    
    // Undo rescaling
    r = ifthen(x0 < RV(R(FP::min_exponent)), RV(0.0), scalbn(r, iround_x));
    
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
#if 0
    r = exp(x) - RV(1.0);
    return ifthen(r == RV(0.0), x, r);
#endif
  }
  
}; // namespace vecmathlib

#endif  // #ifndef MATHFUNCS_EXP_H
