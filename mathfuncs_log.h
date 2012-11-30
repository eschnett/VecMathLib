// -*-C++-*-

#ifndef MATHFUNCS_LOG_H
#define MATHFUNCS_LOG_H

#include "mathfuncs_base.h"

#include <cassert>
#include <cmath>



namespace vecmathlib {
  
  template<typename realvec_t>
  realvec_t mathfuncs<realvec_t>::vml_log2(realvec_t x)
  {
    // Rescale
    assert(all(x > RV(0.0)));
    intvec_t ilogb_x = ilogb(x);
    x = scalbn(x, -ilogb_x);
    assert(all(x >= RV(1.0) && x < RV(2.0)));
    
    // Approximate
    assert(all(x >= RV(1.0) && x < RV(2.0)));
    // log(x) = Sum[n=1,nmax,n%2==1] 2/n ((x-1) / (x+1))^n
    int const nmax = 30;
    x *= RV(M_SQRT1_2);         // shift range to increase accuracy
    realvec_t xm1_xp1 = (x - RV(1.0)) / (x + RV(1.0));
    realvec_t xm1_xp1_2 = xm1_xp1 * xm1_xp1;
    realvec_t y = RV(M_LOG2E * 2.0) * xm1_xp1;
    realvec_t r = y;
    for (int n=3; n<nmax; n+=2) {
      y *= xm1_xp1_2;
      r += y * RV(R(1.0) / R(n));
    }
    r += RV(0.5);               // correct result for range shift
    
    // Undo rescaling
    r += convert_float(ilogb_x);
    
    return r;
  }
  
  
  
  template<typename realvec_t>
  inline
  realvec_t mathfuncs<realvec_t>::vml_log(realvec_t x)
  {
    return log2(x) * RV(M_LN2);
  }

  template<typename realvec_t>
  inline
  realvec_t mathfuncs<realvec_t>::vml_log10(realvec_t x)
  {
    return log(x) * RV(M_LOG10E);
  }

  template<typename realvec_t>
  inline
  realvec_t mathfuncs<realvec_t>::vml_log1p(realvec_t x)
  {
    return log(RV(1.0) + x);
  }
  
}; // namespace vecmathlib

#endif  // #ifndef MATHFUNCS_LOG_H
