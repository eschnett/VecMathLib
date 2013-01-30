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
    // for |x|>0.01: (*)  log(x) = Sum[n=1,nmax,n%2==1] 2/n ((x-1) / (x+1))^n
    // else:         (**) log(x) = Sum[n=1,nmax] (-1)^(n+1) 1/n (x-1)^n
    assert(all(x >= RV(1.0) && x < RV(2.0)));
    
    // nmax   max_error of (*)
    //    5   5.9e-5
    //    7   1.3e-6
    //    9   2.9e-8
    //   15   4.4e-13
    //   17   1.1e-14
    //   19   3.0e-16
    int const nmax = sizeof(real_t)==4 ? 9 : 17;
    x *= RV(M_SQRT1_2);         // shift range to increase accuracy
    
    realvec_t xm1 = x - RV(1.0);
    boolvec_t near1 = fabs(xm1) < RV(0.0001); // epsilon^(1/niters)
    
    // for (*)
    realvec_t xm1_xp1 = xm1 / (x + RV(1.0));
    realvec_t xm1_xp1_2 = xm1_xp1 * xm1_xp1;
    
    // for (**)
    realvec_t mxm1 = - xm1;
    
    realvec_t y  = ifthen(near1, xm1,  RV(2.0) * xm1_xp1);
    realvec_t yf = ifthen(near1, mxm1, xm1_xp1_2);
    y *= RV(M_LOG2E);
    
    realvec_t r = y;
    for (int n=3, nn=2; n<nmax; n+=2, ++nn) {
      y *= yf;
      r += y * ifthen(near1, RV(R(1.0) / R(nn)), RV(R(1.0) / R(n)));
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
#if 0
    // Goldberg, theorem 4
    realvec_t x1 = RV(1.0) + x;
    x1.barrier();
    return ifthen(x1 == x, x, x * log(x1) / (x1 - RV(1.0)));
#endif
  }
  
}; // namespace vecmathlib

#endif  // #ifndef MATHFUNCS_LOG_H
