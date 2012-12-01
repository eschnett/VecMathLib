// -*-C++-*-

#ifndef MATHFUNCS_POW_H
#define MATHFUNCS_POW_H

#include "mathfuncs_base.h"

#include <cassert>
#include <cmath>



namespace vecmathlib {
  
  template<typename realvec_t>
  realvec_t mathfuncs<realvec_t>::vml_pow(realvec_t x, realvec_t y)
  {
    realvec_t r = exp(log(fabs(x)) * y);
    // The result is negative if x<0 and if y is integer and odd
    realvec_t sign = fmod(x, RV(2.0)) + RV(1.0);
    return ifthen(x == RV(0.0), RV(0.0), copysign(r, sign));
  }
  
}; // namespace vecmathlib

#endif  // #ifndef MATHFUNCS_POW_H
