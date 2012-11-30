// -*-C++-*-

#ifndef MATHFUNCS_RCP_H
#define MATHFUNCS_RCP_H

#include "mathfuncs_base.h"

#include <cassert>
#include <cmath>



namespace vecmathlib {
  
  template<typename realvec_t>
  realvec_t mathfuncs<realvec_t>::vml_rcp(realvec_t x)
  {
    // Handle negative values
    realvec_t x0 = x;
    x = fabs(x);
    
    // Initial guess
    assert(all(x > RV(0.0)));
    intvec_t ilogb_x = ilogb(x);
    // For stability, choose a starting value that is below the result
    realvec_t r = scalbn(RV(0.5), -ilogb_x);
    
    // Iterate
    int const nmax = 7;
    for (int n=1; n<nmax; ++n) {
      // Step
      assert(all(x > RV(0.0)));
      // Newton method:
      // Solve   f(r) = 0   for   f(r) = x - 1/r
      //    r <- r - f(r) / f'(r)
      //    r <- 2 r - r^2 x
      r *= RV(2.0) - r * x;
    }
    
    // Handle negative values
    r = copysign(r, x0);
    
    return r;
  }
  
}; // namespace vecmathlib

#endif  // #ifndef MATHFUNCS_RCP_H
