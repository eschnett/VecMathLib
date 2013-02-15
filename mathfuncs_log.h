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
    
    assert(all(x >= RV(1.0) && x < RV(2.0)));
    
    realvec_t y = (x - RV(1.0)) / (x + RV(1.0));
    realvec_t y2 = y*y;
    
    realvec_t r;
    switch (sizeof(real_t)) {
    case 4:
      // float, error=5.98355642684398209498469870525e-9
      r = RV(0.410981538282433293325329456838);
      r = fma(r, y2, RV(0.402155483172044562892705980539));
      r = fma(r, y2, RV(0.57755014627178237959721643293));
      r = fma(r, y2, RV(0.96178780600659929206930296869));
      r = fma(r, y2, RV(2.88539012786343587248965772685));
      break;
    case 8:
      // double, error=9.45037202901655672811489051683e-17
      r = RV(0.259935726478127940817401224248);
      r = fma(r, y2, RV(0.140676370079882918464564658472));
      r = fma(r, y2, RV(0.196513478841924000569879320851));
      r = fma(r, y2, RV(0.221596471338300882039273355617));
      r = fma(r, y2, RV(0.262327298560598641020007602127));
      r = fma(r, y2, RV(0.320598261015170101859472461613));
      r = fma(r, y2, RV(0.412198595799726905825871956187));
      r = fma(r, y2, RV(0.57707801621733949207376840932));
      r = fma(r, y2, RV(0.96179669392666302667713134701));
      r = fma(r, y2, RV(2.88539008177792581277410991327));
      break;
    default:
      __builtin_unreachable();
    }
    r *= y;
    
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
