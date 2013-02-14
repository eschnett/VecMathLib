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
    // Rescale input
    x *= RV(1.0/(2.0*M_PI));
    
    // Reduce range: sin(x) = sin(x + 2pi)
    x -= round(x);
    VML_ASSERT(all(x >= RV(-0.5) && x <= RV(+0.5)));
    
    // Reduce range: sin(x) = -sin(-x)
    realvec_t sign = x;
    x = fabs(x);
    VML_ASSERT(all(x >= RV(0.0) && x <= RV(0.5)));
    
    // Reduce range: sin(x) = sin(pi - x)
    x = fmin(x, RV(0.5)-x);
    VML_ASSERT(all(x >= RV(0.0) && x <= RV(0.25)));
    
    // Polynomial expansion
    realvec_t x2 = x*x;
    realvec_t r;
    switch (sizeof(real_t)) {
    case 4:
      // float, error=9.63199952895688527467513447362e-9
      r = RV(+39.6507114184650539767708004425);
      r = fma(r, x2, RV(-76.564420291497946510082570076));
      r = fma(r, x2, RV(+81.601622084783653475948476293)); 
      r = fma(r, x2, RV(-41.3416645219970165527584702233));
      r = fma(r, x2, RV(+6.2831851984651420764055506572));
      break;
    case 8:
      // double, error=5.7727879332805926019333251099475986171199401686760527068393e-19
      x*=RV(4.0);
      x2*=RV(16.0);
      r = RV(+5.86762142331525639276611151945042252894600631928465887664953969158e-12);
      r = fma(r, x2, RV(-6.6841794237726945077103730471119313610778831290978967922327437499159e-10));
      r = fma(r, x2, RV(+5.69213209691935763964803333647156625719770522961368e-8));
      r = fma(r, x2, RV(-3.598842978568308467678745880795603458651065024316994e-6));
      r = fma(r, x2, RV(+0.000160441184690020871131343239138998965356942754919962001));
      r = fma(r, x2, RV(-0.0046817541352970526059161105189653385731022565875740295038));
      r = fma(r, x2, RV(+0.0796926262461644478160657529010052899007494340748685273131));
      r = fma(r, x2, RV(-0.6459640975062461124233071963502627016003287601029162984167));
      r = fma(r, x2, RV(+1.5707963267948966169881546168049607518586729643258729850326));
      break;
    default:
      __builtin_unreachable();
    }
    r *= x;
    
    // Undo range reduction
    r = copysign(r, sign);
    
    return r;
  }

  template<typename realvec_t>
  realvec_t mathfuncs<realvec_t>::vml_cos(realvec_t x)
  {
    return vml_sin(x + RV(M_PI_2));
  }
  
  template<typename realvec_t>
  realvec_t mathfuncs<realvec_t>::vml_tan(realvec_t x)
  {
    return sin(x) / cos(x);
  }
  
}; // namespace vecmathlib

#endif  // #ifndef MATHFUNCS_SIN_H
