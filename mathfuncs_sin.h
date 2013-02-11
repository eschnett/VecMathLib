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
    return vml_sin(x + RV(M_PI_2));
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

  template<typename realvec_t>
  realvec_t mathfuncs<realvec_t>::vml_sin_chebyshev_single(realvec_t x)
  {
    typename realvec_t::scalar_t const eps __attribute__((__unused__)) = 0.000001;
    realvec_t r, s, t;

    // Reduce range: sin(x) = sin(x + 2PI)
    t = x + RV(M_PI);
    t = t * RV(1.0 / (2.0 * M_PI));
    t = floor(t);
    t = x - (t * RV(2.0 * M_PI));
    assert(all(t >= RV(-(1.0 + eps) * M_PI) && t <= RV((1.0 + eps) * M_PI)));

    // Reduce range: -sin(x) = -sin(-x)
    s = copysign(RV(1.0), t);
    t = fabs(t);

    // Reduce range: sin(x) = sin(pi-x)
    x = ifthen(t > RV(M_PI_2), RV(M_PI) - t, t);
    assert(all(x >= RV(0.0 - eps) && x <= RV(M_PI_2)));

    // Evaluate Chebyshev polynomial expansion r = s*T(x)
    // where T(x) = x + c1*x^3 + c2*x^5 + c3*x^7 + ...
    t = x * x;
    x = x * s;
    s = t * x;
    r = fma(RV(-2.2636462059494963e-09), t, RV(1.9761453001926614e-08));
    r = fma(r, t, RV(-9.1621870444937360e-08));
    r = fma(r, t, RV(2.8673509297532770e-06));
    r = fma(r, t, RV(-1.9850777498822076e-04));
    r = fma(r, t, RV(8.3333708144611930e-03));
    r = fma(r, t, RV(-1.6666667163372040e-01));
    r = fma(r, s, x);
    return r;
  }

  template<typename realvec_t>
  realvec_t mathfuncs<realvec_t>::vml_sin_chebyshev_double(realvec_t x)
  {
    typename realvec_t::scalar_t const eps __attribute__((__unused__)) = 0.000001;
    realvec_t r, s, t;

    // Reduce range: sin(x) = sin(x + 2PI)
    t = x + RV(M_PI);
    t = t * RV(1.0 / (2.0 * M_PI));
    t = floor(t);
    t = x - (t * RV(2.0 * M_PI));
    assert(all(t >= RV(-(1.0 + eps) * M_PI) && t <= RV((1.0 + eps) * M_PI)));

    // Reduce range: -sin(x) = -sin(-x)
    s = copysign(RV(1.0), t);
    t = fabs(t);

    // Reduce range: sin(x) = sin(pi-x)
    x = ifthen(t > RV(M_PI_2), RV(M_PI) - t, t);
    assert(all(x >= RV(0.0 - eps) && x <= RV(M_PI_2)));

    // Evaluate Chebyshev polynomial expansion r = s*T(x)
    // where T(x) = x + c1*x^3 + c2*x^5 + c3*x^7 + ...
    t = x * x;
    x = x * s;
    s = t * x;
    r = fma(RV(1.1888871779171205e-23), t, RV(-1.6213346583948200e-21));
    r = fma(r, t, RV(9.4674830124704450e-19));
    r = fma(r, t, RV(-3.6586864533554100e-16));
    r = fma(r, t, RV(-7.5815669263036780e-13));
    r = fma(r, t, RV(1.6058175109947732e-10));
    r = fma(r, t, RV(-2.5052101017560582e-08));
    r = fma(r, t, RV(2.7557319185249400e-06));
    r = fma(r, t, RV(-1.9841269841152493e-04));
    r = fma(r, t, RV(8.3333333333331560e-03));
    r = fma(r, t, RV(-1.6666666666666666e-01));
    r = fma(r, s, x);
    return r;
  }

  template<typename realvec_t>
  realvec_t mathfuncs<realvec_t>::vml_cos_chebyshev_single(realvec_t x)
  {
    typename realvec_t::boolvec_t m;
    realvec_t r, s, t;

    // Reduce range: cos(x) = cos(x + 2PI)
    t = x + RV(M_PI);
    t = t * RV(1.0 / (2.0 * M_PI));
    t = floor(t);
    t = x - (t * RV(2.0 * M_PI));

    // Reduce range: cos(x) = sin(PI/2 - x)
    t = fabs(t);
    m = t > RV(M_PI_2);
    x = RV(M_PI) - t;
    s = ifthen(m, RV(-1.0), RV(1.0));
    x = ifthen(m, x, t);

    // Evaluate Chebyshev polynomial s*T(x) where
    // T(x) = 1 + c1*x^2 + c2*x^4 + c3*x^6 + ...
    t = x * x;
    x = t * s;
    r = fma(RV(-1.0986536764839979e-11), t, RV(2.0856702467901100e-09));
    r = fma(r, t, RV(-2.7556891974788950e-07));
    r = fma(r, t, RV(2.4801582421938170e-05));
    r = fma(r, t, RV(-1.3888888860984269e-03));
    r = fma(r, t, RV(4.1666666666056330e-02));
    r = fma(r, t, RV(-5.0000000000000000e-01));
    r = fma(r, x, s);
    return r;
  }

  template<typename realvec_t>
  realvec_t mathfuncs<realvec_t>::vml_cos_chebyshev_double(realvec_t x)
  {
    typename realvec_t::boolvec_t m;
    realvec_t r, s, t;

    // Reduce range: cos(x) = cos(x + 2PI)
    t = x + RV(M_PI);
    t = t * RV(1.0 / (2.0 * M_PI));
    t = floor(t);
    t = x - (t * RV(2.0 * M_PI));

    // Reduce range: cos(x) = sin(PI/2 - x)
    t = fabs(t);
    m = t > RV(M_PI_2);
    x = RV(M_PI) - t;
    s = ifthen(m, RV(-1.0), RV(1.0));
    x = ifthen(m, x, t);

    // Evaluate Chebyshev polynomial s*T(x) where
    // T(x) = 1 + c1*x^2 + c2*x^4 + c3*x^6 + ...
    t = x * x;
    x = t * s;
    r = fma(RV(-8.6512994843471700e-22), t, RV(4.1086675770914360e-19));
    r = fma(r, t, RV(-1.5619143199049570e-16));
    r = fma(r, t, RV(4.7794771764282040e-14));
    r = fma(r, t, RV(-1.1470745595224050e-11));
    r = fma(r, t, RV(2.0876756987841530e-09));
    r = fma(r, t, RV(-2.7557319223985710e-07));
    r = fma(r, t, RV(2.4801587301587300e-05));
    r = fma(r, t, RV(-1.3888888888888890e-03));
    r = fma(r, t, RV(4.1666666666666664e-02));
    r = fma(r, t, RV(-5.0000000000000000e-01));
    r = fma(r, x, s);
    return r;
  }

}; // namespace vecmathlib

#endif  // #ifndef MATHFUNCS_SIN_H
