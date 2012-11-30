// -*-C++-*-

#ifndef VEC_BASE_H
#define VEC_BASE_H

#include <iostream>

namespace vecmathlib {
  
  template<typename real_t, int size>
  struct boolvec {
  };
  
  template<typename real_t, int size>
  struct intvec {
  };
  
  template<typename real_t, int size>
  struct realvec {
  };
  

  
  // boolvec wrappers
  
  template<typename real_t, int size>
  inline intvec<real_t, size> as_int(boolvec<real_t, size> x)
  {
    return x.as_int();
  }
  
  template<typename real_t, int size>
  inline intvec<real_t, size> convert_int(boolvec<real_t, size> x)
  {
    return x.convert_int();
  }
  
  template<typename real_t, int size>
  inline bool all(boolvec<real_t, size> x) { return x.all(); }
  
  template<typename real_t, int size>
  inline bool any(boolvec<real_t, size> x) { return x.any(); }
  
  template<typename real_t, int size>
  inline
  intvec<real_t, size> ifthen(boolvec<real_t, size> c,
                               intvec<real_t, size> x,
                               intvec<real_t, size> y)
  {
    return c.ifthen(x, y);
  }
  
  template<typename real_t, int size>
  inline
  realvec<real_t, size> ifthen(boolvec<real_t, size> c,
                                 realvec<real_t, size> x,
                                 realvec<real_t, size> y)
  {
    return c.ifthen(x, y);
  }
  
  
  
  // intvec wrappers
  
  template<typename real_t, int size>
  inline boolvec<real_t, size> as_bool(intvec<real_t, size> x)
  {
    return x.as_bool();
  }
  
  template<typename real_t, int size>
  inline boolvec<real_t, size> convert_bool(intvec<real_t, size> x)
  {
    return x.convert_bool();
  }
  
  template<typename real_t, int size>
  inline realvec<real_t, size> as_float(intvec<real_t, size> x)
  {
    return x.as_float();
  }
  
  template<typename real_t, int size>
  inline realvec<real_t, size> convert_float(intvec<real_t, size> x)
  {
    return x.convert_float();
  }
  
  template<typename real_t, int size>
  inline intvec<real_t, size> lsr(intvec<real_t, size> x,
                                   typename intvec<real_t, size>::int_t n)
  {
    return x.lsr(n);
  }
  
  template<typename real_t, int size>
  inline intvec<real_t, size> lsr(intvec<real_t, size> x,
                                   intvec<real_t, size> n)
  {
    return x.lsr(n);
  }
  
  
  
  // realvec wrappers
  
  template<typename real_t, int size>
  inline intvec<real_t, size> as_int(realvec<real_t, size> x)
  {
    return x.as_int();
  }
  
  template<typename real_t, int size>
  inline intvec<real_t, size> convert_int(realvec<real_t, size> x)
  {
    return x.convert_int();
  }
  
  template<typename real_t, int size>
  inline
  typename realvec<real_t, size>::real_t prod(realvec<real_t, size> x)
  {
    return x.prod();
  }
  
  template<typename real_t, int size>
  inline
  typename realvec<real_t, size>::real_t sum(realvec<real_t, size> x)
  {
    return x.sum();
  }
  
  
  
  template<typename real_t, int size>
  inline realvec<real_t, size> copysign(realvec<real_t, size> x,
                                          realvec<real_t, size> y)
  {
    return x.copysign(y);
  }
  
  template<typename real_t, int size>
  inline realvec<real_t, size> fabs(realvec<real_t, size> x)
  {
    return x.fabs();
  }
  
  template<typename real_t, int size>
  inline intvec<real_t, size> ilogb(realvec<real_t, size> x)
  {
    return x.ilogb();
  }
  
  template<typename real_t, int size>
  inline
  realvec<real_t, size> scalbn(realvec<real_t, size> x,
                                 intvec<real_t, size> n)
  {
    return x.scalbn(n);
  }
  
  template<typename real_t, int size>
  inline boolvec<real_t, size> signbit(realvec<real_t, size> x)
  {
    return x.signbit();
  }
  
  
  
  template<typename real_t, int size>
  inline realvec<real_t, size> acos(realvec<real_t, size> x)
  {
    return x.acos();
  }
  
  template<typename real_t, int size>
  inline realvec<real_t, size> asin(realvec<real_t, size> x)
  {
    return x.asin();
  }
  
  template<typename real_t, int size>
  inline realvec<real_t, size> atan(realvec<real_t, size> x)
  {
    return x.atan();
  }
  
  template<typename real_t, int size>
  inline realvec<real_t, size> floor(realvec<real_t, size> x)
  {
    return x.floor();
  }
  
  template<typename real_t, int size>
  inline realvec<real_t, size> log(realvec<real_t, size> x)
  {
    return x.log();
  }
  
  template<typename real_t, int size>
  inline realvec<real_t, size> log10(realvec<real_t, size> x)
  {
    return x.log10();
  }
  
  template<typename real_t, int size>
  inline realvec<real_t, size> log1p(realvec<real_t, size> x)
  {
    return x.log1p();
  }
  
  template<typename real_t, int size>
  inline realvec<real_t, size> log2(realvec<real_t, size> x)
  {
    return x.log2();
  }
  
  template<typename real_t, int size>
  inline realvec<real_t, size> rcp(realvec<real_t, size> x)
  {
    return x.rcp();
  }
  
  template<typename real_t, int size>
  inline realvec<real_t, size> rsqrt(realvec<real_t, size> x)
  {
    return x.rsqrt();
  }
  
  template<typename real_t, int size>
  inline realvec<real_t, size> sqrt(realvec<real_t, size> x)
  {
    return x.sqrt();
  }
  
  
  
  template<typename real_t, int size>
  std::ostream& operator<<(std::ostream& os, boolvec<real_t, size> const& x)
  {
    os << "[";
    for (int i=0; i<size; ++i) {
      if (i!=0) os << ",";
      os << x[i];
    }
    os << "]";
    return os;
  }
  
  template<typename real_t, int size>
  std::ostream& operator<<(std::ostream& os, intvec<real_t, size> const& x)
  {
    os << "[";
    for (int i=0; i<size; ++i) {
      if (i!=0) os << ",";
      os << x[i];
    }
    os << "]";
    return os;
  }
  
  template<typename real_t, int size>
  std::ostream& operator<<(std::ostream& os, realvec<real_t, size> const& x)
  {
    os << "[";
    for (int i=0; i<size; ++i) {
      if (i!=0) os << ",";
      os << x[i];
    }
    os << "]";
    return os;
  }
  
} // namespace vecmathlib

#endif  // #ifndef VEC_BASE_H
