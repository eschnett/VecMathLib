// -*-C++-*-

#ifndef VEC_MASK_H
#define VEC_MASK_H

#include <cstdlib>



namespace vecmathlib {
  
  template<typename realvec_t>
  class mask_t {
    
    friend realvec_t;
    
    typedef typename realvec_t::boolvec_t boolvec_t;
    typedef typename realvec_t::intvec_t intvec_t;
    static int const size = realvec_t::size;
    
    std::ptrdiff_t imin, imax;
    std::ptrdiff_t i;
    boolvec_t m;
    bool all_m;
    
  public:
    mask_t(boolvec_t m_): m(m_), all_m(all(m)) {}
    mask_t(std::ptrdiff_t imin_, std::ptrdiff_t imax_, std::ptrdiff_t ioff):
      imin(imin_), imax(imax_),
      i(imin - (ioff + imin) % size)
    {
      all_m = i>=imin && i<=imax-size;
      if (__builtin_expect(all_m, true)) {
        m = true;
      } else {
        m = (intvec_t(i) >= intvec_t(imin     ) - intvec_t::iota() &&
             intvec_t(i) <= intvec_t(imax-size) - intvec_t::iota());
      }
    }
    std::ptrdiff_t index() const { return i; }
    operator bool() const { return i<imax; }
    void operator++()
    {
      i += size;
      all_m = i<=imax-size;
      if (__builtin_expect(all_m, true)) {
        m = true;
      } else {
        m = intvec_t(i) <= intvec_t(imax-size) - intvec_t::iota();
      }
    }
  };
  
} // namespace vecmathlib

#endif  // #ifndef VEC_MASK_H
