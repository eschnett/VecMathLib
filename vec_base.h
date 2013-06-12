// -*-C++-*-

#ifndef VEC_BASE_H
#define VEC_BASE_H

#include <iostream>

#include "vec_mask.h"



namespace vecmathlib {
  
  
  
  template<typename T, int N>
  struct vec_types {
  };
  
  
  
  template<typename VT>
  typename VT::boolvec_t make_boolvec(typename VT::bool_t b);
  template<typename VT>
  typename VT::boolvec_t make_boolvec(const typename VT::bool_t* restrict bs);
  
  template<typename VT>
  typename VT::intvec_t make_intvec(typename VT::int_t b);
  template<typename VT>
  typename VT::intvec_t make_intvec(const typename VT::int_t* restrict is);
  template<typename VT>
  typename VT::intvec_t make_iota();
  
  template<typename VT>
  typename VT::realvec_t make_realvec(typename VT::real_t b);
  template<typename VT>
  typename VT::realvec_t make_realvec(const typename VT::real_t* restrict rs);
  
  
  
  template<typename VT>
  std::ostream& operator<<(std::ostream& os, typename VT::boolvec_t x)
  {
    os << "[";
    for (int i=0; i<VT::size; ++i) {
      if (i!=0) os << ",";
      os << x[i];
    }
    os << "]";
    return os;
  }
  
  template<typename VT>
  std::ostream& operator<<(std::ostream& os, typename VT::intvec_t x)
  {
    os << "[";
    for (int i=0; i<VT::size; ++i) {
      if (i!=0) os << ",";
      os << x[i];
    }
    os << "]";
    return os;
  }
  
  template<typename VT>
  std::ostream& operator<<(std::ostream& os, typename VT::realvec_t x)
  {
    os << "[";
    for (int i=0; i<VT::size; ++i) {
      if (i!=0) os << ",";
      os << x[i];
    }
    os << "]";
    return os;
  }
  
} // namespace vecmathlib

#endif  // #ifndef VEC_BASE_H
