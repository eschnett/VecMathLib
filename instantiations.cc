// Instantiante some functions to be able to inspect the generated
// machine code

#define VML_NODEBUG

#include "vecmathlib.h"



namespace vecmathlib {
  
  template<typename realvec_t, int n>
  typename realvec_t::real_t get_elt(realvec_t x)
  {
    return x[n];
  }
  template<typename realvec_t, int n>
  realvec_t set_elt(realvec_t x, typename realvec_t::real_t a)
  {
    return x.set_elt(n, a);
  }
  
  // template realbuiltinvec<float,1> fabs(realbuiltinvec<float,1> x);
  // template realbuiltinvec<float,1> fmin(realbuiltinvec<float,1> x, realbuiltinvec<float,1> y);
  // template intbuiltinvec<float,1> lsr(intbuiltinvec<float,1> x, intbuiltinvec<float,1>::int_t n);
  // template intbuiltinvec<double,1> lsr(intbuiltinvec<double,1> x, intbuiltinvec<double,1>::int_t n);
  // template intbuiltinvec<double,2> lsr(intbuiltinvec<double,2> x, intbuiltinvec<double,2>::int_t n);
  // template intbuiltinvec<double,2> lsr(intbuiltinvec<double,2> x, intbuiltinvec<double,2> n);
  // template realbuiltinvec<float,1> ifthen(realbuiltinvec<float,1>::boolvec_t c, realbuiltinvec<float,1> x, realbuiltinvec<float,1> y);
  // template realbuiltinvec<double,1> ifthen(realbuiltinvec<double,1>::boolvec_t c, realbuiltinvec<double,1> x, realbuiltinvec<double,1> y);
  // template realbuiltinvec<float,4> ifthen(realbuiltinvec<float,4>::boolvec_t c, realbuiltinvec<float,4> x, realbuiltinvec<float,4> y);
  // template realbuiltinvec<double,2> ifthen(realbuiltinvec<double,2>::boolvec_t c, realbuiltinvec<double,2> x, realbuiltinvec<double,2> y);
  
#ifdef VECMATHLIB_HAVE_VEC_FLOAT_1
  template realvec<float,1> round(realvec<float,1> x);
#endif
  
#ifdef VECMATHLIB_HAVE_VEC_DOUBLE_1
  template realvec<double,1> exp(realvec<double,1> x);
  template realvec<double,1> log(realvec<double,1> x);
  template realvec<double,1> sin(realvec<double,1> x);
  template realvec<double,1> sqrt(realvec<double,1> x);
  template realvec<double,1>::real_t get_elt<realvec<double,1>,0>(realvec<double,1> x);
  template realvec<double,1> set_elt<realvec<double,1>,0>(realvec<double,1> x, realvec<double,1>::real_t a);
#endif
  
#ifdef VECMATHLIB_HAVE_VEC_DOUBLE_2
  template realvec<double,2> exp(realvec<double,2> x);
  template realvec<double,2> log(realvec<double,2> x);
  template realvec<double,2> sin(realvec<double,2> x);
  template realvec<double,2> sqrt(realvec<double,2> x);
  template realvec<double,2>::real_t get_elt<realvec<double,2>,0>(realvec<double,2>);
  template realvec<double,2>::real_t get_elt<realvec<double,2>,1>(realvec<double,2>);
  template realvec<double,2> set_elt<realvec<double,2>,0>(realvec<double,2> x, realvec<double,2>::real_t a);
  template realvec<double,2> set_elt<realvec<double,2>,1>(realvec<double,2> x, realvec<double,2>::real_t a);
#endif
  
#ifdef VECMATHLIB_HAVE_VEC_DOUBLE_4
  template realvec<double,4> exp(realvec<double,4> x);
  template realvec<double,4> log(realvec<double,4> x);
  template realvec<double,4> sin(realvec<double,4> x);
  template realvec<double,4> sqrt(realvec<double,4> x);
  template realvec<double,4>::real_t get_elt<realvec<double,4>,0>(realvec<double,4>);
  template realvec<double,4>::real_t get_elt<realvec<double,4>,1>(realvec<double,4>);
  template realvec<double,4>::real_t get_elt<realvec<double,4>,2>(realvec<double,4>);
  template realvec<double,4>::real_t get_elt<realvec<double,4>,3>(realvec<double,4>);
  template realvec<double,4> set_elt<realvec<double,4>,0>(realvec<double,4> x, realvec<double,4>::real_t a);
  template realvec<double,4> set_elt<realvec<double,4>,1>(realvec<double,4> x, realvec<double,4>::real_t a);
  template realvec<double,4> set_elt<realvec<double,4>,2>(realvec<double,4> x, realvec<double,4>::real_t a);
  template realvec<double,4> set_elt<realvec<double,4>,3>(realvec<double,4> x, realvec<double,4>::real_t a);
#endif
  
}
