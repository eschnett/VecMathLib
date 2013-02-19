// Instantiante some functions to be able to inspect the generated
// machine code

#define VML_NODEBUG

#include "vecmathlib.h"

using namespace std;

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
  template struct realvec<double,2>;
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
  template struct realvec<double,4>;
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
