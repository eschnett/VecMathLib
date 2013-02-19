// Instantiante some functions to be able to inspect the generated
// machine code

#define VML_NODEBUG

#include "vecmathlib.h"

using namespace std;

namespace vecmathlib {
  
#ifdef VECMATHLIB_HAVE_VEC_DOUBLE_1
  template realvec<double,1> exp(realvec<double,1> x);
  template realvec<double,1> log(realvec<double,1> x);
  template realvec<double,1> sin(realvec<double,1> x);
  template realvec<double,1> sqrt(realvec<double,1> x);
#endif
  
#ifdef VECMATHLIB_HAVE_VEC_DOUBLE_2
  template realvec<double,2> exp(realvec<double,2> x);
  template realvec<double,2> log(realvec<double,2> x);
  template realvec<double,2> sin(realvec<double,2> x);
  template realvec<double,2> sqrt(realvec<double,2> x);
#endif
  
#ifdef VECMATHLIB_HAVE_VEC_DOUBLE_4
  template realvec<double,4> exp(realvec<double,4> x);
  template realvec<double,4> log(realvec<double,4> x);
  template realvec<double,4> sin(realvec<double,4> x);
  template realvec<double,4> sqrt(realvec<double,4> x);
#endif
  
}
