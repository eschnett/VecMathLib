// -*-C++-*-

#include "vecmathlib.h"

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <vector>

using namespace std;
using namespace vecmathlib;



// Assuming that xptr and yptr are aligned, but ldm can be arbitrary
template<typename realvec_t>
void smoothu(typename realvec_t::real_t const *xptr,
             typename realvec_t::real_t *yptr,
             ptrdiff_t m, ptrdiff_t ldm, ptrdiff_t n)
{
  typedef typename realvec_t::real_t real_t;
  typedef typename realvec_t::mask_t mask_t;
  for (ptrdiff_t j=1; j<n-1; ++j) {
    // Desired loop bounds
    const ptrdiff_t imin = 1;
    const ptrdiff_t imax = m-1;
    // Align actual loop iterations with vector size
    const ptrdiff_t ioff = ldm*j;
    for (mask_t mask(imin, imax, ioff); mask; ++mask) {
      const ptrdiff_t i = mask.index();
      const ptrdiff_t ij = ioff + i;
      const realvec_t x   = realvec_t::loada(xptr+ij);
      const realvec_t xil = realvec_t::loadu(xptr+ij, -1);
      const realvec_t xir = realvec_t::loadu(xptr+ij, +1);
      const realvec_t xjl = realvec_t::loadu(xptr+ij-ldm);
      const realvec_t xjr = realvec_t::loadu(xptr+ij+ldm);
      const realvec_t y =
        realvec_t(real_t(0.5)) * x +
        realvec_t(real_t(0.125)) * (xil + xir + xjl + xjr);
      y.storea(yptr+ij, mask);
    }
  }
}



// Assuming that xptr and yptr are aligned, and ldm is a multiple of
// the vector size
template<typename realvec_t>
void smootha(typename realvec_t::real_t const *xptr,
             typename realvec_t::real_t *yptr,
             ptrdiff_t m, ptrdiff_t ldm, ptrdiff_t n)
{
  typedef typename realvec_t::real_t real_t;
  typedef typename realvec_t::mask_t mask_t;
  assert(ldm % realvec_t::size == 0);
  for (ptrdiff_t j=1; j<n-1; ++j) {
    // Desired loop bounds
    const ptrdiff_t imin = 1;
    const ptrdiff_t imax = m-1;
    // Align actual loop iterations with vector size
    const ptrdiff_t ioff = ldm*j;
    for (mask_t mask(imin, imax, ioff); mask; ++mask) {
      const ptrdiff_t i = mask.index();
      const ptrdiff_t ij = ioff + i;
      const realvec_t x   = realvec_t::loada(xptr+ij);
      const realvec_t xil = realvec_t::loadu(xptr+ij, -1);
      const realvec_t xir = realvec_t::loadu(xptr+ij, +1);
      const realvec_t xjl = realvec_t::loada(xptr+ij-ldm);
      const realvec_t xjr = realvec_t::loada(xptr+ij+ldm);
      const realvec_t y =
        realvec_t(real_t(0.5)) * x +
        realvec_t(real_t(0.125)) * (xil + xir + xjl + xjr);
      y.storea(yptr+ij, mask);
    }
  }
}



static size_t align_up(size_t i, size_t size)
{
  return (i + size - 1) / size * size;
}

int main(int argc, char** argv)
{
  const ptrdiff_t m = 100;
  const ptrdiff_t n = 100;
  
#if defined VECMATHLIB_HAVE_VEC_DOUBLE_4
  typedef realvec<double,4> realvec_t;
#elif defined VECMATHLIB_HAVE_VEC_DOUBLE_2
  typedef realvec<double,2> realvec_t;
#else
  typedef realpseudovec<double,1> realvec_t;
#endif
  
  const ptrdiff_t ldm = align_up(m, realvec_t::size);
  typedef realvec_t::real_t real_t;
  vector<real_t> x(ldm*n, 0.0), y(ldm*n, 0.0);
  
  smoothu<realvec_t>(&x[0], &y[0], m, ldm, n);
  smootha<realvec_t>(&x[0], &y[0], m, ldm, n);
  
  return 0;
}
