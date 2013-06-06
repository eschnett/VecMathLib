// -*-C++-*-

#include "vecmathlib.h"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
#include <sstream>

using namespace std;


  
int num_errors = 0;



template<typename realvec_t>
struct vecmathlib_test {
  
  typedef typename realvec_t::boolvec_t boolvec_t;
  typedef typename realvec_t::intvec_t intvec_t;
  
  typedef typename realvec_t::int_t int_t;
  typedef typename realvec_t::uint_t uint_t;
  typedef typename realvec_t::real_t real_t;
  
  // Short names for type casts
  typedef real_t R;
  typedef int_t I;
  typedef uint_t U;
  typedef realvec_t RV;
  typedef intvec_t IV;
  typedef boolvec_t BV;
  
  typedef vecmathlib::floatprops<real_t> FP;
  typedef vecmathlib::mathfuncs<realvec_t> MF;
  
  
  
  // Test each function with this many random values
  static const int imax = 10000;
  // Require (arbitrarily) that 3/4 of the digits are correct
  static real_t accuracy() { return pow(realvec_t::epsilon(), R(0.75)); }
  
  
  
  static realvec_t random(const real_t xmin, const real_t xmax)
  {
    realvec_t x;
    for (int i=0; i<realvec_t::size; ++i) {
      const real_t r =
        (xmax - xmin) *FP::convert_float(rand()) / FP::convert_float(RAND_MAX);
      x.set_elt(i, xmin + r);
    }
    return x;
  }
  
  static intvec_t random(const int_t nmin, const int_t nmax)
  {
    intvec_t n;
    for (int i=0; i<intvec_t::size; ++i) {
      const real_t r =
        R(nmax - nmin + 1) *
        R(rand()) / (R(RAND_MAX) + R(1.0));
      n.set_elt(i, nmin + FP::convert_int(floor(r)));
    }
    return n;
  }
  
  
  
  static bool is_big_endian()
  {
    const int i = 1;
    unsigned char cs[sizeof i];
    memcpy(cs, &i, sizeof i);
    return cs[0]==0;
  }
  
  template<typename T>
  static string hex(const T x)
  {
    unsigned char cs[sizeof x];
    memcpy(cs, &x, sizeof x);
    ostringstream buf;
    buf << "0x";
    const char* const hexdigits = "0123456789abcdef";
    const int n0 = is_big_endian() ? 0 : sizeof x - 1;
    const int dn = is_big_endian() ? +1 : -1;
    const int n1 = n0 + sizeof x * dn;
    for (int n=n0; n!=n1; n+=dn) {
      buf << hexdigits[cs[n]>>4] << hexdigits[cs[n]&15];
    }
    return buf.str();
  }
  
  
  
  static boolvec_t supported(realvec_t x)
  {
    return x==RV(0.0) || MF::vml_ieee_isnormal(x)
#ifdef VML_HAVE_DENORMALS
      || MF::vml_ieee_isfinite(x)
#endif
#ifdef VML_HAVE_INF
      || MF::vml_ieee_isinf(x)
#endif
#ifdef VML_HAVE_NAN
      || MF::vml_ieee_isnan(x)
#endif
      ;
  }
  
  static boolvec_t supported(intvec_t x)
  {
    return true;
  }
  
  static boolvec_t supported(boolvec_t x)
  {
    return true;
  }
  
  
  
  // Check load memory access
  static void check_mem(const char* const func,
                        const realvec_t x,
                        const real_t* const p,
                        const realvec_t xold,
                        const int mval)
  {
    realvec_t xwant;
    for (int i=0; i<realvec_t::size; ++i) {
      xwant.set_elt(i, mval & (1<<i) ? p[i] : xold[i]);
    }
    const boolvec_t isbad = x != xwant;
    if (any(isbad)) {
      ++ num_errors;
      cout << setprecision(realvec_t::digits10+2)
           << "Error in " << func << ":\n"
           << "   found=" << x << " [" << hex(x) << "]\n"
           << "   expected=" << xwant << " [" << hex(xwant) << "]\n"
           << "   isbad=" << isbad << "\n"
           << flush;
    }
  }
  
  // Check store memory access
  static void check_mem(const char* const func,
                        const real_t* const p,
                        const realvec_t x,
                        const real_t* const pold,
                        const int mval)
  {
    realvec_t pv, pvwant;
    for (int i=0; i<realvec_t::size; ++i) {
      pv.set_elt(i, p[i]);
      pvwant.set_elt(i, mval & (1<<i) ? x[i] : pold[i]);
    }
    const boolvec_t isbad = pv != pvwant;
    if (any(isbad)) {
      ++ num_errors;
      cout << setprecision(realvec_t::digits10+2)
           << "Error in " << func << ":\n"
           << "   found=" << pv << " [" << hex(pv) << "]\n"
           << "   expected=" << pvwant << " [" << hex(pvwant) << "]\n"
           << "   isbad=" << isbad << "\n"
           << flush;
    }
  }
  
  static void check_bool(const char* const func,
                         const bool rstd, const bool rvml)
  {
    const bool dr = rstd ^ rvml;
    const bool isbad = dr;
    if (isbad) {
      ++ num_errors;
      cout << "Error in " << func << ":\n"
           << "   fstd()=" << rstd << " [" << hex(rstd) << "]\n"
           << "   fvml()=" << rvml << " [" << hex(rvml) << "]\n"
           << "   isbad()=" << isbad << "\n"
           << flush;
    }
  }
  
  template<typename A>
  static void check_bool(const char* const func,
                         const bool rstd, const bool rvml, const A x)
  {
    const bool dr = rstd ^ rvml;
    const bool isbad = dr;
    if (isbad) {
      ++ num_errors;
      cout << "Error in " << func << ":\n"
           << "   x=" << x << " [" << hex(x) << "]\n"
           << "   fstd(x)=" << rstd << " [" << hex(rstd) << "]\n"
           << "   fvml(x)=" << rvml << " [" << hex(rvml) << "]\n"
           << "   isbad(x)=" << isbad << "\n"
           << flush;
    }
  }
  
  template<typename A>
  static void check_bool(const char* const func,
                         const boolvec_t rstd, const boolvec_t rvml,
                         const A x)
  {
    boolvec_t dr;
    bool isbad = false;
    for (int i=0; i<realvec_t::size; ++i) {
      dr.set_elt(i, rstd[i] ^ rvml[i]);
      isbad |= dr[i];
    }
    if (isbad) {
      ++ num_errors;
      cout << "Error in " << func << ":\n"
           << "   x=" << x << " [" << hex(x) << "]\n"
           << "   fstd(x)=" << rstd << " [" << hex(rstd) << "]\n"
           << "   fvml(x)=" << rvml << " [" << hex(rvml) << "]\n"
           << "   error(x)=" << dr << " [" << hex(rvml) << "]\n"
           << "   isbad(x)=" << isbad << "\n"
           << flush;
    }
  }
  
  template<typename A, typename B>
  static void check_bool(const char* const func,
                         const boolvec_t rstd, const boolvec_t rvml,
                         const A x, const B y)
  {
    boolvec_t dr;
    bool isbad = false;
    for (int i=0; i<realvec_t::size; ++i) {
      dr.set_elt(i, rstd[i] ^ rvml[i]);
      isbad |= dr[i];
    }
    if (isbad) {
      ++ num_errors;
      cout << "Error in " << func << ":\n"
           << "   x=" << x << " [" << hex(x) << "]\n"
           << "   y=" << y << " [" << hex(y) << "]\n"
           << "   fstd(x,y)=" << rstd << " [" << hex(rstd) << "]\n"
           << "   fvml(x,y)=" << rvml << " [" << hex(rvml) << "]\n"
           << "   error(x,y)=" << dr << " [" << hex(rvml) << "]\n"
           << "   isbad(x,y)=" << isbad << "\n"
           << flush;
    }
  }
  
  template<typename A>
  static void check_bool(const char* const func,
                         bool fstd(typename A::scalar_t x),
                         boolvec_t fvml(A x),
                         const A x)
  {
    boolvec_t rstd;
    for (int i=0; i<boolvec_t::size; ++i) {
      rstd.set_elt(i, fstd(x[i]));
    }
    const boolvec_t rvml = fvml(x);
    const boolvec_t dr = rstd != rvml;
    const boolvec_t isbad = supported(x) && supported(rstd) && dr;
    if (any(isbad)) {
      ++ num_errors;
      cout << setprecision(realvec_t::digits10+2)
           << "Error in " << func << ":\n"
           << "   x=" << x << " [" << hex(x) << "]\n"
           << "   fstd(x)=" << rstd << " [" << hex(rstd) << "]\n"
           << "   fvml(x)=" << rvml << " [" << hex(rvml) << "]\n"
	   << "   error(x)=" << dr << " [" << hex(dr) << "]\n"
           << "   isbad(x)=" << isbad << "\n"
           << flush;
    }
  }
  
  template<typename A, typename B>
  static void check_bool(const char* const func,
                         bool fstd(typename A::scalar_t x,
                                   typename B::scalar_t y),
                         boolvec_t fvml(A x, B y),
                         const A x, const B y)
  {
    boolvec_t rstd;
    for (int i=0; i<boolvec_t::size; ++i) {
      rstd.set_elt(i, fstd(x[i], y[i]));
    }
    const boolvec_t rvml = fvml(x, y);
    const boolvec_t dr = rstd != rvml;
    const boolvec_t isbad = supported(x) && supported(rstd) && dr;
    if (any(isbad)) {
      ++ num_errors;
      cout << setprecision(realvec_t::digits10+2)
           << "Error in " << func << ":\n"
           << "   x=" << x << " [" << hex(x) << "]\n"
           << "   y=" << y << " [" << hex(y) << "]\n"
           << "   fstd(x,y)=" << rstd << " [" << hex(rstd) << "]\n"
           << "   fvml(x,y)=" << rvml << " [" << hex(rvml) << "]\n"
	   << "   error(x,y)=" << dr << " [" << hex(dr) << "]\n"
           << "   isbad(x,y)=" << isbad << "\n"
           << flush;
    }
  }
  
  template<typename A>
  static void check_int(const char* const func,
                        int_t fstd(typename A::scalar_t x),
                        intvec_t fvml(A x),
                        const A x)
  {
    intvec_t rstd;
    for (int i=0; i<intvec_t::size; ++i) {
      rstd.set_elt(i, fstd(x[i]));
    }
    const intvec_t rvml = fvml(x);
    const intvec_t dr = rstd - rvml;
    const boolvec_t isbad = supported(x) && supported(rstd) && convert_bool(dr);
    if (any(isbad)) {
      ++ num_errors;
      cout << setprecision(realvec_t::digits10+2)
           << "Error in " << func << ":\n"
           << "   x=" << x << " [" << hex(x) << "]\n"
           << "   fstd(x)=" << rstd << " [" << hex(rstd) << "]\n"
           << "   fvml(x)=" << rvml << " [" << hex(rvml) << "]\n"
	   << "   error(x)=" << dr << " [" << hex(dr) << "]\n"
           << "   isbad(x)=" << isbad << "\n"
           << flush;
    }
  }
  
  template<typename A, typename B>
  static void check_int(const char* const func,
                        int_t fstd(typename A::scalar_t x, B y),
                        intvec_t fvml(A x, B y),
                        const A x, const B y)
  {
    intvec_t rstd;
    for (int i=0; i<intvec_t::size; ++i) {
      rstd.set_elt(i, fstd(x[i], y));
    }
    const intvec_t rvml = fvml(x, y);
    const intvec_t dr = rstd - rvml;
    const boolvec_t isbad = supported(x) && supported(rstd) && convert_bool(dr);
    if (any(isbad)) {
      ++ num_errors;
      cout << setprecision(realvec_t::digits10+2)
           << "Error in " << func << ":\n"
           << "   x=" << x << " [" << hex(x) << "]\n"
           << "   y=" << y << " [" << hex(y) << "]\n"
           << "   fstd(x,y)=" << rstd << " [" << hex(rstd) << "]\n"
           << "   fvml(x,y)=" << rvml << " [" << hex(rvml) << "]\n"
	   << "   error(x,y)=" << dr << " [" << hex(dr) << "]\n"
           << "   isbad(x,y)=" << isbad << "\n"
           << flush;
    }
  }
  
  template<typename A, typename B>
  static void check_int(const char* const func,
                        int_t fstd(typename A::scalar_t x,
                                   typename B::scalar_t y),
                        intvec_t fvml(A x, B y),
                        const A x, const B y)
  {
    intvec_t rstd;
    for (int i=0; i<intvec_t::size; ++i) {
      rstd.set_elt(i, fstd(x[i], y[i]));
    }
    const intvec_t rvml = fvml(x, y);
    const intvec_t dr = rstd - rvml;
    const boolvec_t isbad =
      supported(x) && supported(y) && supported(rstd) && convert_bool(dr);
    if (any(isbad)) {
      ++ num_errors;
      cout << setprecision(realvec_t::digits10+2)
           << "Error in " << func << ":\n"
           << "   x=" << x << " [" << hex(x) << "]\n"
           << "   y=" << y << " [" << hex(y) << "]\n"
           << "   fstd(x,y)=" << rstd << " [" << hex(rstd) << "]\n"
           << "   fvml(x,y)=" << rvml << " [" << hex(rvml) << "]\n"
	   << "   error(x,y)=" << dr << " [" << hex(dr) << "]\n"
           << "   isbad(x,y)=" << isbad << "\n"
           << flush;
    }
  }
  
  template<typename A, typename B, typename C>
  static void check_int(const char* const func,
                        int_t fstd(typename A::scalar_t x,
                                   typename B::scalar_t y,
                                   typename C::scalar_t z),
                        intvec_t fvml(A x, B y, C z),
                        const A x, const B y, const C z)
  {
    intvec_t rstd;
    for (int i=0; i<intvec_t::size; ++i) {
      rstd.set_elt(i, fstd(x[i], y[i], z[i]));
    }
    const intvec_t rvml = fvml(x, y, z);
    const intvec_t dr = rstd - rvml;
    const boolvec_t isbad =
      supported(x) && supported(y) && supported(z) && supported(rstd) &&
      convert_bool(dr);
    if (any(isbad)) {
      ++ num_errors;
      cout << setprecision(realvec_t::digits10+2)
           << "Error in " << func << ":\n"
           << "   x=" << x << " [" << hex(x) << "]\n"
           << "   y=" << y << " [" << hex(y) << "]\n"
           << "   z=" << z << " [" << hex(z) << "]\n"
           << "   fstd(x,y,z)=" << rstd << " [" << hex(rstd) << "]\n"
           << "   fvml(x,y,z)=" << rvml << " [" << hex(rvml) << "]\n"
	   << "   error(x,y,z)=" << dr << " [" << hex(dr) << "]\n"
           << "   isbad(x,y,z)=" << isbad << "\n"
           << flush;
    }
  }
  
  template<typename A>
  static void check_real(const char* const func,
                         real_t fstd(typename A::scalar_t x),
                         realvec_t fvml(A x),
                         const A x,
                         const real_t accuracy)
  {
    realvec_t rstd;
    for (int i=0; i<realvec_t::size; ++i) {
      rstd.set_elt(i, fstd(x[i]));
    }
    const realvec_t rvml = fvml(x);
    const realvec_t dr = rstd - rvml;
    const realvec_t scale = fabs(rstd) + fabs(rvml) + realvec_t(1.0);
    const boolvec_t isbad =
      supported(x) && supported(rstd) &&
      fabs(dr) > realvec_t(accuracy) * scale;
    if (any(isbad)) {
      ++ num_errors;
      cout << setprecision(realvec_t::digits10+2)
           << "Error in " << func << "(" << x << " "
           << "[" << hex(x) << "]):\n"
           << "   fstd(x)=" << rstd << " [" << hex(rstd) << "]\n"
           << "   fvml(x)=" << rvml << " [" << hex(rvml) << "]\n"
           << "   abs-error(x)=" << fabs(dr) << "\n"
           << "   rel-error(x)=" << fabs(dr) / scale << "\n"
           << "   isbad(x)=" << isbad << "\n"
           << "   accuracy=" << accuracy << "\n"
           << flush;
    }
  }
  
  template<typename A, typename B>
  static void check_real(const char* const func,
                         real_t fstd(typename A::scalar_t x,
                                     typename B::scalar_t y),
                         realvec_t fvml(A x, B y),
                         const A x, const B y,
                         const real_t accuracy)
  {
    realvec_t rstd;
    for (int i=0; i<realvec_t::size; ++i) {
      rstd.set_elt(i, fstd(x[i], y[i]));
    }
    const realvec_t rvml = fvml(x, y);
    const realvec_t dr = rstd - rvml;
    const realvec_t scale = fabs(rstd) + fabs(rvml) + realvec_t(1.0);
    const boolvec_t isbad =
      supported(x) && supported(y) && supported(rstd) &&
      fabs(dr) > realvec_t(accuracy) * scale;
    if (any(isbad)) {
      ++ num_errors;
      cout << setprecision(realvec_t::digits10+2)
           << "Error in " << func << "(" << x << "," << y << " "
           << "[" << hex(x) << "],[" << hex(y) << "]):\n"
           << "   fstd(x,y)=" << rstd << " [" << hex(rstd) << "]\n"
           << "   fvml(x,y)=" << rvml << " [" << hex(rvml) << "]\n"
           << "   abs-error(x,y)=" << fabs(dr) << "\n"
           << "   rel-error(x,y)=" << fabs(dr) / scale << "\n"
           << "   isbad(x,y)=" << isbad << "\n"
           << "   accuracy=" << accuracy << "\n"
           << flush;
    }
  }
  
  template<typename A, typename B, typename C>
  static void check_real(const char* const func,
                         real_t fstd(typename A::scalar_t x,
                                     typename B::scalar_t y,
                                     typename C::scalar_t z),
                         realvec_t fvml(A x, B y, C z),
                         const A x, const B y, C const z,
                         const real_t accuracy)
  {
    realvec_t rstd;
    for (int i=0; i<realvec_t::size; ++i) {
      rstd.set_elt(i, fstd(x[i], y[i], z[i]));
    }
    const realvec_t rvml = fvml(x, y, z);
    const realvec_t dr = rstd - rvml;
    const realvec_t scale = fabs(rstd) + fabs(rvml) + realvec_t(1.0);
    const boolvec_t isbad =
      supported(x) && supported(y) && supported(z) && supported(rstd) &&
      fabs(dr) > realvec_t(accuracy) * scale;
    if (any(isbad)) {
      ++ num_errors;
      cout << setprecision(realvec_t::digits10+2)
           << "Error in " << func << "(" << x << "," << y<< "," << z << " "
           << "[" << hex(x) << "," << hex(y) << "," << hex(z) << "]):\n"
           << "   fstd(x,y,z)=" << rstd << " [" << hex(rstd) << "]\n"
           << "   fvml(x,y,z)=" << rvml << " [" << hex(rvml) << "]\n"
           << "   abs-error(x,y,z)=" << fabs(dr) << "\n"
           << "   rel-error(x,y,z)=" << fabs(dr) / scale << "\n"
           << "   isbad(x,y,z)=" << isbad << "\n"
           << "   accuracy=" << accuracy << "\n"
           << flush;
    }
  }
  
  
  
  static real_t* align_mem(real_t* p)
  {
    const ptrdiff_t alignment = sizeof(realvec_t);
    p = (real_t*)((intptr_t(p) + alignment-1) & -alignment);
    assert(intptr_t(p) % alignment == 0);
    return p;
  }
  static string add_suffix(const char* str, int i)
  {
    ostringstream buf;
    buf << str << "." << i;
    return buf.str();
  }
  static void test_mem()
  {
    cout << "   testing loada loadu storea storeu (errors may lead to segfaults)...\n" << flush;
    const int n = 4;
    const int sz = realvec_t::size;
    const int nbytes = n*sz*sizeof(real_t);
    real_t* const x = align_mem(new real_t[(n+1)*sz]);
    real_t* const xnew = align_mem(new real_t[(n+1)*sz]);
    for (int i=0; i<n; ++i) {
      realvec_t xv = random(R(-10.0), R(+10.0));
      memcpy(&x[i*sz], &xv, sizeof xv);
    }
    const realvec_t z = random(R(-10.0), R(+10.0));
    
    // loada
    {
      const real_t *p = &x[sz];
      realvec_t y = realvec_t::loada(p);
      check_mem("loada", y, p, z, ~0);
    }
    
    // loadu
    for (ptrdiff_t i=0; i<realvec_t::size; ++i) {
      const real_t *p = &x[sz];
      realvec_t y = realvec_t::loadu(p+i);
      check_mem(add_suffix("loadu", i).c_str(), y, p+i, z, ~0);
    }
    
    // loadu(ioff)
    for (ptrdiff_t ioff=0; ioff<realvec_t::size; ++ioff) {
      const real_t *p = &x[sz];
      realvec_t y = realvec_t::loadu(p, ioff);
      check_mem(add_suffix("loadu(ioff)", ioff).c_str(), y, p+ioff, z, ~0);
    }
    
    // storea
    {
      memcpy(xnew, x, nbytes);
      real_t *p = &xnew[sz];
      storea(z, p);
      check_mem("storea", p, z, &x[sz], ~0);
    }
    
    // storeu
    for (ptrdiff_t i=0; i<realvec_t::size; ++i) {
      memcpy(xnew, x, nbytes);
      real_t *p = &xnew[sz];
      storeu(z, p+i);
      check_mem(add_suffix("storeu", i).c_str(), p+i, z, &x[sz+i], ~0);
    }
    
    // storeu
    for (ptrdiff_t ioff=0; ioff<realvec_t::size; ++ioff) {
      memcpy(xnew, x, nbytes);
      real_t *p = &xnew[sz];
      storeu(z, p, ioff);
      check_mem(add_suffix("storeu(ioff)", ioff).c_str(),
                p+ioff, z, &x[sz+ioff], ~0);
    }
    
    for (int mval=0; mval<(1<<realvec_t::size); ++mval) {
      boolvec_t mbool;
      for (int i=0; i<realvec_t::size; ++i) mbool.set_elt(i, mval & (1<<i));
      typename realvec_t::mask_t mask(mbool);
      
      // loada(mask)
      {
        const real_t *p = &x[sz];
        realvec_t y = loada(p, z, mask);
        check_mem("loada(mask)", y, p, z, mval);
      }
      
      // loadu(mask)
      for (ptrdiff_t i=0; i<realvec_t::size; ++i) {
        const real_t *p = &x[sz];
        realvec_t y = loadu(p+i, z, mask);
        check_mem("loadu(mask)", y, p+i, z, mval);
      }
      
      // loadu(ioff, mask)
      for (ptrdiff_t ioff=0; ioff<realvec_t::size; ++ioff) {
        const real_t *p = &x[sz];
        realvec_t y = loadu(p, ioff, z, mask);
        check_mem("loadu(ioff,mask)", y, p+ioff, z, mval);
      }
      
      // storea
      {
	memcpy(xnew, x, nbytes);
        real_t *p = &xnew[sz];
        storea(z, p, mask);
        check_mem("storea(mask)", p, z, &x[sz], mval);
      }
      
      // storeu
      for (ptrdiff_t i=0; i<realvec_t::size; ++i) {
	memcpy(xnew, x, nbytes);
        real_t *p = &xnew[sz];
        storeu(z, p+i, mask);
        check_mem("storeu(mask)", p+i, z, &x[sz+i], mval);
      }
      
      // storeu
      for (ptrdiff_t ioff=0; ioff<realvec_t::size; ++ioff) {
	memcpy(xnew, x, nbytes);
        real_t *p = &xnew[sz];
        storeu(z, p, ioff, mask);
        check_mem("storeu(ioff,mask)", p+ioff, z, &x[sz+ioff], mval);
      }
      
    } // for mval
  }
  
  
  
  template<typename T>
  static T local_ifthen(bool b, T x, T y) { return b ? x : y; }
  static void test_bool()
  {
    cout << "   testing boolean operations...\n" << flush;
    
    const boolvec_t bf = boolvec_t(false);
    const boolvec_t bt = boolvec_t(true);
    for (int i=0; i<realvec_t::size; ++i) {
      check_bool("false", false, bf[i]);
      check_bool("true", true, bt[i]);
    }
    check_bool("all", false, all(bf), false);
    check_bool("all", true, all(bt), true);
    check_bool("any", false, any(bf), false);
    check_bool("any", true, any(bt), true);
    
    boolvec_t b0 = bt;
    boolvec_t b1 = bf;
    for (int n=0; n<realvec_t::size; ++n) {
      b0.set_elt(n, false);
      b1.set_elt(n, true);
      for (int i=0; i<realvec_t::size; ++i) {
        check_bool("set_elt", i<=n ? false : true, b0[i], false);
        check_bool("set_elt", i<=n ? true : false, b1[i], true);
      }
    }
    
    for (int n=0; n<(1<<realvec_t::size); ++n) {
      boolvec_t x;
      for (int i=0; i<realvec_t::size; ++i) {
        x.set_elt(i, n & (1<<i));
      }
      
      {
        boolvec_t rstd;
        for (int i=0; i<realvec_t::size; ++i) {
          rstd.set_elt(i, !x[i]);
        }
        boolvec_t rvml = !x;
        check_bool("!", rstd, rvml, x);
      }
      {
        bool rstd = true;
        for (int i=0; i<realvec_t::size; ++i) {
          rstd &= x[i];
        }
        bool rvml = all(x);
        check_bool("all", rstd, rvml, x);
      }
      {
        bool rstd = false;
        for (int i=0; i<realvec_t::size; ++i) {
          rstd |= x[i];
        }
        bool rvml = any(x);
        check_bool("any", rstd, rvml, x);
      }
      check_int("ifthen(int)",
                local_ifthen<int_t>,
                (intvec_t(*)(boolvec_t,intvec_t,intvec_t))vecmathlib::ifthen,
                x, IV(I(1)), IV(I(2)));
      check_real("ifthen(real)",
                 local_ifthen<real_t>,
                 ((realvec_t(*)(boolvec_t,realvec_t,realvec_t))
                  vecmathlib::ifthen),
                x, RV(1.0), RV(2.0), R(0.0));
    }
    
    for (int n=0; n<(1<<realvec_t::size); ++n) {
      for (int m=0; m<(1<<realvec_t::size); ++m) {
        boolvec_t x, y;
        for (int i=0; i<realvec_t::size; ++i) {
          x.set_elt(i, n & (1<<i));
          y.set_elt(i, m & (1<<i));
        }
        
        {
          boolvec_t rstd;
          for (int i=0; i<realvec_t::size; ++i) {
            rstd.set_elt(i, x[i] && y[i]);
          }
          boolvec_t rvml = x && y;
          check_bool("&&", rstd, rvml, x, y);
        }
        {
          boolvec_t rstd;
          for (int i=0; i<realvec_t::size; ++i) {
            rstd.set_elt(i, x[i] || y[i]);
          }
          boolvec_t rvml = x || y;
          check_bool("||", rstd, rvml, x, y);
        }
        {
          boolvec_t rstd;
          for (int i=0; i<realvec_t::size; ++i) {
            rstd.set_elt(i, x[i] == y[i]);
          }
          boolvec_t rvml = x == y;
          check_bool("==", rstd, rvml, x, y);
        }
        {
          boolvec_t rstd;
          for (int i=0; i<realvec_t::size; ++i) {
            rstd.set_elt(i, x[i] != y[i]);
          }
          boolvec_t rvml = x != y;
          check_bool("!=", rstd, rvml, x, y);
        }
      }
    }
  }
  
  
  
  static bool local_convert_bool(int_t x) { return x; }
  static int_t local_convert_int(bool x) { return x; }
  template<typename T> static T local_pos(T x) { return +x; }
  template<typename T> static T local_neg(T x) { return -x; }
  template<typename T> static T local_not(T x) { return ~x; }
  template<typename T> static T local_add(T x, T y) { return x+y; }
  template<typename T> static T local_sub(T x, T y) { return x-y; }
  template<typename T> static T local_and(T x, T y) { return x&y; }
  template<typename T> static T local_or(T x, T y) { return x|y; }
  template<typename T> static T local_xor(T x, T y) { return x^y; }
  
  static int_t local_lsr(int_t x, int_t y) { return uint_t(x)>>uint_t(y); }
  template<typename T> static T local_srs(T x, typename T::scalar_t y)
  {
    return x>>y;
  }
  template<typename T> static T local_sls(T x, typename T::scalar_t y)
  {
    return x<<y;
  }
  template<typename T> static T local_sr(T x, T y) { return x>>y; }
  template<typename T> static T local_sl(T x, T y) { return x<<y; }

  template<typename T> static bool local_eq(T x, T y) { return x==y; }
  template<typename T> static bool local_ne(T x, T y) { return x!=y; }
  template<typename T> static bool local_lt(T x, T y) { return x<y; }
  template<typename T> static bool local_le(T x, T y) { return x<=y; }
  template<typename T> static bool local_gt(T x, T y) { return x>y; }
  template<typename T> static bool local_ge(T x, T y) { return x>=y; }
  template<typename T> static boolvec_t local_veq(T x, T y) { return x==y; }
  template<typename T> static boolvec_t local_vne(T x, T y) { return x!=y; }
  template<typename T> static boolvec_t local_vlt(T x, T y) { return x<y; }
  template<typename T> static boolvec_t local_vle(T x, T y) { return x<=y; }
  template<typename T> static boolvec_t local_vgt(T x, T y) { return x>y; }
  template<typename T> static boolvec_t local_vge(T x, T y) { return x>=y; }
  static void test_int()
  {
    cout << "   testing integer operations...\n" << flush;
    
    const int_t int_min = std::numeric_limits<int_t>::min();
    const int_t int_max = std::numeric_limits<int_t>::max();
    const int_t values[] = {
      0, 1, 2, 3, -1, -2, -3,
      int_min, int_min+1, int_min+2, int_min+3,
      int_max, int_max-1, int_max-2, int_max-3,
    };
    const int nvalues = sizeof values / sizeof *values;

    for (int i=0; i<nvalues*nvalues+2*imax; ++i) {
      intvec_t x, y;
      if (i<nvalues*nvalues) {
        x = values[i%nvalues];
        y = values[i/nvalues];
      } else if (i<nvalues*nvalues+imax) {
        x = random(I(-100), I(+100));
        y = random(I(-100), I(+100));
      } else {
        x = random(int_min/2, int_max/2);
        y = random(int_min/2, int_max/2);
      }
      boolvec_t b = convert_bool(random(I(0), I(1)));
      
      check_bool("convert_bool(int)",
                 local_convert_bool,
                 (boolvec_t(*)(intvec_t))vecmathlib::convert_bool,
                  x);
      check_int("convert_int(bool)",
                local_convert_int,
                (intvec_t(*)(boolvec_t))vecmathlib::convert_int,
                b);
      
      check_int("+", local_pos<int_t>, local_pos<intvec_t>, x);
      check_int("-", local_neg<int_t>, local_neg<intvec_t>, x);
      check_int("~", local_not<int_t>, local_not<intvec_t>, x);
      
      check_int("+", local_add<int_t>, local_add<intvec_t>, x, y);
      check_int("-", local_sub<int_t>, local_sub<intvec_t>, x, y);
      check_int("&", local_and<int_t>, local_and<intvec_t>, x, y);
      check_int("|", local_or<int_t>, local_or<intvec_t>, x, y);
      check_int("^", local_xor<int_t>, local_xor<intvec_t>, x, y);
      
      const int_t bits = 8*sizeof(int_t);
      check_int("lsr",
                local_lsr, (intvec_t(*)(intvec_t,int_t))vecmathlib::lsr,
                x, y[0] & (bits-1));
      check_int(">>", local_sr<int_t>, local_srs<intvec_t>, x, y[0] & (bits-1));
      check_int("<<", local_sl<int_t>, local_sls<intvec_t>, x, y[0] & (bits-1));
      check_int("lsr", local_lsr, vecmathlib::lsr, x, y & IV(bits-1));
      check_int(">>", local_sr<int_t>, local_sr<intvec_t>, x, y & IV(bits-1));
      check_int("<<", local_sl<int_t>, local_sl<intvec_t>, x, y & IV(bits-1));
      
      check_bool("==", local_eq<int_t>, local_veq<intvec_t>, x, y);
      check_bool("!=", local_ne<int_t>, local_vne<intvec_t>, x, y);
      check_bool("<", local_lt<int_t>, local_vlt<intvec_t>, x, y);
      check_bool("<=", local_le<int_t>, local_vle<intvec_t>, x, y);
      check_bool(">", local_gt<int_t>, local_vgt<intvec_t>, x, y);
      check_bool(">=", local_ge<int_t>, local_vge<intvec_t>, x, y);
    }
  }
  
  // Change signature: "int" -> "int_t"
  static int_t local_ilogb(real_t x)
  {
    int r = std::ilogb(x);
    if (r==FP_ILOGB0) return std::numeric_limits<int_t>::min();
    if (r==FP_ILOGBNAN) return std::numeric_limits<int_t>::max();
    return r;
  }
  static real_t local_ldexp(real_t x, int_t n) { return ldexp(x, n); }
  static void test_fabs()
  {
    cout << "   testing copysign fabs fdim fma fmax fmin ilogb isfinite isinf isnan isnormal ldexp signbit...\n" << flush;
    
    const real_t eps = FP::epsilon();
    const real_t int_min = R(std::numeric_limits<int_t>::min());
    const real_t int_max = R(std::numeric_limits<int_t>::max());
    const real_t uint_min = R(std::numeric_limits<uint_t>::min());
    const real_t uint_max = R(std::numeric_limits<uint_t>::max());
    const real_t values[] = {
      R(+0.0), R(+0.1), R(+0.9), R(+1.0), R(+1.1),
      R(-0.0), R(-0.1), R(-0.9), R(-1.0), R(-1.1),
      R(+0.0)+eps, R(+0.1)+eps, R(+0.9)+eps, R(+1.0)+eps, R(+1.1)+eps,
      R(-0.0)+eps, R(-0.1)+eps, R(-0.9)+eps, R(-1.0)+eps, R(-1.1)+eps,
      R(+0.0)-eps, R(+0.1)-eps, R(+0.9)-eps, R(+1.0)-eps, R(+1.1)-eps,
      R(-0.0)-eps, R(-0.1)-eps, R(-0.9)-eps, R(-1.0)-eps, R(-1.1)-eps,
#ifdef VML_HAVE_DENORMALS
      +FP::min(), +FP::min()*(R(1.0)+eps), +FP::min()*R(2.0),
      -FP::min(), -FP::min()*(R(1.0)+eps), -FP::min()*R(2.0),
#endif
      +FP::max(), +FP::max()*(R(1.0)-eps), +FP::max()*(R(1.0)-R(2.0)*eps),
      -FP::max(), -FP::max()*(R(1.0)-eps), -FP::max()*(R(1.0)-R(2.0)*eps),
      +R(0.5)*FP::max(), +R(0.5)*FP::max()*(R(1.0)+eps),
      -R(0.5)*FP::max(), -R(0.5)*FP::max()*(R(1.0)+eps),
#ifdef VML_HAVE_INF
      +R(1.0/0.0),              // +FP::infinity()
      -R(1.0/0.0),              // -FP::infinity()
#endif
#ifdef VML_HAVE_NAN
      R(0.0/0.0),               // FP::quite_NaN()
#endif
      +int_min, +int_max, +uint_min, +uint_max,
      -int_min, -int_max, -uint_min, -uint_max,
      +int_min+R(0.1), +int_max+R(0.1), +uint_min+R(0.1), +uint_max+R(0.1),
      -int_min+R(0.1), -int_max+R(0.1), -uint_min+R(0.1), -uint_max+R(0.1),
      +int_min-R(0.1), +int_max-R(0.1), +uint_min-R(0.1), +uint_max-R(0.1),
      -int_min-R(0.1), -int_max-R(0.1), -uint_min-R(0.1), -uint_max-R(0.1),
      +int_min+R(1.0), +int_max+R(1.0), +uint_min+R(1.0), +uint_max+R(1.0),
      -int_min+R(1.0), -int_max+R(1.0), -uint_min+R(1.0), -uint_max+R(1.0),
      +int_min-R(1.0), +int_max-R(1.0), +uint_min-R(1.0), +uint_max-R(1.0),
      -int_min-R(1.0), -int_max-R(1.0), -uint_min-R(1.0), -uint_max-R(1.0),
      -R(443.9999425),
    };
    const int nvalues = sizeof values / sizeof *values;
    
    for (int i=0; i<8*nvalues+imax; ++i) {
      const realvec_t x =
        i<8*nvalues && i&1 ? RV(values[i/8]) : random(R(-10.0), R(+10.0));
      const realvec_t y =
        i<8*nvalues && i&2 ? RV(values[i/8]) : random(R(-10.0), R(+10.0));
      const realvec_t z =
        i<8*nvalues && i&4 ? RV(values[i/8]) : random(R(-10.0), R(+10.0));
      const intvec_t n = random(int_t(-10), int_t(+10));
      check_real("copysign",
                 (real_t(*)(real_t,real_t))std::copysign,
                 (realvec_t(*)(realvec_t,realvec_t))vecmathlib::copysign,
                 x, y, R(0.0));
      check_real("fabs",
                 (real_t(*)(real_t))std::fabs,
                 (realvec_t(*)(realvec_t))vecmathlib::fabs,
                 x, 0.0);
      check_real<realvec_t,realvec_t>("fdim", std::fdim, vecmathlib::fdim, x, y, accuracy());
      check_real<realvec_t,realvec_t,realvec_t>("fma", std::fma, vecmathlib::fma, x, y, z, R(2.0)*accuracy());
      check_real<realvec_t,realvec_t>("fmax", std::fmax, vecmathlib::fmax, x, y, 0.0);
      check_real<realvec_t,realvec_t>("fmin", std::fmin, vecmathlib::fmin, x, y, 0.0);
      check_int("ilogb", local_ilogb, (intvec_t(*)(realvec_t))vecmathlib::ilogb, x);
#if defined VML_HAVE_INF && defined VML_HAVE_NAN
      check_bool("isfinite", std::isfinite, vecmathlib::isfinite, x);
#endif
#ifdef VML_HAVE_INF
      check_bool("isinf", std::isinf, vecmathlib::isinf, x);
#endif
#ifdef VML_HAVE_NAN
      check_bool("isnan", std::isnan, vecmathlib::isnan, x);
#endif
#ifdef VML_HAVE_DENORMALS
      check_bool("isnormal", std::isnormal, vecmathlib::isnormal, x);
#endif
      check_real<realvec_t,intvec_t>("ldexp", local_ldexp, vecmathlib::ldexp, x, n, 0.0);
      check_bool("signbit",
                 (bool(*)(real_t))std::signbit,
                 (boolvec_t(*)(realvec_t))vecmathlib::signbit,
                 x);
    }
  }
  
  static void test_convert()
  {
    cout << "   testing ceil convert_float convert_int floor rint round trunc...\n"
         << flush;
    
    const real_t eps = FP::epsilon();
    const real_t int_min = R(std::numeric_limits<int_t>::min());
    const real_t int_max = R(std::numeric_limits<int_t>::max());
    const real_t uint_min = R(std::numeric_limits<uint_t>::min());
    const real_t uint_max = R(std::numeric_limits<uint_t>::max());
    const real_t mantissa_max = (U(1) << (FP::mantissa_bits+1)) - U(1);
    const real_t real_max =
      (((U(1) << (FP::mantissa_bits+1)) - U(1)) << (FP::exponent_bits-1)) +
      (U(1) << (FP::exponent_bits-1)) - U(1);
    const real_t values[] = {
      R(+0.0), R(+0.1), R(+0.9), R(+1.0), R(+1.1),
      R(-0.0), R(-0.1), R(-0.9), R(-1.0), R(-1.1),
      R(+0.0)+eps, R(+0.1)+eps, R(+0.9)+eps, R(+1.0)+eps, R(+1.1)+eps,
      R(-0.0)+eps, R(-0.1)+eps, R(-0.9)+eps, R(-1.0)+eps, R(-1.1)+eps,
      R(+0.0)-eps, R(+0.1)-eps, R(+0.9)-eps, R(+1.0)-eps, R(+1.1)-eps,
      R(-0.0)-eps, R(-0.1)-eps, R(-0.9)-eps, R(-1.0)-eps, R(-1.1)-eps,
#ifdef VML_HAVE_DENORMALS
      +FP::min(), +FP::min()*(R(1.0)+eps), +FP::min()*R(2.0),
      -FP::min(), -FP::min()*(R(1.0)+eps), -FP::min()*R(2.0),
#endif
      +FP::max(), +FP::max()*(R(1.0)-eps), +FP::max()*(R(1.0)-R(2.0)*eps),
      -FP::max(), -FP::max()*(R(1.0)-eps), -FP::max()*(R(1.0)-R(2.0)*eps),
      +R(0.5)*FP::max(), +R(0.5)*FP::max()*(R(1.0)+eps),
      -R(0.5)*FP::max(), -R(0.5)*FP::max()*(R(1.0)+eps),
#ifdef VML_HAVE_INF
      +R(1.0/0.0),              // +FP::infinity()
      -R(1.0/0.0),              // -FP::infinity()
#endif
#ifdef VML_HAVE_NAN
      R(0.0/0.0),               // FP::quite_NaN()
#endif
      +int_min, +int_max, +uint_min, +uint_max,
      -int_min, -int_max, -uint_min, -uint_max,
      +int_min+R(0.1), +int_max+R(0.1), +uint_min+R(0.1), +uint_max+R(0.1),
      -int_min+R(0.1), -int_max+R(0.1), -uint_min+R(0.1), -uint_max+R(0.1),
      +int_min-R(0.1), +int_max-R(0.1), +uint_min-R(0.1), +uint_max-R(0.1),
      -int_min-R(0.1), -int_max-R(0.1), -uint_min-R(0.1), -uint_max-R(0.1),
      +int_min+R(1.0), +int_max+R(1.0), +uint_min+R(1.0), +uint_max+R(1.0),
      -int_min+R(1.0), -int_max+R(1.0), -uint_min+R(1.0), -uint_max+R(1.0),
      +int_min-R(1.0), +int_max-R(1.0), +uint_min-R(1.0), +uint_max-R(1.0),
      -int_min-R(1.0), -int_max-R(1.0), -uint_min-R(1.0), -uint_max-R(1.0),
      +mantissa_max, +mantissa_max-R(1.0), +mantissa_max+R(1.0),
      -mantissa_max, -mantissa_max-R(1.0), -mantissa_max+R(1.0),
      +real_max, +real_max-R(1.0), +real_max+R(1.0),
      -real_max, -real_max-R(1.0), -real_max+R(1.0),
      -R(443.9999425),
    };
    const int nvalues = sizeof values / sizeof *values;
    
    for (int i=0; i<nvalues+imax; ++i) {
      const realvec_t x =
        i<nvalues ? RV(values[i]) : random(R(-1.0e+10), R(+1.0e+10));
      const intvec_t n1 = random(int_t(-100), int_t(+100));
      //const intvec_t n2 = random(int_t(-1000000000), int_t(+1000000000));
      const intvec_t n2 =
        random(std::numeric_limits<int_t>::min() / 2, // avoid overflow
               std::numeric_limits<int_t>::max() / 2);
      const realvec_t fn1 = vecmathlib::convert_float(n1);
      const realvec_t fn2 = vecmathlib::convert_float(n2);
      const realvec_t fn1h = vecmathlib::convert_float(n1) * RV(0.25);
      const realvec_t fn2h = vecmathlib::convert_float(n2) * RV(0.25);
      check_real<intvec_t>("convert_float",
            FP::convert_float, vecmathlib::convert_float, n1, accuracy());
      check_real<intvec_t>("convert_float",
            FP::convert_float, vecmathlib::convert_float, n2, accuracy());
      // Note: RV(int_max) > int_max due to rounding
      if (all(x >= RV(int_min) && x < RV(int_max))) {
        check_int("convert_int",
                  (int_t(*)(real_t))FP::convert_int,
                  (intvec_t(*)(realvec_t))vecmathlib::convert_int, x);
      }
      check_real<realvec_t>("ceil", std::ceil, vecmathlib::ceil, x, accuracy());
      check_real<realvec_t>("ceil", std::ceil, vecmathlib::ceil, fn1, accuracy());
      check_real<realvec_t>("ceil", std::ceil, vecmathlib::ceil, fn2, accuracy());
      check_real<realvec_t>("ceil", std::ceil, vecmathlib::ceil, fn1h, accuracy());
      check_real<realvec_t>("ceil", std::ceil, vecmathlib::ceil, fn2h, accuracy());
      check_real<realvec_t>("floor", std::floor, vecmathlib::floor, x, accuracy());
      check_real<realvec_t>("floor", std::floor, vecmathlib::floor, fn1, accuracy());
      check_real<realvec_t>("floor", std::floor, vecmathlib::floor, fn2, accuracy());
      check_real<realvec_t>("floor", std::floor, vecmathlib::floor, fn1h, accuracy());
      check_real<realvec_t>("floor", std::floor, vecmathlib::floor, fn2h, accuracy());
      check_real<realvec_t>("rint", std::rint, vecmathlib::rint, x, accuracy());
      check_real<realvec_t>("rint", std::rint, vecmathlib::rint, fn1, accuracy());
      check_real<realvec_t>("rint", std::rint, vecmathlib::rint, fn2, accuracy());
      check_real<realvec_t>("rint", std::rint, vecmathlib::rint, fn1h, accuracy());
      check_real<realvec_t>("rint", std::rint, vecmathlib::rint, fn2h, accuracy());
      check_real<realvec_t>("round", std::round, vecmathlib::round, x, accuracy());
      check_real<realvec_t>("round", std::round, vecmathlib::round, fn1, accuracy());
      check_real<realvec_t>("round", std::round, vecmathlib::round, fn2, accuracy());
      check_real<realvec_t>("round", std::round, vecmathlib::round, fn1h, accuracy());
      check_real<realvec_t>("round", std::round, vecmathlib::round, fn2h, accuracy());
      check_real<realvec_t>("trunc", std::trunc, vecmathlib::trunc, x, accuracy());
      check_real<realvec_t>("trunc", std::trunc, vecmathlib::trunc, fn1, accuracy());
      check_real<realvec_t>("trunc", std::trunc, vecmathlib::trunc, fn2, accuracy());
      check_real<realvec_t>("trunc", std::trunc, vecmathlib::trunc, fn1h, accuracy());
      check_real<realvec_t>("trunc", std::trunc, vecmathlib::trunc, fn2h, accuracy());
    }
  }
  
  
  
  static void test_asin()
  {
    cout << "   testing asin acos atan...\n" << flush;
    for (int i=0; i<imax; ++i) {
      const realvec_t x = random(R(-1.0), R(+1.0));
      check_real<realvec_t>("asin", std::asin, vecmathlib::asin, x, accuracy());
      check_real<realvec_t>("acos", std::acos, vecmathlib::acos, x, accuracy());
    }
    for (int i=0; i<imax; ++i) {
      const realvec_t x = random(R(-100.0), R(+100.0));
      check_real<realvec_t>("atan", std::atan, vecmathlib::atan, x, accuracy());
    }
  }
  
  static void test_asinh()
  {
    cout << "   testing asinh acosh atanh...\n" << flush;
    for (int i=0; i<imax; ++i) {
      const realvec_t x = random(R(-1000.0), R(+1000.0));
      check_real<realvec_t>("asinh", std::asinh, vecmathlib::asinh, x, accuracy());
    }
    for (int i=0; i<imax; ++i) {
      const realvec_t x = random(R(1.0), R(1000.0));
      check_real<realvec_t>("acosh", std::acosh, vecmathlib::acosh, x, accuracy());
    }
    for (int i=0; i<imax; ++i) {
      const realvec_t x = random(R(-1.0), R(+1.0));
      check_real<realvec_t>("atanh", std::atanh, vecmathlib::atanh, x, accuracy());
    }
  }
  
  static real_t local_exp10(real_t x) { return pow(R(10.0), x); }
  static void test_exp()
  {
    cout << "   testing exp exp10 exp2 expm1...\n" << flush;
    for (int i=0; i<imax; ++i) {
      const realvec_t x = random(R(-100.0), R(+100.0));
      check_real<realvec_t>("exp", std::exp, vecmathlib::exp, x, accuracy());
      check_real<realvec_t>("exp10", local_exp10, vecmathlib::exp10, x, accuracy());
      check_real<realvec_t>("exp2", std::exp2, vecmathlib::exp2, x, accuracy());
      check_real<realvec_t>("expm1", std::expm1, vecmathlib::expm1, x, accuracy());
    }
  }
  
  static void test_log()
  {
    cout << "   testing log log10 log1p log2...\n" << flush;
    for (int i=0; i<imax; ++i) {
      const realvec_t x = random(R(1.0e-10), R(1.0e+10));
      check_real<realvec_t>("log", std::log, vecmathlib::log, x, accuracy());
      check_real<realvec_t>("log10", std::log10, vecmathlib::log10, x, accuracy());
      check_real<realvec_t>("log1p", std::log1p, vecmathlib::log1p, x, accuracy());
      check_real<realvec_t>("log2", std::log2, vecmathlib::log2, x, accuracy());
    }
  }
  
  static void test_pow()
  {
    cout << "   testing pow...\n" << flush;
    for (int i=0; i<imax; ++i) {
      const realvec_t x = random(R(0.001), R(1000.0));
      const realvec_t y = random(R(-10.0), R(+10.0));
      const realvec_t ya = fabs(y);
      const intvec_t n = random(I(-10), I(+10));
      const realvec_t fn = vecmathlib::convert_float(n);
      check_real<realvec_t,realvec_t>("pow(0,y)", std::pow, vecmathlib::pow, RV(0.0), ya, accuracy());
      check_real<realvec_t,realvec_t>("pow(x,0)", std::pow, vecmathlib::pow, x, RV(0.0), accuracy());
      // just to check
      check_real<realvec_t>("log(x)", std::log, vecmathlib::log, x, accuracy());
      check_real<realvec_t,realvec_t>("pow(x,y)", std::pow, vecmathlib::pow, x, y, accuracy());
      check_real<realvec_t,realvec_t>("pow(-x,n)", std::pow, vecmathlib::pow, -x, fn, accuracy());
    }
  }
  
  static real_t local_rcp(real_t x) { return R(1.0)/x; }
  static void test_rcp()
  {
    cout << "   testing fmod rcp remainder...\n" << flush;
    for (int i=0; i<imax; ++i) {
      const realvec_t x = random(R(-100.0), R(+100.0));
      const realvec_t y = random(R(-100.0), R(+100.0));
      const intvec_t n = random(I(-100), I(+100));
      const intvec_t m = random(I(-100), I(+100));
      const realvec_t fn = vecmathlib::convert_float(n);
      const realvec_t fm = vecmathlib::convert_float
        (m + vecmathlib::convert_int(m == intvec_t(I(0))));
      check_real<realvec_t>("rcp", local_rcp, vecmathlib::rcp, x, accuracy());
      check_real<realvec_t,realvec_t>("fmod(x,y)", std::fmod, vecmathlib::fmod, x, y, 2.0*accuracy());
      check_real<realvec_t,realvec_t>("fmod(x,m)", std::fmod, vecmathlib::fmod, x, fm, 2.0*accuracy());
      check_real<realvec_t,realvec_t>("fmod(n,y)", std::fmod, vecmathlib::fmod, fn, y, 2.0*accuracy());
      check_real<realvec_t,realvec_t>("remainder(x,y)",
				 std::remainder, vecmathlib::remainder, x, y, R(2.0)*accuracy());
      check_real<realvec_t,realvec_t>("remainder(x,m)",
				 std::remainder, vecmathlib::remainder, x, fm, R(2.0)*accuracy());
      check_real<realvec_t,realvec_t>("remainder(n,y)",
				 std::remainder, vecmathlib::remainder, fn, y, R(2.0)*accuracy());
    }
  }
  
  static void test_sin()
  {
    cout << "   testing cos sin tan...\n" << flush;
    for (int i=0; i<imax; ++i) {
      const realvec_t x = random(R(-10.0), R(+10.0));
      check_real<realvec_t>("sin", std::sin, vecmathlib::sin, x, accuracy());
      check_real<realvec_t>("cos", std::cos, vecmathlib::cos, x, accuracy());
    }
    for (int i=0; i<imax; ++i) {
      const realvec_t x0 = random(R(-1.55), R(+1.55));
      const intvec_t n = random(I(-10), I(+10));
      const realvec_t x = x0 + vecmathlib::convert_float(n) * RV(M_PI);
      // tan loses accuracy near pi/2
      // (by definition, not by implementation?)
      check_real<realvec_t>("tan", std::tan, vecmathlib::tan, x, R(100.0)*accuracy());
    }
  }
  
  static void test_sinh()
  {
    cout << "   testing cosh sinh tanh...\n" << flush;
    for (int i=0; i<imax; ++i) {
      const realvec_t x = random(R(-10.0), R(+10.0));
      check_real<realvec_t>("sinh", std::sinh, vecmathlib::sinh, x, accuracy());
      check_real<realvec_t>("cosh", std::cosh, vecmathlib::cosh, x, accuracy());
      check_real<realvec_t>("tanh", std::tanh, vecmathlib::tanh, x, accuracy());
    }
  }
  
  static real_t local_rsqrt(real_t x) { return R(1.0)/sqrt(x); }
  static void test_sqrt()
  {
    cout << "   testing cbrt hypot rsqrt sqrt...\n" << flush;
    for (int i=0; i<imax; ++i) {
      const realvec_t x = random(R(0.0), R(1.0e+3));
      const realvec_t y = random(-R(1.0e+3), R(1.0e+3));
      const realvec_t z = random(-R(1.0e+3), R(1.0e+3));
      check_real<realvec_t>("cbrt", std::cbrt, vecmathlib::cbrt, x, accuracy());
      check_real<realvec_t,realvec_t>("hypot", std::hypot, vecmathlib::hypot, y, z, accuracy());
      check_real<realvec_t>("rsqrt", local_rsqrt, vecmathlib::rsqrt, x, accuracy());
      check_real<realvec_t>("sqrt", std::sqrt, vecmathlib::sqrt, x, accuracy());
    }
  }
  
  
  
  static void test()
  {
    cout << "\n"
         << "Testing math functions for type " << realvec_t::name() << ":\n";
    
    test_bool();
    test_int();
    // test_float();
    
    test_mem();
    
    // Test "basic" functions first
    test_fabs();
    test_convert();
#warning "TODO"
    // test_rcp();
    // test_sqrt();
    // test_exp();
    // test_log();
    // test_pow();
    // test_sin();
    // test_sinh();
    // test_asin();
    // test_asinh();
  }
};



int main(int argc, char** argv)
{
  using namespace vecmathlib;

  cout << "Testing math functions:\n" << flush;
  
  cout << "Vecmathlib configuration: [conf"
#ifdef VML_DEBUG
    "-DEBUG"
#endif
#ifdef __ALTIVEC__
    "-Altivec"
#endif
#ifdef __VSX__
    "-VSX"
#endif
#if defined __bgq__ && defined __VECTOR4DOUBLE__
    "-QPX"
#endif
#ifdef __SSE2__
    "-SSE2"
#endif
#ifdef __SSE3__
    "-SSE3"
#endif
#ifdef __SSE4_1__
    "-SSE4.1"
#endif
#ifdef __SSE4a__
    "-SSE4a"
#endif
#ifdef __AVX__
    "-AVX"
#endif
    "]\n";
  
#warning "TODO"
//   vecmathlib_test<realpseudovec<float,1>>::test();
//   // vecmathlib_test<realbuiltinvec<float,1>>::test();
//   vecmathlib_test<realtestvec<float,1>>::test();
// #ifdef VECMATHLIB_HAVE_VEC_FLOAT_1
//   vecmathlib_test<realvec<float,1>>::test();
// #endif
//   vecmathlib_test<realpseudovec<float,4>>::test();
//   // vecmathlib_test<realbuiltinvec<float,4>>::test();
//   vecmathlib_test<realtestvec<float,4>>::test();
// #ifdef VECMATHLIB_HAVE_VEC_FLOAT_4
//   vecmathlib_test<realvec<float,4>>::test();
// #endif
// #ifdef VECMATHLIB_HAVE_VEC_FLOAT_8
//   vecmathlib_test<realpseudovec<float,8>>::test();
//   // vecmathlib_test<realbuiltinvec<float,8>>::test();
//   vecmathlib_test<realtestvec<float,8>>::test();
//   vecmathlib_test<realvec<float,8>>::test();
// #endif
//   
//   vecmathlib_test<realpseudovec<double,1>>::test();
//   // vecmathlib_test<realbuiltinvec<double,1>>::test();
//   vecmathlib_test<realtestvec<double,1>>::test();
// #ifdef VECMATHLIB_HAVE_VEC_DOUBLE_1
//   vecmathlib_test<realvec<double,1>>::test();
// #endif
//   vecmathlib_test<realpseudovec<double,2>>::test();
//   // vecmathlib_test<realbuiltinvec<double,2>>::test();
//   vecmathlib_test<realtestvec<double,2>>::test();
// #ifdef VECMATHLIB_HAVE_VEC_DOUBLE_2
//   vecmathlib_test<realvec<double,2>>::test();
// #endif
#ifdef VECMATHLIB_HAVE_VEC_DOUBLE_4
//   vecmathlib_test<realpseudovec<double,4>>::test();
//   // vecmathlib_test<realbuiltinvec<double,4>>::test();
//   vecmathlib_test<realtestvec<double,4>>::test();
  vecmathlib_test<realvec<double,4>>::test();
#endif
  
  cout << "\n";
  if (num_errors == 0) {
    cout << "SUCCESS";
  } else {
    cout << "FAILURE";
  }
  cout << ": " << num_errors << " errors found\n" << flush;
  
  return num_errors == 0 ? 0 : 1;
}
