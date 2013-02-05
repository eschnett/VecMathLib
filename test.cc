// -*-C++-*-

#include "vecmathlib.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <typeinfo>

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
  
  
  
  // Test each function with this many random values
  static int const imax = 10000;
  // Require (arbitrarily) that 3/4 of the digits are correct
  static real_t accuracy() { return pow(realvec_t::epsilon(), R(0.75)); }
  
  
  
  static realvec_t random(real_t const xmin, real_t const xmax)
  {
    realvec_t x;
    for (int i=0; i<realvec_t::size; ++i) {
      real_t r =
        (xmax - xmin) *FP::convert_float(rand()) / FP::convert_float(RAND_MAX);
      x.set_elt(i, xmin + r);
    }
    return x;
  }
  
  static intvec_t random(int_t const nmin, int_t const nmax)
  {
    intvec_t n;
    for (int i=0; i<intvec_t::size; ++i) {
      real_t r =
        R(nmax - nmin + 1) *
        R(rand()) / (R(RAND_MAX) + R(1.0));
      n.set_elt(i, nmin + FP::convert_int(std::floor(r)));
    }
    return n;
  }
  
  
  
  static void check_mem(char const* const func,
                        real_t const* p,
                        realvec_t x,
                        realvec_t xorig,
                        int mval)
  {
    realvec_t y;
    for (int i=0; i<realvec_t::size; ++i) {
      y.set_elt(i, mval & (1<<i) ? p[i] : xorig[i]);
    }
    boolvec_t isbad = x != y;
    if (any(isbad)) {
      ++ num_errors;
      cout << setprecision(realvec_t::digits10+2)
           << "Error in " << func << ":\n"
           << "   found=" << x << "\n"
           << "   expected=" << y << "\n"
           << "   isbad=" << isbad << "\n"
           << flush;
    }
  }
  
  static void check_mem(char const* const func,
                        real_t const* p,
                        realvec_t x,
                        real_t const* porig,
                        int mval)
  {
    realvec_t pvec, y;
    for (int i=0; i<realvec_t::size; ++i) {
      pvec.set_elt(i, p[i]);
      y.set_elt(i, mval & (1<<i) ? x[i] : porig[i]);
    }
    boolvec_t isbad = pvec != y;
    if (any(isbad)) {
      ++ num_errors;
      cout << setprecision(realvec_t::digits10+2)
           << "Error in " << func << ":\n"
           << "   found=" << pvec << "\n"
           << "   expected=" << y << "\n"
           << "   isbad=" << isbad << "\n"
           << flush;
    }
  }
  
  template<typename A>
  static void check(char const* const func,
                    real_t fstd(typename A::scalar_t), realvec_t fvml(A),
                    A const x,
                    real_t const accuracy)
  {
    realvec_t rstd;
    for (int i=0; i<realvec_t::size; ++i) {
      rstd.set_elt(i, fstd(x[i]));
    }
    realvec_t const rvml = fvml(x);
    realvec_t const dr = rstd - rvml;
    realvec_t const scale = fabs(rstd) + fabs(rvml) + realvec_t(1.0);
    boolvec_t const isbad = fabs(dr) > realvec_t(accuracy) * scale;
    if (any(isbad)) {
      ++ num_errors;
      cout << setprecision(realvec_t::digits10+2)
           << "Error in " << func << "(" << x << "):\n"
           << "   fstd(x)=" << rstd << "\n"
           << "   fvml(x)=" << rvml << "\n"
           << "   rel-error(x)=" << fabs(dr) / scale << "\n"
           << "   isbad(x)=" << isbad << "\n"
           << flush;
    }
  }
  
  template<typename A, typename B>
  static void check(char const* const func,
                    real_t fstd(typename A::scalar_t, typename B::scalar_t),
                    realvec_t fvml(A, B),
                    A const x, B const y,
                    real_t const accuracy)
  {
    realvec_t rstd;
    for (int i=0; i<realvec_t::size; ++i) {
      rstd.set_elt(i, fstd(x[i], y[i]));
    }
    realvec_t const rvml = fvml(x, y);
    realvec_t const dr = rstd - rvml;
    realvec_t const scale = fabs(rstd) + fabs(rvml) + realvec_t(1.0);
    boolvec_t const isbad = fabs(dr) > realvec_t(accuracy) * scale;
    if (any(isbad)) {
      ++ num_errors;
      cout << setprecision(realvec_t::digits10+1)
           << "Error in " << func << "(" << x << "," << y << "):\n"
           << "   fstd(x,y)=" << rstd << "\n"
           << "   fvml(x,y)=" << rvml << "\n"
           << "   rel-error(x,y)=" << fabs(dr) / scale << "\n"
           << "   isbad(x,y)=" << isbad << "\n"
           << flush;
    }
  }
  
  template<typename A, typename B, typename C>
  static void check(char const* const func,
                    real_t fstd(typename A::scalar_t, typename B::scalar_t,
                                typename C::scalar_t),
                    realvec_t fvml(A, B, C),
                    A const x, B const y, C const z,
                    real_t const accuracy)
  {
    realvec_t rstd;
    for (int i=0; i<realvec_t::size; ++i) {
      rstd.set_elt(i, fstd(x[i], y[i], z[i]));
    }
    realvec_t const rvml = fvml(x, y, z);
    realvec_t const dr = rstd - rvml;
    if (any(fabs(dr) >
            realvec_t(accuracy) * (fabs(rstd) + fabs(rvml) + realvec_t(1.0))))
    {
      ++ num_errors;
      cout << setprecision(realvec_t::digits10+1)
           << "Error in " << func << "(" << x << "," << y<< "," << z << "):\n"
           << "   fstd(x,y,z)=" << rstd << "\n"
           << "   fvml(x,y,z)=" << rvml << "\n"
           << flush;
    }
  }
  
  template<typename A>
  static void check(char const* const func,
                    int_t fstd(typename A::scalar_t), intvec_t fvml(A),
                    A const x)
  {
    intvec_t rstd;
    for (int i=0; i<intvec_t::size; ++i) {
      rstd.set_elt(i, fstd(x[i]));
    }
    intvec_t const rvml = fvml(x);
    intvec_t const dr = rstd - rvml;
    if (any(convert_bool(dr))) {
      cout << setprecision(realvec_t::digits10+1)
           << "Error in " << func << "(" << x << "):\n"
           << "   fstd(x)=" << rstd << "\n"
           << "   fvml(x)=" << rvml << "\n"
           << flush;
    }
  }
  
  template<typename A, typename B>
  static void check(char const* const func,
                    int fstd(typename A::scalar_t, typename B::scalar_t),
                    intvec_t fvml(A, B),
                    A const x, B const y)
  {
    intvec_t rstd;
    for (int i=0; i<intvec_t::size; ++i) {
      rstd.set_elt(i, fstd(x[i], y[i]));
    }
    intvec_t const rvml = fvml(x, y);
    intvec_t const dr = rstd - rvml;
    if (any(convert_bool(dr))) {
      ++ num_errors;
      cout << setprecision(realvec_t::digits10+1)
           << "Error in " << func << "(" << x << "," << y << "):\n"
           << "   fstd(x,y)=" << rstd << "\n"
           << "   fvml(x,y)=" << rvml << "\n"
           << flush;
    }
  }
  
  template<typename A>
  static void check(char const* const func,
                    bool fstd(typename A::scalar_t), boolvec_t fvml(A),
                    A const x)
  {
    boolvec_t rstd;
    for (int i=0; i<boolvec_t::size; ++i) {
      rstd.set_elt(i, fstd(x[i]));
    }
    boolvec_t const rvml = fvml(x);
    boolvec_t const dr = rstd != rvml;
    if (any(dr)) {
      cout << setprecision(realvec_t::digits10+1)
           << "Error in " << func << "(" << x << "):\n"
           << "   fstd(x)=" << rstd << "\n"
           << "   fvml(x)=" << rvml << "\n"
           << flush;
    }
  }
  
  
  
  static void test_mem()
  {
    cout << "   testing loada loadu storea storeu (errors may lead to segfaults)...\n" << flush;
    int const n = 6;
    realvec_t x[n], xnew[n];
    for (int i=0; i<n; ++i) {
      x[i] = random(R(-10.0), R(+10.0));
    }
    realvec_t const z = random(R(-10.0), R(+10.0));
    
    // loada
    {
      real_t *p = (real_t*)&x[1];
      realvec_t y = realvec_t::loada(p);
      check_mem("loada", p, y, z, ~0);
    }
    
    // loadu
    for (ptrdiff_t i=0; i<realvec_t::size; ++i) {
      real_t *p = (real_t*)&x[1];
      realvec_t y = realvec_t::loadu(p+i);
      check_mem("loadu", p+i, y, z, ~0);
    }
    
    // loadu(ioff)
    for (ptrdiff_t ioff=0; ioff<realvec_t::size; ++ioff) {
      real_t *p = (real_t*)&x[1];
      realvec_t y = realvec_t::loadu(p, ioff);
      check_mem("loadu(ioff)", p+ioff, y, z, ~0);
    }
    
    // storea
    {
      memcpy(xnew, x, n*sizeof *xnew);
      real_t *p = (real_t*)&xnew[1];
      storea(z, p);
      check_mem("storea", p, z, (real_t*)&x[1], ~0);
    }
    
    // storeu
    for (ptrdiff_t i=0; i<realvec_t::size; ++i) {
      memcpy(xnew, x, n*sizeof *xnew);
      real_t *p = (real_t*)&xnew[1];
      storeu(z, p+i);
      check_mem("storeu", p+i, z, (real_t*)&x[1]+i, ~0);
    }
    
    // storeu
    for (ptrdiff_t ioff=0; ioff<realvec_t::size; ++ioff) {
      memcpy(xnew, x, n*sizeof *xnew);
      real_t *p = (real_t*)&xnew[1];
      storeu(z, p, ioff);
      check_mem("storeu(ioff)", p+ioff, z, (real_t*)&x[1]+ioff, ~0);
    }
    
    for (int mval=0; mval<(1<<realvec_t::size); ++mval) {
      boolvec_t mbool;
      for (int i=0; i<realvec_t::size; ++i) mbool.set_elt(i, mval & (1<<i));
      typename realvec_t::mask_t mask(mbool);
      
      // loada(mask)
      {
        real_t *p = (real_t*)&x[1];
        realvec_t y = loada(p, z, mask);
        check_mem("loada(mask)", p, y, z, mval);
      }
      
      // loadu(mask)
      for (ptrdiff_t i=0; i<realvec_t::size; ++i) {
        real_t *p = (real_t*)&x[1];
        realvec_t y = loadu(p+i, z, mask);
        check_mem("loadu(mask)", p+i, y, z, mval);
      }
      
      // loadu(ioff, mask)
      for (ptrdiff_t ioff=0; ioff<realvec_t::size; ++ioff) {
        real_t *p = (real_t*)&x[1];
        realvec_t y = loadu(p, ioff, z, mask);
        check_mem("loadu(ioff,mask)", p+ioff, y, z, mval);
      }
      
      // storea
      {
        memcpy(xnew, x, n*sizeof *xnew);
        real_t *p = (real_t*)&xnew[1];
        storea(z, p, mask);
        check_mem("storea(mask)", p, z, (real_t*)&x[1], mval);
      }
      
      // storeu
      for (ptrdiff_t i=0; i<realvec_t::size; ++i) {
        memcpy(xnew, x, n*sizeof *xnew);
        real_t *p = (real_t*)&xnew[1];
        storeu(z, p+i, mask);
        check_mem("storeu(mask)", p+i, z, (real_t*)&x[1]+i, mval);
      }
      
      // storeu
      for (ptrdiff_t ioff=0; ioff<realvec_t::size; ++ioff) {
        memcpy(xnew, x, n*sizeof *xnew);
        real_t *p = (real_t*)&xnew[1];
        storeu(z, p, ioff, mask);
        check_mem("storeu(ioff,mask)", p+ioff, z, (real_t*)&x[1]+ioff, mval);
      }
      
    } // for mval
  }
  
  // Change signature: "int" -> "int_t"
  static int_t ilogb(real_t x) { return std::ilogb(x); }
  static real_t scalbn(real_t x, int_t n) { return std::scalbn(x, n); }
  static void test_fabs()
  {
    cout << "   testing copysign fabs fdim fma fmax fmin ilogb isfinite isinf isnan isnormal scalbn signbit...\n" << flush;
    for (int i=0; i<imax; ++i) {
      realvec_t const x = random(R(-10.0), R(+10.0));
      realvec_t const y = random(R(-10.0), R(+10.0));
      realvec_t const z = random(R(-10.0), R(+10.0));
      intvec_t const n = random(int_t(-10), int_t(+10));
      check("copysign", copysign, vecmathlib::copysign, x, y, 0.0);
      check("fabs", fabs, vecmathlib::fabs, x, 0.0);
      check("fdim", fdim, vecmathlib::fdim, x, y, accuracy());
      check("fma", fma, vecmathlib::fma, x, y, z, accuracy());
      check("fmax", fmax, vecmathlib::fmax, x, y, 0.0);
      check("fmin", fmin, vecmathlib::fmin, x, y, 0.0);
      check("ilogb", ilogb, vecmathlib::ilogb, x);
      check("isfinite", isfinite, vecmathlib::isfinite, x);
      check("isinf", isinf, vecmathlib::isinf, x);
      check("isnan", isnan, vecmathlib::isnan, x);
      check("isnormal", isnormal, vecmathlib::isnormal, x);
      check("scalbn", scalbn, vecmathlib::scalbn, x, n, 0.0);
      check("signbit", signbit, vecmathlib::signbit, x);
    }
  }
  
  static void test_convert()
  {
    cout << "   testing ceil convert_float convert_int floor round...\n"
         << flush;
    
    for (int i=0; i<imax; ++i) {
      realvec_t const x = random(R(-1.0e+10), R(+1.0e+10));
      intvec_t const n = random(int_t(-1000000), int_t(+1000000));
      realvec_t const fn = vecmathlib::convert_float(n);
      check("convert_float",
            FP::convert_float, vecmathlib::convert_float, n, accuracy());
      check("convert_int", FP::convert_int, vecmathlib::convert_int, x);
      check("ceil", ceil, vecmathlib::ceil, x, accuracy());
      check("ceil", ceil, vecmathlib::ceil, fn, accuracy());
      check("floor", floor, vecmathlib::floor, x, accuracy());
      check("floor", floor, vecmathlib::floor, fn, accuracy());
      check("round", round, vecmathlib::round, x, accuracy());
      check("round", round, vecmathlib::round, fn, accuracy());
    }
  }
  
  
  
  static void test_asin()
  {
    cout << "   testing asin acos atan...\n" << flush;
    for (int i=0; i<imax; ++i) {
      realvec_t const x = random(R(-1.0), R(+1.0));
      check("asin", asin, vecmathlib::asin, x, accuracy());
      check("acos", acos, vecmathlib::acos, x, accuracy());
    }
    for (int i=0; i<imax; ++i) {
      realvec_t const x = random(R(-100.0), R(+100.0));
      check("atan", atan, vecmathlib::atan, x, accuracy());
    }
  }
  
  static void test_asinh()
  {
    cout << "   testing asinh acosh atanh...\n" << flush;
    for (int i=0; i<imax; ++i) {
      realvec_t const x = random(R(-1000.0), R(+1000.0));
      check("asinh", asinh, vecmathlib::asinh, x, accuracy());
    }
    for (int i=0; i<imax; ++i) {
      realvec_t const x = random(R(1.0), R(1000.0));
      check("acosh", acosh, vecmathlib::acosh, x, accuracy());
    }
    for (int i=0; i<imax; ++i) {
      realvec_t const x = random(R(-1.0), R(+1.0));
      check("atanh", atanh, vecmathlib::atanh, x, accuracy());
    }
  }
  
  static real_t exp10(real_t x) { return pow(R(10.0), x); }
  static void test_exp()
  {
    cout << "   testing exp exp10 exp2 expm1...\n" << flush;
    for (int i=0; i<imax; ++i) {
      realvec_t const x = random(R(-100.0), R(+100.0));
      check("exp", exp, vecmathlib::exp, x, accuracy());
      check("exp10", exp10, vecmathlib::exp10, x, accuracy());
      check("exp2", exp2, vecmathlib::exp2, x, accuracy());
      check("expm1", expm1, vecmathlib::expm1, x, accuracy());
    }
  }
  
  static void test_log()
  {
    cout << "   testing log log10 log1p log2...\n" << flush;
    for (int i=0; i<imax; ++i) {
      realvec_t const x = random(R(1.0e-10), R(1.0e+10));
      check("log", log, vecmathlib::log, x, accuracy());
      check("log10", log10, vecmathlib::log10, x, accuracy());
      check("log1p", log1p, vecmathlib::log1p, x, accuracy());
      check("log2", log2, vecmathlib::log2, x, accuracy());
    }
  }
  
  static void test_pow()
  {
    cout << "   testing pow...\n" << flush;
    for (int i=0; i<imax; ++i) {
      realvec_t const x = random(R(0.001), R(1000.0));
      realvec_t const y = random(R(-10.0), R(+10.0));
      realvec_t const ya = fabs(y);
      intvec_t const n = random(I(-10), I(+10));
      realvec_t const fn = vecmathlib::convert_float(n);
      check("pow(0,y)", pow, vecmathlib::pow, RV(0.0), ya, accuracy());
      check("pow(x,0)", pow, vecmathlib::pow, x, RV(0.0), accuracy());
      // just to check
      check("log(x)", log, vecmathlib::log, x, accuracy());
      check("pow(x,y)", pow, vecmathlib::pow, x, y, accuracy());
      check("pow(-x,n)", pow, vecmathlib::pow, -x, fn, accuracy());
    }
  }
  
  static real_t rcp(real_t x) { return R(1.0)/x; }
  static void test_rcp()
  {
    cout << "   testing fmod rcp remainder...\n" << flush;
    for (int i=0; i<imax; ++i) {
      realvec_t const x = random(R(-100.0), R(+100.0));
      realvec_t const y = random(R(-100.0), R(+100.0));
      intvec_t const n = random(I(-100), I(+100));
      intvec_t const m = random(I(-100), I(+100));
      realvec_t const fn = vecmathlib::convert_float(n);
      realvec_t const fm = vecmathlib::convert_float
        (m + vecmathlib::convert_int(m == intvec_t(I(0))));
      check("rcp", rcp, vecmathlib::rcp, x, accuracy());
      check("fmod(x,y)", fmod, vecmathlib::fmod, x, y, accuracy());
      check("fmod(x,m)", fmod, vecmathlib::fmod, x, fm, accuracy());
      check("fmod(n,y)", fmod, vecmathlib::fmod, fn, y, accuracy());
      check("remainder(x,y)",
            remainder, vecmathlib::remainder, x, y, accuracy());
      check("remainder(x,m)",
            remainder, vecmathlib::remainder, x, fm, accuracy());
      check("remainder(n,y)",
            remainder, vecmathlib::remainder, fn, y, accuracy());
    }
  }
  
  static void test_sin()
  {
    cout << "   testing cos sin tan...\n" << flush;
    for (int i=0; i<imax; ++i) {
      realvec_t const x = random(R(-10.0), R(+10.0));
      check("sin", sin, vecmathlib::sin, x, accuracy());
      check("cos", cos, vecmathlib::cos, x, accuracy());
    }
    for (int i=0; i<imax; ++i) {
      realvec_t const x0 = random(R(-1.55), R(+1.55));
      intvec_t const n = random(I(-10), I(+10));
      realvec_t const x = x0 + vecmathlib::convert_float(n) * RV(M_PI);
      // tan loses accuracy near pi/2
      // (by definition, not by implementation?)
      check("tan", tan, vecmathlib::tan, x, R(100.0)*accuracy());
    }
  }
  
  static void test_sinh()
  {
    cout << "   testing cosh sinh tanh...\n" << flush;
    for (int i=0; i<imax; ++i) {
      realvec_t const x = random(R(-10.0), R(+10.0));
      check("sinh", sinh, vecmathlib::sinh, x, accuracy());
      check("cosh", cosh, vecmathlib::cosh, x, accuracy());
      check("tanh", tanh, vecmathlib::tanh, x, accuracy());
    }
  }
  
  static real_t rsqrt(real_t x) { return R(1.0)/sqrt(x); }
  static void test_sqrt()
  {
    cout << "   testing rsqrt sqrt...\n" << flush;
    for (int i=0; i<imax; ++i) {
      realvec_t const x = random(R(0.0), R(1.0e+3));
      check("rsqrt", rsqrt, vecmathlib::rsqrt, x, accuracy());
      check("sqrt", sqrt, vecmathlib::sqrt, x, accuracy());
    }
  }
  
  
  
  static void test()
  {
    cout << "\n"
         << "Testing math functions for type " << realvec_t::name() << ":\n";
    
    test_mem();
    
    test_fabs();
    test_convert();
    
    test_asin();
    test_asinh();
    test_exp();
    test_log();
    test_pow();
    test_rcp();
    test_sin();
    test_sinh();
    test_sqrt();
  }
};



int main(int argc, char** argv)
{
  using namespace vecmathlib;

  cout << "Testing math functions:\n" << flush;
  
#ifdef VECMATHLIB_HAVE_VEC_FLOAT_1
  vecmathlib_test<realpseudovec<float,1>>::test();
  vecmathlib_test<realvec<float,1>>::test();
#endif
#ifdef VECMATHLIB_HAVE_VEC_FLOAT_4
  vecmathlib_test<realpseudovec<float,4>>::test();
  vecmathlib_test<realvec<float,4>>::test();
#endif
#ifdef VECMATHLIB_HAVE_VEC_FLOAT_8
  vecmathlib_test<realpseudovec<float,8>>::test();
  vecmathlib_test<realvec<float,8>>::test();
#endif
  
#ifdef VECMATHLIB_HAVE_VEC_DOUBLE_1
  vecmathlib_test<realpseudovec<double,1>>::test();
  vecmathlib_test<realvec<double,1>>::test();
#endif
#ifdef VECMATHLIB_HAVE_VEC_DOUBLE_2
  vecmathlib_test<realpseudovec<double,2>>::test();
  vecmathlib_test<realvec<double,2>>::test();
#endif
#ifdef VECMATHLIB_HAVE_VEC_DOUBLE_4
  vecmathlib_test<realpseudovec<double,4>>::test();
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
