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
  typedef typename realvec_t::real_t real_t;
  
  typedef vecmathlib::floatprops<real_t> FP;
  
  
  
  static int const imax = 1000;
  static real_t constexpr accuracy = pow(realvec_t::epsilon(), real_t(0.75));
  
  
  
  static realvec_t random(real_t const xmin, real_t const xmax)
  {
    realvec_t x;
    for (int i=0; i<realvec_t::size; ++i) {
      x.set_elt(i, xmin + (xmax - xmin) * real_t(rand()) / real_t(RAND_MAX));
    }
    return x;
  }
  
  static intvec_t random(int_t const nmin, int_t const nmax)
  {
    intvec_t n;
    for (int i=0; i<intvec_t::size; ++i) {
      real_t x =
        real_t(nmax - nmin + 1) *
        real_t(rand()) / (real_t(RAND_MAX) + real_t(1.0));
      n.set_elt(i, nmin + FP::convert_int(std::floor(x)));
    }
    return n;
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
    if (any(fabs(dr) >
            realvec_t(accuracy) * (fabs(rstd) + fabs(rvml) + realvec_t(1.0))))
    {
      ++ num_errors;
      cout << setprecision(realvec_t::digits10+1)
           << "Error in " << func << "(" << x << "):\n"
           << "   fstd(x)=" << rstd << "\n"
           << "   fvml(x)=" << rvml << "\n";
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
    if (any(fabs(dr) >
            realvec_t(accuracy) * (fabs(rstd) + fabs(rvml) + realvec_t(1.0))))
    {
      ++ num_errors;
      cout << setprecision(realvec_t::digits10+1)
           << "Error in " << func << "(" << x << "," << y << "):\n"
           << "   fstd(x,y)=" << rstd << "\n"
           << "   fvml(x,y)=" << rvml << "\n";
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
           << "   fvml(x)=" << rvml << "\n";
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
           << "   fvml(x,y)=" << rvml << "\n";
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
           << "   fvml(x)=" << rvml << "\n";
    }
  }
  
  
  
  static int_t ilogb(real_t x) { return std::ilogb(x); }
  static real_t scalbn(real_t x, int_t n) { return std::scalbn(x, n); }
  static void test_fabs()
  {
    cout << "   testing copysign fabs ilogb scalbn signbit...\n";
    for (int i=0; i<imax; ++i) {
      realvec_t const x = random(real_t(-10.0), real_t(+10.0));
      realvec_t const y = random(real_t(-10.0), real_t(+10.0));
      intvec_t const n = random(int_t(-10), int_t(+10));
      check("copysign", copysign, vecmathlib::copysign, x, y, 0.0);
      check("fabs", fabs, vecmathlib::fabs, x, 0.0);
      check("ilogb", ilogb, vecmathlib::ilogb, x);
      check("scalbn", scalbn, vecmathlib::scalbn, x, n, 0.0);
      check("signbit", signbit, vecmathlib::signbit, x);
    }
  }
  
  static void test_convert()
  {
    cout << "   testing ceil convert_float convert_int floor round...\n";
    for (int i=0; i<imax; ++i) {
      realvec_t const x = random(real_t(-10.0), real_t(+10.0));
      intvec_t const n = random(int_t(-10), int_t(+10));
      check("ceil", ceil, vecmathlib::ceil, x, accuracy);
      check("convert_float",
            FP::convert_float, vecmathlib::convert_float, n, accuracy);
      check("convert_int", FP::convert_int, vecmathlib::convert_int, x);
      check("floor", floor, vecmathlib::floor, x, accuracy);
      check("round", round, vecmathlib::round, x, accuracy);
    }
  }
  
  
  
  static void test_asin()
  {
    cout << "   testing asin acos atan...\n";
    for (int i=0; i<imax; ++i) {
      realvec_t const x = random(real_t(-1.0), real_t(+1.0));
      check("asin", asin, vecmathlib::asin, x, accuracy);
      check("acos", acos, vecmathlib::acos, x, accuracy);
    }
    for (int i=0; i<imax; ++i) {
      realvec_t const x = random(real_t(-100.0), real_t(+100.0));
      check("atan", atan, vecmathlib::atan, x, accuracy);
    }
  }
  
  static void test_asinh()
  {
    cout << "   testing asinh acosh atanh...\n";
    for (int i=0; i<imax; ++i) {
      realvec_t const x = random(real_t(-1000.0), real_t(+1000.0));
      // asinh loses accuracy (by definition, not by implementation?)
      check("asinh", asinh, vecmathlib::asinh, x, real_t(1.0e+3)*accuracy);
    }
    for (int i=0; i<imax; ++i) {
      realvec_t const x = random(real_t(1.0), real_t(1000.0));
      check("acosh", acosh, vecmathlib::acosh, x, accuracy);
    }
    for (int i=0; i<imax; ++i) {
      realvec_t const x = random(real_t(-1.0), real_t(+1.0));
      check("atanh", atanh, vecmathlib::atanh, x, accuracy);
    }
  }
  
  static real_t exp10(real_t x) { return pow(real_t(10.0), x); }
  static void test_exp()
  {
    cout << "   testing exp exp10 exp2 expm1...\n";
    for (int i=0; i<imax; ++i) {
      realvec_t const x = random(real_t(-100.0), real_t(+100.0));
      check("exp", exp, vecmathlib::exp, x, accuracy);
      check("exp10", exp10, vecmathlib::exp10, x, accuracy);
      check("exp2", exp2, vecmathlib::exp2, x, accuracy);
      check("expm1", expm1, vecmathlib::expm1, x, accuracy);
    }
  }
  
  static void test_log()
  {
    cout << "   testing log log10 log1p log2...\n";
    for (int i=0; i<imax; ++i) {
      realvec_t const x = random(real_t(1.0e-10), real_t(1.0e+10));
      check("log", log, vecmathlib::log, x, accuracy);
      check("log10", log10, vecmathlib::log10, x, accuracy);
      check("log1p", log1p, vecmathlib::log1p, x, accuracy);
      check("log2", log2, vecmathlib::log2, x, accuracy);
    }
  }
  
  static void test_pow()
  {
    cout << "   testing pow...\n";
    for (int i=0; i<imax; ++i) {
      realvec_t const x = random(real_t(0.001), real_t(1000.0));
      realvec_t const y = random(real_t(-10.0), real_t(+10.0));
      check("pow", pow, vecmathlib::pow, x, y, accuracy);
    }
  }
  
  static real_t rcp(real_t x) { return real_t(1.0)/x; }
  static void test_rcp()
  {
    cout << "   testing fmod rcp remainder...\n";
    for (int i=0; i<imax; ++i) {
      realvec_t const x = random(real_t(-10.0), real_t(+10.0));
      realvec_t const y = random(real_t(-10.0), real_t(+10.0));
      check("fmod", fmod, vecmathlib::fmod, x, y, accuracy);
      check("rcp", rcp, vecmathlib::rcp, x, accuracy);
      check("remainder", remainder, vecmathlib::remainder, x, y, accuracy);
    }
  }
  
  // void test_sin()
  // {
  //   for (int i=0; i<imax; ++i) {
  //     real_t const x = random(-10.0, 10.0);
  //     check("sin", sin, vecmathlib::sin, x, accuracy);
  //     check("cos", cos, vecmathlib::cos, x, accuracy);
  //     check("tan", tan, vecmathlib::tan, x, accuracy);
  //   }
  // }
  
  static void test_sinh()
  {
    cout << "   testing cosh sinh tanh...\n";
    for (int i=0; i<imax; ++i) {
      realvec_t const x = random(real_t(-10.0), real_t(+10.0));
      check("sinh", sinh, vecmathlib::sinh, x, accuracy);
      check("cosh", cosh, vecmathlib::cosh, x, accuracy);
      check("tanh", tanh, vecmathlib::tanh, x, accuracy);
    }
  }
  
  static real_t rsqrt(real_t x) { return real_t(1.0)/sqrt(x); }
  static void test_sqrt()
  {
    cout << "   testing rsqrt sqrt...\n";
    for (int i=0; i<imax; ++i) {
      realvec_t const x = random(real_t(0.0), real_t(1.0e+3));
      check("rsqrt", rsqrt, vecmathlib::rsqrt, x, accuracy);
      check("sqrt", sqrt, vecmathlib::sqrt, x, accuracy);
    }
  }
  
  
  
  static void test()
  {
    cout << "\n"
         << "Testing math functions for type "
      // << typeid(realvec_t).name()
         << realvec_t::name
         << ":\n";
    
    test_fabs();
    test_convert();
    
    test_asin();
    test_asinh();
    test_exp();
    test_log();
    test_pow();
    test_rcp();
    // test_sin();
    test_sinh();
    test_sqrt();
  }
};



int main(int argc, char** argv)
{
  using namespace vecmathlib;

  cout << "Testing math functions:\n";
  
  vecmathlib_test<realvec<float,1>>::test();
  vecmathlib_test<realvec<double,4>>::test();
  
  cout << "\n";
  if (num_errors == 0) {
    cout << "SUCCESS";
  } else {
    cout << "FAILURE";
  }
  cout << ": " << num_errors << " errors found\n";
  return num_errors == 0 ? 0 : 1;
}
