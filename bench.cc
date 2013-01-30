// -*-C++-*-

#define NDEBUG

#include "vecmathlib.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <typeinfo>

#include <sys/time.h>

using namespace std;



typedef unsigned long long ticks;
inline ticks getticks()
{
  ticks a, d;
  asm volatile("rdtsc" : "=a" (a), "=d" (d)); 
  return a | (d << 32); 
}
inline double elapsed(ticks t1, ticks t0)
{
  return t1-t0;
}

double get_sys_time()
{
  timeval tp;
  gettimeofday(&tp, NULL);
  return tp.tv_sec + 1.0e-6 * tp.tv_usec;
}

double measure_tick()
{
  ticks const rstart = getticks();
  double const wstart = get_sys_time();
  while (get_sys_time() - wstart < 0.1) {
    // do nothing, just wait
  }
  ticks const rend = getticks();
  double const wend = get_sys_time();
  assert(wend-wstart >= 0.09);
  return (wend - wstart) / elapsed(rend, rstart);
}



template<typename T, int N>
struct pseudovec {
  T v[N];
  static int const size = N;
  static string name()
  {
    string base;
    if (typeid(T) == typeid(float)) base = "float";
    else if (typeid(T) == typeid(double)) base = "double";
    else base = typeid(T).name();
    return string("<")+to_string(N)+"*"+base+">";
  }
  pseudovec() {}
  pseudovec(T const& w) { for (int i=0; i<N; ++i) v[i]=w; }
  pseudovec& set_elt(int i, T const& w) { v[i]=w; return *this; }
  T operator[] (int i) const { return v[i]; }
  pseudovec& operator+=(pseudovec const& x)
  {
    for (int i=0; i<N; ++i) v[i]+=x.v[i]; return *this;
  }
};
template<typename T, int N>
pseudovec<T,N> map(T f(T), pseudovec<T,N> const& x)
{
  pseudovec<T,N> r;
  for (int i=0; i<N; ++i) r.set_elt(i, f(x[i]));
  return r;
}
template<typename T, int N>
pseudovec<T,N> atan(pseudovec<T,N> const& x) { return map(std::atan, x); }
template<typename T, int N>
pseudovec<T,N> cos(pseudovec<T,N> const& x) { return map(std::cos, x); }
template<typename T, int N>
pseudovec<T,N> exp(pseudovec<T,N> const& x) { return map(std::exp, x); }
template<typename T, int N>
pseudovec<T,N> log(pseudovec<T,N> const& x) { return map(std::log, x); }
template<typename T, int N>
pseudovec<T,N> sin(pseudovec<T,N> const& x) { return map(std::sin, x); }
template<typename T, int N>
pseudovec<T,N> sqrt(pseudovec<T,N> const& x) { return map(std::sqrt, x); }



double global_result = 0.0;
template<typename realvec_t>
void save_result(realvec_t result)
{
  for (int i=0; i<realvec_t::size; ++i) {
    global_result += result[i];
  }
}



template<typename T> inline T identity(T x) { return x; }

#define DECLARE_FUNCTOR(func)                   \
  template<typename T>                          \
  struct functor_##func {                       \
    T operator()(T x) { return func(x); }       \
  }

DECLARE_FUNCTOR(identity);
DECLARE_FUNCTOR(sqrt);
DECLARE_FUNCTOR(exp);
DECLARE_FUNCTOR(log);
DECLARE_FUNCTOR(sin);
DECLARE_FUNCTOR(cos);
DECLARE_FUNCTOR(atan);



template<typename realvec_t, typename func_t>
double bench_func()
{
  realvec_t x0, dx;
  for (int i=0; i<realvec_t::size; ++i) {
    x0.set_elt(i, 1.0f + float(i));
    dx.set_elt(i, 1.0e-6f);
  }
  realvec_t x, y;
  ticks t0, t1;
  double const cycles_per_tick = 1.0; // measure_tick();
  int const numiters = 10000000;
  
  func_t func;
  t0 = getticks();
  x = y = x0;
  for (int n=0; n<numiters; ++n) {
    y += func(x);
    x += dx;
  }
  t1 = getticks();
  save_result(y);
  
  return cycles_per_tick * elapsed(t1,t0) * realvec_t::size / numiters;
}

template<typename realvec_t>
void bench()
{
  cout << "identity(" << realvec_t::name() << "):" << flush;
  double const cycles_identity =
    bench_func<realvec_t, functor_identity<realvec_t>>();
  cout << "   " << cycles_identity << " cycles\n";
  
  cout << "sqrt(" << realvec_t::name() << "):" << flush;
  double const cycles_sqrt =
    bench_func<realvec_t, functor_sqrt<realvec_t>>();
  cout << "   " << cycles_sqrt - cycles_identity << " cycles\n";
  
  cout << "exp(" << realvec_t::name() << "):" << flush;
  double const cycles_exp =
    bench_func<realvec_t, functor_exp<realvec_t>>();
  cout << "   " << cycles_exp - cycles_identity << " cycles\n";
  
  cout << "log(" << realvec_t::name() << "):" << flush;
  double const cycles_log =
    bench_func<realvec_t, functor_log<realvec_t>>();
  cout << "   " << cycles_log - cycles_identity << " cycles\n";
  
  cout << "sin(" << realvec_t::name() << "):" << flush;
  double const cycles_sin =
    bench_func<realvec_t, functor_sin<realvec_t>>();
  cout << "   " << cycles_sin - cycles_identity << " cycles\n";
  
  cout << "cos(" << realvec_t::name() << "):" << flush;
  double const cycles_cos =
    bench_func<realvec_t, functor_cos<realvec_t>>();
  cout << "   " << cycles_cos - cycles_identity << " cycles\n";
  
  cout << "atan(" << realvec_t::name() << "):" << flush;
  double const cycles_atan =
    bench_func<realvec_t, functor_atan<realvec_t>>();
  cout << "   " << cycles_atan - cycles_identity << " cycles\n";
}



int main(int argc, char** argv)
{
  using namespace vecmathlib;

  cout << "Benchmarking math functions:\n"
       << "\n";
  
  bench<pseudovec<float,1>>();
#ifdef VECMATHLIB_HAVE_VEC_FLOAT_1
  bench<realvec<float,1>>();
#endif
#ifdef VECMATHLIB_HAVE_VEC_FLOAT_4
  bench<pseudovec<float,4>>();
  bench<realvec<float,4>>();
#endif
#ifdef VECMATHLIB_HAVE_VEC_FLOAT_8
  bench<pseudovec<float,8>>();
  bench<realvec<float,8>>();
#endif
  
  cout << "\n";
  
  bench<pseudovec<double,1>>();
#ifdef VECMATHLIB_HAVE_VEC_DOUBLE_1
  bench<realvec<double,1>>();
#endif
#ifdef VECMATHLIB_HAVE_VEC_DOUBLE_2
  bench<pseudovec<double,2>>();
  bench<realvec<double,2>>();
#endif
#ifdef VECMATHLIB_HAVE_VEC_DOUBLE_4
  bench<pseudovec<double,4>>();
  bench<realvec<double,4>>();
#endif
  
  cout << "\n"
       << "Outputting global result to prevent optimisation: "
       << global_result << "\n";
  
  return 0;
}
