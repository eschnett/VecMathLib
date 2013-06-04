// -*-C++-*-

#ifndef VEC_DOUBLE_QPX_H
#define VEC_DOUBLE_QPX_H

#include "floatprops.h"
#include "mathfuncs.h"
#include "vec_base.h"

#include <cmath>

// QPX intrinsics
#include <builtins.h>
#include <mass_simd.h>
// #include <qpxintin.h>



namespace vecmathlib {
  
#define VECMATHLIB_HAVE_VEC_DOUBLE_4
  template<> struct boolvec<double,4>;
  template<> struct intvec<double,4>;
  template<> struct realvec<double,4>;
  
  
  
  template<>
  struct boolvec<double,4>: floatprops<double>
  {
    static int const size = 4;
    typedef bool scalar_t;
    typedef vector4double bvector_t;
    static int const alignment = sizeof(bvector_t);
    
    static_assert(size * sizeof(real_t) == sizeof(bvector_t),
                  "vector size is wrong");
    
  private:
    // canonical true is +1.0, canonical false is -1.0
    // >=0 is true, -0 is true, nan is false
    static real_t from_bool(bool a) { return a ? +1.0 : -1.0; }
    static bool to_bool(real_t a) { return a>=0.0; }
  public:
    
    typedef boolvec boolvec_t;
    typedef intvec<real_t, size> intvec_t;
    typedef realvec<real_t, size> realvec_t;
    
    // Short names for type casts
    typedef real_t R;
    typedef int_t I;
    typedef uint_t U;
    typedef realvec_t RV;
    typedef intvec_t IV;
    typedef boolvec_t BV;
    typedef floatprops<real_t> FP;
    typedef mathfuncs<realvec_t> MF;
    
    
    
    bvector_t v;
    
    boolvec() {}
    // Can't have a non-trivial copy constructor; if so, objects won't
    // be passed in registers
    // boolvec(boolvec const& x): v(x.v) {}
    // boolvec& operator=(boolvec const& x) { return v=x.v, *this; }
    boolvec(bvector_t x): v(x) {}
    boolvec(bool a): v(vec_splats(from_bool(a))) {}
    boolvec(bool const* as)
    {
      for (int d=0; d<size; ++d) set_elt(d, from_bool(as[d]));
    }
    
    operator bvector_t() const { return v; }
    bool operator[](int n) const { return to_bool(((real_t const*)&v)[n]); }
    boolvec& set_elt(int n, bool a)
    {
      return ((real_t*)&v)[n]=from_bool(a), *this;
    }
    
    
    
    intvec_t as_int() const;      // defined after intvec
    intvec_t convert_int() const; // defined after intvec
    
    
    
    boolvec operator!() const { return vec_not(v); }
    
    boolvec operator&&(boolvec x) const { return vec_and(v, x.v); }
    boolvec operator||(boolvec x) const { return vec_or(v, x.v); }
    boolvec operator==(boolvec x) const { return vec_logical(v, x.v, 0b1001); }
    boolvec operator!=(boolvec x) const { return vec_xor(v, x.v); }
    
    bool all() const
    {
      return (*this)[0] && (*this)[1] && (*this)[2] && (*this)[3];
    }
    bool any() const
    {
      return (*this)[0] || (*this)[1] || (*this)[2] || (*this)[3];
    }
    
    
    
    // ifthen(condition, then-value, else-value)
    intvec_t ifthen(intvec_t x, intvec_t y) const; // defined after intvec
    realvec_t ifthen(realvec_t x, realvec_t y) const; // defined after realvec
  };
  
  
  
  template<>
  struct intvec<double,4>: floatprops<double>
  {
    static int const size = 4;
    typedef int_t scalar_t;
    typedef vector4double ivector_t;
    static int const alignment = sizeof(ivector_t);
    
    static_assert(size * sizeof(real_t) == sizeof(ivector_t),
                  "vector size is wrong");
    
    typedef boolvec<real_t, size> boolvec_t;
    typedef intvec intvec_t;
    typedef realvec<real_t, size> realvec_t;
    
    // Short names for type casts
    typedef real_t R;
    typedef int_t I;
    typedef uint_t U;
    typedef realvec_t RV;
    typedef intvec_t IV;
    typedef boolvec_t BV;
    typedef floatprops<real_t> FP;
    typedef mathfuncs<realvec_t> MF;
    
    
    
    ivector_t v;
    
    intvec() {}
    // Can't have a non-trivial copy constructor; if so, objects won't
    // be passed in registers
    // intvec(intvec const& x): v(x.v) {}
    // intvec& operator=(intvec const& x) { return v=x.v, *this; }
    intvec(ivector_t x): v(x) {}
    intvec(int_t a): v(vec_splats(as_float(a))) {}
    intvec(int_t const* as)
    {
      for (int d=0; d<size; ++d) set_elt(d, as[d]);
    }
    static intvec iota() { return {0, 1, 2, 3}; }
    
    operator ivector_t() const { return v; }
    int_t operator[](int n) const { return ((int_t const*)&v)[n]; }
    intvec& set_elt(int n, int_t a) { return ((int_t*)&v)[n]=a, *this; }
    
    
    
    // Vector casts do not change the bit battern
    boolvec_t as_bool() const { return v; }
    realvec_t as_float() const;      // defined after realvec
    realvec_t convert_float() const; // defined after realvec
  };
  
  
  
  template<>
  struct realvec<double,4>: floatprops<double>
  {
    static int const size = 4;
    typedef real_t scalar_t;
    typedef __vector double vector_t;
    static int const alignment = sizeof(vector_t);
    
    static char const* name() { return "<QPX:4*double>"; }
    void barrier() { __asm__("": "+v" (v)); }
    
    static_assert(size * sizeof(real_t) == sizeof(vector_t),
                  "vector size is wrong");
    
    typedef boolvec<real_t, size> boolvec_t;
    typedef intvec<real_t, size> intvec_t;
    typedef realvec realvec_t;
    
    // Short names for type casts
    typedef real_t R;
    typedef int_t I;
    typedef uint_t U;
    typedef realvec_t RV;
    typedef intvec_t IV;
    typedef boolvec_t BV;
    typedef floatprops<real_t> FP;
    typedef mathfuncs<realvec_t> MF;
    
    
    
    vector_t v;
    
    realvec() {}
    // Can't have a non-trivial copy constructor; if so, objects won't
    // be passed in registers
    // realvec(realvec const& x): v(x.v) {}
    // realvec& operator=(realvec const& x) { return v=x.v, *this; }
    realvec(vector_t x): v(x) {}
    realvec(real_t a): v(vec_splats(a)) {}
    realvec(real_t const* as)
    {
      for (int d=0; d<size; ++d) set_elt(d, as[d]);
    }
    
    operator vector_t() const { return v; }
    real_t operator[](int n) const { return ((real_t const*)&v)[n]; }
    realvec& set_elt(int n, real_t a) { return ((real_t*)&v)[n]=a, *this; }
    
    
    
    typedef vecmathlib::mask_t<realvec_t> mask_t;
    
    static realvec_t loada(real_t const* p)
    {
      VML_ASSERT(intptr_t(p) % alignment == 0);
      return vec_ld(0, p);
    }
    static realvec_t loadu(real_t const* p)
    {
      realvec_t v0 = vec_ld(0, p);
      realvec_t v1 = vec_ld(15, p);
      return vec_perm(v0.v, v1.v, vec_lvsl(0, p));
    }
    static realvec_t loadu(real_t const* p, std::ptrdiff_t ioff)
    {
      VML_ASSERT(intptr_t(p) % alignment == 0);
      if (ioff % realvec::size == 0) return loada(p+ioff);
      return loadu(p+ioff);
    }
    realvec_t loada(real_t const* p, mask_t const& m) const
    {
      VML_ASSERT(intptr_t(p) % alignment == 0);
      if (__builtin_expect(all(m.m), true)) {
        return loada(p);
      } else {
        return m.m.ifthen(loada(p), *this);
      }
    }
    realvec_t loadu(real_t const* p, mask_t const& m) const
    {
      if (__builtin_expect(m.all_m, true)) {
        return loadu(p);
      } else {
        return m.m.ifthen(loadu(p), *this);
      }
    }
    realvec_t loadu(real_t const* p, std::ptrdiff_t ioff, mask_t const& m) const
    {
      VML_ASSERT(intptr_t(p) % alignment == 0);
      if (ioff % realvec::size == 0) return loada(p+ioff, m);
      return loadu(p+ioff, m);
    }
    
    void storea(real_t* p) const
    {
      VML_ASSERT(intptr_t(p) % alignment == 0);
      vec_st(v, 0, p);
    }
    void storeu(real_t* p) const
    {
      // Vector stores would require vector loads, which would need to
      // be atomic
      // TODO: see <https://developer.apple.com/hardwaredrivers/ve/alignment.html> for good ideas
      p[0] = (*this)[0];
      p[1] = (*this)[1];
      p[2] = (*this)[2];
      p[3] = (*this)[3];
    }
    void storeu(real_t* p, std::ptrdiff_t ioff) const
    {
      VML_ASSERT(intptr_t(p) % alignment == 0);
      if (ioff % realvec::size == 0) return storea(p+ioff);
      storeu(p+ioff);
    }
    void storea(real_t* p, mask_t const& m) const
    {
      VML_ASSERT(intptr_t(p) % alignment == 0);
      if (__builtin_expect(m.all_m, true)) {
        storea(p);
      } else {
        if (m.m[0]) p[0] = (*this)[0];
        if (m.m[1]) p[1] = (*this)[1];
        if (m.m[0]) p[2] = (*this)[2];
        if (m.m[1]) p[3] = (*this)[3];
      }
    }
    void storeu(real_t* p, mask_t const& m) const
    {
      if (__builtin_expect(m.all_m, true)) {
        storeu(p);
      } else {
        if (m.m[0]) p[0] = (*this)[0];
        if (m.m[1]) p[1] = (*this)[1];
        if (m.m[0]) p[2] = (*this)[2];
        if (m.m[1]) p[3] = (*this)[3];
      }
    }
    void storeu(real_t* p, std::ptrdiff_t ioff, mask_t const& m) const
    {
      VML_ASSERT(intptr_t(p) % alignment == 0);
      if (ioff % realvec::size == 0) return storea(p+ioff, m);
      storeu(p+ioff, m);
    }
    
    
    
    intvec_t as_int() const { return v; }
    intvec_t convert_int() const { return vec_ctid(v); }
    
    
    
    realvec operator+() const { return *this; }
    realvec operator-() const { return vec_neg(v); }
    
    realvec operator+(realvec x) const { return vec_add(v, x.v); }
    realvec operator-(realvec x) const { return vec_sub(v, x.v); }
    realvec operator*(realvec x) const { return vec_mul(v, x.v); }
    realvec operator/(realvec x) const { return vec_swdiv_nochk(v, x.v); }
    
    realvec& operator+=(realvec const& x) { return *this=*this+x; }
    realvec& operator-=(realvec const& x) { return *this=*this-x; }
    realvec& operator*=(realvec const& x) { return *this=*this*x; }
    realvec& operator/=(realvec const& x) { return *this=*this/x; }
    
    real_t prod() const
    {
      return (*this)[0] * (*this)[1] * (*this)[2] * (*this)[3];
    }
    real_t sum() const
    {
      return (*this)[0] + (*this)[1] + (*this)[2] + (*this)[3];
    }
    
    
    
    boolvec_t operator==(realvec const& x) const { return vec_cmpeq(v, x.v); }
    boolvec_t operator!=(realvec const& x) const { return ! (*this == x); }
    boolvec_t operator<(realvec const& x) const { return vec_cmplt(v, x.v); }
    boolvec_t operator<=(realvec const& x) const { return ! (*this > x); }
    boolvec_t operator>(realvec const& x) const { return vec_cmpgt(v, x.v); }
    boolvec_t operator>=(realvec const& x) const { return ! (*this < x); }
    
    
    
    realvec acos() const { return MF::vml_acos(*this); }
    realvec acosh() const { return MF::vml_acosh(*this); }
    realvec asin() const { return MF::vml_asin(*this); }
    realvec asinh() const { return MF::vml_asinh(*this); }
    realvec atan() const { return MF::vml_atan(*this); }
    realvec atanh() const { return MF::vml_atanh(*this); }
    realvec cbrt() const { return MF::vml_cbrt(*this); }
    realvec ceil() const { return vec_ceil(v); }
    realvec copysign(realvec y) const { return vec_cpsgn(v, y.v); }
    realvec cos() const { return MF::vml_cos(*this); }
    realvec cosh() const { return MF::vml_cosh(*this); }
    realvec exp() const { return MF::vml_exp(*this); }
    realvec exp10() const { return MF::vml_exp10(*this); }
    realvec exp2() const { return MF::vml_exp2(*this); }
    realvec expm1() const { return MF::vml_expm1(*this); }
    realvec fabs() const { return vec_abs(v); }
    realvec fdim(realvec y) const { return MF::vml_fdim(*this, y); }
    realvec floor() const { return vec_floor(v); }
    realvec fma(realvec y, realvec z) const { return vec_madd(v, y.v, z.v); }
    realvec fmax(realvec y) const { return MF::vml_fmax(v, y.v); }
    realvec fmin(realvec y) const { return MV::vml_fmin(v, y.v); }
    realvec fmod(realvec y) const { return MF::vml_fmod(*this, y); }
    realvec hypot(realvec y) const { return MF::vml_hypot(*this, y); }
    intvec_t ilogb() const { return MF::vml_ilogb(*this); }
    boolvec_t isfinite() const { return MF::vml_isfinite(*this); }
    boolvec_t isinf() const { return MF::vml_isinf(*this); }
    boolvec_t isnan() const { return vec_testnan(v); }
    boolvec_t isnormal() const { return MF::vml_isnormal(*this); }
    realvec ldexp(int_t n) const { return MF::vml_ldexp(*this, n); }
    realvec ldexp(intvec_t n) const { return MF::vml_ldexp(*this, n); }
    realvec log() const { return MF::vml_log(*this); }
    realvec log10() const { return MF::vml_log10(*this); }
    realvec log1p() const { return MF::vml_log1p(*this); }
    realvec log2() const { return MF::vml_log2(*this); }
    realvec pow(realvec y) const { return MF::vml_pow(*this, y); }
    realvec rcp() const { return RV(1.0)/(*this); }
    realvec remainder(realvec y) const { return MF::vml_remainder(*this, y); }
    realvec rint() const { return MF::vml_rint(*this); }
    realvec round() const { return vec_round(v); }
    realvec rsqrt() const
    {
      realvec x = *this;
      realvec r = vec_rsqrte(x.v); // this is only an approximation
      // TODO: use fma
      // one Newton iteration (see vml_rsqrt)
      r += RV(0.5)*r * (RV(1.0) - x * r*r);
      return r;
    }
    boolvec_t signbit() const { return MF::vml_signbit(*this); }
    realvec sin() const { return MF::vml_sin(*this); }
    realvec sinh() const { return MF::vml_sinh(*this); }
    realvec sqrt() const { return vec_sqrtsw_nochk(v); }
    realvec tan() const { return MF::vml_tan(*this); }
    realvec tanh() const { return MF::vml_tanh(*this); }
    realvec trunc() const { return vec_trunc(v); }
  };
  
  
  
  // boolvec definitions
  
  inline
  auto boolvec<double,4>::as_int() const -> intvec_t
  {
    return v;
  }
  
  inline
  auto boolvec<double,4>::convert_int() const -> intvec_t
  {
    return ifthen(IV(I(1)), IV(I(0)));
  }
  
  inline
  auto boolvec<double,4>::ifthen(intvec_t x, intvec_t y) const -> intvec_t
  {
    return ifthen(x.as_float(), y.as_float());
  }
  
  inline
  auto boolvec<double,4>::ifthen(realvec_t x, realvec_t y) const -> realvec_t
  {
    return vec_sel(y.v, x.v, v);
  }
  
  
  
  // intvec definitions
  
  inline auto intvec<double,4>::as_float() const -> realvec_t
  {
    return v;
  }
  
  inline auto intvec<double,4>::convert_float() const -> realvec_t
  {
    return vec_cfid(v);
  }
  
} // namespace vecmathlib

#endif  // #ifndef VEC_DOUBLE_QPX_H