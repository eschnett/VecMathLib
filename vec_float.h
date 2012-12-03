// -*-C++-*-

#ifndef VEC_FLOAT_H
#define VEC_FLOAT_H

#include "floatprops.h"
#include "mathfuncs.h"
#include "vec_base.h"

#include <cmath>



namespace vecmathlib {
  
  template<> struct boolvec<float,1>;
  template<> struct intvec<float,1>;
  template<> struct realvec<float,1>;
  
  
  
  template<>
  struct boolvec<float,1>: floatprops<float>
  {
    static int const size = 1;
    typedef bool scalar_t;
    typedef bool bvector_t;
    
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
    //boolvec(vector_t x): v(x) {}
    boolvec(bool a): v(a) {}
    boolvec(bool const* as): v(as[0]) {}
    
    operator bvector_t() const { return v; }
    bool operator[](int n) const { return v; }
    boolvec& set_elt(int n, bool a) { return v=a, *this; }
    
    
    
    auto as_int() const -> intvec_t;      // defined after intvec
    auto convert_int() const -> intvec_t; // defined after intvec
    
    
    
    boolvec operator!() const { return !v; }
    
    boolvec operator&&(boolvec x) const { return v&&x.v; }
    boolvec operator||(boolvec x) const { return v||x.v; }
    boolvec operator==(boolvec x) const { return v==x.v; }
    boolvec operator!=(boolvec x) const { return v!=x.v; }
    
    bool all() const { return v; }
    bool any() const { return v; }
    
    
    
    // ifthen(condition, then-value, else-value)
    auto ifthen(intvec_t x,
                intvec_t y) const -> intvec_t; // defined after intvec
    auto ifthen(realvec_t x,
                realvec_t y) const -> realvec_t; // defined after realvec
    
  };
  
  
  
  template<>
  struct intvec<float,1>: floatprops<float>
  {
    static int const size = 1;
    typedef int_t scalar_t;
    typedef int_t ivector_t;
    
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
    //intvec(vector_t x): v(x) {}
    intvec(int_t a): v(a) {}
    intvec(int_t const* as): v(as[0]) {}
    
    operator ivector_t() const { return v; }
    int_t operator[](int n) const { return v; }
    intvec& set_elt(int n, int_t a) { return v=a, *this; }
    
    
    
    auto as_bool() const -> boolvec_t { return v; }
    auto convert_bool() const -> boolvec_t { return v; }
    auto as_float() const -> realvec_t; // defined after realvec
    auto convert_float() const -> realvec_t; // defined after realvec
    
    
    
    // Note: not all arithmetic operations are supported!
    
    intvec operator+() const { return +v; }
    intvec operator-() const { return -v; }
    
    intvec operator+(intvec x) const { return v+x.v; }
    intvec operator-(intvec x) const { return v-x.v; }
    
    intvec& operator+=(intvec const& x) { return *this=*this+x; }
    intvec& operator-=(intvec const& x) { return *this=*this-x; }


    
    intvec operator~() const { return ~v; }
    
    intvec operator&(intvec x) const { return v&x.v; }
    intvec operator|(intvec x) const { return v|x.v; }
    intvec operator^(intvec x) const { return v^x.v; }
    
    intvec& operator&=(intvec const& x) { return *this=*this&x; }
    intvec& operator|=(intvec const& x) { return *this=*this|x; }
    intvec& operator^=(intvec const& x) { return *this=*this^x; }
    
    
    
    intvec lsr(int_t n) const { return U(v)>>n; }
    intvec operator>>(int_t n) const { return v>>n; }
    intvec operator<<(int_t n) const { return v<<n; }
    intvec& operator>>=(int_t n) { return *this=*this>>n; }
    intvec& operator<<=(int_t n) { return *this=*this<<n; }
    
    intvec lsr(intvec n) const { return U(v)>>n; }
    intvec operator>>(intvec n) const { return v>>n; }
    intvec operator<<(intvec n) const { return v<<n; }
    intvec& operator>>=(intvec n) { return *this=*this>>n; }
    intvec& operator<<=(intvec n) { return *this=*this<<n; }
  };
  
  
  
  template<>
  struct realvec<float,1>: floatprops<float>
  {
    static int const size = 1;
    typedef real_t scalar_t;
    typedef real_t vector_t;
    
    static constexpr char const* const name = "<1*float>";
    inline void barrier() { asm("": "+x" (v)); }
    
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
    //realvec(vector_t x): v(x) {}
    realvec(real_t a): v(a) {}
    realvec(real_t const* as): v(as[0]) {}
    
    operator vector_t() const { return v; }
    real_t operator[](int n) const { return v; }
    realvec& set_elt(int n, real_t a) { return v=a, *this; }
    
    
    
    intvec_t as_int() const { return FP::as_int(v); }
    intvec_t convert_int() const { return MF::vml_convert_int(v); }
    
    
    
    realvec operator+() const { return +v; }
    realvec operator-() const { return -v; }
    
    realvec operator+(realvec x) const { return v+x.v; }
    realvec operator-(realvec x) const { return v-x.v; }
    realvec operator*(realvec x) const { return v*x.v; }
    realvec operator/(realvec x) const { return v/x.v; }
    
    realvec& operator+=(realvec const& x) { return *this=*this+x; }
    realvec& operator-=(realvec const& x) { return *this=*this-x; }
    realvec& operator*=(realvec const& x) { return *this=*this*x; }
    realvec& operator/=(realvec const& x) { return *this=*this/x; }
    
    real_t prod() const { return v; }
    real_t sum() const { return v; }
    
    
    
    boolvec_t operator==(realvec const& x) const { return v==x.v; }
    boolvec_t operator!=(realvec const& x) const { return v!=x.v; }
    boolvec_t operator<(realvec const& x) const { return v<x.v; }
    boolvec_t operator<=(realvec const& x) const { return v<=x.v; }
    boolvec_t operator>(realvec const& x) const { return v>x.v; }
    boolvec_t operator>=(realvec const& x) const { return v>=x.v; }
    
    
    
    realvec acos() const { return MF::vml_acos(*this); }
    realvec acosh() const { return MF::vml_acosh(*this); }
    realvec asin() const { return MF::vml_asin(*this); }
    realvec asinh() const { return MF::vml_asinh(*this); }
    realvec atan() const { return MF::vml_atan(*this); }
    realvec atanh() const { return MF::vml_atanh(*this); }
    realvec ceil() const { return MF::vml_ceil(*this); }
    realvec copysign(realvec y) const { return MF::vml_copysign(*this, y); }
    realvec cos() const { return MF::vml_cos(*this); }
    realvec cosh() const { return MF::vml_cosh(*this); }
    realvec exp() const { return MF::vml_exp(*this); }
    realvec exp10() const { return MF::vml_exp10(*this); }
    realvec exp2() const { return MF::vml_exp2(*this); }
    realvec expm1() const { return MF::vml_expm1(*this); }
    realvec fabs() const { return MF::vml_fabs(*this); }
    realvec floor() const { return MF::vml_floor(*this); }
    realvec fmod(realvec y) const { return MF::vml_fmod(*this, y); }
    intvec_t ilogb() const { return MF::vml_ilogb(*this); }
    realvec log() const { return MF::vml_log(*this); }
    realvec log10() const { return MF::vml_log10(*this); }
    realvec log1p() const { return MF::vml_log1p(*this); }
    realvec log2() const { return MF::vml_log2(*this); }
    realvec pow(realvec y) const { return MF::vml_pow(*this, y); }
    realvec rcp() const { return MF::vml_rcp(*this); }
    realvec remainder(realvec y) const { return MF::vml_remainder(*this, y); }
    realvec round() const { return MF::vml_round(*this); }
    realvec rsqrt() const { return MF::vml_rsqrt(*this); }
    realvec scalbn(intvec_t n) const { return MF::vml_scalbn(*this, n); }
    boolvec_t signbit() const { return MF::vml_signbit(*this); }
    realvec sin() const { return MF::vml_sin(*this); }
    realvec sinh() const { return MF::vml_sinh(*this); }
    realvec sqrt() const { return MF::vml_sqrt(*this); }
    realvec tan() const { return MF::vml_tan(*this); }
    realvec tanh() const { return MF::vml_tanh(*this); }
  };
  
  
  
  // boolvec definitions
  
  inline
  auto boolvec<float,1>::as_int() const -> intvec_t
  {
    return v;
  }
  
  inline
  auto boolvec<float,1>::convert_int() const -> intvec_t
  {
    return v;
  }
  
  inline
  auto boolvec<float,1>::ifthen(intvec_t x, intvec_t y) const -> intvec_t
  {
    return v ? x.v : y.v;
  }
  
  inline
  auto boolvec<float,1>::ifthen(realvec_t x, realvec_t y) const -> realvec_t
  {
    return v ? x.v : y.v;
  }

  
  
  // intvec definitions
  
  inline auto intvec<float,1>::as_float() const -> realvec_t
  {
    return FP::as_float(v);
  }
  
  inline auto intvec<float,1>::convert_float() const -> realvec_t
  {
    return MF::vml_convert_float(v);
  }
  
} // namespace vecmathlib

#endif  // #ifndef VEC_FLOAT_H
