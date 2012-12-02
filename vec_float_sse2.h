// -*-C++-*-

#ifndef VEC_FLOAT_SSE2_H
#define VEC_FLOAT_SSE2_H

#include "floatprops.h"
#include "mathfuncs.h"
#include "vec_base.h"

#include <cmath>

// SSE2 intrinsics
#include <xmmintrin.h>
#if defined __SSE4_1__          // Intel's SSE 4.1
#  include <smmintrin.h>
#endif
#if defined __SSE4A__           // AMD's SSE 4a
#  include <ammintrin.h>
#endif



namespace vecmathlib {
  
  template<> struct boolvec<float,4>;
  template<> struct intvec<float,4>;
  template<> struct realvec<float,4>;
  
  
  
  template<>
  struct boolvec<float,4>: floatprops<float>
  {
    static int const size = 4;
    typedef bool scalar_t;
    typedef __m128 bvector_t;
    
    static_assert(size * sizeof(real_t) == sizeof(bvector_t),
                  "vector size is wrong");
    
  private:
    // true values have the sign bit set, false values have it unset
    static uint_t from_bool(bool a) { return - uint_t(a); }
    static bool to_bool(uint_t a) { return int_t(a) < int_t(0); }
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
    boolvec(boolvec const& x): v(x.v) {}
    boolvec& operator=(boolvec const& x) { return v=x.v, *this; }
    boolvec(bvector_t x): v(x) {}
    boolvec(bool a):
    v(_mm_castsi128_ps(_mm_set1_epi32(from_bool(a)))) {}
    boolvec(bool const* as):
    v(_mm_castsi128_ps(_mm_set_epi32(from_bool(as[3]),
                                     from_bool(as[2]),
                                     from_bool(as[1]),
                                     from_bool(as[0])))) {}
    
    operator bvector_t() const { return v; }
    bool operator[](int n) const { return to_bool(((uint_t const*)&v)[n]); }
    boolvec& set_elt(int n, bool a)
    {
      return ((int_t*)&v)[n] = from_bool(a), *this;
    }
    
    
    
    auto as_int() const -> intvec_t;      // defined after intvec
    auto convert_int() const -> intvec_t; // defined after intvec
    
    
    
    boolvec operator!() const { return _mm_xor_ps(boolvec(true), v); }
    
    boolvec operator&&(boolvec x) const { return _mm_and_ps(v, x.v); }
    boolvec operator||(boolvec x) const { return _mm_or_ps(v, x.v); }
    boolvec operator==(boolvec x) const { return !(*this==x); }
    boolvec operator!=(boolvec x) const { return _mm_xor_ps(v, x.v); }
    
    bool all() const
    {
      return (*this)[0] && (*this)[1] && (*this)[2] && (*this)[3];
    }
    bool any() const
    {
      return (*this)[0] || (*this)[1] || (*this)[2] || (*this)[3];
    }
    
    
    
    // ifthen(condition, then-value, else-value)
    auto ifthen(intvec_t x,
                intvec_t y) const -> intvec_t; // defined after intvec
    auto ifthen(realvec_t x,
                realvec_t y) const -> realvec_t; // defined after realvec
    
  };
  
  
  
  template<>
  struct intvec<float,4>: floatprops<float>
  {
    static int const size = 4;
    typedef int_t scalar_t;
    typedef __m128i ivector_t;
    
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
    intvec(intvec const& x): v(x.v) {}
    intvec& operator=(intvec const& x) { return v=x.v, *this; }
    intvec(ivector_t x): v(x) {}
    intvec(int_t a): v(_mm_set1_epi32(a)) {}
    intvec(int_t const* as): v(_mm_set_epi32(as[3], as[2], as[1], as[0])) {}
    
    operator ivector_t() const { return v; }
    int_t operator[](int n) const { return ((int_t const*)&v)[n]; }
    intvec& set_elt(int n, int_t a) { return ((int_t*)&v)[n]=a, *this; }
    
    
    
    auto as_bool() const -> boolvec_t { return _mm_castsi128_ps(v); }
    auto convert_bool() const -> boolvec_t
    {
      // Result: convert_bool(0)=false, convert_bool(else)=true
      return ! IV(_mm_cmpeq_epi32(v, IV(0))).as_bool();
    }
    auto as_float() const -> realvec_t; // defined after realvec
    auto convert_float() const -> realvec_t; // defined after realvec
    
    
    
    // Note: not all arithmetic operations are supported!
    
    intvec operator+() const { return *this; }
    intvec operator-() const { return IV(0) - *this; }
    
    intvec operator+(intvec x) const { return _mm_add_epi32(v, x.v); }
    intvec operator-(intvec x) const { return _mm_sub_epi32(v, x.v); }
    
    intvec& operator+=(intvec const& x) { return *this=*this+x; }
    intvec& operator-=(intvec const& x) { return *this=*this-x; }
    
    
    
    intvec operator~() const { return IV(~U(0)) ^ *this; }
    
    intvec operator&(intvec x) const
    {
      return _mm_castps_si128(_mm_and_ps(_mm_castsi128_ps(v),
                                         _mm_castsi128_ps(x.v)));
    }
    intvec operator|(intvec x) const
    {
      return _mm_castps_si128(_mm_or_ps(_mm_castsi128_ps(v),
                                        _mm_castsi128_ps(x.v)));
    }
    intvec operator^(intvec x) const
    {
      return _mm_castps_si128(_mm_xor_ps(_mm_castsi128_ps(v),
                                         _mm_castsi128_ps(x.v)));
    }
    
    intvec& operator&=(intvec const& x) { return *this=*this&x; }
    intvec& operator|=(intvec const& x) { return *this=*this|x; }
    intvec& operator^=(intvec const& x) { return *this=*this^x; }
    
    
    
    intvec lsr(int_t n) const { return _mm_srli_epi32(v, n); }
    intvec operator>>(int_t n) const { return _mm_srai_epi32(v, n); }
    intvec operator<<(int_t n) const { return _mm_slli_epi32(v, n); }
    intvec& operator>>=(int_t n) { return *this=*this>>n; }
    intvec& operator<<=(int_t n) { return *this=*this<<n; }
    
    intvec lsr(intvec n) const { return _mm_srl_epi32(v, n.v); }
    intvec operator>>(intvec n) const { return _mm_sra_epi32(v, n.v); }
    intvec operator<<(intvec n) const { return _mm_sll_epi32(v, n.v); }
    intvec& operator>>=(intvec n) { return *this=*this>>n; }
    intvec& operator<<=(intvec n) { return *this=*this<<n; }
  };
  
  
  
  template<>
  struct realvec<float,4>: floatprops<float>
  {
    static int const size = 4;
    typedef real_t scalar_t;
    typedef __m128 vector_t;
    
    static constexpr char const* const name = "<SSE2:4*float>";
    inline void barrier() { asm("": "+x" (v)); }
    
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
    realvec(realvec const& x): v(x.v) {}
    realvec& operator=(realvec const& x) { return v=x.v, *this; }
    realvec(vector_t x): v(x) {}
    realvec(real_t a): v(_mm_set1_ps(a)) {}
    realvec(real_t const* as): v(_mm_set_ps(as[3], as[2], as[1], as[0])) {}
    
    operator vector_t() const { return v; }
    real_t operator[](int n) const { return ((real_t const*)&v)[n]; }
    realvec& set_elt(int n, real_t a) { return ((real_t*)&v)[n]=a, *this; }
    
    
    
    intvec_t as_int() const { return _mm_castps_si128(v); }
    intvec_t convert_int() const { return _mm_cvtps_epi32(v); }
    
    
    
    realvec operator+() const { return *this; }
    realvec operator-() const { return RV(0.0) - *this; }
    
    realvec operator+(realvec x) const { return _mm_add_ps(v, x.v); }
    realvec operator-(realvec x) const { return _mm_sub_ps(v, x.v); }
    realvec operator*(realvec x) const { return _mm_mul_ps(v, x.v); }
    realvec operator/(realvec x) const { return _mm_div_ps(v, x.v); }
    
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
    
    
    
    boolvec_t operator==(realvec const& x) const
    {
      return _mm_cmpeq_ps(v, x.v);
    }
    boolvec_t operator!=(realvec const& x) const
    {
      return _mm_cmpneq_ps(v, x.v);
    }
    boolvec_t operator<(realvec const& x) const
    {
      return _mm_cmplt_ps(v, x.v);
    }
    boolvec_t operator<=(realvec const& x) const
    {
      return _mm_cmple_ps(v, x.v);
    }
    boolvec_t operator>(realvec const& x) const
    {
      return _mm_cmpgt_ps(v, x.v);
    }
    boolvec_t operator>=(realvec const& x) const
    {
      return _mm_cmpge_ps(v, x.v);
    }
    
    
    
    realvec acos() const { return MF::vml_acos(*this); }
    realvec acosh() const { return MF::vml_acosh(*this); }
    realvec asin() const { return MF::vml_asin(*this); }
    realvec asinh() const { return MF::vml_asinh(*this); }
    realvec atan() const { return MF::vml_atan(*this); }
    realvec atanh() const { return MF::vml_atanh(*this); }
    realvec ceil() const { return _mm_ceil_ps(v); }
    realvec copysign(realvec y) const { return MF::vml_copysign(*this, y); }
    realvec cos() const { return MF::vml_cos(*this); }
    realvec cosh() const { return MF::vml_cosh(*this); }
    realvec exp() const { return MF::vml_exp(*this); }
    realvec exp10() const { return MF::vml_exp10(*this); }
    realvec exp2() const { return MF::vml_exp2(*this); }
    realvec expm1() const { return MF::vml_expm1(*this); }
    realvec fabs() const { return MF::vml_fabs(*this); }
    realvec floor() const { return _mm_floor_ps(v); }
    realvec fmod(realvec y) const { return MF::vml_fmod(*this, y); }
    intvec_t ilogb() const { return MF::vml_ilogb(*this); }
    realvec log() const { return MF::vml_log(*this); }
    realvec log10() const { return MF::vml_log10(*this); }
    realvec log1p() const { return MF::vml_log1p(*this); }
    realvec log2() const { return MF::vml_log2(*this); }
    realvec pow(realvec y) const { return MF::vml_pow(*this, y); }
    realvec rcp() const
    {
      realvec x = *this;
      realvec r = _mm_rcp_ps(x); // this is only an approximation
      r *= RV(2.0) - r*x;        // one Newton iteration (see vml_rcp)
      return r;
    }
    realvec remainder(realvec y) const { return MF::vml_remainder(*this, y); }
    realvec round() const { return _mm_round_ps(v, _MM_FROUND_NINT); }
    realvec rsqrt() const
    {
      realvec x = *this;
      realvec r = _mm_rsqrt_ps(x);    // this is only an approximation
      r *= RV(1.5) - RV(0.5)*x * r*r; // one Newton iteration (see vml_rsqrt)
      return r;
    }
    realvec scalbn(intvec_t n) const { return MF::vml_scalbn(*this, n); }
    boolvec_t signbit() const { return v; }
    realvec sin() const { return MF::vml_sin(*this); }
    realvec sinh() const { return MF::vml_sinh(*this); }
    realvec sqrt() const { return _mm_sqrt_ps(v); }
    realvec tan() const { return MF::vml_tan(*this); }
    realvec tanh() const { return MF::vml_tanh(*this); }
  };
  
  
  
  // boolvec definitions
  
  inline
  auto boolvec<float,4>::as_int() const -> intvec_t
  {
    return _mm_castps_si128(v);
  }
  
  inline
  auto boolvec<float,4>::convert_int() const -> intvec_t
  {
    return lsr(as_int(), bits-1);
  }
  
  inline
  auto boolvec<float,4>::ifthen(intvec_t x, intvec_t y) const -> intvec_t
  {
    return ifthen(x.as_float(), y.as_float()).as_int();
  }
  
  inline
  auto boolvec<float,4>::ifthen(realvec_t x, realvec_t y) const -> realvec_t
  {
    return _mm_blendv_ps(y.v, x.v, v);
  }

  
  
  // intvec definitions
  
  inline auto intvec<float,4>::as_float() const -> realvec_t
  {
    return _mm_castsi128_ps(v);
  }
  
  inline auto intvec<float,4>::convert_float() const -> realvec_t
  {
    return _mm_cvtepi32_ps(v);
  }
  
} // namespace vecmathlib

#endif  // #ifndef VEC_FLOAT_SSE2_H
