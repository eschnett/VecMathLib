// -*-C++-*-

#ifndef VEC_DOUBLE_SSE2_H
#define VEC_DOUBLE_SSE2_H

#include "floatprops.h"
#include "mathfuncs.h"
#include "vec_base.h"

#include <cmath>

// SSE2 intrinsics
#include <emmintrin.h>
#ifdef __SSE4_1__               // Intel's SSE 4.1
#  include <smmintrin.h>
#endif
#ifdef __SSE4A__                // AMD's SSE 4a
#  include <ammintrin.h>
#endif



namespace vecmathlib {
  
#define VECMATHLIB_HAVE_VEC_DOUBLE_2
  template<> struct boolvec<double,2>;
  template<> struct intvec<double,2>;
  template<> struct realvec<double,2>;
  
  
  
  template<>
  struct boolvec<double,2>: floatprops<double>
  {
    static int const size = 2;
    typedef bool scalar_t;
    typedef __m128d bvector_t;
    
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
    // Can't have a non-trivial copy constructor; if so, objects won't
    // be passed in registers
    // boolvec(boolvec const& x): v(x.v) {}
    // boolvec& operator=(boolvec const& x) { return v=x.v, *this; }
    boolvec(bvector_t x): v(x) {}
    boolvec(bool a):
    v(_mm_castsi128_pd(_mm_set1_epi64x(from_bool(a)))) {}
    boolvec(bool const* as):
    v(_mm_castsi128_pd(_mm_set_epi64x(from_bool(as[1]), from_bool(as[0])))) {}
    
    operator bvector_t() const { return v; }
    bool operator[](int n) const
    {
      // return to_bool(((uint_t const*)&v)[n]);
      boolvec x = *this;
      switch (n){
      case 0: /* do nothing */ break;
      case 1: x = _mm_shuffle_pd(x.v, x.v, _MM_SHUFFLE2(0,1)); break;
      default: assert(0);
      }
      // return to_bool(FP::as_int(_mm_cvtsd_f64(x.v)));
      return to_bool(_mm_cvtsi128_si64(_mm_castpd_si128(x.v)));
    }
    boolvec& set_elt(int n, bool a)
    {
      return ((int_t*)&v)[n] = from_bool(a), *this;
    }
    
    
    
    intvec_t as_int() const;      // defined after intvec
    intvec_t convert_int() const; // defined after intvec
    
    
    
    boolvec operator!() const { return _mm_xor_pd(boolvec(true), v); }
    
    boolvec operator&&(boolvec x) const { return _mm_and_pd(v, x.v); }
    boolvec operator||(boolvec x) const { return _mm_or_pd(v, x.v); }
    boolvec operator==(boolvec x) const { return !(*this==x); }
    boolvec operator!=(boolvec x) const { return _mm_xor_pd(v, x.v); }
    
    bool all() const;
    bool any() const;
    
    
    
    // ifthen(condition, then-value, else-value)
    intvec_t ifthen(intvec_t x, intvec_t y) const; // defined after intvec
    realvec_t ifthen(realvec_t x, realvec_t y) const; // defined after realvec
  };
  
  
  
  template<>
  struct intvec<double,2>: floatprops<double>
  {
    static int const size = 2;
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
    // Can't have a non-trivial copy constructor; if so, objects won't
    // be passed in registers
    // intvec(intvec const& x): v(x.v) {}
    // intvec& operator=(intvec const& x) { return v=x.v, *this; }
    intvec(ivector_t x): v(x) {}
    intvec(int_t a): v(_mm_set1_epi64x(a)) {}
    intvec(int_t const* as): v(_mm_set_epi64x(as[1], as[0])) {}
    
    operator ivector_t() const { return v; }
    int_t operator[](int n) const { return ((int_t const*)&v)[n]; }
    intvec& set_elt(int n, int_t a) { return ((int_t*)&v)[n]=a, *this; }
    
    
    
    boolvec_t as_bool() const { return _mm_castsi128_pd(v); }
    boolvec_t convert_bool() const
    {
      // Result: convert_bool(0)=false, convert_bool(else)=true
      // There is no intrinsic to compare with zero. Instead, we check
      // whether x is positive and x-1 is negative.
      intvec x = *this;
      intvec xm1 = x - 1;
      // We know that boolvec values depend only on the sign bit
      // return (~xm1 | x).as_bool();
      return x.as_bool() || !xm1.as_bool();
    }
    realvec_t as_float() const;      // defined after realvec
    realvec_t convert_float() const; // defined after realvec
    
    
    
    // Note: not all arithmetic operations are supported!
    
    intvec operator+() const { return *this; }
    intvec operator-() const { return IV(I(0)) - *this; }
    
    intvec operator+(intvec x) const { return _mm_add_epi64(v, x.v); }
    intvec operator-(intvec x) const { return _mm_sub_epi64(v, x.v); }
    
    intvec& operator+=(intvec const& x) { return *this=*this+x; }
    intvec& operator-=(intvec const& x) { return *this=*this-x; }
    
    
    
    intvec operator~() const { return IV(~U(0)) ^ *this; }
    
    intvec operator&(intvec x) const
    {
      return _mm_castpd_si128(_mm_and_pd(_mm_castsi128_pd(v),
                                         _mm_castsi128_pd(x.v)));
    }
    intvec operator|(intvec x) const
    {
      return _mm_castpd_si128(_mm_or_pd(_mm_castsi128_pd(v),
                                        _mm_castsi128_pd(x.v)));
    }
    intvec operator^(intvec x) const
    {
      return _mm_castpd_si128(_mm_xor_pd(_mm_castsi128_pd(v),
                                         _mm_castsi128_pd(x.v)));
    }
    
    intvec& operator&=(intvec const& x) { return *this=*this&x; }
    intvec& operator|=(intvec const& x) { return *this=*this|x; }
    intvec& operator^=(intvec const& x) { return *this=*this^x; }
    
    
    
    intvec lsr(int_t n) const { return _mm_srli_epi64(v, n); }
    intvec operator>>(int_t n) const
    {
      // There is no _mm_srai_epi64. To emulate it, add 0x80000000
      // before shifting, and subtract the shifted 0x80000000 after
      // shifting
      intvec x = *this;
      // Convert signed to unsiged
      x += U(1) << (bits-1);
      // Shift
      x = x.lsr(x);
      // Undo conversion
      x -= U(1) << (bits-n);
      return x;
    }
    intvec operator<<(int_t n) const { return _mm_slli_epi64(v, n); }
    intvec& operator>>=(int_t n) { return *this=*this>>n; }
    intvec& operator<<=(int_t n) { return *this=*this<<n; }
    
    intvec lsr(intvec n) const { return _mm_srl_epi64(v, n.v); }
    intvec operator>>(intvec n) const
    {
      // There is no _mm_srai_epi64. To emulate it, add 0x80000000
      // before shifting, and subtract the shifted 0x80000000 after
      // shifting
      intvec x = *this;
      // Convert signed to unsiged
      x += U(1) << (bits-1);
      // Shift
      x = x.lsr(n);
      // Undo conversion
      x -= IV(1) << (IV(bits) - n);
      return x;
    }
    intvec operator<<(intvec n) const { return _mm_sll_epi64(v, n.v); }
    intvec& operator>>=(intvec n) { return *this=*this>>n; }
    intvec& operator<<=(intvec n) { return *this=*this<<n; }
  };
  
  
  
  template<>
  struct realvec<double,2>: floatprops<double>
  {
    static int const size = 2;
    typedef real_t scalar_t;
    typedef __m128d vector_t;
    
    static constexpr char const* const name = "<SSE2:2*double>";
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
    // Can't have a non-trivial copy constructor; if so, objects won't
    // be passed in registers
    // realvec(realvec const& x): v(x.v) {}
    // realvec& operator=(realvec const& x) { return v=x.v, *this; }
    realvec(vector_t x): v(x) {}
    realvec(real_t a): v(_mm_set1_pd(a)) {}
    realvec(real_t const* as): v(_mm_set_pd(as[1], as[0])) {}
    
    operator vector_t() const { return v; }
    real_t operator[](int n) const
    {
      // return ((real_t const*)&v)[n];
      switch (n){
      case 0: return _mm_cvtsd_f64(v);
      case 1: return _mm_cvtsd_f64(_mm_shuffle_pd(v, v, _MM_SHUFFLE2(0,1)));
      }
      assert(0);
    }
    realvec& set_elt(int n, real_t a) { return ((real_t*)&v)[n]=a, *this; }
    
    
    
    intvec_t as_int() const { return _mm_castpd_si128(v); }
    intvec_t convert_int() const { return MF::vml_convert_int(*this); }
    
    
    
    realvec operator+() const { return *this; }
    realvec operator-() const { return RV(0.0) - *this; }
    
    realvec operator+(realvec x) const { return _mm_add_pd(v, x.v); }
    realvec operator-(realvec x) const { return _mm_sub_pd(v, x.v); }
    realvec operator*(realvec x) const { return _mm_mul_pd(v, x.v); }
    realvec operator/(realvec x) const { return _mm_div_pd(v, x.v); }
    
    realvec& operator+=(realvec const& x) { return *this=*this+x; }
    realvec& operator-=(realvec const& x) { return *this=*this-x; }
    realvec& operator*=(realvec const& x) { return *this=*this*x; }
    realvec& operator/=(realvec const& x) { return *this=*this/x; }
    
    real_t prod() const
    {
      return (*this)[0] * (*this)[1];
    }
    real_t sum() const
    {
      // return (*this)[0] + (*this)[1];
      return _mm_cvtsd_f64(_mm_hadd_pd(v, v));
    }
    
    
    
    boolvec_t operator==(realvec const& x) const
    {
      return _mm_cmpeq_pd(v, x.v);
    }
    boolvec_t operator!=(realvec const& x) const
    {
      return _mm_cmpneq_pd(v, x.v);
    }
    boolvec_t operator<(realvec const& x) const
    {
      return _mm_cmplt_pd(v, x.v);
    }
    boolvec_t operator<=(realvec const& x) const
    {
      return _mm_cmple_pd(v, x.v);
    }
    boolvec_t operator>(realvec const& x) const
    {
      return _mm_cmpgt_pd(v, x.v);
    }
    boolvec_t operator>=(realvec const& x) const
    {
      return _mm_cmpge_pd(v, x.v);
    }
    
    
    
    realvec acos() const { return MF::vml_acos(*this); }
    realvec acosh() const { return MF::vml_acosh(*this); }
    realvec asin() const { return MF::vml_asin(*this); }
    realvec asinh() const { return MF::vml_asinh(*this); }
    realvec atan() const { return MF::vml_atan(*this); }
    realvec atanh() const { return MF::vml_atanh(*this); }
    realvec ceil() const { return _mm_ceil_pd(v); }
    realvec copysign(realvec y) const { return MF::vml_copysign(*this, y); }
    realvec cos() const { return MF::vml_cos(*this); }
    realvec cosh() const { return MF::vml_cosh(*this); }
    realvec exp() const { return MF::vml_exp(*this); }
    realvec exp10() const { return MF::vml_exp10(*this); }
    realvec exp2() const { return MF::vml_exp2(*this); }
    realvec expm1() const { return MF::vml_expm1(*this); }
    realvec fabs() const { return MF::vml_fabs(*this); }
    realvec floor() const { return _mm_floor_pd(v); }
    realvec fmod(realvec y) const { return MF::vml_fmod(*this, y); }
    intvec_t ilogb() const { return MF::vml_ilogb(*this); }
    realvec log() const { return MF::vml_log(*this); }
    realvec log10() const { return MF::vml_log10(*this); }
    realvec log1p() const { return MF::vml_log1p(*this); }
    realvec log2() const { return MF::vml_log2(*this); }
    realvec pow(realvec y) const { return MF::vml_pow(*this, y); }
    realvec rcp() const { return _mm_div_pd(_mm_set1_pd(1.0), v); }
    realvec remainder(realvec y) const { return MF::vml_remainder(*this, y); }
    realvec round() const { return _mm_round_pd(v, _MM_FROUND_NINT); }
    realvec rsqrt() const { return MF::vml_rsqrt(*this); }
    realvec scalbn(intvec_t n) const { return MF::vml_scalbn(*this, n); }
    boolvec_t signbit() const { return v; }
    realvec sin() const { return MF::vml_sin(*this); }
    realvec sinh() const { return MF::vml_sinh(*this); }
    realvec sqrt() const { return _mm_sqrt_pd(v); }
    realvec tan() const { return MF::vml_tan(*this); }
    realvec tanh() const { return MF::vml_tanh(*this); }
  };
  
  
  
  // boolvec definitions
  
  inline
  auto boolvec<double,2>::as_int() const -> intvec_t
  {
    return _mm_castpd_si128(v);
  }
  
  inline
  auto boolvec<double,2>::convert_int() const -> intvec_t
  {
    //return ifthen(v, U(1), U(0));
    return lsr(as_int(), bits-1);
  }
  
  inline
  bool boolvec<double,2>::all() const
  {
    // return (*this)[0] && (*this)[1];
    boolvec x = *this;
    x = x || _mm_shuffle_pd(x.v, x.v, _MM_SHUFFLE2(0,1));
    return x[0];
  }
  
  inline
  bool boolvec<double,2>::any() const
  {
    // return (*this)[0] || (*this)[1];
    boolvec x = *this;
    x = x && _mm_shuffle_pd(x.v, x.v, _MM_SHUFFLE2(0,1));
    return x[0];
  }
  
  inline
  auto boolvec<double,2>::ifthen(intvec_t x, intvec_t y) const -> intvec_t
  {
    return ifthen(x.as_float(), y.as_float()).as_int();
  }
  
  inline
  auto boolvec<double,2>::ifthen(realvec_t x, realvec_t y) const -> realvec_t
  {
    return _mm_blendv_pd(y.v, x.v, v);
  }
  
  
  
  // intvec definitions
  
  inline auto intvec<double,2>::as_float() const -> realvec_t
  {
    return _mm_castsi128_pd(v);
  }
  
  inline auto intvec<double,2>::convert_float() const -> realvec_t
  {
    return MF::vml_convert_float(*this);
  }
  
} // namespace vecmathlib

#endif  // #ifndef VEC_DOUBLE_SSE2_H
