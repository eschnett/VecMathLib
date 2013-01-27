// -*-C++-*-

#ifndef VEC_FLOAT_AVX_H
#define VEC_FLOAT_AVX_H

#include "floatprops.h"
#include "mathfuncs.h"
#include "vec_base.h"

#include <cmath>

// AVX intrinsics
#include <immintrin.h>



namespace vecmathlib {
  
#define VECMATHLIB_HAVE_VEC_FLOAT_8
  template<> struct boolvec<float,8>;
  template<> struct intvec<float,8>;
  template<> struct realvec<float,8>;
  
  
  
  template<>
  struct boolvec<float,8>: floatprops<float>
  {
    static int const size = 8;
    typedef bool scalar_t;
    typedef __m256 bvector_t;
    
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
    v(_mm256_castsi256_ps(_mm256_set1_epi32(from_bool(a)))) {}
    boolvec(bool const* as):
    v(_mm256_castsi256_ps(_mm256_set_epi32(from_bool(as[7]),
                                           from_bool(as[6]),
                                           from_bool(as[5]),
                                           from_bool(as[4]),
                                           from_bool(as[3]),
                                           from_bool(as[2]),
                                           from_bool(as[1]),
                                           from_bool(as[0])))) {}
    
    operator bvector_t() const { return v; }
    bool operator[](int n) const
    {
      // return to_bool(((uint_t const*)&v)[n]);
      __m128 x =
        n & 4 ? _mm256_extractf128_ps(v, 1) : _mm256_castps256_ps128(v);
      switch (n & 3){
      case 0: /* do nothing */ break;
      case 1: x = _mm_shuffle_ps(x, x, _MM_SHUFFLE(2,3,0,1)); break;
      case 2: x = _mm_shuffle_ps(x, x, _MM_SHUFFLE(1,0,3,2)); break;
      case 3: x = _mm_shuffle_ps(x, x, _MM_SHUFFLE(0,1,2,3)); break;
      default: assert(0);
      }
      // return to_bool(FP::as_int(_mm_cvtss_f32(x)));
      return to_bool(_mm_cvtsi128_si32(_mm_castps_si128(x)));
    }
    boolvec& set_elt(int n, bool a)
    {
      return ((int_t*)&v)[n] = from_bool(a), *this;
    }
    
    
    
    intvec_t as_int() const;      // defined after intvec
    intvec_t convert_int() const; // defined after intvec
    
    
    
    boolvec operator!() const { return _mm256_xor_ps(boolvec(true), v); }
    
    boolvec operator&&(boolvec x) const { return _mm256_and_ps(v, x.v); }
    boolvec operator||(boolvec x) const { return _mm256_or_ps(v, x.v); }
    boolvec operator==(boolvec x) const { return !(*this==x); }
    boolvec operator!=(boolvec x) const { return _mm256_xor_ps(v, x.v); }
    
    bool all() const;
    bool any() const;
    
    
    
    // ifthen(condition, then-value, else-value)
    intvec_t ifthen(intvec_t x, intvec_t y) const; // defined after intvec
    realvec_t ifthen(realvec_t x, realvec_t y) const; // defined after realvec
  };
  
  
  
  template<>
  struct intvec<float,8>: floatprops<float>
  {
    static int const size = 8;
    typedef int_t scalar_t;
    typedef __m256i ivector_t;
    
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
    intvec(int_t a): v(_mm256_set1_epi32(a)) {}
    intvec(int_t const* as): v(_mm256_set_epi32(as[7], as[6], as[5], as[4],
                                                as[3], as[2], as[1], as[0])) {}
    
    operator ivector_t() const { return v; }
    int_t operator[](int n) const
    {
      // return ((int_t const*)&v)[n];
      __m128i x =
        n & 4 ? _mm256_extractf128_si256(v, 1) : _mm256_castsi256_si128(v);
      switch (n & 3){
      case 0: /* do nothing */ break;
      case 1: x = _mm_shuffle_epi32(x, _MM_SHUFFLE(2,3,0,1)); break;
      case 2: x = _mm_shuffle_epi32(x, _MM_SHUFFLE(1,0,3,2)); break;
      case 3: x = _mm_shuffle_epi32(x, _MM_SHUFFLE(0,1,2,3)); break;
      default: assert(0);
      }
      return _mm_cvtsi128_si32(x);
    }
    intvec& set_elt(int n, int_t a) { return ((int_t*)&v)[n]=a, *this; }
    
    
    
    boolvec_t as_bool() const { return _mm256_castsi256_ps(v); }
    boolvec_t convert_bool() const
    {
      // Result: convert_bool(0)=false, convert_bool(else)=true
      // There is no intrinsic to compare with zero. Instead, we check
      // whether x is positive and x-1 is negative.
      intvec x = *this;
      // We know that boolvec values depend only on the sign bit
      // return (~(x-1) | x).as_bool();
      // return x.as_bool() || !(x-1).as_bool();
      return x.as_bool() || (x + (FP::signbit_mask - 1)).as_bool();
    }
    realvec_t as_float() const;      // defined after realvec
    realvec_t convert_float() const; // defined after realvec
    
    
    
    // Note: not all arithmetic operations are supported!
    
    intvec operator+() const { return *this; }
    intvec operator-() const { return IV(0) - *this; }
    
    intvec operator+(intvec x) const
    {
      __m128i vlo = _mm256_castsi256_si128(v);
      __m128i vhi = _mm256_extractf128_si256(v, 1);
      __m128i xvlo = _mm256_castsi256_si128(x.v);
      __m128i xvhi = _mm256_extractf128_si256(x.v, 1);
      vlo = _mm_add_epi32(vlo, xvlo);
      vhi = _mm_add_epi32(vhi, xvhi);
      return _mm256_insertf128_si256(_mm256_castsi128_si256(vlo), vhi, 1);
    }
    intvec operator-(intvec x) const
    {
      __m128i vlo = _mm256_castsi256_si128(v);
      __m128i vhi = _mm256_extractf128_si256(v, 1);
      __m128i xvlo = _mm256_castsi256_si128(x.v);
      __m128i xvhi = _mm256_extractf128_si256(x.v, 1);
      vlo = _mm_sub_epi32(vlo, xvlo);
      vhi = _mm_sub_epi32(vhi, xvhi);
      return _mm256_insertf128_si256(_mm256_castsi128_si256(vlo), vhi, 1);
    }
    
    intvec& operator+=(intvec const& x) { return *this=*this+x; }
    intvec& operator-=(intvec const& x) { return *this=*this-x; }
    
    
    
    intvec operator~() const { return IV(~U(0)) ^ *this; }
    
    intvec operator&(intvec x) const
    {
      return _mm256_castps_si256(_mm256_and_ps(_mm256_castsi256_ps(v),
                                               _mm256_castsi256_ps(x.v)));
    }
    intvec operator|(intvec x) const
    {
      return _mm256_castps_si256(_mm256_or_ps(_mm256_castsi256_ps(v),
                                              _mm256_castsi256_ps(x.v)));
    }
    intvec operator^(intvec x) const
    {
      return _mm256_castps_si256(_mm256_xor_ps(_mm256_castsi256_ps(v),
                                               _mm256_castsi256_ps(x.v)));
    }
    
    intvec& operator&=(intvec const& x) { return *this=*this&x; }
    intvec& operator|=(intvec const& x) { return *this=*this|x; }
    intvec& operator^=(intvec const& x) { return *this=*this^x; }
    
    
    
    intvec lsr(int_t n) const
    {
      __m128i vlo = _mm256_castsi256_si128(v);
      __m128i vhi = _mm256_extractf128_si256(v, 1);
      vlo = _mm_srli_epi32(vlo, n);
      vhi = _mm_srli_epi32(vhi, n);
      return _mm256_insertf128_si256(_mm256_castsi128_si256(vlo), vhi, 1);
    }
    intvec operator>>(int_t n) const
    {
      __m128i vlo = _mm256_castsi256_si128(v);
      __m128i vhi = _mm256_extractf128_si256(v, 1);
      // There is no _mm_srai_epi32. To emulate it, add 0x80000000
      // before shifting, and subtract the shifted 0x80000000 after
      // shifting
      // Convert signed to unsiged
      vlo += _mm_set1_epi32(U(1) << (bits-1));
      vhi += _mm_set1_epi32(U(1) << (bits-1));
      // Shift
      vlo = _mm_srli_epi32(vlo, n);
      vhi = _mm_srli_epi32(vhi, n);
      // Undo conversion
      vlo -= _mm_set1_epi32(U(1) << (bits-n));
      vhi -= _mm_set1_epi32(U(1) << (bits-n));
      return _mm256_insertf128_si256(_mm256_castsi128_si256(vlo), vhi, 1);
    }
    intvec operator<<(int_t n) const
    {
      __m128i vlo = _mm256_castsi256_si128(v);
      __m128i vhi = _mm256_extractf128_si256(v, 1);
      vlo = _mm_slli_epi32(vlo, n);
      vhi = _mm_slli_epi32(vhi, n);
      return _mm256_insertf128_si256(_mm256_castsi128_si256(vlo), vhi, 1);
    }
    intvec& operator>>=(int_t n) { return *this=*this>>n; }
    intvec& operator<<=(int_t n) { return *this=*this<<n; }
    
    intvec lsr(intvec n) const
    {
      __m128i vlo = _mm256_castsi256_si128(v);
      __m128i vhi = _mm256_extractf128_si256(v, 1);
      __m128i nvlo = _mm256_castsi256_si128(n.v);
      __m128i nvhi = _mm256_extractf128_si256(n.v, 1);
      vlo = _mm_srl_epi32(vlo, nvlo);
      vhi = _mm_srl_epi32(vhi, nvhi);
      return _mm256_insertf128_si256(_mm256_castsi128_si256(vlo), vhi, 1);
    }
    intvec operator>>(intvec n) const
    {
      __m128i vlo = _mm256_castsi256_si128(v);
      __m128i vhi = _mm256_extractf128_si256(v, 1);
      __m128i nvlo = _mm256_castsi256_si128(n.v);
      __m128i nvhi = _mm256_extractf128_si256(n.v, 1);
      // Convert signed to unsiged
      vlo += _mm_set1_epi32(U(1) << (bits-1));
      vhi += _mm_set1_epi32(U(1) << (bits-1));
      // Shift
      vlo = _mm_srl_epi32(vlo, nvlo);
      vhi = _mm_srl_epi32(vhi, nvhi);
      // Undo conversion
      vlo -= _mm_sll_epi32(_mm_set1_epi32(1),
                           _mm_sub_epi32(_mm_set1_epi32(bits), nvlo));
      vhi -= _mm_sll_epi32(_mm_set1_epi32(1),
                           _mm_sub_epi32(_mm_set1_epi32(bits), nvhi));
      return _mm256_insertf128_si256(_mm256_castsi128_si256(vlo), vhi, 1);
    }
    intvec operator<<(intvec n) const
    {
      __m128i vlo = _mm256_castsi256_si128(v);
      __m128i vhi = _mm256_extractf128_si256(v, 1);
      __m128i nvlo = _mm256_castsi256_si128(n.v);
      __m128i nvhi = _mm256_extractf128_si256(n.v, 1);
      vlo = _mm_sll_epi32(vlo, nvlo);
      vhi = _mm_sll_epi32(vhi, nvhi);
      return _mm256_insertf128_si256(_mm256_castsi128_si256(vlo), vhi, 1);
    }
    intvec& operator>>=(intvec n) { return *this=*this>>n; }
    intvec& operator<<=(intvec n) { return *this=*this<<n; }
    
    
    
    boolvec_t operator==(intvec const& x) const
    {
      return ! (*this != x);
    }
    boolvec_t operator!=(intvec const& x) const
    {
      return (*this ^ x).convert_bool();
    }
  };
  
  
  
  template<>
  struct realvec<float,8>: floatprops<float>
  {
    static int const size = 8;
    typedef real_t scalar_t;
    typedef __m256 vector_t;
    
    static constexpr char const* const name = "<AVX:8*float>";
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
    realvec(real_t a): v(_mm256_set1_ps(a)) {}
    realvec(real_t const* as): v(_mm256_set_ps(as[7], as[6], as[5], as[4],
                                               as[3], as[2], as[1], as[0])) {}
    
    operator vector_t() const { return v; }
    real_t operator[](int n) const
    {
      // return ((real_t const*)&v)[n];
      __m128 x =
        n & 4 ? _mm256_extractf128_ps(v, 1) : _mm256_castps256_ps128(v);
      switch (n & 3){
      case 0: /* do nothing */ break;
      case 1: x = _mm_shuffle_ps(x, x, _MM_SHUFFLE(2,3,0,1)); break;
      case 2: x = _mm_shuffle_ps(x, x, _MM_SHUFFLE(1,0,3,2)); break;
      case 3: x = _mm_shuffle_ps(x, x, _MM_SHUFFLE(0,1,2,3)); break;
      default: assert(0);
      }
      return _mm_cvtss_f32(x);
    }
    realvec& set_elt(int n, real_t a) { return ((real_t*)&v)[n]=a, *this; }
    
    
    
    intvec_t as_int() const { return _mm256_castps_si256(v); }
    intvec_t convert_int() const { return _mm256_cvtps_epi32(v); }
    
    
    
    realvec operator+() const { return *this; }
    realvec operator-() const { return RV(0.0) - *this; }
    
    realvec operator+(realvec x) const { return _mm256_add_ps(v, x.v); }
    realvec operator-(realvec x) const { return _mm256_sub_ps(v, x.v); }
    realvec operator*(realvec x) const { return _mm256_mul_ps(v, x.v); }
    realvec operator/(realvec x) const { return _mm256_div_ps(v, x.v); }
    
    realvec& operator+=(realvec const& x) { return *this=*this+x; }
    realvec& operator-=(realvec const& x) { return *this=*this-x; }
    realvec& operator*=(realvec const& x) { return *this=*this*x; }
    realvec& operator/=(realvec const& x) { return *this=*this/x; }
    
    real_t prod() const
    {
      return
        (*this)[0] * (*this)[1] * (*this)[2] * (*this)[3] *
        (*this)[4] * (*this)[5] * (*this)[6] * (*this)[7];
    }
    real_t sum() const
    {
      // return
      //   (*this)[0] + (*this)[1] + (*this)[2] + (*this)[3] +
      //   (*this)[4] + (*this)[5] + (*this)[6] + (*this)[7];
      // _m256 x = vhaddps(v, v);
      // x = vhaddps(x, x);
      // __m128 xlo = _mm256_extractf128_ps(x, 0);
      // __m128 xhi = _mm256_extractf128_ps(x, 1);
      // return _mm_cvtsd_f64(xlo) + _mm_cvtsd_f64(xhi);
      realvec x = *this;
      x = _mm256_hadd_ps(x.v, x.v);
      x = _mm256_hadd_ps(x.v, x.v);
      return x[0] + x[4];
    }
    
    
    
    boolvec_t operator==(realvec const& x) const
    {
      return _mm256_cmp_ps(v, x.v, _CMP_EQ_OQ);
    }
    boolvec_t operator!=(realvec const& x) const
    {
      return _mm256_cmp_ps(v, x.v, _CMP_NEQ_OQ);
    }
    boolvec_t operator<(realvec const& x) const
    {
      return _mm256_cmp_ps(v, x.v, _CMP_LT_OQ);
    }
    boolvec_t operator<=(realvec const& x) const
    {
      return _mm256_cmp_ps(v, x.v, _CMP_LE_OQ);
    }
    boolvec_t operator>(realvec const& x) const
    {
      return _mm256_cmp_ps(v, x.v, _CMP_GT_OQ);
    }
    boolvec_t operator>=(realvec const& x) const
    {
      return _mm256_cmp_ps(v, x.v, _CMP_GE_OQ);
    }
    
    
    
    realvec acos() const { return MF::vml_acos(*this); }
    realvec acosh() const { return MF::vml_acosh(*this); }
    realvec asin() const { return MF::vml_asin(*this); }
    realvec asinh() const { return MF::vml_asinh(*this); }
    realvec atan() const { return MF::vml_atan(*this); }
    realvec atanh() const { return MF::vml_atanh(*this); }
    realvec ceil() const { return _mm256_ceil_ps(v); }
    realvec copysign(realvec y) const { return MF::vml_copysign(*this, y); }
    realvec cos() const { return MF::vml_cos(*this); }
    realvec cosh() const { return MF::vml_cosh(*this); }
    realvec exp() const { return MF::vml_exp(*this); }
    realvec exp10() const { return MF::vml_exp10(*this); }
    realvec exp2() const { return MF::vml_exp2(*this); }
    realvec expm1() const { return MF::vml_expm1(*this); }
    realvec fabs() const { return MF::vml_fabs(*this); }
    realvec fdim(realvec y) const { return MF::vml_fdim(*this, y); }
    realvec floor() const { return _mm256_floor_ps(v); }
    realvec fma(realvec y, realvec z) const { return MF::vml_fma(*this, y, z); }
    realvec fmax(realvec y) const { return _mm256_max_ps(v, y.v); }
    realvec fmin(realvec y) const { return _mm256_min_ps(v, y.v); }
    realvec fmod(realvec y) const { return MF::vml_fmod(*this, y); }
    intvec_t ilogb() const { return MF::vml_ilogb(*this); }
    boolvec_t isfinite() const { return MF::vml_isfinite(*this); }
    boolvec_t isinf() const { return MF::vml_isinf(*this); }
    boolvec_t isnan() const { return MF::vml_isnan(*this); }
    boolvec_t isnormal() const { return MF::vml_isnormal(*this); }
    realvec log() const { return MF::vml_log(*this); }
    realvec log10() const { return MF::vml_log10(*this); }
    realvec log1p() const { return MF::vml_log1p(*this); }
    realvec log2() const { return MF::vml_log2(*this); }
    realvec pow(realvec y) const { return MF::vml_pow(*this, y); }
    realvec rcp() const
    {
      realvec x = *this;
      realvec r = _mm256_rcp_ps(x); // this is only an approximation
      r *= RV(2.0) - r*x;        // one Newton iteration (see vml_rcp)
      return r;
    }
    realvec remainder(realvec y) const { return MF::vml_remainder(*this, y); }
    realvec round() const { return _mm256_round_ps(v, _MM_FROUND_NINT); }
    realvec rsqrt() const
    {
      realvec x = *this;
      realvec r = _mm256_rsqrt_ps(x);    // this is only an approximation
      r *= RV(1.5) - RV(0.5)*x * r*r; // one Newton iteration (see vml_rsqrt)
      return r;
    }
    realvec scalbn(intvec_t n) const { return MF::vml_scalbn(*this, n); }
    boolvec_t signbit() const { return v; }
    realvec sin() const { return MF::vml_sin(*this); }
    realvec sinh() const { return MF::vml_sinh(*this); }
    realvec sqrt() const { return _mm256_sqrt_ps(v); }
    realvec tan() const { return MF::vml_tan(*this); }
    realvec tanh() const { return MF::vml_tanh(*this); }
  };
  
  
  
  // boolvec definitions
  
  inline
  auto boolvec<float,8>::as_int() const -> intvec_t
  {
    return _mm256_castps_si256(v);
  }
  
  inline
  auto boolvec<float,8>::convert_int() const -> intvec_t
  {
    return lsr(as_int(), bits-1);
  }
  
  inline
  bool boolvec<float,8>::all() const
  {
    // return
    //   (*this)[0] && (*this)[1] && (*this)[2] && (*this)[3] &&
    //   (*this)[4] && (*this)[5] && (*this)[6] && (*this)[7];
    boolvec x = *this;
    x = x && _mm256_shuffle_ps(x.v, x.v, _MM_SHUFFLE(1,0,3,2));
    x = x && _mm256_shuffle_ps(x.v, x.v, _MM_SHUFFLE(2,3,0,1));
    __m128 y = _mm_and_ps(_mm256_castps256_ps128(x.v),
                          _mm256_extractf128_ps(x.v, 1));
    return to_bool(_mm_cvtsi128_si32(_mm_castps_si128(y)));
  }
  
  inline
  bool boolvec<float,8>::any() const
  {
    // return
    //   (*this)[0] || (*this)[1] || (*this)[2] || (*this)[3] ||
    //   (*this)[4] || (*this)[5] || (*this)[6] || (*this)[7];
    boolvec x = *this;
    x = x || _mm256_shuffle_ps(x.v, x.v, _MM_SHUFFLE(1,0,3,2));
    x = x || _mm256_shuffle_ps(x.v, x.v, _MM_SHUFFLE(2,3,0,1));
    __m128 y = _mm_or_ps(_mm256_castps256_ps128(x.v),
                         _mm256_extractf128_ps(x.v, 1));
    return to_bool(_mm_cvtsi128_si32(_mm_castps_si128(y)));
  }
  
  inline
  auto boolvec<float,8>::ifthen(intvec_t x, intvec_t y) const -> intvec_t
  {
    return ifthen(x.as_float(), y.as_float()).as_int();
  }
  
  inline
  auto boolvec<float,8>::ifthen(realvec_t x, realvec_t y) const -> realvec_t
  {
    return _mm256_blendv_ps(y.v, x.v, v);
  }

  
  
  // intvec definitions
  
  inline auto intvec<float,8>::as_float() const -> realvec_t
  {
    return _mm256_castsi256_ps(v);
  }
  
  inline auto intvec<float,8>::convert_float() const -> realvec_t
  {
    return _mm256_cvtepi32_ps(v);
  }
  
} // namespace vecmathlib

#endif  // #ifndef VEC_FLOAT_AVX_H