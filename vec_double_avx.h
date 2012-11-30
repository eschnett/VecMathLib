// -*-C++-*-

#ifndef VEC_DOUBLE_AVX_H
#define VEC_DOUBLE_AVX_H

#include "floatprops.h"
#include "mathfuncs.h"
#include "vec_base.h"

#include <cmath>

// AVX intrinsics
#include <immintrin.h>



namespace vecmathlib {
  
  template<> struct boolvec<double,4>;
  template<> struct intvec<double,4>;
  template<> struct realvec<double,4>;
  
  
  
  template<>
  struct boolvec<double,4>: floatprops<double>
  {
    static const int size = 4;
    typedef bool scalar_t;
    typedef __m256d bvector_t;
    
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
    v(_mm256_castsi256_pd(_mm256_set1_epi64x(from_bool(a)))) {}
    boolvec(bool const* as):
    v(_mm256_castsi256_pd(_mm256_set_epi64x(from_bool(as[3]),
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
    
    
    
    boolvec operator!() const { return _mm256_xor_pd(boolvec(true), v); }
    
    boolvec operator&&(boolvec x) const { return _mm256_and_pd(v, x.v); }
    boolvec operator||(boolvec x) const { return _mm256_or_pd(v, x.v); }
    boolvec operator==(boolvec x) const { return !(*this==x); }
    boolvec operator!=(boolvec x) const { return _mm256_xor_pd(v, x.v); }
    
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
  struct intvec<double,4>: floatprops<double>
  {
    static const int size = 4;
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
    intvec(intvec const& x): v(x.v) {}
    intvec& operator=(intvec const& x) { return v=x.v, *this; }
    intvec(ivector_t x): v(x) {}
    intvec(int_t a): v(_mm256_set1_epi64x(a)) {}
    intvec(int_t const* as): v(_mm256_set_epi64x(as[3], as[2], as[1], as[0])) {}
    
    operator ivector_t() const { return v; }
    int_t operator[](int n) const { return ((int_t const*)&v)[n]; }
    intvec& set_elt(int n, int_t a) { return ((int_t*)&v)[n]=a, *this; }
    
    
    
    auto as_bool() const -> boolvec_t { return _mm256_castsi256_pd(v); }
    auto convert_bool() const -> boolvec_t
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
    auto as_float() const -> realvec_t; // defined after realvec
    auto convert_float() const -> realvec_t; // defined after realvec
    
    
    
    // Note: not all arithmetic operations are supported!
    
    intvec operator+() const { return *this; }
    intvec operator-() const
    {
      __m128i v01 = _mm256_castsi256_si128(v);
      __m128i v23 = _mm256_extractf128_si256(v, 1);
      v01 = _mm_sub_epi64(_mm_set1_epi64x(0), v01);
      v23 = _mm_sub_epi64(_mm_set1_epi64x(0), v23);
      return _mm256_insertf128_si256(_mm256_castsi128_si256(v01), v23, 1);
    }
    
    intvec operator+(intvec x) const
    {
      __m128i v01 = _mm256_castsi256_si128(v);
      __m128i v23 = _mm256_extractf128_si256(v, 1);
      __m128i xv01 = _mm256_castsi256_si128(x.v);
      __m128i xv23 = _mm256_extractf128_si256(x.v, 1);
      v01 = _mm_add_epi64(v01, xv01);
      v23 = _mm_add_epi64(v23, xv23);
      return _mm256_insertf128_si256(_mm256_castsi128_si256(v01), v23, 1);
    }
    intvec operator-(intvec x) const
    {
      __m128i v01 = _mm256_castsi256_si128(v);
      __m128i v23 = _mm256_extractf128_si256(v, 1);
      __m128i xv01 = _mm256_castsi256_si128(x.v);
      __m128i xv23 = _mm256_extractf128_si256(x.v, 1);
      v01 = _mm_sub_epi64(v01, xv01);
      v23 = _mm_sub_epi64(v23, xv23);
      return _mm256_insertf128_si256(_mm256_castsi128_si256(v01), v23, 1);
    }
    
    intvec& operator+=(intvec const& x) { return *this=*this+x; }
    intvec& operator-=(intvec const& x) { return *this=*this-x; }
    
    
    
    intvec operator~() const
    {
      return _mm256_castpd_si256(_mm256_xor_pd(_mm256_castsi256_pd(IV(~U(0))),
                                               _mm256_castsi256_pd(v)));
    }
    
    intvec operator&(intvec x) const
    {
      return _mm256_castpd_si256(_mm256_and_pd(_mm256_castsi256_pd(v),
                                               _mm256_castsi256_pd(x.v)));
    }
    intvec operator|(intvec x) const
    {
      return _mm256_castpd_si256(_mm256_or_pd(_mm256_castsi256_pd(v),
                                              _mm256_castsi256_pd(x.v)));
    }
    intvec operator^(intvec x) const
    {
      return _mm256_castpd_si256(_mm256_xor_pd(_mm256_castsi256_pd(v),
                                               _mm256_castsi256_pd(x.v)));
    }
    
    intvec& operator&=(intvec const& x) { return *this=*this&x; }
    intvec& operator|=(intvec const& x) { return *this=*this|x; }
    intvec& operator^=(intvec const& x) { return *this=*this^x; }
    
    
    
    // SSE2 shift instructions potentially relevant for 64-bit
    // operations:
    //    _mm_slli_epi64
    //    _mm_srli_epi64
    //    _mm_sll_epi64
    //    _mm_srl_epi64
    //    _mm_srai_epi32
    //    _mm_sra_epi32
    //    _mm_srli_si128
    //    _mm_slli_si128
    intvec lsr(int_t n) const
    {
      __m128i v01 = _mm256_castsi256_si128(v);
      __m128i v23 = _mm256_extractf128_si256(v, 1);
      v01 = _mm_srli_epi64(v01, n);
      v23 = _mm_srli_epi64(v23, n);
      return _mm256_insertf128_si256(_mm256_castsi128_si256(v01), v23, 1);
    }
    intvec operator>>(int_t n) const
    {
      __m128i v01 = _mm256_castsi256_si128(v);
      __m128i v23 = _mm256_extractf128_si256(v, 1);
      // There is no _mm_srai_epi64. To emulate it, add 0x8000000
      // before shifting, and subtracted the shifted 0x80000000 after shifting
#if 0
      __m128i signmask01 = _mm_sub_epi64(_mm_set1_epi64x(0),
                                         _mm_srli_epi64(v01, 63));
      __m128i signmask23 = _mm_sub_epi64(_mm_set1_epi64x(0),
                                         _mm_srli_epi64(v23, 63));
      v01 = _mm_xor_si128(signmask01, v01);
      v23 = _mm_xor_si128(signmask23, v23);
      v01 = _mm_srli_epi64(v01, n);
      v23 = _mm_srli_epi64(v23, n);
      v01 = _mm_xor_si128(signmask01, v01);
      v23 = _mm_xor_si128(signmask23, v23);
#else
      // Convert signed to unsiged
      v01 += _mm_set1_epi64x(U(1) << (bits-1));
      v23 += _mm_set1_epi64x(U(1) << (bits-1));
      // Shift
      v01 = _mm_srli_epi64(v01, n);
      v23 = _mm_srli_epi64(v23, n);
      // Undo conversion
      v01 -= _mm_set1_epi64x(U(1) << (bits-n));
      v23 -= _mm_set1_epi64x(U(1) << (bits-n));
#endif
      return _mm256_insertf128_si256(_mm256_castsi128_si256(v01), v23, 1);
    }
    intvec operator<<(int_t n) const
    {
      __m128i v01 = _mm256_castsi256_si128(v);
      __m128i v23 = _mm256_extractf128_si256(v, 1);
      v01 = _mm_slli_epi64(v01, n);
      v23 = _mm_slli_epi64(v23, n);
      return _mm256_insertf128_si256(_mm256_castsi128_si256(v01), v23, 1);
    }
    intvec& operator>>=(int_t n) { return *this=*this>>n; }
    intvec& operator<<=(int_t n) { return *this=*this<<n; }
    
    intvec lsr(intvec n) const
    {
      __m128i v01 = _mm256_castsi256_si128(v);
      __m128i v23 = _mm256_extractf128_si256(v, 1);
      __m128i nv01 = _mm256_castsi256_si128(n.v);
      __m128i nv23 = _mm256_extractf128_si256(n.v, 1);
      v01 = _mm_srl_epi64(v01, nv01);
      v23 = _mm_srl_epi64(v23, nv23);
      return _mm256_insertf128_si256(_mm256_castsi128_si256(v01), v23, 1);
    }
    intvec operator>>(intvec n) const
    {
      __m128i v01 = _mm256_castsi256_si128(v);
      __m128i v23 = _mm256_extractf128_si256(v, 1);
      __m128i nv01 = _mm256_castsi256_si128(n.v);
      __m128i nv23 = _mm256_extractf128_si256(n.v, 1);
#if 0
      // There is no _mm_srai_epi64. To emulate it, invert all bits
      // before and after shifting if the sign bit is set.
      __m128i signmask01 = _mm_sub_epi64(_mm_set1_epi64x(0),
                                         _mm_srli_epi64(v01, 63));
      __m128i signmask23 = _mm_sub_epi64(_mm_set1_epi64x(0),
                                         _mm_srli_epi64(v23, 63));
      v01 = _mm_xor_si128(signmask01, v01);
      v23 = _mm_xor_si128(signmask23, v23);
      v01 = _mm_srl_epi64(v01, nv01);
      v23 = _mm_srl_epi64(v23, nv23);
      v01 = _mm_xor_si128(signmask01, v01);
      v23 = _mm_xor_si128(signmask23, v23);
#else
      // Convert signed to unsiged
      v01 += _mm_set1_epi64x(U(1) << (bits-1));
      v23 += _mm_set1_epi64x(U(1) << (bits-1));
      // Shift
      v01 = _mm_srl_epi64(v01, nv01);
      v23 = _mm_srl_epi64(v23, nv23);
      // Undo conversion
      v01 -= _mm_sll_epi64(_mm_set1_epi64x(1),
                           _mm_sub_epi64(_mm_set1_epi64x(bits), nv01));
      v23 -= _mm_sll_epi64(_mm_set1_epi64x(1),
                           _mm_sub_epi64(_mm_set1_epi64x(bits), nv23));
#endif
      return _mm256_insertf128_si256(_mm256_castsi128_si256(v01), v23, 1);
    }
    intvec operator<<(intvec n) const
    {
      __m128i v01 = _mm256_castsi256_si128(v);
      __m128i v23 = _mm256_extractf128_si256(v, 1);
      __m128i nv01 = _mm256_castsi256_si128(n.v);
      __m128i nv23 = _mm256_extractf128_si256(n.v, 1);
      v01 = _mm_sll_epi64(v01, nv01);
      v23 = _mm_sll_epi64(v23, nv23);
      return _mm256_insertf128_si256(_mm256_castsi128_si256(v01), v23, 1);
    }
    intvec& operator>>=(intvec n) { return *this=*this>>n; }
    intvec& operator<<=(intvec n) { return *this=*this<<n; }
  };
  
  
  
  template<>
  struct realvec<double,4>: floatprops<double>
  {
    static const int size = 4;
    typedef real_t scalar_t;
    typedef __m256d vector_t;
    
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
    realvec(real_t a): v(_mm256_set1_pd(a)) {}
    realvec(real_t const* as): v(_mm256_set_pd(as[3], as[2], as[1], as[0])) {}
    
    operator vector_t() const { return v; }
    real_t operator[](int n) const { return ((real_t const*)&v)[n]; }
    realvec& set_elt(int n, real_t a) { return ((real_t*)&v)[n]=a, *this; }
    
    
    
    intvec_t as_int() const { return _mm256_castpd_si256(v); }
    intvec_t convert_int() const
    {
#if 0
      __m128i iv0123 = _mm256_cvtpd_epi32(v);
      __m128i iv2301 = _mm_shuffle_ps(iv0123, iv0123, 0b10110001);
      __m256i iv01232301 =
        _mm256_insertf128_si256(_mm256_castsi128_si256(iv0123), iv2301, 1);
      __m256i zero = _mm256_setzero_ps();
      return _mm256_unpacklo_ps(iv01232301, zero);
#else
      realvec x = _mm256_floor_pd(v);
      intvec_t ix = x.as_int();
      boolvec_t sign = ix.as_bool();
      intvec_t exponent = (ix & exponent_mask) >> mantissa_bits;
      ix &= mantissa_mask;
      ix |= U(1) << mantissa_bits;  // add hidden bit
      ix <<= exponent - exponent_offset + 52; // ???
      ix = ifthen(sign, -ix, ix);
      return ix;
#endif
    }
    
    
    
    realvec operator+() const { return *this; }
    realvec operator-() const { return _mm256_sub_pd(_mm256_set1_pd(0.0), v); }
    
    realvec operator+(realvec x) const { return _mm256_add_pd(v, x.v); }
    realvec operator-(realvec x) const { return _mm256_sub_pd(v, x.v); }
    realvec operator*(realvec x) const { return _mm256_mul_pd(v, x.v); }
    realvec operator/(realvec x) const { return _mm256_div_pd(v, x.v); }
    
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
      return _mm256_cmp_pd(v, x.v, _CMP_EQ_OQ);
    }
    boolvec_t operator!=(realvec const& x) const
    {
      return _mm256_cmp_pd(v, x.v, _CMP_NEQ_OQ);
    }
    boolvec_t operator<(realvec const& x) const
    {
      return _mm256_cmp_pd(v, x.v, _CMP_LT_OQ);
    }
    boolvec_t operator<=(realvec const& x) const
    {
      return _mm256_cmp_pd(v, x.v, _CMP_LE_OQ);
    }
    boolvec_t operator>(realvec const& x) const
    {
      return _mm256_cmp_pd(v, x.v, _CMP_GT_OQ);
    }
    boolvec_t operator>=(realvec const& x) const
    {
      return _mm256_cmp_pd(v, x.v, _CMP_GE_OQ);
    }
    
    
    
#if 0
    realvec copysign(realvec y) const
    {
      uint_t signmask = U(1) << (bits-1);
      intvec_t value = as_int() & IV(~signmask);
      intvec_t sign = y.as_int() & IV(signmask);
      return (sign | value).as_float();
    }
    
    realvec fabs() const
    {
      uint_t signmask = U(1) << (bits-1);
      return (as_int() & IV(~signmask)).as_float();
    }
    
    intvec_t ilogb() const
    {
      intvec_t exponent_mask =
        ((U(1) << exponent_bits) - U(1)) << mantissa_bits;
      return lsr(as_int() & exponent_mask, mantissa_bits) - IV(exponent_offset);
    }
    
    realvec scalbn(intvec_t n) const
    {
      return *this * ((n + exponent_offset) << mantissa_bits).as_float();
    }
    
    boolvec_t signbit() const
    {
      return v;
    }
#endif
    
    realvec copysign(realvec y) const { return MF::vml_copysign(*this, y); }
    realvec fabs() const { return MF::vml_fabs(*this); }
    intvec_t ilogb() const { return MF::vml_ilogb(*this); }
    realvec scalbn(intvec_t n) const { return MF::vml_scalbn(*this, n); }
    boolvec_t signbit() const { return v; }
    
    
    
    realvec acos() const { return MF::vml_acos(*this); }
    realvec asin() const { return MF::vml_asin(*this); }
    realvec atan() const { return MF::vml_atan(*this); }
    realvec floor() const { return _mm256_floor_pd(v); }
    realvec log() const { return MF::vml_log(*this); }
    realvec log10() const { return MF::vml_log10(*this); }
    realvec log1p() const { return MF::vml_log1p(*this); }
    realvec log2() const { return MF::vml_log2(*this); }
    realvec rcp() const { return _mm256_div_pd(_mm256_set1_pd(1.0), v); }
    // realvec rcp() const { return MF::vml_rcp(*this); }
    realvec rsqrt() const { return MF::vml_rsqrt(*this); }
    realvec sqrt() const { return _mm256_sqrt_pd(v); }
    // realvec sqrt() const { return MF::vml_sqrt(*this); }
  };
  
  
  
  // boolvec definitions
  
  inline
  auto boolvec<double,4>::as_int() const -> intvec_t
  {
    return _mm256_castpd_si256(v);
  }
  
  inline
  auto boolvec<double,4>::convert_int() const -> intvec_t
  {
    //return ifthen(v, U(1), U(0));
    return lsr(as_int(), bits-1);
  }
  
  inline
  auto boolvec<double,4>::ifthen(intvec_t x, intvec_t y) const -> intvec_t
  {
    return ifthen(x.as_float(), y.as_float()).as_int();
  }
  
  inline
  auto boolvec<double,4>::ifthen(realvec_t x, realvec_t y) const -> realvec_t
  {
    return _mm256_blendv_pd(y.v, x.v, v);
  }

  
  
  // intvec definitions
  
  inline auto intvec<double,4>::as_float() const -> realvec_t
  {
    return _mm256_castsi256_pd(v);
  }
  
  inline auto intvec<double,4>::convert_float() const -> realvec_t
  {
    intvec x = v;
    uint_t signbitmask = U(1) << (bits-1);
    // make unsigned
    x += signbitmask;
    // convert lower 52 bits
    intvec xlo = x & IV((U(1) << mantissa_bits) - 1);
    int_t exponent_0 = exponent_offset + mantissa_bits;
    xlo = xlo | exponent_0;
    realvec_t flo = xlo.as_float() - FP::as_float(exponent_0);
    // convert upper 22 bits
    intvec xhi = x.lsr(U(mantissa_bits));
    int_t exponent_52 = exponent_0 + mantissa_bits;
    xhi = xhi | exponent_52;
    realvec_t fhi = xhi.as_float() - FP::as_float(exponent_52);
    return flo + fhi - R(signbitmask);
  }
  
} // namespace vecmathlib

#endif  // #ifndef VEC_DOUBLE_AVX_H
