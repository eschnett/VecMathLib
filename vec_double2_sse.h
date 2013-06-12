// -*-C++-*-

#ifndef VEC_DOUBLE_SSE2_H
#define VEC_DOUBLE_SSE2_H

#include "mathfuncs.h"
#include "real_properties.h"
#include "vec_base.h"

#include <cmath>

// SSE2 intrinsics
#include <emmintrin.h>
#ifdef __SSE3__                 // Intel's SSE 3
#  include <pmmintrin.h>
#endif
#ifdef __SSE4_1__               // Intel's SSE 4.1
#  include <smmintrin.h>
#endif
#ifdef __SSE4A__                // AMD's SSE 4a
#  include <ammintrin.h>
#endif
#if defined __AVX__             // Intel's AVX
#  include <immintrin.h>
#endif



namespace vecmathlib {
  
  
  
  struct vec_double2_sse;
  
  template<>
  struct vec_types<double,2> {
    typedef vec_double2_sse impl;
  };
#define VECMATHLIB_HAVE_DOUBLE2
  
  
  
  struct vec_double2_sse {
    // typedef __m128d boolvec_t;
    struct boolvec_t { __m128d v; };
    typedef __m128i  intvec_t;
    typedef __m128d  realvec_t;
    typedef bool     bool_t;
    typedef double   bool_impl_t;
    typedef int64_t  int_t;
    typedef uint64_t uint_t;
    typedef double   real_t;
    
    static const int size = sizeof(realvec_t) / sizeof(real_t);
    static_assert(size * sizeof(bool_impl_t) == sizeof(boolvec_t),
                  "sizeof(boolvec_t) is inconsistent");
    static_assert(size * sizeof(int_t) == sizeof(intvec_t),
                  "sizeof(intvec_t) is inconsistent");
    static_assert(size * sizeof(real_t) == sizeof(realvec_t),
                  "sizeof(realvec_t) is inconsistent");
    
    typedef real_properties<real_t> RP;
    
    // true values have the sign bit set, false values have it unset
    static bool_impl_t from_bool(bool_t a)
    {
      return RP::as_real(-int_t(a));
    }
    static bool_t to_bool(bool_impl_t a)
    {
      return RP::as_int(a) < int_t(0);
    }
  };
#define VT vec_double2_sse
  
  
  
  VT::boolvec_t& set_elt(VT::boolvec_t& x, int n, VT::bool_t b);
  VT::bool_t elt(const VT::boolvec_t& x, int n);
  VT::boolvec_t as_bool(VT::boolvec_t x);
  VT::intvec_t as_int(VT::boolvec_t x);
  VT::realvec_t as_real(VT::boolvec_t x);
  VT::boolvec_t operator!=(VT::boolvec_t x, VT::boolvec_t y);
  VT::boolvec_t ifthen(VT::boolvec_t x, VT::boolvec_t y, VT::boolvec_t z);
  VT::intvec_t ifthen(VT::boolvec_t x, VT::intvec_t y, VT::intvec_t z);
  VT::realvec_t ifthen(VT::boolvec_t x, VT::realvec_t y, VT::realvec_t z);
  
  VT::intvec_t& set_elt(VT::intvec_t& x, int n, VT::int_t i);
  VT::int_t elt(const VT::intvec_t& x, int n);
  VT::boolvec_t as_bool(VT::intvec_t x);
  VT::intvec_t as_int(VT::intvec_t x);
  VT::realvec_t as_real(VT::intvec_t x);
  VT::intvec_t lsr(VT::intvec_t x, VT::int_t y);
  VT::intvec_t lsr(VT::intvec_t x, VT::intvec_t y);
    
  VT::realvec_t& set_elt(VT::realvec_t& x, int n, VT::real_t r);
  VT::real_t elt(const VT::realvec_t& x, int n);
  VT::boolvec_t as_bool(VT::realvec_t x);
  VT::intvec_t as_int(VT::realvec_t x);
  VT::realvec_t as_real(VT::realvec_t x);
  
  
  
  // bool operations
  
  template<>
  VT::boolvec_t make_boolvec<VT>(VT::bool_t b)
  {
    return as_bool(make_realvec<VT>(VT::from_bool(b)));
  }
  
  template<>
  VT::boolvec_t make_boolvec<VT>(const VT::bool_t* restrict bs)
  {
    VT::real_t rs[] = { VT::from_bool(bs[0]), VT::from_bool(bs[1]) };
    return as_bool(make_realvec<VT>(rs));
  }
  
  VT::boolvec_t& set_elt(VT::boolvec_t& x, int n, VT::bool_t b)
  {
    return set_elt(x.v, n, VT::from_bool(b)), x;
  }
  
  VT::bool_t elt(const VT::boolvec_t& x, int n)
  {
    return VT::to_bool(elt(as_real(x), n));
  }
  
  VT::boolvec_t as_bool(VT::boolvec_t x)
  {
    return x;
  }
  
  VT::intvec_t as_int(VT::boolvec_t x)
  {
    return as_int(as_real(x));
  }
  
  VT::realvec_t as_real(VT::boolvec_t x)
  {
    return x.v;
  }
  
  VT::boolvec_t convert_bool(VT::boolvec_t x)
  {
    return x;
  }
  
  VT::intvec_t convert_int(VT::boolvec_t x)
  {
    return lsr(as_int(x), VT::RP::bits-1);
  }
  
  VT::boolvec_t operator!(VT::boolvec_t x)
  {
    return as_bool(_mm_xor_pd(as_real(make_boolvec<VT>(true)), as_real(x)));
  }
  
  VT::boolvec_t operator&&(VT::boolvec_t x, VT::boolvec_t y)
  {
    return as_bool(_mm_and_pd(as_real(x), as_real(y)));
  }
  
  VT::boolvec_t operator||(VT::boolvec_t x, VT::boolvec_t y)
  {
    return as_bool(_mm_or_pd(as_real(x), as_real(y)));
  }
  
  VT::boolvec_t operator==(VT::boolvec_t x, VT::boolvec_t y)
  {
    return ! (x!=y);
  }
  
  VT::boolvec_t operator!=(VT::boolvec_t x, VT::boolvec_t y)
  {
    return as_bool(_mm_xor_pd(as_real(x), as_real(y)));
  }
  
  VT::bool_t all(VT::boolvec_t x)
  {
    x = x && as_bool(_mm_shuffle_pd(x.v, x.v, _MM_SHUFFLE2(0,1)));
    return elt(x, 0);
  }
  
  VT::bool_t any(VT::boolvec_t x)
  {
    x = x || as_bool(_mm_shuffle_pd(x.v, x.v, _MM_SHUFFLE2(0,1)));
    return elt(x, 0);
  }
  
  VT::boolvec_t ifthen(VT::boolvec_t x, VT::boolvec_t y, VT::boolvec_t z)
  {
    return as_bool(ifthen(x, as_int(y), as_int(z)));
  }
  
  VT::intvec_t ifthen(VT::boolvec_t x, VT::intvec_t y, VT::intvec_t z)
  {
    return as_int(ifthen(x, as_real(y), as_real(z)));
  }
  
  VT::realvec_t ifthen(VT::boolvec_t x, VT::realvec_t y, VT::realvec_t z)
  {
#ifdef __SSE4_1__
    return _mm_blendv_pd(z, y, x.v);
#else
    return as_real(( -convert_int(x) & as_int(y)) |
                   (~-convert_int(x) & as_int(z)));
#endif
  }
  
  
  
  // int operations
  
  template<>
  VT::intvec_t make_intvec<VT>(VT::int_t i)
  {
    return _mm_set1_epi64x(i);
  }
  
  template<>
  VT::intvec_t make_intvec<VT>(const VT::int_t* restrict is)
  {
    return _mm_set_epi64x(is[1], is[0]);
  }
  
  template<>
  VT::intvec_t make_iota<VT>()
  {
    const VT::int_t is[] = { 0, 1 };
    return make_intvec<VT>(is);
  }
  
  VT::intvec_t& set_elt(VT::intvec_t& x, int n, VT::int_t i)
  {
    return ((VT::int_t*)&x)[n] = i, x;
  }
  
  VT::int_t elt(const VT::intvec_t& x, int n)
  {
    return ((const VT::int_t*)&x)[n];
  }
  
  VT::boolvec_t as_bool(VT::intvec_t x)
  {
    return as_bool(as_real(x));
  }
  
  VT::intvec_t as_int(VT::intvec_t x)
  {
    return x;
  }
  
  VT::realvec_t as_real(VT::intvec_t x)
  {
    return _mm_castsi128_pd(x);
  }
  
  VT::intvec_t convert_int(VT::intvec_t x)
  {
    return x;
  }
  
  VT::realvec_t convert_real(VT::intvec_t x)
  {
    return vml_convert_real(x);
  }
  
  VT::intvec_t operator+(VT::intvec_t x)
  {
    return x;
  }
  
  VT::intvec_t operator-(VT::intvec_t x)
  {
    return make_intvec<VT>(0) - x;
  }
  
  VT::intvec_t operator+(VT::intvec_t x, VT::intvec_t y)
  {
    return _mm_add_epi64(x, y);
  }
  
  VT::intvec_t operator-(VT::intvec_t x, VT::intvec_t y)
  {
    return _mm_sub_epi64(x, y);
  }
  
  VT::intvec_t operator*(VT::intvec_t x, VT::intvec_t y)
  {
    VT::intvec_t r;
    for (int n=0; n<VT::size; ++n) {
      set_elt(r, n, x[n] * y[n]);
    }
    return r;
  }
  
  VT::intvec_t operator/(VT::intvec_t x, VT::intvec_t y)
  {
    VT::intvec_t r;
    for (int n=0; n<VT::size; ++n) {
      set_elt(r, n, x[n] / y[n]);
    }
    return r;
  }
  
  VT::intvec_t operator%(VT::intvec_t x, VT::intvec_t y)
  {
    VT::intvec_t r;
    for (int n=0; n<VT::size; ++n) {
      set_elt(r, n, x[n] % y[n]);
    }
    return r;
  }
  
  VT::intvec_t& operator+=(VT::intvec_t& x, VT::intvec_t y)
  {
    return x = x + y;
  }
  
  VT::intvec_t& operator-=(VT::intvec_t& x, VT::intvec_t y)
  {
    return x = x - y;
  }
  
  VT::intvec_t& operator*=(VT::intvec_t& x, VT::intvec_t y)
  {
    return x = x * y;
  }
  
  VT::intvec_t& operator/=(VT::intvec_t& x, VT::intvec_t y)
  {
    return x = x / y;
  }
  
  VT::intvec_t& operator%=(VT::intvec_t& x, VT::intvec_t y)
  {
    return x = x % y;
  }
  
  VT::intvec_t operator~(VT::intvec_t x)
  {
    return make_intvec<VT>(~0) ^ x;
  }
  
  VT::intvec_t operator&(VT::intvec_t x, VT::intvec_t y)
  {
    return as_int(_mm_and_pd(as_real(x), as_real(y));
  }
  
  VT::intvec_t operator|(VT::intvec_t x, VT::intvec_t y)
  {
    return as_int(_mm_or_pd(as_real(x), as_real(y));
  }
  
  VT::intvec_t operator^(VT::intvec_t x, VT::intvec_t y)
  {
    return as_int(_mm_xor_pd(as_real(x), as_real(y));
  }
  
  VT::intvec_t& operator&=(VT::intvec_t& x, VT::intvec_t y)
  {
    return x = x & y;
  }
  
  VT::intvec_t& operator|=(VT::intvec_t& x, VT::intvec_t y)
  {
    return x = x | y;
  }
  
  VT::intvec_t& operator^=(VT::intvec_t& x, VT::intvec_t y)
  {
    return x = x ^ y;
  }
  
  VT::intvec_t lsr(VT::intvec_t x, VT::int_t y)
  {
    return _mm_srli_epi64(x, y);
  }
  
  VT::intvec_t operator>>(VT::intvec_t x, VT::int_t y)
  {
    // There is no _mm_srai_epi64. To emulate it, add 0x80000000
    // before shifting, and subtract the shifted 0x80000000 after
    // shifting.
    // Convert signed to unsiged
    x += make_intvec<VT>(int_t(1) << (RP::bits-1));
    // Shift
    x = x.lsr(y);
    // Undo conversion
    x -= make_intvec<VT>(int_t(1) << (RP::bits-1-y));
    return x;
  }
  
  VT::intvec_t operator<<(VT::intvec_t x, VT::int_t y)
  {
    return _mm_srli_epi64(x, y);
  }
  
  VT::intvec_t& operator>>=(VT::intvec_t& x, VT::int_t y)
  {
    return x = x >> y;
  }
  
  VT::intvec_t& operator<<=(VT::intvec_t& x, VT::int_t y)
  {
    return x = x << y;
  }
  
  VT::intvec_t lsr(VT::intvec_t x, VT::intvec_t y)
  {
    VT::intvec_t r;
    for (int n=0; n<VT::size; ++n) {
      set_elt(r, n, uint_t(x[n]) << uint_t(y[n]));
    }
    return r;
  }
  
  VT::intvec_t operator>>(VT::intvec_t x, VT::intvec_t y)
  {
    VT::intvec_t r;
    for (int n=0; n<VT::size; ++n) {
      set_elt(r, n, x[n] >> y[n]);
    }
    return r;
  }
  
  VT::intvec_t operator<<(VT::intvec_t x, VT::intvec_t y)
  {
    VT::intvec_t r;
    for (int n=0; n<VT::size; ++n) {
      set_elt(r, n, x[n] << y[n]);
    }
    return r;
  }
  
  VT::intvec_t& operator>>=(VT::intvec_t& x, VT::intvec_t y)
  {
    return x = x >> y;
  }
  
  VT::intvec_t& operator<<=(VT::intvec_t& x, VT::intvec_t y)
  {
    return x = x << y;
  }
  
  VT::boolvec_t operator==(VT::intvec_t x, VT::intvec_t y)
  {
    return ! (x!=y);
  }
  
  VT::boolvec_t operator!=(VT::intvec_t x, VT::intvec_t y)
  {
    return convert_bool(x ^ y);
  }
  
  VT::boolvec_t operator<(VT::intvec_t x, VT::intvec_t y)
  {
    VT::intvec_t r;
    for (int n=0; n<VT::size; ++n) {
      set_elt(r, n, x[n] < y[n]);
    }
    return r;
  }
  
  VT::boolvec_t operator<=(VT::intvec_t x, VT::intvec_t y)
  {
    return ! (x>y);
  }
  
  VT::boolvec_t operator>(VT::intvec_t x, VT::intvec_t y)
  {
    return y<x;
  }
  
  VT::boolvec_t operator>=(VT::intvec_t x, VT::intvec_t y)
  {
    return ! (x<y);
  }
  
  
  
  // real operations
  
  template<>
    VT::realvec_t make_realvec<VT>(real_t r)
  {
    return _mm_set1_pd(r);
  }
  
  template<>
  VT::realvec_t make_realvec<VT>(const real_t* restrict rs)
  {
    return _mm_set_pd(rs[1], rs[0]);
  }
  
  VT::realvec_t& set_elt(VT::intvec_t& x, int n, real_t r)
  {
    return ((real_t*)&x)[n] = r, x;
  }
  
  VT::real_t elt(const VT::realvec_t& x, int n)
  {
    return ((const real_t*)&x)[n];
  }
  
  VT::boolvec_t as_bool(VT::realvec_t x)
  {
    return VT::boolvec_t({x});
  }
  
  VT::intvec_t as_int(VT::realvec_t x)
  {
    return _mm_castpd_si128(x);
  }
  
  VT::realvec_t as_real(VT::realvec_t x)
  {
    return x;
  }
  
  VT::intvec_t convert_int(VT::intvec_t x)
  {
    return vml_convert_int(x);
  }
  
  VT::realvec_t convert_real(VT::intvec_t x)
  {
    return x;
  }
  
  VT::realvec_t operator+(VT::realvec_t x)
  {
    return x;
  }
  
  VT::realvec_t operator-(VT::realvec_t x)
  {
    return make_realvec<VT>(0.0) - x;
  }
  
  VT::realvec_t operator+(VT::realvec_t x, VT::realvec_t y)
  {
    return _mm_add_pd(x, y);
  }
  
  VT::realvec_t operator-(VT::realvec_t x, VT::realvec_t y)
  {
    return _mm_sub_pd(x, y);
  }
  
  VT::realvec_t operator*(VT::realvec_t x, VT::realvec_t y)
  {
    return _mm_mul_pd(x, y);
  }
  
  VT::realvec_t operator/(VT::realvec_t x, VT::realvec_t y)
  {
    return _mm_div_pd(x, y);
  }
  
  VT::realvec_t& operator+=(VT::realvec_t& x, VT::realvec_t y)
  {
    return x = x + y;
  }
  
  VT::realvec_t& operator-=(VT::realvec_t& x, VT::realvec_t y)
  {
    return x = x - y;
  }
  
  VT::realvec_t& operator*=(VT::realvec_t& x, VT::realvec_t y)
  {
    return x = x * y;
  }
  
  VT::realvec_t& operator/=(VT::realvec_t& x, VT::realvec_t y)
  {
    return x = x / y;
  }
  
  VT::real_t prod(VT::realvec_t x)
  {
    return x[0] * x[1];
  }
  
  VT::real_t sum(VT::realvec_t x)
  {
#ifdef __SSE3__
    x = _mm_hadd_pd(x, x);
    return x[0];
#else
    return x[0] + x[1];
#endif
  }
  
  VT::boolvec_t operator==(VT::realvec_t x, VT::realvec_t y)
  {
    return as_bool(_mm_cmpeq_pd(x, y));
  }
  
  VT::boolvec_t operator!=(VT::realvec_t x, VT::realvec_t y)
  {
    return as_bool(_mm_cmpneq_pd(x, y));
  }
  
  VT::boolvec_t operator<(VT::realvec_t x, VT::realvec_t y)
  {
    return as_bool(_mm_cmplt_pd(x, y));
  }
  
  VT::boolvec_t operator<=(VT::realvec_t x, VT::realvec_t y)
  {
    return as_bool(_mm_cmple_pd(x, y));
  }
  
  VT::boolvec_t operator>(VT::realvec_t x, VT::realvec_t y)
  {
    return as_bool(_mm_cmpgt_pd(x, y));
  }
  
  VT::boolvec_t operator>=(VT::realvec_t x, VT::realvec_t y)
  {
    return as_bool(_mm_cmpge_pd(x, y));
  }
  
  
  
//TODO     realvec acos() const { return MF::vml_acos(*this); }
//TODO     realvec acosh() const { return MF::vml_acosh(*this); }
//TODO     realvec asin() const { return MF::vml_asin(*this); }
//TODO     realvec asinh() const { return MF::vml_asinh(*this); }
//TODO     realvec atan() const { return MF::vml_atan(*this); }
//TODO     realvec atan2(realvec y) const { return MF::vml_atan2(*this, y); }
//TODO     realvec atanh() const { return MF::vml_atanh(*this); }
//TODO     realvec cbrt() const { return MF::vml_cbrt(*this); }
//TODO     realvec ceil() const
//TODO     {
//TODO #ifdef __SSE4_1__
//TODO       return _mm_ceil_pd(v);
//TODO #else
//TODO       return MF::vml_ceil(*this);
//TODO #endif
//TODO  }
//TODO     realvec copysign(realvec y) const { return MF::vml_copysign(*this, y); }
//TODO     realvec cos() const { return MF::vml_cos(*this); }
//TODO     realvec cosh() const { return MF::vml_cosh(*this); }
//TODO     realvec exp() const { return MF::vml_exp(*this); }
//TODO     realvec exp10() const { return MF::vml_exp10(*this); }
//TODO     realvec exp2() const { return MF::vml_exp2(*this); }
//TODO     realvec expm1() const { return MF::vml_expm1(*this); }
//TODO     realvec fabs() const { return MF::vml_fabs(*this); }
//TODO     realvec fdim(realvec y) const { return MF::vml_fdim(*this, y); }
//TODO     realvec floor() const
//TODO     {
//TODO #ifdef __SSE4_1__
//TODO       return _mm_floor_pd(v);
//TODO #else
//TODO       return MF::vml_floor(*this);
//TODO #endif
//TODO  }
//TODO     realvec fma(realvec y, realvec z) const { return MF::vml_fma(*this, y, z); }
//TODO     realvec fmax(realvec y) const { return _mm_max_pd(v, y.v); }
//TODO     realvec fmin(realvec y) const { return _mm_min_pd(v, y.v); }
//TODO     realvec fmod(realvec y) const { return MF::vml_fmod(*this, y); }
//TODO     realvec hypot(realvec y) const { return MF::vml_hypot(*this, y); }
//TODO     intvec_t ilogb() const { return MF::vml_ilogb(*this); }
//TODO     boolvec_t isfinite() const { return MF::vml_isfinite(*this); }
//TODO     boolvec_t isinf() const { return MF::vml_isinf(*this); }
//TODO     boolvec_t isnan() const { return _mm_cmpunord_pd(v, v);; }
//TODO     boolvec_t isnormal() const { return MF::vml_isnormal(*this); }
//TODO     realvec ldexp(int_t n) const { return MF::vml_ldexp(*this, n); }
//TODO     realvec ldexp(intvec_t n) const { return MF::vml_ldexp(*this, n); }
//TODO     realvec log() const { return MF::vml_log(*this); }
//TODO     realvec log10() const { return MF::vml_log10(*this); }
//TODO     realvec log1p() const { return MF::vml_log1p(*this); }
//TODO     realvec log2() const { return MF::vml_log2(*this); }
//TODO     realvec nextafter(realvec y) const { return MF::vml_nextafter(*this, y); }
//TODO     realvec pow(realvec y) const { return MF::vml_pow(*this, y); }
//TODO     realvec rcp() const { return _mm_div_pd(_mm_set1_pd(1.0), v); }
//TODO     realvec remainder(realvec y) const { return MF::vml_remainder(*this, y); }
//TODO     realvec rint() const
//TODO     {
//TODO #ifdef __SSE4_1__
//TODO       return _mm_round_pd(v, _MM_FROUND_TO_NEAREST_INT);
//TODO #else
//TODO       return MF::vml_rint(*this);
//TODO #endif
//TODO     }
//TODO     realvec round() const { return MF::vml_round(*this); }
//TODO     realvec rsqrt() const { return MF::vml_rsqrt(*this); }
//TODO     boolvec_t signbit() const { return v; }
//TODO     realvec sin() const { return MF::vml_sin(*this); }
//TODO     realvec sinh() const { return MF::vml_sinh(*this); }
//TODO     realvec sqrt() const { return _mm_sqrt_pd(v); }
//TODO     realvec tan() const { return MF::vml_tan(*this); }
//TODO     realvec tanh() const { return MF::vml_tanh(*this); }
//TODO     realvec trunc() const
//TODO     {
//TODO #ifdef __SSE4_1__
//TODO       return _mm_round_pd(v, _MM_FROUND_TO_ZERO);
//TODO #else
//TODO       return MF::vml_trunc(*this);
//TODO #endif
//TODO  }
  
  
  
#undef VT
  
  
  
} // namespace vecmathlib

#endif  // #ifndef VEC_DOUBLE_SSE2_H
