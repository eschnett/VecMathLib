# See file "BUILD" for instructions

set (GENERATED_SRCS
acos.cc acosh.cc asin.cc asinh.cc atan.cc atanh.cc cbrt.cc ceil.cc copysign.cc cos.cc cosh.cc exp.cc exp2.cc exp10.cc expm1.cc fabs.cc fdim.cc floor.cc fma.cc fmax.cc fmin.cc fmod.cc hypot.cc ilogb_.cc ldexp_.cc log.cc log2.cc log10.cc log1p.cc pow.cc remainder.cc rint.cc round.cc rsqrt.cc sin.cc sinh.cc sqrt.cc tan.cc tanh.cc trunc.cc isfinite.cc isinf.cc isnan.cc isnormal.cc signbit.cc acospi.cl asinpi.cl atanpi.cl atan2.cl atan2pi.cl cospi.cl fmax.cl fmin.cl fract.cl frexp.cl ilogb.cl ldexp.cl logb.cl mad.cl maxmag.cl minmag.cl modf.cl nan.cl pown.cl powr.cl remquo.cl rootn.cl sincos.cl sinpi.cl tanpi.cl half_cos.cl half_divide.cl half_exp.cl half_exp2.cl half_exp10.cl half_log.cl half_log2.cl half_log10.cl half_powr.cl half_recip.cl half_rsqrt.cl half_sin.cl half_sqrt.cl half_tan.cl native_cos.cl native_divide.cl native_exp.cl native_exp2.cl native_exp10.cl native_log.cl native_log2.cl native_log10.cl native_powr.cl native_recip.cl native_rsqrt.cl native_sin.cl native_sqrt.cl native_tan.cl clamp.cl degrees.cl max.cl min.cl mix.cl radians.cl step.cl smoothstep.cl sign.cl isequal.cl isnotequal.cl isgreater.cl isgreaterequal.cl isless.cl islessequal.cl islessgreater.cl isordered.cl isunordered.cl
  )

set (SRCS
  cross.cl dot.cl distance.cl length.cl normalize.cl fast_distance.cl
  fast_length.cl fast_normalize.cl
  ${GENERATED_SRCS}
)

add_custom_command(
  OUTPUT ${GENERATED_SRCS} kernel-vecmathlib.h
  COMMAND ./generate-files.py
  DEPENDS generate-files.py)

add_library (kernel-vecmathlib ${SRCS})
