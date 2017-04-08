# See file "BUILD" for instructions

cmake_minimum_required (VERSION 2.6)

project (vecmathlib)

add_executable (example example.cc)
add_executable (example_float example_float.cc)
add_executable (loop loop.cc)
add_executable (interp interp.cc)
add_executable (bench bench.cc)

add_library (instantiations instantiations.cc)

# Tests

add_executable (selftest selftest.cc)
enable_testing()
add_test(NAME selftest COMMAND ./selftest)

# Install

set(HEADERS
  floatbuiltins.h
  floatprops.h
  floattypes.h
  mathfuncs.h
  mathfuncs_asin.h
  mathfuncs_asinh.h
  mathfuncs_base.h
  mathfuncs_convert.h
  mathfuncs_exp.h
  mathfuncs_fabs.h
  mathfuncs_int.h
  mathfuncs_log.h
  mathfuncs_pow.h
  mathfuncs_rcp.h
  mathfuncs_sin.h
  mathfuncs_sinh.h
  mathfuncs_sqrt.h
  vec_altivec_float4.h
  vec_avx_double4.h
  vec_avx_float8.h
  vec_avx_fp16_16.h
  vec_avx_fp8_32.h
  vec_base.h
  vec_builtin.h
  vec_mask.h
  vec_mic_double8.h
  vec_neon_float2.h
  vec_neon_float4.h
  vec_pseudo.h
  vec_qpx_double4.h
  vec_sse_double1.h
  vec_sse_double2.h
  vec_sse_float1.h
  vec_sse_float4.h
  vec_test.h
  vec_vsx_double2.h
  vecmathlib.h
  )

install(FILES ${HEADERS} DESTINATION include/vecmathlib)

set(PKG_CONFIG_REQUIRES "")
set(PKG_CONFIG_INCLUDEDIR "\${prefix}/include/vecmathlib")
set(PKG_CONFIG_LIBDIR "\${prefix}/lib")
set(PKG_CONFIG_CFLAGS "-I\${includedir}")
set(PKG_CONFIG_LIBS "")

configure_file(
  "${PROJECT_SOURCE_DIR}/pkg-config.pc.in"
  "${PROJECT_BINARY_DIR}/${PROJECT_NAME}.pc"
)

install(FILES "${CMAKE_BINARY_DIR}/${PROJECT_NAME}.pc" DESTINATION lib/pkgconfig)
