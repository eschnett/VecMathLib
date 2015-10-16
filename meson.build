# See file "BUILD" for instructions

project('vecmathlib', 'cpp')

selftest_exe = executable('selftest', 'selftest.cc')
bench_exe = executable('bench', 'bench.cc')

example_exe = executable('example', 'example.cc')
example_float_exe = executable('example_float', 'example_float.cc')
loop_exe = executable('loop', 'loop.cc')
interp_exe = executable('interp', 'interp.cc')

instantiations_lib = static_library('instantiations', 'instantiations.cc')

test('selftest', selftest_exe)
