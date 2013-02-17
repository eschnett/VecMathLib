#! /usr/bin/env python

import re



# Types:
SI = "SI"                       # int/long
SJ = "SJ"                       # int (even for double)
VB = "VB"                       # boolN
VF = "VF"                       # floatN/doubleN
VI = "VI"                       # intN/longN
VJ = "VJ"                       # intN/longN (except int1 for double1)
VK = "VK"                       # intN (even for doubleN)

# Each function is described by a tuple with the following entries:
#    1. name
#    2. external argument types (see above)
#    3. external return type
#    4. vecmathlib argument types (see above)
#    5. vecmathlib return type
# This allows generating externally visible functions with different
# signatures, e.g. to support OpenCL.
funcs = [
    ("acos"     , [VF        ], VF, [VF        ], VF),
    ("acosh"    , [VF        ], VF, [VF        ], VF),
    ("asin"     , [VF        ], VF, [VF        ], VF),
    ("asinh"    , [VF        ], VF, [VF        ], VF),
    ("atan"     , [VF        ], VF, [VF        ], VF),
    ("atanh"    , [VF        ], VF, [VF        ], VF),
    ("cbrt"     , [VF        ], VF, [VF        ], VF),
    ("ceil"     , [VF        ], VF, [VF        ], VF),
    ("copysign" , [VF, VF    ], VF, [VF, VF    ], VF),
    ("cos"      , [VF        ], VF, [VF        ], VF),
    ("cosh"     , [VF        ], VF, [VF        ], VF),
    ("exp"      , [VF        ], VF, [VF        ], VF),
    ("exp2"     , [VF        ], VF, [VF        ], VF),
    ("exp10"    , [VF        ], VF, [VF        ], VF),
    ("expm1"    , [VF        ], VF, [VF        ], VF),
    ("fabs"     , [VF        ], VF, [VF        ], VF),
    ("fdim"     , [VF, VF    ], VF, [VF, VF    ], VF),
    ("floor"    , [VF        ], VF, [VF        ], VF),
    ("fma"      , [VF, VF, VF], VF, [VF, VF, VF], VF),
    ("fmax"     , [VF, VF    ], VF, [VF, VF    ], VF),
    ("fmin"     , [VF, VF    ], VF, [VF, VF    ], VF),
    ("fmod"     , [VF, VF    ], VF, [VF, VF    ], VF),
    ("hypot"    , [VF, VF    ], VF, [VF, VF    ], VF),
    ("ilogb"    , [VF        ], VJ, [VF        ], VI), # should return VK
    ("isfinite" , [VF        ], VJ, [VF        ], VB),
    ("isinf"    , [VF        ], VJ, [VF        ], VB),
    ("isnan"    , [VF        ], VJ, [VF        ], VB),
    ("isnormal" , [VF        ], VJ, [VF        ], VB),
    ("ldexp"    , [VF, VJ    ], VF, [VF, VI    ], VF), # should take VK
    ("ldexp"    , [VF, SJ    ], VF, [VF, SI    ], VF),
    ("log"      , [VF        ], VF, [VF        ], VF),
    ("log2"     , [VF        ], VF, [VF        ], VF),
    ("log10"    , [VF        ], VF, [VF        ], VF),
    ("log1p"    , [VF        ], VF, [VF        ], VF),
    ("pow"      , [VF, VF    ], VF, [VF, VF    ], VF),
    ("remainder", [VF, VF    ], VF, [VF, VF    ], VF),
    ("round"    , [VF        ], VF, [VF        ], VF),
    ("rsqrt"    , [VF        ], VF, [VF        ], VF),
    ("signbit"  , [VF        ], VJ, [VF        ], VB),
    ("sin"      , [VF        ], VF, [VF        ], VF),
    ("sinh"     , [VF        ], VF, [VF        ], VF),
    ("sqrt"     , [VF        ], VF, [VF        ], VF),
    ("tan"      , [VF        ], VF, [VF        ], VF),
    ("tanh"     , [VF        ], VF, [VF        ], VF),
    ("trunc"    , [VF        ], VF, [VF        ], VF),
    ]



def output_prelude():
    print """\
// Instantiante functions to create a library that can be called from
// elsewhere

// Note: This file has been generated automatically.
//       Do not modify directly!
//       Your changes would be lost.

// Note: We use a prefix __vml_ for all functions to avoid namespace
//       collisions with libm.



// Make things go fast (and debugging difficult...)
#define VML_NODEBUG
#include "vecmathlib.h"

#include <algorithm>
#include <cstdint>
#include <cstring>



#define cl_khr_fp64
#define cles_khr_int64



// Generic conversion function
template<typename A, typename B>
static B bitcast(A a)
{
  B b;
  std::memcpy(&b, &a, std::min(sizeof a, sizeof b));
  if (sizeof b > sizeof a) {
    std::memset((char*)&b + sizeof a, 0, sizeof b - sizeof a);
  }
  return b;
}



// Define vector types

using std::int32_t;
#define int int32_t
typedef int int2  __attribute__((__ext_vector_type__( 2)));
typedef int int3  __attribute__((__ext_vector_type__( 3)));
typedef int int4  __attribute__((__ext_vector_type__( 4)));
typedef int int8  __attribute__((__ext_vector_type__( 8)));
typedef int int16 __attribute__((__ext_vector_type__(16)));

#ifdef cles_khr_int64
using std::int64_t;
#define long int64_t
typedef long long2  __attribute__((__ext_vector_type__( 2)));
typedef long long3  __attribute__((__ext_vector_type__( 3)));
typedef long long4  __attribute__((__ext_vector_type__( 4)));
typedef long long8  __attribute__((__ext_vector_type__( 8)));
typedef long long16 __attribute__((__ext_vector_type__(16)));
#endif

typedef float float2  __attribute__((__ext_vector_type__( 2)));
typedef float float3  __attribute__((__ext_vector_type__( 3)));
typedef float float4  __attribute__((__ext_vector_type__( 4)));
typedef float float8  __attribute__((__ext_vector_type__( 8)));
typedef float float16 __attribute__((__ext_vector_type__(16)));

#ifdef cl_khr_fp64
typedef double double2  __attribute__((__ext_vector_type__( 2)));
typedef double double3  __attribute__((__ext_vector_type__( 3)));
typedef double double4  __attribute__((__ext_vector_type__( 4)));
typedef double double8  __attribute__((__ext_vector_type__( 8)));
typedef double double16 __attribute__((__ext_vector_type__(16)));
#endif



"""



def mktype(tp, vectype):
    (basetype, size) = re.match("([A-Za-z]+)([0-9]*)", vectype).groups()
    size = 1 if size=="" else int(size)
    if tp==SJ:
        if size==1: return "int"
        return "int" if basetype=="float" else "long"
    if tp==VI:
        ibasetype = "int" if basetype=="float" else "long"
        return "%s%s" % (ibasetype, "" if size==1 else str(size))
    if tp==VJ:
        if size==1: return "int"
        ibasetype = "int" if basetype=="float" else "long"
        return "%s%d" % (ibasetype, size)
    if tp==VK:
        if size==1: return "int"
        return "int%d" % size
    if tp==VF:
        return vectype
    raise "unreachable"

def mkvmltype(tp, vectype):
    if tp==SI: return vectype+"::int_t"
    if tp==VB: return vectype+"::boolvec_t"
    if tp==VI: return vectype+"::intvec_t"
    if tp==VF: return vectype
    raise "unreachable"



def output_func_vml(func, vectype):
    (name, args, ret, vmlargs, vmlret) = func
    print "// Implement %s by calling vecmathlib" % name
    (basetype, size) = re.match("([A-Za-z]+)([0-9]*)", vectype).groups()
    size = 1 if size=="" else int(size)
    vmltype = "vecmathlib::realvec<%s,%d>" % (basetype, size)
    vmlinttype = "%s::intvec_t" % vmltype
    vmlbooltype = "%s::boolvec_t" % vmltype
    funcargstr = ", ".join(map(lambda (n, arg):
                                   "%s x%d" % (mktype(arg, vectype), n),
                               zip(range(0, 100), args)))
    funcretstr = mktype(ret, vectype)
    print "%s __vml_%s(%s)" % (funcretstr, name, funcargstr)
    print "{"
    for (n, arg, vmlarg) in zip(range(0, 100), args, vmlargs):
        print ("  %s y%d = bitcast<%s,%s>(x%d);" %
               (mkvmltype(vmlarg, vmltype), n,
                mktype(arg, vectype), mkvmltype(vmlarg, vmltype), n))
    callargstr = ", ".join(map(lambda (n, arg): "y%d" % n,
                               zip(range(0, 100), args)))
    callretstr = mkvmltype(vmlret, vmltype)
    print "  %s r = vecmathlib::%s(%s);" % (callretstr, name, callargstr)
    # We may need to convert from the VML type to the OpenCL type
    # before bitcasting. This may be a real conversion, e.g. bool to
    # int. This may also involve a change in size (e.g. long to int),
    # but only if the type is scalar. These conversions are applied
    # before bitcasting.
    # convfunc: conversion function to call
    # convtype: result type of conversion, also input to bitcast
    # bitcasttype: output of bitcast; may differ from function result
    #              if a size change is needed
    if vmlret==ret:
        convfunc    = ""
        convtype    = callretstr
        bitcasttype = funcretstr
    else:
        if vmlret==VI and ret in (VJ,VK):
            convfunc    = ""
            convtype    = callretstr
        elif vmlret==VB and ret in (VJ,VK):
            convfunc = "vecmathlib::convert_int"
            convtype = vmlinttype
        else:
            raise "missing"
        if ret in (VJ,VK):
            bitcasttype = mktype(VI, vectype)
        else:
            raise "missing"
    print "  return bitcast<%s,%s>(%s(r));" % (convtype, bitcasttype, convfunc)
    print "}"

def output_func_libm(func, vectype):
    (name, args, ret, vmlargs, vmlret) = func
    print "// Implement %s by calling libm" % name
    (basetype, size) = re.match("([A-Za-z]+)([0-9]*)", vectype).groups()
    size = 1 if size=="" else int(size)
    othertype = "vecmathlib::realpseudovec<%s,%d>" % (basetype, size)
    otherinttype = "%s::intvec_t" % othertype
    funcargstr = ", ".join(map(lambda (n, arg):
                                   "%s x%d" % (mktype(arg, vectype), n),
                               zip(range(0, 100), args)))
    print "%s __vml_%s(%s)" % (vectype, name, funcargstr)
    print "{"
    for (n, arg) in zip(range(0, 100), args):
        print ("  %s y%d = x%d;" % (othertype, n, n))
    callargstr = ", ".join(map(lambda (n, arg): "y%d" % n,
                               zip(range(0, 100), args)))
    callretstr = othertype if ret==VF else otherinttype
    print "  %s r = vecmathlib::%s(%s);" % (callretstr, name, callargstr)
    print "  return r[0];"
    print "}"

def output_func_upcast(func, vectype):
    (name, args, ret, vmlargs, vmlret) = func
    print "// Implement %s by using a larger vector size" % name
    (basetype, size) = re.match("([A-Za-z]+)([0-9]*)", vectype).groups()
    size = 1 if size=="" else int(size)
    size2 = 4 if size==3 else size*2 # next power of 2
    size2 = "" if size2==1 else str(size2)
    othertype = "%s%s" % (basetype, size2)
    declargstr = ", ".join(map(lambda (n, arg): "%s" % mktype(arg, othertype),
                               zip(range(0, 100), args)))
    print "%s __vml_%s(%s);" % (mktype(ret, othertype), name, declargstr)
    funcargstr = ", ".join(map(lambda (n, arg):
                                   "%s x%d" % (mktype(arg, vectype), n),
                               zip(range(0, 100), args)))
    print "%s __vml_%s(%s)" % (mktype(ret, vectype), name, funcargstr)
    print "{"
    for (n, arg) in zip(range(0, 100), args):
        print ("  %s y%d = bitcast<%s,%s>(x%d);" %
               (mktype(arg, othertype), n,
                mktype(arg, vectype), mktype(arg, othertype), n))
    callargstr = ", ".join(map(lambda (n, arg): "y%d" % n,
                               zip(range(0, 100), args)))
    print "  %s r = __vml_%s(%s);" % (mktype(ret, othertype), name, callargstr)
    print ("  return bitcast<%s,%s>(r);" %
           (mktype(ret, othertype), mktype(ret, vectype)))
    print "}"

def output_func_split(func, vectype):
    (name, args, ret, vmlargs, vmlret) = func
    print "// Implement %s by splitting into a smaller vector size" % name
    (basetype, size) = re.match("([A-Za-z]+)([0-9]*)", vectype).groups()
    size = 1 if size=="" else int(size)
    size2 = (size+1) / 2        # divide by 2, rounding up
    size2 = "" if size2==1 else str(size2)
    othertype = "%s%s" % (basetype, size2)
    declargstr = ", ".join(map(lambda (n, arg): "%s" % mktype(arg, othertype),
                               zip(range(0, 100), args)))
    print "%s __vml_%s(%s);" % (mktype(ret, othertype), name, declargstr)
    funcargstr = ", ".join(map(lambda (n, arg):
                                   "%s x%d" % (mktype(arg, vectype), n),
                               zip(range(0, 100), args)))
    print "%s __vml_%s(%s)" % (mktype(ret, vectype), name, funcargstr)
    print "{"
    print "  struct pair { %s lo, hi; };" % othertype
    for (n, arg) in zip(range(0, 100), args):
        print ("  %s y%d = bitcast<%s,%s>(x%d);" %
               (mktype(arg, othertype), n,
                mktype(arg, vectype), mktype(arg, othertype), n))
    callargstr = ", ".join(map(lambda (n, arg): "y%d" % n,
                               zip(range(0, 100), args)))
    print "  %s r = __vml_%s(%s);" % (mktype(ret, othertype), name, callargstr)
    print ("  return bitcast<%s,%s>(r);" %
           (mktype(ret, othertype), mktype(ret, vectype)))
    print "}"



def output_func(func):
    (name, args, ret, vmlargs, vmlret) = func
    print ""
    print "// %s: %s -> %s" % (name, args, ret)
    for basetype in ["float", "double"]:
        if basetype=="double":
            print "#ifdef cl_khr_fp64"
        for size in [1, 2, 3, 4, 8, 16]:
            # Ignore this prototype for size==1 if there are any
            # scalar arguments; this prevents duplicate definitions
            if size==1 and any(map(lambda arg: arg==SI or arg==SJ, args)):
                continue
            sizename = '' if size==1 else str(size)
            type = basetype + sizename
            # always use vecmathlib if available
            print ("#if defined VECMATHLIB_HAVE_VEC_%s_%d" %
                   (basetype.upper(), size))
            output_func_vml(func, type)
            if size==1:
                # a scalar type: use libm
                print "#else"
                output_func_libm(func, type)
            else:
                # a vector type: try upcasting to next power of 2
                size2 = 4 if size==3 else size*2
                print ("#elif defined VECMATHLIB_HAVE_VEC_%s_%d" %
                       (basetype.upper(), size2))
                output_func_upcast(func, type)
                # a vector type: split into smaller vector type
                print "#else"
                output_func_split(func, type)
            print "#endif"
        if basetype=="double":
            print "#endif // #ifdef cl_khr_fp64"



output_prelude()
map(output_func, funcs)
