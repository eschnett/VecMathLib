#! /usr/bin/env python

import re



# Types:
SI = "SI"                       # int/long
SJ = "SJ"                       # int (even for double)
SF = "SF"                       # float/double
VB = "VB"                       # boolN
VF = "VF"                       # floatN/doubleN
VI = "VI"                       # intN/longN
VJ = "VJ"                       # intN/longN (except int1 for double1)
VK = "VK"                       # intN (even for doubleN)
VU = "VU"                       # uintN/ulongN

# Each function is described by a tuple with the following entries:
#    1. name
#    2. external argument types (see above)
#    3. external return type
#    4. vecmathlib argument types (see above)
#    5. vecmathlib return type
# This allows generating externally visible functions with different
# signatures, e.g. to support OpenCL.
vmlfuncs = [
    ("acos"     , [VF        ], VF, [VF        ], VF), # 6.12.2
    ("acosh"    , [VF        ], VF, [VF        ], VF), # 6.12.2
    ("asin"     , [VF        ], VF, [VF        ], VF), # 6.12.2
    ("asinh"    , [VF        ], VF, [VF        ], VF), # 6.12.2
    ("atan"     , [VF        ], VF, [VF        ], VF), # 6.12.2
    ("atanh"    , [VF        ], VF, [VF        ], VF), # 6.12.2
    ("cbrt"     , [VF        ], VF, [VF        ], VF), # 6.12.2
    ("ceil"     , [VF        ], VF, [VF        ], VF), # 6.12.2
    ("copysign" , [VF, VF    ], VF, [VF, VF    ], VF), # 6.12.2
    ("cos"      , [VF        ], VF, [VF        ], VF), # 6.12.2
    ("cosh"     , [VF        ], VF, [VF        ], VF), # 6.12.2
    ("exp"      , [VF        ], VF, [VF        ], VF), # 6.12.2
    ("exp2"     , [VF        ], VF, [VF        ], VF), # 6.12.2
    ("exp10"    , [VF        ], VF, [VF        ], VF), # 6.12.2
    ("expm1"    , [VF        ], VF, [VF        ], VF), # 6.12.2
    ("fabs"     , [VF        ], VF, [VF        ], VF), # 6.12.2
    ("fdim"     , [VF, VF    ], VF, [VF, VF    ], VF), # 6.12.2
    ("floor"    , [VF        ], VF, [VF        ], VF), # 6.12.2
    ("fma"      , [VF, VF, VF], VF, [VF, VF, VF], VF), # 6.12.2
    ("fmax"     , [VF, VF    ], VF, [VF, VF    ], VF), # 6.12.2
    ("fmin"     , [VF, VF    ], VF, [VF, VF    ], VF), # 6.12.2
    ("fmod"     , [VF, VF    ], VF, [VF, VF    ], VF), # 6.12.2
    ("hypot"    , [VF, VF    ], VF, [VF, VF    ], VF), # 6.12.2
    ("ilogb"    , [VF        ], VJ, [VF        ], VI), # 6.12.2 (but should return VK)
    ("ldexp"    , [VF, VJ    ], VF, [VF, VI    ], VF), # 6.12.2 (but should take VK)
    ("ldexp"    , [VF, SJ    ], VF, [VF, SI    ], VF), # 6.12.2
    ("log"      , [VF        ], VF, [VF        ], VF), # 6.12.2
    ("log2"     , [VF        ], VF, [VF        ], VF), # 6.12.2
    ("log10"    , [VF        ], VF, [VF        ], VF), # 6.12.2
    ("log1p"    , [VF        ], VF, [VF        ], VF), # 6.12.2
    ("pow"      , [VF, VF    ], VF, [VF, VF    ], VF), # 6.12.2
    ("remainder", [VF, VF    ], VF, [VF, VF    ], VF), # 6.12.2
    ("round"    , [VF        ], VF, [VF        ], VF), # 6.12.2
    ("rsqrt"    , [VF        ], VF, [VF        ], VF), # 6.12.2
    ("sin"      , [VF        ], VF, [VF        ], VF), # 6.12.2
    ("sinh"     , [VF        ], VF, [VF        ], VF), # 6.12.2
    ("sqrt"     , [VF        ], VF, [VF        ], VF), # 6.12.2
    ("tan"      , [VF        ], VF, [VF        ], VF), # 6.12.2
    ("tanh"     , [VF        ], VF, [VF        ], VF), # 6.12.2
    ("trunc"    , [VF        ], VF, [VF        ], VF), # 6.12.2
    
    ("isfinite" , [VF        ], VJ, [VF        ], VB), # 6.12.6
    ("isinf"    , [VF        ], VJ, [VF        ], VB), # 6.12.6
    ("isnan"    , [VF        ], VJ, [VF        ], VB), # 6.12.6
    ("isnormal" , [VF        ], VJ, [VF        ], VB), # 6.12.6
    ("signbit"  , [VF        ], VJ, [VF        ], VB), # 6.12.6
    ]



directfuncs = [
    ("acospi"       , [VF        ], VF, "acos(x0)/(scalar_t)M_PI"),     # 6.12.2
    ("asinpi"       , [VF        ], VF, "asin(x0)/(scalar_t)M_PI"),     # 6.12.2
    ("atanpi"       , [VF        ], VF, "atan(x0)/(scalar_t)M_PI"),     # 6.12.2
    ("atan2pi"      , [VF, VF    ], VF, "atan2(x0,x1)/(scalar_t)M_PI"), # 6.12.2
    ("cospi"        , [VF        ], VF, "cos((scalar_t)M_PI*x0)"),      # 6.12.2
    ("fmax"         , [VF, SF    ], VF, "fmax(x0,(vector_t)x1)"),       # 6.12.2
    ("fmin"         , [VF, SF    ], VF, "fmin(x0,(vector_t)x1)"),       # 6.12.2
    ("mad"          , [VF, VF, VF], VF, "fma(x0,x1,x2)"),               # 6.12.2
    ("maxmag"       , [VF, VF    ], VF, "fabs(x0)>fabs(x1) ? x0 : fabs(x1)>fabs(x0) ? x1 : fmax(x0,x1)"), # 6.12.2
    ("minmag"       , [VF, VF    ], VF, "fabs(x0)<fabs(x1) ? x0 : fabs(x1)<fabs(x0) ? x1 : fmin(x0,x1)"), # 6.12.2
    ("nan"          , [VU        ], VF, "(scalar_t)0.0/(scalar_t)0.0"), # 6.12.2
    ("pown"         , [VF, VK    ], VF, "pow(x0,convert_vector_t(x1))"), # 6.12.2
    ("powr"         , [VF, VF    ], VF, "pow(x0,x1)"), # 6.12.2
    ("rint"         , [VF        ], VF, "round(x0)"),  # 6.12.2
    ("rootn"        , [VF, VK    ], VF, "pow(x0,(scalar_t)1.0/convert_vector_t(x1))"), # 6.12.2
    ("sinpi"        , [VF        ], VF, "sin((scalar_t)M_PI*x0)"), # 6.12.2
    ("tanpi"        , [VF        ], VF, "tan((scalar_t)M_PI*x0)"), # 6.12.2
    
    ("half_cos"     , [VF        ], VF, "cos(x0)"),          # 6.12.2
    ("half_divide"  , [VF, VF    ], VF, "x0/x1"),            # 6.12.2
    ("half_exp"     , [VF        ], VF, "exp(x0)"),          # 6.12.2
    ("half_exp2"    , [VF        ], VF, "exp2(x0)"),         # 6.12.2
    ("half_exp10"   , [VF        ], VF, "exp10(x0)"),        # 6.12.2
    ("half_log"     , [VF        ], VF, "log(x0)"),          # 6.12.2
    ("half_log2"    , [VF        ], VF, "log2(x0)"),         # 6.12.2
    ("half_log10"   , [VF        ], VF, "log10(x0)"),        # 6.12.2
    ("half_powr"    , [VF, VF    ], VF, "powr(x0,x1)"),      # 6.12.2
    ("half_recip"   , [VF        ], VF, "(scalar_t)1.0/x0"), # 6.12.2
    ("half_rsqrt"   , [VF        ], VF, "rsqrt(x0)"),        # 6.12.2
    ("half_sin"     , [VF        ], VF, "sin(x0)"),          # 6.12.2
    ("half_sqrt"    , [VF        ], VF, "sqrt(x0)"),         # 6.12.2
    ("half_tan"     , [VF        ], VF, "tan(x0)"),          # 6.12.2
    
    ("native_cos"   , [VF        ], VF, "cos(x0)"),          # 6.12.2
    ("native_divide", [VF, VF    ], VF, "x0/x1"),            # 6.12.2
    ("native_exp"   , [VF        ], VF, "exp(x0)"),          # 6.12.2
    ("native_exp2"  , [VF        ], VF, "exp2(x0)"),         # 6.12.2
    ("native_exp10" , [VF        ], VF, "exp10(x0)"),        # 6.12.2
    ("native_log"   , [VF        ], VF, "log(x0)"),          # 6.12.2
    ("native_log2"  , [VF        ], VF, "log2(x0)"),         # 6.12.2
    ("native_log10" , [VF        ], VF, "log10(x0)"),        # 6.12.2
    ("native_powr"  , [VF, VF    ], VF, "powr(x0,x1)"),      # 6.12.2
    ("native_recip" , [VF        ], VF, "(scalar_t)1.0/x0"), # 6.12.2
    ("native_rsqrt" , [VF        ], VF, "rsqrt(x0)"),        # 6.12.2
    ("native_sin"   , [VF        ], VF, "sin(x0)"),          # 6.12.2
    ("native_sqrt"  , [VF        ], VF, "sqrt(x0)"),         # 6.12.2
    ("native_tan"   , [VF        ], VF, "tan(x0)"),          # 6.12.2
    
    ("clamp"        , [VF, VF, VF], VF, "fmin(fmax(x0,x1),x2)"), # 6.12.4
    ("clamp"        , [VF, SF, SF], VF, "fmin(fmax(x0,x1),x2)"), # 6.12.4
    ("degrees"      , [VF        ], VF, "(scalar_t)(180.0/M_PI)*x0"), # 6.12.4
    ("max"          , [VF, VF    ], VF, "fmax(x0,x1)"),   # 6.12.4
    ("max"          , [VF, SF    ], VF, "fmax(x0,x1)"),   # 6.12.4
    ("min"          , [VF, VF    ], VF, "fmin(x0,x1)"),   # 6.12.4
    ("min"          , [VF, SF    ], VF, "fmin(x0,x1)"),   # 6.12.4
    ("mix"          , [VF, VF, VF], VF, "x0+(x1-x0)*x2"), # 6.12.4
    ("mix"          , [VF, VF, SF], VF, "x0+(x1-x0)*x2"), # 6.12.4
    ("radians"      , [VF        ], VF, "(scalar_t)(M_PI/180.0)*x0"), # 6.12.4
    ("step"         , [VF, VF    ], VF, "x1<x0 ? (vector_t)(scalar_t)0.0 : (vector_t)(scalar_t)1.0"), # 6.12.4
    ("step"         , [SF, VF    ], VF, "x1<x0 ? (vector_t)(scalar_t)0.0 : (vector_t)(scalar_t)1.0"), # 6.12.4
    ("smoothstep"   , [VF, VF, VF], VF, "({ vector_t t = clamp((x2-x0)/(x1-x0), (scalar_t)0.0, (scalar_t)1.0); t*t*((scalar_t)3.0-(scalar_t)2.0*t); })"), # 6.12.4
    ("smoothstep"   , [SF, SF, VF], VF, "({ vector_t t = clamp((x2-x0)/(x1-x0), (scalar_t)0.0, (scalar_t)1.0); t*t*((scalar_t)3.0-(scalar_t)2.0*t); })"), # 6.12.4
    ("sign"         , [VF        ], VF, "copysign(x0!=(scalar_t)0.0 ? (vector_t)(scalar_t)1.0 : (vector_t)(scalar_t)0.0,x0)"), # 6.12.4
    
    ("isequal"       , [VF, VF    ], VJ, "x0==x1"),                 # 6.12.6
    ("isnotequal"    , [VF, VF    ], VJ, "x0!=x1"),                 # 6.12.6
    ("isgreater"     , [VF, VF    ], VJ, "x0>x1"),                  # 6.12.6
    ("isgreaterequal", [VF, VF    ], VJ, "x0>=x1"),                 # 6.12.6
    ("isless"        , [VF, VF    ], VJ, "x0<x1"),                  # 6.12.6
    ("islessequal"   , [VF, VF    ], VJ, "x0<=x1"),                 # 6.12.6
    ("islessgreater" , [VF, VF    ], VJ, "x0<x1 || x0>x1"),         # 6.12.6
    ("isordered"     , [VF, VF    ], VJ, "!isunordered(x0,x1)"),    # 6.12.6
    ("isunordered"   , [VF, VF    ], VJ, "isnan(x0) || isnan(x1)"), # 6.12.6
]

# Missing functions from 6.12.2: atan2, erfc, erf, fract, frexp,
# lgamma, lgamma_r, logb, modf, nextafter, remquo, sincos, tgamma

# Unchecked: 6.12.3 (integer functions)

# Missing functions from 6.12.6 (relational functions): any, all,
# bitselect, select

# Unchecked: 6.12.7 (vector data load and store functions)

# Unchecked: 6.12.12 (miscellaneous vector functions)



outfile = None
outfile_did_truncate = set()
def out(str): outfile.write("%s\n" % str)
def out_open(name):
    global outfile
    global outfile_did_truncate
    if outfile: raise "file already open"
    is_first_open = name not in outfile_did_truncate
    if is_first_open:
        outfile = open(name, "w")
        outfile.close()
        outfile_did_truncate.add(name)
        print name
    outfile = open(name, "a")
    return is_first_open
def out_close():
    global outfile
    outfile.close()
    outfile = None

declfile = None
def decl(str):
    if str=="" or str.startswith("//") or str.startswith("#"):
        declfile.write("%s\n" % str)
    else:
        declfile.write("__attribute__((__overloadable__)) %s;\n" % str)
def decl_open(name):
    global declfile
    declfile = open(name, "w")
def decl_close():
    global declfile
    declfile.close()
    declfile = None



def mktype(tp, vectype):
    (basetype, size) = re.match("([A-Za-z]+)([0-9]*)", vectype).groups()
    size = 1 if size=="" else int(size)
    if tp==SJ:
        if size==1: return "int"
        return "int" if basetype=="float" else "long"
    if tp==SF:
        return basetype
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
    if tp==VU:
        ibasetype = "uint" if basetype=="float" else "ulong"
        return "%s%s" % (ibasetype, "" if size==1 else str(size))
    if tp==VF:
        return vectype
    raise "unreachable"

def mkvmltype(tp, vectype):
    if tp==SI: return vectype+"::int_t"
    if tp==SF: return vectype+"::real_t"
    if tp==VB: return vectype+"::boolvec_t"
    if tp==VI: return vectype+"::intvec_t"
    if tp==VF: return vectype
    raise "unreachable"



def output_vmlfunc_vml(func, vectype):
    (name, args, ret, vmlargs, vmlret) = func
    out("// Implement %s by calling vecmathlib" % name)
    (basetype, size) = re.match("([A-Za-z]+)([0-9]*)", vectype).groups()
    size = 1 if size=="" else int(size)
    vmltype = "vecmathlib::realvec<%s,%d>" % (basetype, size)
    vmlinttype = "%s::intvec_t" % vmltype
    vmlbooltype = "%s::boolvec_t" % vmltype
    funcargstr = ", ".join(map(lambda (n, arg):
                                   "%s x%d" % (mktype(arg, vectype), n),
                               zip(range(0, 100), args)))
    funcretstr = mktype(ret, vectype)
    decl("%s __vml_%s(%s)" % (funcretstr, name, funcargstr))
    out("%s __vml_%s(%s)" % (funcretstr, name, funcargstr))
    out("{")
    for (n, arg, vmlarg) in zip(range(0, 100), args, vmlargs):
        out("  %s y%d = bitcast<%s,%s>(x%d);" %
            (mkvmltype(vmlarg, vmltype), n,
             mktype(arg, vectype), mkvmltype(vmlarg, vmltype), n))
    callargstr = ", ".join(map(lambda (n, arg): "y%d" % n,
                               zip(range(0, 100), args)))
    callretstr = mkvmltype(vmlret, vmltype)
    out("  %s r = vecmathlib::%s(%s);" % (callretstr, name, callargstr))
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
    out("  return bitcast<%s,%s>(%s(r));" % (convtype, bitcasttype, convfunc))
    out("}")

def output_vmlfunc_libm(func, vectype):
    (name, args, ret, vmlargs, vmlret) = func
    out("// Implement %s by calling libm" % name)
    (basetype, size) = re.match("([A-Za-z]+)([0-9]*)", vectype).groups()
    size = 1 if size=="" else int(size)
    othertype = "vecmathlib::realpseudovec<%s,%d>" % (basetype, size)
    otherinttype = "%s::intvec_t" % othertype
    funcargstr = ", ".join(map(lambda (n, arg):
                                   "%s x%d" % (mktype(arg, vectype), n),
                               zip(range(0, 100), args)))
    decl("%s __vml_%s(%s)" % (vectype, name, funcargstr))
    out("%s __vml_%s(%s)" % (vectype, name, funcargstr))
    out("{")
    for (n, arg) in zip(range(0, 100), args):
        out("  %s y%d = x%d;" % (othertype, n, n))
    callargstr = ", ".join(map(lambda (n, arg): "y%d" % n,
                               zip(range(0, 100), args)))
    callretstr = othertype if ret==VF else otherinttype
    out("  %s r = vecmathlib::%s(%s);" % (callretstr, name, callargstr))
    out("  return r[0];")
    out("}")

def output_vmlfunc_upcast(func, vectype):
    (name, args, ret, vmlargs, vmlret) = func
    out("// Implement %s by using a larger vector size" % name)
    (basetype, size) = re.match("([A-Za-z]+)([0-9]*)", vectype).groups()
    size = 1 if size=="" else int(size)
    size2 = 4 if size==3 else size*2 # next power of 2
    size2 = "" if size2==1 else str(size2)
    othertype = "%s%s" % (basetype, size2)
    declargstr = ", ".join(map(lambda (n, arg): "%s" % mktype(arg, othertype),
                               zip(range(0, 100), args)))
    out("%s __vml_%s(%s);" % (mktype(ret, othertype), name, declargstr))
    funcargstr = ", ".join(map(lambda (n, arg):
                                   "%s x%d" % (mktype(arg, vectype), n),
                               zip(range(0, 100), args)))
    decl("%s __vml_%s(%s)" % (mktype(ret, vectype), name, funcargstr))
    out("%s __vml_%s(%s)" % (mktype(ret, vectype), name, funcargstr))
    out("{")
    for (n, arg) in zip(range(0, 100), args):
        out("  %s y%d = bitcast<%s,%s>(x%d);" %
            (mktype(arg, othertype), n,
             mktype(arg, vectype), mktype(arg, othertype), n))
    callargstr = ", ".join(map(lambda (n, arg): "y%d" % n,
                               zip(range(0, 100), args)))
    out("  %s r = __vml_%s(%s);" % (mktype(ret, othertype), name, callargstr))
    out("  return bitcast<%s,%s>(r);" %
        (mktype(ret, othertype), mktype(ret, vectype)))
    out("}")

def output_vmlfunc_split(func, vectype):
    (name, args, ret, vmlargs, vmlret) = func
    out("// Implement %s by splitting into a smaller vector size" % name)
    (basetype, size) = re.match("([A-Za-z]+)([0-9]*)", vectype).groups()
    size = 1 if size=="" else int(size)
    size2 = (size+1) / 2        # divide by 2, rounding up
    size2 = "" if size2==1 else str(size2)
    othertype = "%s%s" % (basetype, size2)
    declargstr = ", ".join(map(lambda (n, arg): "%s" % mktype(arg, othertype),
                               zip(range(0, 100), args)))
    out("%s __vml_%s(%s);" % (mktype(ret, othertype), name, declargstr))
    funcargstr = ", ".join(map(lambda (n, arg):
                                   "%s x%d" % (mktype(arg, vectype), n),
                               zip(range(0, 100), args)))
    decl("%s __vml_%s(%s)" % (mktype(ret, vectype), name, funcargstr))
    out("%s __vml_%s(%s)" % (mktype(ret, vectype), name, funcargstr))
    out("{")
    out("  struct pair { %s lo, hi; };" % othertype)
    for (n, arg) in zip(range(0, 100), args):
        out("  %s y%d = bitcast<%s,%s>(x%d);" %
            (mktype(arg, othertype), n,
             mktype(arg, vectype), mktype(arg, othertype), n))
    callargstr = ", ".join(map(lambda (n, arg): "y%d" % n,
                               zip(range(0, 100), args)))
    out("  %s r = __vml_%s(%s);" % (mktype(ret, othertype), name, callargstr))
    out("  return bitcast<%s,%s>(r);" %
        (mktype(ret, othertype), mktype(ret, vectype)))
    out("}")



def output_directfunc_direct(func, vectype):
    (name, args, ret, impl) = func
    out("// Implement %s directly" % name)
    (basetype, size) = re.match("([A-Za-z]+)([0-9]*)", vectype).groups()
    size = 1 if size=="" else int(size)
    funcargstr = ", ".join(map(lambda (n, arg):
                                   "%s x%d" % (mktype(arg, vectype), n),
                               zip(range(0, 100), args)))
    funcretstr = mktype(ret, vectype)
    decl("%s __vml_%s(%s)" % (funcretstr, name, funcargstr))
    out("__attribute__((__overloadable__))");
    out("%s __vml_%s(%s)" % (funcretstr, name, funcargstr))
    out("{")
    out("  typedef %s scalar_t;" % basetype)
    out("  typedef %s vector_t;" % vectype)
    out("#define convert_vector_t convert_%s" % vectype)
    out("  return %s;" % impl)
    out("#undef convert_vector_t")
    out("}")



def output_vmlfunc(func):
    (name, args, ret, vmlargs, vmlret) = func
    is_first_open = out_open("%s.cc" % name)
    if is_first_open:
        out("// Note: This file has been automatically generated. Do not modify.")
        out("")
        out("#include \"pocl-compat.h\"")
        out("")
    else:
        out("")
        out("")
        out("")
    decl("")
    decl("// %s: %s -> %s" % (name, args, ret))
    decl("#undef %s" % name)
    decl("#define %s __vml_%s" % (name, name))
    out("// %s: %s -> %s" % (name, args, ret))
    for basetype in ["float", "double"]:
        if basetype=="double":
            out("")
            out("#ifdef cl_khr_fp64")
        for size in [1, 2, 3, 4, 8, 16]:
            # Ignore this prototype for size==1 if there are any
            # scalar arguments; this prevents duplicate definitions
            if size==1 and any(map(lambda arg: arg in (SI, SJ, SF), args)):
                continue
            sizename = '' if size==1 else str(size)
            vectype = basetype + sizename
            # always use vecmathlib if available
            out("")
            out("// %s: VF=%s" % (name, vectype))
            out("#if defined VECMATHLIB_HAVE_VEC_%s_%d" %
                (basetype.upper(), size))
            output_vmlfunc_vml(func, vectype)
            if size==1:
                # a scalar type: use libm
                out("#else")
                output_vmlfunc_libm(func, vectype)
            else:
                # a vector type: try upcasting to next power of 2
                size2 = 4 if size==3 else size*2
                out("#elif defined VECMATHLIB_HAVE_VEC_%s_%d" %
                    (basetype.upper(), size2))
                output_vmlfunc_upcast(func, vectype)
                # a vector type: split into smaller vector type
                out("#else")
                output_vmlfunc_split(func, vectype)
            out("#endif")
        if basetype=="double":
            out("")
            out("#endif // #ifdef cl_khr_fp64")
    out_close()



def output_directfunc(func):
    (name, args, ret, impl) = func
    is_first_open = out_open("%s.cl" % name)
    if is_first_open:
        out("// Note: This file has been automatically generated. Do not modify.")
        out("")
    else:
        out("")
        out("")
        out("")
    decl("")
    decl("// %s: %s -> %s" % (name, args, ret))
    decl("#undef %s" % name)
    decl("#define %s __vml_%s" % (name, name))
    out("// %s: %s -> %s" % (name, args, ret))
    for basetype in ["float", "double"]:
        if ((name.startswith("half_") or name.startswith("native_")) and
            basetype=="double"):
            continue
        if basetype=="double":
            out("")
            out("#ifdef cl_khr_fp64")
        for size in [1, 2, 3, 4, 8, 16]:
            # Ignore this prototype for size==1 if there are any
            # scalar arguments; this prevents duplicate definitions
            if size==1 and any(map(lambda arg: arg in (SI, SJ, SF), args)):
                continue
            sizename = '' if size==1 else str(size)
            vectype = basetype + sizename
            # always use vecmathlib if available
            out("")
            out("// %s: VF=%s" % (name, vectype))
            output_directfunc_direct(func, vectype)
        if basetype=="double":
            out("")
            out("#endif // #ifdef cl_khr_fp64")
    out_close()



decl_open("kernel-vecmathlib.h")
decl("// Note: This file has been automatically generated. Do not modify.")
decl("#ifndef KERNEL_VECMATHLIB_H")
decl("#define KERNEL_VECMATHLIB_H 1")
map(output_vmlfunc, vmlfuncs)
map(output_directfunc, directfuncs)
decl("")
decl("#endif // #ifndef KERNEL_VECMATHLIB_H")
decl_close()
