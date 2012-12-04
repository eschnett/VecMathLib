// -*-C++-*-

#include "vecmathlib.h"

#include <iostream>

using namespace std;

float arg_f1(float x) { return -vecmathlib::realvec<float,1>(x)[0]; }
float arg_f4(__m128 x) { return -vecmathlib::realvec<float,4>(x)[0]; }
float arg_f8(__m256 x) { return -vecmathlib::realvec<float,8>(x)[0]; }

int size_f1(float x) { return sizeof(vecmathlib::realvec<float,1>); }
int size_f4(__m128 x) { return sizeof(vecmathlib::realvec<float,4>); }
int size_f8(__m256 x) { return sizeof(vecmathlib::realvec<float,8>); }

float sum_f1(vecmathlib::realvec<float,1> x) { return sum(x); }
float sum_f4(vecmathlib::realvec<float,4> x) { return sum(x); }
float sum_f8(vecmathlib::realvec<float,8> x) { return sum(x); }
double sum_d1(vecmathlib::realvec<double,1> x) { return sum(x); }
double sum_d2(vecmathlib::realvec<double,2> x) { return sum(x); }
double sum_d4(vecmathlib::realvec<double,4> x) { return sum(x); }

bool all_f1(vecmathlib::boolvec<float,1> x) { return all(x); }
bool all_f4(vecmathlib::boolvec<float,4> x) { return all(x); }
bool all_f8(vecmathlib::boolvec<float,8> x) { return all(x); }
bool all_d1(vecmathlib::boolvec<double,1> x) { return all(x); }
bool all_d2(vecmathlib::boolvec<double,2> x) { return all(x); }
bool all_d4(vecmathlib::boolvec<double,4> x) { return all(x); }

float elt0_f1(vecmathlib::realvec<float,1> x) { return x[0]; }
float elt0_f4(vecmathlib::realvec<float,4> x) { return x[0]; }
float elt1_f4(vecmathlib::realvec<float,4> x) { return x[1]; }
float elt2_f4(vecmathlib::realvec<float,4> x) { return x[2]; }
float elt3_f4(vecmathlib::realvec<float,4> x) { return x[3]; }
float elt0_f8(vecmathlib::realvec<float,8> x) { return x[0]; }
float elt1_f8(vecmathlib::realvec<float,8> x) { return x[1]; }
float elt2_f8(vecmathlib::realvec<float,8> x) { return x[2]; }
float elt3_f8(vecmathlib::realvec<float,8> x) { return x[3]; }
float elt4_f8(vecmathlib::realvec<float,8> x) { return x[4]; }
float elt5_f8(vecmathlib::realvec<float,8> x) { return x[5]; }
float elt6_f8(vecmathlib::realvec<float,8> x) { return x[6]; }
float elt7_f8(vecmathlib::realvec<float,8> x) { return x[7]; }
double elt0_d1(vecmathlib::realvec<double,1> x) { return x[0]; }
double elt0_d2(vecmathlib::realvec<double,2> x) { return x[0]; }
double elt1_d2(vecmathlib::realvec<double,2> x) { return x[1]; }
double elt0_d4(vecmathlib::realvec<double,4> x) { return x[0]; }
double elt1_d4(vecmathlib::realvec<double,4> x) { return x[1]; }
double elt2_d4(vecmathlib::realvec<double,4> x) { return x[2]; }
double elt3_d4(vecmathlib::realvec<double,4> x) { return x[3]; }

int main(int argc, char** argv)
{
  using namespace vecmathlib;
  typedef realvec<float,1> realvec_t;
  // typedef realvec<double,4> realvec_t;
  typedef realvec_t::boolvec_t boolvec_t;
  typedef realvec_t::intvec_t intvec_t;
  
  realvec_t x = 1.0;
  realvec_t y = x + realvec_t(1.0);
  y = sqrt(y);
  realvec_t z = log(y);
  boolvec_t b = x < y;
  intvec_t i = convert_int(y);
  
  cout << "x=" << x << "\n";
  cout << "y=" << y << "\n";
  cout << "z=" << z << "\n";
  cout << "b=" << b << "\n";
  cout << "i=" << i << "\n";
  
  return 0;
}
