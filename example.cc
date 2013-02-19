// -*-C++-*-

#include "vecmathlib.h"

#include <iostream>

using namespace std;



int main(int argc, char** argv)
{
  using namespace vecmathlib;
  // typedef realvec<float,1> realvec_t;
  typedef realvec<double,2> realvec_t;
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
