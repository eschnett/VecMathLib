// -*-C++-*-

#ifndef FLOATTYPESS_H
#define FLOATTYPESS_H

#include <cstdint>

namespace vecmathlib {
  
  struct fp8 {
    // 1 bit sign, 4 bits exponent, 3 bits mantissa
    std::uint8_t val;
  };
  
  struct fp16 {
    // 1 bit sign, 5 bits exponent, 10 bits mantissa
    std::uint16_t val;
  };
  
} // namespace vecmathlib

#endif  // #ifndef FLOATTYPES_H
