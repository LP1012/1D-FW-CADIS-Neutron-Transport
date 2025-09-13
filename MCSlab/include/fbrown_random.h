//
// The following code has been adapted from mcnp_random.h
// written by Forrest Brown. This routine adaps the original
// C version to C++.
//

#pragma once

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// define aliases for common types
using LONG = long long;
using ULONG = unsigned long long;
using REAL = double;

class RNG {
public:
  RNG(ULONG seed = 1ULL); // class constructor
  REAL randomNumber();    // random number generator

private:
  ULONG RN_SEED;

  const int RN_INDEX = 1;
  const ULONG RN_MULT = 9219741426499971445ULL;
  const ULONG RN_ADD = 1ULL;
  const int RN_BITS = 63;
  const int RN_SHIFT = 10;
  const ULONG RN_STRIDE = 152917ULL;
  const ULONG RN_SEED0 = 1ULL;
  const ULONG RN_MOD = 1ULL << 63;
  const ULONG RN_MASK = (1ULL << 63) - 1ULL;
  const ULONG RN_PERIOD = 1ULL << 63;
  const double RN_NORM = 1.0 / (double)(1ULL << 53);
};

//----------------------------------------------------------------------
// reference data:  seeds for case of init.seed = 2,
//                  seed numbers for index 1-5, 123456-123460
//----------------------------------------------------------------------
inline constexpr ULONG RN_CHECK[10] = {
    // ***** 2 *****
    9219741426499971446ULL, 666764808255707375ULL,  4935109208453540924ULL,
    7076815037777023853ULL, 5594070487082964434ULL, 7069484152921594561ULL,
    8424485724631982902ULL, 19322398608391599ULL,   8639759691969673212ULL,
    8181315819375227437ULL};
