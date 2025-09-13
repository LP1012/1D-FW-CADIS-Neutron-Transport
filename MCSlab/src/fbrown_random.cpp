#include "fbrown_random.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

RNG::RNG(ULONG seed) { RN_SEED = seed; };

REAL RNG::randomNumber() {
  ULONG I53;
  RN_SEED = (RN_MULT * RN_SEED + RN_ADD) & RN_MASK;
  I53 = RN_SEED >> RN_SHIFT;

  if (!I53)
    I53++;

  return (REAL)(I53 * RN_NORM);
};
