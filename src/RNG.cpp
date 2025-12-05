#include "RNG.h"

UniformRNG::UniformRNG(const double lower, const double upper, const unsigned int seed)
  : _generator(seed), _distribution(lower, upper){};

double
UniformRNG::generateRN()
{
  return _distribution(_generator);
};