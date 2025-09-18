#include "RNG.h"

UniformRNG::UniformRNG(const double lower, const double upper)
    : _distribution(lower, upper){};

double UniformRNG::generateRN() { return _distribution(_generator); };