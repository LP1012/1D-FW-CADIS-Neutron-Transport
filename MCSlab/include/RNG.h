#pragma once

#include <random>

class UniformRNG {
public:
  UniformRNG(const double lower = 0, const double upper = 1.0);
  double generateRN();

private:
  std::default_random_engine _generator;
  std::uniform_real_distribution<double> _distribution;
};