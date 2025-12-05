#pragma once

#include <random>

class UniformRNG {
public:
  UniformRNG(const double lower = 0, const double upper = 1.0,
             const unsigned int seed = std::random_device{}());

  double generateRN();

private:
  std::mt19937 _generator; // choose rng to be universally applicable
  std::uniform_real_distribution<double> _distribution;
};