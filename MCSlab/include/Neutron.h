#pragma once

#include "Point.h"
#include "RNG.h"

class Neutron {
public:
  Neutron(Point position);

  // move neutron to new position
  void move(const double distance);

  // define getter functions
  Point pos() const { return _pos; }
  double ang() const { return _ang; }
  double mu() const { return _mu; }

private:
  Point _pos;
  double _ang;
  double _mu;

  // initialize RNG
  UniformRNG _rng;

  void randomIsoAngle();
};
