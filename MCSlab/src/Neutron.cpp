#include "Neutron.h"
#include "Point.h"
#include <cmath>
#include <random>

Neutron::Neutron(Point position) : _pos(position), _rng() {
  Neutron::randomIsoAngle();
};

void Neutron::move(const double distance) {
  _pos.moveX(distance * _mu); // CHECK MATH
};

void Neutron::randomIsoAngle() {
  _mu = 2 * _rng.generateRN() - 1; // sample between -1 and 1
  _ang = std::acos(_mu);
};
