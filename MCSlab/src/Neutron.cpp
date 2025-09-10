#include "Neutron.h"
#include "mcnp_random.h"

#include <cmath>

Neutron::Neutron(Point position) : _pos(position)
{
  Neutron::randomIsoAngle();
};

void
Neutron::move(const double distance)
{
  _pos.x = distance * _mu; // CHECK
}


void
Neutron::randomIsoAngle()
{
  _mu = 2*rang()-1;// sample between -1 and 1
  _ang = std::acos(_mu);
}
