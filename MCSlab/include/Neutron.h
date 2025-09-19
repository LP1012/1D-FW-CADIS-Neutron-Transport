#pragma once

#include "Point.h"
#include "RNG.h"
#include "Region.h"
#include <vector>

class Neutron {
public:
  Neutron(double position, std::vector<Region> &regions);

  // move neutron to new position
  void movePosition(const double new_position);

  // define getter functions
  double pos() const { return _pos; }
  double ang() const { return _ang; }
  double mu() const { return _mu; }

private:
  double _pos;
  double _ang;
  double _mu;
  unsigned int _region_id;
  const std::vector<Region> &_regions;

  // initialize RNG
  UniformRNG _rng;

  void randomIsoAngle();
  void setRegionID(const std::vector<Region> &regions);
};
