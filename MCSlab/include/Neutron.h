#pragma once

#include "Point.h"
#include "RNG.h"
#include "Region.h"
#include <vector>

class Neutron {
public:
  Neutron(double position);

  // move neutron to new position
  void movePositionAndRegion(const double new_position,
                             const std::vector<Region> &regions);

  /// @brief Set start position of neutron randomly in fissionable regionss
  /// @param fissionable_regions
  void setRandomStartPosition(const std::vector<Region> &fissionable_regions);

  double distanceToCollision();
  double distanceToEdge();

  // define getter functions
  double pos() const { return _pos; }
  double ang() const { return _ang; }
  double mu() const { return _mu; }

private:
  double _pos;     // x-position
  double _ang;     // angle (theta) in radians
  double _mu;      // cosine of angle
  Region *_region; // region currently located in
  bool _is_alive;  // is neutron still being tracked?

  // initialize RNG
  UniformRNG _rng;

  /// @brief Set cosine of angle randomly assuming isotropic distribution
  void randomIsoAngle();

  void setRegion(const std::vector<Region> &regions);

  void kill();
};
