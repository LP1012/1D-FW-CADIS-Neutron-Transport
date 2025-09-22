#pragma once

#include "Point.h"
#include "RNG.h"
#include "Region.h"
#include <vector>

class Neutron {
public:
  Neutron(double position, const std::vector<Region> &regions);

  // move neutron to new position
  void movePosition(const double new_position);

  /// @brief Set start position of neutron randomly in fissionable regionss
  /// @param fissionable_regions
  void setRandomStartPosition(const std::vector<Region> &fissionable_regions);

  double distanceToCollision(const double mfp);
  double distanceToEdge();

  // define getter functions
  double pos() const { return _pos; }
  double ang() const { return _ang; }
  double mu() const { return _mu; }
  unsigned int regionID() const { return _region_id; }

private:
  double _pos;                        // x-position
  double _ang;                        // angle (theta) in radians
  double _mu;                         // cosine of angle
  unsigned int _region_id;            // ID of region currently located in
  Region *_region;                    // region currently located in
  bool _is_alive;                     // is neutron still being tracked?
  const std::vector<Region> _regions; // all regions available

  // initialize RNG
  UniformRNG _rng;

  /// @brief Set cosine of angle randomly assuming isotropic distribution
  void randomIsoAngle();
  /// @brief Set region ID based on current location and regions available
  /// @param regions
  void setRegionID();

  void setRegion();

  void kill();
};
