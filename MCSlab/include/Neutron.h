#pragma once

#include "Point.h"
#include "RNG.h"
#include "Region.h"
#include <optional>
#include <vector>
#include <iostream>

class Neutron
{
  friend class NeutronTest;

public:
  Neutron(double position,
          std::vector<Region> & regions,
          std::optional<unsigned int> seed = std::nullopt);

  // move neutron to new position
  void movePositionAndRegion(const double new_position, std::vector<Region> & regions);

  void movePositionWithinRegion(const double new_position);

  /// @brief Set start position of neutron randomly in fissionable regionss
  /// @param fissionable_regions
  void setRandomStartPosition(const std::vector<Region> & fissionable_regions,
                              std::vector<Region> & regions);

  /// @brief Calculates distance to collision
  double distanceToCollision();

  /// @brief Calculates distance to edge along current direction
  double distanceToEdge();

  /// @brief Terminate particle
  void kill();

  /// @brief Set cosine of angle randomly assuming isotropic distribution
  void randomIsoAngle();

  // define getter functions
  double pos() const { return _pos; }
  double mu() const { return _mu; }
  bool isAlive() const { return _is_alive; }
  // Region region() const { return *_region; }

  const Region & region() const { return *_region; }

private:
  double _pos;      // x-position
  double _mu;       // cosine of angle
  Region * _region; // region currently located in
  bool _is_alive;   // is neutron still being tracked?

  // initialize RNG
  UniformRNG _rng;

  /// @brief Set region based on current location. Will not work on edge.
  /// @param regions
  void setRegion(std::vector<Region> & regions);
};
