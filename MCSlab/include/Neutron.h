#pragma once

#include "Point.h"
#include "RNG.h"
#include "Region.h"
#include <optional>
#include <vector>

class Neutron {
public:
  Neutron(double position, const std::vector<Region> &regions,
          std::optional<unsigned int> seed = std::nullopt);

  // move neutron to new position
  void movePositionAndRegion(const double new_position,
                             const std::vector<Region> &regions);

  /// @brief Set start position of neutron randomly in fissionable regionss
  /// @param fissionable_regions
  void setRandomStartPosition(const std::vector<Region> &fissionable_regions);

  /// @brief Calculates distance to collision
  double distanceToCollision();

  /// @brief Calculates distance to edge along current direction
  double distanceToEdge();

  /// @brief Set position in region on boundary (region depends on direction)
  /// @param new_position
  /// @param new_region
  void setPositionOnBoundary(const double new_position, Region &new_region);

  /// @brief Terminate particle
  void kill();

  /// @brief Set cosine of angle randomly assuming isotropic distribution
  void randomIsoAngle();

  // define getter functions
  double pos() const { return _pos; }
  double mu() const { return _mu; }
  bool isAlive() const { return _is_alive; }
  Region region() const { return *_region; }

private:
  double _pos;     // x-position
  double _mu;      // cosine of angle
  Region *_region; // region currently located in
  bool _is_alive;  // is neutron still being tracked?

  // initialize RNG
  UniformRNG _rng;

  /// @brief Set region based on current location. Will not work on edge.
  /// @param regions
  void setRegion(const std::vector<Region> &regions);
};
