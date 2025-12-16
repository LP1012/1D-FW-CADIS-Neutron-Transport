#pragma once

#include "Point.h"
#include "RNG.h"
#include "Cell.h"
#include <optional>
#include <vector>
#include <iostream>
#include <deque>

class Neutron
{
  friend class NeutronTest;

public:
  Neutron(double position, double mu, Cell * cell, std::optional<unsigned int> seed = std::nullopt);

  // move neutron to new position
  void movePositionAndCell(const double new_position, std::vector<Cell> & cells);

  void movePositionWithinCell(const double new_position);

  /// @brief Set start position of neutron randomly in fissionable regionss
  /// @param fissionable_regions
  void setRandomStartPosition(const std::vector<Cell> & fissionable_cells,
                              std::vector<Cell> & cells);

  /// @brief Calculates distance to collision
  double distanceToCollision();

  /// @brief Calculates distance to edge along current direction
  double distanceToEdge();

  /// @brief Terminate particle
  void kill();

  bool weightIsOkay();
  void roulette();
  void split(std::deque<Neutron> & split_bank);

  /// @brief Set cosine of angle randomly assuming isotropic distribution
  static double randomIsoAngle(UniformRNG rng);

  /// @brief Change weight oe neutron
  /// @param weight
  void changeWeight(const double new_weight);

  // define getter functions
  double pos() const { return _pos; }
  double mu() const { return _mu; }
  bool isAlive() const { return _is_alive; }
  const double weight() const { return _weight; }
  const Cell & cell() const { return *_cell; }

private:
  double _pos;    // x-position
  double _mu;     // cosine of angle
  Cell * _cell;   // region currently located in
  bool _is_alive; // is neutron still being tracked?
  double _weight; // weight of neutron

  // initialize RNG
  UniformRNG _rng;

  /// @brief Set region based on current location. Will not work on edge.
  /// @param regions
  void setCell(std::vector<Cell> & cells);
  void checkNeutron();
};
