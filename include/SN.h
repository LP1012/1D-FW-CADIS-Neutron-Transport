#pragma once

#include "Cell.h"
#include "SNCell.h"
#include <vector>

template <std::size_t N>
class SN
{
public:
  explicit SN(const std::vector<Cell> & cells);

  std::vector<double> forwardFlux();

protected:
  const unsigned int _num_cells;
  std::vector<double> _mus;
  std::vector<SNCell> _sn_cells;
  std::vector<double> _source_vector;
  double _k;

  void populateSNCells(const std::vector<Cell> & cells);

  void normalizeSources();
  void sweepRight(const unsigned int mu_index);
  void sweepLeft(const unsigned int mu_index);

  double computeAngularFlux(const SNCell & cell, const double cell_flux, const double mu);
  double computeCellFlux(const double cell_centered_flux, const double known_flux);
  double computeRightFlux(const int cell_index);
  double computeLeftFlux(const int cell_index);
};