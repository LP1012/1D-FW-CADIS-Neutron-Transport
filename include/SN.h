#pragma once

#include "Cell.h"
#include "SNCell.h"
#include <vector>

template <std::size_t N>
class SN
{
public:
  explicit SN(const std::vector<Cell> & cells);

  void run();
  std::vector<double> getScalarFlux();

protected:
  const unsigned int _num_cells;
  std::vector<double> _mus;
  std::vector<SNCell<N>> _sn_cells;
  double _k;

  void populateSNCells(const std::vector<Cell> & cells);

  void normalizeSources();
  void sweepRight(const unsigned int mu_index);
  void sweepLeft(const unsigned int mu_index);
  void computeScalarFluxAll();

  double computeAngularFlux(const SNCell<N> & cell, const double cell_flux, const double mu);
  double computeCellFlux(const double cell_centered_flux, const double known_flux);

  void updateK();
  double integrateFissionSource(const std::vector<SNCell<N>> & cells);

  void updateSource();
};