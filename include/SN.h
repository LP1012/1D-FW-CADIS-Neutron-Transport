#pragma once

#include "Cell.h"
#include "SNCell.h"
#include <vector>
#include "CellQuadrature.h"

class SN
{
public:
  SN(const std::vector<Cell> & cells, const unsigned int GQ_order);

  void run();
  std::vector<double> getScalarFlux();
  std::vector<double> getSources();
  double L2Norm(const std::vector<double> & vector);

protected:
  const unsigned int _num_cells;
  const unsigned int _gq_order;
  discreteQuadrature::GaussLegendreRule _gauss_legendre_rule;

  bool _is_converged;
  std::vector<double> _mus;
  std::vector<SNCell> _sn_cells;
  double _k;

  void populateSNCells(const std::vector<Cell> & cells);

  void normalizeSources();
  void sweepRight(const unsigned int mu_index);
  void sweepLeft(const unsigned int mu_index);
  void computeScalarFluxAll();

  double computeAngularFlux(const SNCell & cell, const double cell_flux, const double mu);
  double computeCellFlux(const double cell_centered_flux, const double known_flux);

  void updateK();
  double integrateFissionSource(const std::vector<SNCell> & cells);

  void updateSource();
  bool isConverged(const std::vector<double> & old_flux,
                   const std::vector<double> & new_flux,
                   const double old_k,
                   const double new_k);

  void exportToCsv();
};