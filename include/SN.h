#pragma once

#include "Cell.h"
#include "SNCell.h"
#include <vector>

class SN
{
public:
  SN(const std::vector<Cell> & cells, const unsigned int GQ_order);

  std::vector<double> forwardFlux();

protected:
  const unsigned int _gq_order;
  std::vector<SNCell> _sn_cells;
  std::vector<double> _source_vector;
  double _k;

  void populateSNCells(const std::vector<Cell> & cells);

  double computeRightFlux(const int cell_index);
  double computeLeftFlux(const int cell_index);
};