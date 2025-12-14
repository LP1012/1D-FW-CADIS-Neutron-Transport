#pragma once

#include "Cell.h"
#include <vector>

class SN
{
public:
  SN(const std::vector<Cell> & cells);

  std::vector<double> forwardFlux();

protected:
  std::vector<double> _source_vector;
  double _k;

  double computeRightFlux(const int cell_index);
  double computeLeftFlux(const int cell_index);
};